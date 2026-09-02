# 7. Optimization Pipeline — End-to-End Reference

This document describes the full training and evaluation pipeline implemented in
`opt_bin_resp/src/`. It maps every conceptual step to the exact file and function
that implements it, and explains the non-obvious numerical tricks at each stage.

**Theory cross-references:** §02 (MWC model), §03 (latent space & affinities), §06 (memory limits).

---

## Overview

```
RunConfig / SingleRunConfig (config.py)
        │
        ▼
LigandEnvironment  ←──────────────────── 1. World Building
        │
        │  sample_batch()
        ▼
(E_open, concs, masks)                  ─── 2. Batch Sampling
        │
        │  BinaryReceptor.forward()
        ▼
activity  (B, R)                        ─── 3. Receptor Physics
        │
        │  CellReadout  (cell mode only — otherwise this stage is absent)
        ▼
activity  (B, C)                        ─── 3b. Cell Readout
        │
        │  DiscreteExactLoss / MI losses
        ▼
scalar loss                             ─── 4. Loss Computation
        │
        │  Adam + temperature schedule
        ▼
trained env + physics                   ─── 5. Training Loop
        │
        │  _eval_stats()
        ▼
metrics dict  →  ExperimentLogger       ─── 6. Evaluation & Logging
```

---

## Stage 1 — World Building (`environment.py::LigandEnvironment.__init__`)

### 1.1 Family Prototype Centers

`n_families` prototype centers are placed inside a uniform N-ball whose radius is
calibrated so that the **expected pairwise distance equals `avg_family_distance`**.
The calibration uses a Monte Carlo estimate of the average distance in a unit N-ball
(`_generate_family_centers`), then scales. Family centers are **frozen buffers** — they
do not receive gradients.

`SymmetricLigandEnvironment` overrides `_generate_family_centers` to use exact geometric
arrangements (polygon in 2D, tetrahedron/octahedron in 3D) for controlled experiments.

### 1.2 Ligand Pool (fixed)

`n_ligands` specific ligands are drawn once and stored as `ligand_latent` buffers:
```
v_ℓ ~ N(v_f, family_spread)   [gaussian mode]
v_ℓ ~ UniformNBall(v_f, family_spread)   [uniform mode]
```
`UniformNBall` uses the **direction + radius trick**: sample a uniform direction on
the sphere, then scale by `r = radius * U^(1/D)` to get uniform density in the ball
(without the concentration-of-measure bias that plagues high-D Gaussian sampling).

### 1.3 Learnable Unit Parameters

The set of learnable parameters depends on the model mode (controlled by `use_interface_model`
in `SingleRunConfig` / `RunConfig`).

**Classic model** (`use_interface_model=False`, default):

| Parameter | Shape | Init | Meaning |
|---|---|---|---|
| `unit_latent` | (U, D) | `N(0,1)` | Position of the unit in chemical latent space |
| `base_energy_u` | (U,) | `E[ln c]` | log(EC50) at the perfectly matching ligand |
| `max_energy_u_raw` | (U,) | `softplus⁻¹(10)` | Raw param; `softplus(·)` gives E_max (selectivity ceiling) |

**Interface model** (`use_interface_model=True`):

Each unit has two faces reflecting the two interfaces it participates in.  Binding
pockets sit at the interface between adjacent subunits in the pentameric ring.

| Parameter | Shape | Init | Meaning |
|---|---|---|---|
| `unit_latent_plus` | (U, D) | `N(0,1)` | + face embedding (contributed to interface with the *next* unit) |
| `unit_latent_minus` | (U, D) | `N(0,1)` | − face embedding (contributed to interface with the *previous* unit) |
| `base_energy_u_plus/minus` | (U,) | `E[ln c]` | Base energy of each face |
| `max_energy_u_raw_plus/minus` | (U,) | `softplus⁻¹(10)` | Selectivity ceiling of each face |

Pocket embedding for interface *i* of receptor *r* = (u₀, …, u_{k−1}):
```
v_pocket_i = 0.5 * (unit_latent_plus[u_i] + unit_latent_minus[u_{i+1 mod k}])
E_base_i   = 0.5 * (base_energy_u_plus[u_i]  + base_energy_u_minus[u_{i+1 mod k}])
E_max_i    = 0.5 * (softplus(max_raw_plus[u_i]) + softplus(max_raw_minus[u_{i+1 mod k}]))
E_o_i      = E_base_i + E_max_i * (1 − exp(−‖v_pocket_i − v_obs‖² / λ²))
```

Note: the kernel is applied to the **averaged** pocket embedding, not to each face
separately — `f(0.5*(a+b)) ≠ 0.5*(f(a)+f(b))` in general.

**Init rationale (both models):** `base_energy = E[ln c]` ensures EC50 ≈ typical ligand
concentration at iteration 0.  `E_max = 10` places the saturation far above any training
concentration, freeing the optimizer to tune selectivity from an unconstrained start.

### 1.4 Affinity Kernel

Two kernel choices (§03):
- **gaussian** (default): `E_o = E_base + softplus(max_raw) * (1 − exp(−d²/λ²))`
  — saturates at large distance; biophysically correct.
- **quadratic** (legacy): `E_o = E_base + softplus(slope_raw) * d²`
  — grows without bound; kept for backward compatibility.

---

## Stage 2 — Batch Sampling (`environment.py::LigandEnvironment.sample_batch`)

`sample_batch(batch_size, receptor_indices=None)` is composed of
three private helpers so each concern is testable in isolation:

| Helper | Returns | Notes |
|---|---|---|
| `_sample_masks(B)` | `((B,L) mask, (B,s_upper) sparse_idx)` | hierarchical draw; `sparse_idx` carries selected ligand indices |
| `_sample_noisy_ligands(B, sparse_idx)` | `(B, s_upper, D)` | noisy coords for present ligands only; dummy slots zeroed |
| `_compute_energies(v_ligands, receptor_indices)` | `(B,s_upper,U)` or `(B,s_upper,R,k)` | distance trick + kernel; second dim is now s_upper, not L |

`s_upper` is a static bound computed once in `__init__`:
`s_upper = min(L, (⌊μ_src + 4√μ_src⌋ + 1) × max_m)` covering the ~99.99th percentile of
Poisson(mu_sources) × max_m without a per-batch GPU→CPU sync.
For the default config (mu_sources=4, max_m=10): s_upper = 130 vs L = 2000, a **15× reduction**
in the dominant einsum and its backward.  Padding slots carry global index L with concentration 0,
contributing exp(−27) ≈ 0 to the logsumexp in physics — numerically invisible.

### 2.1 Presence sampling — hierarchical source→ligand model

Controlled by `mu_sources`, `mu_ligands_per_source`, and `n_presence_blocks` in
`SingleRunConfig`/`RunConfig`.

**Block partition** — `n_ligands` ligands are divided into `K = n_presence_blocks`
source blocks at environment construction.  Block membership is assigned by a seeded
random permutation with seed `n_ligands × 131071 + n_presence_blocks` (deterministic
from those two values alone), independently of `ligand_family_assignments`.
Stored as buffer `presence_block_id (L,)`.  Helper buffers `_block_members (K, max_m)`,
`_block_valid (K, max_m)`, and `_block_sizes (K,)` index into the partition.
Three additional precomputed buffers enable fully vectorized mask sampling:
`_safe_block_members (K, max_m)` — invalid slots redirected to dummy index L;
`_log_pmf_source (K,)` — normalized ZTP log-PMF for the n_src draw;
`_log_pmf_ligand (K, max_m)` — per-block normalized ZTP log-PMF for n_lig draws.

**Per-batch draw** (implemented in `_sample_masks`) —

1. **Active sources.** Draw `n_src ~ ZTP(mu_sources, K)` (zero-truncated Poisson
   with rate `mu_sources`, support `[1, K]`). Select `n_src` distinct source blocks
   uniformly without replacement via Gumbel-top-k.

2. **Ligands per active source.** For each active block `k` (size `m_k`), draw
   `n_lig ~ ZTP(mu_ligands_per_source, m_k)`. Select `n_lig` distinct ligands within
   block `k` uniformly without replacement via Gumbel-top-k.

3. Set `M[b, picked] = 1`, all else 0.

Zero-truncated Poisson sampling on `[1, n_max]`: unnormalized log-PMF
`s·ln(λ) − ln(s!)` for `s = 1..n_max` is normalized and stored as a buffer at init.
At sample time the **Gumbel-max trick** is used — `argmax(log_pmf + Gumbel(0,1))` gives
an exact categorical draw without `torch.multinomial`, enabling full vectorization across
all K blocks simultaneously with no Python loop and no GPU→CPU synchronization.

Every row has `S ≥ 1` by construction (step 1 forces at least one active source, step
2 forces at least one ligand per source).  No rejection loop is needed; `sample_batch`
asserts `masks.sum(-1).min() >= 1` as a safeguard.

**Knobs:**

| Parameter | Effect at small values | Effect at large values |
|---|---|---|
| `mu_sources` | ~1 active source (sparse source mixture) | many sources active |
| `mu_ligands_per_source` | ~1 ligand per source (pure singleton → pure source) | many ligands per source |

Corner cases: `K = 1` ("one source, mu_ligands_per_source controls mixture size");
`K = n_ligands` (each ligand its own source, `mu_ligands_per_source` then forced to 1).

**Block-shared concentration mean** (`block_shared_conc_mean = True`, gated on
`n_presence_blocks > 1`): within each source block all ligands inherit the same
concentration mean (the average of their individually configured means).  Per-sniff
concentrations still draw independently around that shared mean.  This matches the
turbulent-transport assumption: co-occurrence survives but concentration ratios do not.

### 2.2 Sampling guarantee

The hierarchical sampler guarantees `S = sum(M, dim=-1) >= 1` for every row by
construction (step 1 always selects ≥ 1 source; step 2 always selects ≥ 1 ligand
per source).  There is no empty-mixture rejection loop.  `sample_batch` checks this
invariant with a single `assert masks.sum(-1).min() >= 1`.

**Energy computation** — the identity trick avoids allocating an `(B, L, *, D)` tensor:
```
dist_sq = ||a||² + ||b||² - 2 <a, b>
```

**Classic model** returns `E_open: (B, L, U)` — one energy per unit per ligand observation.

**Interface model** (`use_interface_model=True`, requires `receptor_indices`): computes pocket
embeddings for every interface of every receptor, then returns `E_open: (B, L, R, k_sub)` —
already gathered per receptor.  The peak intermediate tensor is `(B, L, R·k_sub)`, avoiding
a `(B, L, U, U)` all-pairs matrix.

`BinaryReceptor.forward` accepts a `pre_gathered=True` flag and skips the index-gather step
when the interface model's output is passed directly.

---

## Stage 3 — Receptor Physics (`physics.py::BinaryReceptor.p_open`)

The BinaryReceptor implements the EC50-threshold limit of the MWC model (§02.4).

**Classic model:**
```
ln EC50^(r,ℓ) = (1/k_sub) Σ_u E_open^(u,ℓ)        [geometric mean of subunit affinities]

ln_sum_terms = logsumexp_ℓ [ ln(c_ℓ) − ln EC50^(r,ℓ) ]   [mixture aggregation]

p(active) = sigmoid( ln_sum_terms / T )
```

**Why logsumexp?** The mixture binding polynomial is a sum of concentration/affinity
ratios. Taking the log and using logsumexp is numerically stable when individual terms
span many orders of magnitude.

**Temperature calibration** (`compute_initial_temperature`): T_init is set to the
empirical std of `ln_sum_terms` over a calibration batch. This gives the sigmoid
argument unit standard deviation at the start of training, avoiding both the saturated
regime (all 0/1, vanishing gradient) and the mushy regime (all 0.5, non-discriminating).

**Interface model:** `ln EC50^(r,ℓ) = (1/k_sub) Σ_i E_pocket^(i,ℓ)` where the sum runs
over interfaces rather than subunits.  `BinaryReceptor.p_open` is unchanged; it always
averages over the last dimension of `energies_k`, regardless of whether that dimension
indexes subunits (classic) or pocket interfaces (interface model).

**Temperature calibration** (`compute_initial_temperature`): T_init is set to the
empirical std of `ln_sum_terms` over a calibration batch. This gives the sigmoid
argument unit standard deviation at the start of training, avoiding both the saturated
regime (all 0/1, vanishing gradient) and the mushy regime (all 0.5, non-discriminating).
In interface mode, the calibration batch is drawn via `sample_batch(..., receptor_indices)`.

**Quadrature for dose-response** (`get_dose_response`): when `distribution_type='gaussian'`,
the function integrates over the observation noise using Gauss-Hermite quadrature
(§06.2). If the grid would exceed 100,000 points (i.e., `quadrature_degree^latent_dim > 1e5`),
**or when `use_interface_model=True`** (pocket embeddings are pair-wise, making the
quadrature grid impractical), it falls back to the mean-energy approximation.

---

## Stage 3b — Cell Readout (`cells.py`) — optional

Present only when the config declares cells (`cell_gene_sets` or `n_cells` is set;
`SingleRunConfig.is_cell_mode()`). It shifts the sensory unit from the RECEPTOR to
the CELL: instead of listing R receptors, one lists C cells, each defined by the set
of subunit genes it expresses. A cell assembles **every** receptor its genes allow.

### 3b.1 Repertoire expansion

A cell $c$ expressing the gene set $G_c$ with $g = |G_c|$ assembles

- **standard model** — the unordered multisets of size $k_{sub}$ over $G_c$:
  $\binom{g + k_{sub} - 1}{k_{sub}}$ receptors.
  $G_c = \{0,2\}$, $k_{sub}=5$ gives the 6 receptors
  $(0,0,0,0,0)$, $(0,0,0,0,2)$, $(0,0,0,2,2)$, $(0,0,2,2,2)$, $(0,2,2,2,2)$, $(2,2,2,2,2)$.
- **interface model** (`use_interface_model=True`) — the canonical cyclic ring
  arrangements, rotations identified and reflections distinct (same convention as
  `geometry.generate_ordered_receptor_indices`, §"Receptor Sampling Strategies"):
  $\frac{1}{k_{sub}} \sum_{d \mid k_{sub}} \varphi(d)\, g^{k_{sub}/d}$ receptors,
  where $\varphi$ is Euler's totient function. For $g=2$, $k_{sub}=5$: 8 arrangements.

`CellArray` expands every cell, then takes the **deduplicated union** as
`receptor_indices` of size $R_{pool}$ — a receptor shared by several cells is
simulated once — and records the abundances in a matrix $W$ of shape $(C, R_{pool})$
with $W_{cr} = 0$ when cell $c$ cannot assemble receptor $r$.

### 3b.2 Abundance ($W$), `cell_stoichiometry`

Rows of $W$ are normalised: $\sum_r W_{cr} = 1$. A cell has a fixed total number of
receptors regardless of how many genes it expresses; only their diversity changes.

- `'multinomial'` (default) — **random assembly**: all subunits are produced in equal
  amounts and assemble independently into the $k_{sub}$ slots, so each of the $g^{k_{sub}}$
  ordered words is equally likely and a receptor's abundance is proportional to how
  many words collapse onto it. For a receptor with $n_i$ copies of gene $i$:
  $$ W_{cr} \;\propto\; \frac{k_{sub}!}{\prod_i n_i!} \qquad\text{(standard model)} $$
  $$ W_{cr} \;\propto\; \#\{\text{distinct rotations of the ring word}\} \qquad\text{(interface model)} $$
  For $g=2$, $k_{sub}=5$ this gives $1,5,10,10,5,1$ out of 32 — homomers are strongly
  suppressed relative to mixed receptors, with no free parameter.
- `'uniform'` — every receptor type the cell can make is present in equal amount,
  $W_{cr} = 1 / |{\rm repertoire}|$. Implies per-type assembly control rather than
  random mixing.

### 3b.3 Activation, `cell_readout`

The abundance-weighted open fraction (the **drive**) is
$$ S_{bc} \;=\; \sum_{r} W_{cr}\, p_{br} $$
with $p_{br}$ the open probability of receptor $r$ on sniff $b$ from Stage 3.

- `'threshold'` (default) — the cell's ionic current is proportional to the number of
  its open receptors, and it fires above a threshold:
  $$ A_{bc} \;=\; \sigma\!\left( \frac{S_{bc} - \theta_c}{T_{cell}} \right) $$
  $A_{bc}$ is a genuine firing probability, which is what the Bernoulli-mixture
  entropy estimators of Stage 4 assume.
- `'noisy_or'` — the cell is active if any of its receptors opens:
  $A_{bc} = 1 - \exp\!\big(\sum_r k_{sub} W_{cr} \ln(1 - p_{br})\big)$. No threshold
  parameter, but saturates toward 1 for cells with large repertoires.
- `'mean'` — $A_{bc} = S_{bc}$. Diagnostic only: a mean of probabilities is not a
  firing probability, so the Stage 4 estimators over-read it.

**Calibration** (`calibrate_cell_readout`, threshold mode only). Same rationale as
`compute_initial_temperature` in Stage 3: a threshold off the support of $S$ leaves
every cell permanently silent or saturated and the array carries zero entropy.
- $\theta \leftarrow$ `median_threshold`$(S)$ — ONE scalar shared by all cells,
  pinned to the median of the pooled drive so the population fires ~50% of the time.
  Sharing it is deliberate: a per-cell threshold forces every cell to the same firing
  rate, erasing the heterogeneity the cell picture exists to study.
- $T_{cell} \leftarrow$ `drive_scale`$(S, \theta)$ — unit spread in the sigmoid argument.

Both helpers live in `cells.py` and both exist to dodge a specific way the sigmoid
degenerates. They are described in §3b.7.

$\theta$ is **not** learnable by default (`cell_threshold_learnable=False`): it is
pinned to the data, not fitted, so it costs no free parameter. It is **re-pinned** every
`cell_recalibrate_every` epochs, because the drive distribution moves as the chemistry
trains and a median measured at epoch 0 goes stale. Re-pinning also keeps the sigmoid's
transition band on the densest part of the drive, which is where the phase-2 gradient
comes from. A learnable $\theta$ collapses to the fair-coin degeneracy of §3b.6 and is
kept only as an ablation.

**Both $T_{cell}$ endpoints are MULTIPLES of the live drive spread**, never absolute
values: $S$ is a weighted mean of probabilities, so its scale is set by the environment.
The endpoints are recomputed whenever the spread is re-measured.

### 3b.3b Two-phase sharpening schedule

Annealing the receptor and the cell together makes the cell's operating point chase a
drive distribution that is still moving underneath it. In cell mode the schedule is
therefore split at `cell_phase_split` (default 0.5):

| | receptor sharpness $T$ | cell sharpness $T_{cell}$ |
|---|---|---|
| **Phase 1** | anneals $T_{init} \to T_{final}$, reaching $T_{final}$ on the LAST epoch of the phase | held soft at $1\times$ the live spread |
| **Phase 2** | held at $T_{final}$ | anneals down to `cell_temperature` $\times$ spread |

Phase 1 keeps gradients live everywhere while the chemistry arranges the drive around
the threshold. Phase 2 hardens the readout. `cell_temperature` defaults to 0.01, leaving
~0.8% of $(b,c)$ pairs inside the transition band — cells are then effectively
deterministic, which is the condition under which the entropy objective is a valid proxy
for information (§3b.6). Receptor-only runs keep the original 80%-of-training schedule.

The per-epoch measurement always evaluates at the FINAL sharpness, so the convergence
curve during phase 1 already reports the honest hard-readout number.

### 3b.6 Why the cell must be deterministic

The estimators of Stage 4 read every activity as the probability that a coin lands
heads, and charge its binary entropy into the total. Decompose the reported entropy as
$H = H(s \mid \text{sniff}) + I(s; \text{sniff})$: only the second term is information.
A cell sitting at $A = 0.5$ contributes a full bit of the first term and nothing to the
second — the KT self-test `[3b]` in `bin_loss.py` shows this exactly, reporting $H = R$
for an array of fair coins that responds to nothing.

This never bit the receptor picture because the receptor temperature anneals until
$p \in \{0, 1\}$, so $H(s \mid \text{sniff}) \approx 0$ and entropy $\approx$ information.
Cells add a *second* softness, and if $T_{cell}$ is comparable to the spread of $S$ the
reported entropy is almost entirely noise. Hence the phase-2 hardening above.

Two consequences for the readout menu of §3b.3:
- `'mean'` cannot escape this. $S$ is an average over the repertoire, and averages
  concentrate, so it never approaches 0 or 1. Unusable with a binary estimator.
- `'noisy_or'` saturates toward 1 for large repertoires — nearly binary, but stuck ON.

### 3b.7 Two ways the threshold degenerates, and the two fixes

Both come from the same root cause. The drive $S$ is a weighted mean of receptor open
probabilities, and a *sharp* receptor is off unless its ligand is present, so those
probabilities are sigmoids saturated deep into their tails. The resulting drive is
wildly non-uniform: measured on a toy at the final receptor sharpness, **26% of
$(b,c)$ pairs were exactly zero, the median was $\approx 10^{-30}$, and the maximum was
$\approx 1$ — thirty orders of magnitude apart.** Placing a sharp sigmoid on a
distribution shaped like that needs care in two separate places.

**(a) The threshold must not sit on a point mass — `median_threshold`.**
When more than half the drives share one value (typically exactly 0, the sparse-code
case), the plain median IS that value. A sharp sigmoid centred on a point mass returns
$\sigma(0) = 0.5$ for every member of it: the whole mass becomes fair coins and the
array reports maximum entropy while carrying nothing (§3b.6).

The rule is one line: **$\theta$ is the midpoint between the median and the next
distinct value above it.** For a continuous drive this changes nothing of substance —
`torch.median` returns the lower of the two middle order statistics, and the midpoint of
those two is the textbook median, so we merely pick the interior of the median interval
instead of its lower endpoint. For a drive with a point mass at the bottom, the same
rule steps $\theta$ *strictly above* the mass, so those sniffs read cleanly OFF and the
code is sparse but honest. The trap is narrow: any $\theta$ above the mass works, only
$\theta$ exactly on it produces coins.

**(b) The sharpness must be measured where the threshold is — `drive_scale`.**
$T_{cell}$ is a fraction of "the spread of the drive", but *which* spread matters
enormously. A standard deviation of the drive above is $\approx 4 \times 10^{-2}$ — set
entirely by the handful of strongly-firing sniffs, and saying nothing about the
$10^{-30}$ region the threshold actually occupies. Using it gives a $T_{cell}$ that
swamps the middle of the distribution, and every sniff there evaluates to
$\sigma(\approx 0) = 0.5$.

Measured, on the same drive with a correctly-placed $\theta$:

| spread measure | value | $T_{cell}$ | fair coins | firing |
|---|---|---|---|---|
| standard deviation | 3.7e-02 | 3.7e-04 | **91%** | 14% |
| interquartile range | 1.1e-18 | 1.1e-20 | **67%** | 43% |
| **median absolute deviation** | 1.1e-30 | 1.1e-32 | **0%** | 50% |

The IQR is better but still fails, because the upper quartile sits twelve orders of
magnitude above the median. The **MAD** is the deviation of the *typical* point from the
median, so it tracks the middle by construction. `drive_scale` returns the MAD.

A corollary worth stating: **no absolute floor may be applied to $T_{cell}$.** A guard
like `max(T, 1e-6)` looks harmless but reinstates the whole failure when the drive lives
at $10^{-30}$. The floor in `run.py` is $10^{-45}$, just off zero.

Degenerate case: if every drive is identical, no threshold can separate anything.
`median_threshold` flags it (`on_point_mass` in the returned diagnostics, surfaced in the
training log) rather than silently patching, because when the drive has zero spread
$T_{cell}$ collapses too and no nudge can outrun it.

End to end, on a 40-epoch toy run: before these two fixes, 98.5% of activities were fair
coins; after, **0%, with 49% of sniffs firing** — the design target.

### 3b.4 Memory: the receptor pool is the bottleneck

The physics tensors are $(B, s_{upper}, R_{pool}[, k_{sub}])$, and $R_{pool}$ — the
union of every cell's repertoire, capped at $\binom{N_{genes}+k_{sub}-1}{k_{sub}}$ —
grows far faster than $C$. Because all three readouts reduce over the pool through a
term **linear** in the per-receptor contribution, `cell_activity` walks the pool in
slices of `cell_pool_chunk` receptors and sums the $(B, C)$ accumulator, so the peak
forward tensor is $O(B \cdot s_{upper} \cdot \text{chunk})$ instead of
$O(B \cdot s_{upper} \cdot R_{pool})$. During training the graph is still retained
unless `recompute_backward=True` gradient-checkpoints each chunk; without it, chunking
only helps under `no_grad` (evaluation).

`resolve_batch_sizes` is called with `n_receptors=C` and `n_physics_receptors=R_pool`:
the estimator caps (which are exponential or quadratic in the number of channels) size
on the C cells, while only the physics cap sees $R_{pool}$.

### 3b.5 Sampling gene sets

`build_cell_array` takes explicit `gene_sets`, or samples `n_cells` of them:
- `'bernoulli'` — gene $u$ expressed independently with probability `cell_gene_probs[u]`;
  empty cells rejected, `cell_max_genes` caps the repertoire size.
- `'size_pmf'` — two stage: draw the *number* of expressed genes from `cell_size_pmf`
  (entry $i$ = probability of $i+1$ genes), then draw that many genes uniformly without
  replacement. This is the hook for an arbitrarily complex expression model — a
  distribution conditioned on how many genes are already expressed is exactly a choice
  of `cell_size_pmf`.

Sampling happens **once**, in `SingleRunConfig.__post_init__`, which writes the
resolved gene sets back into `cell_gene_sets` and the derived pool into
`receptor_indices`. Both are persisted in `config.json`, so a run is exactly
reproducible from it; `SimulationRunner._initialize` rebuilds $W$ from the stored gene
sets and asserts the pool matches.

---

## Stage 4 — Loss Computation

Loss modules are selected by `cfg.entropy` in `run.py::_build_loss`.

Every loss and every measurement consumes a `(B, N)` activity tensor and is agnostic
to whether `N` counts receptors (Stage 3) or cells (Stage 3b) — the two pictures
differ only at `SimulationRunner._activity`. Read `R` below as "number of channels
in the array": it is `C` in cell mode.

### 4a. `DiscreteExactLoss` (`bin_loss.py`) — maximize joint array entropy

**Objective:** `min −H(A)` where A is the binary array activity (or `min C` for collision).

Five entropy estimators — pick based on array size:

| `entropy_type` | Complexity | When to use |
|---|---|---|
| `'shannon'` | O(B · 2^R) | R < ~15; exact but exponential |
| `'collision'` | O(B² · R) | Scalable lower bound; training loss minimises C directly |
| `'kt'` | O(B² · R) time, O(chunk² · R) mem | Certified Shannon lower bound (Kolchinsky–Tracey Bhattacharyya); tight when components separate. Ceiling log2(total B) via diagonal anchoring |
| `'blocked'` | O(B · 2^block_size · R/block) | Upper bound, captures within-block correlations |
| `'blocked_corrected'` | O(B · 2^block_size + R²) | Tighter: blocked minus cross-block pairwise MI |
| `'proxy'` | O(B · R²) | Fastest; pairwise covariance/repulsion penalty |

**Collision trick:** computes log P(collision) = Sigma_r log P_r(collision) in log-space
via logsumexp, then exponentiates once. Training loss returns C = exp(log_mean_coll_prob)
directly (no log), removing the 1/C gradient blow-up at low collision probability.
Measurement returns H2 = -log2(C) in bits. For B > collision_chunk_size, the batch is
split into ceil(B / collision_chunk_size) chunks; cross-chunk log-probabilities are
concatenated and reduced via a single logsumexp so the full batch is consumed (no
samples discarded). The collision_chunk_size is sized to GPU memory at init time by
``resolve_batch_sizes`` (see §06.4 for the formula); it defaults to 2048 when not
configured. Chunk size (quadratic memory) raises the bit ceiling; chunk count (linear)
only reduces variance.

**KT (Kolchinsky–Tracey):** certified lower bound on the *same* Shannon H(s) that
`shannon` computes exactly, via pairwise Bhattacharyya affinities.
H_KT = H_cond - (1/B) Sigma_i log2( (1/B) Sigma_j BC(i,j) ), with
H_cond = mean_b Sigma_r h2(A_br) and log BC(i,j) = Sigma_r log(sqrt(A_i A_j) +
sqrt((1-A_i)(1-A_j))). ALL pairs, **diagonal kept**, normalised by **total** B: since
BC(i,i)=1, the self-term anchors the inner sum at 1/B and the resolvable ceiling is
log2(total B) — not log2(chunk) as for collision. Two nested chunk loops (outer i-chunks,
inner j-chunks over the whole batch, cat + logsumexp) give O(B²·R) time (quadratic in
total batch) and O(chunk²·R) memory; it reuses `collision_chunk_size` (identical memory
scaling). Dispatched via the generic `-compute_entropy(activity)` path (no collision-style
C-return trick). The config flag `recompute_backward=True` gradient-checkpoints the inner
loop so the training batch is bounded by compute rather than retained memory (see
Batch-size auto-scaling below).

**Blocked-corrected:** H_blocked - Sigma_{cross-block (i,j)} I(A_i; A_j). The pairwise
MI matrix is computed via the batch Gram matrix A^T A / B (vectorised, differentiable).
See §5.5 of `05_optimization.md`.

**Blocked trick — correlation-aware partitioning:** groups receptors by absolute Pearson
correlation affinity via greedy clustering (see §5.4 of `05_optimization.md`).  Each
block captures the most correlated receptors together, giving a tighter upper bound on
the true joint entropy than random partitioning. The partition is computed on detached
activity (stop-gradient) and cached for `block_refresh_interval=50` training steps to
avoid gradient instability. Evaluation calls (`use_cache=False`) always build a fresh
partition from the eval batch. Shared by `blocked` and `blocked_corrected`.

### 4b. `AnnealedEntropyLoss` (`annealed_loss.py`) — annealed blocked -> collision

**Objective:** `min -[(1 - lam) H_blocked + lam H_collision]` where `lam = epoch / epochs`.

Combines the blocked Shannon estimator and collision H2 via linear interpolation.
Selected via `cfg.entropy = 'annealed'`.

### 4b'. `BlockedToCorrectedLoss` (`annealed_loss.py`) — annealed blocked -> blocked-corrected

**Objective:** `min -[(1 - lam) H_blocked + lam H_blocked_corrected]` where `lam = epoch / epochs`.

Stays within the Shannon family while tightening the bound over training. A
`lam_override=1.0` gives pure blocked-corrected from epoch 0. Selected via
`cfg.entropy = 'blocked_to_corrected'`.

Both annealed losses receive `(activity, epoch, epochs)` in `_train` and provide
`compute_entropy(activity, entropy_type)` for measurement helpers.

### 4c. `MaximizeMutualInformationLigandLoss` (`family_mi_loss.py`) — MI(array; mixture)

**Objective:** `max I(A ; M) = H(A) − H(A | M)`

H(A) is computed on the full batch. H(A | M) conditions on the exact mixture identity:
mixture masks are hashed to integer IDs using binary powers
`id = Σ_ℓ M_ℓ · 2^ℓ`, then the batch is grouped by ID and entropy computed per group.

### 4d. `MaximizeMutualInformationConcentrationLoss` (`concentration_mi_loss.py`)

**Objective:** `max I(A ; C) = H(A) − H(A | C)`

Sorts the batch by total concentration, divides into `n_c_bins` quantile bins,
computes entropy per bin.

---

## Stage 5 — Training Loop (`run.py::SimulationRunner._train`)

```
for epoch in range(epochs):
    1. anneal temperature   frac = min(1, epoch / (0.8*epochs))
                            T = T_end + (T_start − T_end) * (1 − frac)
                            # reaches T_end at 80% of epochs, holds it for the last 20%
    2. sample batch         E, concs, masks = env.sample_batch(batch_size)
    3. compute activity     activity = physics(E, concs, receptor_indices)
    4. compute loss         loss = loss_fn(activity, ...)     # see dispatch below
    5. backward + step      loss.backward(); optimizer.step()
    6. log every 1%         per_epoch_measure ? _eval_stats(...) : {loss, train_entropy=-loss}
```

Steps 2–5 are wrapped in `record_function("prof:…")` labels (sample+physics_fwd / loss_fwd /
backward) — inert unless a `torch.profiler` run is active (`tasks/profiling/`), used to group ops
in the profiler trace. Step 6 logs the free training objective instead of a test measurement when
`per_epoch_measure=False` (see Batch-size auto-scaling).

**Loss dispatch** in step 4 depends on the loss module type:
- `DiscreteExactLoss`: `loss_fn(activity)`
- `AnnealedEntropyLoss`: `loss_fn(activity, epoch, self.config.epochs)` — needs
  training progress to compute the interpolation parameter λ.
- `MaximizeMutualInformationLigandLoss`: `loss_fn(activity, mixture_masks=masks)`
- `MaximizeMutualInformationConcentrationLoss`: `loss_fn(activity, concs=concs)`

**Temperature annealing:** starts at `T_init` and decreases linearly to the configured
`temperature`, reaching it at **80% of epochs** and then holding it for the final 20%.
The hold matters because `_eval_stats` always measures at `temperature` (the sharp T):
without it the training T equals the eval T only at the last epoch, so the sharp-T
objective is under-optimized and the logged entropy spikes at the buzzer. A high T keeps
the sigmoid soft early in training (smooth gradients); a low T sharpens to binary
decisions later. `T_init` is controlled by
`initial_temperature` in the config: set to `"auto"` (default) to calibrate it as the
empirical std of pre-sigmoid terms via `compute_initial_temperature`, or to an explicit
float to override (useful when few ligands make the auto-calibration underestimate the
scale needed for exploration).

**Warm-starting** (`SweepRunner`): `warm_start: bool` in `RunConfig` controls whether
chain warm-starting is applied along the trajectory.  At each step after the first,
`SweepRunner.execute` applies:

| Condition | Action |
|---|---|
| `warm_start=False` or first step | **Cold start:** environment built from scratch. |
| `curr.n_genes > prev.n_genes` | **Chain warm-start:** pass `prev_env` forward; `_initialize` calls `clone_with_extra_units` to expand the gene pool. LR is damped 10×. |
| `curr.n_genes <= prev.n_genes` | **Cold start:** n_genes decreased → new env group boundary, reset. |

**Typical usage:**
```python
# Gene-growth sweep — chain warm-start at each n_genes step.
RunConfig(n_genes=[3,4,5,6,7,8], conc_mean=(...), warm_start=True, ...)

# Fixed n_genes, vary other parameters — always cold.
RunConfig(n_genes=5, entropy=["collision","shannon"], conc_mean=(...), warm_start=False, ...)
```

`_initialize` always builds `receptor_indices` fresh from `self.config.receptor_indices`,
which is auto-generated in `SingleRunConfig.__post_init__` from `(n_genes, n_receptors,
receptor_sampling_strategy, receptor_sampling_seed)` when `receptor_indices is None`.
The LR is damped 10× on warm-start to preserve learned representations.

**Optional cosine LR scheduler:** wraps Adam with `CosineAnnealingLR` when
`use_scheduler=True`, annealing from `lr` to `1e-5` over the full run.

---

## Stage 6 — Evaluation & Metrics (`run.py::SimulationRunner._eval_stats`)

Evaluation runs under `torch.no_grad()` with `test_batch_size` samples.
`_eval_stats` draws a single **mixture batch** (natural Bernoulli masks) before the
measurement loop.  All metrics — including conditional-entropy and MI measurements —
operate on this same batch.

When `eval_chunk_size < test_batch_size`, soft metrics use one chunk and hard codeword
metrics are accumulated across all chunks (avoiding CUDA OOM on large eval budgets).

`full_array_entropy` calls `loss_fn.compute_entropy(act, entropy_type=..., use_cache=False)`
so evaluation always builds a fresh correlation-aware partition from the eval batch,
never reading or writing the training cache.

Available metrics (add to `measurement_fns` in config):

| Key | Function | What it measures |
|---|---|---|
| `full_array_entropy` | `analysis_helper.full_array_entropy` | The loss's NATIVE joint-entropy estimator only (`loss_fn.compute_entropy` default: collision for a collision loss, blocked for annealed, blocked_corrected for blocked_to_corrected, kt for a kt loss). Logs a single `full_array_entropy` column. |
| `entropy_collision` / `entropy_blocked` / `entropy_blocked_corrected` / `entropy_kt` / `entropy_kt_upper` | `analysis_helper.entropy_*` | Opt-in extra estimators, each logging one column (`full_array_entropy_collision` / `_blocked` / `_blocked_corrected` / `_kt` / `_kt_upper`). Add to `measurement_fns` to record estimators other than the loss's own (e.g. the `optimizer` task requests all to compare losses on a common footing). Each is a full extra evaluation. `entropy_kt` / `entropy_kt_upper` are computed on the soft assignment (not via `compute_entropy`, which `AnnealedEntropyLoss` rejects for `'kt'`). KT is all-pairs O(B²·R) with resolvable-entropy ceiling log₂(B); it is measured on the **full eval batch** `test_batch_size` (see below), generated in sub-batches and tiled internally so peak memory is bounded by the tile, not B — **no sample cap; only time grows (O(B²))**. `entropy_kt` is the Bhattacharyya **lower** bound, `entropy_kt_upper` the KL-divergence **upper** bound (clamped to R bits); together they bracket the true joint entropy H(s) — certified, and tight in the well-separated (optimized) regime. |
| `codeword_entropy` | `analysis_helper.codeword_entropy` | Hard plug-in + Miller-Madow entropy of binary codewords |
| `mean_receptor_distance` | `analysis_helper.mean_receptor_distance` | Average pairwise latent-space distance between receptors |
| `receptor_distances` | `analysis_helper.receptor_distances` | Full (R, R) pairwise distance matrix |
| `conditional_entropy_ligand` | `analysis_helper.conditional_entropy_ligand` | (1/N_l) Σ H(A \| L_l) |
| `mutual_information_ligand` | `analysis_helper.mutual_information_ligand` | (1/N_l) Σ I(A ; L_l) |
| `conditional_entropy_concentration` | `analysis_helper.conditional_entropy_concentration` | mean over present ligands of H(A \| C_l, l present) — dense per-ligand conc |
| `mutual_information_concentration` | `analysis_helper.mutual_information_concentration` | mean over present ligands of I(A ; C_l \| l present) — concentration (level) coding |
| `concentration_channel` | `analysis_helper.concentration_channel` | H(A \| M) — condition on the full presence pattern; concentration channel I(A;c\|M), total-comparable |
| `identity_channel` | `analysis_helper.identity_channel` | I(A ; M) = H(A) − H(A \| M) — composition/identity channel, total-comparable (identity_channel + concentration_channel = H(A)) |
| `conditional_entropy_family` | `analysis_helper.conditional_entropy_family` | (1/N_f) Σ H(A \| F_f) — latent-space families |
| `mutual_information_family` | `analysis_helper.mutual_information_family` | (1/N_f) Σ I(A ; F_f) |
| `conditional_entropy_block` | `analysis_helper.conditional_entropy_block` | (1/N_b) Σ H(A \| B_b) — source blocks |
| `mutual_information_block` | `analysis_helper.mutual_information_block` | (1/N_b) Σ I(A ; B_b) |
| `rank_ordered_distances` | `analysis_helper.rank_ordered_distances` | Rank-ordered energy gap from preferred ligand |
| `mean_specialization_index` | `analysis_helper.mean_specialization_index` | S_r = (A_max − A_bg)/(A_max + A_bg) |
| `receptor_conditioned_entropy` | `analysis_helper.receptor_conditioned_entropy` | H(M \| a_r > 0.5) — mixture uncertainty when receptor fires |

**Family labels:** `conditional_entropy_family` and `mutual_information_family` require
`family_labels: (B, n_families)` bool, derived lazily in `_eval_stats` from
`env.ligand_family_assignments` and `mixture_masks`:
```python
family_labels[b, f] = True  iff  any ligand from family f is present in sample b
```
No extra batch sample is drawn; the tensor is computed from the already-sampled `masks`.

**Block labels:** `conditional_entropy_block` and `mutual_information_block` require
`block_labels: (B, n_presence_blocks)` bool, derived lazily in `_eval_stats` from
`env.presence_block_id` and `mixture_masks`:
```python
block_labels[b, k] = True  iff  any ligand from source block k is present in sample b
```
Computed identically to `family_labels` but via `env.presence_block_id`.  Since the
block partition is orthogonal to the family partition by construction, comparing
`mutual_information_block` with `mutual_information_family` disentangles the receptor's
sensitivity to *source co-occurrence* (blocks) from its sensitivity to *chemical
similarity* (families).  With the hierarchical sampler, `block_labels[b, k] = True`
exactly when source block `k` was drawn as an active source in sniff `b`.

**`conditional_entropy_family` / `mutual_information_family` — marginal conditioning:**
Each family f is treated as an **independent binary variable** F_f ∈ {0,1}.  For
each f the batch is split into present/absent groups and:

```
H(A | F_f) = P(f=1)·H(A | f present) + P(f=0)·H(A | f absent)
```

The function returns `(1/N_f) Σ_f H(A | F_f)`, so `mutual_information_family`
returns the average pairwise MI: `(1/N_f) Σ_f I(A ; F_f)`.

This avoids the combinatorial explosion of the joint presence pattern (2^N_f
possible groups), which would require enormous batches to estimate reliably.

**`conditional_entropy_ligand` / `mutual_information_ligand` — marginal conditioning:**
Each ligand l is treated as an **independent binary variable** L_l ∈ {0,1}.  For
each l the batch is split into present/absent groups and:

```
H(A | L_l) = P(l=1)·H(A | l present) + P(l=0)·H(A | l absent)
```

The function returns `(1/N_l) Σ_l H(A | L_l)`, so `mutual_information_ligand`
returns the average pairwise MI: `(1/N_l) Σ_l I(A ; L_l)`.

This avoids the combinatorial explosion of the joint mixture pattern (2^N_l
possible groups), which would require enormous batches to estimate reliably.

**`conditional_entropy_concentration` / `mutual_information_concentration` — per-ligand,
present-only, quantile binning:**
Each ligand l is scored on a DENSE `(B, L)` per-ligand concentration `concs_dense`
(identity-aligned, 0 for absent), built by `sample_batch(return_dense_conc=True)` from the
same `sparse_idx`/`concs` that produced the activity. For ligand l we keep only the samples
where it is present (`mixture_masks[:, l] == 1`), sort those by `c_l`, split into `n_c_bins`
equal-count bins, and average:

```
H(A | C_l, l present) ≈ Σ_k (n_k/n_present) · H(A | bin_k)
I(A ; C_l | l present)  = H(A | l present) − H(A | C_l, l present)
```

Ligands present in `< 2·n_c_bins` samples are skipped; the metric is the mean over scored
ligands. The present-only restriction removes the presence/absence confound, so this measures
genuine concentration (level) coding — unlike the old version, which iterated over the sparse
present-*slot* columns (arbitrary order + padding zeros) and was meaningless. This is still a
per-ligand MARGINAL average, NOT on the joint `full_array_entropy` scale.

**`identity_channel` / `concentration_channel` — joint, total-comparable:**
Conditioning on the *full* presence pattern M (samples grouped by identical mask via
`torch.unique`) gives the exact chain-rule split `H(A) = I(A;M) + H(A|M)`: `identity_channel`
is the composition channel I(A;M), `concentration_channel` is the concentration residual
H(A|M) = I(A;c|M), and they sum to `full_array_entropy`. Use the reliance fractions
`identity_channel / full_array_entropy` and `concentration_channel / full_array_entropy`
(which sum to 1) to compare how much an architecture leans on identity vs concentration coding.
Reliable only when patterns repeat (single-ligand / low-`mu_ligands_per_source` sniffs, where M
is the categorical ligand id); in dense mixtures rows are nearly all unique → H(A|M) collapses
to 0 and I(A;M) → H(A) spuriously.

**Miller-Madow correction:** the plug-in entropy estimator is biased downward for finite
batch sizes. The Miller-Madow correction adds `(K_hat − 1) / (2·B·ln2)` where K_hat is
the number of distinct observed codewords, partially correcting this bias.

**Post-hoc measurement (outside the pipeline):** `analysis_helper` exposes reusable
primitives for re-measuring a **reloaded** checkpoint (via `plotlib.load_model`) without
retraining — `sample_activity(env, physics, ri, n_samples)` (estimator-agnostic, chunked
forward so memory is bounded), `kt_bracket(...)` (its KT lower/upper convenience), and
`eval_batch_cap(free_bytes)` (the `(tile, B)` memory cap). `src.testscaling` orchestrates on
top of these (walk sweeps → pick runs → choose test sizes → write `<sweep>/test_scaling.csv`);
the per-task `tasks/*/scripts/test_scaling.py` are thin wrappers (`build_parser` +
`sizes_from_args` + `run`) that just set defaults. Size strategy: `--test_sizes` (absolute),
`--mult` (per-run multiples of the train batch B, e.g. `1 2 4 8 16`), or an auto ×4 ladder —
to check the entropy-vs-samples curve has plateaued.

---

## Sweep Architecture (`config.py::RunConfig + run.py::SweepRunner`)

`RunConfig` accepts scalar or list values for every parameter field.
List-valued fields are **zip-iterated** (not crossed): all axis lists must share
the same length L, producing exactly L steps in a single trajectory.

Fields whose values are inherently arrays (`conc_mean`, `conc_std`,
`kernel_params`, `measurement_fns`, `cell_gene_probs`, `cell_size_pmf`) use a
`tuple` when fixed and a `List[tuple]` when iterated, so `isinstance(val, list)`
uniformly identifies every axis without special-casing. `cell_gene_sets` is nested
one level deeper (a tuple of gene tuples) and so has its own round-trip handling
(`_NESTED_TUPLE_FIELDS`).

Concentration parameters are supplied **directly** as tuples in the config (no
range-based RNG sampling).  To obtain multiple statistically independent runs,
supply multiple entries in the list or run the sweep script multiple times.

**`warm_start: bool` (replaces `warm_start_axis`):**  
When `True`, steps are sorted by `(env_group, n_genes)` — env groups are inferred
from the input order (a new group begins each time n_genes does not strictly
increase, i.e. the sweep restarts).  Within each group steps are sorted by
n_genes ascending and chain warm-started; at group boundaries n_genes decreases,
which triggers a cold reset.  When `False`, steps run in natural order and every
step starts cold.

**Folder layout** (written by `IO.py::_run_rel_path`):
```
{sweep_root}/{scalar_axis_1}_{val}/.../run_{YYYYMMDD_HHMMSS}/
```
Only scalar-valued axes appear as directory components (array-typed axes like
`conc_mean` are recovered from `config.json`).  The timestamp leaf guarantees
uniqueness when identical parameters are run more than once.  The execution
timestamp is also stored as `run_timestamp` in `config.json`.

**`SweepLoader.iter_run_dirs()`** crawls the sweep root for `config.json` files
rather than regenerating paths from the sweep config, making it robust to
interrupted or partial sweeps.

**Fields for heteromer sweeps:**

| Field | Type | Meaning |
|---|---|---|
| `n_receptors` | `Optional[int]` | Target receptor count; triggers `build_heteromer_array` in `__post_init__` |
| `receptor_sampling_strategy` | `str` | `"cascading"` (default) or `"uniform_random"` |
| `receptor_sampling_seed` | `Optional[int]` | RNG seed; same args → same receptor set |
| `use_interface_model` | `bool` | Forwarded to `build_heteromer_array`; routes to ordered-ring variants when `True` |

**Fields for cell sweeps** (§3b). Setting `cell_gene_sets` **or** `n_cells` switches
the array from receptors to cells; `receptor_indices` is then derived and
`n_receptors` / `receptor_sampling_*` are ignored.

| Field | Type | Meaning |
|---|---|---|
| `cell_gene_sets` | `Optional[Tuple[Tuple[int,...],...]]` | Explicit gene set per cell |
| `n_cells` | `Optional[int]` | Sample this many gene sets instead |
| `cell_sampling_strategy` | `str` | `"bernoulli"` (default) or `"size_pmf"` |
| `cell_gene_probs` | `Optional[Tuple[float,...]]` | bernoulli: P(gene *u* expressed); default `2/n_genes` each |
| `cell_size_pmf` | `Optional[Tuple[float,...]]` | size_pmf: entry *i* = P(cell expresses *i*+1 genes) |
| `cell_max_genes` | `Optional[int]` | bernoulli: reject cells above this (caps $R_\text{pool}$) |
| `cell_sampling_seed` | `Optional[int]` | RNG seed; same args → same gene sets |
| `cell_stoichiometry` | `str` | `"multinomial"` (default, random assembly) or `"uniform"` |
| `cell_readout` | `str` | `"threshold"` (default), `"noisy_or"`, `"mean"` |
| `cell_threshold` | `float \| "auto"` | `"auto"` → per-cell median drive at calibration |
| `cell_temperature` | `float` | $T_\text{cell}$ at the end of annealing |
| `cell_initial_temperature` | `float \| "auto"` | `"auto"` → std of the calibrated drive |
| `cell_pool_chunk` | `Optional[int]` | Receptors per pool chunk in the fused physics+readout pass |

**Batch-size auto-scaling** (`run.py::resolve_batch_sizes`): pass `batch_size="auto"` and/or
`test_batch_size="auto"` in `SingleRunConfig` / `RunConfig` to have sizes resolved at init
time based on array size R and entropy estimator:

- **Shannon**: `B_train = max(512, 2^R)` — one sample per histogram bin for good coverage.
  Memory cap: the soft-assignment tensor has shape `(B, 2^R)` float32; budget is
  `B × 2^R ≤ 2^35` floats (~128 GiB), yielding `B_max = 2^(35−R)` (~10^6 at R = 15 on A100).
  The cap binds for R ≥ 18, gracefully reducing B back toward the minimum.
- **Collision**: cost scales as O(m²·R) per chunk, not O(B·2^R).
  `B_train = max(512, 16 · 2^(R/2))`, rounded up to a multiple of the adaptive collision
  chunk size `m_max`. `m_max = max(512, floor(√(budget / (R · 4 · 3))))` — the (R, m, m)
  binding matrix in float32 with 3× safety for backward. The full batch is consumed via
  `ceil(B / m_max)` cross-chunk pairs (no `max_chunks` cap). `resolve_batch_sizes` returns
  `collision_chunk_size = m_max` alongside the batch sizes; it is threaded to
  `DiscreteExactLoss` and on to `compute_collision_entropy`.
- **KT (Kolchinsky–Tracey)**: the double loop retains every `(m,n,R)` block for backward,
  so the retained graph is `B²·R·4` **independent of chunk** — the batch is a single tile at
  `B = floor(√(mem_budget / (R·4·3)))`, with `collision_chunk_size = min(B, KT_TRAIN_TILE)`
  (`KT_TRAIN_TILE = 4096`) only to keep the per-tile `torch.log` transient bounded. The tile
  size does **not** change the retained `B²·R` (all tiles are kept for backward); a bigger tile
  just means fewer, larger (fused) kernels → fewer launches, for an `m²·R` transient that rides
  under the already-present `B²·R`. `compile_kt=True` `torch.compile`s the per-tile kernel. Config flag **`recompute_backward=True`**
  gradient-checkpoints each i-chunk's contribution (`_kt_row_contribution`) so the blocks are
  recomputed in backward, dropping the retained graph to `O(B·R)`. Memory then stops binding
  and the **compute** does: `B = floor(√(KT_COMPUTE_BUDGET / R))` (same √(1/R) shape),
  floored by the physics cap. `KT_COMPUTE_BUDGET` (module constant, `run.py`) defaults to 4×
  the memory-bound work `B²R ≈ 5.3e9` ⇒ ~2× batch, ~4× step time; exact value + gradient.
  See §06.5.
- **Blocked Shannon**: builds `(B, 2^block_size)` histograms (default `block_size=15` →
  `2^15`), *not* `(B, 2^R)` — so it is **not** subject to the Shannon `2^R` cap (that
  misclassification pinned B to the floor of 512). `ceil(R/block_size)·n_partitions` such
  histograms are retained for backward: `B_blocked = budget / (2^block_size · ceil(R/block_size)
  · n_partitions · 4 · 4)`. Independent of R, so the physics cap below typically binds.
  Config flag **`recompute_backward=True`** gradient-checkpoints the histogram (recomputed in
  backward, not retained): the trailing `· 4` safety drops to `· 2`, ~doubling B for ~+20% step
  time, exact gradients. Applies to `blocked` / `blocked_corrected` / `blocked_to_corrected`
  (all histogram-only); **not** `annealed`, whose un-checkpointed collision `(R,B,B)` block would
  then bind. See §06.6.
- **proxy / mi_***: O(B·R²) or O(B·R) — no exponential or `B²` term; physics-bound.
- **Physics bottleneck cap (all estimators)**: in the **interface model** the
  forward+backward holds *many* `(B, n_ligands, R·k_sub)` float32 tensors simultaneously —
  in `_compute_energies` (`ab`, `dist_sq`, `exp(·)`, `E_open`) and again in `p_open`
  (`log_terms_open/closed`), several retained for backward plus their gradients. The
  retained energy graph also coexists with the collision `(R,m,m)` chunk during loss/backward,
  so the factor must leave headroom for that term too. A **16× safety factor** is applied:
  `B_physics = free_mem × 0.8 / (16 · n_ligands · R · k_sub · 4)` (80% of free GPU memory,
  queried via `torch.cuda.mem_get_info()`). `B_train = min(B_train, B_physics)`.
  - The old 4× factor assumed a *single* such tensor and OOM'd for the interface model.
    It survived for the **classic model** only by accident: there the energy tensor is
    `(B, n_ligands, U)` with `U = n_genes` (no `k_sub` axis), so charging `R·k_sub` over-
    estimates by `k_sub` — a hidden ~5× cushion that vanishes once `use_interface_model=True`.
  - Note `s_upper ≤ n_ligands`, so `n_ligands` is a conservative bound on the present-ligand
    axis of the energy tensor. When the presence sampler is dense (e.g. `n_presence_blocks=1`,
    large `mu_ligands_per_source`) `s_upper` saturates to `n_ligands` and the bound is tight.
- **Measurement batches** (when `test_batch_size="auto"`): two sizes are returned.
  `test_perepoch = 4 · batch_size` drives the per-epoch convergence curve (cheap, so KT's
  O(B²) does not dominate every logged epoch); `test_final = min(2^R, memory)` is used once for
  the closing measurement (as many samples as fit). Config field **`test_max_batch`** caps
  `test_final` (`test_final = max(batch_size, min(test_final, test_max_batch))`) to bound the
  O(B²) KT cost of the final test at high R **without** shrinking training or the per-epoch curve
  An explicit `test_batch_size` (non-`"auto"`) is used for both. Evaluation is **chunked**
  (`_eval_stats`): soft metrics run on a single `chunk_size = train batch_size` pass; the batch
  drives multi-pass accumulation and the full-batch KT.
- **`per_epoch_measure=False`** (light mode): `_train` skips the per-epoch `_eval_stats` entirely and
  instead logs the training objective already computed for the gradient step (`loss`, and
  `train_entropy = -loss` — the native entropy for entropy-maximising losses) — **free**, no extra
  sampling. The single closing test then runs at **4×train** (regardless of `test_batch_size`/
  `test_max_batch`). Used by `fig1/het_casc.py`: the impact_of_heteromerization figure needs only the
  final KT bracket, and more samples can be obtained post-hoc from the saved `best_model.pt`.

---

## SQLite Run Index (`src/db.py`)

`runs.db` is a derived lookup table kept in `base_folder`.  It does **not** change
the folder structure or data files — `config.json` remains ground truth.  Delete
and rebuild at any time with `backfill`.

### Schema

| Column | Type | Source |
|---|---|---|
| `path` | TEXT PK | Relative path from `base_folder` |
| `sweep_name` | TEXT | Prefix of sweep dir name (before `_YYYYMMDD_HHMMSS`) |
| `sweep_date` | TEXT | Timestamp regex `(\d{8}_\d{6})` from sweep dir |
| `status` | TEXT | `complete` / `partial` / `missing` |
| `run_mtime` | REAL | `os.path.getmtime` of run dir |
| `git_hash` | TEXT | Short HEAD hash at index time, or NULL |
| `created` / `modified` | TEXT | ISO timestamps (UTC); `created` is immutable after first insert |
| *(scalar config fields)* | varies | All scalar `SingleRunConfig` fields; list-valued fields (`conc_mean`, `conc_std`, `kernel_params`, `measurement_fns`, `receptor_indices`) are skipped → NULL |
| `{metric}_mean` | REAL | Mean of each `test_results.json` metric list; columns added dynamically via `ALTER TABLE` when new metrics appear |

UNIQUE key on `path`; `INSERT OR REPLACE` (UPSERT) makes re-indexing idempotent.
WAL mode + exponential-backoff retry handle parallel sweeps writing simultaneously.

### Automatic hook

`SweepRunner.execute()` calls `db.add_run(run_dir, db_path)` after every
`test_results.json` is written.  The call is wrapped in a bare `try/except` — a DB
failure never aborts a sweep.  The DB is not created automatically; run `init` once
to activate indexing.

### SweepLoader DB integration

`SweepLoader` auto-detects `runs.db` by looking one directory above `sweep_root`
(i.e. at `base_folder/runs.db`).  When present it is used transparently:

| Method | DB present | DB absent |
|---|---|---|
| `load_all_test_results()` | single SQL query (`GLOB sweep_subdir/*`) | disk crawl (config.json + test_results.json) |
| `find_run_dir(**filters)` | SQL `WHERE` on scalar config cols + GLOB | `iter_run_dirs()` + `getattr` filter |
| `iter_run_dirs()` | always disk crawl (needs full config) | disk crawl |
| `load_all_histories()` | always disk crawl (stats.csv not in DB) | disk crawl |

`_mean`-suffixed metric columns from the DB are renamed to bare metric names
(e.g. `full_array_entropy_mean` → `full_array_entropy`) so the DataFrame
schema matches the disk-crawl convention and analysis code is unchanged.

### CLI

```
python -m src.db init       runs.db
python -m src.db backfill   runs.db
python -m src.db add-run    runs.db  path/to/run_dir
python -m src.db sync       runs.db
python -m src.db reconcile  runs.db  [--dry-run]
python -m src.db delete     runs.db  relative/path  [--dry-run]
python -m src.db move       runs.db  old/path  new/path
python -m src.db alter      runs.db  add-col    col_name  TYPE
python -m src.db alter      runs.db  remove-col col_name  [--dry-run]
python -m src.db query      runs.db  [--where EXPR] [--cols c1,c2] [--limit N]
```

`--dry-run` is supported on `reconcile`, `delete`, and `alter remove-col`.

---

## Receptor Sampling Strategies (`geometry.py`)

### Unified entry-point: `build_heteromer_array`

```python
build_heteromer_array(n_genes, k_sub, R_target, strategy="uniform_random", seed=None,
                      use_interface_model=False)
```

Returns `(R_target, k_sub)` long tensor. Called automatically from
`SingleRunConfig.__post_init__` when `n_receptors` is set and `receptor_indices` is None,
with `use_interface_model` forwarded from the config.

| `strategy` | Classic model | Interface model |
|---|---|---|
| `"cascading"` | `generate_cascading_receptors` — homomers first, then 2-mers, … | `generate_cascading_ordered_receptors` — same tier order, cyclic pool |
| `"uniform_random"` | Reservoir-sampling from `combinations_with_replacement` | `generate_ordered_receptor_indices` — sample from canonical cyclic pool |

When `R_target` exceeds the pool size the full pool is returned with a warning.
Determinism contract: same `(n_genes, k_sub, R_target, strategy, seed, use_interface_model)` → identical tensor.

### Pool size: `count_receptor_combinations`

```python
count_receptor_combinations(n_genes, k_sub, use_interface_model=False) -> int
```

Returns the number of distinct receptor types without enumerating them.

- **Classic model** (unordered multisets): $\binom{n\_genes + k\_sub - 1}{k\_sub}$ (stars and bars).
- **Interface model** (cyclic arrangements, rotation-identified, reflection-distinct):
  Burnside's lemma for the cyclic group $C_{k_{sub}}$:

$$\frac{1}{k_{sub}} \sum_{d \mid k_{sub}} \varphi(d)\, n_{genes}^{k_{sub}/d}$$

where $\varphi$ is Euler's totient function. Example: $n_{genes}=3,\, k_{sub}=5$ gives 21 (classic) vs 51 (interface).

### Classic model — unordered compositions (lower-level functions)

| Function | Strategy |
|---|---|
| `generate_receptor_indices` | Random sample from all combinations_with_replacement |
| `generate_targeted_receptors` | Explicit counts per complexity level (n_unique_subunits) |
| `generate_cascading_receptors` | Fill quota by complexity: homomers first, then 2-mers, etc.  Accepts optional `seed` arg. |
| `generate_exp_distributed_receptors` | Draw complexity from exponential distribution |
| `generate_bernoulli_receptors` | Each gene present via Bernoulli(gene_probs); fill slots proportionally |

### Interface model — ordered cyclic arrangements

The interface model requires **ordered** receptor tuples representing the ring layout.
Two tuples are the same receptor iff they are **cyclic rotations** of each other;
reflections are **distinct** because the +/− face asymmetry breaks mirror symmetry.

Canonical representative: the **lexicographically minimum cyclic rotation** of the tuple
(computed by `_canonical_rotation`).

| Function | Strategy |
|---|---|
| `generate_ordered_receptor_indices` | Random sample from all canonical cyclic arrangements; accepts optional `seed` |
| `generate_targeted_ordered_receptors` | Explicit counts per complexity level; accepts optional `seed` |
| `generate_cascading_ordered_receptors` | Fill quota by complexity: homomers first; **lazy per-tier** (see below); accepts optional `seed` |

**Cost note:** the random/targeted variants enumerate the full pool: at most
`|combos| × k_sub!` candidate permutations, then deduplicate via canonical form.
For n_genes=26, k_sub=5 this is ~17 M operations — but it scales steeply with
n_genes (n_genes=35 → ~85 s).

**Cascading is lazy.** `generate_cascading_ordered_receptors` builds one
complexity tier at a time via `_tier_canonical_forms(n_genes, k_sub, n_unique)`
(genes chosen by `C(n_genes, n_unique)`, stoichiometries by
`_positive_compositions`), ascending from homomers, and stops as soon as the
`n_sensors` quota is filled. High-complexity tiers (4-mers, 5-mers) — which
dominate the pool — are never enumerated when `n_sensors` is small. For
n_genes=35, k_sub=5, n_sensors≤49 only tiers 1–2 are built (~50 ms vs ~85 s for
full enumeration). The seeded intra-tier shuffle is applied to a `sorted()`
tier list, so output is reproducible across runs/platforms.
