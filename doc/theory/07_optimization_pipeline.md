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

## Stage 4 — Loss Computation

Loss modules are selected by `cfg.entropy` in `run.py::_build_loss`.

### 4a. `DiscreteExactLoss` (`bin_loss.py`) — maximize joint array entropy

**Objective:** `min −H(A)` where A is the binary array activity (or `min C` for collision).

Five entropy estimators — pick based on array size:

| `entropy_type` | Complexity | When to use |
|---|---|---|
| `'shannon'` | O(B · 2^R) | R < ~15; exact but exponential |
| `'collision'` | O(B² · R) | Scalable lower bound; training loss minimises C directly |
| `'blocked'` | O(B · 2^block_size · R/block) | Upper bound, captures within-block correlations |
| `'blocked_corrected'` | O(B · 2^block_size + R²) | Tighter: blocked minus cross-block pairwise MI |
| `'proxy'` | O(B · R²) | Fastest; pairwise covariance/repulsion penalty |

**Collision trick:** computes log P(collision) = Sigma_r log P_r(collision) in log-space
via logsumexp, then exponentiates once. Training loss returns C = exp(log_mean_coll_prob)
directly (no log), removing the 1/C gradient blow-up at low collision probability.
Measurement returns H2 = -log2(C) in bits. For B > 2048, cross-chunk evaluation across
8 sub-batches avoids the (B, B) matrix.

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
    1. anneal temperature   T = T_end + (T_start − T_end) * (1 − epoch/epochs)
    2. sample batch         E, concs, masks = env.sample_batch(batch_size)
    3. compute activity     activity = physics(E, concs, receptor_indices)
    4. compute loss         loss = loss_fn(activity, ...)     # see dispatch below
    5. backward + step      loss.backward(); optimizer.step()
    6. eval every 1%        _eval_stats(...)
```

**Loss dispatch** in step 4 depends on the loss module type:
- `DiscreteExactLoss`: `loss_fn(activity)`
- `AnnealedEntropyLoss`: `loss_fn(activity, epoch, self.config.epochs)` — needs
  training progress to compute the interpolation parameter λ.
- `MaximizeMutualInformationLigandLoss`: `loss_fn(activity, mixture_masks=masks)`
- `MaximizeMutualInformationConcentrationLoss`: `loss_fn(activity, concs=concs)`

**Temperature annealing:** starts at `T_init` and decreases linearly to the configured
`temperature`. A high T keeps the sigmoid soft early in training (smooth gradients); a
low T sharpens to binary decisions later. `T_init` is controlled by
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
| `full_array_entropy` | `analysis_helper.full_array_entropy` | Collision H2 + blocked + blocked-corrected Shannon of joint activity |
| `codeword_entropy` | `analysis_helper.codeword_entropy` | Hard plug-in + Miller-Madow entropy of binary codewords |
| `mean_receptor_distance` | `analysis_helper.mean_receptor_distance` | Average pairwise latent-space distance between receptors |
| `receptor_distances` | `analysis_helper.receptor_distances` | Full (R, R) pairwise distance matrix |
| `conditional_entropy_ligand` | `analysis_helper.conditional_entropy_ligand` | (1/N_l) Σ H(A \| L_l) |
| `mutual_information_ligand` | `analysis_helper.mutual_information_ligand` | (1/N_l) Σ I(A ; L_l) |
| `conditional_entropy_concentration` | `analysis_helper.conditional_entropy_concentration` | (1/N_l) Σ H(A \| C_l) |
| `mutual_information_concentration` | `analysis_helper.mutual_information_concentration` | (1/N_l) Σ I(A ; C_l) |
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

**`conditional_entropy_concentration` / `mutual_information_concentration` — per-ligand marginal, quantile binning:**
Each ligand l is treated independently. Samples are sorted by c_l and split into
`n_c_bins` equal-quantile bins. The conditional entropy for ligand l is:

```
H(A | C_l) ≈ Σ_k (n_k/B) · H(A | C_l ∈ bin_k)
```

The function returns `(1/N_l) Σ_l H(A | C_l)`, so `mutual_information_concentration`
returns `(1/N_l) Σ_l I(A ; C_l)` — the mean pairwise MI between the receptor array
and each individual ligand's concentration.

**Miller-Madow correction:** the plug-in entropy estimator is biased downward for finite
batch sizes. The Miller-Madow correction adds `(K_hat − 1) / (2·B·ln2)` where K_hat is
the number of distinct observed codewords, partially correcting this bias.

---

## Sweep Architecture (`config.py::RunConfig + run.py::SweepRunner`)

`RunConfig` accepts scalar or list values for every parameter field.
List-valued fields are **zip-iterated** (not crossed): all axis lists must share
the same length L, producing exactly L steps in a single trajectory.

Fields whose values are inherently arrays (`conc_mean`, `conc_std`,
`kernel_params`, `measurement_fns`) use a `tuple` when fixed and a `List[tuple]`
when iterated, so `isinstance(val, list)` uniformly identifies every axis without
special-casing.

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

**Batch-size auto-scaling** (`run.py::resolve_batch_sizes`): pass `batch_size="auto"` and/or
`test_batch_size="auto"` in `SingleRunConfig` / `RunConfig` to have sizes resolved at init
time based on array size R and entropy estimator:

- **Shannon**: `B_train = max(512, 2^R)` — one sample per histogram bin for good coverage.
  Memory cap: the soft-assignment tensor has shape `(B, 2^R)` float32; budget is
  `B × 2^R ≤ 2^35` floats (~128 GiB), yielding `B_max = 2^(35−R)` (~10^6 at R = 15 on A100).
  The cap binds for R ≥ 18, gracefully reducing B back toward the minimum.
- **Collision**: cost scales as O(B²·R), not O(B·2^R), so a smaller B suffices.
  `B_train = max(512, 16 · 2^(R/2))`. This heuristic explodes for large R, so it is also
  capped by the **collision memory bound**: the estimator materialises a `(B, B)` collision
  matrix (`B²·4` bytes), giving `B_collision = √(budget / (4 · 4))` (4× safety). Without this
  cap the `(B,B)` term is bounded only by the physics cap below, which ignores it.
- **Blocked Shannon**: builds `(B, 2^block_size)` histograms (default `block_size=15` →
  `2^15`), *not* `(B, 2^R)` — so it is **not** subject to the Shannon `2^R` cap (that
  misclassification pinned B to the floor of 512). `ceil(R/block_size)·n_partitions` such
  histograms are retained for backward: `B_blocked = budget / (2^block_size · ceil(R/block_size)
  · n_partitions · 4 · 4)`. Independent of R, so the physics cap below typically binds.
- **proxy / mi_***: O(B·R²) or O(B·R) — no exponential or `B²` term; physics-bound.
- **Physics bottleneck cap (all estimators)**: in the **interface model** the
  forward+backward holds *many* `(B, n_ligands, R·k_sub)` float32 tensors simultaneously —
  in `_compute_energies` (`ab`, `dist_sq`, `exp(·)`, `E_open`) and again in `p_open`
  (`log_terms_open/closed`), several retained for backward plus their gradients. The
  retained energy graph also coexists with the collision `(B,B)` matrix during loss/backward,
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
- `test_batch_size = 4 · batch_size` in all cases. Safe despite the larger value because
  evaluation is **chunked** (`_eval_stats`): soft metrics (incl. collision `(B,B)`) run on a
  single `chunk_size = train batch_size` pass; `test_batch_size` only drives multi-pass
  accumulation of hard-codeword metrics.

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
