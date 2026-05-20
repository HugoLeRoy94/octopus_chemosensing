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

Each training step generates a fresh batch of `batch_size` sensory "sniffs".
Steps 1–3 are shared between the two model modes; step 4 diverges.

1. **Mixture mask** — independent Bernoulli draws per ligand: `M_ℓ ~ Bernoulli(p_presence_ℓ)`.
2. **Concentrations** — sampled from `ConcentrationModel`, then zero-masked for absent ligands.
3. **Observation noise** — `v_obs,ℓ = v_ℓ + N(0, σ_noise)` simulates docking variance.
4. **Distances & energies** — the identity trick avoids allocating an (B, L, *, D) tensor:

```python
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

Three loss modules exist, selected by `cfg.entropy` in `run.py::_build_loss`.

### 4a. `DiscreteExactLoss` (`bin_loss.py`) — maximize joint array entropy

**Objective:** `min −H(A)` where A is the binary array activity.

Four entropy estimators — pick based on array size:

| `entropy_type` | Complexity | When to use |
|---|---|---|
| `'shannon'` | O(B · 2^R) | R < ~15; exact but exponential |
| `'renyi'` | O(B² · R) | Default scalable choice; exact Rényi H2 |
| `'blocked'` | O(B · 2^block_size · R/block) | Captures higher-order terms, R up to ~100 |
| `'proxy'` | O(B · R²) | Fastest; pairwise covariance/repulsion penalty |

**Rényi trick:** avoids both the 2^R state space and float underflow.
Computes log P(collision) = Σ_r log P_r(collision), sums in log-space via logsumexp,
then exponentiates only once at the end. Diagonal (self-collision) entries are masked out
before averaging. For large batches (B > 2048), cross-chunk evaluation across 8 random
sub-batches avoids allocating a (B, B) matrix.

**Blocked trick:** randomly partitions the R receptors into blocks of `block_size=12`,
computes exact Shannon entropy within each block, and sums under the independence assumption.
Averaging over `n_partitions=4` random partitions reduces bias from any particular blocking.

### 4b. `MaximizeMutualInformationLigandLoss` (`family_mi_loss.py`)

**Objective:** `max I(A ; M) = H(A) − H(A | M)`

H(A) is computed on the full batch. H(A | M) conditions on the exact mixture identity:
mixture masks are hashed to integer IDs using binary powers
`id = Σ_ℓ M_ℓ · 2^ℓ`, then the batch is grouped by ID and entropy computed per group.

### 4c. `MaximizeMutualInformationConcentrationLoss` (`concentration_mi_loss.py`)

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
    4. compute loss         loss = loss_fn(activity, ...)
    5. backward + step      loss.backward(); optimizer.step()
    6. eval every 1%        _eval_stats(...)
```

**Temperature annealing:** starts at `T_init` (calibrated, ~std of pre-sigmoid terms)
and decreases linearly to the configured `temperature`. A high T keeps the sigmoid soft
early in training (smooth gradients); a low T sharpens to binary decisions later.

**Warm-starting** (`SweepRunner`): the environment from step N is passed to step N+1.
`warm_start_axis` in `RunConfig` controls which axis (or axes) forms the trajectory:

- **Single string** (e.g. `"n_genes"`): existing behaviour — sorted ascending, environment
  warm-started via `clone_with_extra_units` when `n_genes` grows.
- **List of strings** (e.g. `["n_genes", "n_receptors"]`): the listed axes are **zipped**
  together rather than forming a Cartesian product, sorted by the first axis.  Use this to
  grow n_genes and n_receptors jointly (the heteromer case).

`_initialize` always builds `receptor_indices` from `self.config.receptor_indices`, which
is auto-generated in `SingleRunConfig.__post_init__` from `(n_genes, n_receptors,
receptor_sampling_strategy, receptor_sampling_seed)` when `receptor_indices is None`.
The LR is damped 10× on warm-start to preserve learned representations.

**Optional cosine LR scheduler:** wraps Adam with `CosineAnnealingLR` when
`use_scheduler=True`, annealing from `lr` to `1e-5` over the full run.

---

## Stage 6 — Evaluation & Metrics (`run.py::SimulationRunner._eval_stats`)

Evaluation runs under `torch.no_grad()` with `test_batch_size` samples.
When `eval_chunk_size < test_batch_size`, soft metrics use one chunk and hard codeword
metrics are accumulated across all chunks (avoiding CUDA OOM on large eval budgets).

Available metrics (add to `measurement_fns` in config):

| Key | Function | What it measures |
|---|---|---|
| `full_array_entropy` | `analysis_helper.full_array_entropy` | Rényi H2 + blocked Shannon of joint activity |
| `codeword_entropy` | `analysis_helper.codeword_entropy` | Hard plug-in + Miller-Madow entropy of binary codewords |
| `mean_receptor_distance` | `analysis_helper.mean_receptor_distance` | Average pairwise latent-space distance between receptors |
| `receptor_distances` | `analysis_helper.receptor_distances` | Full (R, R) pairwise distance matrix |
| `conditional_entropy_ligand` | `analysis_helper.conditional_entropy_ligand` | H(A \| M) — conditioned on exact ligand identity |
| `mutual_information_ligand` | `analysis_helper.mutual_information_ligand` | I(A ; M) — identity channel |
| `conditional_entropy_concentration` | `analysis_helper.conditional_entropy_concentration` | H(A \| C) |
| `mutual_information_concentration` | `analysis_helper.mutual_information_concentration` | I(A ; C) — concentration channel |
| `conditional_entropy_family` | `analysis_helper.conditional_entropy_family` | H(A \| F) — conditioned on family presence mask |
| `mutual_information_family` | `analysis_helper.mutual_information_family` | I(A ; F) — family channel (coarser than identity) |
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

**Miller-Madow correction:** the plug-in entropy estimator is biased downward for finite
batch sizes. The Miller-Madow correction adds `(K_hat − 1) / (2·B·ln2)` where K_hat is
the number of distinct observed codewords, partially correcting this bias.

---

## Sweep Architecture (`config.py::RunConfig + run.py::SweepRunner`)

`RunConfig` accepts `Union[T, List[T]]` for any parameter. List-valued fields
become sweep axes; a Cartesian product is taken over all independent axes.
The warm-start axis (or axes) is extracted and run sequentially within each trajectory.

Per-trajectory concentration draws (`conc_mean`, `conc_std`, `p_presence`) are
deterministically sampled from a seeded NumPy RNG, so sweeps are fully reproducible
and re-loadable from the saved config.

**New fields for heteromer sweeps:**

| Field | Type | Meaning |
|---|---|---|
| `n_receptors` | `Optional[int]` | Target receptor count; triggers `build_heteromer_array` in `__post_init__` |
| `receptor_sampling_strategy` | `str` | `"cascading"` (default) or `"uniform_random"` |
| `receptor_sampling_seed` | `Optional[int]` | RNG seed; same args → same receptor set |

**Batch-size auto-scaling** (`run.py::resolve_batch_sizes`): pass `batch_size="auto"` and/or
`test_batch_size="auto"` in `SingleRunConfig` / `RunConfig` to have sizes resolved at init
time based on array size R and entropy estimator:

- **Shannon**: `B_train = max(512, 2^R)` — one sample per histogram bin for good coverage.
  Memory cap: the soft-assignment tensor has shape `(B, 2^R)` float32; budget is
  `B × 2^R ≤ 2^35` floats (~128 GiB), yielding `B_max = 2^(35−R)` (~10^6 at R = 15 on A100).
  The cap binds for R ≥ 18, gracefully reducing B back toward the minimum.
- **Rényi-2**: cost scales as O(B²·R), not O(B·2^R), so a smaller B suffices.
  `B_train = max(512, 16 · 2^(R/2))`.
- `test_batch_size = 4 · batch_size` in both cases.

---

## Receptor Sampling Strategies (`geometry.py`)

### Unified entry-point: `build_heteromer_array`

```python
build_heteromer_array(n_genes, k_sub, R_target, strategy="uniform_random", seed=None)
```

Returns `(R_target, k_sub)` long tensor. Called automatically from
`SingleRunConfig.__post_init__` when `n_receptors` is set and `receptor_indices` is None.

| `strategy` | Description |
|---|---|
| `"cascading"` | Homomers first, then 2-mers, … until quota reached.  Biologically motivated — simpler complexes fold more reliably. |
| `"uniform_random"` | Reservoir-sampling from the full `combinations_with_replacement` pool (O(R_target) memory; never materialises the full pool). |

When `R_target` exceeds the pool size the full pool is returned with a warning.
Determinism contract: same `(n_genes, k_sub, R_target, strategy, seed)` → identical tensor.

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
| `generate_ordered_receptor_indices` | Random sample from all canonical cyclic arrangements |
| `generate_targeted_ordered_receptors` | Explicit counts per complexity level |
| `generate_cascading_ordered_receptors` | Fill quota by complexity: homomers first |

**Cost note:** enumeration generates at most `|combos| × k_sub!` candidate permutations,
then deduplicates via canonical form.  For n_genes=26, k_sub=5 this is ~17 M operations
— acceptable at init time but not on the hot path.
