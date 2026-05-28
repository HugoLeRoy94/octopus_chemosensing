# Implementation Spec — Block-Correlated Ligand Presence via Gaussian Copula

## Goal

Replace the independent per-ligand Bernoulli presence draw in the environment with a
**Gaussian-copula** draw that introduces *block-structured correlation* between
ligands ("sources" that emit co-occurring bouquets), while leaving each ligand's
marginal presence probability exactly unchanged.

The feature must be **opt-in** and must **reduce exactly to the current behaviour**
when correlation is off, so existing experiments are untouched.

Reference: del Castillo et al., PNAS 2026, §1.1 — they use a block-covariance
Gaussian copula for odorant presence. We follow the same construction.

---

## Background — what a Gaussian copula does (for the implementer)

To draw a correlated binary presence vector `M ∈ {0,1}^L` with (a) per-ligand
marginal probabilities `p_ℓ` and (b) block-structured correlation:

1. Build an `L × L` correlation matrix `Σ` — block-diagonal. Within a block,
   off-diagonal entries are `rho_block`; across blocks, 0; diagonal 1.
2. Cholesky-factor once: `Σ = R^T R`.
3. Per batch: draw independent standard normals `eps` of shape `(B, L)`,
   correlate them `z = eps @ R`, so each row ~ `N(0, Σ)`.
4. Threshold: `M[b,ℓ] = 1 if z[b,ℓ] < tau_ℓ else 0`, where
   `tau_ℓ = Phi_inv(p_ℓ)` is the standard-normal quantile of the target marginal.

The marginal is hit exactly because `Phi(z_ℓ)` is uniform; the correlation in `Σ`
survives the monotone threshold. `rho_block = 0` → `Σ = I` → independent
Bernoulli, i.e. the current behaviour.

---

## Change 1 — Block assignment at environment construction

In `LigandEnvironment.__init__` (`environment.py`):

- Add config-driven partition of the `n_ligands` ligands into `n_presence_blocks`
  source blocks.
- **Block membership MUST be assigned independently of `ligand_family_assignments`
  and independently of `ligand_latent` positions.** Presence blocks (co-occurrence /
  "sources") and latent families (chemical similarity) are orthogonal structures.
  Do NOT reuse the family partition. Assign blocks by a separate seeded shuffle.
- Store block membership as a buffer (e.g. `presence_block_id`, shape `(n_ligands,)`).
- Block sizes: by default partition as evenly as possible (`n_ligands //
  n_presence_blocks`, distributing the remainder). Keep it simple.

## Change 2 — Build and store the copula factor

Still in `__init__`, when correlated presence is enabled:

- Construct the `L × L` correlation matrix `Sigma`:
  - diagonal = 1
  - entry `(i,j)` for `i != j` = `rho_block` if ligands `i,j` share a block, else 0
- Numerical safety: a block matrix with constant off-diagonal `rho` is positive
  definite for `rho ∈ (−1/(m−1), 1)` where `m` is block size; since we use
  `rho ∈ [0,1)` and `m ≥ 1` this is fine, but still add a tiny jitter
  (`Sigma += 1e-6 * I`) before factorising to be safe.
- Cholesky-factor: `R = torch.linalg.cholesky(Sigma)` (store the factor as a
  buffer; do NOT recompute per batch).
- Precompute thresholds `tau = Phi_inv(p_presence)` once and store as a buffer.
  Use `torch.distributions.Normal(0,1).icdf(p_presence)` or
  `math.sqrt(2) * torch.erfinv(2*p - 1)`.

## Change 3 — Correlated draw in `sample_batch`

In `LigandEnvironment.sample_batch` (`environment.py`):

- Current code draws `M_ℓ ~ Bernoulli(p_presence_ℓ)` independently.
- Replace with: if correlated presence is enabled,
  ```
  eps = torch.randn(batch_size, L, device=...)
  z   = eps @ R.T            # (B, L), each row ~ N(0, Sigma)
  M   = (z < tau).float()    # (B, L) multi-hot mask, broadcasting tau over batch
  ```
  Verify the matrix orientation against the stored Cholesky convention
  (`torch.linalg.cholesky` returns lower-triangular `L` with `Sigma = L @ L.T`;
  then `z = eps @ L.T` gives `cov(z) = L L.T = Sigma`). Add a one-line comment
  fixing the convention so it can't silently transpose-bug.
- If correlated presence is disabled, keep the existing independent Bernoulli path
  untouched.
- Everything downstream of `M` (concentration sampling, zero-masking absent
  ligands, energy computation) is UNCHANGED.

## Change 4 — Concentration: block-shared mean only

Optional but recommended, controlled by its own config flag:

- When block correlation is on, give each *block* a shared concentration mean
  (draw one `conc_mean` per block, then each ligand in the block inherits it).
- **Do NOT couple per-sample concentrations.** Per-sniff concentrations still draw
  independently around the (now block-shared) mean. Rationale (del Castillo et al.,
  §1.1): turbulent transport scrambles concentration ratios even when source
  co-occurrence survives — they preserve presence statistics but explicitly do not
  preserve concentration ratios across samples.
- If this flag is off, concentration draws are exactly as they are now.

## Change 5 — Config plumbing

In `config.py` (`SingleRunConfig` and `RunConfig`):

- New fields, no default value:
  - `n_presence_blocks: int`
  - `rho_block: float`
- `RunConfig` must accept `Union[T, List[T]]` for `n_presence_blocks` and
  `rho_block` so they become sweep axes.

## Change 6 — Determinism & docs

- All new random draws (block partition) must be seeded and reproducible:
  same `(n_ligands, n_presence_blocks)` → identical partition.
- The per-batch `torch.randn` follows the existing batch-sampling RNG convention —
  do not introduce a separate uncontrolled generator.
- Update the relevant docstring / theory-doc comment in `environment.py` to note
  the copula construction and cite del Castillo et al. 2026 §1.1.

---

## Validation checklist (please verify before declaring done)

2. **Marginals preserved under correlation:** with `rho_block = 0.6`, per-ligand
   empirical presence frequency still ≈ `p_presence` (the copula must not shift
   marginals — only correlation).
3. **Correlation appears in the right place:** with `rho_block = 0.6`, the
   empirical binary correlation matrix is block-diagonal — positive within blocks,
   ~0 across blocks. (Note: realized *binary* correlation is monotonically related
   to `rho_block` but smaller — the tetrachoric compression — this is expected, do
   not "fix" it.)
3. **Orthogonality:** confirm the block partition is uncorrelated with
   `ligand_family_assignments` (e.g. mutual information / contingency between the
   two partitions ≈ chance).
4. **PD safety:** Cholesky succeeds for `rho_block` up to at least 0.95 and for
   the largest block size used.
5. **Determinism:** two environments built with the same seeds have identical
   `presence_block_id` and identical `R` factor.

## Scope discipline

Keep changes minimal and localized to `environment.py` and `config.py`. Do not
refactor unrelated sampling code. Do not change the receptor physics, the loss
functions, or the measurement code. The independent-Bernoulli path must remain
the default and must remain a literal code path, not an emulated special case, if
that is cleaner — but the copula path with `rho=0` must also be verified to match.
