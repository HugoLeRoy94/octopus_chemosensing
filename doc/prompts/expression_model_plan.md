# Gene Expression Models — Implementation Plan

## Context

Data is loaded as in `individual_to_collectif_genes_expression.py` and
`correlation_matrix.py`:

- `A`: binary matrix of shape `(n_genes, n_cells)` from `binaryTable_allGenes`
- `CRnames`: gene names (length `n_genes`, typically 24)
- `prob_expr = A.mean(axis=1)`: per-gene expression probabilities `p_i`
- `coexp = A.sum(axis=0)`: number of genes expressed per cell
- `corr_matrix`, `co_expr`: pairwise correlation and co-expression matrices

The goal is to fit successively richer models of `P(σ)` (the joint distribution
over gene expression patterns) and check what each one predicts about the data.

Create a new script `expression_models.py` next to the existing files. Reuse
`_to_str_list` and the loading boilerplate.

---

## Step 1 — Independent baseline (Poisson-binomial)

**Model:** every gene independent, $P(\sigma) = \prod_i p_i^{\sigma_i}(1-p_i)^{1-\sigma_i}$.

**Implement:**

1. Compute `p_i = A.mean(axis=1)`.
2. Compute the Poisson-binomial distribution of $L = \sum_i \sigma_i$ exactly by
   expanding the generating function $G(z) = \prod_i [(1-p_i) + p_i z]$.
   Just iteratively convolve: start with `pmf = [1.0]`, then for each `p_i`
   do `pmf = np.convolve(pmf, [1-p_i, p_i])`. Returns array of length `n_genes+1`.
3. Compute three headline numbers:
   - `mu_data = coexp.mean()`
   - `var_data = coexp.var()`
   - `var_indep = np.sum(p_i * (1 - p_i))`
   Print them and the ratio `var_data / var_indep`. Expect ratio ≫ 1.
4. Refit the geometric/exponential properly (the current fit isn't normalized
   on integers). Use $P(L=l) = (1-q) q^{l-1}$ with $q = e^{-\beta}$; fit by
   linear regression of $\log P(L=l)$ vs $l$ on the nonzero bins.

**Output figures:**

- Single panel comparing three curves on log-y axis: empirical `P(L)`,
  Poisson-binomial prediction, geometric fit. Title should display the
  variance ratio.

---

## Step 2 — Bernoulli mixture (latent cell types)

**Model:** $P(\sigma) = \sum_k \pi_k \prod_i q_{i,k}^{\sigma_i}(1-q_{i,k})^{1-\sigma_i}$.

**Implement:**

1. Use `sklearn.mixture.BernoulliMixture` if available, otherwise hand-roll EM:
   - E-step: posterior $\gamma_{c,k} \propto \pi_k \prod_i q_{i,k}^{A_{i,c}}(1-q_{i,k})^{1-A_{i,c}}$,
     done in log-space.
   - M-step: $\pi_k = \langle \gamma_{c,k}\rangle_c$, $q_{i,k} = \sum_c \gamma_{c,k} A_{i,c} / \sum_c \gamma_{c,k}$.
   - Initialize $q_{i,k}$ with small random perturbations around $p_i$.
   - Stop on log-likelihood plateau (rel tol 1e-5) or 200 iterations.
2. Fit for `K = 1, 2, ..., 8`. For each, compute:
   - train log-likelihood
   - BIC = $-2\log L + (K \cdot n_\text{genes} + K - 1) \log n_\text{cells}$
   - 5-fold held-out log-likelihood per cell
3. Pick best `K` by held-out LL (BIC as sanity check).
4. For the chosen `K`, derive predicted `P(L)` by sampling $10^5$ patterns from
   the fitted mixture and binning, or analytically by convolving per-component
   Poisson-binomials weighted by $\pi_k$.

**Output figures:**

- BIC and held-out LL vs `K`.
- Heatmap of $q_{i,k}$ for the chosen `K` (genes on y-axis, components on x),
  reordered with the same `reordered_idx` as in `correlation_matrix.py` so blocks
  are visible.
- Overlay predicted `P(L)` from the mixture on the Step 1 figure.

---

## Step 3 — Pairwise MaxEnt (Ising)

**Model:** $P(\sigma) \propto \exp\left(\sum_i h_i \sigma_i + \sum_{i<j} J_{ij}\sigma_i\sigma_j\right)$.

**Implement:**

1. Fit by **pseudolikelihood** (don't try exact partition function at N=24 unless
   the run is fast — $2^{24} \approx 1.7 \times 10^7$ states is borderline; pseudolikelihood
   is robust and fast).
   - For each gene $i$, do logistic regression of $A_i$ on all other genes $A_{j\neq i}$.
     Coefficients give $h_i$ (intercept) and $J_{ij}$ (slope on $A_j$).
   - Symmetrize: $J_{ij} \leftarrow (J_{ij} + J_{ji})/2$.
   - Use `sklearn.linear_model.LogisticRegression` with no regularization or
     small L2 (`C=100` or so).
2. Sanity checks the model must pass (these are the actual tests):
   - Predicted $\langle \sigma_i \rangle$ matches data $p_i$. Compute by Gibbs sampling
     ($10^5$ samples after $10^4$ burn-in) and overlay scatter.
   - Predicted $\langle \sigma_i \sigma_j \rangle$ matches data. Scatter plot.
3. **The real test**: predicted `P(L)` from the Ising samples vs data.

**Output figures:**

- Heatmap of $J_{ij}$ reordered with `reordered_idx`. Compare visually to the
  correlation heatmap.
- Two scatter plots: predicted vs observed $\langle \sigma_i \rangle$ and
  $\langle \sigma_i \sigma_j \rangle$.
- Overlay predicted `P(L)` from Ising on the Step 1/2 figure (now four curves:
  data, Poisson-binomial, mixture, Ising).

---

## Step 4 — Triplet test (the falsification step)

**Implement:**

1. Compute observed triplet means $\langle \sigma_i \sigma_j \sigma_k \rangle$ for all triplets
   from data (only triplets where at least, say, 20 cells have all three on, to
   avoid noise — record the count too).
2. Compute the same triplets from Gibbs samples of the Ising model.
3. Compute the same triplets from samples of the Bernoulli mixture.
4. Compute the independent prediction $p_i p_j p_k$.

**Output figure:**

- Scatter plot, observed vs predicted triplet probability, one panel per model
  (independent, mixture, Ising). Add the y=x line. The model that hugs the
  diagonal best is the winner.

---

## Step 5 — Save artifacts

Persist these for downstream use:

- `p_i`, the Poisson-binomial pmf
- mixture parameters $\pi_k, q_{i,k}$ for the chosen `K`
- Ising parameters $h_i, J_{ij}$
- All `P(L)` predictions on the same grid

Use a single `np.savez` with descriptive keys. Print the headline numbers
(variance ratio, chosen K, mean abs error on $\langle\sigma_i\sigma_j\rangle$, triplet
correlation per model) at the end as a summary block.

---

## Notes / gotchas

- All sampling diagnostics should use the **same** evaluation batch size where
  comparable, otherwise the noise floors differ.
- Don't compute the Ising partition function explicitly. Pseudolikelihood for
  fitting, Gibbs sampling for everything else.
- If `n_genes` ≤ 20, you *can* enumerate $2^N$ states for exact Ising — feel free
  to add that as a check, but pseudolikelihood + Gibbs is the default.
- Keep the existing `reordered_idx` ordering throughout for visual continuity
  with the correlation matrix figure.
