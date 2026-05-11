# Topic 2 — Alternatives to the Quadratic Expansion from the Ideal Ligand

## Context

The current latent-space model (see `03_latent_environment.md`) writes the open-state energy as:

$$\tilde{E}_o^{(u,\ell)} = E_\text{base}^{(u)} + \|\textbf{v}_u - \textbf{v}_{\text{obs},\ell}\|^2$$

This is justified as a harmonic (quadratic) expansion of the binding energy landscape around the optimal ligand for unit $u$.

## Current critique

### The harmonic approximation is local
- A Taylor expansion of a binding landscape around its minimum is quadratic *near the minimum*. For ligands close to the unit's optimum, this is reasonable.
- For ligands far from the optimum (the majority case in a diverse environment), the real landscape is effectively flat: the ligand simply doesn't dock, so $\Delta E \approx 0$, not $\Delta E \to \infty$.
- The current model grows quadratically without bound: a "very wrong" ligand has a hugely unfavorable $EC_{50}$ rather than no binding at all. This changes the optimization geometry — receptors get pushed to *repel* non-targets as strongly as they're *attracted* to targets, which isn't biophysical.

### The metric is isotropic
- Squared Euclidean distance is isotropic in latent space. Real binding pockets are anisotropic: a steric clash in "volume" costs differently from a polarity mismatch.
- The optional weight vector $w_u$ provides axis-aligned anisotropy, but full anisotropy would require a Mahalanobis form $(\textbf{v}_u - \textbf{v}_\ell)^\top \Sigma_u^{-1} (\textbf{v}_u - \textbf{v}_\ell)$, at the cost of $O(D^2)$ parameters per unit.

### $E_\text{base}^{(u)}$ does too much
- It absorbs the optimal binding energy, the closed-state conformational energy $\epsilon_u$, and the "lock-and-key" floor.
- This recreates the sloppy-mode degeneracy already noted in the JPCB 2017 paper (Section 2.3 of `02_biophysics_mwc.md`): the optimization can drift along directions that mix "intrinsically leaky" with "binds tightly to optimal ligand."

## Candidate alternatives

### A. Bounded energy via a saturating function
$$\tilde{E}_o^{(u,\ell)} = E_\text{base}^{(u)} + E_\text{max}^{(u)} \cdot \left(1 - \exp(-\|\textbf{v}_u - \textbf{v}_\ell\|^2 / \lambda_u^2)\right)$$
- Bounds the energy of non-binders at $E_\text{base}^{(u)} + E_\text{max}^{(u)}$.
- Introduces a width parameter $\lambda_u$ controlling receptor breadth (broad vs narrow tuning — exactly what experimentalists call "selectivity").
- Cost: one or two extra parameters per unit.
- Locally quadratic (Taylor expansion around $v_u$) — preserves the harmonic justification near the optimum while fixing the far-from-optimum pathology.

### B. Gaussian affinity (equivalent reparameterization)
Define an *affinity* $A^{(u,\ell)} = \exp(-\|v_u - v_\ell\|^2 / \lambda_u^2)$ taking values in $(0, 1]$.
Then $\ln K_o^{(u,\ell)} = E_\text{base}^{(u)} - \ln(A^{(u,\ell)} + \kappa)$ for some floor $\kappa > 0$.
- Same biophysical content as (A), different parameterization.
- More natural for thinking about "fraction of receptors bound" — affinity directly maps to docking probability.
- The floor $\kappa$ represents non-specific binding.

### C. Bilinear (inner-product) energy
$$\tilde{E}_o^{(u,\ell)} = E_\text{base}^{(u)} - \textbf{v}_u^\top \textbf{v}_\ell$$
- Drops the geometric metaphor entirely. Each unit and ligand is a vector; binding energy is their alignment.
- Naturally anisotropic.
- Low-rank factorization of the affinity matrix $K_o^{(u,\ell)}$ — analogous to assuming a small number of independent chemical descriptors.
- Cost: energy unbounded below (perfect alignment $\to -\infty$); requires norm constraints or weight decay on $v_u, v_\ell$.
- Loses the "ideal ligand" geometric picture.

### D. Negative-distance / cosine
$$\tilde{E}_o^{(u,\ell)} = E_\text{base}^{(u)} - \alpha_u \cos(\textbf{v}_u, \textbf{v}_\ell)$$
- Bounded above and below.
- Pure orientation-based binding (magnitude of $v_\ell$ doesn't matter).
- Cleanest if you want $D$ to represent a "type signature" rather than a "physical descriptor."

### E. Keep quadratic but cap
$$\tilde{E}_o^{(u,\ell)} = E_\text{base}^{(u)} + \min\left(\|\textbf{v}_u - \textbf{v}_\ell\|^2, \; E_\text{cap}^{(u)}\right)$$
- Minimal departure from the current model.
- Hard cap is non-differentiable at the boundary; softplus or log-sum-exp smoothing would be needed in practice.
- Preserves the local harmonic interpretation.

## Open questions for the next discussion

1. **Which pathology actually bites in practice?**
   - Test: during training, look at the distribution of $\|v_u - v_\ell\|^2$ across the batch. If the optimization is exploiting the quadratic blow-up, you'll see receptors with very large $\|v_u - v_\ell\|^2$ values for the "least preferred" ligands, growing without bound during training. If not, the quadratic model is fine in practice even if biophysically unrealistic.

2. **What does the answer change?**
   - If the quadratic blow-up is benign, the simplest defense is "quadratic is a local harmonic expansion and the optimization happens to stay near the optimum." Keep the current model.
   - If the blow-up is exploited, option (A) is the smallest change that fixes it while preserving the geometric picture.

3. **Does the choice interact with the family structure?**
   - Bounded energies (A, B, D) imply that receptors have a *finite range* — they're blind to ligands sufficiently far away. This is more compatible with a multi-family environment where each receptor specializes in one family.
   - Unbounded energies (C, current quadratic) make every receptor sensitive to every ligand to some degree. This favors a continuous coding scheme.

4. **Does the anisotropy question (Mahalanobis / $w_u$) need to be answered separately, or is it bundled with the choice of energy form?**

5. **Is the $E_\text{base}^{(u)}$ degeneracy a real problem for optimization, or just an interpretability concern?**
   - If just interpretability, can be fixed post-hoc by gauge-fixing (e.g., constrain $\sum_u E_\text{base}^{(u)} = 0$ or set one unit's $E_\text{base}$ to zero by convention).
   - If it causes optimization drift, would need to be addressed structurally.

## What "simpler" looks like

The user's stated goal: simplify the environment so that explaining the model's behavior doesn't grow proportionally with environment complexity.

- **Simplest possible:** Option (C), bilinear, with $E_\text{base}^{(u)} = 0$ for all $u$. Pure matrix factorization of affinities. Loses a lot of structure but is mathematically transparent.
- **Simplest defensible:** Keep current quadratic, but explicitly argue that for the capacity claims we operate in the local-harmonic regime and check empirically that energies don't blow up.
- **Cleanest physical story:** Option (A), saturating quadratic. Local harmonic + finite range. Costs one parameter per unit.

The decision should be driven by which experiment is being run, not by which model is "more correct" in isolation.
