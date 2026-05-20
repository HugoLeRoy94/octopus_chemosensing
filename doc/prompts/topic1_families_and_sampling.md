# Topic 1 — Baked-in Families vs Emergent Clustering, and the Two-Tier Sampling Question

## Context

The current environment (see `03_latent_environment.md`) generates ligands in two structured ways:

1. **Family hierarchy.** Each ligand is assigned to a family $f$, with a family center $v_f$, and its base coordinate is drawn as $v_\ell \sim \mathcal{N}(v_f, \sigma_\text{shape})$.
2. **Two-tier sampling.** A fixed pool of $N_\text{ligands}$ is instantiated once with frozen base coordinates $v_\ell$. During each batch, an observation noise $v_{\text{obs},\ell} \sim \mathcal{N}(v_\ell, \sigma_\text{noise})$ is added on top, with $\sigma_\text{noise} \ll \sigma_\text{shape}$.

Two design choices are entangled here:

- **Are families baked into the data-generating process, or should they emerge from optimization?**
- **Is the ligand pool fixed (with observation jitter), or freshly sampled every batch?**

## Current critique

### On baked-in families
- Hard-coding $v_\ell \sim \mathcal{N}(v_f, \sigma_\text{shape})$ injects the cluster structure that the array is then asked to "discover." The mutual information $I(\mathcal{A}; F)$ from Section 6.2 partly measures whether the array exploited a prior we built in, not whether family-tuned receptors emerged from environmental statistics alone.
- The user's counter-argument: removing baked-in families and relying on post-hoc clustering yields "weak" families — clustering algorithms always find clusters, but they're not crisp. Baked-in families also enable richer environment designs (controlled within-family vs between-family separation, controlled number of families, controlled overlap).
- This counter-argument is fair: baked-in families are a *controllable knob*, and the alternative (emergent families via clustering) is statistically squishy.

### On the two-tier sampling
- The fixed pool of $N_\text{ligands}$ is computationally convenient but not biophysically motivated. Real chemical space is effectively infinite.
- $\sigma_\text{noise}$ is described as "docking variance / thermal fluctuation," but with $\sigma_\text{noise} \ll \sigma_\text{shape}$ it contributes negligibly to the energy variance; it functions mostly as regularization.
- The Bernoulli mixture mask treats ligand presence as independent, which is unrealistic (ripe fruit emits dozens of co-occurring volatiles) but the user has flagged this as separately addressable.

## Open questions for the next discussion

1. **Should families remain baked in, given the user's research goal of controllable environment complexity?** Or is there a hybrid: families baked in for studies of $I(\mathcal{A}; F)$, but turned off (single isotropic Gaussian) for the capacity-saturation phase-diagram experiments where families would be a confound?

2. **Should the fixed pool of ligands be replaced by fresh-sampling per batch?**
   - *Argument for fresh sampling:* simpler model (one fewer hyperparameter, $\sigma_\text{noise}$ is gone), no arbitrary pool size $N_\text{ligands}$, ligands become a true continuous distribution. The vocabulary bound $H(\mathcal{A}) \le \log_2[M(N+1)]$ from the Cover analysis becomes purely a function of batch size $B$, which is a knob you'd vary anyway.
   - *Argument for fixed pool:* enables ligand-specific properties (presence probability $p_{\text{presence}, \ell}$, ligand-specific concentration distribution $\mu_{c,\ell}$, $\sigma_{c,\ell}$). These encode "some odors are rare, some are common, some always come at high concentration" — a feature of real olfactory environments. Fresh sampling kills this.

3. **If we keep the fixed pool, is $\sigma_\text{noise}$ doing useful work, or can it be set to zero?** If the only role of $\sigma_\text{noise}$ is to provide gradient signal / regularization, that role might be better played by the concentration noise alone, or by an explicit regularizer on receptor parameters.

4. **Is there a clean separation: "structural" environment parameters (families, $D$, $M$) fixed once, vs "statistical" parameters (concentrations, mixture composition) sampled per batch?** This might let us keep the fixed pool's expressivity while removing $\sigma_\text{noise}$ as a separate axis.

## Decision matrix to populate

| Configuration | Fresh ligands per batch | Fixed pool + $\sigma_\text{noise}$ | Fixed pool, no noise |
|---|---|---|---|
| Capacity experiments (no families) | ? | ? | ? |
| Structural experiments ($I(\mathcal{A};F)$) | ? | ? | ? |
| Real-environment-statistics experiments (rare/common odors) | ? | ? | ? |

The goal of the next discussion is to fill in this table with a principled choice for each row.
