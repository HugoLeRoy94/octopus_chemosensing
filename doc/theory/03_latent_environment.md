## 3. Generating Ligand/Unit interactions : The Environment

In previous iterations, the interaction energies between receptors and ligand families were modeled using independent parameter matrices ($\bar{\Sigma}$ and $\bar{M}$). The current implementation abstracts these interactions into a **Chemical Latent Space**, which naturally captures the cross-correlations between families and units by relying on geometric distances.

### 3.1 Chemical Latent Space and Affinities

We assume that the receptor array is designed for sensing a large number of ligands categorized into broad chemical families. Instead of building independent microscopic interactions, each unit $u$ and specific ligand $\ell$ is assigned a vector embedding in a $D$-dimensional latent space:

* $\mathbf{v}_u \in \mathbb{R}^D$: The learnable latent coordinate of a protein unit.
* $\mathbf{v}_\ell \in \mathbb{R}^D$: The fixed base coordinate of a specific ligand $\ell$.

The fundamental assumption is that a smaller geometric distance in this latent space corresponds to stronger chemical binding (lower open-state energy).

**Open-State Energy — Saturating Affinity Kernel:**
The instantaneous open-state interaction energy between unit $u$ and observed ligand $\ell$ is:

$$\boxed{E_o^{(u,\ell)} = E_{\text{base}}^{(u)} + E_{\text{max}}^{(u)} \cdot \left(1 - \exp\!\left(-\frac{\|\mathbf{v}_u - \mathbf{v}_{\text{obs},\ell}\|^2}{\lambda^2}\right)\right)}$$

Three parameters shape this expression:

| Parameter | Type | Physical meaning |
|---|---|---|
| $E_{\text{base}}^{(u)}$ | per-unit, **learned** | Open-state energy at the **optimal** ligand ($\mathbf{v}_\ell = \mathbf{v}_u$, i.e. zero distance). Equals $\Delta E_o^{(u,\text{opt})} - \epsilon_u$ from the MWC derivation (Section 2.5). The EC50 at the perfectly matched ligand is $\exp[E_{\text{base}}^{(u)}]$. |
| $E_{\text{max}}^{(u)}$ | per-unit, **learned** (constrained $> 0$ via softplus) | Maximum extra energy cost for a fully mismatched ligand. As $\|\mathbf{v}_u - \mathbf{v}_\ell\| \to \infty$, $E_o \to E_{\text{base}}^{(u)} + E_{\text{max}}^{(u)}$. Controls how selective each unit is. |
| $\lambda$ | global, **fixed hyperparameter** | Affinity length scale — the characteristic distance beyond which a ligand is "fully rejected". Smaller $\lambda$ means sharper selectivity. |

The EC50 of a heteromeric receptor in the saturating model is:

$$\ln EC_{50}^{(r,\ell)} = \frac{1}{k_\text{sub}} \sum_{u=1}^{k_\text{sub}} \left[ E_{\text{base}}^{(u)} + E_{\text{max}}^{(u)} \cdot \left(1 - e^{-\|\mathbf{v}_u - \mathbf{v}_{\text{obs},\ell}\|^2 / \lambda^2}\right) \right]$$

**Why the saturating form?**

The previous model used a pure quadratic $E_o = E_{\text{base}} + \|\mathbf{v}_u - \mathbf{v}_\ell\|^2$, which grows without bound. A ligand that is far off in chemical space would accumulate an arbitrarily large energy — biophysically wrong, since a non-docking ligand simply contributes ~0 extra cost (the pocket never engages it). The saturating kernel has four key properties:

1. **Correct limiting behavior:** as distance $\to \infty$, the energy saturates at $E_{\text{base}} + E_{\text{max}}$ rather than diverging.
2. **Local quadratic recovery:** near the optimum ($\|\mathbf{v}_u - \mathbf{v}_\ell\| \ll \lambda$), a Taylor expansion gives $E_o \approx E_{\text{base}} + \frac{E_{\text{max}}}{\lambda^2}\|\mathbf{v}_u - \mathbf{v}_\ell\|^2$, recovering the harmonic approximation. Gradients and local geometry are therefore unchanged.
3. **Clean parameter meanings:** $E_{\text{base}}$ is unambiguously the EC50 at the optimal ligand; $E_{\text{max}}$ is unambiguously the selectivity ceiling. They enter additively and have distinct physical meanings, avoiding the sloppy-mode degeneracy identified in the JPCB 2017 MWC paper (supplementary, p. S25).
4. **$\lambda$ decouples breadth from peak depth:** $\lambda$ appears only in the exponent alongside the distance, never multiplying or adding to $E_{\text{base}}$. This is consistent with the Cleland et al. 2021 (odorsampling) convention where the analogous $\sigma$ parameter is normalized out.

**Parameter initialization:**

* $E_{\text{base}}^{(u)}$ is initialized to $\mathbb{E}[\ln c]$, the expected log-concentration averaged over all ligands. This ensures that each receptor's EC50 matches the expected ligand concentration at the start of optimization.
* $E_{\text{max}}^{(u)}$ is initialized to $10$ (in log-concentration units), placing the saturation ceiling well above any typical concentration in the batch, so the optimizer can freely tune selectivity downward.

*(Crucial Note: Because we operate in the threshold limit defined in Section 2, the receptor activation is strictly governed by $EC_{50}$. This means we can completely bypass the calculation of closed-state energies in the current loss functions, and evaluate the array strictly on the open-state geometry).*

### 3.2 The Fixed World vs. Dynamic Mixtures

Instead of generating completely random new chemicals at every step, the environment operates on a two-tier paradigm: a fixed "world" of specific ligands initialized once, and dynamic, noisy "sniffs" (mixtures) of those ligands sampled during optimization.

#### 3.2.1 Fixed Initialization (The Ligand Pool)
During the environment setup, $N_{\text{ligands}}$ specific ligands are instantiated. For each ligand $\ell$, the following parameters are drawn once and permanently fixed:
* **Family Assignment:** The ligand is assigned to a broad family $f$.
* **Base Coordinate ($v_\ell$):** Drawn from the family's spatial distribution: $v_\ell \sim \mathcal{N}(v_f, \sigma_{\text{shape}})$.
* **Concentration Parameters:** A unique mean and variance ($\mu_{c,\ell}, \sigma_{c,\ell}$) governing its specific concentration distribution.
* **Presence Probability ($p_{\text{presence},\ell}$):** The independent probability that this specific ligand will appear in any given environmental sample.

#### 3.2.2 Dynamic Sampling (Batch Optimization)
Every time a batch of sensory events is sampled, the environment dynamically generates a unique *mixture* of these existing ligands:

1. **The Mixture Mask:** A multi-hot binary vector $M$ is generated by drawing from a Bernoulli distribution for each ligand based on its $p_{\text{presence},\ell}$. ($M_\ell = 1$ if present, $0$ if absent).
2. **Concentrations:** A specific concentration $c_\ell$ is sampled for every ligand from its assigned distribution. It is then multiplied by $M_\ell$ so absent ligands have a concentration of $0$.
3. **Observation Noise:** While the ligands have fixed base coordinates, they don't dock in the exact mathematical orientation every time. To simulate docking variance or thermal fluctuations, a small Gaussian noise is added to the base coordinate to create the final observed coordinate:
   $$v_{\text{obs},\ell} \sim \mathcal{N}(v_\ell, \sigma_{\text{noise}})$$
   *(Note: $\sigma_{\text{noise}} \ll \sigma_{\text{shape}}$)*

The instantaneous interaction energies $E_o^{(u,\ell)}$ are then computed using these noisy observed coordinates for all ligands currently present in the mixture.

### 3.3 Concentration Models

The environmental concentration $c_\ell$ of a specific ligand $\ell$ is handled via independent concentration models. Based on biophysical assumptions, the environment can be modeled in two ways:

* **LogNormal Concentration:** Concentration spans orders of magnitude, modeled as $\log_{10}(c_\ell) \sim \mathcal{N}(\mu_{c,\ell}, \sigma_{c,\ell})$.
* **Normal Concentration:** A simple Gaussian assumption, modeled as $c_\ell \sim \mathcal{N}(\mu_{c,\ell}, \sigma_{c,\ell})$, clamped at $c > 0$ to satisfy physical constraints.