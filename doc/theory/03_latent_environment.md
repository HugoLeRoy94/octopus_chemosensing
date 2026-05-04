## 3. Generating Ligand/Unit interactions : The Environment

In previous iterations, the interaction energies between receptors and ligand families were modeled using independent parameter matrices ($\bar{\Sigma}$ and $\bar{M}$). The current implementation abstracts these interactions into a **Chemical Latent Space**, which naturally captures the cross-correlations between families and units by relying on geometric distances.

### 3.1 Chemical Latent Space and Affinities

We assume that the receptor array is designed for sensing a large number of ligands categorized into families. Instead of building independent microscopic interactions, each unit $u$ and each ligand family $f$ is assigned a vector embedding in a $D$-dimensional latent space:

* $v_u \in \mathbb{R}^D$: The learnable latent representation of a protein unit.
* $v_f \in \mathbb{R}^D$: The fixed latent representation of a ligand family, normalized to a unit hypersphere.

The fundamental assumption is that a smaller geometric distance in this latent space corresponds to stronger chemical binding (lower energy).

Additionally, we allow every units to "specialize" in a specific type of interaction, thus adding a weight vector $w_u$ to allow units to be more or less sensitive to the distance in a specific dimension.

**1. Open State Energy:**
The mean interaction energy for the open state, $\mu_o^{(u,f)}$, is derived from the squared Euclidean distance between the unit and the family embeddings:


$$\mu_o^{(u,f)} = E_{\text{base}}^{(u)} + \|\textbf{v}_u - \textbf{v}_f\|^2$$


where $E_{\text{base}}^{(u)}$ is a baseline energy parameter for unit $u$. *(Note: In the code implementation, an optional weight vector $w_u$ may be used to allow units to specialize in specific dimensions, but this is abstracted mathematically).*


**2. Closed State Energy (Thermodynamic Constraint):**
As established, a functional ion channel must have a greater dissociation constant when closed ($K_c > K_o$), which requires $\Delta E_c > \Delta E_o$. To strictly enforce this physical constraint, the closed-state mean energy is modeled as the open-state energy plus a strictly positive shift:

<br>
*(Crucial Note: Because we operate in the threshold limit defined in Section 2, the receptor activation is strictly governed by $EC_{50}$. This means we can completely bypass the calculation of closed-state energies in the current loss functions, and evaluate the array strictly on the open-state geometry).*


### 3.2 Stochastic Energy Sampling

While the mean energies ($\mu_o, \mu_c$) are defined by the latent geometry, ligands within the same family still exhibit natural variation. When a ligand from family $f$ is encountered, its specific interaction energy with unit $u$ in state $\alpha \in \{o, c\}$ is drawn from a Normal distribution:


$$E_\alpha^{(u,f)} \sim \mathcal{N}(\mu_\alpha^{(u,f)}, \sigma_\alpha^{(u,f)})$$


where $\sigma_\alpha^{(u,f)}$ is a learnable standard deviation specific to that unit-family pair.

### 3.3 Concentration Models

The environmental concentration $c$ of a ligand family is handled via independent concentration models. Based on biophysical assumptions, the environment can be modeled in two ways:

* **LogNormal Concentration:** Concentration spans orders of magnitude, modeled as $\log_{10}(c) \sim \mathcal{N}(\mu_c^f, \sigma_c^f)$.
* **Normal Concentration:** A simple Gaussian assumption, modeled as $c \sim \mathcal{N}(\mu_c^f, \sigma_c^f)$, clamped at $c > 0$ to satisfy physical constraints.

---

We assume that the receptor array is designed for sensing a large number of ligand. Any ligand is defined by 

1. its interaction with every units: $\{E_\alpha^{u,\ell}\}_{\alpha,u}$
2. A probability distribution of concentration: $p_\ell(c)$

To model the random encounter of a ligand, we categorize them into families. All the ligands of a given familly is assumed to have similar properties. As a consequence, the random encounter of a ligand of a given familly lies within a fairly narrow distribution (Normal) specific to each sub-unit.
For a given familly $f$:
$$
P_\alpha^{u,f}( E) = P(E | \alpha,u,f) = \mathcal{N}(\mu_\alpha^{u,f},\sigma_\alpha^{u,f}) = \sigma_\alpha^{u,f} \mathcal{N}(0,1) + \mu_\alpha^{u,f}
$$
Where $\alpha \in \{ o, c\}$. The goal is to optimize the mean and standard deviation of the probability distribution of each individual sub-unit to each familly, in order to maximize the mutual information between a ligand's identity and/or its concentration in a given environment.

In total, the number of parameters to optimize for a single receptor is:
$$
N_\text{param} = \underbrace{n_\text{unit}}_\text{\# of different units} \times ( \underbrace{\underbrace{2}_{\mu,\sigma} \times \underbrace{2}_\text{o,c} \times \underbrace{N_\mathcal{F}}_\text{\# of families}}_{interaction} + \underbrace{n_\epsilon}_ \text{\# parameter for $\epsilon$}) + \underbrace{N_\mathcal{F}\times n_c}_\text{\# parameter for the concentration}