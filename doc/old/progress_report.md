# Mutual Information of Heteromerized Ion-Channel
In this work, we investigates the generation of functional diversity in octopus chemosensing channels. These channels are pentamers ($k_\text{sub}=5$) formed from a pool of proteins expressed by $n_\text{genes} = 26$ different genes (e.g., 'CR518', 'CR918', etc.).
In this work, we quantify the coding capability of such an array for the identification of the identity and concentration of a ligand in the environment. To do so, we first introduce our model for individual ion-channel, and derive the relation between the response of a homo-pentamer and hetero-pentamers. We then introduce a model of the environment that allows us to generate an arbitrary number of ligands, with arbitrary complexity. Finally, we compute the mutual information between the environment and an array of receptors. This work allow us to compare the efficiency of hetero-pentamers vs homo-pentamers.


## I. MWC model for ion channel

Each receptor is made of $k_\text{sub}$ sub-units among $n_\text{genes}$ possible units.

The response of a receptor made of the $\mathcal{U}$ sub-units to a ligand $\ell$ at concentration $c$ can be written using the MWC model for ion channel opening:
$$
p_o^r(c,\ell) = \frac{\prod_{u\in \mathcal{U}} (1+c/K_o^{(u,\ell)})}{\prod_{u\in \mathcal{U}} (1+c/K_o^{(u,\ell)}) + e^{-\epsilon_\mathcal{U}}\prod_{u\in \mathcal{U}}(1+c/K_c^{(u,\ell)})}
$$
Where the $K_o$ and $K_c$ are respectively the dissociation constant of the ligand with the channel in the open and closed state.

in [1], they show that this function is over-parametrized, meaning that several combination of values for $K_c$, $K_o$ and $\epsilon$ can lead to the same curve. As a result, they find a simplified model by setting : $\tilde{K}_o = K_o e^{\beta \epsilon }$:
$$
p_o^{(r,l)}(c) = \frac{\prod_u \left(\frac{c}{\tilde{K}_o^{u,l}}\right)}{\prod_u\left(\frac{c}{\tilde{K}_o^{(u,l)}} \right) + \prod_u\left(1+\frac{c}{K_c^{(u,l)}} \right)}.
$$

In many sensory systems, downstream neural processing relies on highly thresholded, binned, or even binary signals (e.g., a neuron firing an action potential or remaining silent).
We work in the limit where the activation of a heteromeric receptor is strictly governed by its $EC_{50}$ : the concentration at which the opening probability is equal to $0.5$.
We find that the half-activation concentration for a heteromer is the geometric mean of its individual subunit affinities:

$$EC_{50}^{(r,\ell)} = \left( \prod_{u=1}^{k_\text{sub}} \tilde{K}_o^{(u,\ell)} \right)^{1/k_\text{sub}}$$

As a result, the response of a receptor to a ligand is given by a single set of parameters that are the effective dissociation constant in the open states, of the individual sub-units: $\tilde{K}_o^{u,\ell}$. In the following we drop the $\sim$.



## II Dissociation constants

Here we derive the dissociation constants from microscopic interactions. We consider the chemical reactions:
$$
U_o + L \rightarrow U_o L
$$
where $U_o$ denote the sub-units of a channel in the open state, $L$ is the ligand, and $U_oL$ is the ligand-channel dimers.
The the dissociation constants reads :
$$
K_o = \frac{[U_o][L]}{[U_oL]}
$$
Where the brackets denotes concentration. From a statistical mechanics point of view, we can write the concentration of any specie in the state $\alpha$ (open/close or bound/unbound) as:
$$
[\alpha] \propto p_\text{eq}(\alpha) \propto e^{\beta E_\alpha},
$$
Where, $p_\text{eq}$ is the equilibrium probability of being in state $\alpha$. As a result, we write the dissociation constants as:
$$
K_o^{(u,\ell)} = \exp\left[ \Delta E_o^{u,\ell} \right]
$$
with:
$$
\Delta E_o^{(u,\ell)} = E(U_oL) - E(U_o) - E(L) 
$$

## III. Model of the environment
### III.1 Chemical Latent Space

Following previous work, showing that there exist a abstract space wherein close distance translate into similar odors, and thus similar receptors affinity. 
We consider a $D$-dimensional chemical latent space, and denote the position of any chemical species in this space by a vector $\textbf{v}$.
In our approach, a ligand is defined by its dissociation constant or affinity energy defined in the previous section.
As a result, we write the energy difference of a ligand between the bound and the unbound state in term of chemical distance:
$$
\Delta E_o^{(u,\ell)} = E^{(u)} + \|\textbf{v}_u - \textbf{v}_\ell\|^2
$$
In this equation, the base energy: $E^{(u)}$, and the chemical properties $\textbf{v}_u$ of a unit defines efficiency of an array, whereas $\textbf{v}_\ell$ is a random variable modeling the random generation of a ligand.

### III.2 Stochastic Energy Sampling

To model the random encounter of a ligand, we categorize them into families. Each familly is defined by a $D$-dimensional distribution in the latent space with its own characteristics. For instance for a familly defined by a normal distribution:
$$
\textbf{v}_\ell \sim \mathcal{N}^D(\boldsymbol{\mu}_\ell, \boldsymbol{\sigma}_\ell),
$$
where $\boldsymbol{\mu}$ and $\boldsymbol{\sigma}$ are two $D$-dimensional vectors defining the mean and standard deviation of the distribution.

### III.3 Concentration Model
The environmental concentration $c$ of a ligand family is handled via independent concentration models. Based on biophysical assumptions, the environment is modeled by a Log-Normal Concentration: 

$$
\log_{10}(c) \sim \mathcal{N}(\mu_c^f, \sigma_c^f).
$$


## IV. Information Produced by an Array of Receptors

We now consider an array of $N$ receptors, and write the activity of this array as $\mathcal{A} = (a_0,a_1, \cdots , a_N)$. To compute the maximum coding capability of such an array, we compute the mutual information between the activity of the array and ligand's identity and/or concentration for a given environment:
$$
I(\mathcal{A};(c,\ell)) = H(\mathcal{A}) - H(\mathcal{A} | (c,\ell)),
$$
Where 
$$
H(\mathcal{A}) = \sum_\mathcal{A} p(\mathcal{A}) \log[p(\mathcal{A})],
$$
is the entropy of the array.
We consider the mapping between affinities, and activity as deterministic, noiseless: $H(A | (c,\ell)) = 0$. Therefore, the mutual information between the array and the environment is simply given by the entropy of the array $H(\mathcal{A})$.



A fundamental challenge in optimizing an array of chemosensors is that the mutual information maximized by the array, $I(\mathcal{A};(c,\ell)) = H(\mathcal{A}) - H(\mathcal{A} | (c,\ell))$, is strictly upper-bounded by the inherent complexity of the environment.
However, if the environment lacks sufficient variance or dimensionality, adding more receptors (or increasing combinatorial diversity through heteromerization) will yield diminishing returns. The array's entropy will aggressively plateau as the sensors are forced to extract redundant information from a limited probability space.


## Results

### Optimization Results

![This is my caption for the SVG.](1fam_5sensor_homo.svg)
### References 

[1]: Einav, Tal, and Rob Phillips. “Monod-Wyman-Changeux Analysis of Ligand-Gated Ion Channel Mutants.” The Journal of Physical Chemistry B 121, no. 15 (2017): 3813–24. https://doi.org/10.1021/acs.jpcb.6b12672.

[2]: Lee, Brian K., Emily J. Mayhew, Benjamin Sanchez-Lengeling, et al. “A Principal Odor Map Unifies Diverse Tasks in Olfactory Perception.” Science 381, no. 6661 (2023): 999–1006. https://doi.org/10.1126/science.ade4401.