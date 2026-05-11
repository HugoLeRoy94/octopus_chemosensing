## 2. Single Receptor Characteristics

### 2.1 MWC Model

Each receptor is made of $k_\text{sub}$ sub-units among $n_\text{unit}$ possible units.
We denote $\textbf{r} = (r_0 \cdots r_{k_\text{sub}}) \in \mathbb{R}^{n_\text{unit} \times k_\text{sub}}$ the identity of a receptor characterized by the identity of each of the units that composes it. 
The response of a receptor $\textbf{r}$ to a ligand $\ell$ at concentration $c$ can be written using the MWC model for ion channel opening:
$$
p_o^r(c,\ell) = \frac{\prod_{u\in \mathcal{U}} (1+c/K_o^{(u,\ell)})}{\prod_{u\in \mathcal{U}} (1+c/K_o^{(u,\ell)}) + e^{-\epsilon_\mathcal{U}}\prod_{u\in \mathcal{U}}(1+c/K_c^{(u,\ell)})}
$$

Where the $K$'s are called the affinities.

### 2.2 Microscopic model of a receptor
To make the link between the dissociation constant, we consider the chemical reactions:
$$
\begin{aligned}
&C_o + L \rightarrow C_o L\\
&C_c + L \rightarrow C_c L,
\end{aligned}
$$
where $C_o$ and $C_c$ are respectively the channel in the open and closed state, $L$ is the ligand, and $C_oL$ and $C_cL$ are the ligand-channel dimers.
The the dissociation constants reads :
$$
\begin{aligned}
K_o &= \frac{[C_o][L]}{[C_oL]}\\
K_c &= \frac{[C_c][L]}{[C_cL]},
\end{aligned}
$$
where $[C_o]$ and $[C_c]$ are respectively the concentration of open and closed channel, whereas, $[L]$ is the concentration of ligand in solution. From a statistical mechanics point of view, we can write for any specie $\alpha$ $[\alpha] \propto p(\alpha) = 1/Z e^{\beta E_\alpha}$. Where $Z$ is the partition function of the system. As a result, we write the dissociation constants as:
$$
K_o^{(r,\ell)} = \exp\left[ \Delta E_o^{r,\ell} \right]
$$
$$
K_c^{(r,\ell)} = \exp\left[ \Delta E_c^{r,\ell} \right],
$$
with:
$$
\begin{aligned}
\Delta E_o^{(r,\ell)} &= E(C_oL) - E(C_o) - E(L) \\
\Delta E_c^{(r,\ell)} &= E(C_cL) - E(C_c) - E(L)
\end{aligned}
$$
Notice that a good ion channel has $K_c > K_o$ (greater dissociation when closed). With that convention that means $\Delta E_c > \Delta E_o$.

Each affinity correspond to the response of a sub-unit to a ligand.
$\epsilon_\mathcal{U}$ is the energy difference between the closed and opened state of the channel when there are no ligands, it is thus independent of $\ell$.
We write it as the simple sum of the contribution of each individual unit:
$$
\epsilon_\mathcal{U} = \sum_u \epsilon_u
$$
$$
\epsilon_u = E_u(\text{closed}) - E_u(\text{open})
$$
### 2.3 Parameters Degeneracy

in [http://dx.doi.org/10.1021/acs.jpcb.6b12672] supplementary material, p.S25, they look at the hessian matrix:
$$
H_{i,j} = \frac{\partial p_o}{\partial x_i \partial x_j}
$$,
where $x = \beta (\epsilon, \Delta E_o, \Delta E_c)$. They find that the Hessian has three eigenvalues : $\{2\times 10^{-2}, 9\times  10^{-4}, 2\times 10^{-7} \}$. The last eigen value has a eigen vector $(2,1,0)$, and is sensitively smaller than the two others, highlighting a sloopy mode. in the $(2,1,0)$ the parameter $e^{-\beta \epsilon/2}K_o$ is constant. As a result, we can rewrite the opening probability assuming that:
$\left(e^{\beta \epsilon/2} + e^{\beta \epsilon/2} \frac{c}{K_o}\right)^2+\left(1+\frac{c}{K_c}\right)^2\approx\left(e^{\beta \epsilon/2}\frac{c}{K_o}\right)^2$
$$
p_o^{(r,l)}(c) = \frac{\prod_u \left(\frac{c}{\tilde{K}_o^{u,l}}\right)}{\prod_u\left(\frac{c}{\tilde{K}_o^{(u,l)}} \right) + \prod_u\left(1+\frac{c}{K_c^{(u,l)}} \right)},
$$
where $\tilde{K}_o^{(u,\ell)} = K_o^{(u,\ell)} e^{- \beta \epsilon_\mathcal{U} / k_\text{sub}}$.

### 2.4 Half-Activation Concentration ($EC_{50}$) Limit

In many sensory systems, downstream neural processing relies on highly thresholded signals (e.g., binned or binary). We focus on the limit where the activation of a heteromeric receptor is strictly governed by its $EC_{50}$, which is the ligand concentration that yields a 50% opening probability. 

Remarkably, the half-activation concentration for a heteromer is simply the geometric mean of its individual subunits' effective affinities:
$$EC_{50}^{(r,\ell)} = \left( \prod_{u=1}^{k_\text{sub}} \tilde{K}_o^{(u,\ell)} \right)^{1/k_\text{sub}}$$

As a result, a receptor's response to a ligand is fully defined by a single set of parameters: the effective open-state dissociation constants of its individual subunits. This effectively allows us to mathematically bypass modeling the closed states ($K_c$) in downstream algorithms.

### 2.5 Interaction energies

In term of interaction energy, we can write the dissociation constant as :
$$
\Delta E_o^{(u,\ell)} = E(u_oL) - E(u_o) - E(L)
$$
with $\epsilon_u = E_u(\text{closed}) - E_u(\text{open})$, we get
$$
\begin{aligned}
EC_{50}^{(r,\ell)} &= \exp \left[ \frac{1}{k_\text{sub}} \sum_{u=1}^{k_\text{sub}}E(u_oL) - E(u_o) - E(L) - E(u_c) + E(u_o) \right]\\
&= \exp \left[ \frac{1}{k_\text{sub}} \sum_{u=1}^{k_\text{sub}}\tilde{E}_o^{(u,\ell)} \right ]
\end{aligned}
$$
where $\tilde{E}_o^{(u,\ell)} = E(u_oL) - E(L) - E(u_c) = \Delta E_o^{(u,\ell)} - \epsilon_u$ is the effective open-state energy. The specific functional form used to parameterize $\tilde{E}_o^{(u,\ell)}$ in the latent space is described in Section 3.1; in the current implementation it takes the form $E_{\text{base}}^{(u)} + E_{\text{max}}^{(u)}(1 - e^{-\|\mathbf{v}_u - \mathbf{v}_\ell\|^2/\lambda^2})$, where $E_{\text{base}}^{(u)} = \tilde{E}_o^{(u,\ell_\text{opt})}$ is the value at the optimally matched ligand.