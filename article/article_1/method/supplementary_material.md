# Supplementary Material

This document details the five points the main text defers to supplementary
material: the Monod-Wyman-Changeux (MWC) extension to heteromeric receptors,
the saturating affinity kernel, the generation of the morphochemical
environment, the "perfect array" regime, and the bounds used to bracket the
joint Shannon entropy $H(A)$ in Fig.~2.

## S1. MWC model for heteromeric receptors

A receptor $r$ is a pentamer ($k_\text{sub}=5$ sub-units) built from a pool of
$n_\text{genes}$ possible unit types. Each unit $u$ in the ring can bind the
ligand $\ell$ in the channel's open ($o$) or closed ($c$) conformation, with
dissociation constants $K_o^{(u,\ell)}$ and $K_c^{(u,\ell)}$ set by the
binding free-energy difference between the ligand-bound and unbound states,
$K_\alpha^{(u,\ell)} = \exp[\Delta E_\alpha^{(u,\ell)}]$, $\alpha\in\{o,c\}$.
Extending the homomeric MWC model to a heteromer of arbitrary composition
gives the opening probability
$$
p_o^{(r,\ell)}(c) = \frac{\prod_{u\in r} (1+c/K_o^{(u,\ell)})}{\prod_{u\in r} (1+c/K_o^{(u,\ell)}) + e^{-\epsilon_r}\prod_{u\in r}(1+c/K_c^{(u,\ell)})},
$$
where $\epsilon_r=\sum_{u\in r}\epsilon_u$ is the sum of the per-unit
closed-open energy gaps (leakiness) and the product over $u\in r$ runs over
the $k_\text{sub}$ units composing $r$. This is the sense in which the
opening probability of a heteromer is set by the affinities of its
constituent subunits.

**EC50 threshold limit.** Because an allosteric channel behaves as a
logarithmic sensor, downstream readout is well approximated by a single
threshold, the half-activation concentration $EC_{50}$ ($p_o=1/2$), rather
than by the full dose-response curve. For a heteromer this reduces to the
geometric mean of its subunits' effective open-state affinities
$\tilde K_o^{(u,\ell)}$,
$$
EC_{50}^{(r,\ell)} = \Bigl(\prod_{u\in r} \tilde K_o^{(u,\ell)}\Bigr)^{1/k_\text{sub}}
\;\;\Longleftrightarrow\;\;
\ln EC_{50}^{(r,\ell)} = \frac{1}{k_\text{sub}}\sum_{u\in r} E_o^{(u,\ell)},
$$
i.e. the activation threshold of a heteromer, in log-concentration space, is
the arithmetic mean of its individual units' open-state energies
$E_o^{(u,\ell)} = \ln \tilde K_o^{(u,\ell)}$. Closed-state energies no longer
enter this expression, so the array can be evaluated from open-state
energies alone. We write $c^* \equiv EC_{50}$ and take the receptor to
activate once $c$ exceeds this threshold, $A_r = \Theta(\ln c - \ln c^*_{r,\ell})$.

## S2. Saturating affinity kernel (Gaussian radial basis function)

Each unit $u$ and ligand $\ell$ carry coordinates $\mathbf{v}_u,\mathbf{v}_\ell \in \mathbb{R}^D$
in the $D$-dimensional morphochemical latent space, with $d=\|\mathbf{v}_u-\mathbf{v}_\ell\|$
the affinity-determining distance. The open-state energy is a Gaussian radial
basis function (RBF) of $d$,
$$
E_o^{(u,\ell)} = E_\text{base}^{(u)} + E_\text{max}^{(u)}\left(1-\exp\!\left(-\frac{d^2}{\lambda^2}\right)\right),
$$
so that, per unit, $c^*$ rises smoothly from $c^*_\text{min}=\exp(E_\text{base}^{(u)})$
at perfect match ($d=0$) and saturates at $c^*_\text{max}=\exp(E_\text{base}^{(u)}+E_\text{max}^{(u)})$
for a fully mismatched ligand ($d\to\infty$); $\lambda$ is the affinity length
scale setting the characteristic distance beyond which a ligand is fully
rejected. Combined with S1, the heteromer threshold $c^*_{r,\ell}$ is the
geometric mean, over its $k_\text{sub}$ units, of these saturating per-unit
thresholds. The saturating form (rather than an unbounded quadratic $d^2$) is
required because a non-docking ligand must contribute a bounded extra cost,
consistent with binding governed by subsite complementarity, and it reduces
to the quadratic (harmonic) form for $d\ll\lambda$.

## S3. Generation of the environment and the simulation loop

The ligand pool is built once and re-sampled every sensing event.

**Fixed pool.** $n_\text{families}$ family prototypes are placed in latent
space with a target average inter-family distance $d_\text{fam}$; $L$ ligands
are then drawn once around their family's prototype, $\mathbf{v}_\ell \sim
\mathcal{N}(\mathbf{v}_f,\sigma_\text{shape})$, together with a fixed presence
probability $p_\ell$ and concentration parameters $(\mu_{c,\ell},\sigma_{c,\ell})$.

**Sensing event (mixture draw).** For each event a Bernoulli presence mask
$n_\ell\sim\text{Bernoulli}(p_\ell)$ selects the ligands in the mixture, their
concentrations are drawn $c_\ell\sim\text{LogNormal}(\mu_{c,\ell},\sigma_{c,\ell})$,
and a small observation noise $\mathbf{v}_{\text{obs},\ell}\sim\mathcal{N}(\mathbf{v}_\ell,\sigma_\text{noise})$
($\sigma_\text{noise}\ll\sigma_\text{shape}$) perturbs the ligand's coordinate
to represent docking/thermal fluctuations.

**What is fixed vs. optimized.** The environment above — family centers,
$\{\mathbf{v}_\ell\}$, $\{p_\ell\}$, concentration parameters, and the kernel
length scale $\lambda$ (S2) — is sampled once and held fixed; it is not
learned. The only learned quantities are the receptor chemistry: the per-face
latent coordinates $\mathbf{v}_u^{+},\mathbf{v}_u^{-}$ and energy parameters
$E_\text{base}^{u,\pm}, E_\text{max}^{u,\pm}$ of every unit $u$ (S2), updated
by gradient descent to maximize the array's mutual information (defined in
the main text).

**From single-ligand threshold to mixture drive.** S1 defines activation for
a single ligand as a hard threshold, $A_r=\Theta(\ln c-\ln c^*_{r,\ell})$. In
a mixture, ligands compete for the same receptor, and each contributes to
activation in proportion to how far its concentration sits above its own
threshold, $c_\ell/EC_{50}^{(r,\ell)}$. Summing these occupancy ratios over
the ligands present defines the receptor's **drive**,
$$
\text{drive}_r(\{c_\ell\}) = \sum_\ell \frac{c_\ell}{EC_{50}^{(r,\ell)}},
$$
which generalizes the single-ligand rule: half-activation is reached exactly
when $\text{drive}_r=1$, consistently reducing to $c=EC_{50}$ for a single
ligand.

**Soft binning.** The hard threshold $\Theta$ has zero gradient almost
everywhere and cannot be optimized by backpropagation. We replace it by a
temperature-scaled sigmoid of the log-drive,
$$
p(r\mid\{c_\ell\}) = \sigma\!\left(\frac{\ln\,\text{drive}_r(\{c_\ell\})}{T}\right)
= \sigma\!\left(\frac{\displaystyle\ln\sum_\ell \exp\bigl(\ln c_\ell - \ln EC_{50}^{(r,\ell)}\bigr)}{T}\right),
$$
where the second form (a logsumexp) is the numerically stable implementation
used in practice, since concentrations span many orders of magnitude. $T$ is
annealed during training, from an initial value that keeps gradients alive
across the whole population of receptors down to a small target value; as
$T\to0$, $p(r\mid\{c_\ell\})$ recovers the hard threshold of S1 and $p\to
A_r\in\{0,1\}$, so the discrete activation pattern $A=(A_1,\dots,A_R)$ used
to define $H(A)$ in the main text is the $T\to0$ limit of this same relaxed
quantity, not a separate model.

## S4. The perfect-array regime

$H(A)$ can in principle depend on every parameter of the environment
(Table~S1 of Sec.~S3) as well as on $R$ itself, so a direct comparison
between homomer and heteromer arrays would be confounded by these choices
unless we first identify a regime where the comparison does not depend on
them. We therefore look for a window of parameter values inside which
$H(A)$ stops changing with the precise value of each parameter: moving
anywhere inside the window leaves the result essentially the same. For the
homomer array, the value that $H(A)$ settles to inside this window is its
architectural ceiling, $H(A)=R$: we call an array that reaches this ceiling
*perfect*. Crucially, in this same window, with the same environment and the
same fitting procedure, heteromer arrays are not perfect: their $H(A)$ falls
short of $R$ (main text). This is the reason the window matters. Because the
homomer array already reaches $R$ throughout it, neither the environment nor
the optimization can be blamed for that shortfall, so any gap measured
between heteromer and homomer arrays of the same $R$ in this window must be
a property of the heteromeric array itself, its subunit-sharing correlations
(Sec.~S1), not an artifact of the environment or of how it is fit.

Outside this window, one of three distinct mechanisms prevents $H(A)$ from
reaching $R$, even for the homomer array.

**1. The mixture itself is too simple (statistical).** $A$ is a
deterministic, noiseless function of the mixture, so $H(A)\le H(\text{mixture})$:
a deterministic map cannot create entropy. If the mixture statistics
themselves carry little entropy (mixtures too sparse, too few ligands),
$H(\text{mixture})$ binds before $R$ does, and $H(A)$ is capped below $R$ no
matter how the receptors are placed.

**2. The receptor geometry lacks the capacity to realize $R$ independent
dichotomies (geometric).** Even when $H(\text{mixture})$ is large, the array
can only reach $H(A)=R$ if $R$ receptors can be placed in the latent space so
as to realize $R$ mutually decorrelated, half-active dichotomies of the
ligand pool. Whether they can is a question about the expressive capacity of
the receptor family in $D$ dimensions, not about $H(\text{mixture})$.

**3. The optimizer does not converge (numerical).** The first two mechanisms
are properties of the physical model. Separately, the gradient-based
procedure of Sec.~S3 used to reach and measure $H(A)=R$ can itself fail,
typically for reasons of computational cost rather than of principle.

| Parameter | Insensitive range | Justification | Mechanism |
|---|---|---|---|
| $\mu_S$, $L$ | not too small | Too sparse or too few ligands makes $H(\text{mixture})$ itself small, capping $H(A)$ below $R$ regardless of receptor geometry. | 1, statistical |
| $\rho=\sigma_\text{shape}\sqrt{D}/\lambda$ | $\sim[0.2,1]$ | Below: ligands within a family all sit inside one receptor's ball, only family identity is recoverable. Above: each ball covers at most one ligand, collapsing the realizable dichotomies. | 2, geometric |
| $d_\text{fam}/\lambda$ | $\sim[0.5,1.5]$ | Same mechanism as $\rho$: below, neighbouring families' balls overlap; above, no ball can bridge between families. | 2, geometric |
| $D$ | as small as possible while satisfying the $\rho,d_\text{fam}/\lambda$ windows above | Larger $D$ never hurts geometric capacity, but costs more memory and more parameters to optimize per unit, so we use the smallest $D$ consistent with the windows above. | 3, numerical |

**Geometric capacity, qualitatively.** Mechanism 2 is a question about how
many independent dichotomies of the ligand pool $R$ receptors can jointly
realize, a question about the capacity of a family of geometric classifiers,
classically studied for linear classifiers (hyperplanes)
\cite{cover_geometrical_1965}: the VC dimension of a hyperplane in
$\mathbb{R}^D$ is $D+1$, meaning any $D+1$ points in general position can be
split into either of their $2^{D+1}$ possible dichotomies, but this
guarantee is lost beyond $D+1$ points. Our receptors are not hyperplanes.
Through the Gaussian kernel of Sec.~S2, each is a "ball" classifier,
responding to ligands within roughly $\lambda$ of its position rather than
on one side of a plane. Kernel classifiers of this kind are generally less
constrained by the ambient dimension than a hyperplane is, since a nonlinear
kernel implicitly maps the ligands into a much higher, for a Gaussian kernel
effectively infinite, dimensional feature space in which separating them is
easier, the same reasoning Cover used to motivate feature expansion. The
hyperplane bound above should therefore be read as a conservative reference
point, not a tight limit for our receptors. What we do keep from the
comparison is the direction of the effect: increasing $D$ can only make the
classifier at least as capable, never less, which is why the smallest $D$
consistent with the $\rho$ and $d_\text{fam}/\lambda$ windows above costs
nothing in capacity.

## S5. Bracketing the joint Shannon entropy

Exact evaluation of $H(A) = -\sum_A p(A)\log_2 p(A)$ requires enumerating all
$2^R$ activation patterns, which is intractable beyond $R\sim15$. We instead
bracket $H(A)$ with two certified, matched pairwise-distance estimators of
Kolchinsky and Tracey (KT), built from the per-receptor activation
probabilities $A_{i,r}\in[0,1]$ of $B$ sampled events, with mean per-event
conditional entropy $H_\text{cond} = \frac{1}{B}\sum_i \sum_r h_2(A_{i,r})$,
$h_2(a)=-a\log_2 a-(1-a)\log_2(1-a)$.

**Lower bound**, from the Bhattacharyya affinity between events $i,j$,
$$
\log \text{BC}(i,j) = \sum_r \log\!\Bigl(\sqrt{A_iA_j}+\sqrt{(1-A_i)(1-A_j)}\Bigr),
\qquad
H_\text{KT}^{\text{low}} = H_\text{cond} - \frac{1}{B}\sum_i \log_2\!\Bigl(\frac{1}{B}\sum_j \text{BC}(i,j)\Bigr).
$$

**Upper bound**, replacing the Bhattacharyya affinity by minus the
Kullback-Leibler (KL) divergence between the same product-Bernoulli
components,
$$
-\text{KL}(i\|j) = \sum_r \Bigl[A_i\log A_j + (1-A_i)\log(1-A_j)\Bigr] - \sum_r\Bigl[A_i\log A_i+(1-A_i)\log(1-A_i)\Bigr],
$$
$$
H_\text{KT}^{\text{up}} = H_\text{cond} - \frac{1}{B}\sum_i \log_2\!\Bigl(\frac{1}{B}\sum_j e^{-\text{KL}(i\|j)}\Bigr).
$$

Both sums keep the diagonal term ($i=j$: $\text{BC}=1$, $\text{KL}=0$) and
normalize by the full sample $B$, which anchors the resolvable-entropy
ceiling of both estimators at $\log_2 B$. $H_\text{KT}^{\text{low}} \le H(A)
\le H_\text{KT}^{\text{up}}$ is a certified bracket for any $B$, and the
bracket tightens as the activation patterns become well separated — which is
precisely the regime the optimization of $\{\mathbf{v}_u\}$ drives the array
toward. In Fig.~2, solid lines report $H_\text{KT}^{\text{up}}$ and dotted
lines $H_\text{KT}^{\text{low}}$; the gap between them, which widens with
$R$, is shown as the error bars.
