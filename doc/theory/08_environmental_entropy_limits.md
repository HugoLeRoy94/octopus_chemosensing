# Environmental Limits on the Array's Achievable Entropy

This document analyses how each parameter of the environment bounds the maximum
information the receptor array can carry. It is the conceptual backbone for
introducing the model: before comparing architectures, we must understand what
sets the ceiling that *any* architecture is competing against.

Notation follows the project convention: $R$ receptors, $L$ ligands, $D$ latent
dimensions, $\lambda$ affinity length scale, $\sigma_{\text{shape}}$ intra-family
spread, $d_{\text{fam}}$ inter-family distance, $p$ per-ligand presence
probability, $\sigma_{\ln c}$ log-concentration spread, $T$ activation
temperature.

---

## 1. The Two Ceilings

The array activity $\mathcal{A}$ is a **deterministic, noiseless** function of the
environmental input (the mixture $\{c_\ell\}$). Two hard bounds follow immediately.

**Architectural ceiling.** With $R$ binary receptors there are at most $2^R$
distinct activity patterns, so
$$H(\mathcal{A}) \le R \text{ bits.}$$

**Input-entropy ceiling.** A deterministic function cannot create entropy — it can
only preserve or destroy it. Therefore
$$H(\mathcal{A}) \le H(\text{input}).$$

The central message of this document: **for most of the biologically relevant
regime, it is the input-entropy ceiling — not the architectural one — that binds.**
The environment, through how it generates mixtures, decides how much information
*exists to be encoded*. The array can only ever be as informative as its input is
variable. Understanding the model therefore begins with understanding $H(\text{input})$
and how each parameter shapes it.

The input entropy splits into two channels that the array can exploit:
$$H(\text{input}) = \underbrace{H(M)}_{\text{composition / identity}} + \underbrace{H(c \mid M)}_{\text{concentration}}$$
where $M \in \{0,1\}^L$ is the binary presence pattern and $c$ the concentrations
of present ligands. These two channels respond *oppositely* to several parameters,
which is the source of most of the interesting structure.

---

## 2. The Composition Channel and the Partition Limit

### 2.1 Presence entropy

For independent presence with uniform probability $p$, the presence vector is $L$
independent Bernoulli coins:
$$H(M) = L \cdot H_b(p), \qquad H_b(p) = -p\log_2 p - (1-p)\log_2(1-p).$$
$H_b$ is maximised at $p = 0.5$ (one full bit per ligand) and collapses as
$p \to 0$:

| $p$ | $H_b(p)$ (bits/ligand) |
|---|---|
| 0.5 | 1.00 |
| 0.3 | 0.88 |
| 0.1 | 0.47 |
| 0.02 | 0.14 |
| 0.01 | 0.08 |

So lowering $p$ to make sparse mixtures directly bleeds composition entropy out of
the input. But the total $H(M) = L H_b(p)$ can still be large if $L$ is large — the
collapse per ligand is compensated by having more ligands. This is the first hint of
the $p$–$L$ trade-off.

### 2.2 The partition limit — why receptors saturate at $R_{\text{eff}} \approx 2\langle S\rangle$

Total presence entropy is not the whole story. Even when $H(M)$ is far above $R$,
the array may be unable to *use* all its receptors, because of a structural
constraint that mirrors Zwicker et al. (2016).

An optimal array (Zwicker's two principles) wants each receptor to fire about half
the time, $\langle a_n\rangle \approx \tfrac{1}{2}$, with receptors mutually
uncorrelated. To fire half the time, a receptor must "cover" a bundle of ligands
whose presence probabilities sum to $\sim \tfrac{1}{2}$. With uniform $p$, that is
about $\tfrac{1}{2p}$ ligands per receptor. The array therefore poses a **partition
problem**: divide the $L$ ligands into $R$ bundles, each summing to $\sim\tfrac{1}{2}$,
with minimal overlap.

The total coverage demand is $R \cdot \tfrac{1}{2p}$ ligand-slots against a supply of
$L$ ligands. Overlap becomes unavoidable — and with it, receptor–receptor
correlations that *destroy* information (Zwicker's quadratic covariance penalty) —
once demand exceeds supply:
$$\frac{R}{2p} > L \quad\Longleftrightarrow\quad R > 2pL = 2\langle S\rangle.$$

This defines an **effective receptor budget**
$$\boxed{R_{\text{eff}} \approx 2\langle S\rangle = 2pL}$$
the number of receptors the environment can usefully support through the composition
channel. Beyond $R_{\text{eff}}$, additional receptors are forced to overlap, become
correlated, and add little information. The array's composition entropy saturates
*below* $R$ bits — set by the mixture statistics, not by the receptor count.

This is the mechanism behind the empirical observation that **reducing $p$ to get
~2 groups per sniff reduces the achievable MI**: small $\langle S\rangle$ means small
$R_{\text{eff}}$, and the array's useful capacity shrinks with it.

### 2.3 Why increasing $L$ and the number of blocks compensates

$R_{\text{eff}} \approx 2pL$ depends on the *product* $pL$. Holding $\langle S\rangle$
fixed while raising $L$ (and lowering $p$ to match) keeps the mean mixture size the
same but raises the total presence entropy $H(M) = L H_b(p)$ — more independent coins,
each quieter. More ligands means more independent dimensions of variation, hence more
room to give each receptor a distinct, non-overlapping bundle, hence a higher ceiling.

With block correlation (see §4), the effective number of independent "presence atoms"
is closer to the number of blocks $k$ than to $L$: in the strong-correlation limit a
block fires as a unit, so the partition problem runs over $k$ super-coins, not $L$
ligands. Few large blocks → few atoms → low ceiling. Many blocks → many atoms → high
ceiling. This is why increasing the number of groups raises the achievable MI.

The unifying statement:
$$R_{\text{eff}} \approx 2 \times (\text{effective number of independent presence atoms}),$$
where the atom count is $\langle S\rangle$ for independent presence and drops toward
the active-block count as correlation strengthens.

---

## 3. The Concentration Channel — A Parallel Axis That Escapes Sparsity

The composition channel collapses as mixtures become sparse, but the **concentration
channel does not**. Even a single-ligand sniff carries concentration entropy
$\sim \log_2(\sigma_{\ln c}/\delta_c)$, where $\delta_c$ is the smallest resolvable
concentration difference. This is independent of $p$ and survives the sparse limit.

This is the high-dimensional generalisation of Zwicker's narrow-vs-wide concentration
trade-off, governed there by the sensitivity-distribution width $\lambda$ (their
notation) and here by $\sigma_{\ln c}$ relative to the activation sharpness $T$:

- $\sigma_{\ln c} \ll T$: concentration variation is invisible — every present ligand
  is read as "on", the channel carries nothing.
- $\sigma_{\ln c} \sim T$: a receptor sitting near a ligand's typical concentration
  fires about half the time *because of concentration variation alone* — maximal
  marginal entropy from this channel.
- $\sigma_{\ln c} \gg$ (inter-ligand mean spread): per-ligand concentration noise
  drowns any systematic concentration-ladder structure across ligands.

**Consequence for the architecture comparison.** As $p$ falls and the composition
channel shrinks, the optimal array shifts toward concentration coding — reading the
*level* of the few present ligands rather than their combinatorial identity.
Heteromers, with their geometric-mean EC50 ladder (a built-in spread of activation
thresholds), are structurally suited to exactly this. The predicted signature: as
$p$ decreases, $I(\mathcal{A};F)$ (identity) falls while $I(\mathcal{A};c)$
(concentration) holds up or relatively dominates. If the decomposition shows this,
it is the same physics confirming itself, and it sharpens the biological story —
*sparse environments push the array toward concentration coding, which is where the
heteromer ladder pays off.*

---

## 4. Parameter-by-Parameter Influence on the Ceiling

The table below summarises how each environment parameter bounds $H(\mathcal{A})$,
distinguishing the *composition* and *concentration* channels where they differ.

| Parameter | Effect on the ceiling | Mechanism |
|---|---|---|
| $p$ (presence prob.) | ↓ as $p\to 0$ (composition); ~flat (concentration) | Sets $\langle S\rangle = pL$, hence $R_{\text{eff}} \approx 2pL$. Low $p$ shrinks usable receptor budget. |
| $L$ (n_ligands) | ↑ with $L$ | Raises total presence entropy $L H_b(p)$ and the number of independent atoms; raises $R_{\text{eff}}$. |
| $\langle S\rangle = pL$ | sets $R_{\text{eff}} \approx 2\langle S\rangle$ | The single most direct control on the composition ceiling. |
| $k$ (n_blocks) | ↑ with $k$ | Sets the number of *independent* presence atoms under correlation; few blocks → low ceiling. |
| $\rho_{\text{block}}$ (correlation) | ↓ as $\rho\to 1$ at fixed $k$ | Strong correlation collapses $L$ ligands toward $k$ atoms, reducing effective composition entropy; also inflates mixture-size variance. |
| $\sigma_{\ln c}$ | ↑ then ~flat (concentration) | Below $T$: channel dead. Near $T$: channel maximal. Far above inter-ligand spread: ladder washed out. |
| $T$ (final temperature) | sets resolution of concentration channel | Smaller $T$ → sharper threshold → finer concentration discrimination, until noise/gradient limits bite. |
| $\rho = \sigma_{\text{shape}}\sqrt{D}/\lambda$ | ↓ outside $[\,\sim0.2, \sim1\,]$ | Controls whether ligands sit in the gradient-rich kernel regime. $\rho\gtrsim1$: saturation, energies collapse, both channels suffer. $\rho\ll0.2$: only family-level identity recoverable. |
| $d_{\text{fam}}/\lambda$ | ↓ outside $[\,\sim0.5,\sim1.5\,]$ | Too small: families overlap, "family" stops being a distinguishable variable. Too large: family centres mutually invisible across saturated tails. |
| $D$ (latent_dim) | ↑ with $D$ (capacity), but couples to $\rho$ | More dimensions → more linearly-separable dichotomies (Cover), raising the geometric ceiling; but at fixed $\sigma_{\text{shape}}$, larger $D$ inflates $\rho$ toward saturation. Hold $\rho$ fixed by scaling $\sigma_{\text{shape}}\propto 1/\sqrt D$. |
| $n_{\text{families}}$ | mild ↑; caps $I(\mathcal{A};F)$ | $I(\mathcal{A};F) \le \log_2 n_{\text{families}}$. Weak effect on total $H(\mathcal{A})$ once $\gtrsim 5$. |

### Reading the table for the model introduction

The headline points for introducing the model:

1. The array competes against an **input-entropy ceiling set entirely by the
   environment**. Neither homomers nor heteromers can exceed it.
2. The composition ceiling is governed by $R_{\text{eff}} \approx 2\langle S\rangle$
   — a partition limit, not a receptor-count limit.
3. The concentration channel is a parallel, sparsity-robust axis that becomes
   dominant exactly where composition fails — and is where the heteromer ladder is
   expected to matter most.
4. The kernel-geometry parameters ($\rho$, $d_{\text{fam}}/\lambda$) gate whether the
   array can *access* the input entropy at all; outside their good windows the ceiling
   is irrelevant because the array cannot reach it.

---

## 5. The Capacity-vs-Efficiency Distinction (consequence for the main comparison)

Because both architectures share the same input-entropy ceiling, the gene-vs-heteromer
comparison takes one of two characters depending on where the run sits relative to
$H(\text{input})$:

- **Below saturation** ($R < R_{\text{eff}}$, $H(\mathcal{A}) < H(\text{input})$): the
  comparison measures *capacity* — who extracts more bits from a not-yet-exhausted
  environment.
- **At saturation** ($R \gtrsim R_{\text{eff}}$, $H(\mathcal{A}) \to H(\text{input})$):
  the comparison measures *efficiency* — who approaches the shared ceiling with fewer
  genes.

These are different claims and should not be conflated. For each environment in the
prior, $H(M)$ (and an estimate of the concentration-channel contribution) should be
computed so that every run can be labelled capacity-limited or efficiency-limited.
A run where both arms have saturated against $H(\text{input})$ speaks to efficiency;
a run below saturation speaks to capacity.

---

## 6. Parameters Not Yet Characterised — Expected Influence

These have not been pinned down analytically or empirically; listed with the
influence I'd expect, to guide which to probe next.

**Concentration-channel structure**

- **Per-ligand vs block-shared concentration means.** If ligands within a source
  share a concentration mean (emission intensity), the concentration channel gains
  *source-identifying* structure even when presence alone is ambiguous. Expected:
  raises $I(\mathcal{A};\text{source})$, partially recovering identity information
  through the concentration axis. Worth testing once block correlation is in.
- **Correlation between concentration and presence.** Currently presence is correlated
  (copula) but concentration is drawn independently (the transport argument). If a
  source's ligands are not just co-present but co-*concentrated*, that is extra
  structure. Expected: modest increase in source-identity MI; possibly fragile to the
  same turbulence argument that justifies decoupling them. Low priority.

**Mixture-size distribution shape**

- **Variance of mixture size at fixed mean.** Correlation inflates $\text{Var}(S)$ at
  fixed $\langle S\rangle$. The array sees a more dispersed range of mixture
  complexities. Expected: a wider $S$-distribution may *help* (the array must handle
  both sparse and dense sniffs, exercising more of its dynamic range) or *hurt*
  (very dense sniffs saturate receptors, very sparse ones underuse them). Sign
  unclear — worth an explicit sweep of $\rho_{\text{block}}$ at fixed $\langle S\rangle$.
- **Non-uniform marginals (Gamma-distributed frequencies, à la del Castillo).** A few
  common ligands plus many rare ones. Expected: prevents empty mixtures naturally,
  and changes the partition problem (common ligands dominate receptor bundles). Likely
  shifts the array toward dedicating receptors to common ligands — the "frequency →
  dynamic range" effect del Castillo report. Could interact with the heteromer ladder
  in an interesting way. Medium priority.

**Receptor / geometry parameters**

- **$k_{\text{sub}}$ (subunits per receptor).** Fixed at 5 throughout. Larger
  $k_{\text{sub}}$ means the geometric-mean EC50 averages over more units → finer,
  more numerous achievable thresholds (bigger combinatorial palette) but each
  individual heteromer is more "averaged" and less distinctive. Expected: raises the
  heteromer ceiling at fixed $n_{\text{genes}}$ up to a point, then diminishing
  returns as averaging washes out specificity. This is arguably a *core* axis for the
  heteromer story and deserves its own analysis, not just robustness.
- **Observation noise $\sigma_{\text{noise}}$ (docking variance).** Currently small.
  Expected: acts like an effective floor on concentration resolution and blurs the
  kernel; mild reduction of both channels. Probably genuinely minor as long as
  $\sigma_{\text{noise}} \ll \sigma_{\text{shape}}$, but worth one confirmatory point.
- **Interface model vs classic model.** The pocket-at-interface variant changes how
  units combine into receptor affinities (averaged pocket embeddings rather than
  averaged subunit energies). Expected: changes the *geometry* of the achievable
  threshold palette, possibly the heteromer advantage itself. This is a modelling
  choice, not an environment parameter, but it sits upstream of the ceiling and should
  be checked for consistency of the headline.

**Family / latent structure**

- **Anisotropic or hierarchical family covariance.** Real chemical families are
  pancake-shaped and nested. Expected: anisotropy lowers the effective within-family
  dimensionality, which interacts with the Cover bound; hierarchy creates
  sub-family structure that a fine-grained array could exploit. The "cherry on top"
  already flagged — likely a robustness/extension item, not core.
- **Overlap geometry of families (beyond $d_{\text{fam}}/\lambda$).** Number of
  families relative to $D$: with many families in low $D$, family centres themselves
  crowd and the family-identity channel caps out below $\log_2 n_{\text{families}}$.
  Expected: a joint $n_{\text{families}}$–$D$ effect not captured by either alone.
  Low-to-medium priority.

**Temporal structure (out of current scope)**

- **Sniff-to-sniff dynamics / temporal correlation.** Both reference papers explicitly
  exclude this. Expected: a whole additional information axis (plume dynamics carry
  identity and location). Genuinely out of scope, but worth one sentence in the
  discussion acknowledging the ceiling computed here is the *static* ceiling.

### Suggested priority order

1. $k_{\text{sub}}$ — core to the heteromer mechanism; not just robustness.
2. Block-shared concentration means — cheap, directly tests the source-identity story.
3. Mixture-size variance at fixed mean ($\rho_{\text{block}}$ sweep) — disentangles
   correlation-structure effects from mean-size effects.
4. Gamma-distributed frequencies — realism + natural empty-mixture prevention.
5. Everything else — robustness / discussion items.
