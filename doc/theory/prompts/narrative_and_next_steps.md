# Narrative Strategy & Next Steps — Gene Expansion vs Heteromerization

## 1. The Headline

**Biological claim:** gene expansion and heteromerization are two competing evolutionary
strategies for generating olfactory-receptor diversity. Most species (mouse, dog, rat)
went the expansion route, with 1000+ functional OR genes plus the regulatory and
developmental machinery that comes with it. Octopuses appear to have gone the
heteromer route. The model quantifies the trade-off in *information bits per gene*.

**Sentence the abstract should aim at:**
> "An array of $k$-mer heteromers built from $n_{\text{genes}} \approx 25$ genes reaches
> $X\%$ of the mutual information attainable by $n_{\text{genes}} = 1000$ homomeric
> receptors, demonstrating that heteromerization is a viable evolutionary alternative
> to gene-family expansion."

The actual value of $X$ is the deliverable of the simulation campaign. Whether the
final framing is "competitive" (X ~ 50–80%), "complementary" (different regimes), or
"dominant" (X > 90%) follows the data, not the other way round.

## 2. Main-Text Figure Design

### The Three Ratios — Be Explicit About Which One Where

The comparison can be quantified by three different ratios, each answering a
different question. The figure should make clear which one each panel is computing.

| Ratio | Interpretation | Where it appears |
|---|---|---|
| $I / n_{\text{genes}}$ | "Information per gene" — the genomic-saving / efficiency-per-gene story | Panel B (main panel) |
| $I_{\text{hetero}}(n_g{=}25, R) / I_{\text{homo}}(n_g{=}R)$ | "Fraction of capacity retained at matched receptor count" — what heteromerization gives up by reusing genes | Panel C or annotated on Panel A |
| $I_{\text{env}_i} / I_{\text{env}_j}$ | Across-environment ratio — robustness check | SI; visualized as band widths in main figure |

Panel A must remain in absolute units because the x-axis spans different $R$ ranges
for the two arms (homomers only exist at $R = n_{\text{genes}}$) — a ratio needs a
common denominator, which only emerges once you fix one axis.

### Panels

**Panel A — raw capacity vs receptor count.**
- x-axis: $R$ (number of receptors in the array), log scale.
- y-axis: mutual information $I(\mathcal{A};\cdot)$ in bits (absolute).
- Curves:
  - Homomer reference: $R = n_{\text{genes}}$, so this is one curve parameterised by $n_{\text{genes}}$.
  - Heteromer curves: fixed $n_{\text{genes}} \in \{5, 8\}$ for the first shot, $R$ varied by subsampling.
- Reading: heteromer curves sit *below* the homomer curve at matched $R$ (heteromers
  individually are less informative per receptor), but each heteromer curve extends
  far beyond its own $n_{\text{genes}}$ on the x-axis.

**Panel B — efficiency per gene.**
- Same curves, but y-axis is $I / n_{\text{genes}}$.
- Heteromer curves cross *above* the homomer curve.
- This is the panel that carries the headline.

**Optional Panel C — matched-$R$ ratio.**
- For each $R$ where both arms exist, plot $I_{\text{hetero}}(n_g{=}25, R) / I_{\text{homo}}(n_g{=}R)$.
- Quantifies what fraction of "matched-budget" capacity is retained when reusing
  genes via heteromerization.
- Possibly redundant with Panel A's visual story; decide after the first shot.

### Aggregation Strategy

Each curve is the **median across an ensemble of environments** sampled from a
defended prior (specified in Section 4 below). Shaded bands show 10–90 percentile.

Rationale: evolution doesn't get to pick a canonical environment. Aggregating bakes
the robustness check into the main panel and makes the visual claim a statement
about medians across plausible worlds, not about a single hand-picked configuration.
The band widths *are* the across-environment robustness ratio, visualized rather
than computed as a separate number.

### Calibration Zone

For the first-shot figure (Section 6), the entire plot lives in the calibrated
regime: $R \leq 15$, exact Shannon entropy computable, Rényi-2 evaluable alongside
for comparison. Later extensions beyond $R = 15$ will need a visual marker (vertical
line or shaded background) labelled **"calibrated regime"** to demarcate where the
comparison is exact from where it relies on the symmetric-degradation argument.

### Heteromer Subsampling

Each heteromer curve is shown for **two subsampling strategies** as paired curves or
overlapping bands:

1. **Uniform random** over all combinations-with-replacement of $k_{\text{sub}}$ units
   from $n_{\text{genes}}$ genes. This is the maximum-entropy null.
2. **Cascading** (`generate_cascading_receptors`-style): homomers first, then 2-mers,
   then 3-mers, etc. Biologically more defensible since simpler complexes likely fold
   more reliably.

If both tell the same story → conclusion is robust. If they diverge → divergence is
biologically informative and worth a dedicated subsection.

## 3. Auxiliary Measurements Run Alongside the Main Sweep

The concentration-vs-identity decomposition is not a separate study — it is a set of
additional measurement functions invoked on every run in the main sweep. The
expectation is that heteromers gain most of their advantage from the natural EC50
ladder (geometric mean of subunit affinities), which should make them
disproportionately good at concentration coding.

To compute alongside the main entropy on every run:

- `full_array_entropy` — $H(\mathcal{A})$, the main signal.
- `mutual_information_ligand` — $I(\mathcal{A}; F)$, identity channel.
- `mutual_information_concentration` — $I(\mathcal{A}; c)$, concentration channel.
- `codeword_entropy` (plug-in + Miller-Madow) — for bias diagnostics.
- `mean_receptor_distance`, `mean_specialization_index` — geometric diagnostics.

The concentration-vs-identity story emerges from post-hoc analysis of these
measurements, presented as a sub-panel or one extra SI figure. No additional
simulation runs required.

## 4. Defended Environmental Prior

### Ligand Latent-Space Distribution

Ligand positions around a family center are drawn from a **multivariate Gaussian**:
$\mathbf{v}_\ell \sim \mathcal{N}(\mathbf{v}_f, \sigma_{\text{shape}}^2\, I_D)$.

This is the maximum-entropy distribution at fixed mean and covariance — the most
defensible "we are not assuming structure beyond the first two moments" choice. Real
chemical-feature distributions (Mordred descriptors, Morgan fingerprints, learned
embeddings) have heavier tails and anisotropic covariance, but a Gaussian captures
the dominant geometric behaviour and keeps the dimensional analysis clean. The
methods section should acknowledge this in one sentence; heavier-tailed or
anisotropic variants are out of scope for the planned campaign.

**Consequence for the dimensionless kernel argument.** A $D$-dimensional isotropic
Gaussian places its mass in a thin shell at radius $\sigma_{\text{shape}}\sqrt{D}$
(concentration of measure). The typical squared distance between a ligand and its
family center is $D\,\sigma_{\text{shape}}^2$. The kernel ratio that controls whether
ligands sit in the gradient-rich regime or in the saturation plateau is therefore:

$$\rho \equiv \frac{\sigma_{\text{shape}}\sqrt{D}}{\lambda}$$

The good regime is $\rho \lesssim 1$. With $\sigma = 0.1$, $\lambda = 1$, $D = 10$:
$\rho \approx 0.32$ (works — quadratic regime). With $\sigma = 1$, $\lambda = 1$,
$D = 10$: $\rho \approx 3.2$ (deep in saturation — the regime where empirical
failure was observed). $D$ is load-bearing in this ratio: it cannot be increased
indefinitely without inflating $\rho$, which is why the prior below samples $\rho$
directly and derives $\sigma_{\text{shape}}$ from $\rho$ and $D$.

### Parameter Bands

Parameters fall into three categories: (i) bound by understood failure modes, with a
plausible band, (ii) bound by sampling/optimization budget, (iii) effectively free.

#### (i) Bounded by Failure Modes — sample uniformly within the band

| Parameter | Plausible band | Lower-end failure mode | Upper-end failure mode |
|---|---|---|---|
| $\rho = \sigma_{\text{shape}}\sqrt{D}/\lambda$ | $0.2 \ldots 0.8$ | $\rho \ll 0.2$: intra-family ligands too close to family center; only identity-by-family is recoverable | $\rho \gtrsim 1$: kernel saturates, energies collapse, gradients vanish |
| $d_{\text{fam}} / \lambda$ | $0.5 \ldots 1.5$ | too small: families overlap as clusters; "family identity" stops existing | too large: family centers mutually invisible across saturated tails |
| $p_{\text{presence}}$ | $0.1 \ldots 0.5$ (per ligand) | $p \to 0$: most mixtures empty/singletons; mask entropy crashes | $p \to 1$ with large $L$: every mixture is the full bouquet; logsumexp washes out identity |
| $\sigma_{\ln c}$ | $0.5 \ldots 2$ (in nat units) | $\sigma_{\ln c} \ll T_{\text{final}}$: concentration variability invisible to sigmoid | $\sigma_{\ln c} \gg$ inter-ligand mean spread: per-ligand noise drowns ladder coding |
| $L$ (n_ligands) | $50 \ldots 200$ | $L \ll 2^R$: codeword space undersampled by mask alphabet | $L$ very large with $p \sim 0.5$: typical mixture too dense, logsumexp loses identity |
| $D$ (latent_dim) | $10 \ldots 20$ | low $D$: linear-separability bound (Cover) caps achievable entropy | high $D$: at fixed $\sigma_{\text{shape}}$, drags the kernel toward saturation via $\rho$ — control by scaling $\sigma_{\text{shape}} \propto 1/\sqrt{D}$ to hold $\rho$ fixed |

#### (ii) Bounded by Compute / Estimator

| Parameter | Constraint | Recommendation |
|---|---|---|
| $B$ (batch size) | Rényi-2 collision count must be statistically meaningful | $B \gtrsim 8 \cdot 2^{H_2/2}$; auto-scaled with $R$ (see §5.4) |
| $T_{\text{final}}$ | Should be small enough for crisp activation, large enough to avoid vanishing gradient | Annealed schedule handles training; eval at the final value |
| Exact Shannon vs Rényi vs blocked | Memory $\mathcal{O}(B \cdot 2^R)$ for Shannon | Exact Shannon up to $R \approx 15$; Rényi or blocked beyond |

#### (iii) Effectively Free / Pinned by Convention

- `distribution_type`: multivariate Gaussian (`"gaussian"`) — the existing default.
  Alternative distributions (heavier-tailed, anisotropic) are interesting but out of
  scope; possibly a future "cherry on top" once the core campaign is done.
- `n_families`: weak effect once $n_{\text{families}} \gtrsim 5$. Pin to a defensible
  value matching the chemistry literature.
- $\lambda$ scalar vs per-unit learnable: pin to global learnable (or fixed) for the
  main figure; per-unit learnable can be revisited later, but invites identifiability
  concerns (sloppy mode with $E_{\text{max}}$).

### Sampling Prescription for Aggregated Bands

```
For each environment sample i = 1 ... N_env (target N_env = 30):
    Draw L ~ Uniform{50, 100, 200}
    Draw D ~ Uniform{10, 15, 20}
    Draw rho ~ Uniform(0.2, 0.8)        # → set sigma_shape = rho * lambda / sqrt(D)
    Draw d_fam/lambda ~ Uniform(0.5, 1.5)
    Draw p_presence ~ Uniform(0.1, 0.5)
    Fix sigma_lnc at a defensible literature-matched value
    Fix n_families at a defensible value
    Fix lambda = 1 (units of latent space)
    Fix distribution_type = "gaussian"

For each environment: run 3 seeds, both homomer and heteromer arms,
both subsampling strategies for heteromers.
```

The key move is sampling $\rho$ directly and deriving $\sigma_{\text{shape}}$ from
$\rho$ and $D$ — this guarantees every environment sits in the gradient-rich regime
regardless of the drawn $D$. Sampling $\sigma_{\text{shape}}$ and $D$ independently
would let the prior wander into saturation whenever both happen to be drawn at the
upper end of their respective bands.

This is the prior the Methods section will defend, one sentence at a time, with one
citation per row of the table above where available.

## 5. Code Work Required Before Running

### 5.1 Heteromer Subsampling — Function to Build

A new function (or family) in `geometry.py` that, given:
- `n_genes` (size of available gene pool)
- `k_sub` (subunits per receptor, typically 5)
- `R_target` (number of receptors in the output array)
- `strategy` ∈ {"uniform_random", "cascading", possibly others}
- `seed` (for reproducibility)

returns a `(R_target, k_sub)` long tensor of unit indices suitable for passing as
`receptor_indices`.

Key requirements:
- **Uniform random**: sample without replacement from `combinations_with_replacement(range(n_genes), k_sub)`. When `R_target` exceeds the total combinatorial pool, this needs an explicit guard.
- **Cascading**: fill the array by complexity level — all $n_{\text{genes}}$ homomers
  first, then all 2-mers (unordered multisets), then 3-mers, until `R_target` is hit.
  Within each complexity level, order is either deterministic (sorted) or randomised
  (seeded). Existing `generate_cascading_receptors` can probably be adapted.
- **Determinism contract**: same seed + same parameters → same receptor set. Critical
  for warm-starting and for sweep reproducibility.

### 5.2 Warm-Starting Along R (Not n_genes)

The current `clone_with_extra_units` warm-starts along `n_genes`. For the main figure,
we need to warm-start along $R$ at fixed `n_genes`: extend the receptor index tensor
while preserving the trained `unit_latent`, `base_energy_u`, `max_energy_u_raw`,
`ligand_latent`, family centers, and concentration model.

The cleanest implementation is probably an alternative warm-start path in `run.py`
that hands forward only the *environment* (which is per-gene) and rebuilds the
`receptor_indices` tensor at each step. Worth checking whether the existing pipeline
covers this case before writing new code.

### 5.3 Concentration-vs-Identity Measurements (no new code)

Already implemented in `analysis_helper.py`. Action item: ensure
`mutual_information_ligand` and `mutual_information_concentration` are in
`measurement_fns` for every run in the main sweep. Zero new code, just config.

### 5.4 Auto-Adjusting Batch Size

A helper that, given the array size $R$ and the entropy estimator in use, returns
appropriate batch sizes for training and evaluation. Rationale: sweeping over $R$
with a single fixed batch size either wastes compute at small $R$ or undersamples
at large $R$.

Proposed scaling (Rényi-2):
- Train batch: $B_{\text{train}} = \max(B_{\min},\, c_{\text{train}} \cdot 2^{R/2})$
  with $B_{\min} \approx 512$ and $c_{\text{train}} \approx 8$ (gives ≥30 collisions
  in expectation when the array reaches its capacity).
- Eval batch: $B_{\text{eval}} = c_{\text{eval}} \cdot B_{\text{train}}$ with
  $c_{\text{eval}} \approx 4$ — bigger because eval is no-grad and we want tight
  Miller-Madow corrections.

Proposed scaling (exact Shannon, $R \leq 15$):
- Train batch: same formula as Rényi, capped by memory: $B \cdot 2^R \lesssim$ GPU
  budget (a few hundred MB for the joint-state tensor).
- Effectively: $B_{\text{train}}$ small to moderate, since the $2^R$ term dominates.

Integration: a `resolve_batch_sizes(cfg)` function called by `SimulationRunner` at
init time, overriding `cfg.batch_size` and `cfg.test_batch_size` only when they are
left unset (sentinel value, e.g. `"auto"`). Manually-set batch sizes always win, so
existing experiments aren't disturbed.

## 6. First-Shot Simulation Campaign (R ≤ 15)

The initial pass is the "calibration zone" reframed as the first attempt at the main
figure. Everything is in the exactly-tractable regime — exact Shannon entropy
computable, no estimator anxiety, theoretical ceiling reachable.

### 6.1 What This First Shot Establishes

1. **The figure layout works.** Panel A, B, (C) on real data; check that heteromer
   curves visually cross the homomer reference where expected.
2. **The optimizer is calibrated.** Compare exact Shannon $H(\mathcal{A})$ from
   `DiscreteExactLoss(entropy_type='shannon')` against the brute-force theoretical
   ceiling at $R = 10, 12, 15$.
3. **Rényi tracks Shannon.** Both estimators are logged on every run. The gap
   $H_{\text{Shannon}} - H_{\text{Rényi}}$ is the proxy bias; if it's small and
   consistent across configurations, Rényi is trustworthy for the later extension to
   larger $R$.
4. **Concentration-vs-identity story is visible.** With $I(\mathcal{A};F)$ and
   $I(\mathcal{A};c)$ logged alongside total $H(\mathcal{A})$, we can already see
   whether heteromers preferentially boost the concentration channel.
5. **Robustness bands close.** Across the 30-environment ensemble, band widths at
   each $R$ should be narrow relative to the homomer–heteromer gap.

### 6.2 Sweep Definition

| Axis | Values |
|---|---|
| $n_{\text{genes}}$ (homomers) | $\{3, 5, 8, 10, 12, 15\}$ — homomer arm runs at $R = n_g$ |
| $n_{\text{genes}}$ (heteromers) | $\{5, 8\}$ — two curves |
| $R$ (heteromers, per curve) | up to 15, log-spaced from $n_g$ |
| Subsampling strategy | uniform_random, cascading |
| Entropy estimator | shannon (training loss); both shannon and renyi logged as measurements |
| Environment sample | $i = 1 \ldots 30$ from prior in §4 |
| Seed | 3 per (env, condition) |

Batch sizes resolved via §5.4 — auto-adjusted to $R$.

Estimated total runs: roughly $(6 + 2 \times 8) \times 2 \times 30 \times 3 \approx 4000$
runs. Each is small ($R \leq 15$, batch sizes modest), so total wall time should be
manageable.

### 6.3 Robustness SI Pilot (in parallel, lightweight)

One single-axis-at-a-time sweep at fixed $R = 12$ (or wherever the main signal is
cleanest in 6.2), varying each parameter from §4 individually around the prior
centre. 5 axes × 5 values × 3 seeds ≈ 75 runs. Confirms that the prior choices are
sensible and identifies any outlier parameter axes before the post-$R$-15 campaign.

### 6.4 Later Extension (after first shot lands)

Once the first-shot figure is committed and the Rényi-vs-Shannon validation passes,
extend the homomer arm to $n_{\text{genes}} \in \{30, 50, 100, 300, 1000\}$ and the
heteromer arm to $R$ up to ~1000 (subsampled). Rényi only at that point. This is
where the headline number ("25 heteromer genes match X% of 1000 homomer genes")
comes from.

## 7. Diagnostics to Monitor Throughout

Things that need to be true for the narrative to hold, and which should be checked
incrementally rather than discovered at submission time.

### Optimization Health

- **Kernel-regime diagnostic.** Histogram of $\|v_u - v_\ell\|^2 / \lambda^2$ at end
  of training. Most mass should sit in $[0, 4]$ where gradients are non-trivial. Mass
  piling up above ~5 means saturation is killing the run; the corresponding environment
  is outside the regime and should be flagged or rejected from the prior.
- **`softplus(max_energy_u_raw)` post-training.** Should remain bounded above zero
  across units. If the optimizer is driving this toward zero across the board, the
  kernel is collapsing into uniformity — sign of a bad regime or bad init.
- **Pre-sigmoid spread.** Standard deviation of `ln_sum_terms` across the batch.
  Should be of order $T_{\text{final}}$, not orders of magnitude away. If it's much
  smaller, every receptor is firing or silent identically; if much larger, temperature
  annealing isn't catching up.

### Entropy Estimator Health

- **Shannon-vs-Rényi gap.** Logged on every first-shot run. A consistent gap across
  configurations means Rényi is a reliable proxy; a config-dependent gap means
  Rényi can't be trusted for the extension to larger $R$.
- **Rényi-2 collision count.** For each evaluation batch, log the number of pairwise
  collisions actually observed. If it drops below ~30, the Rényi estimate is
  unreliable — auto-batch-size should catch this, but it's worth verifying.
- **Miller-Madow gap.** $\hat{H}_{\text{MM}} - \hat{H}_{\text{plugin}}$ at evaluation
  time. Small means trustworthy; large means undersampled. Compare to $\log_2 B$ as
  the trivial ceiling.
- **$\hat{K} / 2^R$.** Fraction of binary codeword alphabet observed. Distinguishes
  "array is producing few patterns" (low fraction, intrinsic to the array) from
  "we didn't sample enough" (low fraction relative to $\log_2 B$).

### Comparison Validity

- **Theoretical-ceiling check.** At $R = 10, 12, 15$, brute-force the maximum
  possible Shannon entropy for the canonical environment and confirm the trained
  array reaches it (within sampling error). This is the existence proof for the
  calibration claim.
- **Warm-start gain.** Compare warm-started vs cold-started runs at matched $R$.
  Warm-starting should improve final MI for both arms; if the gain is asymmetric
  between homomer and heteromer, that's an optimizer bias worth understanding before
  it affects the headline.
- **Subsampling sensitivity.** Uniform-random vs cascading curves should be close.
  Persistent divergence (say > 1 bit) at multiple $R$ means the conclusion depends on
  the subsampling choice, which is a real finding but reframes the paper.

### Robustness Bands

- **Band width.** Across the 30 sampled environments at fixed $R$, the spread of MI
  should be narrow relative to the homomer-heteromer gap. If bands are wider than the
  gap, the comparison isn't environment-robust and the prior needs tightening (or the
  story needs rewriting).
- **Outlier environments.** Individual environment samples that produce wildly
  different curves are worth investigating before the figure is finalized — they
  often reveal a parameter combination on the edge of a failure mode.

## 8. Decision Points

Before committing to the post-$R$-15 extension, three things need to be confirmed by
the first-shot data:

1. **Does the heteromer/homomer gap exist at small $R$?** If at $R = 15$ with
   $n_g^{\text{hetero}} = 5$ heteromers don't show any per-gene efficiency advantage,
   the headline needs reconsidering before scaling up. The first-shot figure itself
   answers this.

2. **Does Shannon ≈ Rényi within the calibrated zone?** If the gap between
   estimators is small and config-independent, the extension to larger $R$ using
   Rényi alone is well-founded. If not, blocked Shannon may be needed and compute
   budgets need re-evaluating.

3. **Do the bands close?** If 10–90 percentile bands at $R = 15$ are wider than the
   homomer-heteromer gap, the environmental prior is too loose. Tighten before
   running the larger campaign.

Once these three confirm, the post-$R$-15 extension in §6.4 is the long pole and
everything else is verification.
