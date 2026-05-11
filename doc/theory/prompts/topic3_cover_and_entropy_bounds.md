# Topic 3 — Cover's Theorem Bound and the Real Maximal Entropy of the Environment

## Context

The user observes that even when latent dimension $D$ is large, their array of $N$ binary receptors often fails to reach the information-theoretic optimum of $N$ bits. The question is whether this is an *optimization* failure or a *capacity / environment* failure.

Three distinct ceilings on $H(\mathcal{A})$ have been identified:

1. **Coding ceiling:** $H(\mathcal{A}) \le N$ bits, trivially.
2. **Geometric capacity ceiling (Cover):** the number of distinct binary codes the array can produce is bounded by the number of cells in an arrangement of $N$ decision surfaces in the input space.
3. **Vocabulary ceiling:** if the environment exposes only $M$ distinct ligands, the array sees at most $M(N+1)$ distinct inputs and thus at most $M(N+1)$ distinct codes.

The user correctly pointed out that the differential entropy of the environment is **not** a valid upper bound on the discrete entropy of the array.

## What was established

### Nats vs bits
Pure unit convention: $\ln$ gives nats, $\log_2$ gives bits, $1\,\text{nat} = \log_2(e) \approx 1.443\,\text{bits}$. Use nats for calculus, bits for codes/channels.

### Why differential entropy doesn't bound discrete entropy
1. Differential entropy $h(X) = -\int p(x) \ln p(x)\,dx$ is **not unit-invariant**: $X \to \alpha X$ gives $h \to h + \ln|\alpha|$.
2. It can be **negative**: $h(\mathcal{N}(0,\sigma^2)) = \tfrac{1}{2}\ln(2\pi e \sigma^2) < 0$ for $\sigma < 1/\sqrt{2\pi e}$.
3. The discrete entropy of a quantization with bin width $\Delta$ is $H \approx h(X) - \log_2 \Delta$, which diverges as $\Delta \to 0$. There is no finite quantization-independent bound.

### The vocabulary bound (counting argument)
For any fixed ligand $\ell$, the array's response as $\ln c$ varies is monotone in $\ln c$ for each receptor (since $a_r = 1 \iff \ln c > \theta_r(\ell)$ for some threshold $\theta_r(\ell)$). Sorting the $N$ thresholds, the code is piecewise constant on $N+1$ intervals.

Therefore the image of the deterministic map $(\ell, c) \to \mathcal{A}$ has size at most $M(N+1)$, giving:
$$H(\mathcal{A}) \le \log_2[M(N+1)]$$

For $N = 26$, full saturation requires $M \ge 2^N/(N+1) \approx 2.5 \times 10^6$ distinct ligands. The concentration axis only adds $\log_2(N+1) \approx 4.7$ bits — a small additive contribution.

### Role of $\sigma_c$
$\sigma_c$ is a **support condition**, not an entropy contribution. The $N+1$ codes per ligand are only realized if $\ln c$ has enough spread to cross all $N$ thresholds with non-negligible probability. Increasing $\sigma_c$ beyond the threshold range gives no further benefit.

### The geometric capacity bound (Cover-style)
For each receptor, the firing rule $\ln c > \frac{1}{k_\text{sub}}\sum_u (E_\text{base}^{(u)} + \|v_u - v_\ell\|^2)$ defines a quadratic surface in $(v_\ell, \ln c)$-space. By lifting (replacing $\|v_\ell\|^2$ with an auxiliary coordinate), this becomes a hyperplane in $\mathbb{R}^{D+2}$.

The number of cells in an arrangement of $N$ hyperplanes in $\mathbb{R}^{D+2}$ is bounded by:
$$\#\text{cells} \le \sum_{k=0}^{D+2}\binom{N}{k}$$

For $\#\text{cells} \ge 2^N$, we need $D+2$ comparable to $N/2$, not the originally-claimed $N/\log_2 N$. Specifically:

- $D \ll N/2$: capacity $\approx (D+2)\log_2 N$ bits — polynomial in $N$, far below saturation.
- $D \approx N/2$: capacity $\approx N-1$ bits — near saturation.
- $D \ge N$: full shattering, all $2^N$ codes achievable.

For $N=26$, full saturation requires $D \approx 12$–$13$, not $D \approx 5$.

### The combined bound
$$H(\mathcal{A}) \le \log_2 \min\left\{ 2^N,\; \sum_{k=0}^{D+2}\binom{N}{k},\; M(N+1) \right\}$$

## Open questions for the next discussion

### 1. Sharpness of the Cover bound
The Cover bound is achieved only by hyperplane arrangements *in general position*. Real receptor arrangements are far from general position (units are shared between receptors via combinatorial subunit sharing in the heteromer construction). What is the *effective* dimension after accounting for this constraint?

Concretely: with $n_\text{unit}$ proteins and $k_\text{sub} = 5$, there are at most $\binom{n_\text{unit} + k_\text{sub} - 1}{k_\text{sub}}$ distinct heteromers, but their thresholds are *correlated* because they share sub-units. This correlation reduces the effective number of independent decision surfaces below $N$.

### 2. The "general position" question for quadratic surfaces
The lifting argument $\|v_\ell\|^2 \to z$ replaces the quadratic surface with a hyperplane. But this hyperplane lives in a $(D+1)$-dim subspace of $\mathbb{R}^{D+2}$ (the lifted coordinate is determined by the original ones). Does Cover's theorem still give the same count, or is there a tighter bound?

### 3. Phase-diagram experimental design
- Vary $D$ at fixed $M$ (large): expect saturation transition near $D \approx N/2$.
- Vary $M$ at fixed $D$ (large): expect saturation transition near $M \approx 2^N/N$.
- Vary both simultaneously: expect a 2D phase diagram with two "walls."
- $\sigma_c$ should be set "large enough" once and not varied.

### 4. Diagnostic for optimization vs capacity failure
At the end of training, count the number of distinct binary codes the array produces over a large evaluation set. Compare to:
- $2^N$ (full saturation).
- $\sum_{k=0}^{D+2}\binom{N}{k}$ (geometric capacity).
- $M(N+1)$ (vocabulary bound).

If observed $\ll$ all three, optimization-limited. If observed $\approx \min$ of the bounds, capacity-limited.

### 5. Does the bound change with continuous (non-binary) receptors?
The Soft Histogram in Section 5.3 of `04_discrete_information.md` uses $K$ bins per receptor. The vocabulary bound becomes $M \cdot (\text{some function of } K, N)$ — likely $M \cdot (KN + 1)$ or similar, depending on how monotonic the activity is in $\ln c$. Worth working out.

### 6. Mixture environments
Everything above assumes a single ligand per sniff. For mixtures of $\le m$ ligands from a vocabulary of $M$, the effective vocabulary becomes $\sum_{k=0}^{m}\binom{M}{k}$ (presence patterns) — combinatorially larger, easily exceeding $2^N/N$ for moderate $m$ and $M$. This may explain why the user's current simulations work at all: the Bernoulli mixture mask implicitly enlarges the effective ligand vocabulary.

### 7. Is the bound $M(N+1)$ tight, or can it be improved with knowledge of the concentration distribution?
If concentrations are correlated across ligands (e.g., ripe-fruit volatiles co-occur at correlated concentrations), the effective vocabulary is smaller than $M$ multiplied independently by the concentration axis. The bound is "honest" in the worst case but possibly loose in practice.

## Sanity-check experiments to propose

1. **Pure vocabulary test:** Fix $D$ huge (e.g., $D = 2N$), vary $M \in \{2^k : k = 5, \ldots, N+5\}$. Verify saturation transition at $M \approx 2^N/N$.

2. **Pure capacity test:** Fix $M$ huge, vary $D \in \{1, 2, \ldots, N\}$. Verify saturation transition at $D \approx N/2$.

3. **Diagnostic: count distinct codes** in the final array's response to a large evaluation batch, and compare to the analytic bounds.

4. **Single-ligand vs mixture comparison:** confirm whether the mixture mask is what allows the current simulations to reach high $H(\mathcal{A})$ despite small $M$.
