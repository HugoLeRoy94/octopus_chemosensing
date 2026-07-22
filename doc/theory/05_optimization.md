# 5. Discrete Optimization & Entropy Estimators

To computationally optimize the receptor array, our objective is to tune the chemical affinities of each genetic subunit. Three per-unit parameters are jointly learned by gradient descent:

* $E_{\text{base}}^{(u)}$: the open-state energy at the optimal (zero-distance) ligand, which sets each unit's baseline EC50.
* $E_{\text{max}}^{(u)}$: the maximum extra energy cost for a fully mismatched ligand, which controls selectivity breadth (constrained $> 0$).
* $\mathbf{v}_u \in \mathbb{R}^D$: the coordinate of the unit in the chemical latent space, which determines which ligands are "close" (high affinity) versus "far" (low affinity).

The global affinity length scale $\lambda$ is a fixed hyperparameter, not learned. Together, these parameters define the saturating affinity kernel described in Section 3.1. We want the full array to produce the widest possible variety of distinct activation patterns, effectively maximizing the mutual information between the environment and the array's responses.

## 5.1 The Optimization Challenge
Because downstream neural processing relies on highly thresholded, discrete signals (e.g., firing vs. silent), we want to optimize for discrete diversity. However, gradient descent cannot optimize sharp, discrete step-functions (Heaviside functions) because their gradients are zero almost everywhere. 

To bypass this, we use **Continuous Relaxation** (or "Soft Binning"). We approximate the discrete probability that a receptor $r$ fires by using a temperature-scaled Sigmoid centered around the receptor's $EC_{50}$ threshold:

$$p_r = \sigma\left( \frac{\ln(c) - \ln(EC_{50}^{(r,\ell)})}{T} \right)$$

As $T \to 0$, this provides the smooth, continuous gradients necessary for backpropagation while faithfully representing discrete states.

## 5.2 The True Shannon Entropy Limit
Ideally, we want to maximize the Shannon Entropy of the array's joint probability distribution:
$$H(\mathcal{A}) = -\sum_{\mathcal{A}} P(\mathcal{A}) \log_2 P(\mathcal{A})$$

However, computing this requires evaluating the probability of all $2^R$ possible discrete states (for an array of $R$ binary receptors). For an array of 26 receptors, this means computing 67 million states per gradient step, which is computationally impossible and causes Out Of Memory (OOM) errors on modern GPUs.

## 5.3 The Collision Entropy

To solve the dimensionality explosion, we rely on a computationally tractable surrogate: the **collision entropy** (Renyi entropy of order 2, H2).

$$H_2(\mathcal{A}) = -\log_2 C, \qquad C = \sum_{\mathcal{A}} P(\mathcal{A})^2$$

where $C$ is the **collision probability** — the probability that two independent draws from the array produce the exact same activation pattern.

The computational advantage is that $C$ can be estimated from pairwise comparisons across a batch of $B$ ligand exposures ($\mathcal{O}(B^2 \cdot R)$), without enumerating the $2^R$ state space.

### Training loss: minimise $C$ directly

The measurement estimator returns $H_2 = -\log_2 C$ (bits). For training, however, the loss minimises $C$ itself rather than $-H_2$:

$$\mathcal{L}_{\text{collision}} = C = \exp(\text{log\_mean\_coll\_prob})$$

This is equivalent ($\arg\min C = \arg\max H_2$) but removes the $1/C$ factor from the gradient of $-\log C$, which blows up in late training when $C$ becomes small. The reported measurement value is unchanged.

## 5.3b Kolchinsky–Tracey Bhattacharyya Lower Bound (`kt`)

The collision estimator is a Rényi-$H_2$ surrogate; the **Kolchinsky–Tracey (KT)** estimator instead gives a certified **lower bound on the true Shannon** joint entropy $H(\mathcal{A})$ of the same uniform mixture $p(s) = \tfrac{1}{B}\sum_b \prod_r \text{Bernoulli}(A_{br})$, $s\in\{0,1\}^R$. It is tight in the well-separated / noiseless regime (distinct near-deterministic codewords). Reference: Kolchinsky & Tracey, *Estimating Mixture Entropy with Pairwise Distances*, Entropy 2017.

Using the Bhattacharyya affinity $\text{BC}(i,j)$ between two Bernoulli-product components:

$$H_{\text{KT}} = H_{\text{cond}} - \frac{1}{B}\sum_i \log_2\!\Big(\frac{1}{B}\sum_j \text{BC}(i,j)\Big)$$

$$H_{\text{cond}} = \frac{1}{B}\sum_b \sum_r h_2(A_{br}),\qquad h_2(a) = -a\log_2 a - (1-a)\log_2(1-a)$$

$$\log \text{BC}(i,j) = \sum_r \log\!\Big(\sqrt{A_i A_j} + \sqrt{(1-A_i)(1-A_j)}\Big)\quad(\text{natural log})$$

$H_{\text{cond}}$ is the mean per-component (conditional) entropy; the second term is the inter-component overlap. **All pairs are summed, the diagonal is kept**, and the inner sum is normalised by the **total** batch $B$. Because $\text{BC}(i,i)=1$, the self-term anchors each inner sum at $1/B$, so the resolvable-entropy ceiling is $\log_2 B$ (total batch) — contrast the collision estimator, whose adjacent-chunk / diagonal-masked scheme caps its ceiling at $\log_2(\texttt{chunk\_size})$. The two limits: $B$ distinct near-deterministic codewords give $H_{\text{KT}} \approx \log_2 B$ (+ small $H_{\text{cond}}$); $B$ identical sniffs at $A=0.5$ give a zero inter-component term and $H_{\text{KT}} = R$.

Complexity: $\mathcal{O}(B^2 R / \texttt{chunk})$ time (quadratic in total batch), $\mathcal{O}(\texttt{chunk}^2 R)$ peak memory (same scaling as the collision $(R,m,m)$ block). Implemented as two nested chunk loops (outer $i$-chunks, inner $j$-chunks spanning the whole batch) with `cat` + `logsumexp` over the full row, so gradients flow through every chunk.

## 5.4 Correlation-Aware Blocked Entropy

The blocked Shannon estimator partitions the $R$ receptors into blocks of at most `block_size` receptors, computes exact Shannon entropy within each block, and sums under a between-block independence assumption:

$$H_{\text{blocked}} = \sum_{k} H_{\text{exact}}(\text{block}_k)$$

This is always an **upper bound** on the true joint entropy $H(\mathcal{A})$, since the between-block independence assumption can only inflate the estimate. The bound is tight when between-block correlations are small.

To minimize this gap, receptors are grouped by **correlation-aware greedy clustering** rather than random partitioning. The algorithm computes the absolute Pearson correlation matrix $|\rho_{ij}|$ over the current batch (a single `x.T @ x` matmul on centred activity, with the diagonal zeroed). It then greedily seeds each block with the highest-affinity remaining pair, and grows it by repeatedly adding the receptor with the maximum summed affinity to the current block members, until the block reaches `block_size`. The procedure repeats until all receptors are assigned.

Because the partition depends on the current batch statistics, it is computed on **detached** activity (stop-gradient): the grouping does not contribute to the computational graph, and gradients flow only through the within-block Shannon entropy terms. To avoid gradient instability from abrupt partition changes, the partition is **cached** and refreshed only every `block_refresh_interval` training steps (default 50). Evaluation calls always compute a fresh partition (`use_cache=False`).

Complexity: $\mathcal{O}(R^2)$ for the correlation matrix plus $\mathcal{O}(B \cdot 2^{\text{block\_size}})$ per block for exact Shannon — unchanged from the random-partition version, with the $R^2$ term negligible in practice.

## 5.5 Blocked-Corrected Entropy

The plain blocked estimator ignores cross-block correlations entirely. The **blocked-corrected** estimator subtracts the pairwise mutual information (MI) between all receptor pairs in *different* blocks:

$$H_{\text{blocked\_corrected}} = H_{\text{blocked}} - \sum_{\substack{(i,j) \\ \text{cross-block}}} I(A_i; A_j)$$

The binary pairwise MI $I(A_i; A_j)$ is computed from the batch Gram matrix:
- $P(A_i=1, A_j=1) = (\mathbf{A}^\top \mathbf{A} / B)_{ij}$
- Marginals from column means.
- Standard four-term MI formula, vectorised over all $(R, R)$ pairs in a single pass.

Within-block joint distributions are already exact in $H_{\text{blocked}}$; the MI correction only applies to cross-block pairs (upper triangle of a boolean mask built from block assignments).

Properties:
- **Tighter than blocked:** $H_{\text{blocked\_corrected}} \leq H_{\text{blocked}}$ always.
- **Pairwise limit:** only pairwise cross-block dependencies are corrected. Higher-order cross-block structure (e.g., XOR of 3 variables across blocks) escapes the correction.
- **Differentiable:** gradients flow through both the blocked term (via soft_assign) and the MI correction (via the activity).
- **Selectable as a first-class training loss:** `entropy='blocked_corrected'` in the config.

## 5.6 Annealed Blocked -> Collision Schedule

Both the blocked Shannon estimator and the collision H2 estimator have complementary strengths: blocked Shannon provides strong, unbiased gradients early in training (fast basin-finding), while collision H2 gives a faithful lower bound on joint entropy that cannot be inflated by correlated receptors. The **annealed** loss (`entropy='annealed'`) combines both via a linear interpolation controlled by training progress:

$$\mathcal{L}_{\text{annealed}} = -\bigl[(1 - \lambda)\, H_{\text{blocked}} + \lambda\, H_{\text{collision}}\bigr], \qquad \lambda = \frac{\text{epoch}}{\text{epochs}}.$$

- **$\lambda = 0$ (epoch 0):** pure blocked Shannon — fast basin-finding.
- **$\lambda = 1$ (final epoch):** pure collision H2 — faithful lower bound.

## 5.7 Blocked-to-Corrected Schedule

An alternative annealing that stays within the Shannon family: the **blocked-to-corrected** loss (`entropy='blocked_to_corrected'`) ramps from plain blocked to blocked-corrected:

$$\mathcal{L}_{\text{b2c}} = -\bigl[(1 - \lambda)\, H_{\text{blocked}} + \lambda\, H_{\text{blocked\_corrected}}\bigr], \qquad \lambda = \frac{\text{epoch}}{\text{epochs}}.$$

This avoids the collision estimator entirely while still tightening the bound over training. A `lam_override=1.0` configuration makes it pure blocked-corrected from epoch 0.

## Naming summary

| Config string | Estimator (bits) | Training loss |
|---|---|---|
| `shannon` | Exact Shannon H | $-H$ |
| `collision` | Collision H2 = $-\log_2 C$ | $C$ (collision probability, no log) |
| `kt` | KT Bhattacharyya lower bound on Shannon $H$ | $-H_{\text{KT}}$ |
| `blocked` | Blocked Shannon (upper bound) | $-H_{\text{blocked}}$ |
| `blocked_corrected` | Blocked - cross-block MI (point estimate) | $-H_{\text{blocked\_corrected}}$ |
| `annealed` | Blocked (measurement) | $-[(1-\lambda) H_{\text{blocked}} + \lambda H_{\text{collision}}]$ |
| `blocked_to_corrected` | Blocked-corrected (measurement) | $-[(1-\lambda) H_{\text{blocked}} + \lambda H_{\text{blocked\_corrected}}]$ |

Measurement brackets: `collision` (certified lower bound) and `blocked` (certified upper bound); `blocked_corrected` is the point estimate between them.
