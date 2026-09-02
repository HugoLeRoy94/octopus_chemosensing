Now I can give you a precise answer instead of the hand-wave. Both estimators have the identical skeleton (see bin_loss.py:244 and bin_loss.py:306):

$$H \approx H_{\text{cond}} - \frac{1}{B}\sum_i \log_2 m_i, \qquad m_i = \frac{1}{B}\sum_j K(i,j)$$

where $A_{ir}$ = probability receptor $r$ is active given sample $i$, $p_i$ = the product-Bernoulli component for sample $i$, $B$ = batch size, $R$ = number of receptors, $H_{\text{cond}}$ = mean per-component entropy, and $K(i,j)$ is a pairwise affinity kernel with $K(i,i)=1$. The two bounds differ only in $K$.

1. Why one kernel is bounded and the other is not
Lower uses the Bhattacharyya coefficient:
$$K(i,j)=\prod_r\big(\sqrt{A_{ir}A_{jr}}+\sqrt{(1-A_{ir})(1-A_{jr})}\big)$$
Each receptor factor is the overlap of two Bernoulli distributions. By Cauchy-Schwarz each factor lies in $[0,1]$, equal to 1 only when $A_{ir}=A_{jr}$. A product of numbers in $[0,1]$ stays in $[0,1]$. So it is a similarity: 1 for identical components, 0 for perfectly distinguishable ones. Bounded by construction.

Upper uses $K(i,j)=e^{-\mathrm{KL}(p_i|p_j)}$ with
$$\mathrm{KL}(i|j)=\sum_r\Big[A_{ir}\log\tfrac{A_{ir}}{A_{jr}}+(1-A_{ir})\log\tfrac{1-A_{ir}}{1-A_{jr}}\Big]$$
KL is a log-ratio, in $[0,\infty)$. It diverges when $p_j$ puts near-zero probability where $p_i$ puts mass ($A_{jr}\to0$ with $A_{ir}>0$ sends a term to $+\infty$). So $e^{-\mathrm{KL}}$ can be arbitrarily close to 0.

The root difference is: an overlap integral (a bounded inner product) versus a log-ratio (diverges on support mismatch). That one fact drives everything below.

2. What "tight" actually guarantees
I overstated it. The only hard theorem (Kolchinsky-Tracey 2017) is the sandwich: $\text{KT}{\text{lower}} \le H(s) \le \text{KT}{\text{upper}}$, for any $B$, exact up to the eps clamp. There is no standalone guarantee that the lower bound equals $H$. Both bounds provably coincide only in the well-separated regime (every pair of components either near-identical or near-disjoint), where the bracket collapses and "tight" is literal. In the partially-overlapping regime the honest statement is just the bracket. What makes the lower bound the one we report is not certified tightness but that it converges (next point) and every value it returns is a certified floor on $H$.

3. Why raising B pushes the upper bound up but not the lower
The key: $K(i,i)=1$, so $m_i \ge 1/B$ always. The diagonal alone gives $m_i=1/B$, i.e. $-\log_2 m_i = \log_2 B$. That is the trivial ceiling: every sample counted as its own distinct symbol.

Raising $B$ does not change $\mathrm{KL}(p_i|p_j)$ between two fixed components. It only adds points. The question is whether the off-diagonal neighbors keep $m_i$ above $1/B$.

Lower: for a typical distinct pair, each receptor factor is $<1$ but $O(0.8\text{–}0.99)$, and the product over $R$ stays well above zero. Point $i$ keeps many soft neighbors, so $m_i \to \mathbb{E}_z[\text{BC}(p_i,p_z)]$, a fixed constant in $(0,1]$. $-\log_2 m_i$ converges. The bound plateaus and more samples only shrink variance (law of large numbers on a bounded kernel).

Upper: $\mathrm{KL}$ is a sum over $R$ receptors of positive terms, so between two genuinely different components it is large and $e^{-\mathrm{KL}}$ underflows to $\approx 0$. Point $i$ then has no effective neighbor but itself, $m_i \approx 1/B$, and $-\log_2 m_i \approx \log_2 B$. As $B$ grows this climbs without bound, so $H_{\text{upper}} \approx H_{\text{cond}} + \log_2 B$, until it hits the clamp at $R$ (bin_loss.py:306). That clamp is the flat top you see.

So it is not that KL grows with $B$. The KL kernel is razor-sharp (its effective bandwidth shrinks to zero in $R$-dimensional space), so it sees every other sample as infinitely far and the estimator degenerates to the "$B$ distinguishable points" ceiling $\log_2 B \to R$. The Bhattacharyya kernel is broad enough to keep real neighbors, so it measures the true entropy and saturates.

Takeaway for the figure: report the lower bound (converges, certified floor). The upper is a certified ceiling but degrades to $\log_2 B \to R$ as $B$ grows, so it is uninformative at large $B$. The widening bracket is the upper bound degrading, not added uncertainty. One caption line saying so keeps a reader from over-reading the upper's climb.

Here is the single argument, focused on degradation.

Fix the target. Let $H^*$ be the true mixture entropy, the entropy of the population of components $p(z)$. It is a fixed number, independent of $B$, and bounded by $R$ (=number of receptors, so at most $R$ bits). The KT bounds do not bracket $H^*$ directly, they bracket $H_B$, the entropy of the equal-weight mixture of the $B$ components actually in the batch. Standard fact: $H_B \to H^*$ as $B$ grows, and once $B$ exceeds the effective number of distinguishable states $2^{H^*}$, $H_B$ plateaus at $H^*$. Your observed lower-bound plateau is the direct evidence that $H_B$ has already saturated.

Now watch the overshoot. Define the upper bound's excess over the true target:
$$\text{slack}{\uparrow}(B) = \text{KT}{\text{upper}}(B) - H_B.$$
A good upper bound has small slack. From the mechanism, in the distinguishable regime $\text{KT}{\text{upper}}(B) \approx H{\text{cond}} + \log_2 B$ (each point's neighbor mass collapses to the diagonal $1/B$, so $-\log_2 m_i \approx \log_2 B$). Therefore
$$\text{slack}{\uparrow}(B) \approx \big(H{\text{cond}} + \log_2 B\big) - H^* ;\xrightarrow{B\to\infty}; \infty,$$
capped only by the clamp at $R$. The target $H_B$ is flat, but the bound climbs, so the excess diverges. That is the degradation, stated quantitatively: the amount by which the upper bound overshoots the true entropy grows like $\log_2 B$.

Why this is degradation and not just "a different number": an upper bound is only as useful as it is tight. In the limit the upper bound converges to $R = $ "$H(s) \le R$ bits", which is the trivial bound you could write without looking at any data (R binary receptors carry at most R bits). So adding samples drives the bound from an informative value toward the content-free maximum. More data, less information from this bound. By contrast $\text{slack}{\downarrow}(B) = H_B - \text{KT}{\text{lower}}(B) \to H^* - L^*$, a constant (both $H_B$ and the lower bound plateau), so the lower bound's error stays bounded. Only the upper bound degrades.