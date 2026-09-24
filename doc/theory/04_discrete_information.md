## 4. Information Produced by an Array of Receptors

We now consider an array of $N$ receptors, and write the activity of this array as $\mathcal{A} = (a_0,a_1, \cdots , a_N)$. The goal of this section is to design an algorithm to optimize the the mutual information between the activity of the array and ligand's identity and/or concentration for a given environment:
$$
I(\mathcal{A};(c,\ell)) = H(\mathcal{A}) - H(\mathcal{A} | (c,\ell)),
$$
For a deterministic, noiseless binary mapping, $H(A | (c,\ell)) = 0$.
Only in that limit can we instead maximize
$$
H(\mathcal{A}) = -\sum_\mathcal{A} p(\mathcal{A}) \log[p(\mathcal{A})] = -\mathbb{E}\left[ \log[p(\mathcal{A})] \right].
$$

### Stochastic binary responses (implemented with `entropy='kt_mi'`)

Write $X$ for the complete sampled input and $Y\in\{0,1\}^C$ for the output
array. `activity[b,c]` is $a_c(x_b)=P(Y_c=1\mid X=x_b)$, not an observed analog
current. Assuming independent output noise **conditional on X**,

$$P(Y=y\mid X=x)=\prod_c a_c(x)^{y_c}[1-a_c(x)]^{1-y_c},$$
$$H(Y\mid X)=\mathbb E_X\sum_c h_2(a_c(X)),\qquad
I(Y;X)=H(Y)-H(Y\mid X).$$

The conditional term is analytic for each input; no repeated response draws are
needed to train. `KTMutualInformationLoss` optimizes the KT lower bound on this MI
directly (§05). Soft sigmoids are allowed at the endpoint. Identical response
probabilities across all inputs give zero MI, including probabilities of 0.5.
This removes the entropy-only incentive to produce noise, but does not guarantee
faster convergence or prevent saturation and vanishing gradients.

Here X includes concentration and any sampled coordinate/observation-noise
realization. If the intended target is a **clean** stimulus S before that noise,
the required channel is $P(Y\mid S)=\mathbb E_{X\mid S}P(Y\mid X)$.
Shared input noise can correlate cells conditional on S, so subtracting the sum
of binary entropies at each noisy input would give $I(Y;X)$, not $I(Y;S)$.
Likewise, MI about the full sniff is distinct from MI about ligand identity alone.

### Final evaluation by output counting

`mutual_information_counting` draws independent $Y_c\sim\mathrm{Bernoulli}(a_c)$
for each new input, counts the observed **joint output vectors**, and reports
plug-in and Miller–Madow entropy estimates minus the analytic batch-average
$H(Y\mid X)$. This avoids enumerating all $2^C$ states and avoids KT's pairwise
cost. Both paths use the same probability clamp `[1e-6, 1-1e-6]`, including in
the conditional term; near hard endpoints this represents a tiny numerical noise.

Counting estimates are not certified bounds or confidence intervals. They can
be negative after subtraction because of finite-sample bias; values are retained
without clipping. The observed distinct-code fraction and sample count are saved
to help assess coverage. Miller–Madow cannot repair a mostly unseen output alphabet.
Increase the evaluation budget and check stability, especially when response
noise spreads mass over many codes. KT bounds apply to the empirical input
mixture, not automatically to the population channel.

`codeword_entropy` remains a separate diagnostic: it thresholds $a_c>0.5$ and
estimates the entropy of that **deterministic** map. Subtracting the stochastic
response entropy from this hard-code entropy would mix different channels.

### Exact grouping of identical cells (`entropy='grouped_mi'`)

Let group $j$ contain $n_j$ cells with the same abundance row $W_{c,:}$ and
shared readout parameters. Their probabilities $p_j(x)$ are identical for every
input, and their binary draws remain independent conditional on $X$. Grouping
is structural, not a clustering of approximately similar sampled responses.
Define $K_j$ as the number of active cells in group $j$:

$$P(K_j=k\mid x)=\binom{n_j}{k}p_j(x)^k[1-p_j(x)]^{n_j-k},$$
$$P(\mathbf K=\mathbf k\mid x)=\prod_j P(K_j=k_j\mid x).$$

Groups need not be independent after averaging over inputs. The estimator
therefore enumerates the **joint** count alphabet, of size
$S=\prod_j(n_j+1)$, and averages its probabilities over the input batch.
It computes $I(X;\mathbf K)=H(\mathbf K)-H(\mathbf K\mid X)$ analytically;
there are no sampled output bits and no Miller–Madow correction. This is exact
for the empirical input mixture, not an exact integration over the population
of environmental inputs. Like any $B$-component empirical-input MI, it cannot
exceed $\log_2 B$.

**Entropy is not invariant under grouping; MI is.** Given a count vector, all
$\prod_j\binom{n_j}{K_j}$ labeled response patterns are equally likely, independent
of $X$. With
$$D=\mathbb E_{\mathbf K}\sum_j\log_2\binom{n_j}{K_j},$$
$$H(Y)=H(\mathbf K)+D,\qquad H(Y\mid X)=H(\mathbf K\mid X)+D,$$
$$I(X;Y)=I(X;\mathbf K).$$

`GroupedCellMutualInformationLoss.compute_entropy` returns reconstructed **full
labeled-response entropy** $H(Y)$; its forward loss is $-I(X;\mathbf K)$.
`grouped_count_entropy` and `grouped_count_conditional_entropy` explicitly refer
to counts. `response_entropy_grouped` and
`conditional_entropy_response_grouped` include $D$. No existing receptor entropy
API is redefined, and the runner rejects grouped mode for receptor-only configs.
The same `[1e-6, 1-1e-6]` probability clamp as KT/counting is used.

Grouping itself does not require enumeration. `src/response_groups.py::CellGrouping`
provides the shared structural partition and group probabilities for exact, KT,
and sampled-count estimators. The mixture is over inputs:
$P(K)=\mathbb E_X\prod_j\mathrm{Binomial}(K_j;n_j,p_j(X))$; groups are
independent conditional on the full input, not generally marginally independent.

For two input-conditioned binomials with the same multiplicity $n$, the binomial
theorem gives their Bhattacharyya coefficient as
$[\sqrt{pq}+\sqrt{(1-p)(1-q)}]^n$. Their KL divergence is also $n$ times the
Bernoulli KL. Thus weighted group columns give exactly the same KT MI bounds as
the full binary array. `entropy='grouped_kt_mi'` uses this identity; its native
entropy remains the KT lower bound on **labeled** $H(Y)$.

`grouped_counting` instead draws $K_j\sim\mathrm{Binomial}(n_j,p_j(x))$ and
counts observed integer vectors. It estimates $H(K)$ by plug-in and Miller–Madow
frequencies and subtracts analytical $\mathbb E_X\sum_j H(K_j\mid X)$.
Individual binomial supports have $\sum_j(n_j+1)=C+J$ entries; the joint alphabet
is never allocated. Reconstruction adds
$D=H(Y\mid X)-H(K\mid X)$ to count entropy. Count and labeled estimates have
distinct metric names (§09). Negative MI estimates are retained; Miller–Madow
is a bias correction, not a guarantee when many outcomes remain unseen.

For two identical cells and equally likely stimuli giving $p=0$ and $p=1/2$,
the unclamped labeled probabilities are $(5/8,1/8,1/8,1/8)$ and count probabilities
are $(5/8,1/4,1/8)$. Thus $H(Y)=1.548795$, $H(\mathbf K)=1.298795$,
$H(Y\mid X)=1$, $H(\mathbf K\mid X)=0.75$, and both MI values equal
$0.548795$ bits. The discarded $D=0.25$ bits describe which equivalent cell fired.
The numerical clamp introduces a tiny endpoint difference in the implementation.

See §08.7 for cell entropy bounds and §09.12 for configuration and measurements.

## 5. Threshold-Based Activation and Discrete Information

In many sensory systems, downstream neural processing relies on highly thresholded, binned, or even binary signals (e.g., a neuron firing an action potential or remaining silent). Optimizing a receptor array for continuous outputs (differential entropy) yields a fundamentally different geometry than optimizing for discrete, thresholded outputs (Shannon entropy).

### 5.1 Simplifying the Heteromer $EC_{50}$

Building on the adimensional activation curve from Section 2.4, we operate in the limit where the activation of a heteromeric receptor is strictly governed by its $EC_{50}$.
The half-activation concentration for a heteromer is the geometric mean of its individual subunit affinities:

$$EC_{50}^{(r,\ell)} = \left( \prod_{u=1}^{k_\text{sub}} \tilde{K}_o^{(u,\ell)} \right)^{1/k_\text{sub}}$$

By taking the natural logarithm and substituting our energy formulation $\tilde{K}_o = \exp(E_o)$, we find that the activation threshold of a heteromer in energy space is simply the arithmetic average of its sub-units' open-state energies:

$$\ln(EC_{50}^{(r,\ell)}) = \frac{1}{k_\text{sub}} \sum_{u=1}^{k_\text{sub}} E_o^{(u,\ell)}$$

This vastly simplifies the biophysics: we can completely bypass the calculation of closed-state energies ($E_c$) and evaluate the array strictly on the open-state geometry.

### 5.2 Continuous Relaxation of Discrete Activation

A true discrete sensor acts as a Heaviside step function: it fires if $\ln(c) > \ln(EC_{50})$ and remains silent otherwise. Because step functions have zero gradients, we cannot optimize them directly using backpropagation. Instead, we use a temperature-scaled Sigmoid as a continuous relaxation:

$$p_1 = \sigma\left( \frac{\ln(c) - \ln(EC_{50})}{T} \right)$$

As $T \to 0$, this function approaches a perfect step function. Crucially, we change our interpretation of this value: $p_1$ is no longer a "continuous physical current", but rather the **probability that the thresholded receptor fires**.

### 5.3 Differentiable Binning (Soft Histogram)

To compute the exact discrete Shannon entropy for an arbitrary number of downstream activation bins $K$ (where $K=2$ is the binary case), we use a "Soft Histogram" trick.
For a given receptor, we define $K$ bin centers evenly spaced between $0.0$ and $1.0$. The probability of the receptor's activity $a$ falling into bin $k$ with center $c_k$ is computed using a temperature-scaled Softmax over the squared Euclidean distances:

$$P(\text{bin}_k) = \frac{\exp\left(-(a - c_k)^2 / T_{bin}\right)}{\sum_{j=1}^K \exp\left(-(a - c_j)^2 / T_{bin}\right)}$$

This differentiable assignment allows us to analytically compute the exact marginal Shannon entropy for each receptor, normalized to a maximum of $1.0$ unit of information by using base-$K$ logarithms:

$$H(a^r) = -\sum_{k=1}^K P(\text{bin}_k) \log_K P(\text{bin}_k)$$

By maximizing these exact 1D marginal entropies while minimizing the continuous linear covariance between pairs of receptors ($\text{Cov}(a^r, a^{r'}) \to 0$), we can efficiently optimize the full array's joint entropy without the computational explosion of calculating the full $K^N$ state space.

## 6. Exact Joint Entropy and Disentangling Information

While the covariance penalty method efficiently guides gradients, accurately evaluating the true capacity of the array requires computing the exact joint entropy.

### 6.1 Exact State Enumeration
Because individual receptors respond independently to a given ligand $\ell$, the probability of the entire array expressing a specific discrete state $\mathcal{A} = (s_1, \dots, s_N)$ is the product of the individual receptor probabilities:
$$P(\mathcal{A}) = \frac{1}{B} \sum_{\ell=1}^B \prod_{r=1}^N P(a_{r,\ell} \in \text{bin}_{s_r})$$
From this exact probability distribution, the joint Shannon entropy is explicitly calculated as $H(\mathcal{A}) = -\sum_{\mathcal{A}} P(\mathcal{A}) \log_2 P(\mathcal{A})$.
*(Note: Because this requires evaluating all $2^R$ states, it scales exponentially and strictly exceeds GPU memory for arrays where $R > 15$).*

### 6.2 Disentangling Identity vs. Concentration
A biological system may not optimize its sensory apparatus to maximize total information blindly. We decompose the total mutual information to independently measure the array's ability to encode ligand identity versus concentration.

**1. Ligand Identity (Family $F$):**
$$I(\mathcal{A}; F) = H(\mathcal{A}) - H(\mathcal{A} | F)$$
Here, $H(\mathcal{A} | F)$ represents the remaining uncertainty in the array's activity when the ligand family is known. Computationally, this is evaluated by masking the continuous relaxation assignments to isolate specific families, computing the entropy for each, and averaging across all unique families in the evaluation batch.

**2. Ligand Concentration ($c$):**
$$I(\mathcal{A}; c) = H(\mathcal{A}) - H(\mathcal{A} | c)$$
Because the concentration $c$ is drawn from a continuous Log-Normal distribution, calculating the conditional entropy $H(\mathcal{A} | c)$ requires discretization. We achieve this by sorting the environmental batch by concentration and partitioning it into $N_c$ equally sized bins. The conditional entropy is approximated as the average joint entropy of the array's activity within each discrete concentration bin.
