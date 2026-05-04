## 4. Information Produced by an Array of Receptors

We now consider an array of $N$ receptors, and write the activity of this array as $\mathcal{A} = (a_0,a_1, \cdots , a_N)$. The goal of this section is to design an algorithm to optimize the the mutual information between the activity of the array and ligand's identity and/or concentration for a given environment:
$$
I(\mathcal{A};(c,\ell)) = H(\mathcal{A}) - H(\mathcal{A} | (c,\ell)),
$$
We consider the mapping between affinities, and activity as deterministic, noiseless: $H(A | (c,\ell)) = 0$.
Thus we want to maximize 
$$
H(\mathcal{A}) = -\sum_\mathcal{A} p(\mathcal{A}) \log[p(\mathcal{A})] = \mathbb{E}\left[ \log[p(\mathcal{A})] \right].
$$

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