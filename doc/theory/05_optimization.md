# 5. Discrete Optimization & The Rényi Entropy Proxy

To computationally optimize the receptor array, our objective is to tune the chemical affinities of each genetic subunit—specifically their baseline energies $E_{\text{base}}^{(u)}$ and their coordinate vectors $\textbf{v}_u$ in the latent space. We want the full array to produce the widest possible variety of distinct activation patterns, effectively maximizing the mutual information between the environment and the array's responses.

## 5.1 The Optimization Challenge
Because downstream neural processing relies on highly thresholded, discrete signals (e.g., firing vs. silent), we want to optimize for discrete diversity. However, gradient descent cannot optimize sharp, discrete step-functions (Heaviside functions) because their gradients are zero almost everywhere. 

To bypass this, we use **Continuous Relaxation** (or "Soft Binning"). We approximate the discrete probability that a receptor $r$ fires by using a temperature-scaled Sigmoid centered around the receptor's $EC_{50}$ threshold:

$$p_r = \sigma\left( \frac{\ln(c) - \ln(EC_{50}^{(r,\ell)})}{T} \right)$$

As $T \to 0$, this provides the smooth, continuous gradients necessary for backpropagation while faithfully representing discrete states.

## 5.2 The True Shannon Entropy Limit
Ideally, we want to maximize the Shannon Entropy of the array's joint probability distribution:
$$H(\mathcal{A}) = -\sum_{\mathcal{A}} P(\mathcal{A}) \log_2 P(\mathcal{A})$$

However, computing this requires evaluating the probability of all $2^R$ possible discrete states (for an array of $R$ binary receptors). For an array of 26 receptors, this means computing 67 million states per gradient step, which is computationally impossible and causes Out Of Memory (OOM) errors on modern GPUs.

## 5.3 The Rényi Entropy Proxy
To solve the dimensionality explosion, we rely on a computationally tractable surrogate: the **Rényi Entropy** (specifically of order 2). 

The Rényi entropy of order 2 is defined as:
$$H_2(\mathcal{A}) = -\log_2 \sum_{\mathcal{A}} P(\mathcal{A})^2$$

The immense computational advantage of Rényi entropy is that it can be computed directly using the pairwise interactions (the Information Potential) of the array's activity across a batch of ligands, rather than requiring the full joint state distribution.

In our simulation pipeline, we seamlessly choose between optimizing the **Exact Entropy** (Shannon) or the **Proxy Entropy** (Rényi) depending on the configuration. By directly optimizing this Rényi proxy, we push the array towards maximum diversity without relying on split heuristic penalties. Because evaluating the information potential scales quadratically $\mathcal{O}(R^2)$ rather than exponentially, it makes the optimization of large arrays computationally feasible while still effectively maximizing the overall information capacity.