## training loop

The training loop consist in the following step, repeted for a number of epochs fixed in the config file.

### 1. temperature annealing
Considering the relation between the concentration and the probability of activation for a single ligand:

$$p_1 = \sigma\left( \frac{\ln(c) - \ln(EC_{50})}{T} \right)$$

In a mixture, the contributions of all ligands are combined via a logsumexp before passing through the sigmoid:

$$p = \sigma\left( \frac{\ln\!\sum_\ell \exp\!\left(\ln c_\ell - \ln EC_{50,\ell}\right)}{T} \right)$$

We see that temperature makes the sigmoid more or less sharp. Increasing $T$ makes it less sharp, and tends to increase the entropy because the receptor becomes sensitive to a wider diversity of ligands in a non-deterministic way. The main reason we don't use infinitely sharp $T$ is because it leads to gradients taking meaningful values only near the threshold, and vanishing everywhere else. As a result, we use annealing, where the start temperature is hard-coded to 1 and annealed to the final temperature set by the config.

-> Is it good to hard code the initial temperature to 1. ? probably not.
-> What should be the initial/final temperature as a function of the family spread. For instance, if I have 100 ligands, in a high dimensional space, with a single, but wide family, their distance is so big, that any receptor is either activated or not activated and the gradient is mostly 0 everywhere.
    -> What would be the good measurement to check if we are in the right regime ?

### 2. sample batches

At each timestep batch_size (fixed in the config file) mixtures are drawn, by pulling every ligands with a probability set by the p_presence array, and their concentration is drawn independantly in each mixture.

**what is not sampled each time, which means what is fixed :**

- The latent vector of the ligands

- the number of families, their position, and spread

- $\lambda$: the affinity length scale, which controls how quickly the binding energy saturates with distance in latent space: $E = E_\text{base} + E_\text{max}\left(1 - e^{-d^2/\lambda^2}\right)$. The slope of the activation sigmoid is controlled by $T$, not $\lambda$.


**what is sampled:**

- The mixtures array.

- The "binding-noise": the additional noise when computing the specific EC50

- The activation array : which goes through the activation probability $p_1$


**what is optimized:**

- the latent vector of the receptors

- Emax of the receptors

- Emin of the receptors

In the limit of low binding-noise:
-> If there is not enough ligands, the training entropy only tries to optimize over the concentration ladder.

### 3. compute activity

Given a set of parameters — energy, concentration, and mask for each ligand — the activity is computed using the sigmoid function above. For each sample (a full mixture), the activity is a single probability between 0 and 1. The soft binning of this probability into discrete bins is deferred to the loss function in step 4.

### 4. loss function

Given a probability of assignement provided by the "physics" of the system, the loss function first transform the single soft bining number into actual bin (70% chance of being activated means $P(0) = 30%$ in the first bin and $P(1) = 70%$ in the second). Followed by the computation of the entropy.

The loss then computes entropy of the full array, either the Rényi entropy, or the shannon entropy, which can be computed in block, or as full.

**Shannon entropy**:

To understand the numerical operations, we can view this function as a transformation from local marginal probabilities to a global joint distribution. The core operation is a recursive **Kronecker-like product** (outer product) that assumes independence between receptors during the construction phase.

Given $R$ receptors, each with $K$ possible states, the function constructs a joint probability tensor. If we denote the probability of receptor $r$ being in state $k$ for a sample $b$ as $p_{b,r}(k)$, the function iteratively computes:

$$P_b(x_1, x_2, \dots, x_R) = \prod_{r=1}^{R} p_{b,r}(x_r)$$

In the code, this is handled by the `unsqueeze` and `view` operations, which effectively perform:

* **Step 1:** Start with $\mathbf{P}_1 \in \mathbb{R}^{B \times K}$.
* **Step $r$:** Compute the outer product of the current accumulated joint states with the next receptor:

$$\mathbf{P}_{joint} = \mathbf{P}_{accumulated} \otimes \mathbf{P}_{next}$$


* **Result:** A flattened tensor of shape $(B, K^R)$ representing every possible combination of receptor activities.

**Batch normalization**

The function then estimates the "true" probability of each global state $a$ by averaging over the batch $B$:

$$\bar{p}_a = \frac{1}{B} \sum_{b=1}^{B} P_b(a)$$

This step collapses the batch dimension, resulting in a single vector $\mathbf{p} \in \mathbb{R}^{K^R}$ where $\sum \bar{p}_a = 1$.


Finally, the Shannon entropy $H$ is calculated using the marginalized probabilities. To ensure numerical stability, a floor $\epsilon$ is applied to the log term:

$$H(X) = - \sum_{a=1}^{K^R} \bar{p}_a \log_2(\max(\bar{p}_a, \epsilon))$$

**The Block Decomposition**

If the batch size is too large, the algorithm assumes the total system $X$ can be partitioned into $k$ disjoint subsets (blocks), $X = \{S_1, S_2, \dots, S_k\}$. Under a **block-independence assumption**, the joint entropy is approximated as the sum of the exact entropies of each block:

$$H_{blocked} \approx \sum_{i=1}^{k} H_{exact}(S_i)$$

This captures **all high-order interactions** within each block exactly, while ignoring correlations between receptors in different blocks.

To mitigate the bias introduced by ignoring cross-block correlations, the function employs a **multiple-shuffling strategy**:

1. **Random Permutation:** For each partition $j$, receptors are randomly reordered.
2. **Sequential Chunking:** Receptors are grouped into blocks of size $M$ (e.g., 12).
3. **Ensemble Mean:** The results from $N$ different random partitions are averaged:

$$\hat{H} = \frac{1}{N_{partitions}} \sum_{j=1}^{n} H_{partition, j}$$


> **Complexity Shift:** Memory scales as $O(B \cdot 2^{block\_size})$ rather than $O(B \cdot 2^R)$. For $R=60$ and $block\_size=12$, this reduces the state space from $2^{60}$ (impossible) to 5 parallel calculations of $2^{12}$ (trivial).

**Rényi entropy**:

The entropy is estimated via the pairwise collision probability — the probability that two random ligands produce the exact same array state. This avoids the $O(2^R)$ enumeration but requires large batches for a reliable estimate.

### Bound Relationship

$$H_{Rényi} \leq H_{Shannon} \leq H_{Blocked}$$

| Method | Boundary | Logic |
| --- | --- | --- |
| **Rényi ($H_2$)** | **Lower Bound** | Collision-based; biased toward the most probable "peaks" of the distribution. |
| **Shannon (Exact)** | **Ground Truth** | The actual system uncertainty (intractable for large $R$). |
| **Blocked** | **Upper Bound** | Summing independent parts always overestimates total uncertainty if hidden correlations exist. |