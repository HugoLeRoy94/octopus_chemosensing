# objective

Our main objective is to compare the relative efficiency of heteromerization to gene expansion. To do so, we first want to charaterize an array of homomers systematically, which means independently of the environment's parameters. The current belief, is that there exist a regime where the environment is sufficiently complex to generate enough diversity of response so that the array of homomers is optimal in the information theoric point of view : 1 bit per receptors. This lower-bound of complexity compete with an upper bound due to optimization difficulties.

## Cover's theorem

There is a theorem that give the number of sets that are linearly separable as a function of the dimension of the space. This seems to apply to our problem, because our array is essentially a classifier with hyperplans. This is not quite linear, but to some extend it should work. As a result, there is a lower bound of the dimension of the latent space below which the entropy of the array is limited by geometrical factors.

## Impact of the number of ligands



## Impact of the spread of families

As the spread of families $\sigma_f$ increase, the distance between ligands increas as $\sigma_d^d$ where $d$ is the dimension of the latent space. The typical minimal distance between a receptor and a ligand leading to a, at least partial activation is given by 

[To complete]

Which means that: the more spread are the families, the higher the initial temperature of the activation profile must be. Additionnally, the $\lambda$ could be higher, but, not sure what it would do compare to adapting $E_\text{max}$.

## What about the LR Scheduler?

The temperature annealing described above is always active. Separately, a cosine annealing learning-rate scheduler can be enabled via config. In general it seems to perform worse than a constant learning rate. Why this is the case is not fully understood.


## Evaluating the entropy estimate: sampling diagnostics

The receptor array produces binary codewords — R-bit patterns encoding which receptors fired for a given mixture. With $B$ evaluation samples, the joint entropy is estimated by counting how often each pattern appears. This count-based estimate is systematically biased downward: if $B$ is not much larger than $2^R$, many codewords are never observed, and their probability mass is missing from the estimate.

Several diagnostics are logged at evaluation time to make this bias visible.

### The plug-in estimator

The naive estimator counts each codeword's frequency and plugs into Shannon's formula:

$$\hat{H}_\text{plugin} = -\sum_k \hat{p}_k \log_2 \hat{p}_k$$

This underestimates the true entropy whenever the sample budget $B$ is insufficient to cover the full alphabet of $2^R$ codewords.

### Miller-Madow correction

A simple bias correction adds a term proportional to the number of distinct codewords observed:

$$\hat{H}_\text{MM} = \hat{H}_\text{plugin} + \frac{\hat{K} - 1}{2\, B\, \ln 2}$$

where $\hat{K}$ is the number of distinct codewords seen in the $B$ samples. This correction is small when sampling is dense relative to the alphabet, and large when the estimate is severely undersampled.

### Sampling regime diagnostics

Three additional quantities are logged to characterize the sampling regime:

- $\hat{K}$: number of distinct codewords observed. If $\hat{K} \approx B$, the sample budget is exhausted — every sample is unique and the estimate is unreliable.
- $\hat{K}/2^R$: fraction of the full binary alphabet that was observed. A value close to 1 means dense coverage; a small value means the array has low entropy (few patterns are ever produced).
- $\log_2 B$: the trivial upper bound on what any plug-in estimator can report. With $B$ samples, you cannot distinguish more than $B$ distinct probabilities, so $\hat{H}_\text{plugin} \leq \log_2 B$ always. If the reported entropy is close to this ceiling, the reading is dominated by the sample limit, not by the true array entropy.

### How to read the diagnostics in practice

| Situation | Signature |
|---|---|
| Trustworthy estimate | $\hat{H}_\text{MM} \approx \hat{H}_\text{plugin}$, $\hat{K} \ll B$, $\hat{H}_\text{plugin} \ll \log_2 B$ |
| Undersampling artifact | $\hat{H}_\text{MM} \gg \hat{H}_\text{plugin}$, $\hat{K} \approx B$, $\hat{H}_\text{plugin} \approx \log_2 B$ |
| Genuinely low entropy | $\hat{K}/2^R \ll 1$ regardless of $B$ — the array itself produces few distinct codewords |

For example, with $R = 20$ receptors and $B = 2^{12} = 4096$ evaluation samples, $\log_2 B = 12$ bits, far below the 20-bit maximum. Any entropy reading above 12 bits is impossible, so undersampling is guaranteed. With $B = 2^{20}$ samples, $\log_2 B = 20$ bits matches the maximum, and the estimate can in principle cover the full range.