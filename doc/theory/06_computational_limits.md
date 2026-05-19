# 6. Computational Limits and Architectural Decisions

Optimizing combinatorial receptor arrays in high-dimensional environments frequently encounters the "Curse of Dimensionality." This document outlines the physical GPU memory limits of the simulation and the architectural fallbacks implemented to avoid Out Of Memory (OOM) errors.

## 6.1 Scaling & Memory Footprints
The simulation code is written in PyTorch and runs on NVIDIA A100 GPUs. The memory footprint scales according to these primary operations:

1. **Ligand Generation & Distances:** $\mathcal{O}(B \cdot U \cdot D)$
   - Batch Size ($B$), Number of Units ($U$), Latent Dimension ($D$). 
   - *Impact:* Extremely lightweight. A tensor of shape $(2000, 26, 10)$ takes $\approx 2$ MB. The latent dimension ($D$) and Number of Families ($F$) are virtually "free" to scale.

2. **Receptor Physics (Combinatorics):** $\mathcal{O}(B \cdot R \cdot k_{sub})$
   - Number of Receptors ($R$), Subunits per receptor ($k_{sub}=5$). 
   - *Impact:* Scales linearly. Even for $10,000$ combinatorial receptors, the energy tensor takes $\approx 400$ MB.

### The Loss Function Memory Bottleneck
The choice between Exact Shannon Entropy and the Rényi Proxy fundamentally changes the underlying PyTorch tensor operations, serving as the primary constraint on array scaling:

3. **Exact Loss Function Bottleneck (`DiscreteExactLoss`):** $\mathcal{O}(B \cdot 2^R)$
   - **The Math:** To compute the exact Shannon entropy, the continuous relaxation outputs an activation probability tensor $\mathbf{A}$ of shape $(B, R)$. We must compute the joint probability of all $2^R$ possible discrete states across the array. This is performed via a tensor product across the receptor dimension, physically allocating a massive tensor of shape $(B, 2^R)$. We then average across the batch to yield the empirical probability of each joint state, and sum the log-probabilities.
   - *Impact:* Scales exponentially. For $R=10$ receptors and $B=2000$, the tensor is $2000 \times 1024$ ($\approx 8$ MB). However, for $R=26$, the tensor becomes $2000 \times 67,108,864$, which requires hundreds of gigabytes. Evaluating exact Shannon entropy is strictly impossible for large arrays (typically failing for $R > 15$).

4. **Proxy Loss Function Bottleneck (`DiscreteProxyLoss`):** $\mathcal{O}(R^2)$
   - **The Math:** To bypass the exponential state space, the Rényi entropy proxy evaluates the diversity of the array using pairwise interactions. Given the activation matrix $\mathbf{A}$ of shape $(B, R)$, the core operation is a matrix multiplication $\mathbf{A}^T \mathbf{A}$. This multiplies an $(R, B)$ tensor by a $(B, R)$ tensor, yielding a dense pairwise covariance/repulsion matrix of shape $(R, R)$. We then penalize the off-diagonal elements of this matrix to forcefully orthogonalize the receptors.
   - *Impact:* **This is the primary memory bottleneck during standard training.** Because it operates in $\mathcal{O}(R^2)$ space, it completely avoids the $2^R$ explosion. If $R=10,000$, the resulting $(R, R)$ pairwise matrix takes $\approx 400$ MB, making the optimization of massive arrays highly manageable.

## 6.2 Key Algorithmic Decisions

### Dynamic Quadrature Fallback
When calculating the expected activation of a receptor over a Gaussian ligand family, high accuracy traditionally requires a dense Gauss-Hermite quadrature grid. However, grid integration scales exponentially with the latent dimension $\mathcal{O}(grid^{D})$. 
To prevent OOM errors in high-dimensional latent spaces, the `physics.py` module tracks the grid size dynamically. Whenever the integration grid would exceed $100,000$ points, the simulation automatically falls back to a mean-energy approximation. 

## 6.3 Rule of Thumb for Scaling
If the simulation encounters a CUDA Out Of Memory error, adjust parameters in this order:
1. **Reduce Batch Size ($B$)**: Directly impacts almost all tensor allocations.
2. **Reduce Number of Receptors ($R$)**: Relieves the $\mathcal{O}(R^2)$ penalty bottleneck or the $\mathcal{O}(B \cdot 2^R)$ exact entropy bottleneck.
3. *Note:* You do not need to reduce the latent space dimension ($D$) or the number of families ($F$), as they contribute negligibly to the training memory footprint.


## 6.4 Example

Given the parameters: $B = 2^{20}$, $L = 100$, $U = 26$, $R = 20$, $D = 20$, $K = 2$ (binary). All in fp32.

| Tensor | Shape | Size | Note |
|---|---|---|---|
| Ligand coords | $(B, L, D)$ | 8.4 GB | per-sample positions |
| Concentrations | $(B, L)$ | 4.2 GB | per-sample ligand concentrations |
| Presence mask | $(B, L)$ | 4.2 GB | which ligands present |
| **`diff` broadcast** | $(B, L, U, D)$ | **218 GB** | the OOM site |
| Squared distances $\|v_u - v_\ell\|^2$ | $(B, L, U)$ | 10.9 GB | after $\sum_D$ |
| $\ln EC_{50}^{(u,\ell)}$ | $(B, L, U)$ | 10.9 GB | per-unit, per-ligand |
| Receptor energy | $(B, L, R)$ | 8.4 GB | after $\frac{1}{k_\text{sub}}\sum$ over the $k_\text{sub}$ units in each receptor |
| Logsumexp score | $(B, R)$ | 83.9 MB | after mixture aggregation |
| Activation $p$ | $(B, R)$ | 83.9 MB | sigmoid output |
| Soft-bin assignment | $(B, R, K)$ | 167.8 MB | for $K=2$ |
| Activation grad (backward) | $\sim 2\times$ forward | $\sim 25$ GB | autograd retains for backward |

Total *peak* forward-only memory: roughly $10.9 + 10.9 + 8.4 \approx 30$ GB.  Add the autograd graph (which keeps activations for backward) raising to $\sim 60$–$80$ GB.