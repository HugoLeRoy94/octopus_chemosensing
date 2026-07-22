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
The choice of entropy estimator fundamentally changes the underlying PyTorch tensor operations, serving as the primary constraint on array scaling:

3. **Exact Loss Function Bottleneck (`DiscreteExactLoss`):** $\mathcal{O}(B \cdot 2^R)$
   - **The Math:** To compute the exact Shannon entropy, the continuous relaxation outputs an activation probability tensor $\mathbf{A}$ of shape $(B, R)$. We must compute the joint probability of all $2^R$ possible discrete states across the array. This is performed via a tensor product across the receptor dimension, physically allocating a massive tensor of shape $(B, 2^R)$. We then average across the batch to yield the empirical probability of each joint state, and sum the log-probabilities.
   - *Impact:* Scales exponentially. For $R=10$ receptors and $B=2000$, the tensor is $2000 \times 1024$ ($\approx 8$ MB). However, for $R=26$, the tensor becomes $2000 \times 67,108,864$, which requires hundreds of gigabytes. Evaluating exact Shannon entropy is strictly impossible for large arrays (typically failing for $R > 15$).

4. **Collision Entropy (`'collision'`):** $\mathcal{O}(B^2 \cdot R)$ per chunk
   - **The Math:** The collision (Rényi H2) estimator computes $\log P(\text{collision}) = \sum_r \log P_r(\text{collision})$ via pairwise matching across samples. The binding matrix has shape $(R, m, m)$ where $m$ is the collision chunk size — the number of samples that can mutually collide within one pass. Memory cost is **quadratic** in $m$: $(R \cdot m^2 \cdot 4)$ bytes for float32.
   - **Retained-graph memory (training):** The $(R,m,m)$ block is the input to `torch.log`, so autograd keeps **one block per chunk alive** until backward. With $n_{\text{chunks}} = \lceil B/m\rceil$ the training peak is
     $$\text{peak} \approx n_{\text{chunks}} \cdot R \cdot m^2 \cdot 4 = B \cdot R \cdot m \cdot 4 \quad(\text{LINEAR in } m).$$
     So the chunk size cannot be sized to fill the whole budget with one block (the old $m_{\max}=\sqrt{\text{mem\_budget}/(R\cdot4\cdot3)}$ did exactly that, then stacked several chunks on top → OOM). Instead $m$ is sized for a fixed number of coexisting blocks:
     $$m_{\max} = \max\!\bigl(512,\; \lfloor\sqrt{\text{mem\_budget} / (\text{TARGET\_CHUNKS}\cdot R \cdot 4 \cdot 3)}\rfloor\bigr),\qquad B = \text{TARGET\_CHUNKS}\cdot m_{\max}.$$
     The factor 3 covers the transient `torch.log` output plus backward buffers. The batch is split into $\lceil B/m_{\max}\rceil$ chunks; cross-chunk log-probabilities are concatenated and reduced via a single logsumexp. No samples are discarded. (Under `no_grad` evaluation each block is freed after use, so eval batches can be far larger; this bound is the *training* limit.)
   - **Honest ceiling & the knob:** The theoretical ceiling is $\log_2(m_{\max})$; chunk *size* (quadratic memory) raises it, chunk *count* only lowers variance. At fixed peak these trade off — `COLLISION_TARGET_CHUNKS` (module constant in `run.py`) is that knob: raising it shrinks $m$ (lower ceiling) for more pairs / a bigger batch $B = \text{COLLISION\_TARGET\_CHUNKS}\cdot m_{\max}$, at constant memory $\approx \text{mem\_budget}/3$. It is the lever to turn when a larger collision batch is needed. Default `= 4`. Example (63 GiB budget, $R=15$): $m_{\max}\approx 9.3\text{k}$, ceiling $\approx 13.2$ bits, $B\approx 37$k, peak $\approx 20$ GiB.

5. **Kolchinsky–Tracey Lower Bound (`'kt'`):** $\mathcal{O}(B^2 \cdot R)$ time, $\mathcal{O}(\text{chunk}^2 \cdot R)$ memory
   - **The Math:** Certified **lower bound** on the true Shannon joint entropy via pairwise Bhattacharyya affinities: $H_{\text{KT}} = H_{\text{cond}} - \tfrac{1}{B}\sum_i \log_2\!\big(\tfrac{1}{B}\sum_j \text{BC}(i,j)\big)$. Two nested chunk loops (outer $i$-chunks of size $m$, inner $j$-chunks spanning the whole batch) build $(m, n)$ blocks of $\log\text{BC}$ by broadcasting $A_i\,(m,1,R)$ against $A_j\,(1,n,R)$ and summing over $R$; the $j$-blocks are `cat`ed into the full $(m, B)$ row and reduced by `logsumexp`.
   - **Retained-graph memory (training):** each $(m,n,R)$ block is fed to `torch.log` and kept for backward. The double loop touches **all** $n_{\text{chunks}}^2$ chunk pairs, so the retained total is
     $$\text{peak} = n_{\text{chunks}}^2 \cdot R \cdot m^2 \cdot 4 = B^2 \cdot R \cdot 4 \quad(\text{INDEPENDENT of } m):$$
     the per-block $m^2$ and the $n_{\text{chunks}}^2$ count cancel, so **chunking saves no training memory** — the only lever is $B$. The training batch is therefore a single chunk at the largest $B$ that fits, $B = \lfloor\sqrt{\text{mem\_budget}/(R\cdot4\cdot3)}\rfloor$ (this differs from collision, which reduces $m$ to fit several chunks). Under `no_grad` eval, blocks are freed per pair and $B$ can be much larger.
   - **Big-batch KT via `recompute_backward` (enabled):** setting the config flag **`recompute_backward=True`** wraps each **$i$-chunk's contribution** (`_kt_row_contribution`: its inner $j$-loop, $(m,n,R)$ blocks and $(m,B)$ row) in **gradient checkpointing** (`torch.utils.checkpoint`), so those tensors are recomputed in backward instead of retained. The retained graph drops from $B^2\!\cdot\!R$ to $\mathcal{O}(B\!\cdot\!R)$ (just the $(B,R)$ $\sqrt{A}$ inputs), so **memory no longer binds** — the unchanged $\mathcal{O}(B^2 R)$ **compute** does. The sizer then caps the *work*, not the RAM: $B = \lfloor\sqrt{\texttt{KT\_COMPUTE\_BUDGET}/R}\rfloor$ (same $\sqrt{1/R}$ shape as the memory cap), still floored by the physics/forward cap. `KT_COMPUTE_BUDGET` defaults to $4\times$ the memory-bound work ($B^2R\approx 5.3\mathrm{e}9$ on an 80 GiB A100) $\Rightarrow \sim\!2\times$ batch, $\sim\!4\times$ step time. Exact (identical value + gradient); costs one extra inner-loop forward in backward. Not applied to the KL **upper** bound (measurement-only, `no_grad` ⇒ checkpoint is a no-op).
   - **Ceiling via diagonal anchoring:** Unlike collision, KT keeps the **diagonal** ($\text{BC}(i,i)=1$) and normalises by the **total** batch $B$. The self-term anchors each inner sum at $1/B$, so the resolvable-entropy ceiling is $\log_2 B$ (total batch), **not** $\log_2(\text{chunk})$ as for collision. Since training runs one chunk ($\text{chunk}=B$), the ceiling is $\log_2 B$ — for the same peak this is *higher* than collision's multi-chunk $\log_2 m$. Collapsing the inner loop to adjacent chunks or masking the diagonal would silently re-cap the ceiling and is wrong here.
   - *Impact:* Time is $\mathcal{O}(B^2 R)$ — **quadratic in the total batch** (unlike collision, whose adjacent-chunk scheme is linear in chunk count). Very large eval batches therefore cost quadratically even when they fit in memory.
   - **Launch overhead & `compile_kt` (enabled):** the tiled double loop is **launch-bound** — each $(m,n,R)$ tile fires several tiny elementwise/reduce kernels (`mul, mul, add, log, sum`), so profiling shows ~60% of CPU in `cudaLaunchKernel`. Setting the config flag **`compile_kt=True`** wraps the per-tile kernel (`_kt_logbc_tile`: `Σ_R log(√(A_iA_j)+√((1-A_i)(1-A_j)))`) in **`torch.compile`** (`dynamic=True`, cached module-level), so TorchInductor fuses it into ~1 kernel and can skip materialising the $(m,n,R)$ intermediate — cutting both launch count and memory traffic. Exact up to fp reordering (value identical, gradient to ~1e-8). One-time compile cost on first call.

6. **Blocked Shannon (`'blocked'`):** $\mathcal{O}(R^2 + B \cdot 2^{\text{block\_size}} \cdot \lceil R/\text{block\_size}\rceil)$
   - **The Math:** Receptors are partitioned into blocks of at most `block_size` (default 15) by correlation-aware greedy clustering. The clustering step computes the absolute Pearson correlation matrix (one `x.T @ x` matmul on centred activity, shape $(R, R)$), then iterates over $R$ receptors in a Python loop — $\mathcal{O}(R^2)$ total. Within each block, exact Shannon entropy is computed on a $(B, 2^{\text{block\_size}})$ histogram. The partition is cached and refreshed every `block_refresh_interval` (default 50) training steps.
   - *Impact:* The $R^2$ correlation step is negligible vs. the per-block histograms for typical $R$. Memory scales as $\mathcal{O}(B \cdot 2^{\text{block\_size}})$ per block — independent of $R$, making it practical for large arrays. For `block_size=15` and $B=5000$, each block histogram is ~640 MB.
   - **Batch cap & `recompute_backward`:** the $(B, 2^{\text{block\_size}})$ histogram is the joint tensor `compute_shannon_joint_entropy` builds per block; it is retained for backward, so the auto batch is capped at $B \le \text{mem\_budget}/(2^{\text{block\_size}}\cdot n_{\text{blk}}\cdot 4\cdot \text{SAFETY})$ with `SAFETY=4` (forward matrix + backward graph). Setting the config flag **`recompute_backward=True`** gradient-checkpoints each block's Shannon term (`torch.utils.checkpoint`), so the histogram is recomputed in backward instead of retained — `SAFETY` drops to $\approx 2$, roughly **doubling** the feasible batch, for one extra forward of the block (~+20% step time). Checkpointing is exact (identical value and gradients). It is **not** applied to `'annealed'`: that loss also holds an un-checkpointed collision $(R,B,B)$ block which would then dominate and OOM at the larger $B$. Conservative estimate — further gains need chunked-over-$B$ histogram accumulation (drops the $n_{\text{blk}}$ coexistence too).

7. **Proxy Loss Function Bottleneck (`DiscreteProxyLoss`):** $\mathcal{O}(R^2)$
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