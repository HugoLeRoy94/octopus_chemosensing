# %%
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat


def _to_str_list(mat_field):
    """Convert MATLAB cell/char arrays to a flat Python list of strings."""
    x = np.squeeze(mat_field)
    if x.dtype.kind in ("U", "S"):
        if x.ndim == 1:
            return [str(s) for s in x.tolist()]
        else:
            return ["".join(row).strip() for row in x.tolist()]
    out = []
    for elem in np.ravel(x):
        if isinstance(elem, np.ndarray) and elem.dtype.kind in ("U", "S"):
            s = "".join(np.atleast_1d(elem).tolist()).strip()
        elif isinstance(elem, np.ndarray) and elem.dtype.kind in ("i", "f"):
            s = str(elem.item()) if elem.size == 1 else "".join(map(str, elem.tolist()))
        else:
            s = str(elem)
        out.append(s)
    return out


# %%
directory = "/mnt/hcleroy/PostDoc2/octopus_smelling/data/octopus_informationCoding/"
mat_path = directory + "/20260302_HiPlexResults.mat"
md = loadmat(mat_path, squeeze_me=True, struct_as_record=False)
if "binaryTable_allGenes" not in md or "geneList" not in md:
    raise KeyError("Expected fields 'binaryTable_allGenes' and 'geneList' in the .mat file.")

A = np.array(md["binaryTable_allGenes"]).astype(bool).astype(int)  # (nR, nC)
CRnames = np.array(_to_str_list(md["geneList"]))
nR, nC = A.shape
prob_expr = A.mean(axis=1)  # (nR,): marginal expression probability per gene

# %%
# ---------------------------------------------------------
# Empirical: number of unique patterns vs. cumulative batch size
# ---------------------------------------------------------
M_cells = A.T  # (nC, nR): each row is one cell's binary expression pattern
batch_sizes = np.arange(10, nC + 1, 10)

unique_empirical = np.array([
    np.unique(M_cells[:n, :], axis=0).shape[0]
    for n in batch_sizes
])

# %%
# ---------------------------------------------------------
# Model under gene independence:
#   P(sigma) = prod_i  p_i^sigma_i * (1-p_i)^(1-sigma_i)
#   E[# unique | N cells] = sum_{sigma in {0,1}^nR} [1 - (1 - P(sigma))^N]
#
# Feasible only for nR <= 25 (~33 M patterns, ~256 MB float64 peak).
# ---------------------------------------------------------
MAX_NR = 25

if nR <= MAX_NR:
    n_patterns = 2**nR
    indices = np.arange(n_patterns, dtype=np.int32)  # (n_patterns,)

    # log P(sigma) = sum_b log(1-p_b) + sum_b sigma_b * [log(p_b) - log(1-p_b)]
    log_p = np.log(np.clip(prob_expr, 1e-10, 1.0))
    log_q = np.log(np.clip(1.0 - prob_expr, 1e-10, 1.0))
    base = np.sum(log_q)   # contribution when all genes are off
    diff = log_p - log_q   # (nR,): log-odds per gene

    log_probs = np.full(n_patterns, base, dtype=np.float64)
    for b in range(nR):
        log_probs += ((indices >> b) & 1).astype(np.float64) * diff[b]

    # (1 - P(sigma))^N in log space for numerical stability
    one_minus_p = np.clip(1.0 - np.exp(log_probs), 0.0, 1.0)
    log_one_minus_p = np.log(np.where(one_minus_p > 0, one_minus_p, 1e-300))

    unique_model = np.array([
        np.sum(1.0 - np.exp(N * log_one_minus_p))
        for N in batch_sizes
    ])
else:
    unique_model = None
    print(f"nR={nR} > {MAX_NR}: skipping exhaustive model computation.")

# %%
# ---------------------------------------------------------
# Plot
# ---------------------------------------------------------
fig, ax = plt.subplots(figsize=(10, 5))

ax.plot(batch_sizes, unique_empirical, label='Empirical', color='steelblue', lw=2)
if unique_model is not None:
    ax.plot(batch_sizes, unique_model, label='Model (independence)', color='tomato',
            lw=2, linestyle='--')

ax.set_xlabel('Number of cells')
ax.set_ylabel('Number of unique expression patterns')
ax.set_title('Gene expression pattern diversity vs. sample size')
ax.legend()
plt.tight_layout()
plt.savefig('unique_patterns_vs_batch.png', dpi=150)
plt.show()

# %%
