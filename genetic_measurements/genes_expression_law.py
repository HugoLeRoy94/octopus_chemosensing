# %%
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
from scipy.optimize import curve_fit


def _to_str_list(mat_field):
    x = np.squeeze(mat_field)
    if x.dtype.kind in ("U", "S"):
        if x.ndim == 1:
            return [str(s) for s in x.tolist()]
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
# =============================================================================
# DATA LOADING
# binaryTable_allGenes: matrix X of shape (n_genes, n_cells), X[g, c] = 1 if
# gene g is detected in cell c, 0 otherwise.
# =============================================================================
directory = "/mnt/hcleroy/PostDoc2/octopus_smelling/data/octopus_informationCoding/"
mat_path = directory + "/20260302_HiPlexResults.mat"
md = loadmat(mat_path, squeeze_me=True, struct_as_record=False)

X = np.array(md["binaryTable_allGenes"]).astype(bool).astype(int)
gene_names = np.array(_to_str_list(md["geneList"]))
n_genes, n_cells = X.shape

# k_c = number of genes expressed in cell c. Cells with k_c = 0 are discarded.
k_all = X.sum(axis=0)
mask_nz = k_all > 0
X_nz = X[:, mask_nz]
k = k_all[mask_nz]
n_cells_nz = X_nz.shape[1]

print(f"n_genes = {n_genes}, n_cells = {n_cells}, cells with k>=1 = {n_cells_nz}")

# %%
# =============================================================================
# OBSERVABLE 1: individual gene expression frequency
# p_g = fraction of cells (with k>=1) in which gene g is expressed.
# Genes are sorted by decreasing p_g; r = 0, 1, ... is the rank.
# =============================================================================
p_g = X_nz.mean(axis=1)
id_sort = np.argsort(-p_g)
p_sorted = p_g[id_sort]
rank = np.arange(n_genes, dtype=float)

# %%
# =============================================================================
# OBSERVABLE 2: probability of k genes expressed per cell
# P(k) estimated by the normalized histogram over k = 1 ... k_max.
# =============================================================================
k_max = int(k.max())
k_values = np.arange(1, k_max + 1)
counts = np.bincount(k, minlength=k_max + 1)[1:]
P_k = counts / counts.sum()

# %%
# =============================================================================
# EXPONENTIAL FITS
# p(r)  = a * exp(-gamma * r)      gamma = decay rate per rank unit
# P(k)  = b * exp(-beta * k)       beta  = decay rate per expressed gene
# Both fitted in log space (least squares on log p), which weights the small
# values as much as the large ones.
# =============================================================================
def exp_law(x, log_amp, rate):
    return log_amp - rate * x


mask_p = p_sorted > 0
popt_p, _ = curve_fit(exp_law, rank[mask_p], np.log(p_sorted[mask_p]), p0=[0.0, 0.1])
log_a, gamma = popt_p

mask_k = P_k > 0
popt_k, _ = curve_fit(exp_law, k_values[mask_k], np.log(P_k[mask_k]), p0=[0.0, 0.5])
log_b, beta = popt_k

print(f"gene frequency : p(r) = {np.exp(log_a):.4f} * exp(-{gamma:.4f} r)")
print(f"count law      : P(k) = {np.exp(log_b):.4f} * exp(-{beta:.4f} k)")

# %%
# =============================================================================
# AVERAGE NUMBER OF GENES EXPRESSED PER CELL
# <k>_data : empirical mean over cells with k>=1.
# <k>_sum  : sum_g p_g, identical to <k>_data (each cell counted once per gene).
# <k>_fit  : mean of the fitted exponential law, normalized on k = 1 ... n_genes.
# =============================================================================
k_mean_data = k.mean()
k_mean_sum = p_g.sum()

k_grid = np.arange(1, n_genes + 1)
w = np.exp(-beta * k_grid)
w /= w.sum()
k_mean_fit = (k_grid * w).sum()

print(f"<k> data        = {k_mean_data:.3f}")
print(f"<k> sum_g p_g   = {k_mean_sum:.3f}")
print(f"<k> exp fit     = {k_mean_fit:.3f}")
print(f"Var(k)/<k>      = {k.var() / k.mean():.3f}")

# %%
# =============================================================================
# FIGURE
# =============================================================================
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

ax = axes[0]
ax.semilogy(rank, p_sorted, "ko", ms=5, label="data")
ax.semilogy(rank, np.exp(exp_law(rank, log_a, gamma)), "r-", lw=2,
            label=f"$p(r)$ = {np.exp(log_a):.3f} exp(-{gamma:.3f} r)")
ax.set_xlabel("gene rank $r$ (most to least expressed)")
ax.set_ylabel("$p_g$  (expression frequency)")
ax.set_title("individual gene expression frequency")
ax.legend(fontsize=8)

ax = axes[1]
ax.semilogy(k_values, P_k, "ko", ms=5, label="data")
ax.semilogy(k_values, np.exp(exp_law(k_values, log_b, beta)), "r-", lw=2,
            label=f"$P(k)$ = {np.exp(log_b):.3f} exp(-{beta:.3f} k)")
ax.axvline(k_mean_data, color="gray", ls="--", lw=1,
           label=f"$\\langle k \\rangle$ = {k_mean_data:.2f}")
ax.set_xlabel("$k$  (genes expressed per cell)")
ax.set_ylabel("$P(k \\mid k \\geq 1)$")
ax.set_title("genes expressed per cell")
ax.legend(fontsize=8)

plt.tight_layout()
plt.show()

# %%
