# %%
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
from sklearn.linear_model import LogisticRegression

def _to_str_list(mat_field):
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
mat_path  = directory + "/20260302_HiPlexResults.mat"
md = loadmat(mat_path, squeeze_me=True, struct_as_record=False)
if "binaryTable_allGenes" not in md or "geneList" not in md:
    raise KeyError("Expected 'binaryTable_allGenes' and 'geneList'.")

A        = np.array(md["binaryTable_allGenes"]).astype(bool).astype(int)  # (nR, nC)
CRnames  = np.array(_to_str_list(md["geneList"]))
nR, nC   = A.shape
p_i      = A.mean(axis=1)
co_expr  = (A @ A.T) / nC

# %%
# Load Step 1-2 artifacts produced by expression_models.py
_artifacts_path = directory + "model_artifacts.npz"
try:
    _saved = np.load(_artifacts_path)
except FileNotFoundError:
    raise FileNotFoundError(
        f"Run expression_models.py first to generate {_artifacts_path}"
    )

pmf_indep     = _saved["pmf_indep"]
pmf_mixture   = _saved["pmf_mixture"]
P_emp         = _saved["P_emp"]
centers       = _saved["centers"]
geom_pmf      = _saved["geom_pmf"]
q_fit         = float(_saved["q_fit"])
var_data      = float(_saved["var_data"])
var_indep     = float(_saved["var_indep"])
pi_mix        = _saved["pi_mix"]
Q_mix         = _saved["Q_mix"]          # (nR, K)
best_K        = int(_saved["best_K"])
reordered_idx = _saved["reordered_idx"].astype(int)

# %%
# =============================================================================
# STEP 3 — Pairwise MaxEnt (Ising)
# =============================================================================

# --- 3a. Pseudolikelihood fitting --------------------------------------------
# For each gene i: logistic regression of σ_i on all σ_{j≠i}.
# intercept → h_i,  coeff on σ_j → J_{ij}.  Then symmetrize.

X        = A.T.astype(float)    # (nC, nR)  — cells × genes
h_ising  = np.zeros(nR)
J_ising  = np.zeros((nR, nR))

print("Fitting Ising model via pseudolikelihood …")
for i in range(nR):
    others    = np.delete(X, i, axis=1)          # (nC, nR-1)
    lr        = LogisticRegression(C=100, fit_intercept=True, max_iter=500, solver="lbfgs")
    lr.fit(others, X[:, i])
    h_ising[i] = lr.intercept_[0]
    J_ising[i, np.delete(np.arange(nR), i)] = lr.coef_[0]

J_ising = (J_ising + J_ising.T) / 2             # enforce symmetry
np.fill_diagonal(J_ising, 0.0)
print("  done.")

# --- 3b. Gibbs sampling ------------------------------------------------------
# n_chains chains updated in parallel (vectorized); gene sweep is sequential.

def _gibbs_sample(h, J, n_burn=3_000, n_collect=20_000, n_chains=8, seed=0):
    """
    Sequential-within-sweep Gibbs sampler, vectorized over n_chains.
    Returns (n_collect * n_chains, nR) binary array.
    """
    rng    = np.random.default_rng(seed)
    nR_    = len(h)
    states = rng.integers(0, 2, (n_chains, nR_)).astype(float)
    collected = []
    for step in range(n_burn + n_collect):
        for i in range(nR_):
            field      = h[i] + states @ J[i]          # (n_chains,)
            prob       = 1.0 / (1.0 + np.exp(-field))
            states[:, i] = (rng.random(n_chains) < prob).astype(float)
        if step >= n_burn:
            collected.append(states.copy())
    return np.vstack(collected)                          # (n_collect*n_chains, nR)

print("Running Gibbs sampling …")
samples_ising = _gibbs_sample(h_ising, J_ising)
n_samp        = samples_ising.shape[0]
print(f"  {n_samp} samples collected.")

# --- 3c. Statistics from Gibbs samples ---------------------------------------
mean_si_ising   = samples_ising.mean(axis=0)
mean_sisj_ising = samples_ising.T @ samples_ising / n_samp

mae_si   = np.abs(mean_si_ising - p_i).mean()
mae_sisj = np.abs(mean_sisj_ising - co_expr).mean()
print(f"  MAE ⟨σ_i⟩:      {mae_si:.4f}   (should be ~0)")
print(f"  MAE ⟨σ_iσ_j⟩:  {mae_sisj:.4f}  (should be ~0 by construction)")

# --- 3d. P(L) from Gibbs samples ---------------------------------------------
L_ising  = samples_ising.sum(axis=1).astype(int)
pmf_ising = np.bincount(L_ising, minlength=nR + 1) / n_samp

# %%
# --- Figure 1: J_ij heatmap vs raw correlation -------------------------------
corr_matrix = np.corrcoef(A)
J_max       = np.abs(J_ising).max()

fig, axes = plt.subplots(1, 2, figsize=(13, 5))

im0 = axes[0].imshow(
    J_ising[reordered_idx, :][:, reordered_idx],
    cmap="RdBu_r", vmin=-J_max, vmax=J_max,
)
axes[0].set_xticks(np.arange(nR))
axes[0].set_xticklabels(CRnames[reordered_idx], rotation=90, fontsize=7)
axes[0].set_yticks(np.arange(nR))
axes[0].set_yticklabels(CRnames[reordered_idx], fontsize=7)
axes[0].set_title("Ising couplings  $J_{ij}$  (pseudolikelihood)")
plt.colorbar(im0, ax=axes[0], label="$J_{ij}$")

im1 = axes[1].imshow(
    corr_matrix[reordered_idx, :][:, reordered_idx],
    cmap="RdBu_r", vmin=-1, vmax=1,
)
axes[1].set_xticks(np.arange(nR))
axes[1].set_xticklabels(CRnames[reordered_idx], rotation=90, fontsize=7)
axes[1].set_yticks(np.arange(nR))
axes[1].set_yticklabels(CRnames[reordered_idx], fontsize=7)
axes[1].set_title("Pearson correlation  $r_{ij}$  (raw)")
plt.colorbar(im1, ax=axes[1], label="$r_{ij}$")

plt.tight_layout()
plt.show()

# %%
# --- Figure 2: Sanity scatter plots (marginals + pairwise) -------------------
mask    = np.triu(np.ones((nR, nR), dtype=bool), k=1)
co_flat = co_expr[mask]
ij_flat = mean_sisj_ising[mask]

fig, axes = plt.subplots(1, 2, figsize=(9, 4))

axes[0].scatter(p_i, mean_si_ising, s=40, alpha=0.8)
lim0 = (0, 1)
axes[0].plot(lim0, lim0, "k--", lw=1)
axes[0].set_xlabel("Data $\\langle\\sigma_i\\rangle$")
axes[0].set_ylabel("Ising (Gibbs) $\\langle\\sigma_i\\rangle$")
axes[0].set_title(f"Marginals   MAE = {mae_si:.4f}")

lim1 = (0, max(co_flat.max(), ij_flat.max()) * 1.05)
axes[1].scatter(co_flat, ij_flat, s=15, alpha=0.5)
axes[1].plot(lim1, lim1, "k--", lw=1)
axes[1].set_xlabel("Data $\\langle\\sigma_i\\sigma_j\\rangle$")
axes[1].set_ylabel("Ising (Gibbs) $\\langle\\sigma_i\\sigma_j\\rangle$")
axes[1].set_title(f"Pairwise coexpressions   MAE = {mae_sisj:.4f}")

plt.tight_layout()
plt.show()

# %%
# --- Figure 3: P(L) — all models on the same axes ---------------------------
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(centers,            P_emp,        "ko-", ms=5, label="data")
ax.plot(np.arange(nR + 1), pmf_indep,    "b--",       label="Poisson-Binomial (indep.)")
ax.plot(centers,            geom_pmf,     "r:",        label=f"Geometric (q={q_fit:.3f})")
ax.plot(np.arange(nR + 1), pmf_mixture,  "g-",        label=f"Bernoulli mixture (K={best_K})")
ax.plot(np.arange(nR + 1), pmf_ising,    "m-.",       label="Ising (MaxEnt)")
ax.set_yscale("log")
ax.set_ylim(1e-4, 1)
ax.set_xlabel("Number of co-expressed genes  $L$")
ax.set_ylabel("$P(L)$")
ax.set_title(f"P(L): data vs all models   (var ratio = {var_data/var_indep:.2f})")
ax.legend()
plt.tight_layout()
plt.show()

# %%
# =============================================================================
# STEP 4 — Triplet test (falsification)
# =============================================================================

def _triplet_means(S):
    """
    ⟨σ_i σ_j σ_k⟩ for all (i,j,k) from binary sample matrix S: (n, nR).
    Uses a loop over i (24 iterations) so intermediate arrays stay at (n, nR).
    T[i,j,k] = mean_n S[n,i]*S[n,j]*S[n,k]
              = (S * S[:,i:i+1]).T @ S / n
    """
    n, nR_ = S.shape
    T = np.zeros((nR_, nR_, nR_))
    for i in range(nR_):
        T[i] = (S * S[:, i:i+1]).T @ S / n
    return T

print("\nComputing triplet statistics …")

# Observed (data)
T_obs   = _triplet_means(A.T.astype(float))

# Ising (Gibbs samples)
T_ising = _triplet_means(samples_ising)

# Independent: p_i * p_j * p_k
T_indep = p_i[:, None, None] * p_i[None, :, None] * p_i[None, None, :]

# Mixture: analytical  Σ_k π_k * q_{i,k} * q_{j,k} * q_{l,k}
T_mix = np.einsum("k,ik,jk,lk->ijl", pi_mix, Q_mix, Q_mix, Q_mix)

# --- Filter triplets with at least min_count cells with all three on --------
min_count = 20
tri_i, tri_j, tri_k = np.where(
    (np.arange(nR)[:, None, None] < np.arange(nR)[None, :, None]) &
    (np.arange(nR)[None, :, None] < np.arange(nR)[None, None, :])
)
all3_count = (T_obs * nC).astype(int)
valid  = all3_count[tri_i, tri_j, tri_k] >= min_count
tri_i, tri_j, tri_k = tri_i[valid], tri_j[valid], tri_k[valid]
n_trip = len(tri_i)
print(f"  Triplets with ≥{min_count} cells all-on: {n_trip} / {nR*(nR-1)*(nR-2)//6}")

obs_vals   = T_obs[tri_i,   tri_j, tri_k]
ising_vals = T_ising[tri_i, tri_j, tri_k]
indep_vals = T_indep[tri_i, tri_j, tri_k]
mix_vals   = T_mix[tri_i,   tri_j, tri_k]

# %%
# --- Figure 4: Triplet scatter — one panel per model -------------------------
fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharey=True)
all_preds = [indep_vals, mix_vals, ising_vals]
lim_max   = max(obs_vals.max(), max(v.max() for v in all_preds)) * 1.1
lim_trip  = (0, lim_max)

for ax, pred, label, color in zip(
    axes,
    all_preds,
    [
        "Independent  ($p_i p_j p_k$)",
        f"Bernoulli mixture  (K={best_K})",
        "Ising  (MaxEnt)",
    ],
    ["steelblue", "forestgreen", "darkorchid"],
):
    r   = np.corrcoef(obs_vals, pred)[0, 1]
    mae = np.abs(obs_vals - pred).mean()
    ax.scatter(obs_vals, pred, s=20, alpha=0.6, color=color)
    ax.plot(lim_trip, lim_trip, "k--", lw=1)
    ax.set_xlim(lim_trip)
    ax.set_ylim(lim_trip)
    ax.set_xlabel("Observed $\\langle\\sigma_i\\sigma_j\\sigma_k\\rangle$")
    ax.set_ylabel("Predicted")
    ax.set_title(f"{label}\n$r$={r:.3f}   MAE={mae:.4f}")

plt.suptitle(
    f"Triplet co-expression test   (n={n_trip} triplets, min count={min_count})",
    y=1.02,
)
plt.tight_layout()
plt.show()

# %%
# =============================================================================
# STEP 5 — Save artifacts + print summary
# =============================================================================
out_path = directory + "model_artifacts.npz"
np.savez(
    out_path,
    # Step 1 — independent baseline
    p_i           = p_i,
    pmf_indep     = pmf_indep,
    P_emp         = P_emp,
    centers       = centers,
    geom_pmf      = geom_pmf,
    q_fit         = np.array(q_fit),
    var_data      = np.array(var_data),
    var_indep     = np.array(var_indep),
    # Step 2 — Bernoulli mixture
    pi_mix        = pi_mix,
    Q_mix         = Q_mix,
    best_K        = np.array(best_K),
    pmf_mixture   = pmf_mixture,
    # Step 3 — Ising
    h_ising       = h_ising,
    J_ising       = J_ising,
    pmf_ising     = pmf_ising,
    # Shared
    reordered_idx = reordered_idx,
)
print(f"All artifacts saved → {out_path}")

# --- Summary block -----------------------------------------------------------
r_indep = np.corrcoef(obs_vals, indep_vals)[0, 1]
r_mix   = np.corrcoef(obs_vals, mix_vals)[0, 1]
r_ising = np.corrcoef(obs_vals, ising_vals)[0, 1]

print("\n" + "=" * 56)
print("  SUMMARY")
print("=" * 56)
print(f"  n_genes = {nR},  n_cells = {nC}")
print(f"  Variance ratio (data/indep):  {var_data/var_indep:.3f}")
print(f"  Best mixture K:               {best_K}")
print(f"  Ising sanity — MAE ⟨σ_i⟩:    {mae_si:.4f}")
print(f"  Ising sanity — MAE ⟨σ_iσ_j⟩: {mae_sisj:.4f}")
print(f"  Triplet r  (independent):     {r_indep:.3f}")
print(f"  Triplet r  (mixture K={best_K}):   {r_mix:.3f}")
print(f"  Triplet r  (Ising):           {r_ising:.3f}")
print("=" * 56)

# %%
