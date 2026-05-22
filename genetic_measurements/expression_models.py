# %%
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat

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
mat_path = directory + "/20260302_HiPlexResults.mat"
md = loadmat(mat_path, squeeze_me=True, struct_as_record=False)
if "binaryTable_allGenes" not in md or "geneList" not in md:
    raise KeyError("Expected 'binaryTable_allGenes' and 'geneList'.")

A = np.array(md["binaryTable_allGenes"]).astype(bool).astype(int)  # (n_genes, n_cells)
CRnames = np.array(_to_str_list(md["geneList"]))
nR, nC = A.shape

# %%
# =============================================================================
# STEP 1 — Independent baseline (Poisson-Binomial)
# =============================================================================

p_i = A.mean(axis=1)                  # per-gene marginal probabilities
coexp = A.sum(axis=0)                  # number of genes expressed per cell
coexp_nz = coexp[coexp > 0]

# --- Poisson-Binomial pmf via generating-function convolution ----------------
# G(z) = prod_i [(1-p_i) + p_i * z]
# Iteratively convolve: pmf[l] = P(L = l)
pmf_indep = np.array([1.0])
for p in p_i:
    pmf_indep = np.convolve(pmf_indep, [1.0 - p, p])
# pmf_indep has length nR+1; index l gives P(L=l)

# --- Variance diagnostics ----------------------------------------------------
mu_data  = coexp.mean()
var_data = coexp.var()
var_indep = np.sum(p_i * (1.0 - p_i))   # Var under independence

print("=" * 50)
print(f"  n_genes = {nR},  n_cells = {nC}")
print(f"  mean coexp  (data)   = {mu_data:.3f}")
print(f"  var  coexp  (data)   = {var_data:.3f}")
print(f"  var  coexp  (indep)  = {var_indep:.3f}")
print(f"  variance ratio       = {var_data / var_indep:.3f}  (expect >> 1)")
print("=" * 50)

# --- Geometric fit via log-linear regression ---------------------------------
# P(L=l) = (1-q) * q^(l-1)  =>  log P = log(1-q) + (l-1)*log(q)
# equivalently  log P = [log(1-q) - log(q)] + l * log(q)
# => slope = log(q),  q = exp(slope)
max_k = int(coexp_nz.max())
centers = np.arange(1, max_k + 1)
counts  = np.bincount(coexp_nz, minlength=max_k + 1)[1:]
P_emp   = counts / counts.sum()

# Only fit on bins with nonzero empirical counts to avoid log(0)
mask = P_emp > 0
log_P = np.log(P_emp[mask])
l_vals = centers[mask].astype(float)

# Linear regression: y = a + b*l
b, a = np.polyfit(l_vals, log_P, 1)   # slope b = log(q), intercept a
q_fit = np.exp(b)
geom_pmf = (1.0 - q_fit) * q_fit ** (centers - 1)

print(f"\n  Geometric fit: q = {q_fit:.4f}  (mean = 1/(1-q) = {1/(1-q_fit):.2f})")
print(f"  Poisson-Binomial mean = {(pmf_indep * np.arange(nR + 1)).sum():.3f}")

# %%
# --- Figure: P(L) comparison on log-y scale ----------------------------------
fig, ax = plt.subplots(figsize=(6, 4))

ax.plot(centers, P_emp,                         "ko-", ms=5, label="data")
ax.plot(np.arange(nR + 1), pmf_indep,           "b--",       label="Poisson-Binomial (indep.)")
ax.plot(centers, geom_pmf,                      "r:",        label=f"Geometric fit (q={q_fit:.3f})")

ax.set_yscale("log")
ax.set_xlabel("Number of co-expressed genes  $L$")
ax.set_ylabel("$P(L)$")
ax.set_title(
    f"P(L): data vs independent model\n"
    f"Var ratio = {var_data/var_indep:.2f}  "
    f"(data var = {var_data:.1f},  indep var = {var_indep:.1f})"
)
ax.legend()
plt.tight_layout()
plt.show()

# %%
# =============================================================================
# STEP 2 — Bernoulli Mixture (latent cell types)
# =============================================================================
from scipy.special import logsumexp
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as ssd

# --- Reordering index (from co-expression hierarchical clustering) -----------
co_expr    = (A @ A.T) / nC
dist_matrix = co_expr.max() - co_expr
np.fill_diagonal(dist_matrix, 0)
dist_matrix = (dist_matrix + dist_matrix.T) / 2          # enforce exact symmetry
condensed_dist = ssd.squareform(dist_matrix)
reordered_idx  = sch.leaves_list(sch.linkage(condensed_dist, method="average"))

# --- EM helpers --------------------------------------------------------------
def _em_bernoulli(X, K, rng, max_iter=200, tol=1e-5):
    """Single EM run. X: (n_cells, n_genes). Returns (pi, Q, ll)."""
    p_mean = X.mean(axis=0)
    Q      = np.clip(p_mean[:, None] + rng.normal(0, 0.1, (X.shape[1], K)), 1e-3, 1 - 1e-3)
    pi     = np.ones(K) / K
    prev_ll = -np.inf
    for _ in range(max_iter):
        log_lik  = X @ np.log(Q) + (1 - X) @ np.log(1 - Q)    # (n_cells, K)
        log_resp = log_lik + np.log(pi)
        log_norm = logsumexp(log_resp, axis=1, keepdims=True)
        log_resp -= log_norm
        ll   = log_norm.sum()
        resp = np.exp(log_resp)
        N_k  = resp.sum(axis=0)
        pi   = N_k / X.shape[0]
        Q    = np.clip(X.T @ resp / N_k, 1e-6, 1 - 1e-6)
        if abs(ll - prev_ll) / (abs(prev_ll) + 1e-10) < tol:
            break
        prev_ll = ll
    return pi, Q, ll


def _fit_mixture(X, K, n_restarts, seed):
    """Best-of-n_restarts EM. Returns (pi, Q, ll)."""
    rng  = np.random.default_rng(seed)
    best = (-np.inf, None, None)
    for _ in range(n_restarts):
        pi, Q, ll = _em_bernoulli(X, K, rng)
        if ll > best[0]:
            best = (ll, pi.copy(), Q.copy())
    return best[1], best[2], best[0]


def _eval_mixture_ll(X, pi, Q):
    Q_c     = np.clip(Q, 1e-10, 1 - 1e-10)
    log_lik = X @ np.log(Q_c) + (1 - X) @ np.log(1 - Q_c)
    return logsumexp(log_lik + np.log(pi), axis=1).sum()


def _pb_pmf(p):
    """Poisson-Binomial pmf via sequential convolution (exact)."""
    pmf = np.array([1.0])
    for pi in p:
        pmf = np.convolve(pmf, [1.0 - pi, pi])
    return pmf


def _mixture_pl_pmf(pi, Q):
    """P(L) for a Bernoulli mixture: π_k-weighted sum of per-component PB pmfs."""
    pmf = np.zeros(Q.shape[0] + 1)
    for k in range(len(pi)):
        pmf += pi[k] * _pb_pmf(Q[:, k])
    return pmf


# --- Fit K = 1 … 8 with 5-fold CV ------------------------------------------
K_values   = np.arange(1, 9)
X          = A.T.astype(float)                                    # (n_cells, n_genes)
fold_edges = np.array_split(np.random.default_rng(42).permutation(nC), 5)

bic_scores   = np.zeros(len(K_values))
cv_ll_scores = np.zeros(len(K_values))
fitted_params = {}

print("\nFitting Bernoulli mixtures …")
for ki, K in enumerate(K_values):
    n_params = K * nR + (K - 1)
    pi, Q, ll = _fit_mixture(X, K, n_restarts=5, seed=ki)
    bic_scores[ki]   = -2 * ll + n_params * np.log(nC)
    fitted_params[K] = (pi, Q, ll)

    fold_ll = 0.0
    for f in range(5):
        val   = fold_edges[f]
        train = np.concatenate([fold_edges[j] for j in range(5) if j != f])
        pi_f, Q_f, _ = _fit_mixture(X[train], K, n_restarts=3, seed=ki * 10 + f)
        fold_ll += _eval_mixture_ll(X[val], pi_f, Q_f)
    cv_ll_scores[ki] = fold_ll / nC

    print(f"  K={K}:  ll={ll:.1f}  BIC={bic_scores[ki]:.1f}  CV LL/cell={cv_ll_scores[ki]:.4f}")

best_K = int(K_values[np.argmax(cv_ll_scores)])
pi_best, Q_best, _ = fitted_params[best_K]
print(f"\n  Best K = {best_K}  (by 5-fold CV LL per cell)")

# %%
# --- Figure A: Model selection — BIC and CV LL vs K -------------------------
fig, axs = plt.subplots(1, 2, figsize=(9, 4))

axs[0].plot(K_values, bic_scores, "o-")
axs[0].axvline(best_K, color="red", linestyle="--", label=f"best K={best_K}")
axs[0].set_xlabel("K")
axs[0].set_ylabel("BIC")
axs[0].set_title("BIC vs K  (lower = better)")
axs[0].legend()

axs[1].plot(K_values, cv_ll_scores, "o-")
axs[1].axvline(best_K, color="red", linestyle="--", label=f"best K={best_K}")
axs[1].set_xlabel("K")
axs[1].set_ylabel("Held-out LL / cell")
axs[1].set_title("5-fold CV log-likelihood vs K")
axs[1].legend()

plt.tight_layout()
plt.show()

# %%
# --- Figure B: Heatmap of q_{i,k} for best K ---------------------------------
fig, ax = plt.subplots(figsize=(max(3, best_K * 1.5), 7))
im = ax.imshow(Q_best[reordered_idx, :], aspect="auto", vmin=0, vmax=1, cmap="viridis")
ax.set_yticks(np.arange(nR))
ax.set_yticklabels(CRnames[reordered_idx], fontsize=7)
ax.set_xticks(np.arange(best_K))
ax.set_xticklabels([f"C{k+1}\n$\\pi$={pi_best[k]:.2f}" for k in range(best_K)], fontsize=9)
ax.set_title(f"Bernoulli mixture: $q_{{i,k}}$  (K={best_K})")
plt.colorbar(im, label="Expression probability $q_{i,k}$")
plt.tight_layout()
plt.show()

# %%
# --- Figure C: P(L) comparison — data + Poisson-Binomial + mixture ----------
pmf_mixture = _mixture_pl_pmf(pi_best, Q_best)

fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(centers,             P_emp,        "ko-", ms=5, label="data")
ax.plot(np.arange(nR + 1),  pmf_indep,    "b--",       label="Poisson-Binomial (indep.)")
ax.plot(centers,             geom_pmf,     "r:",        label=f"Geometric fit (q={q_fit:.3f})")
ax.plot(np.arange(nR + 1),  pmf_mixture,  "g-",        label=f"Bernoulli mixture (K={best_K})")
ax.set_yscale("log")
ax.set_xlabel("Number of co-expressed genes  $L$")
ax.set_ylabel("$P(L)$")
ax.set_title(
    f"P(L): data vs models\n"
    f"Var ratio = {var_data/var_indep:.2f}  "
    f"(data var = {var_data:.1f},  indep var = {var_indep:.1f})"
)
ax.legend()
plt.tight_layout()
ax.set_ylim(10**-4,1)
plt.show()

# %%
# --- Save all Step 1-2 artifacts for expression_models_ising.py -------------
_save_path = directory + "model_artifacts.npz"
np.savez(
    _save_path,
    p_i           = p_i,
    pmf_indep     = pmf_indep,
    P_emp         = P_emp,
    centers       = centers,
    geom_pmf      = geom_pmf,
    q_fit         = np.array(q_fit),
    var_data      = np.array(var_data),
    var_indep     = np.array(var_indep),
    pi_mix        = pi_best,
    Q_mix         = Q_best,
    best_K        = np.array(best_K),
    pmf_mixture   = pmf_mixture,
    reordered_idx = reordered_idx,
)
print(f"Artifacts saved → {_save_path}")

# %%
