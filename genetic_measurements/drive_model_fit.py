# %%
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
from scipy.optimize import minimize, curve_fit
from scipy.special import comb as binom_coeff, expit, logsumexp

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
# =============================================================================
# DATA LOADING
# =============================================================================
directory = "/mnt/hcleroy/PostDoc2/octopus_smelling/data/octopus_informationCoding/"
mat_path = directory + "/20260302_HiPlexResults.mat"
md = loadmat(mat_path, squeeze_me=True, struct_as_record=False)
if "binaryTable_allGenes" not in md or "geneList" not in md:
    raise KeyError("Expected 'binaryTable_allGenes' and 'geneList'.")

A = np.array(md["binaryTable_allGenes"]).astype(bool).astype(int)  # (nR, nC)
CRnames = np.array(_to_str_list(md["geneList"]))
nR, nC = A.shape

# Exclude cells with k=0 (no gene expressed)
coexp = A.sum(axis=0)
mask_nz = coexp > 0
A_nz = A[:, mask_nz]
nC_nz = A_nz.shape[1]
coexp_nz = coexp[mask_nz]

print(f"n_genes = {nR},  n_cells = {nC},  cells with k>=1 = {nC_nz}")

# %%
# =============================================================================
# RANK GENES BY EMPIRICAL MARGINAL (descending)
# =============================================================================
p_emp = A_nz.mean(axis=1)          # per-gene marginal P(x_g=1 | k>=1)
id_sort = np.argsort(-p_emp)       # most expressed first
ranks = np.empty(nR, dtype=int)
ranks[id_sort] = np.arange(nR)     # r_g = 0 for most expressed

# %%
# =============================================================================
# GAUSS-HERMITE QUADRATURE SETUP
# =============================================================================
N_QUAD = 40
gh_nodes, gh_weights = np.polynomial.hermite.hermgauss(N_QUAD)
# s_i = sqrt(2) * sigma_s * t_i,  omega_i = w_i / sqrt(pi)
# We factor out sigma_s at call time.
gh_omega = gh_weights / np.sqrt(np.pi)  # (N_QUAD,)

# %%
# =============================================================================
# LOG-LIKELIHOOD (vectorized over cells and quadrature nodes)
# =============================================================================
# A_nz: (nR, nC_nz),  ranks: (nR,)

def _nll(params, A_data, ranks_vec, return_components=False):
    """
    Negative log-likelihood of the drive model, conditional on k>=1.

    Parameters
    ----------
    params : (c, gamma, log_sigma_s)
        log_sigma_s = log(sigma_s) to allow unconstrained optimisation;
        sigma_s = exp(log_sigma_s).
    A_data : (nR, nC_nz) binary array.
    ranks_vec : (nR,) integer ranks.
    return_components : if True, also return per-cell log-likelihoods.

    Returns
    -------
    nll : scalar negative log-likelihood.
    """
    c, gamma, log_sigma = params
    sigma_s = np.exp(log_sigma)

    h = c - gamma * ranks_vec                       # (nR,)
    s_nodes = np.sqrt(2.0) * sigma_s * gh_nodes     # (Q,)

    # log p_g(s_i) and log(1-p_g(s_i)) for every gene × quad node
    # logit_gi = h_g + s_i,  shape (nR, Q)
    logit = h[:, None] + s_nodes[None, :]            # (nR, Q)

    # Numerically stable log-sigmoid: log p = -softplus(-logit),
    #                                  log(1-p) = -softplus(logit)
    log_p = -np.logaddexp(0.0, -logit)               # (nR, Q)
    log_1mp = -np.logaddexp(0.0, logit)               # (nR, Q)

    # --- Per-cell, per-node log-likelihood ---
    # log L_i(c) = sum_g [ x_cg * log_p_g(s_i) + (1-x_cg) * log(1-p_g(s_i)) ]
    # A_data: (nR, nC_nz),  log_p: (nR, Q)
    # Result: (nC_nz, Q)
    log_L = A_data.T @ log_p + (1.0 - A_data).T @ log_1mp   # (nC_nz, Q)

    # --- log P(x_c) = logsumexp over nodes of [log(omega_i) + log L_i(c)] ---
    log_omega = np.log(gh_omega)                              # (Q,)
    log_Px = logsumexp(log_omega[None, :] + log_L, axis=1)    # (nC_nz,)

    # --- Zero-truncation: P(k>=1) = 1 - P(k=0) ---
    # P(k=0 | s_i) = prod_g (1-p_g(s_i)) = exp(sum_g log(1-p_g(s_i)))
    log_P0_per_node = log_1mp.sum(axis=0)                     # (Q,)
    log_P0 = logsumexp(log_omega + log_P0_per_node)            # scalar
    log_P_ge1 = np.log1p(-np.exp(log_P0))                     # log(1 - P(k=0))

    # --- NLL ---
    log_lik = log_Px - log_P_ge1                               # (nC_nz,)
    nll = -log_lik.sum()

    if return_components:
        return nll, log_lik
    return nll


def _nll_no_drive(params2, A_data, ranks_vec):
    """NLL with sigma_s fixed to 0 (independent model, no drive)."""
    c, gamma = params2
    return _nll(np.array([c, gamma, -50.0]), A_data, ranks_vec)

# %%
# =============================================================================
# INITIALISATION from linear fit of log(p_g) vs rank
# =============================================================================
p_sorted = p_emp[id_sort]
mask_pos = p_sorted > 0
log_p_sorted = np.log(p_sorted[mask_pos])
r_sorted = np.arange(nR)[mask_pos].astype(float)

slope, intercept = np.polyfit(r_sorted, log_p_sorted, 1)
gamma_init = -slope           # gamma > 0 since slope < 0
c_init = intercept            # intercept ≈ logit(p_0) roughly

# Fano factor for sigma_s init
fano = coexp_nz.var() / coexp_nz.mean()
sigma_init = max(np.sqrt(max(fano - 1.0, 0.0)), 0.1)
log_sigma_init = np.log(sigma_init)

print(f"\nInitialisation: c={c_init:.3f}, gamma={gamma_init:.4f}, "
      f"sigma_s={sigma_init:.3f}  (Fano={fano:.3f})")

# %%
# =============================================================================
# FIT FULL MODEL (c, gamma, sigma_s)
# =============================================================================
print("\nFitting full drive model (c, gamma, sigma_s) ...")
res_full = minimize(
    _nll, x0=[c_init, gamma_init, log_sigma_init],
    args=(A_nz, ranks),
    method="L-BFGS-B",
    bounds=[(None, None), (0.0, None), (None, None)],
    options={"maxiter": 500, "ftol": 1e-12},
)
c_fit, gamma_fit, log_sigma_fit = res_full.x
sigma_fit = np.exp(log_sigma_fit)
nll_full = res_full.fun

print(f"  c     = {c_fit:.4f}")
print(f"  gamma = {gamma_fit:.4f}")
print(f"  sigma = {sigma_fit:.4f}")
print(f"  NLL   = {nll_full:.2f}  (converged={res_full.success})")

# %%
# =============================================================================
# FIT INDEPENDENT MODEL (sigma_s = 0)
# =============================================================================
print("\nFitting independent model (sigma_s = 0) ...")
res_indep = minimize(
    _nll_no_drive, x0=[c_init, gamma_init],
    args=(A_nz, ranks),
    method="L-BFGS-B",
    bounds=[(None, None), (0.0, None)],
    options={"maxiter": 500, "ftol": 1e-12},
)
c_indep, gamma_indep = res_indep.x
nll_indep = res_indep.fun

print(f"  c     = {c_indep:.4f}")
print(f"  gamma = {gamma_indep:.4f}")
print(f"  NLL   = {nll_indep:.2f}  (converged={res_indep.success})")
print(f"  delta NLL (indep - full) = {nll_indep - nll_full:.2f}")

# %%
# =============================================================================
# EXISTING 1-PARAM FITS (exponential, combinatorial) for comparison
# =============================================================================
max_k = int(coexp_nz.max())
centers = np.arange(1, max_k + 1)
counts = np.bincount(coexp_nz, minlength=max_k + 1)[1:]
P_emp_k = counts / counts.sum()

def exp_distrib(l, beta):
    return beta * np.exp(-beta * l)

def combinatorial_distrib(k, beta):
    unnorm = binom_coeff(nR, k) * np.exp(-beta * k)
    Z_trunc = (1.0 + np.exp(-beta))**nR - 1.0
    return unnorm / Z_trunc

popt_exp, _ = curve_fit(exp_distrib, centers, P_emp_k, p0=1, bounds=(0.1, 1))
popt_comb, _ = curve_fit(combinatorial_distrib, centers, P_emp_k, p0=1.0, bounds=(0, np.inf))

# %%
# =============================================================================
# MODEL PREDICTIONS: p_g, P(k), correlations, Fano
# =============================================================================

def _model_pg(c_val, gamma_val, sigma_val):
    """Predicted marginal p_g = E_s[sigmoid(h_g + s)] via quadrature."""
    h = c_val - gamma_val * np.arange(nR, dtype=float)  # rank-ordered
    s_nodes = np.sqrt(2.0) * sigma_val * gh_nodes
    logit = h[:, None] + s_nodes[None, :]
    p_gs = expit(logit)                                  # (nR, Q)
    return (p_gs * gh_omega[None, :]).sum(axis=1)         # (nR,)


def _model_Pk_fft(c_val, gamma_val, sigma_val):
    """
    Predicted P(k | k>=1) via Poisson-Binomial (FFT characteristic function)
    averaged over quadrature nodes.

    The characteristic function (CF) of a Poisson-Binomial distribution is
    CF(t_m) = prod_g [ (1-p_g) + p_g * exp(2*pi*i*m/(N+1)) ],  m = 0..N.
    With this +i convention CF is the *inverse* DFT of the pmf, so the pmf
    is recovered by the forward transform: pmf = Re(fft(CF)) / (N+1).
    (Using ifft here would return the pmf reversed in k.)
    """
    N = nR
    h = c_val - gamma_val * np.arange(N, dtype=float)
    s_nodes = np.sqrt(2.0) * sigma_val * gh_nodes       # (Q,)

    pmf_accum = np.zeros(N + 1)
    for qi in range(N_QUAD):
        p_g = expit(h + s_nodes[qi])                     # (N,)
        # FFT-based Poisson-Binomial pmf
        t = 2.0 * np.pi * np.arange(N + 1) / (N + 1)    # (N+1,)
        eit = np.exp(1j * t)                              # (N+1,)
        # CF = prod_g [(1-p_g) + p_g * exp(i*t)]
        cf = np.prod((1.0 - p_g)[:, None] + p_g[:, None] * eit[None, :], axis=0)
        pmf_node = np.real(np.fft.fft(cf, n=N + 1)) / (N + 1)
        pmf_node = np.maximum(pmf_node, 0.0)
        pmf_accum += gh_omega[qi] * pmf_node

    # Zero-truncate
    P0 = pmf_accum[0]
    pmf_trunc = pmf_accum[1:] / (1.0 - P0)
    return pmf_trunc  # length N, index i -> P(k=i+1 | k>=1)


def _model_corr_matrix(c_val, gamma_val, sigma_val):
    """
    Predicted pairwise correlation matrix.

    Cov(x_g, x_h) = E_s[p_g(s) p_h(s)] - E_s[p_g(s)] E_s[p_h(s)].
    The drive s induces all pairwise correlations (gene-agnostic coupling).
    """
    N = nR
    h = c_val - gamma_val * np.arange(N, dtype=float)
    s_nodes = np.sqrt(2.0) * sigma_val * gh_nodes

    p_all = expit(h[:, None] + s_nodes[None, :])         # (N, Q)

    # E[p_g] and E[p_g * p_h]
    Ep = (p_all * gh_omega[None, :]).sum(axis=1)          # (N,)
    Epp = (p_all * gh_omega[None, :]) @ p_all.T           # (N, N): sum_i omega_i p_g(s_i) p_h(s_i)

    cov = Epp - np.outer(Ep, Ep)
    var = Ep * (1.0 - Ep)
    std = np.sqrt(np.maximum(var, 1e-30))
    corr = cov / np.outer(std, std)
    np.fill_diagonal(corr, 1.0)
    return corr


def _model_fano(c_val, gamma_val, sigma_val):
    """Predicted Fano factor Var(k)/E[k], conditional on k>=1."""
    pmf = _model_Pk_fft(c_val, gamma_val, sigma_val)
    ks = np.arange(1, nR + 1, dtype=float)
    mu = (pmf * ks).sum()
    var = (pmf * ks**2).sum() - mu**2
    return var / mu

# --- Full model predictions ---
pg_full = _model_pg(c_fit, gamma_fit, sigma_fit)
Pk_full = _model_Pk_fft(c_fit, gamma_fit, sigma_fit)
corr_full = _model_corr_matrix(c_fit, gamma_fit, sigma_fit)
fano_full = _model_fano(c_fit, gamma_fit, sigma_fit)

# --- Independent model predictions ---
pg_indep = _model_pg(c_indep, gamma_indep, 1e-15)
Pk_indep = _model_Pk_fft(c_indep, gamma_indep, 1e-15)
corr_indep = _model_corr_matrix(c_indep, gamma_indep, 1e-15)
fano_indep = _model_fano(c_indep, gamma_indep, 1e-15)

# --- Empirical correlation matrix (in rank order) ---
A_ranked = A_nz[id_sort, :]     # reorder genes by rank
corr_data = np.corrcoef(A_ranked)

# %%
# =============================================================================
# PRINT RESULTS
# =============================================================================
fano_data = coexp_nz.var() / coexp_nz.mean()

# Mean off-diagonal correlation
mask_upper = np.triu(np.ones((nR, nR), dtype=bool), k=1)
mean_corr_data = corr_data[mask_upper].mean()
mean_corr_full = corr_full[mask_upper].mean()
mean_corr_indep = corr_indep[mask_upper].mean()

print("\n" + "=" * 60)
print("  RESULTS")
print("=" * 60)
print(f"  Full model:   c={c_fit:.4f}, gamma={gamma_fit:.4f}, sigma_s={sigma_fit:.4f}")
print(f"  Indep model:  c={c_indep:.4f}, gamma={gamma_indep:.4f}, sigma_s=0")
print(f"  NLL full      = {nll_full:.2f}")
print(f"  NLL indep     = {nll_indep:.2f}")
print(f"  delta NLL     = {nll_indep - nll_full:.2f}")
print()
print(f"  Fano factor (Var(k)/E[k]):")
print(f"    data        = {fano_data:.4f}  {'(>1: over-dispersed)' if fano_data > 1 else '(<1: sub-dispersed)'}")
print(f"    full model  = {fano_full:.4f}")
print(f"    indep model = {fano_indep:.4f}")
print()
print(f"  Mean off-diagonal pairwise correlation:")
print(f"    data        = {mean_corr_data:.4f}")
print(f"    full model  = {mean_corr_full:.4f}")
print(f"    indep model = {mean_corr_indep:.4f}")
print()
# The drive is gene-agnostic -> the model predicts no specific pairing.
# Quantify: per-pair corr(data, model) measures how much who-with-whom the
# model captures (expected ~0); spread comparison shows the model misses the
# pair-to-pair heterogeneity (only its mean is a genuine prediction).
_cd = corr_data[mask_upper]
_cf = corr_full[mask_upper]
r_pair = np.corrcoef(_cd, _cf)[0, 1]
print(f"  Per-pair corr(data, model)   = {r_pair:.4f}  (~0: no who-with-whom info)")
print(f"  Std of pairwise correlation: data={_cd.std():.4f}, model={_cf.std():.4f}")
print("=" * 60)

# %%
# =============================================================================
# FIGURE (2×2)
# =============================================================================
fig, axes = plt.subplots(2, 2, figsize=(11, 9))

# --- (a) p_g vs rank ---
ax = axes[0, 0]
p_data_ranked = p_emp[id_sort]
ax.semilogy(np.arange(nR), p_data_ranked, "ko", ms=5, label="data")
ax.semilogy(np.arange(nR), pg_full, "r-", lw=2, label=f"full ($\\gamma$={gamma_fit:.3f})")
ax.semilogy(np.arange(nR), pg_indep, "b--", lw=1.5, label=f"indep ($\\gamma$={gamma_indep:.3f})")
ax.set_xlabel("Gene rank (most → least expressed)")
ax.set_ylabel("$p_g$  (marginal expression probability)")
ax.set_title("(a)  Per-gene marginal vs rank")
ax.legend(fontsize=8)

# --- (b) P(k | k>=1) ---
ax = axes[0, 1]
ax.semilogy(centers, P_emp_k, "ko", ms=5, label="data")
ax.semilogy(centers, Pk_full[:max_k], "r-", lw=2, label="full model")
ax.semilogy(centers, Pk_indep[:max_k], "b--", lw=1.5, label="indep ($\\sigma_s$=0)")
ax.semilogy(centers, exp_distrib(centers, popt_exp[0]), "g:", lw=1.5,
            label=f"exp fit ($\\beta$={popt_exp[0]:.3f})")
ax.semilogy(centers, combinatorial_distrib(centers, popt_comb[0]), "m-.",
            lw=1.5, label=f"$C(N,k)e^{{-\\beta k}}/Z$")
ax.set_xlabel("$k$  (genes expressed per cell)")
ax.set_ylabel("$P(k \\mid k \\geq 1)$")
ax.set_title("(b)  Count distribution")
ax.set_ylim(1e-4, 1)
ax.legend(fontsize=7)

# --- (c) Distribution of pairwise correlations ---
# The drive is gene-agnostic: predicted corr_gh depends ONLY on the two
# marginals (corr_gh ~ sigma_s^2 * sqrt(v_g v_h),  v=p(1-p)). It carries NO
# "who-correlates-with-whom" information and is invariant to relabeling genes.
# So only the DISTRIBUTION (its mean, and the marginal-driven spread) is a real
# prediction — a per-pair scatter would imply an identity claim the model never
# makes. We compare distributions, not pairs.
ax = axes[1, 0]
corr_data_flat = corr_data[mask_upper]
corr_full_flat = corr_full[mask_upper]
bins = np.linspace(min(corr_data_flat.min(), corr_full_flat.min()),
                   corr_data_flat.max(), 40)
ax.hist(corr_data_flat, bins=bins, alpha=0.5, color="black", label="data")
ax.hist(corr_full_flat, bins=bins, alpha=0.5, color="red", label="full model")
ax.axvline(mean_corr_data, color="black", ls="--", lw=1.5,
           label=f"data mean = {mean_corr_data:.3f}")
ax.axvline(mean_corr_full, color="red", ls="--", lw=1.5,
           label=f"model mean = {mean_corr_full:.3f}")
ax.axvline(0.0, color="blue", ls=":", lw=1.5, label="indep (all 0)")
ax.set_xlabel("Pairwise correlation")
ax.set_ylabel("# gene pairs")
ax.set_title("(c)  Correlation distribution\n(only the mean is a real prediction)")
ax.legend(fontsize=7)

# --- (d) Fano factor: bar chart ---
ax = axes[1, 1]
labels = ["data", "full model", "indep model"]
fanos = [fano_data, fano_full, fano_indep]
colors = ["black", "red", "blue"]
bars = ax.bar(labels, fanos, color=colors, alpha=0.7)
ax.axhline(1.0, color="gray", ls="--", lw=1, label="Poisson limit (Fano=1)")
ax.set_ylabel("Fano factor  $\\mathrm{Var}(k) / \\langle k \\rangle$")
ax.set_title("(d)  Fano factor comparison")
ax.legend(fontsize=8)
for bar, val in zip(bars, fanos):
    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02,
            f"{val:.3f}", ha="center", va="bottom", fontsize=9)

plt.tight_layout()
plt.savefig(directory + "drive_model_fit.png", dpi=150)
plt.show()

# %%
# =============================================================================
# CONCLUSION
# =============================================================================
print("\n" + "=" * 60)
print("  CONCLUSION")
print("=" * 60)
if fano_data > 1.0:
    print(f"  The data ARE over-dispersed: Fano = {fano_data:.3f} > 1.")
    if sigma_fit > 0.01 and (nll_indep - nll_full) > 1.0:
        print(f"  sigma_s = {sigma_fit:.4f} > 0: the drive is active.")
        print(f"  The full model beats sigma_s=0 by delta_NLL = {nll_indep - nll_full:.1f}.")
        print(f"  The drive reproduces the heavy tail of P(k) and generates")
        print(f"  positive pairwise correlations (mean corr: data={mean_corr_data:.4f},")
        print(f"  full={mean_corr_full:.4f}, indep={mean_corr_indep:.4f}).")
    else:
        print(f"  However sigma_s ≈ 0 or NLL improvement is negligible;")
        print(f"  the heterogeneous p_g alone may suffice.")
else:
    print(f"  The data are NOT over-dispersed: Fano = {fano_data:.3f} < 1.")
    print(f"  The heterogeneous p_g (exponential in rank) alone accounts")
    print(f"  for the near-exponential P(k). No drive needed (sigma_s -> 0).")
    print(f"  sigma_s fitted = {sigma_fit:.4f}.")
print("=" * 60)

# %%
