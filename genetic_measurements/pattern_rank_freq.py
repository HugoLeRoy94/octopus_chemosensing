# %%
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
import re


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
suckerIDX = np.array(_to_str_list(md["suckerIDX"]))
nR, nC = A.shape

sucker_idx_int = np.array([
    int(re.search(r'\d+', s).group()) if re.search(r'\d+', s) else -1
    for s in suckerIDX
])
M_cells = A.T  # (nC, nR): each row is one cell's binary pattern

# %%
# ---------------------------------------------------------
# Global ranked frequency distribution
# ---------------------------------------------------------
global_patterns, global_counts = np.unique(M_cells, axis=0, return_counts=True)
global_rank_order = np.argsort(-global_counts)
global_patterns_sorted = global_patterns[global_rank_order]   # (n_unique, nR)
global_freq = global_counts[global_rank_order] / nC           # descending, normalized

# Map each pattern (as bytes) to its global rank index for fast lookup
pattern_to_rank = {row.tobytes(): i for i, row in enumerate(global_patterns_sorted)}
n_global_patterns = len(global_patterns_sorted)

# Per-sucker distributions (own ranking + global-rank projection)
ordered_suckers = np.sort(np.unique(sucker_idx_int[sucker_idx_int > 0]))
num_suckers = len(ordered_suckers)

sucker_freqs = {}        # own-rank ordering
sucker_global_freqs = {} # projected onto global rank axis (normalized)
sucker_obs_counts = {}   # projected onto global rank axis (raw integer counts)
sucker_sizes = {}

for sucker in ordered_suckers:
    mask = sucker_idx_int == sucker
    if mask.sum() < 2:
        continue
    n_s = mask.sum()
    patterns_s, counts_s = np.unique(M_cells[mask], axis=0, return_counts=True)

    sucker_sizes[sucker] = n_s

    # Own-rank ordering
    sucker_freqs[sucker] = np.sort(counts_s)[::-1] / n_s

    # Project onto global rank axis
    obs = np.zeros(n_global_patterns, dtype=int)
    for pat, cnt in zip(patterns_s, counts_s):
        rank = pattern_to_rank.get(pat.tobytes())
        if rank is not None:
            obs[rank] = cnt
    sucker_obs_counts[sucker] = obs
    sucker_global_freqs[sucker] = obs / n_s

# %%
# ---------------------------------------------------------
# Plot
# ---------------------------------------------------------
fig, ax = plt.subplots(figsize=(10, 6))

cmap = plt.get_cmap('plasma')
norm = plt.Normalize(vmin=ordered_suckers.min(), vmax=ordered_suckers.max())

for sucker, freq in sucker_freqs.items():
    ax.plot(
        np.arange(1, len(freq) + 1), freq,
        color=cmap(norm(sucker)), alpha=0.6, lw=1.0,
    )

# Global curve on top
ax.plot(
    np.arange(1, len(global_freq) + 1), global_freq,
    color='black', lw=2.5, label='Global', zorder=10,
)

sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
plt.colorbar(sm, ax=ax, label='Sucker ID')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Pattern rank')
ax.set_ylabel('Frequency (fraction of cells)')
ax.set_title('Ranked frequency of gene expression patterns\n(global vs. per sucker)')
ax.legend()
plt.tight_layout()
plt.savefig('pattern_rank_freq.png', dpi=150)
plt.show()

# %%
# ---------------------------------------------------------
# Same plot but all curves share the global rank ordering
# ---------------------------------------------------------
global_ranks = np.arange(1, n_global_patterns + 1)

fig, ax = plt.subplots(figsize=(10, 6))

for sucker, freq_global in sucker_global_freqs.items():
    nonzero = freq_global > 0
    ax.scatter(
        global_ranks[nonzero], freq_global[nonzero],
        color=cmap(norm(sucker)), alpha=0.6, lw=1.0,
    )

ax.plot(
    global_ranks, global_freq,
    color='black', lw=2.5, label='Global', zorder=10,
)

sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
plt.colorbar(sm, ax=ax, label='Sucker ID')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Global pattern rank')
ax.set_ylabel('Frequency (fraction of cells in sucker)')
ax.set_title('Ranked frequency of gene expression patterns\n(all curves on global rank axis)')
ax.legend()
plt.tight_layout()
plt.savefig('pattern_rank_freq_global_axis.png', dpi=150)
plt.show()

# %%
# ---------------------------------------------------------
# Chi-squared goodness-of-fit: does each sucker draw from the global distribution?
# H0: sucker s is a random subsample of the global pattern pool.
# Patterns with expected count < MIN_EXPECTED are pooled into an "other" bin.
# ---------------------------------------------------------
from scipy.stats import chi2 as chi2_dist

MIN_EXPECTED = 5

def bh_correct(p_values):
    n = len(p_values)
    order = np.argsort(p_values)
    p_sorted = p_values[order]
    p_adj_sorted = np.minimum(1.0, p_sorted * n / np.arange(1, n + 1))
    p_adj_sorted = np.minimum.accumulate(p_adj_sorted[::-1])[::-1]
    p_adj = np.empty(n)
    p_adj[order] = p_adj_sorted
    return p_adj

results = []
for sucker in ordered_suckers:
    if sucker not in sucker_obs_counts:
        continue
    n_s = sucker_sizes[sucker]
    observed = sucker_obs_counts[sucker].astype(float)   # (n_global_patterns,)
    expected = global_freq * n_s                          # (n_global_patterns,)

    # Pool rare bins (E_k < MIN_EXPECTED) into a single "other" bin
    keep = expected >= MIN_EXPECTED
    obs_arr = np.append(observed[keep], observed[~keep].sum())
    exp_arr = np.append(expected[keep], expected[~keep].sum())

    n_bins = len(obs_arr)
    dof = n_bins - 1
    if dof < 1:
        continue

    chi2_stat = float(np.sum((obs_arr - exp_arr) ** 2 / exp_arr))
    p_val = float(chi2_dist.sf(chi2_stat, dof))

    results.append({
        'sucker': int(sucker),
        'n_cells': n_s,
        'n_bins': n_bins,
        'dof': dof,
        'chi2': chi2_stat,
        'p_value': p_val,
    })

p_raw = np.array([r['p_value'] for r in results])
p_adj = bh_correct(p_raw)
alpha = 0.05
for r, pa in zip(results, p_adj):
    r['p_adj'] = pa
    r['significant'] = pa < alpha

print(f"{'Sucker':>8} {'N cells':>8} {'bins':>6} {'dof':>5} {'chi2':>10} {'p-value':>12} {'p-adj(BH)':>12} {'sig':>4}")
print("-" * 72)
for r in results:
    print(
        f"{r['sucker']:>8} {r['n_cells']:>8} {r['n_bins']:>6} {r['dof']:>5} "
        f"{r['chi2']:>10.2f} {r['p_value']:>12.3e} {r['p_adj']:>12.3e} "
        f"{'*' if r['significant'] else '':>4}"
    )

# %%
# ---------------------------------------------------------
# Plot: -log10(p_adj) per sucker
# ---------------------------------------------------------
suckers_plot = np.array([r['sucker'] for r in results])
neg_log_p = -np.log10(np.clip(p_adj, 1e-300, 1.0))
sig_flags = np.array([r['significant'] for r in results])

fig, ax = plt.subplots(figsize=(10, 4))

colors = ['tomato' if s else 'steelblue' for s in sig_flags]
ax.bar(np.arange(len(suckers_plot)), neg_log_p, color=colors, edgecolor='black', lw=0.5)
ax.axhline(-np.log10(alpha), color='black', linestyle='--', lw=1.5,
           label=f'BH threshold (α={alpha})')

ax.set_xticks(np.arange(len(suckers_plot)))
ax.set_xticklabels(suckers_plot, fontsize=9)
ax.set_xlabel('Sucker ID')
ax.set_ylabel('$-\log_{10}(p_\mathrm{adj})$')
ax.set_title('Chi-squared test: does each sucker draw from the global pattern distribution?\n(red = significant after BH correction)')
ax.legend()
plt.tight_layout()
plt.savefig('pattern_chisq_pvalues.png', dpi=150)
plt.show()

# %%
