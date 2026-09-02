# %%
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
import pandas as pd
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as ssd
from scipy.special import comb
from scipy.optimize import linear_sum_assignment
from matplotlib.colors import CenteredNorm, TwoSlopeNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable

def _to_str_list(mat_field):
    """
    Convert MATLAB cell/char arrays to a flat Python list of strings.
    Handles cases like: cell array of strings, char matrix, or object arrays.
    """
    x = np.squeeze(mat_field)
    if x.dtype.kind in ("U", "S"):  # already a (numpy) string array
        # char matrix (n x m) -> single string, else array of strings
        if x.ndim == 1:
            return [str(s) for s in x.tolist()]
        else:
            return ["".join(row).strip() for row in x.tolist()]
    # likely an object array (cell array)
    out = []
    for elem in np.ravel(x):
        s = elem
        # char arrays come as ndim>=2 with dtype 'U'/'S'
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
# --- Load .mat ---
md = loadmat(mat_path, squeeze_me=True, struct_as_record=False)
if "binaryTable_allGenes" not in md or "geneList" not in md:
    raise KeyError("Expected fields 'binaryTable_Stack' and 'CRnames' in the .mat file.")

A = np.array(md["binaryTable_allGenes"]).astype(bool).astype(int)  # shape: (n_receptors, n_cells)
CRnames = np.array(_to_str_list(md["geneList"]))
suckerIDX = np.array(_to_str_list(md["suckerIDX"]))
nR, nC = A.shape

plt.rcParams.update({
    'font.size': 6,
    'font.family': 'serif',
    'mathtext.fontset': 'cm',
    'svg.fonttype': 'path',   # draw glyphs as paths: mathtext renders identically everywhere
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.major.size': 2,
    'ytick.major.size': 2,
    'xtick.major.pad': 1.5,
    'ytick.major.pad': 1.5,
    'axes.labelpad': .5,
    'axes.linewidth': 0.5,
})
TICKFS = 5   # tick-label size, used everywhere

def mesh(ax, M, cmap, norm, x=None, y=None, aspect='equal', bad='white'):
    """
    imshow replacement that stays vectorial in the svg.

    imshow embeds an AxesImage, which every vector backend writes out as an
    embedded raster bitmap. pcolormesh emits one filled path per cell instead.
    Row 0 is put on top and the cells are centred on integer coordinates, so the
    result is indistinguishable from imshow with the default (origin='upper').
    """
    ny, nx = M.shape
    x = np.arange(nx + 1) - 0.5 if x is None else np.asarray(x)
    y = np.arange(ny + 1) - 0.5 if y is None else np.asarray(y)
    cm = plt.get_cmap(cmap).copy()
    cm.set_bad(color=bad)
    ax.set_facecolor(bad)   # masked/NaN cells are simply not drawn by pcolormesh
    qm = ax.pcolormesh(x, y, np.ma.masked_invalid(M), cmap=cm, norm=norm,
                       shading='flat', linewidth=0, edgecolors='face', antialiased=False)
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(y[-1], y[0])
    ax.set_aspect(aspect)
    return qm

# %%
# ---------------------------------------------------------
# 1. Single-gene expression + number of genes per cell, both as % of cells
# ---------------------------------------------------------
prob_expr = A.mean(axis=1)                 # P(sigma_i = 1)
id_sort = np.argsort(-prob_expr)

total_expressed = A.sum(axis=0)            # (n_cells,): number of genes on per cell
k_values = np.arange(1, 6, 1)
k_mask = (total_expressed[np.newaxis, :] == k_values[:, np.newaxis])
k_counts = k_mask.sum(axis=1)

FIG1_W = 3.2   # wide enough for the 22 gene labels not to overlap
fig = plt.figure(figsize=(FIG1_W, 2.2))
# right column left empty so the axes width matches the colorbar-bearing figures
gs = fig.add_gridspec(2, 2, height_ratios=[2, 1], width_ratios=[20, 1],
                      hspace=0.7, wspace=0.05,
                      left=0.22, right=0.92, top=0.95, bottom=0.14)
ax0 = fig.add_subplot(gs[0, 0])
ax1 = fig.add_subplot(gs[1, 0])
fig.add_subplot(gs[0, 1]).axis('off')
fig.add_subplot(gs[1, 1]).axis('off')

ax0.bar(np.arange(nR), 100 * prob_expr[id_sort], color='steelblue',
        edgecolor='black', linewidth=0.5)
ax0.set_xlim(-0.5, nR - 0.5)
ax0.set_ylim(0, 100 * prob_expr.max() * 1.1)
ax0.set_xticks(np.arange(nR))
ax0.set_xticklabels(CRnames[id_sort], rotation=45, fontsize=TICKFS,
                    ha='right', rotation_mode='anchor')
ax0.set_ylabel('% cells expressing')

pct_k = 100 * k_counts / nC
ax1.bar(k_values, pct_k, color='steelblue', edgecolor='black', linewidth=0.5)
ax1.set_ylim(0, pct_k.max() * 1.15)
ax1.set_xticks(k_values)
ax1.set_xlabel('Expressed genes per cell, k')
ax1.set_ylabel('% cells')

plt.savefig('single_gene_expression.svg', bbox_inches='tight', pad_inches=0.01)
plt.show()

# %%
# ---------------------------------------------------------
# 2. Gene expression probability conditioned on total expressed genes
# ---------------------------------------------------------
# cond_prob[i, k] = P(gene i expressed | total expressed == k)
safe_counts = np.where(k_counts > 0, k_counts, np.nan)
cond_prob = (A @ k_mask.T) / safe_counts[np.newaxis, :]  # (nR, len(k_values))

# enrichment[i, k] = P(i | L=k) / (k/N)  — ratio to flat expectation
flat_expectation = np.where(k_values > 0, k_values / nR, np.nan)
enrichment = cond_prob / flat_expectation[np.newaxis, :]

norm_enrich = CenteredNorm(vcenter=1, halfrange=2)
cmap_enrich = plt.get_cmap('RdBu_r').copy()
cmap_enrich.set_bad(color='lightgray')

fig = plt.figure(figsize=(1.5, 1.9))
# no layout engine: margins and hspace set by hand so hspace can go negative
gs = fig.add_gridspec(2, 2, height_ratios=[1, 3], width_ratios=[20, 1],
                      hspace=0.05, wspace=0.05,
                      left=0.22, right=0.92, top=1.0, bottom=0.10)
ax0 = fig.add_subplot(gs[0, 0])
ax1 = fig.add_subplot(gs[1, 0], sharex=ax0)
cax = fig.add_subplot(gs[1, 1])
axs = [ax0, ax1]

axs[0].bar(k_values, k_counts, color='steelblue', edgecolor='black', linewidth=0.5)
axs[0].set_ylabel('Cell count')
# shared x axis: set_xticks on ax0 is overwritten by ax1, hide the labels/ticks instead
axs[0].tick_params(axis='x', labelbottom=False, bottom=False, top=False)

im = mesh(axs[1], enrichment[id_sort, :], cmap_enrich, norm_enrich,
          x=np.arange(k_values[0], k_values[-1] + 2) - 0.5, aspect='auto', bad='lightgray')
axs[1].set_yticks(np.arange(nR))
axs[1].set_yticklabels(CRnames[id_sort], fontsize=TICKFS)
axs[1].set_xticks(k_values)
axs[1].set_xlim(k_values[0] - 0.5, k_values[-1] + 0.5)
axs[1].set_xlabel('Expressed genes per cell (k)')
axs[1].set_ylabel('Gene')

cb = plt.colorbar(im, cax=cax, label='Enrichment')
cb.outline.set_linewidth(0.5)
cb.solids.set_rasterized(False)
plt.savefig('P_conditionned_k5.svg', bbox_inches='tight', pad_inches=0.01)
plt.show()

# %%
# ---------------------------------------------------------
# 3. Matrices, distances and clusterings
# ---------------------------------------------------------
N_CLUST = 4          # number of flat clusters -> block boundaries drawn on the matrices
LINK_METHOD = 'average'

joint = (A @ A.T) / nC                    # P(sigma_i=1, sigma_j=1); diagonal = P(sigma_i=1)
cond = joint / prob_expr[np.newaxis, :]   # P(sigma_i=1 | sigma_j=1); diagonal = 1, NOT symmetric
dice = 2 * joint / (prob_expr[:, None] + prob_expr[None, :])  # symmetrised cond., diagonal = 1
corr_matrix = np.corrcoef(A)
corr_matrix = (corr_matrix + corr_matrix.T) / 2  # kill floating-point asymmetry

def cluster(dist_matrix, n_clust=None):
    """Hierarchical clustering of a symmetric distance matrix, with optimal leaf ordering."""
    n_clust = N_CLUST if n_clust is None else n_clust
    d = dist_matrix.copy()
    np.fill_diagonal(d, 0)
    condensed = ssd.squareform((d + d.T) / 2, checks=False)
    Z = sch.linkage(condensed, method=LINK_METHOD)
    Z = sch.optimal_leaf_ordering(Z, condensed)   # best leaf order allowed by the tree
    order = sch.leaves_list(Z)
    labels = sch.fcluster(Z, n_clust, criterion='maxclust')
    coph = sch.cophenet(Z, condensed)[0]
    return Z, order, labels, coph

# co-expression, joint reading: high joint probability -> small distance
Z_joint, order_joint, lab_joint, coph_joint = cluster(joint.max() - joint)
# co-expression, conditional reading: Dice similarity in [0,1], distance = 1 - similarity
Z_dice, order_dice, lab_dice, coph_dice = cluster(1.0 - dice)
# correlation: high correlation -> small distance
Z_corr, order_corr, lab_corr, coph_corr = cluster(1.0 - corr_matrix)

print(f"cophenetic corr   joint: {coph_joint:.3f}   dice: {coph_dice:.3f}   corr: {coph_corr:.3f}")

def mask_diag(M):
    """Copy of M with the diagonal set to NaN (drawn with the colormap 'bad' colour)."""
    Mm = M.astype(float).copy()
    np.fill_diagonal(Mm, np.nan)
    return Mm

def draw_blocks(ax, labels_ordered):
    """Outline the diagonal blocks of consecutive equal cluster labels."""
    cuts = np.where(np.diff(labels_ordered) != 0)[0] + 0.5
    edges = np.concatenate(([-0.5], cuts, [labels_ordered.size - 0.5]))
    for s, e in zip(edges[:-1], edges[1:]):
        ax.hlines([s, e], s, e, color='black', lw=0.7)
        ax.vlines([s, e], s, e, color='black', lw=0.7)

def show_matrix(ax, cax, M, order, labels, title, cbar_label,
                cmap='viridis', norm=None, bad='white', blocks=True, extend='neither',
                xlabels=True):
    im = mesh(ax, M[order, :][:, order], cmap, norm, bad=bad)
    ax.set_xticks(np.arange(nR))
    ax.set_xticklabels(CRnames[order] if xlabels else [], rotation=90, fontsize=TICKFS)
    if not xlabels:
        ax.tick_params(axis='x', bottom=False)
    ax.set_yticks(np.arange(nR))
    ax.set_yticklabels(CRnames[order], fontsize=TICKFS)
    ax.set_title(title, fontsize=6)
    if blocks:
        draw_blocks(ax, labels[order])
    cb = plt.colorbar(im, cax=cax, label=cbar_label, extend=extend)
    cb.outline.set_linewidth(0.5)
    cb.ax.tick_params(labelsize=4)
    cb.solids.set_rasterized(False)   # matplotlib rasterizes the gradient by default
    return im

def two_panel_fig():
    fig = plt.figure(figsize=(4.3, 2.0))
    gs = fig.add_gridspec(1, 5, width_ratios=[20, 1, 9, 20, 1], wspace=0.08,
                          left=0.09, right=0.95, top=0.90, bottom=0.16)
    axL, caxL = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    axR, caxR = fig.add_subplot(gs[0, 3]), fig.add_subplot(gs[0, 4])
    fig.add_subplot(gs[0, 2]).axis('off')
    return fig, axL, caxL, axR, caxR

# %%
# --- co-expression: joint vs conditional reading, each with its own clustering ---
fig, axL, caxL, axR, caxR = two_panel_fig()

off = mask_diag(joint)
show_matrix(axL, caxL, off, order_joint, lab_joint,
            r'joint  $P(\sigma_i{=}1,\sigma_j{=}1)$', 'probability',
            norm=plt.Normalize(np.nanmin(off), np.nanmax(off)))

off = mask_diag(dice)
show_matrix(axR, caxR, off, order_dice, lab_dice,
            r'conditional  $2P_{ij}/(P_i{+}P_j)$', 'similarity',
            norm=plt.Normalize(np.nanmin(off), np.nanmax(off)))

plt.show()

# %%
# --- correlation: full range (diagonal kept) vs masked diagonal, rescaled ---
fig, axL, caxL, axR, caxR = two_panel_fig()

show_matrix(axL, caxL, corr_matrix, order_corr, lab_corr,
            'diagonal kept', r'$\rho_{ij}$',
            cmap='RdBu_r', norm=CenteredNorm(vcenter=0, halfrange=1.0))

CLIP_PCT = 100      # < 100 saturates the strongest pairs and boosts the contrast further
off = mask_diag(corr_matrix)
hr = np.nanpercentile(np.abs(off), CLIP_PCT)
show_matrix(axR, caxR, off, order_corr, lab_corr,
            rf'diagonal masked, $|\rho|\leq{hr:.2f}$', r'$\rho_{ij}$',
            cmap='RdBu_r', norm=CenteredNorm(vcenter=0, halfrange=hr))

plt.show()

# %%
# ---------------------------------------------------------
# FINAL FIGURE: correlation matrix, diagonal kept and saturated
# ---------------------------------------------------------
# Colour range set by the off-diagonal values only, so rho_ii = 1 clips to the top
# colour. TwoSlopeNorm keeps white exactly at rho = 0 while allowing the two
# half-ranges to differ, i.e. the lower bound is the smallest off-diagonal value.
off = mask_diag(corr_matrix)
vmin, vmax = np.nanmin(off), np.nanmax(off)
norm_final = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

fig, ax = plt.subplots(figsize=(2.1, 1.9))
# divider: the colorbar follows the height of the square matrix, no vertical gap
cax = make_axes_locatable(ax).append_axes('right', size='4%', pad=0.04)

show_matrix(ax, cax, corr_matrix, order_corr, lab_corr,
            'Gene-gene correlation', r'$\rho_{ij}$',
            cmap='RdBu_r', norm=norm_final, blocks=False, extend='max', xlabels=False)

plt.savefig('correlation_matrix_final.svg', bbox_inches='tight', pad_inches=0.01)
plt.show()

# %%
# ---------------------------------------------------------
# 4. Consistency between the clusterings
# ---------------------------------------------------------
def adjusted_rand(a, b):
    """Adjusted Rand Index: counts co-clustered PAIRS, so it ignores the label names."""
    ct = pd.crosstab(a, b).values.astype(float)
    n = ct.sum()
    sum_ij = comb(ct, 2).sum()
    sum_i = comb(ct.sum(axis=1), 2).sum()
    sum_j = comb(ct.sum(axis=0), 2).sum()
    exp = sum_i * sum_j / comb(n, 2)
    return (sum_ij - exp) / (0.5 * (sum_i + sum_j) - exp)

def matched_crosstab(a, b):
    """Contingency table with the columns permuted to maximise the diagonal (Hungarian)."""
    ct = pd.crosstab(a, b)
    _, c = linear_sum_assignment(-ct.values)
    rest = [j for j in range(ct.shape[1]) if j not in c]
    return ct.iloc[:, list(c) + rest]

pairs = [('joint', lab_joint, Z_joint, 'corr', lab_corr, Z_corr),
         ('cond', lab_dice, Z_dice, 'corr', lab_corr, Z_corr),
         ('joint', lab_joint, Z_joint, 'cond', lab_dice, Z_dice)]

k_range = np.arange(2, min(11, nR))
fig, axs = plt.subplots(1, 3, figsize=(5.4, 1.8))

for na, la, Za, nb, lb, Zb in pairs:
    ari_k = [adjusted_rand(sch.fcluster(Za, k, 'maxclust'),
                           sch.fcluster(Zb, k, 'maxclust')) for k in k_range]
    axs[0].plot(k_range, ari_k, '-o', lw=0.8, ms=2, label=f'{na} vs {nb}')
    print(f"ARI {na:>5} vs {nb:<5} at k={N_CLUST}: {adjusted_rand(la, lb):.3f}")

axs[0].axhline(0, color='gray', ls='--', lw=0.5)
axs[0].set_xlabel('number of clusters k')
axs[0].set_ylabel('ARI')
axs[0].set_ylim(-0.35, 1.05)
axs[0].legend(fontsize=5, frameon=False)

ct = matched_crosstab(lab_dice, lab_corr)
im = axs[1].imshow(ct.values, cmap='Blues')
axs[1].set_xticks(np.arange(ct.shape[1])); axs[1].set_xticklabels(ct.columns, fontsize=TICKFS)
axs[1].set_yticks(np.arange(ct.shape[0])); axs[1].set_yticklabels(ct.index, fontsize=TICKFS)
axs[1].set_xlabel('correlation cluster')
axs[1].set_ylabel('co-expression cluster')
axs[1].set_title(f'shared genes (k={N_CLUST}, matched)', fontsize=6)
for i in range(ct.shape[0]):
    for j in range(ct.shape[1]):
        axs[1].text(j, i, ct.values[i, j], ha='center', va='center', fontsize=TICKFS,
                    color='white' if ct.values[i, j] > ct.values.max() / 2 else 'black')

def partner_jaccard(a, b):
    """Per gene: Jaccard overlap of its cluster mates under labelling a and under b."""
    same_a = (a[:, None] == a[None, :]); np.fill_diagonal(same_a, False)
    same_b = (b[:, None] == b[None, :]); np.fill_diagonal(same_b, False)
    inter = (same_a & same_b).sum(axis=1)
    union = (same_a | same_b).sum(axis=1)
    return np.where(union > 0, inter / np.maximum(union, 1), 1.0)

jac = partner_jaccard(lab_dice, lab_corr)
axs[2].bar(np.arange(nR), jac[order_dice], color='steelblue', edgecolor='black', lw=0.3)
axs[2].set_xticks(np.arange(nR))
axs[2].set_xticklabels(CRnames[order_dice], rotation=90, fontsize=4)
axs[2].set_ylabel('cluster-mate Jaccard')
axs[2].set_title('per-gene agreement', fontsize=6)

plt.tight_layout()
plt.show()

print(pd.DataFrame({'gene': CRnames, 'P': prob_expr.round(3),
                    'c_joint': lab_joint, 'c_cond': lab_dice, 'c_corr': lab_corr}
                   ).sort_values(['c_cond', 'c_corr']).to_string(index=False))

# %%
