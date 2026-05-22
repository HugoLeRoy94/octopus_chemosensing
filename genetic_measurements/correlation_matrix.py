# %%
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
from scipy.optimize import curve_fit
import itertools
import pandas as pd
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as ssd

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

#directory = "/mnt/hcleroy/PostDoc2/octopus_smelling/data/matlab_files/v2"
#mat_path = directory+"/20260120_HiPlexResults.mat"
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

# %%
# ---------------------------------------------------------
# 1. Plot the probability expression of each individual gene
# ---------------------------------------------------------
# Calculate the mean across the columns (cells) for each row (gene)
prob_expr = A.mean(axis=1)
id_sort = np.argsort(-prob_expr)

plt.plot(range(1,id_sort.shape[0]+1),prob_expr[id_sort])
plt.xticks(np.arange(nR), CRnames[id_sort],rotation=45)
plt.ylabel('probability of expression')
plt.show()

# %%

# ---------------------------------------------------------
# 2. Find the first few probable gene patterns
# ---------------------------------------------------------
# Transpose A so that each row is a cell and each column is a gene
# Then convert it into a pandas DataFrame for easy frequency counting
df = pd.DataFrame(A.T)

# Count the occurrences of each unique pattern (row)
pattern_counts = df.value_counts().reset_index(name='count')

# Calculate the probability of each pattern
pattern_counts['probability'] = pattern_counts['count'] / A.shape[1]

# Display the top 10 most frequent gene patterns
print("Top 10 most probable gene patterns:")
top_10_patterns = pattern_counts.head(10)
print(top_10_patterns.to_string(index=False))

# %%

# ---------------------------------------------------------
# 3. Gene expression probability conditioned on total expressed genes
# ---------------------------------------------------------
total_expressed = A.sum(axis=0)  # (n_cells,): number of genes on per cell
k_values = np.arange(nR + 1)

# k_mask[k, c] = 1 iff cell c has exactly k genes expressed
k_mask = (total_expressed[np.newaxis, :] == k_values[:, np.newaxis])  # (nR+1, n_cells)
k_counts = k_mask.sum(axis=1)                                         # (nR+1,)

# cond_prob[i, k] = P(gene i expressed | total expressed == k)
safe_counts = np.where(k_counts > 0, k_counts, np.nan)
cond_prob = (A @ k_mask.T) / safe_counts[np.newaxis, :]  # (nR, nR+1)

# enrichment[i, k] = P(i | L=k) / (k/N)  — ratio to flat expectation; k=0 → nan
flat_expectation = np.where(k_values > 0, k_values / nR, np.nan)  # (nR+1,)
enrichment = cond_prob / flat_expectation[np.newaxis, :]           # (nR, nR+1)

from matplotlib.colors import TwoSlopeNorm,CenteredNorm
enrich_max = np.nanmax(np.abs(enrichment - 1)) + 1
norm_enrich = CenteredNorm(vcenter=1,halfrange=2)#TwoSlopeNorm(vcenter=1, vmin=0, vmax=enrich_max)

cmap_enrich = plt.get_cmap('RdBu_r').copy()
cmap_enrich.set_bad(color='lightgray')

fig, axs = plt.subplots(2, 1, figsize=(14, 8), gridspec_kw={'height_ratios': [1, 3]}, sharex=True)

axs[0].bar(k_values, k_counts, color='steelblue', edgecolor='black')
axs[0].set_ylabel('Cell count')
axs[0].set_title('Distribution of total expressed genes per cell')

im = axs[1].imshow(enrichment[id_sort, :], aspect='auto', cmap=cmap_enrich, norm=norm_enrich)
axs[1].set_yticks(np.arange(nR))
axs[1].set_yticklabels(CRnames[id_sort], fontsize=7)
axs[1].set_xticks(k_values)
axs[1].set_xlabel('Number of expressed genes in cell (k)')
axs[1].set_ylabel('Gene')
axs[1].set_title('Enrichment  P(gene | L=k) / (k/N)  — white = flat expectation')

plt.colorbar(im, ax=axs[1], label='Enrichment (1 is no enrichement)')
plt.tight_layout()
plt.savefig('P_conditionned.png')

# %%

# E[L | σ_i = 1]: expected total number of expressed genes, given gene i is on
# Vectorized: sum of total_expressed over cells where gene i is on, divided by count
gene_on_counts = A.sum(axis=1).astype(float)          # (nR,): number of cells where gene i is on
mean_L_given_on = (A @ total_expressed) / gene_on_counts  # (nR,)

fig, ax = plt.subplots(figsize=(12, 4))
ax.bar(np.arange(nR), mean_L_given_on[id_sort], color='steelblue', edgecolor='black')
ax.axhline(total_expressed.mean(), color='red', linestyle='--', lw=1.5, label=f'global ⟨L⟩ = {total_expressed.mean():.2f}')
ax.set_xticks(np.arange(nR))
ax.set_xticklabels(CRnames[id_sort], rotation=45, fontsize=7)
ax.set_ylabel('E[L | σᵢ = 1]')
ax.set_title('Expected number of expressed genes conditioned on each gene being on')
ax.legend()
plt.tight_layout()
plt.savefig('E_L_given_on.png')
plt.show()

# %%

# Line plot of conditional probabilities, one curve per k, colored by k
valid_k = np.where(k_counts > 0)[0]
cmap_k = plt.get_cmap('plasma')
norm_k = plt.Normalize(vmin=valid_k.min(), vmax=valid_k.max())

fig, ax = plt.subplots(figsize=(14, 4))
for k in range(10):#valid_k:
    ax.plot(range(nR), cond_prob[id_sort, k], color=cmap_k(norm_k(k)), alpha=0.8, lw=1.5)

sm = plt.cm.ScalarMappable(cmap=cmap_k, norm=norm_k)
plt.colorbar(sm, ax=ax, label='k (total expressed genes)')

ax.set_xticks(np.arange(nR))
ax.set_xticklabels(CRnames[id_sort], rotation=45, fontsize=7)
ax.set_ylabel('P(gene expressed | total expressed = k)')
ax.set_title('Conditional gene expression probability per k')
plt.tight_layout()
plt.show()

# %%

# 1. Pairwise Correlation Matrix (24x24)
# Values near 0 indicate independence between that pair of genes.
corr_matrix = np.corrcoef(A)

# 2. Observed vs Expected Co-expression (Example for Gene 0 and Gene 1)
gene_A_idx, gene_B_idx = 0, 1

# Expected joint probability if strictly independent
expected_joint_prob = prob_expr[gene_A_idx] * prob_expr[gene_B_idx]

# Actual observed joint probability
observed_joint_prob = np.mean((A[gene_A_idx, :] == 1) & (A[gene_B_idx, :] == 1))

print(f"Expected: {expected_joint_prob:.4f}, Observed: {observed_joint_prob:.4f}")

# %%


fig,ax = plt.subplots(ncols=3,figsize=(15,5))

############################GROUP BY CORRELATION######################
# Force strict symmetry to fix floating-point errors
corr_matrix = (corr_matrix + corr_matrix.T) / 2
#############################LOOK AT THE COOEXPRESSION PROBA###############""
# Compute individual probabilities
prob_expr = A.mean(axis=1)
num_cells = A.shape[1]
co_expr = (A @ A.T) / num_cells


# 1. Convert correlation to a distance metric
#dist_matrix = 1 - corr_matrix
dist_matrix = max(co_expr.flatten()) - co_expr
np.fill_diagonal(dist_matrix, 0)
# Now squareform will accept it
condensed_dist = ssd.squareform(dist_matrix)
# 2. Compute hierarchical clustering
linkage_matrix = sch.linkage(condensed_dist, method='average')
# 3. Get the new ordered indices
reordered_idx = sch.leaves_list(linkage_matrix)


#reordered_idx[14], reordered_idx[19] = reordered_idx[19], reordered_idx[14]
#reordered_idx[15], reordered_idx[20] = reordered_idx[20], reordered_idx[15]
#reordered_idx[16], reordered_idx[21] = reordered_idx[21], reordered_idx[16]
#reordered_idx[11], reordered_idx[8] = reordered_idx[8], reordered_idx[11]

####################### PLOT ############################
im = ax[1].imshow(corr_matrix[reordered_idx, :][:, reordered_idx])
cbar = plt.colorbar(im)
ax[1].set_xticks(np.arange(nR), CRnames[reordered_idx], rotation=90,fontsize=7)
ax[1].set_yticks(np.arange(nR), CRnames[reordered_idx],fontsize=7)
ax[1].set_title('correlation matrix')


# Color plot
im2 = ax[0].imshow(co_expr[reordered_idx, :][:, reordered_idx])
cbar = plt.colorbar(im2)
ax[0].set_xticks(np.arange(nR), CRnames[reordered_idx], rotation=90,fontsize=7)
ax[0].set_yticks(np.arange(nR), CRnames[reordered_idx],fontsize=7)
ax[0].set_title('co-expression level')

ax[2].plot(range(1,reordered_idx.shape[0]+1),prob_expr[reordered_idx])
ax[2].set_xticks(np.arange(nR), CRnames[reordered_idx],rotation=45,fontsize=7)
ax[2].set_title('single gene expression proba')

#starts = np.array([8.,12.,21.])
#ends = np.array([14.,21,23.])
starts = np.array([8.,9.,12.,14.,15.,19.])
ends = np.array([9.,14.,14.,21,19,21.])
starts+=0.5
ends+=0.5
for i in range(2):
    for s,e in zip(starts,ends):
        ax[i].hlines([s,e],s,e,color='black')
        ax[i].vlines([s,e],s,e,color='black')
ax[2].axvline([9.],color='black')
ax[2].axvline([15.],color='black')
ax[2].axvline([13.],color='green')
ax[2].axvline([22.-0.1],color='green')
ax[2].axvline([22.+0.1],color='red')
ax[2].axvline([24.],color='red')

plt.show()
# %%

# ---------------------------------------------------------
# 3. Gene expression probability across different suckers
# ---------------------------------------------------------
from scipy.stats import chi2_contingency
import re

# Extract integer IDs from suckerIDX (handles cases like '1', 'sucker_1', etc.)
sucker_idx_int = np.array([int(re.search(r'\d+', s).group()) if re.search(r'\d+', s) else -1 for s in suckerIDX])

ordered_suckers = np.arange(1, 28)
num_suckers = len(ordered_suckers)

sucker_counts = np.zeros(num_suckers, dtype=int)
expr_by_sucker_raw = np.zeros((nR, num_suckers))
expr_by_sucker_rel = np.zeros((nR, num_suckers))

# Global probability of expression for each gene (baseline)
global_prob = A.mean(axis=1)

print("Chi-square test for independence between Gene Expression and Sucker:")
for i in range(nR):
    contingency_table = np.zeros((2, num_suckers))
    for j, sucker in enumerate(ordered_suckers):
        idx_sucker = (sucker_idx_int == sucker)
        sucker_counts[j] = np.sum(idx_sucker)
        
        if sucker_counts[j] > 0:
            raw_prob = A[i, idx_sucker].mean()
            expr_by_sucker_raw[i, j] = raw_prob
            expr_by_sucker_rel[i, j] = raw_prob - global_prob[i]
            
            # Populate contingency table for statistical test
            contingency_table[1, j] = np.sum(A[i, idx_sucker] == 1)
            contingency_table[0, j] = np.sum(A[i, idx_sucker] == 0)
        else:
            # No cells in this sucker
            expr_by_sucker_raw[i, j] = np.nan
            expr_by_sucker_rel[i, j] = np.nan
    
    # Test if gene expression significantly varies across suckers
    # We drop empty columns (suckers with 0 cells) to avoid ValueError
    valid_cols = sucker_counts > 0
    if np.sum(valid_cols) > 1:
        try:
            chi2, p, dof, ex = chi2_contingency(contingency_table[:, valid_cols])
            if p < 0.05:
                print(f"  - Gene {CRnames[i]:>8}: p-value = {p:.4e} (Significant)")
        except ValueError:
            pass

from mpl_toolkits.axes_grid1 import make_axes_locatable

fig, axs = plt.subplots(3, 1, figsize=(12, 10), gridspec_kw={'height_ratios': [1, 2, 2]}, sharex=True)

# Top plot: Sucker cell counts
axs[0].bar(np.arange(num_suckers), sucker_counts, color='skyblue', edgecolor='black')
axs[0].set_ylabel('Cell Count')
axs[0].set_title('Cell Count per Sucker')
for i, count in enumerate(sucker_counts):
    axs[0].text(i, count + max(sucker_counts)*0.02, str(count), ha='center', va='bottom', fontsize=8)

divider0 = make_axes_locatable(axs[0])
cax0 = divider0.append_axes("right", size="2%", pad=0.1)
cax0.axis('off') # Invisible colorbar axis to align X axis perfectly with the heatmaps below

# Middle plot: Raw Probability Heatmap
cmap_raw = plt.get_cmap('viridis').copy()
cmap_raw.set_bad(color='lightgray')
im1 = axs[1].imshow(expr_by_sucker_raw[reordered_idx, :], aspect='auto', cmap=cmap_raw, vmin=0, vmax=1)

axs[1].set_yticks(np.arange(nR))
axs[1].set_yticklabels(CRnames[reordered_idx], fontsize=7)
axs[1].set_ylabel('Gene')
axs[1].set_title('Probability of Expression')

divider1 = make_axes_locatable(axs[1])
cax1 = divider1.append_axes("right", size="2%", pad=0.1)
cbar1 = plt.colorbar(im1, cax=cax1)
cbar1.set_label('Prob')

# Bottom plot: Relative Gene Expression Heatmap
vmax = np.nanmax(np.abs(expr_by_sucker_rel))
if np.isnan(vmax) or vmax == 0:
    vmax = 1.0 # Fallback safety

cmap = plt.get_cmap('RdBu_r').copy()
cmap.set_bad(color='lightgray')
im2 = axs[2].imshow(expr_by_sucker_rel[reordered_idx, :], aspect='auto', cmap=cmap, vmin=-vmax, vmax=vmax)

axs[2].set_xticks(np.arange(num_suckers))
axs[2].set_xticklabels(ordered_suckers, fontsize=10)
axs[2].set_yticks(np.arange(nsR))
axs[2].set_yticklabels(CRnames[reordered_idx], fontsize=7)
axs[2].set_xlabel('Sucker IDX')
axs[2].set_ylabel('Gene')
axs[2].set_title('Relative Expression (Sucker Prob - Global Prob)')

divider2 = make_axes_locatable(axs[2])
cax2 = divider2.append_axes("right", size="2%", pad=0.1)
cbar2 = plt.colorbar(im2, cax=cax2)
cbar2.set_label('Relative Prob')

plt.tight_layout()
plt.show()
# %%
[[s,c] for s, c in zip(unique_suckers, sucker_counts)]