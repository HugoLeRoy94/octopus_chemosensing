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
    raise KeyError("Expected fields 'binaryTable_allGenes' and 'geneList' in the .mat file.")

A = np.array(md["binaryTable_allGenes"]).astype(bool).astype(int)  # shape: (n_receptors, n_cells)
CRnames = np.array(_to_str_list(md["geneList"]))
suckerIDX = np.array(_to_str_list(md["suckerIDX"]))
nR, nC = A.shape
# %%
# =========================================
# HISTOGRAM: NUMBER OF COEXPRESSED SUBUNITS
# compute the probability of expressing L different genes
# =========================================
coexp = A.sum(axis=0)  # per cell
coexp_nz = coexp[coexp != 0]
# integer bins: use bincount for exact integer histogram
max_k = int(coexp_nz.max())
counts = np.bincount(coexp_nz, minlength=max_k + 1)[1:]  # skip zero bin
centers = np.arange(1, max_k + 1)
Proba_L = counts/np.sum(counts)

# Fit with an exponential proba
def exp_distrib(l,beta):
    return beta*np.exp(-beta * l )

popt, pcov = curve_fit(exp_distrib, centers, Proba_L, p0=1, bounds=(0.1, 1))
print(popt)

# %%
fig,ax = plt.subplots(ncols=2,figsize=(6,3))
ax[0].plot(centers,exp_distrib(centers,popt[0]))
ax[0].plot(centers,Proba_L)

ax[1].plot(centers,exp_distrib(centers,popt[0]))
ax[1].plot(centers,Proba_L)
#ax[1].scatter(centers,Proba_L)
plt.yscale('log')

for i in range(2):
    ax[i].set_xlabel('Number of coexpressed receptors')
    ax[i].set_ylabel('Probability')

# %%
print(CRnames[id_sort])

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
# ===========================
# BATCH ANALYSIS - MOVING WINDOWS (cumulative windows like MATLAB)
# ===========================
M_cells = A.T
nCells = M_cells.shape[0]
batchEdge = np.arange(10, nCells + 1, 10)  # 10:10:nCells
nBatches = len(batchEdge)
uniquePerBatch = np.zeros(nBatches, dtype=int)
startIdx = 0  # Python 0-based

for b, edge in enumerate(batchEdge):
    stopIdx = min(startIdx + edge, nCells)  # cumulative from start
    if startIdx >= nCells:
        break
    sub = M_cells[startIdx:stopIdx, :]
    # unique rows count
    uniquePerBatch[b] = np.unique(sub, axis=0).shape[0]

plt.figure()
plt.plot(batchEdge, uniquePerBatch)
plt.xlabel("number of cells")
plt.ylabel("Unique # of patterns in batch")
plt.title("Overlapping batches")



# %%
