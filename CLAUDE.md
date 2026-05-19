# Octopus Chemosensing Project - AI Assistant Rules

## Project Overview
This project simulates the chemo-sensory receptor array of octopuses. We optimize the mutual information between an array of combinatorial ion channels (hetero-pentamers) and a high-dimensional chemical latent space representing ligand families.

## Tech Stack & Coding Constraints
- **Language:** Python 3.10+
- **Deep Learning:** PyTorch.
- **Performance:** Strictly use vectorized tensor operations. Avoid Python `for` loops for mathematical operations.
- **Environment:** The code runs inside a headless Docker container on a remote A100 GPU cluster.

## Documentation Map
If you need specific mathematical or biophysical context, refer to the files in `/mnt/hcleroy/PostDoc2/octopus_smelling/doc/theory/`. Do NOT hallucinate equations:
- **Nomenclature & Vars:** `01_nomenclature.md` (Definitions, Vocab, Experimentalist mapping).
- **Receptor Physics:** `02_biophysics_mwc.md` (MWC model equations, microscopic interactions).
- **Environment & Latent Space:** `03_latent_environment.md` (Embeddings, distances, energy sampling).
- **Discrete Info Theory:** `04_discrete_information.md` (Joint entropy, ligand/concentration MI decomposition).
- **Discrete Optimization:** `05_optimization.md` (Continuous relaxation, Rényi entropy proxy).
- **Computational Limits:** `06_computational_limits.md` (Memory footprint scaling, tensor bottlenecks, algorithmic fallbacks).
- **Academic Context:** `/mnt/hcleroy/PostDoc2/octopus_smelling/doc/tex/prgs_rprt.tex`

## Documentation Maintenance Rule
Whenever you modify any file in `opt_bin_resp/src/`, you MUST:
1. Check whether `doc/theory/07_optimization_pipeline.md` needs updating for the changed logic.
2. Check whether any referenced theory file (`01`–`06`) is affected.
3. Update the relevant doc(s) before reporting the task complete.
Do not skip this even for "small" changes — a renamed parameter or new estimator invalidates the pipeline doc.

## Mathematical Definitions Reference
- $k_{sub}$: Number of sub-units in a receptor (typically 5).
- $EC_{50}$: The concentration at half activation. For a heteromer, this is the geometric mean of its individual subunit affinities.
- $T$: Temperature parameter used in the continuous relaxation of the step-function (Soft Histogram).