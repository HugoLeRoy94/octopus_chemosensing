# Octopus Chemosensing Project - AI Assistant Rules

## Git: NEVER commit unless explicitly told to

Do NOT run `git commit`, `git push`, or any history-rewriting git command unless
the user asks for it **in that message**. A previous "commit this" covers that one
batch only and expires immediately; it is never standing permission.

Finish the work, run the tests, report what changed, and leave it as unstaged
changes in the working tree. Offering to commit is fine. Acting on the offer
before being told is not.

## Figures: NEVER save to disk unless explicitly told to

Do NOT call `plt.savefig` / `fig.savefig`, and do not write image files, unless the
user asks for it **in that message**. Build the figure and show it. If a script
needs a save path at all, put it behind a flag that defaults to off.

The user decides what lands in the repo. Writing files they did not ask for
clutters the tree and forces them to clean up after you.

## Environments: NEVER fix the environment to make a comparison clean

Do NOT propose, and do not implement, matching or sharing the environment seed
across conditions, replicates, or sweeps. Not to make a comparison "paired", not
to isolate one variable, not as a control.

The scientific claim is about a **category** of environments, not about any one
of them. A result that holds only on a shared ligand world is a property of that
world, not of the array, and is worthless. Results must be statistically reliable
across the category and independent of its details.

When two conditions disagree and a fixed environment would "explain" it, the
correct response is **more replicates across more environments**, and a
comparison of distributions rather than of two numbers. Report the spread across
independent optimizations, never the spread across test repeats of one model.

This is permanent. Do not raise it again.

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
- **Cell Arrays:** `09_cell_arrays.md` (`src/cells.py` reference: gene-set sampling, stoichiometry, readout and calibration, and the justification for each choice).
- **Academic Context:** `/mnt/hcleroy/PostDoc2/octopus_smelling/doc/tex/prgs_rprt.tex`

## Data Map
Questions about WHERE A NUMBER CAME FROM (not about physics) are answered in
`/mnt/hcleroy/PostDoc2/octopus_smelling/opt_bin_resp/doc/`:
- **Disk to plot:** `data_pipeline.md` — READ THIS FIRST for anything about stored
  data or analysis. Covers the goal/sweep/run directory layout, what each of
  `sweep_config.json`, `config.json`, `stats.csv`, `test_results.json` and
  `experiment.json` contains and who writes it, what `runs.db` is and which
  analyses actually read it, cluster synchronization, and a traced worked example
  from a file on disk to one plotted point with an error bar.
- **Curate & sync commands:** `curation_and_sync.md` (`manage_data.py curate` /
  `sync`, `curation.csv`, filtering an analysis by curation state).

Use it when asked: what is being plotted, where does this value come from, why do
two sweeps disagree, what do the 10 values in `test_results.json` mean (evaluation
noise on ONE model, never replicates), how is data pulled from the cluster.
Do NOT infer the data layout from directory names alone; `data_pipeline.md`
records which files are ground truth and which are regenerable caches.

## Documentation Maintenance Rule
Whenever you modify any file in `opt_bin_resp/src/`, you MUST:
1. Check whether `doc/theory/07_optimization_pipeline.md` needs updating for the changed logic.
2. Check whether any referenced theory file (`01`–`06`) is affected.
3. If you touched `IO.py`, `db.py`, `plotlib.py`, `manage_data.py`, or anything under
   `opt_bin_resp/tasks/*/analysis/`, check `opt_bin_resp/doc/data_pipeline.md`: it
   names specific functions and quotes real numbers, so it goes stale silently.
4. Update the relevant doc(s) before reporting the task complete.
Do not skip this even for "small" changes — a renamed parameter or new estimator invalidates the pipeline doc.

## Mathematical Definitions Reference
- $k_{sub}$: Number of sub-units in a receptor (typically 5).
- $EC_{50}$: The concentration at half activation. For a heteromer, this is the geometric mean of its individual subunit affinities.
- $T$: Temperature parameter used in the continuous relaxation of the step-function (Soft Histogram).