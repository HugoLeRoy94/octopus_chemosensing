# %%
"""G = 5, every design on one axis: fixed expression vs independent expression laws.

SELF-CONTAINED. Everything from "which folders am I reading" to "what is on each
axis" is in this file, so changing this figure cannot change any other.

Six curves, all at G = 5 genes, C = 30 cells:

  fixed size, complete coverage   heteromers / homomers   tasks/cells/gene_expression
      every cell expresses EXACTLY g genes, and every gene is guaranteed to
      appear somewhere in the array. Requires cells to coordinate.

  independent, uniform law        heteromers / homomers   tasks/cells/expression_law
  independent, exponential law    heteromers / homomers   tasks/cells/expression_law
      every cell draws its own gene set alone: first how many genes (size law),
      then which ones (gene law). Nothing guarantees a gene is expressed at all.

## The one caveat to state out loud when presenting

The two tasks were measured by DIFFERENT estimators, because independent
sampling leaves too few identical cells for exact enumeration:

    gene_expression   exact count enumeration, unbiased by construction
    expression_law    sampled counting over 16,777,216 inputs

tasks/cells/gene_expression/scripts/estimator_crosscheck.py measured what that
costs. Sampled counting under-reads MI by 3.34*(B/2^H(K))^-0.92 bits, which at
this budget is under 0.02 bits for every G=5 point plotted here. The bottom
right panel shows the Good-Turing missing mass f1/B that certifies it. So the
comparison is sound, but it rests on that measurement rather than on assumption.
"""
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def _project_root() -> Path:
    """Walk upward until opt_bin_resp/src is found.

    Not Path(__file__).parents[N]: a fixed hop count breaks in an interactive
    kernel, where __file__ is relative and resolves against the kernel's working
    directory rather than this file's location. Searching upward works whether the
    script is run directly, imported, or pasted into a cell.
    """
    here = Path(__file__).resolve().parent if "__file__" in globals() else Path.cwd()
    for candidate in (here, *here.parents):
        if (candidate / "opt_bin_resp" / "src" / "IO.py").is_file():
            return candidate
    raise SystemExit(
        f"Cannot find opt_bin_resp/src above {here}. Run this from inside the "
        "octopus_smelling tree, or set ROOT by hand below.")


ROOT = _project_root()
sys.path.insert(0, str(ROOT / "opt_bin_resp"))
from src.IO import SweepLoader

DATA = ROOT / "opt_bin_resp" / "data"
HERE = ROOT / "pres" / "group_meeting"
SAVE_FIGURES = False
N_GENES = 5


# %%
# ════════════════════════════════════════════════════════════════════════════
# 1. WHICH SWEEPS
# ════════════════════════════════════════════════════════════════════════════
# Named explicitly, with the estimator each one used. "exact" runs report
# mutual_information_grouped; "counting" runs report the plug-in of the sampled
# estimator and additionally carry the missing mass that qualifies it.
SWEEPS = [
    # goal, folder, design label, condition, estimator
    ("gene_expression", "replicates_heteromers_complete_20260923_112537",
     "fixed size, complete coverage", "heteromers", "exact"),
    ("gene_expression", "replicates_homomers_complete_20260923_112624",
     "fixed size, complete coverage", "homomers", "exact"),
    ("expression_law", "expression_law_heteromers_uniform_uniform_20260924_151241",
     "independent, uniform law", "heteromers", "counting"),
    ("expression_law", "expression_law_homomers_uniform_uniform_20260924_151313",
     "independent, uniform law", "homomers", "counting"),
    ("expression_law", "expression_law_heteromers_exponential_exponential_20260924_190805",
     "independent, exponential law", "heteromers", "counting"),
    ("expression_law", "expression_law_homomers_exponential_exponential_20260924_190829",
     "independent, exponential law", "homomers", "counting"),
]

MI_KEY = {"exact": "mutual_information_grouped",
          "counting": "mutual_information_grouped_counting_plugin"}
NOISE_KEY = {"exact": "conditional_entropy_response_grouped",
             "counting": "conditional_entropy_response_grouped_counting"}

missing = [f"{goal}/{name}" for goal, name, *_ in SWEEPS
           if not (DATA / goal / name / "experiment.json").is_file()]
if missing:
    raise SystemExit("Not on disk:\n  " + "\n  ".join(missing))


# %%
# ════════════════════════════════════════════════════════════════════════════
# 2. ONE ROW PER RUN   (first average: the final measurement -> 1 number)
# ════════════════════════════════════════════════════════════════════════════
records = []
for goal, name, design, condition, estimator in SWEEPS:
    sweep_dir = DATA / goal / name
    # The REQUESTED expression level, keyed by the seed that produced the array.
    # Replicates of one design point land on slightly different realised means, so
    # grouping by the realised value would shatter each point into singletons with
    # no error bar. gene_expression has no such field: its target IS the integer.
    manifest = json.loads((sweep_dir / "experiment.json").read_text())
    target = {r["cell_sampling_seed"]: r.get("target_mean_genes")
              for r in manifest["rows"]}
    kept = 0
    for cfg, run_dir in SweepLoader(str(sweep_dir)).iter_run_dirs():
        if cfg.n_genes != N_GENES:
            continue                      # this figure is G=5 only
        results_path = Path(run_dir) / "test_results.json"
        if not results_path.exists():
            continue
        results = json.loads(results_path.read_text())
        sizes = [len(genes) for genes in cfg.cell_gene_sets]
        records.append({
            "design": design, "condition": condition, "estimator": estimator,
            "C": cfg.n_cells,
            # genes per cell, measured from the array actually drawn. The fixed
            # design has one value for every cell; the independent laws do not.
            "genes_per_cell": float(np.mean(sizes)),
            "spread": float(np.std(sizes)),
            "seed": cfg.cell_sampling_seed,
            "target": target.get(cfg.cell_sampling_seed) or float(np.mean(sizes)),
            "mi": float(np.mean(results[MI_KEY[estimator]])),
            "noise": float(np.mean(results[NOISE_KEY[estimator]])),
            "cell_types": len({tuple(g) for g in cfg.cell_gene_sets}),
            "genes_seen": len({g for cell in cfg.cell_gene_sets for g in cell}),
            "missing_mass": float(np.mean(results["grouped_counting_missing_mass"]))
            if "grouped_counting_missing_mass" in results else np.nan,
        })
        kept += 1
    print(f"{name}: {kept} runs at G={N_GENES}")

runs = pd.DataFrame(records)
if runs.empty:
    raise SystemExit(f"No completed G={N_GENES} runs found.")


# %%
# ════════════════════════════════════════════════════════════════════════════
# 3. ONE ROW PER POINT   (second average: independent optimizations -> 1 number)
# ════════════════════════════════════════════════════════════════════════════
# Both tasks draw their gene sets once and hand the SAME arrays to the heteromer
# and homomer arms, so a design point is matched across conditions and the
# advantage below is paired. Rounding only collapses float noise in the realised
# mean; it does not merge distinct design points.
CURVE = ["design", "condition"]

points = runs.groupby(CURVE + ["target"]).agg(
    n=("mi", "size"), mi_mean=("mi", "mean"), mi_sd=("mi", "std"),
    # plotted against what the arrays actually came out at, averaged over replicates
    x=("genes_per_cell", "mean"),
    noise=("noise", "mean"), types=("cell_types", "mean"), C=("C", "max"),
    genes_seen=("genes_seen", "mean"), missing_mass=("missing_mass", "mean"),
    spread=("spread", "mean"),
).reset_index()
points["mi_sem"] = points.mi_sd / np.sqrt(points.n)

baseline = points[points.target == 1.0].set_index(CURVE).mi_mean
points["retained"] = points.mi_mean / pd.Index(
    list(map(tuple, points[CURVE].values))).map(baseline)

# Heteromer advantage, paired within a design at the same expression level.
het = points[points.condition == "heteromers"].set_index(["design", "target"])
hom = points[points.condition == "homomers"].set_index(["design", "target"])
shared = het.index.intersection(hom.index)
advantage = pd.DataFrame({
    "advantage": het.loc[shared].mi_mean - hom.loc[shared].mi_mean,
    "advantage_sem": np.sqrt(het.loc[shared].mi_sem ** 2 + hom.loc[shared].mi_sem ** 2),
    "x": het.loc[shared].x,
}).reset_index()


# %%
# ════════════════════════════════════════════════════════════════════════════
# 4. WHAT IS ABOUT TO BE DRAWN
# ════════════════════════════════════════════════════════════════════════════
for (design, condition), curve in points.groupby(CURVE, sort=False):
    curve = curve.sort_values("x")
    print(f"\n{condition:>10} | {design}")
    print(f"    genes/cell : {[round(v, 2) for v in curve.x]}")
    print(f"    MI [bits]  : {[round(v, 3) for v in curve.mi_mean]}")
    print(f"    +- sem     : {[round(v, 3) for v in curve.mi_sem]}")
    print(f"    n runs     : {list(curve.n)}")
    print(f"    cell types : {[round(v, 1) for v in curve.types]} of {int(curve.C.iloc[0])}")
    if curve.missing_mass.notna().any():
        print(f"    f1/B       : {[round(v, 5) for v in curve.missing_mass]}")

print("\n\nHeteromer advantage, paired at the same expression level:")
for design, part in advantage.groupby("design", sort=False):
    part = part.sort_values("x")
    print(f"  {design}")
    print(f"    x    : {[round(v, 2) for v in part.x]}")
    print(f"    bits : {[round(v, 3) for v in part.advantage]}"
          f"  +- {[round(v, 3) for v in part.advantage_sem]}")


# %%
# %%
# ════════════════════════════════════════════════════════════════════════════
# 5. THREE SLIDE FIGURES, BUILT UP IN SEQUENCE
# ════════════════════════════════════════════════════════════════════════════
# Same size, same font, and the SAME AXIS LIMITS in all three, so that clicking
# through them in a deck adds curves without the axes jumping underneath.
plt.rcParams.update({
    "font.size": 16, "axes.titlesize": 16, "axes.labelsize": 16,
    "xtick.labelsize": 16, "ytick.labelsize": 16, "legend.fontsize": 14,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": False, "ytick.right": False,
    "xtick.major.size": 7, "ytick.major.size": 7,
    #"xtick.minor.size": 4, "ytick.minor.size": 4,
    "axes.linewidth": 1.2, "lines.linewidth": 2.2, "lines.markersize": 8,
})
FIGSIZE = (8, 5.5)

COLOR = {"heteromers": "tab:blue", "homomers": "tab:orange"}
DASH = {"fixed size, complete coverage": "-",
        "independent, uniform law": "--",
        "independent, exponential law": ":"}
MARK = {"fixed size, complete coverage": "o",
        "independent, uniform law": "s",
        "independent, exponential law": "^"}
LABEL = {"fixed size, complete coverage": "fixed size, every gene covered",
         "independent, uniform law": "independent, uniform law",
         "independent, exponential law": "independent, exponential law"}

# Common limits, from every point any of the three figures can show.
_pad_x = 0.05 * (points.x.max() - points.x.min())
_lo = (points.mi_mean - points.mi_sem.fillna(0)).min()
_hi = (points.mi_mean + points.mi_sem.fillna(0)).max()
_pad_y = 0.08 * (_hi - _lo)
XLIM = (points.x.min() - _pad_x, points.x.max() + _pad_x)
YLIM = (_lo - _pad_y, _hi + _pad_y)


def slide(selection, title):
    """One presentation panel. `selection` is a list of (design, condition)."""
    fig, ax = plt.subplots(figsize=FIGSIZE)
    # When every curve is the same condition, saying so on each legend entry just
    # eats width. The TITLE must then name the condition, or the slide loses it.
    one_condition = len({condition for _, condition in selection}) == 1
    for design, condition in selection:
        curve = points[(points.design == design) & (points.condition == condition)]
        curve = curve.sort_values("x")
        label = LABEL[design] if one_condition else f"{condition}, {LABEL[design]}"
        ax.errorbar(curve.x, curve.mi_mean, yerr=curve.mi_sem,
                    color=COLOR[condition], linestyle=DASH[design],
                    marker=MARK[design], capsize=4, label=label)
    ax.set_xlabel("genes expressed per cell")
    ax.set_ylabel("mutual information [bits]")
    ax.set_xlim(*XLIM)
    ax.set_ylim(*YLIM)
    #ax.minorticks_on()
    ax.set_title(title)
    ax.legend(frameon=False, loc="lower left")
    fig.tight_layout()
    return fig


FIXED = "fixed size, complete coverage"
UNIFORM = "independent, uniform law"
EXPONENTIAL = "independent, exponential law"

# %%
# 1 of 3 -- the baseline: one receptor type per gene.
fig = slide([(FIXED, "homomers")], "homomers, G = 5 genes, C = 30 cells")
fig.savefig(HERE / 'homomers.svg')

# %%
# 2 of 3 -- add heteromers on the same arrays. The gap is the combinatorial gain.
fig = slide([(FIXED, "homomers"), (FIXED, "heteromers")], "G = 5 genes, C = 30 cells")
fig.savefig(HERE / 'homomers_heteromers.svg')

# %%
# 3 of 3 -- heteromers only, now letting each cell choose its genes alone.
# Drop the FIXED entry from this list if the reference curve is not wanted.
fig = slide([(FIXED, "heteromers"), (UNIFORM, "heteromers"), (EXPONENTIAL, "heteromers")],
      "heteromers: fixed vs independent expression")
fig.savefig(HERE / 'heteromers_law.svg')


# %%
# Not for the deck: the measurement check behind the independent-law curves.
SHOW_DIAGNOSTICS = True
if SHOW_DIAGNOSTICS:
    fig, ax = plt.subplots(figsize=FIGSIZE)
    for (design, condition), curve in points.groupby(CURVE, sort=False):
        curve = curve.sort_values("x")
        if curve.missing_mass.isna().all():
            continue
        ax.plot(curve.x, curve.missing_mass, color=COLOR[condition],
                linestyle=DASH[design], marker=MARK[design],
                label=f"{condition}, {LABEL[design]}")
    ax.set(xlabel="genes expressed per cell",
           ylabel="Good-Turing missing mass  f1/B", yscale="log")
    ax.minorticks_on()
    ax.legend(frameon=False, fontsize=11)
    fig.tight_layout()
    plt.show()

# %%
