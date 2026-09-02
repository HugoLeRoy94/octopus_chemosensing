# The Cell Picture — status, findings, open questions

Working note, written 2026-08-06. Summarises the shift from a receptor array to a
cell array: what was decided, what is implemented, what broke, and what to run next.

Code lives in `opt_bin_resp/src/cells.py` (new) plus edits to `physics.py`,
`environment.py`, `config.py`, `run.py`, `IO.py`, `analysis_helper.py`.
Theory write-up is in `doc/theory/07_optimization_pipeline.md` §3b, with the
scaling analysis in `doc/theory/06_computational_limits.md` and vocabulary in
`doc/theory/01_nomenclature.md`.

---

## 1. What changed conceptually

**Before:** the array was a list of receptors. Each receptor was 5 gene indices. The
physics gave, per sniff and per receptor, a probability of being open. Entropy was
computed over that table.

**Now:** the array is a list of **cells**. A cell is defined only by *which genes it
expresses*, and it assembles **every receptor those genes allow**.

The output table has the same shape — sniffs × channels — so every loss and every
measurement downstream works unchanged. The two pictures differ at exactly one point
in the code: `SimulationRunner._activity`.

### Why this is worth doing

45 *cells* is a much richer system than 45 *receptors*, because those cells draw on a
pool of several thousand distinct receptors underneath. The entropy estimator only
ever sees the cells, so the combinatorial richness is free on the estimator side.

---

## 2. The pieces

### 2.1 Gene set → receptor repertoire

A cell expressing genes {0, 2} with 5 subunit slots builds 6 receptors:

```
[0,0,0,0,0]  [0,0,0,0,2]  [0,0,0,2,2]  [0,0,2,2,2]  [0,2,2,2,2]  [2,2,2,2,2]
```

Growth in the number of expressed genes `g`:

| g | receptors (standard) | receptors (interface model) |
|---|---|---|
| 1 | 1 | 1 |
| 2 | 6 | 8 |
| 3 | 21 | 51 |
| 5 | 126 | 629 |
| 8 | 792 | 6560 |

The interface model counts ring arrangements (rotations identified, reflections
distinct), reusing the existing `_canonical_rotation` from `geometry.py`.

### 2.2 The pool

Every receptor any cell can build, **deduplicated**. A receptor shared by two cells has
the same physics either way, so it is simulated once. That deduplicated list is what
gets fed to the existing physics code, which has no idea cells exist.

Measured pool sizes with every cell expressing exactly 5 genes:

| available genes | 10 cells | 20 cells | 50 cells | 100 cells |
|---|---|---|---|---|
| 10 | 836 | 1,118 | 1,547 | 1,759 |
| 20 | 1,137 | 1,973 | 4,123 | 6,826 |
| 26 | 1,206 | 2,242 | 4,928 | 8,414 |

The pool **saturates** — it can never exceed the total number of receptors buildable
from the gene set. At 10 genes and 100 cells, deduplication saves 86%.

### 2.3 Abundance (`cell_stoichiometry`) — selectable

Weights always sum to 1: a cell has a fixed number of receptor molecules regardless of
how many genes it expresses. More genes buys **variety**, not quantity.

- **`multinomial`** (default) — subunits produced in equal amounts, assembling at
  random. Receptor abundance is proportional to how many orderings collapse onto it.
  For g=2: weights 1, 5, 10, 10, 5, 1 out of 32 (Pascal's triangle). No free parameter.
- **`uniform`** — every receptor type in equal amount. Implies per-type assembly control.

**Consequence worth remembering:** under random assembly with g=5, homomers are
5 out of 5^5 = **0.16%** of what the cell builds. 99.84% is heteromeric. The claim
"heteromers dominate when multiple genes are expressed" is not an assumption — it
falls out of equal subunit production.

### 2.4 Cell activation (`cell_readout`)

The **drive** is the abundance-weighted fraction of the cell's receptors that are open:

```
drive = sum over receptors of (abundance × open probability)
```

Then `threshold` (default): `activity = sigmoid((drive − theta) / T_cell)`.

Also implemented: `noisy_or` (fires if any receptor opens; no threshold, saturates for
large repertoires) and `mean` (raw drive — diagnostic only, see §4.2 for why it is
wrong to feed a mean-of-probabilities to the entropy estimators).

### 2.5 The threshold

**One scalar shared by all cells, learnable, initialised at the calibrated median.**

Sharing it is deliberate. A per-cell threshold forces every cell to fire exactly 50% of
the time, which erases the heterogeneity we want: a narrowly-tuned cell *should* be
sparse, a broad one *should* fire often.

Calibration (before training, at the starting receptor temperature) sets theta to the
median of the pooled drive so the population starts at ~50% firing. Without it a
hardcoded theta=0.5 leaves every cell permanently silent — measured thresholds land
around 0.04–0.13, nowhere near 0.5.

Because the threshold is shared, some cells can start silent or saturated. A
diagnostic reports how many, at init.

---

## 3. Performance work (all verified exact)

Three optimisations, prompted by the pool being ~200× larger than the old receptor count.

### 3.1 Composition matrix

A receptor's threshold-concentration is the average of its 5 subunit energies. The old
code fetched 5 numbers per receptor into a table and averaged it away — 5× larger than
the answer, and retained by autograd for the backward pass.

Replaced by a counts matrix (rows = genes, columns = receptors, entries = copy number)
and one matmul. Same arithmetic, no 5-wide intermediate, runs on tensor cores.

### 3.2 Interface-model pocket deduplication — the big one

A pocket is fully determined by *(gene on the plus face, gene on the minus face)*. Not
by which receptor, not by ring position. So there are at most `n_genes²` distinct
pockets — 676 for 26 genes. The code was computing one per (receptor, slot):

| genes | cells | pool | pockets computed | actually distinct | redundancy |
|---|---|---|---|---|---|
| 10 | 20 | 7,540 | 37,700 | 100 | **377×** |
| 26 | 50 | 28,706 | 143,530 | 562 | **255×** |
| 26 | 100 | 53,456 | 267,280 | 644 | **415×** |

Fixed the same way: evaluate the distinct pockets once, contract to receptors with a
composition matrix over pocket-pairs.

### 3.3 bfloat16 autocast

Scoped to **training only**; evaluation stays float32 because bf16 moves the reported
entropy by ~0.0035 bits and individual activities by up to 0.05. Off by default
(`use_amp`).

### Measured results

| model | energy tensor elements | fwd+bwd speed |
|---|---|---|
| classic | 61,440 → 61,440 | **4.1× faster** |
| interface | 238,233,600 → **706,560** | **9.7× faster** |

Equivalence to the old path (both models):

```
max |p_gather − p_composition|      3.6e-07 (classic)   3.0e-07 (interface)
cell activity, chunked vs whole     3.0e-08              1.5e-08
cell activity, composition vs old   3.0e-08              2.2e-08
T_init calibration                  identical to 6 decimals
```

Float32 noise. Old gather path retained behind `use_composition=False`, so MWC and
external callers are untouched.

### Why none of this was worth doing before

Not an oversight — the cost structure inverted. The 5-wide transient at the batch each
array size can afford:

| receptors | batch | transient | |
|---|---|---|---|
| 20 | 32,557 | 0.10 GB | negligible |
| 45 | 21,705 | 0.15 GB | negligible |
| 4,900 | 14,560 | 10.63 GB | now it hurts |

At 20–45 receptors this was a tenth of a gigabyte. The bottleneck was entirely the
entropy estimator, which is where the engineering correctly went. Cell mode raises the
receptor count ~200× without changing the channel count.

---

## 4. Findings — including a serious one

### 4.1 The learnable threshold collapses to a degenerate solution

First test, 20 cells:

| | threshold | reported entropy |
|---|---|---|
| frozen | 0.055 → 0.055 | 16.30 → 16.94 bits |
| learnable | 0.042 → **0.000003** | 17.64 → **20.000 bits** |

20.000 out of a 20-cell maximum looks like a triumph. Opening the trained model:

```
receptor open probability → 0.0000      drive → 0.000000
cell activity = 0.5000 everywhere
variation across sniffs = 0.00002        <-- responds to nothing
99.99% of outputs within 0.001 of exactly 0.5
```

Every cell became a fair coin, independent of the stimulus. The existing KT self-test
in `bin_loss.py` predicts this exactly: *`[3b] identical A=0.5: H_kt = R`*. A perfectly
uninformative array reports the maximum possible entropy.

### 4.2 The deeper problem: entropy is only a proxy for information

Splitting the reported entropy into its parts (H = total output entropy,
h_cond = H(response | sniff) = noise, I = H − h_cond = actual information):

| | H | h_cond | **I** |
|---|---|---|---|
| frozen | 16.90 | 15.94 | **0.965 bits** |
| learnable | 20.00 | 20.00 | **0.000 bits** |

**Both runs were nearly worthless.** Maximising output entropy is a good proxy for
information only when the response is near-**deterministic** given the stimulus.

**Why this never bit the receptor picture:** the receptor temperature anneals until the
receptor is a hard switch, so h_cond ≈ 0 and entropy ≈ information. The two coincided.

**Why cells broke it:** cells add a *second* softness, `T_cell`, and its default was
chosen badly — an absolute 0.05 against a drive whose measured spread was 0.0955. The
sigmoid sat in its mushy middle and nearly all the "entropy" was that mush.

**Fix applied:** `cell_temperature` is now a *fraction of the calibrated drive spread*,
not an absolute value. Default 0.05 gives a sigmoid argument spread of ~20, i.e.
near-deterministic cells. Same logic as `compute_initial_temperature` already uses for
receptors; it simply had not been applied to the final temperature.

**Superseded 2026-09-01.** The relative-`T_cell` fix was necessary but not sufficient.
Two further degeneracies were found and fixed; see `doc/theory/07` §3b.7. In short: the
drive spans ~30 orders of magnitude, so (a) the threshold can land on the point mass of
zeros a sharp receptor produces, and (b) a standard deviation measures only the drive's
tail, not the region the threshold occupies. Fixes are `median_threshold` and
`drive_scale` (a MAD). End to end on a toy: 98.5% fair coins before, **0% after, with
49% of sniffs firing**.

> **STILL PENDING:** a trained run confirming information I > 0 on a realistic setup.
> The coin degeneracy is gone, which was the blocker; the positive confirmation has not
> been produced.

### 4.3 Why the threshold runs away

The threshold is one scalar receiving gradient from every (sniff, cell) pair —
measured at −138 from just 32 sniffs and 4 cells. The environment spreads its learning
over hundreds of parameters that act indirectly. So the threshold reaches its optimum
against a *random* environment long before the chemistry adapts, and the chemistry then
has to make the best of a threshold chosen for an array that no longer exists.

---

## 5. Remaining challenges

**Blocking, in order:**

1. **Verify the `T_cell` fix restores I > 0.** Nothing else matters until this is
   confirmed. (§4.2)
2. **Log `I = H − h_cond` as a measurement.** Without it, a run that has collapsed to
   noise reports a *perfect* score and the entropy curve looks great. Not implemented.
3. **Decide how to handle the threshold timescale gap.** Options, cheapest first:
   freeze theta for the first ~30% of training; give it its own smaller learning rate;
   re-calibrate periodically as the chemistry drifts; or reparameterise it as a
   *quantile* of the drive so it cannot run to a boundary.

**Not blocking:**

4. Confirm the bf16 speedup on the actual GPU — never measured, no GPU on the dev box.
5. Whether a shared threshold leaves too many cells silent at scale. Diagnostic exists;
   never checked beyond toy sizes.
6. Ten pre-existing failures in `unit_test/test_hierarchical_presence.py`. Verified
   present before any of this work (identical count on the stashed tree). Unrelated,
   but they are red.

---

## 6. The science question this was built to answer

The hypothesis, in the user's framing:

- A **homomer** is a broad, floppy detector.
- A **heteromer** is selective — opening needs all five subunits favourably engaged,
  so it behaves like an AND.
- A cell fires if *any* of its receptors fires — an OR across its repertoire.
- OR-ing a few broad detectors gives mush.
- OR-ing many narrow detectors gives a genuinely informative cell.
- Therefore multi-gene expression should **help** with heteromers and **hurt** with
  homomers only.

Two structural facts already support this without extra assumptions: heteromers
outnumber homomers combinatorially, and under random assembly they carry 99.84% of the
abundance in a 5-gene cell (§2.3).

### Simulations to run, and what each would show

**A. `multinomial` vs `uniform` stoichiometry.** Same cells, same genes, same
environment; the only difference is how much heteromer character each cell has
(99.84% vs 96%). *If the hypothesis holds, multinomial wins and the gap widens with
more genes per cell.* This is the cleanest single test, and the switch already exists.

**B. Sweep genes-per-cell (1 → 5) at fixed cell count.** The weighted-average readout
means a promiscuous cell averages over a wider repertoire and has a *narrower* dynamic
range. *Expect a genuine optimum* — too few genes and the cell sees little, too many
and it averages itself into blandness. A raw sum would have made more genes always win,
which is why the normalisation matters.

**C. Homomer-only cells vs full-repertoire cells at matched cell count.** Directly
isolates the AND-vs-OR claim by forcing the repertoire to homomers.

**D. Cells vs receptors at matched channel count.** Does 45 cells drawing on a
4,000-receptor pool beat 45 hand-picked receptors? This is the argument for the cell
picture being more than biological realism.

**E. Interface model vs standard model in cell mode.** The interface model is where the
pocket-between-subunits biology lives, and it is the model the AND argument actually
depends on. Now affordable (§3.2); previously ~400× too expensive to combine with cells.

**Constraint for all of the above:** cells capped at 5 expressed genes. That puts the
pool at ~1,200–8,400 (§2.2) and the memory at roughly 1.3 GB — comfortable, and it
means pool chunking (`cell_pool_chunk`) and a custom backward pass are both
unnecessary. They exist for a later "what if cells expressed 8 genes" question.

**Every one of these is invalid until item 1 in §5 is confirmed.** With cells too soft,
all of them would compare noise against noise.

---

## 7. Config quick reference

```python
RunConfig(
    # --- define the cells ---
    cell_gene_sets=((0,2), (0,1), (1,2,3), (4,)),   # explicit, OR:
    n_cells=20, cell_sampling_strategy="size_pmf",  # sampled
    cell_size_pmf=(0,0,0,0,1.0),                    # entry i = P(i+1 genes); this = always 5
    cell_sampling_seed=7,

    # --- the model choices that matter ---
    cell_stoichiometry="multinomial",   # or "uniform"      <-- experiment A
    cell_readout="threshold",           # or "noisy_or", "mean"
    cell_threshold="auto",              # calibrated median
    cell_threshold_learnable=True,      # shared scalar, in the optimizer
    cell_temperature=0.05,              # FRACTION of the calibrated drive spread

    # --- performance ---
    use_composition=True,               # exact; on by default
    use_amp=False,                      # bf16, training only; unverified on GPU
    cell_pool_chunk=None,               # not needed at 5 genes/cell
)
```

`cell_sampling_strategy="size_pmf"` is the hook for an arbitrarily complex expression
model: a distribution conditioned on how many genes are already expressed is exactly a
choice of `cell_size_pmf`.

Sampling happens once, in `SingleRunConfig.__post_init__`, which writes the resolved
gene sets and the derived pool into the config — so a run is exactly reproducible from
its `config.json`.
