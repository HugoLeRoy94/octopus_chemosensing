# 9. Cell Arrays — `src/cells.py`

Reference for the cell picture: what a cell is, how gene sets are sampled, how a cell
fires, and **why each choice was made rather than an equally plausible alternative**.

Cross-references: §01 (notation), §06.1 (the pool memory bottleneck), §07 stage 3b
(where this sits in the pipeline), §07 §3b.6–3b.7 (the degeneracies the readout avoids).

---

## 9.1 The shift, in one paragraph

In the receptor picture the array is a list of $R$ receptors, each 5 gene indices. In the
cell picture the array is a list of $C$ **cells**, each defined only by the set of subunit
genes $G_c$ it expresses — and a cell assembles **every receptor those genes allow**.

Both produce a $(\text{sniffs} \times \text{channels})$ table, and every loss and
measurement downstream consumes that table without asking what a channel is. The two
pictures therefore differ at exactly one point in the code:
`SimulationRunner._activity`. Nothing in `bin_loss.py`, `environment.py` or the
measurement registry needed to change.

**Why bother.** 45 cells is a far richer system than 45 receptors, because those cells
draw on a pool of several thousand distinct receptors underneath. The entropy estimator
only ever sees $C$, so the combinatorial richness is free on the estimator side — it is
paid for on the physics side instead (§9.6).

### Switching a config from receptors to cells

**Set exactly one of three fields. That is the whole switch.**

```python
RunConfig(
    ...,                      # everything else unchanged
    n_cells = 20,             # <- this line turns the array into cells
    cell_sampling_strategy = "size_pmf",
    cell_size_pmf = (0, 0, 0, 0, 1.0),
)
```

| to say | set | shape |
|---|---|---|
| "sample me $C$ cells" | `n_cells` | `int` |
| "these exact gene sets" | `cell_gene_sets` | tuple of gene tuples |
| "these exact repertoires" | `cell_receptors` | tuple of (tuple of receptor tuples) |

Precedence is `cell_receptors` > `cell_gene_sets` > `n_cells`. Any of them makes
`SingleRunConfig.is_cell_mode()` true, and then:

* `receptor_indices` is **derived** (the deduplicated pool) — do not set it yourself;
* `n_receptors` and `receptor_sampling_*` are **ignored**;
* the entropy is computed over the $C$ cells, while the physics still runs over the
  $R_\text{pool}$ receptors underneath (§9.6);
* every `cell_*` parameter below becomes live. All have working defaults, so a config
  that sets only `n_cells` is already valid.

**Going back** is just as simple: remove all three and the array is receptors again.
Nothing else in the config is cell-specific.

Note that the values are **tuples, not lists**. `RunConfig` reads any list-valued field
as a sweep axis, so a list here would be silently mistaken for one.

### The cell parameters at a glance

| parameter | default | what it does | § |
|---|---|---|---|
| `n_cells` | `None` | how many cells to sample | 9.4 |
| `cell_gene_sets` | `None` | exact gene set per cell | 9.4 |
| `cell_receptors` | `None` | exact repertoire per cell, bypassing gene expansion | 9.4 |
| `cell_sampling_strategy` | `"bernoulli"` | `"bernoulli"` (per-gene probability) or `"size_pmf"` (per-cell gene count) | 9.4 |
| `cell_gene_probs` | `None` | bernoulli: $P(\text{gene } u \text{ expressed})$; defaults to $2/n_\text{genes}$ | 9.4 |
| `cell_size_pmf` | `None` | size_pmf: entry $i$ = $P(i{+}1$ genes$)$ | 9.4 |
| `cell_max_genes` | `None` | bernoulli: reject cells above this — the main lever on $R_\text{pool}$ | 9.4 |
| `cell_sampling_seed` | `None` | reproducible gene sets | 9.4 |
| `cell_stoichiometry` | `"multinomial"` | receptor abundance within a cell: random assembly, or `"uniform"` | 9.5 |
| `cell_readout` | `"threshold"` | how the drive becomes an activity: also `"noisy_or"`, `"mean"` | 9.7 |
| `cell_threshold` | `"auto"` | `"auto"` pins $\theta$ to the median drive; a float fixes it | 9.8 |
| `cell_threshold_learnable` | `False` | let Adam move $\theta$ — an ablation, see the warning in 9.8 | 9.8 |
| `cell_temperature` | `0.01` | final $T_{cell}$, as a FRACTION of the live drive spread — never absolute | 9.8, 9.9 |
| `cell_initial_temperature` | `"auto"` | starting $T_{cell}$, same units | 9.9 |
| `cell_phase_split` | `0.5` | fraction of epochs in phase 1 (receptor sharpens, cell held soft) | 9.9 |
| `cell_recalibrate_every` | `25` | epochs between re-pinning $\theta$ to the moving drive | 9.8 |
| `cell_n_molecules` | `1e4` | receptor MOLECULES per cell; floors $\theta$ at $1/N$ (one open channel) | 9.8.1 |
| `cell_pool_chunk` | `None` | receptors per chunk in the fused physics+readout pass | 9.6 |

The two that most often need attention: **`cell_max_genes`** (or a `cell_size_pmf` that
caps the size) because $R_\text{pool}$ is superlinear in genes-per-cell, and
**`cell_temperature`**, which is a *ratio* — an absolute value there reintroduces the
degeneracy of §9.8.

---

## 9.2 Gene set → repertoire (`expand_gene_set`)

A cell expressing $G_c = \{0, 2\}$ with $k_{sub} = 5$ assembles:

```
[0,0,0,0,0]  [0,0,0,0,2]  [0,0,0,2,2]  [0,0,2,2,2]  [0,2,2,2,2]  [2,2,2,2,2]
```

* **Standard model** — unordered multisets: $\binom{g + k_{sub} - 1}{k_{sub}}$ receptors.
* **Interface model** — canonical cyclic ring arrangements, rotations identified and
  reflections DISTINCT (the $\pm$ face asymmetry breaks mirror symmetry). By Burnside:
  $\tfrac{1}{k_{sub}}\sum_{d \mid k_{sub}} \varphi(d)\, g^{k_{sub}/d}$. Reuses
  `geometry._canonical_rotation`, so the convention matches `build_heteromer_array`.

| $g$ | standard | interface |
|---|---|---|
| 1 | 1 | 1 |
| 2 | 6 | 8 |
| 3 | 21 | 51 |
| 5 | 126 | 629 |
| 8 | 792 | 6560 |

The list is returned **sorted**, so the pool order is deterministic across runs — which
is what lets `SingleRunConfig.__post_init__` store the derived pool in `config.json` and
`SimulationRunner._initialize` rebuild $W$ against it with an assertion rather than a hope.

---

## 9.3 The pool (`CellArray`)

`CellArray` expands every cell, then takes the **deduplicated union** as
`receptor_indices` ($R_\text{pool}$) and records abundances in $W$ of shape
$(C, R_\text{pool})$, with $W_{cr} = 0$ where cell $c$ cannot assemble $r$.

Deduplication is not a micro-optimisation. A receptor's open probability is a property of
the receptor, not of the cell holding it, so a receptor shared by several cells must be
simulated **once**. Cells drawn from overlapping gene sets pool very cheaply, and the pool
**saturates** — it can never exceed the full combinatorial pool for the gene set:

| available genes | 10 cells | 20 cells | 50 cells | 100 cells |
|---|---|---|---|---|
| 10 | 836 | 1,118 | 1,547 | 1,759 |
| 20 | 1,137 | 1,973 | 4,123 | 6,826 |
| 26 | 1,206 | 2,242 | 4,928 | 8,414 |

(every cell expressing exactly 5 genes). At 10 genes and 100 cells, deduplication saves
86% of the work.

---

## 9.4 How cells are sampled

Two samplers, both seeded and reproducible. Sampling happens **once**, in
`SingleRunConfig.__post_init__`, which writes the resolved gene sets back into
`cell_gene_sets` and the derived pool into `receptor_indices` — so a run is exactly
reproducible from its own `config.json`, with no RNG replay needed.

### `bernoulli` — `sample_gene_sets_bernoulli`

Gene $u$ is expressed independently with probability `cell_gene_probs[u]` (default
$2/n_\text{genes}$, i.e. two genes per cell on average). Cells expressing nothing are
rejected and redrawn. `cell_max_genes` rejects cells above a size cap.

Use it when you want expression to be a property of the **gene** — some genes common,
some rare — for instance to model a promoter-strength distribution. The number of genes
per cell is then Poisson-binomial, not something you control directly.

### `size_pmf` — `sample_gene_sets_by_size`

Two stages: draw the **number** of expressed genes from `cell_size_pmf` (entry $i$ is the
unnormalised probability of $i+1$ genes), then draw that many genes uniformly without
replacement.

Use it when you want to control the promiscuity distribution directly — which is the
usual experimental question ("what happens as cells express more genes?"). Setting
`cell_size_pmf=(0,0,0,0,1.0)` gives every cell exactly 5 genes.

**This is the hook for an arbitrarily complex expression model.** A distribution over
"how many genes does this cell express", conditioned on anything you like, is exactly a
choice of `cell_size_pmf`. It was designed this way so the sampler never has to be
rewritten to ask a new question — including the conditional-on-already-expressed models
that motivated the cell picture in the first place.

### Why a size cap matters

$R_\text{pool}$ is the binding cost of the whole pipeline (§9.6), and the repertoire
table in §9.2 is superlinear in $g$. Capping genes per cell is therefore the
highest-leverage knob available, well above chunking or precision tricks.

### Explicit repertoires (`cell_receptors`)

Both samplers choose GENE SETS, and the repertoire follows from them. Sometimes the
repertoire itself is what you want to control, and gene sets cannot express it:

> A cell containing **exactly one heteromer** is not a gene set. Expressing the genes of
> `[0,0,0,1,1]` necessarily also produces every other combination of genes 0 and 1. One
> gene per cell does give exactly one receptor — but only ever the homomer.

`cell_receptors` states the receptor list per cell outright, bypassing the expansion.
Weights are uniform over the receptors listed (you have specified the repertoire, so
`cell_stoichiometry` does not apply), and `cell_gene_sets` is derived from them for
reporting. Receptors are canonicalised on the way in — sorted as a multiset, or
minimum-rotation under the interface model — so an explicitly listed receptor is the
same object the gene-set path would have produced.

One receptor per cell is `tuple((r,) for r in RECEPTORS)`. With the list **sorted**, the
pool is that same list and $W$ is exactly the identity — which is what
`tasks/cells/equivalence` uses to check that cell mode reduces to the receptor model.

---

## 9.5 Abundance (`repertoire_weights`, `cell_stoichiometry`)

Rows of $W$ **always sum to 1**. This is a modelling commitment, not bookkeeping: a cell
has a fixed number of receptor molecules in its membrane regardless of how many genes it
expresses. Expressing more genes buys **variety**, not quantity.

The consequence is the interesting one. Because the drive is an average, a promiscuous
cell averages over a wider repertoire, and averages concentrate — so a promiscuous cell
has a *narrower* dynamic range. There should therefore be a genuine optimum in
genes-per-cell: too few and the cell sees little, too many and it averages itself into
blandness. **A raw sum would have made more genes monotonically better**, which would
have been an artefact of the readout rather than biology.

### `multinomial` (default) — random assembly

All subunits are produced in equal amounts and assemble independently into the $k_{sub}$
slots, so each of the $g^{k_{sub}}$ ordered words is equally likely and a receptor's
abundance is proportional to how many words collapse onto it:

$$w_r \propto \frac{k_{sub}!}{\prod_i n_i!} \quad\text{(standard)}, \qquad
  w_r \propto \#\{\text{distinct rotations}\} \quad\text{(interface)}$$

For $g = 2$: $1, 5, 10, 10, 5, 1$ out of 32 — Pascal's triangle. **No free parameter**;
it follows from "equal subunit production" alone, which is what makes it the default.

### `uniform`

Every receptor type present in equal amount, $w_r = 1/|\text{repertoire}|$. Simpler, but
it implies the cell actively regulates assembly per receptor TYPE rather than just
producing subunits.

### The consequence worth remembering

Under random assembly with $g = 5$, homomers are $5$ of $5^5 = 3125$ ordered words:

> **0.16% of what the cell builds. 99.84% is heteromeric.**

So "heteromers dominate once several genes are expressed" is not an assumption bolted on
— it falls out of equal subunit production. This matters for the AND-vs-OR argument: a
heteromer is selective (all five subunits must be favourably engaged), a homomer is
broad, and a cell ORs over its repertoire. `multinomial` vs `uniform` (homomers at 4%) is
therefore a clean two-line experiment on exactly that question.

---

## 9.6 Memory: the pool, not the cells

The physics tensors are $(B, L, R_\text{pool}[, k_{sub}])$, and $R_\text{pool}$ grows far
faster than $C$. Full treatment in §06.1; the three levers, in the order to reach for them:

1. **Deduplicate** — automatic (§9.3).
2. **Composition binding** — `use_composition=True`, default. Removes the $k_{sub}$ factor
   in the standard model and the ~400× pocket redundancy in the interface model (§06.1).
3. **Chunk the pool** — `cell_pool_chunk`. Every readout reduces over $r$ through a term
   LINEAR in the per-receptor contribution, so `cell_activity` walks the pool in slices
   and sums the $(B, C)$ accumulator: peak forward tensor $O(B \cdot L \cdot \text{chunk})$.
   Pair with `recompute_backward=True` or the saving holds only under `no_grad`.
4. **Cap genes per cell** — `cell_max_genes` (§9.4).

At 5 genes per cell, $R_\text{pool}$ lands at ~1,200–8,400 and the transient at ~1.3 GB.
**Levers 3 and 4 are then unnecessary** — they exist for the "what if cells expressed 8
genes" question.

`resolve_batch_sizes` receives `n_receptors=C` and `n_physics_receptors=R_pool`, so the
estimator caps size on the cells while only the physics cap sees the pool.

---

## 9.7 Firing (`CellReadout`) — and the choices behind it

The **drive** is the abundance-weighted open fraction:
$S_{bc} = \sum_r W_{cr}\, p_{br}$ — plainly, *what fraction of this cell's receptors are
open*, hence what fraction of its maximum current is flowing.

`threshold` (default): $A_{bc} = \sigma\big((S_{bc} - \theta)/T_{cell}\big)$.

### Why a threshold at all, when `mean` is simpler

`mean` feeds $S$ straight to the estimators. That is not wrong so much as a **different
model**: it says the cell fires stochastically at a rate linear in its open fraction. The
estimators then read the coin correctly — but the model now contains noise, and this
project treats the output as deterministic.

The defensible justification for the threshold is **not** biophysical sharpness (real f–I
curves are often threshold-*linear*, not step-like). It is that a **binary code is a
decision, and a decision has a threshold**. What biology does supply is the threshold
itself: below rheobase a neuron emits nothing, however long you wait. If the output were
genuinely graded, the right code would be multi-level and these estimators would be the
wrong tool.

### Why not `noisy_or`

$A_{bc} = 1 - \exp(\sum_r k_{sub} W_{cr}\ln(1 - p_{br}))$. Two objections, both fatal:
a cell whose repertoire is 99.84% heteromeric would fire off its few broad, unselective
homomers; and with a large repertoire it **saturates toward 1** — nearly binary, but
stuck ON.

### Why `mean` is unusable here specifically

$S$ is an average over the repertoire, and averages concentrate, so it never approaches
0 or 1. Every activity lands mid-range, every cell reads as a coin, and the reported
entropy is almost entirely noise. Kept as a diagnostic only.

---

## 9.8 Calibration — the two rules that keep the sigmoid honest

`calibrate_cell_readout` sets $\theta$ and $T_{cell}$ from the empirical drive, and is
re-run every `cell_recalibrate_every` epochs because the drive distribution moves as the
chemistry trains. Both helpers exist to dodge a specific degeneracy; full treatment with
measured numbers in §07 §3b.7.

The root cause of both: the drive spans a huge dynamic range. Measured at the final
receptor sharpness, **26% of pairs exactly zero, median $\approx 10^{-30}$, maximum
$\approx 1$ — thirty orders of magnitude.**

**`median_threshold`** — $\theta$ is the midpoint between the median and the next distinct
value above it, **floored at $1/N$**. For a continuous drive the midpoint rule is the
textbook median. For a drive with a point mass at zero (the sparse-code case) it steps
$\theta$ strictly ABOVE the mass, because a sharp sigmoid centred on a point mass returns
exactly 0.5 for every member of it — turning the whole mass into fair coins. The floor is
what keeps $\theta$ inside physics at all, and it is important enough to have its own
section (§9.8.1).

**`drive_scale`** — the spread is a **MAD, not a standard deviation**, because $T_{cell}$
is a fraction of it and the spread must be measured WHERE THE THRESHOLD SITS. A standard
deviation of that drive is set entirely by the few strongly-firing sniffs: it leaves 91%
of activities as coins. IQR leaves 67%. The MAD leaves **0%**, at 50% firing.

Corollary: **no absolute floor may be put on $T_{cell}$.** A guard like `max(T, 1e-6)`
looks harmless and reinstates the entire failure when the drive lives at $10^{-30}$. The
floor in `run.py` is $10^{-45}$, just off zero.

### 9.8.1 $N$: the receptor copy number, and why $\theta$ is floored at $1/N$

**What $N$ is.** The number of receptor **molecules** a single cell carries in its
membrane — its copy number.

**What $N$ is not**, since both are easy to reach for and neither is this quantity:

* *not* the number of receptor **types** the cell can assemble. That is its repertoire,
  $\binom{g+k_{sub}-1}{k_{sub}}$, which is a count of distinct species, not of molecules.
  A cell with a 126-type repertoire still carries $N$ molecules in total, distributed
  across those types by $W$ (§9.5).
* *not* the number of **cells** in a sucker or an animal. That is the array size $C$.
  The drive is an average over the receptors *inside one cell*, so nothing about how
  many cells sit alongside it enters the quantity.

**Why it matters.** The drive $S$ is a *fraction* of the cell's receptors that are open.
With $N$ molecules, the number open is an integer count $n = N S$, so:

> $1/N$ is the resolution of the drive. Below one open channel there is no state for the
> cell to be in.

**What goes wrong without it.** The median is a pure **rank** statistic: it splits the
sample in half and knows nothing about what the numbers mean. In a sparse code most
drives are numerically tiny — and sharpening the receptor makes this *worse*, not better.
For an OFF receptor $p = \sigma(\ln\text{sum}/T) \approx e^{\ln\text{sum}/T}$, so shrinking
$T$ divides $\ln\text{sum}$ by $T$ and spreads the OFF values over **more** decades.
Measured, sweeping the receptor temperature $1.0 \to 0.003$:

| receptor $T$ | 1.0 | 0.3 | 0.1 | 0.03 | 0.01 | 0.003 |
|---|---|---|---|---|---|---|
| median-pinned $\theta$ | 1e-01 | 1e-03 | 1e-09 | 1e-29 | 1e-87 | 1e-290 |

At $N = 10^4$ that last value is $10^{-286}$ open channels. It is not a threshold; it is a
rank statistic that has lost contact with its units. And it is splitting a cloud of values
that do not physically exist — **at $N = 10^5$, 83% of drives sit below $1/N$** — instead
of separating OFF from ON.

**What the floor restores.** Exactly the behaviour a cell should have: no receptors open
$\to$ OFF, receptors open $\to$ ON. Measured at $N = 10^4$ against a 15.8% underlying
receptor firing rate:

| | $\theta$ | agrees with the receptor code |
|---|---|---|
| unfloored | 1.5e-39 | 86.9% |
| **floored at $1/N$** | 1.0e-04 | **98.9%** |

with 17.0% firing and 0% fair coins.

**It fixes $T_{cell}$ for free.** No second floor is needed: `drive_scale` is a MAD
measured *around* $\theta$, so once most drives sit far below $\theta$ the typical
deviation is $\approx \theta$ itself. The MAD lands on $1/N$ and $T_{cell}$ inherits the
same physical scale — verified across $N = 10^3 \ldots 10^6$, where $T_{cell} \cdot N$ stayed
at 0.010 throughout. $T_{cell}$ therefore remains a *ratio* that anneals, and the
two-phase schedule of §9.9 is untouched.

**Why $n = 1$ and not "channels needed to reach rheobase".** A threshold set from the
firing current, $\theta = n_\text{thresh}/N$, is tempting and defensible — but the result
is insensitive to it: $n_\text{thresh} = 1$ vs $100$ moved the firing rate 17.0% $\to$
16.4%. A second constant that changes nothing earns nothing, so the floor is simply one
open channel.

**Two properties to keep in mind.**

* The floor **only bites in the sparse regime**. When the code is dense enough that the
  median sits above $1/N$, $\theta$ is the median and behaviour is unchanged. The
  calibration reports `theta_floored` so you can tell which regime a run is in.
* Where it bites, "the population fires ~50%" becomes **conditional**: a sparse world
  yields a sparse code (17%, not 50%). That is the honest answer rather than a
  regression, but it is a real behavioural change from an unfloored median.

`cell_n_molecules` defaults to $10^4$ — an order-of-magnitude placeholder. Replace it with
a measured copy number; the results above suggest the answer is not sensitive to getting
it exactly right, only to getting it roughly right. Setting it to `None` disables the
floor and restores the pre-floor behaviour, degeneracy included.

### Why $\theta$ is one shared scalar, and not learned

**Shared, not per-cell.** A per-cell threshold forces every cell to fire exactly 50% of
the time, erasing the heterogeneity the cell picture exists to study — a narrowly-tuned
cell *should* be sparse and a broad one *should* fire often. The cost is that a cell whose
drive never reaches $\theta$ is silent, so calibration reports `n_silent` / `n_saturated`.

**Pinned, not learned** (`cell_threshold_learnable=False` by default). $\theta$ is a
single scalar receiving gradient from every $(b, c)$ pair — measured at $-138$ from 32
sniffs and 4 cells — while the environment spreads its learning over hundreds of
parameters acting indirectly. It therefore reaches its optimum against a *random*
environment long before the chemistry adapts. Measured consequence: it collapsed to
$3\times10^{-6}$, parked every cell at activity 0.5, and reported a **perfect 20.000 bits
out of 20 while carrying exactly zero information**. Pinning it to the median costs no
free parameter and is near-optimal anyway (exactly optimal for a single cell; the shortfall
for the joint code is second-order against what the chemistry can do).

The learnable flag is kept for that ablation, and is now safe to use *provided* the
sharpness schedule of §9.9 runs — a deterministic cell has no 0.5 to park at.

---

## 9.9 The two-phase sharpness schedule

Annealing the receptor and the cell together makes the cell's operating point chase a
drive distribution that is still moving underneath it. In cell mode the schedule splits
at `cell_phase_split` (default 0.5):

| | receptor $T$ | cell $T_{cell}$ |
|---|---|---|
| **Phase 1** | anneals to final, reaching it on the LAST epoch of the phase | held soft at $1\times$ the live spread |
| **Phase 2** | held | anneals to `cell_temperature` $\times$ spread |

Phase 1 keeps gradients live everywhere while the chemistry arranges the drive around the
threshold; phase 2 hardens the readout. `cell_temperature` defaults to 0.01, leaving
~0.8% of $(b,c)$ pairs inside the transition band.

**Deterministic cells are not a nicety — they are the condition under which the entropy
objective is a valid proxy for information** (§07 §3b.6). This never bit the receptor
picture because the receptor temperature already anneals until $p \in \{0,1\}$; cells add
a *second* softness, and if $T_{cell}$ is comparable to the drive spread the reported
entropy is almost all noise.

The known risk is phase-2 gradient starvation: once the cell is sharp, only samples inside
the thin transition band pass gradient. Pinning $\theta$ to the median mitigates it — for
a unimodal drive the median is the densest point, so the band sits where the most samples
are. The two choices reinforce each other.

---

## 9.10 Config quick reference

```python
RunConfig(
    # --- define the cells (explicit repertoires, explicit gene sets, or sampled) ---
    cell_receptors=tuple((r,) for r in RECEPTORS),  # exact repertoire per cell, OR:
    cell_gene_sets=((0,2), (0,1), (1,2,3), (4,)),   # explicit gene sets, OR:
    n_cells=20, cell_sampling_strategy="size_pmf",  # sampled
    cell_size_pmf=(0,0,0,0,1.0),                    # entry i = P(i+1 genes)
    cell_gene_probs=None,                           # bernoulli strategy instead
    cell_max_genes=None,                            # caps R_pool; see 9.4
    cell_sampling_seed=7,

    # --- model choices ---
    cell_stoichiometry="multinomial",   # or "uniform"   <- the AND/OR experiment
    cell_readout="threshold",           # or "noisy_or", "mean"
    cell_threshold="auto",              # pinned to the median drive
    cell_threshold_learnable=False,     # True only as an ablation; see 9.8
    cell_temperature=0.01,              # FRACTION of the live drive spread
    cell_phase_split=0.5,
    cell_recalibrate_every=25,
    cell_n_molecules=1e4,               # receptor copy number; floors theta at 1/N

    # --- performance ---
    use_composition=True,               # exact; default on
    cell_pool_chunk=None,               # unnecessary at <= 5 genes/cell
)
```

Setting `cell_receptors`, `cell_gene_sets` **or** `n_cells` switches the array to cell
mode (`cell_receptors` wins over `cell_gene_sets` if both are given);
`n_receptors` / `receptor_sampling_*` are then ignored and `receptor_indices` is derived.
