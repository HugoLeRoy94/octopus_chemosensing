# Cell model: scientific and computational review

Review date: 2026-09-10. This is an assessment of the current checkout, not a record of fixes.

Follow-up implementation: explicit-repertoire reconstruction, the cell-temperature
endpoint, sparse-mixture overflow, and the presence tests' return-type mismatch have
subsequently been fixed. Findings below describe the reviewed pre-fix state. The
thresholded readout and qualitative equivalence/convergence experiments are retained;
the tasks now train with `entropy='kt_mi'` and report an MI bracket, response
conditional entropy, and final stochastic-output counting. The implementation is
documented in theory §04, §05 and §07. Historical findings below describe entropy
training. Finite-copy channel noise remains a proposed extension (theory §09.10).

The approach is useful for testing whether combinatorial assembly increases information at a fixed gene and cell budget. The separation between a shared receptor pool and cell readout is particularly useful. However, the current equivalence task cannot initialize, several claimed bounds are not bounds, and the physical meaning of the two sigmoid temperatures needs to be fixed before interpreting cell information biologically. I would repair the validation layer before launching a broad simulation campaign.

Scope: detailed inspection of the active `opt_bin_resp` cell/physics/environment/configuration/training/entropy/measurement paths, both cell tasks, theory documents, and project context in the expression analyses and progress report. I also inspected the March expression dataset. This is not a line-by-line certification of every historical notebook, fitting script, manuscript, or saved simulation in the repository. No simulation source files were changed. No A100 benchmark or 5,000-epoch training was run. CPU checks used the existing `/home/hugo/miniconda3/envs/work/bin/python` environment, PyTorch 2.11.0+cpu.

## 1. What limits the bits?

Use separate symbols for cells C, receptor types R_pool, molecules per cell N, ligand vocabulary L, and independent evaluation sniffs B. They control different ceilings.

### Hard output and stimulus bounds

For a binary array, always:

$$
I(Y;X)\le H(Y)\le C.
$$

For a deterministic cell readout, duplicate rows of the abundance matrix W with the same threshold produce identical outputs. If there are U_c distinct rows, replace C by U_c. Adding copies of the same noiseless cell adds no information. Independent molecular noise could make repeat cells useful, but that is a different model.

With one of L equiprobable ligands, exactly fixed concentration, no observation noise, and a deterministic readout:

$$
H(Y)=I(Y;\ell)\le\min(C,\log_2 L).
$$

Eight ligands therefore give at most 3 bits, regardless of the number of receptor types or latent dimensions. This is an upper bound; attaining it requires that the allowed architecture can separate all eight ligands.

With one ligand at a time and continuously varying concentration, the current positive-weight readout is monotone in concentration for every cell. Each cell crosses its fixed threshold at most once. Consequently each ligand generates at most C+1 hard codes:

$$
H(Y)\le\min\{C,\log_2[L(C+1)]\}.
$$

This bound assumes fixed ligand coordinates, no observation noise, and frozen thresholds during evaluation. Examples: L=8, C=6 gives a loose ceiling of 5.807 bits; L=100, C=100 gives 13.302 bits. Independent concentrations in mixtures do not obey this one-dimensional counting argument. A common concentration multiplier for a fixed mixture composition does.

For fixed-concentration mixtures drawn from a finite set of masks M, use

$$
I(Y;M)\le H(M),\qquad H(Y)\le\min(C,H(M))
$$

in the deterministic case. If exactly s of L ligands are selected uniformly, H(M)=log2 choose(L,s). With variable mixture size S and uniform selection conditional on S, H(M)=H(S)+E[log2 choose(L,S)].

The actual hierarchical sampler admits an analogous exact composition-entropy calculation. Let A be the number of selected source blocks among K blocks; block k contains m_k ligands and has a zero-truncated, upper-truncated Poisson count N_k. Define

$$
h_k=H(N_k)+\mathbb E[\log_2 {m_k\choose N_k}].
$$

Because selected blocks are uniform and each selected block is nonempty, the mask identifies its active blocks, so the implemented selection law before sparse truncation gives

$$
H(M)=H(A)+\mathbb E[\log_2 {K\choose A}]
+\frac{\mathbb E[A]}{K}\sum_{k=1}^K h_k.
$$

This is a useful inexpensive preflight number. `mu_sources` and `mu_ligands_per_source` are Poisson parameters, not the post-truncation mean counts. In particular, their small-parameter limits approach one selected item, not zero.

### Bounds that should not guide the campaign

`doc/theory/08_environmental_entropy_limits.md` treats a sum of discrete composition entropy and continuous concentration entropy as a hard output-entropy ceiling. Differential entropy cannot bound discrete output entropy this way: changing concentration units changes differential entropy without changing the code. Your `doc/prompts/topic3_cover_and_entropy_bounds.md` already identifies this correctly. Use a discrete task variable, an explicit concentration resolution/noise model, or support-counting bounds.

The proposed R_eff approximately 2 times mean mixture size is a coverage heuristic, not a general information bound. For example, a categorical single-ligand input can be mapped to many binary identification bits by overlapping tuning sets. Disjoint coverage is an extra restriction. The motivation from balanced, weakly correlated responses is useful, but does not establish that hard ceiling. [Zwicker, Murugan and Brenner](https://arxiv.org/abs/1602.02974) provide the underlying optimal-array framework.

Likewise, the hyperplane arrangement argument in the Cover discussion applies to an appropriate quadratic/lifted model. It cannot be carried unchanged to the saturating Gaussian kernel followed by a weighted sum of receptor sigmoids. A loose upper bound being well above a result does not prove optimization failure: attainability under shared genes, weights, and thresholds must also be demonstrated. In the hyperplane count itself, a dimension near C/2 gives nearly C bits, not all 2^C regions; exact saturation and near-saturation differ.

Family identity also needs a precise definition. A categorical family F has H(F)<=log2 n_families. A mixture's family-presence vector is a different variable, potentially carrying up to n_families bits. Current family measurements largely condition on individual binary indicators.

## 2. Estimator limits are likely to bind before cell count

KT is a valid way to bound the entropy of the empirical mixture of product Bernoulli distributions. The Bhattacharyya and KL bounds concern that mixture, not a confidence interval for population entropy. A tight bracket can coexist with severe undersampling. The original work and its corrected version distinguish these mixture bounds from statistical estimation. [Kolchinsky and Tracey](https://arxiv.org/abs/1706.02419).

For activities a_bc interpreted as conditionally independent Bernoulli probabilities, define

$$
H_{cond}=\frac1B\sum_b\sum_c h_2(a_{bc}).
$$

The KT formula in `bin_loss.py:171` has the form H_KT=H_cond+an overlap term. Its overlap term is at most log2 B because the self-pair is retained. Thus it is the signal term, not the total soft entropy, that has the sample ceiling. Setting all activities to 0.5 yields C reported bits and zero information about the stimulus. The existing self-test reproduces this.

Log H_cond and H_KT-H_cond separately on every run. Subtracting H_cond gives a lower bound on information about the empirical component index under the product-Bernoulli model. If docking noise varies across sniffs, this index includes the noise realization: it still does not directly measure information about the clean chemical stimulus. For that, repeat the same clean stimulus with independent observation noise and estimate the appropriate conditional entropy.

| Independent samples B | Maximum empirical hard-code entropy / KT signal ceiling |
|---:|---:|
| 2,048 | 11 bits |
| 4,096 | 12 bits |
| 16,384 | 14 bits |
| 65,536 | 16 bits |
| 1,048,576 | 20 bits |

These are ceilings, not accuracy guarantees. For approximately K equally likely hard codes in the well-sampled regime, the leading plug-in bias is approximately (K-1)/(2 B ln 2). Ten samples per code correspond to about 0.072 bits of leading bias; 100 samples per code to about 0.0072 bits. This approximation is not reliable in the undersampled regime or a poorly sampled long tail. Miller–Madow cannot reconstruct a large unseen alphabet: its implemented correction is at most approximately 0.721 bits because K_hat<=B.

KT costs O(B² C) arithmetic. Increasing B from 4,096 to 65,536 costs about 256 times as many pair-cell operations at fixed C. Increasing to 1,048,576 costs 65,536 times as many. Tiling changes memory and launch overhead, not that operation count; the O(B² C/chunk) time expression in several docstrings is misleading.

For sharp outputs, use streamed hard-code counting on much larger evaluation samples, then check stability under sample doubling and report unseen-code diagnostics. Keep KT for optimization and modest-size soft-output checks. Do not infer 40–100 resolved information bits from a few thousand samples just because there are 40–100 cells. Conversely, a 100-cell array encoding only eight deterministic categories is easy to evaluate accurately: required sampling follows the distribution of codes, not simply 2^C.

The `collision` estimator needs a separate correction to its documentation. Its diagonal is removed, so it does NOT have the advertised hard log2(chunk) ceiling. With all 16 distinct four-bit codes, the current function returned 134.784 bits, while KT returned 3.989 bits. No off-diagonal hard collisions exist; numerical smoothing makes the collision estimate tiny. True population H2<=H is valid, but this finite-sample estimate is not a certified lower bound on population Shannon entropy. Expected collision count, roughly number_of_tested_pairs times collision_probability, governs reliability.

## 3. Receptor repertoire limits and what is actually simulated

For pentamers, unordered compositions number choose(g+4,5); cyclic rings with rotations identified and reflections retained number (g^5+4g)/5.

| Expressed genes g | Standard repertoire | Interface ring repertoire |
|---:|---:|---:|
| 1 | 1 | 1 |
| 2 | 6 | 8 |
| 3 | 21 | 51 |
| 5 | 126 | 629 |
| 8 | 792 | 6,560 |
| 10 | 2,002 | 20,008 |
| 22 | 65,780 | 1,030,744 |
| 26 | 142,506 | 2,376,296 |

These are structural counts. They are not counts of independent tunable sensors or distinguishable response functions.

In the standard model, each receptor's energy vector is a convex combination of gene energy vectors. The energy-table rank is at most the number of genes. This restricts realizable tuning; it does not imply information is bounded by log2(number of genes), or by the energy-table rank in bits.

In the interface model there are at most n_genes² distinct directed pockets. Moreover, two different cyclic rings with the same directed-pocket counts have identical response functions in the current additive model. A direct check found 51 structural three-gene pentamers but only 48 distinct pocket-count vectors. Functional deduplication could therefore reduce cost further, provided abundance weights are summed when merging receptor types.

The active runner always constructs `BinaryReceptor` (`run.py:468`). It supports homomers, unordered heteromers, and directed-interface rings under the threshold approximation. The presence of `MWCReceptor` in `physics.py` does not make full MWC a selectable, validated end-to-end cell model: it expects open/closed energies, does not support composition binding, and its current formula lacks an explicit adjustable unliganded gating-energy term. Antagonism, inhibitory currents, unequal per-type conductance, assembly preferences, finite molecule draws, temporal firing, and spatially different inputs across cells are not implemented in the main cell path.

The geometric-mean EC50 is a model limit, not an identity for arbitrary MWC parameters. In the reduced MWC equation in §02, deriving it additionally requires the closed-state binding factor to be approximately one near half activation. The additive mixture drive is also a separate approximation; it need not equal the product of subunit binding polynomials in a general mixed-ligand MWC model.

### An important boundary on the AND/OR interpretation

An arithmetic mean of subunit energies is not a strict logical AND. A favorable contribution can compensate an unfavorable one. Nor does a larger repertoire guarantee narrower tuning or larger information.

There is a stronger result for the standard model. Let d_u=sum_l c_l exp(-E_ul) be the homomer drive and alpha_u the subunit fractions of a heteromer. Its drive is

$$
d_r=\sum_l c_l\prod_u e^{-\alpha_u E_{ul}}
\le\sum_u\alpha_u d_u\le\max_u d_u,
$$

by weighted arithmetic–geometric mean, since alpha_u>=0 and sum alpha_u=1. Therefore, in the hard receptor limit, if every homomer is OFF, every heteromer is OFF. Because the full repertoire includes the homomers, a literal OR over the complete standard repertoire has exactly the same ON region as an OR over its homomers. Heteromer benefits must arise through abundance weighting, graded/thresholded pooling, or a different pocket model—not merely adding extra receptors to that literal OR. Interface pockets can introduce new affinities, so this particular proof does not transfer to them.

The implemented thresholded weighted drive is not a literal OR. Its information comparison remains interesting, but should be described accordingly. Also, the implemented `noisy_or` uses exponents k_sub W, whose sum is always k_sub. For identical receptor probabilities p it gives 1-(1-p)^k_sub, independent of repertoire size. Claims that it necessarily saturates just because the repertoire grows are therefore unsupported by its formula.

## 4. Finite molecule count is a substantive model boundary

`cell_n_molecules=10000` currently only floors calibration of theta at 1/N. It does not draw finite receptor populations or stochastic openings. The drive is an EXPECTED open fraction, which can legitimately lie below 1/N; only an instantaneous realized count is integer-valued. The floor can be a useful modeling regularizer, but is not by itself a finite-copy-number simulation or a proof of physical resolution.

Under equal independent assembly, a particular homomer has abundance g^-5. At N=10000 its expected molecule count is:

| g | Expected copies of a particular homomer |
|---:|---:|
| 3 | 41.15 |
| 5 | 3.2 |
| 8 | 0.305 |
| 10 | 0.1 |
| 26 | 0.000842 |

At g=8, the probability that this homomer exists in a finite random population is approximately 1-exp(-0.305)=26.3%. The deterministic W model nevertheless includes its fractional contribution in every cell. This is acceptable as an ensemble mean, but matters when interpreting rare-type sensitivity or thresholds near a single channel.

For independent channels with fixed counts N_r, a simple variance model is Var(S|X)=sum_r N_r p_r(1-p_r)/N². Shared noise or stochastic assembly modifies it. Even without implementing this immediately, compare the separation of drives from theta against a stated noise scale.

Separate receptor opening sharpness from numerical cell binarization. Currently annealing the receptor changes the physical drive that cells pool. Choose whether its final temperature represents an assumed physical dose-response slope or an ideal binary-receptor limit; these are different scientific experiments. A deterministic threshold on an expected current should be reported as an idealized upper-performance model until noise robustness is measured.

Per-cell thresholds also do not inherently force all cells to fire 50% of the time: per-cell median calibration would. Distinguish the parameterization from the calibration policy. Shared median recalibration introduces adaptive homeostasis and changes when the cell population changes. For a nested cell-count capacity comparison, freeze an existing readout while adding cells; otherwise changes in the shared threshold confound the effect of adding outputs.

## 5. Concrete memory boundaries

Let S_pad be `env.s_upper`, the allocated ligand axis, and P be the number of energy sources (genes or directed pockets). With composition binding, principal float32 tensors include:

| Tensor | Bytes |
|---|---:|
| Source energies | 4 B S_pad P |
| Full receptor log-EC50 / binding terms, each | 4 B S_pad R_pool |
| Chunked receptor terms, each | 4 B S_pad r_chunk |
| Dense abundance matrix W | 4 C R_pool |
| Dense composition matrix | 4 P R_pool |
| KT pair-cell intermediates, uncheckpointed aggregate scale | O(4 B² C) |
| KT inference tile and row buffers | O(4 m² C + 4 m B) |

These are individual allocations/scalings, not peak memory forecasts. Autograd, temporaries, allocator behavior, and coexisting computations add overhead. `cell_pool_chunk` plus checkpointing reduces retained receptor intermediates; it does not remove dense W/composition, source-energy work, or total receptor computation. In the unbound interface path, the full per-ring energy tensor is already built before cell chunking.

For 100 cells sampled with exactly g genes, n_genes=26, interface model, seed 0, the actual pool builder gives:

| g | Distinct sampled gene sets | R_pool | One B=4096, S_pad=8 receptor tensor |
|---:|---:|---:|---:|
| 2 | 88 | 554 | 0.068 GiB |
| 3 | 99 | 4,190 | 0.511 GiB |
| 5 | 100 | 53,456 | 6.525 GiB |

Multiply the final column by S_pad/8 for another ligand-axis size. A single presence block gives S_pad=L even if only about one ligand is actually selected. Sparse mixture statistics therefore do not necessarily imply sparse allocated physics tensors. The docs' suggestion that chunking is unnecessary at five genes per cell is not general, especially in interface mode.

An additional bottleneck precedes training: `compute_initial_temperature` materializes the full B_calibration × S_pad × R_pool contraction and subsequent terms. It does not use `cell_pool_chunk`. Making training chunked therefore does not guarantee initialization fits. Cell calibration does chunk receptor physics, but the initial receptor-temperature calibration does not.

`resolve_batch_sizes` still budgets using full L × R_pool × k_sub rather than the actual composition/chunked path, and its minimum batch of 512 can exceed a genuinely smaller memory-safe value. It is a heuristic, not an OOM guarantee. Four A100s also do not automatically provide one combined memory pool: the runner uses one device; current task scripts distribute independent runs.

Exact entropy is a separate wall. At B=4096, the final B × 2^C float32 array alone occupies:

| C | Allocation |
|---:|---:|
| 6 | 1 MiB |
| 12 | 64 MiB |
| 16 | 1 GiB |
| 20 | 16 GiB |
| 26 | 1 TiB |

There is no universal hard cutoff at 15 cells; it depends on B and implementation. But exact enumeration clearly cannot support an unrestricted 26–100-cell campaign. The current identity/concentration metrics silently select it for KT runs.

## 6. What the expression data implies

Directly read `data/octopus_informationCoding/20260302_HiPlexResults.mat`, using the same Boolean interpretation as the expression analysis:

- 22 measured genes and 2,545 nonempty cells.
- 27 `suckerIDX` groups, with 36–203 cells each; median 79.
- 577 distinct gene-expression sets across the dataset.
- 9.82% of cells express more than five measured genes; 2.59% express more than eight.
- Two cells express all 22 measured genes.

These are properties of this file, not estimates of all genes or all sensory cells in an animal. The intended simultaneous array should be specified at a sucker or other justified spatial scale; pooling all recorded cells into one array is another assumption.

If all recorded cells are pooled and every permitted receptor is enumerated, either all-22-gene cell forces the full pool: 65,780 standard types or 1,030,744 interface rings. For the interface case, W alone at C=2545 occupies 9.77 GiB, and a 484 × R_pool composition matrix another 1.86 GiB, before any sniff tensor. Construction also enumerates 22^5 ordered words per all-gene cell in Python and is performed again by the runner.

Identical gene sets can be merged for deterministic information calculations. However, preserve their multiplicities for population-weighted threshold calibration and population comparisons. Dropping duplicates before median calibration changes the model unless those multiplicities are retained.

The current `size_pmf` sampler is a useful exchangeable baseline, but it does not reproduce empirical gene-specific frequencies or coexpression clusters: it draws uniformly conditional on gene count. Your Bernoulli-mixture/Ising/drive analyses address structures that a size PMF alone cannot represent. Start with explicit empirical gene sets or compare them with both a matched-count null and a marginal-frequency null. A five-gene cap is a computational/modeling choice that removes about 10% of the observed cells in this dataset; it should be an explicit sensitivity analysis.

One adjacent documentation issue: `genetic_measurements/drive_model.md` states exact equivalence between a quadratic count MaxEnt model and a Gaussian latent drive with conditionally independent logistic genes. Integrating normalized logistic probabilities under a fixed Gaussian prior includes a product of normalization factors depending on the drive. The elementary Hubbard–Stratonovich transform of the unnormalized quadratic model instead induces a reweighted latent prior. The stated equivalence is not generally exact; nor are logistic-normal marginals exactly exponential in gene rank. This does not invalidate fitting the implemented latent-drive model, but matters when transferring its interpretation into cell sampling.

## 7. Equivalence and convergence: current findings

### Equivalence has an initialization bug

`config.py:200` honors explicit `cell_receptors` and derives the eight-receptor pool correctly. `run.py:365` rebuilds `CellArray` from `cell_gene_sets` only. It expands the genes into all permitted receptors and ignores the explicit repertoire. For the actual equivalence list, the rebuilt pool has 16 receptors rather than eight, and `_initialize()` raises:

```text
cell pool mismatch: CellArray rebuilt a different receptor pool than the one stored in config.receptor_indices.
```

Preserve the explicit-repertoire path when rebuilding the array. This also blocks the clean homomer-only-cell control whenever its expressed genes would otherwise expand into heteromers.

### The intended equivalence needs two separate tests

For W=I and `mean`, forward outputs equal receptor probabilities. That is a genuine implementation identity; compare gradients as well as outputs on the same frozen inputs. The existing module self-test verifies the forward case, but misses the runner reconstruction above.

For `threshold`, hard receptor codes use p>0.5 while hard cell codes use p>theta. They coincide for arbitrary p only if theta=0.5, or if sampled p avoids the interval between thresholds. Flooring theta at 1e-4 does not establish equivalence. A direct example p=[0.01,0.1,0.4,0.6] yields receptor codes [0,0,0,1] and cell codes [1,1,1,1] at that threshold.

At receptor temperature T, the cell changes the effective log-drive threshold to T log(theta/(1-theta)). At theta=1e-4 and T=0.1, this is approximately -0.921 instead of zero, corresponding to about 0.40 times the concentration threshold for a single ligand. The difference tends to vanish in a suitable hard-receptor limit, but is material at finite temperature.

The two optimization runs also do not receive identical sniffs simply because the initial torch seed matches: threshold-cell calibration and recalibration consume extra random batches. Their annealing denominators differ as well. Compare frozen functions on identical inputs first; treat independently optimized endpoints as a separate statistical comparison across seeds.

Equal entropy and equal distinct-code counts do not prove equality of codewords. And the assertion that receptor KT must exceed cell KT is not general: thresholding probability vectors is not the same as deterministically processing a sampled Bernoulli output. After separate optimization there is even less reason for such an ordering. Remove that as an invariant.

### Convergence is a useful benchmark, but not a guaranteed attainable optimum

Six output bits can in principle distinguish eight inputs. That does not prove the sampled six three-gene cells, shared threshold, fixed repertoire weights, and shared interface parameters can realize the required mapping in this particular geometry. Construct an explicit realizable fixture before labeling every shortfall “under-trained.” Multiple initialization seeds and a fixed-temperature training tail can distinguish some optimization failures from structural constraints.

The current input is only approximately the stated finite world. `mu_ligands_per_source=1e-6` still assigns nonzero probability to multiple ligands; `conc_std=1e-4` is nonzero continuous variation. Tiny concentration variance is not a mathematical guarantee of zero concentration information if a sharp threshold can resolve it. For the unit fixture, enumerate eight singleton inputs at exactly fixed concentration and zero noise; inject them directly into the physics if the concentration-distribution constructor requires positive scale.

`identity_channel` and `concentration_channel` (`analysis_helper.py:858`) use exact Shannon enumeration except in collision mode. For the six-cell fixture this is appropriate and cheap. Their sum is exact empirical soft-mixture H, not necessarily the KT lower bound displayed alongside it. H(Y|M) includes readout softness and observation noise as well as concentration variation. For general mixtures, singleton mask groups are skipped, which can spuriously reduce conditional entropy and inflate reported identity information.

Require exactly eight distinct hard codewords on the eight deterministic fixture inputs, with a one-to-one label-to-code map and a minimum drive margin. The current `K_hat>=n_lig-1` criterion is permissive and does not test that map. The analysis prints PASS/FAIL but does not fail a test process. Its advice to increase epochs is a hypothesis, not a diagnosis.

### Final-temperature and sampler issues

The cell annealing denominator in `run.py:796` means the final training epoch never reaches the requested cell temperature. `run()` saves and tests the final training state without explicitly setting that endpoint; periodic measurements temporarily use the endpoint. For 300 epochs and phase split 0.8, the last multiplier is 0.0265 rather than 0.01, a factor 2.65. For the 5000-epoch convergence script it is 0.010396, about 4% high. This can make final metrics disagree with the plotted endpoint. Include a fixed-final-temperature tail and explicitly set/save the evaluated state.

The sparse ligand budget is a quantile heuristic rather than a support bound. `_sample_masks` retains only `s_upper` selected ligands in its sparse index but returns the full mask. In a supported configuration with L=K=100 and both Poisson parameters 1, s_upper=6. With torch seed 123 and 100,000 draws, 14 masks contained more than six ligands, up to eight. Those rows lose active ligands in physics while their identity labels retain them. Address overflow before using precise mixture-information comparisons. The present K=1 cell tasks do not hit this truncation case.

Also mask absent/padding concentrations exactly rather than relying on `log(c+1e-12)`: with sufficiently favorable learned energies, nominally absent ligands can contribute. This is a model-consistency issue at extreme affinities, not necessarily the dominant error in the current toy.

## 8. Validation completed and recommended next steps

Executed existing checks:

| Check | Result |
|---|---|
| `python -m src.cells` | All module self-tests passed |
| `python -m src.bin_loss` | All KT self-tests and gradcheck passed |
| `pytest unit_test/test_hierarchical_presence.py -q -p no:cacheprovider` | 10 failed, 1 passed |
| Actual equivalence configuration, runner initialization only | Reproduced 8-versus-16 receptor pool assertion |
| Interface functional count, three genes | 51 rings, 48 pocket-count vectors |
| Off-diagonal collision counterexample | 134.784 reported bits for 16 distinct four-bit inputs |
| Sparse-mask/sparse-index consistency probe | Reproduced 14 truncated rows in 100,000 |

The ten presence-suite failures occur because tests expect `_sample_masks()` to return a tensor; it now returns `(mask, sparse_idx)`. They do not themselves establish that the statistical sampler is wrong. Repair that test interface and include consistency checks between masks, sparse indices, and dense reconstructed concentrations. No locally stored results under the current `data/equivalence` or `data/convergence` task paths were found, so narrative claims about those runs were not independently verified against saved checkpoints.

Recommended order:

1. Repair explicit repertoire reconstruction, test return-type drift, final-temperature handling, and sparse overflow. Add frozen forward/gradient equality checks for composition versus gather and chunked versus unchunked paths in both receptor models. Keep these tests tiny.
2. Establish the exact eight-input fixture, with direct code mapping and exact MI. Test duplicates, constant output, and all-0.5 output as known counterexamples. Add a constructive achievable target and seed the world and optimizer initialization independently.
3. Log output entropy, conditional Bernoulli entropy, their difference, hard-code entropy, marginal firing rates, distinct cell repertoires, response margins, and evaluation sample counts. Check sample doubling on a frozen model.
4. Profile a few bounded jobs before sweeps: C=6–12 and g=1–3 first; then C=20,50,100 with g=2,3; g=5 only with measured pool cost and bounded initialization. Keep B around 2048–4096 initially. These are candidate starting points based on allocations, not measured A100 throughput guarantees. Explicitly budget KT evaluation and repeat count; `_test` performs ten evaluations.
5. Compare full repertoires against homomer-only repertoires using the SAME gene sets, ligand world, number of cells, molecule budget, and threshold policy. Run both frozen-parameter ablations and separately optimized architectures: they answer different questions. Standard versus interface also changes trainable parameter count, so do not attribute its whole gain solely to receptor ordering.
6. Introduce empirical per-sucker gene sets, compare against expression null models, and quantify the influence of high-gene-count cells. At that point decide whether finite assembly/current noise is required for the biological claim.

The first campaign should establish reliable information differences and sample convergence, not aim to fill C bits. The likely practical regime is tens to roughly a hundred simulated cells, small expressed-gene sets, and much lower reliably resolved joint entropy than the raw cell count. The code can represent more cells, but the receptor expansion, quadratic KT work, empirical alphabet coverage, and missing molecular noise impose separate limits that must be reported separately.

Ok, So these are large review, let's take section by section.
1. The Hard bound limits are well understood. OK. The Bounds that shouldn't guide the campgin. You are right that they are empirical bound, and not mathematically proven. But I think that it's fine. We have a vague rational for them, we cannot do more. It's fine.
2. I didn't exactly understood what limits the KT bound. If I understand correctly what you said, you derived some scaling of the number of bits as a function of the batch size which is the main limiting factor. Previously the limit I coud reach was around 20bits. I don't think I can beat that really. Although, you proposed to use the output counting (it has a name no ?) Do you think it's sharper than the KT bracketing ? Maybe because it's less sensitive to batch size ? Concerning the collision, I don't think I need it anymore, it's less sharp, and harder to use than the KT lower-bound anyway.
3. Right so the receptor repertoir is the main limiting factor in te case of simulating cells I understand. I think that I can get away with it by  first reducing the total number of genes expressed. I can demonstrate the implication of multiple genes expression and heteromerization by simulating 5 different genes I think. Concerning the AND/OR interpretation: You've said that they are not strictly true, ok, it's just an approximation to draw expectations. I feel that in your analysis, what you miss is the weight of the receptors, when expressing a finite number of genes in a cell: most receptors are heteromers. so the response of heteromers dominate the behavior. Then concerning the activation of the cell, this is an element that I'm still unsure about. The noisy_or seems really odd to me, but I'm still kinda hesitating between a thresholded activation, or a mean. The mean makes the output stochastic, whereas my previous model were deterministic, but that shouldn't be too hard to fix, isn't it ?
4. so you'r correct on the finite molecule count. It's more of a regularizer than a true physics idea. The point being that in the limit of hard thresholded activation, I want to avoid the optimizer to put the activation in the queue of the sigmoid to get a fair coin for each cell, thus exploiting the optimizer, to get 0 information. Maybe I should change my optimizer from the entropy being equal to the mutual information only if the output is deterministic. I've figured to work in the deterministic limit and only optimize the entropy, but maybe I could actually maximize the information, what do you think ?
7.
