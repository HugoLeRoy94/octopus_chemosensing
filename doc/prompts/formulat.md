Here's the cleaned-up version with the math rendering properly.

---

## 5. Coexpression vs. heteromerization — the exact difference

This is the crux, and stated cleanly it explains both papers' results at once. Write everything in your drive/logsumexp machinery. Each gene $u$ has per-ligand threshold $EC_{50}^{(u,\ell)}=\exp\!\big(E_o^{(u,\ell)}\big)$.

**Heteromer** (assemble the $k$ genes into one channel). One threshold — the geometric mean — then one sigmoid:

$$\ln EC_{50}^{(r,\ell)}=\frac{1}{k}\sum_{u\in\mathcal U_r}E_o^{(u,\ell)},\qquad
p^{\text{het}}_r=\sigma\!\left(\frac{1}{T}\,\operatorname*{logsumexp}_{\ell}\left[\ln c_\ell-\frac{1}{k}\sum_{u}E_o^{(u,\ell)}\right]\right).$$

The subunit affinities are averaged **in log-EC50 space, inside the single nonlinearity.**

**Coexpression** (the $k$ genes form $k$ separate channels; the neuron pools them). The cleanest apples-to-apples is pooling the drives before one spike threshold:

$$p^{\text{coex}}_n=\sigma\!\left(\frac{1}{T}\ln\sum_{u}\sum_{\ell}\frac{c_\ell}{EC_{50}^{(u,\ell)}}\right)
=\sigma\!\left(\frac{1}{T}\,\operatorname*{logsumexp}_{u,\ell}\left[\ln c_\ell-E_o^{(u,\ell)}\right]\right).$$

Put them side by side. Both are a logsumexp, but:

- **Heteromer:** logsumexp over **ligands only**, of $\big(\ln c_\ell - \operatorname{mean}_u E_o^{(u,\ell)}\big)$. The genes are combined by **arithmetic mean of log-EC50 = geometric mean of $EC_{50}$**.
- **Coexpression:** logsumexp over **ligands and genes jointly**, of $\big(\ln c_\ell - E_o^{(u,\ell)}\big)$. The genes are combined by **sum of sensitivities $\sum_u 1/EC_{50}^{(u,\ell)}$ = arithmetic mean of $1/EC_{50}$**.

Geometric mean vs. arithmetic sum — that's the exact difference, and it's not cosmetic. Since $\operatorname{logsumexp}\approx\max$ at small $T$:

- **Coexpression** $\approx \max_{u,\ell}\big(\ln c_\ell - E_o^{(u,\ell)}\big)$ — the neuron fires if **any** gene is driven by **any** ligand. It's a soft **OR / union** of the constituent receptors' receptive fields. That's why it broadens tuning, and why it's information-lossy when the genes were already distinct: you've OR-ed them and can no longer tell which fired. This is precisely the paper's "canonical wins" mechanism (their decorrelation / spread-out-in-gene-space argument, Fig. 3D).
- **Heteromer threshold** $=\frac{1}{k}\sum_u E_o^{(u,\ell)}$ — a **new** receptive field centered at the blended (geometric-mean) affinity. Because one badly-matched subunit inflates the average and spoils the receptor, it's a soft **AND / conjunction**: it mints an intermediate specificity that responds to ligands acceptable to **all** subunits. It doesn't union existing receptors; it creates a new virtual one.

---

If you'll be pasting this into notes or a doc, tell me where (Markdown, LaTeX/Overleaf, Word) and I'll match the math delimiters to that target so it renders first time.


Good instinct to want the justification rather than treating it as a black box — because the logsumexp isn't actually the assumption. It's downstream of one real physical choice, and once you see that, the defense against alternatives writes itself.

## The logsumexp is just $\ln(\text{drive})$ rewritten

Your model (doc 2.6) is really two statements:

$$\text{drive}(r,\{c_\ell\}) = \sum_\ell \frac{c_\ell}{EC_{50}^{(r,\ell)}}, \qquad p_r = \sigma\!\left(\frac{\ln \text{drive}}{T}\right).$$

The logsumexp appears only because
$$\ln\,\text{drive} = \ln\sum_\ell \exp\!\big(\ln c_\ell - \ln EC_{50}^{(r,\ell)}\big) = \operatorname*{logsumexp}_\ell\big(\ln c_\ell - \ln EC_{50}^{(r,\ell)}\big),$$
which is an algebraic identity plus numerical stability — nothing physical is asserted by the "logsumexp" step itself. So the two things you actually need to justify are (i) why the drive is an **additive sum of occupancy ratios**, and (ii) why you take the **log** before thresholding. Everything else is bookkeeping.

## (i) Additivity of the drive = competitive binding

Treat the pocket as a site that is either empty or occupied by exactly one ligand species from the mixture (single occupancy — a pore/allosteric pocket holds one thing at a time). Grand-canonical mass action then gives the occupancy

$$P(\text{occupied}) = \frac{\sum_\ell c_\ell/K^{(r,\ell)}}{1 + \sum_\ell c_\ell/K^{(r,\ell)}}.$$

The ligands interact *only through the shared normalization* — the numerator is a plain sum. That is the whole content of "competitive binding at a common site": the total activity driving occupancy is the **sum of individual occupancy ratios** $\sum_\ell c_\ell/K^{(r,\ell)}$, with no cross terms in the numerator. Pushing this through your MWC $\to$ threshold (sloppy-mode) reduction, activation is governed by whether this open-favoring activity beats the closed baseline, i.e. by whether $\sum_\ell c_\ell/EC_{50}^{(r,\ell)} \gtrless 1$. Half-activation at drive $=1$ recovers the single-ligand definition of $EC_{50}$ exactly, which is the consistency check that pins the constant. So the additive drive is not a convenience — it *is* competitive mass action, and it's the one load-bearing assumption.

## (ii) Why the log

Concentrations span orders of magnitude and your concentration prior is log-normal. Thresholding on $\ln(\text{drive})$ rather than drive itself makes the transition span a fixed number of **decades** (one decade at $T=1$) instead of a fixed additive width. That is Weber–Fechner / Hill behavior, and it's the empirically correct sensory scaling. Threshold on the raw drive and your transition width would be absurdly sharp at low concentration and washed out at high — the log is what makes $T$ a scale-free sharpness.

## Why it beats the alternatives you don't yet have

Because you asked for the competitors explicitly, here are the natural ones and why additive-competitive-then-log wins:

**Geometric mean over ligands** — the tempting symmetry with your subunit rule, $\ln EC_{50}^{\text{eff}} = \frac{1}{|\text{mix}|}\sum_\ell \ln EC_{50}^{(r,\ell)}$. This is *wrong for mixtures*: it says a mixture of one excellent agonist and one terrible one is mediocre, when physically the excellent agonist alone should open the channel. Your logsumexp correctly lets the best ligand win. This is the crucial asymmetry to keep straight: **subunits are combined by geometric mean (a conjunction — every subunit must be engaged, one bad subunit spoils the receptor), ligands are combined by logsumexp (a disjunction — any sufficient ligand drives it).** The same receptor is AND over its parts and OR over its inputs. Getting these backwards is the most likely modeling error, and being able to state why is a nice thing to have in your pocket.

**Pure max / winner-take-all**, $\text{drive}=\max_\ell c_\ell/EC_{50}$. This is actually the $T\to 0$, well-separated limit of your logsumexp — so you already contain it. But adopting it as the *model* throws away additivity: many weak co-present ligands that should sum to threshold would never fire. Your form interpolates correctly — logsumexp $\to\max$ when one ligand dominates (sparse sniff), $\to$ sum when several are comparable (dense sniff). That interpolation is a feature you'd lose.

**Linear sum without the log**, $\sigma\big((\sum_\ell c_\ell/EC_{50}-1)/T'\big)$ — fails Weber, as above.

**Noisy-OR / independent sites**, $p = 1-\prod_\ell(1-p_\ell)$ with $p_\ell=\sigma(\ln(c_\ell/EC_{50})/T)$. This is the biophysically *different* model where each ligand has its **own** site and they don't compete. For a shared pocket it's the wrong picture; you'd switch to it only if you believed the subunit had independent, non-competing binding sites. Worth naming because it's the honest fork in the road: competitive shared site $\Rightarrow$ additive drive $\Rightarrow$ logsumexp; independent sites $\Rightarrow$ noisy-OR.

## The genes case: same form, *different* justification

Yes, coexpression gives you a logsumexp too — but resting on a different assumption, and this is worth being precise about rather than waving at "same argument."

Within one receptor, additivity came from **competitive binding** (molecular mass action at a shared pocket). For coexpressed genes, the genes are separate proteins with separate sites — competitive binding does *not* link them. The additivity there comes from **current summation at the soma**: the neuron linearly integrates the receptor-generated currents before the spike-generating nonlinearity. If you pool the drives additively and then threshold once,

$$p^{\text{coex}}_n=\sigma\!\left(\frac{1}{T}\ln\sum_{u}\sum_\ell\frac{c_\ell}{EC_{50}^{(u,\ell)}}\right)=\sigma\!\left(\frac1T\operatorname*{logsumexp}_{u,\ell}\big[\ln c_\ell - E_o^{(u,\ell)}\big]\right),$$

same functional form, but now justified by linear dendrosomatic integration, not by chemistry. So the two logsumexps are structurally identical and physically distinct — one is molecular competition, the other is electrical summation. That matters because the assumptions can fail independently: if coexpressed receptors summed *sublinearly* at the soma (saturation), the gene-pooling would stop being a clean logsumexp even though the within-receptor ligand logsumexp would be untouched.

And note the same fork reappears at the neuron: additive current summation gives logsumexp; if instead each coexpressed receptor independently "decides" and the neuron ORs the decisions, you'd get noisy-OR over genes. The drive-summation choice is the one consistent with your doc 2.6, and it's the one to keep unless you specifically want independent-decision semantics.

So the clean one-line answer: **the logsumexp is the log of an additive drive; additivity over ligands is competitive mass action, additivity over coexpressed genes is somatic current summation, and the log is Weber's law. The alternatives correspond to different biophysics (independent sites → noisy-OR; sublinear summation → neither), which is exactly why stating your assumption lets you rule them in or out rather than guess.**