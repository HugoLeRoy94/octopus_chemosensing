
- We have $N$ genes. I index them by $g$, and I order them from most-expressed to least-expressed, so that gene 1 is the most frequently expressed. I write $r_g$ for the **rank** of gene $g$ in that ordering, using $r_g = g-1$ so that the top gene has rank 0.
- A single cell either expresses or does not express each gene. I write $x_g = 1$ if gene $g$ is on in that cell, and $x_g = 0$ if it is off.
- The total number of genes a cell expresses is $k = \sum_{g} x_g$.
- A **pattern** $\pi$ is the specific set of genes a cell expresses (for example, "genes 2, 5, and 9"). Its size is $|\pi| = k$.
- $q_\pi$ is the probability that a randomly chosen cell expresses exactly pattern $\pi$. This is what I have been calling the **pattern frequency**.
- $p_g$ is the probability that gene $g$ is expressed, regardless of what else is on. It is obtained from the pattern frequencies by summing over every pattern that contains $g$: $p_g = \sum_{\pi \ni g} q_\pi$. This is the **per-gene marginal**.
- $P(k)$ is the probability that a cell expresses exactly $k$ genes, obtained by summing the frequencies of all patterns of that size: $P(k) = \sum_{|\pi| = k} q_\pi$.

Those six objects ($x_g$, $k$, $\pi$, $q_\pi$, $p_g$, $P(k)$) are all we need. I'll introduce two multipliers ($\beta$ and $\gamma$) when they come up, and define them there.

---


### The first constraint: control the number of expressed genes

Your "first half" is the statement that the system does not want to express too many genes at once. The clean way to encode that, and nothing more, is to fix the average number of expressed genes, $\langle k \rangle$. So the quantity being constrained is $k_\pi = |\pi|$, and its multiplier I call $\beta$.

Here I have to be careful about something I previously glossed over, and it matters for your question.

There are two different maximum-entropy problems depending on what you treat as the elementary outcome.

If the elementary outcome is **the integer $k$ itself** — that is, you only care about the distribution $P(k)$ and treat all cells with the same $k$ as identical — then maximizing $-\sum_k P(k)\ln P(k)$ at fixed $\langle k\rangle$ gives a clean exponential, $P(k) \propto e^{-\beta k}$. This is the result you have been calling the first half, and it is correct **in that setting**.

If instead the elementary outcome is **the pattern $\pi$** — which is unavoidable the moment you also want to talk about $p_g$, because $p_g$ is a property of patterns, not of the bare count — then the same constraint gives $q_\pi \propto e^{-\beta k}$. Summing over all patterns of a given size, and recognizing that there are $\binom{N}{k}$ of them,

$$P(k) = \sum_{|\pi|=k} q_\pi = \frac{1}{Z}\binom{N}{k} e^{-\beta k},$$

where $Z$ is the normalization. This is **not** a pure exponential — it is a binomial distribution, peaked at an interior value of $k$, because the combinatorial factor $\binom{N}{k}$ counts how many ways there are to choose which genes are on.

So the honest statement is: fixing $\langle k\rangle$ gives an exponential $P(k)$ only if you work purely at the level of the count and ignore gene identity. As soon as you move to the pattern level — which you must, to discuss $p_g$ — fixing $\langle k\rangle$ alone produces a binomial in $k$, not an exponential. The difference between the two is exactly the entropy of "which $k$ genes are chosen," the $\ln\binom{N}{k}$ term. I should have flagged this much earlier; it is the root of several of the confusions.

### The second constraint: the rank information

The first constraint says nothing about *which* genes are expressed — every gene enters $k$ identically. To get a non-trivial $p_g$, we need a second constraint that distinguishes genes. The minimal one that does not invent gene-by-gene structure is to fix the average **total rank** of the expressed genes,

$$S_\pi = \sum_{g \in \pi} r_g, \qquad \text{constrain } \langle S \rangle, \quad \text{multiplier } \gamma .$$

In words: on average, how "deep" into the rank list does a cell reach? Fixing $\langle S\rangle$ penalizes expressing low-ranked (rare) genes, in a way that grows linearly with rank. It treats all genes by a single rule (their rank), so it introduces no pair-specific coupling.

### Solving the Lagrangian

With both constraints, the function we extremize is

$$\mathcal{L} = -\sum_\pi q_\pi \ln q_\pi - \lambda\Big(\sum_\pi q_\pi - 1\Big) - \beta\Big(\sum_\pi q_\pi k_\pi - \langle k\rangle\Big) - \gamma\Big(\sum_\pi q_\pi S_\pi - \langle S\rangle\Big),$$

where $\lambda$ enforces normalization. Setting the derivative with respect to a single $q_\pi$ to zero, $\partial\mathcal{L}/\partial q_\pi = -\ln q_\pi - 1 - \lambda - \beta k_\pi - \gamma S_\pi = 0$, and solving,

$$q_\pi = \frac{1}{Z}\, e^{-\beta k_\pi - \gamma S_\pi}.$$

Now substitute $k_\pi = \sum_{g\in\pi} 1$ and $S_\pi = \sum_{g\in\pi} r_g$, so that both exponents are sums over the genes in the pattern. The exponential of a sum is a product, so

$$q_\pi = \frac{1}{Z} \prod_{g \in \pi} e^{-(\beta + \gamma r_g)} .$$

This is the key structural fact: **the pattern probability is a product of independent per-gene factors.** There is no term that couples two specific genes. Each gene contributes a factor $e^{-(\beta + \gamma r_g)}$ if it is on, and a factor 1 if it is off.

### Reading off $p_g$

Because the distribution factorizes, the genes are statistically independent, and each is an independent yes/no variable with

$$p_g = \frac{e^{-(\beta + \gamma r_g)}}{1 + e^{-(\beta + \gamma r_g)}}.$$

When expression is sparse (the typical $p_g$ is small, so the denominator is close to 1), this simplifies to

$$p_g \approx e^{-\beta}\, e^{-\gamma r_g}.$$

This **is** the exponential-in-rank law for $p_g$, and it comes out exactly. The slope of $\ln p_g$ against rank is $-\gamma$. So the second constraint, fixing the average total rank, is precisely what makes the per-gene expression exponential. That part of my earlier claim is solid.

### Reading off $P(k)$, carefully

Now sum the factorized $q_\pi$ over all patterns of fixed size $k$:

$$P(k) = \frac{1}{Z}\, e^{-\beta k} \sum_{|\pi| = k}\ \prod_{g \in \pi} e^{-\gamma r_g}.$$

The remaining sum runs over all ways to pick $k$ genes out of $N$, multiplying their rank-factors. That sum has a name in algebra: it is the **elementary symmetric polynomial of degree $k$** in the variables $\{e^{-\gamma r_g}\}$, written $e_k$. The name is just shorthand for "sum, over all size-$k$ subsets, of the product of the chosen variables." It appeared in my earlier message only as this definition; I did not introduce any new physics with it, and I should have said in plain words what it was.

If we now use the specific rank ladder $r_g = g-1$, then $e^{-\gamma r_g} = (e^{-\gamma})^{g-1}$ is a geometric sequence $1, q, q^2, \dots$ with ratio $q \equiv e^{-\gamma}$. For a geometric sequence there is a classical closed form for the elementary symmetric polynomial,

$$e_k(1, q, q^2, \dots, q^{N-1}) = q^{\binom{k}{2}} \binom{N}{k}_{\!q},$$

and I need to define that last symbol, which I never did. The expression $\binom{N}{k}_{q}$ is the **Gaussian binomial coefficient** (also called the $q$-binomial). The subscript $q$ does not mean "raise to the power $q$"; it labels a generalized version of the ordinary binomial coefficient. It is defined as

$$\binom{N}{k}_{q} = \frac{(1-q^N)(1-q^{N-1})\cdots(1-q^{N-k+1})}{(1-q)(1-q^2)\cdots(1-q^k)},$$

and it has the property that when $q \to 1$ it reduces to the ordinary $\binom{N}{k}$. It is a polynomial in $q$ that counts the same subsets as $\binom{N}{k}$ but weighted by their rank content. I introduced it only because it makes the dependence on $\gamma$ explicit; it is not essential, and if it is more confusing than helpful we can keep $e_k$ unexpanded.

Putting it together,

$$P(k) = \frac{1}{Z}\, e^{-\beta k}\; q^{\binom{k}{2}} \binom{N}{k}_{\!q}, \qquad q = e^{-\gamma}.$$

The thing to take away, without any of the algebra, is the **shape** of $\ln P(k)$. Expanding, $q^{\binom{k}{2}} = e^{-\gamma k(k-1)/2}$, so

$$\ln P(k) = \text{const} - \beta k - \tfrac{\gamma}{2}\,k(k-1) + \ln\binom{N}{k}_{\!q}.$$

There is the linear term $-\beta k$ you wanted — but there is also a **quadratic term** $-\tfrac{\gamma}{2}k^2$, carrying the *same* $\gamma$ that set the slope of $p_g$. So $P(k)$ is a pure exponential only in the regime where that quadratic term is negligible compared to the linear one, that is, when $\gamma k \ll \beta$ — small $\gamma$, or only over a range of small $k$. Outside that regime $P(k)$ is exponential bent by a Gaussian factor.

