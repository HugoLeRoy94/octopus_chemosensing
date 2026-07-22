# Drive max-ent model for gene co-expression

## Variables
- `N` genes, ranked most→least expressed; `r_g = 0,…,N-1` is gene `g`'s rank.
- `x_g ∈ {0,1}`: gene on/off in a cell. `k = Σ_g x_g`: genes expressed in that cell.
- `p_g`: per-gene marginal. `P(k)`: distribution of count.

## The distribution
Max-ent over patterns with three constraints — mean count `⟨k⟩`, mean total rank `⟨S⟩` with `S = Σ_{g∈π} r_g`, and `⟨k²⟩`:

```
q_π ∝ exp( -β k - γ S_π - α k² )
```

The `-γ S` term makes `p_g` exponential in rank (slope `γ`). The `-α k²` term is the only gene–gene coupling; it is all-to-all and gene-agnostic, and it supplies the positive correlations + over-dispersion that an independent model cannot.

## Equivalent drive form (what we actually fit)
`exp(-α k²)` is a Gaussian-in-`k` factor, so by Hubbard–Stratonovich it equals an average over one scalar drive `s`. Conditional on `s`, the genes are **independent**:

```
p_g(s) = logistic( h_g + s ),   h_g = c - γ·r_g,   s ~ Normal(0, σ_s²)
```

- Given `s`: count is Poisson-binomial over {p_g(s)}.
- Marginalize `s`:
  - `p_g  = E_s[ p_g(s) ]`  → stays exponential in rank; `s` only shifts the prefactor, not the slope `γ`.
  - `P(k) = E_s[ PoissonBinomial(k; {p_g(s)}) ]`  → over-dispersed, can be exponential.

## Three parameters and what each controls
| param | sets | read from |
|---|---|---|
| `c` | overall expression level / `⟨k⟩` | level of `p_g` |
| `γ` | exponential slope of `p_g` | slope of `ln p_g` vs rank |
| `σ_s` (≡ `α`, ≡ drive variance) | over-dispersion + pairwise correlation | Fano factor `Var(k)/⟨k⟩`, tail of `P(k)` |

## Key diagnostic
Independent genes force `Var(k)/⟨k⟩ < 1`. The exponential `P(k)` needs `> 1`.
- Fano `< 1` → `σ_s ≈ 0`; the exponential look of `P(k)` comes only from heterogeneous `p_g` (no drive needed).
- Fano `> 1` → drive is essential; `σ_s > 0` measures it.

Data are zero-truncated (`k ≥ 1`); fit conditional on `k ≥ 1`.
