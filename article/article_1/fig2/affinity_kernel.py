"""
EC50 as a function of latent ligand-unit distance.

Saturating affinity kernel (doc 03_latent_environment.md, sec 3.1):

    E_o^(u,l) = E_base^(u) + E_max^(u) * (1 - exp(-||v_u - v_l||^2 / lambda^2))

For a single unit (homomer) the half-activation concentration is

    EC50(d) = exp[ E_base + E_max * (1 - exp(-d^2 / lambda^2)) ]

For a k_sub-heteromer it is the geometric mean over subunits, i.e. the
arithmetic mean of the per-unit open-state energies in log space:

    ln EC50(d) = (1/k_sub) * sum_u [ E_base^(u) + E_max^(u)*(1 - exp(-d^2/lambda^2)) ]

Limits:  d -> 0    (perfect match)     EC50 -> exp(E_base)             (most sensitive)
         d -> inf  (full mismatch)     EC50 -> exp(E_base + E_max)     (saturates)
"""

import numpy as np
import matplotlib.pyplot as plt


def ln_ec50(d, E_base, E_max, lam):
    """log EC50 for a single unit as a function of distance d."""
    return E_base + E_max * (1.0 - np.exp(-(d / lam) ** 2))


def ec50(d, E_base, E_max, lam):
    return np.exp(ln_ec50(d, E_base, E_max, lam))


# ---- parameters (edit these) ------------------------------------------------
E_base = 0.0        # log-EC50 at the optimal ligand (d = 0)
E_max  = 6.0        # selectivity ceiling: extra log-units for a full mismatch
lam    = 1.0        # affinity length scale (hyperparameter)
d_star = 0.9        # a distance to highlight (e.g. the arrow in the UMAP)

d = np.linspace(0.0, 2.5 * lam, 400)

# ---- plot -------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(4.2, 3.2))

ax.plot(d, ec50(d, E_base, E_max, lam), color="#1f6f9e", lw=2.2, zorder=3)

# asymptotes
ax.axhline(np.exp(E_base),         ls=":", lw=1.0, color="#777")
ax.axhline(np.exp(E_base + E_max), ls=":", lw=1.0, color="#777")
ax.text(3. * lam, np.exp(E_base) * 1.15, r"$c^*_\text{min}$: complete match",
        ha="right", va="bottom", fontsize=9, color="#555")
ax.text(3. * lam, np.exp(E_base + E_max) * 0.62,
        r"$c^*_\text{max}$: complete mismatch",
        ha="right", va="top", fontsize=9, color="#555")

# operating point at d_star
#y_star = ec50(d_star, E_base, E_max, lam)
#ax.plot([d_star, d_star], [ec50(0, E_base, E_max, lam), y_star],
#        ls="--", lw=0.9, color="#c0392b")
#ax.plot([0, d_star], [y_star, y_star], ls="--", lw=0.9, color="#c0392b")
#ax.plot(d_star, y_star, "o", color="#c0392b", ms=6, zorder=4)
#ax.annotate(r"$EC_{50}(d)$", (d_star, y_star), (d_star + 0.15, y_star * 0.55),
#            fontsize=9, color="#c0392b")

# lambda marker on the x-axis
#ax.axvline(lam, ls="-", lw=0.6, color="#bbb", zorder=0)
#ax.text(lam, ec50(0, E_base, E_max, lam) * 1.05, r"$\lambda$",
#        ha="center", va="bottom", fontsize=9, color="#999")

ax.set_yscale("log")
ax.set_xlabel(r"latent distance  $d = \|\mathbf{v}_u - \mathbf{v}_\ell\|$",fontsize=9)
ax.set_ylabel(r"$c^*$  (log scale)",fontsize=9)
ax.set_xlim(0, 2.5 * lam)
ax.margins(y=0.15)
ax.spines[["top", "right"]].set_visible(False)
ax.tick_params(direction="in",which="both")
ax.set_xticks([]); ax.set_yticks([])
fig.tight_layout()

#fig.savefig("ec50_vs_distance.png", dpi=200)
fig.savefig("ec50_vs_distance.svg")   # vector version for the figure
plt.show()