#!/usr/bin/env python3
"""
Pentamer of half-coloured discs, with a circular notch bitten out near each
contact.

Geometry
--------
  N subunit centres on a circle of radius R about the centre of mass O.
  Discs of radius r = R sin(pi/N)  ->  neighbours touch exactly at the
  midpoint of each polygon edge, at distance Rc = R cos(pi/N) from O.
  Each disc is split by the diameter collinear with O -> centre, so the
  colour boundary is a radius of the whole assembly.
  The clipping circle of radius rho sits on the same radial ray, centred at
  distance Rc + rho from O, so that its *inner edge* passes through the
  contact point.  It is subtracted from both discs meeting there.

All arcs are exact SVG arcs (no polygonal approximation), so the paths stay
editable as arcs in Inkscape.
"""
import math

# ---------------------------------------------------------------- parameters
N          = 5          # subunits
R          = 100.0      # centre of mass -> subunit centre
RHO_FRAC   = 0.4       # notch radius / subunit radius
OFFSET     = 1.0        # notch centre at Rc + OFFSET*rho  (1.0 -> edge on contact)
COL_A      = "#3d6fa8"  # leading half
COL_B      = "#e2a13c"  # trailing half
STROKE     = "#17191c"
STROKE_W   = 1.8
GUIDES     = True       # pentagon + radii on a separate layer
DEBUG      = False       # overlay the clipping circles, to verify the cut
MARGIN     = 14
OUT        = "pentamer.svg"
# ---------------------------------------------------------------------------

r   = R * math.sin(math.pi / N)          # disc radius (touching condition)
rho = RHO_FRAC * r                       # clipping-circle radius
Rc  = R * math.cos(math.pi / N)          # contact point, distance from O


def pt(angle, rad=1.0, c=(0.0, 0.0)):
    return (c[0] + rad * math.cos(angle), c[1] + rad * math.sin(angle))


def f(p):
    return f"{p[0]:.4f},{p[1]:.4f}"


def wrap(x):
    return (x + math.pi) % (2 * math.pi) - math.pi


theta   = [-math.pi / 2 + 2 * math.pi * i / N for i in range(N)]   # first on top
phi     = [t + math.pi / N for t in theta]                         # contact rays
C       = [pt(t, R) for t in theta]                                # disc centres
CONTACT = [pt(p, Rc) for p in phi]                                 # touching points
NOTCH   = [pt(p, Rc + OFFSET * rho) for p in phi]                  # clip centres


def half_path(i, s):
    """One coloured half of disc i.  s = +-1 selects the side of the
    colour-splitting diameter; the notch on that side is cut out."""
    c, phiB = C[i], theta[i]
    A = pt(phiB + math.pi, r, c)     # inner end of the diameter
    B = pt(phiB, r, c)               # outer end

    # the clipping circle lying on side s of the diameter
    M = next(m for m in (NOTCH[i], NOTCH[(i - 1) % N])
             if wrap(math.atan2(m[1] - c[1], m[0] - c[0]) - phiB) * s > 0)

    # circle(c, r) x circle(M, rho), general case
    dx, dy = M[0] - c[0], M[1] - c[1]
    d      = math.hypot(dx, dy)
    if not abs(d - rho) < r < d + rho:
        raise ValueError("clipping circle does not cut the rim cleanly")
    phiM  = math.atan2(dy, dx)
    a     = (d * d + r * r - rho * rho) / (2 * d)   # along c -> M
    h     = math.sqrt(r * r - a * a)                # perpendicular
    Delta = math.atan2(h, a)                        # half-angle seen from c

    P1 = pt(phiM - s * Delta, r, c)
    P2 = pt(phiM + s * Delta, r, c)

    sw   = 1 if s > 0 else 0        # rim: direction of travel
    nsw  = 1 - sw                   # notch: opposite (concave bite)
    # the bite follows the arc of the small circle nearest to c; that arc is
    # the major one only if the chord lies beyond M, i.e. if a > d
    nlarge = 1 if a > d else 0

    return (f"M {f(A)} L {f(B)} "
            f"A {r:.4f},{r:.4f} 0 0 {sw} {f(P1)} "
            f"A {rho:.4f},{rho:.4f} 0 {nlarge} {nsw} {f(P2)} "
            f"A {r:.4f},{r:.4f} 0 0 {sw} {f(A)} Z")


ext  = R + r + MARGIN
body = [f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape" '
        f'width="{2*ext:.2f}mm" height="{2*ext:.2f}mm" '
        f'viewBox="{-ext:.2f} {-ext:.2f} {2*ext:.2f} {2*ext:.2f}">']

if GUIDES:
    body.append('<g inkscape:groupmode="layer" inkscape:label="guides" '
                'fill="none" stroke="#b0b0b0" stroke-width="0.7" '
                'stroke-dasharray="3 3">')
    body.append('<path d="M ' + ' L '.join(f(c) for c in C) + ' Z"/>')
    for m in NOTCH:
        body.append(f'<path d="M 0,0 L {f(m)}"/>')
    body.append('<circle cx="0" cy="0" r="1.4" fill="#b0b0b0" stroke="none"/>')
    body.append('</g>')

body.append('<g inkscape:groupmode="layer" inkscape:label="subunits" '
            f'stroke="{STROKE}" stroke-width="{STROKE_W}" '
            'stroke-linejoin="round">')
for i in range(N):
    body.append(f'  <g inkscape:label="subunit_{i+1}">')
    body.append(f'    <path fill="{COL_A}" d="{half_path(i, +1)}"/>')
    body.append(f'    <path fill="{COL_B}" d="{half_path(i, -1)}"/>')
    body.append('  </g>')
body.append('</g>')

if DEBUG:
    body.append('<g inkscape:groupmode="layer" inkscape:label="debug" '
                'fill="none" stroke="#d02020" stroke-width="0.8">')
    for m, k in zip(NOTCH, CONTACT):
        body.append(f'<circle cx="{m[0]:.4f}" cy="{m[1]:.4f}" r="{rho:.4f}"/>')
        body.append(f'<circle cx="{k[0]:.4f}" cy="{k[1]:.4f}" r="1.2" '
                    f'fill="#d02020"/>')
    body.append('</g>')

body.append('</svg>')

with open(OUT, "w") as fh:
    fh.write("\n".join(body))

print(f"r    = {r:.4f}   disc radius")
print(f"rho  = {rho:.4f}   clipping radius")
print(f"Rc   = {Rc:.4f}   contact point from O")
print(f"Rn   = {Rc + OFFSET*rho:.4f}   clipping-circle centre from O")
print(f"-> {OUT}")
