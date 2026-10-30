"""Measure the minimizer's behavior at each breakpoint from a weil_breakpoints.py --dump file.

usage: python3 bp_jumps.py bp_lam25_d20.json

The piecewise basis imposes no continuity at the cuts, so the converged eigenvector
can have a jump in value, in slope, or in curvature at each cut.  For every interior
breakpoint this prints the left and right limits of xi, xi', xi'' (from the adjacent
piece polynomials), the three jumps, and each jump divided by the Euler-Lagrange
weight Lambda(p^m) p^{-m/2} * v, where v = xi(0) is the endpoint value and the
breakpoint is identified with log p^m or L - log p^m.  Breakpoints that are not at a
prime-power position (second-order ones, or displaced controls) get weight '-'.

Reading the output.  The basis is even under y -> L - y by construction, so the row for
L - log q is the row for log q with the one-sided values swapped: its jumps are exactly the
negatives and carry no independent information.  When L/2 is itself a prime position
(log 2 at lambda = 2, log 3 at lambda = 3) the value jump there is identically zero by the
same symmetry.  The '/wv' column is jump / (Lambda(q) q^{-1/2} * v) with v = xi(0) the raw
(L^2-unit, sign-arbitrary) eigenvector value at y = 0; v cancels in ratios between cuts.
The note (section 6) uses the value-jump '/wv' only, across degrees: a resolved step would
make it degree-independent, a kink would send it to zero; the measured ratios do neither.
"""
import json, sys
from mpmath import mp, mpf, log, sqrt
import weil_step1 as W

path = sys.argv[1]
d = json.load(open(path))
mp.dps = max(30, d.get("dps", 30))
lam = mpf(d["lam"]); L = mpf(d["L"]); deg = d["d"]
B = [mpf(b) for b in d["B"]]
pieces = [(mpf(a), mpf(b)) for a, b in d["pieces"]]
funcs = d["funcs"]                       # list of [[piece, k, weight], ...]
pnorm = [[mpf(x) for x in row] for row in d["pnorm"]]
xi = [mpf(x) for x in d["xi"]]
npc = len(pieces)

# by_piece: piece -> list of (function index, k, weight, sign)
by_piece = [[] for _ in range(npc)]
for fi, parts in enumerate(funcs):
    for j, (pi, k, w) in enumerate(parts):
        sign = (-1) ** k if (len(parts) == 2 and j == 1) else 1   # mirror copy, as in evals()
        by_piece[pi].append((fi, k, mpf(w), sign))

# Legendre endpoint values and derivatives in the local variable z in [-1, 1]
def P(k, z):            # z = +1 or -1
    return mpf(1) if z > 0 else mpf((-1) ** k)
def dP(k, z):
    return mpf(k * (k + 1)) / 2 * (1 if z > 0 else (-1) ** (k + 1))
def ddP(k, z):
    return mpf((k - 1) * k * (k + 1) * (k + 2)) / 8 * (1 if z > 0 else (-1) ** k)

def side(i, z):
    """(xi, xi', xi'') on piece i at its left end (z=-1) or right end (z=+1), in y units."""
    a, b = pieces[i]; h = (b - a) / 2                    # dz/dy = 1/h
    v0 = v1 = v2 = mpf(0)
    for fi, k, w, s in by_piece[i]:
        c = w * s * pnorm[i][k] * xi[fi]
        v0 += c * P(k, z); v1 += c * dP(k, z) / h; v2 += c * ddP(k, z) / h ** 2
    return v0, v1, v2

v = side(0, -1)[0]

# prime powers up to lambda^2 and the EL weight Lambda(p^m) p^{-m/2}, enumerated by the
# same function as the matrix codes  [code audit B3]
weights = {q: lp / sqrt(mpf(q)) for q, lp in W.prime_powers_upto(lam ** 2)}

def label(y):
    for q, w in weights.items():
        if abs(y - log(mpf(q))) < 1e-6: return f"log {q}", w
        if abs(y - (L - log(mpf(q)))) < 1e-6: return f"L-log {q}", w
    return "-", None

print(f"lambda = {mp.nstr(lam, 6)}, L = {mp.nstr(L, 8)}, degree {deg}, {npc} pieces, "
      f"eps = {d['eps'][:12]}..., endpoint value v = xi(0) = {mp.nstr(v, 8)}")
print(f"EL weights Lambda(q) q^(-1/2): " + ", ".join(f"{q}: {mp.nstr(w, 5)}" for q, w in sorted(weights.items())))
print()
names = [("xi-", "xi+"), ("xi'-", "xi'+"), ("xi''-", "xi''+")]
hdr = f"{'cut':>8} {'at':>9} " + "".join(f" {a:>11} {b:>11} {'jump':>11} {'/wv':>9} |" for a, b in names)
print(hdr.rstrip(" |"))
for j in range(1, len(B) - 1):
    y = B[j]
    lv = side(j - 1, +1); rv = side(j, -1)
    name, w = label(y)
    cells = [f"{mp.nstr(y, 5):>8} {name:>9} "]
    for m in range(3):
        jump = rv[m] - lv[m]
        ratio = mp.nstr(jump / (w * v), 5) if w is not None else "-"
        cells.append(f" {mp.nstr(lv[m], 5):>11} {mp.nstr(rv[m], 5):>11} {mp.nstr(jump, 5):>11} {ratio:>9} |")
    print("".join(cells).rstrip(" |") )
print()
print("note: rows at L-log q are the mirror images of the rows at log q (even basis) and are not independent;")
print("      compare the xi-block '/wv' across degrees and across q (section 6 of the note).")
