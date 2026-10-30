"""Window test for the Legendre tail residual.

usage: python3 leg_window.py leg_lam25_M120.json 100 [shift]

For the tail (degrees >= K0) of the minimizer in the JSON dump, reports the
largest |resid| inside a window of +-0.02 around every prime-power position
log p^m and L - log p^m, and, as controls, around the midpoints between
consecutive positions and around the two endpoints.  An optional third
argument shifts every prime-power position by that amount (the wrong-position
control): with shift=0.05 the "prime" windows should look like the midpoint
windows if the bursts are real.

Also writes the full residual on a 4001-point grid to <stem>_resid_K<K0>.csv
so it can be plotted.
"""
import json, sys, math
from mpmath import mp, mpf
import weil_legendre as WL
import weil_step1 as W

mp.dps = 30
path, K0 = sys.argv[1], int(sys.argv[2])
shift = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0

d = json.load(open(path))
lam = float(d["lam"]); L = mpf(d["L"]); M = d["M"]; ks = d["ks"]
c = [mpf(x) for x in d["c"]]
norm = [mp.sqrt(mpf(2 * k + 1) / L) for k in ks]
Lf = float(L)

ng = 4000
ys = [Lf * i / ng for i in range(ng + 1)]
res = []
for y in ys:
    P = WL.leg_all(M, 2 * mpf(y) / L - 1)
    res.append(float(abs(sum(c[i] * norm[i] * P[k] for i, k in enumerate(ks) if k >= K0))))

# prime powers up to lambda^2, from the same enumeration as the matrix codes  [code audit B3]
pp = sorted(q for q, _ in W.prime_powers_upto(mpf(lam) ** 2))

pos = {}
for q in pp:
    pos[f"log {q}"] = math.log(q) + shift
    pos[f"L-log {q}"] = Lf - math.log(q) - shift

w = 0.02
def wmax(y0):
    return max(r for y, r in zip(ys, res) if abs(y - y0) <= w)

print(f"lambda={lam}, M={M}, K0={K0}, shift={shift:+.3f}, window +-{w}")
print("prime-power positions:")
for name, y0 in sorted(pos.items(), key=lambda kv: kv[1]):
    print(f"  {name:>9} {y0:7.4f}   max|resid| = {wmax(y0):.3e}")

ps = sorted(pos.values())
print("midpoint controls:  (a midpoint closer than 2w = %.2f to a prime position shares grid points with its window and is marked *)" % (2 * w))
for a, b in zip(ps[:-1], ps[1:]):
    y0 = (a + b) / 2
    flag = "*" if min(abs(y0 - p) for p in ps) < 2 * w else " "
    print(f"  {'mid':>9} {y0:7.4f}   max|resid| = {wmax(y0):.3e} {flag}")

print("endpoint windows:")
for y0 in (w, Lf - w):
    print(f"  {'end':>9} {y0:7.4f}   max|resid| = {wmax(y0):.3e}")

out = path.replace(".json", f"_resid_K{K0}.csv")
with open(out, "w") as f:
    for y, r in zip(ys, res):
        f.write(f"{y:.6f},{r:.6e}\n")
print("full residual written to", out)
