"""Per-piece Legendre coefficient decay from a weil_breakpoints.py --dump file.

usage: python3 bp_piece_decay.py bp_lam22_d24.json [bp_lam22_d20.json ...]

For each piece prints |c_k|, k = 0..d, where c_k is the eigenvector's coefficient
on the orthonormal Legendre polynomial of degree k on that piece (paired mirror
functions carry weight 1/sqrt2 each; the center piece has even k only), and the
least-squares slope of log|c_k| against log k over the upper half of the degrees,
even and odd k fitted separately (they can differ when one parity is suppressed by
symmetry).  Rough guide for a singularity at a piece END: step -> coefficients
reach the noise floor (exactly representable); log|y-y0| -> slope about -2;
(y-y0)log|y-y0| -> about -3; |y-y0|^a -> about -(2a+1).  Slopes need k >= ~8
to mean anything; with d = 24 the fit uses 6-7 points per parity, so read them as
indicative.  Several dumps can be given to see whether the slopes are stable in d.
"""
import json, sys, math

def analyze(path):
    d = json.load(open(path))
    deg = d["d"]; L = float(d["L"]); lam = float(d["lam"])
    pieces = [(float(a), float(b)) for a, b in d["pieces"]]
    xi = [float(x) for x in d["xi"]]
    funcs = d["funcs"]
    coef = [dict() for _ in pieces]
    for fi, parts in enumerate(funcs):
        for pi, k, w in parts:
            coef[pi][k] = abs(xi[fi] * float(w))
    print(f"\n=== {path}: lambda = {lam}, degree {deg}, {len(pieces)} pieces, eps = {d['eps'][:10]}")
    half = (len(pieces) + 1) // 2
    for pi in range(half):            # left half (plus the center piece if there is one); mirrors are identical
        a, b = pieces[pi]
        ks = sorted(coef[pi])
        print(f"\npiece {pi}: [{a:.4f}, {b:.4f}]  width {b-a:.4f}")
        print("   k   |c_k|")
        for k in ks:
            print(f"  {k:2d}   {coef[pi][k]:.3e}")
        for parity, name in ((0, "even"), (1, "odd")):
            pts = [(k, coef[pi][k]) for k in ks if k % 2 == parity and k >= deg // 2 and coef[pi][k] > 0]
            if len(pts) < 3:
                continue
            xs = [math.log(k) for k, c in pts]; ys = [math.log(c) for k, c in pts]
            n = len(xs); mx = sum(xs) / n; my = sum(ys) / n
            sxx = sum((x - mx) ** 2 for x in xs); sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
            slope = sxy / sxx
            print(f"   log-log slope, {name} k in [{pts[0][0]}, {pts[-1][0]}] ({n} pts): {slope:+.2f}")

for p in sys.argv[1:]:
    analyze(p)
