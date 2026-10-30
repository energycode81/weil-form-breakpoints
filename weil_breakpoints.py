#!/usr/bin/env python3
"""
weil_breakpoints.py -- the truncated Weil form in a piecewise-polynomial basis with
breakpoints at y = log p^m and L - log p^m (prime powers p^m <= lambda^2), even under
y -> L - y.  Section 6 of the note: the minimizer is non-smooth at these positions, and a
basis that places its cuts there reaches a lower Galerkin value at a given dimension than
the Fourier or plain Legendre bases, or than the same cuts displaced (--no-breaks with
--extra-breaks).

    python3 weil_breakpoints.py --lam 2.2 --d 6,10,14 --dps 40

Basis: on each piece I, orthonormal Legendre polynomials of degree <= d; mirror pieces are
paired into even combinations; a piece symmetric about L/2 uses even degrees only.
q(f,g)(y) = int_0^{L-y} [f(x) g(x+y) + f(x+y) g(x)] dx is integrated exactly (Gauss-Legendre on
the sub-pieces cut by the breakpoints of f and of the shifted g).  The y-integrals of D are
split at all pairwise differences of breakpoints, where q is non-smooth.
"""
import argparse, time, json, bisect
from mpmath import mp, mpf, sqrt, exp, expm1, sinh, cosh, log, matrix, eigsy
from mpmath import gauss_quadrature
import weil_step1 as W
from weil_legendre import leg_all


class BreakpointWeil:
    def __init__(self, lam, d, ygl=None, max_prime=None, no_primes=False, extra_breaks=()):
        # y-integrand contains q(y), a polynomial of degree 2d+1: the Gauss rule must resolve it
        if ygl is None:
            ygl = max(24, 2 * d + 8)
        self.lam = mpf(lam)
        self.D = W.WeilDistribution(lam, max_prime=max_prime, no_primes=no_primes)
        L = self.L = self.D.L
        self.d = d
        # ---- breakpoints, symmetric under y -> L - y.  Points closer than tol are merged,
        # so that e.g. log 2 and L/2 at lambda = 2 (equal in exact arithmetic) cannot produce
        # a sliver piece through a last-bit difference.  [code audit B2]
        tol = L * mpf(10) ** (-(mp.dps // 2))
        B = [mpf(0), L]
        def add(v):
            v = mpf(v)
            if v < 0 or v > L:
                return
            if all(abs(v - b) > tol for b in B):
                B.append(v)
        for q, _ in self.D.pp:
            add(log(mpf(q))); add(L - log(mpf(q)))
        for b in extra_breaks:
            add(b); add(L - mpf(b))
        self.B = sorted(B)
        self.pieces = [(self.B[i], self.B[i + 1]) for i in range(len(self.B) - 1)]
        # ---- even basis: list of (piece index, degree, partner piece index or None, sign)
        # function index -> list of (piece, degree, weight)
        self.funcs = []
        npc = len(self.pieces)
        done = set()
        for i, (a, b) in enumerate(self.pieces):
            if i in done:
                continue
            mirror = npc - 1 - i          # B symmetric => pieces mirror in reverse order
            if mirror == i:
                for k in range(0, d + 1, 2):
                    self.funcs.append([(i, k, mpf(1))])
            else:
                for k in range(0, d + 1):
                    self.funcs.append([(i, k, 1 / sqrt(mpf(2))), (mirror, k, 1 / sqrt(mpf(2)))])
                done.add(mirror)
            done.add(i)
        self.n = len(self.funcs)
        # map piece -> list of (function index, degree, weight); whether an entry is the mirror
        # copy of a paired function is decided in evals() from funcs[fi]
        self.by_piece = [[] for _ in range(npc)]
        for fi, parts in enumerate(self.funcs):
            for (pi, k, wgt) in parts:
                self.by_piece[pi].append((fi, k, wgt))
        self.pnorm = [[sqrt((2 * k + 1) / (b - a)) for k in range(d + 1)] for (a, b) in self.pieces]
        self.gx = gauss_quadrature(d + 2, "legendre")       # exact for degree 2d+2
        self.gy = gauss_quadrature(ygl, "legendre")

    # basis values at a point: sparse list [(func index, value)]
    def evals(self, x):
        L = self.L
        if x < 0 or x > L:
            return []
        i = bisect.bisect_right(self.B, x) - 1
        if i >= len(self.pieces):
            i = len(self.pieces) - 1
        a, b = self.pieces[i]
        # mirrored piece j > center uses P_k(L - y) mapped: for a paired function the mirror
        # piece's polynomial is p_{I,k}(L - y); since P_k(-z) = (-1)^k P_k(z) that equals
        # (-1)^k p on the mirror piece in its own local coordinate.  Handle by local coord + sign.
        z = 2 * (x - a) / (b - a) - 1
        P = leg_all(self.d, z)
        out = []
        for (fi, k, wgt) in self.by_piece[i]:
            parts = self.funcs[fi]
            sign = 1
            if len(parts) == 2 and parts[1][0] == i:       # this is the mirror copy
                sign = (-1) ** k
            out.append((fi, wgt * sign * self.pnorm[i][k] * P[k]))
        return out

    def q_matrix(self, y):
        """Q[i][j] = q(phi_i, phi_j)(y)."""
        L, n = self.L, self.n
        Q = [[mpf(0)] * n for _ in range(n)]
        ell = L - y
        if ell <= 0:
            return Q
        cuts = {mpf(0), ell}
        for b in self.B:
            if 0 < b < ell: cuts.add(b)
            if 0 < b - y < ell: cuts.add(b - y)
        cuts = sorted(cuts)
        zn, wn = self.gx
        for a, b in zip(cuts[:-1], cuts[1:]):
            if b - a <= 0:
                continue
            for z, w in zip(zn, wn):
                x = a + (b - a) * (z + 1) / 2
                wt = (b - a) * w / 2
                A = self.evals(x); Bv = self.evals(x + y)
                for (i, ai) in A:
                    for (j, bj) in Bv:
                        v = wt * ai * bj
                        Q[i][j] += v
                        Q[j][i] += v
        return Q

    def build(self, verbose=True):
        t0 = time.time()
        n, L = self.n, self.L
        # y-cuts: breakpoints and all pairwise differences
        ycuts = {mpf(0), L}
        for b in self.B:
            ycuts.add(b)
            for c in self.B:
                dlt = abs(b - c)
                if 0 < dlt < L: ycuts.add(dlt)
        ycuts = sorted(ycuts)
        T = [[mpf(0)] * n for _ in range(n)]
        Q0 = self.q_matrix(mpf(0))
        cR = self.D.cR
        zn, wn = self.gy
        nseg = len(ycuts) - 1
        for s, (a, b) in enumerate(zip(ycuts[:-1], ycuts[1:])):
            if b - a <= 0:
                continue
            for z, w in zip(zn, wn):
                y = a + (b - a) * (z + 1) / 2
                wt = (b - a) * w / 2
                Q = self.q_matrix(y)
                w02 = wt * 2 * cosh(y / 2)
                e = expm1(y / 2); den = 2 * sinh(y)
                for i in range(n):
                    Ti, Qi, Q0i = T[i], Q[i], Q0[i]
                    for j in range(i, n):
                        Ti[j] += w02 * Qi[j] - wt * (e * Qi[j] + (Qi[j] - Q0i[j])) / den
            if verbose and s % 5 == 0:
                print(f"    y-segment {s}/{nseg}  [{time.time()-t0:.0f}s]", flush=True)
        for i in range(n):
            for j in range(i, n):
                T[i][j] -= cR * Q0[i][j]
        for q, lp in self.D.pp:
            Q = self.q_matrix(log(mpf(q)))
            wgt = lp / sqrt(mpf(q))
            for i in range(n):
                for j in range(i, n):
                    T[i][j] -= wgt * Q[i][j]
        Tm = matrix(n, n)
        for i in range(n):
            for j in range(i, n):
                Tm[i, j] = T[i][j]; Tm[j, i] = T[i][j]
        self.T = Tm
        if verbose:
            print(f"    matrix built ({n}x{n}, {len(self.pieces)} pieces, degree {self.d}) [{time.time()-t0:.0f}s]", flush=True)
        return Tm

    def function(self, coeffs):
        def f(y):
            return sum(c * coeffs[fi] for (fi, c) in self.evals(y))
        return f


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lam", type=str, default="2.2")
    ap.add_argument("--d", type=str, default="6,10,14")
    ap.add_argument("--dps", type=int, default=40)
    ap.add_argument("--ygl", type=int, default=None,
                    help="Gauss points per y-segment (default max(24, 2d+8)); the runs in the note used 24 or 32, "
                         "checked against the default to 15-16 digits at d = 20 and d = 32")
    ap.add_argument("--no-breaks", action="store_true",
                    help="control: drop the prime-position breakpoints (the form keeps every prime); "
                         "with --extra-breaks this gives a basis with cuts at the extra positions only")
    ap.add_argument("--fourier-dump", type=str, default=None)
    ap.add_argument("--extra-breaks", type=str, default="", help="comma list of extra breakpoints y (mirrors added)")
    ap.add_argument("--dump", type=str, default=None, help="write the lowest eigenvector (piece polynomials) to this JSON file")
    args = ap.parse_args()
    mp.dps = args.dps
    lam = mpf(args.lam)
    print(f"lambda = {mp.nstr(lam, 8)}, L = {mp.nstr(2*log(lam), 8)}, {mp.dps} digits, breakpoint basis"
          + ("  [CONTROL: prime-position breakpoints removed; the form keeps all primes]" if args.no_breaks else ""))
    rows = []
    for d in [int(x) for x in args.d.split(",")]:
        xb = [mpf(x) for x in args.extra_breaks.split(",") if x]
        if args.no_breaks:
            BW = BreakpointWeil(lam, d, ygl=args.ygl, no_primes=True, extra_breaks=xb)
            BW.D = W.WeilDistribution(lam)        # keep the primes in D, drop only the breakpoints
        else:
            BW = BreakpointWeil(lam, d, ygl=args.ygl, extra_breaks=xb)
        print(f"\n[d = {d}]  breakpoints: {[mp.nstr(b, 5) for b in BW.B]}")
        T = BW.build()
        E, Q = eigsy(T)
        o = sorted(range(BW.n), key=lambda i: E[i])
        e0, e1 = E[o[0]], E[o[1]]
        xi = [Q[i, o[0]] for i in range(BW.n)]
        if max(xi, key=abs) < 0:                  # sign convention: largest coefficient positive  [code audit B5]
            xi = [-x for x in xi]
        f = BW.function(xi)
        v, pk = f(mpf(0)), f(BW.L / 2)
        print(f"    eps_d = {mp.nstr(e0, 18)}   next = {mp.nstr(e1, 8)}   log eps_d = {mp.nstr(log(e0), 12) if e0 > 0 else 'n/a'}")
        print(f"    endpoint/peak = {mp.nstr(v/pk, 8)}  (log {mp.nstr(log(abs(v/pk)), 8)})")
        # per-piece coefficient magnitude of the top two degrees: shows convergence per piece
        tops = []
        for pi in range(len(BW.pieces)):
            vals = {k: abs(xi[fi]) for (fi, k, wgt) in BW.by_piece[pi]}
            ks = sorted(vals)
            tops.append((mp.nstr(vals[ks[-2]], 2) if len(ks) > 1 else '-', mp.nstr(vals[ks[-1]], 2)))
        print("    |coeff| of top two degrees per piece: " + "  ".join(f"({a},{b})" for a, b in tops))
        if args.fourier_dump:
            dd = json.load(open(args.fourier_dump)); N = dd["N"]; Lf = mpf(dd["L"]); lamf = mpf(dd["lam"])
            xf = [mpf(x) for x in dd["xi"]]
            g = W.xi_function(xf, N, Lf, lamf)
            # project g onto the basis by quadrature on each piece
            c = [mpf(0)] * BW.n; nrm = mpf(0)
            zn, wn = gauss_quadrature(4 * d + 20, "legendre")
            for (a, b) in BW.pieces:
                for z, w in zip(zn, wn):
                    y = a + (b - a) * (z + 1) / 2; wt = (b - a) * w / 2
                    gy = g(exp(y) / lamf); nrm += wt * gy * gy
                    for (fi, val) in BW.evals(y):
                        c[fi] += wt * gy * val
            val = sum(c[i] * T[i, j] * c[j] for i in range(BW.n) for j in range(BW.n))
            print(f"    Fourier eigenvector projected: Q = {mp.nstr(val, 10)}, norm captured {mp.nstr(sum(x*x for x in c)/nrm, 10)}")
        if args.dump:
            S = lambda x: mp.nstr(x, mp.dps)
            json.dump({"lam": S(lam), "L": S(BW.L), "dps": mp.dps, "d": d,
                       "eps": S(e0), "eps1": S(e1),
                       "B": [S(b) for b in BW.B],
                       "pieces": [[S(a), S(b)] for (a, b) in BW.pieces],
                       "funcs": [[[pi_, k, S(w)] for (pi_, k, w) in parts] for parts in BW.funcs],
                       "pnorm": [[S(x) for x in row] for row in BW.pnorm],
                       "xi": [S(x) for x in xi]}, open(args.dump, "w"))
            print(f"    eigenvector written to {args.dump}")
        rows.append((d, BW.n, e0, v / pk))
    print("\nsummary:   d    dim   eps_d                  log eps_d     log(v/peak)")
    for d, n, e0, r in rows:
        print(f"        {d:3d}   {n:4d}   {mp.nstr(e0, 15):>20}   {mp.nstr(log(e0), 10) if e0 > 0 else 'n/a':>12}   {mp.nstr(log(abs(r)), 8):>10}")


if __name__ == "__main__":
    main()
