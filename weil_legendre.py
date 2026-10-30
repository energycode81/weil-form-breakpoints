#!/usr/bin/env python3
"""
weil_legendre.py -- the truncated Weil quadratic form in a basis of shifted Legendre
polynomials on [0, L], L = 2 log lambda, restricted to functions even under y -> L - y.

This is an independent check of weil_step1.py: same distribution D (imported), but
  * a polynomial basis instead of exp(2 pi i n y / L),
  * q(f, g)(y) = int_0^{L-y} [f(x) g(x+y) + f(x+y) g(x)] dx computed exactly by Gauss-Legendre
    in x (the integrand is a polynomial),
  * D applied to q by fixed high-order Gauss-Legendre in y (no oscillation, no subtraction
    problems: q is smooth and q(y) - q(0) is formed at nodes bounded away from 0).

The Hilbert space is L^2([0, L], dy), so eigenvalues are directly comparable with the
Fourier-basis eps_N.  Usage:

    python3 weil_legendre.py --lam 2.5 --M 30,40,50,60 --dps 80

M is the maximal polynomial degree; only even degrees are used (even functions), so the
matrix is (M/2 + 1) x (M/2 + 1).
"""
import argparse, time, json
from mpmath import mp, mpf, sqrt, exp, expm1, sinh, cosh, log, pi, cos, matrix, eigsy, legendre
from mpmath import gauss_quadrature
import weil_step1 as W


def leg_all(kmax, z):
    """P_0(z) .. P_kmax(z) by recurrence."""
    out = [mpf(1), z]
    for k in range(1, kmax):
        out.append(((2 * k + 1) * z * out[k] - k * out[k - 1]) / (k + 1))
    return out[:kmax + 1]


class LegendreWeil:
    def __init__(self, lam, M, ny_sub=12, ny_deg=40, max_prime=None, no_primes=False):
        self.lam = mpf(lam)
        self.D = W.WeilDistribution(lam, max_prime=max_prime, no_primes=no_primes)
        self.L = self.D.L
        self.M = M
        self.ks = list(range(0, M + 1, 2))          # even degrees only
        self.n = len(self.ks)
        L = self.L
        # orthonormal shifted Legendre: phi_k(y) = sqrt((2k+1)/L) P_k(2y/L - 1)
        self.norm = [sqrt(mpf(2 * k + 1) / L) for k in self.ks]

        # --- y-grid for D: Gauss-Legendre on ny_sub subintervals of [0, L]
        zn, wn = gauss_quadrature(ny_deg, "legendre")
        self.ys, self.wy = [], []
        for s in range(ny_sub):
            a, b = L * s / ny_sub, L * (s + 1) / ny_sub
            for z, w in zip(zn, wn):
                self.ys.append(a + (b - a) * (z + 1) / 2)
                self.wy.append((b - a) * w / 2)
        # --- x-grid: exact for polynomials of degree <= 2M+1
        self.zx, self.wx = gauss_quadrature(M + 2, "legendre")
        # --- prime points
        self.pts_primes = [(log(mpf(q)), lp / sqrt(mpf(q))) for q, lp in self.D.pp]

    def phi_matrix(self, xs):
        """Matrix Phi[m][i] = phi_{k_i}(xs[m]) for the even basis."""
        L = self.L
        out = []
        for x in xs:
            P = leg_all(self.M, 2 * x / L - 1)
            out.append([self.norm[i] * P[k] for i, k in enumerate(self.ks)])
        return out

    def q_matrix(self, y):
        """Q[i][j] = q(phi_i, phi_j)(y) for all pairs, y in [0, L]."""
        L, n = self.L, self.n
        ell = L - y
        xs = [ell * (z + 1) / 2 for z in self.zx]
        ws = [ell * w / 2 for w in self.wx]
        A = self.phi_matrix(xs)                       # phi(x)
        B = self.phi_matrix([x + y for x in xs])      # phi(x + y)
        Q = [[mpf(0)] * n for _ in range(n)]
        for m in range(len(xs)):
            w, a, b = ws[m], A[m], B[m]
            for i in range(n):
                wai, wbi = w * a[i], w * b[i]
                row = Q[i]
                for j in range(n):
                    row[j] += wai * b[j] + wbi * a[j]
        return Q

    def build(self, verbose=True):
        t0 = time.time()
        n = self.n
        L = self.L
        T = [[mpf(0)] * n for _ in range(n)]
        Q0 = self.q_matrix(mpf(0))                   # q(0) = 2 <phi_i, phi_j> = 2 delta_ij
        cR = self.D.cR
        # accumulate  W02 - WR  over the y-grid
        for idx, (y, w) in enumerate(zip(self.ys, self.wy)):
            Q = self.q_matrix(y)
            w02 = w * 2 * cosh(y / 2)
            e = expm1(y / 2)
            den = 2 * sinh(y)
            for i in range(n):
                for j in range(i, n):
                    qv = Q[i][j]
                    T[i][j] += w02 * qv - w * (e * qv + (qv - Q0[i][j])) / den
            if verbose and idx % 120 == 0:
                print(f"    y-node {idx}/{len(self.ys)}  [{time.time()-t0:.0f}s]", flush=True)
        # constant part of W_R:  -cR * q(0)
        for i in range(n):
            for j in range(i, n):
                T[i][j] -= cR * Q0[i][j]
        # primes:  - sum Lambda(k) k^{-1/2} q(log k)
        for (ly, wgt) in self.pts_primes:
            Q = self.q_matrix(ly)
            for i in range(n):
                for j in range(i, n):
                    T[i][j] -= wgt * Q[i][j]
        Tm = matrix(n, n)
        for i in range(n):
            for j in range(i, n):
                Tm[i, j] = T[i][j]
                Tm[j, i] = T[i][j]
        self.T = Tm
        if verbose:
            print(f"    matrix built ({n}x{n}) [{time.time()-t0:.0f}s]", flush=True)
        return Tm

    def quad_form_of_function(self, f):
        """Rayleigh quotient of a function f(y) on [0, L] projected onto the basis."""
        zn, wn = gauss_quadrature(4 * self.M + 40, "legendre")
        c = [mpf(0)] * self.n
        nrm = mpf(0)
        for z, w in zip(zn, wn):
            y = self.L * (z + 1) / 2
            wy = self.L * w / 2
            fy = f(y)
            nrm += wy * fy * fy
            P = leg_all(self.M, z)
            for i, k in enumerate(self.ks):
                c[i] += wy * fy * self.norm[i] * P[k]
        val = sum(c[i] * self.T[i, j] * c[j] for i in range(self.n) for j in range(self.n))
        proj = sum(x * x for x in c)
        return val, proj, nrm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lam", type=str, default="2.5")
    ap.add_argument("--M", type=str, default="20,30,40")
    ap.add_argument("--dps", type=int, default=80)
    ap.add_argument("--no-primes", action="store_true")
    ap.add_argument("--max-prime", type=int, default=None)
    ap.add_argument("--fourier-dump", type=str, default=None,
                    help="optional weil_step1 --dump file: cross-check the Fourier eigenvector's Rayleigh quotient")
    ap.add_argument("--dump", type=str, default=None,
                    help="write the minimizer's coefficients to NAME_M{M}.json (the _M{M} suffix is added)")
    args = ap.parse_args()
    mp.dps = args.dps
    lam = mpf(args.lam)
    Ms = [int(x) for x in args.M.split(",")]
    print(f"lambda = {mp.nstr(lam, 8)}, L = {mp.nstr(2*log(lam), 8)}, {mp.dps} digits, even Legendre basis")
    results = []
    for M in Ms:
        print(f"\n[M = {M}]  even degrees 0..{M}, matrix {M//2+1}x{M//2+1}")
        LW = LegendreWeil(lam, M, max_prime=args.max_prime, no_primes=args.no_primes)
        print(f"    prime powers: {[q for q,_ in LW.D.pp]}")
        T = LW.build()
        E, Q = eigsy(T)
        order = sorted(range(LW.n), key=lambda i: E[i])
        e0, e1 = E[order[0]], E[order[1]]
        xi = [Q[i, order[0]] for i in range(LW.n)]
        if max(xi, key=abs) < 0:
            xi = [-x for x in xi]
        print(f"    eps_M = {mp.nstr(e0, 20)}   next = {mp.nstr(e1, 10)}   gap/eps = {mp.nstr(e1/e0 - 1, 5) if e0 != 0 else 'n/a'}")
        print(f"    log eps_M = {mp.nstr(log(e0), 12) if e0 > 0 else 'n/a'}")
        print("    Legendre coefficients |c_k| of the minimizer, k = 0, 2, 4, ...:")
        print("      " + "  ".join(mp.nstr(abs(x), 3) for x in xi))
        # endpoint value and peak:  xi(0) = sum c_i phi_i(0),  xi(L/2) = sum c_i phi_i(L/2)
        P0 = leg_all(M, mpf(-1)); Pm = leg_all(M, mpf(0))
        v = sum(xi[i] * LW.norm[i] * P0[k] for i, k in enumerate(LW.ks))
        pk = sum(xi[i] * LW.norm[i] * Pm[k] for i, k in enumerate(LW.ks))
        print(f"    endpoint xi(0) = {mp.nstr(v, 8)},  xi(L/2) = {mp.nstr(pk, 8)},  v/peak = {mp.nstr(v/pk, 8)}"
              f"  (log = {mp.nstr(log(abs(v/pk)), 8)})")
        # cross-check: quadratic form of cos(2 pi y / L)  vs Fourier a_1 + b_1 (times L/2)
        w1 = 2 * pi / LW.L
        val, proj, nrm = LW.quad_form_of_function(lambda y: cos(w1 * y))
        a1, b1 = LW.D.a(1), LW.D.b(1)
        print(f"    check Q(cos 2pi y/L): Legendre {mp.nstr(val, 15)}  vs Fourier (L/2)(a1+b1) = {mp.nstr((a1+b1)*LW.L/2, 15)}"
              f"   (projection captures {mp.nstr(proj/nrm, 10)} of the norm)")
        if args.fourier_dump:
            d = json.load(open(args.fourier_dump)); N = d["N"]; Lf = mpf(d["L"]); lamf = mpf(d["lam"])
            xf = [mpf(x) for x in d["xi"]]
            f = W.xi_function(xf, N, Lf, lamf)
            val, proj, nrm = LW.quad_form_of_function(lambda y: f(exp(y) / lamf))
            print(f"    Fourier eigenvector (N={N}, eps_N={d['eps'][:10]}) in this basis: Q = {mp.nstr(val, 10)},"
                  f"  norm captured = {mp.nstr(proj/nrm, 12)}")
        results.append((M, e0, e1, v / pk))
        if args.dump:
            with open(args.dump.replace(".json", f"_M{M}.json"), "w") as fh:
                json.dump({"lam": mp.nstr(lam, mp.dps), "M": M, "L": mp.nstr(LW.L, mp.dps),
                           "eps": mp.nstr(e0, mp.dps), "ks": LW.ks, "c": [mp.nstr(x, mp.dps) for x in xi]}, fh)
    print("\nsummary:  M    eps_M                    log eps_M      log(v/peak)")
    for M, e0, e1, r in results:
        print(f"        {M:3d}   {mp.nstr(e0, 15):>22}   {mp.nstr(log(e0), 10) if e0 > 0 else 'n/a':>12}   {mp.nstr(log(abs(r)), 8):>10}")


if __name__ == "__main__":
    main()
