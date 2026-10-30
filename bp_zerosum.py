#!/usr/bin/env python3
"""
bp_zerosum.py -- basis-free check of a breakpoint-basis eigenvalue.

Reads a --dump file written by weil_breakpoints.py, reconstructs the eigenvector
xi(y) as a polynomial on each piece, computes its Fourier transform in closed form,
and evaluates the explicit-formula zero sum

        QW(xi, xi) = sum_rho |xi_hat(t_rho)|^2 = 2 * sum_{k>=1} |xi_hat(t_k)|^2

over the first M zeros of zeta.  If the breakpoint eigenvalue eps_d is genuine,
2 S_M + tail must come out equal to eps_d.  No matrix of either basis is used.

    python3 bp_zerosum.py bp_lam2_d20.json --M 3000 --dps 40

Fourier transform of a Legendre polynomial on a piece [a, b], h = b - a:
    int_a^b P_k(z(y)) e^{-i t y} dy = (h/2) e^{-i t (a + h/2)} * 2 (-i)^k j_k(t h / 2),
with j_k the spherical Bessel function, j_k(w) = sqrt(pi / (2 w)) J_{k+1/2}(w).
"""
import argparse, json, math, time
from mpmath import mp, mpf, mpc, pi, sqrt, exp, log, besselj, quad, zetazero, legendre


def load(path):
    d = json.load(open(path))
    mp.dps = max(mp.dps, d["dps"])
    D = {
        "lam": mpf(d["lam"]), "L": mpf(d["L"]), "d": d["d"],
        "eps": mpf(d["eps"]),
        "pieces": [(mpf(a), mpf(b)) for a, b in d["pieces"]],
        "xi": [mpf(x) for x in d["xi"]],
        "funcs": [[(int(pi_), int(k), mpf(w)) for pi_, k, w in parts] for parts in d["funcs"]],
        "pnorm": [[mpf(x) for x in row] for row in d["pnorm"]],
    }
    return D


def piece_coefficients(D):
    """c[i][k] = coefficient of P_k(z) on piece i (z the local coordinate in [-1, 1]),
    reproducing BreakpointWeil.evals exactly: weight * sign * pnorm * xi."""
    npc = len(D["pieces"]); dmax = D["d"]
    c = [[mpf(0)] * (dmax + 1) for _ in range(npc)]
    for fi, parts in enumerate(D["funcs"]):
        for (pi_, k, wgt) in parts:
            sign = 1
            if len(parts) == 2 and parts[1][0] == pi_:      # mirror copy of a paired function
                sign = (-1) ** k
            c[pi_][k] += wgt * sign * D["pnorm"][pi_][k] * D["xi"][fi]
    return c


def xi_of_y(D, c):
    """xi(y) evaluated from the piece polynomials (for checks)."""
    def f(y):
        for i, (a, b) in enumerate(D["pieces"]):
            if a <= y <= b:
                z = 2 * (y - a) / (b - a) - 1
                return sum(c[i][k] * legendre(k, z) for k in range(D["d"] + 1))
        return mpf(0)
    return f


def sph_bessel_all(kmax, w):
    """[j_0(w), ..., j_kmax(w)] via besselj of half-integer order; j_k(-w) = (-1)^k j_k(w)
    so negative arguments are handled by parity rather than by principal branches.  [code audit B6]"""
    if w == 0:
        return [mpf(1)] + [mpf(0)] * kmax
    aw = abs(w)
    pref = sqrt(pi / (2 * aw))
    js = [pref * besselj(k + mpf(1) / 2, aw) for k in range(kmax + 1)]
    if w < 0:
        js = [j if k % 2 == 0 else -j for k, j in enumerate(js)]
    return js


def xi_hat_factory(D, c):
    """xi_hat(t) = int_0^L xi(y) e^{-i t y} dy, in closed form."""
    pieces, dmax = D["pieces"], D["d"]
    ipow = [(-1j) ** k for k in range(dmax + 1)]
    def f(t):
        t = mpf(t)
        total = mpc(0)
        for i, (a, b) in enumerate(pieces):
            h = b - a
            w = t * h / 2
            js = sph_bessel_all(dmax, w)
            s = sum(c[i][k] * ipow[k] * js[k] for k in range(dmax + 1))
            total += (h / 2) * exp(-1j * t * (a + h / 2)) * 2 * s
        return total
    return f


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dump")
    ap.add_argument("--M", type=int, default=3000, help="number of zeta zeros")
    ap.add_argument("--dps", type=int, default=40, help="digits for xi_hat evaluation")
    ap.add_argument("--zeros-cache", type=str, default=None, help="JSON file to cache/load zeros")
    args = ap.parse_args()
    mp.dps = args.dps
    D = load(args.dump)
    c = piece_coefficients(D)
    L, eps = D["L"], D["eps"]
    xh = xi_hat_factory(D, c)
    xif = xi_of_y(D, c)
    print(f"breakpoint zero-sum check: lambda={mp.nstr(D['lam'],6)}, L={mp.nstr(L,8)}, degree {D['d']}, "
          f"{len(D['pieces'])} pieces, dim {len(D['xi'])}, eps_d={mp.nstr(eps,12)}")

    # ---- sanity 1: ||xi||^2 = 1 from the piece polynomials (basis is orthonormal)
    nrm = sum(quad(lambda y: xif(y) ** 2, [a, b]) for (a, b) in D["pieces"])
    print(f"  ||xi||^2 from piece polynomials      = {mp.nstr(nrm, 15)}   (expect 1)")
    # ---- sanity 2: closed-form transform vs direct quadrature at one t
    t0 = mpf(14.134725141734693790457251983562)
    direct = sum(quad(lambda y: xif(y) * exp(-1j * t0 * y), [a, b]) for (a, b) in D["pieces"])
    print(f"  |xi_hat(t_1)|: closed form {mp.nstr(abs(xh(t0)), 12)}   quadrature {mp.nstr(abs(direct), 12)}")
    # ---- sanity 3: Parseval,  (1/2 pi) int |xi_hat|^2 dt = 1  (finite range; tail ~ 1/t)
    par = quad(lambda t: abs(xh(t)) ** 2, [-400, -100, -30, 0, 30, 100, 400]) / (2 * pi)
    print(f"  Parseval (1/2pi) int_-400^400 |xi_hat|^2 = {mp.nstr(par, 8)}   (expect ~1; tail beyond 400 missing)")

    # ---- zeros: known to ~12 digits beyond sqrt(eps), since xi_hat vanishes at a zero with
    #      O(1) slope and an error dt in the zero gives a spurious |xi_hat|^2 ~ dt^2
    zdps = max(15, int(-math.log10(float(eps)) / 2) + 12)
    ts = None
    if args.zeros_cache:
        try:
            zc = json.load(open(args.zeros_cache))
            if zc["dps"] >= zdps and len(zc["t"]) >= args.M:
                ts = [mpf(x) for x in zc["t"][:args.M]]
                print(f"  loaded {args.M} zeros from {args.zeros_cache}")
        except FileNotFoundError:
            pass
    if ts is None:
        print(f"  computing {args.M} zeros at {zdps} digits ...", flush=True)
        t1 = time.time()
        with mp.workdps(zdps):
            ts = [mpf(zetazero(k).imag) for k in range(1, args.M + 1)]
        print(f"  zeros done [{time.time()-t1:.0f}s]", flush=True)
        if args.zeros_cache:
            with mp.workdps(zdps):
                json.dump({"dps": zdps, "t": [mp.nstr(t, zdps) for t in ts]}, open(args.zeros_cache, "w"))

    # ---- the sum
    t1 = time.time()
    vals = [abs(xh(t)) ** 2 for t in ts]
    print(f"  xi_hat evaluated at {len(ts)} zeros [{time.time()-t1:.0f}s]")
    S = mpf(0); partial = []
    for v in vals:
        S += v; partial.append(S)
    print("\n  partial sums  2 S_M / eps_d:")
    print("     M       t_M        2 S_M / eps_d")
    for k in [1, 2, 3, 5, 10, 20, 50, 100, 200, 300, 500, 1000, 1500, 2000, 3000, 5000, args.M]:
        if k <= args.M:
            print(f"   {k:5d}  {mp.nstr(ts[k-1], 8):>10}   {mp.nstr(2 * partial[k-1] / eps, 10):>12}")

    print("\n  contribution by frequency window (2 x sum over zeros in window):")
    edges = [0, 30, 50, 77, 110, 150, 200, 300, 400, 600, 800, 1200, 1600, 2400, 3200, 1e9]
    for lo, hi in zip(edges[:-1], edges[1:]):
        idx = [i for i, t in enumerate(ts) if lo <= t < hi]
        if not idx:
            continue
        cc = 2 * sum(vals[i] for i in idx)
        print(f"     t in [{lo:5.0f},{min(hi, float(ts[-1])):5.0f})  {len(idx):5d} zeros   {mp.nstr(cc / eps, 6):>10} of eps_d")

    # ---- tail beyond t_M: window means of |xi_hat|^2 ~ A t^-p over the upper half, integrated
    #      against the zero density (1/2 pi) log(t / 2 pi)
    half = args.M // 2
    wins = [(ts[i], vals[i]) for i in range(half, args.M)]
    nb = 8
    xs_, ys_ = [], []
    for b in range(nb):
        chunk = wins[b * len(wins) // nb:(b + 1) * len(wins) // nb]
        xs_.append(math.log(float(sum(t for t, _ in chunk) / len(chunk))))
        ys_.append(math.log(float(sum(v for _, v in chunk) / len(chunk))))
    n = len(xs_); mx, my = sum(xs_) / n, sum(ys_) / n
    p = -sum((x - mx) * (y - my) for x, y in zip(xs_, ys_)) / sum((x - mx) ** 2 for x in xs_)
    A = math.exp(my + p * mx)
    T = float(ts[-1])
    if p < 1.2:
        last = wins[-len(wins) // nb:]
        tm = float(sum(t for t, _ in last) / len(last)); vm = float(sum(v for _, v in last) / len(last))
        A, p = vm * tm ** 2, 2.0
        print(f"\n  tail exponent fit unreliable; using p = 2 with last-window amplitude")
    tail = quad(lambda t: A * t ** (-p) * log(t / (2 * pi)) / (2 * pi), [T, 10 * T, 100 * T, mp.inf])
    print(f"\n  tail beyond t_M = {T:.1f}:  <|xi_hat|^2> ~ {A:.3e} t^(-{p:.2f});  2 x tail = {mp.nstr(2 * tail / eps, 5)} of eps_d")
    total = 2 * partial[-1] + 2 * tail
    print(f"\n  2 S_M + tail = {mp.nstr(total / eps, 10)} x eps_d      (deviation from 1: {mp.nstr(total / eps - 1, 4)})")
    print(f"  endpoint value implied by the 1/t^2 tail (2 v^2 = A): v = {math.sqrt(A / 2):.4e};  xi(0) from the polynomials = {mp.nstr(xif(mpf(0)), 6)}")


if __name__ == "__main__":
    main()
