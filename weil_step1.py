#!/usr/bin/env python3
"""
Step 1 of the numerical program: the truncated Weil quadratic form QW_lambda^N
following Connes-Consani-Moscovici, "Zeta Spectral Triples" (arXiv:2511.22755).

What it does, for given lambda and N:

  (A) Builds the (2N+1)x(2N+1) matrix tau of QW_lambda in the basis V_n, |n|<=N.
      tau_{nn} = a_n,  tau_{nm} = (b_n - b_m)/(n - m)          [Lemma 5.1]
      where  b_n = -(1/pi) D(sin(2 pi n y/L)),  a_n = 2 D((1-y/L) cos(2 pi n y/L))
      and D = W02 - WR - sum_p Wp is the distribution Psi^# pulled back to y=log x
      on [0, L], L = 2 log lambda                               [eqs (3.13)-(3.16)]

  (B) Diagonalizes tau: lowest eigenvalue eps_N, gap to the next, evenness of the
      eigenvector (xi_{-n} = xi_n), decay of its coefficients (convergence check).

  (C) Builds the perturbed operator D' = D - |D xi><eta|, computes its spectrum
      (which by Thm 5.10 is real and equals the zeros of the Fourier transform of xi),
      scales by 2 pi / L and compares with the zeros of zeta(1/2 + i s).

  (D) Builds the prolate "educated guess" k_lambda = E(h_lambda), h_lambda the
      zero-integral combination of prolate functions h_{0,lambda}, h_{4,lambda}
      [eq (7.6)], and computes the cosine similarity with xi.

  (E) Computes 1 - chi_4(lambda), the Fourier-truncation eigenvalue of h_{4,lambda},
      and compares log eps_N with log(1 - chi_4)                  [Figure 4]

Usage:
    python3 weil_step1.py --lam 3 --N 20 --dps 60
    python3 weil_step1.py --lam 2.5 --N 30 --dps 80 --no-prolate
    python3 weil_step1.py --lam 2.5 --N 60 --dps 70 --no-prolate --no-spectrum --dump xi.json
    python3 weil_step1.py --zeros xi.json --M 3000        # explicit-formula zero sum (section 7.1 of the note)
    python3 weil_step1.py --endpoint a.json,b.json        # endpoint diagnostic (section 7.2); prints D1, D2, D3
                                                          # with the minimizer's endpoint v_xi; the note's table
                                                          # uses the prolate endpoint v_k: D1_note = D1 - D2,
                                                          # D3_note = D3 + 2 D2.
Ablation flags --max-prime Q (keep prime powers <= Q) and --no-primes give forms that are NOT the
Weil form; the output says so.

Cost grows like N^3 in high-precision arithmetic; N=20..40 at 60-80 digits is
minutes, N=120 at 200 digits (the paper's setting) is hours.
"""

import argparse
import time
from mpmath import mp, mpf, pi, sin, cos, sinh, cosh, exp, expm1, log, sqrt
from mpmath import quad, matrix, eigsy, eig, euler, legendre, zetazero, zeta


# --------------------------------------------------------------------------
#  Number theory helpers
# --------------------------------------------------------------------------

def prime_powers_upto(x, rel_tol=mpf("1e-9")):
    """List of (k, Lambda(k)) for prime powers k <= x (von Mangoldt weights).

    x may be an mpf; the comparison k <= x is made with a relative tolerance so that
    lambda^2 passed to finite precision (e.g. --lam 1.41421356237 -> lambda^2 = 1.99999...)
    includes the prime power sitting exactly at lambda^2.  All diagnostic scripts use this
    function, so the prime-power set is enumerated in one place.  [code audit B3]
    """
    x = int(mp.floor(mpf(x) * (1 + rel_tol)))
    sieve = bytearray([1]) * (x + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(x ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
    out = []
    for p in range(2, x + 1):
        if sieve[p]:
            lp = log(mpf(p))
            q = p
            while q <= x:
                out.append((q, lp))
                q *= p
    return sorted(out)


# --------------------------------------------------------------------------
#  The distribution D on [0, L]   (D = W02 - WR - sum_p Wp)
# --------------------------------------------------------------------------

class WeilDistribution:
    def __init__(self, lam, nsub=8, max_prime=None, no_primes=False):
        self.lam = mpf(lam)
        self.L = 2 * log(self.lam)
        self.pp = prime_powers_upto(self.lam ** 2)      # p^m <= e^L = lambda^2
        # ablation options (the form is then NO LONGER the Weil form):
        if no_primes:
            self.pp = []
        elif max_prime is not None:
            self.pp = [(q, lp) for q, lp in self.pp if q <= max_prime]
        self.nsub = nsub                                 # base subdivisions for quad
        # constant in W_R^#: (1/2)(gamma + log 4 pi) - (1/2) log((e^L+1)/(e^L-1))
        eL = exp(self.L)
        self.cR = (euler + log(4 * pi)) / 2 - log((eL + 1) / (eL - 1)) / 2

    def _pts(self, freq_index):
        """Subdivision points on [0, L] adapted to the oscillation frequency."""
        k = self.nsub + 2 * abs(int(freq_index))
        return [self.L * i / k for i in range(k + 1)]

    # ---- W_{0,2}^#(phi) = int_0^L phi(y) (e^{y/2} + e^{-y/2}) dy
    def W02(self, phi, n):
        f = lambda y: phi(y) * 2 * cosh(y / 2)
        return quad(f, self._pts(n))

    # ---- W_R^#(phi) = cR*phi(0) + int_0^L (e^{y/2} phi(y) - phi(0)) / (2 sinh y) dy
    def WR(self, phi, n, phi0, phi_minus_phi0):
        # phi_minus_phi0(y) must compute phi(y)-phi(0) WITHOUT cancellation.
        def f(y):
            if y == 0:
                return mpf(0)   # never hit by tanh-sinh, but be safe
            num = expm1(y / 2) * phi(y) + phi_minus_phi0(y)
            return num / (2 * sinh(y))
        return self.cR * phi0 + quad(f, self._pts(n))

    # ---- sum_p W_p^#(phi) = sum_{p^m <= lambda^2} Lambda(p^m) p^{-m/2} phi(log p^m)
    def Wp(self, phi):
        return sum(lp * phi(log(mpf(q))) / sqrt(mpf(q)) for q, lp in self.pp)

    def apply(self, phi, n, phi0, phi_minus_phi0):
        return self.W02(phi, n) - self.WR(phi, n, phi0, phi_minus_phi0) - self.Wp(phi)

    # ---- b_n and a_n of Lemma 5.1
    def b(self, n):
        w = 2 * pi * n / self.L
        phi = lambda y: sin(w * y)
        return -self.apply(phi, n, mpf(0), phi) / pi

    def a(self, n):
        w = 2 * pi * n / self.L
        L = self.L
        phi = lambda y: (1 - y / L) * cos(w * y)
        # phi(y) - 1 = -2 (1 - y/L) sin^2(wy/2) - y/L   (no cancellation near y = 0)
        pm1 = lambda y: -2 * (1 - y / L) * sin(w * y / 2) ** 2 - y / L
        return 2 * self.apply(phi, n, mpf(1), pm1)


def build_tau(D, N, verbose=True):
    """The (2N+1)x(2N+1) matrix tau, indices -N..N mapped to 0..2N."""
    t0 = time.time()
    b = {}
    a = {}
    for n in range(N + 1):
        b[n] = D.b(n)
        a[n] = D.a(n)
        b[-n] = -b[n]
        a[-n] = a[n]
        if verbose and n % 5 == 0:
            print(f"    coefficients n={n:3d}   a={mp.nstr(a[n], 8):>14}  "
                  f"b={mp.nstr(b[n], 8):>14}   [{time.time()-t0:.0f}s]", flush=True)
    tau = matrix(2 * N + 1, 2 * N + 1)
    for i in range(-N, N + 1):
        for j in range(-N, N + 1):
            if i == j:
                tau[i + N, j + N] = a[i]
            else:
                tau[i + N, j + N] = (b[i] - b[j]) / (i - j)
    return tau, a, b


# --------------------------------------------------------------------------
#  Sanity check of W02 against the closed form (4.2) of the paper
# --------------------------------------------------------------------------

def check_W02_formula(D, n, m):
    L = D.L
    closed = 32 * L * sinh(L / 4) ** 2 * (L ** 2 - 16 * pi ** 2 * m * n) / (
        (L ** 2 + 16 * pi ** 2 * m ** 2) * (L ** 2 + 16 * pi ** 2 * n ** 2))
    if n == m:
        phi = lambda y: 2 * (1 - y / L) * cos(2 * pi * n * y / L)
    else:
        phi = lambda y: (sin(2 * pi * m * y / L) - sin(2 * pi * n * y / L)) / (pi * (n - m))
    numeric = D.W02(phi, max(abs(n), abs(m)))
    return closed, numeric


# --------------------------------------------------------------------------
#  Spectrum of the perturbed scaling operator  D' = D - |D xi><eta|
# --------------------------------------------------------------------------

def perturbed_spectrum(xi, N, L):
    """Eigenvalues of D' scaled by 2 pi / L.  xi is a list of 2N+1 real coefficients."""
    eta_xi = sum(xi)                     # <eta | xi>,  eta = sum_j e_j
    if abs(eta_xi) < mpf(10) ** (-mp.dps // 2):
        return None
    c = 1 / eta_xi
    Dp = matrix(2 * N + 1, 2 * N + 1)
    for i in range(-N, N + 1):
        for j in range(-N, N + 1):
            Dp[i + N, j + N] = (i if i == j else 0) - i * xi[i + N] * c
    ev = eig(Dp, left=False, right=False)
    ev = sorted(ev, key=lambda z: (mp.re(z)))
    return [z * 2 * pi / L for z in ev]


# --------------------------------------------------------------------------
#  Prolate spheroidal functions  (Legendre / Bouwkamp expansion, even parity)
# --------------------------------------------------------------------------

def prolate_even(c, K):
    """
    Even eigenfunctions of  -d/dz[(1-z^2) y'] + c^2 z^2 y = chi y  on [-1,1],
    expanded as y = sum_{k even} d_k P_k(z).  Returns (chis, coeff_lists) sorted
    by increasing chi, i.e. n = 0, 2, 4, ...   K = number of even Legendre terms.
    """
    ks = [2 * i for i in range(K)]
    B = matrix(K, K)
    c2 = mpf(c) ** 2
    for i, j in enumerate(ks):
        j = mpf(j)
        B[i, i] = j * (j + 1) + c2 * (2 * j * j + 2 * j - 1) / ((2 * j - 1) * (2 * j + 3))
        if i + 1 < K:
            off = c2 * (j + 1) * (j + 2) / ((2 * j + 3) * sqrt((2 * j + 1) * (2 * j + 5)))
            B[i, i + 1] = off
            B[i + 1, i] = off
    E, Q = eigsy(B)
    order = sorted(range(K), key=lambda i: E[i])
    out = []
    for idx in order:
        # back to raw Legendre coefficients: d_k = e_k * sqrt((2k+1)/2)
        d = [Q[i, idx] * sqrt((2 * ks[i] + 1) / mpf(2)) for i in range(K)]
        out.append((E[idx], ks, d))
    return out


def legendre_sum(ks, d, z):
    return sum(dk * legendre(k, z) for k, dk in zip(ks, d))


class ProlateGuess:
    """h_lambda = zero-integral combination of h_{0,lambda}, h_{4,lambda};  k_lambda = E(h_lambda)."""

    def __init__(self, lam, K=None):
        self.lam = mpf(lam)
        self.c = 2 * pi * self.lam ** 2                     # eq (7.9)
        if K is None:
            K = int(self.c) + 40
        # 1 - chi_n ~ e^{-2c}: need roughly 2c/ln(10) extra digits to resolve it
        self.work_dps = max(mp.dps, int(2 * self.c / 2.3) + 25)
        with mp.workdps(self.work_dps):
            modes = prolate_even(self.c, K)
        self.h0 = modes[0]      # n = 0
        self.h4 = modes[2]      # n = 4  (third even mode)
        # integral over [-lam, lam] of sum d_k P_k(x/lam) = 2 lam d_0
        I0 = 2 * self.lam * self.h0[2][0]
        I4 = 2 * self.lam * self.h4[2][0]
        self.alpha, self.beta = I4, -I0      # alpha*h0 + beta*h4 has zero integral

    def h(self, x):
        z = x / self.lam
        if abs(z) > 1:
            return mpf(0)
        return (self.alpha * legendre_sum(self.h0[1], self.h0[2], z)
                + self.beta * legendre_sum(self.h4[1], self.h4[2], z))

    def k(self, u):
        """E(h)(u) = u^{1/2} sum_{n>=1} h(n u),   u in [1/lam, lam]."""
        s = mpf(0)
        n = 1
        while n * u <= self.lam:
            s += self.h(n * u)
            n += 1
        return sqrt(u) * s

    def one_minus_chi4(self):
        """1 - chi_4(lambda):  hat h4(0) = chi_4 h4(0)  =>  chi_4 = int h4 / h4(0)."""
        ks, d = self.h4[1], self.h4[2]
        with mp.workdps(self.work_dps):
            integral = 2 * self.lam * d[0]
            at0 = legendre_sum(ks, d, mpf(0))
            return +(1 - integral / at0)


def xi_function(xi, N, L, lam):
    """The even eigenvector as a function of u in [1/lam, lam]:  sum xi_k V_k(u)."""
    def f(u):
        y = log(lam * u)
        s = xi[N]
        for k in range(1, N + 1):
            s += 2 * xi[k + N] * cos(2 * pi * k * y / L)     # uses xi_{-k} = xi_k
        return s / sqrt(L)
    return f


# --------------------------------------------------------------------------
#  Rayleigh quotient of the prolate guess in the Weil form's own norm
# --------------------------------------------------------------------------

def project_onto_basis(func_of_y, N, L, nodes=240, dps=None):
    """
    Real even coefficients c_n = L^{-1/2} int_0^L f(y) cos(2 pi n y / L) dy, n=0..N,
    for f given as a function of y = log(lam u) on [0, L].  f is symmetrized first
    (f(y) + f(L-y))/2 so that the coefficient vector is exactly even.
    Returns the full vector over -N..N.
    """
    from mpmath import gauss_quadrature
    with mp.workdps(dps or mp.dps):
        xs, ws = gauss_quadrature(nodes, "legendre")
        ys = [L * (x + 1) / 2 for x in xs]
        wts = [L * w / 2 for w in ws]
        fv = [(func_of_y(y) + func_of_y(L - y)) / 2 for y in ys]
        c = []
        for n in range(N + 1):
            c.append(sum(w * f * cos(2 * pi * n * y / L) for w, f, y in zip(wts, fv, ys)) / sqrt(L))
        return [c[abs(n)] for n in range(-N, N + 1)]


def quad_form(tau, v):
    n = len(v)
    return sum(v[i] * tau[i, j] * v[j] for i in range(n) for j in range(n))


# --------------------------------------------------------------------------
#  One complete analysis at fixed (lambda, N)
# --------------------------------------------------------------------------

NSUB = 8
MAX_PRIME = None
NO_PRIMES = False


def analyze(lam, N, nzeros=10, do_spectrum=True, do_prolate=True, verbose=True, dump=None):
    """Run steps (A)-(E).  Returns a dict of the key numbers."""
    R = {"lam": lam, "N": N, "dps": mp.dps}
    say = print if verbose else (lambda *a, **k: None)
    D = WeilDistribution(lam, nsub=NSUB, max_prime=MAX_PRIME, no_primes=NO_PRIMES)
    L = D.L
    say(f"lambda = {mp.nstr(lam, 10)},  L = 2 log lambda = {mp.nstr(L, 10)},  "
        f"N = {N},  {mp.dps} digits")
    say(f"prime powers <= lambda^2 = {mp.nstr(lam**2, 8)}: {[q for q, _ in D.pp]}"
        + ("   [ABLATION: primes removed/truncated - not the Weil form]" if (NO_PRIMES or MAX_PRIME is not None) else ""))
    say(f"matrix size (2N+1)x(2N+1) = {2*N+1}x{2*N+1}")

    if verbose:
        say("\n[check] W02 against closed form (4.2):")
        for (n, m) in [(0, 0), (1, 0), (2, 1)]:
            cf, nu = check_W02_formula(D, n, m)
            say(f"    (n,m)=({n},{m})  closed={mp.nstr(cf, 12)}  quad={mp.nstr(nu, 12)}  "
                f"diff={mp.nstr(cf - nu, 3)}")

    # ---- (A) matrix ----
    say("\n[A] building tau ...")
    tau, a, b = build_tau(D, N, verbose=verbose)

    # ---- (B) diagonalize ----
    say("\n[B] diagonalizing ...")
    E, Q = eigsy(tau)
    order = sorted(range(2 * N + 1), key=lambda i: E[i])
    e0, e1 = E[order[0]], E[order[1]]
    xi = [Q[i, order[0]] for i in range(2 * N + 1)]
    imax = max(range(2 * N + 1), key=lambda i: abs(xi[i]))
    if xi[imax] < 0:
        xi = [-x for x in xi]
    odd_part = sqrt(sum((xi[k + N] - xi[-k + N]) ** 2 for k in range(1, N + 1)) / 2)
    tail = max(abs(xi[0]), abs(xi[2 * N]))        # |xi_{+-N}|  (index 0 <-> n=-N)
    R.update(eps=e0, eps1=e1, gap=e1 - e0, odd=odd_part, tail=tail,
             log_eps=(log(e0) if e0 > 0 else None))
    say(f"    eps_N   (lowest eigenvalue)  = {mp.nstr(e0, 20)}")
    say(f"    next eigenvalue              = {mp.nstr(e1, 20)}")
    say(f"    gap                          = {mp.nstr(e1 - e0, 10)}")
    say(f"    log eps_N                    = "
        f"{mp.nstr(log(e0), 12) if e0 > 0 else 'n/a (eps_N <= 0)'}")
    say(f"    norm of odd part of xi       = {mp.nstr(odd_part, 5)}")
    # eigenvector accuracy is ~ 10^-dps / gap (mixing with the neighboring eigenvector);
    # the endpoint value we want to read is ~ e^{2.5} lam^4.5 e^{-2 pi lam^2} relative.
    import math as _m
    noise = mpf(10) ** (-mp.dps) / (e1 - e0)
    v_expect = _m.exp(2.5) * float(lam) ** 4.5 * _m.exp(-2 * _m.pi * float(lam) ** 2)
    need = int(_m.log10(1 / (v_expect * 1e-3 * float(e1 - e0)))) + 5
    say(f"    eigenvector noise floor ~ 10^-dps/gap = {mp.nstr(noise, 3)};  expected endpoint value ~ {v_expect:.1e} relative;"
        f"  digits needed to resolve it: ~{need}" + ("   <-- INSUFFICIENT PRECISION for endpoint/tail diagnostics" if float(noise) > 1e-3 * v_expect else ""))
    say("    coefficient decay |xi_k| (k = 0, 1, 2, ...):")
    say("      " + "  ".join(mp.nstr(abs(xi[k + N]), 4) for k in range(0, N + 1)))
    say("    eigenvalues (lowest 6): " + ", ".join(mp.nstr(E[order[i]], 8)
                                                  for i in range(min(6, 2 * N + 1))))
    if dump:
        import json
        with open(dump, "w") as fh:
            json.dump({"lam": mp.nstr(lam, mp.dps), "N": N, "dps": mp.dps, "L": mp.nstr(L, mp.dps),
                       "eps": mp.nstr(e0, mp.dps), "eps1": mp.nstr(e1, mp.dps),
                       "primes": [q for q, _ in D.pp],
                       "xi": [mp.nstr(x, mp.dps) for x in xi]}, fh)
        say(f"    eigenvector written to {dump}")

    # ---- (C) spectrum of D' vs zeta zeros ----
    if do_spectrum:
        say("\n[C] spectrum of the perturbed scaling operator vs zeros of zeta(1/2+is):")
        spec = perturbed_spectrum(xi, N, L)
        if spec is None:
            say("    <eta|xi> = 0, cannot normalize; skipping")
        else:
            max_im = max(abs(mp.im(z)) for z in spec)
            say(f"    max |Im| over spectrum (should be ~0): {mp.nstr(max_im, 5)}")
            pos = sorted([mp.re(z) for z in spec if mp.re(z) > mpf('0.5')])
            say("      k    zeta zero          eigenvalue          difference")
            errs = []
            for k in range(1, min(nzeros, len(pos)) + 1):
                zz = zetazero(k).imag
                errs.append(pos[k - 1] - zz)
                say(f"    {k:3d}   {mp.nstr(zz, 14):>16}   {mp.nstr(pos[k-1], 14):>16}   "
                    f"{mp.nstr(pos[k-1] - zz, 4):>10}")
            R["zero_err1"] = errs[0] if errs else None

    # ---- (D),(E) prolate comparison ----
    if do_prolate:
        say("\n[D] prolate guess k_lambda = E(h_lambda):")
        t0 = time.time()
        P = ProlateGuess(lam)
        say(f"    c = 2 pi lambda^2 = {mp.nstr(P.c, 8)},  Legendre terms = {len(P.h0[1])},  "
            f"working digits = {P.work_dps},  "
            f"chi_0 = {mp.nstr(P.h0[0], 10)}, chi_4 = {mp.nstr(P.h4[0], 10)}   "
            f"[{time.time()-t0:.0f}s]")
        xif = xi_function(xi, N, L, lam)

        # (D1) plain L^2 cosine similarity (20 digits is plenty)
        t0 = time.time()
        with mp.workdps(20):
            from mpmath import gauss_quadrature
            nodes, weights = gauss_quadrature(max(160, 6 * N + 60), "legendre")
            ys = [L * (x + 1) / 2 for x in nodes]
            ws = [L * wt / 2 for wt in weights]
            kv = [P.k(exp(y) / lam) for y in ys]
            xv = [xif(exp(y) / lam) for y in ys]
            kk = sum(w_ * k_ * k_ for w_, k_ in zip(ws, kv))
            xk = sum(w_ * k_ * x_ for w_, k_, x_ in zip(ws, kv, xv))
            xx = sum(w_ * x_ * x_ for w_, x_ in zip(ws, xv))
            cosine = xk / sqrt(kk * xx)
        R["one_minus_cos"] = 1 - abs(cosine)
        say(f"    cosine similarity <xi,k>/|xi||k| = {mp.nstr(cosine, 20)}   "
            f"[{time.time()-t0:.0f}s]")
        say(f"    1 - |cos|                    = {mp.nstr(1 - abs(cosine), 6)}")

        # (D2) Rayleigh quotient of k_lambda in the Weil form's own norm.
        #      Project k onto the V_n basis (working precision high enough that the
        #      projection error, squared, is below eps), then evaluate v^T tau v / v^T v.
        t0 = time.time()
        proj_dps = max(40, mp.dps)
        v = project_onto_basis(lambda y: P.k(exp(y) / lam), N, L,
                               nodes=max(240, 6 * N + 60), dps=proj_dps)
        vv = sum(x * x for x in v)
        v = [x / sqrt(vv) for x in v]
        rq = quad_form(tau, v)
        # residual of the projection: how much of ||k||^2 is outside the span of V_n, |n|<=N
        proj_loss = 1 - vv / kk
        # distance to xi in the Q-norm: Q(v - xi)  (both unit vectors, sign-aligned)
        sgn = 1 if sum(x * y for x, y in zip(v, xi)) > 0 else -1
        dvec = [x - sgn * y for x, y in zip(v, xi)]
        q_dist = quad_form(tau, dvec)
        R.update(rq_k=rq, proj_loss=proj_loss, q_dist=q_dist)
        say(f"    Rayleigh quotient Q(k)/|k|^2 = {mp.nstr(rq, 12)}   "
            f"(eps_N = {mp.nstr(e0, 6)},  ratio = {mp.nstr(rq / e0, 8) if e0 > 0 else 'n/a'})   "
            f"[{time.time()-t0:.0f}s]")
        say(f"    fraction of |k|^2 outside span(V_n,|n|<=N) = {mp.nstr(proj_loss, 5)}")
        say(f"    Q(k_hat - xi) (unit vectors)  = {mp.nstr(q_dist, 8)}")
        say("    values at a few u (xi normalized so that xi(lambda)=1, k likewise):")
        xl, kl = xif(lam), P.k(lam)
        for u in [1 / lam, sqrt(1 / lam), mpf(1), sqrt(lam), lam]:
            say(f"      u={mp.nstr(u, 6):>8}   xi={mp.nstr(xif(u)/xl, 10):>16}   "
                f"k={mp.nstr(P.k(u)/kl, 10):>16}")

        say("\n[E] eps_N versus 1 - chi_4(lambda)   (Figure 4 of the paper):")
        omc = P.one_minus_chi4()
        R["log_1mchi4"] = log(omc)
        say(f"    1 - chi_4(lambda)            = {mp.nstr(omc, 15)}")
        say(f"    log(1 - chi_4)               = {mp.nstr(log(omc), 12)}")
        say(f"    Fuchs exponent -4 pi lam^2 + 9 log lam = "
            f"{mp.nstr(-4*pi*lam**2 + 9*log(lam), 12)}   (slope reference only)")
        if e0 > 0:
            say(f"    log eps_N                    = {mp.nstr(log(e0), 12)}")
            say(f"    log eps_N - log(1-chi_4)     = {mp.nstr(log(e0) - log(omc), 12)}")
            say(f"    eps_N / (1 - chi_4)          = {mp.nstr(e0 / omc, 12)}")
        else:
            say("    eps_N <= 0 at this N; increase N (and dps) before comparing.")
    return R



# --------------------------------------------------------------------------
#  Explicit-formula decomposition:  QW(xi,xi) = sum_rho |xi_hat(t_rho)|^2
# --------------------------------------------------------------------------

def xi_hat_factory(xi, N, L):
    """Closed form (5.25): xi_hat(z) = 2 L^{-1/2} sin(zL/2) sum_j xi_j / (z - 2 pi j / L)."""
    js = [(j, xi[j + N]) for j in range(-N, N + 1)]
    def f(z):
        return 2 / sqrt(L) * sin(z * L / 2) * sum(c / (z - 2 * pi * j / L) for j, c in js)
    return f


def zeros_diagnostic(path, M=300, dps=40, mellin=True):
    import json
    d = json.load(open(path))
    mp.dps = dps
    lam, L, N = mpf(d["lam"]), mpf(d["L"]), d["N"]
    eps = mpf(d["eps"])
    xi = [mpf(x) for x in d["xi"]]
    xh = xi_hat_factory(xi, N, L)
    print(f"explicit-formula check: lambda={mp.nstr(lam,6)}, N={N}, eps_N={mp.nstr(eps,10)}")
    # sanity: Parseval,  ||xi||^2 = (1/2pi) int |xi_hat|^2 dt
    par = quad(lambda t: abs(xh(t)) ** 2, [-200, -50, -10, 0, 10, 50, 200]) / (2 * pi)
    print(f"  Parseval (1/2pi) int_-200^200 |xi_hat|^2 = {mp.nstr(par, 8)}   (||xi||^2 = 1; tail beyond 200 missing)")

    # Mellin transform of the prolate guess h_lambda, for the truncation-error comparison
    Mh = None
    if mellin:
        P = ProlateGuess(lam)
        # M(h)(1/2 - it) = int_0^lam v^{1/2 - it} h(v) dv/v   (h even; support [-lam, lam])
        def Mh(t):
            return quad(lambda v: v ** (mpf(1) / 2 - 1j * t) * P.h(v) / v,
                        [mpf(10) ** (-dps // 2)] + [lam * k / 8 for k in range(1, 9)])
        # normalization: match xi and E(h) at u = 1 (both ~ peak); use L2-optimal scale instead
        # scale c so that c * E(h) ~ xi in L2 on [1/lam, lam]
        with mp.workdps(20):
            from mpmath import gauss_quadrature
            xs, ws = gauss_quadrature(160, "legendre")
            ys = [L * (x + 1) / 2 for x in xs]; wt = [L * w / 2 for w in ws]
            xif = xi_function(xi, N, L, lam)
            kv = [P.k(exp(y) / lam) for y in ys]; xv = [xif(exp(y) / lam) for y in ys]
            c_scale = sum(w * a * b for w, a, b in zip(wt, kv, xv)) / sum(w * a * a for w, a in zip(wt, kv))
        print(f"  scale c with xi ~ c E(h_lambda): {mp.nstr(c_scale, 8)}")

    # xi_hat VANISHES at the zeros with O(1) slope, so an error dt in a zero gives a
    # spurious |xi_hat|^2 ~ dt^2.  Need dt << sqrt(eps): zero precision ~ -log10(sqrt eps) + 12.
    import math
    zdps = max(15, int(-math.log10(float(eps)) / 2) + 12)
    # xi_hat(t) is a sum of O(1e-2) terms canceling to ~sqrt(eps): evaluating it needs
    # about -log10(sqrt eps) digits of headroom on top of the digits wanted in the result.
    # [code audit B4]  Raise the working precision if --dps is below that.
    need = int(-math.log10(float(eps)) / 2) + 15
    if dps < need:
        print(f"  --dps {dps} is below the {need} digits needed for xi_hat at this eps; using {need}")
        mp.dps = dps = need
        xh = xi_hat_factory(xi, N, L)
    print(f"\n  computing {M} zeros at {zdps} digits (floor |xi_hat|^2 ~ 1e-{2*zdps} must be << eps),"
          f" evaluating xi_hat at {dps} digits ...", flush=True)
    t0 = time.time()
    with mp.workdps(zdps):
        ts = [mpf(zetazero(k).imag) for k in range(1, M + 1)]
    print(f"  zeros done [{time.time()-t0:.0f}s]", flush=True)
    vals = [abs(xh(mpf(t))) ** 2 for t in ts]
    S = mpf(0)
    partial = []
    for v in vals:
        S += v
        partial.append(S)
    print(f"\n  partial sums, counting each zero pair rho, rho-bar:   2 S_M = 2 sum_{{k<=M}} |xi_hat(t_k)|^2")
    print("     M       t_M        2 S_M / eps_N")
    for k in [1, 2, 3, 5, 10, 20, 50, 100, 200, 300, 500, 1000, 1500, 2000, 3000, 5000, M]:
        if k <= M:
            print(f"   {k:5d}  {mp.nstr(ts[k-1], 8):>10}   {mp.nstr(2 * partial[k-1] / eps, 8):>10}")
    # contribution by window in t
    print("\n  contribution to eps by frequency window (2 x sum over zeros in window):")
    edges = [0, 30, 50, 77, 110, 150, 200, 300, 400, 600, 800, 1200, 1600, 2400, 3200, 1e9]
    for lo, hi in zip(edges[:-1], edges[1:]):
        idx = [i for i, t in enumerate(ts) if lo <= t < hi]
        if not idx:
            continue
        c = 2 * sum(vals[i] for i in idx)
        avg = sum(vals[i] for i in idx) / len(idx)
        print(f"     t in [{lo:5.0f},{min(hi, float(ts[-1])):5.0f})  {len(idx):5d} zeros   "
              f"{mp.nstr(c / eps, 6):>10} of eps   mean|xi_hat|^2 = {mp.nstr(avg, 4)}")
    # tail: fit window means of |xi_hat|^2 over the upper half of the range to A t^-p,
    # then integrate against the zero density (1/2pi) log(t/2pi)
    import math
    half = M // 2
    wins = [(ts[i], vals[i]) for i in range(half, M)]
    nb = 8
    xs_, ys_ = [], []
    for b in range(nb):
        chunk = wins[b * len(wins) // nb:(b + 1) * len(wins) // nb]
        xs_.append(math.log(float(sum(t for t, _ in chunk) / len(chunk))))
        ys_.append(math.log(float(sum(v for _, v in chunk) / len(chunk))))
    n = len(xs_)
    mx, my = sum(xs_) / n, sum(ys_) / n
    p = -sum((x - mx) * (y - my) for x, y in zip(xs_, ys_)) / sum((x - mx) ** 2 for x in xs_)
    A = math.exp(my + p * mx)
    T = float(ts[-1])
    if p < 1.2:
        print(f"\n  tail fit gave <|xi_hat|^2> ~ {A:.3e} t^(-{p:.2f}): exponent too small for a convergent tail;"
              " fit is unreliable (noise or not yet in the 1/t^2 regime).  Using p=2 with the same"
              " last-window amplitude as a rough estimate.")
        last = wins[-len(wins) // nb:]
        tm = float(sum(t for t, _ in last) / len(last)); vm = float(sum(v for _, v in last) / len(last))
        A, p = vm * tm ** 2, 2.0
    tail = quad(lambda t: A * t ** (-p) * log(t / (2 * pi)) / (2 * pi), [T, 10 * T, 100 * T, mp.inf])
    print(f"\n  tail beyond t_M = {T:.1f}:  window means fit  <|xi_hat|^2> ~ {A:.3e} t^(-{p:.2f});"
          f"  2 x integral against zero density = {mp.nstr(2 * tail / eps, 5)} of eps")
    print(f"  endpoint value implied by the tail (jump picture, 2 v^2 = A):  v = {math.sqrt(A/2):.4e}"
          f"  (divide by the peak value xi(1) to compare with --endpoint)")
    total = 2 * partial[-1] + 2 * tail
    print(f"\n  2 S_M + tail = {mp.nstr(total / eps, 8)} eps_N      (deviation from 1: {mp.nstr(total/eps - 1, 4)})")
    if mellin:
        print("\n  between zeros: |xi_hat(t)| versus c |zeta(1/2 - it) M(h)(1/2 - it)|  (should be ~1 away from zeros)")
        print("       t        |xi_hat|        c|zeta M(h)|      ratio")
        for t in [mpf(x) for x in ['2', '5', '8', '12', '17.5', '23', '28', '35']]:
            a = abs(xh(t)); b = abs(c_scale * zeta(mpf(1) / 2 - 1j * t) * Mh(t))
            print(f"   {mp.nstr(t,5):>7}   {mp.nstr(a,6):>12}   {mp.nstr(b,6):>12}   {mp.nstr(a/b,6):>8}")


# --------------------------------------------------------------------------
#  Endpoint diagnostic: v_xi = xi(1/lam) versus v_k = E(h_lam)(1/lam) and sqrt(1-chi_4)
# --------------------------------------------------------------------------

def endpoint_diagnostic(paths, dps=40):
    import json
    rows = []
    print("endpoint values.  Predictions if the mechanism is right: the three log-differences are lambda-independent.")
    print("  lam    N   log eps    log(v_xi/pk)  log(v_k/pk)  .5log(1-chi4)   D1=log(v_xi/pk)-.5log(1-chi4)   D2=log(v_xi/v_k)   D3=log eps-2log(v_xi/pk)")
    for path in paths:
        d = json.load(open(path))
        mp.dps = max(dps, d["dps"])
        lam, L, N = mpf(d["lam"]), mpf(d["L"]), d["N"]
        eps = mpf(d["eps"])
        xi = [mpf(x) for x in d["xi"]]
        xk = xi[N:]
        # xi(y) = L^{-1/2} (xi_0 + 2 sum_{k>=1} xi_k cos(2 pi k y / L));  y=0 -> endpoint, y=L/2 -> u=1
        def xi_at(y, Kmax):
            return (xk[0] + 2 * sum(xk[k] * cos(2 * pi * k * y / L) for k in range(1, Kmax + 1))) / sqrt(L)
        v_full = xi_at(mpf(0), N)
        v_half = xi_at(mpf(0), N // 2)
        pk = xi_at(L / 2, N)
        P = ProlateGuess(lam)
        v_k = P.k(1 / lam)
        pk_k = P.k(mpf(1))
        omc = P.one_minus_chi4()
        r_xi, r_k = v_full / pk, v_k / pk_k
        D1 = log(abs(r_xi)) - log(omc) / 2
        D2 = log(abs(r_xi / r_k))
        D3 = log(eps) - 2 * log(abs(r_xi))
        rows.append((lam, N, D1, D2, D3))
        print(f"  {mp.nstr(lam,4):>4} {N:4d}  {mp.nstr(log(eps),7):>8}   {mp.nstr(log(abs(r_xi)),7):>10}   "
              f"{mp.nstr(log(abs(r_k)),7):>10}   {mp.nstr(log(omc)/2,7):>10}      {mp.nstr(D1,5):>8}"
              f"                     {mp.nstr(D2,5):>8}           {mp.nstr(D3,5):>8}")
        print(f"        v_xi/pk from N terms = {mp.nstr(r_xi,6)},  from N/2 terms = {mp.nstr(v_half/pk,6)}"
              f"   (agreement gauges truncation error of the endpoint value);  sign(v_xi)={'+' if v_full>0 else '-'}, sign(v_k)={'+' if v_k>0 else '-'}")
    return rows

# --------------------------------------------------------------------------
#  Sweep driver
# --------------------------------------------------------------------------

def fmt(x, d=6):
    return "n/a" if x is None else mp.nstr(x, d)


def sweep(lams, Ns, dps_rule, nzeros, csv_path):
    """Run analyze() over every (lambda, N) pair; print a table and write CSV."""
    cols = ["lam", "N", "dps", "eps", "gap", "odd", "tail", "log_eps",
            "log_1mchi4", "one_minus_cos", "rq_k", "proj_loss", "q_dist", "zero_err1"]
    rows = []
    print("   lam    N  dps       eps_N     gap/eps    |xi_N|   log eps  log(1-chi4)  1-|cos|    Q(k)/eps   Q(k-xi)    zero1 err")
    for lam in lams:
        for N in Ns:
            mp.dps = dps_rule(lam, N)
            t0 = time.time()
            R = analyze(mpf(lam), N, nzeros=nzeros, verbose=False)
            rows.append(R)
            ge = (R["gap"] / R["eps"]) if R["eps"] > 0 else None
            rk = (R["rq_k"] / R["eps"]) if R["eps"] > 0 else None
            print(f"{mp.nstr(mpf(lam),5):>6} {N:4d} {mp.dps:4d} {fmt(R['eps'],4):>11} {fmt(ge,3):>9} "
                  f"{fmt(R['tail'],3):>9} {fmt(R['log_eps'],7):>9} {fmt(R.get('log_1mchi4'),7):>12} "
                  f"{fmt(R.get('one_minus_cos'),3):>9} {fmt(rk,3):>10} {fmt(R.get('q_dist'),3):>9} "
                  f"{fmt(R.get('zero_err1'),3):>10}   [{time.time()-t0:.0f}s]", flush=True)
            with open(csv_path, "w") as fh:
                fh.write(",".join(cols) + "\n")
                for r in rows:
                    fh.write(",".join(fmt(r.get(c), 20) if c not in ("N", "dps") else str(r[c])
                                      for c in cols) + "\n")
    print(f"\nwritten {csv_path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lam", type=str, default="3", help="single lambda (ignored if --sweep-lam)")
    ap.add_argument("--N", type=int, default=20, help="single N (ignored if --sweep-N)")
    ap.add_argument("--dps", type=int, default=60,
                    help="working digits. eps_N needs ~ -log10(eps)+20; the endpoint/tail diagnostics need "
                         "rounding << (endpoint value) x (eigenvalue gap): roughly 4.5*lam^2 + 60 (lam=4.2: ~140)")
    ap.add_argument("--nzeros", type=int, default=10)
    ap.add_argument("--no-prolate", action="store_true")
    ap.add_argument("--no-spectrum", action="store_true")
    ap.add_argument("--sweep-lam", type=str, default=None,
                    help="comma list of lambdas, e.g. '1.4142,2,2.5,3,3.5'")
    ap.add_argument("--sweep-N", type=str, default=None,
                    help="comma list of N values, e.g. '10,20,30,40'")
    ap.add_argument("--auto-dps", action="store_true",
                    help="choose digits from lambda: max(--dps, 4.5*lambda^2 + 30), enough to resolve eps")
    ap.add_argument("--csv", type=str, default="sweep.csv")
    ap.add_argument("--nsub", type=int, default=8, help="base quadrature subdivisions (default 8)")
    ap.add_argument("--no-primes", action="store_true", help="ablation: drop all prime terms (archimedean + pole only)")
    ap.add_argument("--max-prime", type=int, default=None, help="ablation: keep only prime powers <= this")
    ap.add_argument("--dump", type=str, default=None, help="write eigenvector xi to this JSON file")
    ap.add_argument("--zeros", type=str, default=None, help="explicit-formula check on a --dump file")
    ap.add_argument("--endpoint", type=str, default=None, help="comma list of --dump files: endpoint-value diagnostic")
    ap.add_argument("--M", type=int, default=300, help="number of zeta zeros for --zeros")
    ap.add_argument("--no-mellin", action="store_true", help="skip the zeta*M(h) comparison in --zeros")
    args = ap.parse_args()
    global NSUB, MAX_PRIME, NO_PRIMES
    NSUB = args.nsub
    MAX_PRIME = args.max_prime
    NO_PRIMES = args.no_primes
    mp.dps = args.dps          # before any mpf() is formed from the command line  [code audit B7]

    if args.zeros:
        zeros_diagnostic(args.zeros, M=args.M, dps=args.dps, mellin=not args.no_mellin)
        return
    if args.endpoint:
        endpoint_diagnostic(args.endpoint.split(","), dps=args.dps)
        return

    if args.sweep_lam or args.sweep_N:
        lams = [mpf(x) for x in args.sweep_lam.split(",")] if args.sweep_lam else [mpf(args.lam)]
        Ns = [int(x) for x in args.sweep_N.split(",")] if args.sweep_N else [args.N]
        if args.auto_dps:
            # eps ~ exp(-~10 lambda^2): need ~4.5 lambda^2 digits plus headroom
            rule = lambda lam, N: max(args.dps, int(4.5 * float(lam) ** 2) + 30)
        else:
            rule = lambda lam, N: args.dps
        sweep(lams, Ns, rule, args.nzeros, args.csv)
        return

    mp.dps = args.dps
    analyze(mpf(args.lam), args.N, nzeros=args.nzeros,
            do_spectrum=not args.no_spectrum, do_prolate=not args.no_prolate, verbose=True,
            dump=args.dump)


if __name__ == "__main__":
    main()
