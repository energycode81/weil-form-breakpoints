# Time-to-accuracy: Fourier vs breakpoint basis

*Every row of the raw-data table has its log in `logs/ladders/` (`tta_<basis>_lam<λ>_<size>.log`) and its wall time in `logs/ladders/tta.csv`; the bracketed `[…s]` stamps inside the logs are matrix-assembly time only and are 10–25% below the wall times, which include the eigensolve. Several rows were run a second time, with identical ε: the three Fourier N = 100 rows (`logs/fourier/lam2_N100.log`, `lam25_N100.log`, `lam3_N100.log`; times 279/310/423 s), the λ = 2.5 d = 20 (ygl 32), d = 24 and d = 32 rows (`logs/breakpoint/lam25_correct_d20.log`, `_d24`, `_d32`), the λ = 3 plain d = 24 row (`logs/breakpoint/lam3_d24_plain.log`, 1139 s on a loaded machine vs 876 s here), and both λ = 3 extra-cut rows (`logs/breakpoint/controls/lam3_d24_predextras.log`, `lam3_d24_wrongextras.log`).*

Same `--dps` for both bases at each λ (50 / 60 / 90). Fourier runs use `--no-prolate --no-spectrum` so both sides time matrix assembly + eigensolve only. Fourier dimension is 2N+1; breakpoint dimension is the actual matrix size printed by `weil_breakpoints.py` (`--ygl 24`). Every value is a Galerkin upper bound on ε_λ, so the smaller of two values is the closer.

## Reference values

| λ | reference ε_λ | source | spread of last two breakpoint entries |
|---|---|---|---|
| 2 | 7.5187e-13 | breakpoint d = 48 (earlier run, ygl scaled) | 1.3e-4 (d = 28 → 48) |
| 2.5 | 4.7414e-24 | breakpoint d = 32 (dump bp_lam25_d32.json); ladder 4.7532, 4.7462, 4.7433, 4.7414 at d = 16, 20, 24, 32; d = 20 identical to 15 digits under --ygl 32 | 4e-4 (d = 24 → 32), still moving |
| 3 | 2.470e-38 | breakpoint d = 24 + three extra cuts in the slowest piece (dim 225; predicted-position and arbitrary-position sets give 2.4702e-38 and 2.4699e-38); d = 28 gives 2.477e-38 | 3e-3 (upper bound; still drifting slowly) |

At λ = 2 the Fourier 1/N² extrapolation (7.48e-13) agrees with the breakpoint reference to 0.5%. At λ = 2.5 the Fourier extrapolations (4.96–5.29e-24, depending on assumed exponent) are all *above* the breakpoint value, and the Fourier local convergence exponent is erratic (0.3–1.6), so the Fourier extrapolation is unreliable there.

## Table: first ladder entry reaching a given relative accuracy

| λ | accuracy | Fourier (dim, s) | breakpoint (dim, s) | speedup |
|---|---|---|---|---|
| 2 | 10% | N = 60 (121, 123 s) | d = 10 (22, 3.1 s) | 40× |
| 2 | 1% | not reached at N = 100 (201, 362 s; 2.6%) | d = 10 (22, 3.1 s) | > 100× |
| 2 | 0.1% (3 digits) | projected N ≈ 400 (≈ 800, ≈ 1.7 h at t ∝ N²) | d = 12 (26, 4.8 s) | ≈ 1300× (projected) |
| 2.5 | 10% | not reached at N = 100 (201, 385 s; 14%) | d = 16 (77, 84 s) | > 4.6× |
| 2.5 | 1% | out of reach | d = 16 (77, 84 s) | — |
| 2.5 | 0.1% (3 digits) | out of reach | d = 20 (95, 146 s) | — |
| 3 | 10% | not reached at N = 100 (201, 516 s; 27% high) | d = 24 (150, 875 s; 0.7% high) | > 0.6× (Fourier never gets there) |
| 3 | 1% | out of reach | d = 24 (150, 875 s) | — |
| 3 | 0.3% | out of reach | d = 28 (174, 1746 s) or d = 24 + three extra cuts (225, 4598 s) | — |

The λ = 2 "0.1% Fourier" row is a projection from the observed N^(−2.3) tail. The λ = 3 rows are measured (runs of 2026-09-12/13). The Fourier 1/N² extrapolation at λ = 3 (2.97e-38) was 20% high, so Fourier ladders cannot be extrapolated reliably there.

## Cost scaling (measured)

- Fourier: t ∝ N^2.0 at all three λ (55.7 s → 362 s for N = 40 → 100 at λ = 2). Assembly-dominated: O(N) archimedean quadratures, O(N²) multiprecision entry arithmetic.
- Breakpoint: t ∝ d^2.45–2.55 at fixed λ (2.48, 2.45, 2.55 from the d = 10 → 20 / 12 → 24 pairs at λ = 2, 2.5, 3); at fixed d = 12 the cost is 4.8 s / 41 s / 149 s for λ = 2 / 2.5 / 3 (dimension 26 / 59 / 78), i.e. it grows with the number of pieces (which grows with the number of prime powers ≤ λ²).

## Convergence character

Fourier: at every λ, a fast pre-asymptotic phase followed by an algebraic tail. At λ = 3 the relative error drops 15 orders of magnitude from N = 10 to N = 40 (local exponents 31, 26, 8.9), then 1.99, 1.88, 0.59 on the steps 40→60→80→100 — slow, with an unstable exponent. This is the reconciliation with Groskin's "Sobolev exponent" s(c): a three-point power-law fit over N ∈ {40, 60, 80} measures the transient when the transient has not yet ended at those N, which is the case at his c = 23; at c = λ² = 9 the transient ends near N = 40 and the tail exponent is ≈ 2, not 46.

Breakpoint: below a λ-dependent threshold degree the error is enormous (10^5 at λ = 2, d = 4; 10^14 at λ = 2.5, d = 4), then it drops by 0.3–2.4 orders of magnitude per unit degree. Threshold degree (first d within a factor ~2 of the converged value): 8–10 at λ = 2, 12–16 at λ = 2.5, 20–24 at λ = 3, roughly 2–2.5·λ². Once the prime-position non-smoothness is resolved, the remaining error is the 10⁻³·v non-analytic component described in Section 2 of the note, and the ladders are reported as data without a rate claim.

## Caveats

1. Resolved: λ = 2.5, d = 20 rerun with `--ygl 32` gives 4.74617992974983e-24, identical to the `--ygl 24` value to 15 digits.
2. The λ = 2.5 breakpoint ladder stalls between d = 10 and d = 12 (factor 4 instead of the ~10³ per two degrees seen elsewhere) before resuming. Unexplained; could be a quadrature or a parity-pairing effect. Worth one look.
3. The last digits at λ = 2 converge slowly from d = 16 on — the 10⁻³·v non-analytic component at the endpoints and prime positions (Section 2 of the note), which no polynomial basis resolves spectrally.

## Raw data

| λ | basis | size | dim | seconds | ε |
|---|---|---|---|---|---|
| 2 | fourier | N=10 | 21 | 5.3 | 1.9496e-12 |
| 2 | fourier | N=15 | 31 | 9.7 | 1.0676e-12 |
| 2 | fourier | N=20 | 41 | 15.8 | 1.0320e-12 |
| 2 | fourier | N=30 | 61 | 32.5 | 9.6369e-13 |
| 2 | fourier | N=40 | 81 | 55.7 | 8.8828e-13 |
| 2 | fourier | N=60 | 121 | 122.9 | 8.1281e-13 |
| 2 | fourier | N=80 | 161 | 222.6 | 7.8417e-13 |
| 2 | fourier | N=100 | 201 | 362.1 | 7.7130e-13 |
| 2 | breakpoint | d=4 | 10 | 0.6 | 5.8897e-08 |
| 2 | breakpoint | d=6 | 14 | 1.1 | 1.2133e-09 |
| 2 | breakpoint | d=8 | 18 | 1.9 | 4.0613e-12 |
| 2 | breakpoint | d=10 | 22 | 3.1 | 7.5802e-13 |
| 2 | breakpoint | d=12 | 26 | 4.8 | 7.5247e-13 |
| 2 | breakpoint | d=16 | 34 | 9.6 | 7.5228e-13 |
| 2 | breakpoint | d=20 | 42 | 17.3 | 7.5212e-13 |
| 2.5 | fourier | N=10 | 21 | 5.8 | 1.0027e-19 |
| 2.5 | fourier | N=15 | 31 | 10.8 | 1.9488e-22 |
| 2.5 | fourier | N=20 | 41 | 17.7 | 1.0823e-23 |
| 2.5 | fourier | N=30 | 61 | 36.3 | 6.7935e-24 |
| 2.5 | fourier | N=40 | 81 | 62.5 | 6.3703e-24 |
| 2.5 | fourier | N=60 | 121 | 136.0 | 5.6001e-24 |
| 2.5 | fourier | N=80 | 161 | 244.4 | 5.5226e-24 |
| 2.5 | fourier | N=100 | 201 | 384.9 | 5.4099e-24 |
| 2.5 | breakpoint | d=4 | 23 | 4.1 | 7.6228e-10 |
| 2.5 | breakpoint | d=6 | 32 | 8.8 | 1.0787e-14 |
| 2.5 | breakpoint | d=8 | 41 | 16.1 | 6.0379e-18 |
| 2.5 | breakpoint | d=10 | 50 | 26.8 | 6.4642e-21 |
| 2.5 | breakpoint | d=12 | 59 | 41.4 | 1.5119e-21 |
| 2.5 | breakpoint | d=16 | 77 | 84.1 | 4.7532e-24 |
| 2.5 | breakpoint | d=20 | 95 | 146.4 | 4.7462e-24 |
| 3 | fourier | N=10 | 21 | 8.0 | 2.9719e-23 |
| 3 | fourier | N=20 | 41 | 23.6 | 1.2662e-32 |
| 3 | fourier | N=30 | 61 | 48.3 | 4.0210e-37 |
| 3 | fourier | N=40 | 81 | 81.8 | 5.3927e-38 |
| 3 | fourier | N=60 | 121 | 182.8 | 3.7742e-38 |
| 3 | fourier | N=80 | 161 | 322.6 | 3.2286e-38 |
| 3 | fourier | N=100 | 201 | 516.2 | 3.1355e-38 |
| 3 | breakpoint | d=4 | 30 | 14.4 | 6.8990e-13 |
| 3 | breakpoint | d=6 | 42 | 31.3 | 3.9053e-16 |
| 3 | breakpoint | d=8 | 54 | 58.5 | 4.1568e-20 |
| 3 | breakpoint | d=10 | 66 | 96.8 | 2.4500e-22 |
| 3 | breakpoint | d=12 | 78 | 149.2 | 7.5039e-26 |
| 3 | breakpoint | d=16 | 102 | 301.0 | 7.3392e-32 |
| 2.5 | breakpoint | d=20, ygl 32 | 95 | 193.4 | 4.7462e-24 |
| 3 | breakpoint | d=20 | 126 | 550.1 | 1.8149e-35 |
| 3 | breakpoint | d=24 | 150 | 875.8 | 2.4883e-38 |
| 3 | breakpoint | d=28, ygl 32 | 174 | 1745.9 | 2.4770e-38 |
| 3 | breakpoint | d=24, ygl 32, extra breaks 0.28768,0.33647,0.47 (predicted second-order positions) | 225 | 4598.0 | 2.4702e-38 |
| 3 | breakpoint | d=24, ygl 32, extra breaks 0.30,0.40,0.50 (control) | 225 | 4360 | 2.4699e-38 |
| 2.5 | breakpoint | d=24, ygl 32 | 113 | 302 | 4.7433e-24 |
| 2.5 | breakpoint | d=32, ygl 32 | 149 | 650 | 4.7414e-24 |
