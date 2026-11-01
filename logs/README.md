# Breakpoints at the prime-power positions: computing the truncated Weil-form minimizer

Code, logs and eigenvectors for the note

> J. Echt, *An Improved Computation of the Truncated Weil-form Minimizer: Breakpoints at the Prime-power Positions*, October 2026 (`weil_form_breakpoints.md` in this repository).

The note (`weil_form_breakpoints.md`) shows that the Fourier basis used in every published computation of the Connes–Consani–Moscovici truncated Weil form (*Zeta Spectral Triples*, arXiv:2511.22755) does not converge to its lowest eigenvalue ε_λ at any practical size (at dimension 201 it is at least 2.6%, 14% and 27% above the minimum at λ = 2, 2.5 and 3), explains why (the minimizer has corners at y = log p^m and L − log p^m, one pair for each prime power p^m ≤ λ²), and gives a piecewise-polynomial basis with cuts at these corners' positions that, at λ = 2, comes within 0.1% of its dimension-130 value already at dimension 26. Everything in the note can be reproduced from this repository on ordinary hardware.

This analysis was produced with the assistance of Anthropic's Claude (Fable 5.1). Every computation was run by the author, who is responsible for all content.

## Contents

| path | what it is |
|---|---|
| `weil_form_breakpoints.md` | the note |
| `weil_step1.py` | Fourier basis: builds the (2N+1)×(2N+1) matrix, finds ε_N and ξ; options for the spectrum of D_log, the prolate comparison, the explicit-formula zero sum (`--zeros`), endpoint diagnostics (`--endpoint`), and ablation (`--max-prime`, `--no-primes`) |
| `weil_legendre.py` | Legendre basis (even polynomials on the whole interval, no cuts); the control basis of Section 2 of the note |
| `weil_breakpoints.py` | breakpoint basis: piecewise polynomials with cuts at the prime-power positions; `--extra-breaks` and `--no-breaks` for the displacement controls |
| `bp_jumps.py`, `bp_zerosum.py`, `bp_piece_decay.py`, `leg_window.py` | diagnostics run on a dumped eigenvector: jump amplitudes at the cuts, the matrix-free zero sum, per-piece coefficient decay, the Legendre residual scan |
| `check1.py` | one matrix entry (a₀ at λ = 2) recomputed from the definition of the form, independently of the other scripts |
| `logs/` | the output of every run quoted in the note, grouped by what produced it: `fourier/`, `fourier/ablation/`, `breakpoint/`, `breakpoint/controls/`, `legendre/`, `diagnostics/`, `validation/`, `ladders/` |
| `logs/INDEX.md` | one line per log: the command that produced it and the numbers in the note it backs |
| `tta_table.md` | time-to-accuracy ladders of the Fourier and breakpoint bases at λ = 2, 2.5, 3 (ε and wall time against matrix size) |
| `fig_convergence_lam2.png`, `make_figure.py` | Figure 1 of the note (the λ = 2 convergence of both bases) and the script that draws it from the logged values |
| `dumps/` | the eigenvectors (JSON) that the diagnostics were run on, plus the zeros cache for `bp_zerosum.py` |

## Running

Every script prints its options with `--help`. The two that matter most:

- `--dps` sets the working precision in decimal digits (50–140 in the runs reported; ε_λ is 10⁻¹³ at λ = 2 and 10⁻³⁸ at λ = 3, so the precision has to be comfortably more than 13 and 38 digits respectively).
- `--ygl` (breakpoint basis only) is the number of Gauss–Legendre quadrature nodes per integration segment; the default is max(24, 2d + 8) for degree d. The matrix at degree 32 is unchanged to 16 digits between 32 and 72 nodes.

The commands behind the headline numbers, each with the log it produced:

```bash
# Fourier values at N = 100 (the 2.6 / 14 / 27 % figures)
python3 weil_step1.py --lam 2   --N 100 --dps 50 --no-prolate --no-spectrum     # logs/fourier/lam2_N100.log
python3 weil_step1.py --lam 2.5 --N 100 --dps 60 --no-prolate --no-spectrum     # logs/fourier/lam25_N100.log
python3 weil_step1.py --lam 3   --N 100 --dps 90 --no-prolate --no-spectrum     # logs/fourier/lam3_N100.log

# best breakpoint values
python3 weil_breakpoints.py --lam 2   --d 48 --dps 50 --dump bp_lam2_d48_v2.json              # logs/breakpoint/lam2_d48.log
python3 weil_breakpoints.py --lam 2.5 --d 32 --dps 60 --ygl 32 --dump bp_lam25_d32.json   # logs/breakpoint/lam25_correct_d32.log
python3 weil_breakpoints.py --lam 3 --d 24 --dps 90 --ygl 32 --extra-breaks 0.28768,0.33647,0.47000   # logs/breakpoint/controls/lam3_d24_predextras.log

# displacement control at lambda = 2.5, degree 20: correct positions, shifted by -0.05, midpoints
python3 weil_breakpoints.py --lam 2.5 --d 20 --dps 60 --ygl 32 --dump bp_lam25_d20.json                                      # logs/breakpoint/lam25_correct_d20.log
python3 weil_breakpoints.py --lam 2.5 --d 20 --dps 60 --ygl 32 --no-breaks --extra-breaks 0.64315,1.04861,1.33629,1.55944    # logs/breakpoint/controls/lam25_shifted_d20.log
python3 weil_breakpoints.py --lam 2.5 --d 20 --dps 60 --ygl 32 --no-breaks --extra-breaks 0.11157,0.33472,0.56972,0.71356    # logs/breakpoint/controls/lam25_midpoint_d20.log

# amplitude at the cuts, and the matrix-free check
python3 bp_jumps.py bp_lam25_d20.json                                     # logs/diagnostics/bp_jumps_bp_lam25_d20.log
python3 bp_zerosum.py bp_lam2_d20.json --zeros-cache zeros.json           # logs/diagnostics/bp_zerosum_lam2_d20.log
```

`logs/INDEX.md` has the command for every other log. Run times depend on the machine; every run in the note was made on ordinary local hardware.

## Reading the logs

- Each log starts with the parameters of the run (λ, L, N or d, digits, the prime powers used) and, for the Fourier basis, three `[check] W02` lines comparing the pole-term quadrature with the closed form (4.2) of arXiv:2511.22755.
- The `[…s]` stamps inside a log are matrix-assembly time only. Wall times for the time-to-accuracy ladders are in `logs/ladders/tta.csv`.
- Every value printed is a Rayleigh–Ritz upper bound on ε_λ (given exact matrix entries), so the smaller of two values from any two bases is the closer. The exception is `fourier/ablation/`: there the primes are removed or truncated, so the form is no longer the Weil form and its eigenvalues bound nothing about ε_λ.
- The ladder logs in `logs/ladders/` were produced by an earlier version of the scripts and print `diagonalising` where the current scripts print `diagonalizing`; the numbers are unaffected. The N = 100 Fourier rows and several breakpoint rows were rerun with the current scripts (see `tta_table.md`) and reproduce the ladder values to every printed digit.

### Log vocabulary

- **ξ (xi)** is the computed minimizer, normalized to unit L² norm on [λ⁻¹, λ]; **eps_N / eps_d / eps_M** is the lowest eigenvalue in the Fourier, breakpoint or Legendre basis. All positions are in the variable y = log λu ∈ [0, L].
- **Breakpoint logs** (`weil_breakpoints.py`): `breakpoints` lists the cut positions in y; `endpoint/peak` is ξ(0)/ξ(L/2), the endpoint value relative to the peak; the "top two degrees per piece" are the two highest-degree Legendre coefficients on each piece, a per-piece convergence indicator.
- **`bp_jumps` logs**: one row per cut. `at` is the cut position with its name (log q or L − log q). `xi−` and `xi+` are the one-sided values of ξ at the cut, evaluated from the polynomial pieces to its left and right; `jump` is their difference; `/wv` is that jump divided by Λ(q) q^(−1/2) times the endpoint value v = ξ(0) printed in the header. The second and third blocks repeat this for the first and second derivatives ξ′ and ξ″. The quantity quoted in Section 2 of the note (the amplitude ratio) is `/wv` in the first (value) block.
- **`bp_zerosum` logs**: partial sums of the explicit-formula zero sum, reported as fractions of eps_d.
- **`leg_window` logs**: the Euler–Lagrange residual of the Legendre minimizer in a small window around each candidate position.
- **`bp_piece_decay` logs**: the decay of each piece's Legendre coefficients across degrees.

## Notation

The note is plain text: `X_n` is a subscript, `X^n` a superscript, `{…}` groups, and `Λ(p^m)` is the von Mangoldt function. Variables: u ∈ [λ⁻¹, λ] is the paper's variable, y = log λu ∈ [0, L] with L = 2 log λ is the variable the breakpoints are stated in, and Suzuki's x = y − a with a = log λ.

## Citing

If you use the basis or the logs, please cite the note (`weil_form_breakpoints.md`) and the release tag of this repository.

## License

MIT
