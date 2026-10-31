# Breakpoints at the prime-power positions: computing the truncated Weil-form minimizer

Code, logs and eigenvectors for the note

> J. Echt, *An Improved Computation of the Truncated Weil-form Minimizer: Breakpoints at the Prime-power Positions*, October 2026 (`weil_form_breakpoints.md` in this repository).

The note (`weil_form_breakpoints.md`) shows that the Fourier basis used in every published computation of the Connes–Consani–Moscovici truncated Weil form does not converge to its lowest eigenvalue ε_λ at any practical size (at dimension 201 it is at least 2.6%, 14% and 27% above the minimum at λ = 2, 2.5 and 3), explains why (the minimizer has corners at y = log p^m and L − log p^m, one pair for each prime power p^m ≤ λ²), and gives a piecewise-polynomial basis with cuts at these corners' positions that reaches within 0.1% of the minimum at dimension 26 at λ = 2. Everything in the note can be reproduced from this repository on ordinary hardware.

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
| `dumps/` | the eigenvectors (JSON) that the diagnostics were run on, plus the zeros cache for `bp_zerosum.py`. |

## Running

Every script prints its options with `--help`. The two that matter most:

- `--dps` sets the working precision in decimal digits (50–140 in the runs reported; ε_λ is 10⁻¹³ at λ = 2 and 10⁻³⁸ at λ = 3, so the precision has to be well above the digits of ε).
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

- Each log starts with the parameters of the run (λ, L, N or d, digits, the prime powers used) and, for the Fourier basis, three `[check] W02` lines comparing the pole-term quadrature with the paper's closed form.
- The `[…s]` stamps inside a log are matrix-assembly time only. Wall times for the time-to-accuracy ladders are in `logs/ladders/tta.csv`.
- Every value printed is a Rayleigh–Ritz upper bound on ε_λ (given exact matrix entries), so the smaller of two values from any two bases is the closer.
- The ladder logs in `logs/ladders/` were produced by an earlier version of the scripts and print `diagonalising` where the current scripts print `diagonalizing`; the numbers are unaffected. The N = 100 Fourier rows and several breakpoint rows were rerun with the current scripts (see `tta_table.md`) and reproduce the ladder values to every printed digit.

## Notation

The note is plain text: `X_n` is a subscript, `X^n` a superscript, `{…}` groups, and `Λ(p^m)` is the von Mangoldt function. Variables: u ∈ [λ⁻¹, λ] is the paper's variable, y = log λu ∈ [0, L] with L = 2 log λ is the variable the breakpoints are stated in, and Suzuki's x = y − a with a = log λ.

## Citing

If you use the basis or the logs, please cite the note (`weil_form_breakpoints.md`) and the release tag of this repository.

## License

MIT
