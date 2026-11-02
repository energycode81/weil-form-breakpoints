# Log index

One line per log file: the command that produced it and what it backs in the note (`weil_form_breakpoints.md`). Section references are to that note: the Abstract, Setting, Section 1 (the Fourier values and Groskin's tables), Section 2 (the corners, the breakpoint basis and its three tests), Section 3 (the Figure 4 comparison), and "What is not claimed". Logs marked *supporting* are runs whose numbers are not quoted in the note but which the quoted numbers were checked against; they are kept so that nothing behind the note is missing.

All runs are `python3 <script> …` with the scripts in the repository root; `tee logs/<path>` captures the output. The eigenvectors that the diagnostics were run on are in `dumps/`. The `[…s]` stamps inside the logs are matrix-assembly time only; wall times for the ladders are in `ladders/tta.csv` and `tta_table.md`.

Logs are grouped by what produced them, not by when they were run:

| directory | contents |
|---|---|
| `fourier/` | `weil_step1.py` runs in the Fourier basis V_n, one per (λ, N) |
| `fourier/ablation/` | the same with primes removed (`--max-prime`, `--no-primes`) |
| `breakpoint/` | `weil_breakpoints.py` runs with cuts at the prime-power positions |
| `breakpoint/controls/` | the same basis with the cuts moved, added or removed |
| `legendre/` | `weil_legendre.py` runs and the `leg_window.py` residual scans on them |
| `diagnostics/` | scripts run on a dumped eigenvector: `bp_jumps.py`, `bp_zerosum.py`, `bp_piece_decay.py`, `weil_step1.py --endpoint` |
| `validation/` | the independent evaluation of a matrix entry |
| `ladders/` | the time-to-accuracy ladders behind `tta_table.md`, both bases, three λ |

## validation/

| file | command | backs |
|---|---|---|
| `check1_lam2.log` | `check1.py` | Setting ("validated the matrix entries against the definition of the form"): a₀ = 0.0697042428372728… from the definition; archimedean 2.0623050802620368, pole 2.8853900817779268, primes 0.7533807586786172 |
| `lam2_N10_check.log` | `weil_step1.py --lam 2 --N 10 --dps 30 --no-prolate --no-spectrum` | Setting: the matrix side, a₀ = 0.069704243; `[check] W02` lines |

## fourier/

| file | command | backs |
|---|---|---|
| `lam2_N100.log` | `weil_step1.py --lam 2 --N 100 --dps 50 --no-prolate --no-spectrum` | Abstract, Section 1: N = 100 value 7.71·10⁻¹³, 2.6% above the breakpoint value |
| `lam25_N100.log` | `weil_step1.py --lam 2.5 --N 100 --dps 60 --no-prolate --no-spectrum` | Abstract, Section 1: 5.41·10⁻²⁴, 14% above; Section 3: Fourier offset 2.26 |
| `lam3_N100.log` | `weil_step1.py --lam 3 --N 100 --dps 90 --no-prolate --no-spectrum` | Abstract, Section 1: 3.14·10⁻³⁸, 27% above |
| `lam2_N60.log` | `weil_step1.py --lam 2 --N 60 --dps 60 --no-prolate --no-spectrum --dump xi_lam2_N60.json` | Section 2 ("below the Fourier value at dimension 121"): 8.13·10⁻¹³ |
| `lam2_N30_full.log` | `weil_step1.py --lam 2 --N 30 --dps 60` | Section 3: Fourier offset 2.03 at λ = 2 |
| `lam3_N60_full.log` | `weil_step1.py --lam 3 --N 60 --dps 100` | Setting: zeta-zero accuracy 2·10⁻³⁴ at λ = 3 with a 121×121 matrix; Section 3: Fourier offset 2.51 |
| `lam35_N60_full.log` | `weil_step1.py --lam 3.5 --N 60 --dps 85` | Section 3: Fourier offset ≤ 2.94 at λ = 3.5 |
| `lam42_N100.log`, `lam42_N120.log` | `weil_step1.py --lam 4.2 --N <100, 120> --dps 140 …` | What is not claimed: the largest λ run (4.2, prime powers ≤ 17); Fourier only |
| `lamsqrt2_N6.log`, `lamsqrt2_N10.log` | `weil_step1.py --lam 1.41421356237 --N <6, 10> --dps 30 …` | supporting: Setting (simple, even ground state at every λ run) |
| `run_lam3_N20.log`, `lam3_N20_full90.log` | `weil_step1.py --lam 3 --N 20 --dps <60, 90>` | supporting: zeta-zero accuracy at smaller N |
| `lam3_N40_full.log`, `lam3_N10_full.log`, `lam3_N30_full.log`, `lam3_N50_full.log` | `weil_step1.py --lam 3 --N <N> --dps <90, 100>` | supporting: the λ = 3 Fourier ladder with spectrum and prolate comparison |
| `lam2_N8_full.log`, `lam25_N<30,40,60,80,100>_full.log` | `weil_step1.py --lam <λ> --N <N> --dps <30, 70>` | supporting: the λ = 2 and 2.5 Fourier ladders with spectrum and prolate comparison |
| `lam25_N60.log`, `lam3_N60_dump.log`, `lam35_N80.log` | `weil_step1.py … --no-prolate --no-spectrum --dump xi_<λ>_N<N>.json` | supporting: eigenvector dumps for the zero sums and `endpoint_all` |
| `lam25_N60_zeros.log`, `lam3_N60_zeros.log`, `lam35_N80_zeros.log` | `weil_step1.py --zeros xi_<λ>_N<N>.json --M 3000 --dps 30 --no-mellin` | supporting: explicit-formula zero sums for the Fourier eigenvectors (0.972 + tail, 0.963 + tail, 0.958 + tail of ε) |

## fourier/ablation/

| file | command | backs |
|---|---|---|
| `lam25_N40_maxprime4.log`, `_maxprime3`, `_noprimes` | `weil_step1.py --lam 2.5 --N 40 --dps 40 --no-prolate --no-spectrum --max-prime <4, 3>` / `--no-primes` | supporting: with primes removed the minimum turns negative (−0.50 when 5 is removed) |
| `lam3_N40_maxprime7.log`, `_maxprime5`, `_noprimes` | `weil_step1.py --lam 3 --N 40 --dps 60 --no-prolate --no-spectrum --max-prime <7, 5>` / `--no-primes` | supporting: the same at λ = 3 (−0.14, −0.67, −1.53) |

## breakpoint/

| file | command | backs |
|---|---|---|
| `lam2_d28.log`, `lam2_d36.log`, `lam2_d48.log` | `weil_breakpoints.py --lam 2 --d <d> --dps 50 --dump bp_lam2_d<d>_v2.json` (`--ygl` scaled with d) | Abstract, Section 1: best λ = 2 value 7.5187·10⁻¹³ at dimension 98; last ladder step 5·10⁻⁵ (7.5191 → 7.5187) |
| `lam25_correct_d32.log` | `weil_breakpoints.py --lam 2.5 --d 32 --dps 60 --ygl 32 --dump bp_lam25_d32.json` | Abstract, Section 1: best λ = 2.5 value 4.7414·10⁻²⁴ at dimension 149 |
| `lam25_correct_d24.log` | same, `--d 24` | Section 1: last ladder step 4·10⁻⁴ (4.7433 → 4.7414) |
| `lam25_correct_d32_ygl72.log` | `weil_breakpoints.py --lam 2.5 --d 32 --dps 60 --dump …` (default 72-point rule) | Section 1: matrix unchanged to 16 digits under quadrature refinement at degree 32 |
| `lam25_correct_d20.log` | `weil_breakpoints.py --lam 2.5 --d 20 --dps 60 --ygl 32 --dump bp_lam25_d20.json` | Section 2 Positions: correct-position value 4.7462·10⁻²⁴ |
| `lam2_d64.log` | `weil_breakpoints.py --lam 2 --d 64 --dps 50 --dump bp_lam2_d64.json` | Section 1: out-of-sample test of the contraction assumption — predicted 7.5183–7.5186·10⁻¹³, run gave 7.51850·10⁻¹³ |
| `lam22_correct_d12.log`, `_d16`, `_d20` | `weil_breakpoints.py --lam 2.2 --d <d> --dps 60 --ygl 32` | Section 2 Positions: reference ladder 5.73261, 5.72700, 5.72584 ×10⁻¹⁷ |
| `lam3_d24_plain.log` | `weil_breakpoints.py --lam 3 --d 24 --dps 90 --ygl 32` | Section 2 (second-order test): reference value 2.48829·10⁻³⁸ |

## breakpoint/controls/

| file | command | backs |
|---|---|---|
| `lam25_shifted_d20.log` | `weil_breakpoints.py --lam 2.5 --d 20 --dps 60 --ygl 32 --no-breaks --extra-breaks 0.64315,1.04861,1.33629,1.55944` | Abstract, Section 2 Positions: 4.8020·10⁻²⁴ (+1.2%) |
| `lam25_midpoint_d20.log` | `… --no-breaks --extra-breaks 0.11157,0.33472,0.56972,0.71356` | Abstract, Section 2 Positions: 4.8337·10⁻²⁴ (+1.8%) |
| `lam22_shifted_d12.log`, `_d16`, `_d20` | `weil_breakpoints.py --lam 2.2 --d <d> --dps 60 --ygl 32 --no-breaks --extra-breaks 0.64315,1.04861,1.33629` | Section 2 Positions: +1.9%, +1.3%, +0.94% |
| `lam22_midpoint_d20.log` | `… --lam 2.2 --d 20 … --no-breaks --extra-breaks 0.09531,0.33446,0.58573` | supporting: midpoints at λ = 2.2 (+1.7%) |
| `lam3_d24_predextras.log` | `weil_breakpoints.py --lam 3 --d 24 --dps 90 --ygl 32 --extra-breaks 0.28768,0.33647,0.47000` | Abstract, Section 1: best λ = 3 value 2.470·10⁻³⁸ at dimension 225; Section 2: cuts at the predicted second-order positions |
| `lam3_d24_wrongextras.log` | same with `--extra-breaks 0.30,0.40,0.50` | Section 2: arbitrary cuts in the same piece give the same improvement (2.4699·10⁻³⁸) |
| `lam2_d12_cut_log32.log`, `lam22_secondorder_cuts.log`, `lam25_center_cut.log` | single extra cuts (`--extra-breaks`) at λ = 2, 2.2, 2.5 | supporting: single-cut tests without a displaced control |

## legendre/

| file | command | backs |
|---|---|---|
| `leg2_M48.log` | `weil_legendre.py --lam 2 --M 48 --dps 50 --dump leg2.json` | Section 2 ("below the Legendre value at degree 48"): 7.730·10⁻¹³ |
| `leg22_window.log` | `leg_window.py leg22_M80.json 40` | Section 2 Positions: residual scan at λ = 2.2, no preference for the prime positions |
| `leg25_window.log` | `leg_window.py leg25_M120.json 60` | Section 2 Positions: the same at λ = 2.5 |
| `leg22_M80.log`, `leg25_M120.log` | `weil_legendre.py --lam <2.2, 2.5> --M <80, 120> --dps 60 --dump …` | supporting: the eigenvectors the two scans were run on |
| `leg25_M80_M100.log` | `weil_legendre.py --lam 2.5 --M 80,100 --dps 60` | supporting: the λ = 2.5 Legendre ladder |

## diagnostics/

| file | command | backs |
|---|---|---|
| `bp_jumps_bp_lam25_d<16,20,24,32>.log` | `bp_jumps.py bp_lam25_d<d>.json` | Section 2 Amplitude: within 10% across prime powers at degree 32 (0.128–0.141), within a factor of 1.5 at degree 16; ratio from 0.16–0.23 to 0.13–0.14 |
| `bp_jumps_bp_lam22_d<16,20,24>.log` | `bp_jumps.py bp_lam22_d<d>.json` | Section 2 Amplitude: within 5–17% at λ = 2.2 |
| `bp_jumps_bp_lam2_d<8…64>.log` (8 files) | `bp_jumps.py bp_lam2_d<d>.json` (d = 8–20, 64), `bp_jumps.py bp_lam2_d<d>_v2.json` (d = 28, 36, 48) | Section 2 Amplitude: ratio 0.19 → 0.119 over degrees 12–64 at λ = 2 |
| `bp_zerosum_lam2_d20.log` | `bp_zerosum.py bp_lam2_d20.json --zeros-cache zeros.json` | Section 2 Matrix-free check: 0.992 from 3000 zeros, 0.9995 with fitted t^{−2.07} tail |
| `lam22_piece_decay.log` | `bp_piece_decay.py bp_lam22_d16.json bp_lam22_d20.json bp_lam22_d24.json` | supporting: Section 2 ("what remains … on the prime-bounded interior pieces") |
| `endpoint_all.log` | `weil_step1.py --endpoint xi_lam2_N60.json,xi_lam25_N60.json,xi_lam3_N60.json,xi_lam35_N80.json,xi_lam42_N100.json` | supporting: endpoint values of the Fourier eigenvectors |

## ladders/

| files | command | backs |
|---|---|---|
| `tta_fourier_lam<λ>_<N>.log` (23 files; λ = 2, 2.5, 3; N = 10 … 100) | `weil_step1.py --lam <λ> --N <N> --dps <50/60/90> --no-prolate --no-spectrum` | Section 1: the fast phase ending near N ≈ 30 and 40, the local exponents 0.3–2 (λ = 3: 1.99, 1.88, 0.59), the 1/N extrapolations; `tta_table.md` |
| `tta_breakpoint_lam<λ>_<d>.log` (25 files; d = 4 … 20, plus 24, 24_x, 28 at λ = 3 and 20_ygl32 at λ = 2.5) | `weil_breakpoints.py --lam <λ> --d <d> --dps <50/60/90> --ygl 24` | Abstract, Section 2: 0.1% at dimension 26 at λ = 2 (d = 12); Section 1: quadrature check at degree 20 (`20` vs `20_ygl32`); λ = 3 last ladder step; `tta_table.md` |
| `tta.sh` | — | the script that ran the 43 ladder rows of `tta.csv` (N = 10 … 100; d = 4 … 20 at λ = 2 and 2.5, d = 4 … 16 at λ = 3) and timed each run |
| `tta.csv` | `tta.sh` | wall times of those 43 rows; the five later ladder runs (`lam2.5_20_ygl32`, `lam3_20`, `lam3_24`, `lam3_24_x`, `lam3_28`) were made individually and their wall times are recorded only in `tta_table.md` |
