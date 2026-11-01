#!/usr/bin/env python3
"""Figure 1 of the note: convergence of the Fourier and breakpoint bases at lambda = 2.

Plots the relative excess of each computed value over the dimension-130 breakpoint
value (logs/breakpoint/lam2_d64.log). Every point is taken from a log in logs/:

  Fourier (weil_step1.py), logs/ladders/tta_fourier_lam2_<N>.log, N = 10..100
  Breakpoint (weil_breakpoints.py), logs/ladders/tta_breakpoint_lam2_<d>.log, d = 4..20,
    then logs/breakpoint/lam2_d28.log, lam2_d36.log, lam2_d48.log, lam2_d64.log

Output: fig_convergence_lam2.png. Requires matplotlib (used only for this figure).
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REF = 7.51850488862100e-13          # dimension 130 (degree 64), logs/breakpoint/lam2_d64.log

# (matrix dimension, eps)
fourier = [                          # dim = 2N + 1
    (21, 1.950e-12), (31, 1.068e-12), (41, 1.032e-12), (61, 0.964e-12),
    (81, 0.888e-12), (121, 0.813e-12), (161, 0.784e-12), (201, 0.771e-12),
]
breakpoint = [                       # dim = 2(d + 1) for the 4 pieces at lambda = 2
    (10, 5.89e-8), (14, 1.21e-9), (18, 4.06e-12), (22, 7.580e-13),
    (26, 7.5247e-13), (34, 7.5228e-13), (42, 7.5212e-13),
    (58, 7.51971e-13), (74, 7.51909e-13), (98, 7.51870e-13),
]

BLUE, VERM = "#0072B2", "#D55E00"
fx, fy = zip(*[(d, v / REF - 1) for d, v in fourier])
bx, by = zip(*[(d, v / REF - 1) for d, v in breakpoint])

fig, ax = plt.subplots(figsize=(7.2, 4.6), dpi=200)
ax.semilogy(fx, fy, "-o", color=VERM, lw=2, ms=6, label="Fourier basis")
ax.semilogy(bx, by, "-s", color=BLUE, lw=2, ms=6, label="breakpoint basis")

ax.annotate("dimension 201: still 2.6% high", (201, fy[-1]), xytext=(128, 0.09),
            fontsize=9, color=VERM,
            arrowprops=dict(arrowstyle="-", color=VERM, lw=0.8))
ax.annotate("dimension 26: 0.08% high", (26, by[4]), xytext=(44, 2.6e-3),
            fontsize=9, color=BLUE,
            arrowprops=dict(arrowstyle="-", color=BLUE, lw=0.8))

ax.set_xlim(0, 212)
ax.set_ylim(1e-5, 12)
ax.set_xlabel("matrix dimension")
ax.set_ylabel("excess over the lowest computed value")
ax.set_title("λ = 2: both bases are upper bounds; only one settles", fontsize=11)
ax.grid(True, which="major", color="0.88", lw=0.6)
ax.grid(True, which="minor", color="0.94", lw=0.4)
ax.legend(frameon=False, loc="upper right")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.text(0.01, 0.012,
         "Every value is an upper bound on ε_λ; each point is its excess over the lowest value computed, the degree-64\n"
         "breakpoint value (7.51850·10⁻¹³, dimension 130). Breakpoint points at dimensions 10 and 14 lie above the frame.\n"
         "Data: logs/ladders/, logs/breakpoint/.",
         fontsize=9, color="0.35")
fig.tight_layout(rect=(0, 0.115, 1, 1))
fig.savefig("fig_convergence_lam2.png")
print("wrote fig_convergence_lam2.png")
