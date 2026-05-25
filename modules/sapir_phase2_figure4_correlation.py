"""
SAPIR-Net Phase 2 — Figure 4 Generator
========================================
Single-panel figure showing the AMOX_CLAV Scenario B capacity-loss
distribution under four correlation values rho in {0, 0.25, 0.50, 0.75}
for the per-iteration ROW Pareto draws on PENI (HS 294110) and OTHAB
(HS 294190), coupled by a Gaussian copula. Generated from the
correlation-extension run outputs in
results/phase2_monte_carlo_correlation.csv and
results/phase2_monte_carlo_correlation_amox_clav_losses.npz.

Design
------
- Top panel: overlaid kernel density estimates of the AMOX_CLAV
  Scenario B capacity-loss distribution at each rho, zoomed to the
  93-99% loss range where the differences are visible.
- Bottom panel: forest-plot-style summary showing per-rho mean
  (point) and P5–P95 interval (horizontal line), with the rho=0
  reference annotated and the P(severe) = 100% annotation across
  all rho values (per the M-1 band-arithmetic-driven structural-
  guarantee finding).
- Severe-shortage threshold (30%) annotated for context, since the
  full plot range otherwise hides it.

Author: Lead Data Architect / Independent Researcher
Version: Phase 2 Figure 4 v1.0, May 2026
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import gaussian_kde

RESULTS_DIR = "results"

CORR_CSV = os.path.join(RESULTS_DIR, "phase2_monte_carlo_correlation.csv")
CORR_NPZ = os.path.join(RESULTS_DIR,
                        "phase2_monte_carlo_correlation_amox_clav_losses.npz")
OUTPUT_PNG = os.path.join(RESULTS_DIR,
                          "phase2_fig4_dual_chokepoint_correlation.png")

# Load results
corr_df = pd.read_csv(CORR_CSV)
amox_clav_df = corr_df[corr_df.Drug_ID == "AMOX_CLAV"].sort_values("Rho")
print("AMOX_CLAV Scenario B summary across rho values:")
print(amox_clav_df.to_string(index=False))

npz = np.load(CORR_NPZ)
print(f"\nLoaded loss arrays: {list(npz.keys())}")

# Color scheme: rho = 0 navy / rho = 0.25 teal / rho = 0.50 amber / rho = 0.75 burnt
colors = {
    0.00: "#1f3a68",
    0.25: "#1d8a8a",
    0.50: "#d18a1c",
    0.75: "#b03a2e",
}

rho_values = sorted(amox_clav_df["Rho"].unique())

fig, (ax_density, ax_forest) = plt.subplots(
    2, 1, figsize=(8.5, 7.5),
    gridspec_kw={"height_ratios": [3, 2], "hspace": 0.42},
)

# ------------------------------------------------------------
# Top panel: density overlays
# ------------------------------------------------------------
loss_grid = np.linspace(0.91, 0.99, 600)
for rho in rho_values:
    arr = npz[f"rho_{rho:.2f}"]
    kde = gaussian_kde(arr, bw_method=0.18)
    density = kde(loss_grid)
    ax_density.plot(loss_grid * 100, density,
                    color=colors[rho], linewidth=2.0,
                    label=f"rho = {rho:.2f}")
    ax_density.fill_between(loss_grid * 100, 0, density,
                            color=colors[rho], alpha=0.10)

ax_density.set_xlim(91, 99)
ax_density.set_xlabel("AMOX_CLAV capacity loss under Scenario B (%)",
                      fontsize=10)
ax_density.set_ylabel("Density (kernel density estimate)", fontsize=10)
ax_density.set_title(
    "Figure 4. AMOX_CLAV Scenario B capacity-loss distribution under\n"
    "partial-correlation sensitivity on the PENI–OTHAB ROW Pareto draws",
    fontsize=11, loc="left", pad=10,
)
ax_density.legend(loc="upper left", fontsize=9, framealpha=0.95,
                  title=r"Gaussian-copula $\rho$", title_fontsize=9)
ax_density.grid(alpha=0.25, linestyle="--", linewidth=0.5)
ax_density.tick_params(axis='both', labelsize=9)

# Subtle annotation about the zoom
ax_density.text(
    0.99, 0.98,
    "Zoomed to 91–99% loss range where the\n"
    "rho-induced shifts are visible. Full support\n"
    "of the loss distribution is approximately\n"
    "[85%, 100%] for all rho values.",
    transform=ax_density.transAxes,
    fontsize=8, verticalalignment="top",
    horizontalalignment="right",
    bbox=dict(facecolor="white", edgecolor="#888", alpha=0.92,
              boxstyle="round,pad=0.4"),
)

# ------------------------------------------------------------
# Bottom panel: forest-plot summary (mean with P5-P95 interval)
# ------------------------------------------------------------
y_positions = np.arange(len(rho_values))[::-1]  # rho=0 at top

for y, rho in zip(y_positions, rho_values):
    row = amox_clav_df[amox_clav_df.Rho == rho].iloc[0]
    mean = row["Mean_Capacity_Loss_pct"]
    p5 = row["P5_Loss_pct"]
    p95 = row["P95_Loss_pct"]

    ax_forest.plot([p5, p95], [y, y],
                   color=colors[rho], linewidth=2.5, solid_capstyle="butt")
    ax_forest.plot([p5, p5], [y - 0.12, y + 0.12],
                   color=colors[rho], linewidth=1.6)
    ax_forest.plot([p95, p95], [y - 0.12, y + 0.12],
                   color=colors[rho], linewidth=1.6)
    ax_forest.scatter([mean], [y], color=colors[rho], s=70,
                      zorder=5, edgecolor="white", linewidth=1.0)
    # Annotation: mean / P5 / P95
    ax_forest.text(p95 + 0.15, y,
                   f"  mean = {mean:.2f}%   P5 = {p5:.2f}%   "
                   f"P95 = {p95:.2f}%   P(severe) = "
                   f"{row['Prob_Severe_Shortage_pct']:.2f}%",
                   fontsize=8.5, verticalalignment="center")

ax_forest.set_yticks(y_positions)
ax_forest.set_yticklabels([f"$\\rho$ = {r:.2f}" for r in rho_values],
                          fontsize=9.5)
ax_forest.set_xlim(93.0, 100.3)
ax_forest.set_xlabel("AMOX_CLAV Scenario B capacity loss (%)", fontsize=10)
ax_forest.grid(alpha=0.25, linestyle="--", linewidth=0.5, axis="x")
ax_forest.tick_params(axis='x', labelsize=9)

ax_forest.set_title(
    "Per-rho mean (point) and P5–P95 interval (line), N = 10,000 iterations. "
    "P(severe shortage) is\n"
    "100% across all rho values — the binary indicator is band-arithmetic-driven "
    "per the\nSection 3.3.1 structural-guarantee result; "
    "rho induces only magnitude-distribution shifts.",
    fontsize=9, loc="left", pad=10,
)

# Severe-shortage threshold reference annotation (outside the zoom range)
ax_forest.axvspan(99.5, 100.0, alpha=0.0)  # placeholder for x-axis padding
ax_forest.text(
    0.99, -0.32,
    "30% severe-shortage threshold sits far below the displayed loss range; "
    "the entire AMOX_CLAV\nScenario B loss distribution is "
    "above the severe-shortage threshold at every rho.",
    transform=ax_forest.transAxes,
    fontsize=8, verticalalignment="top", horizontalalignment="right",
    style="italic", color="#444",
)

# Final adjustments
plt.subplots_adjust(left=0.10, right=0.97, top=0.93, bottom=0.13)

plt.savefig(OUTPUT_PNG, dpi=160, bbox_inches="tight", facecolor="white")
plt.close(fig)

print(f"\nFigure 4 written: {OUTPUT_PNG}")
print("File size:", os.path.getsize(OUTPUT_PNG), "bytes")
