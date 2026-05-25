"""
SAPIR-Net Phase 2 Module 4A: Data Visualization (Antibiotics)
=============================================================
Generates publication-quality figures from the Phase 2 Monte Carlo
results. Mirrors the Phase 1 Module 4A architecture; antibiotic-
specific labels, four-drug schema, and a third figure (federal-policy-
aware sensitivity) added per Brief Section 4.1.

Inputs
------
results/phase2_monte_carlo_primary.csv
results/phase2_monte_carlo_sensitivity.csv
results/phase2_monte_carlo_raw_losses.npz

Outputs (300 dpi PNG)
---------------------
results/phase2_fig1_vulnerability_heatmap.png
results/phase2_fig2_capacity_loss_densities.png
results/phase2_fig3_federal_policy_sensitivity.png

Style notes
-----------
Restrained academic-policy convention: single-hue sequential heat-map
(white -> muted oxblood); muted navy / steel / oxblood / ochre / slate
palette; no rainbow colour ramps; clean white background; spines
hidden top/right; minimal gridlines; Helvetica/Arial preferred fonts.

Dependencies: numpy, scipy, pandas, matplotlib

Author: Lead Data Architect / Independent Researcher
Version: Phase 2 Module 4A v1.0, May 2026
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
from scipy.stats import gaussian_kde


# ============================================================
# PATHS
# ============================================================

RESULTS_DIR = "results"
PRIMARY_CSV = os.path.join(RESULTS_DIR, "phase2_monte_carlo_primary.csv")
SENSITIVITY_CSV = os.path.join(RESULTS_DIR, "phase2_monte_carlo_sensitivity.csv")
RAW_NPZ = os.path.join(RESULTS_DIR, "phase2_monte_carlo_raw_losses.npz")

FIG1_PATH = os.path.join(RESULTS_DIR, "phase2_fig1_vulnerability_heatmap.png")
FIG2_PATH = os.path.join(RESULTS_DIR, "phase2_fig2_capacity_loss_densities.png")
FIG3_PATH = os.path.join(RESULTS_DIR, "phase2_fig3_federal_policy_sensitivity.png")


# ============================================================
# PALETTE
# ============================================================
# Hand-picked academic-policy convention (Brookings / RAND /
# JAMA Health Forum range). Restrained, muted, low-saturation.

NAVY     = "#1B2A4A"   # primary dark text and emphasis
STEEL    = "#3F5E7A"   # secondary fill / axes
SLATE    = "#7A8B99"   # tertiary / faint elements
OXBLOOD  = "#8B2635"   # risk emphasis (used sparingly)
OXBLOOD_LIGHT = "#C97A82"  # mid-range of heatmap ramp
OCHRE    = "#B07A1F"   # accent for lower-risk scenario
SAND     = "#E8DEC3"   # light fill for ranges
NEAR_WHITE = "#FBFAF7"  # off-white background ramp anchor
WHITE    = "#FFFFFF"


# ============================================================
# GLOBAL MATPLOTLIB STYLE
# ============================================================

plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.sans-serif":   ["Helvetica", "Arial", "Liberation Sans",
                           "DejaVu Sans"],
    "font.size":         10.5,
    "axes.titlesize":    12.5,
    "axes.titleweight":  "normal",
    "axes.labelsize":    10.5,
    "axes.facecolor":    WHITE,
    "figure.facecolor":  WHITE,
    "axes.edgecolor":    SLATE,
    "axes.linewidth":    0.8,
    "axes.grid":         False,
    "xtick.color":       NAVY,
    "ytick.color":       NAVY,
    "xtick.labelsize":   9.5,
    "ytick.labelsize":   9.5,
    "text.color":        NAVY,
    "legend.frameon":    False,
    "legend.fontsize":   9.5,
    "figure.dpi":        110,
    "savefig.dpi":       300,
    "savefig.bbox":      "tight",
    "savefig.facecolor": WHITE,
})


# ============================================================
# LOAD RESULTS
# ============================================================

primary_df = pd.read_csv(PRIMARY_CSV)
sensitivity_df = pd.read_csv(SENSITIVITY_CSV)
raw_losses = np.load(RAW_NPZ)

# Canonical drug ordering used throughout
DRUG_ORDER = ["AMOX", "AMOX_CLAV", "AZI", "DOXY"]
DRUG_LABELS = {
    "AMOX":      "Amoxicillin",
    "AMOX_CLAV": "Amoxicillin/clavulanate\n(Augmentin)",
    "AZI":       "Azithromycin",
    "DOXY":      "Doxycycline",
}
DRUG_LABELS_INLINE = {  # single-line variant for panel titles
    "AMOX":      "Amoxicillin",
    "AMOX_CLAV": "Amoxicillin/clavulanate (Augmentin)",
    "AZI":       "Azithromycin",
    "DOXY":      "Doxycycline",
}

SCENARIO_ORDER = [
    "Baseline (No Shock)",
    "A: Direct Ban (CN+IN)",
    "B: Cascading Upstream Shock",
    "C: Logistics Chokepoint",
]
# Short column labels for the heat-map
SCENARIO_SHORT = {
    "Baseline (No Shock)":         "Baseline\n(No Shock)",
    "A: Direct Ban (CN+IN)":       "A - Direct Ban\n(CN + IN)",
    "B: Cascading Upstream Shock": "B - Cascading\nUpstream Shock",
    "C: Logistics Chokepoint":     "C - Logistics\nChokepoint",
}


def _npz_key(scen_name, drug):
    """Reproduce Module 3's npz key sanitisation."""
    scen_key = (scen_name.replace(" ", "_")
                .replace(":", "")
                .replace("(", "")
                .replace(")", "")
                .replace("+", "and"))
    return f"{scen_key}__{drug}"


# ============================================================
# FIGURE 1 - VULNERABILITY HEAT-MAP
# ============================================================
# 4 drugs (rows) x 4 scenarios (columns); cells annotated with
# P(severe shortage) %. Single-hue sequential palette (white ->
# muted oxblood). AMOX_CLAV x Scenario C cell flagged with a
# distinct border as the dual-chokepoint signature.

def make_fig1_heatmap():
    pivot = primary_df.pivot_table(
        index="Drug_ID", columns="Scenario",
        values="Prob_Severe_Shortage_pct",
        aggfunc="first",
    )
    pivot = pivot.reindex(index=DRUG_ORDER, columns=SCENARIO_ORDER)
    data = pivot.values  # shape (4, 4)

    cmap = LinearSegmentedColormap.from_list(
        "phase2_risk",
        [NEAR_WHITE, "#F0E6D4", "#D9B59E", OXBLOOD_LIGHT, OXBLOOD],
        N=256,
    )

    fig, ax = plt.subplots(figsize=(8.6, 4.8))

    im = ax.imshow(data, cmap=cmap, aspect="auto", vmin=0, vmax=100)

    # Cell annotations
    n_rows, n_cols = data.shape
    for i in range(n_rows):
        for j in range(n_cols):
            val = data[i, j]
            text_color = WHITE if val > 55 else NAVY
            weight = "semibold" if val > 30 else "normal"
            if val == 0.0:
                label = "0%"
            elif abs(val - round(val)) < 0.005:
                label = f"{val:.0f}%"
            else:
                label = f"{val:.2f}%"
            ax.text(j, i, label, ha="center", va="center",
                    color=text_color, fontsize=11, fontweight=weight)

    # Highlight the AMOX_CLAV x Scenario C cell (dual-chokepoint
    # signature). Border drawn as a thin solid rectangle.
    amox_clav_row = DRUG_ORDER.index("AMOX_CLAV")
    scen_c_col = SCENARIO_ORDER.index("C: Logistics Chokepoint")
    rect = Rectangle(
        (scen_c_col - 0.5, amox_clav_row - 0.5), 1.0, 1.0,
        linewidth=1.6, edgecolor=NAVY, facecolor="none", zorder=10,
    )
    ax.add_patch(rect)

    # Axis labels
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels([SCENARIO_SHORT[s] for s in SCENARIO_ORDER],
                       fontsize=9.5)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels([DRUG_LABELS[d] for d in DRUG_ORDER],
                       fontsize=10)

    ax.set_title(
        "Probability of severe shortage (>30% capacity loss) by drug and "
        "disruption scenario\n"
        "SAPIR-Net Phase 2 Monte Carlo simulation, N = 10,000 iterations",
        pad=14, fontsize=12, color=NAVY, loc="left",
    )

    cbar = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.025,
                        ticks=[0, 25, 50, 75, 100])
    cbar.set_label("Probability of severe shortage (%)",
                   fontsize=9.5, color=NAVY)
    cbar.ax.tick_params(labelsize=9, color=SLATE)
    cbar.outline.set_edgecolor(SLATE)
    cbar.outline.set_linewidth(0.6)

    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)

    fig.text(
        0.013, -0.02,
        "Highlighted cell (Amoxicillin/clavulanate x Scenario C): "
        "the 65.42% probability of severe shortage reflects the dual "
        "upstream chokepoint (6-APA + clavulanic acid) propagating as a "
        "multiplicative joint loss under uniform L1->L2 capacity reduction. "
        "Single-component drugs return 20.89% under the same scenario.",
        fontsize=8.0, color=STEEL, ha="left", va="top", wrap=True,
    )

    fig.tight_layout()
    fig.savefig(FIG1_PATH, dpi=300, bbox_inches="tight",
                facecolor=WHITE, edgecolor="none")
    print(f"Exported: {FIG1_PATH}")
    plt.close(fig)


# ============================================================
# FIGURE 2 - CAPACITY-LOSS DENSITIES (2 x 2 BY DRUG)
# ============================================================
# Each panel overlays Scenario B (oxblood, KDE) and Scenario C
# (ochre, KDE). Scenario A is shown as a deterministic vertical
# line (every iteration has the same loss for A). Baseline is a
# point mass at 0% and is noted in the legend, not drawn. The 30%
# severe-shortage threshold is a slate dashed vertical line. The
# 2 x 2 layout places AMOX (top-left) and AMOX_CLAV (top-right)
# side-by-side so the dual-chokepoint compounding is visible by
# direct comparison.

def make_fig2_densities():
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.4),
                              sharex=True, sharey=False)

    x_grid = np.linspace(0, 100, 600)

    for idx, drug in enumerate(DRUG_ORDER):
        ax = axes[idx // 2, idx % 2]

        # Scenario B KDE
        b_arr = raw_losses[_npz_key("B: Cascading Upstream Shock", drug)] * 100
        kde_b = gaussian_kde(b_arr)
        b_density = kde_b(x_grid)
        ax.fill_between(x_grid, b_density, 0,
                        color=OXBLOOD, alpha=0.28, linewidth=0)
        ax.plot(x_grid, b_density, color=OXBLOOD, linewidth=1.8)

        # Scenario C KDE
        c_arr = raw_losses[_npz_key("C: Logistics Chokepoint", drug)] * 100
        kde_c = gaussian_kde(c_arr)
        c_density = kde_c(x_grid)
        ax.fill_between(x_grid, c_density, 0,
                        color=OCHRE, alpha=0.28, linewidth=0)
        ax.plot(x_grid, c_density, color=OCHRE, linewidth=1.8)

        # Scenario A deterministic vertical line
        a_arr = raw_losses[_npz_key("A: Direct Ban (CN+IN)", drug)] * 100
        a_val = float(a_arr[0])  # deterministic -> every iter identical

        # Establish y-limit from the two KDEs before plotting vertical lines
        ymax = max(b_density.max(), c_density.max()) * 1.18
        ax.set_ylim(0, ymax)

        ax.axvline(x=a_val, color=NAVY, linewidth=1.6,
                   linestyle="-", zorder=5)

        # 30% severe-shortage threshold
        ax.axvline(x=30, color=SLATE, linewidth=1.0,
                   linestyle="--", zorder=4)

        # Scenario A annotation
        ax.text(a_val + 0.9, ymax * 0.93,
                f"A: {a_val:.2f}%",
                fontsize=8.5, color=NAVY, ha="left", va="top")

        # 30% threshold annotation (only in top-left panel to avoid
        # redundancy; the dashed line is self-explanatory in the others)
        if idx == 0:
            ax.text(31.0, ymax * 0.96,
                    "30% threshold",
                    fontsize=8, color=SLATE, ha="left", va="top",
                    rotation=90)

        # Panel title
        ax.set_title(DRUG_LABELS_INLINE[drug],
                     fontsize=11.5, color=NAVY, loc="left", pad=8)
        ax.set_xlim(0, 100)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(length=3, color=SLATE)
        ax.xaxis.set_major_formatter(mtick.PercentFormatter(decimals=0))

    # Common axis labels
    for ax in axes[1, :]:
        ax.set_xlabel("Capacity loss (%)", fontsize=10)
    for ax in axes[:, 0]:
        ax.set_ylabel("Density", fontsize=10)

    # Shared legend at top
    legend_handles = [
        plt.Line2D([0], [0], color=NAVY, linewidth=1.6,
                   label="Scenario A: Direct Ban (CN+IN) -- deterministic"),
        plt.Line2D([0], [0], color=OXBLOOD, linewidth=2.5,
                   label="Scenario B: Cascading Upstream Shock"),
        plt.Line2D([0], [0], color=OCHRE, linewidth=2.5,
                   label="Scenario C: Logistics Chokepoint"),
        plt.Line2D([0], [0], color=SLATE, linewidth=1.0, linestyle="--",
                   label="30% severe-shortage threshold"),
    ]
    fig.legend(handles=legend_handles, loc="upper center",
               bbox_to_anchor=(0.5, 0.98),
               ncol=2, fontsize=9.5, frameon=False)

    # Suptitle
    fig.suptitle(
        "Capacity-loss distributions by drug and disruption scenario\n"
        "SAPIR-Net Phase 2 Monte Carlo simulation, N = 10,000 iterations",
        fontsize=12, color=NAVY, y=1.05, x=0.013, ha="left",
    )

    # Footnote
    fig.text(
        0.013, -0.02,
        "Baseline (no shock) omitted from densities (point mass at 0%). "
        "Side-by-side comparison of amoxicillin (top-left) and "
        "amoxicillin/clavulanate (top-right) makes the dual-chokepoint "
        "compounding visible: under both Scenarios B and C, the combination "
        "drug's loss distribution sits to the right of the single-component "
        "drug. The compounding is most pronounced under Scenario C, where "
        "the multiplicative joint loss rule for the 6-APA and clavulanic-acid "
        "chokepoints widens the distribution and shifts its mean from "
        "22.57% (amoxicillin) to 38.62% (amoxicillin/clavulanate).",
        fontsize=8.0, color=STEEL, ha="left", va="top", wrap=True,
    )

    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(FIG2_PATH, dpi=300, bbox_inches="tight",
                facecolor=WHITE, edgecolor="none")
    print(f"Exported: {FIG2_PATH}")
    plt.close(fig)


# ============================================================
# FIGURE 3 - FEDERAL-POLICY-AWARE SENSITIVITY (2 PANELS)
# ============================================================
# Left panel: Scenario B mean capacity loss across the
#   row_china_exposure band (S1 lower / primary / S2 upper)
#   per drug. Floating ranges with primary midpoint marked.
#   Annotation: P(severe) = 100% across the entire band for
#   all four drugs.
# Right panel: Azithromycin baseline-window sensitivity.
#   Scenarios A, B, C: primary 2020-2024 vs S3 2023-2025
#   mean capacity loss. Demonstrates A is window-sensitive
#   while B is dominated by ROW degradation; C is band-
#   independent.

def make_fig3_sensitivity():
    prim_b = primary_df[primary_df.Scenario == "B: Cascading Upstream Shock"]
    s1_b = sensitivity_df[
        (sensitivity_df.Sensitivity_Label == "S1_band_lower")
        & (sensitivity_df.Scenario == "B: Cascading Upstream Shock")
    ]
    s2_b = sensitivity_df[
        (sensitivity_df.Sensitivity_Label == "S2_band_upper")
        & (sensitivity_df.Scenario == "B: Cascading Upstream Shock")
    ]

    fig, (ax_left, ax_right) = plt.subplots(
        1, 2, figsize=(13.5, 5.6),
        gridspec_kw={"width_ratios": [1.35, 1]},
    )

    # ------------- LEFT PANEL -------------
    y_positions = list(range(len(DRUG_ORDER)))[::-1]  # top -> bottom

    for y, drug in zip(y_positions, DRUG_ORDER):
        s1_val = float(s1_b[s1_b.Drug_ID == drug].iloc[0]["Mean_Capacity_Loss_pct"])
        prim_val = float(prim_b[prim_b.Drug_ID == drug].iloc[0]["Mean_Capacity_Loss_pct"])
        s2_val = float(s2_b[s2_b.Drug_ID == drug].iloc[0]["Mean_Capacity_Loss_pct"])

        # Floating range bar
        ax_left.plot([s1_val, s2_val], [y, y],
                     color=STEEL, linewidth=4.5, solid_capstyle="butt",
                     zorder=2, alpha=0.85)
        # Endpoint markers
        ax_left.scatter([s1_val, s2_val], [y, y],
                        marker="|", s=320, color=NAVY, zorder=3,
                        linewidth=2.5)
        # Primary mid-band marker
        ax_left.scatter([prim_val], [y],
                        marker="D", s=70, color=OXBLOOD, zorder=4,
                        edgecolor=NAVY, linewidth=0.6)

        # Endpoint value annotations (outside the bar)
        ax_left.text(s1_val - 0.6, y - 0.32, f"{s1_val:.2f}%",
                     fontsize=8.5, color=NAVY, ha="right", va="top")
        ax_left.text(s2_val + 0.6, y - 0.32, f"{s2_val:.2f}%",
                     fontsize=8.5, color=NAVY, ha="left", va="top")
        # Primary value annotation (above the diamond)
        ax_left.text(prim_val, y + 0.30, f"{prim_val:.2f}%",
                     fontsize=8.5, color=OXBLOOD, ha="center", va="bottom",
                     fontweight="semibold")

    ax_left.set_yticks(y_positions)
    ax_left.set_yticklabels(
        [DRUG_LABELS_INLINE[d] for d in DRUG_ORDER],
        fontsize=10,
    )
    ax_left.set_xlim(68, 103)
    ax_left.set_ylim(-0.8, len(DRUG_ORDER) - 0.2)
    ax_left.xaxis.set_major_formatter(mtick.PercentFormatter(decimals=0))
    ax_left.set_xlabel("Scenario B mean capacity loss (%)", fontsize=10)
    ax_left.set_title(
        "Band sensitivity for Scenario B:\n"
        "row_china_exposure lower (S1) -- primary -- upper (S2)",
        fontsize=11, color=NAVY, loc="left", pad=8,
    )

    # Robustness annotation (positioned in lower-right corner of panel)
    ax_left.text(
        102.5, -0.55,
        "P(severe shortage) = 100% across\n"
        "the entire band for all four drugs --\n"
        "including the lower endpoint (S1).",
        fontsize=9, color=NAVY, ha="right", va="bottom",
        bbox=dict(boxstyle="round,pad=0.45", facecolor=SAND,
                  edgecolor=SLATE, linewidth=0.6),
    )

    # Legend
    legend_handles_left = [
        plt.Line2D([0], [0], marker="|", color=NAVY, markersize=12,
                   linestyle="", markeredgewidth=2.5,
                   label="S1 / S2 band endpoint"),
        plt.Line2D([0], [0], marker="D", color=OXBLOOD, markersize=8,
                   linestyle="", markeredgecolor=NAVY, markeredgewidth=0.6,
                   label="Primary (mid-band)"),
        plt.Line2D([0], [0], color=STEEL, linewidth=4,
                   label="S1 -- S2 range"),
    ]
    ax_left.legend(handles=legend_handles_left, loc="lower left",
                   fontsize=8.5, frameon=False)

    ax_left.spines["top"].set_visible(False)
    ax_left.spines["right"].set_visible(False)
    ax_left.tick_params(length=3, color=SLATE)

    # ------------- RIGHT PANEL -------------
    prim_azi = primary_df[primary_df.Drug_ID == "AZI"]
    s3_azi = sensitivity_df[
        (sensitivity_df.Sensitivity_Label == "S3_HS294150_subrecent_2023_2025")
        & (sensitivity_df.Drug_ID == "AZI")
    ]

    scenarios_show = [
        "A: Direct Ban (CN+IN)",
        "B: Cascading Upstream Shock",
        "C: Logistics Chokepoint",
    ]
    scenarios_label = [
        "A - Direct Ban\n(CN+IN)",
        "B - Cascading\nUpstream Shock",
        "C - Logistics\nChokepoint",
    ]

    x_pos = np.arange(len(scenarios_show))
    bar_width = 0.36

    primary_vals, s3_vals = [], []
    for s in scenarios_show:
        p = float(prim_azi[prim_azi.Scenario == s].iloc[0]["Mean_Capacity_Loss_pct"])
        sv = float(s3_azi[s3_azi.Scenario == s].iloc[0]["Mean_Capacity_Loss_pct"])
        primary_vals.append(p)
        s3_vals.append(sv)

    bars1 = ax_right.bar(x_pos - bar_width / 2, primary_vals, bar_width,
                          color=STEEL, edgecolor=NAVY, linewidth=0.6,
                          label="Primary (2020-2024 baseline)")
    bars2 = ax_right.bar(x_pos + bar_width / 2, s3_vals, bar_width,
                          color=OXBLOOD, edgecolor=NAVY, linewidth=0.6,
                          label="S3 (2023-2025 baseline, HS 294150 only)")

    for b, v in zip(bars1, primary_vals):
        ax_right.text(b.get_x() + b.get_width() / 2, v + 1.4,
                      f"{v:.2f}%", ha="center", va="bottom",
                      fontsize=8.5, color=NAVY)
    for b, v in zip(bars2, s3_vals):
        ax_right.text(b.get_x() + b.get_width() / 2, v + 1.4,
                      f"{v:.2f}%", ha="center", va="bottom",
                      fontsize=8.5, color=NAVY)

    # Shift annotations
    for i, (p, sv) in enumerate(zip(primary_vals, s3_vals)):
        shift = sv - p
        ymax_local = max(p, sv)
        if abs(shift) > 0.5:
            ax_right.annotate(
                f"+{shift:.2f} pp" if shift > 0 else f"{shift:.2f} pp",
                xy=(i, ymax_local + 8), ha="center", va="bottom",
                fontsize=9, color=OXBLOOD, fontweight="semibold",
            )
        else:
            ax_right.annotate(
                "no shift\n(band-\nindependent)",
                xy=(i, ymax_local + 4), ha="center", va="bottom",
                fontsize=8, color=SLATE, style="italic",
            )

    ax_right.set_xticks(x_pos)
    ax_right.set_xticklabels(scenarios_label, fontsize=9.5)
    ax_right.set_ylim(0, 118)
    ax_right.yaxis.set_major_formatter(mtick.PercentFormatter(decimals=0))
    ax_right.set_ylabel("Mean capacity loss (%)", fontsize=10)
    ax_right.set_title(
        "Baseline-window sensitivity (Azithromycin):\n"
        "HS 294150 = 2020-2024 (primary) vs 2023-2025 (S3)",
        fontsize=11, color=NAVY, loc="left", pad=8,
    )
    ax_right.legend(loc="upper right", fontsize=9, frameon=False)
    ax_right.spines["top"].set_visible(False)
    ax_right.spines["right"].set_visible(False)
    ax_right.tick_params(length=3, color=SLATE)

    # Suptitle
    fig.suptitle(
        "Federal-policy-aware sensitivity analysis\n"
        "SAPIR-Net Phase 2 Monte Carlo simulation, N = 10,000 iterations",
        fontsize=12, color=NAVY, y=1.04, x=0.012, ha="left",
    )

    # Footnote
    fig.text(
        0.012, -0.06,
        "Left: even at the band lower endpoint (S1), the cascading upstream "
        "shock places every drug above the 30% severe-shortage threshold. "
        "The federal-policy implication: closing the upstream Chinese KSM/API "
        "concentration is required to materially attenuate Scenario B; "
        "domestic FDF capacity ramps and foreign-API-sourced reserves do not "
        "reach this layer. Right: azithromycin under a sub-recent 2023-2025 "
        "baseline window for HS 294150 (where the direct China share rose "
        "from 56.83% in 2020-2024 to 73.06% in 2023-2025) -- Scenario A "
        "jumps +19.83 pp because the direct ban removes a larger share of "
        "supply, while Scenario B moves only +2.93 pp because ROW Pareto "
        "degradation already dominates the loss distribution. Scenario C "
        "is unchanged across all sensitivity runs because the log-normal "
        "uniform-capacity-reduction shock is independent of the per-HS "
        "row_china_exposure band.",
        fontsize=8.0, color=STEEL, ha="left", va="top", wrap=True,
    )

    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(FIG3_PATH, dpi=300, bbox_inches="tight",
                facecolor=WHITE, edgecolor="none")
    print(f"Exported: {FIG3_PATH}")
    plt.close(fig)


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    print("SAPIR-Net Phase 2 Module 4A: generating figures...")
    print()
    make_fig1_heatmap()
    make_fig2_densities()
    make_fig3_sensitivity()
    print()
    print("Module 4A complete. All three figures ready for Day 5 Red Team "
          "visual audit and Phase 2 white paper integration.")
