"""
SAPIR-Net Phase 2 Module 3: Monte Carlo Simulation Loop (Antibiotics)
======================================================================
Phase 2 extension of the SAPIR-Net Monte Carlo simulation framework
(Zenodo DOI 10.5281/zenodo.19549343; code DOI 10.5281/zenodo.19548610)
from oncology APIs to antibiotic APIs. Methodology unchanged from
Phase 1 (N = 10,000 iterations, 30%% severe shortage threshold, three
scenarios per Phase 1 typology); application domain extended.

Module 3 integrates:
  - Module 1: empirical edge weights from
    results/phase2_hhi_aggregate_2020_2024.csv (5-year USD aggregates)
  - Module 2: scenario classes ScenarioA, ScenarioB, ScenarioC and
    L2 -> L3 propagation logic (with DUAL CHOKEPOINT multiplicative
    joint capacity loss for AMOX_CLAV)
  - phase2_calibration_table.csv: per-chemical-class
    row_china_exposure bands feeding Scenario B's per-HS Pareto sampler

Baseline window
---------------
The Monte Carlo baseline uses 5-year aggregate (2020-2024) USD trade
values, matching the Socal et al. 2025 JAMA Health Forum cross-check
window. This window is also used by Module 1's
phase2_hhi_aggregate_2020_2024.csv export. Phase 1's Monte Carlo
baseline used 5-year aggregates by the same convention. The Module 1
graph build also stored 8-year (2018-2025) aggregates on the L1->L2
edges; the 5-year window is preferred for the Monte Carlo for three
reasons: (a) matches the calibration target window of the
row_china_exposure source stack; (b) aligns with Socal et al. 2025's
2020-2024 disaggregated supplier table; (c) preserves Phase 1
methodological continuity.

Sensitivity analyses
--------------------
Three sensitivity dimensions, each producing a parallel results
table with the same scenario x drug x statistic schema:

  S1. Band-lower-bound sensitivity. Each HS code's
      row_china_exposure parameter is pinned to the lower endpoint
      of its chemical-class band (degenerate band, no shock variation
      within band; the Pareto sampler returns the constant low value).

  S2. Band-upper-bound sensitivity. Each HS code's parameter is
      pinned to the upper endpoint.

  S3. HS 294150 sub-recent-years baseline. The BASELINE_TRADE_USD
      entries for HS 294150 are recomputed using the 2023-2025
      3-year aggregate instead of the 2020-2024 5-year aggregate.
      Per Day 2 close item: HS 294150's empirical direct China share
      shifted from 15.31%% (2018) through 50.36%% (2020) to 83.10%%
      (2024), and 60.72%% (2025). The chemical-class band of [0.60,
      0.85] is held constant; the question this sensitivity answers
      is whether the band remains conservative when the empirical
      baseline is shifted toward the more recent trade pattern.

The primary results table is generated with the full per-HS band
Pareto sampling, which is the Phase 1 default. The two band-endpoint
sensitivities (S1, S2) probe the range of the parameter, and S3
probes baseline-window sensitivity for the one HS code with the
strongest temporal trend.

Outputs
-------
results/phase2_monte_carlo_primary.csv
    Per-(scenario, drug) summary statistics for the primary run.
    Schema: Scenario, Drug, HS_Dependency, Mean_Capacity_Loss_pct,
    Median_Capacity_Loss_pct, P5_Loss_pct, P95_Loss_pct,
    Prob_Severe_Shortage_pct (where severe = >30%% capacity loss).
    Four scenarios x four drugs = 16 rows.

results/phase2_monte_carlo_sensitivity.csv
    Per-(sensitivity_label, scenario, drug) summary statistics for
    the three sensitivity runs (S1 lower, S2 upper, S3 294150
    sub-recent). Same schema as primary plus Sensitivity_Label
    column.

results/phase2_monte_carlo_raw_losses.npz
    Compressed numpy archive with the full N x drug capacity-loss
    arrays for each scenario in the primary run. Loaded by Module 4A
    (Day 4) for density curves and heat-maps.

results/phase2_monte_carlo_audit.txt
    Run log with seed, configuration, summary tables, and the
    structural findings (drug vulnerability matrix; band-endpoint
    range; 294150 baseline-shift effect).

Dependencies: numpy, scipy, pandas

Author: Lead Data Architect / Independent Researcher
Version: Phase 2 Module 3 v1.0, May 2026
"""

import os
import csv
import numpy as np
import pandas as pd

# Module 2 import: scenario classes, propagation logic, ROW China
# exposure bands. Module 2 is in the same working directory.
from sapir_phase2_module2_disruption_engine import (
    ScenarioA, ScenarioB, ScenarioC,
    DRUG_DEPENDENCIES, ROW_CHINA_EXPOSURE,
    propagate_l2_to_l3,
    RNG_SEED,
)

# ============================================================
# CONFIGURATION
# ============================================================

DATA_DIR = "data"
RESULTS_DIR = "results"
os.makedirs(RESULTS_DIR, exist_ok=True)

RAW_COMTRADE_PATH = os.path.join(DATA_DIR,
                                  "sapir_phase2_raw_comtrade_20182025.csv")
HHI_AGGREGATE_PATH = os.path.join(RESULTS_DIR,
                                   "phase2_hhi_aggregate_2020_2024.csv")

# Monte Carlo iteration count -- Phase 1 convention
N_ITER = 10_000

# Severe shortage threshold -- Phase 1 convention
SEVERE_SHORTAGE_THRESHOLD = 0.30

# Reset RNG seed at top of Module 3 to ensure Module 3 starts from
# the canonical state regardless of whether Module 2 audit consumed
# RNG draws on import. (Module 2's audit_distributions() is only
# called via its __main__ block, not on import, but this defensive
# reset preserves reproducibility under any import order.)
np.random.seed(RNG_SEED)

PHASE2_HS_CODES = ["294110", "294130", "294150", "294190"]


# ============================================================
# SECTION 1: LOAD BASELINE TRADE FLOWS (5-YEAR 2020-2024 AGGREGATE)
# ============================================================
# Reads Module 1's phase2_hhi_aggregate_2020_2024.csv to populate
# BASELINE_TRADE_USD. The per-(L1_source, HS) edge weights are
# reconstructed from the per-HS aggregate USD totals plus the
# per-HS China_USD_5yr and India_USD_5yr columns; the ROW edge is
# the residual.

def load_baseline_trade_usd(hhi_agg_path: str) -> tuple:
    """
    Returns:
        baseline_edges: dict keyed by (L1_source, HS_code) -> USD
        baseline_totals: dict keyed by HS_code -> total USD across L1
    """
    df = pd.read_csv(hhi_agg_path)
    # Filter to per-HS rows (exclude the ALL_2941xx aggregate row)
    df = df[df.HS_Code != "ALL_2941xx"].copy()
    # HS_Code can be int after read; normalize to str
    df["HS_Code"] = df["HS_Code"].astype(str)

    baseline_edges = {}
    baseline_totals = {}
    for _, row in df.iterrows():
        hs = row["HS_Code"]
        cn = float(row["China_USD_5yr"])
        ind = float(row["India_USD_5yr"])
        rw = float(row["ROW_USD_5yr"])
        baseline_edges[("CN", hs)] = cn
        baseline_edges[("IN", hs)] = ind
        baseline_edges[("ROW", hs)] = rw
        baseline_totals[hs] = cn + ind + rw

    return baseline_edges, baseline_totals


BASELINE_TRADE_USD, BASELINE_TOTAL = load_baseline_trade_usd(HHI_AGGREGATE_PATH)


# ============================================================
# SECTION 2: HS 294150 SUB-RECENT-YEARS BASELINE (2023-2025)
# ============================================================
# Sensitivity S3 per Day 2 close item: recompute HS 294150 baseline
# using the 2023-2025 3-year window. Direct read from the raw
# Comtrade CSV; uses the same Phase 1 csv.reader trim-to-header
# logic as Module 1.

def compute_subrecent_baseline_294150(raw_path: str,
                                       start_year: int = 2023,
                                       end_year: int = 2025) -> tuple:
    """
    Returns per-L1-source aggregated USD trade values for HS 294150
    over [start_year, end_year]. Same parsing as Module 1.
    """
    with open(raw_path, newline="", encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader)
        rows = [r[: len(header)] for r in reader]
    raw = pd.DataFrame(rows, columns=header)

    df = pd.DataFrame({
        "Year": pd.to_numeric(raw["refYear"]),
        "Source_Country": raw["partnerDesc"].str.strip(),
        "HS_Code": raw["cmdCode"].str.strip(),
        "Trade_Value_USD": pd.to_numeric(raw["primaryValue"],
                                          errors="coerce").fillna(0),
    })

    sub = df[
        (df.HS_Code == "294150")
        & (df.Year >= start_year)
        & (df.Year <= end_year)
    ]

    world = sub.loc[sub.Source_Country == "World", "Trade_Value_USD"].sum()
    china = sub.loc[sub.Source_Country == "China", "Trade_Value_USD"].sum()
    india = sub.loc[sub.Source_Country == "India", "Trade_Value_USD"].sum()
    row_val = world - china - india

    return {
        ("CN", "294150"): float(china),
        ("IN", "294150"): float(india),
        ("ROW", "294150"): float(row_val),
    }, float(world)


SUBRECENT_294150_EDGES, SUBRECENT_294150_TOTAL = (
    compute_subrecent_baseline_294150(RAW_COMTRADE_PATH, 2023, 2025)
)


# ============================================================
# SECTION 3: SIMULATION ENGINE
# ============================================================

def compute_remaining_by_hs(edge_weights: dict) -> dict:
    """Sum remaining USD per HS code across all L1 sources."""
    remaining = {hs: 0.0 for hs in PHASE2_HS_CODES}
    for (source, hs), val in edge_weights.items():
        if hs in remaining:
            remaining[hs] += val
    return remaining


def run_scenario_baseline(n_iter: int, baseline_edges: dict,
                          baseline_totals: dict, bands_by_hs: dict):
    """
    Baseline (no shock). Loss = 0 deterministically.
    bands_by_hs is unused but kept for signature consistency.
    """
    losses = {drug: np.zeros(n_iter) for drug in DRUG_DEPENDENCIES}
    return losses


def run_scenario_a(n_iter: int, baseline_edges: dict,
                   baseline_totals: dict, bands_by_hs: dict):
    """
    Scenario A: CN + IN direct exports zeroed across all HS codes.
    Deterministic. bands_by_hs is unused.
    """
    losses = {drug: np.zeros(n_iter) for drug in DRUG_DEPENDENCIES}

    shocked = ScenarioA.apply(baseline_edges)
    remaining = compute_remaining_by_hs(shocked)
    drug_cap = propagate_l2_to_l3(remaining, baseline_totals)

    for drug in DRUG_DEPENDENCIES:
        loss_frac = 1.0 - drug_cap[drug]
        losses[drug][:] = loss_frac

    return losses


def run_scenario_b(n_iter: int, baseline_edges: dict,
                   baseline_totals: dict, bands_by_hs: dict):
    """
    Scenario B: CN zeroed across all HS codes + ROW Pareto degradation
    sampled per HS code per iteration.

    bands_by_hs: dict keyed by HS_code, values are (low, high) tuples.
                 If low == high, the sampler returns a constant
                 (degenerate band, used by S1/S2 sensitivity runs).
    """
    losses = {drug: np.zeros(n_iter) for drug in DRUG_DEPENDENCIES}

    # Pre-sample per-HS ROW degradation arrays (length n_iter each)
    degradation_arrays = {}
    for hs in PHASE2_HS_CODES:
        low, high = bands_by_hs[hs]
        degradation_arrays[hs] = ScenarioB.sample_row_degradation(
            hs, n_samples=n_iter, band_override=(low, high)
        )

    for i in range(n_iter):
        # Snapshot baseline; apply CN ban and ROW degradation
        weights = baseline_edges.copy()
        # CN zeroed (correlated across HS)
        for key in weights:
            if key[0] == "CN":
                weights[key] = 0.0
        # ROW degraded by per-HS Pareto-rescaled fraction (independent
        # across HS per iteration)
        for key in list(weights.keys()):
            source, hs = key
            if source == "ROW":
                weights[key] = weights[key] * (1.0 - degradation_arrays[hs][i])

        remaining = compute_remaining_by_hs(weights)
        drug_cap = propagate_l2_to_l3(remaining, baseline_totals)
        for drug in DRUG_DEPENDENCIES:
            losses[drug][i] = 1.0 - drug_cap[drug]

    return losses


def run_scenario_c(n_iter: int, baseline_edges: dict,
                   baseline_totals: dict, bands_by_hs: dict):
    """
    Scenario C: Log-normal capacity reduction applied uniformly to ALL
    L1 -> L2 edges. bands_by_hs unused.
    """
    losses = {drug: np.zeros(n_iter) for drug in DRUG_DEPENDENCIES}

    reductions = ScenarioC.sample_capacity_reduction(n_iter)

    for i in range(n_iter):
        # Uniform reduction multiplies every edge equally
        r = reductions[i]
        # remaining_HS = (1 - r) * baseline_total_HS for each HS
        # p_HS = (1 - r) for every HS
        # Single-input drugs: capacity = (1 - r)
        # Dual-input drugs (multiplicative): capacity = (1 - r)^2
        p_uniform = 1.0 - r
        for drug, spec in DRUG_DEPENDENCIES.items():
            if spec["rule"] == "single":
                losses[drug][i] = 1.0 - p_uniform
            elif spec["rule"] == "multiplicative_joint":
                k = len(spec["chokepoints"])
                losses[drug][i] = 1.0 - (p_uniform ** k)

    return losses


SCENARIOS = {
    "Baseline (No Shock)":           run_scenario_baseline,
    "A: Direct Ban (CN+IN)":         run_scenario_a,
    "B: Cascading Upstream Shock":   run_scenario_b,
    "C: Logistics Chokepoint":       run_scenario_c,
}


# ============================================================
# SECTION 4: ONE-RUN HELPER
# ============================================================

def run_all_scenarios(baseline_edges: dict, baseline_totals: dict,
                       bands_by_hs: dict, n_iter: int = N_ITER):
    """
    Execute all four scenarios under a given baseline + band configuration.

    The RNG seed is reset at the start of each call so that sensitivity
    runs are exactly comparable: any difference in outputs across
    runs is attributable to the configuration change (band or
    baseline), not to RNG state drift between runs. In particular,
    Scenario C's log-normal draws (which do not depend on the
    row_china_exposure bands) produce identical loss distributions
    across the primary run and all sensitivity runs.

    Returns:
        results_df:   per-(scenario, drug) summary statistics DataFrame
        loss_arrays:  dict keyed by scenario name, values are dicts
                      keyed by drug -> np.ndarray of capacity losses
    """
    np.random.seed(RNG_SEED)
    all_rows = []
    loss_arrays = {}

    for scenario_name, runner in SCENARIOS.items():
        losses = runner(n_iter, baseline_edges, baseline_totals, bands_by_hs)
        loss_arrays[scenario_name] = losses

        for drug in DRUG_DEPENDENCIES:
            arr = losses[drug]
            mean_l = arr.mean()
            med_l = np.median(arr)
            p5 = np.percentile(arr, 5)
            p95 = np.percentile(arr, 95)
            prob_sev = (arr > SEVERE_SHORTAGE_THRESHOLD).mean()
            hs_dep = ", ".join(DRUG_DEPENDENCIES[drug]["chokepoints"])
            all_rows.append({
                "Scenario":                 scenario_name,
                "Drug":                     DRUG_DEPENDENCIES[drug]["label"],
                "Drug_ID":                  drug,
                "HS_Dependency":            hs_dep,
                "Mean_Capacity_Loss_pct":   round(mean_l * 100, 2),
                "Median_Capacity_Loss_pct": round(med_l * 100, 2),
                "P5_Loss_pct":              round(p5 * 100, 2),
                "P95_Loss_pct":             round(p95 * 100, 2),
                "Prob_Severe_Shortage_pct": round(prob_sev * 100, 2),
            })

    return pd.DataFrame(all_rows), loss_arrays


# ============================================================
# SECTION 5: BAND CONFIGURATIONS
# ============================================================

def midpoint_bands() -> dict:
    """Each HS gets its full chemical-class band [low, high]."""
    return {hs: (ROW_CHINA_EXPOSURE[hs]["low"],
                 ROW_CHINA_EXPOSURE[hs]["high"])
            for hs in PHASE2_HS_CODES}


def lower_pinned_bands() -> dict:
    """Each HS pinned to band lower endpoint (degenerate)."""
    return {hs: (ROW_CHINA_EXPOSURE[hs]["low"],
                 ROW_CHINA_EXPOSURE[hs]["low"])
            for hs in PHASE2_HS_CODES}


def upper_pinned_bands() -> dict:
    """Each HS pinned to band upper endpoint (degenerate)."""
    return {hs: (ROW_CHINA_EXPOSURE[hs]["high"],
                 ROW_CHINA_EXPOSURE[hs]["high"])
            for hs in PHASE2_HS_CODES}


# ============================================================
# SECTION 6: MAIN
# ============================================================

if __name__ == "__main__":
    audit_lines = []

    def alog(s=""):
        audit_lines.append(s)
        print(s)

    alog("=" * 75)
    alog("SAPIR-Net Phase 2 Module 3: Monte Carlo Simulation")
    alog(f"N = {N_ITER:,} iterations  |  severe shortage threshold = "
         f"{SEVERE_SHORTAGE_THRESHOLD:.0%} capacity loss")
    alog(f"RNG seed = {RNG_SEED} (Phase 1 reproducibility convention)")
    alog("=" * 75)
    alog("")
    alog("Baseline window: 2020-2024 5-year USD trade aggregate")
    alog("Source: results/phase2_hhi_aggregate_2020_2024.csv (Module 1 output)")
    alog("")
    alog("BASELINE_TRADE_USD (per (L1_source, HS_code)):")
    for hs in PHASE2_HS_CODES:
        alog(f"  HS {hs}:  CN ${BASELINE_TRADE_USD[('CN', hs)]:>14,.0f}  "
             f"|  IN ${BASELINE_TRADE_USD[('IN', hs)]:>14,.0f}  "
             f"|  ROW ${BASELINE_TRADE_USD[('ROW', hs)]:>14,.0f}  "
             f"|  Total ${BASELINE_TOTAL[hs]:>14,.0f}")
    alog("")

    # ============================================================
    # PRIMARY RUN: full per-HS band Pareto sampling
    # ============================================================
    alog("-" * 75)
    alog("PRIMARY RUN: full per-HS band Pareto sampling")
    alog("-" * 75)
    alog("Bands:")
    for hs in PHASE2_HS_CODES:
        b = ROW_CHINA_EXPOSURE[hs]
        alog(f"  HS {hs} ({b['label']:<45s}): "
             f"[{b['low']:.0%}, {b['high']:.0%}]")
    alog("")

    primary_df, primary_losses = run_all_scenarios(
        BASELINE_TRADE_USD, BASELINE_TOTAL,
        midpoint_bands(),  # full band for Pareto sampling
        n_iter=N_ITER,
    )

    alog("Primary run -- full results table:")
    alog(primary_df.to_string(index=False))
    alog("")

    # Drug vulnerability matrix
    alog("L3 DRUG VULNERABILITY MATRIX  (Probability of Severe Shortage, %%)")
    pivot_prob = primary_df.pivot_table(
        index="Drug", columns="Scenario",
        values="Prob_Severe_Shortage_pct", aggfunc="first",
    )
    col_order = [c for c in SCENARIOS.keys() if c in pivot_prob.columns]
    pivot_prob = pivot_prob[col_order]
    alog(pivot_prob.to_string())
    alog("")

    alog("MEAN CAPACITY LOSS MATRIX  (%%)")
    pivot_mean = primary_df.pivot_table(
        index="Drug", columns="Scenario",
        values="Mean_Capacity_Loss_pct", aggfunc="first",
    )
    pivot_mean = pivot_mean[col_order]
    alog(pivot_mean.to_string())
    alog("")

    # ============================================================
    # SENSITIVITY S1: band lower endpoint pinned (all HS)
    # ============================================================
    alog("-" * 75)
    alog("SENSITIVITY S1: row_china_exposure pinned to lower band endpoint")
    alog("-" * 75)
    alog("Each HS band degenerate at its lower endpoint:")
    for hs in PHASE2_HS_CODES:
        alog(f"  HS {hs}: {ROW_CHINA_EXPOSURE[hs]['low']:.0%}")
    alog("")

    s1_df, _ = run_all_scenarios(
        BASELINE_TRADE_USD, BASELINE_TOTAL,
        lower_pinned_bands(),
        n_iter=N_ITER,
    )
    s1_df.insert(0, "Sensitivity_Label", "S1_band_lower")
    alog(s1_df.to_string(index=False))
    alog("")

    # ============================================================
    # SENSITIVITY S2: band upper endpoint pinned (all HS)
    # ============================================================
    alog("-" * 75)
    alog("SENSITIVITY S2: row_china_exposure pinned to upper band endpoint")
    alog("-" * 75)
    alog("Each HS band degenerate at its upper endpoint:")
    for hs in PHASE2_HS_CODES:
        alog(f"  HS {hs}: {ROW_CHINA_EXPOSURE[hs]['high']:.0%}")
    alog("")

    s2_df, _ = run_all_scenarios(
        BASELINE_TRADE_USD, BASELINE_TOTAL,
        upper_pinned_bands(),
        n_iter=N_ITER,
    )
    s2_df.insert(0, "Sensitivity_Label", "S2_band_upper")
    alog(s2_df.to_string(index=False))
    alog("")

    # ============================================================
    # SENSITIVITY S3: HS 294150 sub-recent-years baseline (2023-2025)
    # ============================================================
    alog("-" * 75)
    alog("SENSITIVITY S3: HS 294150 baseline = 2023-2025 (3-year window)")
    alog("-" * 75)
    alog("All other HS codes keep 2020-2024 baseline.")
    alog("All bands remain at chemical-class values (full per-HS Pareto).")
    alog("")
    alog("HS 294150 sub-recent baseline (2023-2025):")
    alog(f"  CN  ${SUBRECENT_294150_EDGES[('CN', '294150')]:>14,.0f}  |  "
         f"IN  ${SUBRECENT_294150_EDGES[('IN', '294150')]:>14,.0f}  |  "
         f"ROW ${SUBRECENT_294150_EDGES[('ROW', '294150')]:>14,.0f}  |  "
         f"Total ${SUBRECENT_294150_TOTAL:>14,.0f}")
    cn_share_subrecent = (SUBRECENT_294150_EDGES[('CN', '294150')] /
                           SUBRECENT_294150_TOTAL * 100)
    cn_share_primary = (BASELINE_TRADE_USD[('CN', '294150')] /
                         BASELINE_TOTAL['294150'] * 100)
    alog(f"  Direct China share (sub-recent 2023-2025): "
         f"{cn_share_subrecent:.2f}%")
    alog(f"  Direct China share (primary 2020-2024):    "
         f"{cn_share_primary:.2f}%")
    alog(f"  Shift toward CN-direct exposure: "
         f"{cn_share_subrecent - cn_share_primary:+.2f} pp")
    alog("")

    # Construct alternative baseline edges + totals: copy primary,
    # override HS 294150 entries with sub-recent values.
    s3_edges = BASELINE_TRADE_USD.copy()
    s3_totals = BASELINE_TOTAL.copy()
    for k, v in SUBRECENT_294150_EDGES.items():
        s3_edges[k] = v
    s3_totals["294150"] = SUBRECENT_294150_TOTAL

    s3_df, _ = run_all_scenarios(
        s3_edges, s3_totals,
        midpoint_bands(),  # full band Pareto (same as primary)
        n_iter=N_ITER,
    )
    s3_df.insert(0, "Sensitivity_Label", "S3_HS294150_subrecent_2023_2025")
    alog(s3_df.to_string(index=False))
    alog("")

    # ============================================================
    # OUTPUTS
    # ============================================================
    primary_path = os.path.join(RESULTS_DIR, "phase2_monte_carlo_primary.csv")
    primary_df.to_csv(primary_path, index=False)

    sensitivity_df = pd.concat([s1_df, s2_df, s3_df], ignore_index=True)
    sensitivity_path = os.path.join(RESULTS_DIR,
                                     "phase2_monte_carlo_sensitivity.csv")
    sensitivity_df.to_csv(sensitivity_path, index=False)

    # Raw loss arrays (for Module 4A density curves / heat-maps)
    raw_path = os.path.join(RESULTS_DIR,
                             "phase2_monte_carlo_raw_losses.npz")
    raw_dict = {}
    for scen_name, drug_losses in primary_losses.items():
        # numpy npz keys cannot contain "(", ")", ":", spaces well-handled;
        # sanitize.
        scen_key = (scen_name.replace(" ", "_")
                    .replace(":", "")
                    .replace("(", "")
                    .replace(")", "")
                    .replace("+", "and"))
        for drug, arr in drug_losses.items():
            key = f"{scen_key}__{drug}"
            raw_dict[key] = arr
    np.savez_compressed(raw_path, **raw_dict)

    alog("-" * 75)
    alog("EXPORTS")
    alog("-" * 75)
    alog(f"  {primary_path}")
    alog(f"  {sensitivity_path}")
    alog(f"  {raw_path}")
    alog("")

    # ============================================================
    # KEY FINDINGS SUMMARY (BLUE TEAM SELF-CHECK ONLY)
    # ============================================================
    # Per Brief Section 2 Red Team protocol: this chat is Blue Team.
    # The block below is the build-completion status line. The Day 5
    # Red Team session in a separate chat will audit Module 3
    # alongside the rest of the build.

    alog("=" * 75)
    alog("MODULE 3 BUILD STATUS")
    alog("=" * 75)
    alog(f"  Primary run rows produced:                {len(primary_df)} "
         f"(4 scenarios x 4 drugs)")
    alog(f"  Sensitivity run rows produced:            {len(sensitivity_df)} "
         f"(3 sensitivities x 4 scenarios x 4 drugs)")
    alog(f"  Raw loss arrays saved (scenarios x drugs): "
         f"{len(raw_dict)} = {len(primary_losses)} x {len(DRUG_DEPENDENCIES)}")
    alog("")
    alog("Scenario B (Cascading Upstream Shock) primary-run probabilities of")
    alog(f"severe shortage (loss > {SEVERE_SHORTAGE_THRESHOLD:.0%}):")
    scenB = primary_df[primary_df.Scenario == "B: Cascading Upstream Shock"]
    for _, row in scenB.iterrows():
        alog(f"  {row['Drug']:<45s}  P(severe) = "
             f"{row['Prob_Severe_Shortage_pct']:6.2f}%  "
             f"|  mean loss = {row['Mean_Capacity_Loss_pct']:6.2f}%")
    alog("")
    alog("Module 3 complete. Outputs ready for Module 4A (Day 4 build:")
    alog("visualizations) and downstream Phase 2 white paper draft.")
    alog("=" * 75)

    # Write audit log
    audit_path = os.path.join(RESULTS_DIR, "phase2_monte_carlo_audit.txt")
    with open(audit_path, "w", encoding="utf-8") as f:
        f.write("\n".join(audit_lines))
    print(f"\nAudit log written: {audit_path}")
