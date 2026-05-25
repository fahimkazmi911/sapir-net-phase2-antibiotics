"""
SAPIR-Net Phase 2 Module 3 — Partial-Correlation Extension Runner
==================================================================
Day 8 Blue Team v1.2 polish addition. Addresses Day 7 SME Peer-Review
Audit finding M-4 (option a): the multiplicative joint capacity-loss
rule applied to AMOX_CLAV in Scenario B assumes statistical
independence of the per-iteration ROW Pareto degradation on the two
upstream chemical inputs PENI (HS 294110) and OTHAB (HS 294190).
The independence assumption is asserted in v1.1 §2.5.3 and §2.8 to
be conservative under the realistic case that the 6-APA and
clavulanic acid fermentation production bases share Chinese-cluster
upstream infrastructure. The present extension quantifies that
conservatism via a Gaussian copula on the per-iteration PENI and
OTHAB Pareto draws, at four correlation values rho in {0, 0.25,
0.50, 0.75}.

Implementation approach
-----------------------
1. Draw bivariate standard normal pairs (Z1, Z2) with Pearson
   correlation rho via Cholesky-decomposition construction
   (Z2 = rho * Z1 + sqrt(1 - rho^2) * Z_aux).
2. Transform to U[0,1] via the standard normal CDF Phi(.).
3. Transform each U to the rescaled-Pareto marginal via the
   inverse CDF y = 1 - (1-U)^(1/alpha), alpha = 2.5. This is
   the same marginal distribution Module 2's
   ScenarioB.sample_row_degradation produces; verification by
   simulation in audit log.
4. Rescale to the chemical-class band [low, high] for PENI and
   OTHAB respectively.
5. Apply the Scenario B propagation logic (CN ban + ROW
   degradation) with the new correlated (PENI, OTHAB) pair, and
   with independent TETRA and ERYTH draws (single-input drugs
   are unaffected by the PENI-OTHAB joint structure).
6. Compute per-drug capacity-loss summary statistics over
   N = 10,000 iterations per rho.

Design choices documented for Day 9 Red Team verification
---------------------------------------------------------
- Module 3 itself is NOT modified. The Day 4 byte-exact close on
  phase2_monte_carlo_primary.csv, phase2_monte_carlo_sensitivity.csv,
  and phase2_monte_carlo_raw_losses.npz is preserved. The partial-
  correlation extension is a separate runner script that imports
  Module 3's state and adds new outputs.

- The new sampler uses numpy's default_rng (PCG64) seeded with 42.
  Module 3's primary-run sampler uses scipy.stats.pareto.rvs against
  the legacy np.random global state seeded with 42. The two RNG
  paths are different; consequently the new sampler at rho = 0 does
  NOT byte-exact reproduce the existing primary-run Scenario B
  values. Instead, the new sampler at rho = 0 is statistically
  equivalent to the existing primary-run (same marginal distributions,
  same independence structure); the comparison is reported in the
  audit log as a sampling-variability check at N = 10,000.

- Only the AMOX_CLAV Scenario B distribution is affected by the
  PENI-OTHAB correlation. The three single-input drugs (AMOX uses
  PENI only, AZI uses ERYTH only, DOXY uses TETRA only) under
  Scenario B have only one chokepoint each and are by construction
  invariant to the PENI-OTHAB joint structure. They appear in the
  partial-correlation CSV with the rho = 0 statistical-equivalence
  values (sampling variability around the primary-run values).

- Per the v1.2 M-1 reframing: the binary P(severe shortage) outcome
  is band-arithmetic-driven and remains 100% across all rho values
  for every drug in Scenario B. The partial-correlation sensitivity
  therefore quantifies the MAGNITUDE distribution (mean, median,
  P5, P95) rather than the binary indicator.

Outputs
-------
results/phase2_monte_carlo_correlation.csv
    Per-(rho, drug) summary statistics for the partial-correlation
    Scenario B run.
    Schema: Rho, Drug, Drug_ID, HS_Dependency,
        Mean_Capacity_Loss_pct, Median_Capacity_Loss_pct,
        P5_Loss_pct, P95_Loss_pct, Prob_Severe_Shortage_pct.
    4 rho values x 4 drugs = 16 rows.

results/phase2_monte_carlo_correlation_audit.txt
    Run log with RNG seed, configuration, marginal-distribution
    validation, joint-correlation validation, the AMOX_CLAV
    magnitude-distribution shifts, and the sampling-variability
    check against the existing primary-run rho = 0 values.

results/phase2_monte_carlo_correlation_amox_clav_losses.npz
    Compressed numpy archive with the AMOX_CLAV capacity-loss
    arrays at each rho value (used by Figure 4).

Author: Lead Data Architect / Independent Researcher
Version: Phase 2 Module 3 Correlation Extension v1.0, May 2026
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import norm, pearsonr

# Import existing Module 2/3 state. Module 3's import populates
# BASELINE_TRADE_USD, BASELINE_TOTAL, etc. by reading the Module 1
# aggregate CSV from results/.
from sapir_phase2_module2_disruption_engine import (
    DRUG_DEPENDENCIES, ROW_CHINA_EXPOSURE, ScenarioB,
    propagate_l2_to_l3, RNG_SEED,
)
from sapir_phase2_module3_monte_carlo import (
    BASELINE_TRADE_USD, BASELINE_TOTAL, N_ITER,
    SEVERE_SHORTAGE_THRESHOLD, PHASE2_HS_CODES,
    compute_remaining_by_hs, RESULTS_DIR,
)


# ============================================================
# CONFIGURATION
# ============================================================

# Correlation values to evaluate. rho = 0 is the independence
# baseline (statistically matching the existing Module 3 primary
# Scenario B); rho ∈ {0.25, 0.50, 0.75} represent progressively
# more realistic shared-infrastructure scenarios.
RHO_VALUES = [0.0, 0.25, 0.50, 0.75]

PARETO_ALPHA = ScenarioB.PARETO_ALPHA  # = 2.5

# HS code roles for AMOX_CLAV multiplicative joint rule
HS_PENI = "294110"
HS_OTHAB = "294190"
HS_TETRA = "294130"
HS_ERYTH = "294150"


# ============================================================
# CORE SAMPLER
# ============================================================

def sample_correlated_rescaled_pareto_pair(low_a: float, high_a: float,
                                            low_b: float, high_b: float,
                                            n_samples: int, rho: float,
                                            rng: np.random.Generator) -> tuple:
    """
    Sample n_samples pairs (deg_a, deg_b) where each marginal matches the
    rescaled-Pareto(alpha=2.5) distribution on its respective band and
    the joint dependence is a Gaussian copula at Pearson correlation rho
    in the underlying bivariate normal.

    Marginal distribution check
    ---------------------------
    The existing Module 2 sampler draws X ~ Pareto(alpha) on [1, infty)
    and computes Y = 1 - 1/X. The CDF of Y is G(y) = 1 - (1-y)^alpha on
    [0, 1) (a Beta(1, alpha) distribution). The inverse CDF is
        G^{-1}(u) = 1 - (1-u)^(1/alpha)
    for u in [0, 1). This routine draws U ~ Uniform[0, 1) via a Gaussian
    copula (correlated standard normals transformed by the normal CDF)
    and applies G^{-1} to each U component. The resulting per-axis
    marginal is identical in distribution to the existing sampler.

    Joint dependence
    ----------------
    The Gaussian copula has Spearman rank correlation rho_S =
    (6/pi) * arcsin(rho/2) ~ rho for moderate rho. Pearson correlation
    in the rescaled-Pareto marginals depends on the marginals; reported
    in the audit log.

    Degenerate-band handling
    ------------------------
    If both bands degenerate (low == high), returns the constants.
    This shouldn't occur in this extension since AMOX_CLAV's PENI and
    OTHAB bands are nondegenerate, but handled for safety.
    """
    # Degenerate handling
    if low_a == high_a and low_b == high_b:
        return np.full(n_samples, low_a), np.full(n_samples, low_b)

    # Bivariate standard normal with Pearson correlation rho:
    #   Z1 ~ N(0,1); Z_aux ~ N(0,1) independent of Z1;
    #   Z2 = rho * Z1 + sqrt(1 - rho^2) * Z_aux
    z1 = rng.standard_normal(n_samples)
    if rho == 0.0:
        z2 = rng.standard_normal(n_samples)
    else:
        z_aux = rng.standard_normal(n_samples)
        z2 = rho * z1 + np.sqrt(1.0 - rho * rho) * z_aux

    # Transform to U[0,1] via standard normal CDF
    u1 = norm.cdf(z1)
    u2 = norm.cdf(z2)

    # Rescaled-Pareto inverse CDF: y = 1 - (1-u)^(1/alpha)
    normalized_1 = 1.0 - (1.0 - u1) ** (1.0 / PARETO_ALPHA)
    normalized_2 = 1.0 - (1.0 - u2) ** (1.0 / PARETO_ALPHA)

    # Rescale to band
    deg_a = low_a + normalized_1 * (high_a - low_a)
    deg_b = low_b + normalized_2 * (high_b - low_b)

    # Hard clamp for safety (matches Module 2 convention)
    deg_a = np.clip(deg_a, low_a, high_a)
    deg_b = np.clip(deg_b, low_b, high_b)

    return deg_a, deg_b


def sample_independent_rescaled_pareto(low: float, high: float,
                                        n_samples: int,
                                        rng: np.random.Generator) -> np.ndarray:
    """Single-axis rescaled-Pareto draw using the same RNG generator as
    the correlated pair sampler (for the single-input HS codes TETRA
    and ERYTH in the partial-correlation Scenario B run).
    """
    if low == high:
        return np.full(n_samples, low)
    z = rng.standard_normal(n_samples)
    u = norm.cdf(z)
    normalized = 1.0 - (1.0 - u) ** (1.0 / PARETO_ALPHA)
    deg = low + normalized * (high - low)
    return np.clip(deg, low, high)


# ============================================================
# SCENARIO B WITH PARTIAL-CORRELATION ON (PENI, OTHAB)
# ============================================================

def run_scenario_b_correlated(baseline_edges: dict,
                               baseline_totals: dict,
                               rho: float,
                               n_iter: int = N_ITER,
                               seed: int = RNG_SEED) -> tuple:
    """
    Scenario B execution with PENI (HS 294110) and OTHAB (HS 294190)
    ROW Pareto draws coupled by a Gaussian copula at correlation rho.
    TETRA (HS 294130) and ERYTH (HS 294150) ROW draws are independent
    (the chemical-class bands are otherwise unchanged from the primary
    Module 3 configuration).

    Returns
    -------
    losses : dict
        Keyed by drug ID (AMOX, AMOX_CLAV, AZI, DOXY); values are
        np.ndarray of capacity losses of length n_iter.
    draws : dict
        Keyed by 'PENI', 'OTHAB', 'TETRA', 'ERYTH'; values are the
        per-HS degradation samples used for the run (length n_iter
        each). Useful for the audit log marginal/joint diagnostics.
    """
    rng = np.random.default_rng(seed)

    band_peni = (ROW_CHINA_EXPOSURE[HS_PENI]["low"],
                 ROW_CHINA_EXPOSURE[HS_PENI]["high"])
    band_othab = (ROW_CHINA_EXPOSURE[HS_OTHAB]["low"],
                  ROW_CHINA_EXPOSURE[HS_OTHAB]["high"])
    band_tetra = (ROW_CHINA_EXPOSURE[HS_TETRA]["low"],
                  ROW_CHINA_EXPOSURE[HS_TETRA]["high"])
    band_eryth = (ROW_CHINA_EXPOSURE[HS_ERYTH]["low"],
                  ROW_CHINA_EXPOSURE[HS_ERYTH]["high"])

    # Correlated PENI + OTHAB draws (the M-4 partial-correlation extension)
    deg_peni, deg_othab = sample_correlated_rescaled_pareto_pair(
        band_peni[0], band_peni[1],
        band_othab[0], band_othab[1],
        n_iter, rho, rng,
    )

    # Independent TETRA + ERYTH draws (unchanged structure for single-
    # input drugs DOXY and AZI)
    deg_tetra = sample_independent_rescaled_pareto(
        band_tetra[0], band_tetra[1], n_iter, rng)
    deg_eryth = sample_independent_rescaled_pareto(
        band_eryth[0], band_eryth[1], n_iter, rng)

    degradation_by_hs = {
        HS_PENI:  deg_peni,
        HS_OTHAB: deg_othab,
        HS_TETRA: deg_tetra,
        HS_ERYTH: deg_eryth,
    }

    losses = {drug: np.zeros(n_iter) for drug in DRUG_DEPENDENCIES}

    for i in range(n_iter):
        weights = baseline_edges.copy()
        # CN ban (correlated across HS by construction; same as primary B)
        for key in weights:
            if key[0] == "CN":
                weights[key] = 0.0
        # ROW degradation per HS using the new correlated/independent draws
        for key in list(weights.keys()):
            source, hs = key
            if source == "ROW":
                weights[key] = weights[key] * (1.0 - degradation_by_hs[hs][i])

        remaining = compute_remaining_by_hs(weights)
        drug_cap = propagate_l2_to_l3(remaining, baseline_totals)
        for drug in DRUG_DEPENDENCIES:
            losses[drug][i] = 1.0 - drug_cap[drug]

    draws = {
        "PENI":  deg_peni,
        "OTHAB": deg_othab,
        "TETRA": deg_tetra,
        "ERYTH": deg_eryth,
    }

    return losses, draws


# ============================================================
# SUMMARY STATISTICS
# ============================================================

def summarize_losses(losses: dict, rho: float) -> pd.DataFrame:
    """Compute per-drug summary statistics for one rho value."""
    rows = []
    for drug, arr in losses.items():
        mean_l = arr.mean()
        med_l = np.median(arr)
        p5 = np.percentile(arr, 5)
        p95 = np.percentile(arr, 95)
        prob_sev = (arr > SEVERE_SHORTAGE_THRESHOLD).mean()
        hs_dep = ", ".join(DRUG_DEPENDENCIES[drug]["chokepoints"])
        rows.append({
            "Rho":                      rho,
            "Drug":                     DRUG_DEPENDENCIES[drug]["label"],
            "Drug_ID":                  drug,
            "HS_Dependency":            hs_dep,
            "Mean_Capacity_Loss_pct":   round(mean_l * 100, 2),
            "Median_Capacity_Loss_pct": round(med_l * 100, 2),
            "P5_Loss_pct":              round(p5 * 100, 2),
            "P95_Loss_pct":             round(p95 * 100, 2),
            "Prob_Severe_Shortage_pct": round(prob_sev * 100, 2),
        })
    return pd.DataFrame(rows)


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    audit_lines = []

    def alog(s=""):
        audit_lines.append(s)
        print(s)

    alog("=" * 75)
    alog("SAPIR-Net Phase 2 Module 3 — Partial-Correlation Extension")
    alog("(Day 8 Blue Team v1.2 polish; addresses Day 7 SME Peer-Review")
    alog(" Audit finding M-4 option (a): partial-correlation sensitivity")
    alog(" on the AMOX_CLAV dual-chokepoint independence assumption in")
    alog(" Scenario B's per-iteration ROW Pareto sampling.)")
    alog("=" * 75)
    alog("")
    alog(f"N = {N_ITER:,} iterations per rho value")
    alog(f"Severe shortage threshold: {SEVERE_SHORTAGE_THRESHOLD:.0%} capacity loss")
    alog(f"RNG seed: {RNG_SEED} (numpy default_rng / PCG64; Day 4 Module 3")
    alog( "primary run uses scipy pareto.rvs against legacy np.random global")
    alog( "state; the two RNG paths are different, so at rho = 0 the new")
    alog( "sampler does NOT byte-exact reproduce the existing primary-run")
    alog( "AMOX_CLAV Scenario B values but should match within sampling")
    alog( "tolerance at N = 10,000.)")
    alog("")
    alog("Baseline window: 2020-2024 5-year USD trade aggregate (same as")
    alog("Module 3 primary run).")
    alog("")
    alog("Chemical-class bands (unchanged from Module 2 primary configuration):")
    for hs in (HS_PENI, HS_OTHAB, HS_TETRA, HS_ERYTH):
        b = ROW_CHINA_EXPOSURE[hs]
        alog(f"  HS {hs} ({b['label']:<45s}): "
             f"[{b['low']:.0%}, {b['high']:.0%}]")
    alog("")
    alog("Partial-correlation values to evaluate:")
    alog(f"  rho in {RHO_VALUES}")
    alog("  (rho = 0 = independence; rho > 0 = positive Gaussian-copula")
    alog("   correlation between per-iteration PENI and OTHAB ROW Pareto draws)")
    alog("")

    # ============================================================
    # MARGINAL DISTRIBUTION VALIDATION at rho = 0
    # ============================================================
    alog("-" * 75)
    alog("MARGINAL DISTRIBUTION VALIDATION (rho = 0)")
    alog("-" * 75)
    alog("Verify that the new copula-based sampler produces marginal")
    alog("PENI and OTHAB distributions matching the Module 2 reference")
    alog("targets (50,000-sample audit; same target values as")
    alog("phase2_module2_distribution_audit.txt).")
    alog("")

    rng_validation = np.random.default_rng(RNG_SEED)
    n_validation = 50_000
    band_peni_low = ROW_CHINA_EXPOSURE[HS_PENI]["low"]
    band_peni_high = ROW_CHINA_EXPOSURE[HS_PENI]["high"]
    band_othab_low = ROW_CHINA_EXPOSURE[HS_OTHAB]["low"]
    band_othab_high = ROW_CHINA_EXPOSURE[HS_OTHAB]["high"]

    val_peni, val_othab = sample_correlated_rescaled_pareto_pair(
        band_peni_low, band_peni_high,
        band_othab_low, band_othab_high,
        n_validation, 0.0, rng_validation,
    )

    alog(f"  HS 294110 (PENI), 50,000-sample validation at rho = 0:")
    alog(f"    mean   = {val_peni.mean()*100:6.2f}%   "
         f"(Module 2 reference target: 75.7%)")
    alog(f"    median = {np.median(val_peni)*100:6.2f}%   "
         f"(target: 74.8%)")
    alog(f"    P5     = {np.percentile(val_peni, 5)*100:6.2f}%   "
         f"(target: 70.4%)")
    alog(f"    P95    = {np.percentile(val_peni, 95)*100:6.2f}%   "
         f"(target: 83.9%)")
    alog(f"  HS 294190 (OTHAB), 50,000-sample validation at rho = 0:")
    alog(f"    mean   = {val_othab.mean()*100:6.2f}%   "
         f"(Module 2 reference target: 62.2%)")
    alog(f"    median = {np.median(val_othab)*100:6.2f}%   "
         f"(target: 61.1%)")
    alog(f"    P5     = {np.percentile(val_othab, 5)*100:6.2f}%   "
         f"(target: 55.5%)")
    alog(f"    P95    = {np.percentile(val_othab, 95)*100:6.2f}%   "
         f"(target: 72.5%)")
    alog("  Empirical Pearson correlation of validation draws "
         f"(rho_intended = 0):  rho_obs = "
         f"{pearsonr(val_peni, val_othab)[0]:+.4f}")
    alog("")

    # ============================================================
    # MAIN RUN: per-rho Scenario B with partial correlation
    # ============================================================
    alog("-" * 75)
    alog("PARTIAL-CORRELATION SCENARIO B RUNS")
    alog("-" * 75)

    all_results = []
    amox_clav_loss_arrays = {}

    for rho in RHO_VALUES:
        alog("")
        alog(f"### rho = {rho:.2f}")

        losses, draws = run_scenario_b_correlated(
            BASELINE_TRADE_USD, BASELINE_TOTAL,
            rho=rho, n_iter=N_ITER, seed=RNG_SEED,
        )

        # Audit: realized correlation between PENI and OTHAB draws
        rho_realized_pearson = pearsonr(draws["PENI"], draws["OTHAB"])[0]
        alog(f"    PENI–OTHAB realized Pearson correlation: "
             f"{rho_realized_pearson:+.4f} "
             f"(rho_intended = {rho:+.2f}; Gaussian copula in normal")
        alog( "     space; Pearson in rescaled-Pareto marginals is slightly")
        alog( "     attenuated relative to the copula rho by the nonlinear")
        alog( "     marginal CDFs — standard copula behavior)")
        alog(f"    PENI marginal: mean = {draws['PENI'].mean()*100:6.2f}%, "
             f"median = {np.median(draws['PENI'])*100:6.2f}%, "
             f"P95 = {np.percentile(draws['PENI'], 95)*100:6.2f}%")
        alog(f"    OTHAB marginal: mean = {draws['OTHAB'].mean()*100:6.2f}%, "
             f"median = {np.median(draws['OTHAB'])*100:6.2f}%, "
             f"P95 = {np.percentile(draws['OTHAB'], 95)*100:6.2f}%")

        df_rho = summarize_losses(losses, rho)
        all_results.append(df_rho)
        amox_clav_loss_arrays[f"rho_{rho:.2f}"] = losses["AMOX_CLAV"]

        alog("")
        alog("    Per-drug Scenario B summary at this rho:")
        for _, row in df_rho.iterrows():
            alog(f"      {row['Drug']:<45s}  "
                 f"mean = {row['Mean_Capacity_Loss_pct']:6.2f}%  "
                 f"median = {row['Median_Capacity_Loss_pct']:6.2f}%  "
                 f"P5 = {row['P5_Loss_pct']:6.2f}%  "
                 f"P95 = {row['P95_Loss_pct']:6.2f}%  "
                 f"P(sev) = {row['Prob_Severe_Shortage_pct']:6.2f}%")

    # ============================================================
    # CROSS-RHO MAGNITUDE SHIFT TABLE — AMOX_CLAV
    # ============================================================
    alog("")
    alog("-" * 75)
    alog("AMOX_CLAV SCENARIO B MAGNITUDE DISTRIBUTION SHIFTS ACROSS rho")
    alog("-" * 75)
    alog("Per the v1.2 M-1 reframing, the binary P(severe shortage)")
    alog("outcome is band-arithmetic-driven and remains 100% across all")
    alog("rho values. The partial-correlation sensitivity therefore")
    alog("characterizes the MAGNITUDE distribution. The conservatism")
    alog("claim in v1.1 §2.5.3 / §2.8 is upheld iff positive correlation")
    alog("between PENI and OTHAB ROW Pareto draws shifts the AMOX_CLAV")
    alog("magnitude distribution rightward (higher mean / higher P95).")
    alog("")
    alog(f"  {'rho':>6s}  {'mean':>8s}  {'median':>8s}  {'P5':>8s}  "
         f"{'P95':>8s}  {'P(sev)':>8s}")
    for df_rho in all_results:
        ac_row = df_rho[df_rho.Drug_ID == "AMOX_CLAV"].iloc[0]
        alog(f"  {ac_row['Rho']:6.2f}  "
             f"{ac_row['Mean_Capacity_Loss_pct']:7.2f}%  "
             f"{ac_row['Median_Capacity_Loss_pct']:7.2f}%  "
             f"{ac_row['P5_Loss_pct']:7.2f}%  "
             f"{ac_row['P95_Loss_pct']:7.2f}%  "
             f"{ac_row['Prob_Severe_Shortage_pct']:7.2f}%")
    alog("")

    # Magnitude-shift summary
    rho0_amox_clav = all_results[0][all_results[0].Drug_ID == "AMOX_CLAV"].iloc[0]
    alog("Magnitude shift from rho = 0 to rho = 0.75 (AMOX_CLAV Scenario B):")
    for df_rho in all_results[1:]:
        rho_val = df_rho.iloc[0]["Rho"]
        ac_row = df_rho[df_rho.Drug_ID == "AMOX_CLAV"].iloc[0]
        dmean = ac_row["Mean_Capacity_Loss_pct"] - rho0_amox_clav["Mean_Capacity_Loss_pct"]
        dp95 = ac_row["P95_Loss_pct"] - rho0_amox_clav["P95_Loss_pct"]
        dp5 = ac_row["P5_Loss_pct"] - rho0_amox_clav["P5_Loss_pct"]
        alog(f"  rho = {rho_val:.2f}: mean shift = {dmean:+.2f} pp, "
             f"P5 shift = {dp5:+.2f} pp, P95 shift = {dp95:+.2f} pp")
    alog("")

    # ============================================================
    # SAMPLING-VARIABILITY CHECK vs EXISTING PRIMARY RUN at rho = 0
    # ============================================================
    alog("-" * 75)
    alog("SAMPLING-VARIABILITY CHECK: new sampler at rho = 0 vs Module 3")
    alog("primary-run AMOX_CLAV Scenario B")
    alog("-" * 75)
    alog("The existing primary run (phase2_monte_carlo_primary.csv,")
    alog("Scenario B / AMOX_CLAV) reports:")
    alog("  mean = 95.16%, median = 95.10%, P5 = 93.73%, P95 = 96.84%,")
    alog("  P(severe) = 100.00%.")
    alog("")
    alog("The new partial-correlation sampler at rho = 0 should match")
    alog("within sampling variability at N = 10,000:")
    new_rho0 = all_results[0][all_results[0].Drug_ID == "AMOX_CLAV"].iloc[0]
    alog(f"  mean = {new_rho0['Mean_Capacity_Loss_pct']:.2f}%, "
         f"median = {new_rho0['Median_Capacity_Loss_pct']:.2f}%, "
         f"P5 = {new_rho0['P5_Loss_pct']:.2f}%, "
         f"P95 = {new_rho0['P95_Loss_pct']:.2f}%, "
         f"P(severe) = {new_rho0['Prob_Severe_Shortage_pct']:.2f}%.")
    alog("")
    alog("Absolute differences:")
    primary_ref = {"Mean_Capacity_Loss_pct": 95.16,
                   "Median_Capacity_Loss_pct": 95.10,
                   "P5_Loss_pct": 93.73,
                   "P95_Loss_pct": 96.84,
                   "Prob_Severe_Shortage_pct": 100.00}
    for k in ("Mean_Capacity_Loss_pct", "Median_Capacity_Loss_pct",
              "P5_Loss_pct", "P95_Loss_pct", "Prob_Severe_Shortage_pct"):
        diff = abs(new_rho0[k] - primary_ref[k])
        alog(f"  |new - primary| {k:<28s} = {diff:5.2f} pp")
    alog("")

    # ============================================================
    # OUTPUTS
    # ============================================================
    results_df = pd.concat(all_results, ignore_index=True)
    out_csv = os.path.join(RESULTS_DIR, "phase2_monte_carlo_correlation.csv")
    results_df.to_csv(out_csv, index=False)

    out_npz = os.path.join(RESULTS_DIR,
                           "phase2_monte_carlo_correlation_amox_clav_losses.npz")
    np.savez_compressed(out_npz, **amox_clav_loss_arrays)

    alog("-" * 75)
    alog("EXPORTS")
    alog("-" * 75)
    alog(f"  {out_csv}")
    alog(f"  {out_npz}")
    alog("")
    alog("=" * 75)
    alog("PARTIAL-CORRELATION EXTENSION RUN COMPLETE")
    alog("=" * 75)
    alog("Day 4 byte-exact close on Module 3 primary and sensitivity outputs")
    alog("is preserved (Module 3 itself was not modified). New outputs:")
    alog(f"  phase2_monte_carlo_correlation.csv          ({len(results_df)} rows)")
    alog(f"  phase2_monte_carlo_correlation_audit.txt    (this file)")
    alog(f"  phase2_monte_carlo_correlation_amox_clav_losses.npz")
    alog(f"    ({len(amox_clav_loss_arrays)} AMOX_CLAV loss arrays for Figure 4)")
    alog("")

    audit_path = os.path.join(RESULTS_DIR,
                              "phase2_monte_carlo_correlation_audit.txt")
    with open(audit_path, "w", encoding="utf-8") as f:
        f.write("\n".join(audit_lines))
    print(f"\nAudit log written: {audit_path}")
