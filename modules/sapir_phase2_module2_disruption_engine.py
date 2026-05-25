"""
SAPIR-Net Phase 2 Module 2: Disruption Probability Engine (Antibiotics)
========================================================================
Phase 2 extension of the SAPIR-Net methodology (Zenodo DOI
10.5281/zenodo.19549343; code DOI 10.5281/zenodo.19548610) from
oncology APIs to antibiotic APIs. Methodology unchanged from Phase 1
(scenario typology, distributional choices, distribution shape
parameters); application domain extended.

Three core scenarios per Phase 1 typology:
  A) Direct Ban: Total loss of CN + IN direct exports to US for all
     four antibiotic HS codes.
  B) Cascading Upstream Shock: CN direct loss + parameterized ROW
     degradation via per-HS-code china_exposure bands, with Pareto-
     distributed severity variance. Each HS code's ROW edge is
     degraded independently per Monte Carlo iteration.
  C) Global Logistics Chokepoint: Log-normal capacity reduction
     applied uniformly across ALL L1 -> L2 edges, regardless of
     origin country or HS code.

Phase 2 structural difference from Phase 1: the L2 -> L3 propagation
step. Phase 1's three drugs (Cisplatin, Carboplatin, Methotrexate)
each depended on a single L2 chemical commodity. Phase 2 includes
Augmentin (amoxicillin/clavulanate) as a DUAL CHOKEPOINT drug,
requiring both PENI (6-APA via HS 294110) and OTHAB (clavulanic
acid via HS 294190 proxy) simultaneously. The L2 -> L3 propagation
for AMOX_CLAV is MULTIPLICATIVE JOINT capacity loss, modeling the
compounded risk of two independent upstream chokepoint disruptions:

    capacity_remaining(AMOX_CLAV) = p_PENI * p_OTHAB

where p_HS = remaining(HS) / baseline_total(HS). This is the
standard series-system reliability model under the assumption of
independent ROW degradation per HS code (which is how the Pareto
samples are drawn in Module 3). The CN ban is perfectly correlated
across HS codes (CN edges zero out simultaneously) and the model
captures that correctly through the addition of independently
sampled ROW degradations on each HS. Single-input drugs (AMOX,
AZI, DOXY) propagate L2 -> L3 directly: capacity_remaining(drug)
= p_HS for the drug's single upstream chemical commodity.

Calibration source: row_china_exposure bands at the chemical-class
level, from phase2_calibration_table.csv produced by Module 1.
Anchor sources for the bands per Execution Brief Section 4.3:
NASEM 2022 (Building Resilience into the Nation's Medical Product
Supply Chains, ch. 3, citing Schondelmeyer et al. 2020);
Schondelmeyer USCC testimony June 5, 2025; USP Quality Matters 2024
("amoxicillin synthesis ultimately depends on four KSMs, each
produced almost entirely in China"); Socal et al. 2025 JAMA Health
Forum (DOI 10.1001/jamahealthforum.2025.3871).

Output
------
results/phase2_module2_distribution_audit.txt
    Distribution validation report for Red Team audit; samples from
    each scenario at 50,000 draws per HS code; reports moments and
    tail behavior.

Dependencies: numpy, scipy, pandas

Author: Lead Data Architect / Independent Researcher
Version: Phase 2 Module 2 v1.0, May 2026
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import pareto, lognorm

# Reproducibility per Phase 1 method
RNG_SEED = 42
np.random.seed(RNG_SEED)

RESULTS_DIR = "results"
os.makedirs(RESULTS_DIR, exist_ok=True)

# ============================================================
# SECTION 1: COMMODITY-SPECIFIC ROW CHINA EXPOSURE PARAMETERS
# ============================================================
# Phase 2 bands at the chemical-class level (per Execution Brief
# Section 4.3 and phase2_calibration_table.csv). These represent
# the estimated fraction of Rest-of-World upstream API supply that
# itself depends on Chinese sub-inputs (key starting materials,
# fermentation intermediates) for each antibiotic chemical class.
#
# Bands feed Scenario B's Pareto-rescaled degradation sampler.
# Module 3's Monte Carlo loop samples within the per-HS band each
# iteration; sensitivity analysis runs at lower bound, midpoint,
# and upper bound separately.

ROW_CHINA_EXPOSURE = {
    "294110": {
        "label": "Penicillin (beta-lactam) / 6-APA",
        "low": 0.70,
        "high": 0.90,
        "rationale": (
            "NASEM 2022 (citing Schondelmeyer et al. 2020); "
            "Schondelmeyer USCC testimony June 5, 2025 ('more than 70% "
            "of global production of APIs key ingredients like "
            "penicillin G and amoxicillin'); USP Quality Matters 2024 "
            "('amoxicillin synthesis ultimately depends on four KSMs, "
            "each produced almost entirely in China'); Socal et al. "
            "2025 JAMA Health Forum (China = 70.1% of US antibiotic "
            "API imports by volume, 2024). Strongest source stack of "
            "the four bands."
        ),
    },
    "294130": {
        "label": "Tetracycline",
        "low": 0.65,
        "high": 0.85,
        "rationale": (
            "NASEM 2022; Schondelmeyer USCC June 2025 (China dominant "
            "for tetracycline and doxycycline specifically; India has "
            "discontinued doxycycline intermediate production in favor "
            "of Chinese supply); Socal et al. 2025 aggregate. The "
            "2013-2014 doxycycline-class shortage following Chinese "
            "tetracycline API manufacturer exit is the prototypical "
            "Scenario B historical precedent."
        ),
    },
    "294150": {
        "label": "Erythromycin / macrolide",
        "low": 0.60,
        "high": 0.85,
        "rationale": (
            "NASEM 2022; Schondelmeyer USCC June 2025 (China dominant "
            "for erythromycin and azithromycin); Socal et al. 2025 "
            "aggregate. Band slightly wider than HS 294110 because "
            "public concentration data is less granular for macrolide "
            "fermentation. Phase 2 empirical direct-China share for "
            "this HS code shifted from 15.31% (2018) to 83.10% (2024); "
            "sub-recent-years sensitivity reported in Module 3."
        ),
    },
    "294190": {
        "label": "Clavulanic acid (beta-lactamase inhibitor)",
        "low": 0.55,
        "high": 0.80,
        "rationale": (
            "Industry reports on clavulanic acid fermentation "
            "concentration; WHO essential medicines context. Public "
            "concentration data sparser than for the other three "
            "classes; band widest accordingly. Applied to OTHAB -> "
            "AMOX_CLAV edge (DUAL CHOKEPOINT). Chemical-level band "
            "applied directly; the HS 294190 bundled empirical share "
            "(20.79% USD, 2020-2024) is not used as the parameter "
            "value because the bundle includes cephalosporins, "
            "glycopeptides, lincosamides, and other non-clavulanic-"
            "acid antibiotics."
        ),
    },
}


# ============================================================
# SECTION 2: SCENARIO DEFINITIONS
# ============================================================

class ScenarioA:
    """
    Direct Export Ban.
    CN and IN direct exports to the US drop to zero across all four
    antibiotic HS codes simultaneously. ROW edges are unaffected.
    Deterministic scenario -- no stochastic component.
    """

    name = "A: Direct Ban (CN + IN)"

    @staticmethod
    def apply(edge_weights: dict) -> dict:
        """
        Args:
            edge_weights: dict keyed by (L1_source, HS_code) tuples,
                          values are baseline trade flow amounts in USD.
        Returns:
            Modified edge weights with CN and IN zeroed out.
        """
        modified = edge_weights.copy()
        for key in modified:
            source, _ = key
            if source in ("CN", "IN"):
                modified[key] = 0.0
        return modified


class ScenarioB:
    """
    Cascading Upstream Shock.
    CN direct exports drop to zero across all four HS codes.
    ROW edges are degraded by a stochastic fraction drawn from a
    Pareto distribution, bounded by the HS-code-specific
    china_exposure range. Each HS code's ROW degradation is sampled
    INDEPENDENTLY per Monte Carlo iteration.

    Pareto distribution rationale (unchanged from Phase 1):
      Supply chain disruptions exhibit fat-tailed severity. Most
      disruptions cause moderate degradation, but a non-trivial
      probability mass exists for near-total upstream collapse.

    Pareto parameterization (unchanged from Phase 1):
      Shape (alpha) = 2.5: finite variance, moderately heavy tail.
        P(severity > 2x median) ~ 18%.
      The raw Pareto sample is rescaled to the [low, high] exposure
      band via the inverse-CDF mapping 1 - 1/x.
    """

    name = "B: Cascading Upstream Shock (CN ban + ROW degradation)"

    PARETO_ALPHA = 2.5

    @staticmethod
    def sample_row_degradation(hs_code: str, n_samples: int = 1,
                                band_override: tuple = None) -> np.ndarray:
        """
        Sample ROW degradation fractions for a given commodity.

        Returns values in [low, high] range, distributed according to
        a rescaled Pareto. Bulk of mass near 'low', fat tail extending
        toward 'high'.

        Args:
            hs_code: one of "294110", "294130", "294150", "294190".
            n_samples: number of samples to draw.
            band_override: optional (low, high) tuple to override the
                ROW_CHINA_EXPOSURE band. Used by sensitivity analysis
                to run at lower/midpoint/upper bound, and by the HS
                294150 sub-recent-years sensitivity scenario.

        Returns:
            Array of degradation fractions in [low, high].
        """
        if band_override is not None:
            low, high = band_override
        else:
            params = ROW_CHINA_EXPOSURE[hs_code]
            low, high = params["low"], params["high"]

        # Edge case: degenerate band (low == high) -- e.g., a sensitivity
        # run pinned to a single point. Return the constant.
        if low == high:
            return np.full(n_samples, low)

        raw = pareto.rvs(b=ScenarioB.PARETO_ALPHA, size=n_samples)
        # Map [1, inf) -> [0, 1) via the inverse-CDF transform
        normalized = 1.0 - (1.0 / raw)
        degradation = low + normalized * (high - low)
        # Hard clamp to band for safety (Pareto tail can otherwise spill)
        return np.clip(degradation, low, high)

    @staticmethod
    def apply(edge_weights: dict, degradation_samples: dict) -> dict:
        """
        Args:
            edge_weights: dict keyed by (L1_source, HS_code) tuples.
            degradation_samples: dict keyed by HS_code, values are
                                 single float degradation fractions
                                 (one realization per HS for one MC iter).
        Returns:
            Modified edge weights.
        """
        modified = edge_weights.copy()
        for key in modified:
            source, commodity = key
            if source == "CN":
                modified[key] = 0.0
            elif source == "ROW" and commodity in degradation_samples:
                frac = degradation_samples[commodity]
                modified[key] = modified[key] * (1.0 - frac)
        return modified


class ScenarioC:
    """
    Global Logistics Chokepoint.
    A systemic capacity reduction (e.g., major port closure, pandemic-
    era shipping disruption) affecting ALL L1 -> L2 edges regardless
    of origin country or HS code.

    Log-normal distribution rationale (unchanged from Phase 1):
      Logistics disruptions are right-skewed: most events cause moderate
      delays (10-30% capacity loss), but extreme events (Suez blockage,
      COVID port shutdowns) can reach 50-70%. Log-normal captures this
      asymmetry naturally.

    Parameterization (unchanged from Phase 1):
      mu (underlying normal mean) = ln(0.20)  -> median loss ~ 20%
      sigma (underlying normal std) = 0.50    -> 95th-pct loss ~ 50-55%
    """

    name = "C: Global Logistics Chokepoint"

    MU = float(np.log(0.20))
    SIGMA = 0.50

    @staticmethod
    def sample_capacity_reduction(n_samples: int = 1) -> np.ndarray:
        """
        Sample global capacity reduction fractions.

        Returns values in (0, 1) representing the fraction of total
        supply capacity lost across all L1 -> L2 edges. Median ~ 0.20,
        95th percentile ~ 0.55.

        Returns:
            Array of capacity reduction fractions, clamped to
            [0.01, 0.95].
        """
        raw = lognorm.rvs(s=ScenarioC.SIGMA, scale=np.exp(ScenarioC.MU),
                          size=n_samples)
        return np.clip(raw, 0.01, 0.95)

    @staticmethod
    def apply(edge_weights: dict, capacity_reduction: float) -> dict:
        """
        Args:
            edge_weights: dict keyed by (L1_source, HS_code) tuples.
            capacity_reduction: single float, fraction of capacity lost.
        Returns:
            Modified edge weights (all edges reduced uniformly).
        """
        modified = edge_weights.copy()
        for key in modified:
            modified[key] = modified[key] * (1.0 - capacity_reduction)
        return modified


# ============================================================
# SECTION 3: L2 -> L3 PROPAGATION
# ============================================================
# Phase 2 includes one DUAL CHOKEPOINT drug (AMOX_CLAV) requiring
# both PENI (HS 294110) and OTHAB (HS 294190 proxy) simultaneously.
# The propagation rule is MULTIPLICATIVE JOINT capacity loss,
# representing compounded series-system reliability under
# independent ROW degradation per HS. Single-input drugs (AMOX, AZI,
# DOXY) propagate as the capacity fraction of their sole upstream
# chemical commodity.

DRUG_DEPENDENCIES = {
    "AMOX":      {"chokepoints": ["294110"],
                  "rule": "single",
                  "label": "Amoxicillin"},
    "AMOX_CLAV": {"chokepoints": ["294110", "294190"],
                  "rule": "multiplicative_joint",
                  "label": "Amoxicillin/clavulanate (Augmentin)"},
    "AZI":       {"chokepoints": ["294150"],
                  "rule": "single",
                  "label": "Azithromycin"},
    "DOXY":      {"chokepoints": ["294130"],
                  "rule": "single",
                  "label": "Doxycycline"},
}


def propagate_l2_to_l3(remaining_by_hs: dict, baseline_by_hs: dict) -> dict:
    """
    Compute L3 capacity remaining for each drug given L2 capacity
    fractions.

    Args:
        remaining_by_hs: dict keyed by HS_code, values are post-shock
                         remaining trade values (USD).
        baseline_by_hs: dict keyed by HS_code, values are pre-shock
                        baseline trade values (USD).

    Returns:
        dict keyed by drug ID (AMOX, AMOX_CLAV, AZI, DOXY), values
        are fraction of baseline capacity remaining in [0, 1].
    """
    drug_capacity = {}
    for drug, spec in DRUG_DEPENDENCIES.items():
        if spec["rule"] == "single":
            hs = spec["chokepoints"][0]
            drug_capacity[drug] = (
                remaining_by_hs[hs] / baseline_by_hs[hs]
                if baseline_by_hs[hs] > 0 else 0.0
            )
        elif spec["rule"] == "multiplicative_joint":
            p = 1.0
            for hs in spec["chokepoints"]:
                p *= (
                    remaining_by_hs[hs] / baseline_by_hs[hs]
                    if baseline_by_hs[hs] > 0 else 0.0
                )
            drug_capacity[drug] = p
        else:
            raise ValueError(f"Unknown propagation rule: {spec['rule']}")
    return drug_capacity


# ============================================================
# SECTION 4: DISTRIBUTION VALIDATION & AUDIT OUTPUT
# ============================================================

def audit_distributions(n_samples: int = 50_000, write_to_file: bool = True):
    """
    Generate and summarize sample distributions for Red Team audit.
    Writes a text report alongside printing to stdout.
    """
    lines = []

    def emit(s: str = ""):
        lines.append(s)
        print(s)

    emit("=" * 75)
    emit("SAPIR-Net Phase 2 Module 2: Disruption Probability Engine")
    emit("Distribution validation audit (n = {:,} samples per HS code)".format(n_samples))
    emit("=" * 75)
    emit(f"RNG seed: {RNG_SEED} (Phase 1 reproducibility convention)")
    emit("")

    # --- ROW China Exposure Parameters (per chemical class) ---
    emit("-" * 75)
    emit("CHEMICAL-CLASS ROW CHINA EXPOSURE BANDS")
    emit("-" * 75)
    for hs, params in ROW_CHINA_EXPOSURE.items():
        emit(f"\n  HS {hs} ({params['label']}):")
        emit(f"    Exposure band: [{params['low']:.0%}, {params['high']:.0%}]")
        emit(f"    Rationale: {params['rationale']}")

    # --- Scenario A ---
    emit("")
    emit("=" * 75)
    emit(f"SCENARIO A: {ScenarioA.name}")
    emit("=" * 75)
    emit("  Type: Deterministic")
    emit("  CN edges: 0% remaining (all four HS codes)")
    emit("  IN edges: 0% remaining (all four HS codes)")
    emit("  ROW edges: 100% remaining (unaffected)")

    # --- Scenario B ---
    emit("")
    emit("=" * 75)
    emit(f"SCENARIO B: {ScenarioB.name}")
    emit("=" * 75)
    emit(f"  Distribution: Pareto (alpha = {ScenarioB.PARETO_ALPHA}), rescaled to per-HS exposure band")
    emit(f"  CN edges: 0% remaining (total ban; all four HS codes)")
    emit(f"  IN edges: 100% remaining (unaffected)")

    for hs in ["294110", "294130", "294150", "294190"]:
        samples = ScenarioB.sample_row_degradation(hs, n_samples)
        label = ROW_CHINA_EXPOSURE[hs]["label"]
        emit(f"\n  HS {hs} ({label}) ROW degradation samples:")
        emit(f"    Band:                   [{ROW_CHINA_EXPOSURE[hs]['low']:.0%}, {ROW_CHINA_EXPOSURE[hs]['high']:.0%}]")
        emit(f"    Mean degradation:       {samples.mean():.4f} ({samples.mean():.1%})")
        emit(f"    Median degradation:     {np.median(samples):.4f} ({np.median(samples):.1%})")
        emit(f"    Std deviation:          {samples.std():.4f}")
        emit(f"    5th percentile:         {np.percentile(samples, 5):.4f} ({np.percentile(samples, 5):.1%})")
        emit(f"    25th percentile:        {np.percentile(samples, 25):.4f} ({np.percentile(samples, 25):.1%})")
        emit(f"    75th percentile:        {np.percentile(samples, 75):.4f} ({np.percentile(samples, 75):.1%})")
        emit(f"    95th percentile:        {np.percentile(samples, 95):.4f} ({np.percentile(samples, 95):.1%})")
        emit(f"    Max observed:           {samples.max():.4f} ({samples.max():.1%})")

    # --- Scenario C ---
    emit("")
    emit("=" * 75)
    emit(f"SCENARIO C: {ScenarioC.name}")
    emit("=" * 75)
    emit(f"  Distribution: Log-normal (mu = ln(0.20) = {ScenarioC.MU:.3f}, sigma = {ScenarioC.SIGMA})")
    emit(f"  Applied to: ALL L1 -> L2 edges uniformly (12 edges; 4 HS x 3 partner buckets)")

    samples_c = ScenarioC.sample_capacity_reduction(n_samples)
    emit(f"\n  Capacity reduction samples:")
    emit(f"    Mean reduction:         {samples_c.mean():.4f} ({samples_c.mean():.1%})")
    emit(f"    Median reduction:       {np.median(samples_c):.4f} ({np.median(samples_c):.1%})")
    emit(f"    Std deviation:          {samples_c.std():.4f}")
    emit(f"    5th percentile:         {np.percentile(samples_c, 5):.4f} ({np.percentile(samples_c, 5):.1%})")
    emit(f"    25th percentile:        {np.percentile(samples_c, 25):.4f} ({np.percentile(samples_c, 25):.1%})")
    emit(f"    75th percentile:        {np.percentile(samples_c, 75):.4f} ({np.percentile(samples_c, 75):.1%})")
    emit(f"    95th percentile:        {np.percentile(samples_c, 95):.4f} ({np.percentile(samples_c, 95):.1%})")
    emit(f"    Max observed:           {samples_c.max():.4f} ({samples_c.max():.1%})")

    # --- L2 -> L3 propagation table ---
    emit("")
    emit("=" * 75)
    emit("L2 -> L3 PROPAGATION RULES")
    emit("=" * 75)
    propagation_table = pd.DataFrame([
        {"Drug": DRUG_DEPENDENCIES[d]["label"],
         "L2_Chokepoints": ", ".join(DRUG_DEPENDENCIES[d]["chokepoints"]),
         "Rule": DRUG_DEPENDENCIES[d]["rule"],
         "Formula": (
             f"capacity = p_{DRUG_DEPENDENCIES[d]['chokepoints'][0]}"
             if DRUG_DEPENDENCIES[d]["rule"] == "single"
             else "capacity = " + " * ".join(
                 f"p_{hs}" for hs in DRUG_DEPENDENCIES[d]["chokepoints"]
             )
         )}
        for d in DRUG_DEPENDENCIES
    ])
    emit(propagation_table.to_string(index=False))
    emit("")
    emit("AMOX_CLAV (Augmentin) is the DUAL CHOKEPOINT drug: penicillin")
    emit("(6-APA via PENI/HS 294110) AND clavulanic acid (via OTHAB/HS 294190")
    emit("proxy). Multiplicative joint capacity loss represents compounded")
    emit("series-system reliability under independent ROW degradation samples")
    emit("per HS code. CN ban remains perfectly correlated across HS (CN edges")
    emit("zero out simultaneously); only the ROW Pareto draws are per-HS-")
    emit("independent.")

    # --- Combined scenario summary table ---
    emit("")
    emit("=" * 75)
    emit("SCENARIO COMPARISON MATRIX")
    emit("=" * 75)
    summary = pd.DataFrame({
        "Scenario": [
            "A: Direct Ban (CN+IN)",
            "B: Cascading Upstream Shock",
            "C: Logistics Chokepoint",
        ],
        "CN_Direct":     ["Zeroed", "Zeroed", "Reduced (stochastic)"],
        "IN_Direct":     ["Zeroed", "Unchanged", "Reduced (stochastic)"],
        "ROW_Impact":    [
            "None",
            "Degraded (Pareto, per-HS band)",
            "Reduced (Log-normal, uniform)",
        ],
        "Distribution":  ["Deterministic",
                          "Pareto (alpha=2.5)",
                          "Log-normal (sigma=0.50)"],
        "Tail_Behavior": ["N/A",
                          "Fat tail (extreme cascades possible)",
                          "Right-skewed (rare severe events)"],
    })
    emit(summary.to_string(index=False))

    emit("")
    emit("=" * 75)
    emit("Module 2 complete. Scenario classes and L2->L3 propagation logic")
    emit("ready for Module 3 (Monte Carlo Simulation Loop).")
    emit("=" * 75)

    if write_to_file:
        out_path = os.path.join(RESULTS_DIR,
                                "phase2_module2_distribution_audit.txt")
        with open(out_path, "w", encoding="utf-8") as f:
            f.write("\n".join(lines))
        print(f"\nDistribution audit written: {out_path}")


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    audit_distributions(n_samples=50_000, write_to_file=True)
