# Network Modeling of U.S. Generic Antibiotic Supply Chain Concentration and Disruption

**Author:** Syed Fahim Abbas Kazmi, BTech

**Affiliation:** Independent Researcher

**Corresponding author:** Syed Fahim Abbas Kazmi; email [to be supplied at submission]; ORCID 0009-0000-5075-9638

**Running head:** SAPIR-Net Phase 2: Antibiotics

**Word count, main text:** 3,299 (JHF Original Investigation target: 3,000)

---

## Key Points

**Question.** Are the 2025–2026 federal pharmaceutical-manufacturing interventions structured to attenuate the cascading-disruption risk that U.S. generic antibiotic supply chains carry?

**Findings.** In a stochastic network simulation of four antibiotic Harmonized System codes, the four-drug scope (amoxicillin, amoxicillin/clavulanate, azithromycin, doxycycline) faced a 100% probability of severe shortage under a cascading upstream shock, structurally guaranteed by source-stack-calibrated chemical-class lower bands. Amoxicillin/clavulanate additionally faced 65.4% probability under a global logistics scenario, against 20.9% for single-input drugs.

**Meaning.** Federal interventions concentrated at the finished-dosage-form layer do not reach the upstream concentration layer at which residual structural risk is located.

---

## Abstract

**Importance.** U.S. generic antibiotic supply chains carry documented upstream geographic concentration associated with recurrent shortages of clinical consequence. The 2025–2026 federal pharmaceutical-manufacturing policy cycle has commenced material interventions whose supply-chain-layer effects remain analytically unaddressed.

**Objective.** To extend the SAPIR-Net network model from oncology to four U.S. generic antibiotic active pharmaceutical ingredient (API) supply chains, quantifying geographic concentration and cascading-disruption vulnerability and evaluating the supply-chain-layer alignment of current federal interventions against the layer at which the model identifies residual structural risk.

**Design.** Stochastic network simulation. A three-layer directed graph linked geopolitical source nodes (China, India, Rest of World), four antibiotic chemical commodity nodes, and four essential-medicine nodes. Three disruption scenarios were applied (deterministic China-plus-India direct export ban; cascading upstream shock with stochastic Rest-of-World degradation; global logistics chokepoint), with Monte Carlo simulation, N = 10,000 iterations per scenario, at fixed pseudo-random seed.

**Data Sources.** UN Comtrade bilateral trade data for U.S. imports under Harmonized System codes 294110, 294130, 294150, and 294190, calendar years 2018–2025. Chemical-class upstream-exposure parameters were calibrated against a convergent federal advisory, U.S.-China Economic and Security Review Commission, U.S. Pharmacopeia, and peer-reviewed source stack.<sup>7,8,9,10</sup>

**Exposures.** The three disruption scenarios above, plus a no-shock baseline.

**Main Outcomes and Measures.** Per drug per scenario: mean, median, 5th and 95th-percentile capacity loss; probability of severe shortage (capacity loss exceeding 30%). Four sensitivity analyses (parameter-band endpoints, baseline-window shift, multiplicative-joint independence).

**Results.** The Herfindahl-Hirschman Index of U.S. import value for the 2020–2024 aggregate exceeded 4,200 on each of the four codes, reaching 6,318 on the tetracycline code (HS 294130). The cascading upstream shock produced a 100% probability of severe shortage for every drug, structurally guaranteed at lower band endpoints (deterministic minimum capacity loss 74.3%–93.2%). Under the global logistics scenario, the three single-input drugs each faced 20.9% probability; amoxicillin/clavulanate, the one dual-chokepoint combination drug, faced 65.4%. Partial-correlation sensitivity left the 100% binary outcome unchanged across ρ ∈ {0, 0.25, 0.50, 0.75}.

**Conclusions and Relevance.** Residual upstream structural risk is not, on current evidence, materially attenuated by federal interventions at the finished-dosage-form layer. Closing the residual requires upstream-layer-aligned federal investment beyond present ReShoreRx-scale awards. Combination antibiotic products warrant separate supply chain risk classification due to multiplicative compounding of upstream chokepoints.

---

## Introduction

Drug shortages remain the most persistent failure mode of the U.S. pharmaceutical infrastructure, and the SAPIR-Net research program addresses the absence of an open, quantitative framework for forecasting them.<sup>1</sup> Phase 1 of SAPIR-Net applied a three-layer directed-graph model to three essential generic oncology drugs and identified a structurally determined divergence in vulnerability profiles driven by upstream chemical concentration.<sup>1</sup> The present study extends the same model architecture to four essential generic antibiotics: amoxicillin, amoxicillin/clavulanate, azithromycin, and doxycycline. SAPIR-Net[^1] is the network-modeling component of the broader integrated forecasting program described in the Discussion. Antibiotic shortages in the United States are recurrent and structurally similar across the past decade.<sup>2</sup> The FDA Drug Shortage Database reported approximately 80 active drug shortages as of February 2026,<sup>3</sup> and the U.S. Pharmacopeia 2024–2025 Vulnerable Medicines List identified 100 medicines as vulnerable on the basis of essentiality and supply chain risk indicators.<sup>4</sup> The four drugs in the present scope are listed on the World Health Organization's 23rd Model List of Essential Medicines<sup>5</sup>; three are AWaRe "Access" antibiotics and one (azithromycin) is "Watch." Each has been the subject of recurrent shortage listings.<sup>3,6</sup>

[^1]: SAPIR-Net is the network model presented in this research program (continuing from Phase 1, Zenodo DOI 10.5281/zenodo.19549343). It is distinct from the federal Strategic Active Pharmaceutical Ingredients Reserve ("SAPIR") referenced in the August 19, 2025 Presidential Executive Order on essential medicines.<sup>17</sup>

For the antibiotic class, upstream geographic concentration is documented as more severe than for any other major therapeutic category. The National Academies reported in 2022 that "China produces almost all of the APIs used to manufacture such drugs as penicillin G."<sup>7</sup> Schondelmeyer testified to the U.S.-China Economic and Security Review Commission on June 5, 2025 that approximately 70%–80% of the world's antibiotic active pharmaceutical ingredients (APIs) are produced in China, with more than 70% of global production of key API ingredients including penicillin G and amoxicillin originating in China.<sup>8</sup> The U.S. Pharmacopeia found that 679 APIs used in U.S.-approved generic drugs depend on China as the sole supplier of at least one key starting material (KSM); amoxicillin specifically depends on four KSMs each produced almost entirely in China.<sup>9</sup> The most directly comparable peer-reviewed analysis (*JAMA Health Forum*, October 2025) found that China accounted for 70.1% of U.S. antibiotic API imports by volume in 2024, while the finished-dosage-form (FDF) layer is comparatively diversified.<sup>10</sup> Convergent peer-reviewed analyses<sup>11,12</sup> document the same pattern.

Two empirical events anchor the present analysis. The 2013–2014 U.S. doxycycline shortage<sup>13</sup> was triggered by the exit of a major Chinese tetracycline-API manufacturer<sup>14</sup>; the price increase for some formulations reached approximately 5,500%,<sup>15</sup> and Heritage Pharmaceuticals re-entered the market only by re-sourcing API through a European third-party supplier.<sup>14</sup> The 2022–2024 amoxicillin oral suspension shortage was resolved at the dispensing layer by May 2025<sup>3</sup>; the structural upstream conditions that produced it have not been.

The federal response has converged on three actions of direct relevance. The 100-day reviews under Executive Order 14017 produced the most authoritative federal acknowledgment of pharmaceutical supply chain concentration in the public record.<sup>16</sup> The August 19, 2025 Presidential Executive Order on essential medicines established a six-month strategic API supply target.<sup>17</sup> The FDA Commissioner's National Priority Voucher (CNPV) pilot granted its first approval on December 9, 2025: USAntibiotics, Bristol, Tennessee, for Augmentin XR.<sup>18</sup> USAntibiotics has stated that the Bristol facility, at full utilization, can meet 100% of U.S. demand for amoxicillin at the FDF layer; the capacity claim characterizes the FDF layer only, and there is no public evidence that the upstream API and KSM supply on which the facility depends operates independently of the predominantly Chinese global supply base.<sup>4,8,10</sup>

The question this analysis addresses is whether current federal interventions materially attenuate the cascading-disruption risk that the antibiotic supply chain carries upstream. The objectives are to quantify per-code concentration, to model three structural disruption scenarios, to evaluate the federal policy landscape against the supply-chain layer at which residual risk is identified, and to position the analysis relative to two adjacent U.S. drug-supply mapping efforts as methodologically complementary.

---

## Methods

### Network Architecture

The SAPIR-Net model is a three-layer directed graph implemented in Python using NetworkX. Layer 1 nodes represent geopolitical sources (China, India, Rest of World); Layer 2 nodes represent antibiotic chemical commodity classes corresponding to four Harmonized System (HS) codes; Layer 3 nodes represent the four essential medicines. The implementation is organized as four independent, auditable modules (graph construction and concentration analysis; disruption probability engine; Monte Carlo simulation; visualization and reporting). The architecture is identical to SAPIR-Net Phase 1,<sup>1</sup> with the single Phase 2–specific methodological extension of the multiplicative joint capacity-loss rule for combination drugs described below; the application domain extends to antibiotics. The pharmaceutical supply chain itself spans five conventionally numbered layers (0: KSMs and fermentation intermediates; 1: geopolitical sources; 2: chemical commodity APIs; 3: finished APIs; 4: finished dosage forms); the model directly observes Layers 1–2 through UN Comtrade bilateral trade data, computes Layer 3 capacity loss as output, and parameterizes Layer 0 risk through a chemical-class exposure variable applied to the Rest-of-World edge of Layer 1. Layer 4 is documented qualitatively. A full layer-to-model-coverage mapping is provided in the eSupplement.

### Data Sources and Concentration Metric

Bilateral trade data were extracted from the UN Comtrade database<sup>19</sup> for U.S. imports (reporter code 842) under four 6-digit HS codes within HS heading 2941: 294110 (penicillins; amoxicillin and the penicillin component of amoxicillin/clavulanate), 294130 (tetracyclines; doxycycline), 294150 (erythromycin; azithromycin), and 294190 ("other antibiotics, n.e.s."; bundle proxy for the clavulanic acid component of amoxicillin/clavulanate). Partner groups were China, India, and World, with Rest-of-World derived by subtraction. Calendar years 2018 through 2025 were extracted; the 2020–2024 five-year aggregate window was the primary baseline. A separate Jordan-as-partner extraction returned zero rows across all four HS codes for all eight years, corroborating that the HS 2941xx codes filter to the API layer and exclude FDF-only suppliers (full discussion in the eSupplement).<sup>4,10</sup>

Geographic concentration of each commodity's import flow was computed as the Herfindahl-Hirschman Index, HHI = Σ(s_i)², where s_i is the U.S. import share of source country i. The Department of Justice and Federal Trade Commission firm-level HHI thresholds (above 2,500 highly concentrated) are reported as a comparative reference for the trade-data application, consistent with prior peer-reviewed and federal use.<sup>7,10,16</sup>

### Edge Weighting and the Multiplicative Joint Rule

Layer 1 → Layer 2 edges are weighted by 2020–2024 five-year USD-aggregate trade values. Layer 2 → Layer 3 edges carry binary obligate-precursor flags. For single-input drugs, capacity_remaining(drug) = p_HS. For amoxicillin/clavulanate, which depends on two independent fermentation chokepoints simultaneously (6-aminopenicillanic acid from HS 294110 and clavulanic acid via HS 294190), capacity_remaining(AMOX_CLAV) = p_PENI × p_OTHAB, generalizing to Π p_i for any combination drug with N independent upstream chokepoints. This multiplicative joint rule is the single Phase 2–specific methodological extension over Phase 1, applying the standard series-system reliability model under independent disruption sampling; the independence assumption is conservative in direction (any shared Chinese upstream infrastructure would increase the compounded effect), and is quantified by the partial-correlation sensitivity in Results.

Indirect Chinese dependency within Rest-of-World supply is parameterized at the chemical-class level using a variable termed `row_china_exposure`, defined as the estimated fraction of Rest-of-World trade flow that depends on Chinese upstream inputs at Layer 0. Primary calibration bands are 0.70–0.90 (penicillin), 0.65–0.85 (tetracycline), 0.60–0.85 (erythromycin/macrolide), and 0.55–0.80 (clavulanic acid), anchored in the convergent NASEM, Schondelmeyer USCC, USP KSM, and Socal et al. source stack;<sup>7,8,9,10</sup> the full calibration source stack is in the eSupplement.

### Disruption Scenarios and Simulation

Three disruption scenarios were applied plus a no-shock baseline. **Scenario A** (Direct Export Ban) is deterministic: direct exports from China and India across all four HS codes are reduced to zero; Rest-of-World flows are unaffected. **Scenario B** (Cascading Upstream Shock) eliminates China's direct exports and degrades Rest-of-World supply per HS code by a stochastic fraction drawn from a rescaled Pareto distribution (shape α = 2.5) bounded by the chemical-class `row_china_exposure` band, with per-HS draws independent across HS codes in the primary specification. **Scenario C** (Global Logistics Chokepoint) applies a uniform log-normal capacity reduction (μ = ln 0.20, σ = 0.50; median 20%, 95th percentile 50%–55%) to all Layer 1 → Layer 2 edges simultaneously. Empirical anchors for Scenarios B and C are the 2013–2014 doxycycline cascade<sup>13,14,15</sup> and the 2021 Suez Canal blockage plus pandemic-era port shutdowns; full distribution parameterization and validation are in the eSupplement.

Each scenario was executed for N = 10,000 independent iterations at a fixed pseudo-random generator seed (42), regenerating byte-exact identical results on rerun. A severe shortage event was defined as per-iteration capacity loss exceeding 30%. Output metrics per drug per scenario include mean, median, 5th and 95th-percentile capacity loss, and the probability of severe shortage. Reporting follows the non-cost-aspects subset of the CHEERS Reporting Guidelines for decision analytical models.

Four sensitivity analyses were executed: **S1** pinned each chemical-class parameter to its band lower endpoint; **S2** to its band upper endpoint; **S3** replaced the 2020–2024 baseline with a 2023–2025 sub-recent baseline for HS 294150 only (motivated by substantial intra-window drift in the empirical direct China share for that code); and **S4** applied a Gaussian-copula correlation structure to the per-iteration Pareto draws on PENI and OTHAB in Scenario B at correlation values ρ ∈ {0, 0.25, 0.50, 0.75}, quantifying the multiplicative-joint independence assumption.

### Reproducibility and Reporting

All code modules and intermediate outputs are archived at Zenodo (DOI [reserved at deposit, to be inserted in the published version]) and at GitHub ([repository URL to be inserted in the published version]). Full module-by-module methodological detail, the Pareto-distribution validation audit, the complete sensitivity-table panel (49 rows), and the partial-correlation Gaussian-copula derivation are provided in the eSupplement.

---

## Results

### Concentration

Each of the four antibiotic HS codes exhibited highly concentrated U.S. import flows over the 2020–2024 baseline window. Per-code HHI ranged from 4,276 (HS 294150) to 6,318 (HS 294130); aggregate HHI across the four codes was 5,289 (Table 1). The direct China share varied substantially across codes, from 20.8% on HS 294190 to 75.8% on HS 294130 (USD basis). The HS 294190 figure characterizes a heterogeneous "other antibiotics, n.e.s." bundle that includes cephalosporins, glycopeptides, lincosamides, and other antibiotic APIs alongside the clavulanic acid trade that the OTHAB node is intended to capture; it is not an estimate of clavulanic-acid-specific concentration, and the model addresses this by applying the chemical-class band [0.55, 0.80] directly to the OTHAB → AMOX_CLAV edge.

The Phase 2 four-code kg-basis direct China share for 2020–2024 (60.0%; 95.2% kg-valid coverage) sits within 2.6 percentage points of the Socal et al. five-year reference (62.6%) and within 6.1 percentage points of their 2024 single-year reference (70.1%).<sup>10</sup> The two analyses use different trade-data sources and different commodity-code classification systems; the convergence supports the methodological soundness of the present extraction. HS 294150 exhibited a substantial intra-window drift (direct China share 15.3% in 2018 → 50.4% in 2020 → 83.1% in 2024, with a partial correction to 60.7% in 2025); sensitivity analysis S3 addresses this directly.

**Table 1.** *U.S. import concentration for antibiotic HS 2941xx codes, 2020–2024 five-year aggregate.*

| HS Code | Description | China Share (USD) | HHI | DOJ/FTC Reference |
|---|---|---|---|---|
| 294110 | Penicillins | 44.8% | 4,900 | Highly concentrated |
| 294130 | Tetracyclines | 75.8% | 6,318 | Highly concentrated |
| 294150 | Erythromycin | 56.8% | 4,276 | Highly concentrated |
| 294190 | Other antibiotics, n.e.s. (bundle) | 20.8% | 5,932 | Highly concentrated |
| Aggregate (4 codes) | — | 28.8% | 5,289 | Highly concentrated |

### Primary Monte Carlo Outputs

Under the no-shock baseline, all four drugs returned zero capacity loss across all iterations (sanity check). Under Scenario A, the deterministic CN+IN direct export ban produced a 100% probability of severe shortage for every drug, with capacity losses reflecting underlying direct-share concentration: amoxicillin 46.2%, amoxicillin/clavulanate 60.2%, azithromycin 70.9%, doxycycline 75.9%. Under Scenario B, the cascading upstream shock produced a 100% probability of severe shortage for every drug, with mean capacity losses 85.5% (amoxicillin), 95.2% (amoxicillin/clavulanate), 76.4% (azithromycin), and 92.8% (doxycycline); the entire capacity-loss distribution sat above the 30% severe-shortage threshold for every drug. Under Scenario C, the three single-input drugs each faced a 20.9% probability of severe shortage (mean capacity loss 22.6%), while amoxicillin/clavulanate — the one combination drug in scope — faced 65.4% (mean capacity loss 38.6%). The probability matrix is visualized in Figure 1; the Scenario B and C capacity-loss density distributions are visualized in Figure 2.

**Table 2.** *Phase 2 primary Monte Carlo results. N = 10,000 iterations per scenario; RNG seed = 42; chemical-class `row_china_exposure` parameters at primary band values; 2020–2024 baseline.*

| Scenario | Drug | Mean Loss | Median | P5 / P95 | P(Severe) |
|---|---|---|---|---|---|
| Baseline | All four | 0.0% | 0.0% | 0.0% / 0.0% | 0.0% |
| A: Direct Ban | Amoxicillin | 46.2% | 46.2% | 46.2% / 46.2% | 100.0% |
| A: Direct Ban | Amoxicillin/clavulanate | 60.2% | 60.2% | 60.2% / 60.2% | 100.0% |
| A: Direct Ban | Azithromycin | 70.9% | 70.9% | 70.9% / 70.9% | 100.0% |
| A: Direct Ban | Doxycycline | 75.9% | 75.9% | 75.9% / 75.9% | 100.0% |
| B: Cascading Upstream | Amoxicillin | 85.5% | 85.0% | 82.7% / 89.9% | 100.0% |
| B: Cascading Upstream | Amoxicillin/clavulanate | 95.2% | 95.1% | 93.7% / 96.8% | 100.0% |
| B: Cascading Upstream | Azithromycin | 76.4% | 76.1% | 74.5% / 79.4% | 100.0% |
| B: Cascading Upstream | Doxycycline | 92.8% | 92.6% | 91.5% / 94.8% | 100.0% |
| C: Logistics | Amoxicillin | 22.6% | 20.0% | 8.8% / 44.9% | 20.9% |
| C: Logistics | Amoxicillin/clavulanate | 38.6% | 35.9% | 16.8% / 69.7% | **65.4%** |
| C: Logistics | Azithromycin | 22.6% | 20.0% | 8.8% / 44.9% | 20.9% |
| C: Logistics | Doxycycline | 22.6% | 20.0% | 8.8% / 44.9% | 20.9% |

![Figure 1](phase2_fig1_vulnerability_heatmap.png)

*Figure 1.* Probability of severe drug shortage (>30% capacity loss) by disruption scenario and drug. Monte Carlo simulation, N = 10,000 iterations per scenario. The amoxicillin/clavulanate Scenario C cell (65.4%) is the dual-chokepoint signature: more than three times the single-input drugs' probability (20.9%) under the same logistics shock.

![Figure 2](phase2_fig2_capacity_loss_densities.png)

*Figure 2.* Capacity-loss density distributions under Scenarios B (cascading upstream shock) and C (logistics chokepoint), with Scenario A deterministic reference lines and the 30% severe-shortage threshold marked. Under both scenarios, the amoxicillin/clavulanate distribution sits to the right of the single-input amoxicillin distribution, visualizing the dual-chokepoint compounding.

### Sensitivity Analyses

The S1 result is the structural-guarantee disclosure. With each chemical-class parameter pinned to its band lower endpoint and the Pareto sampler returning the constant low value, mean capacity losses under Scenario B were 82.4% (amoxicillin), 93.2% (amoxicillin/clavulanate), 74.3% (azithromycin), and 91.4% (doxycycline) — all above the 30% severe-shortage threshold by 44 to 63 percentage points. The S1 collapse is not specific to the S1 sensitivity: it is the deterministic band-arithmetic floor of every iteration of Scenario B in the primary run, because per-iteration Pareto draws sit at or above the lower band endpoint by construction. The 100% probability outcome in the primary results is the binary indicator equivalent of this structural-guarantee floor; the primary Scenario B varies in the magnitude of loss across iterations but never in the binary outcome. Closing the Scenario B vulnerability requires shifting the chemical-class parameter materially below the bottom of its source-stack-calibrated band — to a level lower than the convergent NASEM, USCC, USP, and peer-reviewed evidence currently supports.<sup>7,8,9,10</sup> Federal interventions located downstream of that parameter do not move it in the relevant direction.

S2 (upper band) confirmed robustness across the entire calibrated band; mean Scenario B capacity losses reached 93.2%–98.6% across the four drugs. S3 (HS 294150 sub-recent baseline) established an asymmetric structural-robustness finding: for azithromycin, Scenario A mean capacity loss jumped 19.83 percentage points under the sub-recent 2023–2025 baseline (because the deterministic direct-ban scenario is tied to the empirical direct China share, which shifted upward substantially), while Scenario B moved only 2.93 percentage points and Scenario C was unchanged (band-independent). A direct-export-ban analysis evaluated against the 2020–2024 aggregate would understate the risk that the more recent 2023–2025 trade pattern actually carries; the cascading-upstream-shock result is robust to the same shift (Figure 3). S4 (partial-correlation sensitivity on the multiplicative-joint independence assumption) left the binary 100% probability outcome unchanged across ρ ∈ {0, 0.25, 0.50, 0.75}; the magnitude distribution shifted in the predicted conservative direction (AMOX_CLAV Scenario B 95th-percentile capacity loss +0.45 percentage points at ρ = 0.75 vs ρ = 0). The independence assumption is conservative in direction. Full per-ρ summary statistics, the corresponding density-distribution visualization, and a one-paragraph Pareto α-sensitivity treatment are reported in the eSupplement.

![Figure 3](phase2_fig3_federal_policy_sensitivity.png)

*Figure 3.* Federal-policy-aware sensitivity analysis. Left: Scenario B band sensitivity (S1 lower endpoint / primary band sampling / S2 upper endpoint) for each of the four drugs, with the P(severe shortage) = 100% annotation across all configurations. Right: azithromycin baseline-window sensitivity comparing the primary 2020–2024 baseline against the S3 2023–2025 sub-recent baseline across Scenarios A, B, and C.

---

## Discussion

The principal contribution of the present analysis is the per-chemical-class disaggregation of an aggregated finding established at the class-aggregate level in the peer-reviewed antibiotic-API trade literature.<sup>10</sup> The disaggregation contributes a per-code HHI panel identifying HS 294130 (tetracycline) at the most highly concentrated position within the four-code scope, and the dual-chokepoint signature discussed below — neither of which the class-aggregate analysis generates by construction.

The central policy implication is structural rather than incremental. Each of the three most consequential federal pharmaceutical-supply interventions of the 2025–2026 cycle — the CNPV pilot, the August 19, 2025 Executive Order, and the ASPR ReShoreRx investment portfolio — addresses a documented vulnerability at a layer at which federal investment is comparatively tractable. The model identifies the residual risk at a different layer. The S1 result quantifies the implication: at every chemical-class band lower endpoint within the source-stack-calibrated envelope, the probability of severe shortage under the cascading upstream shock remains 100% for every drug. Closing that residual requires shifting the chemical-class exposure parameter below the band — a shift not produced by FDF-layer capacity additions, by a six-month API reserve drawn from existing predominantly Chinese-sourced supply, or by the present cumulative scale of upstream-layer-aligned ReShoreRx awards.<sup>21,22</sup> The CNPV first approval reduces FDF-layer concentration substantively and is a real step in U.S. pharmaceutical manufacturing policy on its own terms; the analysis is not a critique of that step, but a quantification of what remains. As of May 2026, eighteen CNPV vouchers have been issued and seven approvals granted; none addresses the upstream API or KSM manufacturing layer for an antibiotic or any other pharmaceutical class.<sup>23</sup> A public layer-classification at the time of each subsequent approval would allow the present framework to evaluate the Scenario B probability attenuation produced by each approval directly. Parallel ongoing policy analyses including the Council on Geostrategy April 2026 report<sup>24</sup> reach convergent qualitative conclusions; the present analysis contributes the quantitative network-simulation framework that adjacent policy work does not.

The dual-chokepoint signature is the most novel finding relative to Phase 1. Amoxicillin/clavulanate's 65.4% probability of severe shortage under the global logistics scenario, against 20.9% for the three single-input drugs in the same scope, identifies combination antibiotic products as a structurally compounded risk category hidden in single-API analyses. The mechanism is the multiplicative joint capacity-loss rule applied to two independent fermentation chokepoints (6-aminopenicillanic acid via penicillin G fermentation and clavulanic acid via *Streptomyces clavuligerus* fermentation). The implication for the August 19, 2025 Executive Order's six-month strategic API reserve is concrete: a reserve organized at the API-class level (e.g., six months of amoxicillin API) would by construction fail to protect amoxicillin/clavulanate through any event affecting the clavulanic acid fermentation base. The chemistry of the class compounds this further: the hygroscopic and thermally labile character of 6-aminopenicillanic acid and clavulanic acid intermediates introduces shelf-life and first-in-first-out rotation requirements that tie any reserve's operational viability to ongoing upstream production capacity rather than to one-time stockpile builds. Reserve design for combination antibiotic products needs to identify and stock each independent upstream chemistry separately. The multiplicative rule generalizes to any combination drug with N independent upstream chokepoints; future SAPIR-Net phases on anesthesia drugs (where fixed-dose combinations are common) will likely encounter the same mechanism.

The present analysis is methodologically distinct from two adjacent U.S. drug-supply mapping efforts. The University of Minnesota Resilient Drug Supply Project<sup>8</sup> is a facility-level production map assembled from FDA Drug Master Files and proprietary industry data. The U.S. Pharmacopeia Medicine Supply Map and Vulnerable Medicines List<sup>4</sup> produce an integrated essentiality-and-shortage-risk composite score; a separate USP analysis reports that 67% of antimicrobial-API Drug Master Files map to facilities located in India and China.<sup>20</sup> SAPIR-Net Phase 2 is a public-data Monte Carlo simulation over a directed-graph supply chain model with stochastic shock propagation; the approaches converge on the same structural finding through independent data layers, reinforcing the policy implication and positioning the present framework as methodologically complementary rather than competitive.

### Limitations

Several limitations bear on interpretation. The HS 294190 bundle includes non-clavulanic-acid antibiotic APIs; the model addresses this by applying the chemical-class band directly to the OTHAB → AMOX_CLAV edge. The `row_china_exposure` variable is an informed estimate rather than a directly measured trade flow; the bands are anchored in the convergent NASEM, USCC, USP, and Socal et al. evidence stack,<sup>7,8,9,10</sup> and the structural-guarantee result holds at the lower band endpoint. The binary stoichiometric Layer 2 → Layer 3 model assumes zero inventory buffer, zero substitution elasticity, and instantaneous propagation, representing a worst-case acute-phase scenario; real-world mitigation attenuates observed impact relative to model output. The dual-chokepoint multiplicative rule's independence assumption is conservative in direction per the S4 partial-correlation sensitivity. The model represents acute-phase vulnerability and does not account for month-to-year market and policy adjustment. The four-drug scope does not cover the entire antibiotic class; broader penicillin-class, cephalosporin-class, fluoroquinolone-class, aminoglycoside-class, and carbapenem-class extensions are held for subsequent work. The federal pharmaceutical-supply policy environment continues to evolve; subsequent CNPV approvals and program developments may refine the policy framing<sup>25</sup> without altering the model's structural findings, which are calibrated against trade data and source-stack evidence and are not contingent on the specific federal policy state of any single date.

### Conclusion

Stochastic network simulation calibrated to UN Comtrade bilateral trade data and to the convergent NASEM, USCC, USP, and peer-reviewed source stack identifies a 100% probability of severe shortage for four essential U.S. generic antibiotics under a cascading upstream shock, structurally guaranteed at the source-stack-calibrated chemical-class lower band endpoints. The federal pharmaceutical-manufacturing interventions of the 2025–2026 cycle operate predominantly at layers downstream of the layer at which the model identifies the residual structural risk. The dual-chokepoint signature in amoxicillin/clavulanate (65.4% Scenario C probability against 20.9% for single-input drugs) identifies combination antibiotic products as a structurally compounded risk category. Phase 2 is the antibiotics extension of the integrated SAPIR-Net program, which also includes the National Pharmaceutical Asset Optimization and Resilience Framework addressing facility-level production reliability through hardware-agnostic predictive maintenance,<sup>26</sup> and the forthcoming National Pharmaceutical Intelligence Network integrating macro-geopolitical exposure with facility-level reliability. The integrated program is the endeavor; Phase 2 is the antibiotics extension of it.

---

## Data and Code Availability

All Phase 2 code modules, the UN Comtrade extraction, intermediate analysis outputs, simulation outputs, the partial-correlation sensitivity runner script, the figure-generation scripts, and the full audit logs are archived at Zenodo (DOI [reserved at deposit, to be inserted at publication]) and at GitHub ([repository URL to be inserted at publication]). The fixed pseudo-random generator seed (42) reproduces all reported numerical outputs to byte-exact tolerance.

---

## Acknowledgments

The author thanks the staff of the University of Minnesota PRIME Institute and the Resilient Drug Supply Project, the U.S. Pharmacopeia Medicine Supply Map program, and the Johns Hopkins Bloomberg School of Public Health antibiotic supply chain research group for the public-record peer-reviewed work and federally recorded testimony that anchor the calibration of the model. The author thanks the United Nations Statistics Division for the open Comtrade data infrastructure, and the FDA, CDC, ASPR, and Federal Register publication systems for the open record of federal pharmaceutical-supply policy. None of those organizations bears responsibility for the analysis or for any of its conclusions.

**Author Contributions.** Sole author; all conception, design, data extraction, model implementation, analysis, interpretation, and drafting.

**Conflicts of Interest.** None declared.

**Funding.** Independent research; no funding to disclose.

**Use of AI in Publication and Research.** The author used Anthropic's Claude (claude-opus-4 family, Anthropic, San Francisco, CA; accessed via the Anthropic API and Claude.ai interface) as a drafting, code-review, and audit-protocol assistant during manuscript preparation. AI assistance was used for prose drafting and editing of expository sections, for Python code generation and review for the four analytical modules, and as the drafting partner for an internal Blue-Team / Red-Team / subject-matter peer-review audit protocol. All methodology decisions, data extraction parameters, chemical-class band calibrations, scenario specifications, and policy framing were made by the human author. All citations were verified against primary sources. The author is solely responsible for the accuracy, integrity, and conclusions of the manuscript.

---

## References

1. Kazmi SFA. SAPIR-Net: Supply Chain Analysis for Pharmaceutical Infrastructure Resilience — A Network Model for Quantifying Cascading Vulnerability in the U.S. Generic Oncology Drug Supply Chain. Phase 1: Cisplatin, Carboplatin, and Methotrexate. Zenodo. 2026. doi:10.5281/zenodo.19549343

2. Malta M, Christian M, Cadwallader AB. Oncology drug shortages: impacts, policy reforms, and advocacy imperatives. *Cancer J*. 2025;31:e0784. doi:10.1097/PPO.0000000000000784

3. U.S. Food and Drug Administration. CDER Drug Shortage Database. Accessed May 2026. https://www.accessdata.fda.gov/scripts/drugshortages/

4. United States Pharmacopeia. 2024–2025 Vulnerable Medicines List for the United States: a data-based approach to identify risks and enable interventions to increase reliability of supply. USP Issue Brief. 2024. https://www.usp.org/sites/default/files/usp/document/our-impact/2024-2025-vulnerable-medicines-list.pdf

5. World Health Organization. WHO Model List of Essential Medicines — 23rd list, 2023. WHO/MHP/HPS/EML/2023.02. Geneva: World Health Organization; 2023. https://www.who.int/publications/i/item/WHO-MHP-HPS-EML-2023.02

6. American Society of Health-System Pharmacists. ASHP Drug Shortages Resource Center. Accessed May 2026. https://www.ashp.org/drug-shortages

7. National Academies of Sciences, Engineering, and Medicine. Building Resilience Into the Nation's Medical Product Supply Chains. National Academies Press; 2022. doi:10.17226/26420

8. Schondelmeyer SW. Designing a resilient U.S. drug supply. Statement before the U.S.-China Economic and Security Review Commission, hearing on "Dominance by Design: China Shock 2.0 and the Supply Chain Chokepoints Eroding U.S. Security," Panel II. June 5, 2025. https://www.uscc.gov/sites/default/files/2025-06/Stephen_Schondelmeyer_Testimony.pdf

9. United States Pharmacopeia. Concentrated origins, widespread risk: new USP insights on key starting materials. USP Quality Matters Blog. 2024. https://qualitymatters.usp.org/concentrated-origins-widespread-risk-new-usp-insights-key-starting-materials

10. Socal MP, Sun Y, Ballreich JM, Lambert JD, Dai T, Dada M. US antibiotic importation and supply chain vulnerabilities. *JAMA Health Forum*. 2025;6(10):e253871. doi:10.1001/jamahealthforum.2025.3871

11. Socal MP, Ahn K, Greene JA, Anderson GF. Competition and vulnerabilities in the global supply chain for US generic active pharmaceutical ingredients. *Health Aff (Millwood)*. 2023;42(3):407-415. doi:10.1377/hlthaff.2022.01120

12. Yang Y, Husain L, Huang Y. China's position and competitiveness in the global antibiotic value chain: implications for global health. *Global Health*. 2024;20(1):87. doi:10.1186/s12992-024-01089-x

13. Centers for Disease Control and Prevention. Nationwide shortage of doxycycline: resources for providers and recommendations for patient care. CDC Health Alert Network advisory CDCHAN-00349. June 12, 2013. https://stacks.cdc.gov/view/cdc/25258

14. DrBicuspid.com. Doxycycline, tetracycline shortage eases. DrBicuspid.com [industry trade press]. November 21, 2013. https://www.drbicuspid.com/dental-hygiene/infection-control/antibiotics/article/15368317/doxycycline-tetracycline-shortage-eases

15. Alpern JD, Zhang L, Stauffer WM, Kesselheim AS. Trends in pricing and generic competition within the oral antibiotic drug market in the United States. *Clin Infect Dis*. 2017;65(11):1848-1852. doi:10.1093/cid/cix634

16. The White House. Building Resilient Supply Chains, Revitalizing American Manufacturing, and Fostering Broad-Based Growth: 100-Day Reviews Under Executive Order 14017. June 2021. https://www.whitehouse.gov/wp-content/uploads/2021/06/100-day-supply-chain-review-report.pdf

17. The White House. Executive Order on Essential Medicines and the Strategic Active Pharmaceutical Ingredients Reserve. *Federal Register*. August 19, 2025;90(158):Doc 2025-15823.

18. U.S. Food and Drug Administration. FDA approves first Commissioner's National Priority Voucher application: Augmentin XR for USAntibiotics. FDA Press Announcement. December 9, 2025. https://www.fda.gov/news-events/press-announcements

19. United Nations. UN Comtrade Database. United Nations Department of Economic and Social Affairs, Statistics Division. Accessed May 2026. https://comtradeplus.un.org/

20. United States Pharmacopeia. USP Medicine Supply Map: antimicrobial Drug Master File analysis. Accessed May 2026. https://www.usp.org/supply-chain/medicine-supply-map

21. U.S. Department of Health and Human Services, Administration for Strategic Preparedness and Response. HHS announces $8.3 million investment to expand U.S.-based manufacturing of essential medicines through the API Innovation Center. ASPR Newsroom. March 25, 2026. https://aspr.hhs.gov/newsroom/Pages/APIIC-propofol-metoprolol-March-2026.aspx

22. U.S. Department of Health and Human Services, Administration for Strategic Preparedness and Response. HHS awards $14 million to the API Innovation Center to onshore manufacturing for asthma, hypertension, and anxiety APIs. ASPR Newsroom. September 2024. https://aspr.hhs.gov/newsroom/Pages/APIIC-September-2024.aspx

23. U.S. Food and Drug Administration. Commissioner's National Priority Voucher pilot program; selection announcements. FDA Press Announcements. October 22, 2025 – May 2026. https://www.fda.gov/news-events/press-announcements

24. Council on Geostrategy. Foreign control of antibiotic supply: the systemic vulnerabilities the transatlantic supply chain carries. Council on Geostrategy Policy Paper. April 2026. https://www.geostrategy.org.uk/research/foreign-control-antibiotic-supply

25. U.S. Food and Drug Administration. Commissioner's National Priority Voucher pilot program; public hearing. Notice of public hearing scheduled for June 4, 2026. *Federal Register*. March 20, 2026;91(56):Doc 2026-05573.

26. Kazmi SFA. NPAORF: National Pharmaceutical Asset Optimization & Resilience Framework — a hardware-agnostic, GxP-compliant predictive maintenance platform for mitigating the 2025–2026 U.S. drug shortage crisis. Zenodo. 2026. doi:10.5281/zenodo.19310969
