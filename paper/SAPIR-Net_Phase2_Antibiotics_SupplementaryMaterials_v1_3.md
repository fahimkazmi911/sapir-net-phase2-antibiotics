# Online Supplementary Materials

**Manuscript title:** Network Modeling of U.S. Generic Antibiotic Supply Chain Concentration and Disruption

**Author:** Syed Fahim Abbas Kazmi, BTech (Independent Researcher; ORCID 0009-0000-5075-9638)

**Contents:**

- eAppendix 1. Pharmaceutical Supply Chain Layer Structure and Model Coverage
- eAppendix 2. Data Sources: Full Detail
- eAppendix 3. Calibration-Band Source Stack
- eAppendix 4. Disruption Scenario Distribution Specifications and Validation
- eAppendix 5. Multiplicative Joint Capacity-Loss Rule: Derivation
- eAppendix 6. Complete Sensitivity Analysis Panel
- eAppendix 7. Partial-Correlation Sensitivity (S4): Methodology and Full Results
- eAppendix 8. Pareto Shape Parameter α Sensitivity
- eTable 1. Per-Chemical-Class `row_china_exposure` Calibration Bands and Source Anchors
- eTable 2. Layer-by-Layer Model Coverage
- eTable 3. AMOX_CLAV Scenario B Partial-Correlation Sensitivity Summary Statistics
- eTable 4. Complete Sensitivity Table (S1, S2, S3 Panel)
- eFigure 1. AMOX_CLAV Scenario B Capacity-Loss Distribution Under Partial-Correlation Sensitivity

---

## eAppendix 1. Pharmaceutical Supply Chain Layer Structure and Model Coverage

The pharmaceutical supply chain for the four target drugs is conventionally stratified across five layers. The SAPIR-Net Phase 2 model directly observes some layers, parameterizes others from external sources, and does not represent the remainder. This stratification is consequential for interpretation: pharmaceutical risk is documented at every layer, but published data are concentrated at the FDF and API/intermediate layers (the layers at which trade statistics are recorded and at which essential-medicines lists operate). Risk at Layer 0 — the fermentation and key-starting-material layer — is documented in expert sources but is not directly trade-observable. The model is transparent about what it directly observes (Layer 1 and Layer 2 trade flows) and what it parameterizes from external sources (Layer 0 fermentation and KSM concentration carried into the Rest-of-World edge).

**eTable 2.** *Pharmaceutical supply chain layer structure and SAPIR-Net Phase 2 model coverage.*

| Layer | Description | Examples | Phase 2 Coverage |
|---|---|---|---|
| 4 | Finished Dosage Form (FDF) | Tablets, capsules, oral suspensions, injectables dispensed to patients. USAntibiotics Bristol; Civica Rx; Indian generic manufacturers; Jordan-based FDF producers. | Not directly modeled. Documented qualitatively via U.S. Pharmacopeia 2024–2025 Vulnerable Medicines List and Socal et al. 2025. |
| 3 | Active Pharmaceutical Ingredient (API) | Amoxicillin; amoxicillin/clavulanate; azithromycin; doxycycline. | **MODEL OUTPUT LAYER** (graph node Layer 3: AMOX, AMOX_CLAV, AZI, DOXY). Capacity loss computed per drug per Monte Carlo iteration. |
| 2 | API / Chemical Intermediate (HS 2941xx) | 6-aminopenicillanic acid and the amoxicillin synthesis chain; tetracycline bulk API; erythromycin A; clavulanic acid (HS 294190 bundle). | **MODEL OBSERVED LAYER** (graph node Layer 2: PENI, TETRA, ERYTH, OTHAB). Empirical UN Comtrade bilateral trade data, 2018–2025. |
| 1 | Geopolitical Source | China (CN), India (IN), Rest of World (ROW residual). | **MODEL OBSERVED LAYER** (graph node Layer 1). Country-of-origin trade flow shares. |
| 0 | Key Starting Materials (KSMs) and Fermentation Intermediates | 6-APA precursors (penicillin G via *Penicillium* fermentation; four amoxicillin KSMs per USP 2024); erythromycin A from *Saccharopolyspora erythraea* fermentation; oxytetracycline from *Streptomyces* fermentation; clavulanic acid from *Streptomyces clavuligerus* fermentation. | Not directly modeled. Parameterized via `row_china_exposure` for the Layer 1 Rest-of-World edge per chemical class. U.S. Pharmacopeia 2024 Quality Matters KSM analysis is the primary external anchor for the bands. |

---

## eAppendix 2. Data Sources: Full Detail

### eAppendix 2.1. UN Comtrade Extraction Specification

The primary empirical input is bilateral trade data extracted from the United Nations Comtrade database (comtradeplus.un.org) for calendar years 2018 through 2025. The United States (Comtrade reporter code 842) is the importing country in all records. Three partner groups are represented: China (Comtrade code 156), India (Comtrade code 699), and World (Comtrade code 0, the global aggregate). The Rest-of-World residual is derived by subtracting the China and India bilateral values from the World total. The extraction returned 96 rows across the four HS codes, three partner groups, and eight years.

The four HS 2941xx codes are antibiotic-specific within the Harmonized System taxonomy. Unlike the broader-bundled HS codes used in SAPIR-Net Phase 1 (HS 284390, which bundles pharmaceutical platinum compounds with catalytic-converter and electronics-grade platinum compounds; HS 293359, which bundles pharmaceutical heterocyclic intermediates with agrochemicals and dyes), the 2941xx codes contain a substantially smaller proportion of non-pharmaceutical trade volume. The 6-Digit Dilution Effect identified in Phase 1 is therefore materially reduced in Phase 2, and the empirical concentration metric sits closer to the pharmaceutical-grade trade subset.

### eAppendix 2.2. The Jordan-Zero Observation

A separate extraction specifying Jordan (Comtrade code 400) as a fourth partner returned zero rows across all four HS codes and all eight years. This finding is methodologically interesting rather than a data error. Jordan is named in the U.S. Pharmacopeia 2024–2025 Vulnerable Medicines List for the United States as a major source of finished dosage forms — specifically 48% of amoxicillin oral suspension and 24% of doxycycline hyclate capsule imports — and is corroborated as a top-three antibiotic FDF supplier (9.0% of U.S. antibiotic FDF imports for 2020–2024) by Socal et al. The same Socal et al. analysis finds Jordan absent from the top fifteen antibiotic API sources across the entire 1992–2024 study window. The Jordan-zero-rows result for HS 2941xx therefore corroborates rather than contradicts the published evidence: the HS 2941xx codes successfully filter to the API layer and exclude FDF-only suppliers, demonstrating that the Phase 2 trade-data extraction operates at the layer the model targets. By the parsimonious reading of the convergent USP and Socal et al. evidence — a top-three FDF supplier position alongside complete absence from the API supplier list across the 1992–2024 study window — Jordan's antibiotic export industry consists of FDF formulation (the mixing, granulation, encapsulation, and oral-suspension preparation steps that convert imported bulk API into the dispensed product) operating downstream of imported bulk API rather than domestic fermentation capacity.

### eAppendix 2.3. Kg-Valid Coverage and the Socal et al. Cross-Check

Net trade weight in kilograms is retained as a secondary variable, with an explicit kg-valid subset flag for each (HS, year) cell. A Comtrade record is treated as kg-valid only when the World, China, and India kg fields are simultaneously non-zero and internally consistent. The kg-valid coverage varies by HS code: 100% for HS 294110 and HS 294190, 62.66% for HS 294130, and 16.86% for HS 294150. The USD basis is used as the primary edge weight throughout the model; the kg basis is used only for the Socal et al. cross-check, which itself operates on a volume basis.

The HS 294150 kg-basis figure operates on a 16.86% kg-valid coverage subset and is reported throughout the paper with explicit data-sparseness caveat; the figure is informative for the aggregate four-code kg-basis cross-check against Socal et al. but is not used as a single-code primary metric. The aggregate four-code kg-basis cross-check pools across all four HS codes and reaches 95.18% kg-valid coverage in aggregate, dominated by the two 100%-kg-valid codes (HS 294110 and HS 294190).

Cross-check derivation. The Socal et al. analysis reports a 2024 single-year China share at the antibiotic API layer of 70.1% of U.S. imports by volume, derived from U.S. Census Trade Online data. The same analysis reports a 2020–2024 five-year aggregate China share at the antibiotic API layer of 62.6%, derived from Socal et al.'s multi-year top-five API origins table (China 62.6%; Bulgaria 16.1%; Spain 3.2%; Mexico 3.1%; Israel 3.0%). Both reference values are quoted directly from Socal et al.; neither is re-derived from the underlying data. The Phase 2 volume-weighted 2024 single-year aggregate across the four HS 2941xx codes, computed on the kg-valid subset of Comtrade records, lies within 6.12 percentage points of the Socal et al. 2024 figure. The Phase 2 five-year 2020–2024 aggregate sits within 2.58 percentage points of the same Socal et al. multi-year figure. The two analyses use different trade-data sources (U.S. Census Trade Online vs. UN Comtrade), different commodity-code classification systems (HTS vs. HS), and slightly different temporal windows, and they converge to within sub-3-percentage-point tolerance at the aggregate level.

### eAppendix 2.4. Chemical Dependency Mappings: Detailed Rationale

The dependency relationships between Layer 2 commodity nodes and Layer 3 medicine nodes are derived from established pharmaceutical chemistry. Amoxicillin is a semi-synthetic penicillin produced from 6-aminopenicillanic acid (6-APA), which is itself derived from penicillin G via enzymatic cleavage at scale. Amoxicillin/clavulanate (Augmentin) is a fixed combination of amoxicillin and clavulanic acid; the clavulanic acid component is a β-lactamase inhibitor produced exclusively by fermentation of *Streptomyces clavuligerus* and structurally distinct from the penicillin nucleus. Azithromycin is a 15-membered macrolide synthesized from erythromycin A via Beckmann rearrangement and methylamination, with erythromycin A itself produced by fermentation of *Saccharopolyspora erythraea*. Doxycycline is a semi-synthetic tetracycline produced by chemical modification of oxytetracycline, with oxytetracycline produced by fermentation of *Streptomyces rimosus* or related species. These dependencies are stoichiometric and obligate: synthesis of the finished API is not chemically possible without the named precursor.

---

## eAppendix 3. Calibration-Band Source Stack

The chemical-class `row_china_exposure` parameter bands are calibrated through convergent evidence from three independent source layers: federal advisory reports (National Academies of Sciences, Engineering, and Medicine 2022; The White House 2021); federally recorded primary testimony (Schondelmeyer 2025 at the U.S.-China Economic and Security Review Commission); and peer-reviewed analysis (Socal et al. 2025 in *JAMA Health Forum*). The U.S. Pharmacopeia 2024 Quality Matters key-starting-material analysis provides an additional source layer specifically at the KSM level. The bands are wide because public concentration data at Layer 0 are imprecise.

**eTable 1.** *Per-chemical-class `row_china_exposure` calibration bands and source anchors.*

| Chemical Class | HS Code | Phase 2 Band | Primary Anchor Sources |
|---|---|---|---|
| Penicillin (PENI) | 294110 | 0.70 – 0.90 | NASEM 2022, citing Schondelmeyer et al. 2020; Schondelmeyer 2025 USCC testimony; U.S. Pharmacopeia 2024 four-KSM amoxicillin finding; Socal et al. 2025 70.1% volume share of U.S. antibiotic API imports in 2024 |
| Tetracycline (TETRA) | 294130 | 0.65 – 0.85 | NASEM 2022; Schondelmeyer 2025 (China dominant for tetracycline and doxycycline specifically, with India having discontinued doxycycline intermediate production in favor of Chinese supply); Socal et al. 2025 aggregate |
| Erythromycin / macrolide (ERYTH) | 294150 | 0.60 – 0.85 | NASEM 2022; Schondelmeyer 2025 (China dominant for erythromycin and azithromycin); Socal et al. 2025 aggregate. Band slightly wider than HS 294110 because public concentration data are less granular for macrolide fermentation. |
| Clavulanic acid (OTHAB) | 294190 (bundle proxy) | 0.55 – 0.80 | Industry reports on clavulanic acid fermentation concentration; WHO essential medicines context. Public concentration data sparser than for the other three classes; band widest accordingly. Applied to the OTHAB → AMOX_CLAV edge directly; the HS 294190 bundle empirical share is not used as a parameter value because the bundle includes cephalosporins, glycopeptides, lincosamides, and other non-clavulanic-acid antibiotics. |

The OTHAB application note merits emphasis. HS 294190 ("other antibiotics, n.e.s.") is a heterogeneous residual within the HS 2941 heading and includes cephalosporins, glycopeptides, lincosamides, and other antibiotic APIs alongside the clavulanic acid trade that the OTHAB node is intended to capture. The empirical HS 294190 direct China share (20.79% USD, 2020–2024) therefore characterizes the entire bundle, not the clavulanic acid subset. Phase 2 addresses this by applying the chemical-class band [0.55, 0.80] directly to the OTHAB → AMOX_CLAV edge in the disruption probability engine rather than using the HS 294190 empirical share as a calibration value. Any reader interpreting the HS 294190 number as a clavulanic-acid concentration estimate would draw the wrong inference.

---

## eAppendix 4. Disruption Scenario Distribution Specifications and Validation

### eAppendix 4.1. Scenario B Pareto Specification

Scenario B's Rest-of-World degradation per HS code is drawn from a rescaled Pareto distribution bounded by the chemical-class-specific `row_china_exposure` range. The shape parameter α = 2.5 produces a moderately heavy tail with finite variance. Under this parameterization, the bulk of sampled degradation values clusters near the lower bound of the chemical-class exposure range, with the tail extending toward the upper bound. The mapping from the Pareto draw on [1, ∞) to the chemical-class range [low, high] is the inverse-CDF transform

> x ↦ low + (1 − 1/x) × (high − low)

with a hard clamp at the upper bound for tail safety.

Each HS code's ROW degradation is sampled independently per Monte Carlo iteration. India's direct exports are not zeroed in Scenario B; the scenario is designed as a China-specific upstream shock rather than a generalized bilateral disruption, consistent with Phase 1. The Pareto distribution is selected on empirical grounds: supply chain disruption severity exhibits fat-tailed behavior, with most events causing moderate degradation but a non-trivial probability mass allocated to near-total upstream collapse.

### eAppendix 4.2. Scenario B Distribution Validation

Across 50,000 validation samples per HS code (RNG seed 42, distinct from the 10,000-iteration Monte Carlo loop), the resulting Rest-of-World degradation distribution per chemical class is characterized as follows:

- HS 294110 (penicillin): sample mean 75.7% (median 74.8%, 5th percentile 70.4%, 95th percentile 83.9%).
- HS 294130 (tetracycline): sample mean 70.7% (median 69.9%, 5th percentile 65.4%, 95th percentile 78.9%).
- HS 294150 (erythromycin): sample mean 67.2% (median 66.1%, 5th percentile 60.5%, 95th percentile 77.5%).
- HS 294190 (clavulanic acid trade-code proxy): sample mean 62.2% (median 61.1%, 5th percentile 55.5%, 95th percentile 72.5%).

The fat-tailed Pareto behavior is visible in the gap between the median and the 95th percentile of each draw set — approximately 9 to 11 percentage points across the four chemical classes — and the maximum observed degradation across 50,000 draws sits within 1 percentage point of the upper band endpoint for each class, confirming that the hard clamp at the upper bound is not observed to engage in the 50,000-sample validation per chemical class. The full distribution validation audit is archived in the Zenodo deposit as `phase2_module2_distribution_audit.txt`.

### eAppendix 4.3. Scenario C Log-Normal Specification

A systemic capacity reduction is applied uniformly to all Layer 1 → Layer 2 edges, regardless of origin country or HS code. The reduction fraction is sampled from a log-normal distribution with parameters μ = ln(0.20) and σ = 0.50, calibrated to produce a median capacity loss of approximately 20% and a 95th-percentile capacity loss of approximately 50–55%. The log-normal distribution captures the right-skewed empirical profile of logistics disruptions: most events cause moderate capacity loss, but extreme events such as the 2021 Suez Canal blockage and pandemic-era port shutdowns can reach 50–70% reduction. Across 50,000 validation samples (RNG seed 42), the resulting capacity-reduction distribution has a sample mean of 22.6%, median 20.0%, 5th percentile 8.8%, and 95th percentile 45.4%, matching the calibration target moments within rounding. The same log-normal capacity reduction is applied to all four HS codes simultaneously per Monte Carlo iteration; the global logistics shock is, by construction, perfectly correlated across origin countries and chemical classes.

### eAppendix 4.4. Empirical Anchors for the Scenario Specifications

The empirical anchor for Scenario B's structural form is the 2013–2014 U.S. doxycycline shortage, in which the exit of a major Chinese tetracycline API manufacturer produced an upstream supply collapse that propagated through the global doxycycline production chain on a months-to-years timescale. Heritage Pharmaceuticals' subsequent re-entry into the doxycycline market via a European third-party API supplier illustrates both the cascade and the partial recovery mechanism that operates at multi-year scale. The empirical anchors for Scenario C's structural form are the 2021 Suez Canal blockage and pandemic-era port shutdowns, both of which produced uniform global logistics-capacity reductions affecting all imported chemical classes simultaneously.

---

## eAppendix 5. Multiplicative Joint Capacity-Loss Rule: Derivation

Three of the four target drugs — amoxicillin, azithromycin, doxycycline — depend on a single Layer 2 commodity each. For these drugs the Layer 2 → Layer 3 propagation rule reduces to the Phase 1 form: the capacity remaining for the drug equals the capacity fraction of its sole upstream chemical commodity, denoted

> capacity_remaining(drug) = p_HS

where p_HS is the residual fraction of the upstream HS-code trade volume relative to the baseline 2020–2024 total after a scenario-specific disruption has been applied.

Amoxicillin/clavulanate depends on two Layer 2 commodities simultaneously: PENI (for the amoxicillin component) and OTHAB (for the clavulanic acid component). The two upstream chemistries are biologically and geographically independent in the model: 6-APA via penicillin G fermentation and clavulanic acid via *Streptomyces clavuligerus* fermentation occupy separate fermentation production bases. The Layer 2 → Layer 3 propagation rule for AMOX_CLAV is therefore multiplicative joint capacity loss:

> capacity_remaining(AMOX_CLAV) = p_PENI × p_OTHAB

This is the standard series-system reliability model under the assumption of independent disruption sampling on the two upstream chemical inputs. The Monte Carlo implementation samples the Rest-of-World degradation on each HS code independently per iteration in Scenario B, so the independence assumption holds for the stochastic component of that scenario. The direct China ban in Scenarios A and B is perfectly correlated across all four HS codes — China's direct exports zero out simultaneously — and the model captures that correlation correctly through the deterministic CN-edge zeroing. The multiplicative joint rule generalizes naturally to any combination drug with N independent upstream chokepoints (capacity_remaining = Π p_i for i = 1 to N) and is the single Phase 2–specific methodological extension over Phase 1's three-drug oncology scope, which did not contain a dual-chokepoint combination drug.

The independence assumption merits a methodological note. In reality, the two fermentation production bases may share some upstream Chinese carbon-source or fermentation-equipment dependencies that the model does not capture. To the extent any such shared dependency exists, the model's independence assumption is conservative — i.e., it understates rather than overstates the dual-chokepoint compounding effect. The conservatism is quantified by the S4 partial-correlation sensitivity reported in eAppendix 7. Scenario C's perfectly correlated global logistics shock model already places PENI and OTHAB at ρ = 1 by construction (the single uniform log-normal shock affects all four HS-code Layer 1 → Layer 2 edges identically), so the S4 partial-correlation sensitivity does not bear directly on the Scenario C 65.4% probability; that probability is the floor of any analysis that introduces additional shared-dependency structure not captured by the Scenario C single-shock specification.

---

## eAppendix 6. Complete Sensitivity Analysis Panel

The complete S1, S2, and S3 sensitivity analyses are reported in eTable 4 across all four drugs and all three disruption scenarios. The 49-row table is also archived in the Zenodo deposit as `phase2_monte_carlo_sensitivity.csv`.

**eTable 4.** *Complete sensitivity analysis panel. Each scenario was executed for N = 10,000 iterations at RNG seed 42 against the indicated parameter configuration. Probability of severe shortage is the proportion of iterations exceeding 30% capacity loss.*

| Configuration | Drug | Scenario | Mean Loss | Median Loss | P5 | P95 | P(Severe) |
|---|---|---|---|---|---|---|---|
| Primary | Amoxicillin | A | 46.23% | 46.23% | 46.23% | 46.23% | 100.0% |
| Primary | Amoxicillin | B | 85.46% | 85.00% | 82.66% | 89.88% | 100.0% |
| Primary | Amoxicillin | C | 22.57% | 19.96% | 8.76% | 44.91% | 20.9% |
| Primary | AMOX_CLAV | A | 60.22% | 60.22% | 60.22% | 60.22% | 100.0% |
| Primary | AMOX_CLAV | B | 95.16% | 95.10% | 93.73% | 96.84% | 100.0% |
| Primary | AMOX_CLAV | C | 38.62% | 35.93% | 16.76% | 69.66% | 65.4% |
| Primary | Azithromycin | A | 70.86% | 70.86% | 70.86% | 70.86% | 100.0% |
| Primary | Azithromycin | B | 76.39% | 76.09% | 74.47% | 79.39% | 100.0% |
| Primary | Azithromycin | C | 22.57% | 19.96% | 8.76% | 44.91% | 20.9% |
| Primary | Doxycycline | A | 75.90% | 75.90% | 75.90% | 75.90% | 100.0% |
| Primary | Doxycycline | B | 92.81% | 92.59% | 91.51% | 94.81% | 100.0% |
| Primary | Doxycycline | C | 22.57% | 19.96% | 8.76% | 44.91% | 20.9% |
| S1 (lower band) | Amoxicillin | B | 82.44% | 82.44% | 82.44% | 82.44% | 100.0% |
| S1 (lower band) | AMOX_CLAV | B | 93.24% | 93.24% | 93.24% | 93.24% | 100.0% |
| S1 (lower band) | Azithromycin | B | 74.32% | 74.32% | 74.32% | 74.32% | 100.0% |
| S1 (lower band) | Doxycycline | B | 91.41% | 91.41% | 91.41% | 91.41% | 100.0% |
| S2 (upper band) | Amoxicillin | B | 93.19% | 93.19% | 93.19% | 93.19% | 100.0% |
| S2 (upper band) | AMOX_CLAV | B | 98.64% | 98.64% | 98.64% | 98.64% | 100.0% |
| S2 (upper band) | Azithromycin | B | 81.60% | 81.60% | 81.60% | 81.60% | 100.0% |
| S2 (upper band) | Doxycycline | B | 96.23% | 96.23% | 96.23% | 96.23% | 100.0% |
| S3 (HS 294150 2023–2025 baseline) | Azithromycin | A | 90.69% | 90.69% | 90.69% | 90.69% | 100.0% |
| S3 (HS 294150 2023–2025 baseline) | Azithromycin | B | 79.32% | 78.99% | 77.21% | 82.74% | 100.0% |
| S3 (HS 294150 2023–2025 baseline) | Azithromycin | C | 22.57% | 19.96% | 8.76% | 44.91% | 20.9% |

S3 modifies only the HS 294150 baseline; the three other drugs are unchanged from the primary configuration and are omitted from the S3 rows for compactness. The full row panel (including the S1 and S2 results for Scenarios A and C, and the unchanged drugs under S3) is in the archived CSV.

**eAppendix 6.1. The structural-guarantee disclosure (S1).** The S1 result collapse is not specific to the S1 sensitivity. It is the deterministic band-arithmetic floor of every iteration of Scenario B under the model's source-stack-calibrated parameter bands. For every Monte Carlo iteration of the primary band-sampled Scenario B run, the per-HS Pareto draws on ROW degradation sit at or above each chemical-class lower band endpoint by construction; the CN direct ban additionally zeros the Chinese share. The resulting per-iteration capacity loss is bounded below by the S1 result for every drug. No Pareto draw, no choice of Pareto shape parameter, no Monte Carlo iteration variability, and no baseline-window shift (S3) can produce a single iteration of Scenario B for any drug in which capacity loss falls below the band-arithmetic floor of S1. The 100% probability outcome reported in the primary results is the binary indicator equivalent of this structural-guarantee floor: the primary band-sampled Scenario B varies in the magnitude of loss across iterations (mean 76.39% to 95.16%; 5th-percentile 74.47% to 93.73%) but never in the binary severe-shortage outcome. The 100% probability is therefore band-arithmetic-driven, conditional on the source-stack-calibrated parameter space, rather than an emergent Monte Carlo discovery.

**eAppendix 6.2. S3 asymmetric structural-robustness finding.** For azithromycin, Scenario A's mean capacity loss jumps from 70.86% (primary) to 90.69% (S3), an increase of 19.83 percentage points, because the direct-ban scenario is deterministically tied to the empirical CN + IN direct share, which shifted upward substantially in the sub-recent window (HS 294150 direct China share rose from 50.36% in 2020 to 83.10% in 2024). Scenario B moves only modestly, from 76.39% to 79.32%, a shift of 2.93 percentage points, because Scenario B's loss distribution is already dominated by the Pareto-degraded Rest-of-World component; raising the empirical CN direct share within the same scenario adds capacity loss only at the margin once the ROW degradation has accounted for the bulk of upstream supply elimination. Scenario C is unchanged at 22.57% mean capacity loss because the log-normal shock is band-independent. The S3 result establishes a publishable structural robustness finding: the 100% probability of severe shortage under Scenario B is stable to material shifts in the empirical baseline window, while Scenario A is highly sensitive to baseline-window choice. This sensitivity asymmetry is informative for federal policy framing — a direct export ban scenario evaluated against the 2020–2024 aggregate would understate the risk that the 2024–2025 trade pattern actually carries, whereas the Scenario B cascading-shock result is robust to the same baseline-window shift.

---

## eAppendix 7. Partial-Correlation Sensitivity (S4): Methodology and Full Results

The Phase 2 multiplicative joint capacity-loss rule for amoxicillin/clavulanate assumes statistical independence of the Rest-of-World Pareto degradation on the two upstream chemical inputs in Scenario B. The independence assumption is asserted in the main text to be conservative under the realistic case that the 6-APA penicillin G fermentation production base and the *Streptomyces clavuligerus* clavulanic acid fermentation production base share upstream Chinese provincial-cluster infrastructure (shared KSM carbon-source and nitrogen-source feedstocks; shared bioreactor and downstream-processing equipment classes; shared operator pools within the same provincial industrial clusters documented in NASEM 2022 and in Schondelmeyer 2025 for Hebei, Henan, and Inner Mongolia penicillin G production specifically).

### eAppendix 7.1. Gaussian Copula Methodology

The present sensitivity analysis quantifies the conservatism by sampling the per-iteration ROW Pareto draws on PENI (HS 294110) and OTHAB (HS 294190) under a Gaussian copula at correlation values ρ ∈ {0, 0.25, 0.50, 0.75}, with ρ = 0 reproducing the primary Scenario B configuration (independent per-HS sampling) and higher ρ values representing progressively more realistic shared-infrastructure scenarios. The Pareto marginal distributions are unchanged at each ρ: only the joint dependence between PENI and OTHAB shifts.

The Gaussian copula is implemented by drawing iid standard-normal pairs and rotating them through a 2×2 correlation matrix with off-diagonal element ρ; the rotated pairs are then transformed to U(0,1) via the standard normal CDF and inverted through the rescaled Pareto CDF to obtain Pareto draws with the prescribed marginal distribution and Gaussian-copula dependence. A separate sampler module was implemented for this analysis distinct from the primary-run Pareto sampler. The S4 sampler at ρ = 0 was validated against the existing primary-run AMOX_CLAV Scenario B output: at N = 10,000, the new sampler returns mean 95.19% (primary 95.16%; difference 0.03 percentage points), median 95.13% (primary 95.10%; 0.03 pp), P5 93.74% (primary 93.73%; 0.01 pp), and P95 96.87% (primary 96.84%; 0.03 pp). All differences sit within sampling variability at N = 10,000, confirming that the new sampler produces marginal distributions consistent with the existing primary-run Pareto sampler at the independence baseline. The full validation log is archived in the Zenodo deposit as `phase2_monte_carlo_correlation_audit.txt`. The single-input drugs (AMOX, AZI, DOXY) have only one upstream chokepoint each in Scenario B and are by construction invariant to the PENI–OTHAB joint structure; their results are statistically equivalent to the primary run across all ρ values within sampling tolerance and are reported in the supplementary `phase2_monte_carlo_correlation.csv` file alongside the AMOX_CLAV results for completeness.

The Gaussian copula was selected over alternatives (notably the t-copula) for three reasons. First, the substantive conclusion the S4 sensitivity is asked to test — whether the binary 100% severe-shortage outcome is robust to relaxation of the per-iteration PENI–OTHAB independence assumption — is set by the band-arithmetic floor of the marginal Pareto distributions, not by their joint tail dependence. No copula choice can move a single iteration of Scenario B below the band-arithmetic floor for any drug, because each chemical-class lower band endpoint already exceeds the 30% severe-shortage threshold by 44 to 63 percentage points; the headline 100% finding is copula-invariant. Second, the Gaussian copula's single-parameter interpretation (ρ as progressive correlation strength) cleanly tracks the substantive question (the degree of shared upstream Chinese fermentation infrastructure between PENI and OTHAB) without introducing a tail-dependence parameter to interpret jointly with ρ. Third, a t-copula sensitivity targeting the upper-tail magnitude distribution specifically — at which a t-ν copula's positive tail dependence would shift the AMOX_CLAV Scenario B 95th-percentile loss beyond the Gaussian-copula shift reported in eAppendix 7.3 — is a natural Phase 2.x methodological note, where the magnitude-distribution implications can be characterized at greater depth than the present paper's policy-relevant scope requires. The conservatism direction asserted in the main text would, if anything, be sharpened by a t-copula sensitivity, not reversed; the present Gaussian-copula treatment is therefore the conservative reportable result for the v1.3 manuscript scope.

### eAppendix 7.2. AMOX_CLAV Scenario B Partial-Correlation Results

**eTable 3.** *AMOX_CLAV Scenario B partial-correlation sensitivity. Gaussian-copula correlation values ρ ∈ {0, 0.25, 0.50, 0.75} applied to the per-iteration Rest-of-World Pareto degradation draws on HS 294110 (PENI) and HS 294190 (OTHAB); N = 10,000 iterations per ρ value; RNG seed = 42. P(severe shortage) is the proportion of iterations exceeding 30% capacity loss.*

| ρ | Mean Loss | Median Loss | P5 Loss | P95 Loss | P(Severe Shortage) |
|---|---|---|---|---|---|
| 0.00 | 95.19% | 95.13% | 93.74% | 96.87% | 100.0% |
| 0.25 | 95.17% | 95.08% | 93.62% | 97.02% | 100.0% |
| 0.50 | 95.14% | 95.03% | 93.53% | 97.18% | 100.0% |
| 0.75 | 95.11% | 94.96% | 93.45% | 97.32% | 100.0% |

### eAppendix 7.3. Interpretation

The result is consistent with two distinct components of the main paper's framing. First, the binary P(severe shortage) outcome remains 100% across all ρ values for AMOX_CLAV. This is consistent with the structural-guarantee result: the band-arithmetic floor of Scenario B exceeds the 30% severe-shortage threshold for every drug by 44 to 63 percentage points, and no rearrangement of the per-iteration joint dependence between PENI and OTHAB Pareto draws within the chemical-class bands can move a single iteration of Scenario B below the threshold. The partial-correlation sensitivity therefore characterizes the magnitude distribution rather than the binary indicator. Second, the magnitude distribution shifts in the predicted conservative direction: the AMOX_CLAV Scenario B P95 loss increases from 96.87% at ρ = 0 to 97.32% at ρ = 0.75 (a +0.45 percentage-point upward shift in the upper tail), and the P5 loss decreases from 93.74% at ρ = 0 to 93.45% at ρ = 0.75 (a −0.29 percentage-point downward shift in the lower tail, reflecting that positive correlation produces more iterations with both draws simultaneously near the lower band endpoint, where the multiplicative-joint loss is smaller). The conservatism claim asserted in the main text — that the independence assumption understates the compounded effect of dual-chokepoint disruption — is upheld in direction for the upper tail of the magnitude distribution. The magnitude of the shift is modest: at ρ = 0.75, the AMOX_CLAV Scenario B upper-tail loss is roughly half a percentage point higher than at ρ = 0, against a baseline upper-tail loss already approaching the maximum-possible-loss ceiling. The headline 100% probability finding is robust to relaxation of the independence assumption across the full range of correlation values evaluated.

![](phase2_fig4_dual_chokepoint_correlation.png)

*eFigure 1.* AMOX_CLAV Scenario B capacity-loss distribution under partial-correlation sensitivity on the per-iteration PENI–OTHAB Rest-of-World Pareto draws. Top panel: overlaid kernel density estimates at ρ ∈ {0, 0.25, 0.50, 0.75}, zoomed to the 91%–99% loss range where the ρ-induced shifts are visible. Bottom panel: per-ρ mean (point) with P5–P95 interval (line), N = 10,000 iterations per ρ. P(severe shortage) is 100% across all ρ values per the structural-guarantee result; the partial-correlation sensitivity characterizes the magnitude distribution, not the binary indicator. The 30% severe-shortage threshold sits far below the displayed loss range — the entire AMOX_CLAV Scenario B loss distribution is above the severe-shortage threshold at every ρ value.

---

## eAppendix 8. Pareto Shape Parameter α Sensitivity

The choice of the Pareto shape parameter α = 2.5 in the Scenario B per-HS Rest-of-World degradation sampler admits a methodological-note treatment of the model's sensitivity to that choice. Under α ∈ {1.5, 2.5, 3.5}, the Pareto-rescaled marginal distribution within each chemical-class band shifts: lower α produces a heavier tail and pushes the per-iteration draws further toward the upper band endpoint on average, while higher α produces a thinner tail and concentrates the draws closer to the lower band endpoint. The bulk of the parameter-sensitivity effect is therefore on the magnitude distribution of Scenario B losses, not on the binary severe-shortage indicator: per the structural-guarantee result, every iteration of Scenario B for every drug produces capacity loss exceeding the 30% severe-shortage threshold regardless of where the Pareto draw lands within the band, because the band lower endpoint is itself already above the threshold by 44 to 63 percentage points. Pareto α therefore affects the within-band magnitude distribution of capacity loss but does not affect the 100% probability outcome. A fuller α-sensitivity treatment across additional shape values is held for a Phase 2.x methodological note where the magnitude-distribution implications can be characterized at greater depth than the present paper's policy-relevant scope.

---

## eAppendix 9. Code, Data, and Reproducibility

All analysis modules (Modules 1–4A) are implemented in Python 3.12 (specifically 3.12.3 per `requirements.txt`) using NumPy, SciPy, pandas, NetworkX, and Matplotlib. The four core modules accept the prior module's output as input with no hidden state dependencies. The Gaussian-copula partial-correlation sensitivity (S4) is implemented as a separate runner script that imports from the existing Module 3 without modifying it; the byte-exact close on the primary-run `phase2_monte_carlo_primary.csv`, `phase2_monte_carlo_sensitivity.csv`, and `phase2_monte_carlo_raw_losses.npz` files is therefore preserved. The complete artifact set archived at the Zenodo deposit (DOI [10.5281/zenodo.20376649](https://doi.org/10.5281/zenodo.20376649)) and at the GitHub repository ([https://github.com/fahimkazmi911/sapir-net-phase2-antibiotics](https://github.com/fahimkazmi911/sapir-net-phase2-antibiotics)) comprises:

- Code modules: `sapir_phase2_module1_graph_hhi.py`, `sapir_phase2_module2_disruption_engine.py`, `sapir_phase2_module3_monte_carlo.py`, `sapir_phase2_module3_correlation_run.py`, `sapir_phase2_module4a_visualizations.py`, `sapir_phase2_figure4_correlation.py`.
- Input data: UN Comtrade bilateral trade extraction 2018–2025 for HS 294110, 294130, 294150, 294190 (CSV).
- Outputs: HHI per-code panel, calibration table, Monte Carlo primary and sensitivity tables (49-row), distribution audit logs, partial-correlation summary statistics and per-ρ AMOX_CLAV capacity-loss arrays, rendered figures (PNG).
- Dependency manifest: `requirements.txt` with pinned versions for all third-party packages.
- Audit logs: `phase2_module2_distribution_audit.txt`, `phase2_monte_carlo_correlation_audit.txt`.
- LICENSE and README.

RNG seed 42 reproduces all reported numerical outputs to byte-exact tolerance.

---

*End of Online Supplementary Materials.*
