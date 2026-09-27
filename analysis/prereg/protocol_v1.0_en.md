# Analysis protocol (preregistration) v1.0 — final, for OSF registration

- Drafted: 2026-09-26 (drafted with AI assistance); approved by the author: 2026-09-26 (all items in Section 9 as recommended)
- Status: **final**. This English version is a sentence-by-sentence translation of the author-approved Chinese version (`analysis/analysis_plan.md`, v1.0); both versions are attached to the OSF registration, and the SHA-256 of both files is recorded in `analysis/prereg/HASHES.txt`. None of the analyses in Section 5 will be run before OSF registration.
- Working title: Examination session and fasting duration as sources of variation in leukocyte-ratio inflammation indices: a randomized-session analysis of NHANES III with replication in NHANES 1999–2023

---

## 0. What had been seen before this protocol was locked (transparency statement)

As required for preregistration of secondary data analyses, we declare our prior knowledge of the data. As of this version:

**Already seen (all at the design level; none answers the research questions)**
- Sample sizes, inclusion flow, numbers of deaths, and follow-up time;
- Univariate distributions of the indices in the whole sample (not by group), e.g., GLR median 1.908, SD of ln(GLR) 0.446;
- Balance of pre-assignment covariates between randomized groups (maximum |SMD| 0.039), randomized group versus session actually attended (compliance 95.4% / 96.2%), and randomized group versus fasting duration and clock time of blood draw (manipulation check);
- In continuous NHANES, the relation between the session attended and covariates (maximum |SMD| 0.037);
- Agreement between the Coulter granulocyte-to-lymphocyte ratio (GLR) and the manual-differential NLR (Spearman ρ 0.896, non-random subsample).

**Not seen**
- Any blood-cell index summarized by session, clock time of blood draw, or fasting duration (means, quantiles, or proportions);
- Any association between an index and mortality;
- The scripts of the preliminary feasibility assessment (2026-09-26) were checked one by one; they likewise computed only the design-level quantities listed above.

---

## 1. Background and aims (brief)

Leukocyte-ratio indices such as NLR, SII, SIRI, and PIV are used extensively in secondary analyses of NHANES (SII/NLR/SIRI/PIV × NHANES: 694 articles in total, 225 in the past 12 months), and these analyses generally do not consider the time of blood collection. Leukocyte subsets have well-established diurnal rhythms, and NHANES examinations are conducted in morning, afternoon, and evening sessions. The official NHANES III documentation states that **households were randomly assigned to morning ("standard") or afternoon/evening ("modified") sessions** (`refs/official_docs/lab-acc.txt`, lines 4530–4533), which amounts to a randomized experiment embedded in a national survey. This random assignment has previously been used to evaluate diurnal differences in fasting glucose (Troisi 2000, JAMA, PMID 11135780) and ALT (Ruhl & Everhart 2013, J Clin Gastroenterol, PMID 23164687); it has not been applied to leukocyte-ratio indices.

**Aims**
1. To estimate, using the NHANES III random assignment, the causal effect of examination session on GLR (the three-part-differential approximation of NLR) and related indices, and to judge whether it exceeds the laboratory-medicine "desirable bias" threshold;
2. To translate the effect into the number of people reclassified around commonly used cutoffs;
3. To assess how much session-related measurement variation affects index–mortality HRs;
4. To replicate the findings by actual session in continuous NHANES (1999–2018, 2021–2023).

---

## 2. Data and design

| | NHANES III (1988–1994) | Continuous NHANES 1999–2018 | Continuous NHANES 2021–2023 |
|---|---|---|---|
| Session information | Random assignment (reconstructed from whether WTPFSD6/WTPFMD6 are blank) + actual session + clock time of blood draw | Actual session (morning/afternoon/evening); no clock time, no assignment variable | Actual session (morning / afternoon or evening) |
| Complete blood count | Coulter S-Plus JR three-part differential (lymphocytes, mononuclear cells, granulocytes) | Five-part differential | Five-part differential |
| Mortality follow-up | 2019 public-use linkage, median 26.2 years | 2019 public-use linkage, median 108 months | None |
| Role in this study | **Primary design (randomized)** | Replication (non-randomized) | Cross-sectional replication |

Rule for reconstructing assignment (as described in the official documentation): for persons assigned to the morning session, WTPFSD6 is non-blank (>0 if examined in the assigned session, =0 if examined in another session) and WTPFMD6 is blank; the reverse holds for persons assigned to the afternoon/evening session. In the ITT sample, 0 persons have both fields non-blank, and 0 MEC-examined adults lack an assignment record.

---

## 3. Study population (numbers from `analysis/output/data_notes.md`)

**NHANES III**
- ITT sample (Part 1 primary analysis): aged ≥20 years, MEC-examined, randomly assigned, no chemotherapy in the past 4 weeks, valid leukocyte/lymphocyte/granulocyte/platelet counts → **15,476 persons** (morning group 7,724; afternoon/evening group 7,752).
  - Clock time, fasting duration, and BMI were measured after assignment; smoking and other questionnaire variables are covariates that the ITT primary estimate does not use. None of these is an **ITT inclusion criterion** (the 15,416 persons in the preliminary assessment additionally required these variables to be valid).
- Clock time/fasting analysis (Part 1b): ITT participants with valid clock time of blood draw and fasting duration → 15,459 persons; persons who fasted ≥24 hours (55) are excluded from this part.
- Mortality analysis (Part 2): ITT participants eligible for mortality linkage → 15,466 persons; 6,634 all-cause deaths, 1,449 cancer deaths.

**Continuous NHANES**
- 1999–2018 main sample: aged ≥20 years, MEC-examined, session recorded, complete five-part differential and platelet count, not pregnant, valid fasting duration <24 hours → **48,021 persons**; mortality analysis 47,934 persons (7,654 all-cause deaths, 1,662 cancer deaths).
- 2021–2023: same criteria → 5,672 persons (cross-sectional only).

---

## 4. Variable definitions

**Primary outcome index**
- NHANES III: GLR = granulocyte count / lymphocyte count (Coulter). Granulocytes include neutrophils, eosinophils, and basophils; GLR is approximately 5.6% higher than NLR.
- Continuous NHANES: NLR = neutrophil count / lymphocyte count.

**Secondary indices**: leukocytes, lymphocytes, granulocytes (or neutrophils), mononuclear cells (NHANES III; not equivalent to monocytes)/monocytes, platelets; SII (in NHANES III, the granulocyte version SII_g = platelets × granulocytes / lymphocytes); PLR; in continuous NHANES additionally SIRI, PIV, and (neutrophils + eosinophils + basophils)/lymphocytes, defined as in NHANES III. CRP: in NHANES III, 65.5% of values are at the limit of detection, so only the binary variable "detectable (>0.21 mg/dL)" is used; in continuous NHANES, LBXCRP is used for 1999–2010 and hs-CRP for 2015–2018 and 2021–2023, with the platforms analyzed separately.

All ratios and counts are analyzed on the natural-log scale; effects are reported as the percentage difference of afternoon/evening relative to morning: (e^β − 1) × 100%.

**Exposures**
- NHANES III: randomized group Z (afternoon/evening = 1, morning = 0); actual session D (afternoon or evening = 1).
- Continuous NHANES: actual session (afternoon and evening, each versus morning; for 2021–2023, afternoon/evening versus morning).

**Covariates (core set, defined as similarly as possible in both surveys)**: age (spline with 4 degrees of freedom), sex, race/ethnicity, survey phase or cycle, smoking (never/former/current/missing), BMI (spline with 3 degrees of freedom + missing indicator), education (categories + missing), poverty-income ratio (<1, 1–<2, 2–<4, ≥4, missing); in continuous NHANES, additionally the examination season. Extended set (sensitivity analyses only): additionally self-reported history of cardiovascular disease and of cancer.

**Mortality outcomes**: all-cause death; cancer death (UCOD_LEADING = 002, cause-specific hazard, other causes of death censored); follow-up starts at the MEC examination and ends on 2019-12-31.

---

## 5. Statistical analysis

General settings: all primary estimates use the sample design (NHANES III: strata SDPSTRA6, PSU SDPPSU6, weight WTPFEX6; continuous NHANES 1999–2018: 20-year weights combined according to NCHS rules, stratum numbers prefixed by cycle; 2021–2023: phlebotomy weight WTPH2YR). Variance: Taylor linearization for linear models; for quantities that require comparing two estimates within the same resample (Δ in Section 5.4), the official Fay BRR replicate weights WTPXRP1–52 (ρ = 0.3) are used in NHANES III and the delete-one-PSU jackknife (JKn) in continuous NHANES. Degrees of freedom: 49 for NHANES III. Random seed 20260926 (used only for the exploratory bootstrap).

### 5.1 Part 1 primary analysis: effect of random session assignment (NHANES III, ITT)

- **Primary estimate**: weighted linear regression ln(GLR) ~ Z, without covariate adjustment. The percentage difference and 95% CI are reported.
- **Scales reported alongside** (no separate decision thresholds):
  - standardized shift = β / SD(ln GLR);
  - partial R² (proportion of the variance of ln GLR explained by Z) with 95% CI (BRR);
  - reclassification ("consequence" endpoints): difference between groups in the proportion with GLR ≥3.0; difference between groups in the proportion above the weighted 75th percentile of the whole sample ("top quartile"); the same for SII_g. Additionally, the proportion of afternoon/evening participants classified as high that is attributable to session = (p₁ − p₀)/p₁.
- **CACE** (effect among those who attended their assigned session): Wald estimator β_ITT / [P(D=1|Z=1) − P(D=1|Z=0)], with CI from BRR.
- **Sensitivity analyses**
  - S1 unweighted;
  - S2 adjusted for the core covariates (for precision);
  - S3 per-protocol analysis (only those examined in their assigned session, using the official session weights WTPFSD6/WTPFMD6);
  - S4 the proportion with a valid CBC differs by 0.9 percentage points between groups (93.9% vs 93.0%); the primary estimate is repeated with inverse probability weighting (valid CBC predicted from group and core covariates);
  - S5 excluding persons for whom the examining physician recorded "possible infection" (PEP13C = 2; coding confirmed in the NHANES III examination codebook: 1 no infection, 2 possible infection, 8 blank but applicable).
- **Heterogeneity (exploratory)**: stratified by sex, age (20–39, 40–59, ≥60), race/ethnicity, and smoking; interaction P values are reported, without conclusions.
- **Manual-differential NLR**: this subsample consists of "a random 10% plus examinees with instrument results beyond predetermined limits" (official documentation), i.e., it was selected on outcome values; it is **not used for the ITT analysis** and serves only to describe the agreement between GLR and NLR.

### 5.2 Part 1b clock time of blood draw and fasting duration (NHANES III, non-randomized, exploratory)

- ln(GLR) ~ natural spline(clock time, 4 df) + natural spline(fasting hours, 3 df) + core covariates; the adjusted diurnal curve is plotted, with reference points 09:00 and 12 hours of fasting.
- The model is repeated within each randomized group to reduce confounding.
- Description of "how much of the session effect accompanies fasting duration": the change in β after adding the fasting spline to the ITT model (fasting occurs after assignment, so this is descriptive only and not a causal decomposition).

### 5.3 Part 1c replication in continuous NHANES (non-randomized)

- 1999–2018: weighted linear regression ln(NLR) ~ session (afternoon, evening vs morning) + core covariates + cycle + examination season. Primary replication quantity: percentage difference, afternoon vs morning.
- The same approach is applied to secondary indices, GLR defined as in NHANES III, and CRP (by platform); consequence endpoints: adjusted differences in the proportions with NLR ≥3.0, SII ≥445.22 (the inflection point of the U-shaped curve between SII and cancer mortality reported by Wu 2025, BMC Public Health, PMID 39833806), and top quartile (marginal standardization).
- 2021–2023: afternoon/evening vs morning, same model.
- Secondary model: additionally including a spline of fasting duration.

### 5.4 Part 2 influence of session-related variation on index–mortality HRs

- **Standardized index**: for participants examined in the afternoon/evening session, the session effect is subtracted from the ln index, converting it to the "morning blood draw" scale. In NHANES III, the CACE from Part 1 × actual session is used; in continuous NHANES, the adjusted effects from Part 1c (afternoon and evening separately).
- **Model**: weighted Cox regression, time scale months of follow-up, core covariates; the exposure is log₂(index) (HR per doubling). The observed-index and standardized-index models have exactly the same covariates.
- **Primary quantity**: Δ = (β_standardized − β_observed) / β_observed, i.e., the relative change in the ln HR per doubling.
- **Paired resampling**: within each replicate-weight/jackknife sample, **the session effect is re-estimated**, the index re-standardized, and both Cox models refitted to obtain Δ_r; these give the 95% CI of Δ.
- **Analytical prediction**: attenuation factor λ = 1 − δ²·p(1−p) / Var(ln index) (δ is the session effect, p the afternoon/evening proportion); the predicted Δ ≈ 1/λ − 1 is reported alongside the observed Δ.
- **Outcomes**: all-cause mortality (precision anchor) and cancer mortality (oncological primary outcome), each judged separately. NHANES III is primary; continuous NHANES 1999–2018 is the replication.
- **Secondary (no decision rules)**: HR per SD; HR for the highest vs lowest quartile and the proportion of people who change quartile after standardization; shift of the restricted cubic spline (knots at the 5th/35th/65th/95th percentiles) inflection point before vs after standardization (in percentile units, CI from 500 bootstrap replicates); in NHANES 2005–2018, the change in the SII 445.22 inflection point before vs after standardization; the proportional hazards assumption is checked with Schoenfeld residuals and reported.

### 5.5 Secondary and exploratory analyses

- Triglycerides (NHANES III TGP; continuous NHANES LBXSTR) as a marker of the postprandial state, described cross-sectionally as a mediator, exploratory only;
- A limited-grid specification curve (continuous NHANES): index {NLR, SII, SIRI, PIV} × exposure coding {per doubling, per SD, Q4 vs Q1, dichotomous cutoff (NLR 3.0, SII 445.22; SIRI and PIV have no accepted cutoff, so the weighted median is used)} × covariate set {core, extended} × outcome {all-cause, cancer}, 64 specifications in total; the distribution of Δ caused by session standardization is reported;
- Cancer survivors (1999–2018, 4,469 with a self-reported cancer history) and female breast cancer survivors (675, 73 cancer deaths): only the Part 1c session differences and descriptive HRs are reported, without judgment;
- Unintentional injury death (182 in NHANES III, 229 in continuous NHANES) as a negative control outcome, reported qualitatively only (power is insufficient for testing);
- The serum calcium negative control has been dropped (conclusion of the preliminary assessment).

### 5.6 Missing data, software, and dual implementation

- The ITT primary analysis requires no covariates; in adjusted models, missing covariate values are handled with missing categories.
- Primary engine R 4.6.1 (survey, survival). Key numbers (Part 1 primary estimate, CACE, both βs and Δ in Part 2) are programmed again in Python by an independent sub-agent that has no access to the R code and works only from this protocol; results enter interpretation only when point estimates agree to 4 decimal places and standard errors differ by <5%.
- **Blinded development**: before registration, all code is written and run end-to-end on "scrambled data" — in NHANES III the group labels are randomly permuted within strata, in continuous NHANES the session is permuted within cycle × stratum, and the mortality outcomes are permuted across the whole population; all output is labeled "SCRAMBLED – NOT REAL". After registration, the code is run once on the real data; any subsequent code change is recorded as a protocol deviation with date and reason.

---

## 6. Decision rules (fixed before any results are seen)

### R1: Is the effect of session on GLR clinically meaningful? (Part 1 primary outcome)

Threshold = the laboratory-medicine "desirable bias" = 0.25 × √(CVI² + CVG²). Under a normal distribution, a shift of 0.25 population standard deviations changes the proportions outside the two sides of the 95% reference interval from 2.5% each to approximately 4.4% and 1.4% (about 5.7% in total, versus 5%); the principle of setting analytical quality goals from the effect of bias on the transferability of reference intervals is described by Gowans et al. (1988, PMID 3238321). Biological variation is taken from the whole-blood meta-analysis estimates of the EFLM Biological Variation Database (accessed 2026-09-26; database updated 2026-02-11):

| Measurand | CVI | CVG | Desirable bias |
|---|---|---|---|
| Leukocytes | 11.1% | 16.9% | 5.1% |
| Neutrophils | 12.5% | 25.9% | 7.2% |
| Lymphocytes | 10.5% | 22.8% | 6.3% |
| Monocytes | 14.0% | 23.0% | 6.7% |
| Platelets | 7.3% | 18.6% | 5.0% |
| **NLR / GLR (derived)** | 16.3% | 34.5% | **9.5%** |
| **SII (derived)** | 17.9% | 39.2% | **10.8%** |

No published biological variation estimates exist for NLR and SII; they are derived by adding the component variances on the log scale (assuming independence between components; if neutrophils and lymphocytes are positively correlated between individuals, the derived CVG is overestimated and the threshold too wide, so a "strictest" case is additionally reported: the smallest component desirable bias, 6.3%, as a sensitivity threshold).

| Position of the 95% CI relative to threshold m (percentage difference) | Judgment |
|---|---|
| CI entirely within (−m, +m) | Session effect **within the desirable bias** (not clinically meaningful) |
| CI lower bound > +m, or upper bound < −m | Session effect **exceeds the desirable bias** (clinically meaningful) |
| CI crosses +m or −m | **Indeterminate**: a clinically meaningful shift cannot be excluded |

The primary judgment uses m = 9.5% (GLR); SII_g uses 10.8%; each component uses its own threshold (secondary). The tier (optimal 0.125×, desirable 0.25×, minimum 0.375×) in which the point estimate falls is also reported.

### R2: Is the influence of session standardization on mortality HRs negligible? (Part 2; judged separately for all-cause and cancer mortality)

Threshold: relative change in ln HR of ±10% (the epidemiological "10% change-in-estimate" convention).

| 95% CI of Δ | Judgment |
|---|---|
| Entirely within (−10%, +10%) | HR **robust** to session-related variation |
| Lower bound > +10% or upper bound < −10% | HR **materially affected** by session-related variation |
| Crosses ±10% | **Indeterminate** |

### R3: Does continuous NHANES replicate NHANES III?

| NHANES III (R1 result) | Continuous NHANES 1999–2018, afternoon vs morning (NLR) | Judgment |
|---|---|---|
| Exceeds desirable bias | Same direction and 95% CI excludes 0 | Replicated |
| Exceeds desirable bias | Opposite direction and CI excludes 0 | Inconsistent (differences between the surveys to be discussed) |
| Exceeds desirable bias | CI includes 0 | Not replicated |
| Within desirable bias | CI entirely within ±9.5% | Replicated (consistent null) |
| Within desirable bias | CI entirely outside ±9.5% | Inconsistent |
| Within desirable bias | CI crosses ±9.5% | Indeterminate |
| Indeterminate | — | Direction and magnitude described only; no "replication" conclusion |

Magnitudes are compared using GLR defined identically, descriptively only. Results for 2021–2023 are descriptive only.

---

## 7. Prespecified wording of conclusions (by R1 × R2 combination)

| R1 | R2 | Main conclusion of the article |
|---|---|---|
| Exceeds desirable bias | Robust | Examination session produces a clinically meaningful shift in GLR/NLR and reclassifies X% of people around cutoffs, but has a negligible influence on mortality HRs. The concern lies in the transferability of cutoffs and reference intervals, not in the association estimates themselves. |
| Exceeds desirable bias | Materially affected | Session changes both the index values and the HRs materially; analyses of NHANES leukocyte ratios need to account for session. |
| Within desirable bias | Robust | Under randomized evidence, NHANES examination session has no clinically meaningful influence on GLR/NLR — an adequately powered null result for the existing literature. |
| Within desirable bias | Materially affected | Should not occur according to the attenuation formula; treated as a signal of an analytical error, code and data are checked first, and no conclusion is drawn from it. |
| Indeterminate | Any | The CI is reported, with the wording "a clinically meaningful shift cannot be excluded", not "no influence". |

The only relation that can be calculated in advance is that between R2 and the size of the session effect: by the attenuation formula, Δ is about 4% for a session effect of 20% and reaches 10% only at about 31%. The result of R1 cannot be calculated in advance.

Wording of conclusions: the NHANES III session effect derives from random assignment, so "effect of examination-session assignment" may be written; continuous NHANES is an observational comparison, for which only "differed by session" is written. "Bias" and "predicts" are not used in the title or abstract.

---

## 8. Differences from the preliminary assessment plan (item 19 of the topic-screening supplement, "required changes")

| Item | Preliminary assessment | This protocol | Reason |
|---|---|---|---|
| ITT inclusion criteria | Required valid clock time, fasting duration, BMI, and smoking (15,416 persons) | Only pre-assignment criteria + valid CBC (15,476 persons) | Variables measured after assignment should not determine the ITT sample; covariates are not used in the unadjusted primary estimate |
| Decision rule 1 | Partial R² CI + SD shift + reclassification proportion (no numerical threshold) | The three are reported as planned, with the EFLM desirable bias added as the numerical threshold | A number fixed before seeing results is needed; desirable bias is the standard yardstick in laboratory medicine for judging pre-analytical variation |
| Paired resampling | Bootstrap by PSU within strata | Official BRR replicate weights in NHANES III; JKn in continuous NHANES | Official variance methods; deterministic, reproducible results |
| True NLR | Not addressed | Manual-differential subsample used only to describe agreement | The subsample was selected on outcome values and cannot be used for ITT |
| Exclusion of fasting ≥24 hours | Excluded from the whole sample | Not excluded from the NHANES III ITT sample; excluded from Part 1b and continuous NHANES | Fasting occurs after assignment |

All other elements (all-cause mortality as precision anchor, cancer mortality as oncological primary outcome, dropping the serum calcium control, qualitative reporting of accidental deaths, exploratory triglycerides, descriptive breast cancer subgroup, no "bias" in the title) are unchanged from the preliminary assessment.

---

## 9. Decisions made (approved by the author on 2026-09-26)

1. **R1 threshold**: EFLM desirable bias — GLR/NLR 9.5%, SII 10.8%; the smallest component desirable bias, 6.3%, is additionally reported as a sensitivity threshold.
2. **R2 threshold**: relative change in ln HR of ±10%.
3. **Part 2 exposure scale**: per doubling (log₂) as primary, per SD as secondary.
4. **Covariates**: core set; diabetes and hypertension are not added (DIQ and BPQ files are not additionally downloaded). The extended set adds only self-reported history of cardiovascular disease and of cancer.
5. **S5**: retained; PEP13C coding confirmed against the official codebook.
6. **OSF registration**: submitted by the author's own account as an Open-Ended Registration, with the full Chinese and English protocols attached; an embargo may be set without affecting the registration timestamp.

---

## 10. Locking procedure

1. Author approval of the Chinese version — completed (2026-09-26);
2. Sentence-by-sentence translation into English, attached to the OSF registration;
3. SHA-256 of both files computed and recorded in `analysis/prereg/HASHES.txt` and `PROJECT_MEMORY.md`;
4. The author submits the OSF registration and obtains the registration timestamp and DOI;
5. Code development on scrambled data (may proceed in parallel with steps 1–4);
6. After registration, a single run on the real data, judged according to Section 6 and concluded according to Section 7.

---

## Appendix: sources

- NHANES III laboratory data file documentation: `refs/official_docs/lab-acc.txt` (random session assignment: lines 4524–4578; BRR replicate weights: from line 4656; definition of the manual-differential subsample: lines 4842–4856)
- NHANES III examination data file codebook: `refs/official_docs/exam-acc.txt` (PEP13A and PEP13C coding, positions 1481 and 1483)
- Continuous NHANES fasting questionnaire codebooks: `refs/official_docs/PH.htm`, `FASTQX_J.htm`, `FASTQX_L.htm`
- EFLM Biological Variation Database (https://biologicalvariation.eu/, accessed 2026-09-26); citation format: Aarsand AK, et al. The EFLM Biological Variation Database. Original meta-analysis of haematological parameters: Coşkun A, et al. Clin Chem Lab Med 2019;58:25–32 (PMID 31503541)
- Normal NLR values: Forget P, et al. BMC Res Notes 2017;10:12 (PMID 28057051), healthy adults 0.78–3.53
- SII 445.22: Wu S, et al. BMC Public Health 2025;25:227 (PMID 39833806), inflection point of the U-shaped curve between SII and cancer mortality in NHANES 2005–2018
- Design precedents: Troisi 2000 JAMA (PMID 11135780); Ruhl & Everhart 2013 J Clin Gastroenterol (PMID 23164687)
- Desirable bias and transferability of reference intervals: Gowans EM, et al. Scand J Clin Lab Invest 1988 (PMID 3238321)
- All references above will undergo citation verification at the writing stage (Stage 23).
