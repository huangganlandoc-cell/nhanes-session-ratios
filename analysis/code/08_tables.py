"""
稿件表格（Markdown + CSV）。数字全部读自分析数据集与 results_REAL.json（第 2 次运行），不手抄。
输出：analysis/tables/
"""
import json
import numpy as np
import pandas as pd
from config import PROJECT, DERIVED

R = json.load(open(PROJECT / "analysis/results/real/results_REAL.json", encoding="utf-8"))
TAB = PROJECT / "analysis/tables"
TAB.mkdir(parents=True, exist_ok=True)


def wmean(x, w):
    m = x.notna()
    return np.average(x[m], weights=w[m])


def wsd(x, w):
    m = x.notna(); mu = np.average(x[m], weights=w[m])
    return np.sqrt(np.average((x[m] - mu) ** 2, weights=w[m]))


def smd_cont(a, b, wa, wb):
    return (wmean(a, wa) - wmean(b, wb)) / np.sqrt((wsd(a, wa) ** 2 + wsd(b, wb) ** 2) / 2)


def md(df):
    cols = list(df.columns)
    out = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    out += ["| " + " | ".join(str(v) for v in r) + " |" for r in df.itertuples(index=False)]
    return "\n".join(out)


# ---------------- Table 1：NHANES III 随机分组基线（加权） ----------------
n3 = pd.read_pickle(DERIVED / "nhanes3_analytic.pkl")
itt = n3[n3.in_itt].copy()
itt["women"] = (itt.sex == 2).astype(float)
itt["educ_lt12"] = np.where(itt.educ_years.isna(), np.nan, (itt.educ_years < 12).astype(float))
itt["pir_lt1"] = np.where(itt.pir.isna(), np.nan, (itt.pir < 1).astype(float))
itt["health_fairpoor"] = np.where(itt.self_health.isna(), np.nan, (itt.self_health >= 4).astype(float))
itt["cvd_any"] = itt[["hx_mi", "hx_chf", "hx_stroke"]].max(axis=1, skipna=True)
A, B = itt[itt["arm"] == "AM"], itt[itt["arm"] == "PMEVE"]
rows = [("Participants, n (unweighted)", f"{len(A):,}", f"{len(B):,}", "")]


def cont_row(label, v, fmt="{:.1f}"):
    return (label, f"{fmt.format(wmean(A[v], A.w_mec))} ({fmt.format(wsd(A[v], A.w_mec))})",
            f"{fmt.format(wmean(B[v], B.w_mec))} ({fmt.format(wsd(B[v], B.w_mec))})",
            f"{smd_cont(B[v], A[v], B.w_mec, A.w_mec):+.3f}")


def pct_row(label, series_a, series_b):
    pa, pb = wmean(series_a, A.w_mec), wmean(series_b, B.w_mec)
    s = (pb - pa) / np.sqrt((pa * (1 - pa) + pb * (1 - pb)) / 2)
    return (label, f"{100 * pa:.1f}", f"{100 * pb:.1f}", f"{s:+.3f}")


rows.append(cont_row("Age, y", "age"))
rows.append(pct_row("Women, %", A.women, B.women))
for code, lab in [(1, "Non-Hispanic White"), (2, "Non-Hispanic Black"), (3, "Mexican American"), (4, "Other")]:
    rows.append(pct_row(f"Race and ethnicity: {lab}, %", (A.race_eth == code).astype(float), (B.race_eth == code).astype(float)))
rows.append(pct_row("Survey phase 1 (1988–1991), %", (A.phase == 1).astype(float), (B.phase == 1).astype(float)))
rows.append(pct_row("Education <12 years, %", A.educ_lt12, B.educ_lt12))
rows.append(pct_row("Poverty–income ratio <1, %", A.pir_lt1, B.pir_lt1))
for s_, lab in [("current", "Current smoker"), ("former", "Former smoker")]:
    rows.append(pct_row(f"{lab}, %", (A.smoke == s_).astype(float), (B.smoke == s_).astype(float)))
rows.append(cont_row("Body mass index, kg/m²", "bmi"))
rows.append(pct_row("Self-rated health fair or poor, %", A.health_fairpoor, B.health_fairpoor))
rows.append(pct_row("Diabetes (self-reported), %", A.hx_diabetes, B.hx_diabetes))
rows.append(pct_row("Hypertension (self-reported), %", A.hx_hypertension, B.hx_hypertension))
rows.append(pct_row("Cardiovascular disease (self-reported), %", A.cvd_any, B.cvd_any))
rows.append(pct_row("Cancer other than skin (self-reported), %", A.hx_other_cancer, B.hx_other_cancer))
rows.append(("Examination characteristics (after assignment)", "", "", ""))
rows.append(pct_row("Examined in assigned session, %", (A.session == 1).astype(float), B.session.isin([2, 3]).astype(float)))
for v, lab in [("fast_h", "Fasting duration, h"), ("draw_clock_h", "Clock time of venipuncture, h")]:
    rows.append((f"{lab}, median (IQR)",
                 f"{A[v].median():.1f} ({A[v].quantile(.25):.1f}–{A[v].quantile(.75):.1f})",
                 f"{B[v].median():.1f} ({B[v].quantile(.25):.1f}–{B[v].quantile(.75):.1f})", ""))
t1 = pd.DataFrame(rows, columns=["Characteristic", "Morning assignment", "Afternoon/evening assignment", "SMD"])
t1.to_csv(TAB / "Table1_baseline_by_arm.csv", index=False)
note1 = ("Values are survey-weighted means (SD) or percentages (NHANES III MEC examination weights); participant counts are "
         "unweighted; examination-characteristic medians are unweighted. SMD, standardized mean difference "
         "(afternoon/evening minus morning).")

# ---------------- Table 2：随机分组效应（Part 1）与敏感性分析 ----------------
p1 = R["part1"]
f = lambda d: f"{d['pct']:+.2f} ({d['pct_lo']:+.2f} to {d['pct_hi']:+.2f})"
rows = [("GLR (primary)", f(p1["primary_glr_taylor"]), "±9.5", p1["primary_glr_taylor"]["R1"])]
names = {"wbc": "Leukocytes", "lym": "Lymphocytes", "gran": "Granulocytes", "mono_nuc": "Mononuclear cells", "plt": "Platelets",
         "sii_g": "SII (granulocyte-based)", "plr": "PLR"}
for k, lab in names.items():
    d = p1["secondary"][k]
    rows.append((lab, f(d), "—" if d["threshold"] is None else f"±{d['threshold']:.1f}", d["R1"] or "—"))
c = p1["cace"]
rows.append(("GLR, complier average causal effect", f(c), "±9.5", "—"))
sn = {"S1_unweighted": "S1 Unweighted", "S2_adjusted": "S2 Covariate-adjusted", "S3_per_protocol": "S3 Per protocol",
      "S4_ipw_cbc": "S4 Inverse probability weighting for missing blood count",
      "S5_excl_possible_infection": "S5 Excluding possible infection at examination"}
for k, lab in sn.items():
    rows.append((f"GLR, {lab}", f(p1["sensitivity"][k]), "±9.5", p1["sensitivity"][k]["R1"]))
t2 = pd.DataFrame(rows, columns=["Measure", "Difference, % (95% CI)", "Desirable bias, %", "Prespecified judgment (R1)"])
t2.to_csv(TAB / "Table2_randomized_effects.csv", index=False)
rc = p1["reclassification"]
reclass = pd.DataFrame([(k, f"{100 * v['p0']:.1f}", f"{100 * v['p1']:.1f}", f"{100 * v['rd']:+.1f} ({100 * v['rd_lo']:+.1f} to {100 * v['rd_hi']:+.1f})")
                        for k, v in rc.items()], columns=["Classification", "Morning, %", "Afternoon/evening, %", "Difference, pp (95% CI)"])

# ---------------- Table 3：死亡关联 ----------------
p2 = R["part2"]
rows = []
for s_, lab in [("nhanes3", "NHANES III (GLR)"), ("continuous", "NHANES 1999–2018 (NLR)")]:
    for ev, evl in [("death_all", "All-cause"), ("death_cancer", "Cancer")]:
        d = p2[s_][ev]
        rows.append((lab, evl, f"{d['events']:,}", f"{d['hr_obs_per_doubling']:.3f}", f"{d['hr_std_per_doubling']:.3f}",
                     f"{d['delta_pct']:+.2f} ({d['lo_pct']:+.2f} to {d['hi_pct']:+.2f})", d["R2"]))
t3 = pd.DataFrame(rows, columns=["Cohort (index)", "Outcome", "Deaths", "HR per doubling, observed", "HR per doubling, standardized",
                                 "Change in log HR, % (95% CI)", "Prespecified judgment (R2)"])
t3.to_csv(TAB / "Table3_mortality.csv", index=False)
s3 = p2["nhanes3_secondary"]
hrci = pd.DataFrame([(ev, f"{s3[ev]['x_obs']['hr_per_doubling']:.3f} ({s3[ev]['x_obs']['lo']:.3f}–{s3[ev]['x_obs']['hi']:.3f})",
                      f"{s3[ev]['x_std']['hr_per_doubling']:.3f} ({s3[ev]['x_std']['lo']:.3f}–{s3[ev]['x_std']['hi']:.3f})",
                      f"{s3[ev]['x_obs']['hr_per_sd']:.3f}", f"{s3[ev]['q4_vs_q1']['hr_obs']:.3f} → {s3[ev]['q4_vs_q1']['hr_std']:.3f}",
                      f"{100 * s3[ev]['q4_vs_q1']['share_changing_quartile']:.1f}", f"{s3[ev]['ph_global_p']:.2g}")
                     for ev in ["death_all", "death_cancer"]],
                    columns=["Outcome", "HR (95% CI), observed", "HR (95% CI), standardized", "HR per SD, observed",
                             "Q4 vs Q1 HR, observed → standardized", "Changed quartile, %", "PH global P"])

doc = ["# Tables (draft; generated by analysis/code/08_tables.py)", "",
       "## Table 1. Characteristics of NHANES III participants by randomized examination-session assignment", "", md(t1), "", note1, "",
       "## Table 2. Effect of randomized session assignment (NHANES III, intention-to-treat)", "", md(t2), "",
       "Percentage differences are afternoon/evening versus morning assignment (survey-weighted; Taylor linearization). "
       "Desirable bias from the EFLM Biological Variation Database (GLR and SII derived from component estimates).", "",
       "### Table 2 (continued). Classification above cutoffs by randomized assignment", "", md(reclass), "",
       "## Table 3. Session standardization and associations with mortality", "", md(t3), "",
       "### Table 3 (continued). NHANES III, design-based hazard ratios per doubling of GLR", "", md(hrci), ""]
(TAB / "tables_draft.md").write_text("\n".join(doc), encoding="utf-8")
print("写出", TAB / "tables_draft.md")
