"""
数据探查（Stage 14a）：流程图人数、变量可得性、随机化平衡与依从、单变量分布。

盲态规则：本脚本只输出
  - 人数、缺失、单变量（全样本合并）分布；
  - 随机分组/场次与"分配前协变量"的平衡（设计诊断，不涉及结局）；
  - 分组与实际场次、禁食时长、抽血钟点之间的关系（暴露—暴露，属于操纵检验）。
不输出任何按场次/禁食/钟点分组的血细胞指标，也不输出指标与死亡的关系。

输出：analysis/DATA_TABLE_NOTES.md 的数字部分（analysis/output/data_notes.md），以及 CSV 明细。
"""
import numpy as np
import pandas as pd
from config import DERIVED, OUT

OUTD = OUT / "data_notes"
OUTD.mkdir(parents=True, exist_ok=True)
lines = []
P = lines.append


def smd_cont(x1, x0, w1=None, w0=None):
    m1, m0 = np.average(x1, weights=w1), np.average(x0, weights=w0)
    v1 = np.average((x1 - m1) ** 2, weights=w1)
    v0 = np.average((x0 - m0) ** 2, weights=w0)
    return (m1 - m0) / np.sqrt((v1 + v0) / 2)


def balance(df, group, g1, g0, specs, wcol=None):
    """specs: [(显示名, 列名, 'cont'|'cat')]；返回每行的均值/比例与 SMD（未加权、加权）。"""
    rows = []
    for label, col, kind in specs:
        sub = df[df[col].notna() & (df[col] != "missing")] if kind == "cat" else df[df[col].notna()]
        a, b = sub[sub[group] == g1], sub[sub[group] == g0]
        wa = a[wcol] if wcol else None
        wb = b[wcol] if wcol else None
        if kind == "cont":
            rows.append((label, "", f"{a[col].mean():.2f}", f"{b[col].mean():.2f}",
                         smd_cont(a[col], b[col]), smd_cont(a[col], b[col], wa, wb) if wcol else np.nan,
                         len(a), len(b)))
        else:
            for lev in sorted(sub[col].unique(), key=str):
                xa, xb = (a[col] == lev).astype(float), (b[col] == lev).astype(float)
                rows.append((label, str(lev), f"{100 * xa.mean():.1f}%", f"{100 * xb.mean():.1f}%",
                             smd_cont(xa, xb), smd_cont(xa, xb, wa, wb) if wcol else np.nan, len(a), len(b)))
    return pd.DataFrame(rows, columns=["变量", "水平", g1, g0, "SMD_未加权", "SMD_加权", f"n_{g1}", f"n_{g0}"])


def md_table(df, floatfmt="{:.3f}"):
    cols = list(df.columns)
    out = ["| " + " | ".join(map(str, cols)) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        out.append("| " + " | ".join(floatfmt.format(v) if isinstance(v, float) else str(v) for v in r) + " |")
    return "\n".join(out)


# =========================== NHANES III ===========================
n3 = pd.read_pickle(DERIVED / "nhanes3_analytic.pkl")
P("## A. NHANES III（1988–1994）\n")
P("### A1. 纳入流程（≥20 岁起算）\n")
s = n3
flow = [("实验室文件中 ≥20 岁", len(s))]
for flag, label in [("ok_mec", "MEC 受检（排除家访体检）"), ("ok_assigned", "有场次随机分配记录"),
                    ("ok_not_chemo", "排除近 4 周化疗"), ("ok_cbc_ok", "白细胞/淋巴/粒细胞/血小板有效")]:
    s = s[s[flag]]
    flow.append((label, len(s)))
flow.append(("= ITT 样本（Part 1 主分析）", int(n3.in_itt.sum())))
flow.append(("其中抽血钟点与禁食时长均有效（Part 1b）", int(n3.in_timing.sum())))
itt = n3[n3.in_itt].copy()
mort = itt[(itt.eligstat == 1) & itt.fu_months_exam.notna()]
flow.append(("其中可链接死亡（eligstat=1 且随访时间有效，Part 2）", len(mort)))
P(md_table(pd.DataFrame(flow, columns=["步骤", "人数"])))
P("\n未进入 ITT 的 MEC 受检成人中，没有分配记录者："
  f"{int((n3.ok_mec & n3.ok_age20 & ~n3.ok_assigned).sum())} 人。\n")

P("### A2. 随机分配与实际场次（依从）\n")
ct = pd.crosstab(itt["arm"], itt.session.map({1: "上午", 2: "下午", 3: "晚上"}), margins=True, margins_name="合计")
P(md_table(ct.reset_index().rename(columns={"arm": "分配"})))
for arm, s1 in [("AM", 1)]:
    pass
am = itt[itt["arm"] == "AM"]
pm = itt[itt["arm"] == "PMEVE"]
P(f"\n- 分配上午者实际在上午：{(am.session == 1).mean() * 100:.1f}%；"
  f"分配下午/晚上者实际在下午或晚上：{pm.session.isin([2, 3]).mean() * 100:.1f}%")
P(f"- 分配场次的权重为 0（=未按分配到场）：上午组 {int((am.w_am == 0).sum())} 人，下午/晚上组 {int((pm.w_pm == 0).sum())} 人；"
  "与上面交叉表的不依从人数应一致。")
P(f"- 家庭序号（DMPFSEQ）不同值：{itt.family_seq.nunique()}；设计层 {itt.strata.nunique()} 个，每层 PSU {sorted(itt.psu.unique().tolist())}\n")

P("### A3. 操纵检验：分配组与禁食时长、抽血钟点（暴露—暴露，不涉及结局）\n")
rows = []
for arm, g in itt.groupby("arm"):
    f_ = g.fast_h.dropna()
    c_ = g.draw_clock_h.dropna()
    rows.append((arm, len(g), f"{f_.median():.1f} ({f_.quantile(.25):.1f}–{f_.quantile(.75):.1f})",
                 f"{(f_ >= 8).mean() * 100:.1f}%", f"{c_.min():.2f}–{c_.max():.2f}", f"{c_.median():.2f}"))
P(md_table(pd.DataFrame(rows, columns=["分配", "n", "禁食小时 中位(IQR)", "禁食≥8h", "抽血钟点范围", "钟点中位"])))
P(f"\n禁食 ≥24 小时：{int((itt.fast_h >= 24).sum())} 人（最长 {itt.fast_h.max():.1f} 小时）。\n")

P("### A4. 随机化平衡（分配前协变量；SMD 绝对值 >0.10 视为不平衡）\n")
itt["sex_lbl"] = itt.sex.map({1: "男", 2: "女"})
itt["race_lbl"] = itt.race_eth.map({1: "非西语裔白人", 2: "非西语裔黑人", 3: "墨西哥裔", 4: "其他"})
itt["phase_lbl"] = itt.phase.map({1: "1988–91", 2: "1991–94"})
specs = [("年龄", "age", "cont"), ("性别", "sex_lbl", "cat"), ("种族/族裔", "race_lbl", "cat"),
         ("调查阶段", "phase_lbl", "cat"), ("贫困收入比", "pir", "cont"), ("受教育年限", "educ_years", "cont"),
         ("自评健康(1极好–5差)", "self_health", "cont"), ("吸烟", "smoke", "cat"),
         ("糖尿病史", "hx_diabetes", "cont"), ("高血压史", "hx_hypertension", "cont"), ("心梗史", "hx_mi", "cont"),
         ("心衰史", "hx_chf", "cont"), ("卒中史", "hx_stroke", "cont"), ("其他癌症史", "hx_other_cancer", "cont"),
         ("BMI（MEC 测，分配后）", "bmi", "cont")]
bal = balance(itt, "arm", "PMEVE", "AM", specs, wcol="w_mec")
bal.to_csv(OUTD / "nhanes3_balance.csv", index=False)
P(md_table(bal))
P(f"\n最大 |SMD|（未加权）= {bal.SMD_未加权.abs().max():.3f}；（加权）= {bal.SMD_加权.abs().max():.3f}\n")

P("### A5. ITT 样本中各变量可得性（不分组）\n")
avail = [("白细胞", "wbc"), ("淋巴细胞", "lym"), ("粒细胞", "gran"), ("单个核细胞(Coulter)", "mono_nuc"),
         ("血小板", "plt"), ("GLR", "glr"), ("SII 粒细胞版", "sii_g"), ("人工分类 NLR", "nlr_manual"),
         ("CRP（mg/dL）", "crp"), ("纤维蛋白原", "fibrinogen"), ("甘油三酯", "tg"), ("血清铁", "iron"),
         ("BMI", "bmi"), ("贫困收入比", "pir"), ("受教育年限", "educ_years")]
rows = [(lab, int(itt[c].notna().sum()), f"{itt[c].notna().mean() * 100:.1f}%") for lab, c in avail]
rows.append(("吸烟状态有效", int((itt.smoke != "missing").sum()), f"{(itt.smoke != 'missing').mean() * 100:.1f}%"))
P(md_table(pd.DataFrame(rows, columns=["变量", "有效 n", "占比"])))
crp = itt.crp.dropna()
P(f"\nCRP 处于检测下限 0.21 mg/dL：{(crp <= 0.21).mean() * 100:.1f}%（{int((crp <= 0.21).sum())}/{len(crp)}）\n")

P("### A6. 指标单变量分布（全 ITT 样本合并，不分组）\n")
rows = []
for lab, c in [("GLR", "glr"), ("SII 粒细胞版", "sii_g"), ("PLR", "plr"), ("人工分类 NLR", "nlr_manual"),
               ("白细胞", "wbc"), ("淋巴细胞", "lym"), ("粒细胞", "gran"), ("血小板", "plt")]:
    x = itt[c].dropna()
    x = x[x > 0]
    rows.append((lab, len(x), f"{x.median():.3f}", f"{x.quantile(.25):.3f}–{x.quantile(.75):.3f}",
                 f"{x.quantile(.99):.2f}", f"{np.log(x).std():.3f}"))
P(md_table(pd.DataFrame(rows, columns=["指标", "n", "中位数", "IQR", "P99", "SD(ln)"])))
g = itt.glr.dropna()
P(f"\nGLR 落在 3.0 上下 ±15%（2.55–3.45）内：{((g >= 2.55) & (g <= 3.45)).mean() * 100:.1f}%；GLR ≥3.0：{(g >= 3).mean() * 100:.1f}%\n")
# 测量效度：Coulter GLR 与人工分类 NLR 的一致性（不涉及场次）。
# 注意：人工分类子样本 = 随机 10% + 所有仪器结果超预设界值者（官方文档），按结局值选择，只作描述。
mm = itt[itt.nlr_manual.notna() & (itt.nlr_manual > 0)]
rho = mm[["glr", "nlr_manual"]].corr(method="spearman").iloc[0, 1]
ratio = (mm.glr / mm.nlr_manual).median()
P(f"Coulter GLR 与人工分类 NLR（n={len(mm)}，非随机子样本）：Spearman ρ = {rho:.3f}；GLR/NLR 比值中位数 {ratio:.3f}\n")

P("### A7. 死亡随访（Part 2 样本，不分组）\n")
P(f"- 人数 {len(mort)}；全因死亡 {int(mort.death_all.sum())}；癌症死亡 {int(mort.death_cancer.sum())}；"
  f"意外伤害死亡 {int(mort.death_accident.sum())}")
P(f"- 随访（自体检起）中位 {mort.fu_months_exam.median() / 12:.1f} 年；总人年 {mort.fu_months_exam.sum() / 12:,.0f}")
P(f"- 女性 {int((mort.sex == 2).sum())} 人（癌症死亡 {int(mort[mort.sex == 2].death_cancer.sum())}）；"
  f"自报其他癌症史 {int((mort.hx_other_cancer == 1).sum())} 人\n")

# =========================== 连续 NHANES ===========================
c = pd.read_pickle(DERIVED / "continuous_analytic.pkl")
P("## B. 连续 NHANES（1999–2018 合并；2021–2023 单列）\n")
for per, sub in [("1999–2018", c[c.cycle != "2021-2023"]), ("2021–2023", c[c.cycle == "2021-2023"])]:
    s = sub
    flow = [("≥20 岁 MEC 受检", len(s))]
    for flag, label in [("ok_session", "场次有记录"), ("ok_cbc", "五分类血常规+血小板完整"),
                        ("ok_not_preg", "排除妊娠"), ("ok_fast", "禁食时长有效且 <24 小时")]:
        s = s[s[flag]]
        flow.append((label, len(s)))
    P(f"### B1. 纳入流程：{per}\n")
    P(md_table(pd.DataFrame(flow, columns=["步骤", "人数"])))
    P("")

main = c[c.in_main].copy()
P("### B2. 各周期人数与场次分布（暴露分布，不涉及结局）\n")
lab_s = {0: "上午", 1: "下午", 2: "晚上"}
main["sess_lbl"] = np.where(main.cycle == "2021-2023", main.session.map({0: "上午", 1: "下午/晚上"}),
                            main.session.map(lab_s))
ct = pd.crosstab(main.cycle, main.sess_lbl, margins=True, margins_name="合计")
P(md_table(ct.reset_index()))
P("")

m18 = main[main.cycle != "2021-2023"].copy()
P("### B3. 1999–2018：场次与协变量（非随机，仅描述；SMD 为下午或晚上 vs 上午）\n")
m18["sex_lbl"] = m18.sex.map({1: "男", 2: "女"})
m18["race_lbl"] = m18.race_eth.map({1: "墨西哥裔", 2: "其他西语裔", 3: "非西语裔白人", 4: "非西语裔黑人", 5: "其他"})
m18["pmeve"] = np.where(m18.session == 0, "AM", "PMEVE")
specs_c = [("年龄", "age", "cont"), ("性别", "sex_lbl", "cat"), ("种族/族裔", "race_lbl", "cat"),
           ("贫困收入比", "pir", "cont"), ("教育(1–5)", "educ", "cont"), ("吸烟", "smoke", "cat"),
           ("BMI", "bmi", "cont"), ("心血管病史", "hx_cvd", "cont"), ("癌症史", "hx_cancer", "cont"),
           ("体检季节(1冬/2夏)", "exam_season", "cont")]
balc = balance(m18, "pmeve", "PMEVE", "AM", specs_c, wcol="w_mec_pooled")
balc.to_csv(OUTD / "continuous_session_covariates.csv", index=False)
P(md_table(balc))
P(f"\n最大 |SMD|（加权）= {balc.SMD_加权.abs().max():.3f}\n")

P("### B4. 禁食时长与场次（暴露—暴露）\n")
rows = []
for k_, g in main.groupby("sess_lbl"):
    f_ = g.fast_h
    rows.append((k_, len(g), f"{f_.median():.1f} ({f_.quantile(.25):.1f}–{f_.quantile(.75):.1f})",
                 f"{(f_ >= 8).mean() * 100:.1f}%"))
P(md_table(pd.DataFrame(rows, columns=["场次", "n", "禁食小时 中位(IQR)", "禁食≥8h"])))
P(f"\n上午空腹子样本权重 >0：{int((m18.w_fast_sub > 0).sum())} 人（1999–2018 主样本内）\n")

P("### B5. 可得性与单变量分布（1999–2018 主样本，不分组）\n")
rows = []
for lab, col in [("NLR", "nlr"), ("SII", "sii"), ("SIRI", "siri"), ("PIV", "piv"), ("PLR", "plr"), ("GLR 口径", "glr")]:
    x = m18[col].dropna()
    x = x[x > 0]
    rows.append((lab, len(x), f"{x.median():.3f}", f"{x.quantile(.25):.3f}–{x.quantile(.75):.3f}", f"{np.log(x).std():.3f}"))
P(md_table(pd.DataFrame(rows, columns=["指标", "n", "中位数", "IQR", "SD(ln)"])))
nlr = m18.nlr.dropna()
sii = m18.sii.dropna()
P(f"\n- NLR 3.0 上下 ±15%（2.55–3.45）内：{((nlr >= 2.55) & (nlr <= 3.45)).mean() * 100:.1f}%；NLR ≥3.0：{(nlr >= 3).mean() * 100:.1f}%")
P(f"- SII 445.22 上下 ±15% 内：{((sii >= 445.22 * .85) & (sii <= 445.22 * 1.15)).mean() * 100:.1f}%；"
  f"445.22 位于第 {(sii < 445.22).mean() * 100:.1f} 百分位")
P(f"- CRP（1999–2010，LBXCRP）有效 {int(m18.crp_1999_2010.notna().sum())}；hs-CRP（2015–2018）有效 {int(m18.hscrp.notna().sum())}；"
  f"2021–2023 hs-CRP 有效 {int(main[main.cycle == '2021-2023'].hscrp.notna().sum())}")
P(f"- 甘油三酯：参考法 LBXTR 有效 {int(m18.tg_ref.notna().sum())}；生化 LBXSTR 有效 {int(m18.tg_biochem.notna().sum())}")
P(f"- 协变量缺失：BMI {int(m18.bmi.isna().sum())}；吸烟 {int((m18.smoke == 'missing').sum())}；PIR {int(m18.pir.isna().sum())}；教育 {int(m18.educ.isna().sum())}\n")

P("### B6. 死亡随访（1999–2018 主样本，不分组）\n")
e = m18[(m18.eligstat == 1) & m18.fu_months_exam.notna()]
P(f"- 人数 {len(e)}；全因死亡 {int(e.death_all.sum())}；癌症死亡 {int(e.death_cancer.sum())}；意外伤害死亡 {int(e.death_accident.sum())}")
P(f"- 随访中位 {e.fu_months_exam.median():.0f} 个月；总人年 {e.fu_months_exam.sum() / 12:,.0f}")
byc = e.groupby("cycle").death_cancer.sum().astype(int).to_dict()
P(f"- 各周期癌症死亡：{byc}")
P(f"- 自报癌症史 {int((e.hx_cancer == 1).sum())} 人（癌症死亡 {int(e[e.hx_cancer == 1].death_cancer.sum())}）；"
  f"女性乳腺癌史 {int((e.hx_breast_cancer & (e.sex == 2)).sum())} 人（癌症死亡 {int(e[e.hx_breast_cancer & (e.sex == 2)].death_cancer.sum())}）\n")

(OUT / "data_notes.md").write_text("\n".join(lines), encoding="utf-8")
print("\n".join(lines))
