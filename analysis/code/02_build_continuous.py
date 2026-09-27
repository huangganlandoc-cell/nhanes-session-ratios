"""
构建连续 NHANES（1999–2018 合并 + 2021–2023 单列）分析数据集。

连续 NHANES 的公开文件没有"分配到哪个场次"的变量，只有实际场次（PHDSESN：0 上午, 1 下午, 2 晚上；
2021–2023 为 PHDSESNZ：0 上午, 1 下午/晚上），也没有抽血钟点，因此只作按实际场次的复现分析。

不做什么：不计算任何"场次/禁食 → 指标"或"指标 → 死亡"的关联（方案锁定前保持盲态）。

输出：analysis/derived/continuous_analytic.pkl（全部 ≥20 岁记录，含纳入标记）
"""
import numpy as np
import pandas as pd
from config import XPT, MORT, DERIVED

# 各周期文件名：(周期, DEMO, 血常规, 禁食问卷, 甘油三酯参考法, MCQ, BMX, SMQ, CRP, 生化, 死亡文件后缀)
CYCLES = [
    ("1999-2000", "DEMO", "LAB25", "PH", "LAB13AM", "MCQ", "BMX", "SMQ", "LAB11", "LAB18", "1999_2000"),
    ("2001-2002", "DEMO_B", "L25_B", "PH_B", "L13AM_B", "MCQ_B", "BMX_B", "SMQ_B", "L11_B", "L40_B", "2001_2002"),
    ("2003-2004", "DEMO_C", "L25_C", "PH_C", "L13AM_C", "MCQ_C", "BMX_C", "SMQ_C", "L11_C", "L40_C", "2003_2004"),
    ("2005-2006", "DEMO_D", "CBC_D", "FASTQX_D", "TRIGLY_D", "MCQ_D", "BMX_D", "SMQ_D", "CRP_D", "BIOPRO_D", "2005_2006"),
    ("2007-2008", "DEMO_E", "CBC_E", "FASTQX_E", "TRIGLY_E", "MCQ_E", "BMX_E", "SMQ_E", "CRP_E", "BIOPRO_E", "2007_2008"),
    ("2009-2010", "DEMO_F", "CBC_F", "FASTQX_F", "TRIGLY_F", "MCQ_F", "BMX_F", "SMQ_F", "CRP_F", "BIOPRO_F", "2009_2010"),
    ("2011-2012", "DEMO_G", "CBC_G", "FASTQX_G", "TRIGLY_G", "MCQ_G", "BMX_G", "SMQ_G", None, "BIOPRO_G", "2011_2012"),
    ("2013-2014", "DEMO_H", "CBC_H", "FASTQX_H", "TRIGLY_H", "MCQ_H", "BMX_H", "SMQ_H", None, "BIOPRO_H", "2013_2014"),
    ("2015-2016", "DEMO_I", "CBC_I", "FASTQX_I", "TRIGLY_I", "MCQ_I", "BMX_I", "SMQ_I", "HSCRP_I", "BIOPRO_I", "2015_2016"),
    ("2017-2018", "DEMO_J", "CBC_J", "FASTQX_J", "TRIGLY_J", "MCQ_J", "BMX_J", "SMQ_J", "HSCRP_J", "BIOPRO_J", "2017_2018"),
    ("2021-2023", "DEMO_L", "CBC_L", "FASTQX_L", None, "MCQ_L", "BMX_L", "SMQ_L", "HSCRP_L", None, None),
]

KEEP = {
    "demo": ["RIDSTATR", "RIAGENDR", "RIDAGEYR", "RIDRETH1", "DMDEDUC2", "INDFMPIR", "RIDEXPRG", "RIDEXMON",
             "WTMEC2YR", "WTMEC4YR", "SDMVPSU", "SDMVSTRA"],
    "cbc": ["WTPH2YR", "LBXWBCSI", "LBDLYMNO", "LBDNENO", "LBDMONO", "LBDEONO", "LBDBANO", "LBXPLTSI"],
    "fast": ["PHAFSTHR", "PHAFSTMN", "PHDSESN", "PHDSESNZ"],
    "tg": ["WTSAF2YR", "WTSAF4YR", "LBXTR"],
    "mcq": ["MCQ220", "MCQ230A", "MCQ230B", "MCQ230C", "MCQ230D", "MCQ160B", "MCQ160C", "MCQ160E", "MCQ160F"],
    "bmx": ["BMXBMI", "BMXHT"],
    "smq": ["SMQ020", "SMQ040"],
    "crp": ["LBXCRP", "LBXHSCRP", "LBDHRPLC"],
    "bio": ["LBXSTR"],
}


def rd(name, cols):
    t = pd.read_sas(XPT / f"{name}.xpt", format="xport", encoding="latin-1")
    t = t[["SEQN"] + [c for c in cols if c in t.columns]]
    assert t.SEQN.is_unique, name
    return t


def read_mort(suffix):
    m = pd.read_fwf(MORT / f"NHANES_{suffix}_MORT_2019_PUBLIC.dat",
                    colspecs=[(0, 6), (14, 15), (15, 16), (16, 19), (42, 45), (45, 48)],
                    names=["SEQN", "eligstat", "mortstat", "ucod", "permth_int", "permth_exm"], dtype=str)
    for c in m.columns:
        m[c] = pd.to_numeric(m[c].str.strip(), errors="coerce")
    return m


parts = []
for k, (cyc, demo, cbc, fast, tg, mcq, bmx, smq, crp, bio, mort) in enumerate(CYCLES):
    d = rd(demo, KEEP["demo"])
    for name, key in [(cbc, "cbc"), (fast, "fast"), (tg, "tg"), (mcq, "mcq"), (bmx, "bmx"),
                      (smq, "smq"), (crp, "crp"), (bio, "bio")]:
        if name:
            d = d.merge(rd(name, KEEP[key]), on="SEQN", how="left")
    if mort:
        d = d.merge(read_mort(mort), on="SEQN", how="left")
    d["cycle"] = cyc
    d["cycle_idx"] = k
    parts.append(d)
D = pd.concat(parts, ignore_index=True)

# SAS 传输格式里极小的数（如 5.4e-79）实际是 0
num = D.select_dtypes("number").columns
D[num] = D[num].where(~(D[num].abs() < 1e-10), 0.0)

o = pd.DataFrame({"SEQN": D.SEQN, "cycle": D.cycle, "cycle_idx": D.cycle_idx})
o["examined"] = D.RIDSTATR == 2
o["sex"] = D.RIAGENDR                         # 1 男, 2 女
o["age"] = D.RIDAGEYR                         # 80 岁封顶（部分周期 85）
o["race_eth"] = D.RIDRETH1                    # 1 墨西哥裔, 2 其他西语裔, 3 非西语裔白人, 4 非西语裔黑人, 5 其他
o["educ"] = D.DMDEDUC2.where(D.DMDEDUC2.isin([1, 2, 3, 4, 5]))
o["pir"] = D.INDFMPIR
o["pregnant"] = D.RIDEXPRG == 1
o["exam_season"] = D.RIDEXMON                 # 1 = 11–4 月, 2 = 5–10 月
o["psu"] = D.SDMVPSU
o["strata"] = D.cycle_idx * 1000 + D.SDMVSTRA  # 各周期层号加前缀，保证唯一
# 权重：1999–2018 合并 20 年按 NCHS 规则（1999–2002 用 4 年权重×4/20，其余 2 年权重×2/20）
early = D.cycle.isin(["1999-2000", "2001-2002"])
o["w_mec_pooled"] = np.where(D.cycle == "2021-2023", np.nan,
                             np.where(early, D.WTMEC4YR * 4 / 20, D.WTMEC2YR * 2 / 20))
o["w_mec_2y"] = D.WTMEC2YR
o["w_phleb_2123"] = D.WTPH2YR                 # 仅 2021–2023 的血常规用抽血权重

# 场次与禁食
o["session"] = np.where(D.cycle == "2021-2023", D.PHDSESNZ, D.PHDSESN)   # 0 上午；1 下午（2021–23 为下午/晚上）；2 晚上
o["fast_h"] = D.PHAFSTHR + D.PHAFSTMN / 60
o["fast_topcoded"] = (D.cycle == "2021-2023") & (D.PHAFSTHR == 30)        # 2021–23 "30 小时及以上"

# 血常规（五分类）
o["wbc"] = D.LBXWBCSI
o["lym"] = D.LBDLYMNO
o["neu"] = D.LBDNENO
o["mono"] = D.LBDMONO
o["eos"] = D.LBDEONO
o["baso"] = D.LBDBANO
o["plt"] = D.LBXPLTSI
ok_lym = o.lym > 0
o["nlr"] = (o.neu / o.lym).where(ok_lym)
o["sii"] = (o.plt * o.neu / o.lym).where(ok_lym)
o["siri"] = (o.neu * o.mono / o.lym).where(ok_lym)
o["piv"] = (o.neu * o.plt * o.mono / o.lym).where(ok_lym)
o["plr"] = (o.plt / o.lym).where(ok_lym)
o["glr"] = ((o.neu + o.eos + o.baso) / o.lym).where(ok_lym)   # 与 NHANES III 粒细胞口径对齐的比值

# 其他实验室
o["w_fast_sub"] = np.where(early, D.WTSAF4YR, D.WTSAF2YR)    # 上午空腹子样本权重
o["tg_ref"] = D.LBXTR
o["tg_biochem"] = D.LBXSTR
o["crp_1999_2010"] = D.LBXCRP                                # mg/dL
o["hscrp"] = D.LBXHSCRP                                      # mg/L（2015–2018、2021–2023）
o["hscrp_below_lod"] = D.LBDHRPLC == 1

# 协变量
o["smoke"] = np.select([D.SMQ020 == 2, (D.SMQ020 == 1) & D.SMQ040.isin([1, 2]), (D.SMQ020 == 1) & (D.SMQ040 == 3)],
                       ["never", "current", "former"], default="missing")
o["bmi"] = D.BMXBMI
o["height_cm"] = D.BMXHT
yn = lambda s: np.where(s == 1, 1.0, np.where(s == 2, 0.0, np.nan))
o["hx_cancer"] = yn(D.MCQ220)
o["hx_breast_cancer"] = D[["MCQ230A", "MCQ230B", "MCQ230C", "MCQ230D"]].eq(14).any(axis=1)
cvd = pd.DataFrame({c: yn(D[c]) for c in ["MCQ160B", "MCQ160C", "MCQ160E", "MCQ160F"]})
o["hx_cvd"] = np.where(cvd.eq(1).any(axis=1), 1.0, np.where(cvd.eq(0).all(axis=1), 0.0, np.nan))

# 死亡（2019 公共链接，2021–2023 无）
o["eligstat"] = D.eligstat
o["death_all"] = np.where(D.eligstat == 1, (D.mortstat == 1).astype(float), np.nan)
o["death_cancer"] = np.where(D.eligstat == 1, ((D.mortstat == 1) & (D.ucod == 2)).astype(float), np.nan)
o["death_accident"] = np.where(D.eligstat == 1, ((D.mortstat == 1) & (D.ucod == 4)).astype(float), np.nan)
o["fu_months_exam"] = D.permth_exm

# 纳入标记
o["ok_age20"] = o.age >= 20
o["ok_session"] = o.session.notna()
o["ok_cbc"] = o[["wbc", "lym", "neu", "mono", "plt"]].notna().all(axis=1) & ok_lym
o["ok_not_preg"] = ~o.pregnant
o["ok_fast"] = o.fast_h.notna() & (o.fast_h < 24)          # 预设：禁食 ≥24 小时视为不可信，剔除
o["in_main"] = o.examined & o.ok_age20 & o.ok_session & o.ok_cbc & o.ok_not_preg & o.ok_fast

DERIVED.mkdir(parents=True, exist_ok=True)
keep = o[o.examined & o.ok_age20].copy()
keep.to_pickle(DERIVED / "continuous_analytic.pkl")
m = keep[keep.in_main]
print("≥20 岁受检者", len(keep), "| 主样本", len(m),
      "| 1999–2018", int((m.cycle != "2021-2023").sum()), "| 2021–2023", int((m.cycle == "2021-2023").sum()))
