"""
构建 NHANES III（1988–1994）分析数据集。

做什么：
  1. 从实验室、体检、成人问卷三个定长文件读出所需变量，按 SEQN 合并；
  2. 用 WTPFSD6 / WTPFMD6 两个权重"是否空白"还原每个人被随机分到的场次
     （官方文档：家户随机分到上午 standard 或下午/晚上 modified 场次；
       分到该场次但去了别的场次的人，该场次权重记为 0；另一个场次权重为空白）；
  3. 把 8888 类"应答但空白"编码转成缺失，派生白细胞比值指标；
  4. 合并 2019 版公共死亡链接。

不做什么：不计算任何"场次/时相 → 指标"或"指标 → 死亡"的关联（方案锁定前保持盲态）。

输出：analysis/derived/nhanes3_analytic.pkl（全部 ≥20 岁 MEC 受检者，含纳入标记）
"""
import re
import numpy as np
import pandas as pd
from config import N3, DERIVED
from n3io import read_raw, to_num

# ---------- 1. 读原始变量 ----------
LAB_VARS = ["SEQN", "DMPFSEQ", "DMPSTAT", "HSSEX", "HSAGEIR", "HSAGEU", "DMARETHN", "DMPPIR",
            "SDPPHASE", "SDPPSU6", "SDPSTRA6", "WTPFEX6", "WTPFSD6", "WTPFMD6",
            "MXPSESSR", "MXPTIDW", "PHPCHM2", "PHPFAST", "PHPBEST",
            "WCP", "LMP", "MOP", "GRP", "PLP", "HGP",
            "GRPDIF", "LMPDIF", "MOPDIF", "EOP", "BOP", "BAP", "LAP",
            "TGP", "CRP", "FBP", "FEP", "SGP", "G1P"]
ADULT_VARS = ["SEQN", "HFA8R", "HAB1", "HAD1", "HAE2", "HAF10", "HAC1C", "HAC1D",
              "HAC1N", "HAC1O", "HAC3OS", "HAR1", "HAR3"]
EXAM_VARS = ["SEQN", "BMPBMI", "BMPHT", "PEP13A", "PEP13C"]

REP_VARS = [f"WTPXRP{i}" for i in range(1, 53)]     # Fay BRR 重复权重（MEC 受检样本，官方提供，rho = 0.3）
lab = read_raw("lab.dat", "lab.sas", LAB_VARS + REP_VARS)
adult = read_raw("adult.dat", "adult.sas", ADULT_VARS)
exam = read_raw("exam.dat", "exam.sas", EXAM_VARS)

d = lab.merge(adult, on="SEQN", how="left").merge(exam, on="SEQN", how="left")
assert d.SEQN.is_unique

# ---------- 2. 还原随机分配的场次（看原始字符串，区分"空白"和"0"） ----------
sd_blank = d.WTPFSD6 == ""
md_blank = d.WTPFMD6 == ""
assert not ((~sd_blank) & (~md_blank)).any(), "有人两个场次权重都非空白，与官方文档不符"
d["arm"] = np.select([~sd_blank, ~md_blank], ["AM", "PMEVE"], default="none")

# ---------- 3. 数值化并处理缺失编码 ----------
num = lambda v, codes=(): to_num(d[v], codes)
out = pd.DataFrame({"SEQN": num("SEQN")})
out["family_seq"] = num("DMPFSEQ")
out["exam_status"] = num("DMPSTAT")                # 2 = MEC 受检, 3 = 家访体检
out["sex"] = num("HSSEX")                           # 1 男, 2 女
out["age"] = num("HSAGEIR")                         # 本文件年龄单位全为"年"（HSAGEU=2），90 岁封顶
assert (d.HSAGEU == "2").all()
out["race_eth"] = num("DMARETHN")                   # 1 非西语裔白人, 2 非西语裔黑人, 3 墨西哥裔, 4 其他
out["pir"] = num("DMPPIR", ("888888",))
out["phase"] = num("SDPPHASE")
out["psu"] = num("SDPPSU6")
out["strata"] = num("SDPSTRA6")
out["w_mec"] = num("WTPFEX6")
out["w_am"] = num("WTPFSD6")                        # 空白 -> NaN；0 保留为 0
out["w_pm"] = num("WTPFMD6")
out["arm"] = d["arm"].values                       # 随机分配：AM 上午 / PMEVE 下午或晚上 / none 无分配记录
out["session"] = num("MXPSESSR")                    # 实际场次 1 上午, 2 下午, 3 晚上
out["dow"] = num("MXPTIDW")
out["chemo4wk"] = num("PHPCHM2", ("8",))            # 1 近 4 周化疗（这些人不抽血）


def clock_to_hours(s):
    s = s.strip()
    m = re.fullmatch(r"(\d\d):(\d\d)", s)
    return int(m.group(1)) + int(m.group(2)) / 60 if m else np.nan


out["draw_clock_h"] = d["PHPBEST"].map(clock_to_hours)          # 抽血钟点（小时，07:32–22:02）
out["fast_h"] = num("PHPFAST", ("88888",))                       # 计算的禁食时长（小时）

# 自动血细胞计数（Coulter S-Plus JR，三分类）
out["wbc"] = num("WCP", ("88888",))       # 10^3/uL
out["lym"] = num("LMP", ("88888",))
out["mono_nuc"] = num("MOP", ("8888",))   # 单个核细胞（Coulter 分群），不等同于单核细胞
out["gran"] = num("GRP", ("88888",))      # 粒细胞（含中性、嗜酸、嗜碱）
out["plt"] = num("PLP", ("88888",))       # 10^3/uL
out["hgb"] = num("HGP", ("88888",))
# 人工分类（子样本，100 个细胞百分比）
for v, name in [("GRPDIF", "man_seg_pct"), ("LMPDIF", "man_lym_pct"), ("MOPDIF", "man_mono_pct"),
                ("EOP", "man_eos_pct"), ("BOP", "man_baso_pct"), ("BAP", "man_band_pct"),
                ("LAP", "man_atyp_pct")]:
    out[name] = num(v, ("888", "88"))
# 其他实验室指标
out["tg"] = num("TGP", ("8888",))
out["crp"] = num("CRP", ("88888",))       # mg/dL，检测下限 0.21
out["fibrinogen"] = num("FBP", ("8888",))
out["iron"] = num("FEP", ("888",))
out["glu_serum"] = num("SGP", ("888",))
out["glu_plasma"] = num("G1P", ("88888",))

# 协变量（成人问卷，入户访谈时收集，早于体检与场次）
out["educ_years"] = num("HFA8R", ("88", "99"))
out["self_health"] = num("HAB1", ("8", "9"))          # 1 极好 ... 5 差
for v, name in [("HAD1", "hx_diabetes"), ("HAE2", "hx_hypertension"), ("HAF10", "hx_mi"),
                ("HAC1C", "hx_chf"), ("HAC1D", "hx_stroke"), ("HAC1N", "hx_skin_cancer"),
                ("HAC1O", "hx_other_cancer")]:
    x = num(v, ("8", "9"))
    out[name] = np.where(x == 1, 1.0, np.where(x == 2, 0.0, np.nan))
out["cancer_site_code"] = num("HAC3OS", ("88", "99"))  # 癌症部位编码（未核对码本；本研究分析不使用）
har1 = num("HAR1", ("8",))
har3 = num("HAR3", ("8",))
out["smoke"] = np.select([har1 == 2, (har1 == 1) & (har3 == 1), (har1 == 1) & (har3 == 2)],
                         ["never", "current", "former"], default="missing")
out["bmi"] = num("BMPBMI", ("8888",))
out["height_cm"] = num("BMPHT", ("88888",))
out["pe_health"] = num("PEP13A", ("8",))              # 体检医生评估健康状况：1 极好 … 5 差（exam-acc 码本位置 1481）
out["pe_infection"] = num("PEP13C", ("8",))           # 体检医生记录：1 无感染, 2 可能有感染（exam-acc 码本位置 1483）

# ---------- 4. 派生指标 ----------
# 注意：NHANES III 只有三分类，"NLR"只能用粒细胞/淋巴细胞比（GLR）近似
out["glr"] = out.gran / out.lym
out["sii_g"] = out.plt * out.gran / out.lym             # SII 的粒细胞版近似（×10^3/uL 单位下的数值）
out["plr"] = out.plt / out.lym
# 人工分类得到的真正 NLR（分叶核 + 杆状核）/ 淋巴细胞
out["nlr_manual"] = (out.man_seg_pct + out.man_band_pct.fillna(0)) / out.man_lym_pct
out.loc[out.man_lym_pct <= 0, "nlr_manual"] = np.nan

# 场次合规：实际场次与分配是否一致
out["session_pmeve"] = np.where(out.session.isin([2, 3]), 1.0, np.where(out.session == 1, 0.0, np.nan))
out["arm_pmeve"] = np.where(out["arm"] == "PMEVE", 1.0, np.where(out["arm"] == "AM", 0.0, np.nan))
out["complier"] = np.where(out["arm"] == "none", np.nan,
                           (out.arm_pmeve == out.session_pmeve).astype(float))

# ---------- 5. 死亡链接（2019 公共版） ----------
m = pd.read_fwf(N3 / "NHANES_III_MORT_2019_PUBLIC.dat",
                colspecs=[(0, 6), (14, 15), (15, 16), (16, 19), (42, 45), (45, 48)],
                names=["SEQN", "eligstat", "mortstat", "ucod", "permth_int", "permth_exm"], dtype=str)
for c in m.columns:
    m[c] = pd.to_numeric(m[c].str.strip(), errors="coerce")
assert m.SEQN.is_unique
out = out.merge(m, on="SEQN", how="left")
# 死因（UCOD_LEADING）：1 心脏病, 2 恶性肿瘤, 4 意外伤害, 其余见官方说明
out["death_all"] = np.where(out.eligstat == 1, (out.mortstat == 1).astype(float), np.nan)
out["death_cancer"] = np.where(out.eligstat == 1, ((out.mortstat == 1) & (out.ucod == 2)).astype(float), np.nan)
out["death_accident"] = np.where(out.eligstat == 1, ((out.mortstat == 1) & (out.ucod == 4)).astype(float), np.nan)
out["fu_months_exam"] = out.permth_exm

# ---------- 6. 纳入标记（逐层，便于画流程图） ----------
# ITT 样本只用分配前就确定的条件 + 结局本身可测；钟点、禁食、BMI 都是分配后才测的，不作为 ITT 纳入条件
f = pd.DataFrame(index=out.index)
f["age20"] = out.age >= 20
f["mec"] = out.session.notna()                           # MEC 受检（排除家访体检）
f["assigned"] = out["arm"].isin(["AM", "PMEVE"])         # 有随机分配记录（≥12 岁 MEC 受检者均应有）
f["not_chemo"] = out.chemo4wk != 1
f["cbc_ok"] = out[["wbc", "lym", "gran", "plt"]].notna().all(axis=1) & (out.lym > 0)
f["timing_ok"] = out.draw_clock_h.notna() & out.fast_h.notna()
for c in f.columns:
    out["ok_" + c] = f[c].values
out["in_itt"] = f[["age20", "mec", "assigned", "not_chemo", "cbc_ok"]].all(axis=1).values
out["in_timing"] = out["in_itt"] & f["timing_ok"].values          # 钟点/禁食分析（Part 1b）

rep = pd.DataFrame({v.lower(): to_num(d[v]).values for v in REP_VARS}, index=d.index)
rep["SEQN"] = to_num(d["SEQN"]).values
out = out.merge(rep, on="SEQN", how="left")

DERIVED.mkdir(parents=True, exist_ok=True)
keep = out[out.age >= 20].copy()          # 只存成人，体积小
keep.to_pickle(DERIVED / "nhanes3_analytic.pkl")
print(f"NHANES III 成人记录 {len(keep)}，ITT 样本 {int(keep.in_itt.sum())}，钟点/禁食样本 {int(keep.in_timing.sum())}")
