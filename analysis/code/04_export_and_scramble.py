"""
导出给 R 用的 CSV，并生成"置乱数据"（方案 5.6 盲态开发）。

置乱规则（方案 v1.0 第 5.6 节）：
  - NHANES III：在设计层（strata）内随机置换"分组标签"。为保留依从结构，分组、实际场次、场次权重、
    抽血钟点、禁食时长作为一整块一起置换，置换范围是"≥20 岁、MEC 受检、有分配、非近 4 周化疗"者；
  - 连续 NHANES：在周期 × 层内置换场次（与禁食时长一起）；
  - 死亡结局（状态、死因、随访时间）在人群间整体置换。
置乱后，场次与指标、指标与死亡之间的真实关联全部被打断，只能用于调试代码。

输出：
  derived/nhanes3_REAL.csv.gz、derived/continuous_REAL.csv.gz（真实数据，登记后才允许分析）
  derived/nhanes3_SCRAMBLED.csv.gz、derived/continuous_SCRAMBLED.csv.gz
"""
import numpy as np
import pandas as pd
from config import DERIVED, SEED

rng = np.random.default_rng(SEED)


def permute_block(df, cols, group_cols, mask):
    """在 mask 选中的行里，按 group_cols 分组，把 cols 这一整块行一起随机重排。"""
    out = df.copy()
    idx = df.index[mask]
    for _, g in df.loc[idx].groupby(group_cols, dropna=False):
        src = g.index.to_numpy()
        dst = rng.permutation(src)
        out.loc[src, cols] = df.loc[dst, cols].to_numpy()
    return out


# ---------------- NHANES III ----------------
n3 = pd.read_pickle(DERIVED / "nhanes3_analytic.pkl")
n3.to_csv(DERIVED / "nhanes3_REAL.csv.gz", index=False)

expo3 = ["arm", "arm_pmeve", "w_am", "w_pm", "session", "session_pmeve", "complier", "dow",
         "draw_clock_h", "fast_h", "ok_timing_ok"]
base3 = n3.ok_age20 & n3.ok_mec & n3.ok_assigned & n3.ok_not_chemo
s3 = permute_block(n3, expo3, ["strata"], base3)
mort3 = ["eligstat", "mortstat", "ucod", "permth_int", "permth_exm", "death_all", "death_cancer",
         "death_accident", "fu_months_exam"]
s3 = permute_block(s3, mort3, [pd.Series(0, index=s3.index)], s3.in_itt)
s3["in_timing"] = s3.in_itt & s3.ok_timing_ok.astype(bool)
s3["SCRAMBLED"] = 1
s3.to_csv(DERIVED / "nhanes3_SCRAMBLED.csv.gz", index=False)

# ---------------- 连续 NHANES ----------------
c = pd.read_pickle(DERIVED / "continuous_analytic.pkl")
c.to_csv(DERIVED / "continuous_REAL.csv.gz", index=False)

expoc = ["session", "fast_h", "fast_topcoded", "ok_session", "ok_fast"]
sc = permute_block(c, expoc, ["cycle", "strata"], c.examined & c.ok_age20)
mortc = ["eligstat", "death_all", "death_cancer", "death_accident", "fu_months_exam"]
has_mort = c.cycle != "2021-2023"
sc = permute_block(sc, mortc, [pd.Series(0, index=sc.index)], has_mort & sc.ok_age20)
sc["in_main"] = sc.examined & sc.ok_age20 & sc.ok_session.astype(bool) & sc.ok_cbc & sc.ok_not_preg & sc.ok_fast.astype(bool)
sc["SCRAMBLED"] = 1
sc.to_csv(DERIVED / "continuous_SCRAMBLED.csv.gz", index=False)

# 自检：置乱后分组比例、依从率、样本量应与真实数据基本一致
itt, itt_s = n3[n3.in_itt], s3[s3.in_itt]
print("NHANES III ITT", len(itt), "| 置乱后上午组比例", round((itt_s["arm"] == "AM").mean(), 4),
      "（真实", round((itt["arm"] == "AM").mean(), 4), "）| 置乱后依从率", round(itt_s.complier.mean(), 4),
      "（真实", round(itt.complier.mean(), 4), "）")
print("连续 NHANES 主样本：真实", int(c.in_main.sum()), "置乱", int(sc.in_main.sum()))
print("已写出 REAL 与 SCRAMBLED 两套 CSV")
