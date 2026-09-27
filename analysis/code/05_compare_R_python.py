"""
R 主分析与 Python 独立实现的逐项核对（方案 5.6）：点估计差 <0.00005（小数点后 4 位一致）、标准误相对差 <5%。
用法：python3 05_compare_R_python.py [scrambled|real]
"""
import json
import sys
from config import PROJECT

mode = (sys.argv[1] if len(sys.argv) > 1 else "scrambled").lower()
tag = mode.upper()
R = json.load(open(PROJECT / "analysis/results" / mode / f"results_{tag}.json", encoding="utf-8"))
P = json.load(open(PROJECT / "analysis/py_check" / f"key_estimates_{tag}.json", encoding="utf-8"))

p1, p2 = R["part1"], R.get("part2", {}).get("nhanes3", {})
rows = [
    ("ITT β（ln GLR）", p1["primary_glr_taylor"]["beta"], P["beta_itt"], "est"),
    ("ITT SE（Taylor）", p1["primary_glr_taylor"]["se"], P["se_taylor"], "se"),
    ("ITT SE（BRR）", p1["primary_glr_brr"]["se"], P["se_brr"], "se"),
    ("CACE", p1["cace"]["beta"], P["cace"], "est"),
    ("CACE SE（BRR）", p1["cace"]["se"], P["cace_se_brr"], "se"),
]
for ev, key in [("death_all", "all"), ("death_cancer", "cancer")]:
    if ev in p2:
        rows += [(f"{key}: β_obs", p2[ev]["beta_obs"], P[key]["beta_obs"], "est"),
                 (f"{key}: β_std", p2[ev]["beta_std"], P[key]["beta_std"], "est"),
                 (f"{key}: Δ", p2[ev]["delta"], P[key]["delta"], "est"),
                 (f"{key}: SE(Δ)（BRR）", p2[ev]["se"], P[key]["delta_se_brr"], "se")]
if "analytic" in p2:
    rows.append(("预测 Δ（解析）", p2["analytic"]["predicted_delta_pct"] / 100, P["all"]["predicted_delta"], "est"))

ok_all = True
missing = [ev for ev in ("death_all", "death_cancer") if ev not in p2] + ([] if "analytic" in p2 else ["analytic"])
if missing:
    ok_all = False
    print(f"注意：R 结果缺少 Part 2 项 {missing}，核对不完整")
print(f"{'项目':<22}{'R':>22}{'Python':>22}{'差':>14}  结果")
for name, r, p, kind in rows:
    if kind == "est":
        diff = abs(r - p); ok = diff < 5e-5; shown = f"{diff:.2e}"
    else:
        diff = abs(r - p) / abs(p); ok = diff < 0.05; shown = f"{100 * diff:.2f}%"
    ok_all &= ok
    print(f"{name:<22}{r:>22.10f}{p:>22.10f}{shown:>14}  {'通过' if ok else '不通过'}")
print("\n总判定：", "全部通过，可以进入解读" if ok_all else "有不通过项，先查原因，不进入解读")
(PROJECT / "analysis/results" / mode / f"compare_R_python_{tag}.txt").write_text(
    "\n".join(f"{n}\tR={r!r}\tPython={p!r}" for n, r, p, _ in rows) + f"\n总判定\t{'PASS' if ok_all else 'FAIL'}\n",
    encoding="utf-8")
