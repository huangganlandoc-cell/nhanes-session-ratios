"""
把 R 结果（results_<MODE>.json）写成中文摘要，并按方案 v1.0 第 7 节把 R1 × R2 组合对应到预先写好的结论措辞。
用法：python3 06_write_summary.py [scrambled|real]
"""
import json
import sys
from config import PROJECT

mode = (sys.argv[1] if len(sys.argv) > 1 else "scrambled").lower()
tag = mode.upper()
R = json.load(open(PROJECT / "analysis/results" / mode / f"results_{tag}.json", encoding="utf-8"))
L = []
P = L.append
if tag == "SCRAMBLED":
    P("> **SCRAMBLED – NOT REAL**：以下数字来自置乱数据，只用于检查代码，没有任何实际含义。\n")

f = lambda x, d=2: f"{x:+.{d}f}"
ci = lambda r: f"{f(r['pct'])}%（95% CI {f(r['pct_lo'])} 至 {f(r['pct_hi'])}）"

P(f"# 结果摘要（{tag}）\n")
P(f"运行时间：{R['meta']['run_time']}；方案：{R['meta']['protocol']}；{R['meta']['R']}，survey {R['meta']['survey']}\n")

p1 = R["part1"]
pr = p1["primary_glr_taylor"]
P("## Part 1：场次随机分配对 GLR 的效应（NHANES III，ITT）\n")
P(f"- 样本：{p1['n_itt']} 人（下午/晚上组 {p1['n_z1']}，上午组 {p1['n_z0']}）")
P(f"- **主估计**：下午/晚上 vs 上午 {ci(pr)}；BRR 置信区间 {f(p1['primary_glr_brr']['pct_lo'])} 至 {f(p1['primary_glr_brr']['pct_hi'])}")
P(f"- **R1（界值 9.5%）：{pr['R1']}**；最严界值 6.3%：{pr['R1_strict_6.3']}；点估计所在档：{pr['tier']}")
P(f"- CACE：{ci(p1['cace'])}（依从：P(D=1|Z=1) = {p1['cace']['p1']:.3f}，P(D=1|Z=0) = {p1['cace']['p0']:.3f}）")
P(f"- 标准化位移 {p1['sd_shift']['est']:+.3f} SD；偏 R² {100 * p1['partial_r2']['est']:.3f}%（{100 * p1['partial_r2']['lo']:.3f}–{100 * p1['partial_r2']['hi']:.3f}%）")
P("- 重分类（比例差，百分点）：")
for k, v in p1["reclassification"].items():
    P(f"  - {k}：下午/晚上 {100 * v['p1']:.1f}% vs 上午 {100 * v['p0']:.1f}%，差 {100 * v['rd']:+.1f}（{100 * v['rd_lo']:+.1f} 至 {100 * v['rd_hi']:+.1f}）；"
      f"归因比例 {100 * v['af']:+.1f}%")
P("\n| 敏感性分析 | 效应 | R1 |\n|---|---|---|")
for k, v in p1["sensitivity"].items():
    P(f"| {k} | {ci(v)} | {v['R1']} |")
P("\n| 次要指标 | 效应 | 界值 | R1 |\n|---|---|---|---|")
for k, v in p1["secondary"].items():
    if "pct" in v:
        P(f"| {k} | {ci(v)} | {v['threshold'] if v['threshold'] is not None else '—'} | {v['R1'] if v['R1'] is not None else '—'} |")
c = p1["secondary"]["crp_detectable_rd"]
P(f"\nCRP 可检出比例差：{100 * c['rd']:+.1f} 个百分点（{100 * c['lo']:+.1f} 至 {100 * c['hi']:+.1f}）\n")

if "part1c" in R:
    q = R["part1c"]
    a, e = q["nlr_primary"]["afternoon"], q["nlr_primary"]["evening"]
    P("## Part 1c：连续 NHANES 1999–2018 复现（非随机）\n")
    P(f"- 样本 {q['n_main']}；NLR 下午 vs 上午 {ci(a)}；晚上 vs 上午 {ci(e)}")
    P(f"- **R3：{q['R3']}**")
    P(f"- 2021–2023（n = {q['c2021_2023']['n']}）：下午/晚上 vs 上午 {ci(q['c2021_2023']['pmeve_vs_morning'])}")
    for k, v in q["consequences"].items():
        if isinstance(v, dict):
            P(f"- {k}：上午 {100 * v['p_morning']:.1f}%，下午 {100 * v['p_afternoon']:.1f}%（差 {100 * v['rd_afternoon']:+.1f}，"
              f"{100 * v['rd_afternoon_lo']:+.1f} 至 {100 * v['rd_afternoon_hi']:+.1f}），晚上 {100 * v['p_evening']:.1f}%")
    P("")

if "part2" in R:
    q = R["part2"]
    P("## Part 2：场次标准化对死亡 HR 的影响\n")
    P("| 数据 | 结局 | 事件 | HR/翻倍（观察） | HR/翻倍（标准化） | Δ（ln HR 相对变化） | R2 |\n|---|---|---|---|---|---|---|")
    for s in ["nhanes3", "continuous"]:
        for ev in ["death_all", "death_cancer"]:
            r = q[s][ev]
            P(f"| {s} | {ev} | {r['events']} | {r['hr_obs_per_doubling']:.3f} | {r['hr_std_per_doubling']:.3f} | "
              f"{f(r['delta_pct'])}%（{f(r['lo_pct'])} 至 {f(r['hi_pct'])}） | {r['R2']} |")
    P(f"\n解析预测的 Δ：NHANES III {q['nhanes3']['analytic']['predicted_delta_pct']:+.2f}%；"
      f"连续 NHANES {q['continuous']['analytic']['predicted_delta_pct']:+.2f}%\n")

    # 方案第 7 节：按 R1 × R2（NHANES III；全因与癌症分别）对应预设结论
    concl = {
        ("exceeds desirable bias", "robust"): "采血场次使 GLR/NLR 发生有临床意义的偏移，并使截断值附近的人被重分类；但对死亡 HR 的影响可忽略。风险在于截断值与参考区间的可迁移性，不在关联估计本身。",
        ("exceeds desirable bias", "materially affected"): "场次既改变指标数值，也实质改变 HR；使用 NHANES 白细胞比值的分析需要处理场次。",
        ("within desirable bias", "robust"): "在随机化证据下，NHANES 采血场次对 GLR/NLR 没有临床意义的影响——对现有文献是有检出力的阴性结论。",
        ("within desirable bias", "materially affected"): "按衰减公式不应出现；视为分析错误信号，先查代码与数据，不据此下结论。",
    }
    P("## 按方案第 7 节对应的结论措辞\n")
    for ev in ["death_all", "death_cancer"]:
        r1, r2 = pr["R1"], q["nhanes3"][ev]["R2"]
        txt = concl.get((r1, r2)) or ("报告置信区间，写\"不能排除有临床意义的偏移\"，不写\"无影响\"。" if r1 == "indeterminate"
                                     else "R2 不确定：报告置信区间，不下\"稳健\"或\"受影响\"的结论。")
        P(f"- {ev}：R1 = {r1}，R2 = {r2} → {txt}")

out = PROJECT / "analysis/results" / mode / f"summary_{tag}.md"
out.write_text("\n".join(L) + "\n", encoding="utf-8")
print(out)
