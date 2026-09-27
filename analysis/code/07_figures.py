"""
发表级图（Stage 14d）。所有数字直接读自 analysis/results/real/results_REAL.json（第 2 次运行）与分析数据集，不手抄。
输出：analysis/figures/*.pdf（矢量）与 *.png（300 dpi）
配色：dataviz 参考调色板前 3 槽（蓝 #2a78d6、橙 #eb6834、青 #1baf7a，已用 validate_palette.js 验证）；
     形状作第二编码（圆 = 计数/观察值，方 = 比值/标准化值），黑白打印也可区分。
"""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from config import PROJECT, DERIVED

R = json.load(open(PROJECT / "analysis/results/real/results_REAL.json", encoding="utf-8"))
FIG = PROJECT / "analysis/figures"
FIG.mkdir(parents=True, exist_ok=True)

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID, BAND = "#0b0b0b", "#52514e", "#e6e5e1", "#f0efec"
plt.rcParams.update({
    "font.family": "Arial", "font.size": 7.5, "axes.titlesize": 8, "axes.labelsize": 7.5,
    "xtick.labelsize": 7, "ytick.labelsize": 7.5, "legend.fontsize": 7, "axes.edgecolor": INK2,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0, "xtick.color": INK2,
    "ytick.color": INK, "text.color": INK, "axes.labelcolor": INK, "pdf.fonttype": 42, "svg.fonttype": "none",
})
MM = 1 / 25.4


def save(fig, name):
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("写出", name)


def clean_axes(ax):
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.grid(axis="x", color=GRID, linewidth=0.5)
    ax.set_axisbelow(True)


def forest(ax, rows, xlim, xlabel, band_key="threshold"):
    """rows: list of dict(label, est, lo, hi, color, marker, band) ；自上而下绘制；band 为 ±期望偏倚阴影。"""
    n = len(rows)
    for i, r in enumerate(rows):
        y = n - 1 - i
        if r.get(band_key):
            ax.add_patch(plt.Rectangle((-r[band_key], y - 0.38), 2 * r[band_key], 0.76, color=BAND, lw=0, zorder=0))
        if r.get("est") is None:
            continue
        ax.plot([r["lo"], r["hi"]], [y, y], color=r["color"], lw=1.2, solid_capstyle="round", zorder=2)
        ax.plot(r["est"], y, marker=r["marker"], ms=5, color=r["color"], mec="white", mew=0.8, zorder=3)
    ax.axvline(0, color=INK2, lw=0.6, zorder=1)
    ax.set_yticks(range(n))
    ax.set_yticklabels([r["label"] for r in rows][::-1])
    ax.set_xlim(*xlim)
    ax.set_ylim(-0.6, n - 0.4)
    ax.set_xlabel(xlabel)
    clean_axes(ax)


# ============ 图 1：研究流程 ============
n3 = pd.read_pickle(DERIVED / "nhanes3_analytic.pkl")
c = pd.read_pickle(DERIVED / "continuous_analytic.pkl")
itt = n3[n3.in_itt]
am, pm = itt[itt["arm"] == "AM"], itt[itt["arm"] == "PMEVE"]
am_ok, pm_ok = int((am.session == 1).sum()), int(pm.session.isin([2, 3]).sum())
mort3 = itt[(itt.eligstat == 1) & itt.fu_months_exam.notna()]
m18 = c[c.in_main & (c.cycle != "2021-2023")]
m23 = c[c.in_main & (c.cycle == "2021-2023")]
mort18 = m18[(m18.eligstat == 1) & m18.fu_months_exam.notna()]
base_mec = n3[n3.ok_mec]
excl_chemo = int((base_mec.ok_assigned & ~base_mec.ok_not_chemo).sum())
excl_cbc = int((base_mec.ok_assigned & base_mec.ok_not_chemo & ~base_mec.ok_cbc_ok).sum())

fig, ax = plt.subplots(figsize=(178 * MM, 90 * MM))
ax.set_xlim(0, 178); ax.set_ylim(0, 90); ax.axis("off")
FS = 6.2


def box(x, y, w, h, text, fc="white", ec=INK2, ha="center"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.3,rounding_size=1.2", fc=fc, ec=ec, lw=0.6))
    tx = x + w / 2 if ha == "center" else x + 2.5
    ax.text(tx, y + h / 2, text, ha=ha, va="center", fontsize=FS, linespacing=1.35)


def arrow(x1, y1, x2, y2):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle="-|>", lw=0.6, color=INK2, mutation_scale=6))


n_excl = len(n3) - len(itt)
ax.text(2, 86.5, "A  NHANES III (1988–1994): randomized sessions", fontsize=7.2, fontweight="bold")
box(4, 70, 52, 11, f"Adults aged ≥20 y in the NHANES III\nlaboratory file (n = {len(n3):,})")
box(60, 56, 41, 19, f"Excluded (n = {n_excl:,})\n• Home examination: {len(n3) - int(n3.ok_mec.sum()):,}\n"
                   f"• Chemotherapy <4 wk: {excl_chemo:,}\n• Missing blood count: {excl_cbc:,}", ec=GRID, ha="left")
arrow(30, 70, 30, 59); arrow(30, 65, 60, 65)
box(4, 48, 52, 11, f"Intention-to-treat sample (n = {len(itt):,})\nsession randomized by household")
arrow(18, 48, 18, 42); arrow(42, 48, 42, 42)
box(3, 23, 30, 19, f"Morning\nsession\nn = {len(am):,}\nExamined as assigned\n{am_ok:,} ({100 * am_ok / len(am):.1f}%)",
    fc="#eef4fc", ec=BLUE)
box(35, 23, 30, 19, f"Afternoon/evening\nsession\nn = {len(pm):,}\nExamined as assigned\n{pm_ok:,} ({100 * pm_ok / len(pm):.1f}%)",
    fc="#fdf0ea", ec=ORANGE)
arrow(18, 23, 18, 17); arrow(50, 23, 50, 17)
box(3, 2, 62, 15, f"Linked mortality through 2019 (n = {len(mort3):,})\nAll-cause deaths {int(mort3.death_all.sum()):,}; "
                  f"cancer deaths {int(mort3.death_cancer.sum()):,}\nMedian follow-up {mort3.fu_months_exam.median() / 12:.1f} years")

ax.text(104, 86.5, "B  Continuous NHANES: actual session", fontsize=7.2, fontweight="bold")
box(104, 57, 72, 24, f"NHANES 1999–2018 (n = {len(m18):,}; non-randomized)\nAdults ≥20 y with session, complete differential,\n"
                     f"not pregnant, fasting <24 h\nMorning {int((m18.session == 0).sum()):,} · Afternoon {int((m18.session == 1).sum()):,} · "
                     f"Evening {int((m18.session == 2).sum()):,}")
arrow(140, 57, 140, 50)
box(104, 35, 72, 15, f"Linked mortality through 2019 (n = {len(mort18):,})\nAll-cause deaths {int(mort18.death_all.sum()):,}; "
                     f"cancer deaths {int(mort18.death_cancer.sum()):,}")
box(104, 2, 72, 15, f"NHANES 2021–2023, cross-sectional (n = {len(m23):,})\nMorning {int((m23.session == 0).sum()):,} · "
                    f"Afternoon/evening {int((m23.session == 1).sum()):,}")
save(fig, "Figure1_flow")

# ============ 图 2：随机分配的场次效应（NHANES III ITT） ============
p1 = R["part1"]; s = p1["secondary"]; pr = p1["primary_glr_taylor"]


def row(label, d, color, marker, band):
    return dict(label=label, est=d["pct"], lo=d["pct_lo"], hi=d["pct_hi"], color=color, marker=marker, threshold=band)


rows_counts = [row("Leukocytes", s["wbc"], BLUE, "o", s["wbc"]["threshold"]),
               row("Lymphocytes", s["lym"], BLUE, "o", s["lym"]["threshold"]),
               row("Granulocytes", s["gran"], BLUE, "o", s["gran"]["threshold"]),
               row("Mononuclear cells", s["mono_nuc"], BLUE, "o", s["mono_nuc"]["threshold"]),
               row("Platelets", s["plt"], BLUE, "o", s["plt"]["threshold"])]
rows_ratios = [row("GLR (primary)", pr, ORANGE, "s", 9.5),
               row("SII (granulocyte-based)", s["sii_g"], ORANGE, "s", s["sii_g"]["threshold"]),
               row("PLR", s["plr"], ORANGE, "s", None)]
fig, axes = plt.subplots(2, 1, figsize=(110 * MM, 88 * MM), gridspec_kw=dict(height_ratios=[5, 3], hspace=0.55), sharex=True)
forest(axes[0], rows_counts, (-16, 18), "")
forest(axes[1], rows_ratios, (-16, 18), "Difference, afternoon/evening vs morning assignment (%)")
axes[0].set_title("Cell counts", loc="left", fontweight="bold", color=BLUE)
axes[1].set_title("Leukocyte ratio indices", loc="left", fontweight="bold", color=ORANGE)
for ax, rr in ((axes[0], rows_counts), (axes[1], rows_ratios)):
    n = len(rr)
    for i, r in enumerate(rr):
        ax.text(1.03, n - 1 - i, f"{r['est']:+.1f} ({r['lo']:+.1f} to {r['hi']:+.1f})", transform=ax.get_yaxis_transform(),
                ha="left", va="center", fontsize=6.5, color=INK2)
        ax.text(1.03 + 0.34, n - 1 - i, "–" if not r.get("threshold") else f"±{r['threshold']:.1f}",
                transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=6.5, color=INK2)
axes[0].text(1.03, len(rows_counts) - 0.3, "% (95% CI)", transform=axes[0].get_yaxis_transform(), ha="left", va="bottom",
             fontsize=6.5, fontweight="bold")
axes[0].text(1.37, len(rows_counts) - 0.3, "DB", transform=axes[0].get_yaxis_transform(), ha="left", va="bottom",
             fontsize=6.5, fontweight="bold")
save(fig, "Figure2_randomized_effects")

# ============ 图 3：连续 NHANES 复现 ============
q = R["part1c"]; qs = q["secondary"]
labels = [("Leukocytes", "wbc"), ("Lymphocytes", "lym"), ("Neutrophils", "neu"), ("Monocytes", "mono"), ("Platelets", "plt"),
          ("NLR", None), ("GLR-equivalent", "glr"), ("SII", "sii"), ("SIRI", "siri"), ("PIV", "piv"), ("PLR", "plr")]
fig, (ax, bx) = plt.subplots(1, 2, figsize=(178 * MM, 84 * MM), gridspec_kw=dict(width_ratios=[1.35, 1], wspace=0.35))
n = len(labels)
for i, (lab, key) in enumerate(labels):
    d = q["nlr_primary"] if key is None else qs[key]
    y = n - 1 - i
    band = {"NLR": 9.5, "GLR-equivalent": 9.5, "SII": 10.8, "Leukocytes": 5.1, "Lymphocytes": 6.3, "Neutrophils": 7.2,
            "Monocytes": 6.7, "Platelets": 5.0}.get(lab)
    if band:
        ax.add_patch(plt.Rectangle((-band, y - 0.38), 2 * band, 0.76, color=BAND, lw=0, zorder=0))
    for off, sess, col, mk in ((0.14, "afternoon", BLUE, "o"), (-0.14, "evening", ORANGE, "s")):
        e = d[sess]
        ax.plot([e["pct_lo"], e["pct_hi"]], [y + off, y + off], color=col, lw=1.1, zorder=2)
        ax.plot(e["pct"], y + off, marker=mk, ms=4.5, color=col, mec="white", mew=0.7, zorder=3)
ax.axvline(0, color=INK2, lw=0.6)
ax.axhline(n - 5.5, color=GRID, lw=0.8)
ax.set_yticks(range(n)); ax.set_yticklabels([l for l, _ in labels][::-1])
ax.set_ylim(-0.6, n - 0.4); ax.set_xlim(-20, 26)
ax.set_xlabel("Adjusted difference vs morning session (%)")
clean_axes(ax)
ax.plot([], [], color=BLUE, marker="o", ms=4.5, lw=1.1, label="Afternoon")
ax.plot([], [], color=ORANGE, marker="s", ms=4.5, lw=1.1, label="Evening")
ax.legend(loc="upper left", frameon=False)
ax.set_title("A  Session differences, NHANES 1999–2018", loc="left", fontweight="bold")

cons = q["consequences"]
items = [("NLR ≥3.0", "nlr_ge_3.0"), ("SII ≥445.22", "sii_ge_445.22"), ("NLR in top\nquartile", "nlr_top_quartile")]
x = np.arange(len(items))
for j, (sess, col, mk) in enumerate((("morning", AQUA, "D"), ("afternoon", BLUE, "o"), ("evening", ORANGE, "s"))):
    vals = [100 * cons[k][f"p_{sess}"] for _, k in items]
    bx.plot(x + (j - 1) * 0.22, vals, ls="none", marker=mk, ms=5, color=col, mec="white", mew=0.7, label=sess.capitalize())
for i, (_, k) in enumerate(items):
    d = cons[k]
    bx.text(i, 100 * max(d["p_morning"], d["p_afternoon"], d["p_evening"]) + 3.5,
            f"Aft. {100 * d['rd_afternoon']:+.1f} pp\n({100 * d['rd_afternoon_lo']:+.1f} to {100 * d['rd_afternoon_hi']:+.1f})",
            ha="center", fontsize=6, color=INK2)
bx.set_xticks(x); bx.set_xticklabels([l for l, _ in items])
bx.set_ylabel("Participants above cutoff (%)"); bx.set_ylim(0, 75); bx.set_xlim(-0.6, len(items) - 0.4)
for sp in ("top", "right"):
    bx.spines[sp].set_visible(False)
bx.grid(axis="y", color=GRID, lw=0.5); bx.set_axisbelow(True)
bx.legend(frameon=False, loc="upper left", ncol=3, columnspacing=0.9, handletextpad=0.2, borderaxespad=0.2)
bx.set_title("B  Classification at common cutoffs", loc="left", fontweight="bold")
save(fig, "Figure3_continuous_replication")

# ============ 图 4：死亡 HR 与 Δ ============
p2 = R["part2"]
fig, ax = plt.subplots(figsize=(120 * MM, 58 * MM))
rows = [("NHANES III, GLR — all-cause", p2["nhanes3"]["death_all"]), ("NHANES III, GLR — cancer", p2["nhanes3"]["death_cancer"]),
        ("NHANES 1999–2018, NLR — all-cause", p2["continuous"]["death_all"]),
        ("NHANES 1999–2018, NLR — cancer", p2["continuous"]["death_cancer"])]
n = len(rows)
ax.axvspan(-10, 10, color=BAND, lw=0, zorder=0)
for i, (lab, d) in enumerate(rows):
    y = n - 1 - i
    ax.plot([d["lo_pct"], d["hi_pct"]], [y, y], color=BLUE, lw=1.2, zorder=2)
    ax.plot(d["delta_pct"], y, marker="o", ms=5, color=BLUE, mec="white", mew=0.8, zorder=3)
    ax.text(1.03, y, f"{d['hr_obs_per_doubling']:.3f} → {d['hr_std_per_doubling']:.3f}", transform=ax.get_yaxis_transform(), va="center", fontsize=6.5, color=INK2)
ax.axvline(0, color=INK2, lw=0.6)
ax.set_yticks(range(n)); ax.set_yticklabels([r[0] for r in rows][::-1])
ax.set_xlim(-15, 15); ax.set_ylim(-0.6, n - 0.4)
ax.text(1.03, n - 0.45, "HR per doubling\nobserved → standardized", transform=ax.get_yaxis_transform(), va="bottom", fontsize=6.5, fontweight="bold")
ax.set_xlabel("Change in log hazard ratio per doubling after session standardization (%)")
clean_axes(ax)
save(fig, "Figure4_mortality_robustness")

# ============ 补充图 S1：钟点与禁食（NHANES III，非随机） ============
b = R["part1b"]["all"]
cl = pd.DataFrame(b["clock_vs_0900"]); fa = pd.DataFrame(b["fasting_vs_12h"])
fig, (ax, bx) = plt.subplots(1, 2, figsize=(178 * MM, 62 * MM), gridspec_kw=dict(wspace=0.3))
for a_, d, ref, xl in ((ax, cl, 9, "Clock time of venipuncture (h)"), (bx, fa, 12, "Fasting duration (h)")):
    a_.fill_between(d.x, d.lo, d.hi, color="#cde2fb", lw=0)
    a_.plot(d.x, d.pct, color=BLUE, lw=1.4)
    a_.axhline(0, color=INK2, lw=0.6)
    a_.plot(ref, 0, marker="o", ms=4, color=INK)
    a_.set_xlabel(xl); a_.set_ylabel("Difference in GLR vs reference (%)")
    for sp in ("top", "right"):
        a_.spines[sp].set_visible(False)
    a_.grid(color=GRID, lw=0.5); a_.set_axisbelow(True)
ax.set_title("A  Clock time (reference 09:00)", loc="left", fontweight="bold")
bx.set_title("B  Fasting duration (reference 12 h)", loc="left", fontweight="bold")
save(fig, "FigureS1_clock_fasting")

# ============ 补充图 S2：异质性（探索） ============
het = p1["heterogeneity"]
names = {"sex_f": ("Sex", {"1": "Men", "2": "Women"}), "agegrp_f": ("Age", {"20-39": "20–39 y", "40-59": "40–59 y", "60+": "≥60 y"}),
         "race_f": ("Race/ethnicity", {"1": "Non-Hispanic White", "2": "Non-Hispanic Black", "3": "Mexican American", "4": "Other"}),
         "smoke_f": ("Smoking", {"never": "Never", "former": "Former", "current": "Current"})}
rows = []
for k, (grp, labs) in names.items():
    for lev, lab in labs.items():
        d = het[k]["by_level"][lev]
        rows.append(dict(label=f"{grp}: {lab}", est=d["pct"], lo=d["pct_lo"], hi=d["pct_hi"], color=ORANGE, marker="s", threshold=9.5))
fig, ax = plt.subplots(figsize=(110 * MM, 80 * MM))
forest(ax, rows, (-15, 12), "GLR difference, afternoon/evening vs morning assignment (%)")
yy = len(rows) - 1
for k, (grp, labs) in names.items():
    ax.text(11.8, yy, f"P-int {het[k]['p_interaction']:.3f}", ha="right", va="center", fontsize=6, color=INK2)
    yy -= len(labs)
save(fig, "FigureS2_heterogeneity")

# ============ 补充图 S3：specification curve ============
sp = pd.read_csv(PROJECT / "analysis/results/real/spec_curve.csv").sort_values("delta_pct").reset_index(drop=True)
fig, ax = plt.subplots(figsize=(120 * MM, 55 * MM))
ax.axhspan(-10, 10, color=BAND, lw=0)
q4 = sp.coding == "q4_vs_q1"
ax.plot(sp.index[~q4], sp.delta_pct[~q4], ls="none", marker="o", ms=3.5, color=BLUE, label="Continuous or dichotomous coding")
ax.plot(sp.index[q4], sp.delta_pct[q4], ls="none", marker="s", ms=3.5, color=ORANGE, label="Top vs bottom quartile")
ax.axhline(0, color=INK2, lw=0.6)
ax.set_xlabel("Specification (ranked)"); ax.set_ylabel("Change in log HR after\nstandardization (%)")
for s_ in ("top", "right"):
    ax.spines[s_].set_visible(False)
ax.grid(axis="y", color=GRID, lw=0.5); ax.set_axisbelow(True)
ax.legend(frameon=False, loc="upper left")
ax.set_title("64 specifications: index × coding × covariates × outcome (NHANES 1999–2018)", loc="left", fontsize=7.5)
save(fig, "FigureS3_specification_curve")
