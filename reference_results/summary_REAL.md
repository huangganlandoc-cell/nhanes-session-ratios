# 结果摘要（REAL）

运行时间：2026-09-26 19:14:21 PDT；方案：v1.0 (2026-09-26)；R version 4.6.1 (2026-06-24)，survey 4.5

## Part 1：场次随机分配对 GLR 的效应（NHANES III，ITT）

- 样本：15476 人（下午/晚上组 7752，上午组 7724）
- **主估计**：下午/晚上 vs 上午 -1.56%（95% CI -3.57 至 +0.49）；BRR 置信区间 -3.50 至 +0.42
- **R1（界值 9.5%）：within desirable bias**；最严界值 6.3%：within desirable bias；点估计所在档：within optimal
- CACE：-1.72%（95% CI -3.85 至 +0.47）（依从：P(D=1|Z=1) = 0.959，P(D=1|Z=0) = 0.051）
- 标准化位移 -0.038 SD；偏 R² 0.036%（0.000–0.132%）
- 重分类（比例差，百分点）：
  - glr_ge_3.0：下午/晚上 13.6% vs 上午 15.1%，差 -1.5（-3.0 至 +0.1）；归因比例 -10.9%
  - glr_top_quartile：下午/晚上 24.3% vs 上午 25.7%，差 -1.4（-3.8 至 +0.9）；归因比例 -5.9%
  - siig_top_quartile：下午/晚上 25.1% vs 上午 24.8%，差 +0.3（-1.5 至 +2.2）；归因比例 +1.2%

| 敏感性分析 | 效应 | R1 |
|---|---|---|
| S1_unweighted | -3.03%（95% CI -4.30 至 -1.75） | within desirable bias |
| S2_adjusted | -1.71%（95% CI -3.65 至 +0.27） | within desirable bias |
| S3_per_protocol | -2.26%（95% CI -4.28 至 -0.21） | within desirable bias |
| S4_ipw_cbc | -1.60%（95% CI -3.62 至 +0.46） | within desirable bias |
| S5_excl_possible_infection | -1.31%（95% CI -3.45 至 +0.89） | within desirable bias |

| 次要指标 | 效应 | 界值 | R1 |
|---|---|---|---|
| wbc | +11.97%（95% CI +10.40 至 +13.57） | 5.1 | exceeds desirable bias |
| lym | +13.62%（95% CI +11.83 至 +15.43） | 6.3 | exceeds desirable bias |
| gran | +11.85%（95% CI +9.85 至 +13.88） | 7.2 | exceeds desirable bias |
| mono_nuc | +3.25%（95% CI +0.50 至 +6.08） | 6.7 | within desirable bias |
| plt | +2.26%（95% CI +1.17 至 +3.36） | 5 | within desirable bias |
| sii_g | +0.67%（95% CI -1.77 至 +3.16） | 10.8 | within desirable bias |
| plr | -9.99%（95% CI -11.80 至 -8.15） | — | — |

CRP 可检出比例差：+0.7 个百分点（-1.2 至 +2.7）

## Part 1c：连续 NHANES 1999–2018 复现（非随机）

- 样本 48021；NLR 下午 vs 上午 +4.16%（95% CI +3.03 至 +5.31）；晚上 vs 上午 -2.68%（95% CI -4.03 至 -1.32）
- **R3：replicated (consistent null)**
- 2021–2023（n = 5672）：下午/晚上 vs 上午 +0.06%（95% CI -2.78 至 +2.97）
- nlr_ge_3.0：上午 16.5%，下午 18.0%（差 +1.5，+0.6 至 +2.5），晚上 13.7%
- sii_ge_445.22：上午 55.2%，下午 61.2%（差 +6.0，+4.8 至 +7.2），晚上 55.3%
- nlr_top_quartile：上午 24.8%，下午 26.6%（差 +1.9，+0.7 至 +3.0），晚上 21.6%

## Part 2：场次标准化对死亡 HR 的影响

| 数据 | 结局 | 事件 | HR/翻倍（观察） | HR/翻倍（标准化） | Δ（ln HR 相对变化） | R2 |
|---|---|---|---|---|---|---|
| nhanes3 | death_all | 6634 | 1.208 | 1.208 | +0.04%（-0.82 至 +0.90） | robust |
| nhanes3 | death_cancer | 1449 | 1.163 | 1.162 | -0.73%（-3.08 至 +1.62） | robust |
| continuous | death_all | 7654 | 1.302 | 1.300 | -0.42%（-1.31 至 +0.47） | robust |
| continuous | death_cancer | 1662 | 1.058 | 1.060 | +1.86%（-9.05 至 +12.77） | indeterminate |

解析预测的 Δ：NHANES III +0.04%；连续 NHANES +0.31%

## 按方案第 7 节对应的结论措辞

- death_all：R1 = within desirable bias，R2 = robust → 在随机化证据下，NHANES 采血场次对 GLR/NLR 没有临床意义的影响——对现有文献是有检出力的阴性结论。
- death_cancer：R1 = within desirable bias，R2 = robust → 在随机化证据下，NHANES 采血场次对 GLR/NLR 没有临床意义的影响——对现有文献是有检出力的阴性结论。
