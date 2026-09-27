# Part 1 / 1b：NHANES III 场次随机分配（方案 5.1、5.2；判定规则 R1）

n3 <- read_derived("nhanes3")
ad <- n3[n3$w_mec > 0, ]                                  # MEC 受检成人（设计建在这一层，ITT 作为子域）
ad <- make_covs_n3(ad, ad$in_itt)
ad$lnG <- log(ad$glr)
ad$Z <- ad$arm_pmeve
ad$D <- ad$session_pmeve

# 钟点与禁食样条（节点取自 Part 1b 样本：钟点 4 df，禁食 3 df）
tim_s <- ad$in_timing & !is.na(ad$fast_h) & ad$fast_h < 24
CLK_K <- quantile(ad$draw_clock_h[tim_s], c(.25, .5, .75)); CLK_B <- range(ad$draw_clock_h[tim_s])
FST_K <- quantile(ad$fast_h[tim_s], c(1/3, 2/3)); FST_B <- range(ad$fast_h[tim_s])
basis_clk <- function(x) ns(x, knots = CLK_K, Boundary.knots = CLK_B)
basis_fst <- function(x) ns(x, knots = FST_K, Boundary.knots = FST_B)
ok_t <- !is.na(ad$draw_clock_h) & !is.na(ad$fast_h)
Bc <- matrix(NA_real_, nrow(ad), 4); Bc[ok_t, ] <- basis_clk(ad$draw_clock_h[ok_t]); colnames(Bc) <- paste0("clk_ns", 1:4)
Bf <- matrix(NA_real_, nrow(ad), 3); Bf[ok_t, ] <- basis_fst(pmin(ad$fast_h[ok_t], FST_B[2])); colnames(Bf) <- paste0("fst_ns", 1:3)
ad <- cbind(ad, Bc, Bf)
ad$in_timing24 <- tim_s

des3 <- svydesign(ids = ~psu, strata = ~strata, weights = ~w_mec, nest = TRUE, data = ad)
rep3 <- svrepdesign(data = ad, repweights = "wtpxrp[0-9]+", weights = ~w_mec, type = "Fay", rho = FAY_RHO,
                    combined.weights = TRUE, mse = TRUE)   # 方差以全样本估计为中心（NHANES/WesVar 惯例，与 brr() 一致）
d_itt <- subset(des3, in_itt)
itt <- ad[ad$in_itt, ]
cat("ITT n =", nrow(itt), " Z=1:", sum(itt$Z == 1), " Z=0:", sum(itt$Z == 0), "\n")

# ---------- 5.1 主估计 ----------
m1 <- svyglm(lnG ~ Z, design = d_itt)
b1 <- coef(m1)[["Z"]]; se1 <- SE(m1)[["Z"]]
m1r <- svyglm(lnG ~ Z, design = subset(rep3, in_itt))
se1r <- SE(m1r)[["Z"]]
prim <- ci_pct_list(b1, se1, DF_N3)
prim_brr <- ci_pct_list(b1, se1r, DF_N3)
prim$R1 <- judge_R1(c(prim$pct_lo, prim$pct_hi), M_R1$glr)
prim$R1_strict_6.3 <- judge_R1(c(prim$pct_lo, prim$pct_hi), M_R1$strict)
prim$tier <- tier_of(prim$pct)
RES$part1 <- list(n_itt = nrow(itt), n_z1 = sum(itt$Z == 1), n_z0 = sum(itt$Z == 0),
                  primary_glr_taylor = prim, primary_glr_brr = prim_brr,
                  tiers_glr_pct = as.list(TIERS_GLR))
cat(sprintf("主估计 GLR：%+.2f%% (95%% CI %+.2f to %+.2f)  R1 = %s\n", prim$pct, prim$pct_lo, prim$pct_hi, prim$R1))

# ---------- 同时报告的尺度（BRR） ----------
scale_fun <- function(dat, w) {
  z1 <- dat$Z == 1
  b <- wmean(dat$lnG[z1], w[z1]) - wmean(dat$lnG[!z1], w[!z1])
  p <- wmean(dat$Z, w); v <- wvar(dat$lnG, w)
  c(sd_shift = b / sqrt(v), partial_r2 = b^2 * p * (1 - p) / v)
}
sc <- brr(itt, scale_fun)
RES$part1$sd_shift <- list(est = sc$est[1], lo = sc$lo[1], hi = sc$hi[1])
RES$part1$partial_r2 <- list(est = sc$est[2], lo = max(0, sc$lo[2]), hi = sc$hi[2])

reclass_fun <- function(var, cut = NULL) function(dat, w) {
  x <- dat[[var]]
  hi <- if (is.null(cut)) as.numeric(x > wquantile(x, w, 0.75)) else as.numeric(x >= cut)
  z1 <- dat$Z == 1
  p1 <- wmean(hi[z1], w[z1]); p0 <- wmean(hi[!z1], w[!z1])
  c(p1 = p1, p0 = p0, rd = p1 - p0, af = (p1 - p0) / p1)
}
rc <- list()
for (spec in list(list("glr_ge_3.0", "glr", 3.0), list("glr_top_quartile", "glr", NULL), list("siig_top_quartile", "sii_g", NULL))) {
  r <- brr(itt, reclass_fun(spec[[2]], spec[[3]]))
  rc[[spec[[1]]]] <- list(p1 = r$est[1], p0 = r$est[2], rd = r$est[3], rd_lo = r$lo[3], rd_hi = r$hi[3],
                          af = r$est[4], af_lo = r$lo[4], af_hi = r$hi[4])
}
RES$part1$reclassification <- rc

# ---------- CACE ----------
cace_fun <- function(dat, w) {
  z1 <- dat$Z == 1
  b <- wmean(dat$lnG[z1], w[z1]) - wmean(dat$lnG[!z1], w[!z1])
  p1 <- wmean(dat$D[z1], w[z1]); p0 <- wmean(dat$D[!z1], w[!z1])
  c(cace = b / (p1 - p0), p1 = p1, p0 = p0)
}
cc <- brr(itt, cace_fun)
RES$part1$cace <- c(ci_pct_list(cc$est[1], cc$se[1], DF_N3), list(p1 = cc$est[2], p0 = cc$est[3]))
DELTA_CACE <- cc$est[1]

# ---------- 敏感性分析 S1–S5 ----------
sens <- list()
ad$one <- 1
des_unw <- svydesign(ids = ~psu, strata = ~strata, weights = ~one, nest = TRUE, data = ad)
m <- svyglm(lnG ~ Z, design = subset(des_unw, in_itt))
sens$S1_unweighted <- ci_pct_list(coef(m)[["Z"]], SE(m)[["Z"]], DF_N3)
m <- svyglm(as.formula(paste("lnG ~ Z +", COVS_N3)), design = d_itt)
sens$S2_adjusted <- ci_pct_list(coef(m)[["Z"]], SE(m)[["Z"]], DF_N3)
pp <- itt
pp$w_pp <- ifelse(pp$arm == "AM", pp$w_am, pp$w_pm)
pp <- pp[!is.na(pp$w_pp) & pp$w_pp > 0, ]
des_pp <- svydesign(ids = ~psu, strata = ~strata, weights = ~w_pp, nest = TRUE, data = pp)
m <- svyglm(lnG ~ Z, design = des_pp)
sens$S3_per_protocol <- c(ci_pct_list(coef(m)[["Z"]], SE(m)[["Z"]], DF_N3), list(n = nrow(pp)))
base <- ad$ok_age20 & ad$ok_mec & ad$ok_assigned & ad$ok_not_chemo
ad_b <- make_covs_n3(ad, base)
ad_b$cbc <- as.numeric(ad_b$ok_cbc_ok)
des_b <- svydesign(ids = ~psu, strata = ~strata, weights = ~w_mec, nest = TRUE, data = ad_b)
mm <- svyglm(as.formula(paste("cbc ~ Z +", COVS_N3)), design = subset(des_b, base), family = quasibinomial())
ad$p_cbc <- NA_real_
ad$p_cbc[base] <- predict(mm, newdata = ad_b[base, ], type = "response")
ad$w_ipw <- ad$w_mec / ad$p_cbc
des_ipw <- svydesign(ids = ~psu, strata = ~strata, weights = ~w_ipw, nest = TRUE, data = ad[ad$in_itt, ])
m <- svyglm(lnG ~ Z, design = des_ipw)
sens$S4_ipw_cbc <- ci_pct_list(coef(m)[["Z"]], SE(m)[["Z"]], DF_N3)
m <- svyglm(lnG ~ Z, design = subset(d_itt, is.na(pe_infection) | pe_infection != 2))
sens$S5_excl_possible_infection <- c(ci_pct_list(coef(m)[["Z"]], SE(m)[["Z"]], DF_N3),
                                     list(n_excluded = sum(itt$pe_infection == 2, na.rm = TRUE)))
for (k in names(sens)) sens[[k]]$R1 <- judge_R1(c(sens[[k]]$pct_lo, sens[[k]]$pct_hi), M_R1$glr)
RES$part1$sensitivity <- sens

# ---------- 异质性（探索） ----------
ad$agegrp_f <- cut(ad$age, c(19, 39, 59, Inf), labels = c("20-39", "40-59", "60+"))
des3 <- svydesign(ids = ~psu, strata = ~strata, weights = ~w_mec, nest = TRUE, data = ad)
d_itt <- subset(des3, in_itt)
het <- list()
for (mod in c("sex_f", "agegrp_f", "race_f", "smoke_f")) {
  dd <- if (mod == "smoke_f") subset(d_itt, smoke != "missing") else d_itt
  m <- svyglm(as.formula(paste0("lnG ~ Z * ", mod)), design = dd)
  p_int <- regTermTest(m, as.formula(paste0("~Z:", mod)))$p
  lev <- levels(droplevels(dd$variables[[mod]]))
  strata_eff <- lapply(lev, function(l) {
    ms <- svyglm(lnG ~ Z, design = subset(dd, get(mod) == l))
    ci_pct_list(coef(ms)[["Z"]], SE(ms)[["Z"]], DF_N3)
  })
  names(strata_eff) <- lev
  het[[mod]] <- list(p_interaction = as.numeric(p_int), by_level = strata_eff)
}
RES$part1$heterogeneity <- het

# ---------- 次要指标 ----------
sec_specs <- list(wbc = list("wbc", M_R1$wbc), lym = list("lym", M_R1$lym), gran = list("gran", M_R1$neu),
                  mono_nuc = list("mono_nuc", M_R1$mono), plt = list("plt", M_R1$plt),
                  sii_g = list("sii_g", M_R1$sii), plr = list("plr", NA))
sec <- list()
for (k in names(sec_specs)) {
  v <- sec_specs[[k]][[1]]; mthr <- sec_specs[[k]][[2]]
  ok <- itt[[v]] > 0
  m <- svyglm(as.formula(paste0("log(", v, ") ~ Z")), design = subset(d_itt, get(v) > 0))
  r <- ci_pct_list(coef(m)[["Z"]], SE(m)[["Z"]], DF_N3)
  r$n_nonpositive_excluded <- sum(!ok)
  r$threshold <- mthr
  r$R1 <- if (is.na(mthr)) NA else judge_R1(c(r$pct_lo, r$pct_hi), mthr)
  sec[[k]] <- r
}
ad$crp_detect <- as.numeric(ad$crp > 0.21)
des3 <- svydesign(ids = ~psu, strata = ~strata, weights = ~w_mec, nest = TRUE, data = ad)
d_itt <- subset(des3, in_itt)
m <- svyglm(crp_detect ~ Z, design = subset(d_itt, !is.na(crp_detect)))
sec$crp_detectable_rd <- list(rd = coef(m)[["Z"]], lo = coef(m)[["Z"]] - qt(.975, DF_N3) * SE(m)[["Z"]],
                              hi = coef(m)[["Z"]] + qt(.975, DF_N3) * SE(m)[["Z"]], n = sum(!is.na(itt$crp)))
RES$part1$secondary <- sec

# 人工分类 NLR 与 GLR 的一致性（描述，不分组）
mm_ <- itt[!is.na(itt$nlr_manual) & itt$nlr_manual > 0, ]
RES$part1$manual_nlr_agreement <- list(n = nrow(mm_), spearman = cor(mm_$glr, mm_$nlr_manual, method = "spearman"),
                                       median_ratio = median(mm_$glr / mm_$nlr_manual))

# ---------- 5.2 Part 1b：钟点与禁食（非随机，探索） ----------
d_tim <- subset(des3, in_timing24)
f1b <- as.formula(paste("lnG ~ clk_ns1 + clk_ns2 + clk_ns3 + clk_ns4 + fst_ns1 + fst_ns2 + fst_ns3 +", COVS_N3))
curve_fit <- function(design) {
  m <- svyglm(f1b, design = design)
  V <- vcov(m); bt <- coef(m)
  cc <- paste0("clk_ns", 1:4); fc <- paste0("fst_ns", 1:3)
  grid_c <- seq(8, 21, by = 0.5); ref_c <- basis_clk(9)
  Xc <- sweep(basis_clk(grid_c), 2, ref_c)
  grid_f <- seq(0, 23, by = 1); ref_f <- basis_fst(12)
  Xf <- sweep(basis_fst(grid_f), 2, ref_f)
  mk <- function(X, cols, grid) {
    est <- as.vector(X %*% bt[cols]); se <- sqrt(rowSums((X %*% V[cols, cols]) * X))
    data.frame(x = grid, pct = pct(est), lo = pct(est - qt(.975, DF_N3) * se), hi = pct(est + qt(.975, DF_N3) * se))
  }
  list(clock_vs_0900 = mk(Xc, cc, grid_c), fasting_vs_12h = mk(Xf, fc, grid_f), n = nrow(design$variables))
}
p1b <- list(all = curve_fit(d_tim), within_AM = curve_fit(subset(d_tim, Z == 0)), within_PMEVE = curve_fit(subset(d_tim, Z == 1)))
ma <- svyglm(lnG ~ Z, design = d_tim)
mb <- svyglm(lnG ~ Z + fst_ns1 + fst_ns2 + fst_ns3, design = d_tim)
p1b$itt_without_vs_with_fasting <- list(without = ci_pct_list(coef(ma)[["Z"]], SE(ma)[["Z"]], DF_N3),
                                        with_fasting_spline = ci_pct_list(coef(mb)[["Z"]], SE(mb)[["Z"]], DF_N3))
RES$part1b <- p1b
write.csv(p1b$all$clock_vs_0900, file.path(OUTDIR, "part1b_clock_curve.csv"), row.names = FALSE)
write.csv(p1b$all$fasting_vs_12h, file.path(OUTDIR, "part1b_fasting_curve.csv"), row.names = FALSE)
cat("Part 1 / 1b 完成\n")
