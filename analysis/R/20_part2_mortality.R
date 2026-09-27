# Part 2：场次变异对"指标—死亡"HR 的影响（方案 5.4；判定规则 R2）
# 依赖 10_part1_nhanes3.R（ad、cace_fun）与 12_part1c_continuous.R（c18、main18、WREP、RSC、SCL、DF_C）

library(parallel)
NCORES <- max(1, min(12, detectCores() - 2))

cox_beta <- function(dat, w, event, xvar, covs, extra = "") {
  keep <- w > 0
  f <- as.formula(paste0("Surv(fu_months_exam, ", event, ") ~ ", xvar, " + ", covs, extra))
  fit <- coxph(f, data = dat[keep, ], weights = w[keep], ties = "breslow")
  b <- unname(coef(fit)[xvar])
  if (!is.finite(b) || abs(b) > 20) stop("暴露系数不收敛：", xvar, " = ", b)   # 协变量零事件类别发散无妨，暴露发散必须停
  b
}
rel_ci <- function(est, se, df) list(delta = est, delta_pct = 100 * est, se = se,
                                     lo_pct = 100 * (est - qt(0.975, df) * se), hi_pct = 100 * (est + qt(0.975, df) * se))

# ======================= NHANES III（主） =======================
M_s <- ad$in_itt & ad$eligstat == 1 & !is.na(ad$fu_months_exam)
adm <- make_covs_n3(ad, M_s)                                  # 协变量样条节点取自死亡样本
itt_df <- adm[adm$in_itt, ]
mort <- adm[M_s, ]
mort$x_obs <- log2(mort$glr)
cat("NHANES III 死亡样本 n =", nrow(mort), " 全因死亡 =", sum(mort$death_all), " 癌症死亡 =", sum(mort$death_cancer), "\n")

delta_n3 <- function(w_itt, w_m, event, extra = "") {
  dl <- cace_fun(itt_df, w_itt)[["cace"]]
  d2 <- mort; d2$x_std <- d2$x_obs - (dl / log(2)) * d2$D
  bo <- cox_beta(d2, w_m, event, "x_obs", COVS_N3, extra)
  bs <- cox_beta(d2, w_m, event, "x_std", COVS_N3, extra)
  c(beta_obs = bo, beta_std = bs, delta = (bs - bo) / bo, cace = dl)
}
p2 <- list()
for (ev in c("death_all", "death_cancer")) {
  est <- delta_n3(itt_df$w_mec, mort$w_mec, ev)
  reps <- do.call(rbind, mclapply(1:52, function(r) {
    wc <- paste0("wtpxrp", r); delta_n3(itt_df[[wc]], mort[[wc]], ev) }, mc.cores = NCORES))
  se <- sqrt(colSums(sweep(reps, 2, est)^2) / (52 * (1 - FAY_RHO)^2))
  r <- rel_ci(est[["delta"]], se[["delta"]], DF_N3)
  r$R2 <- judge_R2(c(r$lo_pct, r$hi_pct))
  r$beta_obs <- est[["beta_obs"]]; r$beta_std <- est[["beta_std"]]
  r$hr_obs_per_doubling <- exp(est[["beta_obs"]]); r$hr_std_per_doubling <- exp(est[["beta_std"]])
  r$se_beta_obs_brr <- se[["beta_obs"]]; r$se_beta_std_brr <- se[["beta_std"]]
  r$events <- sum(mort[[ev]])
  p2$nhanes3[[ev]] <- r
}
# 解析预测
DL <- cace_fun(itt_df, itt_df$w_mec)[["cace"]]
pD <- wmean(mort$D, mort$w_mec); vL <- wvar(log(mort$glr), mort$w_mec)
lam <- 1 - DL^2 * pD * (1 - pD) / vL
p2$nhanes3$analytic <- list(delta_session = DL, p_pmeve = pD, var_ln_index = vL, lambda = lam, predicted_delta_pct = 100 * (1 / lam - 1))

# 设计法（Taylor）HR 与比例风险检验；每 SD HR；四分位
mort$x_std <- mort$x_obs - (DL / log(2)) * mort$D
des_m <- svydesign(ids = ~psu, strata = ~strata, weights = ~w_mec, nest = TRUE, data = mort)
sec3 <- list()
for (ev in c("death_all", "death_cancer")) {
  out <- list()
  for (xv in c("x_obs", "x_std")) {
    f <- as.formula(paste0("Surv(fu_months_exam, ", ev, ") ~ ", xv, " + ", COVS_N3))
    sf <- svycoxph(f, design = des_m, method = "breslow")      # svycoxph 用 method= 指定 Breslow（ties= 会被当成变量）
    b <- coef(sf)[[xv]]; s <- SE(sf)[[xv]]
    sdx <- sqrt(wvar(mort[[xv]], mort$w_mec))
    out[[xv]] <- list(hr_per_doubling = exp(b), lo = exp(b - qt(.975, DF_N3) * s), hi = exp(b + qt(.975, DF_N3) * s),
                      hr_per_sd = exp(b * sdx), sd_log2 = sdx)
  }
  fit <- coxph(as.formula(paste0("Surv(fu_months_exam, ", ev, ") ~ x_obs + ", COVS_N3)), data = mort,
               weights = mort$w_mec, ties = "breslow")
  out$ph_global_p <- cox.zph(fit)$table["GLOBAL", "p"]
  qo <- sapply(c(.25, .5, .75), function(p) wquantile(mort$x_obs, mort$w_mec, p))
  qs <- sapply(c(.25, .5, .75), function(p) wquantile(mort$x_std, mort$w_mec, p))
  mort$q_obs <- factor(findInterval(mort$x_obs, qo, left.open = TRUE) + 1)
  mort$q_std <- factor(findInterval(mort$x_std, qs, left.open = TRUE) + 1)
  hq <- sapply(c("q_obs", "q_std"), function(q) {
    fit <- coxph(as.formula(paste0("Surv(fu_months_exam, ", ev, ") ~ ", q, " + ", COVS_N3)), data = mort,
                 weights = mort$w_mec, ties = "breslow")
    exp(coef(fit)[[paste0(q, "4")]]) })
  out$q4_vs_q1 <- list(hr_obs = hq[["q_obs"]], hr_std = hq[["q_std"]],
                       share_changing_quartile = wmean(as.numeric(mort$q_obs != mort$q_std), mort$w_mec))
  ext <- " + hx_cvd_any + hx_cancer_any"
  mort$hx_cvd_any <- factor(ifelse(is.na(mort$hx_mi) & is.na(mort$hx_chf) & is.na(mort$hx_stroke), "missing",
                            ifelse(pmax(mort$hx_mi, mort$hx_chf, mort$hx_stroke, na.rm = TRUE) == 1, "yes", "no")))
  mort$hx_cancer_any <- factor(ifelse(is.na(mort$hx_other_cancer), "missing", ifelse(mort$hx_other_cancer == 1, "yes", "no")))
  bo <- cox_beta(mort, mort$w_mec, ev, "x_obs", COVS_N3, ext); bs <- cox_beta(mort, mort$w_mec, ev, "x_std", COVS_N3, ext)
  out$extended_covariates_delta_pct <- 100 * (bs - bo) / bo
  sec3[[ev]] <- out
}
# 阴性对照结局：意外伤害死亡（只作定性）
acc <- sapply(c("x_obs", "x_std"), function(xv) {
  fit <- coxph(as.formula(paste0("Surv(fu_months_exam, death_accident) ~ ", xv, " + ", COVS_N3)), data = mort,
               weights = mort$w_mec, ties = "breslow"); exp(coef(fit)[[xv]]) })
sec3$accident_negative_control <- list(events = sum(mort$death_accident), hr_obs = acc[["x_obs"]], hr_std = acc[["x_std"]])
p2$nhanes3_secondary <- sec3

# ======================= 连续 NHANES 1999–2018（复现） =======================
Mc_s <- c18$in_main & c18$eligstat == 1 & !is.na(c18$fu_months_exam)
cm <- make_covs_c(c18, Mc_s)[Mc_s, ]                           # Cox 用：样条节点取自死亡样本
cm$x_obs <- log2(cm$nlr)
f_delta <- as.formula(paste("lnN ~ sess_f +", COVS_C, "+ cycle_f"))
idx_main <- which(c18$in_main); idx_mort <- which(Mc_s)
cat("连续 NHANES 死亡样本 n =", nrow(cm), " 全因 =", sum(cm$death_all), " 癌症 =", sum(cm$death_cancer), "\n")

delta_c <- function(w_main, w_m, event, extra = "") {
  keep <- w_main > 0
  dd <- main18[keep, ]; dd$.w <- w_main[keep]          # 权重作为数据列传入（lm 在公式环境里找权重变量）
  fit <- lm(f_delta, data = dd, weights = .w)
  da <- coef(fit)[["sess_f1"]]; de <- coef(fit)[["sess_f2"]]
  d2 <- cm; d2$x_std <- d2$x_obs - (da * (d2$session == 1) + de * (d2$session == 2)) / log(2)
  bo <- cox_beta(d2, w_m, event, "x_obs", paste(COVS_C, "+ cycle_f"), extra)
  bs <- cox_beta(d2, w_m, event, "x_std", paste(COVS_C, "+ cycle_f"), extra)
  c(beta_obs = bo, beta_std = bs, delta = (bs - bo) / bo, d_aft = da, d_eve = de)
}
for (ev in c("death_all", "death_cancer")) {
  est <- delta_c(main18$w_mec_pooled, cm$w_mec_pooled, ev)
  reps <- do.call(rbind, mclapply(seq_len(ncol(WREP)), function(r)
    delta_c(WREP[idx_main, r], WREP[idx_mort, r], ev), mc.cores = NCORES))
  se <- sqrt(SCL * colSums(RSC * sweep(reps, 2, est)^2))
  r <- rel_ci(est[["delta"]], se[["delta"]], DF_C)
  r$R2 <- judge_R2(c(r$lo_pct, r$hi_pct))
  r$beta_obs <- est[["beta_obs"]]; r$beta_std <- est[["beta_std"]]
  r$hr_obs_per_doubling <- exp(est[["beta_obs"]]); r$hr_std_per_doubling <- exp(est[["beta_std"]])
  r$events <- sum(cm[[ev]])
  p2$continuous[[ev]] <- r
}
off <- DELTA_C[["aft"]] * (cm$session == 1) + DELTA_C[["eve"]] * (cm$session == 2)
lamc <- 1 - wvar(off, cm$w_mec_pooled) / wvar(log(cm$nlr), cm$w_mec_pooled)
p2$continuous$analytic <- list(delta_aft = DELTA_C[["aft"]], delta_eve = DELTA_C[["eve"]], lambda = lamc,
                               predicted_delta_pct = 100 * (1 / lamc - 1))
cm$x_std <- cm$x_obs - off / log(2)
accc <- sapply(c("x_obs", "x_std"), function(xv) {
  fit <- coxph(as.formula(paste0("Surv(fu_months_exam, death_accident) ~ ", xv, " + ", COVS_C, " + cycle_f")),
               data = cm, weights = cm$w_mec_pooled, ties = "breslow"); exp(coef(fit)[[xv]]) })
p2$continuous_secondary <- list(accident_negative_control = list(events = sum(cm$death_accident),
                                hr_obs = accc[["x_obs"]], hr_std = accc[["x_std"]]))
RES$part2 <- p2
for (s in c("nhanes3", "continuous")) for (ev in c("death_all", "death_cancer")) {
  r <- p2[[s]][[ev]]
  cat(sprintf("%s %s：Δ = %+.2f%% (%+.2f to %+.2f)  R2 = %s\n", s, ev, r$delta_pct, r$lo_pct, r$hi_pct, r$R2))
}
cat("Part 2 完成\n")
