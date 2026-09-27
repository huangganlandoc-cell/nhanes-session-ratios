# 次要与探索性分析（方案 5.4 次要项、5.5）。均不设判定规则。
# 依赖 10、12、20 中的对象（ad、c18、main18、mort、cm、DELTA_C、cace_fun 等）

library(parallel)
NCORES <- max(1, min(12, detectCores() - 2))
EXP <- list()

# ---------- 1. 限制性立方样条拐点（节点 5/35/65/95 百分位），标准化前后位移 ----------
# 拐点 = 5–95 百分位范围内预测对数风险的最低点；以样本加权百分位表示；bootstrap 500 次（subbootstrap）
nadir <- function(x, dat, w, event, covs, kn) {
  B <- ns(x, knots = kn[2:3], Boundary.knots = kn[c(1, 4)]); colnames(B) <- paste0("rcs", 1:3)
  d2 <- cbind(dat, B); keep <- w > 0
  f <- as.formula(paste0("Surv(fu_months_exam, ", event, ") ~ rcs1 + rcs2 + rcs3 + ", covs))
  fit <- coxph(f, data = d2[keep, ], weights = w[keep], ties = "breslow")
  cr <- coef(fit)[paste0("rcs", 1:3)]
  if (any(!is.finite(cr)) || max(abs(cr)) > 50) return(c(value = NA_real_, pctl = NA_real_, interior = NA_real_))  # 样条发散的重抽样记为缺失
  g <- seq(kn[1], kn[4], length.out = 400)
  lp <- as.vector(ns(g, knots = kn[2:3], Boundary.knots = kn[c(1, 4)]) %*% coef(fit)[paste0("rcs", 1:3)])
  xm <- g[which.min(lp)]
  c(value = xm, pctl = 100 * sum(w[x <= xm]) / sum(w), interior = as.numeric(which.min(lp) > 1 & which.min(lp) < 400))
}
rcs_shift <- function(dat, xobs, event, covs, get_w, get_xstd, R) {
  # get_w(r)、get_xstd(r)：r = 0 为全样本，r ≥ 1 为第 r 个自助重复样本（各行权重均取自同一重复样本）
  kn_o <- quantile(xobs, c(.05, .35, .65, .95))
  one <- function(r) {
    w <- get_w(r); xs <- get_xstd(r); kn_s <- quantile(xs, c(.05, .35, .65, .95))
    o <- nadir(xobs, dat, w, event, covs, kn_o); s <- nadir(xs, dat, w, event, covs, kn_s)
    c(obs_value = o[["value"]], obs_pctl = o[["pctl"]], obs_interior = o[["interior"]],
      std_value = s[["value"]], std_pctl = s[["pctl"]], std_interior = s[["interior"]], shift_pctl = s[["pctl"]] - o[["pctl"]])
  }
  est <- one(0)
  reps <- do.call(rbind, mclapply(seq_len(R), function(r) tryCatch(one(r), error = function(e) rep(NA_real_, 7)),
                                  mc.cores = NCORES))
  ok <- is.finite(reps[, 7])
  c(as.list(est), list(shift_lo = unname(quantile(reps[ok, 7], .025)), shift_hi = unname(quantile(reps[ok, 7], .975)),
                       n_boot_ok = sum(ok)))
}
SEED <- 20260926; set.seed(SEED)
# NHANES III：GLR；自助法重复权重取自 MEC 成人设计（行序与 adm 相同），每个重复样本内重估 CACE
des_ad <- svydesign(ids = ~psu, strata = ~strata, weights = ~w_mec, nest = TRUE, data = adm)
bs3 <- as.svrepdesign(des_ad, type = "subbootstrap", replicates = 500)
W3 <- weights(bs3, "analysis")
im <- which(M_s); ii <- which(adm$in_itt)
w3_m <- function(r) if (r == 0) mort$w_mec else W3[im, r]
x3_std <- function(r) {
  dl <- cace_fun(itt_df, if (r == 0) itt_df$w_mec else W3[ii, r])[["cace"]]
  mort$glr * exp(-dl * mort$D)
}
rcs <- list()
for (ev in c("death_all", "death_cancer"))
  rcs[[paste0("nhanes3_glr_", ev)]] <- rcs_shift(mort, mort$glr, ev, COVS_N3, w3_m, x3_std, ncol(W3))
cat("拐点（NHANES III）完成\n")

# 连续 NHANES：NLR（全因、癌症）；SII 2005–2018（癌症，对应 Wu 2025 的 445.22）
bsc <- as.svrepdesign(desc, type = "subbootstrap", replicates = 500)
WC <- weights(bsc, "analysis")                          # 行序与 c18 相同
imc <- which(Mc_s); imain <- which(c18$in_main)
sess_offset <- function(var, w_main) {
  keep <- w_main > 0 & main18[[var]] > 0
  dd <- main18[keep, ]; dd$.w <- w_main[keep]
  fit <- lm(as.formula(paste0("log(", var, ") ~ sess_f + ", COVS_C, " + cycle_f")), data = dd, weights = .w)
  function(sess) coef(fit)[["sess_f1"]] * (sess == 1) + coef(fit)[["sess_f2"]] * (sess == 2)
}
covc <- paste(COVS_C, "+ cycle_f")
wmain <- function(r) if (r == 0) main18$w_mec_pooled else WC[imain, r]
for (ev in c("death_all", "death_cancer"))
  rcs[[paste0("continuous_nlr_", ev)]] <- rcs_shift(cm, cm$nlr, ev, covc,
    function(r) if (r == 0) cm$w_mec_pooled else WC[imc, r],
    function(r) cm$nlr * exp(-sess_offset("nlr", wmain(r))(cm$session)), ncol(WC))
s0518 <- cm$cycle %in% c("2005-2006", "2007-2008", "2009-2010", "2011-2012", "2013-2014", "2015-2016", "2017-2018")
cm05 <- cm[s0518, ]
rcs$continuous_sii_2005_2018_death_cancer <- rcs_shift(cm05, cm05$sii, "death_cancer", covc,
  function(r) if (r == 0) cm05$w_mec_pooled else WC[imc[s0518], r],
  function(r) cm05$sii * exp(-sess_offset("sii", wmain(r))(cm05$session)), ncol(WC))
rcs$continuous_sii_2005_2018_death_cancer$wu2025_reference_value <- 445.22
EXP$rcs_turning_points <- rcs
cat("拐点（连续 NHANES）完成\n")

# ---------- 2. Specification curve（连续 NHANES，64 种设定，只报点估计） ----------
cm$hx_cvd_f <- factor(ifelse(is.na(cm$hx_cvd), "missing", cm$hx_cvd))
cm$hx_cancer_f <- factor(ifelse(is.na(cm$hx_cancer), "missing", cm$hx_cancer))
spec_rows <- list()
for (ix in c("nlr", "sii", "siri", "piv")) {
  ok <- cm[[ix]] > 0
  dd <- cm[ok, ]; w <- dd$w_mec_pooled
  offf <- sess_offset(ix, main18$w_mec_pooled * (main18[[ix]] > 0))
  xo <- dd[[ix]]; xs <- xo * exp(-offf(dd$session))
  cut_fixed <- switch(ix, nlr = 3.0, sii = 445.22, wquantile(xo, w, 0.5))
  codings <- list(
    per_doubling = function(x) data.frame(e = log2(x)),
    per_sd = function(x) data.frame(e = log(x) / sqrt(wvar(log(x), w))),
    q4_vs_q1 = function(x) { q <- sapply(c(.25, .5, .75), function(p) wquantile(x, w, p))
                             data.frame(e = factor(findInterval(x, q, left.open = TRUE) + 1)) },
    dichotomous = function(x) data.frame(e = as.numeric(x >= cut_fixed)))
  for (cd_ in names(codings)) for (cv in c("core", "extended")) for (ev in c("death_all", "death_cancer")) {
    covs <- paste(covc, if (cv == "extended") "+ hx_cvd_f + hx_cancer_f" else "")
    get_b <- function(x) {
      d2 <- cbind(dd, codings[[cd_]](x))
      fit <- coxph(as.formula(paste0("Surv(fu_months_exam, ", ev, ") ~ e + ", covs)), data = d2, weights = w, ties = "breslow")
      cf <- coef(fit); if (cd_ == "q4_vs_q1") cf[["e4"]] else cf[["e"]]
    }
    bo <- get_b(xo); bs <- get_b(xs)
    spec_rows[[length(spec_rows) + 1]] <- data.frame(index = ix, coding = cd_, covariates = cv, outcome = ev,
                                                     hr_obs = exp(bo), hr_std = exp(bs), delta_pct = 100 * (bs - bo) / bo)
  }
}
spec <- do.call(rbind, spec_rows)
write.csv(spec, file.path(OUTDIR, "spec_curve.csv"), row.names = FALSE)
EXP$spec_curve <- list(n_specs = nrow(spec), delta_pct_median = median(spec$delta_pct),
                       delta_pct_range = range(spec$delta_pct), file = "spec_curve.csv")
cat("Specification curve 完成\n")

# ---------- 3. 甘油三酯（餐后状态标志）横断面描述 ----------
d_tg <- subset(des3, in_itt & !is.na(tg) & tg > 0)
m_a <- svyglm(lnG ~ Z, design = d_tg); m_b <- svyglm(lnG ~ Z + log(tg), design = d_tg)
tg3 <- list(without_tg = ci_pct_list(coef(m_a)[["Z"]], SE(m_a)[["Z"]], DF_N3),
            with_ln_tg = ci_pct_list(coef(m_b)[["Z"]], SE(m_b)[["Z"]], DF_N3))
d_tgc <- subset(dm, !is.na(tg_biochem) & tg_biochem > 0)
m_a <- svyglm(f1c("lnN"), design = d_tgc); m_b <- svyglm(f1c("lnN", "+ log(tg_biochem)"), design = d_tgc)
tgc <- list(without_tg = ci_pct_list(coef(m_a)[["sess_f1"]], SE(m_a)[["sess_f1"]], DF_C),
            with_ln_tg = ci_pct_list(coef(m_b)[["sess_f1"]], SE(m_b)[["sess_f1"]], DF_C))
EXP$triglycerides <- list(nhanes3_itt = tg3, continuous_afternoon = tgc)

# ---------- 4. 癌症幸存者与女性乳腺癌幸存者（描述） ----------
surv <- list()
for (g in c("cancer_survivors", "breast_cancer_survivors_female")) {
  sel <- if (g == "cancer_survivors") quote(!is.na(hx_cancer) & hx_cancer == 1) else quote(hx_breast_cancer & sex == 2)
  dg <- subset(dm, eval(sel))
  ng <- nrow(dg$variables)
  # 单一性别亚组（女性乳腺癌幸存者）去掉性别变量，否则性别只有一个水平而报错
  covs_g <- if (g == "breast_cancer_survivors_female") sub("sex_f \\+ ", "", COVS_C) else COVS_C
  r <- tryCatch({
    m <- svyglm(as.formula(paste("lnN ~ sess_f +", covs_g, "+ cycle_f")), design = dg)
    list(afternoon = ci_pct_list(coef(m)[["sess_f1"]], SE(m)[["sess_f1"]], DF_C),
         evening = ci_pct_list(coef(m)[["sess_f2"]], SE(m)[["sess_f2"]], DF_C), covariates = covs_g)
  }, error = function(e) list(error = conditionMessage(e)))
  cg <- cm[with(cm, eval(sel)), ]
  hr <- sapply(c("death_all", "death_cancer"), function(ev) sapply(c("x_obs", "x_std"), function(xv) {
    if (sum(cg[[ev]]) < 5) return(NA_real_)
    # 小亚组只调整年龄（调整 10 个周期时部分周期零事件导致不收敛）；不收敛记为 NA
    fit <- tryCatch(withCallingHandlers(
             coxph(as.formula(paste0("Surv(fu_months_exam, ", ev, ") ~ ", xv, " + age_ns1 + age_ns2 + age_ns3 + age_ns4")),
                   data = cg, weights = cg$w_mec_pooled, ties = "breslow"),
             warning = function(w) if (grepl("did not converge", conditionMessage(w))) stop("not converged")),
           error = function(e) NULL)
    if (is.null(fit) || !is.finite(coef(fit)[[xv]])) NA_real_ else exp(coef(fit)[[xv]]) }))
  surv[[g]] <- list(n = ng, session_effect_nlr = r, n_mort = nrow(cg), deaths_all = sum(cg$death_all),
                    deaths_cancer = sum(cg$death_cancer), hr_per_doubling_age_adjusted = hr)
}
EXP$survivors_descriptive <- surv
RES$exploratory <- EXP
cat("探索性分析完成\n")
