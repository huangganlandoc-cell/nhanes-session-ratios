# Part 1c：连续 NHANES 按实际场次复现（方案 5.3；判定规则 R3）

library(parallel)
NCORES <- max(1, min(12, detectCores() - 2))

cd <- read_derived("continuous")

# ---------- 1999–2018 ----------
c18 <- cd[cd$cycle != "2021-2023" & !is.na(cd$w_mec_pooled) & cd$w_mec_pooled > 0, ]
c18 <- make_covs_c(c18, c18$in_main)
c18$lnN <- log(c18$nlr)
tim_c <- c18$in_main                                 # 主样本已限定禁食 <24 h
FC_K <- quantile(c18$fast_h[tim_c], c(1/3, 2/3)); FC_B <- range(c18$fast_h[tim_c])
okf <- !is.na(c18$fast_h)
Bf <- matrix(NA_real_, nrow(c18), 3)
Bf[okf, ] <- ns(pmin(pmax(c18$fast_h[okf], FC_B[1]), FC_B[2]), knots = FC_K, Boundary.knots = FC_B)
colnames(Bf) <- paste0("fst_ns", 1:3); c18 <- cbind(c18, Bf)
desc <- svydesign(ids = ~psu, strata = ~strata, weights = ~w_mec_pooled, nest = TRUE, data = c18)
dm <- subset(desc, in_main)
DF_C <- degf(desc)
main18 <- c18[c18$in_main, ]
cat("连续 NHANES 1999–2018 主样本 n =", nrow(main18), " 设计自由度 =", DF_C, "\n")

f1c <- function(y, extra = "") as.formula(paste(y, "~ sess_f +", COVS_C, "+ cycle_f", extra))
fit_sess <- function(y, design, extra = "") {
  m <- svyglm(f1c(y, extra), design = design)
  list(afternoon = ci_pct_list(coef(m)[["sess_f1"]], SE(m)[["sess_f1"]], DF_C),
       evening = ci_pct_list(coef(m)[["sess_f2"]], SE(m)[["sess_f2"]], DF_C), n = nrow(design$variables))
}
p1c <- list(n_main = nrow(main18), df = DF_C)
p1c$nlr_primary <- fit_sess("lnN", dm)
m_nlr <- svyglm(f1c("lnN"), design = dm)
DELTA_C <- c(aft = coef(m_nlr)[["sess_f1"]], eve = coef(m_nlr)[["sess_f2"]])

sec <- list()
for (v in c("sii", "siri", "piv", "plr", "wbc", "lym", "neu", "mono", "plt", "glr")) {
  r <- fit_sess(paste0("log(", v, ")"), subset(dm, get(v) > 0))
  r$n_nonpositive_excluded <- sum(!(main18[[v]] > 0), na.rm = TRUE)
  sec[[v]] <- r
}
sec$crp_1999_2010 <- fit_sess("log(crp_1999_2010)", subset(dm, !is.na(crp_1999_2010) & crp_1999_2010 > 0))
sec$hscrp_2015_2018 <- fit_sess("log(hscrp)", subset(dm, !is.na(hscrp) & hscrp > 0))
p1c$secondary <- sec
p1c$nlr_with_fasting_spline <- fit_sess("lnN", dm, "+ fst_ns1 + fst_ns2 + fst_ns3")

# R3（方案第 6 节）：以 NHANES III 的 R1 结果为条件
r1 <- RES$part1$primary_glr_taylor$R1
b3 <- RES$part1$primary_glr_taylor$pct
a <- p1c$nlr_primary$afternoon
ci <- c(a$pct_lo, a$pct_hi)
p1c$R3 <- if (r1 == "exceeds desirable bias") {
  if (ci[1] > 0 || ci[2] < 0) { if (sign(a$pct) == sign(b3)) "replicated" else "inconsistent" } else "not replicated"
} else if (r1 == "within desirable bias") {
  if (ci[1] > -M_R1$nlr && ci[2] < M_R1$nlr) "replicated (consistent null)"
  else if (ci[1] > M_R1$nlr || ci[2] < -M_R1$nlr) "inconsistent" else "indeterminate"
} else "descriptive only (NHANES III indeterminate)"

# ---------- 后果终点：边际标准化 + JKn 置信区间 ----------
repc <- as.svrepdesign(desc, type = "JKn")
WREP <- weights(repc, "analysis")
RSC <- repc$rscales; SCL <- repc$scale
jk_var <- function(est, reps) SCL * colSums(RSC * sweep(reps, 2, est)^2)
in_m <- c18$in_main
q75_nlr <- wquantile(main18$nlr, main18$w_mec_pooled, 0.75)
outs <- list(nlr_ge_3.0 = function(d) as.numeric(d$nlr >= 3.0),
             sii_ge_445.22 = function(d) as.numeric(d$sii >= 445.22),
             nlr_top_quartile = function(d) as.numeric(d$nlr > q75_nlr))
fb <- as.formula(paste("yy ~ sess_f +", COVS_C, "+ cycle_f"))
marg <- function(w, yy) {
  # 权重缩放到均值 1（与 svyglm 相同；不改变估计值），否则 logistic 初始值被推到 0/1 附近而不收敛
  dat <- main18; dat$yy <- yy; dat$ww <- w / mean(w[w > 0])
  dat <- dat[dat$ww > 0, ]
  fit <- glm(fb, data = dat, weights = ww, family = quasibinomial(), control = glm.control(maxit = 100))
  if (!fit$converged) stop("边际标准化的 logistic 模型不收敛")
  if (abs(wmean(fitted(fit), dat$ww) - wmean(dat$yy, dat$ww)) > 1e-6) stop("加权平均预测值与观察比例不符")
  pr <- sapply(c(0, 1, 2), function(s) { nd <- dat; nd$sess_f <- factor(s, levels = c(0, 1, 2))
                 wmean(predict(fit, newdata = nd, type = "response"), dat$ww) })
  c(pr, pr[2] - pr[1], pr[3] - pr[1])
}
cons <- list()
for (k in names(outs)) {
  yy <- outs[[k]](main18)
  est <- marg(main18$w_mec_pooled, yy)
  reps <- do.call(rbind, mclapply(seq_len(ncol(WREP)), function(r) marg(WREP[in_m, r], yy), mc.cores = NCORES))
  se <- sqrt(jk_var(est, reps))
  z <- qt(0.975, DF_C)
  cons[[k]] <- list(p_morning = est[1], p_afternoon = est[2], p_evening = est[3],
                    rd_afternoon = est[4], rd_afternoon_lo = est[4] - z * se[4], rd_afternoon_hi = est[4] + z * se[4],
                    rd_evening = est[5], rd_evening_lo = est[5] - z * se[5], rd_evening_hi = est[5] + z * se[5])
}
cons$nlr_q75_cutpoint <- q75_nlr
cons$n_jk_replicates <- ncol(WREP)
p1c$consequences <- cons

# ---------- 2021–2023 ----------
c23 <- cd[cd$cycle == "2021-2023" & !is.na(cd$w_phleb_2123) & cd$w_phleb_2123 > 0, ]
c23 <- make_covs_c(c23, c23$in_main)
c23$sess2_f <- factor(c23$session, levels = c(0, 1))
c23$lnN <- log(c23$nlr)
des23 <- svydesign(ids = ~psu, strata = ~strata, weights = ~w_phleb_2123, nest = TRUE, data = c23)
d23 <- subset(des23, in_main)
DF23 <- degf(des23)
m23 <- svyglm(as.formula(paste("lnN ~ sess2_f +", COVS_C)), design = d23)
p1c$c2021_2023 <- list(n = sum(c23$in_main), df = DF23,
                       pmeve_vs_morning = ci_pct_list(coef(m23)[["sess2_f1"]], SE(m23)[["sess2_f1"]], DF23))
m23h <- svyglm(as.formula(paste("log(hscrp) ~ sess2_f +", COVS_C)), design = subset(d23, !is.na(hscrp) & hscrp > 0))
p1c$c2021_2023$hscrp <- ci_pct_list(coef(m23h)[["sess2_f1"]], SE(m23h)[["sess2_f1"]], DF23)

RES$part1c <- p1c
cat(sprintf("Part 1c NLR 下午 vs 上午：%+.2f%% (%+.2f to %+.2f)  R3 = %s\n", a$pct, a$pct_lo, a$pct_hi, p1c$R3))
cat("Part 1c 完成\n")
