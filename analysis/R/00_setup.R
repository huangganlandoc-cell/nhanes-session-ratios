# 公共设置：读数据（带盲态闸门）、抽样设计、协变量构造、重复权重工具函数。
# 方案：analysis/prereg/protocol_v1.0_zh.md（v1.0，2026-09-26 审定）。

suppressPackageStartupMessages({
  library(survey); library(survival); library(splines); library(data.table); library(jsonlite)
})

options(warn = 1)                                     # 警告即时打印，便于在日志里定位是哪个模型
ROOT <- normalizePath(file.path(getwd(), "..", ".."))
stopifnot(dir.exists(file.path(ROOT, "analysis", "R")))          # 须在 analysis/R 目录下运行

# ---------- 盲态闸门 ----------
# 默认只读置乱数据；对真实数据运行需要 MODE=real 且 OSF 登记文件存在并含 osf.io 链接
MODE <- Sys.getenv("MODE", "scrambled")
if (MODE == "real") {
  reg <- file.path(ROOT, "analysis/prereg/OSF_REGISTRATION.txt")
  if (!file.exists(reg) || !any(grepl("osf\\.io", readLines(reg, warn = FALSE))))
    stop("尚未登记 OSF（analysis/prereg/OSF_REGISTRATION.txt 不存在或无 osf.io 链接），拒绝对真实数据运行")
  TAG <- "REAL"
} else {
  MODE <- "scrambled"; TAG <- "SCRAMBLED"
}
OUTDIR <- file.path(ROOT, "analysis/results", tolower(TAG))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
banner <- function() if (TAG == "SCRAMBLED") cat("\n########  SCRAMBLED – NOT REAL  ########\n\n")
banner()

read_derived <- function(stem) {
  f <- file.path(ROOT, "analysis/derived", paste0(stem, "_", TAG, ".csv.gz"))
  d <- fread(cmd = paste("gzip -dc", shQuote(f)))
  for (v in names(d)) if (is.character(d[[v]]) && all(d[[v]] %in% c("True", "False", "")))
    set(d, j = v, value = d[[v]] == "True")
  as.data.frame(d)
}

# ---------- 常量（方案第 6 节） ----------
DF_N3 <- 49                                   # NHANES III 设计自由度（98 PSU − 49 层）
FAY_RHO <- 0.3
BV <- list(wbc = c(11.1, 16.9), neu = c(12.5, 25.9), lym = c(10.5, 22.8), mono = c(14.0, 23.0), plt = c(7.3, 18.6))
db <- function(cvi, cvg, k = 0.25) k * sqrt(cvi^2 + cvg^2)
CV_NLR <- c(sqrt(BV$neu[1]^2 + BV$lym[1]^2), sqrt(BV$neu[2]^2 + BV$lym[2]^2))
CV_SII <- c(sqrt(BV$plt[1]^2 + BV$neu[1]^2 + BV$lym[1]^2), sqrt(BV$plt[2]^2 + BV$neu[2]^2 + BV$lym[2]^2))
M_R1 <- list(glr = 9.5, nlr = 9.5, sii = 10.8, strict = 6.3,       # 登记的界值（取一位小数）
             wbc = 5.1, neu = 7.2, lym = 6.3, mono = 6.7, plt = 5.0)
TIERS_GLR <- c(optimal = db(CV_NLR[1], CV_NLR[2], 0.125), desirable = db(CV_NLR[1], CV_NLR[2], 0.25),
               minimum = db(CV_NLR[1], CV_NLR[2], 0.375))
M_R2 <- 10

pct <- function(b) (exp(b) - 1) * 100

# ---------- 判定规则 ----------
judge_R1 <- function(ci_pct, m) {
  if (ci_pct[1] > -m && ci_pct[2] < m) "within desirable bias"
  else if (ci_pct[1] > m || ci_pct[2] < -m) "exceeds desirable bias"
  else "indeterminate"
}
judge_R2 <- function(ci_rel_pct, m = M_R2) {
  if (ci_rel_pct[1] > -m && ci_rel_pct[2] < m) "robust"
  else if (ci_rel_pct[1] > m || ci_rel_pct[2] < -m) "materially affected"
  else "indeterminate"
}
tier_of <- function(p, tiers = TIERS_GLR) {
  a <- abs(p)
  if (a <= tiers["optimal"]) "within optimal" else if (a <= tiers["desirable"]) "within desirable"
  else if (a <= tiers["minimum"]) "within minimum" else "beyond minimum"
}

# ---------- 加权分位数（加权经验分布函数首次 ≥ p 的值） ----------
wquantile <- function(x, w, p) {
  o <- order(x); x <- x[o]; w <- w[o]
  cw <- cumsum(w) / sum(w)
  x[which(cw >= p)[1]]
}
wmean <- function(x, w) sum(w * x) / sum(w)
wvar <- function(x, w) { m <- wmean(x, w); sum(w * (x - m)^2) / sum(w) }

# ---------- Fay BRR：θ 在全样本权重与 52 个重复权重下各算一次 ----------
brr <- function(data, fun, wcol = "w_mec", R = 52, rho = FAY_RHO, df = DF_N3) {
  est <- fun(data, data[[wcol]])
  reps <- vapply(seq_len(R), function(r) fun(data, data[[paste0("wtpxrp", r)]]), numeric(length(est)))
  reps <- matrix(reps, nrow = length(est))
  se <- sqrt(rowSums((reps - est)^2) / (R * (1 - rho)^2))
  list(est = est, se = se, lo = est - qt(0.975, df) * se, hi = est + qt(0.975, df) * se)
}

# ---------- 协变量（核心集，方案第 4 节）；样条节点取自给定分析样本 ----------
# 年龄：自然样条 4 df（内节点 25/50/75 百分位，边界节点取全距）；BMI：缺失以中位数填补 + 缺失指示，自然样条 3 df
make_covs_n3 <- function(d, in_sample) {
  s <- d[in_sample, ]
  age_k <- quantile(s$age, c(.25, .5, .75)); age_b <- range(s$age)
  bmi_med <- median(s$bmi, na.rm = TRUE)
  d$bmi_missing <- as.integer(is.na(d$bmi))
  d$bmi_filled <- ifelse(is.na(d$bmi), bmi_med, d$bmi)
  bf <- d$bmi_filled[in_sample]
  bmi_k <- quantile(bf, c(1/3, 2/3)); bmi_b <- range(bf)
  A <- ns(d$age, knots = age_k, Boundary.knots = age_b); colnames(A) <- paste0("age_ns", 1:4)
  B <- ns(d$bmi_filled, knots = bmi_k, Boundary.knots = bmi_b); colnames(B) <- paste0("bmi_ns", 1:3)
  d <- cbind(d[, setdiff(names(d), c(colnames(A), colnames(B)))], as.data.frame(unclass(A)), as.data.frame(unclass(B)))
  d$sex_f <- factor(d$sex); d$race_f <- factor(d$race_eth); d$phase_f <- factor(d$phase)
  d$smoke_f <- factor(d$smoke, levels = c("never", "former", "current", "missing"))
  d$educ_f <- factor(ifelse(is.na(d$educ_years), "missing", ifelse(d$educ_years < 9, "<9",
                     ifelse(d$educ_years < 12, "9-11", ifelse(d$educ_years == 12, "12", ">12")))),
                     levels = c("12", "<9", "9-11", ">12", "missing"))
  d$pir_f <- factor(ifelse(is.na(d$pir), "missing", ifelse(d$pir < 1, "<1", ifelse(d$pir < 2, "1-2",
                    ifelse(d$pir < 4, "2-4", ">=4")))), levels = c("2-4", "<1", "1-2", ">=4", "missing"))
  d
}
COVS_N3 <- "age_ns1 + age_ns2 + age_ns3 + age_ns4 + sex_f + race_f + phase_f + smoke_f + bmi_ns1 + bmi_ns2 + bmi_ns3 + bmi_missing + educ_f + pir_f"

make_covs_c <- function(d, in_sample) {
  s <- d[in_sample, ]
  age_k <- quantile(s$age, c(.25, .5, .75)); age_b <- range(s$age)
  bmi_med <- median(s$bmi, na.rm = TRUE)
  d$bmi_missing <- as.integer(is.na(d$bmi))
  d$bmi_filled <- ifelse(is.na(d$bmi), bmi_med, d$bmi)
  bf <- d$bmi_filled[in_sample]
  bmi_k <- quantile(bf, c(1/3, 2/3)); bmi_b <- range(bf)
  A <- ns(d$age, knots = age_k, Boundary.knots = age_b); colnames(A) <- paste0("age_ns", 1:4)
  B <- ns(d$bmi_filled, knots = bmi_k, Boundary.knots = bmi_b); colnames(B) <- paste0("bmi_ns", 1:3)
  d <- cbind(d[, setdiff(names(d), c(colnames(A), colnames(B)))], as.data.frame(unclass(A)), as.data.frame(unclass(B)))
  d$sex_f <- factor(d$sex); d$race_f <- factor(d$race_eth)
  d$smoke_f <- factor(d$smoke, levels = c("never", "former", "current", "missing"))
  d$educ_f <- factor(ifelse(is.na(d$educ), "missing", as.character(d$educ)), levels = c("3", "1", "2", "4", "5", "missing"))
  d$pir_f <- factor(ifelse(is.na(d$pir), "missing", ifelse(d$pir < 1, "<1", ifelse(d$pir < 2, "1-2",
                    ifelse(d$pir < 4, "2-4", ">=4")))), levels = c("2-4", "<1", "1-2", ">=4", "missing"))
  d$cycle_f <- factor(d$cycle)
  d$season_f <- factor(ifelse(is.na(d$exam_season), "missing", d$exam_season))
  d$sess_f <- factor(d$session, levels = c(0, 1, 2))
  d
}
COVS_C <- "age_ns1 + age_ns2 + age_ns3 + age_ns4 + sex_f + race_f + smoke_f + bmi_ns1 + bmi_ns2 + bmi_ns3 + bmi_missing + educ_f + pir_f + season_f"

# ---------- 结果收集 ----------
RES <- list(meta = list(mode = TAG, protocol = "v1.0 (2026-09-26)", run_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
                        R = R.version.string, survey = as.character(packageVersion("survey"))))
ci_pct_list <- function(b, se, df) {
  lo <- b - qt(0.975, df) * se; hi <- b + qt(0.975, df) * se
  list(beta = b, se = se, df = df, pct = pct(b), pct_lo = pct(lo), pct_hi = pct(hi))
}
save_results <- function() {
  write(toJSON(RES, auto_unbox = TRUE, digits = NA, pretty = TRUE, na = "null"),
        file.path(OUTDIR, paste0("results_", TAG, ".json")))
}
