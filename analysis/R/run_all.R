# 用法：cd analysis/R && Rscript run_all.R            （默认置乱数据）
#       MODE=real Rscript run_all.R                  （需 OSF 登记文件，见 00_setup.R）
#       Rscript run_all.R 10 12                      （只跑指定脚本编号）
invisible(Sys.setlocale("LC_ALL", "en_US.UTF-8"))       # 路径含中文，必须 UTF-8
args <- commandArgs(trailingOnly = TRUE)
source("00_setup.R")
scripts <- c("10" = "10_part1_nhanes3.R", "12" = "12_part1c_continuous.R", "20" = "20_part2_mortality.R", "30" = "30_exploratory.R")
todo <- if (length(args)) scripts[args] else scripts
t0 <- Sys.time()
# 开发用缓存：RESUME=1 时从缓存恢复 Part 1/1c 的对象，只跑后面的脚本（仅限置乱数据；揭盲运行一律从头跑）
CACHE <- file.path(Sys.getenv("NHANES_DATA_DIR", file.path(Sys.getenv("HOME"), "data", "nhanes_public")), "cache", paste0("after_part1c_", TAG, ".RData"))
if (Sys.getenv("RESUME") == "1" && TAG == "SCRAMBLED" && file.exists(CACHE)) {
  load(CACHE); t0 <- Sys.time(); todo <- todo[!names(todo) %in% c("10", "12")]; cat("已从缓存恢复 Part 1 / 1c\n")
}
# 循环变量用 .step：子脚本在同一环境里运行，里面的 k 等变量会覆盖主循环变量
for (.step in names(todo)) {
  cat("\n=====", todo[[.step]], "=====\n"); source(todo[[.step]]); save_results()
  if (.step == "12" && TAG == "SCRAMBLED") { dir.create(dirname(CACHE), recursive = TRUE, showWarnings = FALSE); save(list = ls(all.names = TRUE), file = CACHE) }
}
cat("\n完成，用时", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "分钟；结果：", OUTDIR, "\n")
banner()
