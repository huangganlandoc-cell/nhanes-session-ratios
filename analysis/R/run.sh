#!/bin/sh
# 运行全部分析。每个并行子进程只用单线程矩阵运算，避免抢核。
# 用法：sh run.sh              （置乱数据）
#       MODE=real sh run.sh    （需 analysis/prereg/OSF_REGISTRATION.txt，见 00_setup.R）
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 RESUME=${RESUME:-0}
M=${MODE:-scrambled}
mkdir -p ../results/$M
Rscript run_all.R "$@" > ../results/$M/run_log.txt 2>&1
echo "exit=$?"
tail -25 ../results/$M/run_log.txt
