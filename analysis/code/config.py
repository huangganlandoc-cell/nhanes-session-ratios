"""项目路径与全局设置。所有脚本都从这里取路径；原始数据位置可用环境变量 NHANES_DATA_DIR 指定。"""
import os
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[2]          # maosaozi-sci/
RAW = Path(os.environ.get("NHANES_DATA_DIR", Path.home() / "data" / "nhanes_public"))   # 原始公开数据（默认 ~/data/nhanes_public）
N3 = RAW / "nhanes3"                                   # NHANES III 定长文本 + SAS 读入程序
XPT = RAW / "continuous_xpt"                           # 连续 NHANES 1999–2023 的 .xpt
MORT = RAW / "mortality_2019"                          # 连续 NHANES 公共死亡链接（随访至 2019-12-31）
DERIVED = PROJECT / "analysis" / "derived"             # 构建好的分析数据集
OUT = PROJECT / "analysis" / "output"                  # 表、图、数字

SEED = 20260926                                        # 所有随机过程（bootstrap、置换）统一用这个种子
