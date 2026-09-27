"""读取 NHANES III 定长文本文件（.dat），列位置从官方 SAS 读入程序（.sas）的 INPUT 段解析。"""
import re
import numpy as np
import pandas as pd
from config import N3


def layout(sasfile):
    """返回 {变量名: (起始列, 结束列, 是否字符型)}，列号从 1 开始，与官方文档一致。"""
    txt = open(N3 / sasfile, encoding="latin-1").read()
    i = txt.index("INPUT")
    j = txt.index(";", i)
    pos = {}
    for m in re.finditer(r"^\s+([A-Z][A-Z0-9_]*)\s+(\$\s*)?(\d+)(?:\s*-\s*(\d+))?\s*$", txt[i + 5:j], re.M):
        a = int(m.group(3))
        b = int(m.group(4) or m.group(3))
        pos[m.group(1)] = (a, b, bool(m.group(2)))
    return pos


def read_raw(datfile, sasfile, varlist):
    """按列位置读出原始字符串（不做任何转换），保留"空白"与"0"的区别。"""
    pos = layout(sasfile)
    missing = [v for v in varlist if v not in pos]
    if missing:
        raise KeyError(f"{sasfile} 中找不到变量: {missing}")
    colspecs = [(pos[v][0] - 1, pos[v][1]) for v in varlist]
    df = pd.read_fwf(N3 / datfile, colspecs=colspecs, names=varlist, header=None,
                     dtype=str, encoding="latin-1", keep_default_na=False)
    return df.apply(lambda s: s.str.strip())


def to_num(s, blank_codes=()):
    """字符串转数值；blank_codes 里的编码（如 '8888' 表示 blank but applicable）记为缺失。"""
    s = s.where(~s.isin(list(blank_codes)) & (s != ""))
    return pd.to_numeric(s, errors="coerce")
