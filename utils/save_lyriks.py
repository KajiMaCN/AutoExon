#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
from pathlib import Path

# ===== 输入文件 =====
meta_path = Path("/workspace/datasets/metadata/metadata-lyriks407.csv")
psm_path  = Path("/workspace/results/match/unique_psm_sum_by_sample_with_sample_name.csv")

# ===== 输出文件 =====
out_path = psm_path.with_name(psm_path.stem + "_with_meta.csv")

# ===== 读取 =====
df_psm  = pd.read_csv(psm_path)
df_meta = pd.read_csv(meta_path)

# 只保留 metadata 中必要列（若缺列则填充为 NA）
need_cols = ["Polypeptide.Novogene.ID", "period", "gender", "age", "label"]
for c in need_cols:
    if c not in df_meta.columns:
        df_meta[c] = pd.NA
df_meta = df_meta[need_cols].copy()

# ===== 内连接：仅保留匹配到的行 =====
merged = pd.merge(
    df_psm,
    df_meta,
    on="Polypeptide.Novogene.ID",
    how="inner"
)

# ===== label 改名为 origin_label，并基于其生成新列 label（放到最后一列）=====
merged = merged.rename(columns={"label": "origin_label"})
s = merged["origin_label"].astype(str)

# 规则：
# 1) "convert" -> "Convert"
# 2) ["maintain","early_remit","late_remit","relapse","remit"] -> "Non-convert"
s = s.replace({"convert": "Convert"})
s = s.mask(s.isin(["maintain", "early_remit", "late_remit", "relapse", "remit"]), "Non-convert")

merged["label"] = s

# 将 label 放到最后一列
cols = [c for c in merged.columns if c != "label"]
merged = merged[cols + ["label"]]

# ===== 保存 =====
merged.to_csv(out_path, index=False, encoding="utf-8")
print(f"[OK] 匹配并转换完成：{len(merged)} 条记录，已保存到：{out_path}")
