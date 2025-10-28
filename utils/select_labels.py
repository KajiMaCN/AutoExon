#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
from pathlib import Path

# ===== 输入输出路径 =====
in_path  = Path("/workspace/results/match/lyriks/unique_psm_sum_by_sample_with_sample_name_with_meta.csv")
out_path = in_path.with_name(in_path.stem + "_filtered.csv")

# ===== 读取文件 =====
df = pd.read_csv(in_path)

# ===== 过滤 =====
df_filtered = df[df["label"].isin(["Convert", "Non-convert"])].copy()

# ===== 保存 =====
df_filtered.to_csv(out_path, index=False)

print(f"[OK] 已保留 label 为 Convert / Non-convert 的行，共 {len(df_filtered)} 条。")
print(f"[OUT] {out_path}")
