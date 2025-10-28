#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
from pathlib import Path

# ===== 输入路径 =====
meta_path = Path("/workspace/datasets/metadata/experimental-metadata.csv")
psm_path  = Path("/workspace/results/match/unique_psm_sum_by_sample_nopsm0.csv")
out_path  = psm_path.with_name("unique_psm_sum_by_sample_with_sample_name.csv")

# ===== 读取文件 =====
meta_df = pd.read_csv(meta_path)
psm_df  = pd.read_csv(psm_path)

# ===== 检查关键列 =====
if "Polypeptide.Novogene.ID" not in meta_df.columns:
    raise ValueError("experimental-metadata.csv 中缺少列 'Polypeptide.Novogene.ID'")
if "Polypeptide.Novogene.ID" not in psm_df.columns:
    raise ValueError("unique_psm_sum_by_sample_nolabel.csv 中缺少列 'Polypeptide.Novogene.ID'")
if "Sample.Name" not in meta_df.columns:
    raise ValueError("experimental-metadata.csv 中缺少列 'Sample.Name'")

# ===== 合并匹配 =====
merged = pd.merge(
    psm_df,
    meta_df[["Polypeptide.Novogene.ID", "Sample.Name"]],
    on="Polypeptide.Novogene.ID",
    how="left"
)

# ===== 打印未匹配到的 ID =====
unmatched = merged[merged["Sample.Name"].isna()]["Polypeptide.Novogene.ID"].unique()
if len(unmatched) > 0:
    print("⚠️ 未匹配到的 Polypeptide.Novogene.ID：")
    for uid in unmatched:
        print("  -", uid)
    print(f"共 {len(unmatched)} 个未匹配。")
else:
    print("✅ 所有 ID 均成功匹配。")

# ===== 将 Sample.Name 放到第一列 =====
cols = ["Sample.Name"] + [c for c in merged.columns if c != "Sample.Name"]
merged = merged[cols]

# ===== 保存结果 =====
merged.to_csv(out_path, index=False)
print(f"[OK] 已生成：{out_path}")
