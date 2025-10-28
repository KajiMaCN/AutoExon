#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
from pathlib import Path

# ===== 文件路径 =====
csv_path = Path("/workspace/results/match/unique_psm_sum_by_sample_nolabel.csv")
save_path = csv_path.with_name("unique_psm_sum_by_sample_nopsm0.csv")
# ===== 读取 CSV =====
df = pd.read_csv(csv_path)

# ===== 删除 psm_value_sum 为 0 的行 =====
df_filtered = df[df["psm_value_sum"] != 0]

# ===== 保存结果（覆盖原文件）=====
df_filtered.to_csv(save_path, index=False)

print(f"[OK] 已删除 psm_value_sum=0 的行，剩余 {len(df_filtered)} 条记录。")
