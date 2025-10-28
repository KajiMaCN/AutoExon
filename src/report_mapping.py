# -*- coding: utf-8 -*-
from pathlib import Path
from typing import Dict, Set
import pandas as pd
from tqdm.auto import tqdm

def read_sep(path: str) -> str:
    return "\t" if str(path).endswith((".tsv", ".tsv.gz")) else ","

def build_gene2prot_from_report(report_path: Path, genes_col: str, protein_col: str) -> Dict[str, Set[str]]:
    print("[STEP] 构建 Gene→Protein 映射（直接去重，不切分）...")
    if not report_path.exists():
        raise SystemExit(f"[ERR] report 文件不存在：{report_path}")

    sep = read_sep(str(report_path))
    usecols = [genes_col, protein_col, "Global.Q.Value"]
    df = pd.read_csv(report_path, sep=sep, usecols=[c for c in usecols if c is not None],
                     dtype=str, low_memory=False)
    if df.empty:
        raise SystemExit("[ERR] report 为空或缺少必要列。")

    # 统一过滤：去分号 + Q<0.01
    if genes_col in df.columns:
        df = df[~df[genes_col].astype(str).str.contains(";", na=False)]
    if protein_col in df.columns:
        df = df[~df[protein_col].astype(str).str.contains(";", na=False)]
    if "Global.Q.Value" in df.columns:
        q = pd.to_numeric(df["Global.Q.Value"], errors="coerce")
        df = df[q < 0.01]

    df = df.dropna(subset=[genes_col, protein_col]).drop_duplicates(subset=[genes_col, protein_col])

    mapping: Dict[str, Set[str]] = {}
    for g, p in tqdm(df[[genes_col, protein_col]].itertuples(index=False),
                    total=len(df), desc="聚合映射", leave=False):
        mapping.setdefault(g, set()).add(p)

    print(f"[INFO] 基因数：{len(mapping)}")
    return mapping
