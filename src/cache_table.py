# -*- coding: utf-8 -*-
from pathlib import Path
from typing import Dict, Set, List
import pandas as pd

def ensure_cache_table(cache_csv: Path, gene2prot: Dict[str, Set[str]]) -> pd.DataFrame:
    """
    生成/对齐缓存状态表：
    列：Gene, Protein, UniProtDone, PepMapVizDone（全为 bool）
    """
    cache_csv.parent.mkdir(parents=True, exist_ok=True)

    pairs = []
    for g, ps in gene2prot.items():
        for p in ps:
            pairs.append((g, p))
    base_df = pd.DataFrame(pairs, columns=["Gene", "Protein"]).drop_duplicates()
    base_df["UniProtDone"] = False
    base_df["PepMapVizDone"] = False

    if cache_csv.exists():
        old = pd.read_csv(cache_csv)
        for col in ["Gene", "Protein", "UniProtDone", "PepMapVizDone"]:
            if col not in old.columns:
                old[col] = False if col in ["UniProtDone", "PepMapVizDone"] else pd.NA
        merged = base_df.merge(
            old[["Gene", "Protein", "UniProtDone", "PepMapVizDone"]],
            on=["Gene", "Protein"], how="left", suffixes=("", "_old")
        )
        for col in ["UniProtDone", "PepMapVizDone"]:
            merged[col] = merged[f"{col}_old"].fillna(False).astype(bool)
            merged.drop(columns=[f"{col}_old"], inplace=True)
        merged.to_csv(cache_csv, index=False)
        return merged
    else:
        base_df.to_csv(cache_csv, index=False)
        return base_df

def load_genes_from_mapping(gene2prot: Dict[str, Set[str]], gene_arg: str) -> List[str]:
    print("[STEP] 准备基因列表 ...")
    if gene_arg and gene_arg.strip().upper() != "ALL":
        genes = [gene_arg.strip()]
    else:
        genes = sorted(g for g in gene2prot.keys() if g and g.strip())
    if not genes:
        raise SystemExit("[ERR] 从 report 中未解析到任何 Gene。")
    print(f"[INFO] 计划处理 {len(genes)} 个基因。")
    return genes
