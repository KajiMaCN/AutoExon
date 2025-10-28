# -*- coding: utf-8 -*-
from pathlib import Path
from typing import Dict, Set, Tuple, List
import json
import pandas as pd
from tqdm.auto import tqdm

from src.psm import RES_CACHE_BY_GENE, load_gene_match_result

def process_one_gene(gene: str, args, gene2prot: Dict[str, Set[str]]) -> Tuple[pd.DataFrame | None, dict | None]:
    """
    仅统计（不再调用 R）：
      - 读取 out_uniport/{gene}/pepmapviz_positions.csv
      - 读取 out_uniport/{gene}/{gene}_exons.json
      - 以 PepMapViz 返回的肽段坐标为基准，若与任意 unique-exon 区间有交集，value=1，否则 value=0。
      - 不合并相邻/相接的肽段坐标（如 1–17 与 17–32 仍为两条）。
    """
    try:
        print(f"\n[STEP] >>> 处理基因: {gene}")

        # 1) 取 pepmapviz 结果：优先磁盘
        sub = load_gene_match_result(gene, Path(args.out_uniport))
        if sub is None or sub.empty:
            sub = RES_CACHE_BY_GENE.get(gene)
        if sub is None or sub.empty:
            print(f"[INFO] {gene} 无匹配结果（pepmapviz）")
            return None, None

        # 2) 读取 JSON
        json_path = Path(args.out_uniport) / gene / f"{gene}_exons.json"
        if not json_path.exists():
            print(f"[WARN] {gene} 缺少 {json_path.name}，跳过")
            return None, None
        with open(json_path, "r", encoding="utf-8") as f:
            json_data = json.load(f)

        # 3) 计算 canonical
        ids_set = gene2prot.get(gene, set())
        ids_list = list(ids_set)
        canonical_ids = [x for x in ids_list if "-" not in x]
        picked_acc = canonical_ids[0] if canonical_ids else (ids_list[0] if ids_list else "")
        canonical_acc = picked_acc.split("-")[0].strip() if picked_acc else ""

        # 4) unique-exon 区间（合并重叠）
        print("[STEP] 计算 unique-exon 区间（按 isoform 合并重叠） ...")
        iso2uniq: Dict[str, List[tuple[int, int]]] = {}
        for iso_entry in (json_data.get("isoforms") or []):
            acc_i = iso_entry.get("accession")
            uniq_intervals = []
            for ex in (iso_entry.get("exons") or []):
                # 注意：把 start/end 为 0 的也保留（你有此需求）
                ps = ex.get("start"); pe = ex.get("end")
                is_unique = ex.get("unique", False) is True
                if is_unique and isinstance(ps, int) and isinstance(pe, int):
                    if ps > pe: ps, pe = pe, ps
                    uniq_intervals.append((ps, pe))
            if uniq_intervals:
                uniq_intervals.sort()
                merged = []
                for a, b in uniq_intervals:
                    if not merged or a > merged[-1][1]:
                        merged.append([a, b])
                    else:
                        merged[-1][1] = max(merged[-1][1], b)
                iso2uniq[acc_i] = [(x, y) for x, y in merged]

        stats = {
            "gene": gene,
            "canonical": canonical_acc,
            "n_isoforms_with_unique": sum(1 for _k, v in iso2uniq.items() if v),
            "n_unique_intervals_total": sum(len(v) for v in iso2uniq.values())
        }
        if not iso2uniq:
            print("[INFO] 无 unique-exon 区间，跳过命中统计")
            return None, stats

        # 5) 命中统计（以 pepmapviz 坐标为基准；未命中也写出 value=0）
        print("[STEP] 统计 unique-exon 命中（pep 坐标为基准；含 value=0） ...")
        rows = []

        isoforms_with_unique = {iso for iso, ivs in iso2uniq.items() if ivs}

        sub = sub.copy()
        if "isoform" not in sub.columns:
            print("[WARN] pepmapviz 结果缺少 isoform 列，跳过该基因")
            return None, stats

        area_col = args.area_col
        if area_col not in sub.columns:
            sub[area_col] = 0

        for _, r in tqdm(sub.iterrows(), total=len(sub), desc=f"{gene} 肽段命中", leave=False):
            iso_acc = str(r.get("isoform", ""))
            if iso_acc not in isoforms_with_unique:
                continue

            # pepmapviz 区间（基准）
            try:
                s1 = int(r["start_position"]); s2 = int(r["end_position"])
            except Exception:
                continue
            if s1 > s2:
                s1, s2 = s2, s1

            # 是否与任一 unique-exon 区间相交
            hit_any = False
            for a, b in iso2uniq[iso_acc]:
                if not (s2 < a or b < s1):
                    hit_any = True
                    break

            # 面积：命中才累加，否则 0
            try:
                area_val_full = float(r.get(area_col, 0)) if pd.notna(r.get(area_col, None)) else 0.0
            except Exception:
                area_val_full = 0.0
            area_val = area_val_full if hit_any else 0.0

            rows.append({
                "Polypeptide.Novogene.ID": r.get(args.id_col, pd.NA),
                "gene": gene,
                "canonical": canonical_acc,
                "isoform": iso_acc,
                "label": r.get(args.label_col, "All"),
                "start": s1,
                "end": s2,
                "value": 1 if hit_any else 0,                 # 命中=1；未命中=0（亦保存）
                "Precursor.Normalised.sum": area_val,          # 未命中时为 0
                "Global.Q.Value": r.get("Global.Q.Value", pd.NA)
            })

        uniq_df = pd.DataFrame(rows) if rows else None
        return uniq_df, stats

    except Exception as e:
        print(f"[WARN] 任务异常（{gene}）：{e}")
        return None, None
