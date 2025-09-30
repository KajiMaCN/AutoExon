import argparse
import json
from pathlib import Path

import pandas as pd

from src.exon_fetcher import GeneIsoformExonFetcher
from src.psm import IsoformCoverageCalculator


def _load_genes(summary_path: Path, gene_arg: str) -> list[str]:
    """
    根据参数决定处理哪些基因：
      - gene_arg != 'ALL' 时，仅返回该基因（单基因模式）
      - 否则从 summary CSV 里筛选 has_isoforms == True 的所有基因
    """
    if gene_arg and gene_arg.strip().upper() != "ALL":
        return [gene_arg.strip()]

    df = pd.read_csv(summary_path)
    # 兼容 has_isoforms 列是 bool 或字符串 'True'/'False'
    if df["has_isoforms"].dtype == bool:
        mask = df["has_isoforms"]
    else:
        mask = df["has_isoforms"].astype(str).str.lower().eq("true")

    genes = (
        df.loc[mask, "Gene"]
        .dropna()
        .astype(str)
        .unique()
        .tolist()
    )
    return genes


def main():
    ap = argparse.ArgumentParser(
        description="用 GeneIsoformExonFetcher + IsoformCoverageCalculator 端到端计算覆盖；支持批量基因；输出统一到 /workspace/results/{Gene}"
    )
    # 必需/常用参数
    ap.add_argument("--report", default="/workspace/datasets/Match_result-X401SC24071912_Z01_F001_B1_43/report-matched-M0.csv",
                    help="鉴定/定量报告（.tsv/.csv 自动识别分隔符）")
    ap.add_argument("--gene", default="ALL",
                    help="目标基因符号；默认 ALL 表示从 --summary 中读取 has_isoforms==True 的全部基因")
    ap.add_argument("--summary", default="/workspace/data/combined_69plus10_summary.csv",
                    help="基因汇总表（需包含列：Gene, has_isoforms），默认：/workspace/data/combined_69plus10_summary.csv")

    # Fetcher 侧参数
    ap.add_argument("--tax-id", default="9606", help="物种 TaxID，默认 9606（人）")

    # Calculator 侧参数
    ap.add_argument("--isoform",   default="ALL", help="逗号分隔；ALL 表示全部")
    ap.add_argument("--id-col",    default="Polypeptide.Novogene.ID")
    ap.add_argument("--label-col", default="label")
    ap.add_argument("--seq-col",   default="Stripped.Sequence")
    ap.add_argument("--protein-col", default="Protein.Group")
    ap.add_argument("--genes-col", default="Genes")
    ap.add_argument("--area-col",  default="Precursor.Normalised")
    ap.add_argument("--quantify",  default="PSM", choices=["PSM", "Area"])
    ap.add_argument("--agg",       default="mean", choices=["mean", "sum", "median"])

    args = ap.parse_args()

    # ==== 准备基因列表 ====
    summary_path = Path(args.summary)
    if not summary_path.exists():
        raise SystemExit(f"[ERR] summary 文件不存在：{summary_path}")
    genes = _load_genes(summary_path, args.gene)
    if not genes:
        raise SystemExit("[ERR] 没有可处理的基因（检查 has_isoforms 列或 --gene 参数）")

    print(f"[INFO] 计划处理 {len(genes)} 个基因：{', '.join(genes[:10])}{' ...' if len(genes) > 10 else ''}")

    # ==== 初始化工具 ====
    fetcher = GeneIsoformExonFetcher(timeout=30)
    calc = IsoformCoverageCalculator('/workspace/src/pms.r')

    gene2prot = IsoformCoverageCalculator.read_gene_protein_map(args.report, args.genes_col, args.protein_col)
    
    # ==== 逐基因处理 ====
    uniq_psm_all_rows = []   # 存放所有基因的 unique-exon PSM 明细（by_sample 的子集）
    uniq_ranges_stats = []   # 可选：记录每基因/isoform 的 unique 区间数量（用于日志或检查）
    for gene in genes:
        out_dir = Path(f"/workspace/results/{gene}")
        out_dir.mkdir(parents=True, exist_ok=True)

        print(f"\n[INFO] >>> 处理基因: {gene}")

        # 0) 从 gene2prot 拿一个 UniProt accession（可能有多个，优先 canonical）
        ids_set = gene2prot.get(gene)
        if not ids_set:
            print(f"[WARN] 报告中没找到该基因的 UniProt ID：{gene}，跳过")
            continue

        # 可能是 set，安全取一个；并尽量选 canonical（无连字符的）
        ids_list = list(ids_set)
        canonical_ids = [x for x in ids_list if "-" not in x]
        picked_acc = canonical_ids[0] if canonical_ids else ids_list[0]

        # 如果拿到的是 isoform 形式（如 P02768-2），取前缀作为 canonical（P02768）
        canonical_acc = picked_acc.split("-")[0].strip()

        # 1) 在线获取并构建 isoform/exon 信息
        try:
            fetched = fetcher.query_by_uniprot_accession(canonical_acc, args.tax_id, gene_symbol=gene)
        except Exception as e:
            print(f"[WARN] 获取 {gene}/{canonical_acc} 信息失败：{e}，跳过")
            continue

        if not fetched or not isinstance(fetched, dict) or not fetched.get("isoforms"):
            print(f"[WARN] {gene}/{canonical_acc} 未返回 isoforms，跳过")
            continue

        # 2) 保存“矩阵”CSV
        try:
            fetcher.build_tables_for_saving(fetched, out_prefix=str(out_dir / gene))
        except Exception as e:
            print(f"[WARN] 写入矩阵 CSV 失败：{e}")

        # 3) 准备给 Calculator 的 JSON 字典
        try:
            json_path = out_dir / f"{gene}_exons.json"
            json_data = fetcher.save_json(fetched, str(json_path))
        except Exception as e:
            print(f"[WARN] 生成/保存 JSON 失败：{e}，跳过该基因后续计算")
            continue

        # 4) 选择一个“有序列”的 isoform（优先有外显子的）
        isoforms = json_data.get("isoforms") or []
        iso_with_seq = [
            x for x in isoforms
            if isinstance(x.get("sequence"), str) and x["sequence"]
        ]
        if not iso_with_seq:
            print(f"[WARN] {gene}/{canonical_acc} 所有 isoform 均无有效序列，跳过计算")
            continue

        # 优先选择含有至少一个 exon 的 isoform；否则就选第一个有序列的
        def has_exon(it): 
            return bool(it.get("exons")) and any(
                isinstance(ex.get("start"), int) and isinstance(ex.get("end"), int)
                for ex in (it.get("exons") or [])
            )
        iso = next((x for x in iso_with_seq if has_exon(x)), iso_with_seq[0])

        acc = iso.get("accession", "NA")
        seq_len = len(iso.get("sequence", "")) if isinstance(iso.get("sequence"), str) else 0
        print(f"[INFO] 选用 isoform: {acc}  | 序列长度={seq_len} | 外显子数={len(iso.get('exons') or [])}")

        # 5) 运行覆盖计算（若 json_data 内没有任何可用外显子，Calculator 端可能会直接产出空表）
        try:
            out_sample_csv = str(out_dir / f"{gene}_by_sample.csv")
            out_label_csv  = str(out_dir / f"{gene}_by_label.csv")

            by_label, by_sample=calc.run_and_save(
                report_path=args.report,
                gene=gene,
                json_data=json_data,             # 直接传 dict
                out_sample_csv=out_sample_csv,
                out_label_csv=out_label_csv,
                isoform_select=args.isoform,
                id_col=args.id_col,
                label_col=args.label_col,
                seq_col=args.seq_col,
                genes_col=args.genes_col,
                area_col=args.area_col,
                quantify=args.quantify,
                agg=args.agg,
            )
        except Exception as e:
            print(f"[WARN] Calculator 计算失败：{gene} | {e}")
            continue

        # 6) 基于 json_data 的 unique exon 区间，筛出 by_sample 中落入区间的 PSM
        try:
            # 1) 组装 isoform -> [(start,end), ...] 的 unique 区间（蛋白坐标）
            iso2uniq = {}
            for iso_entry in (json_data.get("isoforms") or []):
                acc_i = iso_entry.get("accession")
                uniq_intervals = []
                for ex in (iso_entry.get("exons") or []):
                    # JSON 的 exons 来自 fetcher.save_json（已按蛋白 start/end 填写）
                    if ex.get("unique", False) is True:
                        ps = ex.get("start"); pe = ex.get("end")
                        if isinstance(ps, int) and isinstance(pe, int):
                            if ps > pe: ps, pe = pe, ps
                            uniq_intervals.append((ps, pe))
                if uniq_intervals:
                    # 合并重叠/相邻区间（保障筛选效率与一致性）
                    uniq_intervals.sort()
                    merged = []
                    for a, b in uniq_intervals:
                        if not merged or a > merged[-1][1]:
                            merged.append([a, b])
                        else:
                            merged[-1][1] = max(merged[-1][1], b)
                    iso2uniq[acc_i] = [(x, y) for x, y in merged]

            # 可选：记录统计
            uniq_ranges_stats.append({
                "gene": gene,
                "canonical": canonical_acc,
                "n_isoforms_with_unique": sum(1 for _k, v in iso2uniq.items() if v),
                "n_unique_intervals_total": sum(len(v) for v in iso2uniq.values())
            })

            if iso2uniq:
                # 2) 逐 isoform 过滤 by_sample：Position 落在任一 unique 区间内
                #    by_sample 列：["Character","Position","id","label","value","isoform"]
                sub_all = []
                for iso_acc, intervals in iso2uniq.items():
                    sub_iso = by_sample.loc[by_sample["isoform"] == iso_acc].copy()
                    if sub_iso.empty:
                        continue

                    # 对每个 interval 做筛选后 union（或用向量化做一次性筛选）
                    # 这里用向量化：对每个 interval 构造布尔掩码，之后 OR 起来
                    pos = sub_iso["Position"].astype(int)
                    m = pd.Series(False, index=sub_iso.index)
                    for a, b in intervals:
                        m |= (pos >= a) & (pos <= b)
                    sub_iso = sub_iso.loc[m].copy()
                    if sub_iso.empty:
                        continue

                    # 附加溯源信息
                    sub_iso["gene"] = gene
                    sub_iso["canonical"] = canonical_acc
                    # 列顺序更友好些
                    sub_iso = sub_iso[["gene", "canonical", "isoform", "Character", "Position", "id", "label", "value"]]
                    sub_all.append(sub_iso)

                if sub_all:
                    uniq_psm_all_rows.append(pd.concat(sub_all, ignore_index=True))

        except Exception as e:
            print(f"[WARN] unique-exon PSM 过滤失败：{gene} | {e}")
    
    # ==== 所有基因处理完成后：写出 unique-exon 内的 PSM 明细与按 sample 汇总 ====
    if uniq_psm_all_rows:
        uniq_all_df = pd.concat(uniq_psm_all_rows, ignore_index=True)

        # 明细表
        out_dir = Path("/workspace/results/_summary")
        out_dir.mkdir(parents=True, exist_ok=True)
        out_psm = out_dir / "unique_psm_by_sample_all.csv"
        uniq_all_df.to_csv(out_psm, index=False)
        print(f"[OK] 写出 unique-exon PSM 明细：{out_psm}")

        # 1) 按 sample(label) 汇总相加（跨所有基因）
        sum_df = (
            uniq_all_df
            .groupby(["label"], as_index=False)["value"]
            .sum()
            .rename(columns={"value": "value_sum"})
        )
        out_sum = out_dir / "unique_psm_sum_by_sample.csv"
        sum_df.to_csv(out_sum, index=False)
        print(f"[OK] 写出 unique-exon 按 sample 汇总：{out_sum}")

        # 2) 按 id 聚合（可同时保留 gene/canonical/isoform/label）
        sum_by_id = (
            uniq_all_df
            .groupby(["gene", "canonical", "isoform", "id", "label"], as_index=False)["value"]
            .sum()
            .rename(columns={"value": "value_sum"})
        )
        out_sum_id = out_dir / "unique_psm_sum_by_id.csv"
        sum_by_id.to_csv(out_sum_id, index=False)
        print(f"[OK] 写出 unique-exon 按 id 聚合：{out_sum_id}")

    else:
        print("[INFO] 没有任何 PSM 落入 unique exon 区间，未生成汇总文件。")

if __name__ == "__main__":
    main()
