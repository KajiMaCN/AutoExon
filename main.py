import argparse
import json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd

from src.exon_fetcher import GeneIsoformExonFetcher
from src.psm import IsoformCoverageCalculator


def _load_genes(summary_path: Path, gene_arg: str) -> list[str]:
    if gene_arg and gene_arg.strip().upper() != "ALL":
        return [gene_arg.strip()]
    df = pd.read_csv(summary_path)
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


def process_one_gene(gene: str, args, gene2prot) -> tuple[pd.DataFrame | None, dict | None]:
    """
    单基因任务：返回
      - uniq_df: 该基因 unique-exon 内的 PSM 明细（by_sample 子集），或 None
      - stats:   简单统计信息（可选），或 None
    所有“落盘”的动作仍在各自基因目录下完成。
    """
    try:
        out_dir = Path(f"/workspace/results/{gene}")
        out_dir.mkdir(parents=True, exist_ok=True)
        print(f"\n[INFO] >>> 处理基因: {gene}")

        # 拿一个 UniProt accession（优先 canonical）
        ids_set = gene2prot.get(gene)
        if not ids_set:
            print(f"[WARN] 报告中没找到该基因的 UniProt ID：{gene}，跳过")
            return None, None

        ids_list = list(ids_set)
        canonical_ids = [x for x in ids_list if "-" not in x]
        picked_acc = canonical_ids[0] if canonical_ids else ids_list[0]
        canonical_acc = picked_acc.split("-")[0].strip()

        # 各线程各自实例化，避免共享 Session/状态
        fetcher = GeneIsoformExonFetcher(timeout=30)
        calc = IsoformCoverageCalculator('/workspace/src/pms.r')

        # 1) 获取 isoform/exon
        try:
            fetched = fetcher.query_by_uniprot_accession(canonical_acc, args.tax_id, gene_symbol=gene)
        except Exception as e:
            print(f"[WARN] 获取 {gene}/{canonical_acc} 信息失败：{e}，跳过")
            return None, None

        if not fetched or not isinstance(fetched, dict) or not fetched.get("isoforms"):
            print(f"[WARN] {gene}/{canonical_acc} 未返回 isoforms，跳过")
            return None, None

        # 2) 写矩阵 CSV
        try:
            fetcher.build_tables_for_saving(fetched, out_prefix=str(out_dir / gene))
        except Exception as e:
            print(f"[WARN] 写入矩阵 CSV 失败：{e}")

        # 3) 生成 JSON
        try:
            json_path = out_dir / f"{gene}_exons.json"
            json_data = fetcher.save_json(fetched, str(json_path))
        except Exception as e:
            print(f"[WARN] 生成/保存 JSON 失败：{e}，跳过该基因后续计算")
            return None, None

        # 4) 选择一个“有序列”的 isoform
        isoforms = json_data.get("isoforms") or []
        iso_with_seq = [x for x in isoforms if isinstance(x.get("sequence"), str) and x["sequence"]]
        if not iso_with_seq:
            print(f"[WARN] {gene}/{canonical_acc} 所有 isoform 均无有效序列，跳过计算")
            return None, None

        def has_exon(it):
            return bool(it.get("exons")) and any(
                isinstance(ex.get("start"), int) and isinstance(ex.get("end"), int)
                for ex in (it.get("exons") or [])
            )
        iso = next((x for x in iso_with_seq if has_exon(x)), iso_with_seq[0])

        acc = iso.get("accession", "NA")
        seq_len = len(iso.get("sequence", "")) if isinstance(iso.get("sequence"), str) else 0
        print(f"[INFO] 选用 isoform: {acc}  | 序列长度={seq_len} | 外显子数={len(iso.get('exons') or [])}")

        # 5) 运行覆盖计算
        out_sample_csv = str(out_dir / f"{gene}_by_sample.csv")
        out_label_csv  = str(out_dir / f"{gene}_by_label.csv")
        try:
            by_label, by_sample = calc.run_and_save(
                report_path=args.report,
                gene=gene,
                json_data=json_data,
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
            return None, None

        # 6) 基于 JSON 的 unique 区间筛选 unique-exon PSM
        iso2uniq = {}
        for iso_entry in (json_data.get("isoforms") or []):
            acc_i = iso_entry.get("accession")
            uniq_intervals = []
            for ex in (iso_entry.get("exons") or []):
                if ex.get("unique", False) is True:
                    ps = ex.get("start"); pe = ex.get("end")
                    if isinstance(ps, int) and isinstance(pe, int):
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
            return None, stats

        sub_all = []
        for iso_acc, intervals in iso2uniq.items():
            sub_iso = by_sample.loc[by_sample["isoform"] == iso_acc].copy()
            if sub_iso.empty:
                continue
            pos = sub_iso["Position"].astype(int)
            m = pd.Series(False, index=sub_iso.index)
            for a, b in intervals:
                m |= (pos >= a) & (pos <= b)
            sub_iso = sub_iso.loc[m].copy()
            if sub_iso.empty:
                continue
            sub_iso["gene"] = gene
            sub_iso["canonical"] = canonical_acc
            sub_iso = sub_iso[["gene", "canonical", "isoform", "Character", "Position", "id", "label", "value"]]
            sub_all.append(sub_iso)

        uniq_df = pd.concat(sub_all, ignore_index=True) if sub_all else None
        return uniq_df, stats

    except Exception as e:
        print(f"[WARN] 任务异常（{gene}）：{e}")
        return None, None


def main():
    ap = argparse.ArgumentParser(
        description="用 GeneIsoformExonFetcher + IsoformCoverageCalculator 并发计算覆盖；输出到 /workspace/results/{Gene}"
    )
    ap.add_argument("--report", default="/workspace/datasets/Match_result-X401SC24071912_Z01_F001_B1_43/report-matched-M0.csv")
    ap.add_argument("--gene", default="ALL")
    ap.add_argument("--summary", default="/workspace/data/combined_69plus10_summary.csv")
    ap.add_argument("--tax-id", default="9606")

    ap.add_argument("--isoform",   default="ALL")
    ap.add_argument("--id-col",    default="Polypeptide.Novogene.ID")
    ap.add_argument("--label-col", default="label")
    ap.add_argument("--seq-col",   default="Stripped.Sequence")
    ap.add_argument("--protein-col", default="Protein.Group")
    ap.add_argument("--genes-col", default="Genes")
    ap.add_argument("--area-col",  default="Precursor.Normalised")
    ap.add_argument("--quantify",  default="PSM", choices=["PSM", "Area"])
    ap.add_argument("--agg",       default="mean", choices=["mean", "sum", "median"])

    # 并发度
    ap.add_argument("--workers", type=int, default=64, help="并发线程数")

    args = ap.parse_args()

    # ==== 准备基因列表 ====
    summary_path = Path(args.summary)
    if not summary_path.exists():
        raise SystemExit(f"[ERR] summary 文件不存在：{summary_path}")
    genes = _load_genes(summary_path, args.gene)
    if not genes:
        raise SystemExit("[ERR] 没有可处理的基因（检查 has_isoforms 列或 --gene 参数）")

    print(f"[INFO] 计划处理 {len(genes)} 个基因：{', '.join(genes[:10])}{' ...' if len(genes) > 10 else ''}")

    # 读 gene->protein 映射（一次即可；线程中只读不会修改）
    gene2prot = IsoformCoverageCalculator.read_gene_protein_map(args.report, args.genes_col, args.protein_col)

    # ==== 并发提交 ====
    uniq_psm_all_rows: list[pd.DataFrame] = []
    uniq_ranges_stats: list[dict] = []

    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as ex:
        futures = {ex.submit(process_one_gene, gene, args, gene2prot): gene for gene in genes}
        for fut in as_completed(futures):
            gene = futures[fut]
            try:
                uniq_df, stats = fut.result()
                if uniq_df is not None and not uniq_df.empty:
                    uniq_psm_all_rows.append(uniq_df)
                if stats:
                    uniq_ranges_stats.append(stats)
            except Exception as e:
                print(f"[WARN] {gene} 结果收集失败：{e}")

    # ==== 汇总写盘 ====
    if uniq_psm_all_rows:
        uniq_all_df = pd.concat(uniq_psm_all_rows, ignore_index=True)

        out_dir = Path("/workspace/results/_summary")
        out_dir.mkdir(parents=True, exist_ok=True)

        # 明细
        out_psm = out_dir / "unique_psm_by_sample_all.csv"
        uniq_all_df.to_csv(out_psm, index=False)
        print(f"[OK] 写出 unique-exon PSM 明细：{out_psm}")

        # 按 sample(label) 汇总
        sum_df = (
            uniq_all_df
            .groupby(["label"], as_index=False)["value"]
            .sum()
            .rename(columns={"value": "value_sum"})
        )
        out_sum = out_dir / "unique_psm_sum_by_sample.csv"
        sum_df.to_csv(out_sum, index=False)
        print(f"[OK] 写出 unique-exon 按 sample 汇总：{out_sum}")

        # 按 id 聚合（保留 gene/canonical/isoform/label 维度）
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
