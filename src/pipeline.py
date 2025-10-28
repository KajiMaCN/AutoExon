# -*- coding: utf-8 -*-
from pathlib import Path
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm.auto import tqdm
import pandas as pd

from src.report_mapping import build_gene2prot_from_report
from src.cache_table import ensure_cache_table, load_genes_from_mapping
from src.uniprot_prefetch import prefetch_gene_json_map
from src.pepmapviz_runner import pepmapviz_thread_worker
from src.statistics import process_one_gene

def run_pipeline(args):
    # 1) Gene→Protein
    gene2prot = build_gene2prot_from_report(Path(args.report), args.genes_col, args.protein_col)

    # 2) 缓存状态表
    cache_csv = Path(args.out_cache) / "gene_protein_status.csv"
    cache_df = ensure_cache_table(cache_csv, gene2prot)  # noqa: F841

    # 3) 基因列表
    genes = load_genes_from_mapping(gene2prot, args.gene)
    print(f"[INFO] 示例前 10 个基因：{', '.join(genes[:10])}{' ...' if len(genes) > 10 else ''}")

    # 4) UniProt 预取
    prefetch_gene_json_map(
        genes, gene2prot, args.tax_id, Path(args.out_uniport), cache_csv
    )

    # 5) PepMapViz 线程（等待完成后再统计）
    pep_done = threading.Event()
    t = threading.Thread(
        target=pepmapviz_thread_worker,
        args=(args, genes, cache_csv, pep_done),
        daemon=True
    )
    t.start()
    pep_done.wait()

    # 6) 并发统计（只用缓存与 JSON，不再调用 R）
    print("[STEP] 并发处理各基因（仅统计） ...")
    uniq_psm_all_rows: list[pd.DataFrame] = []
    uniq_ranges_stats: list[dict] = []

    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as ex:
        futures = {ex.submit(process_one_gene, gene, args, gene2prot): gene for gene in genes}
        for fut in tqdm(as_completed(futures), total=len(futures), desc="基因整体进度"):
            gene = futures[fut]
            try:
                uniq_df, stats = fut.result()
                if uniq_df is not None and not uniq_df.empty:
                    uniq_psm_all_rows.append(uniq_df)
                if stats:
                    uniq_ranges_stats.append(stats)
            except Exception as e:
                print(f"[WARN] {gene} 结果收集失败：{e}")

    # 7) 写盘
    print("[STEP] 写出结果 ...")
    out_dir = Path(args.out_match)
    out_dir.mkdir(parents=True, exist_ok=True)

    if uniq_psm_all_rows:
        uniq_all_df = pd.concat(uniq_psm_all_rows, ignore_index=True)

        print("[STEP] 生成明细表 ...")
        merged = uniq_all_df.copy()
        cols_order_detail = [
            "Polypeptide.Novogene.ID",
            "gene", "canonical", "isoform", "label",
            "start", "end", "value",
            "Precursor.Normalised.sum", "Global.Q.Value"
        ]
        for c in cols_order_detail:
            if c not in merged.columns:
                merged[c] = pd.NA
        merged = merged[cols_order_detail]

        out_psm = out_dir / "unique_psm_by_sample_all.csv"
        merged.to_csv(out_psm, index=False)
        print(f"[OK] 写出明细：{out_psm}")

        # 汇总 1：按 id + pep 区间
        print("[STEP] 生成汇总（按 id + pep 区间） ...")
        sum_by_id = (
            merged
            .groupby(
                ["Polypeptide.Novogene.ID", "gene", "canonical", "isoform", "label", "start", "end"],
                as_index=False
            )
            .agg({
                "value": "sum",
                "Precursor.Normalised.sum": "sum"
            })
            .rename(columns={
                "value": "psm_value_sum",
                "Precursor.Normalised.sum": "Precursor.Normalised_sum"
            })
        )
        out_sum_id = out_dir / "unique_psm_sum_by_sample.csv"
        sum_by_id.to_csv(out_sum_id, index=False)
        print(f"[OK] 写出：{out_sum_id}")

        # 汇总 2：按 pep 区间（不区分 id）
        print("[STEP] 生成汇总（按 pep 区间） ...")
        sum_by_pep = (
            merged
            .groupby(
                ["gene", "canonical", "isoform", "label", "start", "end"],
                as_index=False
            )
            .agg({
                "value": "sum",
                "Precursor.Normalised.sum": "sum"
            })
            .rename(columns={
                "value": "psm_value_sum",
                "Precursor.Normalised.sum": "Precursor.Normalised_sum"
            })
        )
        out_sum_pep = out_dir / "unique_psm_sum_by_exon.csv"  # 兼容旧下游命名
        sum_by_pep.to_csv(out_sum_pep, index=False)
        print(f"[OK] 写出：{out_sum_pep}")
    else:
        print("[INFO] 没有任何命中行，未生成汇总文件。")
