#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import argparse
import time
import json
import threading
import pandas as pd

from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,ProcessPoolExecutor, as_completed
from tqdm.auto import tqdm

from src.exon_fetcher import GeneIsoformExonFetcher
from src.psm import (
    IsoformCoverageCalculator, RES_CACHE_BY_GENE,
    save_gene_match_result, load_gene_match_result
)


def _read_sep(path: str) -> str:
    return "\t" if str(path).endswith((".tsv", ".tsv.gz")) else ","


def _build_gene2prot_from_report(report_path: Path, genes_col: str, protein_col: str) -> dict[str, set[str]]:
    print("[STEP] 构建 Gene→Protein 映射（直接去重，不切分）...")
    if not report_path.exists():
        raise SystemExit(f"[ERR] report 文件不存在：{report_path}")

    sep = _read_sep(str(report_path))
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

    mapping: dict[str, set[str]] = {}
    for g, p in tqdm(df[[genes_col, protein_col]].itertuples(index=False),
                    total=len(df), desc="聚合映射", leave=False):
        mapping.setdefault(g, set()).add(p)

    print(f"[INFO] 基因数：{len(mapping)}")
    return mapping


def _ensure_cache_table(cache_csv: Path, gene2prot: dict[str, set[str]]) -> pd.DataFrame:
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
        # 兼容旧列
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


def _load_genes_from_mapping(gene2prot: dict[str, set[str]], gene_arg: str) -> list[str]:
    print("[STEP] 准备基因列表 ...")
    if gene_arg and gene_arg.strip().upper() != "ALL":
        genes = [gene_arg.strip()]
    else:
        genes = sorted(g for g in gene2prot.keys() if g and g.strip())
    if not genes:
        raise SystemExit("[ERR] 从 report 中未解析到任何 Gene。")
    print(f"[INFO] 计划处理 {len(genes)} 个基因。")
    return genes


def _prefetch_gene_json_map(genes: list[str],
                            gene2prot: dict[str, set[str]],
                            tax_id: str,
                            out_uniport: Path,
                            cache_csv: Path) -> None:
    """
    并发预取 isoform/exon JSON 到 out_uniport/{gene}/{gene}_exons.json
    完成后将该 Gene 所有 pair 的 UniProtDone=True
    """
    print("[STEP] 并发预取 isoform/exon 元数据（UniProt，支持断点续跑） ...")

    # 读取最新缓存表
    cache_df = pd.read_csv(cache_csv)
    if "UniProtDone" not in cache_df.columns:
        cache_df["UniProtDone"] = False
    pending_genes = sorted(cache_df.loc[~cache_df["UniProtDone"], "Gene"].unique().tolist())

    def _task(g: str):
        ids_set = gene2prot.get(g)
        if not ids_set:
            return g, False

        ids_list = list(ids_set)
        canonical_ids = [x for x in ids_list if "-" not in x]
        picked_acc = canonical_ids[0] if canonical_ids else ids_list[0]
        canonical_acc = picked_acc.split("-")[0].strip()

        # 如已有 JSON，直接视为完成
        out_dir = out_uniport / g
        json_path = out_dir / f"{g}_exons.json"
        if json_path.exists():
            return g, True

        fetcher = GeneIsoformExonFetcher(timeout=30)
        try:
            fetched = fetcher.query_by_uniprot_accession(canonical_acc, tax_id, gene_symbol=g)
            if not fetched or not isinstance(fetched, dict) or not fetched.get("isoforms"):
                return g, False

            out_dir.mkdir(parents=True, exist_ok=True)
            try:
                fetcher.build_tables_for_saving(fetched, out_prefix=str(out_dir / g))
            except Exception:
                pass
            try:
                fetcher.save_json(fetched, str(json_path))
            except Exception:
                # 兜底
                with open(json_path, "w", encoding="utf-8") as f:
                    json.dump(fetched, f, ensure_ascii=False)
            return g, True
        except Exception:
            return g, False

    with ThreadPoolExecutor(max_workers=32) as ex:
        futs = {ex.submit(_task, g): g for g in pending_genes}
        for f in tqdm(as_completed(futs), total=len(futs), desc="UniProt 预取进度"):
            g, ok = f.result()
            if ok:
                cache_df.loc[cache_df["Gene"].eq(g), "UniProtDone"] = True
                cache_df.to_csv(cache_csv, index=False)

    # 控制台提示
    ok_genes = cache_df.loc[cache_df["UniProtDone"], "Gene"].nunique()
    print(f"[INFO] UniProt 预取完成/已有：{ok_genes}/{len(genes)} 个基因。")

# ========= 顶层：给多进程调用的子任务函数 =========
def _pepmapviz_run_one_gene(task):
    """
    子进程执行：单基因 PepMapViz
    入参 task: dict，字段：
      - gene, out_uniport, r_script, r_timeout, report, id_col, label_col,
        seq_col, genes_col, area_col, isoform
    返回 (gene, status)：
      status ∈ {"ok","cached","empty","no_json","error"}
    """
    try:
        g = task["gene"]
        out_uniport = Path(task["out_uniport"])
        gene_dir = out_uniport / g
        pep_csv = gene_dir / "pepmapviz_positions.csv"
        json_path = gene_dir / f"{g}_exons.json"

        # 1) 先看磁盘缓存
        if pep_csv.exists():
            df_cached = load_gene_match_result(g, out_uniport)
            if df_cached is not None and not df_cached.empty:
                return g, "cached"

        # 2) 无缓存需要 JSON
        if not json_path.exists():
            return g, "no_json"

        with open(json_path, "r", encoding="utf-8") as f:
            json_data = json.load(f)

        # 3) 跑 R
        calc = IsoformCoverageCalculator(task["r_script"], r_timeout=task["r_timeout"])
        df = calc.run_match_positions(
            report_path=task["report"],
            gene=g,
            json_data=json_data,
            isoform_select=task["isoform"],
            id_col=task["id_col"],
            label_col=task["label_col"],
            seq_col=task["seq_col"],
            genes_col=task["genes_col"],
            area_col=task["area_col"],
            normalize_IL=True,
        )
        if df is not None and not df.empty:
            save_gene_match_result(g, df, out_uniport)
            return g, "ok"
        else:
            return g, "empty"

    except Exception as e:
        print(f"[WARN] 子进程 PepMapViz 失败（{task.get('gene','?')}）：{e}")
        return task.get("gene", "?"), "error"


def _pepmapviz_thread_worker(args, genes: list[str], cache_csv: Path, done_event: threading.Event):
    """
    多进程版 PepMapViz（12 进程）：
      - 优先读取 out_uniport/{gene}/pepmapviz_positions.csv
      - 如无则调用 R（单基因），保存至 out_uniport/{gene}/pepmapviz_positions.csv
      - 保存成功或已缓存则把缓存表 PepMapVizDone=True
      - 最终把可用结果装入 RES_CACHE_BY_GENE
    """
    try:
        base_dir = Path(args.out_uniport)
        base_dir.mkdir(parents=True, exist_ok=True)

        # 读取/规范缓存表
        if cache_csv.exists():
            cache_df = pd.read_csv(cache_csv)
        else:
            cache_df = pd.DataFrame(columns=["Gene", "Protein", "Done", "PepMapVizDone"])
        if "PepMapVizDone" not in cache_df.columns:
            cache_df["PepMapVizDone"] = False

        # 预判哪些已经有 pep csv，可直接标记并跳过计算
        already = set()
        for g in genes:
            if (base_dir / g / "pepmapviz_positions.csv").exists():
                cache_df.loc[cache_df["Gene"].eq(g), "PepMapVizDone"] = True
                already.add(g)
        # 写回一次，以便主进程崩溃也保存状态
        cache_df.to_csv(cache_csv, index=False)

        # 待计算列表：没有 pep csv 且（有 JSON 或者 PepMapVizDone 为 False）
        to_run = []
        for g in genes:
            if g in already:
                continue
            json_path = base_dir / g / f"{g}_exons.json"
            if json_path.exists():
                to_run.append(g)

        # 准备任务包
        tasks = []
        for g in to_run:
            tasks.append({
                "gene": g,
                "out_uniport": str(base_dir),
                "r_script": args.r_script,
                "r_timeout": args.r_timeout,
                "report": args.report,
                "id_col": args.id_col,
                "label_col": args.label_col,
                "seq_col": args.seq_col,
                "genes_col": args.genes_col,
                "area_col": args.area_col,
                "isoform": args.isoform,
            })

        results = []
        per_gene_df = {}

        # 启动 12 个进程；限制 BLAS/OpenMP 线程，避免核上加核
        os.environ.setdefault("OMP_NUM_THREADS", "1")
        os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
        os.environ.setdefault("MKL_NUM_THREADS", "1")
        os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
        os.environ.setdefault("R_PARALLEL_NUM_THREADS", "1")

        print(f"[STEP] PepMapViz 多进程启动：{len(tasks)} 个任务，max_workers=12")
        with ProcessPoolExecutor(max_workers=int(args.workers/2)) as ex:
            futs = {ex.submit(_pepmapviz_run_one_gene, t): t["gene"] for t in tasks}
            for fut in tqdm(list(as_completed(futs)), total=len(futs), desc="PepMapViz 并行进度", dynamic_ncols=True):
                g, status = fut.result()
                results.append((g, status))

        # 汇总：把 ok/cached 写回缓存表，并装入内存缓存
        changed = False
        for g, status in results:
            if status in ("ok", "cached"):
                cache_df.loc[cache_df["Gene"].eq(g), "PepMapVizDone"] = True
                changed = True
                # 读取磁盘结果装入内存缓存
                df_cached = load_gene_match_result(g, base_dir)
                if df_cached is not None and not df_cached.empty:
                    per_gene_df[g] = df_cached

        if changed:
            cache_df.to_csv(cache_csv, index=False)

        if per_gene_df:
            # 装入全局缓存，供后续 match 使用
            IsoformCoverageCalculator.install_cache(per_gene_df)

        ok_n = sum(1 for _, s in results if s in ("ok", "cached")) + len(already)
        print(f"[OK] PepMapViz 线程完成：{ok_n}/{len(genes)} 个基因可用（含缓存 {len(already)}）。")

    except Exception as e:
        print(f"[WARN] PepMapViz 线程异常：{e}")
    finally:
        done_event.set()


def process_one_gene(gene: str, args, gene2prot) -> tuple[pd.DataFrame | None, dict | None]:
    """
    仅统计（不再调用 R）：
      - 读取 out_uniport/{gene}/pepmapviz_positions.csv
      - 读取 out_uniport/{gene}/{gene}_exons.json
      - 命中计数 value+=1；Precursor.Normalised 累加
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
        iso2uniq: dict[str, list[tuple[int, int]]] = {}
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
            print("[INFO] 无 unique-exon 区间，跳过命中统计")
            return None, stats

        # 5) 命中统计
        print("[STEP] 统计 unique-exon 命中 ...")
        rows = []
        for iso_acc, intervals in iso2uniq.items():
            sub_iso = sub.loc[sub["isoform"] == iso_acc].copy()
            if sub_iso.empty:
                continue

            area_col = args.area_col
            if area_col not in sub_iso.columns:
                sub_iso[area_col] = 0

            for _, r in tqdm(sub_iso.iterrows(), total=len(sub_iso), desc=f"{gene}:{iso_acc} 肽段命中", leave=False):
                s1 = int(r["start_position"]); s2 = int(r["end_position"])
                try:
                    area_val = float(r.get(area_col, 0)) if pd.notna(r.get(area_col, None)) else 0.0
                except Exception:
                    area_val = 0.0

                for a, b in intervals:
                    if not (s2 < a or b < s1):
                        rows.append({
                            "Polypeptide.Novogene.ID": r[args.id_col],
                            "gene": gene,
                            "canonical": canonical_acc,
                            "isoform": iso_acc,
                            "label": r.get(args.label_col, "All"),
                            "start": a,
                            "end": b,
                            "value": 1,
                            "Precursor.Normalised.sum": area_val,
                            "Global.Q.Value": r.get("Global.Q.Value", pd.NA)
                        })

        uniq_df = pd.DataFrame(rows) if rows else None
        return uniq_df, stats

    except Exception as e:
        print(f"[WARN] 任务异常（{gene}）：{e}")
        return None, None


def main():
    ap = argparse.ArgumentParser(description="计算覆盖 PSM（断点续跑 + 双缓存 + 并行版）")

    # 输入
    ap.add_argument("--report", default="/workspace/datasets/Match_result-X401SC24071912_Z01_F001_B1_43/report-4123692_8.tsv")
    ap.add_argument("--gene", default="ALL")
    ap.add_argument("--tax-id", default="9606")

    ap.add_argument("--isoform",   default="ALL")
    ap.add_argument("--id-col",    default="Polypeptide.Novogene.ID")
    ap.add_argument("--label-col", default="label")
    ap.add_argument("--seq-col",   default="Stripped.Sequence")
    ap.add_argument("--protein-col", default="Protein.Group")
    ap.add_argument("--genes-col", default="Genes")
    ap.add_argument("--area-col",  default="Precursor.Normalised")

    # 目录（PepMapViz 结果也写回 uniport 目录）
    ap.add_argument("--out-uniport",   default="/workspace/results/uniport")
    ap.add_argument("--out-match",     default="/workspace/results/match")
    ap.add_argument("--out-cache",     default="/workspace/results/caches")

    # R 调用配置（使用 psm.r）
    ap.add_argument("--r-script", default="src/psm.r")
    ap.add_argument("--r-timeout", type=int, default=900)

    # 运行控制
    ap.add_argument("--workers", type=int, default=256, help="并发线程数（用于match统计阶段）")

    args = ap.parse_args()
    start_time = time.time()

    # ==== 1) Gene→Protein 映射 ====
    gene2prot = _build_gene2prot_from_report(Path(args.report), args.genes_col, args.protein_col)

    # ==== 2) 缓存状态表（含 UniProtDone & PepMapVizDone） ====
    cache_csv = Path(args.out_cache) / "gene_protein_status.csv"
    cache_df = _ensure_cache_table(cache_csv, gene2prot)

    # ==== 3) 基因列表 ====
    genes = _load_genes_from_mapping(gene2prot, args.gene)
    print(f"[INFO] 示例前 10 个基因：{', '.join(genes[:10])}{' ...' if len(genes) > 10 else ''}")

    # ==== 4) UniProt 预取（断点续跑） ====
    _prefetch_gene_json_map(
        genes, gene2prot, args.tax_id, Path(args.out_uniport), cache_csv
    )

    # ==== 5) PepMapViz 线程（与 UniProt 独立；结果存回 uniport 目录；更新 PepMapVizDone） ====
    pep_done = threading.Event()
    t = threading.Thread(
        target=_pepmapviz_thread_worker,
        args=(args, genes, cache_csv, pep_done),
        daemon=True
    )
    t.start()

    # match 需要 pep 结果 + json，这里等待线程完成
    pep_done.wait()

    # ==== 6) 并发统计（只用缓存与 JSON，不再调用 R） ====
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

    # ==== 7) 汇总写盘（到 match 目录） ====
    print("[STEP] 写出结果 ...")
    out_dir = Path(args.out_match)
    out_dir.mkdir(parents=True, exist_ok=True)

    if uniq_psm_all_rows:
        uniq_all_df = pd.concat(uniq_psm_all_rows, ignore_index=True)

        # 明细表（不再使用 Sample.Name）
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

        # 汇总 1：按 id
        print("[STEP] 生成汇总（按 id） ...")
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

        # 汇总 2：按 exon（不区分 id）
        print("[STEP] 生成汇总（按 exon） ...")
        sum_by_exon = (
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
        out_sum_exon = out_dir / "unique_psm_sum_by_exon.csv"
        sum_by_exon.to_csv(out_sum_exon, index=False)
        print(f"[OK] 写出：{out_sum_exon}")
    else:
        print("[INFO] 没有任何命中行，未生成汇总文件。")

    print(f"[DONE] 全部完成，耗时 {time.time() - start_time:.1f} 秒")


if __name__ == "__main__":
    main()
