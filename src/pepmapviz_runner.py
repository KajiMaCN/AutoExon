# -*- coding: utf-8 -*-
import os
import json
import threading
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm.auto import tqdm
import pandas as pd

from src.psm import (
    IsoformCoverageCalculator, RES_CACHE_BY_GENE,
    save_gene_match_result, load_gene_match_result
)

def _run_one_gene(task):
    """
    子进程执行：单基因 PepMapViz
    返回 (gene, status)：status ∈ {"ok","cached","empty","no_json","error"}
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


def pepmapviz_thread_worker(args, genes, cache_csv: Path, done_event: threading.Event):
    """
    多进程版 PepMapViz：
      - 优先读取 out_uniport/{gene}/pepmapviz_positions.csv
      - 如无则调用 R（单基因），保存至 out_uniport/{gene}/pepmapviz_positions.csv
      - 保存成功或已缓存则把缓存表 PepMapVizDone=True
      - 最终把可用结果装入 RES_CACHE_BY_GENE
    """
    try:
        base_dir = Path(args.out_uniport)
        base_dir.mkdir(parents=True, exist_ok=True)

        if cache_csv.exists():
            cache_df = pd.read_csv(cache_csv)
        else:
            cache_df = pd.DataFrame(columns=["Gene", "Protein", "Done", "PepMapVizDone"])
        if "PepMapVizDone" not in cache_df.columns:
            cache_df["PepMapVizDone"] = False

        # 已有 pep csv 的直接标记
        already = set()
        for g in genes:
            if (base_dir / g / "pepmapviz_positions.csv").exists():
                cache_df.loc[cache_df["Gene"].eq(g), "PepMapVizDone"] = True
                already.add(g)
        cache_df.to_csv(cache_csv, index=False)

        # 待跑任务（有 JSON 且还没 pep）
        to_run = []
        for g in genes:
            if g in already:
                continue
            json_path = base_dir / g / f"{g}_exons.json"
            if json_path.exists():
                to_run.append(g)

        # 任务包
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

        # 限制 BLAS/OpenMP 线程，避免核上加核
        os.environ.setdefault("OMP_NUM_THREADS", "1")
        os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
        os.environ.setdefault("MKL_NUM_THREADS", "1")
        os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
        os.environ.setdefault("R_PARALLEL_NUM_THREADS", "1")

        max_workers = max(1, int(args.workers // 2))
        print(f"[STEP] PepMapViz 多进程启动：{len(tasks)} 个任务，max_workers={max_workers}")
        with ProcessPoolExecutor(max_workers=max_workers) as ex:
            futs = {ex.submit(_run_one_gene, t): t["gene"] for t in tasks}
            for fut in tqdm(list(as_completed(futs)), total=len(futs), desc="PepMapViz 并行进度", dynamic_ncols=True):
                g, status = fut.result()
                results.append((g, status))

        # 汇总：把 ok/cached 写回缓存表，并装入内存缓存
        changed = False
        for g, status in results:
            if status in ("ok", "cached"):
                cache_df.loc[cache_df["Gene"].eq(g), "PepMapVizDone"] = True
                changed = True
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
