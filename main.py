#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
from pathlib import Path
import time

from src.pipeline import run_pipeline

def build_args():
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
    ap.add_argument("--workers", type=int, default=512, help="并发线程数（用于match统计阶段）")
    return ap.parse_args()

def main():
    args = build_args()
    t0 = time.time()
    run_pipeline(args)
    print(f"[DONE] 全部完成，耗时 {time.time() - t0:.1f} 秒")

if __name__ == "__main__":
    main()
