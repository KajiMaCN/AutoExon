# -*- coding: utf-8 -*-
from pathlib import Path
import json
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm.auto import tqdm

from src.exon_fetcher import GeneIsoformExonFetcher

def prefetch_gene_json_map(
    genes,
    gene2prot,
    tax_id: str,
    out_uniport: Path,
    cache_csv: Path,
) -> None:
    """
    并发预取 isoform/exon JSON 到 out_uniport/{gene}/{gene}_exons.json
    完成后将该 Gene 所有 pair 的 UniProtDone=True
    """
    print("[STEP] 并发预取 isoform/exon 元数据（UniProt，支持断点续跑） ...")

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

    ok_genes = cache_df.loc[cache_df["UniProtDone"], "Gene"].nunique()
    print(f"[INFO] UniProt 预取完成/已有：{ok_genes}/{len(genes)} 个基因。")
