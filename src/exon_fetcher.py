# -*- coding: utf-8 -*-
import json
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple
from collections import defaultdict

import requests
import pandas as pd


class GeneIsoformExonFetcher:
    """
    统一封装：
      - （改）直接根据传入的 UniProt canonical accession 进行后续步骤
      - 获取 isoform 列表
      - 通过 EBI Proteins coordinates 获取转录本/外显子坐标（基因组 & 蛋白）
      - 构建“外显子 x isoform”的矩阵，判定 unique 外显子
      - 保存矩阵 CSV、以及精简 JSON
    """

    UNI_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
    UNI_ENTRY  = "https://rest.uniprot.org/uniprotkb/{acc}.json"
    EBI_COORDS = "https://www.ebi.ac.uk/proteins/api/coordinates"

    def __init__(self, timeout: int = 30, session: Optional[requests.Session] = None):
        """
        :param timeout: HTTP 请求超时（秒）
        :param session: 可注入的 requests.Session（用于连接复用/重试策略）
        """
        self.timeout = timeout
        self.session = session or requests.Session()
        self._default_headers = {"Accept": "application/json"}

    # -------------------- 基础工具方法 --------------------

    def _get_json(self, url: str, params: Optional[dict] = None, headers: Optional[dict] = None) -> Any:
        h = dict(self._default_headers)
        if headers:
            h.update(headers)
        r = self.session.get(url, params=params, headers=h, timeout=self.timeout)
        r.raise_for_status()
        return r.json()

    @staticmethod
    def _to_int_maybe(x):
        """把可能的字符串/浮点/字典形式位置统一转 int。"""
        if x is None:
            return None
        if isinstance(x, dict):
            # 兼容 {"position": 123}, {"value": 123}, {"pos": 123}, {"start":..}, {"end":..}
            for k in ("position", "value", "pos", "start", "end"):
                if k in x:
                    try:
                        return int(x[k])
                    except Exception:
                        pass
            return None
        try:
            return int(x)
        except Exception:
            try:
                return int(str(x).strip())
            except Exception:
                return None

    @staticmethod
    def _pick(obj: dict, *candidates):
        """按候选键依次取第一个存在的子对象。"""
        for k in candidates:
            v = obj.get(k)
            if v is not None:
                return v
        return None

    @staticmethod
    def _fmt_thousand(n):
        return f"{n:,}" if isinstance(n, int) else "NA"

    # -------------------- UniProt / EBI 查询 --------------------

    # （保留，但当前流程不再使用基因名查询）
    def uniprot_get_canonical_accession_by_gene(self, gene_symbol: str, tax_id: str = "9606") -> str:
        q = f"gene_exact:{gene_symbol} AND organism_id:{tax_id} AND reviewed:true"
        params = {
            "query": q,
            "fields": "accession,organism_name,reviewed",
            "format": "json",
            "size": 1
        }
        data = self._get_json(self.UNI_SEARCH, params=params)
        results = data.get("results", [])
        if not results:
            raise ValueError(f"未找到基因 {gene_symbol} (tax_id={tax_id}) 的 reviewed 条目。")
        return results[0]["primaryAccession"]  # canonical accession

    def uniprot_get_isoform_accessions(self, canonical_acc: str) -> List[str]:
        """
        从 UniProt 主条目 JSON 的 'ALTERNATIVE PRODUCTS' 注释里抽取 isoform ID，
        例如返回 ['P19827-1', 'P19827-3']；若无同种型，至少返回 canonical 本体。
        """
        entry = self._get_json(self.UNI_ENTRY.format(acc=canonical_acc))
        isoforms: List[str] = []

        for c in entry.get("comments", []):
            if c.get("commentType") == "ALTERNATIVE PRODUCTS":
                for iso in c.get("isoforms", []):
                    ids = iso.get("isoformIds") or []
                    for i in ids:
                        if isinstance(i, str):
                            isoforms.append(i)

        isoforms = sorted(set(isoforms))
        if not isoforms:
            isoforms = [canonical_acc]
        return isoforms

    def ebi_get_coordinates_for_accession(self, iso_acc: str, tax_id: Optional[str] = None) -> Dict[str, Any]:
        """
        调 EBI Proteins coordinates：返回该 isoform 的所有转录本及其外显子（左基因组坐标 / 右蛋白坐标）。
        """
        params = {"accession": iso_acc}
        if tax_id:
            params["taxid"] = tax_id
        data = self._get_json(self.EBI_COORDS, params=params)
        if not isinstance(data, list) or not data:
            return {"accession": iso_acc, "sequence": None, "transcripts": []}

        entry = data[0]
        res: Dict[str, Any] = {
            "accession": entry.get("accession"),
            "sequence": entry.get("sequence"),
            "transcripts": []
        }

        gn_coords = entry.get("gnCoordinate") or []
        if not gn_coords:
            return res

        def _exon_sort_key(ex):
            gl = self._pick(ex, "genomicLocation", "genomeLocation", "location") or {}
            b = self._pick(gl, "begin", "start")
            return self._to_int_maybe(b) or 0

        for trans in gn_coords:
            tx_id = trans.get("ensemblTranscriptId")
            gloc = trans.get("genomicLocation") or {}
            chrom = gloc.get("chromosome")
            strand_val = gloc.get("strand")
            strand = strand_val if isinstance(strand_val, str) else ("+" if strand_val in (1, "+") else "-")
            exons = gloc.get("exon") or []

            tx = {"ensembl_transcript_id": tx_id, "chromosome": chrom, "strand": strand, "exons": []}
            for i, ex in enumerate(sorted(exons, key=_exon_sort_key), start=1):
                # 蛋白坐标
                pl = self._pick(ex, "proteinLocation", "protein_location", "protLocation") or {}
                p_begin = self._to_int_maybe(self._pick(pl, "begin"))
                p_end   = self._to_int_maybe(self._pick(pl, "end"))

                # 基因组坐标
                gl = self._pick(ex, "genomicLocation", "genomeLocation", "location") or {}
                g_begin = self._to_int_maybe(self._pick(gl, "begin", "start"))
                g_end   = self._to_int_maybe(self._pick(gl, "end", "stop"))

                def fmt_num(n):
                    return f"{n:,}" if isinstance(n, int) else "NA"

                genomic_label = (
                    f"{chrom}:{fmt_num(g_begin)} - {fmt_num(g_end)}" if chrom
                    else f"{fmt_num(g_begin)} - {fmt_num(g_end)}"
                )

                tx["exons"].append({
                    "index": i,
                    "exon_id": ex.get("id"),
                    "genomic": {"chromosome": chrom, "start": g_begin, "end": g_end, "strand": strand},
                    "genomic_label": genomic_label,
                    "protein": {"start": p_begin, "end": p_end},
                })
            res["transcripts"].append(tx)
        return res

    # -------------------- 汇总 / 扁平化 / 构表 --------------------

    def query_gene_isoforms_and_exons(self, gene_symbol: str, tax_id: str = "9606") -> Dict[str, Any]:
        """
        （旧接口，保留）基因名 -> canonical accession -> isoforms -> 每个 isoform 的外显子(左+右坐标)。
        """
        canonical = self.uniprot_get_canonical_accession_by_gene(gene_symbol, tax_id)
        return self.query_by_uniprot_accession(canonical_acc=canonical, tax_id=tax_id, gene_symbol=gene_symbol)
    
    def query_by_uniprot_accession(self, canonical_acc: str, tax_id: str = "9606",
                                   gene_symbol: Optional[str] = None) -> Dict[str, Any]:
        """
        新接口：直接从 canonical UniProt accession 出发，获取 isoforms 与外显子坐标。
        :param canonical_acc: 例如 'P02768'
        :param tax_id: 物种（默认人）
        :param gene_symbol: 可选，仅用于输出里给个提示信息；不影响查询逻辑
        """
        iso_accs = self.uniprot_get_isoform_accessions(canonical_acc)
        out = {
            "gene": gene_symbol,
            "tax_id": tax_id,
            "canonical": canonical_acc,
            "isoforms": []
        }
        for iso in iso_accs:
            coords = self.ebi_get_coordinates_for_accession(iso, tax_id)
            out["isoforms"].append({"accession": iso, "coordinates": coords})
        return out

    @staticmethod
    def _flatten_isoform_exons(result_dict: Dict[str, Any]) -> Tuple[Dict[str, Dict[Tuple[str,int,int], Tuple[Optional[int], Optional[int]]]], List[Tuple[str,int,int]]]:
        """
        拍平成:
          exons_by_iso: dict[iso] -> dict[(chrom,start,end)] -> (p_start, p_end)
          all_coords:   list[(chrom,start,end)]  按 chr,start,end 排序
        """
        exons_by_iso = defaultdict(dict)
        all_coords_set = set()

        def _to_int_safe(v):
            try:
                return int(v)
            except Exception:
                return None

        for iso in result_dict.get("isoforms", []):
            iso_acc = iso["accession"]
            coords = iso["coordinates"]
            for tx in coords.get("transcripts", []):
                chrom = tx.get("chromosome")
                for ex in tx.get("exons", []):
                    g = (ex.get("genomic") or {})
                    p = (ex.get("protein") or {})
                    s, e = g.get("start"), g.get("end")
                    if chrom is None or s is None or e is None:
                        continue
                    s, e = int(s), int(e)
                    if s > e:
                        s, e = e, s
                    key = (str(chrom), s, e)
                    all_coords_set.add(key)

                    ps, pe = _to_int_safe(p.get("start")), _to_int_safe(p.get("end"))
                    # 如果同一 iso+坐标出现多次，只保留第一组非空蛋白坐标
                    if key not in exons_by_iso[iso_acc]:
                        exons_by_iso[iso_acc][key] = (ps, pe)
                    else:
                        old_ps, old_pe = exons_by_iso[iso_acc][key]
                        if old_ps is None and ps is not None:
                            exons_by_iso[iso_acc][key] = (ps, pe)

        all_coords = sorted(all_coords_set, key=lambda t: (t[0], t[1], t[2]))
        return exons_by_iso, all_coords

    def build_matrix_and_unique(self, result_dict: Dict[str, Any], canonical_first: bool = True):
        exons_by_iso, all_coords = self._flatten_isoform_exons(result_dict)
        isoforms = [iso["accession"] for iso in result_dict.get("isoforms", [])]

        if canonical_first and result_dict.get("canonical") in isoforms:
            can = result_dict["canonical"]
            isoforms = sorted(isoforms, key=lambda x: (0 if x == can or x.startswith(can + "-1") else 1, x))

        # 用“新规则”计算每行是否 unique
        unique_map = self.compute_unique_map(result_dict)

        header = ["EXON COORDINATES"] + isoforms
        rows = []
        matrix_rows = []

        for chrom, s, e in all_coords:
            is_unique = bool(unique_map.get((chrom, s, e), False))
            left_label = f"{chrom}:{self._fmt_thousand(s)} - {self._fmt_thousand(e)}"
            if is_unique:
                left_label += " *"

            row = [left_label]
            row_struct = {"chrom": chrom, "start": s, "end": e, "unique": is_unique}
            for iso in isoforms:
                ps, pe = exons_by_iso.get(iso, {}).get((chrom, s, e), (None, None))
                cell = f"{ps} - {pe}" if (ps is not None and pe is not None) else "-"
                row.append(cell)
                row_struct[iso] = cell
            rows.append(row)
            matrix_rows.append(row_struct)

        col_widths = [max(len(h), *(len(r[i]) for r in rows)) for i, h in enumerate(header)]
        return header, rows, col_widths, matrix_rows

    def compute_unique_and_reason_map(self, result_dict: Dict[str, Any]) -> Dict[Tuple[str,int,int], Dict[str, Any]]:
        """
        返回：coord -> {'unique': bool, 'reason': str}
        reason 取值：
        - 'single-isoform' : 有效 isoform 仅 1 个（仅该 isoform 真的提供了外显子坐标）
        - 'multi-isoform'  : 该坐标出现在多个 isoform
        - 'overlap_prev'   : 与上一相邻区间重叠
        - 'overlap_next'   : 与下一相邻区间重叠
        - 'overlap_both'   : 同时与前/后相邻区间重叠
        - 'unique'         : 满足唯一且与相邻区间均不重叠
        - 'no-owner'       : 少见（无归属）
        """
        owners = self.build_unique_index(result_dict)  # coord -> set(isoforms)
        _, all_coords = self._flatten_isoform_exons(result_dict)

        # 关键修正：以 owners 的并集，计算“有效 isoform”数量（确实有外显子坐标参与的 isoform）
        effective_isoforms = set()
        for s in owners.values():
            effective_isoforms.update(s)

        # 若有效 isoform <= 1，则全部标记为 single-isoform（不判定 unique）
        if len(effective_isoforms) <= 1:
            return {coord: {'unique': False, 'reason': 'single-isoform'} for coord in all_coords}

        def overlap(a, b) -> bool:
            # 仅同染色体才判断重叠；闭区间重叠判定
            return a[0] == b[0] and max(a[1], b[1]) <= min(a[2], b[2])

        info_map: Dict[Tuple[str,int,int], Dict[str, Any]] = {}
        n = len(all_coords)
        for i, coord in enumerate(all_coords):
            own = owners.get(coord, set())
            if not own:
                info_map[coord] = {'unique': False, 'reason': 'no-owner'}
                continue
            if len(own) != 1:
                info_map[coord] = {'unique': False, 'reason': 'multi-isoform'}
                continue

            # 候选唯一：检查与相邻区间是否重叠
            has_prev = (i - 1 >= 0 and overlap(all_coords[i-1], coord))
            has_next = (i + 1 < n  and overlap(coord, all_coords[i+1]))

            if has_prev and has_next:
                info_map[coord] = {'unique': False, 'reason': 'overlap_both'}
            elif has_prev:
                info_map[coord] = {'unique': False, 'reason': 'overlap_prev'}
            elif has_next:
                info_map[coord] = {'unique': False, 'reason': 'overlap_next'}
            else:
                info_map[coord] = {'unique': True,  'reason': 'unique'}

        return info_map

    @staticmethod
    def print_exon_matrix(header, rows, col_widths):
        """按列宽打印成和截图类似的表格；缺失用 '-'。"""
        def _fmt_row(cols):
            return "  ".join(c.ljust(w) for c, w in zip(cols, col_widths))
        print(_fmt_row(header))
        print(_fmt_row(["-" * w for w in col_widths]))
        for r in rows:
            print(_fmt_row(r))

    # -------------------- 保存：CSV / JSON --------------------

    def build_tables_for_saving(self, result_dict: Dict[str, Any], out_prefix: str, canonical_first: bool = True) -> pd.DataFrame:
        """
        生成并保存：
        {out_prefix}_exon_matrix.csv —— 宽表（矩阵）

        规则：
        - unique 判定使用 compute_unique_and_reason_map
        - 缺失用 "-" 占位
        - 如果 isoform 列整列都是 "-"，则删除该列
        - 列排序：canonical 在前，其后同一前缀的 -k 按数字升序
        - 行排序：若主链为 '+' 则 (chrom, start) 升序；若为 '-' 则 (chrom 升序, start 降序)
        """
        from collections import Counter

        # --- 聚合底层结构 ---
        exons_by_iso, all_coords = self._flatten_isoform_exons(result_dict)
        isoforms_all = [iso["accession"] for iso in result_dict.get("isoforms", [])]

        # canonical 初步靠前
        can = result_dict.get("canonical")
        if canonical_first and can in isoforms_all:
            isoforms_all = sorted(isoforms_all, key=lambda x: (0 if (x == can or x.startswith(can + "-1")) else 1, x))
        else:
            isoforms_all = sorted(isoforms_all)

        # 唯一性与原因
        uniq_info = self.compute_unique_and_reason_map(result_dict)

        # 主链方向（统计 transcripts 的 strand 众数）
        strands = []
        for iso in result_dict.get("isoforms", []):
            coords = iso.get("coordinates") or {}
            for tx in coords.get("transcripts", []):
                s = tx.get("strand")
                if s in ("+", "-"):
                    strands.append(s)
        main_strand = Counter(strands).most_common(1)[0][0] if strands else "+"

        # --- 行构造 ---
        matrix_rows = []
        for chrom, s, e in all_coords:
            info = uniq_info.get((chrom, s, e), {'unique': False, 'reason': 'no-owner'})
            is_unique = bool(info['unique'])
            reason = str(info['reason'])

            overlap_conflict = (not is_unique) and reason.startswith('overlap')

            unique_iso = None
            if is_unique:
                for iso in isoforms_all:
                    if (chrom, s, e) in exons_by_iso.get(iso, {}):
                        unique_iso = iso
                        break

            row = {
                "chrom": chrom,
                "start": s,
                "end": e,
                "unique": is_unique,
                "unique_isoform": unique_iso,
                "overlap_conflict": overlap_conflict,
                "overlap_reason": (reason if overlap_conflict else reason)
            }

            for iso in isoforms_all:
                ps, pe = exons_by_iso.get(iso, {}).get((chrom, s, e), (None, None))
                row[iso] = (f"{ps} - {pe}" if (ps is not None and pe is not None) else "-")

            matrix_rows.append(row)

        wide_df = pd.DataFrame(matrix_rows)

        # 删除整列均为 "-" 的 isoform 列
        for iso in isoforms_all:
            if iso in wide_df.columns and all(wide_df[iso] == "-"):
                wide_df.drop(columns=[iso], inplace=True)

        # 列顺序：基础列 + isoform 列（排序）
        def isoform_sort_key(acc: str):
            if "-" not in acc:
                return (acc, 0, 0)
            prefix, suffix = acc.split("-", 1)
            num_str = "".join(ch for ch in suffix if ch.isdigit())
            try:
                num = int(num_str) if num_str else float("inf")
            except Exception:
                num = float("inf")
            return (prefix, 1, num)

        base_cols = ["chrom", "start", "end", "unique", "unique_isoform", "overlap_conflict", "overlap_reason"]
        iso_cols_sorted = sorted([c for c in wide_df.columns if c not in base_cols], key=isoform_sort_key)
        wide_df = wide_df[base_cols + iso_cols_sorted]

        # 行排序
        def chrom_rank(ch):
            s = str(ch).strip()
            s = s[3:] if s.lower().startswith("chr") else s
            u = s.upper()
            specials = {"X": 23, "Y": 24, "MT": 25, "M": 25}
            if u in specials:
                return (0, specials[u])
            try:
                return (0, int(u))
            except Exception:
                return (1, u)

        asc_start = (main_strand == "+")
        wide_df = wide_df.sort_values(
            by=["chrom", "start"],
            ascending=[True, asc_start],
            key=lambda col: col.map(chrom_rank) if col.name == "chrom" else col
        ).reset_index(drop=True)

        # --- 保存 ---
        out_prefix = Path(out_prefix)
        out_prefix.parent.mkdir(parents=True, exist_ok=True)
        wide_csv = f"{out_prefix}_exon_matrix.csv"
        wide_df.to_csv(wide_csv, index=False)
        print(f"[OK] saved: {wide_csv}")
        return wide_df

    @staticmethod
    def _coord_key_from_exon(ex, default_chrom=None):
        g = (ex.get("genomic") or {})
        chrom = g.get("chromosome", default_chrom)
        s = g.get("start"); e = g.get("end")
        if chrom is None or s is None or e is None:
            return None
        s, e = int(s), int(e)
        if s > e: s, e = e, s
        return (str(chrom), s, e)

    def build_unique_index(self, result_dict: Dict[str, Any]):
        """统计每个 (chrom,start,end) 出现在哪些 isoform，用于判定 unique。"""
        owners = defaultdict(set)  # coord -> set(iso)
        for iso in result_dict.get("isoforms", []):
            iso_acc = iso["accession"]
            coords = iso["coordinates"]
            for tx in coords.get("transcripts", []):
                chrom = tx.get("chromosome")
                for ex in tx.get("exons", []):
                    key = self._coord_key_from_exon(ex, default_chrom=chrom)
                    if key is not None:
                        owners[key].add(iso_acc)
        return owners

    def build_json_dict(self, result_dict: Dict[str, Any]) -> Dict[str, Any]:
        """
        构造 JSON：
        - 每个 isoform.exons[i] 含：exon_id, start, end, unique, reason, exon(序号)
        - unique / reason 与 CSV 完全一致（来自 compute_unique_and_reason_map）
        """
        owners_index = self.build_unique_index(result_dict)          # coord -> set(iso)
        uniq_info    = self.compute_unique_and_reason_map(result_dict)  # coord -> {'unique':bool,'reason':str}

        out = {
            "gene": result_dict.get("gene"),
            "tax_id": result_dict.get("tax_id"),
            "canonical": result_dict.get("canonical"),
            "isoforms": []
        }

        for iso in result_dict.get("isoforms", []):
            iso_acc = iso.get("accession")
            seq = (iso.get("coordinates") or {}).get("sequence")
            items = []
            seen_coords = set()

            for tx in (iso.get("coordinates") or {}).get("transcripts", []):
                chrom_tx = tx.get("chromosome")
                for ex in tx.get("exons", []):
                    g = (ex.get("genomic") or {})
                    p = (ex.get("protein") or {})
                    chrom = g.get("chromosome", chrom_tx)
                    s = g.get("start"); e = g.get("end")
                    ps = p.get("start"); pe = p.get("end")
                    if chrom is None or s is None or e is None or ps is None or pe is None:
                        continue
                    s, e = int(s), int(e)
                    if s > e: s, e = e, s
                    key = (str(chrom), s, e)
                    if key in seen_coords:
                        continue
                    seen_coords.add(key)

                    # 正确读取唯一性与原因
                    info = uniq_info.get(key, {'unique': False, 'reason': 'no-owner'})
                    owners = owners_index.get(key, set())
                    coord_unique = bool(info.get('unique', False)) and (iso_acc in owners)
                    reason = str(info.get('reason', 'no-owner'))

                    items.append({
                        "exon_id": ex.get("exon_id"),
                        "start": int(ps),
                        "end": int(pe),
                        "unique": coord_unique,
                        "reason": reason
                    })

            # 以蛋白坐标升序
            items.sort(key=lambda d: (d["start"], d["end"]))
            for i, it in enumerate(items, start=1):
                it["exon"] = i

            out["isoforms"].append({
                "accession": iso_acc,
                "sequence": seq,
                "exons": items
            })
        return out

    def save_json(self, result_dict: Dict[str, Any], out_json_path: Optional[str]):
        """
        写文件（若提供路径）并返回构造的字典。
        :param out_json_path: 文件路径；传 None 则不写盘，仅返回字典。
        :return: 构造好的 JSON 字典（可直接用于后续计算）
        """
        out = self.build_json_dict(result_dict)

        if out_json_path:
            outp = Path(out_json_path)
            outp.parent.mkdir(parents=True, exist_ok=True)
            with open(outp, "w", encoding="utf-8") as f:
                json.dump(out, f, ensure_ascii=False, indent=2)
            print(f"[OK] saved: {out_json_path}")

        return out


# -------------------- 示例用法 --------------------
if __name__ == "__main__":
    fetcher = GeneIsoformExonFetcher(timeout=30)

    # 直接使用 UniProt canonical accession（例如 ALB 的 canonical: P02768）
    CANONICAL_ACC = "P02768"
    TAX_ID = "9606"

    # （可选）仅用于输出展示用，不参与查询
    GENE_SYMBOL = "ALB"

    out_dir = Path("/workspace/results") / (GENE_SYMBOL or CANONICAL_ACC)
    out_dir.mkdir(parents=True, exist_ok=True)

    # 新入口：由 canonical accession 出发
    data = fetcher.query_by_uniprot_accession(CANONICAL_ACC, TAX_ID, gene_symbol=GENE_SYMBOL)

    # 打印矩阵
    header, rows, col_widths, _ = fetcher.build_matrix_and_unique(data, canonical_first=True)
    print(f"Gene: {data.get('gene')}  TaxID: {data['tax_id']}  Canonical: {data['canonical']}")
    fetcher.print_exon_matrix(header, rows, col_widths)

    # 保存矩阵 CSV
    out_prefix = str(out_dir / (GENE_SYMBOL or CANONICAL_ACC))
    fetcher.build_tables_for_saving(data, out_prefix=out_prefix)

    # 保存“中等简易版”JSON
    fetcher.save_json(data, str(out_dir / f"{GENE_SYMBOL or CANONICAL_ACC}_exons.json"))
