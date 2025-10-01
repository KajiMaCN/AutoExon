import json
import requests
import pandas as pd

from collections import defaultdict
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple


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

            # 注意：EBI 返回的 exon 顺序常已按基因组排序；按需保留原顺序或轻微排序
            tx = {"ensembl_transcript_id": tx_id, "chromosome": chrom, "strand": strand, "exons": []}
            for i, ex in enumerate(sorted(exons, key=_exon_sort_key), start=1):
                # ---- 蛋白坐标：兼容只有 position 的情况 ----
                pl = self._pick(ex, "proteinLocation", "protein_location", "protLocation") or {}
                p_begin = self._to_int_maybe(self._pick(pl, "begin"))
                p_end   = self._to_int_maybe(self._pick(pl, "end"))
                if p_begin is None and p_end is None:
                    # 只有 position 时（可能是 {"position":{"position":123,"status":"certain"}}）
                    pos_one = self._to_int_maybe(self._pick(pl, "position"))
                    if pos_one is not None:
                        p_begin = pos_one
                        p_end   = pos_one
                # -------------------------------------------

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
    def _flatten_isoform_exons(
        result_dict: Dict[str, Any]
    ) -> Tuple[
        Dict[str, Dict[Tuple[str, int, int], Tuple[Optional[int], Optional[int]]]],
        List[Tuple[str, int, int]]
    ]:
        """
        收集 isoform → exons 映射，同时返回全体外显子坐标的排序列表。
        排序规则：
        - 先按 chrom 升序（字符串比较）
        - 同一 chrom 内按 start 升序
        - start 相同则按 end 升序
        """

        from collections import defaultdict as _dd

        exons_by_iso = _dd(dict)
        seen = set()

        def _to_int_safe(v):
            try:
                return int(v)
            except Exception:
                return None

        collected = []

        for iso in (result_dict.get("isoforms") or []):
            iso_acc = iso.get("accession")
            coords = iso.get("coordinates") or {}
            for tx in (coords.get("transcripts") or []):
                chrom = tx.get("chromosome")
                if chrom is None:
                    continue
                ch = str(chrom)
                for ex in (tx.get("exons") or []):
                    g = (ex.get("genomic") or {})
                    p = (ex.get("protein") or {})

                    s, e = g.get("start"), g.get("end")
                    if s is None or e is None:
                        continue
                    try:
                        s, e = int(s), int(e)
                    except Exception:
                        continue
                    if s > e:
                        s, e = e, s
                    key = (ch, s, e)

                    # protein 坐标
                    ps, pe = _to_int_safe(p.get("start")), _to_int_safe(p.get("end"))
                    if ps is None and pe is None:
                        pos_obj = p.get("position") or {}
                        pos = _to_int_safe(pos_obj.get("position"))
                        if pos is not None:
                            ps = pe = pos

                    if key not in exons_by_iso[iso_acc]:
                        exons_by_iso[iso_acc][key] = (ps, pe)
                    else:
                        old_ps, _ = exons_by_iso[iso_acc][key]
                        if old_ps is None and ps is not None:
                            exons_by_iso[iso_acc][key] = (ps, pe)

                    if key not in seen:
                        seen.add(key)
                        collected.append(key)

        # —— 按 (chrom, start, end) 排序 —— #
        ordered_coords = sorted(
            collected,
            key=lambda x: (x[0], x[1], x[2])  # chrom, start, end
        )

        return exons_by_iso, ordered_coords

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
        相邻夹心判断版：
        - 'single-isoform'：有效 isoform 仅 1 个
        - 'multi-isoform' ：该坐标属于多个 isoform
        - 'overlap_prev'  ：curr.start ∈ [prev.start, prev.end]
        - 'overlap_next'  ：curr.end   ∈ [next.start, next.end]
        - 'overlap_both'  ：两侧同时满足
        - 'unique'        ：仅属 1 个 isoform 且两侧均不满足
        - 'no-owner'
        """
        owners = self.build_unique_index(result_dict)  # coord -> set(isoforms)
        _, ordered_coords = self._flatten_isoform_exons(result_dict)  # 已保序 (chrom,start,end)

        # 统计有效 isoform
        effective_isoforms = set()
        for s in owners.values():
            effective_isoforms.update(s)
        if len(effective_isoforms) <= 1:
            return {coord: {'unique': False, 'reason': 'single-isoform'} for coord in ordered_coords}

        # 按染色体分组，保留组内顺序
        per_chrom: Dict[str, List[Tuple[str, int, int, int]]] = {}  # ch -> [(ch,s,e,global_idx), ...]
        for idx, (ch, s, e) in enumerate(ordered_coords):
            per_chrom.setdefault(ch, []).append((ch, s, e, idx))

        # 两侧“夹心”标志（按 global_idx 存）
        start_in_prev = [False] * len(ordered_coords)
        end_in_next   = [False] * len(ordered_coords)

        # 组内按 start,end 稳定排序（ordered_coords 已保序，保险起见再按数值排序）
        for ch, arr in per_chrom.items():
            arr.sort(key=lambda t: (t[1], t[2]))  # (ch,s,e,idx)
            n = len(arr)
            for j, (_, s, e, gi) in enumerate(arr):
                # 前一个
                if j - 1 >= 0:
                    _, ps, pe, _ = arr[j-1]
                    if ps is not None and pe is not None:
                        start_in_prev[gi] = (ps <= s <= pe)
                # 后一个
                if j + 1 < n:
                    _, ns, ne, _ = arr[j+1]
                    if ns is not None and ne is not None:
                        end_in_next[gi] = (ns <= e <= ne)

        info_map: Dict[Tuple[str,int,int], Dict[str, Any]] = {}
        for idx, coord in enumerate(ordered_coords):
            own = owners.get(coord, set())
            if not own:
                info_map[coord] = {'unique': False, 'reason': 'no-owner'}
                continue
            if len(own) != 1:
                info_map[coord] = {'unique': False, 'reason': 'multi-isoform'}
                continue

            has_prev = bool(start_in_prev[idx])
            has_next = bool(end_in_next[idx])

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
        写 {out_prefix}_exon_matrix.csv

        完全保序（不排序）版本：
        - 行顺序：严格使用 _flatten_isoform_exons 返回的 ordered_coords（该函数已按 JSON 遍历顺序去重保序）
        - 列顺序：严格按 result_dict["isoforms"] 出现顺序；仅删除“整列均为 '-'”的 isoform 列
        - 单元格缺失仍用 '-' 占位
        - unique / overlap_reason 来自 compute_unique_and_reason_map，但只做标注，不影响顺序
        """
        # 1) 取 isoform 的原始顺序
        isoforms_in_order = [
            iso.get("accession") for iso in (result_dict.get("isoforms") or [])
            if iso.get("accession")
        ]

        # 2) 底层映射与“按 JSON 顺序的坐标列表”
        exons_by_iso, ordered_coords = self._flatten_isoform_exons(result_dict)

        # 3) 唯一性/重叠原因（仅标注）
        uniq_info = self.compute_unique_and_reason_map(result_dict)

        # 4) 按 ordered_coords 原样构造行
        matrix_rows = []
        for chrom, s, e in ordered_coords:
            info = uniq_info.get((chrom, s, e), {'unique': False, 'reason': 'no-owner'})
            is_unique = bool(info.get('unique', False))
            reason = str(info.get('reason', ''))
            overlap_conflict = (not is_unique) and reason.startswith('overlap')

            # 如果 unique==True，找该坐标归属的 isoform（按原始列序第一个命中的）
            unique_iso = None
            if is_unique:
                for iso_acc in isoforms_in_order:
                    if (chrom, s, e) in exons_by_iso.get(iso_acc, {}):
                        unique_iso = iso_acc
                        break

            row = {
                "chrom": chrom,
                "start": s,
                "end": e,
                "unique": is_unique,
                "unique_isoform": unique_iso,
                "overlap_conflict": overlap_conflict,
                "overlap_reason": reason if (overlap_conflict or not is_unique) else "unique"
            }

            # 填每个 isoform 的蛋白坐标；缺失用 '-'
            for iso_acc in isoforms_in_order:
                ps, pe = exons_by_iso.get(iso_acc, {}).get((chrom, s, e), (None, None))
                row[iso_acc] = (f"{ps} - {pe}" if (ps is not None and pe is not None) else "-")

            matrix_rows.append(row)

        wide_df = pd.DataFrame(matrix_rows)

        # 5) 删除整列均为 '-' 的 isoform 列（不改变其余列顺序）
        base_cols = ["chrom", "start", "end", "unique", "unique_isoform", "overlap_conflict", "overlap_reason"]
        keep_iso_cols = []
        for iso_acc in isoforms_in_order:
            if iso_acc in wide_df.columns:
                col = wide_df[iso_acc]
                if not (col == "-").all():
                    keep_iso_cols.append(iso_acc)
                else:
                    wide_df.drop(columns=[iso_acc], inplace=True)

        # 6) 组装最终列顺序（严格保序）
        final_cols = base_cols + keep_iso_cols
        wide_df = wide_df[final_cols]

        # 7) 直接写盘（不再 sort）
        out_prefix = Path(out_prefix)
        out_prefix.parent.mkdir(parents=True, exist_ok=True)
        out_csv = f"{out_prefix}_exon_matrix.csv"
        wide_df.to_csv(out_csv, index=False)
        print(f"[OK] saved: {out_csv}")
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
