import os
import json
import subprocess
import pandas as pd

from pathlib import Path
from typing import Dict, List, Optional
from threading import Lock
from tqdm.auto import tqdm

# ====== 全局缓存：一次 R 返回，所有线程复用 ======
RES_CACHE_BY_GENE: dict[str, pd.DataFrame] = {}
_RES_CACHE_LOCK = Lock()

# 需要保留精度的列名（统一在这里配置）
_PRESERVE_COLS = ("Global.Q.Value", "Precursor.Normalised")

def _cast_preserve_cols_to_str(df: pd.DataFrame) -> pd.DataFrame:
    """将需要保留位数的列转为字符串（发送给 R 前调用）。"""
    df = df.copy()
    for col in _PRESERVE_COLS:
        if col in df.columns:
            df[col] = df[col].astype(str)
    return df

def _attach_numeric_clones(df: pd.DataFrame) -> pd.DataFrame:
    """
    为保留位数的字符串列增加数值副本（_num 后缀），用于后续计算；
    原始列保持字符串，不改变写盘展示。
    """
    df = df.copy()
    for col in _PRESERVE_COLS:
        if col in df.columns:
            df[f"{col}_num"] = pd.to_numeric(df[col], errors="coerce")
    return df


def install_cache(per_gene_df: Dict[str, pd.DataFrame]):
    """把已有的每基因 DataFrame 装入全局缓存。"""
    with _RES_CACHE_LOCK:
        RES_CACHE_BY_GENE.update({k: v for k, v in per_gene_df.items() if v is not None and not v.empty})


def save_gene_match_result(gene: str, df: pd.DataFrame, out_uniport_dir: Path) -> Path:
    """
    把单基因的 PepMapViz 匹配结果保存到 uniport/{gene}/pepmapviz_positions.csv
    保留小数精度：不做格式化，让 pandas 原样写出（不省略）。
    """
    out_dir = out_uniport_dir / gene
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "pepmapviz_positions.csv"
    # 不使用 float_format，避免改变小数位（按用户要求不省略、不改动位数）
    df.to_csv(out_file, index=False)
    return out_file


def load_gene_match_result(gene: str, out_uniport_dir: Path) -> Optional[pd.DataFrame]:
    """
    从 uniport/{gene}/pepmapviz_positions.csv 读取单基因缓存；不存在则返回 None
    """
    f = out_uniport_dir / gene / "pepmapviz_positions.csv"
    if f.exists():
        try:
            df = pd.read_csv(f)
            return df
        except Exception:
            return None
    return None


def _apply_report_filters(rep: pd.DataFrame, *, genes_col: str, protein_col: str) -> pd.DataFrame:
    """
    统一对 report 执行两条过滤：
      1) 去除 Protein.Group 或 Genes 中包含分号 ';' 的行
      2) 仅保留 Global.Q.Value < 0.01 的行（非数值会被置为 NaN 并剔除）
    """
    df = rep.copy()
    if genes_col in df.columns:
        df = df[~df[genes_col].astype(str).str.contains(";", na=False)]
    if protein_col in df.columns:
        df = df[~df[protein_col].astype(str).str.contains(";", na=False)]
    qcol = "Global.Q.Value"
    if qcol in df.columns:
        q = pd.to_numeric(df[qcol], errors="coerce")
        df = df[q < 0.01]
    return df


class IsoformCoverageCalculator:
    """
    PepMapViz（R）匹配：
    - run_match_positions：逐 gene 逐 isoform（与 psm.r 单基因接口匹配）
    - 另保留批量接口（如需后续扩展）
    """

    def __init__(self, pepmapviz_script: str, rscript_bin: str = "Rscript", r_timeout: int = 900):
        self.pepmapviz_script = pepmapviz_script
        self.rscript_bin = rscript_bin
        self.r_timeout = r_timeout

    # ---------- 工具 ----------
    @staticmethod
    def read_table_auto(path: str, usecols: Optional[List[str]] = None, dtype=None) -> pd.DataFrame:
        sep = "\t" if str(path).endswith((".tsv", ".tsv.gz")) else ","
        header = pd.read_csv(path, sep=sep, nrows=0, low_memory=False)
        cols = header.columns.tolist()
        if usecols is not None:
            usecols = [c for c in usecols if c in cols]
        print(f"[STEP] 读取 report（{sep.strip() or ','} 分隔） ...")
        return pd.read_csv(path, sep=sep, usecols=usecols, dtype=dtype, low_memory=False)

    @staticmethod
    def normalize_label(s: pd.Series) -> pd.Series:
        s = s.astype(str)
        s = s.replace({"convert": "Convert"})
        s = s.mask(s.isin(["maintain", "early_remit", "late_remit", "relapse", "remit"]), "Non-convert")
        return s

    @staticmethod
    def filter_by_gene(df: pd.DataFrame, genes_col: str, target_gene: str) -> pd.DataFrame:
        s = df[genes_col].astype(str).fillna("")
        mask = s.str.split(r"[;,\s]+", regex=True).apply(lambda arr: target_gene in arr)
        return df.loc[mask].copy()

    @staticmethod
    def extract_all_isoforms_from_dict(data: dict, name_hint: str = "unknown") -> List[Dict[str, str]]:
        out: List[Dict[str, str]] = []
        if "sequence" in data and isinstance(data["sequence"], str) and data["sequence"].strip():
            acc = data.get("accession") or data.get("canonical") or (name_hint + "-1")
            out.append({"isoform": str(acc), "sequence": str(data["sequence"]).strip().upper()})
        iso_list = data.get("isoforms", [])
        if isinstance(iso_list, list) and iso_list:
            for i, it in enumerate(iso_list, start=1):
                seq = (it or {}).get("sequence")
                if not isinstance(seq, str) or not seq.strip():
                    continue
                acc = it.get("accession") or it.get("name") or f"isoform_{i}"
                out.append({"isoform": str(acc), "sequence": seq.strip().upper()})
        if not out:
            raise ValueError("JSON 中未找到任何可用的 isoform 序列。")
        return out

    @staticmethod
    def _extract_json_from_text(text: str) -> Optional[dict]:
        last_curly = text.rfind("}")
        if last_curly == -1:
            return None
        opens = [i for i, ch in enumerate(text[: last_curly + 1]) if ch == "{"]
        for start in reversed(opens):
            snippet = text[start : last_curly + 1]
            try:
                return json.loads(snippet)
            except Exception:
                continue
        return None

    # ---------- 单基因（与 psm.r 对应） ----------
    def _run_pepmapviz_match(
        self,
        *,
        gene: str,
        sequence: str,
        peptide_rows: pd.DataFrame,
        id_col: str,
        label_col: str,
        seq_col: str,
        genes_col: str,
        area_col: str,
        column_keep: Optional[List[str]] = None,
        normalize_IL: bool = True,
    ) -> pd.DataFrame:
        print("[STEP] _run_pepmapviz_match 调用 R ...")
        if not self.pepmapviz_script or not os.path.exists(self.pepmapviz_script):
            raise FileNotFoundError(f"找不到 PepMapViz R 脚本：{self.pepmapviz_script}")

        keep = [id_col, label_col]
        if "Global.Q.Value" in peptide_rows.columns:
            keep.append("Global.Q.Value")
        if area_col in peptide_rows.columns:
            keep.append(area_col)
        if column_keep:
            for c in column_keep:
                if c not in keep:
                    keep.append(c)

        # === 关键：发送给 R 前把两列转为字符串，保留原始位数 ===
        peptide_rows_out = _cast_preserve_cols_to_str(peptide_rows)

        payload = {
            "gene": gene,
            "sequence": sequence,
            "peptide_rows": peptide_rows_out.to_dict(orient="records"),
            "seq_col": seq_col,
            "genes_col": genes_col,
            "column_keep": keep,
            "normalize_IL": bool(normalize_IL),
        }

        try:
            proc = subprocess.run(
                [self.rscript_bin, self.pepmapviz_script],
                input=json.dumps(payload, ensure_ascii=False).encode("utf-8"),
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=self.r_timeout,
                check=False,
            )
        except FileNotFoundError as e:
            raise RuntimeError(f"无法执行 Rscript（检查 rscript_bin 参数）：{e}")
        except subprocess.TimeoutExpired:
            raise RuntimeError(f"PepMapViz 调用超时（>{self.r_timeout}s）：{self.pepmapviz_script}")

        if proc.returncode != 0:
            stderr = proc.stderr.decode(errors="ignore")
            stdout = proc.stdout.decode(errors="ignore")
            raise RuntimeError(f"PepMapViz 执行失败，返回码={proc.returncode}\nSTDERR:\n{stderr}\nSTDOUT:\n{stdout}")

        raw = proc.stdout.decode("utf-8", errors="ignore")
        try:
            out_obj = json.loads(raw)
        except Exception:
            out_obj = self._extract_json_from_text(raw)
            if out_obj is None:
                raise RuntimeError(
                    "PepMapViz 输出非 JSON，且无法从日志中提取。\n"
                    f"STDERR(前1KB):\n{proc.stderr.decode(errors='ignore')[:1024]}\n"
                    f"STDOUT(前1KB):\n{raw[:1024]}"
                )

        if out_obj.get("status") != "ok":
            raise RuntimeError(f"PepMapViz 返回错误：{out_obj}")

        df = pd.DataFrame(out_obj.get("result") or [])

        # === 关键：R 返回后附加数值副本列，原始字符串列不动 ===
        df = _attach_numeric_clones(df)

        print(f"[INFO] R 返回匹配行数：{len(df)}")
        return df

    def run_match_positions(
        self,
        *,
        report_path: str,
        gene: str,
        json_data: dict,
        isoform_select: str = "ALL",
        id_col: str = "Polypeptide.Novogene.ID",
        label_col: str = "label",
        seq_col: str = "Stripped.Sequence",
        genes_col: str = "Genes",
        area_col: str = "Precursor.Normalised",
        normalize_IL: bool = True,
    ) -> Optional[pd.DataFrame]:
        print("[STEP] run_match_positions 开始 ...")

        need = [seq_col, genes_col, id_col, label_col, "Global.Q.Value", area_col]
        rep = self.read_table_auto(report_path, usecols=None)
        missing = [c for c in [seq_col, genes_col, id_col] if c not in rep.columns]
        if missing:
            raise SystemExit(f"[ERR] report 缺少列：{', '.join(missing)}")

        rep = _apply_report_filters(rep, genes_col=genes_col, protein_col="Protein.Group")
        cols_needed = [c for c in need if c in rep.columns]
        report = rep[cols_needed].copy()

        if label_col in report.columns:
            report[label_col] = self.normalize_label(report[label_col])
        else:
            report[label_col] = "All"

        print(f"[STEP] 过滤基因 {gene} 的肽段 ...")
        df_gene = self.filter_by_gene(report, genes_col, gene)
        print(f"[INFO] 该基因肽段行数：{len(df_gene)}")
        if df_gene.empty:
            return None

        # === 关键：发送给 R 前，把需要保留位数的列改为字符串 ===
        df_gene = _cast_preserve_cols_to_str(df_gene)

        isoforms = self.extract_all_isoforms_from_dict(json_data, name_hint=gene)
        wanted = None if isoform_select.strip().upper() == "ALL" else {s.strip(): 1 for s in isoform_select.split(",")}
        if wanted is not None:
            isoforms = [it for it in isoforms if it["isoform"] in wanted]
            if not isoforms:
                raise SystemExit(f"[ERR] --isoform 没有匹配任何 JSON 内的 isoform：{isoform_select}")
        print(f"[INFO] 准备匹配的 isoform 数：{len(isoforms)}")

        all_res: List[pd.DataFrame] = []
        for it in tqdm(isoforms, desc=f"{gene} | isoform 匹配", leave=False):
            iso_name = it["isoform"]
            prot_seq = (it["sequence"] or "").strip().upper()
            if not prot_seq:
                continue

            df_res = self._run_pepmapviz_match(
                gene=gene,
                sequence=prot_seq,
                peptide_rows=df_gene,          # 已转为字符串的表
                id_col=id_col,
                label_col=label_col,
                seq_col=seq_col,
                genes_col=genes_col,
                area_col=area_col,
                column_keep=None,
                normalize_IL=normalize_IL,
            )
            if df_res is None or df_res.empty:
                continue

            # 关键列兜底
            for c in [id_col, label_col, genes_col, "start_position", "end_position"]:
                if c not in df_res.columns:
                    df_res[c] = pd.NA

            # 标注 isoform
            df_res["isoform"] = iso_name
            all_res.append(df_res)

        if not all_res:
            print("[INFO] 所有 isoform 均无匹配位置")
            return None

        res_all = pd.concat(all_res, ignore_index=True)

        base_cols = [id_col, label_col, genes_col, "start_position", "end_position", "isoform"]
        extra_cols = [c for c in res_all.columns if c not in base_cols]
        res_all = res_all[base_cols + extra_cols]

        # === 保险：如果上游某些结果没有附加数值列，这里再补一次 ===
        res_all = _attach_numeric_clones(res_all)

        print(f"[INFO] 合并后匹配行数：{len(res_all)}")
        return res_all