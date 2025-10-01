import os
import json
import uuid
import subprocess
import pandas as pd

from typing import Dict, List, Optional, Tuple
from pathlib import Path


class IsoformCoverageCalculator:
    """
    仅保存版：多 isoform 覆盖计算（调用外部 PepMapViz R 脚本）并必定保存 CSV。
    - 只接受已解析 JSON(dict) 的 isoform 数据
    - 无 run_once / 无 dry-run
    - 解析 R 输出具备容错（支持 out_file 或从 stdout 混杂文本中提取 JSON）
    """

    def __init__(self, pepmapviz_script: str, rscript_bin: str = "Rscript", r_timeout: int = 900):
        self.pepmapviz_script = pepmapviz_script
        self.rscript_bin = rscript_bin
        self.r_timeout = r_timeout

    
    # ------- IO 工具 -------
    @staticmethod
    def read_gene_protein_map(path: str,
                          gene_col: str,
                          protein_col: str) -> dict:
        """
        从报告文件提取基因-蛋白映射关系。
        规则：
        - gene_col 与 protein_col 可能是用 ';' 分隔的多值
        - 按分号顺序一一对应展开
        - 去重
        返回:
        { gene: set([protein1, protein2, ...]) }
        """
        sep = "\t" if path.endswith((".tsv", ".tsv.gz")) else ","
        df = pd.read_csv(path, sep=sep, usecols=[gene_col, protein_col], dtype=str, low_memory=False)

        df = df.dropna(subset=[gene_col, protein_col])

        mapping = {}

        for g_str, p_str in zip(df[gene_col], df[protein_col]):
            genes = [g.strip() for g in str(g_str).split(";") if g.strip()]
            prots = [p.strip() for p in str(p_str).split(";") if p.strip()]

            # 如果数量相等，逐一对应
            if len(genes) == len(prots):
                pairs = zip(genes, prots)
            else:
                # 数量不等时，采取保守做法：所有基因都映射到所有蛋白（笛卡尔积）
                pairs = [(g, p) for g in genes for p in prots]

            for g, p in pairs:
                mapping.setdefault(g, set()).add(p)

        return mapping

    @staticmethod
    def read_table_auto(path: str, usecols: Optional[List[str]] = None, dtype=None) -> pd.DataFrame: 
        sep = "\t" if path.endswith((".tsv", ".tsv.gz")) else "," 
        header = pd.read_csv(path, sep=sep, nrows=0, low_memory=False) 
        cols = header.columns.tolist() 
        if usecols is not None: 
            usecols = [c for c in usecols if c in cols] 
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

    # ------- JSON -> isoforms -------
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
            raise ValueError("JSON 数据中未找到任何可用的 isoform 序列。")
        return out

    # ------- 从混杂文本提取最后一个 JSON 片段 -------
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

    # ------- 调用 R 计算 -------
    def _run_pepmapviz(
        self,
        *,
        gene: str,
        sequence: str,
        peptide_rows: pd.DataFrame,
        id_col: str,
        label_col: str,
        seq_col: str,
        genes_col: str,
        quantify: str,
        area_col: str,
        agg: str,
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        if not self.pepmapviz_script or not os.path.exists(self.pepmapviz_script):
            raise FileNotFoundError(f"找不到 PepMapViz R 脚本：{self.pepmapviz_script}")

        # 建议让 R 写文件（若 R 实现了 out_file 就会生效）
        out_dir = Path(f"/workspace/results/{gene}")
        out_dir.mkdir(parents=True, exist_ok=True)
        out_file = out_dir / f"pepmapviz_{uuid.uuid4().hex}.json"

        payload = {
            "gene": gene,
            "sequence": sequence,
            "id_col": id_col,
            "label_col": label_col,
            "seq_col": seq_col,
            "genes_col": genes_col,
            "quantify": quantify,
            "area_col": area_col,
            "agg": agg,
            "peptide_rows": peptide_rows.to_dict(orient="records"),
            # "out_file": str(out_file),  # R 支持则写文件；不支持则忽略
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

        # 1) 优先从 out_file 读（如果 R 真写了文件）
        out_obj = None
        if out_file.exists():
            try:
                out_txt = out_file.read_text(encoding="utf-8")
                out_obj = json.loads(out_txt)
            except Exception:
                out_obj = None

        # 2) 否则从 stdout 读取，并兼容日志噪声
        if out_obj is None:
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

        by_sample = pd.DataFrame(out_obj.get("by_sample") or [], columns=["Character", "Position", "id", "label", "value"])
        by_label = pd.DataFrame(out_obj.get("by_label") or [], columns=["Character", "Position", "label", "value"])
        return by_sample, by_label

    # ------- 仅保存版入口 -------
    def run_and_save(
        self,
        *,
        report_path: str,
        gene: str,
        json_data: dict,
        out_sample_csv: str,
        out_label_csv: str,
        isoform_select: str = "ALL",
        id_col: str = "Polypeptide.Novogene.ID",
        label_col: str = "label",
        seq_col: str = "Stripped.Sequence",
        genes_col: str = "Genes",
        area_col: str = "Precursor.Normalised",
        quantify: str = "PSM",
        agg: str = "mean",
    ) -> None:
        # 1) 读报告
        need = [seq_col, genes_col, id_col, label_col]
        if quantify.upper() == "AREA":
            need.append(area_col)
        report = self.read_table_auto(report_path, usecols=need)

        if label_col in report.columns:
            report[label_col] = self.normalize_label(report[label_col])
        else:
            report[label_col] = "All"

        df_gene = self.filter_by_gene(report, genes_col, gene)
        if df_gene.empty:
            raise SystemExit(f"[ERR] 基因 {gene} 未在 {genes_col} 中找到。")

        # 2) isoforms
        isoforms = self.extract_all_isoforms_from_dict(json_data, name_hint=gene)
        wanted = None if isoform_select.strip().upper() == "ALL" else {s.strip(): 1 for s in isoform_select.split(",")}
        if wanted is not None:
            isoforms = [it for it in isoforms if it["isoform"] in wanted]
            if not isoforms:
                raise SystemExit(f"[ERR] --isoform 没有匹配任何 JSON 内的 isoform：{isoform_select}")

        pep_cols = [seq_col, genes_col, id_col, label_col]
        if quantify.upper() == "AREA" and area_col in df_gene.columns:
            pep_cols.append(area_col)
        pep_rows = df_gene[pep_cols].copy()

        all_by_sample: List[pd.DataFrame] = []
        all_by_label: List[pd.DataFrame] = []

        for it in isoforms:
            iso_name = it["isoform"]
            prot_seq = (it["sequence"] or "").strip().upper()
            if not prot_seq:
                continue

            by_sample_iso, by_label_iso = self._run_pepmapviz(
                gene=gene,
                sequence=prot_seq,
                peptide_rows=pep_rows,
                id_col=id_col,
                label_col=label_col,
                seq_col=seq_col,
                genes_col=genes_col,
                quantify=quantify.upper(),
                area_col=area_col,
                agg=agg.lower(),
            )

            if not by_sample_iso.empty:
                by_sample_iso["isoform"] = iso_name
                by_sample_iso = by_sample_iso[["Character", "Position", "id", "label", "value", "isoform"]]
                all_by_sample.append(by_sample_iso)

            if not by_label_iso.empty:
                by_label_iso["isoform"] = iso_name
                by_label_iso = by_label_iso[["Character", "Position", "label", "value", "isoform"]]
                all_by_label.append(by_label_iso)

        if not all_by_sample:
            raise SystemExit("[ERR] 所有 isoform 上均未得到覆盖。请检查序列/肽段是否匹配。")

        by_sample = pd.concat(all_by_sample, ignore_index=True)
        by_label = (
            pd.concat(all_by_label, ignore_index=True)
            if all_by_label
            else by_sample.groupby(["Character", "Position", "label", "isoform"], as_index=False)["value"].mean()
        )

        # 3) 必定保存（打印绝对路径）
        out_sample_csv = str(Path(out_sample_csv).resolve())
        out_label_csv  = str(Path(out_label_csv).resolve())
        Path(out_sample_csv).parent.mkdir(parents=True, exist_ok=True)
        Path(out_label_csv).parent.mkdir(parents=True, exist_ok=True)

        by_sample.to_csv(out_sample_csv, index=False)
        by_label.to_csv(out_label_csv, index=False)

        print(f"[OK] 已写出：{out_sample_csv}")
        print(f"[OK] 已写出：{out_label_csv}")
        return by_label, by_sample
