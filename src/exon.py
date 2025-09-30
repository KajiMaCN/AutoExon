#!/usr/bin/env python3
import argparse, json, subprocess, sys, os, re, shutil
import pandas as pd
from pathlib import Path

# ---------- 读表（自动识别 .csv/.tsv；usecols 不存在的列会被忽略） ----------
def read_table_auto(path: str, usecols=None, dtype=None) -> pd.DataFrame:
    sep = "\t" if path.endswith((".tsv",".tsv.gz")) else ","
    header = pd.read_csv(path, sep=sep, nrows=0, low_memory=False)
    cols = header.columns.tolist()
    if usecols is not None:
        usecols = [c for c in usecols if c in cols]
    return pd.read_csv(path, sep=sep, usecols=usecols, dtype=dtype, low_memory=False)

# ---------- label 规范（沿用你最初规则） ----------
def normalize_label(s: pd.Series) -> pd.Series:
    s = s.astype(str)
    s = s.replace({"convert": "Convert"})
    s = s.mask(s.isin(["maintain","early_remit","late_remit","relapse","remit"]), "Non-convert")
    return s

# ---------- 基因过滤：支持多值（; , 或空格分隔） ----------
def filter_by_gene(df: pd.DataFrame, genes_col: str, target_gene: str) -> pd.DataFrame:
    s = df[genes_col].astype(str).fillna("")
    mask = s.str.split(r"[;,\s]+", regex=True).apply(lambda arr: target_gene in arr)
    return df.loc[mask].copy()

# ---------- 在蛋白序列中找到肽段所有匹配（1-based，允许重叠） ----------
def find_all(st: str, sub: str):
    res, i, stl, subl = [], 0, st.lower(), sub.lower()
    while True:
        j = stl.find(subl, i)
        if j == -1: break
        res.append(j+1)
        i = j+1
    return res

# ---------- 从多 isoform JSON 选定 isoform，并标准化 exons ----------
def extract_isoform_from_json(json_path: str, accession: str|None=None,
                              index: int|None=None, prefer_canonical: bool=True):
    J = json.loads(Path(json_path).read_text(encoding="utf-8"))
    if "sequence" in J and "exons" in J:
        seq = J["sequence"]; exons = pd.DataFrame(J["exons"])
        return standardize_exons(exons), seq
    isoforms = J.get("isoforms", [])
    if not isoforms:
        raise ValueError("JSON 未包含 sequence/exons，也没有 isoforms 列表。")
    sel = None
    if accession:
        for it in isoforms:
            if it.get("accession") == accession: sel = it; break
    if sel is None and isinstance(index, int) and 1 <= index <= len(isoforms):
        sel = isoforms[index-1]
    if sel is None and prefer_canonical:
        can = J.get("canonical")
        if can:
            cand1 = f"{can}-1"
            for it in isoforms:
                if it.get("accession") == cand1: sel = it; break
            if sel is None:
                for it in isoforms:
                    if it.get("accession") == can: sel = it; break
    if sel is None:
        for it in isoforms:
            if it.get("sequence"): sel = it; break
    if sel is None or not sel.get("sequence"):
        raise ValueError("未找到带 sequence 的 isoform。")
    seq = sel["sequence"]; exons = pd.DataFrame(sel.get("exons", []))
    return standardize_exons(exons), seq

def standardize_exons(exons: pd.DataFrame) -> pd.DataFrame:
    cols = {c.lower(): c for c in exons.columns}
    def pick(cands):
        for c in cands:
            if c in cols: return cols[c]
        return None
    c_start = pick(["start","begin","from","start_pos","pos_start"])
    c_end   = pick(["end","to","stop","end_pos","pos_end"])
    c_name  = pick(["exon","name","id","label","type","exon_id"])
    if not c_start or not c_end:
        raise ValueError("exons 里找不到 start/end 列")
    if not c_name:
        exons["exon"] = [f"Exon_{i}" for i in range(1, len(exons)+1)]
    else:
        exons = exons.rename(columns={c_name: "exon"})
    exons = exons.rename(columns={c_start: "start", c_end: "end"})
    if "domain_color" not in exons.columns:
        if "unique" in exons.columns:
            exons["domain_color"] = exons["unique"].map(lambda v: "#F8766D" if bool(v) else "#B79F00")
        else:
            pal = ["#F8766D","#B79F00"]
            exons["domain_color"] = [pal[i % len(pal)] for i in range(len(exons))]
    exons["start"] = exons["start"].astype(int)
    exons["end"]   = exons["end"].astype(int)
    exons["exon"]  = exons["exon"].astype(str)
    return exons[["start","end","exon","domain_color"]]

# ---------- 从 R 输出中宽松提取 JSON ----------
def extract_json_from_stream(s: str) -> dict:
    start = s.find("{"); end = s.rfind("}")
    if start != -1 and end != -1 and end > start:
        frag = s[start:end+1]
        return json.loads(frag)
    raise ValueError("未在 R 输出中找到 JSON 对象")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", default="/workspace/datasets/Match_result-X401SC24071912_Z01_F001_B1_43/report-matched-M0.csv")
    ap.add_argument("--gene", default="ITIH1")
    ap.add_argument("--json", default="/workspace/results/exon/ITIH1/ITIH1_exons.json")
    ap.add_argument("--isoform", default="P19827-1")
    ap.add_argument("--outfile", default="/workspace/results/plot/peptides-ITIH1-PSM-bySample.pdf")
    ap.add_argument("--title",   default=None)
    ap.add_argument("--quantify", default="PSM", choices=["PSM","Area"])
    ap.add_argument("--area-col", default="Precursor.Normalised")
    ap.add_argument("--id-col",    default="Polypeptide.Novogene.ID")
    ap.add_argument("--label-col", default="label")
    ap.add_argument("--seq-col",   default="Stripped.Sequence")
    ap.add_argument("--genes-col", default="Genes")
    ap.add_argument("--width", type=float, default=12.0)
    ap.add_argument("--height", type=float, default=100)
    ap.add_argument("--agg", default="mean", choices=["mean","sum","median"])
    ap.add_argument("--plot", type=int, choices=[0,1], default=1)
    ap.add_argument("--rscript", default=None)  # 自动探测
    ap.add_argument("--rfile",   default="/workspace/src/pms.r")
    ap.add_argument("--out-sample-csv", default="/workspace/results/psm_by_sample.csv")
    ap.add_argument("--out-label-csv",  default="/workspace/results/psm_by_label.csv")
    ap.add_argument("--plot-level", dest="plot_level",
                default="sample", choices=["label","sample"],
                help="绘图层级：label=按组别聚合绘图；sample=按样本绘图")
    args = ap.parse_args()

    # 读取 CSV（只取需要的列）
    need = [args.seq_col, args.genes_col, args.id_col, args.label_col, args.area_col]
    report = read_table_auto(args.report, usecols=need)
    # 规范 label
    if args.label_col in report.columns:
        report[args.label_col] = normalize_label(report[args.label_col])
    else:
        report[args.label_col] = "All"

    # 按基因过滤（支持多基因字符串）
    df_gene = filter_by_gene(report[[c for c in [args.seq_col,args.genes_col,args.id_col,args.label_col,args.area_col] if c in report.columns]],
                             args.genes_col, args.gene)
    if df_gene.empty:
        sys.exit(f"[ERR] 基因 {args.gene} 在 {args.genes_col} 中未找到。")

    # 从 JSON 选择 isoform & exons
    exons, protein_seq = extract_isoform_from_json(args.json, accession=args.isoform)

    # 构建 residue-level matching_result（如 quantify=Area，会把面积列展开到每个残基）
    cols_for_pep = [args.seq_col, args.genes_col, args.id_col, args.label_col]
    if args.quantify == "Area" and args.area_col in df_gene.columns:
        cols_for_pep.append(args.area_col)
    pep_rows = df_gene[cols_for_pep].copy()
    if pep_rows.empty:
        sys.exit("[ERR] 传给 PepMapViz 的 peptide_rows 为空，请检查基因过滤与列名设置。")
    # 组装 payload 给 R

    payload = {
        "gene": args.gene,
        "sequence": protein_seq,                      # 供 PepMapViz 做匹配
        "peptide_rows": pep_rows.to_dict(orient="records"),  # ← 肽段级表，一条记录=一个样本的一条肽
        "id_col": args.id_col,
        "label_col": args.label_col,
        "seq_col": args.seq_col,
        "genes_col": args.genes_col,
        "quantify": args.quantify,                    # "PSM" / "Area"
        "area_col": args.area_col,                    # 仅 Area 模式会用到
        "agg": args.agg                               # by_label 的聚合方式
        # 不需要任何绘图相关字段
    }

    # Rscript 路径 & 环境（抑制 locale 噪音）
    rscript = args.rscript or shutil.which("Rscript") or "/usr/bin/Rscript"
    if not os.path.exists(rscript):
        sys.exit(f"[ERR] 找不到 Rscript：{rscript}")
    if not os.path.exists(args.rfile):
        sys.exit(f"[ERR] 找不到 R 脚本：{args.rfile}")
    env = os.environ.copy()
    env.update({"LC_ALL":"C.UTF-8","LANG":"C.UTF-8"})

    # 调 R
    proc = subprocess.run(
        [rscript, args.rfile],
        input=json.dumps(payload, ensure_ascii=False).encode("utf-8"),
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env
    )
    raw_stdout = proc.stdout.decode("utf-8", errors="ignore")
    raw_stderr = proc.stderr.decode("utf-8", errors="ignore")
    if not raw_stdout.strip():
        sys.stderr.write(raw_stderr)
        sys.exit(proc.returncode or 1)

    # 宽松解析 R 的 JSON 输出
    try:
        res = extract_json_from_stream(raw_stdout)
    except Exception as e:
        print("----- R STDOUT BEGIN -----")
        print(raw_stdout)
        print("----- R STDOUT END -----")
        print("----- R STDERR BEGIN -----", file=sys.stderr)
        print(raw_stderr, file=sys.stderr)
        print("----- R STDERR END -----", file=sys.stderr)
        sys.stderr.write(f"[ERR] 解析 R 输出失败：{e}\n")
        sys.exit(1)

    if res.get("status") != "ok":
        sys.stderr.write(f"[ERR] R 端失败：{res.get('message')}\n")
        sys.exit(1)

    # 拿到 by_sample / by_label
    df_by_sample = pd.DataFrame(res.get("by_sample") or {})
    df_by_label  = pd.DataFrame(res.get("by_label") or {})
    print(f"by_sample={df_by_sample.shape}  by_label={df_by_label.shape}  metric={res.get('metric')}  agg={res.get('agg')}")

    if args.out_sample_csv:
        df_by_sample.to_csv(args.out_sample_csv, index=False)
        print(f"[OK] 已写出：{args.out_sample_csv}")
    if args.out_label_csv:
        df_by_label.to_csv(args.out_label_csv, index=False)
        print(f"[OK] 已写出：{args.out_label_csv}")

if __name__ == "__main__":
    main()
