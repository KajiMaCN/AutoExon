import argparse, json, subprocess, sys, os, re, shutil
import pandas as pd
from pathlib import Path

from dataloader import load_config, load_dataset

class Exons:
    def __init__(self,accession, is_unique, sequence, exon_list):
        self.accession = accession
        self.is_unique = is_unique
        self.sequence=sequence
        self.exon_list = exon_list
    
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

class ExonPSM:
    def __init__(self,):
        self.exon=None

    def length(self):
        return self.end - self.start + 1

    def load_data(self, filepath):
        load_dataset(filepath)
    
    @staticmethod
    def load_json(file_path):
        with open(file_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
        return data
    
    def create_exons(self, isoforms):
        exon_detial={}
        gene=isoforms.get("gene","unknown")
        tax_id=isoforms.get("tax_id","unknown")
        canonical=isoforms.get("canonical","unknown")
        isoform_list=isoforms.get("isoforms",[])


        if isoform_list is None or len(isoform_list)==0:
            raise ValueError("No isoforms found in the provided JSON data.")
        else:
            for item in isoform_list:
                accession=item.get("accession","unknown")
                sequence=item.get("sequence","")
                exons_data=item.get("exons",[])
                for exon in exons_data:
                    exon_detial['start']=exon.get("start",0)
                    exon_detial['end']=exon.get("end",0)
                    exon_detial['exon']=exon.get("exon","unknown")

                if exons_data is None or len(exons_data)==0:
                    print(f"No exon data found for isoform {accession}. Skipping.")
                    continue
                exons_df=pd.DataFrame(exons_data)
                standardized_exons=Exons.standardize_exons(exons_df)
                print(f"Gene: {gene}, Tax ID: {tax_id}, Isoform: {accession}, Exons:\n{standardized_exons}")



    def normalize_label(s: pd.Series) -> pd.Series:
        s = s.astype(str)
        s = s.replace({"convert": "Convert"})
        s = s.mask(s.isin(["maintain","early_remit","late_remit","relapse","remit"]), "Non-convert")
        return s
    
    def run(self, isoforms):
        self.create_exons(isoforms)
        

if __name__ == "__main__":
    exon_psm=ExonPSM()
    isoforms=exon_psm.load_json("/workspace/results/exon/ITIH1/ITIH1_exons.json")
    exon_psm.run(isoforms)