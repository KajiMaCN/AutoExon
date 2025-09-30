import json
import requests
import time
import pprint
import pandas as pd

def get_isoforms(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    r = requests.get(url)
    if r.status_code != 200:
        print(f"[{uniprot_id}] Failed with status {r.status_code}")
        return None
    data = r.json()
    result = {}
    # Canonical sequence
    result["canonical"] = {
        "sequence": data["sequence"]["value"],
        "length": data["sequence"]["length"]
    }
    # Isoforms
    isoforms = []
    for comment in data.get("comments", []):
        if comment.get("commentType") == "ALTERNATIVE PRODUCTS":
            for iso in comment.get("isoforms", []):
                entry = {
                    "isoform_id": iso.get("isoformIds", [""])[0],
                    "name": iso.get("name", ""),
                    "seq_ids": iso.get("sequenceIds", []),
                }
                isoforms.append(entry)
    result["isoforms"] = isoforms
    return result


def get_exon_coordinates(uniprot_id, tax_id=None):
    """
    获取给定 UniProt 蛋白/异构体的外显子在蛋白序列上的坐标（begin/end 为蛋白位点）。
    - uniprot_id: 可为主号如 'Q8TES7'，也可为异构体如 'Q8TES7-6'
    - tax_id: 可选，如 '9606'
    """
    base_url = "https://www.ebi.ac.uk/proteins/api/coordinates"

    # 拆分主号与异构体号
    if "-" in uniprot_id:
        acc, iso = uniprot_id.split("-", 1)
        iso = f"{acc}-{iso}"
    else:
        acc, iso = uniprot_id, None

    params = {"accession": acc}
    if iso:
        # 明确指定异构体
        params["isoformAccession"] = iso
    if tax_id:
        params["taxid"] = str(tax_id)

    r = requests.get(base_url, params=params, headers={"Accept": "application/json"}, timeout=30)
    r.raise_for_status()
    data = r.json()

    if not isinstance(data, list) or len(data) == 0:
        raise RuntimeError(f"[{uniprot_id}] coordinates API 返回为空：{data!r}")

    # 可能返回多个条目，优先挑有 gnCoordinate/exon/proteinLocation 的
    def has_protein_location(entry):
        gcs = entry.get("gnCoordinate") or []
        for gc in gcs:
            exons = ((gc.get("genomicLocation") or {}).get("exon")) or []
            for ex in exons:
                pl = ex.get("proteinLocation")
                if isinstance(pl, dict):
                    return True
        return False

    # 先按是否有 proteinLocation 过滤，否则退而求其次用第一个
    entry = next((e for e in data if has_protein_location(e)), data[0])

    result = {
        "accession": entry.get("accession"),
        "sequence": entry.get("sequence"),
        "exons": []
    }

    gene_coords = entry.get("gnCoordinate") or []
    if len(gene_coords) == 0:
        # 常见于某些物种或未注释条目
        print(f"[{uniprot_id}] 未找到 gnCoordinate。原始 keys: {list(entry.keys())}")
        return result

    # 也可能不止一个转录本；此处简单选第一个包含 exon 的
    chosen_gc = None
    for gc in gene_coords:
        exons = ((gc.get("genomicLocation") or {}).get("exon")) or []
        if exons:
            chosen_gc = gc
            break
    if chosen_gc is None:
        print(f"[{uniprot_id}] gnCoordinate 存在但无 exon 列表。")
        return result

    exons = ((chosen_gc.get("genomicLocation") or {}).get("exon")) or []

    def _pos(x):
        """安全取得 position 值：x 可能为 None / dict / 直接数值。"""
        if x is None:
            return None
        if isinstance(x, dict):
            return x.get("position")
        # 偶尔服务也可能给数值/字符串
        return x

    exon_info = []
    for i, exon in enumerate(exons, start=1):
        pl = exon.get("proteinLocation")
        if not isinstance(pl, dict):
            # 没有蛋白定位信息，跳过或记录为 None
            # 你也可以选择保留记录但给 None
            # exon_info.append({"exon": i, "exon_id": exon.get("id"), "start": None, "end": None})
            continue

        begin = _pos(pl.get("begin"))
        end = _pos(pl.get("end"))

        # 有些条目只有 begin 或只有 end，或标注未知；这里仅保留两端都有的
        if begin is None or end is None:
            # 如果希望保留也可以改为记录 None
            # exon_info.append({"exon": i, "exon_id": exon.get("id"), "start": begin, "end": end})
            continue

        exon_info.append({
            "exon": i,
            "exon_id": exon.get("id"),
            "start": int(begin),
            "end": int(end),
        })

    result["exons"] = exon_info
    return result

if __name__ == "__main__":

    kng1_1_uid = "Q8TES7-1"
    kng1_2_uid = "Q8TES7-6"

    result = get_isoforms(kng1_1_uid)
    result = get_isoforms("Q8TES7")
    pprint.pprint(result)

    # Retrieve sequence and exon coordinates for ITIH1-1
    # result_itih1_1 = get_exon_coordinates(itih1_1_uid)
    # result_itih1_3 = get_exon_coordinates(itih1_3_uid)

    result = get_exon_coordinates(kng1_2_uid)
    pprint.pprint(result)

    # Save to JSON
    filepath = "/workspace/results/alb-1.json"
    with open(filepath, 'w') as f:
        json.dump(result, f, indent=2)


    filepath = 'tmp/astral/lyriks402/biomarkers/biomarkers-ancova.csv'
    bm_ancova = pd.read_csv(filepath, index_col=0)

    filepath = 'tmp/astral/lyriks402/biomarkers/biomarkers-elasticnet-nestedkfold.csv'
    bm_enet = pd.read_csv(filepath, index_col=0)
    bm_enet.head()

    filepath = 'tmp/astral/lyriks402/biomarkers/mongan-etable5.csv'
    mongan = pd.read_csv(filepath, index_col=0)
    monganq = mongan[mongan.q < 0.05]

    for uid, gene in zip(bm_ancova.index, bm_ancova.Gene):
        result = get_isoforms(uid)
        isoforms = result.get('isoforms', [])
        if isoforms:
            print(uid, gene)
            # pprint.pprint(isoforms)
            print()

    for uid, gene in zip(bm_enet.index, bm_enet.Gene):
        result = get_isoforms(uid)
        isoforms = result.get('isoforms', [])
        if isoforms:
            print(uid, gene)
            # pprint.pprint(isoforms)
            print()

    for uid, gene in zip(monganq.index, monganq['Protein name']):
        result = get_isoforms(uid)
        isoforms = result.get('isoforms', [])
        if isoforms:
            print(uid, gene)
            pprint.pprint(isoforms)
            print()
