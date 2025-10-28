# 🧬 AutoExon

> Automated isoform/exon retrieval and unique-exon–aware PSM coverage pipeline (Python + R).

## Overview

<p align="center">
  🌏 <b>Language:</b> 
  <a href="../README.md">中文</a> | 
  <a href="./README_EN.md">English</a>
</p>

**AutoExon** integrates UniProt & EBI Proteins APIs with DIA-NN (or similar) peptide reports to:

* Retrieve **isoforms** and **exon coordinates** for a gene (UniProt/EBI).
* Run **PepMapViz** (R) to map peptides onto protein sequences and obtain peptide **start/end** positions per isoform.
* Detect **unique exons** (isoform-specific) and evaluate whether each peptide overlaps any unique-exon interval.
* Produce both **detailed** and **aggregated** CSV outputs for downstream quantification.

The pipeline is **parallel**, **resume-safe** (cache table), and **precision-preserving** for critical numeric columns.

---

## Features

| Module                | Purpose                                                                                        |
| --------------------- | ---------------------------------------------------------------------------------------------- |
| `main.py`             | Entry point; parses args and runs the pipeline.                                                |
| `pipeline.py`         | Orchestrates all stages end-to-end.                                                            |
| `report_mapping.py`   | Builds Gene→Protein mapping from the peptide report with strict filtering.                     |
| `cache_table.py`      | Maintains a resume-safe cache table (`UniProtDone`, `PepMapVizDone`).                          |
| `uniprot_prefetch.py` | Concurrent UniProt/EBI isoform & exon retrieval and JSON/CSV saving.                           |
| `exon_fetcher.py`     | Canonical-accession–based API access, exon matrix building, unique detection, CSV/JSON saving. |
| `pepmapviz_runner.py` | Parallel PepMapViz matching (R) per gene; caches per-gene results.                             |
| `psm.py`              | R invocation, precision-safe I/O, gene-level matching, and in-memory cache.                    |
| `statistics.py`       | Unique-exon hit evaluation and summary tables.                                                 |
| `psm.r`               | R side: PepMapViz glue; computes peptide positions per isoform.                                |

---

## Workflow

```
report (DIA-NN) ──► Gene→Protein mapping ──► Cache table
                         │
                         ├──► UniProt/EBI prefetch (isoform + exon JSON/CSV)
                         │
                         ├──► PepMapViz (R) parallel matching per gene
                         │
                         └──► Unique-exon hit evaluation & aggregation
                                      │
                                      └──► CSV outputs (detail & summaries)
```

---

## Output Layout

Default under `/workspace/results/`:

```
results/
├─ uniport/
│  └─ {GENE}/
│     ├─ {GENE}_exons.json
│     ├─ {GENE}_exon_matrix.csv
│     └─ pepmapviz_positions.csv
├─ match/
│  ├─ unique_psm_by_sample_all.csv
│  ├─ unique_psm_sum_by_sample.csv
│  └─ unique_psm_sum_by_exon.csv
└─ caches/
   └─ gene_protein_status.csv
```

---

## Command-line Arguments

Defined in `main.py::build_args()`:

| Arg             | Default                                | Description                                                             |
| --------------- | -------------------------------------- | ----------------------------------------------------------------------- |
| `--report`      | `/workspace/datasets/.../report-*.tsv` | DIA-NN (or similar) report file.                                        |
| `--gene`        | `ALL`                                  | Process a single gene or all genes.                                     |
| `--tax-id`      | `9606`                                 | NCBI taxonomy ID (human = 9606).                                        |
| `--isoform`     | `ALL`                                  | Target isoform(s) or ALL.                                               |
| `--id-col`      | `Polypeptide.Novogene.ID`              | Peptide ID column.                                                      |
| `--label-col`   | `label`                                | Sample/group label.                                                     |
| `--seq-col`     | `Stripped.Sequence`                    | Peptide sequence.                                                       |
| `--protein-col` | `Protein.Group`                        | Protein column.                                                         |
| `--genes-col`   | `Genes`                                | Gene column.                                                            |
| `--area-col`    | `Precursor.Normalised`                 | Quant value column.                                                     |
| `--out-uniport` | `/workspace/results/uniport`           | UniProt/EBI outputs.                                                    |
| `--out-match`   | `/workspace/results/match`             | Final result tables.                                                    |
| `--out-cache`   | `/workspace/results/caches`            | Resume cache table.                                                     |
| `--r-script`    | `src/psm.r`                            | R script (PepMapViz glue).                                              |
| `--r-timeout`   | `900`                                  | Rscript timeout (s).                                                    |
| `--workers`     | `512`                                  | Thread pool size for statistics; PepMapViz uses `workers//2` processes. |

---

## 🧩 Environment Setup

### 🐍 Python Environment

Create a virtual environment (recommended):

```bash
python3 -m venv venv
source venv/bin/activate
```

Install all dependencies:

```bash
pip install -r requirements.txt
```

**requirements.txt** (for reference):

```txt
pandas>=2.0.0
numpy>=1.24.0
tqdm>=4.65.0
requests>=2.31.0
concurrent-log-handler>=0.9.24
```

> 💡 Optional: if you are using a minimal Linux container, also install system packages:
>
> ```bash
> apt-get update && apt-get install -y python3-pip python3-venv
> ```

---

### 🧬 R Environment

Install base R and required packages:

```bash
apt-get update && apt-get install -y r-base r-base-dev
Rscript -e 'install.packages(c("ggplot2"), repos="https://cloud.r-project.org")'
Rscript -e 'install.packages(c("data.table","dplyr","jsonlite","magrittr","remotes"), repos="https://cloud.r-project.org")'
Rscript -e 'remotes::install_github("Genentech/PepMapViz")'
```

> ⚠️ **Notes**
>
> * The **PepMapViz** package from Genentech is required for peptide–sequence matching.
> * Make sure the `Rscript` command is available in your system `PATH`.
> * Test installation with:
>
>   ```bash
>   Rscript -e 'library(PepMapViz); cat("PepMapViz OK\n")'
>   ```

---

### ✅ Verify Installation

Run the following to check that both environments are ready:

```bash
python3 -c "import pandas, requests, tqdm; print('Python OK')"
Rscript -e 'library(PepMapViz); library(data.table); cat("R OK\n")'
```

Once all checks pass, you can launch the full AutoExon pipeline:

```bash
python3 main.py --report path/to/report.tsv --workers 128
```

---

## End-to-End Steps

### 1) Gene→Protein Mapping (`report_mapping.py`)

* Reads `Genes`, `Protein.Group`, `Global.Q.Value`.
* **Filters**: remove rows with `;` in `Genes` or `Protein.Group`; keep `Global.Q.Value < 0.01`.
* Builds a de-duplicated Gene→Protein set.

> Failing to find required columns or empty report will terminate with a clear error.

### 2) Resume Cache Table (`cache_table.py`)

Creates/updates `caches/gene_protein_status.csv`:

| Gene | Protein | UniProtDone | PepMapVizDone |
| ---- | ------- | ----------- | ------------- |

Used for safe resume: completed stages are not re-run.

### 3) UniProt/EBI Prefetch (`uniprot_prefetch.py`)

* For each gene, choose a **canonical** accession (strip isoform suffix).
* Fetch **isoforms** and **exon coordinates** via APIs.
* Save `{gene}_exons.json` and `{gene}_exon_matrix.csv`, then set `UniProtDone=True`.

### 4) PepMapViz Matching in Parallel (`pepmapviz_runner.py` + `psm.py` + `psm.r`)

* One process per gene (up to `workers//2`), invoking `Rscript src/psm.r`.
* Inputs are precision-preserved (critical numeric columns kept as **strings** to avoid rounding).
* Output per gene: `pepmapviz_positions.csv` (+ optional `result.json` for debugging).
* On success, `PepMapVizDone=True`.

### 5) Unique-exon Hit Evaluation (`statistics.py`)

* Load `{gene}_exons.json` and `pepmapviz_positions.csv`.
* Merge overlapping **unique-exon** intervals per isoform.
* For each peptide interval (`start_position`, `end_position`), mark:

  * `value=1` if it intersects any unique-exon interval of its isoform;
  * `value=0` otherwise (kept in detail output).
* Produce **detail** and **aggregations**.

### 6) Outputs (`pipeline.py`)

* **Detail**: `unique_psm_by_sample_all.csv`
* **By ID & interval**: `unique_psm_sum_by_sample.csv`
* **By interval only**: `unique_psm_sum_by_exon.csv`

---

## Notes & Best Practices

1. **R Environment**

   * Requires R with packages:

     ```
     PepMapViz
     data.table
     jsonlite
     stringi
     ```
   * `Rscript` must be on `PATH`.

2. **Network**

   * UniProt & EBI API calls require Internet access.
   * Partial failures are tolerated; resume via cache table.

3. **Precision Preservation**

   * Columns like `Global.Q.Value` and `Precursor.Normalised` are cast to **string** before R I/O, then numeric clones (`*_num`) are added for downstream arithmetic. This prevents unwanted rounding/formatting.

4. **Parallelism**

   * If resources are constrained, lower `--workers`.
   * PepMapViz processes additionally limit BLAS/OpenMP threads to avoid oversubscription.

5. **Filtering Rules (Report)**

   * Rows containing `;` in `Genes` or `Protein.Group` are **removed**.
   * Keep only `Global.Q.Value < 0.01`.

6. **Resume Safety**

   * `caches/gene_protein_status.csv` drives idempotent reruns.
   * Existing `{gene}/pepmapviz_positions.csv` are reused.

---

## Example Run

```bash
python3 main.py \
  --report /workspace/datasets/....tsv \
  --out-uniport /workspace/results/uniport \
  --out-match /workspace/results/match \
  --out-cache /workspace/results/caches \
  --r-script src/psm.r \
  --workers 256
```

---

## Sample Outputs

**Cache table**

```csv
Gene,Protein,UniProtDone,PepMapVizDone
ALB,P02768,True,True
TPM3,P06753,True,False
```

**Per-gene exon matrix (`{GENE}_exon_matrix.csv`)**

| chrom |  start |    end | unique | unique_isoform | overlap_conflict | overlap_reason | P02768-1 | P02768-2 |
| ----- | -----: | -----: | :----: | :------------: | :--------------: | :------------- | :------: | :------: |
| 4     | 734001 | 734078 |  TRUE  |    P02768-1    |       FALSE      | unique         |   20–30  |     -    |

**Aggregated by ID & interval (`unique_psm_sum_by_sample.csv`)**

| Polypeptide.Novogene.ID | gene | canonical | isoform  | label   | start | end | psm_value_sum | Precursor.Normalised_sum |
| ----------------------- | ---- | --------- | -------- | ------- | ----: | --: | ------------: | -----------------------: |
| ...        | ALB  | P02768    | P02768-1 | Convert |    12 |  24 |             1 |                    32.58 |

---

## Developer Notes

* Modular design; each stage can be tested independently.
* JSON/CSV are human-inspectable for debugging (per-gene subfolders).
* Unique-exon logic is derived from unified isoform–exon ownership and adjacency/overlap reasoning.

---

## Summary

AutoExon provides a robust, high-throughput pipeline for **unique-exon–aware** peptide coverage analysis:

* Automated retrieval of isoform/exon metadata,
* Precise peptide-to-isoform mapping,
* Fast parallel processing with resume support,
* Clear detailed and aggregated outputs for downstream quant/visualization.

---

