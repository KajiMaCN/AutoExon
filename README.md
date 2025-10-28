# AutoExon 使用说明文档

## 📘 项目简介

<p align="center">
  🌏 <b>语言:</b> 
  <a href="./README.md">中文</a> | 
  <a href="./documents/README_EN.md">English</a>
</p>

**AutoExon** 是一个自动化的蛋白质同种型（isoform）外显子解析与 PSM（Peptide Spectrum Match）覆盖计算系统。
该系统整合了 UniProt 与 EBI Proteins API 数据，结合 DIA-NN 或其他定量蛋白质组学软件输出的 report 文件，通过 R 脚本（PepMapViz）和 Python 多线程管线，实现：

* 自动获取基因的 **UniProt isoform / exon 坐标信息**；
* 调用 **PepMapViz** 匹配肽段与 isoform 的覆盖坐标；
* 自动识别 **unique-exon**；
* 计算肽段与 unique-exon 的交集关系；
* 输出详细与汇总的 **PSM 覆盖结果表**。

---

## ⚙️ 功能总览

| 模块                    | 功能描述                                    |
| --------------------- | --------------------------------------- |
| `main.py`             | 主入口，负责参数解析与调用完整处理流程。                    |
| `pipeline.py`         | 主管线执行逻辑，按顺序执行各阶段任务。                     |
| `report_mapping.py`   | 从 report 构建 Gene→Protein 对应表。           |
| `cache_table.py`      | 管理断点续跑状态缓存表。                            |
| `uniprot_prefetch.py` | 并发从 UniProt/EBI 获取 isoform 与 exon 坐标数据。 |
| `exon_fetcher.py`     | 提供统一的 exon 查询、矩阵构建与 unique 判定逻辑。        |
| `pepmapviz_runner.py` | 并行执行 PepMapViz 匹配（R 脚本调用）。              |
| `psm.py`              | 调用 Rscript 执行 PSM 匹配，管理缓存与数值精度。         |
| `statistics.py`       | 对已匹配结果进行 unique-exon 命中统计与汇总。           |
| `psm.r`               | 调用 R 端 PepMapViz 计算肽段在 isoform 上的坐标位置。  |

---

## 🔁 核心流程图

```
        ┌─────────────────────────────┐
        │          main.py            │
        │   (构建参数, 调用pipeline)   │
        └──────────────┬──────────────┘
                       │
        ┌──────────────▼────────────────┐
        │         pipeline.py           │
        │ 1. 构建 Gene→Protein 映射      │
        │ 2. 建立缓存状态表              │
        │ 3. UniProt 预取 isoform/exon  │
        │ 4. 调用 R 计算 PepMapViz 匹配  │
        │ 5. 统计 unique exon 命中情况   │
        └──────────────┬────────────────┘
                       │
         ┌─────────────▼────────────────┐
         │        输出结果文件           │
         │ unique_psm_by_sample_all.csv │
         │ unique_psm_sum_by_sample.csv │
         │ unique_psm_sum_by_exon.csv   │
         └──────────────────────────────┘
```

---

## 📂 输出目录结构

默认输出目录在 `/workspace/results/`：

```
results/
├── uniport/
│   ├── GENE1/
│   │   ├── GENE1_exons.json          # EBI exon坐标数据
│   │   ├── GENE1_exon_matrix.csv     # exon矩阵表
│   │   └── pepmapviz_positions.csv   # PepMapViz匹配结果
│   └── GENE2/...
├── match/
│   ├── unique_psm_by_sample_all.csv  # 明细
│   ├── unique_psm_sum_by_sample.csv  # 按样本聚合
│   └── unique_psm_sum_by_exon.csv    # 按外显子聚合
└── caches/
    └── gene_protein_status.csv       # 断点续跑缓存表
```

---

## 🧩 参数配置说明

在 `main.py` 的 `build_args()` 中定义：

| 参数名             | 默认值                                       | 说明                          |
| --------------- | ----------------------------------------- | --------------------------- |
| `--report`      | `/workspace/datasets/Match_result-...tsv` | DIA-NN 或其他蛋白质组学 report 文件路径 |
| `--gene`        | `ALL`                                     | 指定基因名称（或 `ALL` 表示全部）        |
| `--tax-id`      | `9606`                                    | 物种编号（人类为 9606）              |
| `--isoform`     | `ALL`                                     | 指定 isoform 名称（或 ALL 表示全部）   |
| `--id-col`      | `Polypeptide.Novogene.ID`                 | 肽段 ID 列名                    |
| `--label-col`   | `label`                                   | 分组标签列                       |
| `--seq-col`     | `Stripped.Sequence`                       | 肽段序列列                       |
| `--protein-col` | `Protein.Group`                           | 蛋白列                         |
| `--genes-col`   | `Genes`                                   | 基因列                         |
| `--area-col`    | `Precursor.Normalised`                    | 定量信号列                       |
| `--out-uniport` | `/workspace/results/uniport`              | UniProt/EBI 数据保存路径          |
| `--out-match`   | `/workspace/results/match`                | 结果保存路径                      |
| `--out-cache`   | `/workspace/results/caches`               | 缓存状态保存路径                    |
| `--r-script`    | `src/psm.r`                               | R 脚本路径                      |
| `--r-timeout`   | `900`                                     | Rscript 超时时间（秒）             |
| `--workers`     | `512`                                     | 并发线程数                       |

---

## 🧩 环境安装

AutoExon 运行依赖 **Python** 与 **R** 两个环境，其中 R 用于执行 PepMapViz 匹配逻辑，Python 负责主控管线与并行调度。

### 🐍 Python 环境

可使用虚拟环境进行隔离安装：

```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

**requirements.txt**

```txt
pandas>=2.0.0
numpy>=1.24.0
tqdm>=4.65.0
requests>=2.31.0
concurrent-log-handler>=0.9.24
```

> 💡 若在精简版系统中运行，请先安装基础组件：
>
> ```bash
> apt-get update && apt-get install -y python3-pip python3-venv
> ```

---

### 🧬 R 环境

R 端主要用于调用 **PepMapViz** 包执行肽段与蛋白序列的匹配。安装方法如下：

```bash
apt-get update && apt-get install -y r-base r-base-dev
Rscript -e 'install.packages(c("ggplot2"), repos="https://cloud.r-project.org")'
Rscript -e 'install.packages(c("data.table","dplyr","jsonlite","magrittr","remotes"), repos="https://cloud.r-project.org")'
Rscript -e 'remotes::install_github("Genentech/PepMapViz")'
```

> ⚠️ **说明：**
>
> * `PepMapViz` 是 Genentech 官方维护的匹配工具包。
> * 需确保命令行中 `Rscript` 可直接调用。
> * 可通过以下命令验证：
>
>   ```bash
>   Rscript -e 'library(PepMapViz); cat("PepMapViz OK\n")'
>   ```

---

### ✅ 环境检测

安装完成后可执行以下命令验证环境是否可用：

```bash
python3 -c "import pandas, requests, tqdm; print('Python OK')"
Rscript -e 'library(PepMapViz); library(data.table); cat('R OK\n')'
```

当两者均返回 OK，即可运行主程序：

```bash
python3 main.py --report path/to/report.tsv --workers 128
```

---


## 🧠 运行流程详解

### **① 构建 Gene→Protein 映射**

```python
build_gene2prot_from_report()
```

* 从 report 文件读取 `Genes`, `Protein.Group`, `Global.Q.Value`；
* 去除包含 `;` 的行；
* 过滤 `Global.Q.Value < 0.01`；
* 生成去重后的 Gene→Protein 映射字典。

> ✅ 注意：若 report 文件为空或列名错误，将直接终止程序。

---

### **② 初始化缓存状态表**

```python
ensure_cache_table()
```

生成或更新 `gene_protein_status.csv`，字段：

| Gene | Protein | UniProtDone | PepMapVizDone |
| ---- | ------- | ----------- | ------------- |
| ...  | ...     | False/True  | False/True    |

该表用于断点续跑，即程序中断后可从已完成阶段继续。

---

### **③ 并发预取 UniProt / EBI Exon 数据**

```python
prefetch_gene_json_map()
```

* 自动从 **UniProt API** 获取 canonical accession；
* 调用 **EBI Proteins API** 获取 isoform → exon 坐标；
* 生成：

  * `${gene}_exons.json`
  * `${gene}_exon_matrix.csv`

若成功则更新缓存表 `UniProtDone=True`。

> ⚠️ 若 EBI 或 UniProt 网络波动，会自动跳过并继续其他基因。

---

### **④ 并行调用 R 进行 PepMapViz 匹配**

```python
pepmapviz_runner.pepmapviz_thread_worker()
```

* 逐基因调用 R 脚本 `psm.r`；
* 匹配每个 isoform 上的肽段坐标；
* 保存 `pepmapviz_positions.csv`；
* 成功后更新 `PepMapVizDone=True`。

**并发策略**

* 每次最多运行 `workers // 2` 个子进程；
* 内部设置环境变量防止 OpenBLAS/MKL 并行爆核。

---

### **⑤ 并发统计 unique-exon 覆盖情况**

```python
statistics.process_one_gene()
```

* 读取：

  * `pepmapviz_positions.csv`
  * `${gene}_exons.json`
* 从 JSON 中提取 unique 外显子；
* 判断每条肽段坐标是否与任意 unique-exon 区间相交；
* 输出：

  * 命中 value=1
  * 未命中 value=0（仍保留）

---

### **⑥ 汇总与输出结果**

共生成三类输出：

| 文件名                            | 内容说明              |
| ------------------------------ | ----------------- |
| `unique_psm_by_sample_all.csv` | 每条肽段的命中明细（带样本与标签） |
| `unique_psm_sum_by_sample.csv` | 按 ID 与外显子区间聚合求和   |
| `unique_psm_sum_by_exon.csv`   | 按外显子区间聚合求和（不分样本）  |

---

## 🧾 注意事项

1. **R 依赖环境**

   * 需安装 R 以及以下包：

     ```
     PepMapViz
     data.table
     jsonlite
     stringi
     ```
   * Rscript 命令需在 PATH 中可执行。

2. **网络依赖**

   * UniProt 与 EBI 均需外网访问；
   * 若遇断连，可在缓存表中查看 `UniProtDone` 状态后重跑。

3. **断点续跑**

   * 每次运行前自动检测 `gene_protein_status.csv`；
   * 已完成的部分不会重复下载或计算。

4. **文件命名**

   * 目录名统一按基因名；
   * 所有临时或中间文件均会保存，以便人工检查。

5. **性能建议**

   * 推荐使用多核服务器；
   * 若内存或 CPU 不足，可降低 `--workers`。

6. **报错处理**

   * 日志中 `[WARN]` 表示非致命错误（如某基因缺少 JSON）；
   * `[ERR]` 表示严重错误，程序会中止。

---

## 🧪 示例运行

```bash
python3 main.py \
  --report /workspace/datasets/....tsv \
  --out-uniport /workspace/results/uniport \
  --out-match /workspace/results/match \
  --out-cache /workspace/results/caches \
  --r-script src/psm.r \
  --workers 256
```

程序会自动依次执行所有阶段，并在终端打印进度条与日志信息。

---

## 📜 典型输出示例

### 1️⃣ `gene_protein_status.csv`

```csv
Gene,Protein,UniProtDone,PepMapVizDone
ALB,P02768,True,True
TPM3,P06753,True,False
```

### 2️⃣ `ALB_exon_matrix.csv`

| chrom | start  | end    | unique | unique_isoform | overlap_reason | P02768-1 | P02768-2 |
| ----- | ------ | ------ | ------ | -------------- | -------------- | -------- | -------- |
| 4     | 734001 | 734078 | TRUE   | P02768-1       | unique         | 20-30    | -        |

### 3️⃣ `unique_psm_sum_by_sample.csv`

| Polypeptide.Novogene.ID | gene | canonical | isoform  | label   | start | end | psm_value_sum | Precursor.Normalised_sum |
| ----------------------- | ---- | --------- | -------- | ------- | ----- | --- | ------------- | ------------------------ |
| ...       | ALB  | P02768    | P02768-1 | Convert | 12    | 24  | 1             | 32.58                    |

---

## 🧩 开发者备注

* 源码结构遵循模块化设计，可独立运行或组合调用；
* 所有数据流均为显式输入输出，便于调试；
* JSON 与 CSV 输出格式与 PepMapViz 完全兼容。

---

## 🏁 总结

AutoExon 是一个高并发、断点可续跑、跨语言（Python+R）协同的自动化外显子识别与肽段匹配框架。
其主要优势：

* 自动拉取 UniProt/EBI 外显子数据；
* 精确计算肽段覆盖关系；
* 支持大规模并行；
* 提供详细与聚合层次的输出；
* 可灵活扩展到其他物种或分析流程。

