#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  suppressWarnings(library(PepMapViz))
  suppressWarnings(library(data.table))
  suppressWarnings(library(jsonlite))
  suppressWarnings(library(stringi))
})

# 本地化无害，尽量干净
try({Sys.setenv(LANG="C.UTF-8"); suppressWarnings(Sys.setlocale("LC_ALL","C.UTF-8"))}, silent=TRUE)

J    <- function(x) toJSON(x, auto_unbox=TRUE, null="null", dataframe="columns")
fail <- function(msg){ cat(J(list(status="error", message=as.character(msg)))); quit(status=1) }
ok_empty <- function(gene, metric="PSM", agg="mean"){
  out <- list(
    status="ok", gene=gene, metric=metric, agg=agg,
    by_sample=data.frame(Character=character(), Position=integer(), id=character(), label=character(), value=double()),
    by_label =data.frame(Character=character(), Position=integer(), label=character(), value=double())
  )
  return(out)
}
IL <- function(x) chartr("IJ","LL", toupper(x %||% ""))  # I/J -> L
`%||%` <- function(a,b) if (is.null(a)) b else a

main <- function() {
  X <- fromJSON(paste(readLines(file("stdin"), warn=FALSE), collapse="\n"), simplifyVector=TRUE)
  req <- function(nm){ if (is.null(X[[nm]])) fail(paste("缺少字段:", nm)); X[[nm]] }
  opt <- function(nm, dflt){ if (!is.null(X[[nm]])) X[[nm]] else dflt }

  gene       <- req("gene")
  seq_string <- req("sequence")
  id_col     <- opt("id_col", "Polypeptide.Novogene.ID")
  label_col  <- opt("label_col", "label")
  seq_col    <- opt("seq_col", "Stripped.Sequence")
  genes_col  <- opt("genes_col", "Genes")
  quantify   <- toupper(opt("quantify","PSM"))          # "PSM"/"AREA"
  area_col   <- opt("area_col","Precursor.Normalised")
  agg        <- tolower(opt("agg","mean"))              # "mean"/"sum"/"median"
  out_file   <- opt("out_file","")

  # ---- 输入表 ----
  if (is.null(X[["peptide_rows"]])) fail("需要 peptide_rows（肽段级表）。")
  pep_tbl <- as.data.frame(X[["peptide_rows"]], check.names=FALSE, stringsAsFactors=FALSE)

  need <- c(seq_col, genes_col, id_col, label_col)
  miss <- setdiff(need, names(pep_tbl))
  if (length(miss)) fail(paste("peptide_rows 缺列:", paste(miss, collapse=", ")))

  # 轻度清洗：去空串
  pep_tbl[[seq_col]] <- trimws(as.character(pep_tbl[[seq_col]]))
  pep_tbl <- pep_tbl[nzchar(pep_tbl[[seq_col]]), , drop=FALSE]

  # 序列为空直接返回空结果（防崩）
  if (is.null(seq_string) || !nzchar(seq_string) || nrow(pep_tbl) == 0L) {
    out <- ok_empty(gene, metric=if (quantify=="AREA") "Area" else "PSM", agg=agg)
    if (nzchar(out_file)) jsonlite::write_json(out, out_file, auto_unbox=TRUE)
    cat(J(out)); quit(save="no")
  }

  # ---- I/L 等价预处理（提升命中率，仍是严格匹配逻辑）----
  seq_il <- IL(seq_string)
  pep_tbl[[seq_col]] <- IL(pep_tbl[[seq_col]])

  # PepMapViz 需要 whole_seq$Region_Sequence
  whole_seq <- data.frame(Region_Sequence = seq_il, Genes = gene, stringsAsFactors = FALSE)

  # 仅保留需要透传到后续聚合的列
  keep_cols <- c(id_col, label_col)
  if (quantify=="AREA" && area_col %in% names(pep_tbl)) keep_cols <- c(keep_cols, area_col)

  # ---- 匹配定位（使用传入列名；不要硬编码 "Genes"）----
  # 注意：采用“位置参数”的调用方式与 PepMapViz 的函数签名一致
  mr <- tryCatch({
    match_and_calculate_positions(
      pep_tbl,
      column = seq_col,
      whole_seq,
      match_columns = c(genes_col),
      column_keep   = keep_cols
    )
  }, error=function(e) fail(paste("match_and_calculate_positions 失败：", e$message)))

  # 匹配为空：直接返回空结果（防崩）
  if (is.null(mr) || nrow(mr) == 0L) {
    out <- ok_empty(gene, metric=if (quantify=="AREA") "Area" else "PSM", agg=agg)
    if (nzchar(out_file)) jsonlite::write_json(out, out_file, auto_unbox=TRUE)
    cat(J(out)); quit(save="no")
  }

  # 仅在非空时安全加列
  if (nrow(mr) > 0L) {
    mr$reps <- NA_integer_
    mr$PTM_position <- NA
  }

  # ---- 定量与汇总 ----
  qa <- list(
    whole_seq        = whole_seq,
    matching_result  = mr,
    matching_columns = c(genes_col),                    # ← 修正为参数
    distinct_columns = unique(c(id_col, label_col)),
    quantify_method  = quantify,
    with_PTM         = FALSE,
    reps             = FALSE
  )
  if (quantify=="AREA" && area_col %in% names(mr)) qa$area_column <- area_col

  dws <- tryCatch({
    do.call(peptide_quantification, qa)
  }, error=function(e) fail(paste("peptide_quantification 失败：", e$message)))

  metric <- if (quantify=="AREA") "Area" else "PSM"
  if (!metric %in% names(dws)) {
    # 即便缺列也给出空结果而不是报错中断
    out <- ok_empty(gene, metric=metric, agg=agg)
    if (nzchar(out_file)) jsonlite::write_json(out, out_file, auto_unbox=TRUE)
    cat(J(out)); quit(save="no")
  }

  dt <- as.data.table(dws)
  # 容错：缺 id/label 时补空列，避免 data.table 报名不存在
  if (!id_col   %in% names(dt)) dt[, (id_col)   := ""]
  if (!label_col%in% names(dt)) dt[, (label_col):= ""]

  by_sample <- dt[, .(value = get(metric)),
                  by = .(Character, Position, id = get(id_col), label = get(label_col))]

  if (nrow(by_sample) == 0L) {
    out <- ok_empty(gene, metric=metric, agg=agg)
    if (nzchar(out_file)) jsonlite::write_json(out, out_file, auto_unbox=TRUE)
    cat(J(out)); quit(save="no")
  }

  by_label <- switch(
    agg,
    "mean"   = by_sample[, .(value = mean(value, na.rm=TRUE)),   by=.(Character, Position, label)],
    "sum"    = by_sample[, .(value = sum(value, na.rm=TRUE)),    by=.(Character, Position, label)],
    "median" = by_sample[, .(value = median(value, na.rm=TRUE)), by=.(Character, Position, label)],
    fail("agg 仅支持 mean/sum/median")
  )

  out <- list(status="ok", gene=gene, metric=metric, agg=agg,
              by_sample=as.data.frame(by_sample), by_label=as.data.frame(by_label))
  if (nzchar(out_file)) jsonlite::write_json(out, out_file, auto_unbox=TRUE)
  cat(J(out))
}

tryCatch(main(), error=function(e) fail(e))
