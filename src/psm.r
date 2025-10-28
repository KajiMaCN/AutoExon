#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  suppressWarnings(library(PepMapViz))
  suppressWarnings(library(data.table))
  suppressWarnings(library(jsonlite))
  suppressWarnings(library(stringi))
})

# 尽量干净的本地化
try({Sys.setenv(LANG="C.UTF-8"); suppressWarnings(Sys.setlocale("LC_ALL","C.UTF-8"))}, silent=TRUE)

`%||%` <- function(a,b) if (is.null(a)) b else a
J <- function(x) jsonlite::toJSON(x, auto_unbox=TRUE, null="null", dataframe="columns")
fail <- function(msg){ cat(J(list(status="error", message=as.character(msg)))); quit(status=1) }

# I/J -> L，并转大写（可选）
IL <- function(x) chartr("IJ","LL", toupper(x %||% ""))

.run_one_isoform <- function(gene, seq_string, pep_tbl, seq_col, genes_col, column_keep, normalizeIL) {
  if (!nchar(seq_string)) return(data.frame())

  if (normalizeIL) {
    seq_il <- IL(seq_string)
    pep_tbl[[seq_col]] <- IL(pep_tbl[[seq_col]])
  } else {
    seq_il <- seq_string
  }

  whole_seq <- data.frame(
    Region_Sequence = seq_il,
    Genes = gene %||% "",
    stringsAsFactors = FALSE
  )

  mr <- tryCatch({
    match_and_calculate_positions(
      pep_tbl,
      column        = seq_col,
      whole_seq     = whole_seq,
      match_columns = c(genes_col),
      column_keep   = column_keep
    )
  }, error=function(e) fail(paste("match_and_calculate_positions 调用失败：", e$message)))

  if (is.null(mr)) mr <- data.frame()
  as.data.frame(mr, check.names=FALSE, stringsAsFactors=FALSE)
}

.handle_one_gene <- function(item, top_out_pep_dir="") {
  gene        <- item$gene %||% ""
  seq_col     <- item$seq_col %||% "Stripped.Sequence"
  genes_col   <- item$genes_col %||% "Genes"
  column_keep <- item$column_keep %||% c("Polypeptide.Novogene.ID","Global.Q.Value","Precursor.Normalised")
  normalizeIL <- isTRUE(item$normalize_IL %||% TRUE)
  pep_tbl_in  <- item$peptide_rows
  isoforms    <- item$isoforms %||% list()
  out_dir     <- item$out_dir %||% top_out_pep_dir

  if (is.null(pep_tbl_in)) return(list(gene=gene, result=data.frame()))
  pep_tbl <- tryCatch(as.data.frame(pep_tbl_in, check.names=FALSE, stringsAsFactors=FALSE),
                      error=function(e) fail(paste("peptide_rows 需要可转为 data.frame：", e$message)))
  if (!nrow(pep_tbl)) return(list(gene=gene, result=data.frame()))
  if (!seq_col %in% names(pep_tbl))  fail(paste0("找不到序列列 seq_col=`", seq_col, "`"))
  if (!genes_col %in% names(pep_tbl)) fail(paste0("找不到匹配列 genes_col=`", genes_col, "`"))

  pep_tbl[[seq_col]] <- trimws(as.character(pep_tbl[[seq_col]]))
  pep_tbl <- pep_tbl[nzchar(pep_tbl[[seq_col]]), , drop=FALSE]
  if (!nrow(pep_tbl)) return(list(gene=gene, result=data.frame()))

  res_list <- list()
  if (length(isoforms) == 0) {
    # 兼容单序列模式
    seq_str <- item$sequence %||% ""
    df <- .run_one_isoform(gene, seq_str, pep_tbl, seq_col, genes_col, column_keep, normalizeIL)
    if (nrow(df)) df$isoform <- item$isoform %||% (item$accession %||% NA_character_)
    res_list[[length(res_list)+1]] <- df
  } else {
    for (it in isoforms) {
      seq_str <- it$sequence %||% ""
      iso_acc <- it$isoform  %||% NA_character_
      df <- .run_one_isoform(gene, seq_str, pep_tbl, seq_col, genes_col, column_keep, normalizeIL)
      if (nrow(df)) { df$isoform <- iso_acc; res_list[[length(res_list)+1]] <- df }
    }
  }

  out_df <- if (length(res_list)) data.table::rbindlist(res_list, fill=TRUE) else data.frame()

  # 落盘：{out_dir}/{gene}/pepmapviz_positions.csv
  if (nchar(out_dir)) {
    gene_dir <- file.path(out_dir, gene)
    try(dir.create(gene_dir, recursive=TRUE, showWarnings=FALSE), silent=TRUE)
    try(data.table::fwrite(out_df, file.path(gene_dir, "pepmapviz_positions.csv")), silent=TRUE)
    # 可选写 JSON 便于排查
    safe <- list(status="ok", gene=gene, n_match=nrow(out_df), result=out_df)
    try(jsonlite::write_json(safe, file.path(gene_dir, "result.json"), auto_unbox=TRUE, null="null"), silent=TRUE)
  }

  list(gene=gene, result=out_df)
}

main <- function(){
  X <- tryCatch({
    fromJSON(paste(readLines(file("stdin"), warn=FALSE), collapse="\n"), simplifyVector=TRUE)
  }, error=function(e) fail(paste("无法解析输入 JSON：", e$message)))

  # --------- 批模式（推荐） ---------
  if (!is.null(X$batch)) {
    top_out <- X$out_pep_dir %||% ""
    items   <- X$batch
    if (length(items) == 0) { cat(J(list(status="ok", batch_results=list()))); return(invisible(NULL)) }

    results <- lapply(items, function(it) {
      r <- .handle_one_gene(it, top_out_pep_dir = top_out)
      list(gene = r$gene, result = r$result)
    })

    cat(J(list(status="ok", batch_results=results)))
    return(invisible(NULL))
  }

  # --------- 单基因模式（兼容） ---------
  req <- function(nm){ if (is.null(X[[nm]])) fail(paste("缺少字段:", nm)); X[[nm]] }
  opt <- function(nm, dflt) X[[nm]] %||% dflt

  gene        <- opt("gene", "")
  seq_string  <- req("sequence")
  pep_tbl_in  <- req("peptide_rows")
  seq_col     <- opt("seq_col",   "Stripped.Sequence")
  genes_col   <- opt("genes_col", "Genes")
  column_keep <- opt("column_keep",
                     c("Polypeptide.Novogene.ID","Global.Q.Value","Precursor.Normalised"))
  normalizeIL <- isTRUE(opt("normalize_IL", TRUE))
  out_file    <- opt("out_file","")
  out_dir     <- opt("out_dir","")

  pep_tbl <- tryCatch(as.data.frame(pep_tbl_in, check.names=FALSE, stringsAsFactors=FALSE),
                      error=function(e) fail(paste("peptide_rows 需要是可转为 data.frame 的对象：", e$message)))
  if (!nchar(seq_string)) fail("sequence 不能为空。")
  if (!nrow(pep_tbl))     fail("peptide_rows 为空。")
  if (!seq_col %in% names(pep_tbl))  fail(paste0("找不到序列列 seq_col=`", seq_col, "`"))
  if (!genes_col %in% names(pep_tbl)) fail(paste0("找不到匹配列 genes_col=`", genes_col, "`"))

  pep_tbl[[seq_col]] <- trimws(as.character(pep_tbl[[seq_col]]))
  pep_tbl <- pep_tbl[nzchar(pep_tbl[[seq_col]]), , drop=FALSE]
  if (!nrow(pep_tbl)) fail("清洗后 peptide_rows 为空。")

  if (normalizeIL) {
    seq_il <- IL(seq_string)
    pep_tbl[[seq_col]] <- IL(pep_tbl[[seq_col]])
  } else {
    seq_il <- seq_string
  }

  whole_seq <- data.frame(
    Region_Sequence = seq_il,
    Genes = gene %||% "",
    stringsAsFactors = FALSE
  )

  mr <- tryCatch({
    match_and_calculate_positions(
      pep_tbl,
      column        = seq_col,
      whole_seq     = whole_seq,
      match_columns = c(genes_col),
      column_keep   = column_keep
    )
  }, error=function(e) fail(paste("match_and_calculate_positions 调用失败：", e$message)))

  if (is.null(mr)) {
    mr <- data.frame()
  }

  # 落盘：{out_dir}/{gene}/pepmapviz_positions.csv
  if (nchar(out_dir)) {
    gene_dir <- file.path(out_dir, gene)
    try(dir.create(gene_dir, recursive=TRUE, showWarnings=FALSE), silent=TRUE)
    try(data.table::fwrite(mr, file.path(gene_dir, "pepmapviz_positions.csv")), silent=TRUE)
    safe <- list(status="ok", gene=gene, n_match=nrow(mr), result=mr)
    try(jsonlite::write_json(safe, file.path(gene_dir, "result.json"), auto_unbox=TRUE, null="null"), silent=TRUE)
  }

  out <- list(
    status  = "ok",
    gene    = gene,
    n_match = if (!is.null(nrow(mr))) nrow(mr) else 0L,
    result  = as.data.frame(mr, check.names=FALSE, stringsAsFactors=FALSE)
  )

  if (nzchar(out_file)) {
    jsonlite::write_json(out, out_file, auto_unbox=TRUE, null="null")
  }
  cat(J(out))
}

tryCatch(main(), error=function(e) fail(e))
