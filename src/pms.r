#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  suppressWarnings(library(PepMapViz))
  suppressWarnings(library(data.table))
  suppressWarnings(library(jsonlite))
})
try({Sys.setenv(LANG="C.UTF-8"); suppressWarnings(Sys.setlocale("LC_ALL","C.UTF-8"))}, silent=TRUE)
J <- function(x) toJSON(x, auto_unbox=TRUE, null="null", dataframe="columns")
fail <- function(msg){ cat(J(list(status="error", message=as.character(msg)))); quit(status=1) }

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

  whole_seq <- data.frame(Region_Sequence = seq_string, Genes = gene, stringsAsFactors = FALSE)

  if (is.null(X[["peptide_rows"]])) fail("需要 peptide_rows（肽段级表）。")
  pep_tbl <- as.data.frame(X[["peptide_rows"]], check.names=FALSE, stringsAsFactors=FALSE)
  need <- c(seq_col, genes_col, id_col, label_col)
  miss <- setdiff(need, names(pep_tbl)); if (length(miss)) fail(paste("peptide_rows 缺列:", paste(miss, collapse=", ")))

  keep_cols <- c(id_col, label_col)
  if (quantify=="AREA" && area_col %in% names(pep_tbl)) keep_cols <- c(keep_cols, area_col)

  # 关键：按位置传参（不要写 data= / whole_seq=）
  mr <- match_and_calculate_positions(
    pep_tbl,
    column = seq_col,
    whole_seq,
    match_columns = c(genes_col),
    column_keep   = keep_cols
  )
  mr$reps <- NA; mr$PTM_position <- NA

  qa <- list(
    whole_seq        = whole_seq,
    matching_result  = mr,
    matching_columns = c("Genes"),
    distinct_columns = c(id_col, label_col),
    quantify_method  = quantify,
    with_PTM         = FALSE,
    reps             = FALSE
  )
  if (quantify=="AREA" && area_col %in% names(mr)) qa$area_column <- area_col

  dws <- do.call(peptide_quantification, qa)
  metric <- if (quantify=="AREA") "Area" else "PSM"
  if (!metric %in% names(dws)) fail(paste("PepMapViz 返回中未找到列:", metric))

  dt <- as.data.table(dws)
  by_sample <- dt[, .(value = get(metric)),
                  by = .(Character, Position, id = get(id_col), label = get(label_col))]

  if (agg=="mean") {
    by_label <- by_sample[, .(value = mean(value, na.rm=TRUE)),  by=.(Character, Position, label)]
  } else if (agg=="sum") {
    by_label <- by_sample[, .(value = sum(value, na.rm=TRUE)),   by=.(Character, Position, label)]
  } else if (agg=="median") {
    by_label <- by_sample[, .(value = median(value, na.rm=TRUE)),by=.(Character, Position, label)]
  } else fail("agg 仅支持 mean/sum/median")

  cat(J(list(status="ok", gene=gene, metric=metric, agg=agg,
             by_sample=as.data.frame(by_sample), by_label=as.data.frame(by_label))))
}
tryCatch(main(), error=function(e) fail(e))
