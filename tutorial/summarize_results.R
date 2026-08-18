#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if(length(args) >= 1) args[[1]] else file.path("tutorial", "output")
suppressPackageStartupMessages(library(data.table))

read_scores <- function(cell_type){
  path <- file.path(output_dir, paste0("cell_type_", cell_type, "_scores.tsv"))
  if(!file.exists(path)) stop("Missing tutorial output: ", path)
  tab <- data.table::fread(path)
  tab[, cell_type := paste0("Cell_type_", cell_type)]
  tab
}

scores <- data.table::rbindlist(list(read_scores("A"), read_scores("B")), fill = TRUE)
scores[, abs_Final_Heat := abs(Final.Heat)]
data.table::setorder(scores, cell_type, -abs_Final_Heat, SYMBOL)
scores[, rank_by_Final_Heat := seq_len(.N), by = cell_type]

comparison <- scores[, .(
  cell_type,
  SYMBOL,
  TWAS.Z,
  Mean_expression,
  Initial.Heat,
  Final.Heat,
  delta_NESTA,
  rank_by_Final_Heat
)]

data.table::fwrite(comparison,
                   file.path(output_dir, "cell_type_ranking_comparison.tsv"),
                   sep = "\t")

cat("\nTop five genes by absolute Final Heat\n")
print(comparison[rank_by_Final_Heat <= 5])
