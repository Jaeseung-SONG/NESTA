#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: 07_build_final_reports.R <config.yaml>")
source(file.path(dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1]))), "lib_config.R"))
suppressPackageStartupMessages(library(data.table))
cfg <- read_simple_config(args[[1]])
out <- cfg_get(cfg, "output_dir")
tab <- fread(file.path(out, "tables", "negative_control_46_irnt_gene_level_q99_overlap.tsv"))
known_col <- if ("result2_known_marker_count_total" %in% names(tab)) "result2_known_marker_count_total" else "result2_known_marker_count"
report <- c(
  "# Negative-Control GWAS Final Report",
  "",
  "Phenotype: `46_irnt` - Hand grip strength, left.",
  "",
  sprintf("Gene-level selected genes: %s.", tab$selected_gene_count_q99_union[1]),
  sprintf("Result 2 known GD marker overlap: %s/%s.", tab$result2_overlap_count[1], tab[[known_col]][1]),
  sprintf("Enrichment relative to available universe: %.6f.", tab$enrichment_relative_to_available_gene_universe[1]),
  sprintf("Hypergeometric p-value: %.6g.", tab$hypergeometric_p_value[1])
)
writeLines(report, file.path(out, "reports", "negative_control_46_irnt_summary_report.md"))
