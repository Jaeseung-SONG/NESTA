#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Usage: Rscript scripts/03_make_manuscript_tables.R <result_dir>")
result_dir <- args[[1]]
src <- file.path(result_dir, "manuscript_tables", "table_delta_threshold_sensitivity_compact.tsv")
if (!file.exists(src)) stop("Missing compact table: ", src)
cat("Manuscript-ready compact threshold table:\n", src, "\n", sep = "")
