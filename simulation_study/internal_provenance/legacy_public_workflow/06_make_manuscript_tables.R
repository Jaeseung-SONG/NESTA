#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
report_dir <- Sys.getenv("NESTA_REPORT_DIR", "")
if ("--report-dir" %in% args) report_dir <- args[[match("--report-dir", args) + 1]]
if (!nzchar(report_dir) && length(args) == 1) report_dir <- args[[1]]
if (!nzchar(report_dir)) report_dir <- "/tmp/nesta_simulation_study_report"
dir.create(file.path(report_dir, "manuscript_tables"), recursive=TRUE, showWarnings=FALSE)
copy_if <- function(src, dst) if (file.exists(src)) file.copy(src, dst, overwrite=TRUE)
copy_if(file.path(report_dir, "SCENARIO_SUMMARY_TABLE.csv"), file.path(report_dir, "manuscript_tables", "table_simulation_scenario_summary.csv"))
copy_if(file.path(report_dir, "REFERENCE_RESULT_COMPARISON.csv"), file.path(report_dir, "manuscript_tables", "table_reproducibility_checks.csv"))
cat("manuscript tables written\n")
