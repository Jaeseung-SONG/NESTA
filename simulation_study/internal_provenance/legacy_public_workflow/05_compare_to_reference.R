#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
report_dir <- if (length(args)) args[[1]] else Sys.getenv("NESTA_REPORT_DIR", "/tmp/nesta_simulation_study_report")
dir.create(report_dir, recursive=TRUE, showWarnings=FALSE)
source(file.path("/home/js/NESTA/simulation_study", "R", "reproducibility_checks.R"))
res <- verify_references(report_dir)
cat(if (res$passed) "PASS\n" else "FAIL\n")
