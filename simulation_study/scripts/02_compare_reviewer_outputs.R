#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Usage: Rscript scripts/02_compare_reviewer_outputs.R <result_dir>")
result_dir <- args[[1]]
source(file.path("/home/js/NESTA", "simulation_study", "R", "reproducibility_checks.R"))
checks <- verify_delta_workflow(result_dir)
cat(if (checks$passed) "PASS\n" else "FAIL\n")
