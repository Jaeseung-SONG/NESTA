#!/usr/bin/env Rscript

study_runner <- file.path("/home/js/NESTA", "simulation_study", "run_simulation_study.R")
args <- commandArgs(trailingOnly = TRUE)
status <- system2("Rscript", c(study_runner, args))
if (!identical(status, 0L)) quit(status = status)
