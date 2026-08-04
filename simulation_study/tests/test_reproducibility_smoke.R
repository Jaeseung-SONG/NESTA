#!/usr/bin/env Rscript
root <- "/home/js/NESTA/simulation_study"
source(file.path(root, "R", "reproducibility_checks.R"))
fidelity <- code_fidelity_audit()
stopifnot(all(fidelity$passed))
refs <- read.csv(file.path(root, "reference_manifest", "reference_result_manifest.csv"))
stopifnot(nrow(refs) == 4, all(dir.exists(refs$reference_dir)))
config <- readLines(file.path(root, "config", "final_simulation_study.yaml"), warn = FALSE)
stopifnot(any(grepl("data_dir: /home/js/NESTA/simulation_study_data", config, fixed = TRUE)))
stopifnot(any(grepl("local_results_dir: /home/js/NESTA/simulation_study_results", config, fixed = TRUE)))
out_files <- list.files(file.path(root, "outputs"), all.files = FALSE, no.. = TRUE)
stopifnot(setequal(out_files, "README_KEEP_EMPTY.md"))
cat("simulation_study smoke test passed\n")
