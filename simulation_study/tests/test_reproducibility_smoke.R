#!/usr/bin/env Rscript

root <- "/home/js/NESTA/simulation_study"
source(file.path(root, "R", "reproducibility_checks.R"))

fidelity <- code_fidelity_audit(root)
stopifnot(all(fidelity$passed))

required_scripts <- file.path(root, "scripts", c(
  "00_prepare.R",
  "01_run_delta_threshold_sensitivity.R",
  "02_compare_reviewer_outputs.R",
  "03_make_manuscript_tables.R",
  "08_per_gene_export_delta_sensitivity.R"
))
stopifnot(all(file.exists(required_scripts)))

legacy_public <- file.path(root, "internal_provenance", "legacy_public_workflow")
stopifnot(dir.exists(legacy_public))
stopifnot(!file.exists(file.path(root, "R", "run_comparators.R")))

config <- readLines(file.path(root, "config", "final_simulation_study.yaml"), warn = FALSE)
stopifnot(any(grepl("primary_score: Final_Heat", config, fixed = TRUE)))
stopifnot(any(grepl("complementary_score: delta_NESTA", config, fixed = TRUE)))
stopifnot(any(grepl("top_1pct_framing: conservative operating point", config, fixed = TRUE)))

expected <- final_metric_expectations()
stopifnot(nrow(expected) == 10)
threshold_expected <- threshold_compact_expectations()
stopifnot(all(c("Final Heat only", "delta_NESTA only",
                "Final OR delta_NESTA", "Final AND delta_NESTA") %in%
                threshold_expected$selection_rule))

out_files <- list.files(file.path(root, "outputs"), all.files = FALSE, no.. = TRUE)
stopifnot(setequal(out_files, "README_KEEP_EMPTY.md"))

cat("simulation_study smoke test passed\n")
