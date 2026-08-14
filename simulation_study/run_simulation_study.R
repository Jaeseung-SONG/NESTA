#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = "") {
  if (flag %in% args) {
    idx <- match(flag, args)
    if (idx < length(args)) return(args[[idx + 1]])
  }
  default
}

repo_root <- "/home/js/NESTA"
study_dir <- file.path(repo_root, "simulation_study")
timestamp <- arg_value("--timestamp", Sys.getenv("NESTA_DELTA_EXPORT_TIMESTAMP", ""))
if (!nzchar(timestamp)) timestamp <- format(Sys.time(), "%d%m%y_%H%M")
report_dir <- arg_value("--report-dir", file.path(repo_root, "simulation_study_results",
                                                  paste0("reviewer_delta_threshold_sensitivity_", timestamp)))
sync_dropbox <- "--sync-dropbox" %in% args || identical(Sys.getenv("NESTA_SYNC_DROPBOX"), "1")

cmd_args <- c(file.path(study_dir, "scripts", "08_per_gene_export_delta_sensitivity.R"),
              "--timestamp", timestamp,
              "--output-dir", report_dir)
if (sync_dropbox) cmd_args <- c(cmd_args, "--sync-dropbox")

status <- system2("Rscript", cmd_args)
if (!identical(status, 0L)) stop("Delta threshold sensitivity workflow failed with exit status ", status)

source(file.path(study_dir, "R", "reproducibility_checks.R"))
checks <- verify_delta_workflow(report_dir)
fidelity <- code_fidelity_audit(study_dir)
utils::write.table(fidelity, file.path(report_dir, "code_fidelity_audit.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)

compact <- read_tsv_checked(file.path(report_dir, "manuscript_tables",
                                      "table_delta_threshold_sensitivity_compact.tsv"))
preview <- compact[
  compact$threshold_label %in% c("top_0_5pct", "top_1pct", "top_2pct", "top_5pct") &
    compact$selection_rule %in% c("Final Heat only", "delta_NESTA only",
                                  "Final OR delta_NESTA", "Final AND delta_NESTA"),
  c("threshold_label", "selection_rule", "selected_gene_count_mean",
    "true_target_recall_mean", "FPR_mean",
    "risk_target_recall_mean", "protective_target_recall_mean",
    "opposite_sign_decoy_selection_rate_mean",
    "high_score_decoy_selection_rate_mean"),
  drop = FALSE
]

writeLines(c(
  "# Reviewer-facing Simulation Workflow Report",
  "",
  paste("Output directory:", report_dir),
  "",
  sprintf("Final aggregate reproducibility checks passed: %d/%d.",
          sum(checks$final_comparison$required_passed), nrow(checks$final_comparison)),
  sprintf("Compact threshold table checks passed: %d/%d.",
          sum(checks$compact_comparison$passed), nrow(checks$compact_comparison)),
  sprintf("Code fidelity checks passed: %d/%d.", sum(fidelity$passed), nrow(fidelity)),
  "",
  "Primary public outputs:",
  "- `DELTA_THRESHOLD_SENSITIVITY_REPORT.md`",
  "- `delta_threshold_sensitivity_summary.tsv`",
  "- `delta_threshold_sensitivity_by_scenario.tsv`",
  "- `delta_threshold_sensitivity_directional.tsv`",
  "- `delta_threshold_sensitivity_decoy_control.tsv`",
  "- `delta_threshold_sensitivity_selected_count_summary.tsv`",
  "- `manuscript_tables/table_delta_threshold_sensitivity.tsv`",
  "- `manuscript_tables/table_delta_threshold_sensitivity_compact.tsv`",
  "",
  "The top 1% rule is reported as a conservative operating point. Relaxed thresholds are sensitivity analyses for the recovery/specificity tradeoff."
), file.path(report_dir, "SIMULATION_WORKFLOW_REPORT.md"))

utils::write.table(preview, file.path(report_dir, "compact_table_preview.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)

if (!checks$passed || !all(fidelity$passed)) {
  stop("Reviewer-facing workflow completed but verification checks failed. See ", report_dir)
}

cat(report_dir, "\n")
