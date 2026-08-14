#!/usr/bin/env Rscript
source(file.path("/home/js/NESTA/simulation", "R", "utils.R"))
dropbox_report <- read_report_dir()
report_dir <- project_file("results/reports")
safe_dir_create(report_dir)

impl <- c(list.files(project_file("R"), pattern = "[.]R$", full.names = TRUE),
          list.files(project_file("scripts"), pattern = "[.]R$", full.names = TRUE),
          project_file("config/FROZEN_CONFIG.yaml"))
impl <- impl[file.exists(impl)]
manifest <- data.frame(
  relative_path = sub(paste0(project_dir, "/"), "", impl),
  bytes = file.info(impl)$size,
  sha256 = vapply(impl, file_sha256, character(1)),
  stringsAsFactors = FALSE
)
for (dst_dir in c(report_dir, dropbox_report)) {
  p <- file.path(dst_dir, "IMPLEMENTATION_MANIFEST.csv")
  if (file.exists(p)) unlink(p)
  atomic_write_csv(manifest, p)
}

required <- c(
  "STUDY_PLAN_0705_OBSERVED_METRIC_ADAPTIVE_CALIBRATION.md",
  "FINAL_REPORT.md", "STOP_GO_REPORT.md", "RUN_STATUS.md", "CODE_FIDELITY_AUDIT.md",
  "NETWORK_GENERATOR_AUDIT.md", "PATH_STRATIFICATION_AUDIT.md",
  "DIRECTIONAL_SIGNAL_AUDIT.md", "DEGREE_DISTRIBUTION_AUDIT.md",
  "BENCHMARK_IMPLEMENTATION_AUDIT.md",
  "TARGET_INITIAL_SIGNAL_AUDIT.md", "CONTROL_DISRUPTION_AUDIT.md",
  "BRANCH_SPECIFICITY_AUDIT.md", "BRANCH_CONDUCTANCE_AUDIT.md",
  "RELAY_STRUCTURE_AUDIT.md",
  "CALIBRATION_CANDIDATE_AUDIT.csv", "CALIBRATION_ADAPTIVE_TRACE.md", "DIFFUSION_RETENTION_AUDIT.md",
  "CONTROL_SENSITIVITY_DIAGNOSIS.md",
  "NESTA_FINAL_HEAT_BENCHMARK_REPORT.md",
  "NESTA_DELTA_INTERPRETATION_REPORT.md",
  "NESTA_REPRIORITIZATION_REPORT.md",
  "NESTA_DEGREE_BIAS_REPORT.md",
  "RWR_PPR_DIRECTIONAL_BENCHMARK_REPORT.md",
  "DIFFICULTY_GRADIENT_REPORT.md",
  "MANUSCRIPT_READY_SUMMARY.md",
  "PRIMARY_FINAL_HEAT_METRICS.csv", "PRIMARY_FINAL_HEAT_CONTRASTS.csv",
  "DIRECTION_AWARE_METRICS.csv", "DIRECTION_AWARE_CONTRASTS.csv",
  "REPRIORITIZATION_METRICS.csv", "DELTA_INTERPRETATION_METRICS.csv",
  "DEGREE_BIAS_METRICS.csv", "DEGREE_BIAS_CONTRASTS.csv",
  "BRANCH_SPECIFICITY_AUDIT.csv", "BRANCH_CONDUCTANCE_AUDIT.csv",
  "RELAY_STRUCTURE_AUDIT.csv", "DIFFUSION_RETENTION_AUDIT.csv",
  "BENCHMARK_METRICS.csv", "BENCHMARK_CONTRASTS.csv",
  "BENCHMARK_DIRECTION_AWARE_METRICS.csv", "BENCHMARK_DEGREE_BIAS_METRICS.csv",
  "TARGET_INITIAL_SIGNAL_AUDIT.csv", "CONTROL_DISRUPTION_AUDIT.csv",
  "NULL_BIAS_GUARDRAILS.csv", "TEMPLATE_OUTLIER_AUDIT.csv",
  "TEMPLATE_EXCLUSION_AUDIT.csv",
  "CLEANUP_MANIFEST.csv", "IMPLEMENTATION_MANIFEST.csv", "unit_test_results.csv",
  file.path("figure_source_data", "panel_A_simulation_design_counts.csv"),
  file.path("figure_source_data", "panel_B_difficulty_prevalence.csv"),
  file.path("figure_source_data", "panel_C_final_heat_fold_enrichment.csv"),
  file.path("figure_source_data", "panel_D_sign_concordant_recovery.csv"),
  file.path("figure_source_data", "panel_E_decoy_and_degree_bias.csv"),
  file.path("figure_source_data", "panel_F_reprioritized_blank_targets.csv")
)
missing <- required[!file.exists(file.path(dropbox_report, required))]
if (length(missing)) stop("Missing required report files: ", paste(missing, collapse = ", "))

checks <- list.files(dropbox_report, recursive = TRUE, full.names = TRUE)
checks <- checks[file.info(checks)$isdir == FALSE]
checks <- checks[!basename(checks) %in% c("CHECKSUMS.sha256", "EXPORT_CHECKSUMS.sha256")]
lines <- vapply(sort(checks), function(f) paste(file_sha256(f), sub(paste0(dropbox_report, "/"), "", f)), character(1))
for (dst in c(file.path(dropbox_report, "CHECKSUMS.sha256"), file.path(report_dir, "CHECKSUMS.sha256"))) {
  if (file.exists(dst)) unlink(dst)
  atomic_write_lines(lines, dst)
}
cat("Checksum manifest written for ", length(lines), " files\n", sep = "")
