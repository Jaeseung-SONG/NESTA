Sys.setenv(NESTA_BIDIRECTIONAL_SOURCE_ONLY = "1")
source(file.path("/home/js/NESTA/simulation", "R", "study_0708_bidirectional_rescue.R"))

primary_decision_arm <- "P_combined_decoy_suppressed"
primary_rankings <- c("NESTA_two_tail_balanced", "NESTA_two_tail_direction_matched")

decision_repair_stage <- function(topology) {
  run_bidirectional_stage(topology, primary_decision_arm)
}

benchmark_decision_summary <- function(benchmarks) {
  stats::aggregate(cbind(top100_recall, top150_recall, top200_recall,
                         raw_AUPRC, direction_aware_AUPRC,
                         sign_concordant_top100_recall,
                         opposite_sign_decoy_top100_rate,
                         high_degree_decoy_top100_rate,
                         score_degree_spearman) ~ topology_arm + topology_label + rescue_arm + score_name,
                   data = benchmarks, FUN = function(z) mean(z, na.rm = TRUE))
}

decision_summary <- function(metrics) {
  bidirectional_success_summary(metrics)
}

primary_f_pass <- function(summary, audit) {
  rows <- summary[summary$topology_arm == "F" &
                    summary$rescue_arm == primary_decision_arm &
                    summary$ranking_mode %in% primary_rankings, , drop = FALSE]
  raw_ok <- all(audit$raw_TWAS_top100_A2_fraction[audit$topology_arm == "F"] <= 0.05) &&
    all(audit$M2_top100_A2_fraction[audit$topology_arm == "F"] <= 0.10)
  rows$primary_F_pass <- rows$top100_recall >= 0.70 &
    rows$top150_recall >= 0.90 &
    rows$top200_recall >= 0.95 &
    rows$risk_top100_recall >= 0.60 &
    rows$protective_top100_recall >= 0.60 &
    rows$opposite_sign_decoy_top100_rate <= 0.10 &
    rows$high_degree_decoy_top100_rate <= 0.10 &
    raw_ok
  rows
}

robust_h_pass <- function(summary) {
  rows <- summary[summary$topology_arm == "H" &
                    summary$rescue_arm == primary_decision_arm &
                    summary$ranking_mode %in% primary_rankings, , drop = FALSE]
  if (!nrow(rows)) return(rows)
  rows$robustness_H_pass <- rows$top100_recall >= 0.60 &
    rows$top150_recall >= 0.80 &
    rows$risk_top100_recall >= 0.50 &
    rows$protective_top100_recall >= 0.50 &
    rows$opposite_sign_decoy_top100_rate <= 0.15 &
    rows$high_degree_decoy_top100_rate <= 0.15
  rows
}

saturation_rows <- function(summary) {
  out <- summary[, c("topology_arm", "topology_label", "rescue_arm", "ranking_mode",
                     "top150_recall", "top200_recall",
                     "opposite_sign_decoy_top150_rate",
                     "opposite_sign_decoy_top200_rate",
                     "high_degree_decoy_top150_rate",
                     "high_degree_decoy_top200_rate",
                     "score_degree_spearman",
                     "score_strength_spearman"), drop = FALSE]
  out$interpretation <- "Exploratory saturation diagnostic; not a hard stop after top150/top200 A2 saturation."
  out
}

network_path_audit <- function(audit) {
  audit[, c("topology_arm", "topology_label", "rescue_arm", "replicate",
            "mean_seed_to_A2_weight", "mean_seed_to_opposite_decoy_weight",
            "mean_seed_to_high_degree_decoy_weight",
            "median_seed_to_A2_weight", "median_seed_to_opposite_decoy_weight",
            "median_seed_to_high_degree_decoy_weight"), drop = FALSE]
}

run_decision_rule_repair <- function() {
  verify_project_path()
  verify_binding_plan()
  report_dir <- read_report_dir()
  safe_dir_create(report_dir)
  copy_binding_plan_to_report(report_dir)
  set_run_status("IN_PROGRESS", "YES", "NO", "0708 decision-rule repair robustness")

  f <- decision_repair_stage("F")
  metrics <- f$metrics
  benchmarks <- f$benchmarks
  audit <- f$audit
  summary <- decision_summary(metrics)
  f_decision <- primary_f_pass(summary, audit)
  f_pass <- any(f_decision$primary_F_pass, na.rm = TRUE)
  h <- NULL
  if (f_pass) {
    h <- decision_repair_stage("H")
    metrics <- rbind(metrics, h$metrics)
    benchmarks <- rbind(benchmarks, h$benchmarks)
    audit <- rbind(audit, h$audit)
    summary <- decision_summary(metrics)
  }
  h_decision <- robust_h_pass(summary)
  h_pass <- if (is.null(h)) FALSE else any(h_decision$robustness_H_pass, na.rm = TRUE)
  f_decision <- primary_f_pass(summary, audit)
  best_f <- f_decision[order(f_decision$primary_F_pass, f_decision$top100_recall,
                             f_decision$top150_recall, f_decision$top200_recall,
                             decreasing = TRUE), ][1, , drop = FALSE]
  if (f_pass && (is.null(h) || h_pass)) {
    reason <- "decision_rule_repair_passed_confirmatory_ready"
  } else if (f_pass && !h_pass) {
    reason <- "F_passed_H_failed_robustness"
  } else if (any(f_decision$top100_recall < 0.70 | f_decision$top150_recall < 0.90 |
                 f_decision$top200_recall < 0.95, na.rm = TRUE)) {
    reason <- "bidirectional_recovery_failed"
  } else if (any(f_decision$opposite_sign_decoy_top100_rate > 0.10 |
                 f_decision$high_degree_decoy_top100_rate > 0.10, na.rm = TRUE)) {
    reason <- "top100_decoy_guardrail_failed"
  } else if (any(audit$raw_TWAS_top100_A2_fraction > 0.05 |
                 audit$M2_top100_A2_fraction > 0.10, na.rm = TRUE)) {
    reason <- "raw_or_no_diffusion_baseline_too_strong"
  } else {
    reason <- "implementation_error"
  }
  set_run_status("STOPPED", "YES", "NO", reason)

  bench_summary <- benchmark_decision_summary(benchmarks)
  decoy <- summary[, c("topology_arm", "topology_label", "rescue_arm", "ranking_mode",
                       "opposite_sign_decoy_top100_rate", "opposite_sign_decoy_top150_rate",
                       "opposite_sign_decoy_top200_rate", "high_degree_decoy_top100_rate",
                       "high_degree_decoy_top150_rate", "high_degree_decoy_top200_rate",
                       "score_degree_spearman"), drop = FALSE]
  direction <- summary[, c("topology_arm", "topology_label", "rescue_arm", "ranking_mode",
                           "risk_top100_recall", "protective_top100_recall",
                           "risk_top150_recall", "protective_top150_recall",
                           "risk_top200_recall", "protective_top200_recall",
                           "direction_aware_AUPRC", "sign_concordant_top100_recall"), drop = FALSE]
  schema <- data.frame(file = c("PRIMARY_FINAL_HEAT_METRICS.csv", "SATURATION_DIAGNOSTICS.csv",
                                "BENCHMARK_METRICS.csv"),
                       required_columns_present = c(all(c("top100_recall", "risk_top100_recall",
                                                          "protective_top100_recall") %in% names(metrics)),
                                                    all(c("opposite_sign_decoy_top200_rate",
                                                          "high_degree_decoy_top200_rate") %in% names(summary)),
                                                    all(c("direction_aware_AUPRC",
                                                          "sign_concordant_top100_recall") %in% names(benchmarks))),
                       stringsAsFactors = FALSE)

  write_csv_over(metrics, file.path(report_dir, "PRIMARY_FINAL_HEAT_METRICS.csv"))
  write_csv_over(metrics, file.path(report_dir, "BIDIRECTIONAL_RANKING_METRICS.csv"))
  write_csv_over(direction, file.path(report_dir, "DIRECTION_SPECIFIC_RECOVERY_METRICS.csv"))
  write_csv_over(decoy, file.path(report_dir, "DECOY_GUARDRAIL_METRICS.csv"))
  write_csv_over(saturation_rows(summary), file.path(report_dir, "SATURATION_DIAGNOSTICS.csv"))
  write_csv_over(benchmarks, file.path(report_dir, "BENCHMARK_METRICS.csv"))
  write_csv_over(bench_summary, file.path(report_dir, "BENCHMARK_CONTRASTS.csv"))
  write_csv_over(audit, file.path(report_dir, "INITIAL_SIGNAL_FIELD_AUDIT.csv"))
  write_csv_over(network_path_audit(audit), file.path(report_dir, "NETWORK_PATH_AUDIT.csv"))
  write_csv_over(schema, file.path(report_dir, "METRIC_SCHEMA_AUDIT.csv"))

  write_lines_over(c("# Decision Rule Repair Summary", "",
                     paste0("Outcome classification: `", reason, "`."),
                     sprintf("Primary F ranking `%s`: top100/top150/top200 %.4f / %.4f / %.4f.",
                             best_f$ranking_mode, best_f$top100_recall, best_f$top150_recall, best_f$top200_recall),
                     sprintf("Risk/protective top100 %.4f / %.4f.",
                             best_f$risk_top100_recall, best_f$protective_top100_recall),
                     sprintf("Top100 opposite/high-degree decoy %.4f / %.4f.",
                             best_f$opposite_sign_decoy_top100_rate, best_f$high_degree_decoy_top100_rate),
                     sprintf("Top200 opposite/high-degree decoy saturation diagnostics %.4f / %.4f.",
                             best_f$opposite_sign_decoy_top200_rate, best_f$high_degree_decoy_top200_rate),
                     "Top200 decoy metrics are reported as saturation diagnostics, not hard stop criteria."),
                   file.path(report_dir, "DECISION_RULE_REPAIR_SUMMARY.md"))
  if (is.null(h)) {
    h_lines <- c("# Robustness Topology H Summary", "", "Topology H robustness was not run because topology F did not pass the repaired primary gate.")
  } else {
    h_best <- h_decision[order(h_decision$robustness_H_pass, h_decision$top100_recall,
                               h_decision$top150_recall, decreasing = TRUE), ][1, ]
    h_lines <- c("# Robustness Topology H Summary", "",
                 sprintf("Topology H robustness pass: %s.", h_pass),
                 sprintf("Best H ranking `%s`: top100/top150 %.4f / %.4f.",
                         h_best$ranking_mode, h_best$top100_recall, h_best$top150_recall),
                 sprintf("Risk/protective top100 %.4f / %.4f; decoy top100 opposite/high-degree %.4f / %.4f.",
                         h_best$risk_top100_recall, h_best$protective_top100_recall,
                         h_best$opposite_sign_decoy_top100_rate, h_best$high_degree_decoy_top100_rate))
  }
  write_lines_over(h_lines, file.path(report_dir, "ROBUSTNESS_TOPOLOGY_H_SUMMARY.md"))
  write_lines_over(c("# Code Fidelity Audit", "",
                     paste0("Binding plan SHA256: `", binding_plan_sha256, "`."),
                     "Faithful submitted M2 arithmetic, TWAS.P conversion, strict filtering, signed positive/absolute-negative diffusion, zero-weight edge behavior, self-loop behavior, and diffuStats `n.perm` were retained.",
                     "Submitted Final Heat values were not modified; this round repaired the decision rule and robustness evaluation."),
                   file.path(report_dir, "CODE_FIDELITY_AUDIT.md"))
  write_lines_over(c("# Metric Schema Audit", "",
                     sprintf("Primary metrics schema pass: %s.", schema$required_columns_present[1]),
                     sprintf("Saturation diagnostics schema pass: %s.", schema$required_columns_present[2]),
                     sprintf("Benchmark metrics schema pass: %s.", schema$required_columns_present[3])),
                   file.path(report_dir, "METRIC_SCHEMA_AUDIT.md"))
  write_lines_over(c("# STOP/GO Report", "", "STOP.",
                     paste0("Reason: `", reason, "`."),
                     "Diagnostic stage started: YES.",
                     "Confirmatory started: NO.",
                     "Top200 decoy rates were treated as exploratory saturation diagnostics, not hard stop criteria."),
                   file.path(report_dir, "STOP_GO_REPORT.md"))
  write_lines_over(c("# Final Report", "",
                     paste0("Binding plan SHA256: `", binding_plan_sha256, "`."),
                     paste0("Final outcome classification: `", reason, "`."),
                     sprintf("Primary F `%s` with `%s`: top100/top150/top200 %.4f / %.4f / %.4f.",
                             primary_decision_arm, best_f$ranking_mode,
                             best_f$top100_recall, best_f$top150_recall, best_f$top200_recall),
                     sprintf("Risk/protective top100 %.4f / %.4f; top100 decoy opposite/high-degree %.4f / %.4f.",
                             best_f$risk_top100_recall, best_f$protective_top100_recall,
                             best_f$opposite_sign_decoy_top100_rate, best_f$high_degree_decoy_top100_rate),
                     sprintf("Topology H robustness executed: %s; H pass: %s.", !is.null(h), h_pass),
                     "Top200 decoy rates are saturation diagnostics because A2 recovery is already near-complete by top150/top200.",
                     "Confirmatory execution started: NO."),
                   file.path(report_dir, "FINAL_REPORT.md"))
  if (file.exists(project_file("results/reports/summary_tables/unit_test_results.csv"))) {
    file.copy(project_file("results/reports/summary_tables/unit_test_results.csv"),
              file.path(report_dir, "unit_test_results.csv"), overwrite = TRUE)
  }
  if (file.exists(project_file("CLEANUP_MANIFEST.csv"))) {
    file.copy(project_file("CLEANUP_MANIFEST.csv"), file.path(report_dir, "CLEANUP_MANIFEST.csv"), overwrite = TRUE)
  }
  write_csv_over(data.frame(path = c(project_file("R/study_0708_decision_rule_repair.R"),
                                     project_file("R/study_0708_bidirectional_rescue.R"),
                                     project_file("R/study_0708_initial_signal_rescue.R"),
                                     project_file("R/fidelity.R"), project_file("R/utils.R"),
                                     project_file("scripts/run_decision_rule_repair.R")),
                            sha256 = c(sha(project_file("R/study_0708_decision_rule_repair.R")),
                                       sha(project_file("R/study_0708_bidirectional_rescue.R")),
                                       sha(project_file("R/study_0708_initial_signal_rescue.R")),
                                       sha(project_file("R/fidelity.R")),
                                       sha(project_file("R/utils.R")),
                                       if (file.exists(project_file("scripts/run_decision_rule_repair.R"))) sha(project_file("scripts/run_decision_rule_repair.R")) else NA_character_),
                            role = c("decision_rule_repair_runner", "bidirectional_ranking_helpers",
                                     "initial_signal_carryforward", "faithful_nesta_and_benchmarks",
                                     "binding_plan_and_io_guards", "script_entrypoint"),
                            stringsAsFactors = FALSE),
                 file.path(report_dir, "IMPLEMENTATION_MANIFEST.csv"))
  files <- list.files(report_dir, recursive = TRUE, full.names = TRUE)
  files <- files[file.info(files)$isdir == FALSE & basename(files) != "CHECKSUMS.sha256"]
  writeLines(vapply(files, function(f) paste(sha(f), sub(paste0("^", report_dir, "/"), "", f)), character(1)),
             file.path(report_dir, "CHECKSUMS.sha256"))
  invisible(list(reason = reason, f_pass = f_pass, h_pass = h_pass, best = best_f))
}

if (!identical(Sys.getenv("NESTA_DECISION_REPAIR_SOURCE_ONLY"), "1")) run_decision_rule_repair()
