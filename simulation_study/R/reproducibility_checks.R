read_tsv_checked <- function(path) {
  if (!file.exists(path)) stop("Missing required file: ", path)
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

final_metric_expectations <- function() {
  data.frame(
    scenario = c(
      "decision_rule_repair", "decision_rule_repair", "decision_rule_repair",
      "comparator_framed_confirmatory", "comparator_framed_confirmatory",
      "comparator_framed_confirmatory", "comparator_framed_confirmatory",
      "comparator_framed_confirmatory", "comparator_framed_confirmatory",
      "comparator_framed_confirmatory"
    ),
    metric = c(
      "F_top100", "F_top150", "F_top200",
      "F_top100", "F_top150", "F_top200",
      "risk_direction_recovery", "protective_direction_recovery",
      "opposite_sign_decoy_selection", "high_score_decoy_selection"
    ),
    expected = c(0.79125, 0.9975, 1.0, 0.7795, 0.99875, 1.0,
                 0.79625, 0.76275, 0.0, 0.0),
    stringsAsFactors = FALSE
  )
}

threshold_compact_expectations <- function() {
  data.frame(
    threshold_label = c("top_1pct", "top_1pct", "top_5pct", "top_5pct",
                        "top_5pct", "top_5pct"),
    selection_rule = c("Final Heat only", "Final OR delta_NESTA",
                       "Final Heat only", "delta_NESTA only",
                       "Final OR delta_NESTA", "Final AND delta_NESTA"),
    selected_gene_count_mean = c(10.0, 12.4, 50.0, 50.0, 57.6, 42.4),
    true_target_recall_mean = c(0.00019, 0.00050, 0.1304, 0.2191, 0.2248, 0.1246),
    FPR_mean = c(0.0105, 0.0130, 0.0471, 0.0434, 0.0512, 0.0394),
    risk_target_recall_mean = c(0.0044, 0.0068, 0.8074, 0.9251, 0.9270, 0.8055),
    protective_target_recall_mean = c(0.0041, 0.0074, 0.7821, 0.9135, 0.9159, 0.7798),
    opposite_sign_decoy_selection_rate_mean = 0,
    high_score_decoy_selection_rate_mean = 0,
    stringsAsFactors = FALSE
  )
}

compare_final_metrics <- function(out_dir, tolerance = 1e-6) {
  observed <- read_tsv_checked(file.path(out_dir, "reproduced_final_metric_comparison.tsv"))
  required <- final_metric_expectations()
  key <- paste(required$scenario, required$metric, sep = "||")
  got_key <- paste(observed$scenario, observed$metric, sep = "||")
  if (!all(key %in% got_key)) {
    stop("Missing required reproducibility metric(s): ",
         paste(key[!(key %in% got_key)], collapse = ", "))
  }
  observed <- observed[match(key, got_key), , drop = FALSE]
  observed$expected_required <- required$expected
  observed$required_abs_diff <- abs(as.numeric(observed$observed) - observed$expected_required)
  observed$required_tolerance <- tolerance
  observed$required_passed <- observed$required_abs_diff <= tolerance
  observed
}

compare_compact_threshold_table <- function(out_dir, tolerance = 0.02,
                                            count_tolerance = 0.5) {
  path <- file.path(out_dir, "manuscript_tables", "table_delta_threshold_sensitivity_compact.tsv")
  observed <- read_tsv_checked(path)
  expected <- threshold_compact_expectations()
  key <- paste(expected$threshold_label, expected$selection_rule, sep = "||")
  got_key <- paste(observed$threshold_label, observed$selection_rule, sep = "||")
  if (!all(key %in% got_key)) {
    stop("Missing required compact threshold row(s): ",
         paste(key[!(key %in% got_key)], collapse = ", "))
  }
  obs <- observed[match(key, got_key), , drop = FALSE]
  metrics <- setdiff(names(expected), c("threshold_label", "selection_rule"))
  rows <- list()
  i <- 1L
  for (metric in metrics) {
    tol <- if (metric == "selected_gene_count_mean") count_tolerance else tolerance
    diff <- abs(as.numeric(obs[[metric]]) - as.numeric(expected[[metric]]))
    rows[[i]] <- data.frame(
      threshold_label = expected$threshold_label,
      selection_rule = expected$selection_rule,
      metric = metric,
      observed = as.numeric(obs[[metric]]),
      expected = as.numeric(expected[[metric]]),
      abs_diff = diff,
      tolerance = tol,
      passed = diff <= tol,
      stringsAsFactors = FALSE
    )
    i <- i + 1L
  }
  do.call(rbind, rows)
}

verify_delta_workflow <- function(out_dir, report_dir = out_dir) {
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  final_comparison <- compare_final_metrics(out_dir)
  compact_comparison <- compare_compact_threshold_table(out_dir)
  utils::write.table(final_comparison,
                     file.path(report_dir, "reviewer_final_metric_comparison.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(compact_comparison,
                     file.path(report_dir, "reviewer_compact_threshold_table_comparison.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  list(
    final_comparison = final_comparison,
    compact_comparison = compact_comparison,
    passed = all(final_comparison$required_passed) && all(compact_comparison$passed)
  )
}

code_fidelity_audit <- function(root = "/home/js/NESTA/simulation_study") {
  export_script <- file.path(root, "scripts", "08_per_gene_export_delta_sensitivity.R")
  provenance_dir <- file.path(root, "internal_provenance", "finalized_simulation_code")
  src <- if (file.exists(export_script)) readLines(export_script, warn = FALSE) else character()
  data.frame(
    check = c(
      "finalized_provenance_code_present",
      "delta_nesta_exported",
      "public_per_gene_output_uses_clean_initialization_name",
      "old_internal_output_column_removed",
      "dropbox_sync_optional"
    ),
    passed = c(
      dir.exists(file.path(provenance_dir, "R")) &&
        file.exists(file.path(provenance_dir, "config", "FROZEN_CONFIG.yaml")),
      any(grepl("delta_NESTA", src, fixed = TRUE)) &&
        any(grepl("Final_Heat", src, fixed = TRUE)) &&
        any(grepl("TWAS.Z", src, fixed = TRUE)),
      any(grepl("expression_weighted_initialization", src, fixed = TRUE)),
      !any(grepl("M2_no_diffusion_final_heat", src, fixed = TRUE)),
      any(grepl("--sync-dropbox", src, fixed = TRUE)) &&
        any(grepl("NESTA_SYNC_DROPBOX", src, fixed = TRUE))
    ),
    detail = c(
      "Vendored finalized code/config are isolated under internal_provenance for deterministic reproduction.",
      "Public per-gene tables expose Final_Heat, TWAS.Z, and delta_NESTA.",
      "The public initialization column is expression_weighted_initialization.",
      "The public export script does not write the old internal no-diffusion output column.",
      "Dropbox copying is not part of the default public workflow."
    ),
    stringsAsFactors = FALSE
  )
}
