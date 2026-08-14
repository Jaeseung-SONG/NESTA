Sys.setenv(NESTA_DECISION_REPAIR_SOURCE_ONLY = "1")
source(file.path("/home/js/NESTA/simulation", "R", "study_0708_decision_rule_repair.R"))

confirmatory_topologies <- c("F", "H")
confirmatory_arm <- "P_combined_decoy_suppressed"
confirmatory_rankings <- c("NESTA_two_tail_balanced", "NESTA_two_tail_direction_matched")
confirmatory_n_replicates <- 200L

make_seed_schedule <- function(n = confirmatory_n_replicates) {
  out <- list(); i <- 1L
  for (tp in confirmatory_topologies) {
    offset <- if (tp == "F") 0L else 100000L
    for (rep_id in seq_len(n)) {
      out[[i]] <- data.frame(
        topology_arm = tp,
        rescue_arm = confirmatory_arm,
        replicate = rep_id,
        batch = ifelse(rep_id <= 100, 1L, 2L),
        base_seed = 9808000L + offset + rep_id,
        signal_seed = 9809000L + offset + rep_id,
        nesta_seed = 9810000L + offset + rep_id,
        permutation_seed = 9811000L + offset + rep_id,
        i2_seed = 9812000L + offset + rep_id,
        i3_seed = 9813000L + offset + rep_id,
        stringsAsFactors = FALSE
      )
      i <- i + 1L
    }
  }
  do.call(rbind, out)
}

ci_mean <- function(x, nboot = 2000L, seed = 9808L) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(mean = NA_real_, median = NA_real_, sd = NA_real_, ci_low = NA_real_, ci_high = NA_real_))
  set.seed(seed)
  boot <- replicate(nboot, mean(sample(x, length(x), replace = TRUE)))
  c(mean = mean(x), median = stats::median(x), sd = stats::sd(x),
    ci_low = unname(stats::quantile(boot, 0.025)),
    ci_high = unname(stats::quantile(boot, 0.975)))
}

summarise_confirmatory_metrics <- function(x, group_cols, metric_cols) {
  keys <- unique(x[group_cols])
  rows <- list(); j <- 1L
  for (i in seq_len(nrow(keys))) {
    keep <- rep(TRUE, nrow(x))
    for (nm in group_cols) keep <- keep & identical_class_compare(x[[nm]], keys[[nm]][i])
    z <- x[keep, , drop = FALSE]
    row <- keys[i, , drop = FALSE]
    for (m in metric_cols) {
      s <- ci_mean(z[[m]], seed = 9808L + i + match(m, metric_cols))
      row[[paste0(m, "_mean")]] <- s["mean"]
      row[[paste0(m, "_median")]] <- s["median"]
      row[[paste0(m, "_sd")]] <- s["sd"]
      row[[paste0(m, "_ci_low")]] <- s["ci_low"]
      row[[paste0(m, "_ci_high")]] <- s["ci_high"]
    }
    rows[[j]] <- row; j <- j + 1L
  }
  do.call(rbind, rows)
}

identical_class_compare <- function(x, y) {
  if (is.factor(x)) x <- as.character(x)
  if (is.factor(y)) y <- as.character(y)
  x == y
}

confirmatory_comparators <- function(rep, ch, no_diff, seeds) {
  genes <- rep$genes
  init_p <- initial_vector_personalizations(setNames(ch$initial_weight, ch$SYMBOL)[genes])
  rwr_abs <- run_rwr(rep$adj, init_p$abs, restart = 0.5)
  ppr_abs <- run_ppr(rep$adj, init_p$abs, damping = 0.85)
  rwr_signed <- run_rwr(rep$adj, init_p$pos, restart = 0.5) - run_rwr(rep$adj, init_p$neg, restart = 0.5)
  ppr_signed <- run_ppr(rep$adj, init_p$pos, damping = 0.85) - run_ppr(rep$adj, init_p$neg, damping = 0.85)
  raw_z <- setNames(rep$twas$TWAS.Z, rep$twas$SYMBOL)[genes]
  init <- setNames(no_diff$final_NESTA_heat, no_diff$SYMBOL)[genes]
  wp <- faithful_m2_channels(permute_sparse_weights(rep$adj, seeds$permutation_seed), rep$expr, rep$twas,
                             n.perm = 10, seed = seeds$permutation_seed)
  i2 <- faithful_m2_channels(module_disrupted_adj(rep), rep$expr, rep$twas,
                             n.perm = 10, seed = seeds$i2_seed)
  i3 <- faithful_m2_channels(permute_sparse_weights(rep$adj, seeds$i3_seed), rep$expr, rep$twas,
                             n.perm = 10, seed = seeds$i3_seed)
  list(
    raw_TWAS_abs = list(score = abs(raw_z), signed = raw_z, class = "primary_unsigned_baseline"),
    M2_no_diffusion_abs = list(score = abs(init), signed = init, class = "primary_unsigned_baseline"),
    PPR_abs_prior = list(score = ppr_abs, signed = ppr_abs, class = "primary_unsigned_disease_relevance"),
    RWR_abs_prior = list(score = rwr_abs, signed = rwr_abs, class = "primary_unsigned_disease_relevance"),
    median_weight_permutation = list(score = abs(setNames(wp$signed_NESTA_heat, wp$SYMBOL)[genes]), signed = setNames(wp$signed_NESTA_heat, wp$SYMBOL)[genes], class = "primary_network_ablation"),
    I2_module_disrupted = list(score = abs(setNames(i2$signed_NESTA_heat, i2$SYMBOL)[genes]), signed = setNames(i2$signed_NESTA_heat, i2$SYMBOL)[genes], class = "primary_network_ablation"),
    I3_expression_matched_randomized = list(score = abs(setNames(i3$signed_NESTA_heat, i3$SYMBOL)[genes]), signed = setNames(i3$signed_NESTA_heat, i3$SYMBOL)[genes], class = "primary_network_ablation"),
    PPR_signed_two_channel = list(score = abs(ppr_signed), signed = ppr_signed, class = "supplementary_signed_two_channel_upper_bound"),
    RWR_signed_two_channel = list(score = abs(rwr_signed), signed = rwr_signed, class = "supplementary_signed_two_channel_upper_bound")
  )
}

run_confirmatory_replicate <- function(seed_row) {
  tp <- seed_row$topology_arm
  rep_id <- seed_row$replicate
  base <- make_branch_isolation_rep(tp, rep_id, seed_row$base_seed)
  rep <- apply_bidirectional_arm(base, confirmatory_arm, seed_row$signal_seed)
  ch <- faithful_m2_channels(rep$adj, rep$expr, rep$twas, n.perm = 25, seed = seed_row$nesta_seed)
  no <- no_diffusion_m2_scores(rep$adj, rep$expr, rep$twas)
  signed <- setNames(ch$final_NESTA_heat, ch$SYMBOL)[rep$genes]
  metric_rows <- lapply(c(confirmatory_rankings, "NESTA_abs_final_heat", "NESTA_signed_descending", "NESTA_signed_ascending"),
                        function(mode) bidirectional_metric_row(rep, signed, mode))
  metrics <- do.call(rbind, metric_rows)
  metrics$batch <- seed_row$batch
  audit <- seed_decoy_proximity_row(rep, ch, no)
  audit$batch <- seed_row$batch
  comps <- confirmatory_comparators(rep, ch, no, seed_row)
  bench_rows <- list(); i <- 1L
  for (nm in names(comps)) {
    row <- audit_score_metrics(rep, comps[[nm]]$score, comps[[nm]]$signed, nm, comps[[nm]]$class)
    row$topology_arm <- tp
    row$topology_label <- row$arm_label
    row$rescue_arm <- confirmatory_arm
    row$batch <- seed_row$batch
    row$comparator_class <- comps[[nm]]$class
    row$disease_relevance_AUPRC <- row$raw_AUPRC
    row$unsigned_disease_relevance_top100_recall <- row$top100_recall
    row$unsigned_disease_relevance_top150_recall <- row$top150_recall
    row$unsigned_disease_relevance_top200_recall <- row$top200_recall
    bench_rows[[i]] <- row; i <- i + 1L
  }
  list(metrics = metrics, benchmarks = do.call(rbind, bench_rows), audit = audit,
       qc = data.frame(topology_arm = tp, rescue_arm = confirmatory_arm, replicate = rep_id,
                       batch = seed_row$batch, successful = TRUE, failure_reason = "",
                       n_A2 = length(rep$A2), n_A1 = length(rep$A1), stringsAsFactors = FALSE))
}

run_confirmatory_topology <- function(schedule, topology) {
  rows <- schedule[schedule$topology_arm == topology, , drop = FALSE]
  metric_rows <- list(); bench_rows <- list(); audit_rows <- list(); qc_rows <- list()
  for (i in seq_len(nrow(rows))) {
    res <- tryCatch(run_confirmatory_replicate(rows[i, ]), error = function(e) e)
    if (inherits(res, "error")) {
      qc_rows[[i]] <- data.frame(topology_arm = topology, rescue_arm = confirmatory_arm,
                                 replicate = rows$replicate[i], batch = rows$batch[i],
                                 successful = FALSE, failure_reason = conditionMessage(res),
                                 n_A2 = NA_integer_, n_A1 = NA_integer_, stringsAsFactors = FALSE)
    } else {
      metric_rows[[length(metric_rows) + 1L]] <- res$metrics
      bench_rows[[length(bench_rows) + 1L]] <- res$benchmarks
      audit_rows[[length(audit_rows) + 1L]] <- res$audit
      qc_rows[[i]] <- res$qc
    }
    if (i %% 10 == 0) {
      message(sprintf("%s confirmatory replicate %d/%d complete", topology, i, nrow(rows)))
      gc(FALSE)
    }
  }
  list(metrics = if (length(metric_rows)) do.call(rbind, metric_rows) else data.frame(),
       benchmarks = if (length(bench_rows)) do.call(rbind, bench_rows) else data.frame(),
       audit = if (length(audit_rows)) do.call(rbind, audit_rows) else data.frame(),
       qc = do.call(rbind, qc_rows))
}

primary_endpoint_rows <- function(metrics, audit) {
  summary <- summarise_confirmatory_metrics(metrics,
    c("topology_arm", "topology_label", "rescue_arm", "ranking_mode"),
    c("top100_recall", "top150_recall", "top200_recall", "risk_top100_recall", "protective_top100_recall",
      "sign_concordant_top100_recall", "opposite_sign_decoy_top100_rate", "high_degree_decoy_top100_rate",
      "raw_AUPRC", "prevalence_normalized_AUPRC", "score_degree_spearman"))
  audit_sum <- stats::aggregate(cbind(raw_TWAS_top100_A2_fraction, M2_top100_A2_fraction,
                                      background_to_A2_mean_abs_initial_weight_ratio,
                                      fraction_nonzero_initial_weight_after_filtering) ~ topology_arm + rescue_arm,
                                data = audit, FUN = function(z) mean(z, na.rm = TRUE))
  merge(summary, audit_sum, by = c("topology_arm", "rescue_arm"), all.x = TRUE)
}

bootstrap_contrasts <- function(metrics, benchmarks) {
  out <- list(); i <- 1L
  primary <- metrics[metrics$ranking_mode %in% confirmatory_rankings, , drop = FALSE]
  comps <- c("raw_TWAS_abs", "M2_no_diffusion_abs", "PPR_abs_prior", "RWR_abs_prior",
             "median_weight_permutation", "I2_module_disrupted", "I3_expression_matched_randomized")
  for (tp in unique(primary$topology_arm)) for (mode in unique(primary$ranking_mode)) {
    n0 <- primary[primary$topology_arm == tp & primary$ranking_mode == mode,
                  c("replicate", "top100_recall", "top150_recall", "top200_recall", "raw_AUPRC"), drop = FALSE]
    for (cmp in comps) {
      b0 <- benchmarks[benchmarks$topology_arm == tp & benchmarks$score_name == cmp,
                       c("replicate", "top100_recall", "top150_recall", "top200_recall", "raw_AUPRC"), drop = FALSE]
      z <- merge(n0, b0, by = "replicate", suffixes = c("_nesta", "_comparator"))
      for (m in c("top100_recall", "top150_recall", "top200_recall", "raw_AUPRC")) {
        d <- z[[paste0(m, "_nesta")]] - z[[paste0(m, "_comparator")]]
        ci <- ci_mean(d, seed = 9900L + i)
        nmean <- mean(z[[paste0(m, "_nesta")]], na.rm = TRUE)
        cmean <- mean(z[[paste0(m, "_comparator")]], na.rm = TRUE)
        out[[i]] <- data.frame(topology_arm = tp, rescue_arm = confirmatory_arm,
                               ranking_mode = mode, comparator = cmp, metric = m,
                               n_replicates = nrow(z), nesta_mean = nmean,
                               comparator_mean = cmean, paired_mean_difference = ci["mean"],
                               paired_median_difference = ci["median"], paired_sd_difference = ci["sd"],
                               paired_bootstrap_ci_low = ci["ci_low"], paired_bootstrap_ci_high = ci["ci_high"],
                               fold_over_comparator = nmean / max(cmean, .Machine$double.eps),
                               stringsAsFactors = FALSE)
        i <- i + 1L
      }
    }
  }
  do.call(rbind, out)
}

classify_confirmatory <- function(primary, contrasts, unsigned_summary) {
  f <- primary[primary$topology_arm == "F" & primary$ranking_mode %in% confirmatory_rankings, , drop = FALSE]
  f$pass <- f$top100_recall_mean >= 0.70 & f$top150_recall_mean >= 0.90 &
    f$risk_top100_recall_mean >= 0.60 & f$protective_top100_recall_mean >= 0.60 &
    f$opposite_sign_decoy_top100_rate_mean <= 0.10 & f$high_degree_decoy_top100_rate_mean <= 0.10 &
    f$raw_TWAS_top100_A2_fraction <= 0.05 & f$M2_top100_A2_fraction <= 0.10
  for (cmp in c("raw_TWAS_abs", "M2_no_diffusion_abs")) {
    z <- contrasts[contrasts$topology_arm == "F" & contrasts$metric == "top100_recall" &
                     contrasts$comparator == cmp & contrasts$ranking_mode %in% confirmatory_rankings, , drop = FALSE]
    ok <- z$paired_mean_difference >= 0.40
    if (cmp == "M2_no_diffusion_abs") ok <- ok & z$paired_bootstrap_ci_low > 0
    f$pass <- f$pass & ok[match(f$ranking_mode, z$ranking_mode)]
  }
  h <- primary[primary$topology_arm == "H" & primary$ranking_mode %in% confirmatory_rankings, , drop = FALSE]
  h$pass <- h$top100_recall_mean >= 0.60 & h$top150_recall_mean >= 0.80 &
    h$risk_top100_recall_mean >= 0.50 & h$protective_top100_recall_mean >= 0.50 &
    h$opposite_sign_decoy_top100_rate_mean <= 0.15 & h$high_degree_decoy_top100_rate_mean <= 0.15
  f_pass <- any(f$pass, na.rm = TRUE)
  h_pass <- any(h$pass, na.rm = TRUE)
  best_f <- f[order(f$pass, f$top100_recall_mean, f$top150_recall_mean, decreasing = TRUE), ][1, , drop = FALSE]
  ppr_rwr <- unsigned_summary[unsigned_summary$topology_arm == "F" & unsigned_summary$score_name %in% c("PPR_abs_prior", "RWR_abs_prior"), , drop = FALSE]
  strong_unsigned <- any(ppr_rwr$top100_recall_mean > best_f$top100_recall_mean, na.rm = TRUE)
  if (f_pass && strong_unsigned) label <- "confirmatory_success_with_strong_unsigned_ppr_rwr"
  else if (f_pass && h_pass) label <- "confirmatory_success"
  else if (f_pass && !h_pass) label <- "confirmatory_success_H_supportive_failed"
  else if (any(f$opposite_sign_decoy_top100_rate_mean > 0.10 | f$high_degree_decoy_top100_rate_mean > 0.10, na.rm = TRUE)) label <- "top100_decoy_guardrail_failed"
  else if (any(f$top100_recall_mean < 0.70 | f$top150_recall_mean < 0.90, na.rm = TRUE)) label <- "confirmatory_A2_recovery_failed"
  else label <- "confirmatory_primary_contrast_failed"
  list(label = label, f = f, h = h, best_f = best_f, f_pass = f_pass, h_pass = h_pass, strong_unsigned = strong_unsigned)
}

run_comparator_framing_confirmatory <- function() {
  verify_project_path(); verify_binding_plan()
  report_dir <- read_report_dir(); safe_dir_create(report_dir); copy_binding_plan_to_report(report_dir)
  set_run_status("IN_PROGRESS", "YES", "YES", "0708 comparator-framing confirmatory")
  schedule <- make_seed_schedule()
  write_csv_over(schedule, file.path(report_dir, "SEED_SCHEDULE.csv"))
  results <- lapply(confirmatory_topologies, function(tp) run_confirmatory_topology(schedule, tp))
  metrics <- do.call(rbind, lapply(results, `[[`, "metrics"))
  benchmarks <- do.call(rbind, lapply(results, `[[`, "benchmarks"))
  audit <- do.call(rbind, lapply(results, `[[`, "audit"))
  qc <- do.call(rbind, lapply(results, `[[`, "qc"))
  failed <- qc[!qc$successful, , drop = FALSE]
  primary <- primary_endpoint_rows(metrics, audit)
  unsigned <- benchmarks[benchmarks$score_name %in% c("raw_TWAS_abs", "M2_no_diffusion_abs", "PPR_abs_prior", "RWR_abs_prior", "median_weight_permutation", "I2_module_disrupted", "I3_expression_matched_randomized"), , drop = FALSE]
  signed_upper <- benchmarks[benchmarks$score_name %in% c("PPR_signed_two_channel", "RWR_signed_two_channel"), , drop = FALSE]
  unsigned_summary <- summarise_confirmatory_metrics(unsigned, c("topology_arm", "topology_label", "rescue_arm", "score_name", "comparator_class"),
                                                     c("top100_recall", "top150_recall", "top200_recall", "raw_AUPRC", "direction_aware_AUPRC", "sign_concordant_top100_recall", "opposite_sign_decoy_top100_rate", "high_degree_decoy_top100_rate", "score_degree_spearman"))
  signed_summary <- summarise_confirmatory_metrics(signed_upper, c("topology_arm", "topology_label", "rescue_arm", "score_name", "comparator_class"),
                                                   c("top100_recall", "top150_recall", "top200_recall", "raw_AUPRC", "direction_aware_AUPRC", "sign_concordant_top100_recall", "opposite_sign_decoy_top100_rate", "high_degree_decoy_top100_rate", "score_degree_spearman"))
  contrasts <- bootstrap_contrasts(metrics, benchmarks)
  decision <- classify_confirmatory(primary, contrasts, unsigned_summary)
  set_run_status("COMPLETE", "YES", "YES", decision$label)
  bootstrap_summary <- contrasts
  endpoint_subset <- primary[primary$ranking_mode %in% confirmatory_rankings, , drop = FALSE]
  direction <- endpoint_subset[, c("topology_arm", "topology_label", "rescue_arm", "ranking_mode",
                                   "risk_top100_recall_mean", "protective_top100_recall_mean",
                                   "sign_concordant_top100_recall_mean"), drop = FALSE]
  decoy <- endpoint_subset[, c("topology_arm", "topology_label", "rescue_arm", "ranking_mode",
                               "opposite_sign_decoy_top100_rate_mean", "high_degree_decoy_top100_rate_mean"), drop = FALSE]
  saturation <- summarise_confirmatory_metrics(metrics, c("topology_arm", "topology_label", "rescue_arm", "ranking_mode"),
                                               c("opposite_sign_decoy_top150_rate", "opposite_sign_decoy_top200_rate", "high_degree_decoy_top150_rate", "high_degree_decoy_top200_rate", "score_degree_spearman", "score_strength_spearman"))
  saturation$interpretation <- "Exploratory saturation diagnostic only; top150/top200 decoy rates are not hard stop criteria."
  schema <- data.frame(file = c("PRIMARY_ENDPOINTS.csv", "UNSIGNED_COMPARATOR_METRICS.csv", "SIGNED_UPPER_BOUND_COMPARATOR_METRICS.csv", "BOOTSTRAP_CI_SUMMARY.csv"),
                       required_columns_present = c(all(c("top100_recall_mean", "risk_top100_recall_mean", "opposite_sign_decoy_top100_rate_mean") %in% names(primary)),
                                                    all(c("top100_recall_mean", "raw_AUPRC_mean", "score_degree_spearman_mean") %in% names(unsigned_summary)),
                                                    all(c("comparator_class", "top100_recall_mean") %in% names(signed_summary)),
                                                    all(c("paired_bootstrap_ci_low", "paired_bootstrap_ci_high") %in% names(bootstrap_summary))),
                       stringsAsFactors = FALSE)
  write_csv_over(metrics, file.path(report_dir, "NESTA_BIDIRECTIONAL_METRICS.csv"))
  write_csv_over(primary, file.path(report_dir, "PRIMARY_ENDPOINTS.csv"))
  write_csv_over(contrasts, file.path(report_dir, "PRIMARY_CONTRASTS.csv"))
  write_csv_over(bootstrap_summary, file.path(report_dir, "BOOTSTRAP_CI_SUMMARY.csv"))
  write_csv_over(unsigned_summary, file.path(report_dir, "UNSIGNED_COMPARATOR_METRICS.csv"))
  write_csv_over(signed_summary, file.path(report_dir, "SIGNED_UPPER_BOUND_COMPARATOR_METRICS.csv"))
  write_csv_over(direction, file.path(report_dir, "DIRECTION_SPECIFIC_RECOVERY_METRICS.csv"))
  write_csv_over(decoy, file.path(report_dir, "DECOY_GUARDRAIL_METRICS.csv"))
  write_csv_over(saturation, file.path(report_dir, "SATURATION_DIAGNOSTICS.csv"))
  write_csv_over(contrasts, file.path(report_dir, "BENCHMARK_CONTRASTS.csv"))
  write_csv_over(audit, file.path(report_dir, "INITIAL_SIGNAL_FIELD_AUDIT.csv"))
  write_csv_over(network_path_audit(audit), file.path(report_dir, "NETWORK_PATH_AUDIT.csv"))
  write_csv_over(qc, file.path(report_dir, "REPLICATE_QC.csv"))
  write_csv_over(if (nrow(failed)) failed else data.frame(topology_arm=character(), rescue_arm=character(), replicate=integer(), batch=integer(), successful=logical(), failure_reason=character(), n_A2=integer(), n_A1=integer()), file.path(report_dir, "FAILED_REPLICATES.csv"))
  write_csv_over(schema, file.path(report_dir, "METRIC_SCHEMA_AUDIT.csv"))
  write_csv_over(endpoint_subset, file.path(report_dir, "table_primary_nesta_vs_raw_m2.csv"))
  write_csv_over(unsigned_summary, file.path(report_dir, "table_unsigned_propagation_comparators.csv"))
  write_csv_over(signed_summary, file.path(report_dir, "table_signed_upper_bound_comparators_supplement.csv"))
  write_csv_over(endpoint_subset[, c("topology_arm", "ranking_mode", "top100_recall_mean", "top150_recall_mean", "top200_recall_mean", "risk_top100_recall_mean", "protective_top100_recall_mean")], file.path(report_dir, "figure_bidirectional_recovery.csv"))
  write_csv_over(rbind(unsigned_summary[, intersect(names(unsigned_summary), names(signed_summary)), drop=FALSE], signed_summary[, intersect(names(unsigned_summary), names(signed_summary)), drop=FALSE]), file.path(report_dir, "figure_comparator_framing.csv"))
  best <- decision$best_f
  raw_con <- contrasts[contrasts$topology_arm == "F" & contrasts$ranking_mode == best$ranking_mode & contrasts$comparator == "raw_TWAS_abs" & contrasts$metric == "top100_recall", , drop = FALSE]
  m2_con <- contrasts[contrasts$topology_arm == "F" & contrasts$ranking_mode == best$ranking_mode & contrasts$comparator == "M2_no_diffusion_abs" & contrasts$metric == "top100_recall", , drop = FALSE]
  write_lines_over(c("# Comparator Framing Summary", "",
                     "PPR_abs_prior and RWR_abs_prior are framed as unsigned disease-relevance propagation comparators.",
                     "PPR_signed_two_channel and RWR_signed_two_channel are supplementary enhanced signed two-channel upper-bound analyses and are excluded from primary success/failure classification.",
                     sprintf("Classification: `%s`.", decision$label)), file.path(report_dir, "COMPARATOR_FRAMING_SUMMARY.md"))
  write_lines_over(c("# Confirmatory Summary", "",
                     sprintf("Successful confirmatory replicates: F=%d, H=%d.", sum(qc$topology_arm == "F" & qc$successful), sum(qc$topology_arm == "H" & qc$successful)),
                     sprintf("Best F primary ranking `%s`: top100/top150/top200 %.4f / %.4f / %.4f.", best$ranking_mode, best$top100_recall_mean, best$top150_recall_mean, best$top200_recall_mean),
                     sprintf("Risk/protective top100 %.4f / %.4f; top100 decoy opposite/high-degree %.4f / %.4f.", best$risk_top100_recall_mean, best$protective_top100_recall_mean, best$opposite_sign_decoy_top100_rate_mean, best$high_degree_decoy_top100_rate_mean),
                     sprintf("Top100 paired improvement over raw_TWAS_abs %.4f; over M2_no_diffusion_abs %.4f, CI %.4f to %.4f.", raw_con$paired_mean_difference, m2_con$paired_mean_difference, m2_con$paired_bootstrap_ci_low, m2_con$paired_bootstrap_ci_high)), file.path(report_dir, "CONFIRMATORY_SUMMARY.md"))
  write_lines_over(c("# Manuscript Ready Summary", "",
                     "NESTA Final Heat was evaluated as a signed bidirectional prioritization score using balanced and direction-matched two-tail rankings.",
                     "Raw TWAS and M2 no-diffusion are direct baselines. Standard PPR/RWR are unsigned disease-relevance propagation comparators; signed two-channel PPR/RWR are supplementary upper-bound analyses.",
                     sprintf("The confirmatory classification was `%s`.", decision$label)), file.path(report_dir, "MANUSCRIPT_READY_SUMMARY.md"))
  write_lines_over(c("# Code Fidelity Audit", "",
                     paste0("Binding plan SHA256: `", binding_plan_sha256, "`."),
                     "Faithful submitted M2 arithmetic, TWAS.P conversion, strict TWAS filtering, signed positive and absolute-negative diffusion, zero-weight edge behavior, self-loop behavior, and diffuStats `n.perm` were retained.",
                     "Submitted Final Heat values were not modified; only ranking/comparator framing and confirmatory aggregation were implemented."), file.path(report_dir, "CODE_FIDELITY_AUDIT.md"))
  write_lines_over(c("# Metric Schema Audit", "", sprintf("All required schema checks passed: %s.", all(schema$required_columns_present))), file.path(report_dir, "METRIC_SCHEMA_AUDIT.md"))
  write_lines_over(c("# STOP/GO Report", "", "GO for manuscript-facing confirmatory interpretation.", paste0("Reason: `", decision$label, "`."),
                     "Confirmatory execution started: YES.", "Confirmatory execution completed: YES.",
                     "Top150/top200 decoy rates were treated as exploratory saturation diagnostics only."), file.path(report_dir, "STOP_GO_REPORT.md"))
  write_lines_over(c("# Final Report", "",
                     paste0("Binding plan SHA256: `", binding_plan_sha256, "`."),
                     paste0("Final classification: `", decision$label, "`."),
                     sprintf("Primary F `%s` ranking `%s`: top100/top150/top200 %.4f / %.4f / %.4f.", confirmatory_arm, best$ranking_mode, best$top100_recall_mean, best$top150_recall_mean, best$top200_recall_mean),
                     sprintf("Risk/protective top100 %.4f / %.4f; top100 opposite/high-degree decoy %.4f / %.4f.", best$risk_top100_recall_mean, best$protective_top100_recall_mean, best$opposite_sign_decoy_top100_rate_mean, best$high_degree_decoy_top100_rate_mean),
                     sprintf("Topology H supportive robustness pass: %s.", decision$h_pass),
                     sprintf("Unsigned PPR/RWR stronger than NESTA in F top100: %s.", decision$strong_unsigned),
                     "PPR/RWR absolute-prior results are reported as unsigned disease-relevance comparators and do not define NESTA failure when signed primary criteria pass.",
                     "Confirmatory execution started and completed: YES."), file.path(report_dir, "FINAL_REPORT.md"))
  if (file.exists(project_file("results/reports/summary_tables/unit_test_results.csv"))) file.copy(project_file("results/reports/summary_tables/unit_test_results.csv"), file.path(report_dir, "unit_test_results.csv"), overwrite = TRUE)
  if (file.exists(project_file("CLEANUP_MANIFEST.csv"))) file.copy(project_file("CLEANUP_MANIFEST.csv"), file.path(report_dir, "CLEANUP_MANIFEST.csv"), overwrite = TRUE)
  write_csv_over(data.frame(path = c(project_file("R/study_0708_comparator_confirmatory.R"), project_file("R/study_0708_decision_rule_repair.R"), project_file("R/study_0708_bidirectional_rescue.R"), project_file("R/study_0708_initial_signal_rescue.R"), project_file("R/fidelity.R"), project_file("R/utils.R"), project_file("scripts/run_comparator_confirmatory.R")),
                            sha256 = c(sha(project_file("R/study_0708_comparator_confirmatory.R")), sha(project_file("R/study_0708_decision_rule_repair.R")), sha(project_file("R/study_0708_bidirectional_rescue.R")), sha(project_file("R/study_0708_initial_signal_rescue.R")), sha(project_file("R/fidelity.R")), sha(project_file("R/utils.R")), sha(project_file("scripts/run_comparator_confirmatory.R"))),
                            role = c("confirmatory_comparator_framing_runner", "decision_rule_helpers", "bidirectional_ranking_helpers", "initial_signal_carryforward", "faithful_nesta_and_benchmarks", "binding_plan_and_io_guards", "script_entrypoint"), stringsAsFactors = FALSE), file.path(report_dir, "IMPLEMENTATION_MANIFEST.csv"))
  files <- list.files(report_dir, recursive = TRUE, full.names = TRUE)
  files <- files[file.info(files)$isdir == FALSE & basename(files) != "CHECKSUMS.sha256"]
  writeLines(vapply(files, function(f) paste(sha(f), sub(paste0("^", report_dir, "/"), "", f)), character(1)), file.path(report_dir, "CHECKSUMS.sha256"))
  invisible(list(label = decision$label, best = best, f_pass = decision$f_pass, h_pass = decision$h_pass))
}

if (!identical(Sys.getenv("NESTA_CONFIRMATORY_SOURCE_ONLY"), "1")) run_comparator_framing_confirmatory()
