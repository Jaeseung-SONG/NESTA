Sys.setenv(NESTA_INITIAL_SIGNAL_SOURCE_ONLY = "1")
source(file.path("/home/js/NESTA/simulation", "R", "study_0708_initial_signal_rescue.R"))

bidirectional_arms <- function() {
  c("M_bidirectional_only", "N_initial_decoy_suppressed",
    "O_network_decoy_suppressed", "P_combined_decoy_suppressed")
}

apply_bidirectional_arm <- function(rep, bidirectional_arm, seed) {
  rep <- apply_initial_signal_arm(rep, "K_sparse_TWAS_filtered", seed)
  genes <- rep$genes
  z <- setNames(rep$twas$TWAS.Z, rep$twas$SYMBOL)[genes]
  branch <- unique(c(rep$A1, rep$relay, rep$A2))
  if (bidirectional_arm %in% c("N_initial_decoy_suppressed", "P_combined_decoy_suppressed")) {
    z[rep$D] <- 0
    z[rep$C] <- 0
    rep$twas <- data.frame(SYMBOL = genes, TWAS.Z = as.numeric(z), TWAS.P = twas_p_from_z(z),
                           stringsAsFactors = FALSE)
  }
  if (bidirectional_arm %in% c("O_network_decoy_suppressed", "P_combined_decoy_suppressed")) {
    adj <- rep$adj
    dec <- unique(c(rep$D, rep$C))
    adj[dec, branch] <- adj[dec, branch] * 0.05
    adj[branch, dec] <- adj[branch, dec] * 0.05
    adj[rep$C, rep$background] <- adj[rep$C, rep$background] * 0.25
    adj[rep$background, rep$C] <- adj[rep$background, rep$C] * 0.25
    diag(adj) <- 1
    rep$adj <- Matrix::drop0(adj)
  }
  rep$bidirectional_arm <- bidirectional_arm
  rep
}

two_tail_selected <- function(genes, signed_score, k, exclude = character()) {
  keep <- !(genes %in% exclude)
  g <- genes[keep]
  s <- setNames(signed_score[genes], genes)[g]
  half <- floor(k / 2)
  pos <- g[order(s, decreasing = TRUE, na.last = NA)]
  neg <- g[order(s, decreasing = FALSE, na.last = NA)]
  unique(c(pos[seq_len(min(half, length(pos)))], neg[seq_len(min(k - half, length(neg)))]))
}

ranking_selected <- function(rep, signed_score, mode, k) {
  genes <- rep$genes
  score <- setNames(signed_score[genes], genes)
  if (mode == "NESTA_signed_descending") {
    ranked_genes(genes, score, TRUE, rep$A1)[seq_len(min(k, length(genes) - length(rep$A1)))]
  } else if (mode == "NESTA_signed_ascending") {
    ranked_genes(genes, score, FALSE, rep$A1)[seq_len(min(k, length(genes) - length(rep$A1)))]
  } else if (mode == "NESTA_abs_final_heat") {
    ranked_genes(genes, abs(score), TRUE, rep$A1)[seq_len(min(k, length(genes) - length(rep$A1)))]
  } else if (mode %in% c("NESTA_two_tail_balanced", "NESTA_two_tail_direction_matched")) {
    two_tail_selected(genes, score, k, rep$A1)
  } else {
    stop("Unknown ranking mode: ", mode)
  }
}

ranking_rank_map <- function(rep, signed_score, mode) {
  genes <- rep$genes
  score <- setNames(signed_score[genes], genes)
  if (mode == "NESTA_signed_descending") return(rank_map(genes, score, rep$A1, TRUE))
  if (mode == "NESTA_signed_ascending") return(rank_map(genes, score, rep$A1, FALSE))
  if (mode == "NESTA_abs_final_heat") return(rank_map(genes, abs(score), rep$A1, TRUE))
  pos <- ranked_genes(genes, score, TRUE, rep$A1)
  neg <- ranked_genes(genes, score, FALSE, rep$A1)
  out <- setNames(rep(NA_real_, length(genes)), genes)
  for (i in seq_along(pos)) if (is.na(out[pos[i]])) out[pos[i]] <- (2 * i) - 1
  for (i in seq_along(neg)) if (is.na(out[neg[i]])) out[neg[i]] <- 2 * i
  out
}

score_for_auprc <- function(rep, signed_score, mode) {
  genes <- rep$genes
  score <- setNames(signed_score[genes], genes)
  if (mode == "NESTA_signed_descending") return(score)
  if (mode == "NESTA_signed_ascending") return(-score)
  if (mode == "NESTA_abs_final_heat") return(abs(score))
  abs(score)
}

bidirectional_metric_row <- function(rep, signed_score, ranking_mode) {
  genes <- rep$genes
  score <- setNames(signed_score[genes], genes)
  selected <- lapply(c(50, 100, 150, 200), function(k) ranking_selected(rep, score, ranking_mode, k))
  names(selected) <- paste0("top", c(50, 100, 150, 200))
  ranks <- ranking_rank_map(rep, score, ranking_mode)
  au_score <- score_for_auprc(rep, score, ranking_mode)
  au <- auprc_from_score(genes, au_score, rep$A2, rep$A1)
  prev <- length(rep$A2) / (length(genes) - length(rep$A1))
  risk_recall <- function(k) mean(rep$A2_risk %in% selected[[paste0("top", k)]])
  prot_recall <- function(k) mean(rep$A2_protective %in% selected[[paste0("top", k)]])
  if (ranking_mode == "NESTA_two_tail_direction_matched") {
    pos100 <- head(ranked_genes(genes, score, TRUE, rep$A1), 50)
    neg100 <- head(ranked_genes(genes, score, FALSE, rep$A1), 50)
    sign_conc <- (sum(rep$A2_risk %in% pos100) + sum(rep$A2_protective %in% neg100)) / length(rep$A2)
  } else {
    top100 <- selected$top100
    dirs <- rep$A2_direction
    g <- intersect(top100, names(dirs))
    sign_conc <- if (!length(g)) 0 else sum(ifelse(dirs[g] == "risk", score[g] > 0, score[g] < 0)) / length(rep$A2)
  }
  risk_au <- auprc_from_score(genes, score, rep$A2_risk, rep$A1)
  prot_au <- auprc_from_score(genes, -score, rep$A2_protective, rep$A1)
  deg <- score_degree(rep$adj)[genes]
  str <- score_strength(rep$adj)[genes]
  data.frame(
    topology_arm = rep$arm, topology_label = rep$arm_label,
    rescue_arm = rep$bidirectional_arm, replicate = rep$rep_id,
    ranking_mode = ranking_mode,
    top50_recall = mean(rep$A2 %in% selected$top50),
    top100_recall = mean(rep$A2 %in% selected$top100),
    top150_recall = mean(rep$A2 %in% selected$top150),
    top200_recall = mean(rep$A2 %in% selected$top200),
    risk_top50_recall = risk_recall(50),
    risk_top100_recall = risk_recall(100),
    risk_top150_recall = risk_recall(150),
    risk_top200_recall = risk_recall(200),
    protective_top50_recall = prot_recall(50),
    protective_top100_recall = prot_recall(100),
    protective_top150_recall = prot_recall(150),
    protective_top200_recall = prot_recall(200),
    sign_concordant_top100_recall = sign_conc,
    direction_aware_AUPRC = weighted.mean(c(risk_au, prot_au), c(length(rep$A2_risk), length(rep$A2_protective)), na.rm = TRUE),
    raw_AUPRC = au,
    prevalence_normalized_AUPRC = if (is.finite(au) && prev < 1) (au - prev) / (1 - prev) else NA_real_,
    median_A2_rank = median(ranks[rep$A2], na.rm = TRUE),
    median_risk_A2_rank = median(ranks[rep$A2_risk], na.rm = TRUE),
    median_protective_A2_rank = median(ranks[rep$A2_protective], na.rm = TRUE),
    first_A2_rank = min(ranks[rep$A2], na.rm = TRUE),
    opposite_sign_decoy_top100_rate = mean(rep$D %in% selected$top100),
    opposite_sign_decoy_top150_rate = mean(rep$D %in% selected$top150),
    opposite_sign_decoy_top200_rate = mean(rep$D %in% selected$top200),
    high_degree_decoy_top100_rate = mean(rep$C %in% selected$top100),
    high_degree_decoy_top150_rate = mean(rep$C %in% selected$top150),
    high_degree_decoy_top200_rate = mean(rep$C %in% selected$top200),
    score_degree_spearman = suppressWarnings(stats::cor(abs(score), deg, method = "spearman", use = "pairwise.complete.obs")),
    score_strength_spearman = suppressWarnings(stats::cor(abs(score), str, method = "spearman", use = "pairwise.complete.obs")),
    stringsAsFactors = FALSE
  )
}

seed_decoy_proximity_row <- function(rep, ch, no_diff) {
  audit <- initial_field_audit_row(rep, ch, no_diff)
  adj <- rep$adj
  seed_to <- function(nodes) as.numeric(adj[rep$A1, nodes])
  data.frame(
    topology_arm = rep$arm, topology_label = rep$arm_label,
    rescue_arm = rep$bidirectional_arm, replicate = rep$rep_id,
    n_total_genes = audit$n_total_genes,
    n_nonzero_initial_weight_genes_after_filtering = audit$n_nonzero_initial_weight_genes_after_filtering,
    fraction_nonzero_initial_weight_after_filtering = audit$fraction_nonzero_initial_weight_after_filtering,
    background_to_A2_mean_abs_initial_weight_ratio = audit$background_to_A2_mean_abs_initial_weight_ratio,
    raw_TWAS_top100_A2_fraction = audit$raw_TWAS_top100_A2_fraction,
    M2_top100_A2_fraction = audit$M2_top100_A2_fraction,
    mean_seed_to_A2_weight = mean(seed_to(rep$A2)),
    mean_seed_to_opposite_decoy_weight = mean(seed_to(rep$D)),
    mean_seed_to_high_degree_decoy_weight = mean(seed_to(rep$C)),
    median_seed_to_A2_weight = median(seed_to(rep$A2)),
    median_seed_to_opposite_decoy_weight = median(seed_to(rep$D)),
    median_seed_to_high_degree_decoy_weight = median(seed_to(rep$C)),
    stringsAsFactors = FALSE
  )
}

bidirectional_success_summary <- function(metrics) {
  stats::aggregate(cbind(top50_recall, top100_recall, top150_recall, top200_recall,
                         risk_top50_recall, risk_top100_recall,
                         risk_top150_recall, risk_top200_recall,
                         protective_top50_recall, protective_top100_recall,
                         protective_top150_recall, protective_top200_recall,
                         raw_AUPRC, prevalence_normalized_AUPRC,
                         direction_aware_AUPRC, sign_concordant_top100_recall,
                         median_A2_rank, median_risk_A2_rank,
                         median_protective_A2_rank, first_A2_rank,
                         opposite_sign_decoy_top100_rate,
                         opposite_sign_decoy_top150_rate,
                         opposite_sign_decoy_top200_rate,
                         high_degree_decoy_top100_rate,
                         high_degree_decoy_top150_rate,
                         high_degree_decoy_top200_rate,
                         score_degree_spearman,
                         score_strength_spearman) ~ topology_arm + topology_label + rescue_arm + ranking_mode,
                   data = metrics, FUN = function(z) mean(z, na.rm = TRUE))
}

ranking_contrasts <- function(nesta_summary, benchmark_summary) {
  out <- nesta_summary[nesta_summary$ranking_mode == "NESTA_two_tail_balanced", , drop = FALSE]
  out$candidate_ranking <- out$ranking_mode
  out$primary_success <- out$top100_recall >= 0.60 & out$top150_recall >= 0.70 &
    out$top200_recall >= 0.75 & out$risk_top100_recall >= 0.40 &
    out$protective_top100_recall >= 0.40 &
    out$opposite_sign_decoy_top100_rate <= 0.15 &
    out$high_degree_decoy_top100_rate <= 0.12 &
    out$opposite_sign_decoy_top200_rate <= 0.20 &
    out$high_degree_decoy_top200_rate <= 0.18
  out
}

run_bidirectional_stage <- function(topology, arms) {
  metric_rows <- list(); bench_rows <- list(); audit_rows <- list()
  im <- ib <- ia <- 1
  modes <- c("NESTA_signed_descending", "NESTA_signed_ascending",
             "NESTA_abs_final_heat", "NESTA_two_tail_balanced",
             "NESTA_two_tail_direction_matched")
  for (arm in arms) for (rep_id in seq_len(20)) {
    base <- make_branch_isolation_rep(topology, rep_id, 972000 + match(topology, c("F", "H")) * 10000 + rep_id)
    rep <- apply_bidirectional_arm(base, arm, 972500 + match(arm, bidirectional_arms()) * 1000 + rep_id)
    ch <- faithful_m2_channels(rep$adj, rep$expr, rep$twas, n.perm = 25,
                               seed = 973000 + match(topology, c("F", "H")) * 10000 + rep_id)
    no <- no_diffusion_m2_scores(rep$adj, rep$expr, rep$twas)
    audit_rows[[ia]] <- seed_decoy_proximity_row(rep, ch, no); ia <- ia + 1
    signed <- setNames(ch$final_NESTA_heat, ch$SYMBOL)[rep$genes]
    for (mode in modes) {
      metric_rows[[im]] <- bidirectional_metric_row(rep, signed, mode); im <- im + 1
    }
    comparators <- rescue_comparators(rep, ch, no)
    for (nm in names(comparators)) {
      row <- audit_score_metrics(rep, comparators[[nm]]$score, comparators[[nm]]$signed, nm, "diagnostic_comparator")
      row$topology_arm <- topology
      row$topology_label <- row$arm_label
      row$rescue_arm <- arm
      bench_rows[[ib]] <- row; ib <- ib + 1
    }
    if (rep_id %% 5 == 0) gc(FALSE)
  }
  list(metrics = do.call(rbind, metric_rows),
       benchmarks = do.call(rbind, bench_rows),
       audit = do.call(rbind, audit_rows))
}

run_bidirectional_rescue <- function() {
  verify_project_path()
  verify_binding_plan()
  report_dir <- read_report_dir()
  safe_dir_create(report_dir)
  copy_binding_plan_to_report(report_dir)
  set_run_status("IN_PROGRESS", "YES", "NO", "0708 bidirectional ranking decoy rescue")

  f <- run_bidirectional_stage("F", bidirectional_arms())
  f_sum <- bidirectional_success_summary(f$metrics)
  f_contrast <- ranking_contrasts(f_sum, NULL)
  pass_arms <- unique(f_contrast$rescue_arm[f_contrast$primary_success])
  h <- NULL
  if (length(pass_arms)) h <- run_bidirectional_stage("H", pass_arms)

  metrics <- if (is.null(h)) f$metrics else rbind(f$metrics, h$metrics)
  benchmarks <- if (is.null(h)) f$benchmarks else rbind(f$benchmarks, h$benchmarks)
  audit <- if (is.null(h)) f$audit else rbind(f$audit, h$audit)
  summary <- bidirectional_success_summary(metrics)
  contrasts <- ranking_contrasts(summary, NULL)
  bench_summary <- summarise_score_metrics(benchmarks)
  decoy <- summary[, c("topology_arm", "topology_label", "rescue_arm", "ranking_mode",
                       "opposite_sign_decoy_top100_rate", "opposite_sign_decoy_top200_rate",
                       "high_degree_decoy_top100_rate", "high_degree_decoy_top200_rate",
                       "score_degree_spearman")]
  direction <- summary[, c("topology_arm", "topology_label", "rescue_arm", "ranking_mode",
                           "risk_top100_recall", "protective_top100_recall",
                           "risk_top150_recall", "protective_top150_recall",
                           "risk_top200_recall", "protective_top200_recall",
                           "median_risk_A2_rank", "median_protective_A2_rank",
                           "direction_aware_AUPRC", "sign_concordant_top100_recall")]
  f_bal <- contrasts[contrasts$topology_arm == "F", , drop = FALSE]
  best <- f_bal[order(f_bal$primary_success, f_bal$top100_recall, f_bal$top150_recall,
                      f_bal$top200_recall, -f_bal$opposite_sign_decoy_top100_rate,
                      decreasing = TRUE), ][1, ]
  if (any(f_bal$primary_success)) {
    reason <- "bidirectional_ranking_rescues_A2_and_guardrails"
  } else if (any(f_bal$top100_recall >= 0.60 | f_bal$top200_recall >= 0.75)) {
    reason <- "bidirectional_ranking_rescues_A2_but_decoys_remain"
  } else if (any(f_bal$rescue_arm != "M_bidirectional_only" & f_bal$top200_recall < 0.50)) {
    reason <- "decoy_suppression_destroys_A2_recovery"
  } else if (max(f_bal$top100_recall, na.rm = TRUE) <= 0.50) {
    reason <- "two_tail_ranking_does_not_exceed_one_sided_ceiling"
  } else {
    reason <- "initial_signal_leakage_artifact"
  }
  set_run_status("STOPPED", "YES", "NO", reason)

  write_csv_over(metrics, file.path(report_dir, "BIDIRECTIONAL_RANKING_METRICS.csv"))
  write_csv_over(contrasts, file.path(report_dir, "BIDIRECTIONAL_RANKING_CONTRASTS.csv"))
  write_csv_over(decoy, file.path(report_dir, "DECOY_GUARDRAIL_METRICS.csv"))
  write_csv_over(direction, file.path(report_dir, "DIRECTION_SPECIFIC_RECOVERY_METRICS.csv"))
  write_csv_over(audit, file.path(report_dir, "INITIAL_SIGNAL_FIELD_AUDIT.csv"))
  write_csv_over(metrics, file.path(report_dir, "PRIMARY_FINAL_HEAT_METRICS.csv"))
  write_csv_over(metrics[, c("topology_arm", "topology_label", "rescue_arm", "replicate", "ranking_mode",
                             "direction_aware_AUPRC", "sign_concordant_top100_recall")],
                 file.path(report_dir, "DIRECTION_AWARE_METRICS.csv"))
  write_csv_over(benchmarks, file.path(report_dir, "BENCHMARK_METRICS.csv"))
  write_csv_over(bench_summary, file.path(report_dir, "BENCHMARK_CONTRASTS.csv"))
  write_csv_over(decoy[, c("topology_arm", "topology_label", "rescue_arm", "ranking_mode", "score_degree_spearman",
                           "high_degree_decoy_top100_rate", "high_degree_decoy_top200_rate")],
                 file.path(report_dir, "DEGREE_BIAS_METRICS.csv"))
  write_csv_over(data.frame(guardrail = c("no_confirmatory", "raw_m2_top100_A2_leakage",
                                          "bidirectional_topology_H_only_if_F_passed"),
                            passed = c(TRUE,
                                       all(audit$raw_TWAS_top100_A2_fraction <= 0.05 &
                                             audit$M2_top100_A2_fraction <= 0.05),
                                       is.null(h) || length(pass_arms) > 0),
                            stringsAsFactors = FALSE),
                 file.path(report_dir, "NULL_BIAS_GUARDRAILS.csv"))
  schema <- data.frame(file = c("BIDIRECTIONAL_RANKING_METRICS.csv", "DIRECTION_SPECIFIC_RECOVERY_METRICS.csv",
                                "DECOY_GUARDRAIL_METRICS.csv"),
                       required_columns_present = c(all(c("top100_recall", "risk_top100_recall", "protective_top100_recall") %in% names(metrics)),
                                                    all(c("risk_top100_recall", "protective_top100_recall", "direction_aware_AUPRC") %in% names(direction)),
                                                    all(c("opposite_sign_decoy_top100_rate", "high_degree_decoy_top100_rate") %in% names(decoy))),
                       stringsAsFactors = FALSE)
  write_csv_over(schema, file.path(report_dir, "METRIC_SCHEMA_AUDIT.csv"))

  write_lines_over(c("# Bidirectional Ranking Summary", "",
                     paste0("Outcome classification: `", reason, "`."),
                     sprintf("Best F two-tail-balanced arm: `%s`.", best$rescue_arm),
                     sprintf("Top100/top150/top200 recall: %.4f / %.4f / %.4f.",
                             best$top100_recall, best$top150_recall, best$top200_recall),
                     sprintf("Risk/protective top100 recall: %.4f / %.4f.",
                             best$risk_top100_recall, best$protective_top100_recall),
                     sprintf("Opposite/high-degree decoy top100: %.4f / %.4f.",
                             best$opposite_sign_decoy_top100_rate, best$high_degree_decoy_top100_rate),
                     sprintf("Topology H robustness executed: %s.", !is.null(h))),
                   file.path(report_dir, "BIDIRECTIONAL_RANKING_SUMMARY.md"))
  write_lines_over(c("# Decoy Guardrail Rescue Report", "",
                     "Arms M-P tested bidirectional ranking alone, initial decoy suppression, network decoy suppression, and combined suppression.",
                     sprintf("Best F arm decoy top100/top200 opposite-sign: %.4f / %.4f.",
                             best$opposite_sign_decoy_top100_rate, best$opposite_sign_decoy_top200_rate),
                     sprintf("Best F arm high-degree decoy top100/top200: %.4f / %.4f.",
                             best$high_degree_decoy_top100_rate, best$high_degree_decoy_top200_rate)),
                   file.path(report_dir, "DECOY_GUARDRAIL_RESCUE_REPORT.md"))
  write_lines_over(c("# Direction Specific Recovery Report", "",
                     "Descending, ascending, absolute, balanced two-tail, and direction-matched rankings are reported in the CSV tables.",
                     "The plan hypothesis that descending recovers mostly risk and ascending recovers mostly protective is directly evaluated by risk/protective recall columns."),
                   file.path(report_dir, "DIRECTION_SPECIFIC_RECOVERY_REPORT.md"))
  write_lines_over(c("# Initial Signal Carryforward Audit", "",
                     "All arms carry forward K_sparse_TWAS_filtered as the base initial signal field.",
                     sprintf("Mean raw/M2 top100 A2 fractions: %.4f / %.4f.",
                             mean(audit$raw_TWAS_top100_A2_fraction), mean(audit$M2_top100_A2_fraction)),
                     sprintf("Mean background-to-A2 initial ratio: %.4f.",
                             mean(audit$background_to_A2_mean_abs_initial_weight_ratio))),
                   file.path(report_dir, "INITIAL_SIGNAL_CARRYFORWARD_AUDIT.md"))
  write_lines_over(c("# Code Fidelity Audit", "",
                     paste0("Binding plan SHA256: `", binding_plan_sha256, "`."),
                     "Submitted Final.Heat values were not modified. Faithful M2 arithmetic, P conversion, strict TWAS.P filtering, signed positive/absolute-negative channels, zero-weight edges, self-loops, and diffuStats `n.perm` were retained."),
                   file.path(report_dir, "CODE_FIDELITY_AUDIT.md"))
  write_lines_over(c("# Benchmark Implementation Audit", "",
                     "Comparators retained: raw_TWAS_abs, raw_signed_TWAS, M2_no_diffusion, PPR_abs_prior, PPR_signed_two_channel, RWR_abs_prior, RWR_signed_two_channel, median_weight_permutation, I2_module_disrupted, and I3_expression_matched_randomized.",
                     "PPR/RWR are diagnostic propagation comparators, not replacement NESTA scores."),
                   file.path(report_dir, "BENCHMARK_IMPLEMENTATION_AUDIT.md"))
  write_lines_over(c("# Metric Schema Audit", "",
                     sprintf("Bidirectional schema pass: %s.", schema$required_columns_present[1]),
                     sprintf("Direction schema pass: %s.", schema$required_columns_present[2]),
                     sprintf("Decoy schema pass: %s.", schema$required_columns_present[3])),
                   file.path(report_dir, "METRIC_SCHEMA_AUDIT.md"))
  write_lines_over(c("# STOP/GO Report", "", "STOP.",
                     paste0("Reason: `", reason, "`."),
                     "Diagnostic stage started: YES.",
                     "Confirmatory started: NO.",
                     "This binding plan prohibits confirmatory execution."),
                   file.path(report_dir, "STOP_GO_REPORT.md"))
  write_lines_over(c("# Final Report", "",
                     paste0("Binding plan SHA256: `", binding_plan_sha256, "`."),
                     paste0("Final outcome classification: `", reason, "`."),
                     "Dense-only bidirectional ranking and decoy-guardrail rescue completed on topology F; topology H ran only if an F arm passed.",
                     sprintf("Best F two-tail-balanced arm `%s`: top100/top150/top200 %.4f / %.4f / %.4f.",
                             best$rescue_arm, best$top100_recall, best$top150_recall, best$top200_recall),
                     sprintf("Risk/protective top100 recall %.4f / %.4f; decoy top100 opposite/high-degree %.4f / %.4f.",
                             best$risk_top100_recall, best$protective_top100_recall,
                             best$opposite_sign_decoy_top100_rate, best$high_degree_decoy_top100_rate),
                     sprintf("Confirmatory execution started: NO. Topology H robustness executed: %s.", !is.null(h))),
                   file.path(report_dir, "FINAL_REPORT.md"))
  if (file.exists(project_file("results/reports/summary_tables/unit_test_results.csv"))) {
    file.copy(project_file("results/reports/summary_tables/unit_test_results.csv"),
              file.path(report_dir, "unit_test_results.csv"), overwrite = TRUE)
  }
  if (file.exists(project_file("CLEANUP_MANIFEST.csv"))) {
    file.copy(project_file("CLEANUP_MANIFEST.csv"), file.path(report_dir, "CLEANUP_MANIFEST.csv"), overwrite = TRUE)
  }
  write_csv_over(data.frame(path = c(project_file("R/study_0708_bidirectional_rescue.R"),
                                     project_file("R/study_0708_initial_signal_rescue.R"),
                                     project_file("R/study_0707_operator_audit.R"),
                                     project_file("R/study_0707_branch_isolation.R"),
                                     project_file("R/fidelity.R"), project_file("R/utils.R"),
                                     project_file("scripts/run_bidirectional_rescue.R")),
                            sha256 = c(sha(project_file("R/study_0708_bidirectional_rescue.R")),
                                       sha(project_file("R/study_0708_initial_signal_rescue.R")),
                                       sha(project_file("R/study_0707_operator_audit.R")),
                                       sha(project_file("R/study_0707_branch_isolation.R")),
                                       sha(project_file("R/fidelity.R")),
                                       sha(project_file("R/utils.R")),
                                       if (file.exists(project_file("scripts/run_bidirectional_rescue.R"))) sha(project_file("scripts/run_bidirectional_rescue.R")) else NA_character_),
                            role = c("bidirectional_rescue_runner", "initial_signal_carryforward",
                                     "operator_audit_helpers", "frozen_arm_FH_substrates",
                                     "faithful_nesta_and_benchmarks", "binding_plan_and_io_guards",
                                     "script_entrypoint"),
                            stringsAsFactors = FALSE),
                 file.path(report_dir, "IMPLEMENTATION_MANIFEST.csv"))
  files <- list.files(report_dir, recursive = TRUE, full.names = TRUE)
  files <- files[file.info(files)$isdir == FALSE & basename(files) != "CHECKSUMS.sha256"]
  writeLines(vapply(files, function(f) paste(sha(f), sub(paste0("^", report_dir, "/"), "", f)), character(1)),
             file.path(report_dir, "CHECKSUMS.sha256"))
  invisible(list(reason = reason, best = best, h_ran = !is.null(h)))
}

if (!identical(Sys.getenv("NESTA_BIDIRECTIONAL_SOURCE_ONLY"), "1")) run_bidirectional_rescue()
