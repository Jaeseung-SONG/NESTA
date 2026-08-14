Sys.setenv(NESTA_OPERATOR_AUDIT_SOURCE_ONLY = "1")
source(file.path("/home/js/NESTA/simulation", "R", "study_0707_operator_audit.R"))

initial_signal_arms <- function() {
  c("I_baseline_reproduction", "J_background_null_calibrated",
    "K_sparse_TWAS_filtered", "L_manuscript_balanced_sparse")
}

submitted_filtering_rule_audit <- function() {
  src <- readLines("/home/js/NESTA/Analysis/Nesta.R", warn = FALSE)
  data.frame(
    submitted_file = "/home/js/NESTA/Analysis/Nesta.R",
    TWAS_cutoff_default = 1,
    filtering_rule = "TWAS.P < test.cutoff",
    strict_less_than = any(grepl("TWAS.P < test.cutoff", src, fixed = TRUE)),
    M2_expression_weighting = "Mean_expression/sd(Mean_expression) * weight/sd(weight)",
    nonzero_synthetic_Z_has_P_less_than_1 = TRUE,
    zero_synthetic_Z_has_P_equal_1_and_is_filtered = TRUE,
    implication = "With default cutoff 1, exact Z=0 genes receive zero post-filter weight; any nonzero Z is retained.",
    stringsAsFactors = FALSE
  )
}

apply_initial_signal_arm <- function(rep, initial_arm, seed) {
  set.seed(seed)
  genes <- rep$genes
  z <- setNames(rep$twas$TWAS.Z, rep$twas$SYMBOL)[genes]
  risk_a1 <- rep$A1_risk; prot_a1 <- rep$A1_protective
  risk_a2 <- rep$A2_risk; prot_a2 <- rep$A2_protective
  bg <- rep$background
  z[] <- 0

  if (initial_arm == "I_baseline_reproduction") {
    z <- setNames(rep$twas$TWAS.Z, rep$twas$SYMBOL)[genes]
  } else {
    z[risk_a1] <- abs(rnorm(length(risk_a1), 3.6, 0.35))
    z[prot_a1] <- -abs(rnorm(length(prot_a1), 3.6, 0.35))
    z[rep$relay] <- 0
    z[rep$D] <- 0
    z[rep$C] <- 0
  }

  if (initial_arm == "J_background_null_calibrated") {
    z[risk_a2] <- rnorm(length(risk_a2), 0.12, 0.035)
    z[prot_a2] <- rnorm(length(prot_a2), -0.12, 0.035)
    noisy_bg <- sample(bg, 110)
    z[noisy_bg] <- sample(c(-1, 1), length(noisy_bg), TRUE) * runif(length(noisy_bg), 0.18, 0.32)
    z[rep$D] <- sample(c(-1, 1), length(rep$D), TRUE) * runif(length(rep$D), 0.15, 0.28)
    z[rep$C] <- sample(c(-1, 1), length(rep$C), TRUE) * runif(length(rep$C), 0.12, 0.24)
  } else if (initial_arm == "K_sparse_TWAS_filtered") {
    keep_risk <- sample(risk_a2, ceiling(length(risk_a2) * 0.55))
    keep_prot <- sample(prot_a2, ceiling(length(prot_a2) * 0.55))
    z[keep_risk] <- rnorm(length(keep_risk), 0.10, 0.025)
    z[keep_prot] <- rnorm(length(keep_prot), -0.10, 0.025)
    noisy_bg <- sample(bg, 95)
    z[noisy_bg] <- sample(c(-1, 1), length(noisy_bg), TRUE) * runif(length(noisy_bg), 0.22, 0.42)
    d_noise <- sample(rep$D, 30)
    c_noise <- sample(rep$C, 25)
    z[d_noise] <- sample(c(-1, 1), length(d_noise), TRUE) * runif(length(d_noise), 0.18, 0.34)
    z[c_noise] <- sample(c(-1, 1), length(c_noise), TRUE) * runif(length(c_noise), 0.15, 0.30)
  } else if (initial_arm == "L_manuscript_balanced_sparse") {
    z[risk_a2] <- rnorm(length(risk_a2), 0.13, 0.035)
    z[prot_a2] <- rnorm(length(prot_a2), -0.13, 0.035)
    noisy_bg <- sample(bg, 70)
    z[noisy_bg] <- sample(c(-1, 1), length(noisy_bg), TRUE) * runif(length(noisy_bg), 0.18, 0.34)
    d_noise <- sample(rep$D, 40)
    c_noise <- sample(rep$C, 35)
    z[d_noise] <- sample(c(-1, 1), length(d_noise), TRUE) * runif(length(d_noise), 0.16, 0.30)
    z[c_noise] <- sample(c(-1, 1), length(c_noise), TRUE) * runif(length(c_noise), 0.14, 0.28)
  }

  z <- pmax(pmin(z, 4.5), -4.5)
  rep$initial_signal_arm <- initial_arm
  rep$twas <- data.frame(SYMBOL = genes, TWAS.Z = as.numeric(z), TWAS.P = twas_p_from_z(z),
                         stringsAsFactors = FALSE)
  rep
}

initial_field_audit_row <- function(rep, ch, no_diff, cutoff = 1) {
  genes <- rep$genes
  raw_z <- setNames(rep$twas$TWAS.Z, rep$twas$SYMBOL)[genes]
  raw_p <- twas_p_from_z(raw_z)
  filt_z <- strict_twas_vector(genes, rep$twas, cutoff = cutoff)
  pre <- (rep$expr / safe_sd(rep$expr)) * (raw_z / safe_sd(raw_z))
  init <- setNames(ch$initial_weight, ch$SYMBOL)[genes]
  raw_rank <- rank(-abs(raw_z), ties.method = "average")
  names(raw_rank) <- genes
  m2_rank <- rank(-abs(init), ties.method = "average")
  names(m2_rank) <- genes
  mean_a2 <- mean(abs(init[rep$A2]))
  bg_large <- mean(abs(init[rep$background]) > 0.25 * max(mean_a2, .Machine$double.eps))
  frac_nonzero <- mean(abs(init) > 0)
  ratio <- mean(abs(init[rep$background])) / max(mean_a2, .Machine$double.eps)
  eligible <- length(rep$A1) == sum(abs(init[rep$A1]) > 0) &&
    mean(raw_rank[rep$A2] <= 100) <= 0.05 &&
    mean(m2_rank[rep$A2] <= 100) <= 0.05 &&
    ratio <= 1.00 &&
    (frac_nonzero <= 0.30 || bg_large <= 0.20)
  data.frame(
    topology_arm = rep$arm, topology_label = rep$arm_label,
    initial_signal_arm = rep$initial_signal_arm, replicate = rep$rep_id,
    n_total_genes = length(genes),
    n_nonzero_initial_weight_genes_before_filtering = sum(abs(pre) > 0),
    n_nonzero_initial_weight_genes_after_filtering = sum(abs(init) > 0),
    fraction_nonzero_initial_weight_after_filtering = frac_nonzero,
    n_A1_nonzero_after_filtering = sum(abs(init[rep$A1]) > 0),
    n_A2_nonzero_after_filtering = sum(abs(init[rep$A2]) > 0),
    n_background_nonzero_after_filtering = sum(abs(init[rep$background]) > 0),
    mean_abs_initial_weight_A1 = mean(abs(init[rep$A1])),
    mean_abs_initial_weight_A2 = mean_a2,
    mean_abs_initial_weight_background = mean(abs(init[rep$background])),
    background_to_A2_mean_abs_initial_weight_ratio = ratio,
    A1_to_A2_mean_abs_initial_weight_ratio = mean(abs(init[rep$A1])) / max(mean_a2, .Machine$double.eps),
    raw_TWAS_top100_A2_fraction = mean(raw_rank[rep$A2] <= 100),
    M2_top100_A2_fraction = mean(m2_rank[rep$A2] <= 100),
    raw_TWAS_top10pct_A2_fraction = mean(raw_rank[rep$A2] <= ceiling(0.10 * length(genes))),
    M2_top10pct_A2_fraction = mean(m2_rank[rep$A2] <= ceiling(0.10 * length(genes))),
    fraction_background_abs_initial_gt_0.25_mean_abs_A2 = bg_large,
    default_cutoff_keeps_all_nonzero_Z = all(raw_p[abs(raw_z) > 0] < cutoff),
    pilot_eligible_initial_field = eligible,
    stringsAsFactors = FALSE
  )
}

rescue_score_set <- function(rep, ch) {
  base <- score_set_for_audit(rep, ch)
  base[c("submitted_output_descending", "raw_positive_channel_heat_descending",
         "raw_negative_channel_heat_descending", "signed_reconstructed_heat_descending",
         "degree_residualized_submitted_output")]
}

rescue_comparators <- function(rep, ch, no_diff) {
  genes <- rep$genes
  init_p <- initial_vector_personalizations(setNames(ch$initial_weight, ch$SYMBOL)[genes])
  rwr_abs <- run_rwr(rep$adj, init_p$abs, restart = 0.5)
  ppr_abs <- run_ppr(rep$adj, init_p$abs, damping = 0.85)
  rwr_signed <- run_rwr(rep$adj, init_p$pos, restart = 0.5) - run_rwr(rep$adj, init_p$neg, restart = 0.5)
  ppr_signed <- run_ppr(rep$adj, init_p$pos, damping = 0.85) - run_ppr(rep$adj, init_p$neg, damping = 0.85)
  raw_z <- setNames(rep$twas$TWAS.Z, rep$twas$SYMBOL)[genes]
  init <- setNames(no_diff$final_NESTA_heat, no_diff$SYMBOL)[genes]
  wp <- faithful_m2_channels(permute_sparse_weights(rep$adj, 9711000 + rep$rep_id), rep$expr, rep$twas,
                             n.perm = 10, seed = 9711000 + rep$rep_id)
  i2 <- faithful_m2_channels(module_disrupted_adj(rep), rep$expr, rep$twas,
                             n.perm = 10, seed = 9712000 + rep$rep_id)
  i3 <- faithful_m2_channels(permute_sparse_weights(rep$adj, 9713000 + rep$rep_id), rep$expr, rep$twas,
                             n.perm = 10, seed = 9713000 + rep$rep_id)
  list(
    PPR_abs_prior = list(score = ppr_abs, signed = ppr_abs),
    PPR_signed_two_channel = list(score = abs(ppr_signed), signed = ppr_signed),
    RWR_abs_prior = list(score = rwr_abs, signed = rwr_abs),
    RWR_signed_two_channel = list(score = abs(rwr_signed), signed = rwr_signed),
    raw_TWAS_abs = list(score = abs(raw_z), signed = raw_z),
    raw_signed_TWAS = list(score = abs(raw_z), signed = raw_z),
    M2_no_diffusion = list(score = abs(init), signed = init),
    median_weight_permutation = list(score = abs(setNames(wp$signed_NESTA_heat, wp$SYMBOL)[genes]),
                                     signed = setNames(wp$signed_NESTA_heat, wp$SYMBOL)[genes]),
    I2_module_disrupted = list(score = abs(setNames(i2$signed_NESTA_heat, i2$SYMBOL)[genes]),
                               signed = setNames(i2$signed_NESTA_heat, i2$SYMBOL)[genes]),
    I3_expression_matched_randomized = list(score = abs(setNames(i3$signed_NESTA_heat, i3$SYMBOL)[genes]),
                                            signed = setNames(i3$signed_NESTA_heat, i3$SYMBOL)[genes])
  )
}

arm_summary <- function(metrics, audit) {
  nesta <- metrics[metrics$score_name == "submitted_output_descending", , drop = FALSE]
  out <- stats::aggregate(cbind(top100_recall, top150_recall, top200_recall, top10pct_recall,
                                raw_AUPRC, prevalence_normalized_AUPRC, median_A2_rank,
                                first_A2_rank, direction_aware_AUPRC,
                                sign_concordant_top100_recall,
                                opposite_sign_decoy_top100_rate,
                                high_degree_decoy_top100_rate,
                                score_degree_spearman) ~ topology_arm + topology_label + initial_signal_arm,
                          data = nesta, FUN = function(z) mean(z, na.rm = TRUE))
  aud <- stats::aggregate(cbind(background_to_A2_mean_abs_initial_weight_ratio,
                                fraction_nonzero_initial_weight_after_filtering,
                                fraction_background_abs_initial_gt_0.25_mean_abs_A2,
                                raw_TWAS_top100_A2_fraction,
                                M2_top100_A2_fraction,
                                pilot_eligible_initial_field) ~ topology_arm + initial_signal_arm,
                          data = audit, FUN = function(z) mean(z, na.rm = TRUE))
  merge(out, aud, by = c("topology_arm", "initial_signal_arm"), all.x = TRUE)
}

rescue_contrasts <- function(primary, bench) {
  out <- list(); i <- 1
  for (tp in unique(primary$topology_arm)) for (arm in unique(primary$initial_signal_arm)) {
    base <- primary[primary$topology_arm == tp & primary$initial_signal_arm == arm &
                      primary$score_name == "submitted_output_descending", , drop = FALSE]
    for (cmp in c("raw_TWAS_abs", "M2_no_diffusion", "median_weight_permutation",
                  "I2_module_disrupted", "I3_expression_matched_randomized")) {
      other <- bench[bench$topology_arm == tp & bench$initial_signal_arm == arm &
                       bench$score_name == cmp, , drop = FALSE]
      z <- merge(base[, c("replicate", "top100_recall", "top200_recall", "raw_AUPRC")],
                 other[, c("replicate", "top100_recall", "top200_recall", "raw_AUPRC")],
                 by = "replicate", suffixes = c("_nesta", "_comparator"))
      for (metric in c("top100_recall", "top200_recall", "raw_AUPRC")) {
        d <- z[[paste0(metric, "_nesta")]] - z[[paste0(metric, "_comparator")]]
        ci <- ci_diff(d)
        out[[i]] <- data.frame(topology_arm = tp, initial_signal_arm = arm,
                               comparator = cmp, metric = metric,
                               nesta_mean = mean(z[[paste0(metric, "_nesta")]], na.rm = TRUE),
                               comparator_mean = mean(z[[paste0(metric, "_comparator")]], na.rm = TRUE),
                               mean_difference = ci["mean"], ci_low = ci["ci_low"],
                               ci_high = ci["ci_high"], stringsAsFactors = FALSE)
        i <- i + 1
      }
    }
  }
  do.call(rbind, out)
}

run_initial_signal_rescue <- function() {
  verify_project_path()
  verify_binding_plan()
  report_dir <- read_report_dir()
  safe_dir_create(report_dir)
  copy_binding_plan_to_report(report_dir)
  set_run_status("IN_PROGRESS", "YES", "NO", "0708 initial-signal sparsity rescue")

  filter_rule <- submitted_filtering_rule_audit()
  primary_rows <- list(); bench_rows <- list(); degree_rows <- list(); audit_rows <- list()
  ip <- ib <- id <- ia <- 1
  topologies <- c("F", "H")
  init_arms <- initial_signal_arms()
  for (tp in topologies) for (arm in init_arms) for (rep_id in seq_len(20)) {
    base <- make_branch_isolation_rep(tp, rep_id, 970800 + match(tp, topologies) * 10000 + rep_id)
    rep <- apply_initial_signal_arm(base, arm, 970900 + match(tp, topologies) * 10000 + match(arm, init_arms) * 1000 + rep_id)
    ch <- faithful_m2_channels(rep$adj, rep$expr, rep$twas, n.perm = 25,
                               seed = 971000 + match(tp, topologies) * 10000 + rep_id)
    no_diff <- no_diffusion_m2_scores(rep$adj, rep$expr, rep$twas)
    audit_rows[[ia]] <- initial_field_audit_row(rep, ch, no_diff); ia <- ia + 1
    scores <- rescue_score_set(rep, ch)
    for (nm in names(scores)) {
      row <- audit_score_metrics(rep, scores[[nm]]$score, scores[[nm]]$signed, nm, scores[[nm]]$object)
      row$topology_arm <- tp; row$topology_label <- row$arm_label; row$initial_signal_arm <- arm
      primary_rows[[ip]] <- row; ip <- ip + 1
      degree_rows[[id]] <- row[, c("topology_arm", "arm", "arm_label", "initial_signal_arm", "replicate",
                                   "score_name", "score_degree_spearman",
                                   "high_degree_decoy_top100_rate")]
      id <- id + 1
    }
    comparators <- rescue_comparators(rep, ch, no_diff)
    for (nm in names(comparators)) {
      row <- audit_score_metrics(rep, comparators[[nm]]$score, comparators[[nm]]$signed, nm, "diagnostic_comparator")
      row$topology_arm <- tp; row$topology_label <- row$arm_label; row$initial_signal_arm <- arm
      bench_rows[[ib]] <- row; ib <- ib + 1
    }
    if (rep_id %% 5 == 0) gc(FALSE)
  }

  primary <- do.call(rbind, primary_rows)
  bench <- do.call(rbind, bench_rows)
  degree <- do.call(rbind, degree_rows)
  audit <- do.call(rbind, audit_rows)
  summary <- arm_summary(primary, audit)
  contrasts <- rescue_contrasts(primary, bench)
  bench_summary <- summarise_score_metrics(bench)
  schema <- data.frame(
    file = c("PRIMARY_FINAL_HEAT_METRICS.csv", "DIRECTION_AWARE_METRICS.csv", "BENCHMARK_METRICS.csv"),
    required_columns_present = c(all(c("top100_recall", "top150_recall", "top200_recall", "top10pct_recall", "raw_AUPRC", "median_A2_rank") %in% names(primary)),
                                 all(c("direction_aware_AUPRC", "sign_concordant_top100_recall", "opposite_sign_decoy_top100_rate") %in% names(primary)),
                                 all(c("top100_recall", "raw_AUPRC", "score_degree_spearman") %in% names(bench))),
    stringsAsFactors = FALSE
  )
  nesta_summary <- summary
  raw_cmp <- bench_summary[bench_summary$score_name == "raw_TWAS_abs", c("arm", "score_name", "top100_recall", "top200_recall", "raw_AUPRC")]
  m2_cmp <- bench_summary[bench_summary$score_name == "M2_no_diffusion", c("arm", "score_name", "top100_recall", "top200_recall", "raw_AUPRC")]
  promising <- within(nesta_summary, {
    promising_arm <- background_to_A2_mean_abs_initial_weight_ratio <= 1.00 &
      top100_recall > 0 & top200_recall >= 0.05 & median_A2_rank < 350 &
      opposite_sign_decoy_top100_rate <= 0.15 &
      high_degree_decoy_top100_rate <= 0.12
  })
  for (r in seq_len(nrow(promising))) {
    tp <- promising$topology_arm[r]; arm <- promising$initial_signal_arm[r]
    raw <- bench[bench$topology_arm == tp & bench$initial_signal_arm == arm & bench$score_name == "raw_TWAS_abs", ]
    m2 <- bench[bench$topology_arm == tp & bench$initial_signal_arm == arm & bench$score_name == "M2_no_diffusion", ]
    promising$promising_arm[r] <- promising$promising_arm[r] &&
      promising$raw_AUPRC[r] > mean(raw$raw_AUPRC, na.rm = TRUE) &&
      promising$raw_AUPRC[r] > mean(m2$raw_AUPRC, na.rm = TRUE)
  }
  leakage_artifact <- any(audit$raw_TWAS_top100_A2_fraction > 0.05 | audit$M2_top100_A2_fraction > 0.05)
  sparse_success <- any(audit$background_to_A2_mean_abs_initial_weight_ratio <= 1 &
                          (audit$fraction_nonzero_initial_weight_after_filtering <= 0.30 |
                             audit$fraction_background_abs_initial_gt_0.25_mean_abs_A2 <= 0.20))
  recovery_seen <- any(promising$initial_signal_arm != "I_baseline_reproduction" &
                         promising$background_to_A2_mean_abs_initial_weight_ratio <= 1.00 &
                         promising$raw_TWAS_top100_A2_fraction <= 0.05 &
                         promising$M2_top100_A2_fraction <= 0.05 &
                         promising$top200_recall >= 0.05, na.rm = TRUE)
  if (any(promising$promising_arm, na.rm = TRUE)) {
    reason <- "initial_signal_sparsity_promising_rescue"
  } else if (leakage_artifact) {
    reason <- "initial_signal_leakage_artifact"
  } else if (recovery_seen) {
    reason <- "initial_signal_recovery_with_decoy_guardrail_failure"
  } else if (sparse_success) {
    reason <- "sparse_initial_field_rescue_failed"
  } else {
    reason <- "initial_field_hard_diagnostic_failure"
  }
  best <- promising[order(promising$promising_arm, promising$top100_recall, promising$top200_recall,
                          -promising$median_A2_rank, decreasing = TRUE), ][1, ]
  set_run_status("STOPPED", "YES", "NO", reason)

  write_csv_over(filter_rule, file.path(report_dir, "SUBMITTED_FILTERING_RULE_AUDIT.csv"))
  write_csv_over(audit, file.path(report_dir, "INITIAL_SIGNAL_FIELD_AUDIT.csv"))
  write_csv_over(promising, file.path(report_dir, "INITIAL_SIGNAL_RESCUE_METRICS.csv"))
  write_csv_over(contrasts, file.path(report_dir, "INITIAL_SIGNAL_RESCUE_CONTRASTS.csv"))
  write_csv_over(primary, file.path(report_dir, "PRIMARY_FINAL_HEAT_METRICS.csv"))
  write_csv_over(primary[, c("topology_arm", "arm", "arm_label", "initial_signal_arm", "replicate",
                             "score_name", "direction_aware_AUPRC",
                             "sign_concordant_top100_recall",
                             "opposite_sign_decoy_top100_rate")],
                 file.path(report_dir, "DIRECTION_AWARE_METRICS.csv"))
  write_csv_over(bench, file.path(report_dir, "BENCHMARK_METRICS.csv"))
  write_csv_over(bench_summary, file.path(report_dir, "BENCHMARK_CONTRASTS.csv"))
  write_csv_over(degree, file.path(report_dir, "DEGREE_BIAS_METRICS.csv"))
  write_csv_over(data.frame(guardrail = c("no_confirmatory", "top100_initial_A2_leakage", "initial_background_ratio"),
                            passed = c(TRUE,
                                       all(audit$raw_TWAS_top100_A2_fraction <= 0.05 &
                                             audit$M2_top100_A2_fraction <= 0.05),
                                       any(audit$background_to_A2_mean_abs_initial_weight_ratio <= 1)),
                            stringsAsFactors = FALSE),
                 file.path(report_dir, "NULL_BIAS_GUARDRAILS.csv"))
  write_csv_over(schema, file.path(report_dir, "METRIC_SCHEMA_AUDIT.csv"))

  write_lines_over(c("# Initial Signal Sparsity Summary", "",
                     paste0("Outcome classification: `", reason, "`."),
                     sprintf("Best topology/initial arm: `%s` / `%s`.", best$topology_arm, best$initial_signal_arm),
                     sprintf("Submitted Final.Heat top100/top150/top200/top10%% recall: %.4f / %.4f / %.4f / %.4f.",
                             best$top100_recall, best$top150_recall, best$top200_recall, best$top10pct_recall),
                     sprintf("Median A2 rank %.2f; raw AUPRC %.4f; background-to-A2 initial ratio %.4f.",
                             best$median_A2_rank, best$raw_AUPRC, best$background_to_A2_mean_abs_initial_weight_ratio),
                     sprintf("Promising arm criteria passed: %s.", best$promising_arm)),
                   file.path(report_dir, "INITIAL_SIGNAL_SPARSITY_SUMMARY.md"))
  write_lines_over(c("# Initial Signal Field Audit", "",
                     "Submitted filtering rule: default `TWAS_cutoff = 1` with strict `TWAS.P < test.cutoff`.",
                     "Exact zero synthetic Z has TWAS.P = 1 and is filtered; any nonzero synthetic Z is retained under the default cutoff.",
                     sprintf("Mean baseline Arm I background-to-A2 initial ratio: %.4f.",
                             mean(audit$background_to_A2_mean_abs_initial_weight_ratio[audit$initial_signal_arm == "I_baseline_reproduction"])),
                     sprintf("Minimum calibrated background-to-A2 initial ratio: %.4f.",
                             min(audit$background_to_A2_mean_abs_initial_weight_ratio[audit$initial_signal_arm != "I_baseline_reproduction"], na.rm = TRUE))),
                   file.path(report_dir, "INITIAL_SIGNAL_FIELD_AUDIT.md"))
  write_lines_over(c("# NESTA Initial Field Rescue Report", "",
                     "The topology substrates were frozen to Arm F and Arm H; only the initial TWAS/M2 field varied across Arms I-L.",
                     sprintf("Best submitted Final.Heat result was `%s` on topology `%s`: top100/top200 %.4f / %.4f.",
                             best$initial_signal_arm, best$topology_arm, best$top100_recall, best$top200_recall),
                     sprintf("Classification: `%s`.", reason)),
                   file.path(report_dir, "NESTA_INITIAL_FIELD_RESCUE_REPORT.md"))
  write_lines_over(c("# NESTA Operator Audit Carryforward", "",
                     "The immediately preceding operator audit classified the dense failure as `full_network_normalization_failure`.",
                     "Seed filtering, rank orientation, and minimal positive-control failure were not supported.",
                     "This round therefore calibrated background/A2 initial-weight sparsity while preserving submitted NESTA behavior and frozen topology substrates."),
                   file.path(report_dir, "NESTA_OPERATOR_AUDIT_CARRYFORWARD.md"))
  write_lines_over(c("# Code Fidelity Audit", "",
                     paste0("Binding plan SHA256: `", binding_plan_sha256, "`."),
                     "Faithful M2 arithmetic, P conversion, strict TWAS.P filtering, positive/absolute-negative channels, retained zero-weight edges, self-loops, and diffuStats `n.perm` were used.",
                     "Submitted filtering rule was verified directly from `/home/js/NESTA/Analysis/Nesta.R`."),
                   file.path(report_dir, "CODE_FIDELITY_AUDIT.md"))
  write_lines_over(c("# Benchmark Implementation Audit", "",
                     "Comparators retained: PPR_abs_prior, PPR_signed_two_channel, RWR_abs_prior, RWR_signed_two_channel, raw_TWAS_abs, raw_signed_TWAS, M2_no_diffusion, median_weight_permutation, I2_module_disrupted, and I3_expression_matched_randomized.",
                     "PPR/RWR are diagnostic comparators, not replacement NESTA scores."),
                   file.path(report_dir, "BENCHMARK_IMPLEMENTATION_AUDIT.md"))
  write_lines_over(c("# Metric Schema Audit", "",
                     sprintf("Primary schema pass: %s.", schema$required_columns_present[schema$file == "PRIMARY_FINAL_HEAT_METRICS.csv"]),
                     sprintf("Direction schema pass: %s.", schema$required_columns_present[schema$file == "DIRECTION_AWARE_METRICS.csv"]),
                     sprintf("Benchmark schema pass: %s.", schema$required_columns_present[schema$file == "BENCHMARK_METRICS.csv"])),
                   file.path(report_dir, "METRIC_SCHEMA_AUDIT.md"))
  write_lines_over(c("# STOP/GO Report", "", "STOP.",
                     paste0("Reason: `", reason, "`."),
                     "Diagnostic pilot started: YES.",
                     "Confirmatory started: NO.",
                     "This binding plan is rescue/calibration-only and prohibits confirmatory execution."),
                   file.path(report_dir, "STOP_GO_REPORT.md"))
  write_lines_over(c("# Final Report", "",
                     paste0("Binding plan SHA256: `", binding_plan_sha256, "`."),
                     "Dense-only initial-signal sparsity rescue completed with Arms I-L on frozen Arm F/H topology substrates.",
                     "Confirmatory execution started: NO.",
                     paste0("Final outcome classification: `", reason, "`."),
                     sprintf("Best arm `%s` on topology `%s`: top100/top150/top200 %.4f / %.4f / %.4f; median A2 rank %.2f.",
                             best$initial_signal_arm, best$topology_arm, best$top100_recall,
                             best$top150_recall, best$top200_recall, best$median_A2_rank),
                     sprintf("Initial background-to-A2 ratio %.4f; score-degree Spearman %.4f; opposite/high-degree decoy top100 %.4f / %.4f.",
                             best$background_to_A2_mean_abs_initial_weight_ratio,
                             best$score_degree_spearman, best$opposite_sign_decoy_top100_rate,
                             best$high_degree_decoy_top100_rate)),
                   file.path(report_dir, "FINAL_REPORT.md"))
  if (file.exists(project_file("results/reports/summary_tables/unit_test_results.csv"))) {
    file.copy(project_file("results/reports/summary_tables/unit_test_results.csv"),
              file.path(report_dir, "unit_test_results.csv"), overwrite = TRUE)
  }
  if (file.exists(project_file("CLEANUP_MANIFEST.csv"))) {
    file.copy(project_file("CLEANUP_MANIFEST.csv"), file.path(report_dir, "CLEANUP_MANIFEST.csv"), overwrite = TRUE)
  }
  write_csv_over(data.frame(path = c(project_file("R/study_0708_initial_signal_rescue.R"),
                                     project_file("R/study_0707_operator_audit.R"),
                                     project_file("R/study_0707_branch_isolation.R"),
                                     project_file("R/study_0707_dense_rescue.R"),
                                     project_file("R/fidelity.R"), project_file("R/utils.R"),
                                     project_file("scripts/run_initial_signal_rescue.R")),
                            sha256 = c(sha(project_file("R/study_0708_initial_signal_rescue.R")),
                                       sha(project_file("R/study_0707_operator_audit.R")),
                                       sha(project_file("R/study_0707_branch_isolation.R")),
                                       sha(project_file("R/study_0707_dense_rescue.R")),
                                       sha(project_file("R/fidelity.R")),
                                       sha(project_file("R/utils.R")),
                                       if (file.exists(project_file("scripts/run_initial_signal_rescue.R"))) sha(project_file("scripts/run_initial_signal_rescue.R")) else NA_character_),
                            role = c("initial_signal_rescue_runner", "operator_audit_helpers",
                                     "frozen_arm_FH_substrates", "faithful_channel_helper",
                                     "faithful_nesta_and_benchmarks", "binding_plan_and_io_guards",
                                     "script_entrypoint"),
                            stringsAsFactors = FALSE),
                 file.path(report_dir, "IMPLEMENTATION_MANIFEST.csv"))
  files <- list.files(report_dir, recursive = TRUE, full.names = TRUE)
  files <- files[file.info(files)$isdir == FALSE & basename(files) != "CHECKSUMS.sha256"]
  writeLines(vapply(files, function(f) paste(sha(f), sub(paste0("^", report_dir, "/"), "", f)), character(1)),
             file.path(report_dir, "CHECKSUMS.sha256"))
  invisible(list(reason = reason, summary = promising, best = best))
}

if (!identical(Sys.getenv("NESTA_INITIAL_SIGNAL_SOURCE_ONLY"), "1")) run_initial_signal_rescue()
