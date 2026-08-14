Sys.setenv(NESTA_DENSE_RESCUE_SOURCE_ONLY = "1")
source(file.path("/home/js/NESTA/simulation", "R", "study_0707_dense_rescue.R"))

branch_isolation_arms <- function() {
  list(
    E = list(label = "Arm_E_branch_isolated_Arm_D", bg_edges = 1200, branch_bg_edges = 12,
             relay_cluster_p = 0.018, relay_contacts = 4, relay_a2 = c(0.205, 0.225),
             a1_relay = c(0.275, 0.290), a2_block_p = 0.58, a2_w = c(0.13, 0.16)),
    F = list(label = "Arm_F_relay_pass_through_Arm_D", bg_edges = 1500, branch_bg_edges = 10,
             relay_cluster_p = 0.005, relay_contacts = 8, relay_a2 = c(0.215, 0.240),
             a1_relay = c(0.260, 0.280), a2_block_p = 0.58, a2_w = c(0.13, 0.16)),
    G = list(label = "Arm_G_A2_rank_compression_Arm_D", bg_edges = 1000, branch_bg_edges = 8,
             relay_cluster_p = 0.012, relay_contacts = 7, relay_a2 = c(0.210, 0.235),
             a1_relay = c(0.270, 0.285), a2_block_p = 0.72, a2_w = c(0.15, 0.19)),
    H = list(label = "Arm_H_combined_branch_relay_A2_capture", bg_edges = 700, branch_bg_edges = 4,
             relay_cluster_p = 0.002, relay_contacts = 9, relay_a2 = c(0.225, 0.255),
             a1_relay = c(0.260, 0.275), a2_block_p = 0.76, a2_w = c(0.16, 0.20))
  )
}

make_branch_isolation_rep <- function(arm_id, rep_id, seed) {
  arm <- branch_isolation_arms()[[arm_id]]
  set.seed(seed)
  n <- 1000
  genes <- sprintf("G%05d", seq_len(n))
  a1 <- genes[1:10]
  relay <- genes[11:38]
  a2 <- genes[39:78]
  decoy <- genes[79:138]
  cdecoy <- genes[139:188]
  bg <- genes[189:n]
  risk_a1 <- a1[1:5]
  prot_a1 <- a1[6:10]
  risk_relay <- relay[seq(1, length(relay), 2)]
  prot_relay <- setdiff(relay, risk_relay)
  risk_a2 <- a2[seq(1, length(a2), 2)]
  prot_a2 <- setdiff(a2, risk_a2)
  a2_dir <- c(setNames(rep("risk", length(risk_a2)), risk_a2),
              setNames(rep("protective", length(prot_a2)), prot_a2))
  relay_dir <- c(setNames(rep("risk", length(risk_relay)), risk_relay),
                 setNames(rep("protective", length(prot_relay)), prot_relay))

  ii <- integer(); jj <- integer(); xx <- numeric()
  add <- function(a, b, w) {
    a <- rep(a, length.out = max(length(a), length(b), length(w)))
    b <- rep(b, length.out = length(a)); w <- rep(w, length.out = length(a))
    keep <- !is.na(match(a, genes)) & !is.na(match(b, genes)) & a != b & is.finite(w)
    a <- a[keep]; b <- b[keep]; w <- w[keep]
    ia <- match(a, genes); ib <- match(b, genes)
    ii <<- c(ii, ia, ib); jj <<- c(jj, ib, ia); xx <<- c(xx, w, w)
  }
  add_block <- function(nodes, p, lo, hi) {
    if (length(nodes) < 2) return()
    cmb <- combn(nodes, 2)
    keep <- runif(ncol(cmb)) < p
    cmb <- cmb[, keep, drop = FALSE]
    if (ncol(cmb)) for (k in seq_len(ncol(cmb))) add(cmb[1, k], cmb[2, k], runif(1, lo, hi))
  }

  add_block(c(risk_a1, risk_relay, risk_a2), 0.12, 0.12, 0.17)
  add_block(c(prot_a1, prot_relay, prot_a2), 0.12, 0.12, 0.17)
  add_block(risk_a2, arm$a2_block_p, arm$a2_w[1], arm$a2_w[2])
  add_block(prot_a2, arm$a2_block_p, arm$a2_w[1], arm$a2_w[2])
  add_block(relay, arm$relay_cluster_p, 0.08, 0.11)

  for (r in risk_relay) add(r, sample(risk_a1, 2), runif(2, arm$a1_relay[1], arm$a1_relay[2]))
  for (r in prot_relay) add(r, sample(prot_a1, 2), runif(2, arm$a1_relay[1], arm$a1_relay[2]))
  for (g in risk_a2) add(g, sample(risk_relay, arm$relay_contacts), runif(arm$relay_contacts, arm$relay_a2[1], arm$relay_a2[2]))
  for (g in prot_a2) add(g, sample(prot_relay, arm$relay_contacts), runif(arm$relay_contacts, arm$relay_a2[1], arm$relay_a2[2]))
  for (g in sample(risk_a2, 6)) add(g, sample(risk_a1, 1), runif(1, 0.11, 0.14))
  for (g in sample(prot_a2, 6)) add(g, sample(prot_a1, 1), runif(1, 0.11, 0.14))

  add(sample(bg, arm$bg_edges, TRUE), sample(bg, arm$bg_edges, TRUE), runif(arm$bg_edges, 0.001, 0.035))
  add(sample(c(a1, relay, a2), arm$branch_bg_edges, TRUE), sample(bg, arm$branch_bg_edges, TRUE), runif(arm$branch_bg_edges, 0.001, 0.012))
  add(sample(cdecoy, 120, TRUE), sample(bg, 120, TRUE), runif(120, 0.02, 0.06))

  adj <- sparseMatrix(i = ii, j = jj, x = xx, dims = c(n, n), dimnames = list(genes, genes))
  adj <- forceSymmetric(adj, uplo = "U")
  diag(adj) <- 1

  expr <- setNames(rlnorm(n, 0, 0.30), genes)
  expr[relay] <- pmin(expr[relay], stats::median(expr) * 0.90)
  expr[a2] <- stats::median(expr) * runif(length(a2), 0.92, 1.08)
  expr[cdecoy] <- stats::median(expr) * runif(length(cdecoy), 1.00, 1.12)

  z <- setNames(rnorm(n, 0, 1), genes)
  z[risk_a1] <- abs(rnorm(length(risk_a1), 3.5, 0.45))
  z[prot_a1] <- -abs(rnorm(length(prot_a1), 3.5, 0.45))
  z[relay] <- rnorm(length(relay), 0, 0.10)
  z[risk_a2] <- rnorm(length(risk_a2), 0.16, 0.12)
  z[prot_a2] <- rnorm(length(prot_a2), -0.16, 0.12)
  z[sample(risk_a2, 1)] <- runif(1, 1.05, 1.20)
  z[sample(prot_a2, 1)] <- -runif(1, 1.05, 1.20)
  z[decoy[1:30]] <- -abs(rnorm(30, 1.15, 0.35))
  z[decoy[31:60]] <- abs(rnorm(30, 1.15, 0.35))
  z <- pmax(pmin(z, 4.5), -4.5)

  twas <- data.frame(SYMBOL = genes, TWAS.Z = as.numeric(z), TWAS.P = twas_p_from_z(z),
                     stringsAsFactors = FALSE)
  list(condition = "dense_1000", arm = arm_id, arm_label = arm$label, rep_id = rep_id,
       genes = genes, adj = adj, expr = expr, twas = twas, A1 = a1, A1_risk = risk_a1,
       A1_protective = prot_a1, relay = relay, relay_risk = risk_relay,
       relay_protective = prot_relay, relay_direction = relay_dir,
       A2 = a2, A2_risk = risk_a2, A2_protective = prot_a2, A2_direction = a2_dir,
       D = decoy, C = cdecoy, background = bg)
}

branch_leakage_row <- function(rep) {
  adj <- rep$adj
  branch <- unique(c(rep$A1, rep$relay, rep$A2))
  bg <- rep$background
  seed_cut <- sum(adj[rep$A1, bg])
  seed_total <- sum(adj[rep$A1, setdiff(rep$genes, rep$A1)])
  branch_cut <- sum(adj[branch, bg])
  internal <- sum(adj[branch, branch]) - length(branch)
  bg_sum <- sum(adj[bg, bg]) - length(bg)
  bg_n <- sum(adj[bg, bg] != 0) - length(bg)
  data.frame(arm = rep$arm, arm_label = rep$arm_label, replicate = rep$rep_id,
             seed_neighborhood_background_fraction = seed_cut / max(seed_total, .Machine$double.eps),
             A_branch_background_cut_fraction = branch_cut / max(branch_cut + internal, .Machine$double.eps),
             mean_background_weak_edge_weight = bg_sum / max(bg_n, 1),
             branch_leakage_qc_pass = seed_cut / max(seed_total, .Machine$double.eps) <= 0.010 &&
               branch_cut / max(branch_cut + internal, .Machine$double.eps) <= 0.080,
             stringsAsFactors = FALSE)
}

summarise_branch_metrics <- function(primary, direction, terminal, ranks, leakage, degree) {
  nesta <- primary[primary$score_name == "NESTA_unsigned_final_heat", , drop = FALSE]
  signed <- direction[direction$score_name == "NESTA_signed_reconstructed_heat", , drop = FALSE]
  do.call(rbind, lapply(split(nesta, nesta$arm), function(x) {
    arm <- unique(x$arm)
    sx <- signed[signed$arm == arm, , drop = FALSE]
    tx <- terminal[terminal$arm == arm, , drop = FALSE]
    lx <- leakage[leakage$arm == arm, , drop = FALSE]
    dx <- degree[degree$arm == arm & degree$score_name == "NESTA_unsigned_final_heat", , drop = FALSE]
    data.frame(arm = arm, arm_label = unique(x$arm_label),
               mean_top100_recall = mean(x$top100_recall, na.rm = TRUE),
               mean_top150_recall = mean(x$top150_recall, na.rm = TRUE),
               mean_top200_recall = mean(x$top200_recall, na.rm = TRUE),
               mean_top10pct_recall = mean(x$top10pct_recall, na.rm = TRUE),
               mean_raw_AUPRC = mean(x$raw_AUPRC, na.rm = TRUE),
               direction_aware_AUPRC = mean(sx$direction_aware_AUPRC, na.rm = TRUE),
               sign_concordant_top100_recall = mean(sx$sign_concordant_top100_recall, na.rm = TRUE),
               opposite_sign_decoy_top100_rate = mean(sx$opposite_sign_decoy_top100_rate, na.rm = TRUE),
               relay_top100_rate = mean(tx$relay_top100_rate, na.rm = TRUE),
               high_degree_decoy_top100_rate = mean(dx$high_degree_decoy_top100_rate, na.rm = TRUE),
               score_degree_spearman = mean(dx$score_degree_spearman, na.rm = TRUE),
               relay_to_A2_heat_ratio = mean(tx$relay_to_A2_heat_ratio, na.rm = TRUE),
               retained_in_A_branch = mean(tx$fraction_seed_heat_retained_in_A_branch, na.rm = TRUE),
               heat_reaching_A2 = mean(tx$fraction_seed_heat_reaching_A2, na.rm = TRUE),
               background_leakage = mean(tx$fraction_seed_heat_leaking_to_background, na.rm = TRUE),
               median_A2_rank = mean(tx$median_A2_rank, na.rm = TRUE),
               first_A2_rank = mean(ranks$first_A2_rank[ranks$arm == arm & ranks$score_name == "NESTA_unsigned_final_heat"], na.rm = TRUE),
               seed_neighborhood_background_fraction = mean(lx$seed_neighborhood_background_fraction, na.rm = TRUE),
               A_branch_background_cut_fraction = mean(lx$A_branch_background_cut_fraction, na.rm = TRUE),
               promising = mean(tx$fraction_seed_heat_retained_in_A_branch, na.rm = TRUE) >= 0.10 &&
                 mean(tx$fraction_seed_heat_leaking_to_background, na.rm = TRUE) <= 0.70 &&
                 mean(tx$fraction_seed_heat_reaching_A2, na.rm = TRUE) >= 0.020 &&
                 mean(tx$relay_to_A2_heat_ratio, na.rm = TRUE) <= 1.00 &&
                 mean(tx$median_A2_rank, na.rm = TRUE) < 500 &&
                 mean(x$top200_recall, na.rm = TRUE) > 0.05 &&
                 mean(x$top100_recall, na.rm = TRUE) > 0 &&
                 mean(sx$opposite_sign_decoy_top100_rate, na.rm = TRUE) <= 0.15,
               stringsAsFactors = FALSE)
  }))
}

run_branch_isolation_rescue <- function() {
  verify_project_path()
  verify_binding_plan()
  report_dir <- read_report_dir()
  safe_dir_create(report_dir)
  copy_binding_plan_to_report(report_dir)
  set_run_status("IN_PROGRESS", "YES", "NO", "0707 dense branch-isolation relay-sink rescue")

  smoke <- do.call(rbind, lapply(names(branch_isolation_arms()), function(arm) {
    rows <- lapply(seq_len(2), function(r) {
      rep <- make_branch_isolation_rep(arm, r, 977000 + match(arm, names(branch_isolation_arms())) * 100 + r)
      no_diff <- no_diffusion_m2_scores(rep$adj, rep$expr, rep$twas)
      sig <- initial_signal_ok(rep, no_diff)
      leak <- branch_leakage_row(rep)
      data.frame(arm = arm, replicate = r, raw_TWAS_top100_A2_fraction = sig$raw_top100,
                 M2_initial_top100_A2_fraction = sig$m2_top100,
                 seed_neighborhood_background_fraction = leak$seed_neighborhood_background_fraction,
                 A_branch_background_cut_fraction = leak$A_branch_background_cut_fraction,
                 smoke_pass = sig$raw_top100 <= 0.05 && sig$m2_top100 <= 0.05 &&
                   leak$branch_leakage_qc_pass,
                 stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
  }))

  if (!all(smoke$smoke_pass)) {
    write_csv_over(smoke, file.path(report_dir, "BRANCH_LEAKAGE_AUDIT.csv"))
    reason <- "pre_pilot_branch_leakage_or_initial_signal_failure"
    set_run_status("STOPPED", "NO", "NO", reason)
    write_lines_over(c("# STOP/GO Report", "", "STOP.", paste0("Reason: `", reason, "`."),
                       "Pilot started: NO.", "Confirmatory started: NO."),
                     file.path(report_dir, "STOP_GO_REPORT.md"))
    write_lines_over(c("# Final Report", "", paste0("Binding plan SHA256: `", binding_plan_sha256, "`."),
                       paste0("STOP classification: `", reason, "`.")),
                     file.path(report_dir, "FINAL_REPORT.md"))
    return(invisible(reason))
  }

  primary_rows <- list(); direction_rows <- list(); polarity_rows <- list()
  terminal_rows <- list(); rank_rows <- list(); degree_rows <- list(); leakage_rows <- list(); signal_rows <- list()
  ip <- id <- ipo <- it <- ir <- ideg <- il <- is <- 1
  arms <- names(branch_isolation_arms())
  for (arm in arms) for (rep_id in seq_len(20)) {
    rep <- make_branch_isolation_rep(arm, rep_id, 977700 + match(arm, arms) * 1000 + rep_id)
    ch <- faithful_m2_channels(rep$adj, rep$expr, rep$twas, n.perm = 25, seed = 977700 + rep_id)
    no_diff <- no_diffusion_m2_scores(rep$adj, rep$expr, rep$twas)
    sig <- initial_signal_ok(rep, no_diff)
    signal_rows[[is]] <- data.frame(arm = rep$arm, arm_label = rep$arm_label, replicate = rep$rep_id,
                                    raw_TWAS_top100_A2_fraction = sig$raw_top100,
                                    M2_initial_top100_A2_fraction = sig$m2_top100,
                                    raw_TWAS_top10pct_A2_fraction = sig$raw_top10pct,
                                    M2_initial_top10pct_A2_fraction = sig$m2_top10pct,
                                    stringsAsFactors = FALSE); is <- is + 1
    leakage_rows[[il]] <- branch_leakage_row(rep); il <- il + 1
    scores <- benchmarks_for_rep(rep, ch, no_diff)
    for (nm in names(scores)) {
      primary_rows[[ip]] <- non_directional_metrics(rep, scores[[nm]]$score, nm); ip <- ip + 1
      direction_rows[[id]] <- directional_metrics(rep, scores[[nm]]$signed, nm); id <- id + 1
      rank_rows[[ir]] <- a2_rank_row(rep, scores[[nm]]$score, nm); ir <- ir + 1
      degree_rows[[ideg]] <- score_degree_row(rep, scores[[nm]]$score, nm); ideg <- ideg + 1
    }
    polarity_rows[[ipo]] <- channel_polarity_row(rep, ch); ipo <- ipo + 1
    terminal_rows[[it]] <- terminal_capture_row(rep, ch); it <- it + 1
    if (rep_id %% 5 == 0) gc(FALSE)
  }

  primary <- do.call(rbind, primary_rows)
  direction <- do.call(rbind, direction_rows)
  ranks <- do.call(rbind, rank_rows)
  degree <- do.call(rbind, degree_rows)
  polarity <- do.call(rbind, polarity_rows)
  terminal <- do.call(rbind, terminal_rows)
  leakage <- do.call(rbind, leakage_rows)
  signal <- do.call(rbind, signal_rows)
  summary <- summarise_branch_metrics(primary, direction, terminal, ranks, leakage, degree)
  contrasts <- make_contrasts(primary, direction, ranks)
  schema <- metric_schema_audit(primary, direction)

  write_csv_over(summary, file.path(report_dir, "DENSE_BRANCH_ISOLATION_ARM_METRICS.csv"))
  write_csv_over(contrasts, file.path(report_dir, "DENSE_BRANCH_ISOLATION_ARM_CONTRASTS.csv"))
  write_csv_over(primary, file.path(report_dir, "PRIMARY_FINAL_HEAT_METRICS.csv"))
  write_csv_over(direction, file.path(report_dir, "DIRECTION_AWARE_METRICS.csv"))
  write_csv_over(primary, file.path(report_dir, "BENCHMARK_METRICS.csv"))
  write_csv_over(contrasts, file.path(report_dir, "BENCHMARK_CONTRASTS.csv"))
  write_csv_over(terminal, file.path(report_dir, "TERMINAL_CAPTURE_AUDIT.csv"))
  write_csv_over(terminal[, c("arm", "arm_label", "replicate", "relay_heat_mass", "A2_heat_mass",
                              "relay_to_A2_heat_ratio", "relay_top100_rate", "A2_top100_rate",
                              "median_A2_rank", "median_relay_rank")],
                 file.path(report_dir, "RELAY_VS_A2_HEAT_AUDIT.csv"))
  write_csv_over(ranks, file.path(report_dir, "A2_RANK_DISTRIBUTION_AUDIT.csv"))
  write_csv_over(leakage, file.path(report_dir, "BRANCH_LEAKAGE_AUDIT.csv"))
  write_csv_over(polarity, file.path(report_dir, "NESTA_CHANNEL_POLARITY_AUDIT.csv"))
  write_csv_over(degree, file.path(report_dir, "DEGREE_BIAS_METRICS.csv"))
  write_csv_over(schema, file.path(report_dir, "METRIC_SCHEMA_AUDIT.csv"))
  write_csv_over(data.frame(guardrail = c("raw_m2_top100_A2_leakage", "opposite_sign_decoy", "high_degree_decoy"),
                            passed = c(all(signal$raw_TWAS_top100_A2_fraction <= 0.05 &
                                             signal$M2_initial_top100_A2_fraction <= 0.05),
                                       all(summary$opposite_sign_decoy_top100_rate <= 0.15),
                                       all(summary$high_degree_decoy_top100_rate <= 0.12)),
                            stringsAsFactors = FALSE),
                 file.path(report_dir, "NULL_BIAS_GUARDRAILS.csv"))

  if (!any(summary$mean_top200_recall > 0.05)) {
    reason <- "branch_isolation_rescue_insufficient"
  } else if (any(summary$heat_reaching_A2 >= 0.020 & summary$median_A2_rank >= 500)) {
    reason <- "A2_heat_nonzero_but_rank_compression_failure"
  } else {
    reason <- "dense_branch_isolation_diagnostic_completed_no_confirmatory"
  }
  set_run_status("STOPPED", "YES", "NO", reason)
  best <- summary[order(summary$mean_top200_recall, summary$mean_top100_recall,
                        -summary$median_A2_rank, decreasing = TRUE), ][1, ]

  write_lines_over(c("# Dense Branch-Isolation Arm Summary", "",
                     paste0("Best arm by top200/top100/rank ordering: `", best$arm, "` (", best$arm_label, ")."),
                     sprintf("Top100/top150/top200 recall: %.4f / %.4f / %.4f.",
                             best$mean_top100_recall, best$mean_top150_recall, best$mean_top200_recall),
                     sprintf("Median A2 rank: %.2f; first A2 rank mean: %.2f.",
                             best$median_A2_rank, best$first_A2_rank),
                     sprintf("Retained in A / reaching A2 / background leakage: %.4f / %.4f / %.4f.",
                             best$retained_in_A_branch, best$heat_reaching_A2, best$background_leakage)),
                   file.path(report_dir, "DENSE_BRANCH_ISOLATION_ARM_SUMMARY.md"))
  write_lines_over(c("# Terminal Capture Audit", "",
                     sprintf("Mean seed heat retained in A branch: %.4f.", mean(terminal$fraction_seed_heat_retained_in_A_branch)),
                     sprintf("Mean heat reaching relay/A2: %.4f / %.4f.",
                             mean(terminal$fraction_seed_heat_reaching_relay), mean(terminal$fraction_seed_heat_reaching_A2)),
                     sprintf("Mean background leakage: %.4f.", mean(terminal$fraction_seed_heat_leaking_to_background)),
                     sprintf("Mean relay-to-A2 heat ratio: %.4f.", mean(terminal$relay_to_A2_heat_ratio))),
                   file.path(report_dir, "TERMINAL_CAPTURE_AUDIT.md"))
  write_lines_over(c("# Relay Vs A2 Heat Audit", "",
                     "Relay and A2 heat masses, top-k rates, and median ranks are reported in RELAY_VS_A2_HEAT_AUDIT.csv."),
                   file.path(report_dir, "RELAY_VS_A2_HEAT_AUDIT.md"))
  write_lines_over(c("# A2 Rank Distribution Audit", "",
                     sprintf("Best-arm mean median A2 rank: %.2f.", best$median_A2_rank),
                     "Full per-replicate rank distributions are in A2_RANK_DISTRIBUTION_AUDIT.csv."),
                   file.path(report_dir, "A2_RANK_DISTRIBUTION_AUDIT.md"))
  write_lines_over(c("# Branch Leakage Audit", "",
                     sprintf("Mean seed-neighborhood background fraction: %.5f.", mean(leakage$seed_neighborhood_background_fraction)),
                     sprintf("Mean A-branch background cut fraction: %.5f.", mean(leakage$A_branch_background_cut_fraction))),
                   file.path(report_dir, "BRANCH_LEAKAGE_AUDIT.md"))
  write_lines_over(c("# NESTA Channel Polarity Audit", "",
                     sprintf("Risk A2 positive signed heat fraction: %.4f.", mean(polarity$fraction_risk_A2_positive_signed_heat)),
                     sprintf("Protective A2 negative signed heat fraction: %.4f.", mean(polarity$fraction_protective_A2_negative_signed_heat)),
                     sprintf("A2 nonzero final heat fraction: %.4f.", mean(polarity$fraction_A2_nonzero_final_heat))),
                   file.path(report_dir, "NESTA_CHANNEL_POLARITY_AUDIT.md"))
  write_lines_over(c("# Metric Schema Audit", "",
                     sprintf("PRIMARY_FINAL_HEAT_METRICS.csv schema pass: %s.",
                             schema$required_columns_present[schema$file == "PRIMARY_FINAL_HEAT_METRICS.csv"]),
                     sprintf("DIRECTION_AWARE_METRICS.csv schema pass: %s.",
                             schema$required_columns_present[schema$file == "DIRECTION_AWARE_METRICS.csv"])),
                   file.path(report_dir, "METRIC_SCHEMA_AUDIT.md"))
  write_lines_over(c("# Code Fidelity Audit", "",
                     paste0("Binding plan SHA256: `", binding_plan_sha256, "`."),
                     "Faithful M2 arithmetic, P conversion, signed positive/absolute-negative diffusion, retained zero-weight edges, self-loop behavior, and `n.perm` are used."),
                   file.path(report_dir, "CODE_FIDELITY_AUDIT.md"))
  write_lines_over(c("# Benchmark Implementation Audit", "",
                     "Comparators evaluated: raw_TWAS_abs, raw_signed_TWAS, M2_no_diffusion, RWR_abs_prior, PPR_abs_prior, RWR_signed_two_channel, PPR_signed_two_channel, median_weight_permutation, I2_module_disrupted, and I3_expression_matched_randomized.",
                     "PPR_abs_prior and RWR_abs_prior are direction-blind sensitivity comparators only."),
                   file.path(report_dir, "BENCHMARK_IMPLEMENTATION_AUDIT.md"))
  write_lines_over(c("# STOP/GO Report", "", "STOP.", paste0("Reason: `", reason, "`."),
                     "Pilot started: YES.", "Confirmatory started: NO.",
                     "This binding plan is diagnostic-only and prohibits confirmatory execution."),
                   file.path(report_dir, "STOP_GO_REPORT.md"))
  write_lines_over(c("# Final Report", "",
                     paste0("Binding plan SHA256: `", binding_plan_sha256, "`."),
                     "Dense-only branch-isolation relay-sink rescue completed with 20 replicates per arm across Arms E-H.",
                     "Pilot execution started: YES.",
                     "Confirmatory execution started: NO.",
                     paste0("STOP classification: `", reason, "`."),
                     sprintf("Best arm `%s`: top100/top150/top200 %.4f / %.4f / %.4f; median A2 rank %.2f.",
                             best$arm, best$mean_top100_recall, best$mean_top150_recall,
                             best$mean_top200_recall, best$median_A2_rank),
                     sprintf("Best-arm retained A branch %.4f, A2 heat %.4f, background leakage %.4f, relay-to-A2 heat ratio %.4f.",
                             best$retained_in_A_branch, best$heat_reaching_A2,
                             best$background_leakage, best$relay_to_A2_heat_ratio)),
                   file.path(report_dir, "FINAL_REPORT.md"))

  if (file.exists(project_file("results/reports/summary_tables/unit_test_results.csv"))) {
    file.copy(project_file("results/reports/summary_tables/unit_test_results.csv"),
              file.path(report_dir, "unit_test_results.csv"), overwrite = TRUE)
  }
  if (file.exists(project_file("CLEANUP_MANIFEST.csv"))) {
    file.copy(project_file("CLEANUP_MANIFEST.csv"), file.path(report_dir, "CLEANUP_MANIFEST.csv"),
              overwrite = TRUE)
  }
  write_csv_over(data.frame(path = c(project_file("R/study_0707_branch_isolation.R"),
                                     project_file("R/study_0707_dense_rescue.R"),
                                     project_file("R/fidelity.R"), project_file("R/utils.R"),
                                     project_file("scripts/run_branch_isolation_rescue.R")),
                            sha256 = c(sha(project_file("R/study_0707_branch_isolation.R")),
                                       sha(project_file("R/study_0707_dense_rescue.R")),
                                       sha(project_file("R/fidelity.R")),
                                       sha(project_file("R/utils.R")),
                                       if (file.exists(project_file("scripts/run_branch_isolation_rescue.R"))) sha(project_file("scripts/run_branch_isolation_rescue.R")) else NA_character_),
                            role = c("branch_isolation_runner", "dense_rescue_shared_helpers",
                                     "faithful_nesta_and_benchmarks", "binding_plan_and_io_guards",
                                     "script_entrypoint"),
                            stringsAsFactors = FALSE),
                 file.path(report_dir, "IMPLEMENTATION_MANIFEST.csv"))
  files <- list.files(report_dir, recursive = TRUE, full.names = TRUE)
  files <- files[file.info(files)$isdir == FALSE & basename(files) != "CHECKSUMS.sha256"]
  writeLines(vapply(files, function(f) paste(sha(f), sub(paste0("^", report_dir, "/"), "", f)), character(1)),
             file.path(report_dir, "CHECKSUMS.sha256"))
  invisible(list(reason = reason, summary = summary))
}

if (!identical(Sys.getenv("NESTA_BRANCH_RESCUE_SOURCE_ONLY"), "1")) run_branch_isolation_rescue()
