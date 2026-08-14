Sys.setenv(NESTA_CONFIRMATORY_SOURCE_ONLY = "1")
source(file.path("/home/js/NESTA/simulation", "R", "study_0708_comparator_confirmatory.R"))

strict_fractions <- c(`top_1pct` = 0.01, `top_2_5pct` = 0.025, `top_5pct` = 0.05, `top_10pct` = 0.10)
strict_nesta_modes <- c("NESTA_signed_descending", "NESTA_signed_ascending", "NESTA_abs_final_heat", "NESTA_two_tail_balanced", "NESTA_two_tail_direction_matched")
strict_primary_comparators <- c("raw_TWAS_abs", "M2_no_diffusion_abs", "PPR_abs_prior", "RWR_abs_prior", "I2_module_disrupted", "I3_expression_matched_randomized", "median_weight_permutation")
strict_upper_comparators <- c("PPR_signed_two_channel", "RWR_signed_two_channel")

strict_ranked <- function(rep, signed_score, mode) {
  genes <- rep$genes
  score <- setNames(as.numeric(signed_score[genes]), genes)
  if (mode == "NESTA_signed_descending") return(list(ranked = ranked_genes(genes, score, TRUE, rep$A1), score_for_degree = score, pos = NULL, neg = NULL))
  if (mode == "NESTA_signed_ascending") return(list(ranked = ranked_genes(genes, score, FALSE, rep$A1), score_for_degree = -score, pos = NULL, neg = NULL))
  if (mode == "NESTA_abs_final_heat") return(list(ranked = ranked_genes(genes, abs(score), TRUE, rep$A1), score_for_degree = abs(score), pos = NULL, neg = NULL))
  pos <- ranked_genes(genes, score, TRUE, rep$A1)
  neg <- ranked_genes(genes, score, FALSE, rep$A1)
  list(ranked = NULL, score_for_degree = abs(score), pos = pos, neg = neg)
}

strict_selected <- function(rep, signed_score, mode, k) {
  r <- strict_ranked(rep, signed_score, mode)
  if (!is.null(r$ranked)) return(list(selected = head(r$ranked, k), pos_tail = character(), neg_tail = character(), score_for_degree = r$score_for_degree))
  pos_n <- ceiling(k / 2); neg_n <- floor(k / 2)
  pos <- head(r$pos, pos_n); neg <- head(r$neg, neg_n)
  list(selected = unique(c(pos, neg)), pos_tail = pos, neg_tail = neg, score_for_degree = r$score_for_degree)
}

strict_metric_row <- function(rep, signed_score, mode, cutoff_fraction, score_name = mode, score_class = "NESTA_final_heat") {
  genes <- rep$genes
  ranked_universe_n <- length(genes)
  k <- ceiling(cutoff_fraction * ranked_universe_n)
  sel <- strict_selected(rep, signed_score, mode, k)
  selected <- sel$selected
  pos_tail <- sel$pos_tail
  neg_tail <- sel$neg_tail
  if (!length(pos_tail) && !length(neg_tail)) {
    if (mode == "NESTA_signed_descending") pos_tail <- selected
    else if (mode == "NESTA_signed_ascending") neg_tail <- selected
  }
  total_hits <- sum(rep$A2 %in% selected)
  risk_hits <- if (mode == "NESTA_two_tail_direction_matched") sum(rep$A2_risk %in% pos_tail) else sum(rep$A2_risk %in% selected)
  prot_hits <- if (mode == "NESTA_two_tail_direction_matched") sum(rep$A2_protective %in% neg_tail) else sum(rep$A2_protective %in% selected)
  risk_den <- if (length(pos_tail)) length(pos_tail) else k
  prot_den <- if (length(neg_tail)) length(neg_tail) else k
  deg <- score_degree(rep$adj)[genes]
  sc <- sel$score_for_degree[genes]
  au_score <- if (mode %in% c("NESTA_signed_descending")) signed_score[genes] else if (mode == "NESTA_signed_ascending") -signed_score[genes] else abs(signed_score[genes])
  au <- auprc_from_score(genes, au_score, rep$A2, rep$A1)
  random_expected <- k / max(length(setdiff(genes, rep$A1)), 1)
  data.frame(topology_arm = rep$arm, topology_label = rep$arm_label, rescue_arm = rep$bidirectional_arm,
             replicate = rep$rep_id, score_name = score_name, ranking_mode = mode, score_class = score_class,
             cutoff_fraction = cutoff_fraction, cutoff_label = names(strict_fractions)[match(cutoff_fraction, strict_fractions)], cutoff_k = k,
             total_A2_hits = total_hits, total_A2_recall = total_hits / length(rep$A2), total_A2_precision = total_hits / k,
             enrichment_over_random = (total_hits / length(rep$A2)) / random_expected,
             risk_A2_hits = risk_hits, risk_A2_recall = risk_hits / length(rep$A2_risk), risk_A2_precision = risk_hits / risk_den,
             protective_A2_hits = prot_hits, protective_A2_recall = prot_hits / length(rep$A2_protective), protective_A2_precision = prot_hits / prot_den,
             opposite_direction_decoy_hits = sum(rep$D %in% selected), opposite_direction_decoy_rate = mean(rep$D %in% selected),
             high_degree_decoy_hits = sum(rep$C %in% selected), high_degree_decoy_rate = mean(rep$C %in% selected),
             score_degree_spearman = suppressWarnings(stats::cor(sc, deg, method = "spearman", use = "pairwise.complete.obs")),
             full_ranking_AUPRC = au, ranked_universe_n = ranked_universe_n, stringsAsFactors = FALSE)
}

strict_comparator_row <- function(rep, score, signed_score, score_name, score_class, cutoff_fraction) {
  genes <- rep$genes
  score <- setNames(as.numeric(score[genes]), genes)
  signed_score <- setNames(as.numeric(signed_score[genes]), genes)
  k <- ceiling(cutoff_fraction * length(genes))
  ranked <- ranked_genes(genes, score, TRUE, rep$A1)
  selected <- head(ranked, k)
  total_hits <- sum(rep$A2 %in% selected)
  risk_hits <- sum(rep$A2_risk %in% selected)
  prot_hits <- sum(rep$A2_protective %in% selected)
  deg <- score_degree(rep$adj)[genes]
  au <- auprc_from_score(genes, score, rep$A2, rep$A1)
  random_expected <- k / max(length(setdiff(genes, rep$A1)), 1)
  data.frame(topology_arm = rep$arm, topology_label = rep$arm_label, rescue_arm = rep$bidirectional_arm,
             replicate = rep$rep_id, score_name = score_name, score_class = score_class,
             cutoff_fraction = cutoff_fraction, cutoff_label = names(strict_fractions)[match(cutoff_fraction, strict_fractions)], cutoff_k = k,
             total_A2_hits = total_hits, total_A2_recall = total_hits / length(rep$A2), total_A2_precision = total_hits / k,
             enrichment_over_random = (total_hits / length(rep$A2)) / random_expected,
             risk_A2_hits = risk_hits, risk_A2_recall = risk_hits / length(rep$A2_risk), risk_A2_precision = risk_hits / k,
             protective_A2_hits = prot_hits, protective_A2_recall = prot_hits / length(rep$A2_protective), protective_A2_precision = prot_hits / k,
             opposite_direction_decoy_hits = sum(rep$D %in% selected), opposite_direction_decoy_rate = mean(rep$D %in% selected),
             high_degree_decoy_hits = sum(rep$C %in% selected), high_degree_decoy_rate = mean(rep$C %in% selected),
             score_degree_spearman = suppressWarnings(stats::cor(score, deg, method = "spearman", use = "pairwise.complete.obs")),
             full_ranking_AUPRC = au, ranked_universe_n = length(genes), stringsAsFactors = FALSE)
}

strict_summary <- function(x, group_cols) {
  stats::aggregate(cbind(cutoff_k, total_A2_hits, total_A2_recall, total_A2_precision, enrichment_over_random,
                         risk_A2_hits, risk_A2_recall, risk_A2_precision, protective_A2_hits, protective_A2_recall,
                         protective_A2_precision, opposite_direction_decoy_hits, opposite_direction_decoy_rate,
                         high_degree_decoy_hits, high_degree_decoy_rate, score_degree_spearman, full_ranking_AUPRC) ~ .,
                   data = x[, c(group_cols, "cutoff_k", "total_A2_hits", "total_A2_recall", "total_A2_precision", "enrichment_over_random",
                                "risk_A2_hits", "risk_A2_recall", "risk_A2_precision", "protective_A2_hits", "protective_A2_recall",
                                "protective_A2_precision", "opposite_direction_decoy_hits", "opposite_direction_decoy_rate",
                                "high_degree_decoy_hits", "high_degree_decoy_rate", "score_degree_spearman", "full_ranking_AUPRC"), drop = FALSE],
                   FUN = function(z) mean(z, na.rm = TRUE))
}

strict_contrasts <- function(nesta, comps) {
  out <- list(); i <- 1L
  n0 <- nesta[nesta$ranking_mode %in% c("NESTA_two_tail_balanced", "NESTA_two_tail_direction_matched"), , drop = FALSE]
  for (tp in unique(n0$topology_arm)) for (mode in unique(n0$ranking_mode)) for (frac in strict_fractions) {
    base <- n0[n0$topology_arm == tp & n0$ranking_mode == mode & n0$cutoff_fraction == frac, c("replicate", "total_A2_recall", "total_A2_precision"), drop = FALSE]
    for (cmp in unique(comps$score_name)) {
      other <- comps[comps$topology_arm == tp & comps$score_name == cmp & comps$cutoff_fraction == frac, c("replicate", "total_A2_recall", "total_A2_precision"), drop = FALSE]
      z <- merge(base, other, by = "replicate", suffixes = c("_nesta", "_comparator"))
      for (metric in c("total_A2_recall", "total_A2_precision")) {
        d <- z[[paste0(metric, "_nesta")]] - z[[paste0(metric, "_comparator")]]
        ci <- ci_mean(d, seed = 9709L + i)
        out[[i]] <- data.frame(topology_arm = tp, ranking_mode = mode, comparator = cmp, cutoff_fraction = frac,
                               cutoff_label = names(strict_fractions)[match(frac, strict_fractions)], metric = metric,
                               nesta_mean = mean(z[[paste0(metric, "_nesta")]], na.rm = TRUE),
                               comparator_mean = mean(z[[paste0(metric, "_comparator")]], na.rm = TRUE),
                               paired_mean_difference = ci["mean"], paired_ci_low = ci["ci_low"], paired_ci_high = ci["ci_high"],
                               stringsAsFactors = FALSE)
        i <- i + 1L
      }
    }
  }
  do.call(rbind, out)
}

run_strict_replicate <- function(seed_row) {
  tp <- seed_row$topology_arm
  base <- make_branch_isolation_rep(tp, seed_row$replicate, seed_row$base_seed)
  rep <- apply_bidirectional_arm(base, confirmatory_arm, seed_row$signal_seed)
  ch <- faithful_m2_channels(rep$adj, rep$expr, rep$twas, n.perm = 25, seed = seed_row$nesta_seed)
  no <- no_diffusion_m2_scores(rep$adj, rep$expr, rep$twas)
  signed <- setNames(ch$final_NESTA_heat, ch$SYMBOL)[rep$genes]
  nesta_rows <- list(); i <- 1L
  for (mode in strict_nesta_modes) for (frac in strict_fractions) {
    nesta_rows[[i]] <- strict_metric_row(rep, signed, mode, frac); i <- i + 1L
  }
  comps <- confirmatory_comparators(rep, ch, no, seed_row)
  comp_rows <- list(); upper_rows <- list(); ic <- iu <- 1L
  for (nm in strict_primary_comparators) for (frac in strict_fractions) {
    comp_rows[[ic]] <- strict_comparator_row(rep, comps[[nm]]$score, comps[[nm]]$signed, nm, comps[[nm]]$class, frac); ic <- ic + 1L
  }
  for (nm in strict_upper_comparators) for (frac in strict_fractions) {
    upper_rows[[iu]] <- strict_comparator_row(rep, comps[[nm]]$score, comps[[nm]]$signed, nm, comps[[nm]]$class, frac); iu <- iu + 1L
  }
  list(nesta = do.call(rbind, nesta_rows), comps = do.call(rbind, comp_rows), upper = do.call(rbind, upper_rows))
}

run_strict_top_fraction_audit <- function() {
  verify_project_path(); verify_binding_plan()
  report_dir <- read_report_dir(); safe_dir_create(report_dir); copy_binding_plan_to_report(report_dir)
  set_run_status("IN_PROGRESS", "YES", "NO", "0709 strict top-fraction audit")
  previous_dir <- "/home/js_subdir/Dropbox/NESTA_revision/NESTA_simulation_080726_1951"
  posthoc_available <- any(grepl("gene|rank|score", list.files(previous_dir), ignore.case = TRUE)) && any(grepl("per", list.files(previous_dir), ignore.case = TRUE))
  schedule <- make_seed_schedule()
  nesta_rows <- list(); comp_rows <- list(); upper_rows <- list(); idx <- 1L
  for (r in seq_len(nrow(schedule))) {
    res <- run_strict_replicate(schedule[r, ])
    nesta_rows[[r]] <- res$nesta; comp_rows[[r]] <- res$comps; upper_rows[[r]] <- res$upper
    if (r %% 10 == 0) { message(sprintf("strict audit replicate %d/%d complete", r, nrow(schedule))); gc(FALSE) }
  }
  nesta <- do.call(rbind, nesta_rows)
  comps <- do.call(rbind, comp_rows)
  upper <- do.call(rbind, upper_rows)
  nesta_sum <- strict_summary(nesta, c("topology_arm", "topology_label", "rescue_arm", "ranking_mode", "score_name", "score_class", "cutoff_fraction", "cutoff_label"))
  comp_sum <- strict_summary(comps, c("topology_arm", "topology_label", "rescue_arm", "score_name", "score_class", "cutoff_fraction", "cutoff_label"))
  upper_sum <- strict_summary(upper, c("topology_arm", "topology_label", "rescue_arm", "score_name", "score_class", "cutoff_fraction", "cutoff_label"))
  contrasts <- strict_contrasts(nesta, comps)
  dir_sum <- nesta_sum[nesta_sum$ranking_mode %in% c("NESTA_signed_descending", "NESTA_signed_ascending", "NESTA_two_tail_direction_matched"), , drop = FALSE]
  decoys <- rbind(nesta_sum[, c("topology_arm", "rescue_arm", "ranking_mode", "score_name", "cutoff_fraction", "cutoff_label", "opposite_direction_decoy_hits", "opposite_direction_decoy_rate", "high_degree_decoy_hits", "high_degree_decoy_rate"), drop = FALSE],
                  transform(comp_sum[, c("topology_arm", "rescue_arm", "score_name", "cutoff_fraction", "cutoff_label", "opposite_direction_decoy_hits", "opposite_direction_decoy_rate", "high_degree_decoy_hits", "high_degree_decoy_rate"), drop = FALSE], ranking_mode = score_name))
  f5 <- nesta_sum[nesta_sum$topology_arm == "F" & nesta_sum$ranking_mode == "NESTA_two_tail_direction_matched" & abs(nesta_sum$cutoff_fraction - 0.05) < 1e-12, , drop = FALSE]
  p5 <- comp_sum[comp_sum$topology_arm == "F" & comp_sum$score_name %in% c("PPR_abs_prior", "RWR_abs_prior") & abs(comp_sum$cutoff_fraction - 0.05) < 1e-12, , drop = FALSE]
  p1 <- comp_sum[comp_sum$topology_arm == "F" & comp_sum$score_name %in% c("PPR_abs_prior", "RWR_abs_prior") & abs(comp_sum$cutoff_fraction - 0.01) < 1e-12, , drop = FALSE]
  if (nrow(f5) && f5$risk_A2_recall >= 0.60 && f5$protective_A2_recall >= 0.60 && max(p5$total_A2_recall, na.rm = TRUE) >= f5$total_A2_recall) {
    label <- "strict_fraction_audit_completed_ppr_rwr_remain_strong"
  } else if (nrow(f5) && f5$risk_A2_recall >= 0.60 && f5$protective_A2_recall >= 0.60) {
    label <- "strict_fraction_audit_completed_nesta_directional_signal_retained"
  } else if (nrow(p1) && max(p1$total_A2_recall, na.rm = TRUE) < 0.50 && nrow(p5) && max(p5$total_A2_recall, na.rm = TRUE) < 0.80) {
    label <- "strict_fraction_audit_reveals_top100_leniency"
  } else {
    label <- "strict_fraction_audit_completed_nesta_directional_signal_retained"
  }
  write_csv_over(nesta_sum, file.path(report_dir, "STRICT_TOP_FRACTION_METRICS.csv"))
  write_csv_over(dir_sum, file.path(report_dir, "STRICT_DIRECTION_SPECIFIC_METRICS.csv"))
  write_csv_over(comp_sum, file.path(report_dir, "STRICT_COMPARATOR_METRICS.csv"))
  write_csv_over(upper_sum, file.path(report_dir, "STRICT_SIGNED_UPPER_BOUND_METRICS.csv"))
  write_csv_over(decoys, file.path(report_dir, "STRICT_DECOY_GUARDRAILS.csv"))
  write_csv_over(contrasts, file.path(report_dir, "STRICT_FRACTION_CONTRASTS.csv"))
  write_csv_over(nesta_sum[, c("topology_arm", "ranking_mode", "cutoff_fraction", "cutoff_label", "total_A2_recall", "risk_A2_recall", "protective_A2_recall")], file.path(report_dir, "figure_top_fraction_recall.csv"))
  write_csv_over(nesta_sum[, c("topology_arm", "ranking_mode", "cutoff_fraction", "cutoff_label", "total_A2_precision", "risk_A2_precision", "protective_A2_precision")], file.path(report_dir, "figure_top_fraction_precision.csv"))
  write_csv_over(dir_sum, file.path(report_dir, "figure_direction_specific_top_fraction.csv"))
  write_csv_over(nesta_sum, file.path(report_dir, "table_strict_shortlist_results.csv"))
  write_csv_over(comp_sum[comp_sum$score_name %in% c("PPR_abs_prior", "RWR_abs_prior"), ], file.path(report_dir, "table_ppr_rwr_strict_audit.csv"))
  schema <- data.frame(file = c("STRICT_TOP_FRACTION_METRICS.csv", "STRICT_COMPARATOR_METRICS.csv", "STRICT_DIRECTION_SPECIFIC_METRICS.csv"),
                       required_columns_present = c(all(c("cutoff_fraction", "cutoff_k", "total_A2_recall", "total_A2_precision", "risk_A2_recall", "protective_A2_recall") %in% names(nesta_sum)),
                                                    all(c("score_name", "cutoff_fraction", "total_A2_recall", "total_A2_precision") %in% names(comp_sum)),
                                                    all(c("risk_A2_precision", "protective_A2_precision") %in% names(dir_sum))), stringsAsFactors = FALSE)
  write_csv_over(schema, file.path(report_dir, "METRIC_SCHEMA_AUDIT.csv"))
  bal5 <- nesta_sum[nesta_sum$topology_arm == "F" & nesta_sum$ranking_mode == "NESTA_two_tail_balanced" & abs(nesta_sum$cutoff_fraction - 0.05) < 1e-12, , drop = FALSE]
  ppr5 <- comp_sum[comp_sum$topology_arm == "F" & comp_sum$score_name == "PPR_abs_prior" & abs(comp_sum$cutoff_fraction - 0.05) < 1e-12, , drop = FALSE]
  rwr5 <- comp_sum[comp_sum$topology_arm == "F" & comp_sum$score_name == "RWR_abs_prior" & abs(comp_sum$cutoff_fraction - 0.05) < 1e-12, , drop = FALSE]
  write_lines_over(c("# Strict Top-Fraction Summary", "", paste0("Classification: `", label, "`."),
                     sprintf("Post-hoc per-gene rankings available from 080726_1951: %s; frozen-design rerun performed: %s.", posthoc_available, !posthoc_available),
                     sprintf("F NESTA two-tail top5%% recall/precision: %.4f / %.4f.", bal5$total_A2_recall, bal5$total_A2_precision),
                     sprintf("F NESTA direction-matched top5%% risk/protective recall: %.4f / %.4f.", f5$risk_A2_recall, f5$protective_A2_recall),
                     sprintf("F PPR_abs/RWR_abs top5%% recall: %.4f / %.4f.", ppr5$total_A2_recall, rwr5$total_A2_recall)), file.path(report_dir, "STRICT_TOP_FRACTION_SUMMARY.md"))
  write_lines_over(c("# Manuscript Shortlist Recommendation", "", "Use strict top-fraction metrics alongside fixed top100 metrics.",
                     "For signed NESTA, emphasize direction-specific tail recovery when total two-tail shortlists are split across risk and protective directions.",
                     "PPR/RWR absolute-prior results should remain framed as unsigned disease-relevance propagation comparators."), file.path(report_dir, "MANUSCRIPT_SHORTLIST_RECOMMENDATION.md"))
  write_lines_over(c("# Code Fidelity Audit", "", paste0("Binding plan SHA256: `", binding_plan_sha256, "`."),
                     "The frozen successful design was rerun only because prior per-gene ranks were unavailable.",
                     "Faithful submitted M2 arithmetic, TWAS.P conversion, strict TWAS filtering, signed positive/absolute-negative diffusion, zero-weight edges, self-loops, and diffuStats `n.perm` were retained.",
                     "Submitted Final Heat values were not modified."), file.path(report_dir, "CODE_FIDELITY_AUDIT.md"))
  write_lines_over(c("# Metric Schema Audit", "", sprintf("All schema checks passed: %s.", all(schema$required_columns_present))), file.path(report_dir, "METRIC_SCHEMA_AUDIT.md"))
  write_lines_over(c("# Final Report", "", paste0("Binding plan SHA256: `", binding_plan_sha256, "`."), paste0("Final classification: `", label, "`."),
                     sprintf("F NESTA two-tail-balanced top5%% recall/precision %.4f / %.4f.", bal5$total_A2_recall, bal5$total_A2_precision),
                     sprintf("F NESTA direction-matched top5%% risk/protective recall %.4f / %.4f.", f5$risk_A2_recall, f5$protective_A2_recall),
                     sprintf("F PPR_abs_prior top1%%/top5%% recall %.4f / %.4f.", comp_sum$total_A2_recall[comp_sum$topology_arm == "F" & comp_sum$score_name == "PPR_abs_prior" & abs(comp_sum$cutoff_fraction - 0.01) < 1e-12], ppr5$total_A2_recall),
                     sprintf("F RWR_abs_prior top1%%/top5%% recall %.4f / %.4f.", comp_sum$total_A2_recall[comp_sum$topology_arm == "F" & comp_sum$score_name == "RWR_abs_prior" & abs(comp_sum$cutoff_fraction - 0.01) < 1e-12], rwr5$total_A2_recall),
                     "This was an audit round; no adaptive simulation changes were made after strict cutoff results."), file.path(report_dir, "FINAL_REPORT.md"))
  if (file.exists(project_file("CLEANUP_MANIFEST.csv"))) file.copy(project_file("CLEANUP_MANIFEST.csv"), file.path(report_dir, "CLEANUP_MANIFEST.csv"), overwrite = TRUE)
  if (file.exists(project_file("results/reports/summary_tables/unit_test_results.csv"))) file.copy(project_file("results/reports/summary_tables/unit_test_results.csv"), file.path(report_dir, "unit_test_results.csv"), overwrite = TRUE)
  write_csv_over(data.frame(path = c(project_file("R/study_0709_strict_top_fraction_audit.R"), project_file("R/study_0708_comparator_confirmatory.R"), project_file("R/study_0708_bidirectional_rescue.R"), project_file("R/fidelity.R"), project_file("R/utils.R"), project_file("scripts/run_strict_top_fraction_audit.R")),
                            sha256 = c(sha(project_file("R/study_0709_strict_top_fraction_audit.R")), sha(project_file("R/study_0708_comparator_confirmatory.R")), sha(project_file("R/study_0708_bidirectional_rescue.R")), sha(project_file("R/fidelity.R")), sha(project_file("R/utils.R")), sha(project_file("scripts/run_strict_top_fraction_audit.R"))),
                            role = c("strict_top_fraction_audit_runner", "frozen_confirmatory_helpers", "bidirectional_ranking_helpers", "faithful_nesta_and_benchmarks", "binding_plan_and_io_guards", "script_entrypoint"), stringsAsFactors = FALSE), file.path(report_dir, "IMPLEMENTATION_MANIFEST.csv"))
  files <- list.files(report_dir, recursive = TRUE, full.names = TRUE)
  files <- files[file.info(files)$isdir == FALSE & basename(files) != "CHECKSUMS.sha256"]
  writeLines(vapply(files, function(f) paste(sha(f), sub(paste0("^", report_dir, "/"), "", f)), character(1)), file.path(report_dir, "CHECKSUMS.sha256"))
  set_run_status("COMPLETE", "YES", "NO", label)
  invisible(list(label = label, nesta = nesta_sum, comps = comp_sum, upper = upper_sum))
}

if (!identical(Sys.getenv("NESTA_STRICT_FRACTION_SOURCE_ONLY"), "1")) run_strict_top_fraction_audit()
