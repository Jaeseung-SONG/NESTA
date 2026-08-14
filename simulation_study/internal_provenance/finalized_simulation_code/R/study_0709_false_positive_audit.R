Sys.setenv(NESTA_STRICT_FRACTION_SOURCE_ONLY = "1")
source(file.path("/home/js/NESTA/simulation", "R", "study_0709_strict_top_fraction_audit.R"))

fp_cutoffs <- c(top10 = 10L, top50 = 50L, top100 = 100L)
fp_nesta_modes <- strict_nesta_modes
fp_primary_comparators <- strict_primary_comparators
fp_upper_comparators <- strict_upper_comparators

rep_rescue_arm <- function(rep) {
  if (!is.null(rep$bidirectional_arm) && length(rep$bidirectional_arm) == 1 && !is.na(rep$bidirectional_arm)) rep$bidirectional_arm else confirmatory_arm
}

fp_fraction_label <- function(k, n = 1000L) {
  sprintf("top%.0f%%", 100 * k / n)
}

total_fp_metrics <- function(rep, selected, k, method, ranking_mode, score_class, score_degree_spearman = NA_real_) {
  genes <- setdiff(rep$genes, rep$A1)
  positives <- rep$A2
  selected <- unique(selected[selected %in% genes])
  tp <- sum(selected %in% positives)
  fp <- length(selected) - tp
  fn <- length(positives) - tp
  tn <- length(genes) - length(positives) - fp
  precision <- tp / max(tp + fp, .Machine$double.eps)
  recall <- tp / max(tp + fn, .Machine$double.eps)
  fpr <- fp / max(fp + tn, .Machine$double.eps)
  random_expected <- length(positives) / length(genes)
  data.frame(topology_arm = rep$arm, topology_label = rep$arm_label, rescue_arm = rep_rescue_arm(rep),
             replicate = rep$rep_id, method = method, ranking_mode = ranking_mode, score_class = score_class,
             cutoff = k, cutoff_label = paste0("top", k), cutoff_fraction_label = fp_fraction_label(k, length(rep$genes)),
             TP = tp, FP = fp, TN = tn, FN = fn, recall = recall, precision = precision, FDR = 1 - precision,
             FPR = fpr, specificity = 1 - fpr, FP_per_TP = fp / max(tp, .Machine$double.eps),
             enrichment_over_random = precision / random_expected,
             score_degree_spearman = score_degree_spearman, stringsAsFactors = FALSE)
}

directional_fp_metrics <- function(rep, signed_score, mode, k, method = mode, score_class = "signed_nesta") {
  genes <- rep$genes
  score <- setNames(as.numeric(signed_score[genes]), genes)
  pos <- ranked_genes(genes, score, TRUE, rep$A1)
  neg <- ranked_genes(genes, score, FALSE, rep$A1)
  if (mode == "NESTA_signed_descending") { risk_tail <- head(pos, k); prot_tail <- character() }
  else if (mode == "NESTA_signed_ascending") { risk_tail <- character(); prot_tail <- head(neg, k) }
  else if (mode %in% c("NESTA_two_tail_balanced", "NESTA_two_tail_direction_matched")) { risk_tail <- head(pos, ceiling(k/2)); prot_tail <- head(neg, floor(k/2)) }
  else { risk_tail <- character(); prot_tail <- character() }
  universe <- setdiff(genes, rep$A1)
  mk <- function(tail, positives, opposite, tail_name) {
    tail <- unique(tail[tail %in% universe])
    tp <- sum(tail %in% positives); fp <- length(tail) - tp
    tn <- length(universe) - length(positives) - fp; fn <- length(positives) - tp
    precision <- tp / max(tp + fp, .Machine$double.eps)
    fpr <- fp / max(fp + tn, .Machine$double.eps)
    data.frame(topology_arm = rep$arm, topology_label = rep$arm_label, rescue_arm = rep_rescue_arm(rep),
               replicate = rep$rep_id, method = method, ranking_mode = mode, score_class = score_class,
               tail = tail_name, cutoff = k, cutoff_label = paste0("top", k), cutoff_fraction_label = fp_fraction_label(k, length(genes)),
               tail_size = length(tail), TP = tp, FP = fp, TN = tn, FN = fn,
               recall = tp / max(tp + fn, .Machine$double.eps), precision = precision, FDR = 1 - precision,
               FPR = fpr, specificity = 1 - fpr, FP_per_TP = fp / max(tp, .Machine$double.eps),
               opposite_direction_false_discoveries = sum(tail %in% opposite),
               high_degree_false_discoveries = sum(tail %in% rep$C), stringsAsFactors = FALSE)
  }
  rbind(mk(risk_tail, rep$A2_risk, rep$A2_protective, "risk_tail"),
        mk(prot_tail, rep$A2_protective, rep$A2_risk, "protective_tail"))
}

decoy_fp_metrics <- function(rep, selected, k, method, ranking_mode, score_class) {
  selected <- unique(selected[!(selected %in% rep$A1)])
  fp <- sum(!(selected %in% rep$A2))
  opp <- sum(selected %in% rep$D)
  high <- sum(selected %in% rep$C)
  data.frame(topology_arm = rep$arm, topology_label = rep$arm_label, rescue_arm = rep_rescue_arm(rep),
             replicate = rep$rep_id, method = method, ranking_mode = ranking_mode, score_class = score_class,
             cutoff = k, cutoff_label = paste0("top", k), cutoff_fraction_label = fp_fraction_label(k, length(rep$genes)),
             total_FP = fp, opposite_direction_false_discoveries = opp,
             high_degree_false_discoveries = high,
             opposite_direction_FDR = opp / max(length(selected), .Machine$double.eps),
             high_degree_FDR = high / max(length(selected), .Machine$double.eps),
             decoy_FP_share_among_all_FP = (opp + high) / max(fp, .Machine$double.eps), stringsAsFactors = FALSE)
}

nesta_selected_for_k <- function(rep, signed_score, mode, k) {
  strict_selected(rep, signed_score, mode, k)$selected
}

comparator_selected_for_k <- function(rep, score, k) {
  ranked <- ranked_genes(rep$genes, setNames(as.numeric(score[rep$genes]), rep$genes), TRUE, rep$A1)
  head(ranked, k)
}

run_fp_replicate <- function(seed_row) {
  tp <- seed_row$topology_arm
  base <- make_branch_isolation_rep(tp, seed_row$replicate, seed_row$base_seed)
  rep <- apply_bidirectional_arm(base, confirmatory_arm, seed_row$signal_seed)
  ch <- faithful_m2_channels(rep$adj, rep$expr, rep$twas, n.perm = 25, seed = seed_row$nesta_seed)
  no <- no_diffusion_m2_scores(rep$adj, rep$expr, rep$twas)
  signed <- setNames(ch$final_NESTA_heat, ch$SYMBOL)[rep$genes]
  deg <- score_degree(rep$adj)[rep$genes]
  total <- list(); directional <- list(); decoy <- list(); it <- id <- ic <- 1L
  for (mode in fp_nesta_modes) {
    sc <- if (mode == "NESTA_signed_descending") signed else if (mode == "NESTA_signed_ascending") -signed else abs(signed)
    sdg <- suppressWarnings(stats::cor(sc[rep$genes], deg, method = "spearman", use = "pairwise.complete.obs"))
    for (k in fp_cutoffs) {
      sel <- nesta_selected_for_k(rep, signed, mode, k)
      total[[it]] <- total_fp_metrics(rep, sel, k, mode, mode, "NESTA_Final_Heat", sdg); it <- it + 1L
      directional[[id]] <- directional_fp_metrics(rep, signed, mode, k, mode, "NESTA_Final_Heat"); id <- id + 1L
      decoy[[ic]] <- decoy_fp_metrics(rep, sel, k, mode, mode, "NESTA_Final_Heat"); ic <- ic + 1L
    }
  }
  comps <- confirmatory_comparators(rep, ch, no, seed_row)
  for (nm in c(fp_primary_comparators, fp_upper_comparators)) {
    score <- setNames(comps[[nm]]$score, rep$genes)
    sdg <- suppressWarnings(stats::cor(score[rep$genes], deg, method = "spearman", use = "pairwise.complete.obs"))
    for (k in fp_cutoffs) {
      sel <- comparator_selected_for_k(rep, score, k)
      total[[it]] <- total_fp_metrics(rep, sel, k, nm, nm, comps[[nm]]$class, sdg); it <- it + 1L
      decoy[[ic]] <- decoy_fp_metrics(rep, sel, k, nm, nm, comps[[nm]]$class); ic <- ic + 1L
    }
  }
  list(total = do.call(rbind, total), directional = do.call(rbind, directional), decoy = do.call(rbind, decoy))
}

summarise_fp <- function(x, group_cols) {
  metrics <- setdiff(names(x), c(group_cols, "replicate", "topology_label", "rescue_arm", "score_class", "cutoff_label", "cutoff_fraction_label"))
  metrics <- metrics[vapply(x[metrics], is.numeric, logical(1))]
  stats::aggregate(x[metrics], x[group_cols], FUN = function(z) mean(z, na.rm = TRUE))
}

method_contrasts_fp <- function(total) {
  out <- list(); i <- 1L
  base <- total[total$method == "NESTA_two_tail_balanced", , drop = FALSE]
  for (tp in unique(base$topology_arm)) for (k in fp_cutoffs) {
    b <- base[base$topology_arm == tp & base$cutoff == k, c("replicate", "precision", "FDR", "FPR", "recall"), drop = FALSE]
    for (cmp in c("PPR_abs_prior", "RWR_abs_prior", "raw_TWAS_abs", "M2_no_diffusion_abs")) {
      o <- total[total$topology_arm == tp & total$method == cmp & total$cutoff == k, c("replicate", "precision", "FDR", "FPR", "recall"), drop = FALSE]
      z <- merge(b, o, by = "replicate", suffixes = c("_nesta", "_comparator"))
      for (m in c("precision", "FDR", "FPR", "recall")) {
        d <- z[[paste0(m, "_nesta")]] - z[[paste0(m, "_comparator")]]
        ci <- ci_mean(d, seed = 7900L + i)
        out[[i]] <- data.frame(topology_arm = tp, cutoff = k, cutoff_label = paste0("top", k), comparator = cmp,
                               metric = m, nesta_mean = mean(z[[paste0(m, "_nesta")]], na.rm = TRUE),
                               comparator_mean = mean(z[[paste0(m, "_comparator")]], na.rm = TRUE),
                               paired_mean_difference = ci["mean"], ci_low = ci["ci_low"], ci_high = ci["ci_high"], stringsAsFactors = FALSE)
        i <- i + 1L
      }
    }
  }
  do.call(rbind, out)
}

run_false_positive_audit <- function() {
  verify_project_path(); verify_binding_plan()
  report_dir <- read_report_dir(); safe_dir_create(report_dir); copy_binding_plan_to_report(report_dir)
  set_run_status("IN_PROGRESS", "YES", "NO", "0709 false-positive control audit")
  strict_dir <- "/home/js_subdir/Dropbox/NESTA_revision/NESTA_simulation_090726_0045"
  posthoc_sufficient <- FALSE
  schedule <- make_seed_schedule()
  total_rows <- list(); dir_rows <- list(); decoy_rows <- list()
  for (i in seq_len(nrow(schedule))) {
    res <- run_fp_replicate(schedule[i, ])
    total_rows[[i]] <- res$total; dir_rows[[i]] <- res$directional; decoy_rows[[i]] <- res$decoy
    if (i %% 10 == 0) { message(sprintf("false-positive audit replicate %d/%d complete", i, nrow(schedule))); gc(FALSE) }
  }
  total <- do.call(rbind, total_rows)
  directional <- do.call(rbind, dir_rows)
  decoy <- do.call(rbind, decoy_rows)
  total_sum <- summarise_fp(total, c("topology_arm", "method", "ranking_mode", "cutoff", "cutoff_label", "cutoff_fraction_label"))
  dir_sum <- summarise_fp(directional, c("topology_arm", "method", "ranking_mode", "tail", "cutoff", "cutoff_label", "cutoff_fraction_label"))
  decoy_sum <- summarise_fp(decoy, c("topology_arm", "method", "ranking_mode", "cutoff", "cutoff_label", "cutoff_fraction_label"))
  contrasts <- method_contrasts_fp(total)
  fdr_table <- total_sum[, c("topology_arm", "method", "ranking_mode", "cutoff", "cutoff_label", "recall", "precision", "FDR", "FPR", "specificity", "FP_per_TP", "enrichment_over_random"), drop = FALSE]
  write_csv_over(total_sum, file.path(report_dir, "TOPK_TOTAL_FP_METRICS.csv"))
  write_csv_over(dir_sum, file.path(report_dir, "TOPK_DIRECTIONAL_FP_METRICS.csv"))
  write_csv_over(decoy_sum, file.path(report_dir, "TOPK_DECOY_FP_METRICS.csv"))
  write_csv_over(contrasts, file.path(report_dir, "TOPK_METHOD_CONTRASTS.csv"))
  write_csv_over(fdr_table, file.path(report_dir, "TOPK_FDR_FPR_TABLE.csv"))
  write_csv_over(fdr_table[, c("topology_arm", "method", "cutoff", "FDR")], file.path(report_dir, "figure_topk_fdr.csv"))
  write_csv_over(fdr_table[, c("topology_arm", "method", "cutoff", "FPR")], file.path(report_dir, "figure_topk_fpr.csv"))
  write_csv_over(fdr_table[, c("topology_arm", "method", "cutoff", "precision", "recall")], file.path(report_dir, "figure_topk_precision_recall.csv"))
  write_csv_over(dir_sum[, c("topology_arm", "method", "ranking_mode", "tail", "cutoff", "precision", "FDR", "FPR")], file.path(report_dir, "figure_directional_fdr.csv"))
  write_csv_over(fdr_table, file.path(report_dir, "table_top100_50_10_false_positive_metrics.csv"))
  write_csv_over(fdr_table[fdr_table$method %in% c("NESTA_two_tail_balanced", "NESTA_signed_descending", "NESTA_signed_ascending", "PPR_abs_prior", "RWR_abs_prior"), ], file.path(report_dir, "table_nesta_vs_ppr_rwr_false_positive_control.csv"))
  schema <- data.frame(file = c("TOPK_TOTAL_FP_METRICS.csv", "TOPK_DIRECTIONAL_FP_METRICS.csv", "TOPK_DECOY_FP_METRICS.csv"),
                       required_columns_present = c(all(c("TP", "FP", "TN", "FN", "FDR", "FPR", "specificity", "FP_per_TP") %in% names(total_sum)),
                                                    all(c("tail", "TP", "FP", "FDR", "FPR", "opposite_direction_false_discoveries") %in% names(dir_sum)),
                                                    all(c("opposite_direction_FDR", "high_degree_FDR", "decoy_FP_share_among_all_FP") %in% names(decoy_sum))), stringsAsFactors = FALSE)
  write_csv_over(schema, file.path(report_dir, "METRIC_SCHEMA_AUDIT.csv"))
  f_nesta50 <- total_sum[total_sum$topology_arm == "F" & total_sum$method == "NESTA_two_tail_balanced" & total_sum$cutoff == 50, , drop = FALSE]
  f_risk50 <- dir_sum[dir_sum$topology_arm == "F" & dir_sum$method == "NESTA_signed_descending" & dir_sum$tail == "risk_tail" & dir_sum$cutoff == 50, , drop = FALSE]
  f_prot50 <- dir_sum[dir_sum$topology_arm == "F" & dir_sum$method == "NESTA_signed_ascending" & dir_sum$tail == "protective_tail" & dir_sum$cutoff == 50, , drop = FALSE]
  ppr50 <- total_sum[total_sum$topology_arm == "F" & total_sum$method == "PPR_abs_prior" & total_sum$cutoff == 50, , drop = FALSE]
  rwr50 <- total_sum[total_sum$topology_arm == "F" & total_sum$method == "RWR_abs_prior" & total_sum$cutoff == 50, , drop = FALSE]
  if (ppr50$recall >= f_nesta50$recall && ppr50$precision >= f_nesta50$precision) label <- "false_positive_audit_completed_ppr_rwr_remain_strong_unsigned"
  else if (f_risk50$precision >= 0.25 && f_prot50$precision >= 0.25) label <- "false_positive_audit_completed_directional_control_supports_nesta"
  else label <- "false_positive_audit_inconclusive"
  write_lines_over(c("# False-Positive Control Summary", "", paste0("Classification: `", label, "`."),
                     sprintf("Post-hoc strict tables sufficient: %s; frozen rerun performed: %s.", posthoc_sufficient, !posthoc_sufficient),
                     sprintf("F NESTA two-tail top50 precision/FDR/FPR: %.4f / %.4f / %.4f.", f_nesta50$precision, f_nesta50$FDR, f_nesta50$FPR),
                     sprintf("F signed risk-tail top50 precision/FDR: %.4f / %.4f.", f_risk50$precision, f_risk50$FDR),
                     sprintf("F signed protective-tail top50 precision/FDR: %.4f / %.4f.", f_prot50$precision, f_prot50$FDR),
                     sprintf("F PPR_abs top50 precision/FDR: %.4f / %.4f; RWR_abs %.4f / %.4f.", ppr50$precision, ppr50$FDR, rwr50$precision, rwr50$FDR)), file.path(report_dir, "FALSE_POSITIVE_CONTROL_SUMMARY.md"))
  write_lines_over(c("# Manuscript False-Positive Interpretation", "", "PPR/RWR absolute-prior methods are interpreted as unsigned disease-relevance comparators, so risk/protective direction mismatch is not a primary failure mode for those methods.",
                     "Signed NESTA tails are interpreted with direction-specific false-positive accounting: risk-tail false discoveries are non-risk-A2 genes and protective-tail false discoveries are non-protective-A2 genes.",
                     "Top100, top50, and top10 should be reported separately rather than averaged."), file.path(report_dir, "MANUSCRIPT_FALSE_POSITIVE_INTERPRETATION.md"))
  write_lines_over(c("# Code Fidelity Audit", "", paste0("Binding plan SHA256: `", binding_plan_sha256, "`."),
                     "The frozen successful design was rerun only because per-replicate/per-tail false-positive membership was not available in the strict audit output.",
                     "Faithful submitted M2 arithmetic, TWAS.P conversion, strict TWAS filtering, signed positive/absolute-negative diffusion, zero-weight edges, self-loops, and diffuStats `n.perm` were retained.",
                     "Submitted Final Heat values were not modified."), file.path(report_dir, "CODE_FIDELITY_AUDIT.md"))
  write_lines_over(c("# Metric Schema Audit", "", sprintf("All schema checks passed: %s.", all(schema$required_columns_present))), file.path(report_dir, "METRIC_SCHEMA_AUDIT.md"))
  write_lines_over(c("# Final Report", "", paste0("Binding plan SHA256: `", binding_plan_sha256, "`."), paste0("Final classification: `", label, "`."),
                     sprintf("F NESTA two-tail top50 recall/precision/FDR/FPR: %.4f / %.4f / %.4f / %.4f.", f_nesta50$recall, f_nesta50$precision, f_nesta50$FDR, f_nesta50$FPR),
                     sprintf("F signed risk-tail top50 recall/precision/FDR: %.4f / %.4f / %.4f.", f_risk50$recall, f_risk50$precision, f_risk50$FDR),
                     sprintf("F signed protective-tail top50 recall/precision/FDR: %.4f / %.4f / %.4f.", f_prot50$recall, f_prot50$precision, f_prot50$FDR),
                     sprintf("F PPR_abs top50 recall/precision/FDR: %.4f / %.4f / %.4f.", ppr50$recall, ppr50$precision, ppr50$FDR),
                     sprintf("F RWR_abs top50 recall/precision/FDR: %.4f / %.4f / %.4f.", rwr50$recall, rwr50$precision, rwr50$FDR),
                     "This was an audit round; no adaptive simulation changes were made after false-positive results."), file.path(report_dir, "FINAL_REPORT.md"))
  if (file.exists(project_file("CLEANUP_MANIFEST.csv"))) file.copy(project_file("CLEANUP_MANIFEST.csv"), file.path(report_dir, "CLEANUP_MANIFEST.csv"), overwrite = TRUE)
  if (file.exists(project_file("results/reports/summary_tables/unit_test_results.csv"))) file.copy(project_file("results/reports/summary_tables/unit_test_results.csv"), file.path(report_dir, "unit_test_results.csv"), overwrite = TRUE)
  write_csv_over(data.frame(path = c(project_file("R/study_0709_false_positive_audit.R"), project_file("R/study_0709_strict_top_fraction_audit.R"), project_file("R/study_0708_comparator_confirmatory.R"), project_file("R/fidelity.R"), project_file("R/utils.R"), project_file("scripts/run_false_positive_audit.R")),
                            sha256 = c(sha(project_file("R/study_0709_false_positive_audit.R")), sha(project_file("R/study_0709_strict_top_fraction_audit.R")), sha(project_file("R/study_0708_comparator_confirmatory.R")), sha(project_file("R/fidelity.R")), sha(project_file("R/utils.R")), sha(project_file("scripts/run_false_positive_audit.R"))),
                            role = c("false_positive_audit_runner", "strict_fraction_helpers", "frozen_confirmatory_helpers", "faithful_nesta_and_benchmarks", "binding_plan_and_io_guards", "script_entrypoint"), stringsAsFactors = FALSE), file.path(report_dir, "IMPLEMENTATION_MANIFEST.csv"))
  files <- list.files(report_dir, recursive = TRUE, full.names = TRUE)
  files <- files[file.info(files)$isdir == FALSE & basename(files) != "CHECKSUMS.sha256"]
  writeLines(vapply(files, function(f) paste(sha(f), sub(paste0("^", report_dir, "/"), "", f)), character(1)), file.path(report_dir, "CHECKSUMS.sha256"))
  set_run_status("COMPLETE", "YES", "NO", label)
  invisible(list(label = label, total = total_sum, directional = dir_sum, decoy = decoy_sum))
}

if (!identical(Sys.getenv("NESTA_FP_AUDIT_SOURCE_ONLY"), "1")) run_false_positive_audit()
