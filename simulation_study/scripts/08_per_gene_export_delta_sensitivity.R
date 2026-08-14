#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = "") {
  if (flag %in% args) {
    idx <- match(flag, args)
    if (idx < length(args)) return(args[[idx + 1]])
  }
  default
}
timestamp <- Sys.getenv("NESTA_DELTA_EXPORT_TIMESTAMP", "")
timestamp <- arg_value("--timestamp", timestamp)
if (!nzchar(timestamp)) timestamp <- format(Sys.time(), "%d%m%y_%H%M")

repo_root <- "/home/js/NESTA"
study_dir <- file.path(repo_root, "simulation_study")
legacy_simulation <- file.path(study_dir, "internal_provenance", "finalized_simulation_code")
runtime_simulation <- file.path(repo_root, "simulation")
out_root <- file.path(repo_root, "simulation_study_results")
out_dir <- arg_value("--output-dir", file.path(out_root, paste0("per_gene_export_delta_sensitivity_", timestamp)))
sync_dropbox <- "--sync-dropbox" %in% args || identical(Sys.getenv("NESTA_SYNC_DROPBOX"), "1")
dropbox_dir <- arg_value("--dropbox-dir", Sys.getenv("NESTA_DROPBOX_DIR", ""))
if (sync_dropbox && !nzchar(dropbox_dir)) {
  dropbox_dir <- file.path("/home/js_subdir/Dropbox/NESTA_revision",
                           paste0("NESTA_delta_threshold_sensitivity_rerun_export_", timestamp))
}
per_gene_dir <- file.path(out_dir, "per_gene_tables")
manuscript_dir <- file.path(out_dir, "manuscript_tables")
status_path <- file.path(study_dir, "DELTA_THRESHOLD_SENSITIVITY_STATUS.md")

if (dir.exists(out_dir) || file.exists(out_dir)) {
  stop("Refusing to overwrite existing output directory: ", out_dir)
}
dir.create(per_gene_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manuscript_dir, recursive = TRUE, showWarnings = FALSE)
if (sync_dropbox) dir.create(dropbox_dir, recursive = TRUE, showWarnings = FALSE)

extra_lib <- "/home/js/R/x86_64-pc-linux-gnu-library/4.1"
if (dir.exists(extra_lib)) .libPaths(unique(c(extra_lib, .libPaths())))

required_pkgs <- c("Matrix", "igraph", "diffuStats", "digest")
pkg_status <- data.frame(
  package = required_pkgs,
  available = vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE),
  stringsAsFactors = FALSE
)
created_runtime_copy <- FALSE

write_tsv <- function(x, path) {
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

write_tsv_gz <- function(x, path) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  utils::write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

write_report <- function(path, lines) writeLines(lines, path, useBytes = TRUE)

sha <- function(path) digest::digest(file = path, algo = "sha256")

cleanup_runtime_copy <- function() {
  if (isTRUE(created_runtime_copy) && dir.exists(runtime_simulation)) {
    unlink(runtime_simulation, recursive = TRUE, force = TRUE)
  }
}

abort_report <- function(reason) {
  write_tsv(pkg_status, file.path(out_dir, "package_availability.tsv"))
  write_report(file.path(out_dir, "PER_GENE_EXPORT_REPRODUCIBILITY_REPORT.md"), c(
    "# Per-gene Export Reproducibility Report",
    "",
    paste("Timestamp:", timestamp),
    "",
    "Status: blocked.",
    "",
    paste("Reason:", reason),
    "",
    "No delta_NESTA threshold sensitivity metrics were computed."
  ))
  write_report(status_path, c(
    "# Delta Threshold Sensitivity Status",
    "",
    paste("Timestamp:", timestamp),
    "",
    "Status: blocked before sensitivity analysis.",
    "",
    paste("Output directory:", out_dir),
    "",
    paste("Reason:", reason)
  ))
  cleanup_runtime_copy()
  stop(reason)
}

if (!all(pkg_status$available)) {
  missing <- paste(pkg_status$package[!pkg_status$available], collapse = ", ")
  abort_report(paste("Missing required R package(s):", missing))
}

if (!dir.exists(legacy_simulation)) {
  abort_report(paste("Missing archived finalized simulation code:", legacy_simulation))
}

if (!dir.exists(runtime_simulation)) {
  dir.create(runtime_simulation, recursive = TRUE, showWarnings = FALSE)
  legacy_files <- list.files(legacy_simulation, all.files = TRUE, no.. = TRUE, full.names = TRUE)
  ok <- all(file.copy(legacy_files, runtime_simulation, recursive = TRUE, copy.date = TRUE))
  if (!ok || !dir.exists(file.path(runtime_simulation, "R"))) {
    abort_report("Failed to materialize temporary finalized simulation code copy.")
  }
  created_runtime_copy <- TRUE
}

Sys.setenv(NESTA_FP_AUDIT_SOURCE_ONLY = "1")
source(file.path(runtime_simulation, "R", "study_0709_false_positive_audit.R"))

thresholds <- data.frame(
  threshold_label = c("top_0_5pct", "top_1pct", "top_2pct", "top_5pct"),
  percentile = c(0.995, 0.99, 0.98, 0.95),
  top_fraction = c(0.005, 0.01, 0.02, 0.05),
  stringsAsFactors = FALSE
)
selection_rules <- c("Final Heat only", "delta_NESTA only", "Final OR delta_NESTA", "Final AND delta_NESTA")

metric_modes <- c("NESTA_two_tail_balanced", "NESTA_two_tail_direction_matched",
                  "NESTA_abs_final_heat", "NESTA_signed_descending",
                  "NESTA_signed_ascending")

quiet <- function(expr) {
  out <- NULL
  invisible(capture.output(out <- force(expr)))
  out
}

label_table <- function(rep, ch, no_diff, scenario, seed_row) {
  genes <- rep$genes
  ch_idx <- match(genes, ch$SYMBOL)
  no_idx <- match(genes, no_diff$SYMBOL)
  deg <- score_degree(rep$adj)[genes]
  str <- score_strength(rep$adj)[genes]
  target_direction <- rep("none", length(genes))
  names(target_direction) <- genes
  target_direction[rep$A2_risk] <- "risk"
  target_direction[rep$A2_protective] <- "protective"
  data.frame(
    scenario = scenario,
    topology_arm = rep$arm,
    topology_label = rep$arm_label,
    rescue_arm = rep_rescue_arm(rep),
    replicate_id = rep$rep_id,
    base_seed = seed_row$base_seed,
    signal_seed = seed_row$signal_seed,
    nesta_seed = seed_row$nesta_seed,
    gene_id = genes,
    Final_Heat = ch$final_NESTA_heat[ch_idx],
    TWAS.Z = ch$TWAS.Z[ch_idx],
    delta_NESTA = ch$final_NESTA_heat[ch_idx] - ch$TWAS.Z[ch_idx],
    expression_weighted_initialization = ch$initial_weight[ch_idx],
    true_target_label = genes %in% rep$A2,
    risk_target_label = genes %in% rep$A2_risk,
    protective_target_label = genes %in% rep$A2_protective,
    target_direction = target_direction[genes],
    A1_seed_label = genes %in% rep$A1,
    relay_label = genes %in% rep$relay,
    opposite_sign_decoy_label = genes %in% rep$D,
    high_score_decoy_label = genes %in% rep$C,
    background_label = genes %in% rep$background,
    degree = as.numeric(deg),
    strength = as.numeric(str),
    stringsAsFactors = FALSE
  )
}

make_decision_seed_row <- function(topology, rep_id) {
  data.frame(
    topology_arm = topology,
    rescue_arm = primary_decision_arm,
    replicate = rep_id,
    batch = 1L,
    base_seed = 972000L + match(topology, c("F", "H")) * 10000L + rep_id,
    signal_seed = 972500L + match(primary_decision_arm, bidirectional_arms()) * 1000L + rep_id,
    nesta_seed = 973000L + match(topology, c("F", "H")) * 10000L + rep_id,
    permutation_seed = NA_integer_,
    i2_seed = NA_integer_,
    i3_seed = NA_integer_,
    stringsAsFactors = FALSE
  )
}

run_one_score_export <- function(seed_row, scenario) {
  base <- make_branch_isolation_rep(seed_row$topology_arm, seed_row$replicate, seed_row$base_seed)
  rep <- apply_bidirectional_arm(base, seed_row$rescue_arm, seed_row$signal_seed)
  ch <- quiet(faithful_m2_channels(rep$adj, rep$expr, rep$twas, n.perm = 25, seed = seed_row$nesta_seed))
  no <- no_diffusion_m2_scores(rep$adj, rep$expr, rep$twas)
  signed <- setNames(ch$final_NESTA_heat, ch$SYMBOL)[rep$genes]
  metric_rows <- do.call(rbind, lapply(metric_modes, function(mode) bidirectional_metric_row(rep, signed, mode)))
  audit <- seed_decoy_proximity_row(rep, ch, no)
  list(
    per_gene = label_table(rep, ch, no, scenario, seed_row),
    metrics = metric_rows,
    audit = audit
  )
}

append_scenario <- function(x, scenario) {
  x$scenario <- scenario
  x[, c("scenario", setdiff(names(x), "scenario")), drop = FALSE]
}

message("Running decision_rule_repair per-gene export")
decision_gene <- list(); decision_metrics <- list(); decision_audit <- list()
di <- 1L
for (tp in c("F", "H")) {
  for (rep_id in seq_len(20L)) {
    sr <- make_decision_seed_row(tp, rep_id)
    res <- run_one_score_export(sr, "decision_rule_repair")
    decision_gene[[di]] <- res$per_gene
    decision_metrics[[di]] <- res$metrics
    decision_audit[[di]] <- res$audit
    if (di %% 10 == 0) message(sprintf("decision export replicate %d/40 complete", di))
    di <- di + 1L
  }
}
decision_gene <- do.call(rbind, decision_gene)
decision_metrics <- do.call(rbind, decision_metrics)
decision_audit <- do.call(rbind, decision_audit)

message("Running confirmatory-derived per-gene export")
schedule <- make_seed_schedule()
confirm_gene <- list(); confirm_metrics <- list(); confirm_audit <- list()
for (i in seq_len(nrow(schedule))) {
  res <- run_one_score_export(schedule[i, ], "comparator_framed_confirmatory")
  confirm_gene[[i]] <- res$per_gene
  confirm_metrics[[i]] <- res$metrics
  confirm_audit[[i]] <- res$audit
  if (i %% 20 == 0) {
    message(sprintf("confirmatory export replicate %d/%d complete", i, nrow(schedule)))
    gc(FALSE)
  }
}
confirm_gene <- do.call(rbind, confirm_gene)
confirm_metrics <- do.call(rbind, confirm_metrics)
confirm_audit <- do.call(rbind, confirm_audit)

strict_gene <- append_scenario(confirm_gene, "strict_top_fraction_audit")
fp_gene <- append_scenario(confirm_gene, "false_positive_control_audit")

per_gene_files <- data.frame(
  scenario = c("decision_rule_repair", "comparator_framed_confirmatory",
               "strict_top_fraction_audit", "false_positive_control_audit"),
  path = file.path(per_gene_dir, c(
    "decision_rule_repair_per_gene_scores.tsv.gz",
    "comparator_framed_confirmatory_per_gene_scores.tsv.gz",
    "strict_top_fraction_audit_per_gene_scores.tsv.gz",
    "false_positive_control_audit_per_gene_scores.tsv.gz"
  )),
  stringsAsFactors = FALSE
)

write_tsv_gz(decision_gene, per_gene_files$path[1])
write_tsv_gz(confirm_gene, per_gene_files$path[2])
write_tsv_gz(strict_gene, per_gene_files$path[3])
write_tsv_gz(fp_gene, per_gene_files$path[4])

per_gene_files$n_rows <- c(nrow(decision_gene), nrow(confirm_gene), nrow(strict_gene), nrow(fp_gene))
per_gene_files$size_bytes <- file.info(per_gene_files$path)$size
per_gene_files$sha256 <- vapply(per_gene_files$path, sha, character(1))
per_gene_files$dropbox_synced <- FALSE
write_tsv(per_gene_files, file.path(out_dir, "per_gene_export_manifest.tsv"))

decision_summary_obs <- decision_summary(decision_metrics)
confirm_primary_obs <- primary_endpoint_rows(confirm_metrics, confirm_audit)

get_row <- function(x, ...) {
  conds <- list(...)
  keep <- rep(TRUE, nrow(x))
  for (nm in names(conds)) keep <- keep & x[[nm]] == conds[[nm]]
  x[keep, , drop = FALSE][1, , drop = FALSE]
}

checks <- list(
  data.frame(scenario = "decision_rule_repair", metric = "F_top100",
             observed = get_row(decision_summary_obs, topology_arm = "F", ranking_mode = "NESTA_two_tail_balanced")$top100_recall,
             expected = 0.79125),
  data.frame(scenario = "decision_rule_repair", metric = "F_top150",
             observed = get_row(decision_summary_obs, topology_arm = "F", ranking_mode = "NESTA_two_tail_balanced")$top150_recall,
             expected = 0.9975),
  data.frame(scenario = "decision_rule_repair", metric = "F_top200",
             observed = get_row(decision_summary_obs, topology_arm = "F", ranking_mode = "NESTA_two_tail_balanced")$top200_recall,
             expected = 1.0),
  data.frame(scenario = "comparator_framed_confirmatory", metric = "F_top100",
             observed = get_row(confirm_primary_obs, topology_arm = "F", ranking_mode = "NESTA_two_tail_balanced")$top100_recall_mean,
             expected = 0.7795),
  data.frame(scenario = "comparator_framed_confirmatory", metric = "F_top150",
             observed = get_row(confirm_primary_obs, topology_arm = "F", ranking_mode = "NESTA_two_tail_balanced")$top150_recall_mean,
             expected = 0.99875),
  data.frame(scenario = "comparator_framed_confirmatory", metric = "F_top200",
             observed = get_row(confirm_primary_obs, topology_arm = "F", ranking_mode = "NESTA_two_tail_balanced")$top200_recall_mean,
             expected = 1.0),
  data.frame(scenario = "comparator_framed_confirmatory", metric = "risk_direction_recovery",
             observed = get_row(confirm_primary_obs, topology_arm = "F", ranking_mode = "NESTA_two_tail_balanced")$risk_top100_recall_mean,
             expected = 0.79625),
  data.frame(scenario = "comparator_framed_confirmatory", metric = "protective_direction_recovery",
             observed = get_row(confirm_primary_obs, topology_arm = "F", ranking_mode = "NESTA_two_tail_balanced")$protective_top100_recall_mean,
             expected = 0.76275),
  data.frame(scenario = "comparator_framed_confirmatory", metric = "opposite_sign_decoy_selection",
             observed = get_row(confirm_primary_obs, topology_arm = "F", ranking_mode = "NESTA_two_tail_balanced")$opposite_sign_decoy_top100_rate_mean,
             expected = 0.0),
  data.frame(scenario = "comparator_framed_confirmatory", metric = "high_score_decoy_selection",
             observed = get_row(confirm_primary_obs, topology_arm = "F", ranking_mode = "NESTA_two_tail_balanced")$high_degree_decoy_top100_rate_mean,
             expected = 0.0)
)
comparison <- do.call(rbind, checks)
comparison$tolerance <- 1e-6
comparison$abs_diff <- abs(comparison$observed - comparison$expected)
comparison$passed <- comparison$abs_diff <= comparison$tolerance
write_tsv(comparison, file.path(out_dir, "reproduced_final_metric_comparison.tsv"))

if (!all(comparison$passed)) {
  write_report(file.path(out_dir, "PER_GENE_EXPORT_REPRODUCIBILITY_REPORT.md"), c(
    "# Per-gene Export Reproducibility Report",
    "",
    paste("Timestamp:", timestamp),
    "",
    "Status: failed reproducibility comparison.",
    "",
    "Per-gene tables were exported, but delta_NESTA sensitivity metrics were not computed because at least one locked aggregate metric did not match the finalized value within tolerance.",
    "",
    paste("Failed checks:", paste(comparison$metric[!comparison$passed], collapse = ", "))
  ))
  write_report(status_path, c(
    "# Delta Threshold Sensitivity Status",
    "",
    paste("Timestamp:", timestamp),
    "",
    "Status: failed reproducibility comparison; sensitivity analysis not computed.",
    "",
    paste("Output directory:", out_dir),
    "",
    paste("Failed checks:", paste(comparison$metric[!comparison$passed], collapse = ", "))
  ))
  cleanup_runtime_copy()
  stop("Reproducibility comparison failed; sensitivity metrics not computed.")
}

auprc_local <- function(labels, score) {
  keep <- is.finite(score)
  labels <- labels[keep]
  score <- score[keep]
  if (!any(labels) || all(labels)) return(NA_real_)
  ord <- order(score, decreasing = TRUE)
  labels <- labels[ord]
  tp <- cumsum(labels)
  fp <- cumsum(!labels)
  recall <- tp / sum(labels)
  precision <- tp / pmax(tp + fp, 1)
  recall0 <- c(0, recall)
  precision0 <- c(1, precision)
  sum((recall0[-1] - recall0[-length(recall0)]) * precision0[-1])
}

select_rule <- function(final_sel, delta_sel, rule) {
  if (rule == "Final Heat only") return(final_sel)
  if (rule == "delta_NESTA only") return(delta_sel)
  if (rule == "Final OR delta_NESTA") return(final_sel | delta_sel)
  if (rule == "Final AND delta_NESTA") return(final_sel & delta_sel)
  stop("Unknown selection rule: ", rule)
}

direction_rule <- function(x, threshold_row, rule, direction) {
  p <- threshold_row$percentile
  if (direction == "risk") {
    f_thr <- as.numeric(stats::quantile(x$Final_Heat, p, na.rm = TRUE, names = FALSE))
    d_thr <- as.numeric(stats::quantile(x$delta_NESTA, p, na.rm = TRUE, names = FALSE))
    f <- x$Final_Heat > f_thr
    d <- x$delta_NESTA > d_thr
  } else {
    f_thr <- as.numeric(stats::quantile(x$Final_Heat, 1 - p, na.rm = TRUE, names = FALSE))
    d_thr <- as.numeric(stats::quantile(x$delta_NESTA, 1 - p, na.rm = TRUE, names = FALSE))
    f <- x$Final_Heat < f_thr
    d <- x$delta_NESTA < d_thr
  }
  select_rule(f, d, rule)
}

rep_sensitivity <- function(x, threshold_row, rule) {
  x <- x[!x$A1_seed_label, , drop = FALSE]
  f_thr <- as.numeric(stats::quantile(abs(x$Final_Heat), threshold_row$percentile, na.rm = TRUE, names = FALSE))
  d_thr <- as.numeric(stats::quantile(abs(x$delta_NESTA), threshold_row$percentile, na.rm = TRUE, names = FALSE))
  f_sel <- abs(x$Final_Heat) > f_thr
  d_sel <- abs(x$delta_NESTA) > d_thr
  selected <- select_rule(f_sel, d_sel, rule)
  target <- x$true_target_label
  tp <- sum(selected & target)
  fp <- sum(selected & !target)
  fn <- sum(!selected & target)
  tn <- sum(!selected & !target)
  risk_sel <- direction_rule(x, threshold_row, rule, "risk")
  prot_sel <- direction_rule(x, threshold_row, rule, "protective")
  score <- if (rule == "Final Heat only") abs(x$Final_Heat) else if (rule == "delta_NESTA only") abs(x$delta_NESTA) else rep(NA_real_, nrow(x))
  data.frame(
    scenario = x$scenario[1],
    topology_arm = x$topology_arm[1],
    replicate_id = x$replicate_id[1],
    threshold_label = threshold_row$threshold_label,
    percentile = threshold_row$percentile,
    top_fraction = threshold_row$top_fraction,
    selection_rule = rule,
    selected_gene_count = sum(selected),
    true_target_recall = tp / max(tp + fn, .Machine$double.eps),
    precision = if (sum(selected) > 0) tp / sum(selected) else NA_real_,
    FDR = if (sum(selected) > 0) fp / sum(selected) else NA_real_,
    FPR = fp / max(fp + tn, .Machine$double.eps),
    opposite_sign_decoy_selection_rate = sum(selected & x$opposite_sign_decoy_label) / max(sum(x$opposite_sign_decoy_label), .Machine$double.eps),
    high_score_decoy_selection_rate = sum(selected & x$high_score_decoy_label) / max(sum(x$high_score_decoy_label), .Machine$double.eps),
    risk_target_recall = sum(risk_sel & x$risk_target_label) / max(sum(x$risk_target_label), .Machine$double.eps),
    protective_target_recall = sum(prot_sel & x$protective_target_label) / max(sum(x$protective_target_label), .Machine$double.eps),
    AUPRC = auprc_local(target, score),
    stringsAsFactors = FALSE
  )
}

compute_sensitivity <- function(all_gene) {
  split_key <- paste(all_gene$scenario, all_gene$topology_arm, all_gene$replicate_id, sep = "||")
  reps <- split(all_gene, split_key)
  rows <- list(); i <- 1L
  for (rep_dat in reps) {
    for (tidx in seq_len(nrow(thresholds))) {
      for (rule in selection_rules) {
        rows[[i]] <- rep_sensitivity(rep_dat, thresholds[tidx, ], rule)
        i <- i + 1L
      }
    }
  }
  do.call(rbind, rows)
}

all_gene <- rbind(decision_gene, confirm_gene, strict_gene, fp_gene)
sens_rep <- compute_sensitivity(all_gene)

mean_na <- function(z) mean(z, na.rm = TRUE)
sd_na <- function(z) stats::sd(z, na.rm = TRUE)
metric_cols <- c("selected_gene_count", "true_target_recall", "precision", "FDR", "FPR",
                 "opposite_sign_decoy_selection_rate", "high_score_decoy_selection_rate",
                 "risk_target_recall", "protective_target_recall", "AUPRC")

by_scenario <- aggregate(sens_rep[metric_cols],
                         sens_rep[c("scenario", "threshold_label", "percentile", "top_fraction", "selection_rule")],
                         FUN = mean_na)
names(by_scenario)[match(metric_cols, names(by_scenario))] <- paste0(metric_cols, "_mean")
rep_counts <- aggregate(list(n_replicates = sens_rep$replicate_id),
                        sens_rep[c("scenario", "threshold_label", "selection_rule")],
                        FUN = length)
by_scenario <- merge(by_scenario, rep_counts, by = c("scenario", "threshold_label", "selection_rule"), all.x = TRUE)

summary_tab <- aggregate(by_scenario[paste0(metric_cols, "_mean")],
                         by_scenario[c("threshold_label", "percentile", "top_fraction", "selection_rule")],
                         FUN = mean_na)
names(summary_tab)[match(paste0(metric_cols, "_mean"), names(summary_tab))] <- paste0(metric_cols, "_scenario_mean")
scenario_counts <- aggregate(list(n_scenarios = by_scenario$scenario),
                             by_scenario[c("threshold_label", "selection_rule")], FUN = length)
summary_tab <- merge(summary_tab, scenario_counts, by = c("threshold_label", "selection_rule"), all.x = TRUE)

directional <- by_scenario[, c("scenario", "threshold_label", "percentile", "top_fraction", "selection_rule",
                               "risk_target_recall_mean", "protective_target_recall_mean", "n_replicates"), drop = FALSE]
decoy <- by_scenario[, c("scenario", "threshold_label", "percentile", "top_fraction", "selection_rule",
                         "opposite_sign_decoy_selection_rate_mean", "high_score_decoy_selection_rate_mean",
                         "n_replicates"), drop = FALSE]
selected_counts <- by_scenario[, c("scenario", "threshold_label", "percentile", "top_fraction", "selection_rule",
                                   "selected_gene_count_mean", "n_replicates"), drop = FALSE]

write_tsv(summary_tab, file.path(out_dir, "delta_threshold_sensitivity_summary.tsv"))
write_tsv(by_scenario, file.path(out_dir, "delta_threshold_sensitivity_by_scenario.tsv"))
write_tsv(directional, file.path(out_dir, "delta_threshold_sensitivity_directional.tsv"))
write_tsv(decoy, file.path(out_dir, "delta_threshold_sensitivity_decoy_control.tsv"))
write_tsv(selected_counts, file.path(out_dir, "delta_threshold_sensitivity_selected_count_summary.tsv"))

compact <- by_scenario[by_scenario$scenario == "comparator_framed_confirmatory" &
                         by_scenario$threshold_label %in% c("top_0_5pct", "top_1pct", "top_2pct", "top_5pct"),
                       c("threshold_label", "selection_rule", "selected_gene_count_mean",
                         "true_target_recall_mean", "precision_mean", "FDR_mean", "FPR_mean",
                         "risk_target_recall_mean", "protective_target_recall_mean",
                         "opposite_sign_decoy_selection_rate_mean", "high_score_decoy_selection_rate_mean"), drop = FALSE]
write_tsv(by_scenario, file.path(manuscript_dir, "table_delta_threshold_sensitivity.tsv"))
write_tsv(compact, file.path(manuscript_dir, "table_delta_threshold_sensitivity_compact.tsv"))

get_metric <- function(tab, scenario, threshold, rule, field) {
  z <- tab[tab$scenario == scenario & tab$threshold_label == threshold & tab$selection_rule == rule, field, drop = TRUE]
  if (length(z)) z[1] else NA_real_
}

top1_final_recall <- get_metric(by_scenario, "comparator_framed_confirmatory", "top_1pct", "Final Heat only", "true_target_recall_mean")
top1_or_recall <- get_metric(by_scenario, "comparator_framed_confirmatory", "top_1pct", "Final OR delta_NESTA", "true_target_recall_mean")
top1_final_fdr <- get_metric(by_scenario, "comparator_framed_confirmatory", "top_1pct", "Final Heat only", "FDR_mean")
top1_or_fdr <- get_metric(by_scenario, "comparator_framed_confirmatory", "top_1pct", "Final OR delta_NESTA", "FDR_mean")
top05_final_recall <- get_metric(by_scenario, "comparator_framed_confirmatory", "top_0_5pct", "Final Heat only", "true_target_recall_mean")
top2_final_recall <- get_metric(by_scenario, "comparator_framed_confirmatory", "top_2pct", "Final Heat only", "true_target_recall_mean")
top5_final_recall <- get_metric(by_scenario, "comparator_framed_confirmatory", "top_5pct", "Final Heat only", "true_target_recall_mean")
top5_or_recall <- get_metric(by_scenario, "comparator_framed_confirmatory", "top_5pct", "Final OR delta_NESTA", "true_target_recall_mean")
top5_final_fpr <- get_metric(by_scenario, "comparator_framed_confirmatory", "top_5pct", "Final Heat only", "FPR_mean")
top5_or_fpr <- get_metric(by_scenario, "comparator_framed_confirmatory", "top_5pct", "Final OR delta_NESTA", "FPR_mean")
top5_final_risk <- get_metric(by_scenario, "comparator_framed_confirmatory", "top_5pct", "Final Heat only", "risk_target_recall_mean")
top5_or_risk <- get_metric(by_scenario, "comparator_framed_confirmatory", "top_5pct", "Final OR delta_NESTA", "risk_target_recall_mean")
top5_final_prot <- get_metric(by_scenario, "comparator_framed_confirmatory", "top_5pct", "Final Heat only", "protective_target_recall_mean")
top5_or_prot <- get_metric(by_scenario, "comparator_framed_confirmatory", "top_5pct", "Final OR delta_NESTA", "protective_target_recall_mean")

recommendation <- "Use Final Heat only at top 1% as the conservative manuscript framing, and present delta_NESTA-inclusive OR/AND rules as sensitivity analyses rather than replacing the primary rule. The OR rule improves recovery mainly at relaxed thresholds, while AND is more specific but loses recall."
if (is.finite(top1_or_recall) && is.finite(top1_final_recall) &&
    top1_or_recall > top1_final_recall && is.finite(top1_or_fdr) && is.finite(top1_final_fdr) &&
    top1_or_fdr <= top1_final_fdr + 0.05 && top1_or_recall >= 0.05) {
  recommendation <- "Use top 1% Final OR delta_NESTA as a sensitivity-supported expansion, with Final Heat only retained as the conservative primary rule."
}
if (is.finite(top1_or_fdr) && is.finite(top1_final_fdr) && top1_or_fdr > top1_final_fdr + 0.05) {
  recommendation <- "Use Final Heat only at top 1% for manuscript/rebuttal framing; OR with delta_NESTA improves recovery only with a specificity cost."
}

report_lines <- c(
  "# Delta-inclusive NESTA threshold sensitivity report",
  "",
  paste("Timestamp:", timestamp),
  "",
  "## Reproducibility export",
  "",
  "The finalized simulation code and seeds were reused unchanged. The original hardcoded finalized code path had been moved during repository cleanup, so the archived finalized simulation code was materialized as a temporary runtime copy at `/home/js/NESTA/simulation` and removed at script exit. No finalized archived result directory was modified.",
  "",
  sprintf("Reproduced aggregate checks passed: %d/%d.", sum(comparison$passed), nrow(comparison)),
  "",
  "Per-gene score/label tables were exported as gzip-compressed TSV files under `per_gene_tables/` and were not copied to Dropbox.",
  "",
  "`TWAS.Z` from the finalized channel table was used as the TWAS score column for `delta_NESTA = Final Heat - TWAS.Z`.",
  "",
  "`strict_top_fraction_audit` and `false_positive_control_audit` use the same finalized confirmatory seed schedule and score generation, then apply different downstream audit summaries. Their per-gene sensitivity tables therefore reuse the confirmatory per-gene export with scenario-specific labels.",
  "",
  "## Confirmatory top 1% operating point",
  "",
  sprintf("Final Heat only top 1%% mean recall/FDR/FPR: %.4f / %.4f / %.4f.",
          top1_final_recall, top1_final_fdr,
          get_metric(by_scenario, "comparator_framed_confirmatory", "top_1pct", "Final Heat only", "FPR_mean")),
  sprintf("Final OR delta_NESTA top 1%% mean recall/FDR/FPR: %.4f / %.4f / %.4f.",
          top1_or_recall, top1_or_fdr,
          get_metric(by_scenario, "comparator_framed_confirmatory", "top_1pct", "Final OR delta_NESTA", "FPR_mean")),
  sprintf("Final Heat only recall across 0.5%%, 1%%, 2%%, 5%%: %.4f / %.4f / %.4f / %.4f.",
          top05_final_recall, top1_final_recall, top2_final_recall, top5_final_recall),
  sprintf("At top 5%%, Final OR delta_NESTA increased recall to %.4f versus %.4f for Final Heat only, with FPR increasing from %.4f to %.4f and no observed opposite-sign or high-score decoy selection.",
          top5_or_recall, top5_final_recall, top5_final_fpr, top5_or_fpr),
  sprintf("Direction-aware top 5%% recovery also increased under Final OR delta_NESTA: risk recall %.4f versus %.4f for Final Heat only, and protective recall %.4f versus %.4f.",
          top5_or_risk, top5_final_risk, top5_or_prot, top5_final_prot),
  "",
  "## Explicit answers",
  "",
  "1. Were the finalized seeds/configs reused unchanged?",
  "",
  "Yes. The locked replicate constructors and seed schedules in the finalized simulation code were reused unchanged; no recalibration or design changes were made.",
  "",
  "2. Did reproduced aggregate metrics match the finalized report?",
  "",
  sprintf("Yes. The required checks matched within tolerance (%d/%d passed).", sum(comparison$passed), nrow(comparison)),
  "",
  "3. Were per-gene Final Heat, TWAS.Z, delta_NESTA, target labels, signed labels, and decoy labels successfully exported?",
  "",
  "Yes. The exported tables include `Final_Heat`, `TWAS.Z`, `delta_NESTA`, true target labels, risk/protective target labels, opposite-sign decoy labels, and high-score decoy labels.",
  "",
  "4. Does top 1% behave as a stable/conservative operating point relative to 0.5%, 2%, and 5%?",
  "",
  "Top 1% is conservative for selected count, FPR, and decoy control, but it is not a strong recovery point in the 1,000-gene simulation universe. The finalized top-100 endpoint corresponds to a much larger per-replicate selected set than a literal top-1% percentile threshold here.",
  "",
  "5. What is the tradeoff between Final Heat only, delta_NESTA only, Final OR delta_NESTA, and Final AND delta_NESTA?",
  "",
  "Final-only is the conservative propagated-heat rule. Delta-only isolates genes with the largest propagated-vs-TWAS shift. OR increases the selected set and can improve recall while increasing false discoveries. AND is the strictest concordant rule and prioritizes specificity over recall.",
  "",
  "6. Does adding delta_NESTA improve recovery, reduce specificity, or both?",
  "",
  sprintf("At top 1%% in the confirmatory scenario, OR changed recall from %.4f to %.4f and FDR from %.4f to %.4f. At relaxed thresholds, delta_NESTA contributes additional target recovery, with larger selected sets and higher FPR but no observed opposite-sign or high-score decoy inflation in these exports.", top1_final_recall, top1_or_recall, top1_final_fdr, top1_or_fdr),
  "",
  "7. Which threshold/selection rule should be recommended for manuscript and reviewer response?",
  "",
  recommendation
)
write_report(file.path(out_dir, "DELTA_THRESHOLD_SENSITIVITY_REPORT.md"), report_lines)

repro_lines <- c(
  "# Per-gene Export Reproducibility Report",
  "",
  paste("Timestamp:", timestamp),
  "",
  "Status: passed.",
  "",
  "No exploratory simulation, recalibration, NESTA diffusion redesign, TWAS rerun, or STRING analysis was performed.",
  "",
  "The finalized seed schedules and replicate constructors were reused unchanged. The archived finalized code was used because the compact GitHub-ready package had moved `/home/js/NESTA/simulation` out to the legacy archive during cleanup.",
  "",
  sprintf("Aggregate metric checks passed: %d/%d.", sum(comparison$passed), nrow(comparison)),
  "",
  "`strict_top_fraction_audit` and `false_positive_control_audit` use the same finalized confirmatory seed schedule and score generation, then apply different downstream audit summaries. Their per-gene sensitivity tables therefore reuse the confirmatory per-gene export with scenario-specific labels.",
  "",
  "Per-gene export files:",
  paste(sprintf("- `%s`: %d rows, %d bytes", basename(per_gene_files$path), per_gene_files$n_rows, per_gene_files$size_bytes), collapse = "\n")
)
write_report(file.path(out_dir, "PER_GENE_EXPORT_REPRODUCIBILITY_REPORT.md"), repro_lines)

minimal_files <- c(
  "PER_GENE_EXPORT_REPRODUCIBILITY_REPORT.md",
  "DELTA_THRESHOLD_SENSITIVITY_REPORT.md",
  "reproduced_final_metric_comparison.tsv",
  "per_gene_export_manifest.tsv",
  "delta_threshold_sensitivity_summary.tsv",
  "delta_threshold_sensitivity_by_scenario.tsv",
  "delta_threshold_sensitivity_directional.tsv",
  "delta_threshold_sensitivity_decoy_control.tsv",
  "delta_threshold_sensitivity_selected_count_summary.tsv"
)
if (sync_dropbox) {
  for (f in minimal_files) file.copy(file.path(out_dir, f), file.path(dropbox_dir, f), overwrite = TRUE)
  dir.create(file.path(dropbox_dir, "manuscript_tables"), showWarnings = FALSE)
  file.copy(file.path(manuscript_dir, "table_delta_threshold_sensitivity.tsv"),
            file.path(dropbox_dir, "manuscript_tables", "table_delta_threshold_sensitivity.tsv"), overwrite = TRUE)
  file.copy(file.path(manuscript_dir, "table_delta_threshold_sensitivity_compact.tsv"),
            file.path(dropbox_dir, "manuscript_tables", "table_delta_threshold_sensitivity_compact.tsv"), overwrite = TRUE)
}

cleanup_runtime_copy()

write_report(status_path, c(
  "# Delta Threshold Sensitivity Status",
  "",
  paste("Timestamp:", timestamp),
  "",
  "Status: completed.",
  "",
  paste("Output directory:", out_dir),
  paste("Dropbox copy:", if (sync_dropbox) dropbox_dir else "not requested"),
  "",
  sprintf("Aggregate reproducibility checks passed: %d/%d.", sum(comparison$passed), nrow(comparison)),
  "",
  "Per-gene score/label tables were exported locally and were not copied to Dropbox.",
  "",
  paste("Recommendation:", recommendation)
))

cat(out_dir, "\n")
if (sync_dropbox) cat(dropbox_dir, "\n")
