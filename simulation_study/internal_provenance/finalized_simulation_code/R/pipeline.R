source(file.path("/home/js/NESTA/simulation", "R", "utils.R"))
source(file.path("/home/js/NESTA/simulation", "R", "fidelity.R"))
source(file.path("/home/js/NESTA/simulation", "R", "tom_library.R"))
source(file.path("/home/js/NESTA/simulation", "R", "generator.R"))

path_target_sets <- function(modules) {
  list(
    A2_all = modules$A2_all,
    A2_proximal = modules$A2_proximal,
    A2_intermediate = modules$A2_intermediate,
    A2_distal = modules$A2_distal,
    A2_intermediate_degree_capped = modules$A2_intermediate_degree_capped,
    A2_risk = modules$A2_risk,
    A2_protective = modules$A2_protective,
    A2_low_degree = modules$A2_low_degree,
    A2_mid_degree = modules$A2_mid_degree,
    A2_high_degree = modules$A2_high_degree,
    A2_extreme_high_degree = modules$A2_extreme_high_degree,
    A2_intermediate_degree_capped_risk = modules$A2_intermediate_degree_capped_risk,
    A2_intermediate_degree_capped_protective = modules$A2_intermediate_degree_capped_protective,
    D_opposite_sign_decoy = modules$D_opposite_sign_decoy,
    C_high_degree_decoy = modules$C_high_degree_decoy
  )
}

ranked_genes <- function(scores, score_col, exclude = character()) {
  dat <- scores[!(scores$SYMBOL %in% exclude), , drop = FALSE]
  dat <- dat[is.finite(dat[[score_col]]), , drop = FALSE]
  dat$SYMBOL[order(dat[[score_col]], decreasing = TRUE)]
}

path_recovery_metrics <- function(scores, positives, exclude, score_col = "score") {
  ranked <- ranked_genes(scores, score_col, exclude)
  positives <- positives[positives %in% ranked]
  n_ranked <- length(ranked)
  prev <- if (n_ranked) length(positives) / n_ranked else NA_real_
  score <- scores[[score_col]]
  au <- if (length(positives)) auprc_from_score(scores$SYMBOL, score, positives, exclude) else NA_real_
  top100 <- ranked[seq_len(min(100, length(ranked)))]
  pos_ranks <- match(positives, ranked)
  out <- data.frame(
    n_targets = length(positives),
    ranked_n_genes = n_ranked,
    prevalence = prev,
    raw_AUPRC = au,
    partial_AUPRC_top100 = if (length(positives)) {
      auprc_from_score(top100, rev(seq_len(length(top100))), positives)
    } else NA_real_,
    prevalence_normalized_AUPRC = if (is.finite(au) && is.finite(prev) && prev < 1) (au - prev) / (1 - prev) else NA_real_,
    first_target_rank = if (length(pos_ranks)) min(pos_ranks, na.rm = TRUE) else NA_real_,
    mean_reciprocal_rank = if (length(pos_ranks)) mean(1 / pos_ranks, na.rm = TRUE) else NA_real_,
    stringsAsFactors = FALSE
  )
  for (k in c(50, 100)) {
    top <- ranked[seq_len(min(k, length(ranked)))]
    recall <- if (length(positives)) mean(positives %in% top) else NA_real_
    random_expected <- if (n_ranked) min(k, n_ranked) / n_ranked else NA_real_
    out[[paste0("top", k, "_recall")]] <- recall
    out[[paste0("top", k, "_fold_enrichment_over_random")]] <- recall / random_expected
  }
  out
}

path_metric_rows <- function(scores, modules, exclude, score_col, score_name, network_label,
                             method, replicate, template_key, null = FALSE) {
  rows <- list()
  for (target_name in names(path_target_sets(modules))) {
    targets <- path_target_sets(modules)[[target_name]]
    met <- path_recovery_metrics(scores, targets, exclude, score_col)
    met$target_set <- target_name
    met$score_name <- score_name
    met$network_label <- network_label
    met$method <- method
    met$replicate <- replicate
    met$template_key <- template_key
    met$cell_type <- sub("::.*", "", template_key)
    met$null <- null
    rows[[length(rows) + 1]] <- met
  }
  do.call(rbind, rows)
}

signed_auprc <- function(symbols, signed_score, positives, direction, exclude = character()) {
  score <- if (direction == "risk") signed_score else -signed_score
  auprc_from_score(symbols, score, positives, exclude)
}

direction_metric_row <- function(scores, modules, directions, exclude, rank_col, signed_col,
                                 score_name, network_label, method, replicate, template_key,
                                 null = FALSE) {
  dat <- scores[!(scores$SYMBOL %in% exclude), , drop = FALSE]
  dat <- dat[is.finite(dat[[rank_col]]) & is.finite(dat[[signed_col]]), , drop = FALSE]
  ranked <- dat$SYMBOL[order(dat[[rank_col]], decreasing = TRUE)]
  signed <- setNames(dat[[signed_col]], dat$SYMBOL)
  targets <- modules$A2_intermediate_degree_capped
  target_dir <- directions$A2[targets]
  concordant <- function(genes) {
    g <- intersect(genes, targets)
    if (!length(g)) return(logical())
    ifelse(target_dir[g] == "risk", signed[g] > 0, signed[g] < 0)
  }
  top50 <- ranked[seq_len(min(50, length(ranked)))]
  top100 <- ranked[seq_len(min(100, length(ranked)))]
  d <- modules$D_opposite_sign_decoy
  ddir <- directions$D[d]
  dsign <- ifelse(ddir == "risk", 1, -1)
  dsel50 <- intersect(d, top50)
  dsel100 <- intersect(d, top100)
  data.frame(
    score_name = score_name,
    network_label = network_label,
    method = method,
    replicate = replicate,
    template_key = template_key,
    cell_type = sub("::.*", "", template_key),
    null = null,
    sign_concordant_top50_recall = if (length(targets)) sum(concordant(top50), na.rm = TRUE) / length(targets) else NA_real_,
    sign_concordant_top100_recall = if (length(targets)) sum(concordant(top100), na.rm = TRUE) / length(targets) else NA_real_,
    sign_concordant_precision_top50 = if (length(intersect(top50, targets))) mean(concordant(top50), na.rm = TRUE) else NA_real_,
    sign_concordant_precision_top100 = if (length(intersect(top100, targets))) mean(concordant(top100), na.rm = TRUE) else NA_real_,
    signed_AUPRC_risk = signed_auprc(dat$SYMBOL, dat[[signed_col]], modules$A2_intermediate_degree_capped_risk, "risk", character()),
    signed_AUPRC_protective = signed_auprc(dat$SYMBOL, dat[[signed_col]], modules$A2_intermediate_degree_capped_protective, "protective", character()),
    direction_accuracy_among_top50_A2 = if (length(intersect(top50, targets))) mean(concordant(top50), na.rm = TRUE) else NA_real_,
    direction_accuracy_among_top100_A2 = if (length(intersect(top100, targets))) mean(concordant(top100), na.rm = TRUE) else NA_real_,
    opposite_sign_decoy_top50_rate = if (length(d)) length(dsel50) / length(d) else NA_real_,
    opposite_sign_decoy_top100_rate = if (length(d)) length(dsel100) / length(d) else NA_real_,
    D_top50_selection_rate = if (length(d)) length(dsel50) / length(d) else NA_real_,
    D_top100_selection_rate = if (length(d)) length(dsel100) / length(d) else NA_real_,
    D_sign_discordant_selection_rate = if (length(dsel100)) {
      mean(sign(signed[dsel100]) == dsign[dsel100], na.rm = TRUE)
    } else 0,
    D_fold_enrichment_over_random = if (length(d) && length(ranked)) {
      (length(dsel100) / length(d)) / (min(100, length(ranked)) / length(ranked))
    } else NA_real_,
    stringsAsFactors = FALSE
  )
}

degree_bias_metric_row <- function(scores, rep, score_col, score_name, network_label,
                                   method, replicate, template_key, null = FALSE,
                                   threshold = 0.1) {
  dat <- scores[!(scores$SYMBOL %in% rep$modules$A1), , drop = FALSE]
  dat <- dat[is.finite(dat[[score_col]]), , drop = FALSE]
  adj <- rep$networks[[if (network_label %in% names(rep$networks)) network_label else "relevant"]]
  bin <- adj > threshold
  diag(bin) <- FALSE
  degree <- rowSums(bin)
  strength <- rowSums(adj * bin)
  deg_pct <- whole_degree_percentile(adj, threshold = threshold)
  deg_tab <- rep$degree_table
  deg_pct[deg_tab$SYMBOL] <- deg_tab$degree_percentile
  score <- setNames(dat[[score_col]], dat$SYMBOL)
  top50 <- names(sort(score, decreasing = TRUE))[seq_len(min(50, length(score)))]
  top100 <- names(sort(score, decreasing = TRUE))[seq_len(min(100, length(score)))]
  cdecoy <- rep$modules$C_high_degree_decoy
  data.frame(
    score_name = score_name,
    network_label = network_label,
    method = method,
    replicate = replicate,
    template_key = template_key,
    cell_type = sub("::.*", "", template_key),
    null = null,
    n_A2_primary = length(rep$modules$A2_intermediate_degree_capped),
    n_A2_low_degree = length(rep$modules$A2_low_degree),
    n_A2_mid_degree = length(rep$modules$A2_mid_degree),
    n_A2_high_degree = length(rep$modules$A2_high_degree),
    n_A2_extreme_high_degree = length(rep$modules$A2_extreme_high_degree),
    fraction_A2_extreme_high_degree = if (length(rep$modules$A2_intermediate_degree_capped)) length(rep$modules$A2_extreme_high_degree) / length(rep$modules$A2_intermediate_degree_capped) else NA_real_,
    median_A2_degree_percentile = stats::median(rep$degree_table$degree_percentile[match(rep$modules$A2_intermediate_degree_capped, rep$degree_table$SYMBOL)], na.rm = TRUE),
    A2_vs_background_degree_KS = rep$degree_qc$A2_vs_background_degree_KS,
    A2_vs_C_decoy_degree_KS = rep$degree_qc$A2_vs_C_decoy_degree_KS,
    score_degree_spearman = suppressWarnings(stats::cor(score[names(score) %in% names(degree)],
                                                        degree[names(score)[names(score) %in% names(degree)]],
                                                        method = "spearman", use = "complete.obs")),
    score_strength_spearman = suppressWarnings(stats::cor(score[names(score) %in% names(strength)],
                                                          strength[names(score)[names(score) %in% names(strength)]],
                                                          method = "spearman", use = "complete.obs")),
    top50_degree_percentile_median = stats::median(deg_pct[top50], na.rm = TRUE),
    top100_degree_percentile_median = stats::median(deg_pct[top100], na.rm = TRUE),
    C_high_degree_decoy_top50_rate = if (length(cdecoy)) mean(cdecoy %in% top50) else NA_real_,
    C_high_degree_decoy_top100_rate = if (length(cdecoy)) mean(cdecoy %in% top100) else NA_real_,
    stringsAsFactors = FALSE
  )
}

score_payload_rows <- function(rep, rep_id, null = FALSE, n_weight_perm = 10) {
  mean_expression <- rep$universe$mean_expression
  twas <- rep$twas
  modules <- rep$modules
  rows <- list()
  dir_rows <- list()
  degree_rows <- list()
  bench_rows <- list()
  bench_dir_rows <- list()
  bench_degree_rows <- list()

  add_primary <- function(scores, score_col, signed_col, score_name, network_label, method) {
    rows[[length(rows) + 1]] <<- path_metric_rows(scores, modules, modules$A1, score_col, score_name,
                                                  network_label, method, rep_id, rep$template_key, null)
    if (!all(is.na(scores[[signed_col]]))) {
      dir_rows[[length(dir_rows) + 1]] <<- direction_metric_row(scores, modules, rep$directions,
                                                                modules$A1, score_col, signed_col,
                                                                score_name, network_label, method,
                                                                rep_id, rep$template_key, null)
    }
    degree_rows[[length(degree_rows) + 1]] <<- degree_bias_metric_row(scores, rep, score_col, score_name,
                                                                      network_label, method,
                                                                      rep_id, rep$template_key, null)
  }

  rel <- faithful_m2_scores(rep$networks$relevant, mean_expression, twas, method = "ber_p", n.perm = 300)
  rel$raw_abs_Z <- abs(rel$TWAS.Z)
  rel$raw_signed_TWAS <- rel$TWAS.Z
  add_primary(rel, "NESTA_M2_composite", "final_NESTA_heat", "NESTA_M2_composite", "relevant", "NESTA_M2_faithful")
  add_primary(rel, "NESTA_M2_abs_delta_NESTA", "delta_NESTA", "legacy_delta_only", "relevant", "NESTA_M2_faithful")
  add_primary(rel, "NESTA_M2_abs_final_heat", "final_NESTA_heat", "heat_only", "relevant", "NESTA_M2_faithful")
  add_primary(rel, "raw_abs_Z", "raw_signed_TWAS", "raw_TWAS_abs", "relevant", "raw_TWAS")

  nd <- no_diffusion_m2_scores(rep$networks$relevant, mean_expression, twas)
  add_primary(nd, "NESTA_M2_composite", "final_NESTA_heat", "no_diffusion_M2", "relevant", "M2_no_diffusion")

  for (nm in c("I2", "I3")) {
    s <- faithful_m2_scores(rep$networks[[nm]], mean_expression, twas, method = "ber_p", n.perm = 300)
    add_primary(s, "NESTA_M2_composite", "final_NESTA_heat", "NESTA_M2_composite", nm, "NESTA_M2_faithful")
  }

  if (!null && n_weight_perm > 0) {
    for (i in seq_len(n_weight_perm)) {
      adj <- permute_weights(rep$networks$relevant, seed = 604000 + rep_id * 100 + i)
      s <- faithful_m2_scores(adj, mean_expression, twas, method = "ber_p", n.perm = 300)
      add_primary(s, "NESTA_M2_composite", "final_NESTA_heat", "NESTA_M2_composite",
                  paste0("weight_perm_", i), "NESTA_M2_faithful")
    }
  }

  if (!null) {
    for (net_nm in c("relevant", "I2", "I3")) {
      b <- benchmark_scores(rep$networks[[net_nm]], twas, include_sensitivity = TRUE)
      for (score_nm in names(b)) {
        tab <- b[[score_nm]]
        bench_rows[[length(bench_rows) + 1]] <- path_metric_rows(tab, modules, modules$A1, "score",
                                                                 score_nm, net_nm, score_nm,
                                                                 rep_id, rep$template_key, null)
        if ("score_signed" %in% names(tab) && !all(is.na(tab$score_signed))) {
          bench_dir_rows[[length(bench_dir_rows) + 1]] <- direction_metric_row(tab, modules, rep$directions,
                                                                               modules$A1, "score", "score_signed",
                                                                               score_nm, net_nm, score_nm,
                                                                               rep_id, rep$template_key, null)
        }
        bench_degree_rows[[length(bench_degree_rows) + 1]] <- degree_bias_metric_row(tab, rep, "score",
                                                                                    score_nm, net_nm, score_nm,
                                                                                    rep_id, rep$template_key, null)
      }
    }
  }

  list(primary_path = do.call(rbind, rows),
       primary_direction = do.call(rbind, dir_rows),
       primary_degree = do.call(rbind, degree_rows),
       benchmark_path = if (length(bench_rows)) do.call(rbind, bench_rows) else data.frame(),
       benchmark_signed = if (length(bench_dir_rows)) do.call(rbind, bench_dir_rows) else data.frame(),
       benchmark_degree = if (length(bench_degree_rows)) do.call(rbind, bench_degree_rows) else data.frame())
}

aggregate_weight_permutation <- function(metrics) {
  wp <- metrics[grepl("^weight_perm_", metrics$network_label), , drop = FALSE]
  if (!nrow(wp)) return(metrics)
  num <- names(wp)[vapply(wp, is.numeric, logical(1))]
  keys <- c("replicate", "template_key", "cell_type", "target_set", "null")
  form <- stats::as.formula(paste("cbind(", paste(num, collapse = ","), ") ~ ",
                                  paste(keys, collapse = "+")))
  med <- aggregate(form, wp, stats::median, na.rm = TRUE)
  med$score_name <- "median_weight_permutation"
  med$network_label <- "median_weight_permutation"
  med$method <- "weight_permutation_control"
  for (nm in setdiff(names(metrics), names(med))) med[[nm]] <- NA
  rbind(metrics[!grepl("^weight_perm_", metrics$network_label), , drop = FALSE],
        med[, names(metrics), drop = FALSE])
}

wide_metric <- function(metrics, rowspec, metric) {
  x <- metrics
  for (nm in names(rowspec)) x <- x[x[[nm]] == rowspec[[nm]], , drop = FALSE]
  out <- x[, c("replicate", metric), drop = FALSE]
  names(out)[2] <- "value"
  out
}

paired_contrast_table <- function(metrics, specs, metrics_to_compare, target_set = "A2_intermediate_degree_capped") {
  rows <- list()
  m0 <- metrics[metrics$target_set == target_set, , drop = FALSE]
  for (metric in metrics_to_compare) {
    base <- wide_metric(m0, specs$base, metric)
    names(base)[2] <- "base"
    for (nm in setdiff(names(specs), "base")) {
      cmp <- wide_metric(m0, specs[[nm]], metric)
      names(cmp)[2] <- "cmp"
      tab <- merge(base, cmp, by = "replicate", all = FALSE)
      diff <- tab$base - tab$cmp
      ci <- paired_bootstrap_ci(diff)
      rows[[length(rows) + 1]] <- data.frame(
        target_set = target_set,
        contrast = nm,
        metric = metric,
        mean = ci["mean"],
        median = ci["median"],
        ci_low = ci["lo"],
        ci_high = ci["hi"],
        improved_fraction = mean(diff > 0, na.rm = TRUE),
        unchanged_fraction = mean(diff == 0, na.rm = TRUE),
        deteriorated_fraction = mean(diff < 0, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

direction_contrast_table <- function(metrics) {
  specs <- list(
    base = list(score_name = "NESTA_M2_composite", network_label = "relevant", method = "NESTA_M2_faithful"),
    raw_signed_TWAS = list(score_name = "raw_TWAS_abs", network_label = "relevant", method = "raw_TWAS"),
    unsigned_RWR_abs_prior = list(score_name = "RWR_abs_prior", network_label = "relevant", method = "RWR_abs_prior"),
    signed_two_channel_RWR = list(score_name = "RWR_signed_two_channel", network_label = "relevant", method = "RWR_signed_two_channel"),
    unsigned_PPR_abs_prior = list(score_name = "PPR_abs_prior", network_label = "relevant", method = "PPR_abs_prior"),
    signed_two_channel_PPR = list(score_name = "PPR_signed_two_channel", network_label = "relevant", method = "PPR_signed_two_channel")
  )
  paired_contrast_table(transform(metrics, target_set = "A2_intermediate_degree_capped"), specs,
                        c("sign_concordant_top100_recall", "direction_accuracy_among_top100_A2",
                          "opposite_sign_decoy_top100_rate"), "A2_intermediate_degree_capped")
}

degree_contrast_table <- function(metrics) {
  specs <- list(
    base = list(score_name = "NESTA_M2_composite", network_label = "relevant", method = "NESTA_M2_faithful"),
    RWR_abs_prior = list(score_name = "RWR_abs_prior", network_label = "relevant", method = "RWR_abs_prior"),
    PPR_abs_prior = list(score_name = "PPR_abs_prior", network_label = "relevant", method = "PPR_abs_prior"),
    RWR_signed_two_channel = list(score_name = "RWR_signed_two_channel", network_label = "relevant", method = "RWR_signed_two_channel"),
    PPR_signed_two_channel = list(score_name = "PPR_signed_two_channel", network_label = "relevant", method = "PPR_signed_two_channel")
  )
  paired_contrast_table(transform(metrics, target_set = "A2_intermediate_degree_capped"), specs,
                        c("C_high_degree_decoy_top100_rate", "score_degree_spearman",
                          "score_strength_spearman", "top100_degree_percentile_median"),
                        "A2_intermediate_degree_capped")
}

summarize_primary_outputs <- function(primary_path, primary_direction, primary_degree, out_dir) {
  safe_dir_create(out_dir)
  primary_path <- aggregate_weight_permutation(primary_path)
  atomic_write_csv(primary_path, file.path(out_dir, "PRIMARY_PATH_STRATIFIED_METRICS.csv"))
  atomic_write_csv(primary_direction, file.path(out_dir, "PRIMARY_DIRECTION_AWARE_METRICS.csv"))
  atomic_write_csv(primary_degree, file.path(out_dir, "PRIMARY_DEGREE_BIAS_METRICS.csv"))
  specs <- list(
    base = list(score_name = "NESTA_M2_composite", network_label = "relevant", method = "NESTA_M2_faithful"),
    raw_TWAS_abs = list(score_name = "raw_TWAS_abs", network_label = "relevant", method = "raw_TWAS"),
    legacy_delta_only = list(score_name = "legacy_delta_only", network_label = "relevant", method = "NESTA_M2_faithful"),
    no_diffusion_M2 = list(score_name = "no_diffusion_M2", network_label = "relevant", method = "M2_no_diffusion"),
    I2_module_disrupted = list(score_name = "NESTA_M2_composite", network_label = "I2", method = "NESTA_M2_faithful"),
    I3_expression_matched_randomized = list(score_name = "NESTA_M2_composite", network_label = "I3", method = "NESTA_M2_faithful"),
    median_weight_permutation = list(score_name = "median_weight_permutation", network_label = "median_weight_permutation", method = "weight_permutation_control"),
    heat_only_relevant = list(score_name = "heat_only", network_label = "relevant", method = "NESTA_M2_faithful")
  )
  contrasts <- paired_contrast_table(primary_path, specs,
                                     c("top50_recall", "top100_recall",
                                       "top50_fold_enrichment_over_random",
                                       "top100_fold_enrichment_over_random",
                                       "raw_AUPRC", "partial_AUPRC_top100",
                                       "prevalence_normalized_AUPRC"))
  atomic_write_csv(contrasts, file.path(out_dir, "PRIMARY_PATH_STRATIFIED_CONTRASTS.csv"))
  d_contrasts <- direction_contrast_table(rbind(primary_direction, data.frame()))
  atomic_write_csv(d_contrasts, file.path(out_dir, "PRIMARY_DIRECTION_AWARE_CONTRASTS.csv"))
  degree_contrasts <- degree_contrast_table(primary_degree)
  atomic_write_csv(degree_contrasts, file.path(out_dir, "PRIMARY_DEGREE_BIAS_CONTRASTS.csv"))
  list(path_metrics = primary_path, path_contrasts = contrasts,
       direction_metrics = primary_direction, direction_contrasts = d_contrasts,
       degree_metrics = primary_degree, degree_contrasts = degree_contrasts)
}

summarize_benchmark_outputs <- function(benchmark_path, benchmark_signed, benchmark_degree,
                                        primary_path, primary_direction, primary_degree, out_dir) {
  atomic_write_csv(benchmark_path, file.path(out_dir, "BENCHMARK_PATH_STRATIFIED_METRICS.csv"))
  atomic_write_csv(benchmark_signed, file.path(out_dir, "BENCHMARK_SIGNED_METRICS.csv"))
  atomic_write_csv(benchmark_degree, file.path(out_dir, "BENCHMARK_DEGREE_BIAS_METRICS.csv"))
  nesta <- primary_path[primary_path$score_name == "NESTA_M2_composite" &
                          primary_path$network_label == "relevant" &
                          primary_path$method == "NESTA_M2_faithful", , drop = FALSE]
  bench_rel <- benchmark_path[benchmark_path$network_label == "relevant" &
                                benchmark_path$score_name %in% c("RWR_abs_prior", "PPR_abs_prior",
                                                                  "RWR_signed_two_channel", "PPR_signed_two_channel",
                                                                  "NESTA_common_prior"), , drop = FALSE]
  combo <- rbind(nesta[, names(bench_rel), drop = FALSE], bench_rel)
  specs <- list(
    base = list(score_name = "NESTA_M2_composite", network_label = "relevant"),
    RWR_abs_prior = list(score_name = "RWR_abs_prior", network_label = "relevant"),
    PPR_abs_prior = list(score_name = "PPR_abs_prior", network_label = "relevant"),
    RWR_signed_two_channel = list(score_name = "RWR_signed_two_channel", network_label = "relevant"),
    PPR_signed_two_channel = list(score_name = "PPR_signed_two_channel", network_label = "relevant"),
    NESTA_common_prior = list(score_name = "NESTA_common_prior", network_label = "relevant")
  )
  b_contrasts <- paired_contrast_table(combo, specs,
                                       c("top50_recall", "top100_recall",
                                         "raw_AUPRC", "prevalence_normalized_AUPRC"))
  atomic_write_csv(b_contrasts, file.path(out_dir, "BENCHMARK_PATH_STRATIFIED_CONTRASTS.csv"))
  s_contrasts <- if (nrow(benchmark_signed)) direction_contrast_table(rbind(
    primary_direction[, names(benchmark_signed), drop = FALSE], benchmark_signed
  )) else data.frame()
  atomic_write_csv(s_contrasts, file.path(out_dir, "BENCHMARK_SIGNED_CONTRASTS.csv"))
  degree_contrasts <- if (nrow(benchmark_degree)) degree_contrast_table(rbind(
    primary_degree[, names(benchmark_degree), drop = FALSE], benchmark_degree
  )) else data.frame()
  atomic_write_csv(degree_contrasts, file.path(out_dir, "BENCHMARK_DEGREE_BIAS_CONTRASTS.csv"))
  list(path_contrasts = b_contrasts, signed_contrasts = s_contrasts,
       degree_contrasts = degree_contrasts)
}

directional_qc_row <- function(rep) {
  dirs <- rep$directions$A2
  risk <- rep$modules$A2_risk
  prot <- rep$modules$A2_protective
  data.frame(
    replicate = rep$rep_id,
    template_key = rep$template_key,
    n_A1_risk = length(rep$modules$A1_risk),
    n_A1_protective = length(rep$modules$A1_protective),
    n_A2_risk = length(risk),
    n_A2_protective = length(prot),
    n_A2_intermediate_degree_capped_risk = length(rep$modules$A2_intermediate_degree_capped_risk),
    n_A2_intermediate_degree_capped_protective = length(rep$modules$A2_intermediate_degree_capped_protective),
    n_D_opposite_sign_decoy = length(rep$modules$D_opposite_sign_decoy),
    extreme_high_degree_A2_count = rep$degree_qc$extreme_high_degree_A2_count,
    fraction_A2_extreme_high_degree = rep$degree_qc$fraction_A2_extreme_high_degree,
    median_A2_degree_percentile = rep$degree_qc$median_A2_degree_percentile,
    same_direction_connectivity_fraction = if (length(rep$modules$A2_intermediate_degree_capped)) {
      mean(rep$modules$A2_intermediate_degree_capped %in% names(dirs))
    } else 0,
    opposite_direction_shortcut_fraction = NA_real_,
    directional_qc_pass = length(rep$modules$A1_risk) >= 2 &&
      length(rep$modules$A1_protective) >= 2 &&
      length(rep$modules$A2_intermediate_degree_capped_risk) > 0 &&
      length(rep$modules$A2_intermediate_degree_capped_protective) > 0 &&
      length(rep$modules$D_opposite_sign_decoy) >= 30,
    stringsAsFactors = FALSE
  )
}

run_topology_qc <- function(lib, n_reps = 40, out_dir = project_file("results/topology_qc"),
                            seed_base = 2026062900) {
  safe_dir_create(out_dir)
  templates <- paste(lib$templates$cell_type, lib$templates$module, sep = "::")
  if (length(templates) < 20 || length(unique(lib$templates$cell_type)) < 3) {
    reason <- "insufficient_empirical_templates"
    atomic_write_csv(data.frame(pilot_go = FALSE, reason = reason), file.path(out_dir, "topology_qc_decision.csv"))
    return(list(pass = FALSE, reason = reason))
  }
  reps <- vector("list", n_reps)
  rows <- list()
  match_rows <- list()
  path_rows <- list()
  direction_rows <- list()
  degree_rows <- list()
  for (i in seq_len(n_reps)) {
    key <- templates[((i - 1) %% length(templates)) + 1]
    rep <- make_replicate(lib, key, i, seed = seed_base + i)
    reps[[i]] <- rep
    path_rows[[i]] <- transform(rep$path_qc, replicate = i, template_key = key,
                                n_A2_all = length(rep$modules$A2_all),
                                n_A2_intermediate_degree_capped = length(rep$modules$A2_intermediate_degree_capped),
                                n_thinned_edges = length(rep$thinning_log))
    direction_rows[[i]] <- directional_qc_row(rep)
    degree_rows[[i]] <- transform(rep$degree_qc, replicate = i, template_key = key)
    for (nm in names(rep$networks)) {
      m <- network_metrics(rep$networks[[nm]], rep$modules, rep$universe$mean_expression)
      m$replicate <- i
      m$template_key <- key
      m$network_label <- nm
      m$hard_global_pass <- with(m,
        density >= 0.0008 & density <= 0.22 &
          median_degree > 0 & median_strength > 0 &
          largest_connected_component_fraction >= 0.04 &
          module_a_connected_component_fraction >= 0.75 &
          isolated_fraction <= 0.97 &
          finite_edge_weights & !duplicated_gene_names)
      m$module_a_local_pass <- with(m,
        length(rep$modules$A) >= 30 & length(rep$modules$A) <= 160 &
          length(rep$modules$A2_intermediate_degree_capped) >= 12 &
          within_between_a_ratio >= 1.40 &
          is.finite(median_a1_a2_path) & median_a1_a2_path >= 2 & median_a1_a2_path <= 4 &
          a1_a2_same_component_fraction >= 0.75 &
          direct_1hop_fraction <= 0.60 &
          two_hop_fraction >= 0.40 &
          clustering_a > background_clustering_p40)
      rows[[length(rows) + 1]] <- m
    }
    for (nm in c("I1", "I2", "I3")) {
      q <- matching_qc(rep$networks$relevant, rep$networks[[nm]], rep$universe$mean_expression, rep$modules)
      q$replicate <- i
      q$template_key <- key
      q$network_label <- nm
      match_rows[[length(match_rows) + 1]] <- q
    }
  }
  qc <- do.call(rbind, rows)
  matching <- do.call(rbind, match_rows)
  path_qc <- do.call(rbind, path_rows)
  direction_qc <- do.call(rbind, direction_rows)
  degree_qc <- do.call(rbind, degree_rows)
  matching$strict_matched_control <- matching$network_label %in% c("I2", "I3")
  matching$matching_pass <- with(matching,
    ifelse(strict_matched_control,
           density_diff <= 0.025 &
             median_strength_ratio >= 0.45 &
             median_strength_ratio <= 2.20 &
             expression_ks <= 0.22 &
             identical_twas &
             module_disrupted,
           TRUE))
  relevant <- qc$network_label == "relevant"
  matched <- qc$network_label %in% c("relevant", "I2", "I3")
  hard_global_pass_fraction <- mean(qc$hard_global_pass[matched])
  module_a_local_pass_fraction <- mean(qc$module_a_local_pass[relevant])
  i2_i3_matching_pass_fraction <- mean(matching$matching_pass[matching$strict_matched_control])
  topology_pass <- all(qc$hard_global_pass[matched]) &&
    module_a_local_pass_fraction >= 0.95 &&
    i2_i3_matching_pass_fraction >= 0.85
  path_pass <- all(path_qc$path_stratification_pass)
  directional_pass <- all(direction_qc$directional_qc_pass)
  degree_pass <- all(degree_qc$degree_distribution_qc_pass)
  pass <- topology_pass && path_pass && directional_pass && degree_pass
  reason <- if (pass) "topology_path_degree_direction_qc_passed" else "topology_path_degree_direction_qc_failed"
  atomic_write_csv(qc, file.path(out_dir, "topology_qc_metrics.csv"))
  atomic_write_csv(matching, file.path(out_dir, "relevant_irrelevant_matching_qc.csv"))
  atomic_write_csv(path_qc, file.path(out_dir, "PATH_STRATIFICATION_AUDIT.csv"))
  atomic_write_csv(direction_qc, file.path(out_dir, "DIRECTIONAL_SIGNAL_AUDIT.csv"))
  atomic_write_csv(degree_qc, file.path(out_dir, "DEGREE_DISTRIBUTION_AUDIT.csv"))
  atomic_write_csv(data.frame(pilot_go = pass, reason = reason,
                              topology_qc_pass = topology_pass,
                              path_stratification_qc_pass = path_pass,
                              degree_distribution_qc_pass = degree_pass,
                              directional_qc_pass = directional_pass,
                              template_count = length(templates),
                              template_cell_types = length(unique(lib$templates$cell_type)),
                              hard_global_pass_fraction = hard_global_pass_fraction,
                              module_a_local_pass_fraction = module_a_local_pass_fraction,
                              i2_i3_matching_pass_fraction = i2_i3_matching_pass_fraction,
                              path_fallback_used_fraction = mean(path_qc$path_fallback_used),
                              i1_reported_nonblocking = TRUE),
                   file.path(out_dir, "topology_qc_decision.csv"))
  atomic_save_rds(reps, file.path(out_dir, "frozen_replicate_designs.rds"))
  list(pass = pass, reason = reason, topology_qc_pass = topology_pass,
       path_qc_pass = path_pass, degree_qc_pass = degree_pass,
       directional_qc_pass = directional_pass, reps = reps)
}

run_pilot <- function(lib, reps, out_dir = project_file("results/pilot"),
                      null_dir = project_file("results/null"),
                      cores = max(1, min(16, parallel::detectCores(logical = FALSE))),
                      null_seed_base = 2026063000) {
  safe_dir_create(out_dir)
  safe_dir_create(null_dir)
  run_one <- function(i) {
    path <- file.path(out_dir, sprintf("pilot_replicate_%03d_metrics.rds", i))
    if (file.exists(path)) return(readRDS(path))
    x <- score_payload_rows(reps[[i]], i, null = FALSE, n_weight_perm = 10)
    atomic_save_rds(x, path)
    x
  }
  rows <- parallel::mclapply(seq_along(reps), run_one, mc.cores = cores, mc.preschedule = FALSE)
  primary_path <- do.call(rbind, lapply(rows, `[[`, "primary_path"))
  primary_direction <- do.call(rbind, lapply(rows, `[[`, "primary_direction"))
  primary_degree <- do.call(rbind, lapply(rows, `[[`, "primary_degree"))
  benchmark_path <- do.call(rbind, lapply(rows, `[[`, "benchmark_path"))
  benchmark_signed <- do.call(rbind, lapply(rows, `[[`, "benchmark_signed"))
  benchmark_degree <- do.call(rbind, lapply(rows, `[[`, "benchmark_degree"))
  primary <- summarize_primary_outputs(primary_path, primary_direction, primary_degree, out_dir)
  bench <- summarize_benchmark_outputs(benchmark_path, benchmark_signed, benchmark_degree,
                                       primary$path_metrics, primary$direction_metrics,
                                       primary$degree_metrics, out_dir)

  null_reps <- lapply(seq_along(reps), function(i) make_replicate(lib, reps[[i]]$template_key, i,
                                                                  seed = null_seed_base + i, null = TRUE))
  run_null_one <- function(i) {
    path <- file.path(null_dir, sprintf("null_replicate_%03d_metrics.rds", i))
    if (file.exists(path)) return(readRDS(path))
    x <- score_payload_rows(null_reps[[i]], i, null = TRUE, n_weight_perm = 0)
    atomic_save_rds(x, path)
    x
  }
  null_rows <- parallel::mclapply(seq_along(null_reps), run_null_one, mc.cores = cores, mc.preschedule = FALSE)
  null_primary <- do.call(rbind, lapply(null_rows, `[[`, "primary_path"))
  null_direction <- do.call(rbind, lapply(null_rows, `[[`, "primary_direction"))
  guard <- summarize_null_guardrails(null_reps, null_primary, null_direction, null_dir)
  list(path_metrics = primary$path_metrics,
       path_contrasts = primary$path_contrasts,
       direction_metrics = primary$direction_metrics,
       direction_contrasts = primary$direction_contrasts,
       degree_metrics = primary$degree_metrics,
       degree_contrasts = primary$degree_contrasts,
       benchmark_path_metrics = benchmark_path,
       benchmark_path_contrasts = bench$path_contrasts,
       benchmark_signed_metrics = benchmark_signed,
       benchmark_signed_contrasts = bench$signed_contrasts,
       benchmark_degree_metrics = benchmark_degree,
       benchmark_degree_contrasts = bench$degree_contrasts,
       null_guardrails = guard)
}

summarize_null_guardrails <- function(reps, null_primary, null_direction, out_dir) {
  rows <- list()
  for (i in seq_along(reps)) {
    rep <- reps[[i]]
    m2 <- faithful_m2_scores(rep$networks$relevant, rep$universe$mean_expression, rep$twas)
    ranked <- m2$SYMBOL[order(m2$NESTA_M2_composite, decreasing = TRUE)]
    top50 <- ranked[seq_len(min(50, length(ranked)))]
    top100 <- ranked[seq_len(min(100, length(ranked)))]
    prev_b <- length(rep$modules$B) / (length(rep$universe$genes) - length(rep$modules$A1))
    prev_c <- length(rep$modules$C) / (length(rep$universe$genes) - length(rep$modules$A1))
    bridge <- rep$modules$A_bridge
    d <- rep$modules$D_opposite_sign_decoy
    rows[[length(rows) + 1]] <- data.frame(
      replicate = i,
      module_b_representation_ratio = mean(rep$modules$B %in% top50) / prev_b,
      module_c_representation_ratio = mean(rep$modules$C %in% top50) / prev_c,
      proximal_like_background_top100_rate = mean(setdiff(rep$modules$A2_proximal, rep$modules$A2_intermediate_degree_capped) %in% top100),
      bridge_like_background_top100_rate = if (length(bridge)) mean(bridge %in% top100) else 0,
      opposite_sign_decoy_top100_rate = if (length(d)) mean(d %in% top100) else 0,
      C_high_degree_decoy_top100_rate = if (length(rep$modules$C_high_degree_decoy)) mean(rep$modules$C_high_degree_decoy %in% top100) else 0,
      score_degree_spearman = degree_bias_metric_row(m2, rep, "NESTA_M2_composite", "NESTA_M2_composite",
                                                     "relevant", "NESTA_M2_faithful", i, rep$template_key, TRUE)$score_degree_spearman,
      sign_discordant_top100_false_positive_rate = NA_real_,
      max_composite = max(m2$NESTA_M2_composite),
      stringsAsFactors = FALSE
    )
  }
  tab <- do.call(rbind, rows)
  guard <- data.frame(
    mean_b_ratio = mean(tab$module_b_representation_ratio, na.rm = TRUE),
    mean_c_ratio = mean(tab$module_c_representation_ratio, na.rm = TRUE),
    mean_bridge_like_top100_rate = mean(tab$bridge_like_background_top100_rate, na.rm = TRUE),
    mean_opposite_sign_decoy_top100_rate = mean(tab$opposite_sign_decoy_top100_rate, na.rm = TRUE),
    mean_C_high_degree_decoy_top100_rate = mean(tab$C_high_degree_decoy_top100_rate, na.rm = TRUE),
    mean_score_degree_spearman = mean(tab$score_degree_spearman, na.rm = TRUE),
    empirical_fwe = mean(tab$max_composite > stats::quantile(tab$max_composite, 0.95), na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  guard$guardrail_pass <- with(guard, mean_b_ratio <= 2 & mean_c_ratio <= 2 &
                                mean_bridge_like_top100_rate <= 0.35 &
                                mean_opposite_sign_decoy_top100_rate <= 0.35 &
                                mean_C_high_degree_decoy_top100_rate <= 0.35 &
                                mean_score_degree_spearman <= 0.50 &
                                empirical_fwe <= 0.075)
  atomic_write_csv(tab, file.path(out_dir, "null_bias_replicate_metrics.csv"))
  atomic_write_csv(guard, file.path(out_dir, "NULL_BIAS_GUARDRAILS.csv"))
  guard
}

go_decision <- function(topology, pilot = NULL) {
  if (!isTRUE(topology$pass)) {
    return(list(status = "STOPPED", reason = topology$reason,
                pilot_started = "NO", confirmatory_started = "NO", stop_go = "STOP"))
  }
  if (is.null(pilot)) {
    return(list(status = "STOPPED", reason = "pilot_not_run",
                pilot_started = "NO", confirmatory_started = "NO", stop_go = "STOP"))
  }
  pc <- pilot$path_contrasts
  top100 <- pc[pc$target_set == "A2_intermediate_degree_capped" & pc$metric == "top100_recall", , drop = FALSE]
  required <- c("raw_TWAS_abs", "legacy_delta_only", "no_diffusion_M2",
                "I2_module_disrupted", "I3_expression_matched_randomized",
                "median_weight_permutation")
  idx <- match(required, top100$contrast)
  ok_primary <- all(!is.na(idx)) && all(top100$mean[idx] > 0, na.rm = TRUE)
  ok_raw_ci <- top100$ci_low[top100$contrast == "raw_TWAS_abs"] > 0
  ok_delta_ci <- top100$ci_low[top100$contrast == "legacy_delta_only"] > 0
  ok_network_ci <- all(top100$ci_low[top100$contrast %in% c("no_diffusion_M2", "I2_module_disrupted",
                                                            "I3_expression_matched_randomized",
                                                            "median_weight_permutation")] >= -0.02)
  dm <- pilot$direction_metrics[pilot$direction_metrics$score_name == "NESTA_M2_composite" &
                                  pilot$direction_metrics$network_label == "relevant", , drop = FALSE]
  ok_direction <- mean(dm$sign_concordant_top100_recall, na.rm = TRUE) > 0 &&
    mean(dm$direction_accuracy_among_top100_A2, na.rm = TRUE) > 0.5 &&
    mean(dm$opposite_sign_decoy_top100_rate, na.rm = TRUE) <= 0.35
  ok_null <- isTRUE(pilot$null_guardrails$guardrail_pass[1])
  deg <- pilot$degree_metrics[pilot$degree_metrics$score_name == "NESTA_M2_composite" &
                                pilot$degree_metrics$network_label == "relevant", , drop = FALSE]
  ok_degree <- nrow(deg) > 0 &&
    all(deg$n_A2_extreme_high_degree <= pmin(7, ceiling(0.25 * pmax(deg$n_A2_primary, 1))), na.rm = TRUE) &&
    mean(deg$C_high_degree_decoy_top100_rate, na.rm = TRUE) <= 0.35 &&
    abs(mean(deg$score_degree_spearman, na.rm = TRUE)) <= 0.20
  bc <- pilot$benchmark_path_contrasts
  bench <- bc[bc$target_set == "A2_intermediate_degree_capped" & bc$metric == "top100_recall" &
                bc$contrast %in% c("RWR_abs_prior", "PPR_abs_prior",
                                   "RWR_signed_two_channel", "PPR_signed_two_channel"), , drop = FALSE]
  rwr_ppr_dominate <- nrow(bench) >= 4 && all(bench$mean < 0, na.rm = TRUE)
  completion_fraction <- length(unique(pilot$path_metrics$replicate)) / 40
  audit_pass <- completion_fraction >= 0.95 && isTRUE(topology$topology_qc_pass) &&
    isTRUE(topology$path_qc_pass) && isTRUE(topology$degree_qc_pass) &&
    isTRUE(topology$directional_qc_pass)
  if (audit_pass && ok_null && ok_primary && isTRUE(ok_raw_ci) && isTRUE(ok_delta_ci) && ok_network_ci &&
      ok_degree && ok_direction && !rwr_ppr_dominate) {
    list(status = "GO", reason = "pilot_go_criteria_passed_confirmatory_required",
         pilot_started = "YES", confirmatory_started = "NO", stop_go = "GO")
  } else {
    list(status = "STOPPED", reason = "pilot_go_criteria_failed",
         pilot_started = "YES", confirmatory_started = "NO", stop_go = "STOP")
  }
}
