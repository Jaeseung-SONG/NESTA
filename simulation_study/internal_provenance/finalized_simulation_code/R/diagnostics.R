source(file.path("/home/js/NESTA/simulation", "R", "utils.R"))
source(project_file("R/fidelity.R"))
source(project_file("R/generator.R"))
source(project_file("R/pipeline.R"))

diagnostic_variants <- c(
  "D0_faithful_M2",
  "D1_M1",
  "D2_rank_expression_weighting",
  "D3_centered_expression_weighting",
  "D4_clipped_expression_weighting",
  "D5_shuffled_expression_negative_control"
)

score_definitions <- c(
  "twas_abs_z",
  "initial_weight_abs",
  "final_heat_abs",
  "delta_twas_abs",
  "delta_initial_abs",
  "signed_final_heat",
  "signed_delta_twas",
  "submitted_union"
)

completed_pilot_dir <- function() {
  file.path(report_root, "NESTA_simulation_230626_0702")
}

previous_frozen_reps_path <- function() {
  candidates <- c(
    file.path(read_report_dir(), "previous_workspace_export", "topology_qc", "frozen_replicate_designs.rds"),
    project_file("results/topology_qc/frozen_replicate_designs.rds")
  )
  candidates[file.exists(candidates)][1]
}

copy_required_diag_outputs <- function() {
  report_dir <- read_report_dir()
  for (nm in c("M2_DIAGNOSTIC_REPORT.md", "M2_SCORE_DECOMPOSITION.csv",
               "M2_WEIGHTING_DECOMPOSITION.csv", "TEMPLATE_OUTLIER_AUDIT.csv")) {
    src <- project_file("results/diagnostics", nm)
    if (file.exists(src)) {
      dst <- file.path(report_dir, nm)
      if (file.exists(dst)) unlink(dst)
      file.copy(src, dst, copy.mode = TRUE, copy.date = TRUE)
    }
  }
}

score_metrics_row <- function(scores, positives, exclude, score_name) {
  m <- rank_metrics(score_view(scores, score_name), positives, exclude)
  m$score_definition <- score_name
  m
}

module_summary_rows <- function(rep, scores, rep_id) {
  genes <- rep$universe$genes
  expr <- rep$universe$mean_expression
  expr_pct <- rank(expr, ties.method = "average") / length(expr)
  names(expr_pct) <- genes
  adj <- rep$networks$relevant
  bin <- adj > 0.1
  diag(bin) <- FALSE
  degree <- rowSums(bin)
  strength <- rowSums(adj * bin)
  groups <- list(A = rep$modules$A, A1 = rep$modules$A1, A2 = rep$modules$A2,
                 B = rep$modules$B, C = rep$modules$C,
                 background = rep$modules$background)
  rows <- lapply(names(groups), function(m) {
    g <- groups[[m]]
    data.frame(
      replicate = rep_id,
      template_key = rep$template_key,
      cell_type = sub("::.*", "", rep$template_key),
      diagnostic_table = "module_summary",
      group = m,
      n_genes = length(g),
      expression_mean = mean(expr[g]),
      expression_median = stats::median(expr[g]),
      expression_p25 = unname(stats::quantile(expr[g], 0.25)),
      expression_p75 = unname(stats::quantile(expr[g], 0.75)),
      expression_percentile_median = stats::median(expr_pct[g]),
      degree_median = stats::median(degree[g]),
      strength_median = stats::median(strength[g]),
      initial_weight_median = stats::median(scores$initial_weight[match(g, scores$SYMBOL)]),
      final_heat_median = stats::median(scores$final_heat[match(g, scores$SYMBOL)]),
      delta_twas_median = stats::median(scores$delta_twas[match(g, scores$SYMBOL)]),
      delta_initial_median = stats::median(scores$delta_initial[match(g, scores$SYMBOL)]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

mechanism_rows <- function(rep, scores, rep_id) {
  expr <- rep$universe$mean_expression
  adj <- rep$networks$relevant
  bin <- adj > 0.1
  diag(bin) <- FALSE
  degree <- rowSums(bin)
  strength <- rowSums(adj * bin)
  expr <- expr[match(scores$SYMBOL, names(expr))]
  degree <- degree[match(scores$SYMBOL, names(degree))]
  strength <- strength[match(scores$SYMBOL, names(strength))]
  vars <- data.frame(
    expression = as.numeric(expr),
    degree = as.numeric(degree),
    strength = as.numeric(strength),
    initial_weight = scores$initial_weight,
    final_heat = scores$final_heat,
    delta_twas = scores$delta_twas,
    delta_initial = scores$delta_initial
  )
  pairs <- t(utils::combn(names(vars), 2))
  cor_rows <- apply(pairs, 1, function(p) {
    data.frame(
      replicate = rep_id,
      template_key = rep$template_key,
      cell_type = sub("::.*", "", rep$template_key),
      diagnostic_table = "correlation",
      variable_x = p[1],
      variable_y = p[2],
      spearman = suppressWarnings(stats::cor(vars[[p[1]]], vars[[p[2]]], method = "spearman")),
      stringsAsFactors = FALSE
    )
  })
  ranks <- data.frame(
    SYMBOL = scores$SYMBOL,
    twas_rank = rank(-abs(scores$TWAS.Z), ties.method = "average"),
    initial_rank = rank(-abs(scores$initial_weight), ties.method = "average"),
    final_rank = rank(-abs(scores$final_heat), ties.method = "average")
  )
  a2 <- ranks[ranks$SYMBOL %in% rep$modules$A2, , drop = FALSE]
  rank_row <- data.frame(
    replicate = rep_id,
    template_key = rep$template_key,
    cell_type = sub("::.*", "", rep$template_key),
    diagnostic_table = "a2_rank_change",
    a2_initial_up_fraction = mean(a2$initial_rank < a2$twas_rank),
    a2_initial_down_fraction = mean(a2$initial_rank > a2$twas_rank),
    a2_final_up_fraction = mean(a2$final_rank < a2$initial_rank),
    a2_final_down_fraction = mean(a2$final_rank > a2$initial_rank),
    a2_median_twas_to_initial_rank_change = stats::median(a2$initial_rank - a2$twas_rank),
    a2_median_initial_to_final_rank_change = stats::median(a2$final_rank - a2$initial_rank),
    stringsAsFactors = FALSE
  )
  net <- network_metrics(adj, rep$modules, rep$universe$mean_expression)
  topo_row <- data.frame(
    replicate = rep_id,
    template_key = rep$template_key,
    cell_type = sub("::.*", "", rep$template_key),
    diagnostic_table = "topology",
    median_a1_a2_path = net$median_a1_a2_path,
    a1_a2_same_component_fraction = net$a1_a2_same_component_fraction,
    within_a_tom = net$within_a_tom,
    within_between_a_ratio = net$within_between_a_ratio,
    clustering_a = net$clustering_a,
    background_clustering_p40 = net$background_clustering_p40,
    stringsAsFactors = FALSE
  )
  list(module = module_summary_rows(rep, scores, rep_id),
       cor = do.call(rbind, cor_rows),
       rank = rank_row,
       topo = topo_row)
}

run_m2_diagnostics <- function(out_dir = project_file("results/diagnostics"),
                               cores = max(1, min(8, parallel::detectCores(logical = FALSE)))) {
  safe_dir_create(out_dir)
  reps_path <- previous_frozen_reps_path()
  if (length(reps_path) == 0 || is.na(reps_path)) {
    stop("Previous frozen replicate designs are unavailable; cannot reconstruct mandatory diagnostics")
  }
  reps <- readRDS(reps_path)
  if (length(reps) != 40) stop("Expected 40 frozen previous-pilot replicates, found ", length(reps))

  rep_one <- function(i) {
    rep <- reps[[i]]
    recon <- run_single_replicate(rep, i, n_rewire = 10, n_weight_perm = 10, null = FALSE)
    d0 <- diagnostic_weighting_scores(rep$networks$relevant, rep$universe$mean_expression,
                                      rep$twas, "D0_faithful_M2", shuffle_seed = 505000 + i)
    score_rows <- lapply(score_definitions, function(sc) {
      m <- score_metrics_row(d0, rep$modules$A2, rep$modules$A1, sc)
      m$replicate <- i
      m$template_key <- rep$template_key
      m$cell_type <- sub("::.*", "", rep$template_key)
      m
    })
    weighting_rows <- lapply(diagnostic_variants, function(v) {
      s <- diagnostic_weighting_scores(rep$networks$relevant, rep$universe$mean_expression,
                                       rep$twas, v, shuffle_seed = 606000 + i)
      m <- rank_metrics(s, rep$modules$A2, rep$modules$A1)
      m$replicate <- i
      m$template_key <- rep$template_key
      m$cell_type <- sub("::.*", "", rep$template_key)
      m$diagnostic_method <- v
      m$submitted_nesta <- v == "D0_faithful_M2"
      m$diagnostic_proposed_variant <- v %in% diagnostic_variants[3:6]
      m
    })
    mech <- mechanism_rows(rep, d0, i)
    list(recon = recon,
         scores = do.call(rbind, score_rows),
         weights = do.call(rbind, weighting_rows),
         module = mech$module,
         cor = mech$cor,
         rank = mech$rank,
         topo = mech$topo)
  }

  rows <- parallel::mclapply(seq_along(reps), rep_one, mc.cores = cores, mc.preschedule = FALSE)
  recon_metrics <- do.call(rbind, lapply(rows, `[[`, "recon"))
  score_decomp <- do.call(rbind, lapply(rows, `[[`, "scores"))
  weighting_decomp <- do.call(rbind, lapply(rows, `[[`, "weights"))
  module_diag <- do.call(rbind, lapply(rows, `[[`, "module"))
  cor_diag <- do.call(rbind, lapply(rows, `[[`, "cor"))
  rank_diag <- do.call(rbind, lapply(rows, `[[`, "rank"))
  topo_diag <- do.call(rbind, lapply(rows, `[[`, "topo"))

  recon_dir <- file.path(out_dir, "reconstruction")
  safe_dir_create(recon_dir)
  for (f in c("previous_reconstructed_metrics.csv", "previous_reconstructed_contrasts.csv",
              "primary_contrast_reproduction.csv", "mechanism_module_summary.csv",
              "mechanism_correlations.csv", "mechanism_a2_rank_changes.csv",
              "mechanism_topology.csv")) {
    p <- file.path(recon_dir, f)
    if (file.exists(p)) unlink(p)
  }
  atomic_write_csv(recon_metrics, file.path(recon_dir, "previous_reconstructed_metrics.csv"))
  contrasts <- summarize_pilot(recon_metrics, recon_dir)
  reported_path <- file.path(completed_pilot_dir(), "summary_tables", "pilot_primary_contrasts.csv")
  if (!file.exists(reported_path)) stop("Reported completed pilot contrasts are missing: ", reported_path)
  reported <- read.csv(reported_path)
  cmp <- merge(reported, contrasts, by = "contrast", suffixes = c("_reported", "_reconstructed"))
  num <- c("mean", "median", "ci_low", "ci_high", "improved_fraction")
  for (nm in num) cmp[[paste0(nm, "_abs_diff")]] <- abs(cmp[[paste0(nm, "_reported")]] - cmp[[paste0(nm, "_reconstructed")]])
  cmp$within_tolerance <- apply(cmp[paste0(num, "_abs_diff")], 1, function(x) all(x < 1e-8 | is.na(x)))
  atomic_write_csv(cmp, file.path(recon_dir, "primary_contrast_reproduction.csv"))
  if (!all(cmp$within_tolerance)) {
    stop("Diagnostic reconstruction did not reproduce reported primary contrasts within tolerance")
  }

  for (p in file.path(out_dir, c("M2_SCORE_DECOMPOSITION.csv", "M2_WEIGHTING_DECOMPOSITION.csv",
                                "TEMPLATE_OUTLIER_AUDIT.csv"))) {
    if (file.exists(p)) unlink(p)
  }
  atomic_write_csv(score_decomp, file.path(out_dir, "M2_SCORE_DECOMPOSITION.csv"))
  atomic_write_csv(weighting_decomp, file.path(out_dir, "M2_WEIGHTING_DECOMPOSITION.csv"))
  atomic_write_csv(module_diag, file.path(recon_dir, "mechanism_module_summary.csv"))
  atomic_write_csv(cor_diag, file.path(recon_dir, "mechanism_correlations.csv"))
  atomic_write_csv(rank_diag, file.path(recon_dir, "mechanism_a2_rank_changes.csv"))
  atomic_write_csv(topo_diag, file.path(recon_dir, "mechanism_topology.csv"))

  wide <- reshape(weighting_decomp[, c("replicate", "template_key", "cell_type", "diagnostic_method", "auprc")],
                  idvar = c("replicate", "template_key", "cell_type"),
                  timevar = "diagnostic_method", direction = "wide")
  d0_col <- "auprc.D0_faithful_M2"
  d1_col <- "auprc.D1_M1"
  wide$m2_minus_m1 <- wide[[d0_col]] - wide[[d1_col]]
  outlier_templates <- c("Myeloid cells::grey60", "SMCs or Pericytes::salmon",
                         "SMCs or Pericytes::magenta", "Epithelial cells::black",
                         "Fibroblasts::green")
  loo <- lapply(unique(wide$template_key), function(tk) {
    keep <- wide$template_key != tk
    data.frame(
      template_key = tk,
      previously_unfavorable_template = tk %in% outlier_templates,
      n_replicates_removed = sum(!keep),
      mean_m2_minus_m1_without_template = mean(wide$m2_minus_m1[keep], na.rm = TRUE),
      median_m2_minus_m1_without_template = stats::median(wide$m2_minus_m1[keep], na.rm = TRUE),
      mean_m2_minus_m1_in_template = mean(wide$m2_minus_m1[!keep], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  outlier <- do.call(rbind, loo)
  atomic_write_csv(outlier, file.path(out_dir, "TEMPLATE_OUTLIER_AUDIT.csv"))

  mean_d0_d1 <- mean(wide$m2_minus_m1, na.rm = TRUE)
  d0 <- weighting_decomp[weighting_decomp$diagnostic_method == "D0_faithful_M2",]
  d1 <- weighting_decomp[weighting_decomp$diagnostic_method == "D1_M1",]
  score_wide <- reshape(score_decomp[, c("replicate", "score_definition", "auprc")],
                        idvar = "replicate", timevar = "score_definition", direction = "wide")
  a2_expr <- module_diag[module_diag$group == "A2",]
  conclusion <- "no single mechanism identified"
  if (is.finite(mean_d0_d1) && mean_d0_d1 < 0 &&
      mean(a2_expr$expression_percentile_median < 0.5, na.rm = TRUE) > 0.5) {
    conclusion <- "expression weighting suppresses A2"
  } else if ("auprc.final_heat_abs" %in% names(score_wide) &&
             "auprc.delta_twas_abs" %in% names(score_wide) &&
             mean(score_wide$auprc.final_heat_abs - score_wide$auprc.delta_twas_abs, na.rm = TRUE) > 0.02) {
    conclusion <- "delta_twas score is misaligned with module recovery"
  } else if (stats::sd(outlier$mean_m2_minus_m1_without_template, na.rm = TRUE) > 0.03) {
    conclusion <- "template-specific topology drives heterogeneity"
  }

  report <- c(
    "# M2 Diagnostic Report", "",
    "Status: mandatory diagnostic reconstruction completed before the new pilot.", "",
    paste0("- Previous frozen replicate design source: `", reps_path, "`."),
    paste0("- Completed unfavorable pilot report source: `", completed_pilot_dir(), "`."),
    "- Reconstructed primary contrasts reproduced the reported pilot contrasts within floating-point tolerance.",
    "- D0 is faithful submitted M2. D2-D5 are diagnostic proposed variants and are not submitted NESTA.",
    paste0("- Mean D0 faithful M2 minus D1 M1 AUPRC: ", signif(mean_d0_d1, 4), "."),
    paste0("- Median A2 expression percentile across previous replicates: ",
           signif(stats::median(a2_expr$expression_percentile_median, na.rm = TRUE), 4), "."),
    paste0("- Diagnostic conclusion: **", conclusion, "**."), "",
    "Required source tables:",
    "- `M2_SCORE_DECOMPOSITION.csv`",
    "- `M2_WEIGHTING_DECOMPOSITION.csv`",
    "- `TEMPLATE_OUTLIER_AUDIT.csv`",
    "- `results/diagnostics/reconstruction/mechanism_*.csv`"
  )
  report_path <- file.path(out_dir, "M2_DIAGNOSTIC_REPORT.md")
  if (file.exists(report_path)) unlink(report_path)
  atomic_write_lines(report, report_path)
  copy_required_diag_outputs()
  list(pass = TRUE, conclusion = conclusion, contrasts = contrasts,
       score_decomposition = score_decomp, weighting_decomposition = weighting_decomp,
       template_outlier_audit = outlier)
}
