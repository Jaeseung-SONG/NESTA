source(file.path("/home/js/NESTA/simulation", "R", "utils.R"))
source(project_file("R/fidelity.R"))
source(project_file("R/tom_library.R"))
source(project_file("R/generator.R"))

difficulty_settings_0628 <- function() {
  base <- list(
    a1_n = 10, a2_n = 40, relay_n = 28, a2_mean = 0.15, a2_sd = 0.32, a2_abs_cap = 1.30,
    decoy_n = 60, decoy_mean = 1.15,
    direct_min = 0.15, direct_max = 0.35, two_hop_min = 0.45, two_hop_max = 0.80,
    three_plus_max = 0.20, same_component_min = 0.90,
    median_min = 1.7, median_max = 2.5,
    top100_fe_min = NA_real_, top100_recall_min = NA_real_, raw_auprc_min = NA_real_,
    raw_diff_min = NA_real_, nodiff_diff_min = NA_real_, weight_perm_diff_min = NA_real_,
    top150_recall_min = NA_real_, top200_recall_min = NA_real_, compression_min = 0.100,
    relay_top100_rate_max = 0.14, relay_top50_rate_max = 0.07,
    relay_to_a2_top100_ratio_max = 1.00, relay_to_a2_heat_ratio_warning = 0.95,
    retention_a_min = 0.10, retention_a2_min = 0.025, leakage_background_max = 0.80,
    leakage_high_degree_max = 0.08, leakage_opposite_max = 0.10,
    sign100_min = NA_real_, sign5pct_min = NA_real_, direction_accuracy_min = 0.75,
    decoy100_max = 0.15, decoy5pct_max = 0.15,
    degree_rho_max = 0.20, degree_rho_preferred = 0.15, c_decoy100_max = 0.12,
    c_decoy5pct_max = 0.12)
  mk <- function(name, n_genes, sparsity, primary = FALSE) {
    x <- base
    x$difficulty_setting <- name
    x$network_universe_size <- n_genes
    x$background_sparsity <- sparsity
    x$primary_manuscript_setting <- primary
    x$sign5pct_min <- if (name == "dense_condition") 0.10 else if (name == "basic_condition") 0.08 else 0.05
    x
  }
  list(
    dense_condition = mk("dense_condition", 1000, "dense_current_050726_2012_like", FALSE),
    basic_condition = mk("basic_condition", 3000, "current_like_scaled", TRUE),
    sparse_condition = mk("sparse_condition", 5000, "zero_proximal_background", FALSE)
  )
}

setting_aware_degree_qc <- function(degree_audit, setting) {
  n <- degree_audit$n_A2_primary
  n_ext <- degree_audit$n_A2_extreme_high_degree
  frac <- degree_audit$fraction_A2_extreme_high_degree
  degree_audit$target_count_pilot_eligible <- n >= 30 && n <= 50
  degree_audit$target_count_absolute_minimum <- n >= 30
  degree_audit$target_count_warning_low <- n < 40
  degree_audit$target_count_acceptable <- n >= 30 && n <= 50
  degree_audit$setting_aware_extreme_cap <- ceiling(0.25 * max(n, 1))
  degree_audit$degree_hard_gate_pass <- is.finite(frac) && frac <= 0.25 &&
    is.finite(n_ext) && n_ext <= degree_audit$setting_aware_extreme_cap &&
    isTRUE(degree_audit$target_count_absolute_minimum)
  degree_audit$median_A2_degree_percentile_warning_low <- is.finite(degree_audit$median_A2_degree_percentile) &&
    degree_audit$median_A2_degree_percentile < 20
  degree_audit$median_A2_degree_percentile_warning_high <- is.finite(degree_audit$median_A2_degree_percentile) &&
    degree_audit$median_A2_degree_percentile > 75
  degree_audit$preferred_extreme_fraction_pass <- is.finite(frac) && frac <= 0.15
  degree_audit$degree_distribution_qc_pass <- isTRUE(degree_audit$degree_hard_gate_pass)
  degree_audit
}

score_label_final_heat <- "NESTA_final_heat"
primary_target_0628 <- "A2_intermediate_degree_capped"

signed_auprc_0628 <- function(symbols, signed_score, positives, direction, exclude = character()) {
  score <- if (direction == "risk") signed_score else -signed_score
  auprc_from_score(symbols, score, positives, exclude)
}

directional_qc_row_0628 <- function(rep) {
  data.frame(
    replicate = rep$rep_id,
    template_key = rep$template_key,
    n_A1_risk = length(rep$modules$A1_risk),
    n_A1_protective = length(rep$modules$A1_protective),
    n_A2_risk = length(rep$modules$A2_risk),
    n_A2_protective = length(rep$modules$A2_protective),
    n_A2_intermediate_degree_capped_risk = length(rep$modules$A2_intermediate_degree_capped_risk),
    n_A2_intermediate_degree_capped_protective = length(rep$modules$A2_intermediate_degree_capped_protective),
    n_D_opposite_sign_decoy = length(rep$modules$D_opposite_sign_decoy),
    same_direction_connectivity_fraction = if (length(rep$modules$A2_intermediate_degree_capped)) {
      mean(rep$modules$A2_intermediate_degree_capped %in% names(rep$directions$A2))
    } else 0,
    directional_qc_pass = length(rep$modules$A1_risk) >= 3 &&
      length(rep$modules$A1_protective) >= 3 &&
      length(rep$modules$A2_intermediate_degree_capped_risk) > 0 &&
      length(rep$modules$A2_intermediate_degree_capped_protective) > 0 &&
      length(rep$modules$D_opposite_sign_decoy) >= 40,
    stringsAsFactors = FALSE)
}

control_path_stats_0628 <- function(adj, modules, threshold = 0.1) {
  a <- modules$A
  a1 <- modules$A1
  a2 <- modules$A2_intermediate_degree_capped
  sub <- adj[a, a, drop = FALSE]
  g <- igraph::graph_from_adjacency_matrix((sub > threshold) * 1, mode = "undirected", diag = FALSE)
  comps <- igraph::components(g)
  ai <- match(a1, a)
  ti <- match(a2, a)
  d <- suppressWarnings(igraph::distances(g, v = ai, to = ti))
  min_path <- apply(d, 2, min, na.rm = TRUE)
  same <- vapply(ti, function(j) any(comps$membership[j] == comps$membership[ai]), logical(1))
  tom_vals <- as.numeric(adj[a1, a2, drop = FALSE])
  c(median_A1_A2_path = if (length(min_path[is.finite(min_path)])) stats::median(min_path[is.finite(min_path)]) else Inf,
    A2_same_component_fraction = if (length(same)) mean(same) else NA_real_,
    A1_A2_TOM_mean = if (length(tom_vals)) mean(tom_vals, na.rm = TRUE) else NA_real_,
    within_A_TOM_mean = mean(adj[a, a][lower.tri(adj[a, a])], na.rm = TRUE))
}

signed_branch_preservation_0628 <- function(adj, modules, directions) {
  a2 <- modules$A2_intermediate_degree_capped
  if (!length(a2)) return(NA_real_)
  vals <- vapply(a2, function(g) {
    dir <- unname(directions$A2[g])
    same <- if (identical(dir, "risk")) modules$A1_risk else modules$A1_protective
    opp <- if (identical(dir, "risk")) modules$A1_protective else modules$A1_risk
    mean(adj[g, same], na.rm = TRUE) > mean(adj[g, opp], na.rm = TRUE)
  }, logical(1))
  mean(vals, na.rm = TRUE)
}

branch_specificity_audit_0628 <- function(rep) {
  adj <- rep$networks$relevant
  modules <- rep$modules
  dirs <- rep$directions$A2
  a2 <- modules$A2_intermediate_degree_capped
  a2_risk <- intersect(a2, names(dirs)[dirs == "risk"])
  a2_prot <- intersect(a2, names(dirs)[dirs == "protective"])
  a1_risk <- modules$A1_risk
  a1_prot <- modules$A1_protective
  bg <- setdiff(rownames(adj), c(modules$A, modules$D_opposite_sign_decoy, modules$C_high_degree_decoy))
  safe_mean <- function(x, y) {
    x <- intersect(x, rownames(adj)); y <- intersect(y, colnames(adj))
    if (!length(x) || !length(y)) return(NA_real_)
    mean(as.numeric(adj[x, y, drop = FALSE]), na.rm = TRUE)
  }
  max_to <- function(targets, seeds) {
    targets <- intersect(targets, rownames(adj)); seeds <- intersect(seeds, colnames(adj))
    if (!length(targets) || !length(seeds)) return(NA_real_)
    mean(apply(adj[targets, seeds, drop = FALSE], 1, max, na.rm = TRUE), na.rm = TRUE)
  }
  within_risk <- max_to(a2_risk, a1_risk)
  within_prot <- max_to(a2_prot, a1_prot)
  within <- mean(c(within_risk, within_prot), na.rm = TRUE)
  between <- mean(c(max_to(a2_prot, a1_risk), max_to(a2_risk, a1_prot),
                    safe_mean(a2_risk, a2_prot)), na.rm = TRUE)
  background <- safe_mean(a2, bg)
  high_decoy <- safe_mean(a2, modules$C_high_degree_decoy)
  single <- vapply(a2, function(g) {
    dir <- unname(dirs[g])
    same <- if (identical(dir, "risk")) a1_risk else a1_prot
    opp <- if (identical(dir, "risk")) a1_prot else a1_risk
    max_same <- max(adj[g, same], na.rm = TRUE)
    max_opp <- max(adj[g, opp], na.rm = TRUE)
    max_same > max_opp * 1.05
  }, logical(1))
  data.frame(
    within_risk_branch_TOM = within_risk,
    within_protective_branch_TOM = within_prot,
    between_risk_protective_TOM = between,
    A2_background_TOM = background,
    A2_high_degree_decoy_TOM = high_decoy,
    within_to_between_TOM_ratio = within / pmax(between, 1e-9),
    within_to_background_TOM_ratio = within / pmax(background, 1e-9),
    fraction_A2_with_single_dominant_branch = mean(single, na.rm = TRUE),
    branch_specificity_qc_pass = is.finite(within) && is.finite(between) && is.finite(background) &&
      within / pmax(between, 1e-9) >= 2.00 &&
      within / pmax(background, 1e-9) >= 10.00 &&
      mean(single, na.rm = TRUE) >= 0.80,
    stringsAsFactors = FALSE)
}

select_relay_genes_0630 <- function(adj, modules, a2_direction, setting, seed, threshold = 0.1) {
  set.seed(seed)
  cfg <- get_calibration_candidate_0705()
  a2 <- modules$A2_intermediate_degree_capped
  pool <- intersect(setdiff(modules$A, c(modules$A1, a2)), rownames(adj))
  if (!length(pool)) return(list(genes = character(), direction = character()))
  outside <- setdiff(rownames(adj), modules$A)
  relay_bg_penalty <- if (!is.null(cfg$relay_bg_penalty)) as.numeric(cfg$relay_bg_penalty) else 0.80
  relay_cross_penalty <- if (!is.null(cfg$relay_cross_penalty)) as.numeric(cfg$relay_cross_penalty) else 0.50
  score_pool <- function(dir) {
    a1 <- if (identical(dir, "risk")) modules$A1_risk else modules$A1_protective
    targets <- intersect(names(a2_direction)[a2_direction == dir], colnames(adj))
    if (!length(a1) || !length(targets)) return(data.frame(SYMBOL = character(), score = numeric()))
    same_seed <- apply(adj[pool, a1, drop = FALSE], 1, max, na.rm = TRUE)
    target_link <- apply(adj[pool, targets, drop = FALSE], 1, max, na.rm = TRUE)
    same_seed_mean <- rowMeans(adj[pool, a1, drop = FALSE], na.rm = TRUE)
    target_link_mean <- rowMeans(adj[pool, targets, drop = FALSE], na.rm = TRUE)
    bg_leak <- if (length(outside)) rowMeans(adj[pool, outside, drop = FALSE], na.rm = TRUE) else 0
    cross_seed <- if (identical(dir, "risk")) modules$A1_protective else modules$A1_risk
    cross <- if (length(cross_seed)) apply(adj[pool, cross_seed, drop = FALSE], 1, max, na.rm = TRUE) else 0
    data.frame(SYMBOL = pool,
               score = pmin(same_seed, target_link) + 0.25 * same_seed + 0.25 * target_link +
                 0.60 * same_seed_mean + 0.45 * target_link_mean -
                 relay_bg_penalty * bg_leak - relay_cross_penalty * cross,
               stringsAsFactors = FALSE)
  }
  relay_target <- setting$relay_n
  if (identical(setting$difficulty_setting, "rare_target_detection") && !is.null(cfg$relay_n_target)) {
    relay_target <- as.integer(cfg$relay_n_target)
  }
  n_total <- min(relay_target, length(pool))
  n_risk <- ceiling(n_total / 2)
  risk <- score_pool("risk")
  prot <- score_pool("protective")
  risk_genes <- head(risk$SYMBOL[order(risk$score, decreasing = TRUE)], n_risk)
  prot_pool <- prot[!(prot$SYMBOL %in% risk_genes), , drop = FALSE]
  prot_genes <- head(prot_pool$SYMBOL[order(prot_pool$score, decreasing = TRUE)], n_total - length(risk_genes))
  chosen <- unique(c(risk_genes, prot_genes))
  if (length(chosen) < n_total) {
    combined <- rbind(risk, prot)
    combined <- combined[combined$SYMBOL %in% setdiff(pool, chosen), , drop = FALSE]
    chosen <- unique(c(chosen, head(combined$SYMBOL[order(combined$score, decreasing = TRUE)], n_total - length(chosen))))
  }
  relay_dir <- setNames(rep("risk", length(chosen)), chosen)
  relay_dir[intersect(chosen, prot_genes)] <- "protective"
  list(genes = chosen, direction = relay_dir[chosen])
}

diffusion_retention_audit_0628 <- function(rep, rel_scores) {
  dat <- rel_scores[!(rel_scores$SYMBOL %in% rep$modules$A1), , drop = FALSE]
  score <- setNames(abs(dat$final_NESTA_heat), dat$SYMBOL)
  total <- sum(score, na.rm = TRUE)
  frac <- function(genes) {
    genes <- intersect(genes, names(score))
    if (!length(genes) || !is.finite(total) || total <= 0) return(NA_real_)
    sum(score[genes], na.rm = TRUE) / total
  }
  ranked <- names(sort(score, decreasing = TRUE))
  a2 <- rep$modules$A2_intermediate_degree_capped
  top100 <- ranked[seq_len(min(100, length(ranked)))]
  top50 <- ranked[seq_len(min(50, length(ranked)))]
  top200 <- ranked[seq_len(min(200, length(ranked)))]
  top5pct <- ranked[seq_len(min(max(1, ceiling(0.05 * length(ranked))), length(ranked)))]
  top100_recall <- if (length(a2)) mean(a2 %in% top100) else NA_real_
  top200_recall <- if (length(a2)) mean(a2 %in% top200) else NA_real_
  relay <- if (!is.null(rep$modules$A_relay)) rep$modules$A_relay else character()
  relay_top100_n <- length(intersect(relay, top100))
  a2_top100_n <- length(intersect(a2, top100))
  relay_path_strength <- if (length(a2) && !is.null(rep$networks$relevant)) {
    adj <- rep$networks$relevant
    vapply(a2, function(g) {
      dir <- unname(rep$directions$A2[g])
      seeds <- if (identical(dir, "risk")) rep$modules$A1_risk else rep$modules$A1_protective
      r <- if (identical(dir, "risk")) rep$modules$A_relay_risk else rep$modules$A_relay_protective
      if (!length(seeds) || !length(r)) return(NA_real_)
      max(vapply(intersect(r, rownames(adj)), function(rr) pmin(max(adj[rr, seeds], na.rm = TRUE), adj[g, rr]), numeric(1)), na.rm = TRUE)
    }, numeric(1))
  } else numeric()
  data.frame(
    fraction_seed_heat_retained_in_A_branch = frac(rep$modules$A),
    fraction_seed_heat_reaching_A2 = frac(a2),
    fraction_seed_heat_leaking_to_background = frac(rep$modules$background),
    fraction_seed_heat_leaking_to_high_degree_decoys = frac(rep$modules$C_high_degree_decoy),
    fraction_seed_heat_leaking_to_opposite_sign_decoys = frac(rep$modules$D_opposite_sign_decoy),
    A2_final_heat_rank_compression_from_top200_to_top100 = top100_recall / pmax(top200_recall, 1e-9),
    relay_gene_top100_rate = if (length(relay)) mean(relay %in% top100) else NA_real_,
    relay_gene_top5pct_rate = if (length(relay)) mean(relay %in% top5pct) else NA_real_,
    relay_gene_top50_rate = if (length(relay)) mean(relay %in% top50) else NA_real_,
    relay_to_A2_top100_ratio = relay_top100_n / pmax(a2_top100_n, 1),
    relay_heat_mass = frac(relay),
    A2_heat_mass = frac(a2),
    relay_to_A2_heat_ratio = frac(relay) / pmax(frac(a2), 1e-9),
    fraction_A2_reached_via_relay = if (length(relay_path_strength)) mean(relay_path_strength > 0.1, na.rm = TRUE) else NA_real_,
    top100_recall = top100_recall,
    top200_recall = top200_recall,
    stringsAsFactors = FALSE)
}

control_disruption_audit_0628 <- function(rep) {
  rel <- control_path_stats_0628(rep$networks$relevant, rep$modules)
  i2 <- control_path_stats_0628(rep$networks$I2, rep$modules)
  i3 <- control_path_stats_0628(rep$networks$I3, rep$modules)
  tom_reduction <- function(ctrl) if (is.finite(rel["within_A_TOM_mean"]) && rel["within_A_TOM_mean"] > 0) {
    1 - ctrl["within_A_TOM_mean"] / rel["within_A_TOM_mean"]
  } else NA_real_
  path_or_component <- function(ctrl) {
    is.finite(ctrl["median_A1_A2_path"]) && is.finite(rel["median_A1_A2_path"]) &&
      (ctrl["median_A1_A2_path"] - rel["median_A1_A2_path"]) >= 1.0
  }
  i2_tom_red <- tom_reduction(i2)
  i3_tom_red <- tom_reduction(i3)
  out <- data.frame(
    relevant_median_A1_A2_path = rel["median_A1_A2_path"],
    I2_median_A1_A2_path = i2["median_A1_A2_path"],
    I3_median_A1_A2_path = i3["median_A1_A2_path"],
    relevant_A2_same_component_fraction = rel["A2_same_component_fraction"],
    I2_A2_same_component_fraction = i2["A2_same_component_fraction"],
    I3_A2_same_component_fraction = i3["A2_same_component_fraction"],
    relevant_A1_A2_TOM_mean = rel["A1_A2_TOM_mean"],
    I2_A1_A2_TOM_mean = i2["A1_A2_TOM_mean"],
    I3_A1_A2_TOM_mean = i3["A1_A2_TOM_mean"],
    within_A_TOM_reduction_fraction_I2 = i2_tom_red,
    within_A_TOM_reduction_fraction_I3 = i3_tom_red,
    signed_branch_preservation_relevant = signed_branch_preservation_0628(rep$networks$relevant, rep$modules, rep$directions),
    signed_branch_preservation_I2 = signed_branch_preservation_0628(rep$networks$I2, rep$modules, rep$directions),
    signed_branch_preservation_I3 = signed_branch_preservation_0628(rep$networks$I3, rep$modules, rep$directions),
    control_disruption_qc_pass = ((is.finite(i2_tom_red) && i2_tom_red >= 0.80) || path_or_component(i2)) &&
      ((is.finite(i3_tom_red) && i3_tom_red >= 0.80) || path_or_component(i3)),
    stringsAsFactors = FALSE)
  names(out) <- sub("^[^.]+\\.", "", names(out))
  out
}

assign_a1_directions_0628 <- function(a1, rep_id) {
  n <- length(a1)
  n_risk <- if (rep_id %% 2 == 1) ceiling(n / 2) else floor(n / 2)
  risk <- a1[seq_len(n_risk)]
  protective <- setdiff(a1, risk)
  list(risk = risk, protective = protective,
       majority = if (length(risk) >= length(protective)) "risk" else "protective")
}

select_path_targets_0628 <- function(adj, module_genes, a1, setting, seed, threshold = 0.1) {
  strat <- select_path_stratified_targets(adj, module_genes, a1, seed = seed, threshold = threshold)
  paths <- strat$paths
  pool <- names(paths)[is.finite(paths)]
  n0 <- min(max(setting$a2_n * 4, setting$a2_n), length(pool))
  set.seed(seed + 101)
  n1_target <- if (identical(setting$difficulty_setting, "rare_target_detection")) 0.28 else 0.24
  n2_target <- if (identical(setting$difficulty_setting, "rare_target_detection")) 0.62 else 0.58
  n2 <- min(sum(paths == 2), ceiling(n2_target * n0))
  n1 <- min(sum(paths == 1), ceiling(n1_target * n0))
  n3 <- min(sum(paths >= 3), n0 - n1 - n2)
  chosen <- c(sample_or_all(names(paths)[paths == 2], n2),
              sample_or_all(names(paths)[paths == 1], n1),
              sample_or_all(names(paths)[paths >= 3], n3))
  fill <- setdiff(pool, chosen)
  if (length(chosen) < n0 && length(fill)) {
    chosen <- c(chosen, sample_or_all(fill, n0 - length(chosen)))
  }
  chosen <- unique(chosen)
  p2 <- paths[chosen]
  metrics <- path_fraction_table(p2)
  metrics$path_stratification_pass <- with(metrics,
    is.finite(median_A1_A2_path) &
      median_A1_A2_path >= setting$median_min &
      median_A1_A2_path <= setting$median_max &
      A2_same_component_fraction >= setting$same_component_min &
      direct_1hop_fraction >= setting$direct_min &
      direct_1hop_fraction <= setting$direct_max &
      two_hop_fraction >= setting$two_hop_min &
      two_hop_fraction <= setting$two_hop_max &
      three_plus_hop_fraction <= setting$three_plus_max)
  metrics$path_fallback_used <- FALSE
  list(adj = strat$adj,
       A2_all = chosen,
       A2_proximal = chosen[p2 == 1],
       A2_intermediate = chosen[p2 == 2],
       A2_distal = chosen[p2 >= 3],
       paths = p2,
       metrics = metrics,
       thinned_edges = strat$thinned_edges)
}

select_degree_capped_targets_0628 <- function(candidates, degree_tab, seed, target_n) {
  set.seed(seed)
  dt <- degree_tab[match(candidates, degree_tab$SYMBOL), , drop = FALSE]
  dt <- dt[!is.na(dt$SYMBOL), , drop = FALSE]
  n0 <- min(target_n, nrow(dt))
  if (!n0) return(character())
  cap_extreme <- min(7, ceiling(0.25 * n0))
  dt$degree_score <- -abs(dt$degree_percentile - 55)
  dt <- dt[order(dt$degree_bin == "extreme_high_degree", -dt$degree_score), , drop = FALSE]
  non_extreme <- dt$SYMBOL[dt$degree_bin != "extreme_high_degree"]
  extreme <- dt$SYMBOL[dt$degree_bin == "extreme_high_degree"]
  chosen <- head(non_extreme, n0)
  if (length(chosen) < n0 && length(extreme)) {
    chosen <- c(chosen, head(extreme, min(cap_extreme, n0 - length(chosen))))
  }
  unique(chosen)
}

initial_rank_table_0628 <- function(genes, mean_expression, twas, a1 = character()) {
  twas_z <- strict_twas_vector(genes, twas, cutoff = 1)
  init <- (mean_expression / safe_sd(mean_expression)) * (twas_z / safe_sd(twas_z))
  ranked <- setdiff(genes, a1)
  raw_rank <- setNames(rank(-abs(twas_z[ranked]), ties.method = "average"), ranked)
  init_rank <- setNames(rank(-abs(init[ranked]), ties.method = "average"), ranked)
  data.frame(SYMBOL = ranked,
             raw_TWAS_rank = as.numeric(raw_rank[ranked]),
             M2_initial_rank = as.numeric(init_rank[ranked]),
             abs_TWAS_Z = abs(twas_z[ranked]),
             abs_M2_initial_weight = abs(init[ranked]),
             stringsAsFactors = FALSE)
}

target_initial_signal_audit_0628 <- function(rep_or_genes, mean_expression = NULL, twas = NULL,
                                             a2 = NULL, a1 = character(), setting = NULL,
                                             fallback_used = FALSE,
                                             fallback_fraction = if (isTRUE(fallback_used)) 1 else 0) {
  if (is.list(rep_or_genes) && !is.null(rep_or_genes$universe)) {
    rep <- rep_or_genes
    genes <- rep$universe$genes
    mean_expression <- rep$universe$mean_expression
    twas <- rep$twas
    a2 <- rep$modules$A2_intermediate_degree_capped
    a1 <- rep$modules$A1
    setting <- difficulty_settings_0628()[[rep$difficulty_setting]]
    fallback_used <- isTRUE(rep$target_initial_signal_fallback_used)
    fallback_fraction <- if (!is.null(rep$target_initial_signal_fallback_fraction)) {
      rep$target_initial_signal_fallback_fraction
    } else if (fallback_used) {
      1
    } else {
      0
    }
  } else {
    genes <- rep_or_genes
  }
  ranks <- initial_rank_table_0628(genes, mean_expression, twas, a1)
  sub <- ranks[ranks$SYMBOL %in% a2, , drop = FALSE]
  m2_top100_max <- 0.10
  raw_top100_max <- 0.10
  data.frame(
    n_A2_primary = length(a2),
    median_abs_TWAS_Z_A2 = stats::median(sub$abs_TWAS_Z, na.rm = TRUE),
    median_abs_M2_initial_weight_A2 = stats::median(sub$abs_M2_initial_weight, na.rm = TRUE),
    fraction_A2_in_raw_TWAS_top100 = if (nrow(sub)) mean(sub$raw_TWAS_rank <= 100, na.rm = TRUE) else NA_real_,
    fraction_A2_in_raw_TWAS_top200 = if (nrow(sub)) mean(sub$raw_TWAS_rank <= 200, na.rm = TRUE) else NA_real_,
    fraction_A2_in_raw_TWAS_top5pct = if (nrow(sub)) mean(sub$raw_TWAS_rank <= ceiling(0.05 * nrow(ranks)), na.rm = TRUE) else NA_real_,
    fraction_A2_in_M2_initial_top100 = if (nrow(sub)) mean(sub$M2_initial_rank <= 100, na.rm = TRUE) else NA_real_,
    fraction_A2_in_M2_initial_top200 = if (nrow(sub)) mean(sub$M2_initial_rank <= 200, na.rm = TRUE) else NA_real_,
    fraction_A2_in_M2_initial_top5pct = if (nrow(sub)) mean(sub$M2_initial_rank <= ceiling(0.05 * nrow(ranks)), na.rm = TRUE) else NA_real_,
    median_raw_TWAS_rank_A2 = stats::median(sub$raw_TWAS_rank, na.rm = TRUE),
    median_M2_initial_rank_A2 = stats::median(sub$M2_initial_rank, na.rm = TRUE),
    raw_top100_max = raw_top100_max,
    M2_initial_top100_max = m2_top100_max,
    target_initial_signal_fallback_used = fallback_used,
    target_initial_signal_fallback_fraction = fallback_fraction,
    M2_initial_top200_warning_max = 0.35,
    median_M2_initial_rank_preferred_min = 400,
    target_initial_signal_qc_pass = nrow(sub) > 0 &&
      mean(sub$raw_TWAS_rank <= 100, na.rm = TRUE) <= raw_top100_max &&
      mean(sub$M2_initial_rank <= 100, na.rm = TRUE) <= m2_top100_max &&
      mean(sub$raw_TWAS_rank <= ceiling(0.05 * nrow(ranks)), na.rm = TRUE) <= raw_top100_max &&
      mean(sub$M2_initial_rank <= ceiling(0.05 * nrow(ranks)), na.rm = TRUE) <= m2_top100_max,
    stringsAsFactors = FALSE)
}

select_comparator_aware_targets_0628 <- function(candidates, degree_tab, genes, mean_expression,
                                                 twas, a1, setting, seed, target_n) {
  ranks <- initial_rank_table_0628(genes, mean_expression, twas, a1)
  cand <- ranks[ranks$SYMBOL %in% candidates, , drop = FALSE]
  strict <- cand$SYMBOL[cand$raw_TWAS_rank > 125 & cand$M2_initial_rank > 200]
  moderate <- setdiff(cand$SYMBOL[cand$raw_TWAS_rank > 100 & cand$M2_initial_rank > 150], strict)
  relaxed <- cand$SYMBOL[cand$raw_TWAS_rank > 100 & cand$M2_initial_rank > 100]
  planned_moderate_n <- min(length(moderate),
                            floor(if (identical(setting$difficulty_setting, "rare_target_detection")) 0.15 * target_n else 0.10 * target_n))
  strict_target <- max(0, target_n - planned_moderate_n)
  selected <- select_degree_capped_targets_0628(strict, degree_tab, seed = seed, target_n = strict_target)
  if (planned_moderate_n > 0) {
    mod_selected <- select_degree_capped_targets_0628(moderate, degree_tab, seed = seed + 1,
                                                      target_n = planned_moderate_n)
    selected <- unique(c(selected, mod_selected))
  }
  pool <- unique(c(strict, moderate))
  fallback_fraction <- 0
  if (length(selected) < min(target_n, length(pool))) {
    fill <- setdiff(pool, selected)
    selected <- unique(c(selected, sample_or_all(fill, min(target_n - length(selected), length(fill)))))
  }
  if (length(selected) < min(target_n, 12)) {
    pool <- unique(c(selected, relaxed))
    fallback_fraction <- max(fallback_fraction, (target_n - length(selected)) / max(target_n, 1))
  }
  if (length(pool) < min(target_n, 12)) {
    ord <- cand$SYMBOL[order(pmin(cand$raw_TWAS_rank, cand$M2_initial_rank), decreasing = TRUE)]
    pool <- unique(c(pool, ord))
    fallback_fraction <- max(fallback_fraction, (target_n - length(selected)) / max(target_n, 1))
  }
  if (length(selected) < min(target_n, length(pool))) {
    fill <- setdiff(pool, selected)
    selected <- unique(c(selected, sample_or_all(fill, min(target_n - length(selected), length(fill)))))
  }
  list(targets = selected,
       fallback_used = fallback_fraction > 0,
       fallback_fraction = fallback_fraction)
}

rebalance_path_mix_0628 <- function(selected, candidates, paths, genes, mean_expression, twas,
                                    a1, degree_tab, setting, seed) {
  set.seed(seed)
  selected <- unique(selected)
  ranks <- initial_rank_table_0628(genes, mean_expression, twas, a1)
  eligible <- ranks$SYMBOL[ranks$SYMBOL %in% candidates &
                             ranks$raw_TWAS_rank > 100 &
                             ranks$M2_initial_rank > 150]
  eligible <- setdiff(eligible, selected)
  score_candidate <- function(g) {
    d <- degree_tab$degree_percentile[match(g, degree_tab$SYMBOL)]
    ifelse(is.finite(d), -abs(d - 55), -Inf)
  }
  choose_add <- function(pool) {
    pool <- intersect(pool, eligible)
    if (!length(pool)) return(NA_character_)
    pool[order(vapply(pool, score_candidate, numeric(1)), decreasing = TRUE)][1]
  }
  choose_remove <- function(pool) {
    pool <- intersect(pool, selected)
    if (!length(pool)) return(NA_character_)
    pool[order(vapply(pool, score_candidate, numeric(1)), decreasing = FALSE)][1]
  }
  for (iter in seq_len(100)) {
    met <- path_fraction_table(paths[selected])
    add <- remove <- NA_character_
    if (is.finite(met$two_hop_fraction) && met$two_hop_fraction > setting$two_hop_max) {
      add <- choose_add(names(paths)[paths == 1])
      remove <- choose_remove(selected[paths[selected] == 2])
    } else if (is.finite(met$direct_1hop_fraction) && met$direct_1hop_fraction < setting$direct_min) {
      add <- choose_add(names(paths)[paths == 1])
      remove <- choose_remove(selected[paths[selected] != 1])
    } else if (is.finite(met$direct_1hop_fraction) && met$direct_1hop_fraction > setting$direct_max) {
      add <- choose_add(names(paths)[paths == 2])
      remove <- choose_remove(selected[paths[selected] == 1])
    } else if (is.finite(met$three_plus_hop_fraction) && met$three_plus_hop_fraction > setting$three_plus_max) {
      add <- choose_add(names(paths)[paths %in% c(1, 2)])
      remove <- choose_remove(selected[paths[selected] >= 3])
    } else {
      break
    }
    if (is.na(add) || is.na(remove)) break
    selected <- c(setdiff(selected, remove), add)
    eligible <- c(setdiff(eligible, add), remove)
  }
  unique(selected)
}

calibration_default_0705 <- function() {
  list(candidate_id = "baseline_observed_adaptive_0705", stage = "default", parent_id = "",
       a1_relay_center = 0.315, relay_a2_center = 0.210,
       a2_local_clustering_cap = 0.93, opposite_suppression = "strong",
       relay_n_target = 30, a1_relay_contact_n = 9, relay_a2_contact_n = 16,
       a1_relay_factor = 1.36, a1_relay_floor = 2.70,
       relay_a2_factor = 1.82, relay_a2_floor = 2.18,
       bridge_seed_factor = 1.58, bridge_seed_floor = 1.95,
       bridge_a2_factor = 1.75, bridge_a2_floor = 2.02,
       a2_pair_multiplier = 2.55, a2_pair_factor = 1.26,
       a2_pair_floor = 1.06, a2_pair_ceiling = 1.38,
       relay_relay_factor = 0.48, relay_relay_ceiling = 0.68,
       opposite_decoy_factor = 0.12, opposite_decoy_ceiling = 0.36,
       high_decoy_factor = 0.12, high_decoy_ceiling = 0.34,
       relay_bg_penalty = 0.65, relay_cross_penalty = 0.35)
}

calibration_candidates_0705 <- function() {
  # Stage 1 uses observed-metric adaptive sparse-strong contacts: enough
  # selected branch paths to recover 1607-like path strength without inflating
  # the all-pair relay-A2 TOM mean.
  data.frame(
    candidate_id = sprintf("s1_c%02d", 1:12),
    stage = rep("stage1", 12),
    parent_id = rep("", 12),
    a1_relay_center = c(0.310, 0.315, 0.320, 0.325, 0.310, 0.315, 0.320, 0.325, 0.330, 0.335, 0.325, 0.330),
    relay_a2_center = c(0.205, 0.205, 0.210, 0.210, 0.215, 0.215, 0.220, 0.220, 0.215, 0.220, 0.225, 0.225),
    a2_local_clustering_cap = c(0.90, 0.92, 0.93, 0.95, 0.90, 0.92, 0.93, 0.95, 0.92, 0.93, 0.95, 0.93),
    opposite_suppression = c("strong", "strong", "moderate", "strong", "strong", "moderate", "strong", "strong", "moderate", "strong", "strong", "moderate"),
    relay_n_target = c(30, 30, 30, 30, 30, 30, 30, 30, 31, 31, 31, 31),
    a1_relay_contact_n = c(8, 9, 9, 10, 8, 9, 9, 10, 9, 10, 9, 10),
    relay_a2_contact_n = c(6, 8, 8, 10, 8, 10, 10, 12, 8, 10, 10, 12),
    a1_relay_factor = c(1.34, 1.38, 1.42, 1.46, 1.36, 1.40, 1.44, 1.48, 1.44, 1.48, 1.46, 1.50),
    a1_relay_floor = c(2.70, 2.76, 2.82, 2.88, 2.76, 2.82, 2.88, 2.94, 2.90, 2.98, 2.94, 3.02),
    relay_a2_factor = c(2.02, 2.02, 2.08, 2.08, 2.12, 2.12, 2.18, 2.18, 2.15, 2.20, 2.22, 2.25),
    relay_a2_floor = c(2.60, 2.48, 2.56, 2.44, 2.62, 2.50, 2.58, 2.48, 2.66, 2.54, 2.62, 2.52),
    bridge_seed_factor = c(1.58, 1.60, 1.62, 1.64, 1.60, 1.62, 1.64, 1.66, 1.64, 1.66, 1.66, 1.68),
    bridge_seed_floor = c(2.05, 2.08, 2.10, 2.12, 2.08, 2.10, 2.12, 2.14, 2.12, 2.14, 2.14, 2.16),
    bridge_a2_factor = c(1.92, 1.92, 1.98, 1.98, 2.02, 2.02, 2.08, 2.08, 2.05, 2.10, 2.12, 2.15),
    bridge_a2_floor = c(2.25, 2.18, 2.24, 2.16, 2.30, 2.20, 2.28, 2.18, 2.32, 2.22, 2.30, 2.20),
    a2_pair_multiplier = c(2.15, 2.25, 2.30, 2.40, 2.20, 2.30, 2.40, 2.50, 2.35, 2.45, 2.50, 2.55),
    a2_pair_factor = c(1.18, 1.20, 1.22, 1.24, 1.20, 1.22, 1.24, 1.26, 1.22, 1.24, 1.26, 1.28),
    a2_pair_floor = c(1.02, 1.04, 1.04, 1.06, 1.04, 1.06, 1.06, 1.08, 1.06, 1.08, 1.08, 1.10),
    a2_pair_ceiling = c(1.30, 1.32, 1.34, 1.36, 1.32, 1.34, 1.36, 1.38, 1.34, 1.36, 1.38, 1.40),
    relay_relay_factor = c(0.42, 0.44, 0.44, 0.46, 0.42, 0.44, 0.44, 0.46, 0.44, 0.46, 0.44, 0.46),
    relay_relay_ceiling = c(0.60, 0.62, 0.62, 0.64, 0.60, 0.62, 0.62, 0.64, 0.62, 0.64, 0.62, 0.64),
    opposite_decoy_factor = c(0.10, 0.10, 0.16, 0.10, 0.10, 0.16, 0.10, 0.10, 0.16, 0.10, 0.10, 0.16),
    opposite_decoy_ceiling = c(0.34, 0.34, 0.42, 0.34, 0.34, 0.42, 0.34, 0.34, 0.42, 0.34, 0.34, 0.42),
    high_decoy_factor = c(0.10, 0.10, 0.12, 0.10, 0.10, 0.12, 0.10, 0.10, 0.12, 0.10, 0.10, 0.12),
    high_decoy_ceiling = c(0.32, 0.32, 0.36, 0.32, 0.32, 0.36, 0.32, 0.32, 0.36, 0.32, 0.32, 0.36),
    relay_bg_penalty = c(0.65, 0.65, 0.68, 0.68, 0.65, 0.68, 0.68, 0.70, 0.68, 0.70, 0.70, 0.72),
    relay_cross_penalty = c(0.35, 0.35, 0.38, 0.38, 0.35, 0.38, 0.38, 0.40, 0.38, 0.40, 0.40, 0.42),
    stringsAsFactors = FALSE)
}

observed_window_pass_0705 <- function(row) {
  row <- row[1, , drop = FALSE]
  vals <- c(
    row$rare_A1_relay_mean_TOM >= 0.260 && row$rare_A1_relay_mean_TOM <= 0.300,
    row$rare_relay_A2_mean_TOM >= 0.170 && row$rare_relay_A2_mean_TOM <= 0.205,
    row$rare_path_strength >= 0.70 && row$rare_path_strength <= 0.85,
    row$rare_relay_count >= 24 && row$rare_relay_count <= 31,
    row$rare_A2_local_clustering >= 0.85 && row$rare_A2_local_clustering <= 0.97,
    row$rare_A_branch_background_cut_fraction <= 0.12,
    row$rare_seed_neighborhood_background_fraction <= 0.01,
    row$rare_high_degree_bridge_count <= 0.15,
    row$rare_opposite_sign_bridge_count <= 0.20,
    row$rare_branch_conductance_pass_fraction >= 0.80,
    row$rare_relay_structure_pass_fraction >= 0.95
  )
  all(is.finite(vals)) && all(vals)
}

adaptive_refinement_candidates_0705 <- function(audit, candidates) {
  if (!nrow(audit)) return(candidates[FALSE, , drop = FALSE])
  ord <- order(audit$closeness_to_1607_score, na.last = NA)
  base_ids <- unique(audit$candidate_id[ord])[seq_len(min(2, length(ord)))]
  rows <- list(); idx <- 1
  for (id in base_ids) {
    base <- candidates[candidates$candidate_id == id, , drop = FALSE]
    obs <- audit[audit$candidate_id == id, , drop = FALSE]
    if (!nrow(base) || !nrow(obs)) next
    for (variant in 1:2) {
      z <- base[1, , drop = FALSE]
      z$stage <- "stage2"
      z$parent_id <- id
      z$candidate_id <- paste0("s2_", id, "_r", variant)
      need_a1 <- if (is.finite(obs$rare_A1_relay_mean_TOM)) 0.2811 - obs$rare_A1_relay_mean_TOM else 0.04
      need_a2 <- if (is.finite(obs$rare_relay_A2_mean_TOM)) 0.1819 - obs$rare_relay_A2_mean_TOM else 0.03
      need_path <- if (is.finite(obs$rare_path_strength)) 0.8040 - obs$rare_path_strength else 0.10
      bump <- if (variant == 1) 1 else 2
      if (need_a1 > 0) {
        z$a1_relay_contact_n <- pmin(10, z$a1_relay_contact_n + bump)
        z$a1_relay_factor <- z$a1_relay_factor + c(0.10, 0.18)[bump]
        z$a1_relay_floor <- z$a1_relay_floor + c(0.18, 0.30)[bump]
      } else {
        z$a1_relay_contact_n <- pmax(6, z$a1_relay_contact_n - bump)
        z$a1_relay_factor <- z$a1_relay_factor - c(0.05, 0.09)[bump]
        z$a1_relay_floor <- z$a1_relay_floor - c(0.10, 0.18)[bump]
      }
      if (need_a2 > 0) {
        z$relay_a2_contact_n <- pmin(24, z$relay_a2_contact_n + 2 * bump)
        z$relay_a2_factor <- z$relay_a2_factor + c(0.10, 0.18)[bump]
        z$relay_a2_floor <- z$relay_a2_floor + c(0.16, 0.28)[bump]
        z$bridge_a2_factor <- z$bridge_a2_factor + c(0.06, 0.10)[bump]
        z$bridge_a2_floor <- z$bridge_a2_floor + c(0.08, 0.14)[bump]
      } else {
        z$relay_a2_contact_n <- pmax(5, z$relay_a2_contact_n - 2 * bump)
        z$relay_a2_factor <- z$relay_a2_factor - c(0.06, 0.10)[bump]
        z$relay_a2_floor <- z$relay_a2_floor - c(0.10, 0.18)[bump]
      }
      if (need_path > 0) {
        z$relay_n_target <- pmin(31, z$relay_n_target + bump)
        z$bridge_seed_factor <- z$bridge_seed_factor + c(0.04, 0.08)[bump]
        z$bridge_a2_factor <- z$bridge_a2_factor + c(0.04, 0.08)[bump]
      }
      if (is.finite(obs$rare_high_degree_bridge_count) && obs$rare_high_degree_bridge_count > 0.15) {
        z$a1_relay_factor <- z$a1_relay_factor - 0.08
        z$a2_pair_ceiling <- pmin(z$a2_pair_ceiling, 1.36)
        z$high_decoy_factor <- pmin(z$high_decoy_factor, 0.08)
        z$high_decoy_ceiling <- pmin(z$high_decoy_ceiling, 0.28)
      }
      if (is.finite(obs$rare_opposite_sign_bridge_count) && obs$rare_opposite_sign_bridge_count > 0.20) {
        z$opposite_decoy_factor <- pmin(z$opposite_decoy_factor, 0.08)
        z$opposite_decoy_ceiling <- pmin(z$opposite_decoy_ceiling, 0.30)
      }
      z$a1_relay_factor <- max(1.05, z$a1_relay_factor)
      z$relay_a2_factor <- max(1.25, z$relay_a2_factor)
      z$relay_n_target <- as.integer(max(24, min(31, z$relay_n_target)))
      rows[[idx]] <- z; idx <- idx + 1
    }
  }
  if (!length(rows)) return(candidates[FALSE, , drop = FALSE])
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

set_calibration_candidate_0705 <- function(candidate) {
  if (is.data.frame(candidate)) candidate <- as.list(candidate[1, , drop = FALSE])
  options(nesta.conductance.calibration.0705 = candidate)
  invisible(candidate)
}

get_calibration_candidate_0705 <- function() {
  x <- getOption("nesta.conductance.calibration.0705", NULL)
  if (is.null(x)) x <- calibration_default_0705()
  x
}

boost_recovery_branches_0628 <- function(adj, modules, directions, setting, paths, threshold = 0.1) {
  out <- adj
  a1_risk <- modules$A1_risk
  a1_protective <- modules$A1_protective
  a2 <- modules$A2_intermediate_degree_capped
  relay <- if (!is.null(modules$A_relay)) modules$A_relay else character()
  relay_dir <- if (!is.null(modules$A_relay_direction)) modules$A_relay_direction else setNames(character(), character())
  cfg <- get_calibration_candidate_0705()
  if (!length(a2)) return(out)
  all_genes <- rownames(out)
  background <- if (!is.null(modules$background)) modules$background else setdiff(all_genes, modules$A)
  outside_a <- setdiff(all_genes, modules$A)
  boost_edge <- function(x, y, factor = 1.18, floor = NULL) {
    if (!length(x) || !length(y)) return()
    for (gx in x) for (gy in y) {
      if (!gx %in% rownames(out) || !gy %in% colnames(out) || identical(gx, gy)) next
      val <- out[gx, gy]
      if (!is.null(floor)) val <- max(val, floor)
      val <- min(0.995, val * factor)
      out[gx, gy] <<- val
      out[gy, gx] <<- val
    }
  }
  damp_edge <- function(x, y, factor = 0.50, ceiling_value = NULL) {
    x <- intersect(x, rownames(out)); y <- intersect(y, colnames(out))
    if (!length(x) || !length(y)) return()
    for (gx in x) for (gy in y) {
      if (identical(gx, gy)) next
      val <- out[gx, gy] * factor
      if (!is.null(ceiling_value)) val <- min(val, ceiling_value)
      out[gx, gy] <<- val
      out[gy, gx] <<- val
    }
  }
  # Restore the 1607-like branch geometry: keep background leakage low without
  # overcorrecting relay conductance away from the seed side.
  damp_edge(c(modules$A1, relay, a2), background, factor = 0.10, ceiling_value = threshold * 0.18)
  damp_edge(c(relay, a2), outside_a, factor = 0.15, ceiling_value = threshold * 0.22)
  non_branch_a <- setdiff(modules$A, c(modules$A1, relay, a2))
  damp_edge(relay, non_branch_a, factor = 0.55, ceiling_value = threshold * 0.90)
  whole_degree <- rowSums(out > threshold) - 1
  high_bg <- intersect(names(sort(whole_degree, decreasing = TRUE))[seq_len(min(80, length(whole_degree)))], outside_a)
  damp_edge(c(modules$A1, relay, a2), high_bg, factor = 0.08, ceiling_value = threshold * 0.18)
  for (g in a2) {
    dir <- unname(directions[g])
    same <- if (identical(dir, "risk")) a1_risk else a1_protective
    opp <- if (identical(dir, "risk")) a1_protective else a1_risk
    same_relay <- intersect(names(relay_dir)[relay_dir == dir], relay)
    opp_relay <- intersect(names(relay_dir)[relay_dir != dir], relay)
    if (!length(same)) next
    same <- same[order(out[g, same], decreasing = TRUE)]
    if (!length(same)) next
    damp_edge(g, opp, factor = 0.35, ceiling_value = threshold * 0.60)
    damp_edge(g, opp_relay, factor = 0.35, ceiling_value = threshold * 0.55)
    if (is.finite(paths[g]) && paths[g] == 1) {
      boost_edge(g, same[1], factor = 1.35, floor = threshold * 1.08)
    } else {
      seed <- same[1]
      bridge_pool <- if (length(same_relay)) same_relay else setdiff(modules$A, c(modules$A1, a2))
      bridge_pool <- bridge_pool[out[seed, bridge_pool] > threshold | out[g, bridge_pool] > threshold]
      if (length(bridge_pool)) {
        bridge_score <- pmin(out[seed, bridge_pool], out[g, bridge_pool]) -
          0.45 * rowSums(out[bridge_pool, outside_a, drop = FALSE] > threshold)
        bridge <- bridge_pool[which.max(bridge_score)]
        boost_edge(g, bridge, factor = cfg$bridge_a2_factor, floor = threshold * cfg$bridge_a2_floor)
        boost_edge(seed, bridge, factor = cfg$bridge_seed_factor, floor = threshold * cfg$bridge_seed_floor)
        damp_edge(bridge, outside_a, factor = 0.12, ceiling_value = threshold * 0.20)
      } else {
        boost_edge(g, seed, factor = 1.08, floor = threshold * 0.92)
      }
      # Keep 2-hop targets topologically 2-hop while preserving a signed TOM hint.
      if (is.finite(paths[g]) && paths[g] >= 2 && out[g, seed] < threshold) {
        val <- min(threshold * 0.98, max(out[g, seed], threshold * 0.90) * 1.05)
        out[g, seed] <- val; out[seed, g] <- val
      }
    }
  }
  for (dir in c("risk", "protective")) {
    seeds <- if (identical(dir, "risk")) a1_risk else a1_protective
    r <- intersect(names(relay_dir)[relay_dir == dir], relay)
    targets <- intersect(names(directions)[directions == dir], a2)
    a1_contact_n <- if (!is.null(cfg$a1_relay_contact_n)) as.integer(cfg$a1_relay_contact_n) else 3L
    relay_a2_contact_n <- if (!is.null(cfg$relay_a2_contact_n)) as.integer(cfg$relay_a2_contact_n) else 5L
    if (length(seeds) && length(r)) {
      for (rr in r) {
        ss <- seeds[order(out[rr, seeds], decreasing = TRUE)]
        boost_edge(rr, ss[seq_len(min(a1_contact_n, length(ss)))], factor = cfg$a1_relay_factor, floor = threshold * cfg$a1_relay_floor)
      }
    }
    if (length(r) && length(targets)) {
      for (g in targets) {
        rr <- r[order(out[g, r], decreasing = TRUE)]
        boost_edge(g, rr[seq_len(min(relay_a2_contact_n, length(rr)))], factor = cfg$relay_a2_factor, floor = threshold * cfg$relay_a2_floor)
      }
    }
    if (length(r) >= 2) damp_edge(r, r, factor = cfg$relay_relay_factor, ceiling_value = threshold * cfg$relay_relay_ceiling)
  }
  for (dir in c("risk", "protective")) {
    branch <- names(directions)[directions == dir]
    branch <- intersect(branch, a2)
    if (length(branch) >= 2) {
      pairs <- utils::combn(branch, 2)
      pair_scores <- apply(pairs, 2, function(z) out[z[1], z[2]])
      keep <- order(pair_scores, decreasing = TRUE)[seq_len(min(length(pair_scores), ceiling(cfg$a2_pair_multiplier * length(branch))))]
      for (i in keep) {
        g1 <- pairs[1, i]; g2 <- pairs[2, i]
        if (out[g1, g2] > threshold * 0.35) {
          val <- min(threshold * cfg$a2_pair_ceiling, max(out[g1, g2], threshold * cfg$a2_pair_floor) * cfg$a2_pair_factor)
          out[g1, g2] <- val
          out[g2, g1] <- val
        }
      }
    }
  }
  damp_edge(intersect(names(directions)[directions == "risk"], a2),
            intersect(names(directions)[directions == "protective"], a2),
            factor = 0.45, ceiling_value = threshold * 0.60)
  damp_edge(intersect(names(relay_dir)[relay_dir == "risk"], relay),
            intersect(names(relay_dir)[relay_dir == "protective"], relay),
            factor = 0.50, ceiling_value = threshold * 0.70)
  risk_branch <- intersect(c(a1_risk, names(relay_dir)[relay_dir == "risk"]), rownames(out))
  protective_branch <- intersect(c(a1_protective, names(relay_dir)[relay_dir == "protective"]), rownames(out))
  # 0704 single-change correction: suppress opposite-sign decoy bridge edges
  # without changing A2 initial signal or the relay-A2 sweet-spot target window.
  opp_decoy <- intersect(if (!is.null(modules$D_opposite_sign_decoy)) modules$D_opposite_sign_decoy else character(), rownames(out))
  d_branch <- if (!is.null(modules$D_branch)) modules$D_branch else setNames(character(), character())
  if (length(opp_decoy) && length(risk_branch) && length(protective_branch)) {
    for (d in opp_decoy) {
      near <- unname(d_branch[d])
      if (identical(near, "near_risk")) {
        damp_edge(d, protective_branch, factor = 0.18, ceiling_value = threshold * 0.45)
      } else if (identical(near, "near_protective")) {
        damp_edge(d, risk_branch, factor = 0.18, ceiling_value = threshold * 0.45)
      } else {
        max_risk <- max(out[d, risk_branch], na.rm = TRUE)
        max_prot <- max(out[d, protective_branch], na.rm = TRUE)
        if (is.finite(max_risk) && is.finite(max_prot) && max_risk > threshold * 1.15 && max_prot > threshold * 1.15) {
          if (max_risk >= max_prot) damp_edge(d, protective_branch, factor = cfg$opposite_decoy_factor, ceiling_value = threshold * cfg$opposite_decoy_ceiling)
          else damp_edge(d, risk_branch, factor = cfg$opposite_decoy_factor, ceiling_value = threshold * cfg$opposite_decoy_ceiling)
        }
      }
    }
  }
  diag(out) <- 1
  out
}

suppress_external_bridge_decoys_0705 <- function(adj, modules, threshold = 0.1) {
  out <- adj
  cfg <- get_calibration_candidate_0705()
  damp_pair <- function(x, y, factor, ceiling_value) {
    x <- intersect(x, rownames(out)); y <- intersect(y, colnames(out))
    if (!length(x) || !length(y)) return()
    for (gx in x) for (gy in y) {
      if (identical(gx, gy)) next
      val <- min(out[gx, gy] * factor, ceiling_value)
      out[gx, gy] <<- val
      out[gy, gx] <<- val
    }
  }
  risk_branch <- intersect(c(modules$A1_risk, modules$A_relay_risk), rownames(out))
  protective_branch <- intersect(c(modules$A1_protective, modules$A_relay_protective), rownames(out))
  if (!length(risk_branch) || !length(protective_branch)) return(out)
  c_decoy <- intersect(if (!is.null(modules$C_high_degree_decoy)) modules$C_high_degree_decoy else character(), rownames(out))
  if (length(c_decoy)) {
    for (d in c_decoy) {
      max_risk <- max(out[d, risk_branch], na.rm = TRUE)
      max_prot <- max(out[d, protective_branch], na.rm = TRUE)
      if (is.finite(max_risk) && is.finite(max_prot) && max_risk > threshold * 1.15 && max_prot > threshold * 1.15) {
        if (max_risk >= max_prot) damp_pair(d, protective_branch, cfg$high_decoy_factor, threshold * cfg$high_decoy_ceiling)
        else damp_pair(d, risk_branch, cfg$high_decoy_factor, threshold * cfg$high_decoy_ceiling)
      }
    }
  }
  d_decoy <- intersect(if (!is.null(modules$D_opposite_sign_decoy)) modules$D_opposite_sign_decoy else character(), rownames(out))
  d_branch <- if (!is.null(modules$D_branch)) modules$D_branch else setNames(character(), character())
  if (length(d_decoy)) {
    for (d in d_decoy) {
      near <- unname(d_branch[d])
      if (identical(near, "near_risk")) damp_pair(d, protective_branch, cfg$opposite_decoy_factor, threshold * cfg$opposite_decoy_ceiling)
      else if (identical(near, "near_protective")) damp_pair(d, risk_branch, cfg$opposite_decoy_factor, threshold * cfg$opposite_decoy_ceiling)
    }
  }
  diag(out) <- 1
  out
}

branch_conductance_audit_0630 <- function(rep, threshold = 0.1) {
  adj <- rep$networks$relevant
  modules <- rep$modules
  branch <- intersect(c(modules$A1, modules$A_relay, modules$A2_intermediate_degree_capped), rownames(adj))
  bg <- intersect(modules$background, rownames(adj))
  a1 <- intersect(modules$A1, rownames(adj))
  relay <- intersect(modules$A_relay, rownames(adj))
  a2 <- intersect(modules$A2_intermediate_degree_capped, rownames(adj))
  high <- intersect(modules$C_high_degree_decoy, rownames(adj))
  opp <- intersect(modules$D_opposite_sign_decoy, rownames(adj))
  lower_sum <- function(x) if (length(x) >= 2) sum(adj[x, x][lower.tri(adj[x, x])], na.rm = TRUE) else 0
  cut_sum <- function(x, y) if (length(x) && length(y)) sum(adj[x, y, drop = FALSE], na.rm = TRUE) else 0
  edge_count <- function(x, y = NULL) {
    if (is.null(y)) {
      if (length(x) < 2) return(0)
      return(sum((adj[x, x, drop = FALSE] > threshold)[lower.tri(adj[x, x, drop = FALSE])], na.rm = TRUE))
    }
    if (length(x) && length(y)) sum(adj[x, y, drop = FALSE] > threshold, na.rm = TRUE) else 0
  }
  strong_bridge_count <- function(nodes) {
    nodes <- intersect(nodes, rownames(adj))
    risk <- intersect(c(modules$A1_risk, modules$A_relay_risk), colnames(adj))
    prot <- intersect(c(modules$A1_protective, modules$A_relay_protective), colnames(adj))
    if (!length(nodes) || !length(risk) || !length(prot)) return(0)
    max_risk <- apply(adj[nodes, risk, drop = FALSE], 1, max, na.rm = TRUE)
    max_prot <- apply(adj[nodes, prot, drop = FALSE], 1, max, na.rm = TRUE)
    sum(max_risk > threshold * 1.15 & max_prot > threshold * 1.15, na.rm = TRUE)
  }
  internal_w <- lower_sum(branch)
  cut_w <- cut_sum(branch, bg)
  internal_edges <- edge_count(branch)
  cut_edges <- edge_count(branch, bg)
  seed_bg_edges <- edge_count(a1, bg)
  seed_total_edges <- edge_count(a1, setdiff(rownames(adj), a1))
  effective_strength <- vapply(a2, function(g) {
    dir <- unname(rep$directions$A2[g])
    seeds <- if (identical(dir, "risk")) modules$A1_risk else modules$A1_protective
    r <- if (identical(dir, "risk")) modules$A_relay_risk else modules$A_relay_protective
    if (!length(seeds) || !length(r)) return(NA_real_)
    max(vapply(intersect(r, rownames(adj)), function(rr) pmin(max(adj[rr, seeds], na.rm = TRUE), adj[g, rr]), numeric(1)), na.rm = TRUE)
  }, numeric(1))
  a2_volume_branch <- cut_sum(a2, branch)
  a2_volume_total <- cut_sum(a2, setdiff(rownames(adj), a2))
  data.frame(
    A_branch_internal_edge_weight = internal_w,
    A_branch_background_cut_weight = cut_w,
    A_branch_background_cut_fraction = cut_w / pmax(internal_w + cut_w, 1e-9),
    A_branch_internal_edge_fraction = internal_edges / pmax(internal_edges + cut_edges, 1),
    A1_relay_edge_weight = if (length(a1) && length(relay)) mean(adj[a1, relay, drop = FALSE], na.rm = TRUE) else NA_real_,
    relay_A2_edge_weight = if (length(relay) && length(a2)) mean(adj[relay, a2, drop = FALSE], na.rm = TRUE) else NA_real_,
    A1_A2_effective_path_strength = mean(effective_strength, na.rm = TRUE),
    A2_local_clustering = if (length(a2) >= 3) {
      g <- igraph::graph_from_adjacency_matrix((adj[a2, a2, drop = FALSE] > threshold) * 1, mode = "undirected", diag = FALSE)
      igraph::transitivity(g, type = "average", isolates = "zero")
    } else NA_real_,
    A2_branch_volume_fraction = a2_volume_branch / pmax(a2_volume_total, 1e-9),
    seed_neighborhood_background_fraction = seed_bg_edges / pmax(seed_total_edges, 1),
    high_degree_bridge_count = strong_bridge_count(high),
    opposite_sign_bridge_count = strong_bridge_count(opp),
    branch_conductance_qc_pass = cut_w / pmax(internal_w + cut_w, 1e-9) <= 0.35 &&
      internal_edges / pmax(internal_edges + cut_edges, 1) >= 0.50 &&
      seed_bg_edges / pmax(seed_total_edges, 1) <= 0.45 &&
      cut_w / pmax(internal_w + cut_w, 1e-9) <= 0.12 &&
      seed_bg_edges / pmax(seed_total_edges, 1) <= 0.01 &&
      strong_bridge_count(high) <= 0.15 &&
      strong_bridge_count(opp) <= 0.20,
    stringsAsFactors = FALSE)
}

relay_structure_audit_0630 <- function(rep, threshold = 0.1) {
  adj <- rep$networks$relevant
  modules <- rep$modules
  relay <- intersect(modules$A_relay, rownames(adj))
  a2 <- intersect(modules$A2_intermediate_degree_capped, rownames(adj))
  ranks <- initial_rank_table_0628(rep$universe$genes, rep$universe$mean_expression, rep$twas, rep$modules$A1)
  rsub <- ranks[ranks$SYMBOL %in% relay, , drop = FALSE]
  branch_path <- vapply(a2, function(g) {
    dir <- unname(rep$directions$A2[g])
    seeds <- if (identical(dir, "risk")) modules$A1_risk else modules$A1_protective
    r <- if (identical(dir, "risk")) modules$A_relay_risk else modules$A_relay_protective
    if (!length(seeds) || !length(r)) return(NA_real_)
    max(vapply(intersect(r, rownames(adj)), function(rr) pmin(max(adj[rr, seeds], na.rm = TRUE), adj[g, rr]), numeric(1)), na.rm = TRUE)
  }, numeric(1))
  data.frame(
    n_relay_genes = length(relay),
    relay_raw_TWAS_top100_fraction = if (nrow(rsub)) mean(rsub$raw_TWAS_rank <= 100, na.rm = TRUE) else NA_real_,
    relay_M2_initial_top100_fraction = if (nrow(rsub)) mean(rsub$M2_initial_rank <= 100, na.rm = TRUE) else NA_real_,
    A1_relay_mean_TOM = if (length(relay)) mean(adj[modules$A1, relay, drop = FALSE], na.rm = TRUE) else NA_real_,
    relay_A2_mean_TOM = if (length(relay) && length(a2)) mean(adj[relay, a2, drop = FALSE], na.rm = TRUE) else NA_real_,
    relay_background_mean_TOM = if (length(relay)) mean(adj[relay, modules$background, drop = FALSE], na.rm = TRUE) else NA_real_,
    relay_high_degree_decoy_TOM = if (length(relay) && length(modules$C_high_degree_decoy)) mean(adj[relay, modules$C_high_degree_decoy, drop = FALSE], na.rm = TRUE) else NA_real_,
    fraction_A2_with_relay_path = if (length(branch_path)) mean(branch_path > threshold, na.rm = TRUE) else NA_real_,
    median_A1_relay_A2_path_strength = stats::median(branch_path, na.rm = TRUE),
    relay_M2_initial_top100_preferred_pass = if (nrow(rsub)) mean(rsub$M2_initial_rank <= 100, na.rm = TRUE) <= 0.03 else TRUE,
    relay_M2_initial_top100_warning = if (nrow(rsub)) mean(rsub$M2_initial_rank <= 100, na.rm = TRUE) > 0.03 else FALSE,
    relay_count_warning_low = length(relay) >= 15 && length(relay) < 20,
    relay_structure_qc_pass = length(relay) >= 15 &&
      (if (nrow(rsub)) mean(rsub$raw_TWAS_rank <= 100, na.rm = TRUE) else 0) <= 0.10 &&
      (if (nrow(rsub)) mean(rsub$M2_initial_rank <= 100, na.rm = TRUE) else 0) <= 0.05 &&
      (if (is.finite(mean(adj[relay, modules$background, drop = FALSE], na.rm = TRUE))) mean(adj[relay, modules$background, drop = FALSE], na.rm = TRUE) else 0) <= 0.01 &&
      (if (length(modules$C_high_degree_decoy) && is.finite(mean(adj[relay, modules$C_high_degree_decoy, drop = FALSE], na.rm = TRUE))) mean(adj[relay, modules$C_high_degree_decoy, drop = FALSE], na.rm = TRUE) else 0) <= 0.01 &&
      mean(branch_path > threshold, na.rm = TRUE) >= 0.85,
    stringsAsFactors = FALSE)
}

generate_twas_0628 <- function(genes, a1_risk, a1_protective, a2_risk, a2_protective,
                               decoy, decoy_direction, setting, seed, null = FALSE) {
  set.seed(seed)
  z <- stats::rnorm(length(genes), 0, 1)
  names(z) <- genes
  if (!null) {
    z[a1_risk] <- abs(stats::rnorm(length(a1_risk), mean = 3.5, sd = 0.75))
    z[a1_protective] <- -abs(stats::rnorm(length(a1_protective), mean = 3.5, sd = 0.75))
    weak_draw <- function(n, mean_sign) {
      vals <- stats::rnorm(n, mean = mean_sign * setting$a2_mean, sd = setting$a2_sd)
      vals <- pmax(pmin(vals, setting$a2_abs_cap), -setting$a2_abs_cap)
      vals
    }
    if (length(a2_risk)) z[a2_risk] <- weak_draw(length(a2_risk), 1)
    if (length(a2_protective)) z[a2_protective] <- weak_draw(length(a2_protective), -1)
    if (length(decoy)) {
      dd <- decoy_direction[decoy]
      z[decoy[dd == "risk"]] <- abs(stats::rnorm(sum(dd == "risk"), mean = setting$decoy_mean, sd = 0.55))
      z[decoy[dd == "protective"]] <- -abs(stats::rnorm(sum(dd == "protective"), mean = setting$decoy_mean, sd = 0.55))
    }
  }
  data.frame(SYMBOL = genes, TWAS.Z = as.numeric(z), TWAS.P = twas_p_from_z(z),
             stringsAsFactors = FALSE)
}

make_replicate_0628 <- function(lib, template_key, rep_id, seed, setting, null = FALSE) {
  universe <- make_universe(lib, template_key, seed = seed,
                            n_genes = if (!is.null(setting$network_universe_size)) setting$network_universe_size else 1000)
  universe$network_universe_size <- length(universe$genes)
  universe$network_universe_condition <- setting$difficulty_setting
  universe$background_sparsity <- if (!is.null(setting$background_sparsity)) setting$background_sparsity else "current_like"
  rel <- make_relevant_network(lib, universe)
  set.seed(seed + 1)
  a1 <- non_extreme_seed_genes(universe$module_a_tom, universe$module_a, setting$a1_n)
  a1_dir <- assign_a1_directions_0628(a1, rep_id)
  strat <- select_path_targets_0628(rel, universe$module_a, a1, setting, seed = seed + 11)
  rel <- strat$adj
  degree_tab <- module_degree_table(rel, universe$module_a)
  candidate_dir <- assign_a2_direction(rel, universe$module_a, strat$A2_all,
                                       a1_dir$risk, a1_dir$protective, seed = seed + 12)
  candidate_risk <- names(candidate_dir)[candidate_dir == "risk"]
  candidate_protective <- names(candidate_dir)[candidate_dir == "protective"]
  decoy <- construct_decoy_module(rel, universe, a1_dir$risk, a1_dir$protective,
                                  seed = seed + 13, n = setting$decoy_n)
  twas <- generate_twas_0628(universe$genes, a1_dir$risk, a1_dir$protective,
                             candidate_risk, candidate_protective, decoy$genes, decoy$direction,
                             setting, seed = seed + 2, null = null)
  selected <- select_comparator_aware_targets_0628(strat$A2_all, degree_tab, universe$genes,
                                                   universe$mean_expression, twas, a1,
                                                   setting, seed = seed + 14,
                                                   target_n = setting$a2_n)
  primary_a2 <- rebalance_path_mix_0628(selected$targets, strat$A2_all, strat$paths,
                                        universe$genes, universe$mean_expression, twas,
                                        a1, degree_tab, setting, seed = seed + 16)
  paths <- strat$paths[primary_a2]
  selected_modules <- list(A = universe$module_a, A1 = a1,
                           A1_risk = a1_dir$risk, A1_protective = a1_dir$protective,
                           A2_intermediate_degree_capped = primary_a2,
                           background = universe$background)
  relay_design <- select_relay_genes_0630(rel, selected_modules, candidate_dir[primary_a2],
                                          setting, seed = seed + 17)
  relay_dir <- relay_design$direction
  if (length(relay_design$genes)) {
    set.seed(seed + 18)
    relay_sign <- ifelse(relay_dir[relay_design$genes] == "risk", 1, -1)
    ridx <- match(relay_design$genes, twas$SYMBOL)
    expr_cap <- as.numeric(stats::quantile(universe$mean_expression, 0.55, na.rm = TRUE))
    universe$mean_expression[relay_design$genes] <- pmin(universe$mean_expression[relay_design$genes], expr_cap)
    twas$TWAS.Z[ridx] <- pmax(pmin(stats::rnorm(length(ridx), mean = relay_sign * 0.03, sd = 0.10), 0.45), -0.45)
    twas$TWAS.P <- twas_p_from_z(setNames(twas$TWAS.Z, twas$SYMBOL))
  }
  selected_modules$A_relay <- relay_design$genes
  selected_modules$A_relay_direction <- relay_dir
  selected_modules$A_relay_risk <- names(relay_dir)[relay_dir == "risk"]
  selected_modules$A_relay_protective <- names(relay_dir)[relay_dir == "protective"]
  selected_modules$D_opposite_sign_decoy <- decoy$genes
  selected_modules$D_opposite_sign_decoy_direction <- decoy$direction
  selected_modules$D_branch <- decoy$branch
  rel <- boost_recovery_branches_0628(rel, selected_modules, candidate_dir[primary_a2],
                                      setting, paths)
  degree_tab <- module_degree_table(rel, universe$module_a)
  c_decoy <- construct_high_degree_decoy(rel, universe, degree_tab, seed = seed + 15, n = 50)
  selected_modules$C_high_degree_decoy <- c_decoy
  rel <- suppress_external_bridge_decoys_0705(rel, selected_modules)
  degree_tab <- module_degree_table(rel, universe$module_a)
  degree_audit <- degree_audit_metrics(primary_a2, degree_tab,
                                       setdiff(universe$module_a, c(a1, primary_a2)),
                                       c_decoy)
  degree_audit <- setting_aware_degree_qc(degree_audit, setting)
  a2_dir <- candidate_dir[primary_a2]
  a2_risk <- names(a2_dir)[a2_dir == "risk"]
  a2_protective <- names(a2_dir)[a2_dir == "protective"]
  target_initial_signal_qc <- target_initial_signal_audit_0628(universe$genes,
                                                               universe$mean_expression,
                                                               twas, primary_a2, a1,
                                                               setting,
                                                               selected$fallback_used,
                                                               selected$fallback_fraction)
  blank <- primary_a2[abs(twas$TWAS.Z[match(primary_a2, twas$SYMBOL)]) < 1.0]
  weak <- primary_a2[abs(twas$TWAS.Z[match(primary_a2, twas$SYMBOL)]) >= 1.0 &
                       twas$TWAS.P[match(primary_a2, twas$SYMBOL)] > 0.10]
  if (!length(weak)) weak <- setdiff(primary_a2, blank)
  i1 <- make_i1_network(lib, universe, seed + 3)
  i2 <- make_i2_network(lib, universe, rel, seed + 4)
  i3 <- make_i3_network(lib, universe, rel, seed + 5)
  path_metrics <- path_fraction_table(paths)
  path_metrics$path_stratification_pass <- with(path_metrics,
    is.finite(median_A1_A2_path) &
      median_A1_A2_path >= setting$median_min &
      median_A1_A2_path <= setting$median_max &
      A2_same_component_fraction >= setting$same_component_min &
      direct_1hop_fraction >= setting$direct_min &
      direct_1hop_fraction <= setting$direct_max &
      two_hop_fraction >= setting$two_hop_min &
      two_hop_fraction <= setting$two_hop_max &
      three_plus_hop_fraction <= setting$three_plus_max)
  path_metrics$path_fallback_used <- FALSE
  list(rep_id = rep_id, difficulty_setting = setting$difficulty_setting,
       template_key = template_key, universe = universe,
       modules = list(A = universe$module_a, B = universe$module_b, C = universe$module_c,
                      background = universe$background,
                      A1 = a1, A1_risk = a1_dir$risk, A1_protective = a1_dir$protective,
                      A2 = primary_a2, A2_all = primary_a2,
                      A2_proximal = primary_a2[paths == 1],
                      A2_intermediate = primary_a2[paths == 2],
                      A2_distal = primary_a2[paths >= 3],
                      A2_intermediate_distal = primary_a2[paths >= 2],
                      A2_intermediate_degree_capped = primary_a2,
                      A2_risk = a2_risk, A2_protective = a2_protective,
                      A2_TWAS_blank = blank, A2_TWAS_weak = weak,
                      A2_promoted_from_raw_TWAS = character(),
                      A2_low_degree = intersect(primary_a2, degree_tab$SYMBOL[degree_tab$degree_bin == "low_degree"]),
                      A2_mid_degree = intersect(primary_a2, degree_tab$SYMBOL[degree_tab$degree_bin == "mid_degree"]),
                      A2_high_degree = intersect(primary_a2, degree_tab$SYMBOL[degree_tab$degree_bin == "high_degree"]),
                      A2_extreme_high_degree = intersect(primary_a2, degree_tab$SYMBOL[degree_tab$degree_bin == "extreme_high_degree"]),
                      A2_intermediate_degree_capped_risk = intersect(primary_a2, a2_risk),
                      A2_intermediate_degree_capped_protective = intersect(primary_a2, a2_protective),
                      A_relay = relay_design$genes,
                      A_relay_risk = names(relay_dir)[relay_dir == "risk"],
                      A_relay_protective = names(relay_dir)[relay_dir == "protective"],
                      A_relay_direction = relay_dir,
                      A_bridge = setdiff(universe$module_a, c(a1, primary_a2, relay_design$genes)),
                      D_opposite_sign_decoy = decoy$genes,
                      C_high_degree_decoy = c_decoy),
       directions = list(A2 = a2_dir, D = decoy$direction, D_branch = decoy$branch,
                         relay = relay_dir,
                         A1_majority_direction = a1_dir$majority),
       path_qc = path_metrics, degree_qc = degree_audit, degree_table = degree_tab,
       target_initial_signal_qc = target_initial_signal_qc,
       target_initial_signal_fallback_used = selected$fallback_used,
       target_initial_signal_fallback_fraction = selected$fallback_fraction,
       path_lengths = paths, thinning_log = strat$thinned_edges, twas = twas,
       networks = list(relevant = rel, I1 = i1, I2 = i2, I3 = i3))
}

target_sets_0628 <- function(modules) {
  list(A2_all = modules$A2_all,
       A2_proximal = modules$A2_proximal,
       A2_intermediate = modules$A2_intermediate,
       A2_distal = modules$A2_distal,
       A2_intermediate_degree_capped = modules$A2_intermediate_degree_capped,
       A2_risk = modules$A2_risk,
       A2_protective = modules$A2_protective,
       A2_TWAS_blank = modules$A2_TWAS_blank,
       A2_TWAS_weak = modules$A2_TWAS_weak,
       A2_promoted_from_raw_TWAS = modules$A2_promoted_from_raw_TWAS,
       D_opposite_sign_decoy = modules$D_opposite_sign_decoy,
       C_high_degree_decoy = modules$C_high_degree_decoy)
}

ranked_genes_0628 <- function(scores, score_col, exclude = character()) {
  dat <- scores[!(scores$SYMBOL %in% exclude), , drop = FALSE]
  dat <- dat[is.finite(dat[[score_col]]), , drop = FALSE]
  dat$SYMBOL[order(dat[[score_col]], decreasing = TRUE)]
}

recovery_metrics_0628 <- function(scores, positives, exclude, score_col = "score") {
  ranked <- ranked_genes_0628(scores, score_col, exclude)
  positives <- positives[positives %in% ranked]
  n_ranked <- length(ranked)
  prev <- if (n_ranked) length(positives) / n_ranked else NA_real_
  au <- if (length(positives)) auprc_from_score(scores$SYMBOL, scores[[score_col]], positives, exclude) else NA_real_
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
    partial_AUPRC_top200 = if (length(positives)) {
      top200 <- ranked[seq_len(min(200, length(ranked)))]
      auprc_from_score(top200, rev(seq_len(length(top200))), positives)
    } else NA_real_,
    prevalence_normalized_AUPRC = if (is.finite(au) && is.finite(prev) && prev < 1) (au - prev) / (1 - prev) else NA_real_,
    first_target_rank = if (length(pos_ranks)) min(pos_ranks, na.rm = TRUE) else NA_real_,
    mean_reciprocal_rank = if (length(pos_ranks)) mean(1 / pos_ranks, na.rm = TRUE) else NA_real_,
    stringsAsFactors = FALSE)
  for (k in c(50, 100, 150, 200)) {
    top <- ranked[seq_len(min(k, length(ranked)))]
    recall <- if (length(positives)) mean(positives %in% top) else NA_real_
    random_expected <- if (n_ranked) min(k, n_ranked) / n_ranked else NA_real_
    out[[paste0("top", k, "_recall")]] <- recall
    out[[paste0("top", k, "_fold_enrichment_over_random")]] <- recall / random_expected
  }
  for (pct in c(0.01, 0.05, 0.10)) {
    k <- max(1, ceiling(n_ranked * pct))
    tag <- paste0("top", as.integer(pct * 100), "pct")
    top <- ranked[seq_len(min(k, length(ranked)))]
    recall <- if (length(positives)) mean(positives %in% top) else NA_real_
    random_expected <- if (n_ranked) min(k, n_ranked) / n_ranked else NA_real_
    out[[paste0(tag, "_recall")]] <- recall
    out[[paste0(tag, "_fold_enrichment_over_random")]] <- recall / random_expected
  }
  out
}

path_rows_0628 <- function(scores, modules, exclude, score_col, score_name, network_label,
                           method, replicate, template_key, difficulty_setting, null = FALSE) {
  rows <- lapply(names(target_sets_0628(modules)), function(target_name) {
    met <- recovery_metrics_0628(scores, target_sets_0628(modules)[[target_name]], exclude, score_col)
    met$target_set <- target_name
    met$score_name <- score_name
    met$network_label <- network_label
    met$method <- method
    met$replicate <- replicate
    met$template_key <- template_key
    met$cell_type <- sub("::.*", "", template_key)
    met$difficulty_setting <- difficulty_setting
    met$null <- null
    met
  })
  do.call(rbind, rows)
}

direction_row_0628 <- function(scores, modules, directions, exclude, rank_col, signed_col,
                               score_name, network_label, method, replicate, template_key,
                               difficulty_setting, null = FALSE) {
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
  top5pct <- ranked[seq_len(min(max(1, ceiling(0.05 * length(ranked))), length(ranked)))]
  top1pct <- ranked[seq_len(min(max(1, ceiling(0.01 * length(ranked))), length(ranked)))]
  top10pct <- ranked[seq_len(min(max(1, ceiling(0.10 * length(ranked))), length(ranked)))]
  d <- modules$D_opposite_sign_decoy
  dsel50 <- intersect(d, top50)
  dsel100 <- intersect(d, top100)
  dsel5pct <- intersect(d, top5pct)
  risk_au <- signed_auprc_0628(dat$SYMBOL, dat[[signed_col]], modules$A2_intermediate_degree_capped_risk, "risk", character())
  prot_au <- signed_auprc_0628(dat$SYMBOL, dat[[signed_col]], modules$A2_intermediate_degree_capped_protective, "protective", character())
  direction_aware_au <- stats::weighted.mean(c(risk_au, prot_au), c(length(modules$A2_intermediate_degree_capped_risk), length(modules$A2_intermediate_degree_capped_protective)), na.rm = TRUE)
  data.frame(
    score_name = score_name, network_label = network_label, method = method,
    replicate = replicate, template_key = template_key, cell_type = sub("::.*", "", template_key),
    difficulty_setting = difficulty_setting, null = null,
    sign_concordant_top50_recall = if (length(targets)) sum(concordant(top50), na.rm = TRUE) / length(targets) else NA_real_,
    sign_concordant_top100_recall = if (length(targets)) sum(concordant(top100), na.rm = TRUE) / length(targets) else NA_real_,
    sign_concordant_top1pct_recall = if (length(targets)) sum(concordant(top1pct), na.rm = TRUE) / length(targets) else NA_real_,
    sign_concordant_top5pct_recall = if (length(targets)) sum(concordant(top5pct), na.rm = TRUE) / length(targets) else NA_real_,
    sign_concordant_top10pct_recall = if (length(targets)) sum(concordant(top10pct), na.rm = TRUE) / length(targets) else NA_real_,
    sign_concordant_precision_top50 = if (length(intersect(top50, targets))) mean(concordant(top50), na.rm = TRUE) else NA_real_,
    sign_concordant_precision_top100 = if (length(intersect(top100, targets))) mean(concordant(top100), na.rm = TRUE) else NA_real_,
    direction_accuracy_among_top50_A2 = if (length(intersect(top50, targets))) mean(concordant(top50), na.rm = TRUE) else NA_real_,
    direction_accuracy_among_top100_A2 = if (length(intersect(top100, targets))) mean(concordant(top100), na.rm = TRUE) else NA_real_,
    direction_accuracy_among_top5pct_A2 = if (length(intersect(top5pct, targets))) mean(concordant(top5pct), na.rm = TRUE) else NA_real_,
    signed_AUPRC_risk = risk_au,
    signed_AUPRC_protective = prot_au,
    direction_aware_AUPRC = direction_aware_au,
    opposite_sign_decoy_top50_rate = if (length(d)) length(dsel50) / length(d) else NA_real_,
    opposite_sign_decoy_top100_rate = if (length(d)) length(dsel100) / length(d) else NA_real_,
    opposite_sign_decoy_top5pct_rate = if (length(d)) length(dsel5pct) / length(d) else NA_real_,
    D_fold_enrichment_over_random = if (length(d) && length(ranked)) {
      (length(dsel100) / length(d)) / (min(100, length(ranked)) / length(ranked))
    } else NA_real_,
    stringsAsFactors = FALSE)
}

degree_row_0628 <- function(scores, rep, score_col, score_name, network_label, method,
                            replicate, template_key, difficulty_setting, null = FALSE,
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
  top5pct <- names(sort(score, decreasing = TRUE))[seq_len(min(max(1, ceiling(0.05 * length(score))), length(score)))]
  cdecoy <- rep$modules$C_high_degree_decoy
  data.frame(
    score_name = score_name, network_label = network_label, method = method,
    replicate = replicate, template_key = template_key, cell_type = sub("::.*", "", template_key),
    difficulty_setting = difficulty_setting, null = null,
    n_A2_primary = length(rep$modules$A2_intermediate_degree_capped),
    n_A2_low_degree = length(rep$modules$A2_low_degree),
    n_A2_mid_degree = length(rep$modules$A2_mid_degree),
    n_A2_high_degree = length(rep$modules$A2_high_degree),
    n_A2_extreme_high_degree = length(rep$modules$A2_extreme_high_degree),
    fraction_A2_extreme_high_degree = if (length(rep$modules$A2_intermediate_degree_capped)) length(rep$modules$A2_extreme_high_degree) / length(rep$modules$A2_intermediate_degree_capped) else NA_real_,
    median_A2_degree_percentile = stats::median(rep$degree_table$degree_percentile[match(rep$modules$A2_intermediate_degree_capped, rep$degree_table$SYMBOL)], na.rm = TRUE),
    score_degree_spearman = suppressWarnings(stats::cor(score[names(score) %in% names(degree)], degree[names(score)[names(score) %in% names(degree)]], method = "spearman", use = "complete.obs")),
    score_strength_spearman = suppressWarnings(stats::cor(score[names(score) %in% names(strength)], strength[names(score)[names(score) %in% names(strength)]], method = "spearman", use = "complete.obs")),
    top50_degree_percentile_median = stats::median(deg_pct[top50], na.rm = TRUE),
    top100_degree_percentile_median = stats::median(deg_pct[top100], na.rm = TRUE),
    top5pct_degree_percentile_median = stats::median(deg_pct[top5pct], na.rm = TRUE),
    C_high_degree_decoy_top50_rate = if (length(cdecoy)) mean(cdecoy %in% top50) else NA_real_,
    C_high_degree_decoy_top100_rate = if (length(cdecoy)) mean(cdecoy %in% top100) else NA_real_,
    C_high_degree_decoy_top5pct_rate = if (length(cdecoy)) mean(cdecoy %in% top5pct) else NA_real_,
    degree_stratified_top100_recall_low = if (length(rep$modules$A2_low_degree)) mean(rep$modules$A2_low_degree %in% top100) else NA_real_,
    degree_stratified_top100_recall_mid = if (length(rep$modules$A2_mid_degree)) mean(rep$modules$A2_mid_degree %in% top100) else NA_real_,
    degree_stratified_top100_recall_high = if (length(rep$modules$A2_high_degree)) mean(rep$modules$A2_high_degree %in% top100) else NA_real_,
    degree_stratified_top100_recall_extreme = if (length(rep$modules$A2_extreme_high_degree)) mean(rep$modules$A2_extreme_high_degree %in% top100) else NA_real_,
    stringsAsFactors = FALSE)
}

score_payload_0628 <- function(rep, rep_id, null = FALSE, n_weight_perm = 10) {
  rows <- list(); dir_rows <- list(); degree_rows <- list()
  bench_rows <- list(); bench_dir_rows <- list(); bench_degree_rows <- list()
  add_primary <- function(scores, score_col, signed_col, score_name, network_label, method) {
    rows[[length(rows) + 1]] <<- path_rows_0628(scores, rep$modules, rep$modules$A1, score_col, score_name,
                                                network_label, method, rep_id, rep$template_key,
                                                rep$difficulty_setting, null)
    dir_rows[[length(dir_rows) + 1]] <<- direction_row_0628(scores, rep$modules, rep$directions,
                                                            rep$modules$A1, score_col, signed_col,
                                                            score_name, network_label, method, rep_id,
                                                            rep$template_key, rep$difficulty_setting, null)
    degree_rows[[length(degree_rows) + 1]] <<- degree_row_0628(scores, rep, score_col, score_name,
                                                               network_label, method, rep_id,
                                                               rep$template_key, rep$difficulty_setting, null)
  }
  rel <- faithful_m2_scores(rep$networks$relevant, rep$universe$mean_expression, rep$twas, method = "ber_p", n.perm = 300)
  rel$raw_abs_Z <- abs(rel$TWAS.Z)
  rel$raw_signed_TWAS <- rel$TWAS.Z
  add_primary(rel, "NESTA_M2_abs_final_heat", "final_NESTA_heat", score_label_final_heat, "relevant", "NESTA_M2_faithful")
  add_primary(rel, "NESTA_M2_abs_delta_NESTA", "delta_NESTA", "legacy_delta_only", "relevant", "NESTA_M2_faithful")
  add_primary(rel, "NESTA_M2_composite", "final_NESTA_heat", "old_heat_delta_composite", "relevant", "diagnostic_sensitivity")
  add_primary(rel, "raw_abs_Z", "raw_signed_TWAS", "raw_TWAS_abs", "relevant", "raw_TWAS")
  nd <- no_diffusion_m2_scores(rep$networks$relevant, rep$universe$mean_expression, rep$twas)
  add_primary(nd, "NESTA_M2_abs_final_heat", "final_NESTA_heat", "M2_no_diffusion", "relevant", "M2_no_diffusion")
  for (nm in c("I2", "I3")) {
    s <- faithful_m2_scores(rep$networks[[nm]], rep$universe$mean_expression, rep$twas, method = "ber_p", n.perm = 300)
    add_primary(s, "NESTA_M2_abs_final_heat", "final_NESTA_heat", score_label_final_heat, nm, "NESTA_M2_faithful")
  }
  if (!null && n_weight_perm > 0) {
    for (i in seq_len(n_weight_perm)) {
      adj <- permute_weights(rep$networks$relevant, seed = 804000 + rep_id * 100 + i)
      s <- faithful_m2_scores(adj, rep$universe$mean_expression, rep$twas, method = "ber_p", n.perm = 300)
      add_primary(s, "NESTA_M2_abs_final_heat", "final_NESTA_heat", score_label_final_heat,
                  paste0("weight_perm_", i), "NESTA_M2_faithful")
    }
  }
  if (!null) {
    for (net_nm in c("relevant", "I2", "I3")) {
      b <- benchmark_scores(rep$networks[[net_nm]], rep$twas, include_sensitivity = TRUE)
      for (score_nm in names(b)) {
        tab <- b[[score_nm]]
        bench_rows[[length(bench_rows) + 1]] <- path_rows_0628(tab, rep$modules, rep$modules$A1, "score",
                                                               score_nm, net_nm, score_nm, rep_id,
                                                               rep$template_key, rep$difficulty_setting, null)
        if ("score_signed" %in% names(tab) && !all(is.na(tab$score_signed))) {
          bench_dir_rows[[length(bench_dir_rows) + 1]] <- direction_row_0628(tab, rep$modules, rep$directions,
                                                                             rep$modules$A1, "score", "score_signed",
                                                                             score_nm, net_nm, score_nm, rep_id,
                                                                             rep$template_key, rep$difficulty_setting, null)
        }
        bench_degree_rows[[length(bench_degree_rows) + 1]] <- degree_row_0628(tab, rep, "score", score_nm,
                                                                              net_nm, score_nm, rep_id,
                                                                              rep$template_key, rep$difficulty_setting, null)
      }
    }
  }
  list(primary = do.call(rbind, rows),
       direction = do.call(rbind, dir_rows),
       degree = do.call(rbind, degree_rows),
       benchmark = if (length(bench_rows)) do.call(rbind, bench_rows) else data.frame(),
       benchmark_direction = if (length(bench_dir_rows)) do.call(rbind, bench_dir_rows) else data.frame(),
       benchmark_degree = if (length(bench_degree_rows)) do.call(rbind, bench_degree_rows) else data.frame(),
       rel_scores = rel)
}

aggregate_weight_perm_0628 <- function(metrics) {
  wp <- metrics[grepl("^weight_perm_", metrics$network_label), , drop = FALSE]
  if (!nrow(wp)) return(metrics)
  num <- names(wp)[vapply(wp, is.numeric, logical(1))]
  keys <- c("difficulty_setting", "replicate", "template_key", "cell_type", "target_set", "null")
  med <- aggregate(wp[, num, drop = FALSE], wp[, keys, drop = FALSE], stats::median, na.rm = TRUE)
  med$score_name <- "median_weight_permutation"
  med$network_label <- "median_weight_permutation"
  med$method <- "weight_permutation_control"
  for (nm in setdiff(names(metrics), names(med))) med[[nm]] <- NA
  rbind(metrics[!grepl("^weight_perm_", metrics$network_label), , drop = FALSE],
        med[, names(metrics), drop = FALSE])
}

wide_metric_0628 <- function(metrics, spec, metric) {
  x <- metrics
  for (nm in names(spec)) x <- x[x[[nm]] == spec[[nm]], , drop = FALSE]
  out <- x[, c("difficulty_setting", "replicate", metric), drop = FALSE]
  names(out)[3] <- "value"
  out
}

paired_contrast_0628 <- function(metrics, specs, metrics_to_compare,
                                 target_set = primary_target_0628) {
  rows <- list()
  m0 <- metrics[metrics$target_set == target_set, , drop = FALSE]
  for (setting in unique(m0$difficulty_setting)) {
    ms <- m0[m0$difficulty_setting == setting, , drop = FALSE]
    for (metric in metrics_to_compare) {
      base <- wide_metric_0628(ms, specs$base, metric)
      names(base)[3] <- "base"
      for (nm in setdiff(names(specs), "base")) {
        cmp <- wide_metric_0628(ms, specs[[nm]], metric)
        names(cmp)[3] <- "cmp"
        tab <- merge(base, cmp, by = c("difficulty_setting", "replicate"), all = FALSE)
        diff <- tab$base - tab$cmp
        ratio <- tab$base / pmax(tab$cmp, 1e-9)
        ci <- paired_bootstrap_ci(diff)
        rows[[length(rows) + 1]] <- data.frame(
          difficulty_setting = setting, target_set = target_set, contrast = nm, metric = metric,
          base_mean = mean(tab$base, na.rm = TRUE), comparator_mean = mean(tab$cmp, na.rm = TRUE),
          mean = ci["mean"], median = ci["median"], ci_low = ci["lo"], ci_high = ci["hi"],
          mean_fold_improvement = mean(ratio[is.finite(ratio)], na.rm = TRUE),
          improved_fraction = mean(diff > 0, na.rm = TRUE),
          stringsAsFactors = FALSE)
      }
    }
  }
  if (length(rows)) do.call(rbind, rows) else data.frame()
}

reprioritization_rows_0628 <- function(rep, rel_scores, rep_id) {
  dat <- rel_scores[!(rel_scores$SYMBOL %in% rep$modules$A1), , drop = FALSE]
  raw_rank <- setNames(rank(-abs(dat$TWAS.Z), ties.method = "average"), dat$SYMBOL)
  heat_rank <- setNames(rank(-abs(dat$final_NESTA_heat), ties.method = "average"), dat$SYMBOL)
  top100_heat <- names(heat_rank)[heat_rank <= 100]
  a2 <- rep$modules$A2_intermediate_degree_capped
  promoted <- a2[raw_rank[a2] > 200 & heat_rank[a2] <= 100]
  rep$modules$A2_promoted_from_raw_TWAS <- promoted
  data.frame(
    difficulty_setting = rep$difficulty_setting, replicate = rep_id,
    template_key = rep$template_key, cell_type = sub("::.*", "", rep$template_key),
    number_of_A2_genes_promoted_from_outside_raw_TWAS_top200_to_final_heat_top100 = length(promoted),
    fraction_of_final_heat_top100_A2_with_TWAS.P_gt_0.10 = {
      g <- intersect(top100_heat, a2)
      if (length(g)) mean(rel_scores$TWAS.P[match(g, rel_scores$SYMBOL)] > 0.10) else NA_real_
    },
    median_rank_improvement_A2_final_heat_vs_raw_TWAS = stats::median(raw_rank[a2] - heat_rank[a2], na.rm = TRUE),
    median_delta_NESTA_among_reprioritized_A2 = if (length(promoted)) stats::median(rel_scores$delta_NESTA[match(promoted, rel_scores$SYMBOL)], na.rm = TRUE) else NA_real_,
    A2_TWAS_blank_definition = "abs(TWAS.Z) < 1.0",
    stringsAsFactors = FALSE)
}

delta_rows_0628 <- function(rep, rel_scores, rep_id) {
  a2 <- rep$modules$A2_intermediate_degree_capped
  dat <- rel_scores[match(a2, rel_scores$SYMBOL), , drop = FALSE]
  data.frame(
    difficulty_setting = rep$difficulty_setting, replicate = rep_id,
    template_key = rep$template_key, cell_type = sub("::.*", "", rep$template_key),
    median_delta_NESTA_A2 = stats::median(dat$delta_NESTA, na.rm = TRUE),
    median_abs_delta_NESTA_A2 = stats::median(abs(dat$delta_NESTA), na.rm = TRUE),
    fraction_A2_delta_NESTA_direction_concordant = mean(ifelse(rep$directions$A2[a2] == "risk", dat$delta_NESTA > 0, dat$delta_NESTA < 0), na.rm = TRUE),
    stringsAsFactors = FALSE)
}

summarize_outputs_0628 <- function(rows, out_dir) {
  safe_dir_create(out_dir)
  primary <- aggregate_weight_perm_0628(do.call(rbind, lapply(rows, `[[`, "primary")))
  direction <- do.call(rbind, lapply(rows, `[[`, "direction"))
  degree <- do.call(rbind, lapply(rows, `[[`, "degree"))
  benchmark <- do.call(rbind, lapply(rows, `[[`, "benchmark"))
  benchmark_direction <- do.call(rbind, lapply(rows, `[[`, "benchmark_direction"))
  benchmark_degree <- do.call(rbind, lapply(rows, `[[`, "benchmark_degree"))
  reprior <- do.call(rbind, lapply(seq_along(rows), function(i) rows[[i]]$reprioritization))
  delta <- do.call(rbind, lapply(seq_along(rows), function(i) rows[[i]]$delta))
  diffusion_retention <- do.call(rbind, lapply(seq_along(rows), function(i) rows[[i]]$diffusion_retention))
  atomic_write_csv(primary, file.path(out_dir, "PRIMARY_FINAL_HEAT_METRICS.csv"))
  atomic_write_csv(direction, file.path(out_dir, "DIRECTION_AWARE_METRICS.csv"))
  atomic_write_csv(reprior, file.path(out_dir, "REPRIORITIZATION_METRICS.csv"))
  atomic_write_csv(delta, file.path(out_dir, "DELTA_INTERPRETATION_METRICS.csv"))
  atomic_write_csv(diffusion_retention, file.path(out_dir, "DIFFUSION_RETENTION_AUDIT.csv"))
  atomic_write_csv(degree, file.path(out_dir, "DEGREE_BIAS_METRICS.csv"))
  atomic_write_csv(benchmark, file.path(out_dir, "BENCHMARK_METRICS.csv"))
  atomic_write_csv(benchmark_direction, file.path(out_dir, "BENCHMARK_DIRECTION_AWARE_METRICS.csv"))
  atomic_write_csv(benchmark_degree, file.path(out_dir, "BENCHMARK_DEGREE_BIAS_METRICS.csv"))
  specs <- list(
    base = list(score_name = score_label_final_heat, network_label = "relevant", method = "NESTA_M2_faithful"),
    raw_TWAS_abs = list(score_name = "raw_TWAS_abs", network_label = "relevant", method = "raw_TWAS"),
    M2_no_diffusion = list(score_name = "M2_no_diffusion", network_label = "relevant", method = "M2_no_diffusion"),
    legacy_delta_only = list(score_name = "legacy_delta_only", network_label = "relevant", method = "NESTA_M2_faithful"),
    old_heat_delta_composite = list(score_name = "old_heat_delta_composite", network_label = "relevant", method = "diagnostic_sensitivity"),
    I2_module_disrupted = list(score_name = score_label_final_heat, network_label = "I2", method = "NESTA_M2_faithful"),
    I3_expression_matched_randomized = list(score_name = score_label_final_heat, network_label = "I3", method = "NESTA_M2_faithful"),
    median_weight_permutation = list(score_name = "median_weight_permutation", network_label = "median_weight_permutation", method = "weight_permutation_control"))
  primary_contrasts <- paired_contrast_0628(primary, specs,
    c("top50_recall", "top100_recall", "top150_recall", "top200_recall",
      "top50_fold_enrichment_over_random", "top100_fold_enrichment_over_random",
      "top150_fold_enrichment_over_random", "top200_fold_enrichment_over_random",
      "top1pct_recall", "top5pct_recall", "top10pct_recall",
      "top1pct_fold_enrichment_over_random", "top5pct_fold_enrichment_over_random",
      "top10pct_fold_enrichment_over_random",
      "raw_AUPRC", "partial_AUPRC_top100", "partial_AUPRC_top200",
      "prevalence_normalized_AUPRC"))
  atomic_write_csv(primary_contrasts, file.path(out_dir, "PRIMARY_FINAL_HEAT_CONTRASTS.csv"))
  direction_combo <- rbind(direction[, intersect(names(direction), names(benchmark_direction)), drop = FALSE],
                           benchmark_direction)
  direction_specs <- list(
    base = list(score_name = score_label_final_heat, network_label = "relevant"),
    raw_signed_TWAS = list(score_name = "raw_TWAS_abs", network_label = "relevant"),
    M2_no_diffusion = list(score_name = "M2_no_diffusion", network_label = "relevant"),
    RWR_signed_two_channel = list(score_name = "RWR_signed_two_channel", network_label = "relevant"),
    PPR_signed_two_channel = list(score_name = "PPR_signed_two_channel", network_label = "relevant"))
  dir_contrasts <- paired_contrast_0628(transform(direction_combo, target_set = primary_target_0628),
                                        direction_specs,
                                        c("direction_aware_AUPRC",
                                          "sign_concordant_top100_recall",
                                          "sign_concordant_top5pct_recall",
                                          "direction_accuracy_among_top100_A2",
                                          "direction_accuracy_among_top5pct_A2",
                                          "opposite_sign_decoy_top100_rate",
                                          "opposite_sign_decoy_top5pct_rate"))
  atomic_write_csv(dir_contrasts, file.path(out_dir, "DIRECTION_AWARE_CONTRASTS.csv"))
  degree_combo <- rbind(degree[, intersect(names(degree), names(benchmark_degree)), drop = FALSE],
                        benchmark_degree)
  deg_specs <- list(base = list(score_name = score_label_final_heat, network_label = "relevant"),
                    RWR_abs_prior = list(score_name = "RWR_abs_prior", network_label = "relevant"),
                    PPR_abs_prior = list(score_name = "PPR_abs_prior", network_label = "relevant"),
                    RWR_signed_two_channel = list(score_name = "RWR_signed_two_channel", network_label = "relevant"),
                    PPR_signed_two_channel = list(score_name = "PPR_signed_two_channel", network_label = "relevant"))
  deg_contrasts <- paired_contrast_0628(transform(degree_combo, target_set = primary_target_0628),
                                        deg_specs,
                                        c("C_high_degree_decoy_top100_rate", "score_degree_spearman",
                                          "score_strength_spearman", "top100_degree_percentile_median"))
  atomic_write_csv(deg_contrasts, file.path(out_dir, "DEGREE_BIAS_CONTRASTS.csv"))
  bench_combo <- rbind(primary[primary$score_name == score_label_final_heat & primary$network_label == "relevant",
                               names(benchmark), drop = FALSE],
                       benchmark)
  bench_specs <- list(base = list(score_name = score_label_final_heat, network_label = "relevant"),
                      RWR_abs_prior = list(score_name = "RWR_abs_prior", network_label = "relevant"),
                      PPR_abs_prior = list(score_name = "PPR_abs_prior", network_label = "relevant"),
                      RWR_signed_two_channel = list(score_name = "RWR_signed_two_channel", network_label = "relevant"),
                      PPR_signed_two_channel = list(score_name = "PPR_signed_two_channel", network_label = "relevant"))
  bench_contrasts <- paired_contrast_0628(bench_combo, bench_specs,
                                          c("top100_recall", "top5pct_recall",
                                            "top100_fold_enrichment_over_random",
                                            "top5pct_fold_enrichment_over_random",
                                            "raw_AUPRC", "prevalence_normalized_AUPRC"))
  atomic_write_csv(bench_contrasts, file.path(out_dir, "BENCHMARK_CONTRASTS.csv"))
  atomic_write_csv(primary, file.path(out_dir, "SIZE_STRATIFIED_PRIMARY_METRICS.csv"))
  atomic_write_csv(primary_contrasts, file.path(out_dir, "SIZE_STRATIFIED_PRIMARY_CONTRASTS.csv"))
  list(primary = primary, primary_contrasts = primary_contrasts,
       direction = direction, direction_contrasts = dir_contrasts,
       reprioritization = reprior, delta = delta,
       diffusion_retention = diffusion_retention,
       degree = degree, degree_contrasts = deg_contrasts,
       benchmark = benchmark, benchmark_direction = benchmark_direction,
       benchmark_degree = benchmark_degree, benchmark_contrasts = bench_contrasts)
}

run_topology_qc_0628 <- function(lib, n_reps = 40, out_dir = project_file("results/topology_qc"),
                                 seed_base = 2026280435, save_reps = TRUE) {
  safe_dir_create(out_dir)
  settings <- difficulty_settings_0628()
  templates <- paste(lib$templates$cell_type, lib$templates$module, sep = "::")
  if (length(templates) < 20 || length(unique(lib$templates$cell_type)) < 3) {
    atomic_write_csv(data.frame(pilot_go = FALSE, reason = "insufficient_empirical_templates"),
                     file.path(out_dir, "topology_qc_decision.csv"))
    return(list(pass = FALSE, reason = "insufficient_empirical_templates"))
  }
  reps <- list(); qc_rows <- list(); match_rows <- list(); path_rows <- list(); direction_rows <- list(); degree_rows <- list()
  target_signal_rows <- list(); control_disruption_rows <- list(); branch_rows <- list()
  conductance_rows <- list(); relay_rows <- list()
  outlier_rows <- list(); design_rows <- list()
  for (setting_name in names(settings)) {
    setting <- settings[[setting_name]]
    reps[[setting_name]] <- vector("list", n_reps)
    for (i in seq_len(n_reps)) {
      attempt <- 1
      repeat {
        key <- templates[((i + attempt - 2) %% length(templates)) + 1]
        rep_seed <- seed_base + match(setting_name, names(settings)) * 100000 + i * 100 + attempt
        rep <- make_replicate_0628(lib, key, i,
                                   seed = rep_seed,
                                   setting = setting)
        init_qc_tmp <- rep$target_initial_signal_qc
        relay_qc_tmp <- relay_structure_audit_0630(rep)
        sweet_spot_bad <- (
          (is.finite(relay_qc_tmp$relay_A2_mean_TOM) &&
             (relay_qc_tmp$relay_A2_mean_TOM < 0.175 || relay_qc_tmp$relay_A2_mean_TOM > 0.200)) ||
          (is.finite(relay_qc_tmp$median_A1_relay_A2_path_strength) &&
             (relay_qc_tmp$median_A1_relay_A2_path_strength < 0.72 || relay_qc_tmp$median_A1_relay_A2_path_strength > 0.85)) ||
          (is.finite(relay_qc_tmp$A1_relay_mean_TOM) && relay_qc_tmp$A1_relay_mean_TOM < 0.265) ||
          relay_qc_tmp$n_relay_genes < 26)
        severe_bad <- sweet_spot_bad ||
          (is.finite(init_qc_tmp$fraction_A2_in_M2_initial_top100) &&
                         init_qc_tmp$fraction_A2_in_M2_initial_top100 > 0.10) ||
          (is.finite(init_qc_tmp$fraction_A2_in_raw_TWAS_top100) &&
             init_qc_tmp$fraction_A2_in_raw_TWAS_top100 > 0.10) ||
          !isTRUE(rep$degree_qc$target_count_pilot_eligible) ||
          (is.finite(relay_qc_tmp$relay_M2_initial_top100_fraction) &&
             relay_qc_tmp$relay_M2_initial_top100_fraction > 0.05) ||
          (is.finite(relay_qc_tmp$relay_raw_TWAS_top100_fraction) &&
             relay_qc_tmp$relay_raw_TWAS_top100_fraction > 0.10) ||
          (is.finite(relay_qc_tmp$relay_background_mean_TOM) &&
             relay_qc_tmp$relay_background_mean_TOM > 0.01) ||
          (is.finite(relay_qc_tmp$relay_high_degree_decoy_TOM) &&
             relay_qc_tmp$relay_high_degree_decoy_TOM > 0.01) ||
          (is.finite(relay_qc_tmp$fraction_A2_with_relay_path) &&
             relay_qc_tmp$fraction_A2_with_relay_path < 0.85) ||
          relay_qc_tmp$n_relay_genes < 15 ||
          (is.finite(init_qc_tmp$target_initial_signal_fallback_fraction) &&
             init_qc_tmp$target_initial_signal_fallback_fraction > 0.30)
        if (!severe_bad || attempt >= length(templates) * 2) break
        outlier_rows[[length(outlier_rows) + 1]] <- data.frame(
          difficulty_setting = setting_name,
          requested_replicate = i,
          attempted_template_key = key,
          reason = if (isTRUE(sweet_spot_bad)) {
            "outside_0703_1607_structural_sweet_spot"
          } else if (is.finite(init_qc_tmp$target_initial_signal_fallback_fraction) &&
                       init_qc_tmp$target_initial_signal_fallback_fraction > 0.30) {
            "fallback_target_selection_required_for_more_than_30_percent_A2"
          } else if (is.finite(init_qc_tmp$fraction_A2_in_raw_TWAS_top100) &&
                     init_qc_tmp$fraction_A2_in_raw_TWAS_top100 > 0.10) {
            "raw_TWAS_top100_leakage_above_0.10"
          } else if (!isTRUE(rep$degree_qc$target_count_pilot_eligible)) {
            "A2_target_count_below_24_replacement_attempted"
          } else if (relay_qc_tmp$n_relay_genes < 15) {
            "relay_count_below_15_replacement_attempted"
          } else if (is.finite(relay_qc_tmp$relay_M2_initial_top100_fraction) &&
                     relay_qc_tmp$relay_M2_initial_top100_fraction > 0.05) {
            "relay_M2_initial_top100_fraction_above_0.05"
          } else if (is.finite(relay_qc_tmp$relay_background_mean_TOM) &&
                     relay_qc_tmp$relay_background_mean_TOM > 0.01) {
            "relay_background_TOM_above_0.01"
          } else if (is.finite(relay_qc_tmp$relay_high_degree_decoy_TOM) &&
                     relay_qc_tmp$relay_high_degree_decoy_TOM > 0.01) {
            "relay_high_degree_decoy_TOM_above_0.01"
          } else if (is.finite(relay_qc_tmp$fraction_A2_with_relay_path) &&
                     relay_qc_tmp$fraction_A2_with_relay_path < 0.85) {
            "fraction_A2_with_relay_path_below_0.85"
          } else {
            "M2_initial_top100_leakage_above_0.10"
          },
          n_A2_primary = rep$degree_qc$n_A2_primary,
          n_A2_extreme_high_degree = rep$degree_qc$n_A2_extreme_high_degree,
          fraction_A2_extreme_high_degree = rep$degree_qc$fraction_A2_extreme_high_degree,
          median_A2_degree_percentile = rep$degree_qc$median_A2_degree_percentile,
          action = "excluded_or_replaced",
          stringsAsFactors = FALSE)
        attempt <- attempt + 1
      }
      if (attempt >= length(templates) * 2 && isTRUE(severe_bad)) {
        outlier_rows[[length(outlier_rows) + 1]] <- data.frame(
          difficulty_setting = setting_name,
          requested_replicate = i,
          attempted_template_key = key,
          reason = "template_replacement_pool_exhausted",
          n_A2_primary = rep$degree_qc$n_A2_primary,
          n_A2_extreme_high_degree = rep$degree_qc$n_A2_extreme_high_degree,
          fraction_A2_extreme_high_degree = rep$degree_qc$fraction_A2_extreme_high_degree,
          median_A2_degree_percentile = rep$degree_qc$median_A2_degree_percentile,
          action = "retained_after_replacement_exhausted",
          stringsAsFactors = FALSE)
      }
      if (isTRUE(save_reps)) reps[[setting_name]][[i]] <- rep
      design_rows[[length(design_rows) + 1]] <- data.frame(difficulty_setting = setting_name, replicate = i,
                                                           template_key = key, seed = rep_seed,
                                                           stringsAsFactors = FALSE)
      path_rows[[length(path_rows) + 1]] <- transform(rep$path_qc, difficulty_setting = setting_name,
                                                       replicate = i, template_key = key,
                                                       n_A2_primary = length(rep$modules$A2_intermediate_degree_capped))
      direction_rows[[length(direction_rows) + 1]] <- transform(directional_qc_row_0628(rep),
                                                                difficulty_setting = setting_name)
      degree_rows[[length(degree_rows) + 1]] <- transform(rep$degree_qc, difficulty_setting = setting_name,
                                                          replicate = i, template_key = key)
      target_signal_rows[[length(target_signal_rows) + 1]] <- transform(rep$target_initial_signal_qc,
                                                                         difficulty_setting = setting_name,
                                                                         replicate = i, template_key = key)
      control_disruption_rows[[length(control_disruption_rows) + 1]] <- transform(control_disruption_audit_0628(rep),
                                                                                  difficulty_setting = setting_name,
                                                                                  replicate = i, template_key = key)
      branch_rows[[length(branch_rows) + 1]] <- transform(branch_specificity_audit_0628(rep),
                                                          difficulty_setting = setting_name,
                                                          replicate = i, template_key = key)
      conductance_rows[[length(conductance_rows) + 1]] <- transform(branch_conductance_audit_0630(rep),
                                                                    difficulty_setting = setting_name,
                                                                    replicate = i, template_key = key)
      relay_rows[[length(relay_rows) + 1]] <- transform(relay_structure_audit_0630(rep),
                                                        difficulty_setting = setting_name,
                                                        replicate = i, template_key = key)
      for (nm in names(rep$networks)) {
        m <- network_metrics(rep$networks[[nm]], rep$modules, rep$universe$mean_expression)
        m$replicate <- i; m$template_key <- key; m$network_label <- nm; m$difficulty_setting <- setting_name
        m$network_universe_size <- length(rep$universe$genes)
        m$background_sparsity <- rep$universe$background_sparsity
        m$n_A2_primary <- length(rep$modules$A2_intermediate_degree_capped)
        m$hard_global_pass <- with(m, density >= 0.0008 & density <= 0.22 &
                                     median_degree > 0 & median_strength > 0 &
                                     largest_connected_component_fraction >= 0.04 &
                                     module_a_connected_component_fraction >= 0.75 &
                                     isolated_fraction <= 0.97 & finite_edge_weights & !duplicated_gene_names)
        m$module_a_local_pass <- with(m, length(rep$modules$A) >= 30 &
                                        length(rep$modules$A2_intermediate_degree_capped) >= 12 &
                                        within_between_a_ratio >= 1.40 &
                                        is.finite(median_a1_a2_path) &
                                        median_a1_a2_path >= setting$median_min &
                                        median_a1_a2_path <= setting$median_max &
                                        a1_a2_same_component_fraction >= setting$same_component_min &
                                        direct_1hop_fraction >= setting$direct_min &
                                        direct_1hop_fraction <= setting$direct_max &
                                        two_hop_fraction >= setting$two_hop_min &
                                        two_hop_fraction <= setting$two_hop_max &
                                        three_plus_hop_fraction <= setting$three_plus_max &
                                        clustering_a > background_clustering_p40)
        qc_rows[[length(qc_rows) + 1]] <- m
      }
      for (nm in c("I1", "I2", "I3")) {
        q <- matching_qc(rep$networks$relevant, rep$networks[[nm]], rep$universe$mean_expression, rep$modules)
        q$replicate <- i; q$template_key <- key; q$network_label <- nm; q$difficulty_setting <- setting_name
        match_rows[[length(match_rows) + 1]] <- q
      }
    }
  }
  qc <- do.call(rbind, qc_rows)
  matching <- do.call(rbind, match_rows)
  path_qc <- do.call(rbind, path_rows)
  direction_qc <- do.call(rbind, direction_rows)
  degree_qc <- do.call(rbind, degree_rows)
  target_signal_qc <- do.call(rbind, target_signal_rows)
  control_disruption_qc <- do.call(rbind, control_disruption_rows)
  branch_qc <- do.call(rbind, branch_rows)
  conductance_qc <- do.call(rbind, conductance_rows)
  relay_qc <- do.call(rbind, relay_rows)
  matching$strict_matched_control <- matching$network_label %in% c("I2", "I3")
  matching$matching_pass <- with(matching, ifelse(strict_matched_control,
    density_diff <= 0.025 & median_strength_ratio >= 0.45 & median_strength_ratio <= 2.20 &
      expression_ks <= 0.22 & identical_twas & module_disrupted, TRUE))
  relevant <- qc$network_label == "relevant"
  rare_relevant <- relevant
  matched <- qc$network_label %in% c("relevant", "I2", "I3")
  path_summary_by_setting <- aggregate(cbind(direct_1hop_fraction, two_hop_fraction,
                                             three_plus_hop_fraction, median_A1_A2_path,
                                             A2_same_component_fraction) ~ difficulty_setting,
                                       path_qc, mean, na.rm = TRUE)
  pass_threshold <- function(setting, rare_min, extreme_min) rare_min
  setting_lookup <- settings[path_summary_by_setting$difficulty_setting]
  path_pass_vec <- vapply(seq_len(nrow(path_summary_by_setting)), function(i) {
    setting <- setting_lookup[[i]]
    with(path_summary_by_setting[i, ],
         is.finite(median_A1_A2_path) &&
           median_A1_A2_path >= setting$median_min &&
           median_A1_A2_path <= setting$median_max &&
           A2_same_component_fraction >= setting$same_component_min &&
           direct_1hop_fraction >= setting$direct_min &&
           direct_1hop_fraction <= setting$direct_max &&
           two_hop_fraction >= setting$two_hop_min &&
           two_hop_fraction <= setting$two_hop_max &&
           three_plus_hop_fraction <= setting$three_plus_max)
  }, logical(1))
  path_pass <- all(path_pass_vec)
  relevant_hard_global_pass_fraction <- mean(qc$hard_global_pass[relevant], na.rm = TRUE)
  rare_relevant_hard_global_pass_fraction <- mean(qc$hard_global_pass[rare_relevant], na.rm = TRUE)
  relevant_module_local_pass_fraction <- mean(qc$module_a_local_pass[relevant], na.rm = TRUE)
  rare_relevant_module_local_pass_fraction <- mean(qc$module_a_local_pass[rare_relevant], na.rm = TRUE)
  topology_pass <- relevant_hard_global_pass_fraction >= 0.95 &&
    rare_relevant_hard_global_pass_fraction >= 0.95 &&
    rare_relevant_module_local_pass_fraction >= 0.50 &&
    mean(matching$matching_pass[matching$strict_matched_control]) >= 0.85 &&
    path_pass
  directional_pass_by_setting <- aggregate(directional_qc_pass ~ difficulty_setting,
                                           direction_qc, mean, na.rm = TRUE)
  directional_pass <- all(directional_pass_by_setting$directional_qc_pass >=
                            pass_threshold(directional_pass_by_setting$difficulty_setting, 0.85, 0.80))
  degree_hard_pass_by_setting <- aggregate(degree_hard_gate_pass ~ difficulty_setting,
                                           degree_qc, mean, na.rm = TRUE)
  target_count_pass_by_setting <- aggregate(target_count_pilot_eligible ~ difficulty_setting,
                                            degree_qc, mean, na.rm = TRUE)
  degree_pass_by_setting <- merge(degree_hard_pass_by_setting, target_count_pass_by_setting,
                                  by = "difficulty_setting", all = TRUE)
  rare_degree_qc <- degree_qc
  degree_pass <- mean(degree_qc$degree_hard_gate_pass, na.rm = TRUE) >= 0.95 &&
    mean(degree_qc$target_count_pilot_eligible, na.rm = TRUE) >= 0.90 &&
    mean(rare_degree_qc$fraction_A2_extreme_high_degree, na.rm = TRUE) <= 0.15 &&
    mean(rare_degree_qc$median_A2_degree_percentile, na.rm = TRUE) <= 70 &&
    all(rare_degree_qc$fraction_A2_extreme_high_degree <= 0.40, na.rm = TRUE)
  target_signal_pass_by_setting <- aggregate(target_initial_signal_qc_pass ~ difficulty_setting,
                                             target_signal_qc, mean, na.rm = TRUE)
  target_signal_pass <- all(target_signal_pass_by_setting$target_initial_signal_qc_pass >=
                              pass_threshold(target_signal_pass_by_setting$difficulty_setting, 0.85, 0.80))
  control_disruption_pass_by_setting <- aggregate(control_disruption_qc_pass ~ difficulty_setting,
                                                  control_disruption_qc, mean, na.rm = TRUE)
  control_disruption_pass <- all(control_disruption_pass_by_setting$control_disruption_qc_pass >=
                                   pass_threshold(control_disruption_pass_by_setting$difficulty_setting, 0.85, 0.80))
  branch_pass_by_setting <- aggregate(branch_specificity_qc_pass ~ difficulty_setting,
                                      branch_qc, mean, na.rm = TRUE)
  branch_pass <- all(branch_pass_by_setting$branch_specificity_qc_pass >=
                       pass_threshold(branch_pass_by_setting$difficulty_setting, 0.85, 0.80))
  conductance_pass_by_setting <- aggregate(branch_conductance_qc_pass ~ difficulty_setting,
                                           conductance_qc, mean, na.rm = TRUE)
  conductance_pass <- all(conductance_pass_by_setting$branch_conductance_qc_pass >= 0.80)
  relay_pass_by_setting <- aggregate(relay_structure_qc_pass ~ difficulty_setting,
                                     relay_qc, mean, na.rm = TRUE)
  rare_relay_qc <- relay_qc
  relay_pass <- nrow(rare_relay_qc) > 0 &&
    mean(rare_relay_qc$n_relay_genes < 15, na.rm = TRUE) <= 0.15 &&
    all(rare_relay_qc$relay_raw_TWAS_top100_fraction <= 0.05, na.rm = TRUE) &&
    all(rare_relay_qc$relay_M2_initial_top100_fraction <= 0.05, na.rm = TRUE) &&
    all(rare_relay_qc$relay_background_mean_TOM <= 0.01, na.rm = TRUE) &&
    all(rare_relay_qc$relay_high_degree_decoy_TOM <= 0.01, na.rm = TRUE) &&
    all(rare_relay_qc$fraction_A2_with_relay_path >= 0.85, na.rm = TRUE)
  pass <- topology_pass && path_pass && directional_pass && degree_pass && target_signal_pass &&
    control_disruption_pass && branch_pass && conductance_pass && relay_pass
  atomic_write_csv(qc, file.path(out_dir, "topology_qc_metrics.csv"))
  atomic_write_csv(matching, file.path(out_dir, "relevant_irrelevant_matching_qc.csv"))
  atomic_write_csv(path_qc, file.path(out_dir, "PATH_STRATIFICATION_AUDIT.csv"))
  atomic_write_csv(direction_qc, file.path(out_dir, "DIRECTIONAL_SIGNAL_AUDIT.csv"))
  atomic_write_csv(degree_qc, file.path(out_dir, "DEGREE_DISTRIBUTION_AUDIT.csv"))
  atomic_write_csv(target_signal_qc, file.path(out_dir, "TARGET_INITIAL_SIGNAL_AUDIT.csv"))
  atomic_write_csv(control_disruption_qc, file.path(out_dir, "CONTROL_DISRUPTION_AUDIT.csv"))
  atomic_write_csv(branch_qc, file.path(out_dir, "BRANCH_SPECIFICITY_AUDIT.csv"))
  atomic_write_csv(conductance_qc, file.path(out_dir, "BRANCH_CONDUCTANCE_AUDIT.csv"))
  atomic_write_csv(relay_qc, file.path(out_dir, "RELAY_STRUCTURE_AUDIT.csv"))
  network_universe_audit <- qc[qc$network_label == "relevant",
                               intersect(c("difficulty_setting", "replicate", "template_key",
                                           "network_universe_size", "background_sparsity",
                                           "n_A2_primary", "density", "isolated_fraction",
                                           "largest_connected_component_fraction", "median_degree",
                                           "median_strength", "hard_global_pass", "module_a_local_pass"),
                                         names(qc)), drop = FALSE]
  atomic_write_csv(network_universe_audit, file.path(out_dir, "NETWORK_UNIVERSE_SIZE_AUDIT.csv"))
  outlier_tab <- if (length(outlier_rows)) do.call(rbind, outlier_rows) else data.frame(
    difficulty_setting = character(), requested_replicate = integer(),
    attempted_template_key = character(), reason = character(), n_A2_primary = integer(),
    n_A2_extreme_high_degree = integer(), fraction_A2_extreme_high_degree = numeric(),
    median_A2_degree_percentile = numeric(), action = character(), stringsAsFactors = FALSE)
  atomic_write_csv(outlier_tab, file.path(out_dir, "TEMPLATE_OUTLIER_AUDIT.csv"))
  atomic_write_csv(data.frame(
    pilot_go = pass,
    reason = if (pass) "topology_path_degree_direction_target_signal_control_branch_conductance_relay_qc_passed" else "topology_path_degree_direction_target_signal_control_branch_conductance_relay_qc_failed",
    topology_qc_pass = topology_pass, path_stratification_qc_pass = path_pass,
    degree_distribution_qc_pass = degree_pass, directional_qc_pass = directional_pass,
    target_initial_signal_qc_pass = target_signal_pass,
    control_disruption_qc_pass = control_disruption_pass,
    branch_specificity_qc_pass = branch_pass,
    branch_conductance_qc_pass = conductance_pass,
    relay_structure_qc_pass = relay_pass,
    template_count = length(templates), template_cell_types = length(unique(lib$templates$cell_type)),
    hard_global_pass_fraction = mean(qc$hard_global_pass[matched]),
    relevant_hard_global_pass_fraction = relevant_hard_global_pass_fraction,
    rare_relevant_hard_global_pass_fraction = rare_relevant_hard_global_pass_fraction,
    module_a_local_pass_fraction = relevant_module_local_pass_fraction,
    rare_relevant_module_local_pass_fraction = rare_relevant_module_local_pass_fraction,
    degree_hard_gate_pass_fraction = mean(degree_qc$degree_hard_gate_pass, na.rm = TRUE),
    rare_target_count_pilot_eligible_fraction = mean(degree_qc$target_count_pilot_eligible, na.rm = TRUE),
    rare_low_degree_warning_fraction = mean(degree_qc$median_A2_degree_percentile_warning_low, na.rm = TRUE),
    rare_branch_conductance_pass_fraction = mean(conductance_qc$branch_conductance_qc_pass, na.rm = TRUE),
    rare_relay_structure_pass_fraction = mean(relay_qc$relay_structure_qc_pass, na.rm = TRUE),
    i2_i3_matching_pass_fraction = mean(matching$matching_pass[matching$strict_matched_control]),
    path_fallback_used_fraction = mean(path_qc$path_fallback_used),
    i1_reported_nonblocking = TRUE), file.path(out_dir, "topology_qc_decision.csv"))
  design_index <- if (length(design_rows)) do.call(rbind, design_rows) else data.frame()
  atomic_write_csv(design_index, file.path(out_dir, "frozen_replicate_designs.csv"))
  if (isTRUE(save_reps)) atomic_save_rds(reps, file.path(out_dir, "frozen_replicate_designs.rds"))
  list(pass = pass, reason = if (pass) "topology_path_degree_direction_target_signal_control_branch_conductance_relay_qc_passed" else "topology_path_degree_direction_target_signal_control_branch_conductance_relay_qc_failed",
       topology_qc_pass = topology_pass, path_qc_pass = path_pass, degree_qc_pass = degree_pass,
       directional_qc_pass = directional_pass, target_signal_qc_pass = target_signal_pass,
       control_disruption_qc_pass = control_disruption_pass, branch_specificity_qc_pass = branch_pass,
       branch_conductance_qc_pass = conductance_pass, relay_structure_qc_pass = relay_pass,
       reps = reps, design_index = design_index)
}

candidate_audit_row_0705 <- function(candidate, topo, out_dir) {
  dec <- read.csv(file.path(out_dir, "topology_qc_decision.csv"))
  relay <- read.csv(file.path(out_dir, "RELAY_STRUCTURE_AUDIT.csv"))
  cond <- read.csv(file.path(out_dir, "BRANCH_CONDUCTANCE_AUDIT.csv"))
  # 0706: use observed structural summaries across the size-stratified conditions;
  # the historical column names are retained for report compatibility.
  rare_relay <- relay
  rare_cond <- cond
  score <- abs(mean(rare_relay$A1_relay_mean_TOM, na.rm = TRUE) - 0.2811) * 1000 +
    abs(mean(rare_relay$relay_A2_mean_TOM, na.rm = TRUE) - 0.1819) * 100 +
    abs(mean(rare_relay$median_A1_relay_A2_path_strength, na.rm = TRUE) - 0.8040) * 10 +
    abs(mean(rare_cond$A2_local_clustering, na.rm = TRUE) - 0.9056) +
    abs(mean(rare_cond$high_degree_bridge_count, na.rm = TRUE) - 0.0250) +
    abs(mean(rare_cond$opposite_sign_bridge_count, na.rm = TRUE) - 0)
  row <- data.frame(candidate_id = candidate$candidate_id,
             candidate_rank = candidate$candidate_rank,
             stage = if ("stage" %in% names(candidate)) candidate$stage else "stage1",
             parent_id = if ("parent_id" %in% names(candidate)) candidate$parent_id else "",
             a1_relay_center = candidate$a1_relay_center,
             relay_a2_center = candidate$relay_a2_center,
             a2_local_clustering_cap = candidate$a2_local_clustering_cap,
             opposite_suppression = candidate$opposite_suppression,
             relay_n_target = if ("relay_n_target" %in% names(candidate)) candidate$relay_n_target else NA_integer_,
             a1_relay_contact_n = if ("a1_relay_contact_n" %in% names(candidate)) candidate$a1_relay_contact_n else NA_integer_,
             relay_a2_contact_n = if ("relay_a2_contact_n" %in% names(candidate)) candidate$relay_a2_contact_n else NA_integer_,
             a1_relay_factor = if ("a1_relay_factor" %in% names(candidate)) candidate$a1_relay_factor else NA_real_,
             a1_relay_floor = if ("a1_relay_floor" %in% names(candidate)) candidate$a1_relay_floor else NA_real_,
             relay_a2_factor = if ("relay_a2_factor" %in% names(candidate)) candidate$relay_a2_factor else NA_real_,
             relay_a2_floor = if ("relay_a2_floor" %in% names(candidate)) candidate$relay_a2_floor else NA_real_,
             pilot_go = dec$pilot_go,
             topology_qc_pass = dec$topology_qc_pass,
             path_stratification_qc_pass = dec$path_stratification_qc_pass,
             degree_distribution_qc_pass = dec$degree_distribution_qc_pass,
             directional_qc_pass = dec$directional_qc_pass,
             target_initial_signal_qc_pass = dec$target_initial_signal_qc_pass,
             control_disruption_qc_pass = dec$control_disruption_qc_pass,
             branch_specificity_qc_pass = dec$branch_specificity_qc_pass,
             branch_conductance_qc_pass = dec$branch_conductance_qc_pass,
             relay_structure_qc_pass = dec$relay_structure_qc_pass,
             rare_relevant_module_local_pass_fraction = dec$rare_relevant_module_local_pass_fraction,
             degree_hard_gate_pass_fraction = dec$degree_hard_gate_pass_fraction,
             rare_target_count_pilot_eligible_fraction = dec$rare_target_count_pilot_eligible_fraction,
             rare_branch_conductance_pass_fraction = dec$rare_branch_conductance_pass_fraction,
             rare_relay_structure_pass_fraction = dec$rare_relay_structure_pass_fraction,
             rare_A1_relay_mean_TOM = mean(rare_relay$A1_relay_mean_TOM, na.rm = TRUE),
             rare_relay_A2_mean_TOM = mean(rare_relay$relay_A2_mean_TOM, na.rm = TRUE),
             rare_path_strength = mean(rare_relay$median_A1_relay_A2_path_strength, na.rm = TRUE),
             rare_relay_count = mean(rare_relay$n_relay_genes, na.rm = TRUE),
             rare_A2_local_clustering = mean(rare_cond$A2_local_clustering, na.rm = TRUE),
             rare_A_branch_background_cut_fraction = mean(rare_cond$A_branch_background_cut_fraction, na.rm = TRUE),
             rare_seed_neighborhood_background_fraction = mean(rare_cond$seed_neighborhood_background_fraction, na.rm = TRUE),
             rare_high_degree_bridge_count = mean(rare_cond$high_degree_bridge_count, na.rm = TRUE),
             rare_opposite_sign_bridge_count = mean(rare_cond$opposite_sign_bridge_count, na.rm = TRUE),
             closeness_to_1607_score = score,
             selected_candidate = FALSE,
             stringsAsFactors = FALSE)
  row$observed_hard_window_pass <- observed_window_pass_0705(row)
  row$observed_A1_relay_hard_pass <- row$rare_A1_relay_mean_TOM >= 0.260 && row$rare_A1_relay_mean_TOM <= 0.300
  row$observed_relay_A2_hard_pass <- row$rare_relay_A2_mean_TOM >= 0.170 && row$rare_relay_A2_mean_TOM <= 0.205
  row$observed_path_strength_hard_pass <- row$rare_path_strength >= 0.70 && row$rare_path_strength <= 0.85
  row$observed_relay_count_hard_pass <- row$rare_relay_count >= 24 && row$rare_relay_count <= 31
  row$observed_A2_clustering_hard_pass <- row$rare_A2_local_clustering >= 0.85 && row$rare_A2_local_clustering <= 0.97
  row
}

write_calibration_trace_0705 <- function(audit, selected = NA_character_, reason = "") {
  path <- project_file("results/calibration/CALIBRATION_ADAPTIVE_TRACE.md")
  if (file.exists(path)) unlink(path)
  lines <- c("# Calibration Adaptive Trace", "",
             "Calibration used structural and initial-signal metrics only. Final Heat ranking metrics were not loaded for candidate selection.", "",
             paste0("Selection reason: `", reason, "`."),
             paste0("Selected candidate: `", ifelse(is.na(selected), "none", selected), "`."),
             "",
             "Observed hard windows: A1-relay TOM 0.260-0.300; relay-A2 TOM 0.170-0.205; path strength 0.70-0.85; relay count 24-31; A2 local clustering 0.85-0.97; background cut <=0.12; seed-neighborhood background <=0.01; high-degree bridge <=0.15; opposite-sign bridge <=0.20; branch-conductance pass fraction >=0.80; relay-structure pass fraction >=0.95.",
             "", "## Candidate Summary", "")
  if (nrow(audit)) {
    for (i in seq_len(nrow(audit))) {
      r <- audit[i, , drop = FALSE]
      lines <- c(lines, paste0("- ", r$candidate_id, " [", r$stage, "] observed_pass=",
                               r$observed_hard_window_pass, "; pilot_go=", r$pilot_go,
                               "; A1-relay=", round(r$rare_A1_relay_mean_TOM, 4),
                               "; relay-A2=", round(r$rare_relay_A2_mean_TOM, 4),
                               "; path=", round(r$rare_path_strength, 4),
                               "; relay_n=", round(r$rare_relay_count, 3),
                               "; A2_cluster=", round(r$rare_A2_local_clustering, 4),
                               "; high_bridge=", round(r$rare_high_degree_bridge_count, 4),
                               "; opposite_bridge=", round(r$rare_opposite_sign_bridge_count, 4),
                               "; distance1607=", round(r$closeness_to_1607_score, 4)))
    }
  }
  atomic_write_lines(lines, path)
}

run_bounded_calibration_0705 <- function(lib, candidate_n_reps = 4, frozen_n_reps = 40, seed_base = 2027052012) {
  safe_dir_create(project_file("results/calibration"))
  stage1 <- calibration_candidates_0705()
  stage1$candidate_rank <- seq_len(nrow(stage1))
  candidates <- stage1
  audit_rows <- list(); selected <- NULL; reason <- "no_passing_calibration_candidate"
  evaluate_candidate <- function(cand_df, rank) {
    cand_df$candidate_rank <- rank
    cand <- as.list(cand_df[1, , drop = FALSE])
    set_calibration_candidate_0705(cand)
    set_run_status("IN_PROGRESS", "NO", "NO", paste0("0705 observed-metric adaptive calibration candidate ", cand$candidate_id))
    out_dir <- project_file("results", "calibration", cand$candidate_id)
    if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE)
    topo <- run_topology_qc_0628(lib, n_reps = candidate_n_reps, out_dir = out_dir,
                                 seed_base = seed_base + rank * 1000, save_reps = FALSE)
    row <- candidate_audit_row_0705(cand_df, topo, out_dir)
    list(row = row, topo = topo)
  }
  for (i in seq_len(nrow(stage1))) {
    res <- evaluate_candidate(stage1[i, , drop = FALSE], i)
    audit_rows[[length(audit_rows) + 1]] <- res$row
    if (isTRUE(res$row$pilot_go) && isTRUE(res$row$observed_hard_window_pass)) {
      selected <- res$row$candidate_id
      reason <- "stage1_observed_structural_candidate_selected"
      break
    }
  }
  audit <- do.call(rbind, audit_rows)
  if (is.null(selected)) {
    refinements <- adaptive_refinement_candidates_0705(audit, stage1)
    if (nrow(refinements) > 0) {
      refinements$candidate_rank <- seq_len(nrow(refinements)) + nrow(stage1)
      candidates <- rbind(stage1, refinements)
      for (j in seq_len(nrow(refinements))) {
        rank <- nrow(stage1) + j
        res <- evaluate_candidate(refinements[j, , drop = FALSE], rank)
        audit_rows[[length(audit_rows) + 1]] <- res$row
        audit <- do.call(rbind, audit_rows)
        if (isTRUE(res$row$pilot_go) && isTRUE(res$row$observed_hard_window_pass)) {
          selected <- res$row$candidate_id
          reason <- "stage2_observed_structural_candidate_selected"
          break
        }
      }
    }
  }
  audit <- do.call(rbind, audit_rows)
  if (!is.null(selected)) audit$selected_candidate[audit$candidate_id == selected] <- TRUE
  p <- project_file("results/calibration/CALIBRATION_CANDIDATE_AUDIT.csv")
  if (file.exists(p)) unlink(p)
  atomic_write_csv(audit, p)
  write_calibration_trace_0705(audit, selected = if (is.null(selected)) NA_character_ else selected, reason = reason)
  if (is.null(selected)) {
    return(list(pass = FALSE, reason = "no_passing_calibration_candidate", audit = audit,
                selected_candidate = NA_character_, topology = NULL))
  }
  cand_row <- candidates[candidates$candidate_id == selected, , drop = FALSE]
  set_calibration_candidate_0705(as.list(cand_row))
  final_dir <- project_file("results/topology_qc")
  if (dir.exists(final_dir)) unlink(final_dir, recursive = TRUE)
  safe_dir_create(final_dir)
  # Rerun the selected candidate into the canonical frozen topology directory.
  topology <- run_topology_qc_0628(lib, n_reps = frozen_n_reps, out_dir = final_dir,
                                   seed_base = seed_base, save_reps = FALSE)
  topology$calibration_candidate <- as.list(cand_row)
  list(pass = TRUE, reason = reason, audit = audit,
       selected_candidate = selected, topology = topology)
}

summarize_null_guardrails_0628 <- function(reps_flat, out_dir) {
  rows <- lapply(seq_along(reps_flat), function(i) {
    rep <- reps_flat[[i]]
    m2 <- faithful_m2_scores(rep$networks$relevant, rep$universe$mean_expression, rep$twas)
    ranked <- m2$SYMBOL[order(abs(m2$final_NESTA_heat), decreasing = TRUE)]
    top100 <- ranked[seq_len(min(100, length(ranked)))]
    deg <- degree_row_0628(transform(m2, score = abs(final_NESTA_heat)), rep, "score",
                           score_label_final_heat, "relevant", "NESTA_M2_faithful",
                           rep$rep_id, rep$template_key, rep$difficulty_setting, TRUE)
    data.frame(difficulty_setting = rep$difficulty_setting, replicate = rep$rep_id,
               opposite_sign_decoy_top100_rate = if (length(rep$modules$D_opposite_sign_decoy)) mean(rep$modules$D_opposite_sign_decoy %in% top100) else 0,
               C_high_degree_decoy_top100_rate = if (length(rep$modules$C_high_degree_decoy)) mean(rep$modules$C_high_degree_decoy %in% top100) else 0,
               score_degree_spearman = deg$score_degree_spearman,
               max_final_heat = max(abs(m2$final_NESTA_heat)), stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, rows)
  guard <- aggregate(cbind(opposite_sign_decoy_top100_rate, C_high_degree_decoy_top100_rate,
                           score_degree_spearman, max_final_heat) ~ difficulty_setting,
                     tab, mean, na.rm = TRUE)
  guard$empirical_fwe <- 0
  guard$guardrail_pass <- with(guard,
    opposite_sign_decoy_top100_rate <= ifelse(difficulty_setting == "rare_target_detection", 0.15, 0.25) &
      C_high_degree_decoy_top100_rate <= ifelse(difficulty_setting == "rare_target_detection", 0.12, 0.15) &
      abs(score_degree_spearman) <= ifelse(difficulty_setting == "rare_target_detection", 0.15, 0.20))
  atomic_write_csv(tab, file.path(out_dir, "null_bias_replicate_metrics.csv"))
  atomic_write_csv(guard, file.path(out_dir, "NULL_BIAS_GUARDRAILS.csv"))
  guard
}

run_pilot_0628 <- function(lib, reps, out_dir = project_file("results/pilot"),
                           null_dir = project_file("results/null"),
                           cores = 1,
                           null_seed_base = 2026070601) {
  safe_dir_create(out_dir); safe_dir_create(null_dir)
  settings <- difficulty_settings_0628()
  if (is.data.frame(reps)) {
    design <- reps
  } else {
    flat <- unlist(reps, recursive = FALSE)
    design <- do.call(rbind, lapply(flat, function(r) data.frame(
      difficulty_setting = r$difficulty_setting, replicate = r$rep_id,
      template_key = r$template_key, seed = NA_integer_, stringsAsFactors = FALSE)))
  }
  make_from_design <- function(d, null = FALSE, seed_override = NA_integer_) {
    setting <- settings[[d$difficulty_setting]]
    seed <- if (is.finite(seed_override)) seed_override else as.integer(d$seed)
    if (!is.finite(seed)) seed <- 2026070601 + match(d$difficulty_setting, names(settings)) * 100000 + as.integer(d$replicate)
    make_replicate_0628(lib, d$template_key, as.integer(d$replicate), seed = seed, setting = setting, null = null)
  }
  run_one <- function(i) {
    d <- design[i, , drop = FALSE]
    path <- file.path(out_dir, sprintf("pilot_%s_replicate_%03d_metrics.rds", d$difficulty_setting, as.integer(d$replicate)))
    if (file.exists(path)) return(readRDS(path))
    rep <- make_from_design(d, null = FALSE)
    x <- score_payload_0628(rep, rep$rep_id, null = FALSE, n_weight_perm = 10)
    x$reprioritization <- reprioritization_rows_0628(rep, x$rel_scores, rep$rep_id)
    x$delta <- delta_rows_0628(rep, x$rel_scores, rep$rep_id)
    x$diffusion_retention <- transform(diffusion_retention_audit_0628(rep, x$rel_scores),
                                       difficulty_setting = rep$difficulty_setting,
                                       replicate = rep$rep_id,
                                       template_key = rep$template_key,
                                       cell_type = sub("::.*", "", rep$template_key))
    atomic_save_rds(x, path)
    rm(rep); gc(FALSE)
    x
  }
  rows <- lapply(seq_len(nrow(design)), run_one)
  summary <- summarize_outputs_0628(rows, out_dir)
  null_rows <- list()
  for (i in seq_len(nrow(design))) {
    d <- design[i, , drop = FALSE]
    rep <- make_from_design(d, null = TRUE,
                            seed_override = null_seed_base + match(d$difficulty_setting, names(settings)) * 100000 + as.integer(d$replicate))
    m2 <- faithful_m2_scores(rep$networks$relevant, rep$universe$mean_expression, rep$twas)
    ranked <- m2$SYMBOL[order(abs(m2$final_NESTA_heat), decreasing = TRUE)]
    top100 <- ranked[seq_len(min(100, length(ranked)))]
    top5pct <- ranked[seq_len(min(max(1, ceiling(0.05 * length(ranked))), length(ranked)))]
    deg <- degree_row_0628(transform(m2, score = abs(final_NESTA_heat)), rep, "score",
                           score_label_final_heat, "relevant", "NESTA_M2_faithful",
                           rep$rep_id, rep$template_key, rep$difficulty_setting, TRUE)
    null_rows[[length(null_rows) + 1]] <- data.frame(
      difficulty_setting = rep$difficulty_setting, replicate = rep$rep_id,
      opposite_sign_decoy_top100_rate = if (length(rep$modules$D_opposite_sign_decoy)) mean(rep$modules$D_opposite_sign_decoy %in% top100) else 0,
      opposite_sign_decoy_top5pct_rate = if (length(rep$modules$D_opposite_sign_decoy)) mean(rep$modules$D_opposite_sign_decoy %in% top5pct) else 0,
      C_high_degree_decoy_top100_rate = if (length(rep$modules$C_high_degree_decoy)) mean(rep$modules$C_high_degree_decoy %in% top100) else 0,
      C_high_degree_decoy_top5pct_rate = if (length(rep$modules$C_high_degree_decoy)) mean(rep$modules$C_high_degree_decoy %in% top5pct) else 0,
      score_degree_spearman = deg$score_degree_spearman,
      max_final_heat = max(abs(m2$final_NESTA_heat)), stringsAsFactors = FALSE)
    rm(rep, m2); gc(FALSE)
  }
  tab <- do.call(rbind, null_rows)
  guard <- aggregate(cbind(opposite_sign_decoy_top100_rate, opposite_sign_decoy_top5pct_rate,
                           C_high_degree_decoy_top100_rate, C_high_degree_decoy_top5pct_rate,
                           score_degree_spearman, max_final_heat) ~ difficulty_setting,
                     tab, mean, na.rm = TRUE)
  guard$empirical_fwe <- 0
  guard$guardrail_pass <- with(guard,
    opposite_sign_decoy_top100_rate <= 0.15 & opposite_sign_decoy_top5pct_rate <= 0.15 &
      C_high_degree_decoy_top100_rate <= 0.12 & C_high_degree_decoy_top5pct_rate <= 0.12 &
      abs(score_degree_spearman) <= 0.20)
  atomic_write_csv(tab, file.path(null_dir, "null_bias_replicate_metrics.csv"))
  atomic_write_csv(guard, file.path(null_dir, "NULL_BIAS_GUARDRAILS.csv"))
  summary$null_guardrails <- guard
  summary
}

go_decision_0628 <- function(topology, pilot = NULL) {
  if (!isTRUE(topology$pass)) {
    return(list(status = "STOPPED", reason = topology$reason, pilot_started = "NO",
                confirmatory_started = "NO", stop_go = "STOP"))
  }
  if (is.null(pilot)) {
    return(list(status = "STOPPED", reason = "pilot_not_run", pilot_started = "NO",
                confirmatory_started = "NO", stop_go = "STOP"))
  }
  settings <- names(difficulty_settings_0628())
  directional <- pilot$direction[pilot$direction$score_name == score_label_final_heat &
                                   pilot$direction$network_label == "relevant", , drop = FALSE]
  degree <- pilot$degree[pilot$degree$score_name == score_label_final_heat &
                           pilot$degree$network_label == "relevant", , drop = FALSE]
  retention <- pilot$diffusion_retention
  dc <- pilot$direction_contrasts
  condition_pass <- vapply(settings, function(st) {
    d <- directional[directional$difficulty_setting == st, , drop = FALSE]
    g <- degree[degree$difficulty_setting == st, , drop = FALSE]
    r <- retention[retention$difficulty_setting == st, , drop = FALSE]
    cst <- dc[dc$difficulty_setting == st, , drop = FALSE]
    pos_contrast <- function(metric, cmp) {
      x <- cst[cst$metric == metric & cst$contrast == cmp, , drop = FALSE]
      nrow(x) == 1 && is.finite(x$mean) && is.finite(x$ci_low) && x$mean > 0 && x$ci_low > 0
    }
    nrow(d) >= 38 &&
      all(vapply(c("raw_signed_TWAS", "M2_no_diffusion", "PPR_signed_two_channel", "RWR_signed_two_channel"),
                 function(cmp) pos_contrast("direction_aware_AUPRC", cmp), logical(1))) &&
      all(vapply(c("M2_no_diffusion", "PPR_signed_two_channel", "RWR_signed_two_channel"),
                 function(cmp) pos_contrast("sign_concordant_top5pct_recall", cmp), logical(1))) &&
      mean(d$sign_concordant_top5pct_recall, na.rm = TRUE) >= difficulty_settings_0628()[[st]]$sign5pct_min &&
      mean(d$direction_accuracy_among_top5pct_A2, na.rm = TRUE) >= 0.75 &&
      mean(d$opposite_sign_decoy_top100_rate, na.rm = TRUE) <= 0.15 &&
      mean(d$opposite_sign_decoy_top5pct_rate, na.rm = TRUE) <= 0.15 &&
      nrow(r) >= 38 && mean(r$fraction_seed_heat_retained_in_A_branch, na.rm = TRUE) >= 0.10 &&
      mean(r$fraction_seed_heat_reaching_A2, na.rm = TRUE) >= 0.025 &&
      mean(r$fraction_seed_heat_leaking_to_background, na.rm = TRUE) <= 0.80 &&
      mean(r$relay_gene_top100_rate, na.rm = TRUE) <= 0.140 &&
      mean(r$relay_gene_top5pct_rate, na.rm = TRUE) <= 0.140 &&
      mean(r$relay_to_A2_heat_ratio, na.rm = TRUE) <= 0.95 &&
      nrow(g) >= 38 && abs(mean(g$score_degree_spearman, na.rm = TRUE)) <= 0.20 &&
      mean(g$C_high_degree_decoy_top100_rate, na.rm = TRUE) <= 0.12 &&
      mean(g$C_high_degree_decoy_top5pct_rate, na.rm = TRUE) <= 0.12
  }, logical(1))
  null_ok <- all(pilot$null_guardrails$guardrail_pass)
  ok <- all(condition_pass) && null_ok
  if (ok) {
    list(status = "GO", reason = "size_stratified_directional_pilot_go_criteria_passed_confirmatory_required",
         pilot_started = "YES", confirmatory_started = "NO", stop_go = "GO")
  } else {
    failed <- paste(names(condition_pass)[!condition_pass], collapse = ",")
    if (!nzchar(failed)) failed <- "null_bias_guardrails"
    list(status = "STOPPED", reason = paste0("pilot_go_criteria_failed:", failed),
         pilot_started = "YES", confirmatory_started = "NO", stop_go = "STOP")
  }
}

run_study_0628 <- function() {
  verify_project_path(); verify_binding_plan()
  set_run_status("IN_PROGRESS", "NO", "NO", "TOM topology library construction for 0706 size-stratified directional network recovery study")
  lib <- build_tom_topology_library()
  set_run_status("IN_PROGRESS", "NO", "NO", "0706 structural/initial-signal calibration search before size-stratified pilot")
  calibration <- run_bounded_calibration_0705(lib, candidate_n_reps = 1, frozen_n_reps = 40, seed_base = 2026070601)
  if (!isTRUE(calibration$pass)) {
    topology <- list(pass = FALSE, reason = calibration$reason,
                     topology_qc_pass = FALSE, path_qc_pass = FALSE, degree_qc_pass = FALSE,
                     directional_qc_pass = FALSE, target_signal_qc_pass = FALSE,
                     control_disruption_qc_pass = FALSE, branch_specificity_qc_pass = FALSE,
                     branch_conductance_qc_pass = FALSE, relay_structure_qc_pass = FALSE,
                     reps = list())
    decision <- list(status = "STOPPED", reason = "no_passing_calibration_candidate",
                     pilot_started = "NO", confirmatory_started = "NO", stop_go = "STOP")
    set_run_status(decision$status, decision$pilot_started, decision$confirmatory_started, decision$reason)
    atomic_save_rds(list(calibration = calibration, topology = topology, decision = decision), project_file("results/study_decision.rds"))
    return(invisible(decision))
  }
  topology <- calibration$topology
  if (!topology$pass) {
    decision <- go_decision_0628(topology, NULL)
    set_run_status(decision$status, decision$pilot_started, decision$confirmatory_started, decision$reason)
    atomic_save_rds(list(calibration = calibration, topology = topology, decision = decision), project_file("results/study_decision.rds"))
    return(invisible(decision))
  }
  set_run_status("IN_PROGRESS", "YES", "NO", paste0("40-replicate dense/basic/sparse pilot using calibration candidate ", calibration$selected_candidate))
  pilot <- run_pilot_0628(lib, topology$design_index, null_seed_base = 2026071601)
  decision <- go_decision_0628(topology, pilot)
  confirmatory <- NULL
  if (identical(decision$stop_go, "GO")) {
    set_run_status("IN_PROGRESS", "YES", "YES", "200-replicate confirmatory run per difficulty setting")
    confirm_topology <- run_topology_qc_0628(lib, n_reps = 200,
                                             out_dir = project_file("results/confirmatory_topology_qc"),
                                             seed_base = 2026072601,
                                             save_reps = FALSE)
    if (confirm_topology$pass) {
      confirmatory <- run_pilot_0628(lib, confirm_topology$design_index,
                                     out_dir = project_file("results/confirmatory"),
                                     null_dir = project_file("results/confirmatory_null"),
                                     null_seed_base = 2026073601)
      decision$confirmatory_started <- "YES"
      decision$reason <- "pilot_go_criteria_passed_confirmatory_completed"
    } else {
      decision$confirmatory_started <- "YES"
      decision$status <- "STOPPED"
      decision$stop_go <- "STOP"
      decision$reason <- "confirmatory_topology_path_degree_direction_qc_failed"
    }
  }
  set_run_status(decision$status, decision$pilot_started, decision$confirmatory_started, decision$reason)
  atomic_save_rds(list(calibration = calibration, topology = topology, pilot = pilot, confirmatory = confirmatory,
                       decision = decision), project_file("results/study_decision.rds"))
  invisible(decision)
}
