Sys.setenv(NESTA_DENSE_RESCUE_SOURCE_ONLY = "1",
           NESTA_BRANCH_RESCUE_SOURCE_ONLY = "1")
source(file.path("/home/js/NESTA/simulation", "R", "study_0707_branch_isolation.R"))

normalize_nonnegative <- function(x) {
  x[!is.finite(x)] <- 0
  s <- sum(x)
  if (!is.finite(s) || s <= 0) return(x)
  x / s
}

initial_vector_personalizations <- function(initial_weight) {
  initial_weight <- setNames(as.numeric(initial_weight), names(initial_weight))
  list(
    abs = normalize_nonnegative(abs(initial_weight)),
    pos = normalize_nonnegative(ifelse(initial_weight > 0, initial_weight, 0)),
    neg = normalize_nonnegative(ifelse(initial_weight < 0, abs(initial_weight), 0))
  )
}

score_strength <- function(adj) {
  x <- Matrix::rowSums(adj) - Matrix::diag(adj)
  setNames(as.numeric(x), rownames(adj))
}

score_degree <- function(adj) {
  x <- Matrix::rowSums(adj > 0.1) - 1
  setNames(as.numeric(x), rownames(adj))
}

degree_residualize <- function(score, degree) {
  ok <- is.finite(score) & is.finite(degree)
  out <- setNames(rep(NA_real_, length(score)), names(score))
  if (sum(ok) < 3 || stats::sd(score[ok]) == 0) return(out)
  out[ok] <- stats::resid(stats::lm(score[ok] ~ degree[ok]))
  out
}

audit_score_metrics <- function(rep, score, signed_score, score_name, source_object = "available") {
  genes <- rep$genes
  score <- setNames(as.numeric(score[genes]), genes)
  signed_score <- setNames(as.numeric(signed_score[genes]), genes)
  ranked <- ranked_genes(genes, score, TRUE, rep$A1)
  n_ranked <- length(ranked)
  top <- function(k) ranked[seq_len(min(k, n_ranked))]
  ranks <- rank_map(genes, score, rep$A1, TRUE)
  deg <- score_degree(rep$adj)[genes]
  str <- score_strength(rep$adj)[genes]
  au <- auprc_from_score(genes, score, rep$A2, rep$A1)
  prev <- length(rep$A2) / n_ranked
  dirs <- rep$A2_direction
  sv <- signed_score
  top100 <- top(100)
  concord <- intersect(top100, names(dirs))
  concord_ok <- if (!length(concord)) logical() else ifelse(dirs[concord] == "risk", sv[concord] > 0, sv[concord] < 0)
  risk_au <- auprc_from_score(genes, sv, rep$A2_risk, rep$A1)
  prot_au <- auprc_from_score(genes, -sv, rep$A2_protective, rep$A1)
  data.frame(
    arm = rep$arm, arm_label = rep$arm_label, condition = rep$condition,
    replicate = rep$rep_id, score_name = score_name, source_object = source_object,
    top50_recall = recall_at(ranked, rep$A2, 50),
    top100_recall = recall_at(ranked, rep$A2, 100),
    top150_recall = recall_at(ranked, rep$A2, 150),
    top200_recall = recall_at(ranked, rep$A2, 200),
    top10pct_recall = recall_at(ranked, rep$A2, ceiling(0.10 * n_ranked)),
    raw_AUPRC = au,
    prevalence_normalized_AUPRC = if (is.finite(au) && prev < 1) (au - prev) / (1 - prev) else NA_real_,
    direction_aware_AUPRC = weighted.mean(c(risk_au, prot_au), c(length(rep$A2_risk), length(rep$A2_protective)), na.rm = TRUE),
    sign_concordant_top100_recall = sum(concord_ok, na.rm = TRUE) / length(rep$A2),
    median_A2_rank = median(ranks[rep$A2], na.rm = TRUE),
    first_A2_rank = min(ranks[rep$A2], na.rm = TRUE),
    opposite_sign_decoy_top100_rate = mean(rep$D %in% top100),
    high_degree_decoy_top100_rate = mean(rep$C %in% top100),
    score_degree_spearman = suppressWarnings(stats::cor(score, deg, method = "spearman", use = "pairwise.complete.obs")),
    score_strength_spearman = suppressWarnings(stats::cor(score, str, method = "spearman", use = "pairwise.complete.obs")),
    stringsAsFactors = FALSE
  )
}

score_set_for_audit <- function(rep, ch) {
  genes <- rep$genes
  final <- setNames(ch$final_NESTA_heat, ch$SYMBOL)[genes]
  signed <- setNames(ch$signed_NESTA_heat, ch$SYMBOL)[genes]
  pos_raw <- setNames(ch$positive_channel_heat_raw, ch$SYMBOL)[genes]
  neg_raw <- setNames(ch$negative_channel_heat_raw, ch$SYMBOL)[genes]
  pos_scaled <- setNames(ch$positive_channel_heat, ch$SYMBOL)[genes]
  neg_scaled <- setNames(ch$negative_channel_heat, ch$SYMBOL)[genes]
  deg <- score_degree(rep$adj)[genes]
  raw_diff <- pos_raw - neg_raw
  list(
    submitted_output_descending = list(score = final, signed = final, available = TRUE,
                                       object = "Final.Heat/signed_NESTA_heat descending"),
    submitted_output_ascending = list(score = -final, signed = final, available = TRUE,
                                      object = "Final.Heat ascending"),
    abs_submitted_output_descending = list(score = abs(final), signed = final, available = TRUE,
                                           object = "abs(Final.Heat) descending"),
    raw_positive_channel_heat_descending = list(score = pos_raw, signed = pos_raw, available = TRUE,
                                                object = "raw positive diffusion channel"),
    raw_negative_channel_heat_descending = list(score = neg_raw, signed = -neg_raw, available = TRUE,
                                                object = "raw absolute-negative diffusion channel"),
    signed_reconstructed_heat_descending = list(score = signed, signed = signed, available = TRUE,
                                                object = "positive channel minus negative channel descending"),
    signed_reconstructed_heat_ascending = list(score = -signed, signed = signed, available = TRUE,
                                               object = "positive channel minus negative channel ascending"),
    abs_signed_reconstructed_heat_descending = list(score = abs(signed), signed = signed, available = TRUE,
                                                    object = "abs(reconstructed signed heat) descending"),
    positive_minus_negative_descending = list(score = pos_scaled - neg_scaled, signed = pos_scaled - neg_scaled,
                                             available = TRUE, object = "scaled positive minus scaled negative"),
    positive_plus_abs_negative_descending = list(score = pos_scaled + abs(neg_scaled), signed = pos_scaled - neg_scaled,
                                                available = TRUE, object = "scaled positive plus absolute negative"),
    degree_residualized_submitted_output = list(score = degree_residualize(final, deg), signed = final,
                                                available = TRUE, object = "lm residual Final.Heat ~ degree"),
    degree_residualized_raw_heat = list(score = degree_residualize(raw_diff, deg), signed = raw_diff,
                                        available = TRUE, object = "lm residual raw positive-minus-negative ~ degree")
  )
}

unavailable_score_rows <- function(rep) {
  data.frame(
    arm = rep$arm, arm_label = rep$arm_label, condition = rep$condition,
    replicate = rep$rep_id,
    score_name = c("negative_log10_submitted_p_if_pvalue",
                   "negative_log10_empirical_p_if_available",
                   "diffusion_z_score_if_available"),
    source_object = "not_available",
    reason = c("Submitted non-grid output is Final.Heat, not a P-value object.",
               "diffuStats permutation internals are not returned by the faithful submitted call.",
               "No separate diffusion z-score object is exported; submitted Final.Heat is sample-SD scaled heat."),
    stringsAsFactors = FALSE
  )
}

benchmark_same_initial <- function(rep, ch) {
  genes <- rep$genes
  init <- setNames(ch$initial_weight, ch$SYMBOL)[genes]
  p <- initial_vector_personalizations(init)
  rwr_abs <- run_rwr(rep$adj, p$abs, restart = 0.5)
  ppr_abs <- run_ppr(rep$adj, p$abs, damping = 0.85)
  rwr_signed <- run_rwr(rep$adj, p$pos, restart = 0.5) - run_rwr(rep$adj, p$neg, restart = 0.5)
  ppr_signed <- run_ppr(rep$adj, p$pos, damping = 0.85) - run_ppr(rep$adj, p$neg, damping = 0.85)
  list(
    RWR_abs_prior = list(score = rwr_abs, signed = rwr_abs),
    PPR_abs_prior = list(score = ppr_abs, signed = ppr_abs),
    RWR_signed_two_channel = list(score = abs(rwr_signed), signed = rwr_signed),
    PPR_signed_two_channel = list(score = abs(ppr_signed), signed = ppr_signed)
  )
}

alignment_rows <- function(rep, nesta_scores, benchmarks) {
  genes <- rep$genes
  branch <- unique(c(rep$A1, rep$relay, rep$A2))
  a2relay <- unique(c(rep$A2, rep$relay))
  deg <- score_degree(rep$adj)[genes]
  out <- list(); i <- 1
  for (sn in names(nesta_scores)) {
    s <- setNames(nesta_scores[[sn]]$score[genes], genes)
    s_ranked <- ranked_genes(genes, s, TRUE, rep$A1)
    s_top100 <- s_ranked[seq_len(min(100, length(s_ranked)))]
    s_a2_rank <- rank_map(genes, s, rep$A1)
    s_deg <- suppressWarnings(stats::cor(s, deg, method = "spearman", use = "pairwise.complete.obs"))
    for (bn in names(benchmarks)) {
      b <- setNames(benchmarks[[bn]]$score[genes], genes)
      b_ranked <- ranked_genes(genes, b, TRUE, rep$A1)
      b_top100 <- b_ranked[seq_len(min(100, length(b_ranked)))]
      b_a2_rank <- rank_map(genes, b, rep$A1)
      b_deg <- suppressWarnings(stats::cor(b, deg, method = "spearman", use = "pairwise.complete.obs"))
      out[[i]] <- data.frame(
        arm = rep$arm, arm_label = rep$arm_label, replicate = rep$rep_id,
        nesta_score_name = sn, benchmark_score_name = bn,
        all_gene_spearman = suppressWarnings(stats::cor(s, b, method = "spearman", use = "pairwise.complete.obs")),
        A_branch_spearman = suppressWarnings(stats::cor(s[branch], b[branch], method = "spearman", use = "pairwise.complete.obs")),
        A2_relay_spearman = suppressWarnings(stats::cor(s[a2relay], b[a2relay], method = "spearman", use = "pairwise.complete.obs")),
        top100_overlap = length(intersect(s_top100, b_top100)) / 100,
        median_A2_rank_nesta = median(s_a2_rank[rep$A2], na.rm = TRUE),
        median_A2_rank_benchmark = median(b_a2_rank[rep$A2], na.rm = TRUE),
        first_A2_rank_nesta = min(s_a2_rank[rep$A2], na.rm = TRUE),
        first_A2_rank_benchmark = min(b_a2_rank[rep$A2], na.rm = TRUE),
        degree_spearman_difference = s_deg - b_deg,
        stringsAsFactors = FALSE
      )
      i <- i + 1
    }
  }
  do.call(rbind, out)
}

seed_vector_row <- function(rep, ch, cutoff = 1) {
  genes <- rep$genes
  raw_z <- setNames(rep$twas$TWAS.Z, rep$twas$SYMBOL)[genes]
  raw_p <- twas_p_from_z(raw_z)
  filt_z <- strict_twas_vector(genes, rep$twas, cutoff = cutoff)
  init <- setNames(ch$initial_weight, ch$SYMBOL)[genes]
  initial_pre_filter <- (rep$expr / safe_sd(rep$expr)) * (raw_z / safe_sd(raw_z))
  data.frame(
    arm = rep$arm, arm_label = rep$arm_label, replicate = rep$rep_id,
    total_genes = length(genes),
    A1_seed_genes = length(rep$A1),
    nonzero_initial_weights_before_filtering = sum(abs(initial_pre_filter) > 0),
    nonzero_initial_weights_after_filtering = sum(abs(init) > 0),
    A1_nonzero_after_filtering = sum(abs(init[rep$A1]) > 0),
    A2_nonzero_after_filtering = sum(abs(init[rep$A2]) > 0),
    mean_abs_initial_weight_A1 = mean(abs(init[rep$A1])),
    mean_abs_initial_weight_A2 = mean(abs(init[rep$A2])),
    mean_abs_initial_weight_background = mean(abs(init[rep$background])),
    TWAS_P_threshold_used = cutoff,
    fraction_A1_removed_by_filtering = mean(raw_p[rep$A1] >= cutoff | filt_z[rep$A1] == 0),
    fraction_A2_removed_by_filtering = mean(raw_p[rep$A2] >= cutoff | filt_z[rep$A2] == 0),
    A1_initial_nonzero_pass = sum(abs(init[rep$A1]) > 0) >= ceiling(0.8 * length(rep$A1)),
    A2_initial_weak = mean(rank(-abs(init), ties.method = "average")[match(rep$A2, genes)] <= 100) <= 0.05,
    stringsAsFactors = FALSE
  )
}

operator_row <- function(rep, ch) {
  g <- make_faithful_graph(rep$adj, edge_cutoff = 0.1, retain_zero_edges = TRUE,
                           submitted_diagonal = TRUE)
  data.frame(
    arm = rep$arm, arm_label = rep$arm_label, replicate = rep$rep_id,
    n_vertices = igraph::vcount(g),
    n_edges_including_self_loops_and_zero_weight_retained = igraph::ecount(g),
    n_self_loops = sum(igraph::which_loop(g)),
    n_zero_weight_edges_retained = sum(igraph::E(g)$weight == 0),
    n_positive_initial_seed_weights = sum(ch$initial_weight > 0),
    n_negative_initial_seed_weights = sum(ch$initial_weight < 0),
    final_heat_sd = stats::sd(ch$final_NESTA_heat),
    positive_channel_raw_sd = stats::sd(ch$positive_channel_heat_raw),
    negative_channel_raw_sd = stats::sd(ch$negative_channel_heat_raw),
    stringsAsFactors = FALSE
  )
}

make_positive_control_rep <- function(pc_id, rep_id = 1) {
  set.seed(990000 + rep_id + ifelse(pc_id == "PC1", 0, 1000))
  n <- if (pc_id == "PC1") 200 else 1000
  genes <- sprintf("%s_G%04d", pc_id, seq_len(n))
  a1 <- genes[1:10]
  relay <- genes[11:30]
  a2 <- genes[31:50]
  decoy <- genes[51:80]
  cdecoy <- genes[81:110]
  bg <- genes[111:n]
  risk_a1 <- a1[1:5]; prot_a1 <- a1[6:10]
  risk_relay <- relay[seq(1, length(relay), 2)]
  prot_relay <- setdiff(relay, risk_relay)
  risk_a2 <- a2[seq(1, length(a2), 2)]
  prot_a2 <- setdiff(a2, risk_a2)
  ii <- integer(); jj <- integer(); xx <- numeric()
  add <- function(a, b, w) {
    a <- rep(a, length.out = max(length(a), length(b), length(w)))
    b <- rep(b, length.out = length(a)); w <- rep(w, length.out = length(a))
    keep <- a != b & a %in% genes & b %in% genes
    ia <- match(a[keep], genes); ib <- match(b[keep], genes); ww <- w[keep]
    ii <<- c(ii, ia, ib); jj <<- c(jj, ib, ia); xx <<- c(xx, ww, ww)
  }
  for (r in risk_relay) add(r, sample(risk_a1, 3), runif(3, 0.35, 0.45))
  for (r in prot_relay) add(r, sample(prot_a1, 3), runif(3, 0.35, 0.45))
  for (g in risk_a2) add(g, sample(risk_relay, 5), runif(5, 0.30, 0.42))
  for (g in prot_a2) add(g, sample(prot_relay, 5), runif(5, 0.30, 0.42))
  cmb <- combn(risk_a2, 2); for (k in seq_len(ncol(cmb))) if (runif(1) < 0.35) add(cmb[1, k], cmb[2, k], runif(1, 0.16, 0.24))
  cmb <- combn(prot_a2, 2); for (k in seq_len(ncol(cmb))) if (runif(1) < 0.35) add(cmb[1, k], cmb[2, k], runif(1, 0.16, 0.24))
  if (length(bg) > 1) add(sample(bg, min(500, length(bg) * 2), TRUE), sample(bg, min(500, length(bg) * 2), TRUE), runif(min(500, length(bg) * 2), 0.001, 0.035))
  if (pc_id == "PC2") add(sample(c(a1, relay, a2), 8, TRUE), sample(bg, 8, TRUE), runif(8, 0.001, 0.010))
  adj <- sparseMatrix(i = ii, j = jj, x = xx, dims = c(n, n), dimnames = list(genes, genes))
  adj <- forceSymmetric(adj, uplo = "U"); diag(adj) <- 1
  expr <- setNames(rlnorm(n, 0, 0.25), genes)
  expr[a2] <- stats::median(expr) * runif(length(a2), 0.95, 1.05)
  z <- setNames(rnorm(n, 0, 0.55), genes)
  z[risk_a1] <- abs(rnorm(length(risk_a1), 3.7, 0.30))
  z[prot_a1] <- -abs(rnorm(length(prot_a1), 3.7, 0.30))
  z[relay] <- rnorm(length(relay), 0, 0.08)
  z[risk_a2] <- rnorm(length(risk_a2), 0.12, 0.08)
  z[prot_a2] <- rnorm(length(prot_a2), -0.12, 0.08)
  z[decoy] <- sample(c(-1, 1), length(decoy), TRUE) * abs(rnorm(length(decoy), 0.9, 0.25))
  twas <- data.frame(SYMBOL = genes, TWAS.Z = as.numeric(z), TWAS.P = twas_p_from_z(z),
                     stringsAsFactors = FALSE)
  list(condition = "dense_positive_control", arm = pc_id, arm_label = ifelse(pc_id == "PC1", "tiny_200_minimal_path", "explicit_A_component_1000"),
       rep_id = rep_id, genes = genes, adj = adj, expr = expr, twas = twas,
       A1 = a1, A1_risk = risk_a1, A1_protective = prot_a1,
       relay = relay, relay_risk = risk_relay, relay_protective = prot_relay,
       relay_direction = c(setNames(rep("risk", length(risk_relay)), risk_relay),
                           setNames(rep("protective", length(prot_relay)), prot_relay)),
       A2 = a2, A2_risk = risk_a2, A2_protective = prot_a2,
       A2_direction = c(setNames(rep("risk", length(risk_a2)), risk_a2),
                        setNames(rep("protective", length(prot_a2)), prot_a2)),
       D = decoy, C = cdecoy, background = bg)
}

run_one_operator_replicate <- function(rep, n.perm = 25, seed = 1) {
  ch <- faithful_m2_channels(rep$adj, rep$expr, rep$twas, n.perm = n.perm, seed = seed)
  ss <- score_set_for_audit(rep, ch)
  bm <- benchmark_same_initial(rep, ch)
  metrics <- do.call(rbind, lapply(names(ss), function(nm) audit_score_metrics(rep, ss[[nm]]$score, ss[[nm]]$signed, nm, ss[[nm]]$object)))
  bench_metrics <- do.call(rbind, lapply(names(bm), function(nm) audit_score_metrics(rep, bm[[nm]]$score, bm[[nm]]$signed, nm, "same_initial_vector_benchmark")))
  list(ch = ch, score_metrics = metrics, benchmark_metrics = bench_metrics,
       unavailable = unavailable_score_rows(rep), alignment = alignment_rows(rep, ss, bm),
       seed = seed_vector_row(rep, ch), operator = operator_row(rep, ch))
}

summarise_score_metrics <- function(x) {
  stats::aggregate(cbind(top50_recall, top100_recall, top150_recall, top200_recall,
                         top10pct_recall, raw_AUPRC, prevalence_normalized_AUPRC,
                         direction_aware_AUPRC, sign_concordant_top100_recall,
                         median_A2_rank, first_A2_rank, opposite_sign_decoy_top100_rate,
                         high_degree_decoy_top100_rate, score_degree_spearman,
                         score_strength_spearman) ~ arm + arm_label + score_name,
                   data = x, FUN = function(z) mean(z, na.rm = TRUE))
}

run_operator_audit <- function() {
  verify_project_path()
  verify_binding_plan()
  report_dir <- read_report_dir()
  safe_dir_create(report_dir)
  copy_binding_plan_to_report(report_dir)
  set_run_status("IN_PROGRESS", "YES", "NO", "0707 NESTA operator ranking audit")

  metric_rows <- list(); bench_rows <- list(); unavailable_rows <- list()
  align_rows <- list(); seed_rows <- list(); operator_rows <- list()
  i <- ib <- iu <- ia <- is <- io <- 1
  for (arm in c("F", "H")) for (rep_id in seq_len(20)) {
    rep <- make_branch_isolation_rep(arm, rep_id, 990700 + match(arm, c("F", "H")) * 1000 + rep_id)
    got <- run_one_operator_replicate(rep, n.perm = 25, seed = 990700 + rep_id)
    metric_rows[[i]] <- got$score_metrics; i <- i + 1
    bench_rows[[ib]] <- got$benchmark_metrics; ib <- ib + 1
    unavailable_rows[[iu]] <- got$unavailable; iu <- iu + 1
    align_rows[[ia]] <- got$alignment; ia <- ia + 1
    seed_rows[[is]] <- got$seed; is <- is + 1
    operator_rows[[io]] <- got$operator; io <- io + 1
    if (rep_id %% 5 == 0) gc(FALSE)
  }

  pc_metrics <- list(); pc_i <- 1
  for (pc in c("PC1", "PC2")) {
    rep <- make_positive_control_rep(pc, 1)
    got <- run_one_operator_replicate(rep, n.perm = 25, seed = 992000 + pc_i)
    pc_metrics[[pc_i]] <- rbind(got$score_metrics, got$benchmark_metrics)
    pc_i <- pc_i + 1
  }

  metrics <- do.call(rbind, metric_rows)
  bench <- do.call(rbind, bench_rows)
  unavailable <- do.call(rbind, unavailable_rows)
  align <- do.call(rbind, align_rows)
  seed <- do.call(rbind, seed_rows)
  operator <- do.call(rbind, operator_rows)
  pc <- do.call(rbind, pc_metrics)

  score_summary <- summarise_score_metrics(metrics)
  bench_summary <- summarise_score_metrics(bench)
  pc_summary <- summarise_score_metrics(pc)
  output_object <- data.frame(
    object_name = c("submitted_Final.Heat", "positive_channel_heat_raw", "negative_channel_heat_raw",
                    "positive_channel_heat_scaled", "negative_channel_heat_scaled",
                    "negative_log10_submitted_p_if_pvalue", "empirical_p_if_available", "diffusion_z_score_if_available"),
    availability = c(rep("available", 5), rep("not_available", 3)),
    interpretation = c("Submitted `Final.Heat` equals sample-SD scaled positive-channel heat minus sample-SD scaled absolute-negative-channel heat, scaled again by sample SD.",
                       "Raw output of submitted positive diffusion channel before submitted SD scaling.",
                       "Raw output of submitted absolute-negative diffusion channel before submitted SD scaling.",
                       "Positive channel divided by sample SD as in submitted recombination.",
                       "Negative channel divided by sample SD as in submitted recombination.",
                       "Not produced by non-grid submitted `Nesta.R` output.",
                       "Permutation internals are not returned as a rankable object by the faithful call.",
                       "No separate exported z-score object; `Final.Heat` is an SD-scaled heat score."),
    larger_better = c("signed positive orientation; audit also tests ascending and absolute", "yes for positive branch heat",
                      "yes for protective/negative channel magnitude", "yes for positive branch heat",
                      "yes for protective/negative channel magnitude", "not_available", "not_available", "not_available"),
    stringsAsFactors = FALSE
  )
  permutation_audit <- unique(unavailable[, c("arm", "arm_label", "replicate")])
  permutation_audit$submitted_permutation_score_object <- "not_available"
  permutation_audit$note <- "Submitted non-grid Final.Heat path calls diffuStats but does not expose empirical permutation P-values as a rankable object."

  write_csv_over(output_object, file.path(report_dir, "NESTA_OUTPUT_OBJECT_AUDIT.csv"))
  write_csv_over(metrics, file.path(report_dir, "NESTA_RANK_ORIENTATION_AUDIT.csv"))
  unavailable_metrics <- metrics[0, , drop = FALSE]
  unavailable_metrics[seq_len(nrow(unavailable)), ] <- NA
  for (nm in intersect(names(unavailable), names(unavailable_metrics))) {
    unavailable_metrics[[nm]] <- unavailable[[nm]]
  }
  unavailable_metrics$source_object <- unavailable$source_object
  write_csv_over(rbind(metrics, unavailable_metrics),
                 file.path(report_dir, "NESTA_INTERMEDIATE_SCORE_AUDIT.csv"))
  write_csv_over(seed, file.path(report_dir, "NESTA_SEED_VECTOR_AUDIT.csv"))
  write_csv_over(operator, file.path(report_dir, "NESTA_DIFFUSION_OPERATOR_AUDIT.csv"))
  write_csv_over(align, file.path(report_dir, "NESTA_PPR_ALIGNMENT_AUDIT.csv"))
  write_csv_over(permutation_audit, file.path(report_dir, "NESTA_PERMUTATION_SCORE_AUDIT.csv"))
  write_csv_over(pc, file.path(report_dir, "MINIMAL_POSITIVE_CONTROL_METRICS.csv"))
  write_csv_over(metrics, file.path(report_dir, "PRIMARY_FINAL_HEAT_METRICS.csv"))
  write_csv_over(metrics[, c("arm", "arm_label", "condition", "replicate", "score_name",
                             "direction_aware_AUPRC", "sign_concordant_top100_recall",
                             "opposite_sign_decoy_top100_rate")],
                 file.path(report_dir, "DIRECTION_AWARE_METRICS.csv"))
  write_csv_over(bench, file.path(report_dir, "BENCHMARK_METRICS.csv"))
  write_csv_over(bench_summary, file.path(report_dir, "BENCHMARK_CONTRASTS.csv"))
  schema <- data.frame(
    file = c("PRIMARY_FINAL_HEAT_METRICS.csv", "DIRECTION_AWARE_METRICS.csv", "BENCHMARK_METRICS.csv"),
    required_columns_present = c(all(c("top50_recall", "top100_recall", "top150_recall", "top200_recall", "top10pct_recall", "raw_AUPRC", "prevalence_normalized_AUPRC") %in% names(metrics)),
                                 all(c("direction_aware_AUPRC", "sign_concordant_top100_recall", "opposite_sign_decoy_top100_rate") %in% names(metrics)),
                                 all(c("top100_recall", "raw_AUPRC", "score_degree_spearman") %in% names(bench))),
    stringsAsFactors = FALSE
  )
  write_csv_over(schema, file.path(report_dir, "METRIC_SCHEMA_AUDIT.csv"))

  if (file.exists(project_file("results/reports/summary_tables/unit_test_results.csv"))) {
    file.copy(project_file("results/reports/summary_tables/unit_test_results.csv"),
              file.path(report_dir, "unit_test_results.csv"), overwrite = TRUE)
  }

  seed_failure <- mean(seed$fraction_A1_removed_by_filtering > 0.5, na.rm = TRUE) > 0
  pc_nesta <- pc_summary[!grepl("PPR|RWR", pc_summary$score_name), , drop = FALSE]
  pc_success <- any(pc_nesta$top100_recall >= 0.05 | pc_nesta$top200_recall >= 0.15, na.rm = TRUE)
  arm_success_rows <- score_summary[score_summary$top100_recall >= 0.05 | score_summary$top200_recall >= 0.15, , drop = FALSE]
  promising <- nrow(arm_success_rows) > 0
  best <- score_summary[order(score_summary$top100_recall, score_summary$top200_recall,
                              -score_summary$median_A2_rank, decreasing = TRUE), ][1, ]
  best_bench <- bench_summary[order(bench_summary$top100_recall, bench_summary$top200_recall,
                                    -bench_summary$median_A2_rank, decreasing = TRUE), ][1, ]
  if (seed_failure) {
    reason <- "seed_filtering_failure"
  } else if (!pc_success) {
    reason <- "minimal_positive_control_failure"
  } else if (promising && any(grepl("ascending", arm_success_rows$score_name))) {
    reason <- "ranking_orientation_error_found"
  } else if (promising) {
    reason <- "intermediate_score_rescues_A2"
  } else if (pc_success) {
    reason <- "full_network_normalization_failure"
  } else {
    reason <- "operator_audit_no_rescue"
  }
  set_run_status("STOPPED", "YES", "NO", reason)

  write_lines_over(c("# NESTA Output Object Audit", "",
                     "The submitted non-grid object ranked in prior runs is `Final.Heat` from the submitted score table.",
                     "In `Nesta.R`, `Final.Heat` is positive-channel diffusion heat divided by sample SD minus absolute-negative-channel diffusion heat divided by sample SD, then divided by sample SD again.",
                     "The submitted non-grid output is not an empirical P-value table, adjusted P-value table, or separately exported z-score table.",
                     "Positive and negative channels are recombined as `F.score.Pos/sd(F.score.Pos) - F.score.Neg/sd(F.score.Neg)` followed by another SD scaling."),
                   file.path(report_dir, "NESTA_OUTPUT_OBJECT_AUDIT.md"))
  write_lines_over(c("# NESTA Rank Orientation Audit", "",
                     sprintf("Best NESTA-derived Arm F/H score: `%s` on arm `%s`.", best$score_name, best$arm),
                     sprintf("Mean top100/top150/top200 recall: %.4f / %.4f / %.4f.",
                             best$top100_recall, best$top150_recall, best$top200_recall),
                     sprintf("Median A2 rank: %.2f; first A2 rank: %.2f.",
                             best$median_A2_rank, best$first_A2_rank),
                     sprintf("Promising audit result threshold met: %s.", promising)),
                   file.path(report_dir, "NESTA_RANK_ORIENTATION_AUDIT.md"))
  write_lines_over(c("# NESTA PPR Alignment Audit", "",
                     "PPR/RWR benchmarks were computed on the exact same graph and exact same faithful M2 initial vector used by submitted NESTA.",
                     sprintf("Best benchmark score: `%s` on arm `%s`, top100/top200 recall %.4f / %.4f.",
                             best_bench$score_name, best_bench$arm, best_bench$top100_recall, best_bench$top200_recall),
                     "All-gene, A-branch, A2/relay Spearman correlations, top100 overlap, A2 rank summaries, and degree-dependence differences are in NESTA_PPR_ALIGNMENT_AUDIT.csv."),
                   file.path(report_dir, "NESTA_PPR_ALIGNMENT_AUDIT.md"))
  write_lines_over(c("# NESTA Seed Vector Audit", "",
                     sprintf("Mean A1 removed by TWAS.P filtering: %.4f.", mean(seed$fraction_A1_removed_by_filtering)),
                     sprintf("Mean A1 nonzero after filtering: %.2f of %.2f.",
                             mean(seed$A1_nonzero_after_filtering), mean(seed$A1_seed_genes)),
                     sprintf("Mean A2 nonzero after filtering: %.2f.", mean(seed$A2_nonzero_after_filtering)),
                     "The audit table reports total genes, nonzero initial weights before/after filtering, A1/A2 counts, mean absolute initial weights, TWAS.P threshold, and removal fractions."),
                   file.path(report_dir, "NESTA_SEED_VECTOR_AUDIT.md"))
  write_lines_over(c("# Minimal Positive Control Report", "",
                     sprintf("Minimal positive control success: %s.", pc_success),
                     "PC1 is a 200-gene A1-relay-A2 graph with minimal background.",
                     "PC2 is a 1,000-gene graph with an explicit A-only connected component and tiny controlled leakage.",
                     "Per-score recovery metrics are in MINIMAL_POSITIVE_CONTROL_METRICS.csv."),
                   file.path(report_dir, "MINIMAL_POSITIVE_CONTROL_REPORT.md"))
  write_lines_over(c("# NESTA Operator Audit Summary", "",
                     paste0("Final outcome classification: `", reason, "`."),
                     sprintf("Best NESTA-derived Arm F/H score `%s`: top100/top200 %.4f / %.4f.",
                             best$score_name, best$top100_recall, best$top200_recall),
                     sprintf("Best same-initial benchmark `%s`: top100/top200 %.4f / %.4f.",
                             best_bench$score_name, best_bench$top100_recall, best_bench$top200_recall),
                     sprintf("Promising audit result for next dense validation round: %s.", promising)),
                   file.path(report_dir, "NESTA_OPERATOR_AUDIT_SUMMARY.md"))
  write_lines_over(c("# Code Fidelity Audit", "",
                     paste0("Binding plan SHA256: `", binding_plan_sha256, "`."),
                     "Faithful submitted M2 arithmetic uses SCT row-mean expression, uncentered sample-SD scaling, strict TWAS.P filtering, signed positive and absolute-negative channels, retained zero-weight edge and self-loop behavior, and diffuStats `n.perm`.",
                     "No submitted NESTA method modification was made; this round ranks submitted and directly available intermediate outputs only."),
                   file.path(report_dir, "CODE_FIDELITY_AUDIT.md"))
  write_lines_over(c("# Benchmark Implementation Audit", "",
                     "PPR_abs_prior, PPR_signed_two_channel, RWR_abs_prior, and RWR_signed_two_channel use the same faithful M2 initial vector as NESTA, normalized as nonnegative personalization vectors.",
                     "PPR_abs/RWR_abs are direction-blind sensitivity comparators; signed two-channel PPR/RWR are directional propagation comparators."),
                   file.path(report_dir, "BENCHMARK_IMPLEMENTATION_AUDIT.md"))
  write_lines_over(c("# Metric Schema Audit", "",
                     sprintf("Primary schema pass: %s.", schema$required_columns_present[schema$file == "PRIMARY_FINAL_HEAT_METRICS.csv"]),
                     sprintf("Direction schema pass: %s.", schema$required_columns_present[schema$file == "DIRECTION_AWARE_METRICS.csv"]),
                     sprintf("Benchmark schema pass: %s.", schema$required_columns_present[schema$file == "BENCHMARK_METRICS.csv"])),
                   file.path(report_dir, "METRIC_SCHEMA_AUDIT.md"))
  write_lines_over(c("# STOP/GO Report", "", "STOP.",
                     paste0("Reason: `", reason, "`."),
                     "Pilot started: YES, diagnostic audit only.",
                     "Confirmatory started: NO.",
                     "This binding plan prohibits confirmatory execution."),
                   file.path(report_dir, "STOP_GO_REPORT.md"))
  write_lines_over(c("# Final Report", "",
                     paste0("Binding plan SHA256: `", binding_plan_sha256, "`."),
                     "Dense-only operator and ranking-orientation audit completed on frozen Arm F and Arm H diagnostic substrates.",
                     "Confirmatory execution started: NO.",
                     paste0("Final outcome classification: `", reason, "`."),
                     sprintf("Best NESTA-derived score `%s` on arm `%s`: top100/top150/top200 %.4f / %.4f / %.4f; median A2 rank %.2f.",
                             best$score_name, best$arm, best$top100_recall, best$top150_recall, best$top200_recall, best$median_A2_rank),
                     sprintf("Best same-initial PPR/RWR score `%s` on arm `%s`: top100/top200 %.4f / %.4f.",
                             best_bench$score_name, best_bench$arm, best_bench$top100_recall, best_bench$top200_recall),
                     sprintf("Seed filtering failure: %s. Minimal positive controls passed: %s. Promising NESTA intermediate/orientation result: %s.",
                             seed_failure, pc_success, promising)),
                   file.path(report_dir, "FINAL_REPORT.md"))
  if (file.exists(project_file("CLEANUP_MANIFEST.csv"))) {
    file.copy(project_file("CLEANUP_MANIFEST.csv"), file.path(report_dir, "CLEANUP_MANIFEST.csv"), overwrite = TRUE)
  }
  write_csv_over(data.frame(path = c(project_file("R/study_0707_operator_audit.R"),
                                     project_file("R/study_0707_branch_isolation.R"),
                                     project_file("R/study_0707_dense_rescue.R"),
                                     project_file("R/fidelity.R"), project_file("R/utils.R"),
                                     project_file("scripts/run_operator_audit.R")),
                            sha256 = c(sha(project_file("R/study_0707_operator_audit.R")),
                                       sha(project_file("R/study_0707_branch_isolation.R")),
                                       sha(project_file("R/study_0707_dense_rescue.R")),
                                       sha(project_file("R/fidelity.R")),
                                       sha(project_file("R/utils.R")),
                                       if (file.exists(project_file("scripts/run_operator_audit.R"))) sha(project_file("scripts/run_operator_audit.R")) else NA_character_),
                            role = c("operator_audit_runner", "frozen_arm_FH_substrates",
                                     "faithful_channel_helper", "faithful_nesta_and_benchmarks",
                                     "binding_plan_and_io_guards", "script_entrypoint"),
                            stringsAsFactors = FALSE),
                 file.path(report_dir, "IMPLEMENTATION_MANIFEST.csv"))
  files <- list.files(report_dir, recursive = TRUE, full.names = TRUE)
  files <- files[file.info(files)$isdir == FALSE & basename(files) != "CHECKSUMS.sha256"]
  writeLines(vapply(files, function(f) paste(sha(f), sub(paste0("^", report_dir, "/"), "", f)), character(1)),
             file.path(report_dir, "CHECKSUMS.sha256"))
  invisible(list(reason = reason, best = best, best_benchmark = best_bench,
                 promising = promising, pc_success = pc_success))
}

if (!identical(Sys.getenv("NESTA_OPERATOR_AUDIT_SOURCE_ONLY"), "1")) run_operator_audit()
