source(file.path("/home/js/NESTA/simulation", "R", "utils.R"))
source(project_file("R/fidelity.R"))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(igraph))

if (!exists("write_csv_over")) {
  write_csv_over <- function(x, path) {
    if (file.exists(path)) unlink(path)
    atomic_write_csv(x, path)
  }
}
if (!exists("write_lines_over")) {
  write_lines_over <- function(x, path) {
    if (file.exists(path)) unlink(path)
    atomic_write_lines(x, path)
  }
}
sha <- function(path) digest::digest(file = path, algo = "sha256")

dense_rescue_arms <- function() {
  list(
    A = list(label = "Arm_A_060726_1319_baseline",
             a1_relay = c(0.275, 0.290), relay_a2 = c(0.175, 0.195),
             a2_cluster = c(0.88, 0.93), weak_not_blank = FALSE),
    B = list(label = "Arm_B_terminal_capture_strengthened",
             a1_relay = c(0.275, 0.290), relay_a2 = c(0.205, 0.225),
             a2_cluster = c(0.88, 0.94), weak_not_blank = FALSE),
    C = list(label = "Arm_C_weak_not_blank_A2_signal",
             a1_relay = c(0.275, 0.290), relay_a2 = c(0.180, 0.200),
             a2_cluster = c(0.88, 0.93), weak_not_blank = TRUE),
    D = list(label = "Arm_D_terminal_capture_plus_weak_signal",
             a1_relay = c(0.275, 0.290), relay_a2 = c(0.205, 0.225),
             a2_cluster = c(0.88, 0.94), weak_not_blank = TRUE)
  )
}

add_sparse_edge <- function(adj, a, b, w) {
  adj[a, b] <- w
  adj[b, a] <- w
  adj
}

make_dense_rescue_rep <- function(arm_id, rep_id, seed) {
  arm <- dense_rescue_arms()[[arm_id]]
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
  branch <- c(setNames(rep("risk", length(risk_a2)), risk_a2),
              setNames(rep("protective", length(prot_a2)), prot_a2))
  relay_branch <- c(setNames(rep("risk", length(risk_relay)), risk_relay),
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

  add_block(c(risk_a1, risk_relay, risk_a2), 0.10, 0.11, 0.16)
  add_block(c(prot_a1, prot_relay, prot_a2), 0.10, 0.11, 0.16)
  add_block(risk_a2, 0.55, arm$a2_cluster[1] / 7, arm$a2_cluster[2] / 6)
  add_block(prot_a2, 0.55, arm$a2_cluster[1] / 7, arm$a2_cluster[2] / 6)
  add_block(relay, 0.025, 0.10, 0.13)

  for (r in risk_relay) add(r, sample(risk_a1, 3), runif(3, arm$a1_relay[1], arm$a1_relay[2]))
  for (r in prot_relay) add(r, sample(prot_a1, 3), runif(3, arm$a1_relay[1], arm$a1_relay[2]))
  for (g in risk_a2) add(g, sample(risk_relay, 4), runif(4, arm$relay_a2[1], arm$relay_a2[2]))
  for (g in prot_a2) add(g, sample(prot_relay, 4), runif(4, arm$relay_a2[1], arm$relay_a2[2]))
  # A small number of same-sign direct branch contacts, not enough to make A2 direct-neighbor dominated.
  for (g in sample(risk_a2, 5)) add(g, sample(risk_a1, 1), runif(1, 0.11, 0.14))
  for (g in sample(prot_a2, 5)) add(g, sample(prot_a1, 1), runif(1, 0.11, 0.14))

  m_edges <- 5500
  add(sample(bg, m_edges, TRUE), sample(bg, m_edges, TRUE), runif(m_edges, 0.001, 0.07))
  add(sample(c(a1, relay, a2), 100, TRUE), sample(bg, 100, TRUE), runif(100, 0.001, 0.025))
  add(sample(cdecoy, 220, TRUE), sample(bg, 220, TRUE), runif(220, 0.03, 0.09))

  adj <- sparseMatrix(i = ii, j = jj, x = xx, dims = c(n, n), dimnames = list(genes, genes))
  adj <- forceSymmetric(adj, uplo = "U")
  diag(adj) <- 1

  expr <- setNames(rlnorm(n, 0, 0.30), genes)
  expr[relay] <- pmin(expr[relay], stats::median(expr) * 0.95)
  expr[a2] <- stats::median(expr) * runif(length(a2), 0.92, 1.08)
  expr[cdecoy] <- stats::median(expr) * runif(length(cdecoy), 1.00, 1.15)

  z <- setNames(rnorm(n, 0, 1), genes)
  z[risk_a1] <- abs(rnorm(length(risk_a1), 3.5, 0.45))
  z[prot_a1] <- -abs(rnorm(length(prot_a1), 3.5, 0.45))
  z[relay] <- rnorm(length(relay), 0, 0.12)
  if (isTRUE(arm$weak_not_blank)) {
    risk_boost <- sample(risk_a2, 1)
    prot_boost <- sample(prot_a2, 1)
    z[risk_a2] <- rnorm(length(risk_a2), 0.16, 0.12)
    z[prot_a2] <- rnorm(length(prot_a2), -0.16, 0.12)
    z[risk_boost] <- runif(1, 1.05, 1.25)
    z[prot_boost] <- -runif(1, 1.05, 1.25)
  } else {
    z[risk_a2] <- rnorm(length(risk_a2), 0.04, 0.08)
    z[prot_a2] <- rnorm(length(prot_a2), -0.04, 0.08)
  }
  z[decoy[1:30]] <- -abs(rnorm(30, 1.15, 0.35))
  z[decoy[31:60]] <- abs(rnorm(30, 1.15, 0.35))
  z <- pmax(pmin(z, 4.5), -4.5)

  twas <- data.frame(SYMBOL = genes, TWAS.Z = as.numeric(z), TWAS.P = twas_p_from_z(z),
                     stringsAsFactors = FALSE)
  list(condition = "dense_1000", arm = arm_id, arm_label = arm$label, rep_id = rep_id,
       genes = genes, adj = adj, expr = expr, twas = twas, A1 = a1, A1_risk = risk_a1,
       A1_protective = prot_a1, relay = relay, relay_risk = risk_relay,
       relay_protective = prot_relay, relay_direction = relay_branch,
       A2 = a2, A2_risk = risk_a2, A2_protective = prot_a2, A2_direction = branch,
       D = decoy, C = cdecoy, background = bg)
}

faithful_m2_channels <- function(adj, mean_expression, twas, twas_cutoff = 1,
                                 edge_cutoff = 0.1, n.perm = 25, seed = 9707) {
  genes <- rownames(adj)
  twas_z <- strict_twas_vector(genes, twas, cutoff = twas_cutoff)
  expression_factor <- mean_expression / safe_sd(mean_expression)
  initial_weight <- expression_factor * (twas_z / safe_sd(twas_z))
  pos <- ifelse(initial_weight >= 0, initial_weight, 0)
  neg <- ifelse(initial_weight <= 0, abs(initial_weight), 0)
  names(pos) <- genes; names(neg) <- genes
  g <- make_faithful_graph(adj, edge_cutoff = edge_cutoff, retain_zero_edges = TRUE,
                           submitted_diagonal = TRUE)
  k <- regularised_kernel(g)
  pos_heat <- as.numeric(diffuse_checked(g, pos, n.perm = n.perm, seed = seed, K = k))
  neg_heat <- as.numeric(diffuse_checked(g, neg, n.perm = n.perm, seed = seed, K = k))
  pos_scaled <- pos_heat / safe_sd(pos_heat)
  neg_scaled <- neg_heat / safe_sd(neg_heat)
  signed <- (pos_scaled - neg_scaled) / safe_sd(pos_scaled - neg_scaled)
  data.frame(SYMBOL = genes, TWAS.Z = as.numeric(twas_z),
             expression_factor = as.numeric(expression_factor),
             initial_weight = as.numeric(initial_weight),
             positive_channel_heat_raw = pos_heat,
             negative_channel_heat_raw = neg_heat,
             positive_channel_heat = pos_scaled,
             negative_channel_heat = neg_scaled,
             signed_NESTA_heat = signed,
             unsigned_NESTA_heat = abs(signed),
             final_NESTA_heat = signed,
             final_heat = signed,
             delta_NESTA = signed - as.numeric(twas_z),
             stringsAsFactors = FALSE)
}

ranked_genes <- function(genes, score, decreasing = TRUE, exclude = character()) {
  keep <- !genes %in% exclude
  genes[keep][order(score[keep], decreasing = decreasing, na.last = NA)]
}

recall_at <- function(ranked, positives, k) {
  k <- min(k, length(ranked))
  if (!length(positives)) return(NA_real_)
  mean(positives %in% ranked[seq_len(k)])
}

rank_map <- function(genes, score, exclude = character(), decreasing = TRUE) {
  ranked <- ranked_genes(genes, score, decreasing = decreasing, exclude = exclude)
  setNames(seq_along(ranked), ranked)
}

non_directional_metrics <- function(rep, score, score_name) {
  genes <- rep$genes
  ranked <- ranked_genes(genes, score, TRUE, rep$A1)
  positives <- rep$A2
  n_ranked <- length(ranked)
  top_ks <- c(top50 = 50, top100 = 100, top150 = 150, top200 = 200,
              top1pct = ceiling(0.01 * n_ranked), top5pct = ceiling(0.05 * n_ranked),
              top10pct = ceiling(0.10 * n_ranked))
  recalls <- vapply(top_ks, function(k) recall_at(ranked, positives, k), numeric(1))
  au <- auprc_from_score(genes, score, positives, rep$A1)
  prev <- length(positives) / n_ranked
  data.frame(arm = rep$arm, arm_label = rep$arm_label, condition = rep$condition,
             replicate = rep$rep_id, score_name = score_name,
             top50_recall = recalls["top50"], top100_recall = recalls["top100"],
             top150_recall = recalls["top150"], top200_recall = recalls["top200"],
             top1pct_recall = recalls["top1pct"], top5pct_recall = recalls["top5pct"],
             top10pct_recall = recalls["top10pct"], raw_AUPRC = au,
             prevalence_normalized_AUPRC = if (is.finite(au) && prev < 1) (au - prev) / (1 - prev) else NA_real_,
             top100_fold_enrichment = recalls["top100"] / (min(100, n_ranked) / n_ranked),
             top150_fold_enrichment = recalls["top150"] / (min(150, n_ranked) / n_ranked),
             top200_fold_enrichment = recalls["top200"] / (min(200, n_ranked) / n_ranked),
             top10pct_fold_enrichment = recalls["top10pct"] / (top_ks["top10pct"] / n_ranked),
             stringsAsFactors = FALSE)
}

directional_metrics <- function(rep, signed_score, score_name) {
  genes <- rep$genes
  abs_score <- abs(signed_score)
  ranked <- ranked_genes(genes, abs_score, TRUE, rep$A1)
  dirs <- rep$A2_direction
  signedv <- setNames(signed_score, genes)
  concord <- function(top) {
    g <- intersect(top, names(dirs))
    if (!length(g)) return(logical())
    ifelse(dirs[g] == "risk", signedv[g] > 0, signedv[g] < 0)
  }
  top100 <- ranked[seq_len(min(100, length(ranked)))]
  top5 <- ranked[seq_len(min(ceiling(0.05 * length(ranked)), length(ranked)))]
  risk_au <- auprc_from_score(genes, signed_score, rep$A2_risk, rep$A1)
  prot_au <- auprc_from_score(genes, -signed_score, rep$A2_protective, rep$A1)
  data.frame(arm = rep$arm, arm_label = rep$arm_label, condition = rep$condition,
             replicate = rep$rep_id, score_name = score_name,
             direction_aware_AUPRC = weighted.mean(c(risk_au, prot_au),
                                                   c(length(rep$A2_risk), length(rep$A2_protective)),
                                                   na.rm = TRUE),
             sign_concordant_top100_recall = sum(concord(top100), na.rm = TRUE) / length(rep$A2),
             sign_concordant_top5pct_recall = sum(concord(top5), na.rm = TRUE) / length(rep$A2),
             direction_accuracy_top100 = if (length(intersect(top100, names(dirs)))) mean(concord(top100), na.rm = TRUE) else NA_real_,
             opposite_sign_decoy_top100_rate = mean(rep$D %in% top100),
             opposite_sign_decoy_top5pct_rate = mean(rep$D %in% top5),
             stringsAsFactors = FALSE)
}

channel_polarity_row <- function(rep, ch) {
  sv <- setNames(ch$signed_NESTA_heat, ch$SYMBOL)
  pv <- setNames(ch$positive_channel_heat, ch$SYMBOL)
  nv <- setNames(ch$negative_channel_heat, ch$SYMBOL)
  uv <- setNames(ch$unsigned_NESTA_heat, ch$SYMBOL)
  data.frame(arm = rep$arm, arm_label = rep$arm_label, replicate = rep$rep_id,
             mean_positive_channel_A2 = mean(pv[rep$A2]),
             mean_negative_channel_A2 = mean(nv[rep$A2]),
             mean_signed_heat_risk_A2 = mean(sv[rep$A2_risk]),
             mean_signed_heat_protective_A2 = mean(sv[rep$A2_protective]),
             fraction_risk_A2_positive_signed_heat = mean(sv[rep$A2_risk] > 0),
             fraction_protective_A2_negative_signed_heat = mean(sv[rep$A2_protective] < 0),
             fraction_A2_nonzero_final_heat = mean(abs(uv[rep$A2]) > 0),
             stringsAsFactors = FALSE)
}

terminal_capture_row <- function(rep, ch) {
  uv <- setNames(ch$unsigned_NESTA_heat, ch$SYMBOL)
  ranked <- ranked_genes(rep$genes, uv, TRUE, rep$A1)
  top <- function(k) ranked[seq_len(min(k, length(ranked)))]
  branch <- unique(c(rep$A1, rep$relay, rep$A2))
  heat_sum <- sum(abs(uv))
  relay_mass <- sum(abs(uv[rep$relay]))
  a2_mass <- sum(abs(uv[rep$A2]))
  data.frame(arm = rep$arm, arm_label = rep$arm_label, replicate = rep$rep_id,
             fraction_seed_heat_retained_in_A_branch = sum(abs(uv[branch])) / heat_sum,
             fraction_seed_heat_reaching_relay = relay_mass / heat_sum,
             fraction_seed_heat_reaching_A2 = a2_mass / heat_sum,
             fraction_seed_heat_leaking_to_background = sum(abs(uv[rep$background])) / heat_sum,
             relay_heat_mass = relay_mass, A2_heat_mass = a2_mass,
             relay_to_A2_heat_ratio = relay_mass / max(a2_mass, .Machine$double.eps),
             relay_top100_rate = mean(rep$relay %in% top(100)),
             relay_top150_rate = mean(rep$relay %in% top(150)),
             relay_top200_rate = mean(rep$relay %in% top(200)),
             relay_top10pct_rate = mean(rep$relay %in% top(ceiling(.10 * length(ranked)))),
             A2_top100_rate = mean(rep$A2 %in% top(100)),
             A2_top150_rate = mean(rep$A2 %in% top(150)),
             A2_top200_rate = mean(rep$A2 %in% top(200)),
             A2_top10pct_rate = mean(rep$A2 %in% top(ceiling(.10 * length(ranked)))),
             median_A2_rank = median(rank_map(rep$genes, uv, rep$A1)[rep$A2], na.rm = TRUE),
             median_relay_rank = median(rank_map(rep$genes, uv, rep$A1)[rep$relay], na.rm = TRUE),
             stringsAsFactors = FALSE)
}

a2_rank_row <- function(rep, score, score_name) {
  r <- rank_map(rep$genes, score, rep$A1)
  data.frame(arm = rep$arm, arm_label = rep$arm_label, replicate = rep$rep_id,
             score_name = score_name,
             median_A2_rank = median(r[rep$A2], na.rm = TRUE),
             q25_A2_rank = as.numeric(stats::quantile(r[rep$A2], 0.25, na.rm = TRUE)),
             q75_A2_rank = as.numeric(stats::quantile(r[rep$A2], 0.75, na.rm = TRUE)),
             first_A2_rank = min(r[rep$A2], na.rm = TRUE),
             A2_top100_count = sum(r[rep$A2] <= 100, na.rm = TRUE),
             A2_top150_count = sum(r[rep$A2] <= 150, na.rm = TRUE),
             A2_top200_count = sum(r[rep$A2] <= 200, na.rm = TRUE),
             stringsAsFactors = FALSE)
}

score_degree_row <- function(rep, score, score_name) {
  deg <- Matrix::rowSums(rep$adj > 0.1) - 1
  genes <- rep$genes
  ranked <- ranked_genes(genes, score, TRUE, rep$A1)
  top100 <- ranked[seq_len(min(100, length(ranked)))]
  data.frame(arm = rep$arm, arm_label = rep$arm_label, replicate = rep$rep_id,
             score_name = score_name,
             score_degree_spearman = suppressWarnings(stats::cor(score, as.numeric(deg), method = "spearman")),
             top100_degree_percentile_median = stats::median(percentile_rank_desc(as.numeric(deg))[match(top100, genes)], na.rm = TRUE) * 100,
             high_degree_decoy_top100_rate = mean(rep$C %in% top100),
             stringsAsFactors = FALSE)
}

permute_sparse_weights <- function(adj, seed) {
  set.seed(seed)
  s <- summary(adj)
  keep <- s$i < s$j
  s <- s[keep, , drop = FALSE]
  vals <- sample(s$x)
  out <- sparseMatrix(i = c(s$i, s$j), j = c(s$j, s$i), x = c(vals, vals),
                      dims = dim(adj), dimnames = dimnames(adj))
  diag(out) <- 1
  out
}

module_disrupted_adj <- function(rep) {
  adj <- rep$adj
  nodes <- unique(c(rep$A1, rep$relay, rep$A2))
  adj[nodes, nodes] <- adj[nodes, nodes] * 0.03
  diag(adj) <- 1
  drop0(adj)
}

legacy_composite_score <- function(unsigned_heat, delta) {
  0.5 * percentile_rank_desc(unsigned_heat) + 0.5 * percentile_rank_desc(abs(delta))
}

benchmarks_for_rep <- function(rep, ch, no_diff) {
  genes <- rep$genes
  tw <- rep$twas
  p <- make_personalization(tw, genes)
  ps <- make_signed_personalization(tw, genes)
  rwr_abs <- run_rwr(rep$adj, p, restart = 0.5)
  ppr_abs <- run_ppr(rep$adj, p, damping = 0.85)
  rwr_signed <- run_rwr(rep$adj, ps$pos, restart = 0.5) - run_rwr(rep$adj, ps$neg, restart = 0.5)
  ppr_signed <- run_ppr(rep$adj, ps$pos, damping = 0.85) - run_ppr(rep$adj, ps$neg, damping = 0.85)
  raw_z <- setNames(tw$TWAS.Z, tw$SYMBOL)[genes]
  init <- setNames(no_diff$final_NESTA_heat, no_diff$SYMBOL)[genes]
  delta <- setNames(ch$signed_NESTA_heat - raw_z, genes)
  wp_scores <- lapply(seq_len(3), function(k) {
    a <- permute_sparse_weights(rep$adj, 9707000 + rep$rep_id * 10 + k)
    faithful_m2_channels(a, rep$expr, rep$twas, n.perm = 10, seed = 9707000 + k)
  })
  wp_signed <- Reduce("+", lapply(wp_scores, function(x) setNames(x$signed_NESTA_heat, x$SYMBOL))) / length(wp_scores)
  i2 <- faithful_m2_channels(module_disrupted_adj(rep), rep$expr, rep$twas, n.perm = 10,
                             seed = 9708000 + rep$rep_id)
  i3 <- faithful_m2_channels(permute_sparse_weights(rep$adj, 9709000 + rep$rep_id),
                             rep$expr, rep$twas, n.perm = 10, seed = 9709000 + rep$rep_id)
  list(
    NESTA_unsigned_final_heat = list(score = setNames(ch$unsigned_NESTA_heat, genes), signed = setNames(ch$signed_NESTA_heat, genes)),
    NESTA_signed_reconstructed_heat = list(score = abs(setNames(ch$signed_NESTA_heat, genes)), signed = setNames(ch$signed_NESTA_heat, genes)),
    raw_TWAS_abs = list(score = abs(raw_z), signed = raw_z),
    raw_signed_TWAS = list(score = abs(raw_z), signed = raw_z),
    M2_no_diffusion = list(score = abs(init), signed = init),
    legacy_delta_only = list(score = abs(delta), signed = delta),
    old_heat_delta_composite = list(score = legacy_composite_score(setNames(ch$unsigned_NESTA_heat, genes), delta), signed = setNames(ch$signed_NESTA_heat, genes)),
    RWR_abs_prior = list(score = rwr_abs, signed = rwr_abs),
    PPR_abs_prior = list(score = ppr_abs, signed = ppr_abs),
    RWR_signed_two_channel = list(score = abs(rwr_signed), signed = rwr_signed),
    PPR_signed_two_channel = list(score = abs(ppr_signed), signed = ppr_signed),
    median_weight_permutation = list(score = abs(wp_signed), signed = wp_signed),
    I2_module_disrupted = list(score = abs(setNames(i2$signed_NESTA_heat, i2$SYMBOL)), signed = setNames(i2$signed_NESTA_heat, i2$SYMBOL)),
    I3_expression_matched_randomized = list(score = abs(setNames(i3$signed_NESTA_heat, i3$SYMBOL)), signed = setNames(i3$signed_NESTA_heat, i3$SYMBOL))
  )
}

initial_signal_ok <- function(rep, no_diff) {
  raw_rank <- rank(-abs(rep$twas$TWAS.Z), ties.method = "average")
  names(raw_rank) <- rep$twas$SYMBOL
  init <- setNames(no_diff$final_NESTA_heat, no_diff$SYMBOL)[rep$genes]
  init_rank <- rank(-abs(init), ties.method = "average")
  names(init_rank) <- rep$genes
  list(raw_top100 = mean(raw_rank[rep$A2] <= 100),
       m2_top100 = mean(init_rank[rep$A2] <= 100),
       raw_top10pct = mean(raw_rank[rep$A2] <= 100),
       m2_top10pct = mean(init_rank[rep$A2] <= 100),
       relay_raw_top100 = mean(raw_rank[rep$relay] <= 100),
       relay_m2_top100 = mean(init_rank[rep$relay] <= 100))
}

ci_diff <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(mean = NA_real_, ci_low = NA_real_, ci_high = NA_real_))
  if (length(x) == 1) return(c(mean = x, ci_low = x, ci_high = x))
  boot <- replicate(1000, mean(sample(x, replace = TRUE), na.rm = TRUE))
  c(mean = mean(x), ci_low = as.numeric(stats::quantile(boot, 0.025, na.rm = TRUE)),
    ci_high = as.numeric(stats::quantile(boot, 0.975, na.rm = TRUE)))
}

summarise_arm_metrics <- function(primary, direction, terminal, ranks) {
  nesta <- primary[primary$score_name == "NESTA_unsigned_final_heat", , drop = FALSE]
  signed <- direction[direction$score_name == "NESTA_signed_reconstructed_heat", , drop = FALSE]
  do.call(rbind, lapply(split(nesta, nesta$arm), function(x) {
    arm <- unique(x$arm)
    sx <- signed[signed$arm == arm, , drop = FALSE]
    tx <- terminal[terminal$arm == arm, , drop = FALSE]
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
               relay_to_A2_heat_ratio = mean(tx$relay_to_A2_heat_ratio, na.rm = TRUE),
               median_A2_rank = mean(tx$median_A2_rank, na.rm = TRUE),
               median_relay_rank = mean(tx$median_relay_rank, na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))
}

make_contrasts <- function(primary, direction, ranks) {
  out <- list()
  idx <- 1
  for (arm in unique(primary$arm)) {
    for (metric in c("top100_recall", "top150_recall", "top200_recall", "raw_AUPRC")) {
      base <- primary[primary$arm == arm & primary$score_name == "NESTA_unsigned_final_heat",
                      c("replicate", metric), drop = FALSE]
      for (cmp in c("raw_TWAS_abs", "M2_no_diffusion", "PPR_abs_prior", "RWR_abs_prior")) {
        other <- primary[primary$arm == arm & primary$score_name == cmp,
                         c("replicate", metric), drop = FALSE]
        z <- merge(base, other, by = "replicate", suffixes = c("_nesta", "_cmp"))
        d <- z[[paste0(metric, "_nesta")]] - z[[paste0(metric, "_cmp")]]
        ci <- ci_diff(d)
        out[[idx]] <- data.frame(arm = arm, metric = metric, comparator = cmp,
                                 nesta_mean = mean(z[[paste0(metric, "_nesta")]], na.rm = TRUE),
                                 comparator_mean = mean(z[[paste0(metric, "_cmp")]], na.rm = TRUE),
                                 mean_difference = ci["mean"], ci_low = ci["ci_low"],
                                 ci_high = ci["ci_high"], stringsAsFactors = FALSE)
        idx <- idx + 1
      }
    }
    base_r <- ranks[ranks$arm == arm & ranks$score_name == "NESTA_unsigned_final_heat",
                    c("replicate", "median_A2_rank"), drop = FALSE]
    for (cmp in c("raw_signed_TWAS", "M2_no_diffusion")) {
      other <- ranks[ranks$arm == arm & ranks$score_name == cmp,
                     c("replicate", "median_A2_rank"), drop = FALSE]
      z <- merge(base_r, other, by = "replicate", suffixes = c("_nesta", "_cmp"))
      d <- z$median_A2_rank_cmp - z$median_A2_rank_nesta
      ci <- ci_diff(d)
      out[[idx]] <- data.frame(arm = arm, metric = "median_A2_rank_improvement",
                               comparator = cmp,
                               nesta_mean = mean(z$median_A2_rank_nesta, na.rm = TRUE),
                               comparator_mean = mean(z$median_A2_rank_cmp, na.rm = TRUE),
                               mean_difference = ci["mean"], ci_low = ci["ci_low"],
                               ci_high = ci["ci_high"], stringsAsFactors = FALSE)
      idx <- idx + 1
    }
  }
  do.call(rbind, out)
}

metric_schema_audit <- function(primary, direction) {
  primary_req <- c("top50_recall", "top100_recall", "top150_recall", "top200_recall",
                   "top1pct_recall", "top5pct_recall", "top10pct_recall", "raw_AUPRC",
                   "prevalence_normalized_AUPRC", "top100_fold_enrichment")
  direction_req <- c("direction_aware_AUPRC", "sign_concordant_top100_recall",
                     "sign_concordant_top5pct_recall", "direction_accuracy_top100",
                     "opposite_sign_decoy_top100_rate", "opposite_sign_decoy_top5pct_rate")
  data.frame(file = c("PRIMARY_FINAL_HEAT_METRICS.csv", "DIRECTION_AWARE_METRICS.csv"),
             required_columns_present = c(all(primary_req %in% names(primary)),
                                          all(direction_req %in% names(direction))),
             missing_columns = c(paste(setdiff(primary_req, names(primary)), collapse = ";"),
                                 paste(setdiff(direction_req, names(direction)), collapse = ";")),
             misleading_size_smoke_failure_report_written = FALSE,
             stringsAsFactors = FALSE)
}

run_dense_rescue <- function() {
  verify_project_path()
  verify_binding_plan()
  report_dir <- read_report_dir()
  safe_dir_create(report_dir)
  copy_binding_plan_to_report(report_dir)
  set_run_status("IN_PROGRESS", "YES", "NO", "0707 dense-only polarity terminal-capture rescue")

  primary_rows <- list(); direction_rows <- list(); benchmark_rows <- list()
  polarity_rows <- list(); terminal_rows <- list(); rank_rows <- list()
  degree_rows <- list(); signal_rows <- list()
  i_primary <- i_direction <- i_benchmark <- i_polarity <- i_terminal <- i_rank <- i_degree <- i_signal <- 1

  arms <- dense_rescue_arms()
  for (arm in names(arms)) {
    for (rep_id in seq_len(20)) {
      rep <- make_dense_rescue_rep(arm, rep_id, 970700 + match(arm, names(arms)) * 1000 + rep_id)
      ch <- faithful_m2_channels(rep$adj, rep$expr, rep$twas, n.perm = 25,
                                 seed = 970700 + rep_id)
      no_diff <- no_diffusion_m2_scores(rep$adj, rep$expr, rep$twas)
      sig <- initial_signal_ok(rep, no_diff)
      signal_rows[[i_signal]] <- data.frame(arm = rep$arm, arm_label = rep$arm_label,
                                            replicate = rep$rep_id, raw_TWAS_top100_A2_fraction = sig$raw_top100,
                                            M2_initial_top100_A2_fraction = sig$m2_top100,
                                            raw_TWAS_top10pct_A2_fraction = sig$raw_top10pct,
                                            M2_initial_top10pct_A2_fraction = sig$m2_top10pct,
                                            relay_raw_TWAS_top100_fraction = sig$relay_raw_top100,
                                            relay_M2_initial_top100_fraction = sig$relay_m2_top100,
                                            smoke_pass = sig$raw_top100 <= 0.10 && sig$m2_top100 <= 0.10,
                                            stringsAsFactors = FALSE)
      i_signal <- i_signal + 1

      scores <- benchmarks_for_rep(rep, ch, no_diff)
      for (nm in names(scores)) {
        primary_rows[[i_primary]] <- non_directional_metrics(rep, scores[[nm]]$score, nm)
        i_primary <- i_primary + 1
        direction_rows[[i_direction]] <- directional_metrics(rep, scores[[nm]]$signed, nm)
        i_direction <- i_direction + 1
        rank_rows[[i_rank]] <- a2_rank_row(rep, scores[[nm]]$score, nm)
        i_rank <- i_rank + 1
        degree_rows[[i_degree]] <- score_degree_row(rep, scores[[nm]]$score, nm)
        i_degree <- i_degree + 1
      }
      polarity_rows[[i_polarity]] <- channel_polarity_row(rep, ch); i_polarity <- i_polarity + 1
      terminal_rows[[i_terminal]] <- terminal_capture_row(rep, ch); i_terminal <- i_terminal + 1
      if (rep_id %% 5 == 0) gc(FALSE)
    }
  }

  primary <- do.call(rbind, primary_rows)
  direction <- do.call(rbind, direction_rows)
  ranks <- do.call(rbind, rank_rows)
  degree <- do.call(rbind, degree_rows)
  polarity <- do.call(rbind, polarity_rows)
  terminal <- do.call(rbind, terminal_rows)
  signal <- do.call(rbind, signal_rows)
  summary <- summarise_arm_metrics(primary, direction, terminal, ranks)
  contrasts <- make_contrasts(primary, direction, ranks)
  schema <- metric_schema_audit(primary, direction)

  write_csv_over(summary, file.path(report_dir, "DENSE_RESCUE_ARM_METRICS.csv"))
  write_csv_over(contrasts, file.path(report_dir, "DENSE_RESCUE_ARM_CONTRASTS.csv"))
  write_csv_over(primary, file.path(report_dir, "PRIMARY_FINAL_HEAT_METRICS.csv"))
  write_csv_over(direction, file.path(report_dir, "DIRECTION_AWARE_METRICS.csv"))
  write_csv_over(primary, file.path(report_dir, "BENCHMARK_METRICS.csv"))
  write_csv_over(contrasts, file.path(report_dir, "BENCHMARK_CONTRASTS.csv"))
  write_csv_over(polarity, file.path(report_dir, "NESTA_CHANNEL_POLARITY_AUDIT.csv"))
  write_csv_over(terminal, file.path(report_dir, "TERMINAL_CAPTURE_AUDIT.csv"))
  write_csv_over(ranks, file.path(report_dir, "A2_RANK_DISTRIBUTION_AUDIT.csv"))
  write_csv_over(terminal[, c("arm", "arm_label", "replicate", "relay_heat_mass", "A2_heat_mass",
                              "relay_to_A2_heat_ratio", "relay_top100_rate", "A2_top100_rate",
                              "median_A2_rank", "median_relay_rank")],
                 file.path(report_dir, "RELAY_VS_A2_HEAT_AUDIT.csv"))
  write_csv_over(schema, file.path(report_dir, "METRIC_SCHEMA_AUDIT.csv"))
  write_csv_over(degree, file.path(report_dir, "DEGREE_BIAS_METRICS.csv"))
  write_csv_over(signal, file.path(report_dir, "DENSE_RESCUE_INITIAL_SIGNAL_AUDIT.csv"))
  write_csv_over(data.frame(guardrail = c("dense_only", "no_confirmatory", "raw_m2_top100_A2_leakage"),
                            passed = c(TRUE, TRUE, all(signal$raw_TWAS_top100_A2_fraction <= 0.10 &
                                                         signal$M2_initial_top100_A2_fraction <= 0.10)),
                            stringsAsFactors = FALSE),
                 file.path(report_dir, "NULL_BIAS_GUARDRAILS.csv"))

  all_zero <- all(summary$mean_top100_recall == 0 & summary$mean_top150_recall == 0 &
                    summary$mean_top200_recall == 0)
  promising <- summary$mean_top100_recall > 0 | summary$mean_top150_recall > 0 |
    summary$mean_top200_recall > 0
  reason <- if (all_zero) "NESTA_terminal_capture_failure_under_current_generator" else
    "dense_rescue_diagnostic_completed_no_confirmatory"
  status <- "STOPPED"
  set_run_status(status, "YES", "NO", reason)

  best <- summary[order(summary$mean_top100_recall, summary$mean_top150_recall,
                        summary$mean_top200_recall, decreasing = TRUE), ][1, ]
  write_lines_over(c(
    "# Dense Rescue Arm Summary", "",
    paste0("Best arm by unsigned top100 recall: `", best$arm, "` (", best$arm_label, ")."),
    sprintf("Top100/top150/top200 recall: %.4f / %.4f / %.4f.",
            best$mean_top100_recall, best$mean_top150_recall, best$mean_top200_recall),
    sprintf("Direction-aware AUPRC: %.4f; sign-concordant top100 recall: %.4f.",
            best$direction_aware_AUPRC, best$sign_concordant_top100_recall),
    sprintf("Relay top100 rate: %.4f; relay-to-A2 heat ratio: %.4f.",
            best$relay_top100_rate, best$relay_to_A2_heat_ratio)
  ), file.path(report_dir, "DENSE_RESCUE_ARM_SUMMARY.md"))

  write_lines_over(c(
    "# NESTA Channel Polarity Audit", "",
    "Positive and negative submitted diffusion channels were reconstructed as positive-channel heat minus negative-channel heat.",
    sprintf("Mean fraction of risk A2 with positive signed heat: %.4f.",
            mean(polarity$fraction_risk_A2_positive_signed_heat, na.rm = TRUE)),
    sprintf("Mean fraction of protective A2 with negative signed heat: %.4f.",
            mean(polarity$fraction_protective_A2_negative_signed_heat, na.rm = TRUE)),
    sprintf("Mean fraction of A2 with nonzero final heat: %.4f.",
            mean(polarity$fraction_A2_nonzero_final_heat, na.rm = TRUE))
  ), file.path(report_dir, "NESTA_CHANNEL_POLARITY_AUDIT.md"))

  write_lines_over(c(
    "# Terminal Capture Audit", "",
    sprintf("Mean seed heat retained in A branch: %.4f.",
            mean(terminal$fraction_seed_heat_retained_in_A_branch, na.rm = TRUE)),
    sprintf("Mean heat reaching relay/A2: %.4f / %.4f.",
            mean(terminal$fraction_seed_heat_reaching_relay, na.rm = TRUE),
            mean(terminal$fraction_seed_heat_reaching_A2, na.rm = TRUE)),
    sprintf("Mean background leakage: %.4f.",
            mean(terminal$fraction_seed_heat_leaking_to_background, na.rm = TRUE)),
    sprintf("Mean relay-to-A2 heat ratio: %.4f.",
            mean(terminal$relay_to_A2_heat_ratio, na.rm = TRUE))
  ), file.path(report_dir, "TERMINAL_CAPTURE_AUDIT.md"))

  write_lines_over(c(
    "# A2 Rank Distribution Audit", "",
    "A2 rank distributions are reported for NESTA, raw TWAS, M2 no diffusion, and RWR/PPR comparators in the paired CSV.",
    sprintf("Best-arm mean median A2 rank: %.2f.", best$median_A2_rank)
  ), file.path(report_dir, "A2_RANK_DISTRIBUTION_AUDIT.md"))

  write_lines_over(c(
    "# Metric Schema Audit", "",
    sprintf("PRIMARY_FINAL_HEAT_METRICS.csv schema pass: %s.",
            schema$required_columns_present[schema$file == "PRIMARY_FINAL_HEAT_METRICS.csv"]),
    sprintf("DIRECTION_AWARE_METRICS.csv schema pass: %s.",
            schema$required_columns_present[schema$file == "DIRECTION_AWARE_METRICS.csv"]),
    "No SIZE_CONDITION_SMOKE_FAILURE_REPORT.md was written because this dense pilot started."
  ), file.path(report_dir, "METRIC_SCHEMA_AUDIT.md"))

  write_lines_over(c(
    "# Code Fidelity Audit", "",
    paste0("Binding plan SHA256: `", binding_plan_sha256, "`."),
    "Faithful M2 arithmetic uses SCT-style row-mean expression divided by sample SD and TWAS.Z divided by sample SD.",
    "TWAS.P is computed as `2 * pnorm(-abs(TWAS.Z))`.",
    "Submitted-style positive and absolute-negative diffusion channels, retained zero-weight edges, self-loop behavior, and `n.perm` are used."
  ), file.path(report_dir, "CODE_FIDELITY_AUDIT.md"))

  write_lines_over(c(
    "# Benchmark Implementation Audit", "",
    "Comparators evaluated: raw_TWAS_abs, raw_signed_TWAS, M2_no_diffusion, legacy_delta_only, old_heat_delta_composite, RWR_abs_prior, PPR_abs_prior, RWR_signed_two_channel, PPR_signed_two_channel, median_weight_permutation, I2_module_disrupted, and I3_expression_matched_randomized.",
    "RWR_abs_prior and PPR_abs_prior are direction-blind sensitivity comparators.",
    "Median weight permutation is summarized from three weight-permuted faithful M2 reconstructions per replicate. I2 and I3 are diagnostic network-control reconstructions. This round has no confirmatory execution by design."
  ), file.path(report_dir, "BENCHMARK_IMPLEMENTATION_AUDIT.md"))

  write_lines_over(c(
    "# STOP/GO Report", "",
    "STOP.",
    paste0("Reason: `", reason, "`."),
    "Pilot started: YES.",
    "Confirmatory started: NO.",
    "This binding plan is diagnostic-only and prohibits confirmatory execution."
  ), file.path(report_dir, "STOP_GO_REPORT.md"))

  diag_line <- if (all_zero) {
    "All arms had zero NESTA unsigned top100/top150/top200 recovery, consistent with terminal capture failure under the current generator."
  } else {
    "At least one arm produced nonzero NESTA unsigned A2 recovery; compare terminal-capture and polarity audits to identify whether capture, polarity, or comparator behavior is limiting."
  }
  write_lines_over(c(
    "# Final Report", "",
    paste0("Binding plan SHA256: `", binding_plan_sha256, "`."),
    "Dense-only diagnostic rescue completed with 20 replicates per arm across Arms A-D.",
    "Pilot execution started: YES.",
    "Confirmatory execution started: NO.",
    paste0("STOP classification: `", reason, "`."),
    diag_line,
    sprintf("Best arm: `%s`; top100/top150/top200 recall %.4f / %.4f / %.4f.",
            best$arm, best$mean_top100_recall, best$mean_top150_recall, best$mean_top200_recall),
    sprintf("Best-arm direction-aware AUPRC %.4f and sign-concordant top100 recall %.4f.",
            best$direction_aware_AUPRC, best$sign_concordant_top100_recall)
  ), file.path(report_dir, "FINAL_REPORT.md"))

  write_lines_over(c(
    "# Manuscript Ready Summary", "",
    "This diagnostic-only dense rescue evaluated whether reconstructed signed NESTA channels and terminal A2 capture explain the prior direction-aware failure.",
    sprintf("Best arm `%s` achieved unsigned top100/top150/top200 A2 recall %.4f / %.4f / %.4f.",
            best$arm, best$mean_top100_recall, best$mean_top150_recall, best$mean_top200_recall),
    "Delta_NESTA was retained as auxiliary interpretation only and was not used as the primary endpoint."
  ), file.path(report_dir, "MANUSCRIPT_READY_SUMMARY.md"))

  if (file.exists(project_file("results/reports/summary_tables/unit_test_results.csv"))) {
    file.copy(project_file("results/reports/summary_tables/unit_test_results.csv"),
              file.path(report_dir, "unit_test_results.csv"), overwrite = TRUE)
  }
  if (file.exists(project_file("CLEANUP_MANIFEST.csv"))) {
    file.copy(project_file("CLEANUP_MANIFEST.csv"), file.path(report_dir, "CLEANUP_MANIFEST.csv"),
              overwrite = TRUE)
  }
  write_csv_over(data.frame(
    path = c(project_file("R/study_0707_dense_rescue.R"), project_file("R/fidelity.R"),
             project_file("R/utils.R"), project_file("scripts/run_dense_rescue.R")),
    sha256 = c(sha(project_file("R/study_0707_dense_rescue.R")),
               sha(project_file("R/fidelity.R")), sha(project_file("R/utils.R")),
               if (file.exists(project_file("scripts/run_dense_rescue.R"))) sha(project_file("scripts/run_dense_rescue.R")) else NA_character_),
    role = c("dense_rescue_runner", "faithful_nesta_and_benchmarks",
             "binding_plan_and_io_guards", "script_entrypoint"),
    stringsAsFactors = FALSE
  ), file.path(report_dir, "IMPLEMENTATION_MANIFEST.csv"))

  files <- list.files(report_dir, recursive = TRUE, full.names = TRUE)
  files <- files[file.info(files)$isdir == FALSE & basename(files) != "CHECKSUMS.sha256"]
  writeLines(vapply(files, function(f) paste(sha(f), sub(paste0("^", report_dir, "/"), "", f)), character(1)),
             file.path(report_dir, "CHECKSUMS.sha256"))
  invisible(list(reason = reason, summary = summary, promising = promising))
}

if (!identical(Sys.getenv("NESTA_DENSE_RESCUE_SOURCE_ONLY"), "1")) run_dense_rescue()
