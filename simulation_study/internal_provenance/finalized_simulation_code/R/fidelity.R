source(file.path("/home/js/NESTA/simulation", "R", "utils.R"))

twas_p_from_z <- function(z) {
  2 * stats::pnorm(-abs(z))
}

resolve_duplicate_twas <- function(twas) {
  req <- c("SYMBOL", "TWAS.Z", "TWAS.P")
  if (!all(req %in% names(twas))) stop("TWAS table must contain SYMBOL, TWAS.Z, TWAS.P")
  twas <- twas[order(twas$SYMBOL, twas$TWAS.P), , drop = FALSE]
  twas[!duplicated(twas$SYMBOL), , drop = FALSE]
}

strict_twas_vector <- function(genes, twas, cutoff = 1) {
  twas <- resolve_duplicate_twas(twas)
  twas <- twas[twas$TWAS.P < cutoff, , drop = FALSE]
  out <- setNames(rep(0, length(genes)), genes)
  idx <- match(twas$SYMBOL, genes)
  keep <- !is.na(idx)
  out[idx[keep]] <- twas$TWAS.Z[keep]
  out
}

reject_nperm <- function(args) {
  if ("nperm" %in% names(args)) stop("Unsupported diffuStats argument `nperm`; use `n.perm`")
  TRUE
}

diffuse_checked <- function(graph, scores, method = "ber_p", n.perm = 300, seed = 9703, ...) {
  dots <- list(...)
  reject_nperm(dots)
  diffuStats::diffuse(graph = graph, scores = scores, method = method,
                      n.perm = n.perm, seed = seed, ...)
}

regularised_kernel <- function(graph) {
  getFromNamespace("regularisedLaplacianKernel", "diffuStats")(graph = graph)
}

make_faithful_graph <- function(adj, edge_cutoff = 0.1, retain_zero_edges = TRUE,
                                submitted_diagonal = TRUE) {
  stopifnot(nrow(adj) == ncol(adj), !is.null(rownames(adj)), identical(rownames(adj), colnames(adj)))
  m <- adj
  if (submitted_diagonal) diag(m) <- 1
  if (retain_zero_edges) {
    g <- igraph::graph_from_adjacency_matrix(m, mode = "undirected", weighted = TRUE, diag = submitted_diagonal)
    igraph::E(g)$weight <- ifelse(igraph::E(g)$weight > edge_cutoff, igraph::E(g)$weight, 0)
  } else {
    m[m <= edge_cutoff] <- 0
    g <- igraph::graph_from_adjacency_matrix(m, mode = "undirected", weighted = TRUE, diag = submitted_diagonal)
  }
  g
}

faithful_m2_scores <- function(adj, mean_expression, twas, method = "ber_p",
                               twas_cutoff = 1, edge_cutoff = 0.1,
                               n.perm = 300, seed = 9703) {
  genes <- rownames(adj)
  if (!identical(genes, names(mean_expression))) stop("Mean-expression vector is not aligned to adjacency")
  twas_z <- strict_twas_vector(genes, twas, cutoff = twas_cutoff)
  expression_factor <- mean_expression / safe_sd(mean_expression)
  initial_weight <- expression_factor * (twas_z / safe_sd(twas_z))
  diffuse_weighted_scores(adj, twas_z, initial_weight, method = method,
                          edge_cutoff = edge_cutoff, n.perm = n.perm,
                          seed = seed, method_label = "D0_faithful_M2",
                          expression_factor = expression_factor)
}

diffuse_weighted_scores <- function(adj, twas_z, initial_weight, method = "ber_p",
                                    edge_cutoff = 0.1, n.perm = 300, seed = 9703,
                                    method_label = NA_character_,
                                    expression_factor = rep(NA_real_, length(twas_z))) {
  genes <- rownames(adj)
  if (!identical(genes, names(twas_z))) stop("TWAS vector is not aligned to adjacency")
  if (!identical(genes, names(initial_weight))) stop("Initial-weight vector is not aligned to adjacency")
  pos <- ifelse(initial_weight >= 0, initial_weight, 0)
  neg <- ifelse(initial_weight <= 0, abs(initial_weight), 0)
  names(pos) <- genes
  names(neg) <- genes
  g <- make_faithful_graph(adj, edge_cutoff = edge_cutoff, retain_zero_edges = TRUE, submitted_diagonal = TRUE)
  K <- regularised_kernel(g)
  pos_heat <- diffuse_checked(g, pos, method = method, n.perm = n.perm, seed = seed, K = K)
  neg_heat <- diffuse_checked(g, neg, method = method, n.perm = n.perm, seed = seed, K = K)
  pos_heat <- as.numeric(pos_heat)
  neg_heat <- as.numeric(neg_heat)
  final_heat <- (pos_heat / safe_sd(pos_heat)) - (neg_heat / safe_sd(neg_heat))
  final_heat <- final_heat / safe_sd(final_heat)
  data.frame(
    SYMBOL = genes,
    TWAS.Z = as.numeric(twas_z),
    expression_factor = as.numeric(expression_factor),
    initial_weight = as.numeric(initial_weight),
    final_NESTA_heat = as.numeric(final_heat),
    final_heat = as.numeric(final_heat),
    delta_NESTA = as.numeric(final_heat - twas_z),
    delta_twas = as.numeric(final_heat - twas_z),
    delta_initial = as.numeric(final_heat - initial_weight),
    score_abs_delta_NESTA = abs(as.numeric(final_heat - twas_z)),
    score_abs_delta_twas = abs(as.numeric(final_heat - twas_z)),
    diagnostic_method = method_label,
    stringsAsFactors = FALSE
  ) |> add_nesta_composite()
}

m1_scores <- function(adj, twas, method = "ber_p", twas_cutoff = 1, edge_cutoff = 0.1,
                      n.perm = 300, seed = 9703) {
  genes <- rownames(adj)
  twas_z <- strict_twas_vector(genes, twas, cutoff = twas_cutoff)
  diffuse_weighted_scores(adj, twas_z, twas_z, method = method,
                          edge_cutoff = edge_cutoff, n.perm = n.perm,
                          seed = seed, method_label = "D1_M1",
                          expression_factor = setNames(rep(1, length(genes)), genes))
}

no_diffusion_m2_scores <- function(adj, mean_expression, twas, twas_cutoff = 1) {
  genes <- rownames(adj)
  twas_z <- strict_twas_vector(genes, twas, cutoff = twas_cutoff)
  initial_weight <- (mean_expression / safe_sd(mean_expression)) * (twas_z / safe_sd(twas_z))
  final_heat <- initial_weight / safe_sd(initial_weight)
  data.frame(
    SYMBOL = genes,
    TWAS.Z = as.numeric(twas_z),
    expression_factor = as.numeric(mean_expression / safe_sd(mean_expression)),
    initial_weight = as.numeric(initial_weight),
    final_NESTA_heat = as.numeric(final_heat),
    final_heat = as.numeric(final_heat),
    delta_NESTA = as.numeric(final_heat - twas_z),
    delta_twas = as.numeric(final_heat - twas_z),
    delta_initial = as.numeric(final_heat - initial_weight),
    score_abs_delta_NESTA = abs(as.numeric(final_heat - twas_z)),
    score_abs_delta_twas = abs(as.numeric(final_heat - twas_z)),
    diagnostic_method = "M2_no_diffusion",
    stringsAsFactors = FALSE
  ) |> add_nesta_composite()
}

no_diffusion_m1_scores <- function(adj, twas, twas_cutoff = 1) {
  genes <- rownames(adj)
  twas_z <- strict_twas_vector(genes, twas, cutoff = twas_cutoff)
  final_heat <- twas_z / safe_sd(twas_z)
  data.frame(
    SYMBOL = genes,
    TWAS.Z = as.numeric(twas_z),
    expression_factor = rep(1, length(genes)),
    initial_weight = as.numeric(twas_z),
    final_NESTA_heat = as.numeric(final_heat),
    final_heat = as.numeric(final_heat),
    delta_NESTA = as.numeric(final_heat - twas_z),
    delta_twas = as.numeric(final_heat - twas_z),
    delta_initial = as.numeric(final_heat - twas_z),
    score_abs_delta_NESTA = abs(as.numeric(final_heat - twas_z)),
    score_abs_delta_twas = abs(as.numeric(final_heat - twas_z)),
    diagnostic_method = "M1_no_diffusion",
    stringsAsFactors = FALSE
  ) |> add_nesta_composite()
}

percentile_rank_desc <- function(x) {
  if (length(x) == 1) return(1)
  r <- rank(x, ties.method = "average")
  (r - 1) / (length(x) - 1)
}

add_nesta_composite <- function(scores) {
  scores$NESTA_M2_abs_final_heat <- abs(scores$final_NESTA_heat)
  scores$NESTA_M2_abs_delta_NESTA <- abs(scores$delta_NESTA)
  scores$score_final_percentile <- percentile_rank_desc(scores$NESTA_M2_abs_final_heat)
  scores$score_delta_percentile <- percentile_rank_desc(scores$NESTA_M2_abs_delta_NESTA)
  scores$NESTA_M2_composite <- 0.5 * scores$score_final_percentile + 0.5 * scores$score_delta_percentile
  scores
}

diagnostic_weighting_scores <- function(adj, mean_expression, twas, variant,
                                        method = "ber_p", twas_cutoff = 1,
                                        edge_cutoff = 0.1, n.perm = 300,
                                        seed = 9703, shuffle_seed = 1) {
  genes <- rownames(adj)
  if (!identical(genes, names(mean_expression))) stop("Mean-expression vector is not aligned to adjacency")
  twas_z <- strict_twas_vector(genes, twas, cutoff = twas_cutoff)
  z_scaled <- twas_z / safe_sd(twas_z)
  expr <- mean_expression
  if (variant == "D0_faithful_M2") {
    expr_factor <- expr / safe_sd(expr)
    init <- expr_factor * z_scaled
  } else if (variant == "D1_M1") {
    expr_factor <- setNames(rep(1, length(genes)), genes)
    init <- twas_z
  } else if (variant == "D2_rank_expression_weighting") {
    expr_factor <- rank(expr, ties.method = "average") / length(expr)
    names(expr_factor) <- genes
    init <- expr_factor * twas_z
  } else if (variant == "D3_centered_expression_weighting") {
    expr_factor <- as.numeric(scale(expr))
    names(expr_factor) <- genes
    init <- expr_factor * z_scaled
  } else if (variant == "D4_clipped_expression_weighting") {
    faithful_factor <- expr / safe_sd(expr)
    lim <- stats::quantile(faithful_factor, c(0.05, 0.95), na.rm = TRUE)
    expr_factor <- pmin(pmax(faithful_factor, lim[1]), lim[2])
    names(expr_factor) <- genes
    init <- expr_factor * z_scaled
  } else if (variant == "D5_shuffled_expression_negative_control") {
    set.seed(shuffle_seed)
    shuffled <- sample(expr)
    names(shuffled) <- genes
    expr_factor <- shuffled / safe_sd(shuffled)
    init <- expr_factor * z_scaled
  } else {
    stop("Unknown diagnostic weighting variant: ", variant)
  }
  names(init) <- genes
  diffuse_weighted_scores(adj, twas_z, init, method = method,
                          edge_cutoff = edge_cutoff, n.perm = n.perm,
                          seed = seed, method_label = variant,
                          expression_factor = expr_factor)
}

score_view <- function(scores, score_name) {
  out <- scores
  if (score_name == "twas_abs_z") {
    out$score_abs_delta_twas <- abs(out$TWAS.Z)
  } else if (score_name == "initial_weight_abs") {
    out$score_abs_delta_twas <- abs(out$initial_weight)
  } else if (score_name == "final_heat_abs") {
    out$score_abs_delta_twas <- abs(out$final_heat)
  } else if (score_name == "delta_twas_abs") {
    out$score_abs_delta_twas <- abs(out$delta_twas)
  } else if (score_name == "delta_initial_abs") {
    out$score_abs_delta_twas <- abs(out$delta_initial)
  } else if (score_name == "signed_final_heat") {
    out$score_abs_delta_twas <- out$final_heat
  } else if (score_name == "signed_delta_twas") {
    out$score_abs_delta_twas <- out$delta_twas
  } else if (score_name == "submitted_union") {
    r1 <- rank(-abs(out$final_heat), ties.method = "average")
    r2 <- rank(-abs(out$delta_twas), ties.method = "average")
    out$score_abs_delta_twas <- -pmin(r1, r2)
  } else {
    stop("Unknown score view: ", score_name)
  }
  out
}

score_table_view <- function(scores, score_name) {
  out <- scores
  if (score_name == "NESTA_M2_composite") {
    out$score <- out$NESTA_M2_composite
  } else if (score_name == "NESTA_M2_abs_final_heat") {
    out$score <- abs(out$final_NESTA_heat)
  } else if (score_name == "NESTA_M2_abs_delta_NESTA") {
    out$score <- abs(out$delta_NESTA)
  } else if (score_name == "raw_abs_Z") {
    out$score <- abs(out$TWAS.Z)
  } else if (score_name == "no_diffusion_M2") {
    out$score <- out$NESTA_M2_composite
  } else {
    stop("Unknown score table view: ", score_name)
  }
  out
}

auprc_score <- function(scores, positives, exclude = character()) {
  dat <- scores[!(scores$SYMBOL %in% exclude), , drop = FALSE]
  lab <- dat$SYMBOL %in% positives
  if (!any(lab) || all(lab)) return(NA_real_)
  ord <- order(dat$score_abs_delta_twas, decreasing = TRUE)
  lab <- lab[ord]
  tp <- cumsum(lab)
  fp <- cumsum(!lab)
  recall <- tp / sum(lab)
  precision <- tp / pmax(tp + fp, 1)
  recall0 <- c(0, recall)
  precision0 <- c(1, precision)
  sum((recall0[-1] - recall0[-length(recall0)]) * precision0[-1])
}

rank_metrics <- function(scores, positives, exclude = character(), ks = c(10, 25, 50, 100)) {
  dat <- scores[!(scores$SYMBOL %in% exclude), , drop = FALSE]
  ord <- order(dat$score_abs_delta_twas, decreasing = TRUE)
  ranked <- dat$SYMBOL[ord]
  pos_ranks <- match(positives, ranked)
  pos_ranks <- pos_ranks[!is.na(pos_ranks)]
  out <- data.frame(
    auprc = auprc_score(scores, positives, exclude),
    first_a2_rank = if (length(pos_ranks)) min(pos_ranks) else NA_real_,
    mean_reciprocal_rank = if (length(pos_ranks)) mean(1 / pos_ranks) else NA_real_,
    partial_auprc_top100 = auprc_score(scores[scores$SYMBOL %in% ranked[seq_len(min(100, length(ranked)))], ], positives),
    stringsAsFactors = FALSE
  )
  for (k in ks) {
    out[[paste0("top", k, "_recall")]] <- mean(positives %in% ranked[seq_len(min(k, length(ranked)))])
  }
  out
}

paired_bootstrap_ci <- function(x, nboot = 1000, seed = 20260617) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(mean = NA_real_, median = NA_real_, lo = NA_real_, hi = NA_real_))
  set.seed(seed)
  boot <- replicate(nboot, mean(sample(x, length(x), replace = TRUE)))
  c(mean = mean(x), median = stats::median(x), lo = unname(stats::quantile(boot, 0.025)),
    hi = unname(stats::quantile(boot, 0.975)))
}

auprc_from_score <- function(symbols, score, positives, exclude = character()) {
  keep <- !(symbols %in% exclude)
  symbols <- symbols[keep]
  score <- score[keep]
  lab <- symbols %in% positives
  if (!any(lab) || all(lab)) return(NA_real_)
  ord <- order(score, decreasing = TRUE)
  lab <- lab[ord]
  tp <- cumsum(lab)
  fp <- cumsum(!lab)
  recall <- tp / sum(lab)
  precision <- tp / pmax(tp + fp, 1)
  recall0 <- c(0, recall)
  precision0 <- c(1, precision)
  sum((recall0[-1] - recall0[-length(recall0)]) * precision0[-1])
}

recovery_metrics <- function(scores, positives, exclude = character(), score_col = "score",
                             primary_ks = c(50, 100), sensitivity_ks = c(25, 200)) {
  dat <- scores[!(scores$SYMBOL %in% exclude), , drop = FALSE]
  ord <- order(dat[[score_col]], decreasing = TRUE)
  ranked <- dat$SYMBOL[ord]
  pos <- positives[positives %in% ranked]
  pos_ranks <- match(pos, ranked)
  prevalence <- length(pos) / length(ranked)
  au <- auprc_from_score(dat$SYMBOL, dat[[score_col]], pos)
  out <- data.frame(
    A2_n_targets = length(pos),
    ranked_n_genes = length(ranked),
    A2_prevalence = prevalence,
    A2_AUPRC = au,
    A2_partial_AUPRC_top100 = auprc_from_score(ranked[seq_len(min(100, length(ranked)))],
                                               rev(seq_len(min(100, length(ranked)))),
                                               pos),
    A2_prevalence_normalized_AUPRC = (au - prevalence) / (1 - prevalence),
    first_A2_rank = if (length(pos_ranks)) min(pos_ranks) else NA_real_,
    A2_mean_reciprocal_rank = if (length(pos_ranks)) mean(1 / pos_ranks) else NA_real_,
    stringsAsFactors = FALSE
  )
  for (k in c(primary_ks, sensitivity_ks)) {
    recall <- mean(pos %in% ranked[seq_len(min(k, length(ranked)))])
    random_expected <- min(k, length(ranked)) / length(ranked)
    out[[paste0("A2_top", k, "_recall")]] <- recall
    out[[paste0("A2_top", k, "_fold_enrichment_over_random")]] <- recall / random_expected
  }
  out
}

make_personalization <- function(twas, genes, winsorize_top_1pct = TRUE) {
  tw <- resolve_duplicate_twas(twas)
  tw$TWAS.P <- twas_p_from_z(tw$TWAS.Z)
  p <- setNames(rep(0, length(genes)), genes)
  idx <- match(tw$SYMBOL, genes)
  keep <- !is.na(idx)
  raw <- -log10(tw$TWAS.P[keep])
  finite <- is.finite(raw)
  if (any(!finite) && any(finite)) raw[!finite] <- max(raw[finite])
  raw[!is.finite(raw)] <- 0
  p[idx[keep]] <- raw
  if (winsorize_top_1pct && any(p > 0)) {
    cap <- stats::quantile(p, 0.99, na.rm = TRUE)
    p <- pmin(p, cap)
  }
  p[!is.finite(p)] <- 0
  s <- sum(p)
  if (!is.finite(s) || s <= 0) stop("Personalization vector has non-positive sum")
  p / s
}

make_signed_personalization <- function(twas, genes, winsorize_top_1pct = TRUE) {
  tw <- resolve_duplicate_twas(twas)
  tw$TWAS.P <- twas_p_from_z(tw$TWAS.Z)
  mag <- setNames(rep(0, length(genes)), genes)
  idx <- match(tw$SYMBOL, genes)
  keep <- !is.na(idx)
  raw <- -log10(tw$TWAS.P[keep])
  finite <- is.finite(raw)
  if (any(!finite) && any(finite)) raw[!finite] <- max(raw[finite])
  raw[!is.finite(raw)] <- 0
  mag[idx[keep]] <- raw
  if (winsorize_top_1pct && any(mag > 0)) {
    cap <- stats::quantile(mag, 0.99, na.rm = TRUE)
    mag <- pmin(mag, cap)
  }
  z <- setNames(rep(0, length(genes)), genes)
  z[idx[keep]] <- tw$TWAS.Z[keep]
  pos <- ifelse(z > 0, mag, 0)
  neg <- ifelse(z < 0, mag, 0)
  normalize <- function(x) {
    x[!is.finite(x)] <- 0
    s <- sum(x)
    if (!is.finite(s) || s <= 0) return(x)
    x / s
  }
  list(pos = normalize(pos), neg = normalize(neg))
}

weighted_transition_matrix <- function(adj, edge_cutoff = 0.1, normalize = "column") {
  m <- adj
  m[m <= edge_cutoff] <- 0
  diag(m) <- 0
  if (!any(m > 0)) stop("Weighted graph has no positive edges after cutoff")
  if (normalize == "column") {
    cs <- colSums(m)
    cs[cs == 0] <- 1
    sweep(m, 2, cs, "/")
  } else if (normalize == "row") {
    rs <- rowSums(m)
    rs[rs == 0] <- 1
    sweep(m, 1, rs, "/")
  } else {
    stop("Unknown normalization: ", normalize)
  }
}

run_rwr <- function(adj, p, restart = 0.50, tol = 1e-10, max_iter = 10000,
                    edge_cutoff = 0.1) {
  genes <- rownames(adj)
  if (!identical(genes, names(p))) stop("RWR personalization vector is not aligned")
  W <- weighted_transition_matrix(adj, edge_cutoff = edge_cutoff, normalize = "column")
  s <- p
  for (i in seq_len(max_iter)) {
    nxt <- as.numeric((1 - restart) * (W %*% s) + restart * p)
    names(nxt) <- genes
    if (sum(abs(nxt - s)) < tol) return(nxt)
    s <- nxt
  }
  s
}

run_ppr <- function(adj, p, damping = 0.85, tol = 1e-10, max_iter = 10000,
                    edge_cutoff = 0.1) {
  genes <- rownames(adj)
  if (!identical(genes, names(p))) stop("PPR personalization vector is not aligned")
  W <- weighted_transition_matrix(adj, edge_cutoff = edge_cutoff, normalize = "column")
  s <- p
  restart <- 1 - damping
  for (i in seq_len(max_iter)) {
    nxt <- as.numeric(damping * (W %*% s) + restart * p)
    names(nxt) <- genes
    if (sum(abs(nxt - s)) < tol) return(nxt)
    s <- nxt
  }
  s
}

common_prior_nesta_scores <- function(adj, twas, method = "ber_p", edge_cutoff = 0.1,
                                      n.perm = 300, seed = 9703) {
  genes <- rownames(adj)
  p <- make_personalization(twas, genes)
  g <- make_faithful_graph(adj, edge_cutoff = edge_cutoff, retain_zero_edges = TRUE,
                           submitted_diagonal = TRUE)
  K <- regularised_kernel(g)
  heat <- as.numeric(diffuse_checked(g, p, method = method, n.perm = n.perm, seed = seed, K = K))
  heat <- heat / safe_sd(heat)
  data.frame(SYMBOL = genes, score = heat, stringsAsFactors = FALSE)
}

benchmark_scores <- function(adj, twas, include_sensitivity = TRUE) {
  genes <- rownames(adj)
  p <- make_personalization(twas, genes)
  ps <- make_signed_personalization(twas, genes)
  rwr_signed <- run_rwr(adj, ps$pos, restart = 0.50) - run_rwr(adj, ps$neg, restart = 0.50)
  ppr_signed <- run_ppr(adj, ps$pos, damping = 0.85) - run_ppr(adj, ps$neg, damping = 0.85)
  out <- list(
    raw_minus_log10P_abs_prior = data.frame(SYMBOL = genes, score = as.numeric(p), score_signed = as.numeric(p), stringsAsFactors = FALSE),
    RWR_abs_prior = data.frame(SYMBOL = genes, score = as.numeric(run_rwr(adj, p, restart = 0.50)), score_signed = NA_real_, stringsAsFactors = FALSE),
    PPR_abs_prior = data.frame(SYMBOL = genes, score = as.numeric(run_ppr(adj, p, damping = 0.85)), score_signed = NA_real_, stringsAsFactors = FALSE),
    RWR_signed_two_channel = data.frame(SYMBOL = genes, score = abs(as.numeric(rwr_signed)), score_signed = as.numeric(rwr_signed), stringsAsFactors = FALSE),
    PPR_signed_two_channel = data.frame(SYMBOL = genes, score = abs(as.numeric(ppr_signed)), score_signed = as.numeric(ppr_signed), stringsAsFactors = FALSE),
    NESTA_common_prior = common_prior_nesta_scores(adj, twas)
  )
  out$NESTA_common_prior$score_signed <- out$NESTA_common_prior$score
  if (include_sensitivity) {
    out$RWR_abs_prior_restart_030 <- data.frame(SYMBOL = genes, score = as.numeric(run_rwr(adj, p, restart = 0.30)), score_signed = NA_real_, stringsAsFactors = FALSE)
    out$RWR_abs_prior_restart_070 <- data.frame(SYMBOL = genes, score = as.numeric(run_rwr(adj, p, restart = 0.70)), score_signed = NA_real_, stringsAsFactors = FALSE)
    out$PPR_abs_prior_damping_075 <- data.frame(SYMBOL = genes, score = as.numeric(run_ppr(adj, p, damping = 0.75)), score_signed = NA_real_, stringsAsFactors = FALSE)
    out$PPR_abs_prior_damping_095 <- data.frame(SYMBOL = genes, score = as.numeric(run_ppr(adj, p, damping = 0.95)), score_signed = NA_real_, stringsAsFactors = FALSE)
    r030 <- run_rwr(adj, ps$pos, restart = 0.30) - run_rwr(adj, ps$neg, restart = 0.30)
    r070 <- run_rwr(adj, ps$pos, restart = 0.70) - run_rwr(adj, ps$neg, restart = 0.70)
    d075 <- run_ppr(adj, ps$pos, damping = 0.75) - run_ppr(adj, ps$neg, damping = 0.75)
    d095 <- run_ppr(adj, ps$pos, damping = 0.95) - run_ppr(adj, ps$neg, damping = 0.95)
    out$RWR_signed_two_channel_restart_030 <- data.frame(SYMBOL = genes, score = abs(as.numeric(r030)), score_signed = as.numeric(r030), stringsAsFactors = FALSE)
    out$RWR_signed_two_channel_restart_070 <- data.frame(SYMBOL = genes, score = abs(as.numeric(r070)), score_signed = as.numeric(r070), stringsAsFactors = FALSE)
    out$PPR_signed_two_channel_damping_075 <- data.frame(SYMBOL = genes, score = abs(as.numeric(d075)), score_signed = as.numeric(d075), stringsAsFactors = FALSE)
    out$PPR_signed_two_channel_damping_095 <- data.frame(SYMBOL = genes, score = abs(as.numeric(d095)), score_signed = as.numeric(d095), stringsAsFactors = FALSE)
  }
  out
}
