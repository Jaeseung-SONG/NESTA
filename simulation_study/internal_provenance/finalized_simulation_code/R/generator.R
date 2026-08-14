source(file.path("/home/js/NESTA/simulation", "R", "utils.R"))
source(file.path("/home/js/NESTA/simulation", "R", "fidelity.R"))

sample_weights <- function(pool, n) sample(pool, n, replace = TRUE)

sym_from_pool <- function(genes, pool) {
  n <- length(genes)
  m <- matrix(0, n, n, dimnames = list(genes, genes))
  m[lower.tri(m)] <- sample_weights(pool, n * (n - 1) / 2)
  m <- m + t(m)
  diag(m) <- 1
  m
}

choose_module_a <- function(template) {
  genes <- template$genes
  tom <- template$tom
  strength <- rowSums(tom)
  g <- igraph::graph_from_adjacency_matrix((tom > 0.1) * 1, mode = "undirected", diag = FALSE)
  comps <- igraph::components(g)
  largest <- which.max(comps$csize)
  keep <- which(comps$membership == largest)
  if (length(keep) < 30) {
    stop("Template lacks a 30-gene connected Module A component at TOM threshold 0.1")
  }
  if (length(keep) > 150) {
    subg <- igraph::induced_subgraph(g, keep)
    centrality <- suppressWarnings(igraph::closeness(subg, normalized = TRUE))
    if (all(!is.finite(centrality))) centrality <- rowSums(tom[keep, keep, drop = FALSE])
    center_local <- which.max(centrality)
    d <- as.numeric(igraph::distances(subg, v = center_local)[1, ])
    ord <- order(d, -rowSums(tom[keep, keep, drop = FALSE]), na.last = NA)
    keep <- keep[ord[seq_len(min(150, length(ord)))]]
  }
  genes <- genes[keep]
  tom <- tom[keep, keep, drop = FALSE]
  list(genes = genes, tom = tom, expression = template$expression[match(genes, template$genes)])
}

generate_twas <- function(genes, a1_risk, a1_protective, a2_risk = character(),
                          a2_protective = character(), decoy = character(),
                          decoy_direction = NULL, seed, null = FALSE) {
  set.seed(seed)
  z <- stats::rnorm(length(genes), 0, 1)
  names(z) <- genes
  if (!null) {
    z[a1_risk] <- abs(stats::rnorm(length(a1_risk), mean = 3.5, sd = 0.75))
    z[a1_protective] <- -abs(stats::rnorm(length(a1_protective), mean = 3.5, sd = 0.75))
    if (length(a2_risk)) z[a2_risk] <- stats::rnorm(length(a2_risk), mean = 0.25, sd = 1)
    if (length(a2_protective)) z[a2_protective] <- stats::rnorm(length(a2_protective), mean = -0.25, sd = 1)
    if (length(decoy)) {
      dd <- decoy_direction[decoy]
      z[decoy[dd == "risk"]] <- abs(stats::rnorm(sum(dd == "risk"), mean = 3.0, sd = 0.6))
      z[decoy[dd == "protective"]] <- -abs(stats::rnorm(sum(dd == "protective"), mean = 3.0, sd = 0.6))
    }
  }
  data.frame(SYMBOL = genes, TWAS.Z = as.numeric(z), TWAS.P = twas_p_from_z(z),
             stringsAsFactors = FALSE)
}

non_extreme_seed_genes <- function(module_adj, genes, n = 5) {
  g <- igraph::graph_from_adjacency_matrix((module_adj > 0.1) * 1, mode = "undirected", diag = FALSE)
  comps <- igraph::components(g)
  largest <- which.max(comps$csize)
  genes <- genes[comps$membership == largest]
  module_adj <- module_adj[genes, genes, drop = FALSE]
  deg <- rowSums(module_adj > 0.1) - 1
  close <- suppressWarnings(igraph::closeness(igraph::induced_subgraph(g, which(comps$membership == largest)), normalized = TRUE))
  ok <- deg >= stats::quantile(deg, 0.2) & deg <= stats::quantile(deg, 0.8)
  pool <- genes[ok]
  if (length(pool) < n) pool <- genes
  pool <- pool[order(close[match(pool, genes)], decreasing = TRUE, na.last = NA)]
  if (length(pool) < n) pool <- genes
  sample(pool[seq_len(min(length(pool), max(n, 20)))], n)
}

assign_a1_directions <- function(a1, rep_id) {
  if (rep_id %% 2 == 1) {
    list(risk = a1[seq_len(min(3, length(a1)))],
         protective = setdiff(a1, a1[seq_len(min(3, length(a1)))]),
         majority = "risk")
  } else {
    list(protective = a1[seq_len(min(3, length(a1)))],
         risk = setdiff(a1, a1[seq_len(min(3, length(a1)))]),
         majority = "protective")
  }
}

graph_distances_from_a1 <- function(adj, module_genes, a1, threshold = 0.1) {
  sub <- adj[module_genes, module_genes, drop = FALSE]
  g <- igraph::graph_from_adjacency_matrix((sub > threshold) * 1, mode = "undirected", diag = FALSE)
  d <- suppressWarnings(igraph::distances(g, v = match(a1, module_genes), to = seq_along(module_genes)))
  mind <- apply(d, 2, min, na.rm = TRUE)
  names(mind) <- module_genes
  mind
}

path_fraction_table <- function(paths) {
  paths <- paths[is.finite(paths)]
  n <- length(paths)
  if (!n) {
    return(data.frame(direct_1hop_fraction = NA_real_, two_hop_fraction = NA_real_,
                      three_plus_hop_fraction = NA_real_, median_A1_A2_path = Inf,
                      A2_same_component_fraction = 0))
  }
  data.frame(direct_1hop_fraction = mean(paths == 1),
             two_hop_fraction = mean(paths == 2),
             three_plus_hop_fraction = mean(paths >= 3),
             median_A1_A2_path = stats::median(paths),
             A2_same_component_fraction = mean(is.finite(paths)))
}

thin_a1_shortcuts <- function(adj, a_genes, a1, seed, threshold = 0.1) {
  set.seed(seed)
  out <- adj
  nonseed <- setdiff(a_genes, a1)
  thinned <- character()
  add_ring <- function(m) {
    ring_w <- max(0.12, stats::quantile(m[a_genes, a_genes][lower.tri(m[a_genes, a_genes])], 0.65, na.rm = TRUE))
    for (i in seq_along(a_genes)) {
      j <- if (i == length(a_genes)) 1 else i + 1
      m[a_genes[i], a_genes[j]] <- max(m[a_genes[i], a_genes[j]], ring_w)
      m[a_genes[j], a_genes[i]] <- max(m[a_genes[j], a_genes[i]], ring_w)
    }
    m
  }
  out <- add_ring(out)
  direct <- nonseed[apply(out[a1, nonseed, drop = FALSE] > threshold, 2, any)]
  if (length(direct)) {
    target <- sample(direct, ceiling(length(direct) * 0.65))
    for (g in target) {
      hit <- a1[out[a1, g] > threshold]
      if (length(hit)) {
        out[hit, g] <- pmin(out[hit, g], threshold * 0.95)
        out[g, hit] <- out[hit, g]
        thinned <- c(thinned, paste(hit, g, sep = "--"))
      }
    }
  }
  diag(out) <- 1
  list(adj = out, thinned_edges = unique(thinned))
}

assign_a2_direction <- function(adj, module_genes, a2, a1_risk, a1_protective, seed,
                                threshold = 0.1) {
  set.seed(seed)
  sub <- adj[module_genes, module_genes, drop = FALSE]
  g <- igraph::graph_from_adjacency_matrix((sub > threshold) * 1, mode = "undirected", diag = FALSE)
  dr <- suppressWarnings(igraph::distances(g, v = match(a1_risk, module_genes), to = match(a2, module_genes)))
  dp <- suppressWarnings(igraph::distances(g, v = match(a1_protective, module_genes), to = match(a2, module_genes)))
  mr <- apply(dr, 2, min, na.rm = TRUE)
  mp <- apply(dp, 2, min, na.rm = TRUE)
  dir <- ifelse(mr < mp, "risk", ifelse(mp < mr, "protective", sample(c("risk", "protective"), length(a2), replace = TRUE)))
  names(dir) <- a2
  dir
}

select_path_stratified_targets <- function(adj, module_genes, a1, seed, threshold = 0.1) {
  shaped <- thin_a1_shortcuts(adj, module_genes, a1, seed = seed, threshold = threshold)
  adj2 <- shaped$adj
  nonseed <- setdiff(module_genes, a1)
  paths <- graph_distances_from_a1(adj2, module_genes, a1, threshold = threshold)[nonseed]
  pool1 <- names(paths)[is.finite(paths) & paths == 1]
  pool2 <- names(paths)[is.finite(paths) & paths == 2]
  pool3 <- names(paths)[is.finite(paths) & paths >= 3]
  target_n <- min(length(nonseed), max(30, min(70, floor(length(nonseed) * 0.65))))
  set.seed(seed + 17)
  n2 <- min(length(pool2), ceiling(0.65 * target_n))
  n1 <- min(length(pool1), max(0, ceiling(0.30 * target_n)))
  n3 <- min(length(pool3), max(0, target_n - n1 - n2))
  chosen <- c(sample_or_all(pool1, n1), sample_or_all(pool2, n2), sample_or_all(pool3, n3))
  fill <- setdiff(names(paths)[is.finite(paths)], chosen)
  if (length(chosen) < target_n && length(fill)) chosen <- c(chosen, sample_or_all(fill, target_n - length(chosen)))
  chosen <- unique(chosen)
  paths_chosen <- paths[chosen]
  metrics <- path_fraction_table(paths_chosen)
  metrics$path_stratification_pass <- with(metrics,
    is.finite(median_A1_A2_path) & median_A1_A2_path >= 2 & median_A1_A2_path <= 4 &
      A2_same_component_fraction >= 0.85 &
      ((direct_1hop_fraction <= 0.60 & two_hop_fraction >= 0.50) |
         (direct_1hop_fraction <= 0.55 & two_hop_fraction >= 0.40)))
  metrics$path_fallback_used <- with(metrics, !(direct_1hop_fraction <= 0.60 & two_hop_fraction >= 0.50))
  list(adj = adj2,
       A2_all = chosen,
       A2_proximal = chosen[paths_chosen == 1],
       A2_intermediate = chosen[paths_chosen == 2],
       A2_distal = chosen[paths_chosen >= 3],
       paths = paths_chosen,
       thinned_edges = shaped$thinned_edges,
       metrics = metrics)
}

module_degree_table <- function(adj, module_genes, threshold = 0.1) {
  sub <- adj[module_genes, module_genes, drop = FALSE] > threshold
  diag(sub) <- FALSE
  degree <- rowSums(sub)
  pct <- percentile_rank_desc(degree) * 100
  bin <- ifelse(pct <= 40, "low_degree",
                ifelse(pct <= 80, "mid_degree",
                       ifelse(pct <= 90, "high_degree", "extreme_high_degree")))
  data.frame(SYMBOL = module_genes, module_a_degree = as.numeric(degree),
             degree_percentile = as.numeric(pct), degree_bin = bin,
             stringsAsFactors = FALSE)
}

whole_degree_percentile <- function(adj, threshold = 0.1) {
  bin <- adj > threshold
  diag(bin) <- FALSE
  deg <- rowSums(bin)
  out <- percentile_rank_desc(deg) * 100
  names(out) <- rownames(adj)
  out
}

select_degree_capped_primary <- function(two_hop, degree_tab, seed, target_n = 30) {
  set.seed(seed)
  dt <- degree_tab[match(two_hop, degree_tab$SYMBOL), , drop = FALSE]
  dt <- dt[!is.na(dt$SYMBOL), , drop = FALSE]
  low <- dt$SYMBOL[dt$degree_bin == "low_degree"]
  mid <- dt$SYMBOL[dt$degree_bin == "mid_degree"]
  high <- dt$SYMBOL[dt$degree_bin == "high_degree"]
  extreme <- dt$SYMBOL[dt$degree_bin == "extreme_high_degree"]
  n0 <- min(target_n, nrow(dt))
  if (!n0) return(character())
  n_low <- min(length(low), ceiling(0.30 * n0))
  n_mid <- min(length(mid), ceiling(0.40 * n0))
  cap_extreme <- min(7, ceiling(0.25 * n0))
  n_high <- min(length(high), max(ceiling(0.15 * n0), n0 - n_low - n_mid - cap_extreme))
  n_extreme <- min(length(extreme), cap_extreme, max(0, n0 - n_low - n_mid - n_high))
  chosen <- c(sample_or_all(low, n_low), sample_or_all(mid, n_mid),
              sample_or_all(high, n_high), sample_or_all(extreme, n_extreme))
  fill <- setdiff(c(mid, high, low), chosen)
  if (length(chosen) < n0 && length(fill)) {
    chosen <- c(chosen, sample_or_all(fill, n0 - length(chosen)))
  }
  unique(chosen)
}

degree_audit_metrics <- function(primary, degree_tab, background_genes, c_decoy) {
  dt <- degree_tab[match(primary, degree_tab$SYMBOL), , drop = FALSE]
  bg <- degree_tab[degree_tab$SYMBOL %in% background_genes, , drop = FALSE]
  cd <- degree_tab[degree_tab$SYMBOL %in% c_decoy, , drop = FALSE]
  n_primary <- length(primary)
  n_ext <- sum(dt$degree_bin == "extreme_high_degree", na.rm = TRUE)
  data.frame(
    n_A2_primary = n_primary,
    n_A2_low_degree = sum(dt$degree_bin == "low_degree", na.rm = TRUE),
    n_A2_mid_degree = sum(dt$degree_bin == "mid_degree", na.rm = TRUE),
    n_A2_high_degree = sum(dt$degree_bin == "high_degree", na.rm = TRUE),
    n_A2_extreme_high_degree = n_ext,
    extreme_high_degree_A2_count = n_ext,
    fraction_A2_extreme_high_degree = if (n_primary) n_ext / n_primary else NA_real_,
    median_A2_degree_percentile = stats::median(dt$degree_percentile, na.rm = TRUE),
    median_background_degree_percentile = stats::median(bg$degree_percentile, na.rm = TRUE),
    median_C_high_degree_decoy_percentile = stats::median(cd$degree_percentile, na.rm = TRUE),
    A2_vs_background_degree_KS = if (nrow(dt) && nrow(bg)) ks_stat(dt$degree_percentile, bg$degree_percentile) else NA_real_,
    A2_vs_C_decoy_degree_KS = if (nrow(dt) && nrow(cd)) ks_stat(dt$degree_percentile, cd$degree_percentile) else NA_real_,
    A2_low_mid_fraction = if (n_primary) mean(dt$degree_bin %in% c("low_degree", "mid_degree"), na.rm = TRUE) else NA_real_,
    A2_moderate_connector_fraction = if (n_primary) mean(dt$degree_bin %in% c("mid_degree", "high_degree"), na.rm = TRUE) else NA_real_,
    median_A2_degree_in_target_range = is.finite(stats::median(dt$degree_percentile, na.rm = TRUE)) &&
      stats::median(dt$degree_percentile, na.rm = TRUE) >= 45 &&
      stats::median(dt$degree_percentile, na.rm = TRUE) <= 65,
    median_A2_degree_in_diagnostic_range = is.finite(stats::median(dt$degree_percentile, na.rm = TRUE)) &&
      stats::median(dt$degree_percentile, na.rm = TRUE) >= 35 &&
      stats::median(dt$degree_percentile, na.rm = TRUE) <= 75,
    degree_distribution_qc_pass = n_primary >= 12 &&
      n_ext <= min(7, ceiling(0.25 * max(n_primary, 1))),
    stringsAsFactors = FALSE
  )
}

construct_high_degree_decoy <- function(adj, universe, degree_tab, seed, n = 50, threshold = 0.1) {
  set.seed(seed)
  pct <- whole_degree_percentile(adj, threshold = threshold)
  pool <- setdiff(c(universe$module_c, universe$background, universe$module_b), universe$module_a)
  pool <- pool[pool %in% names(pct)]
  ord <- pool[order(pct[pool], decreasing = TRUE, na.last = NA)]
  sample_or_all(ord, min(n, length(ord)))
}

construct_decoy_module <- function(adj, universe, a1_risk, a1_protective, seed,
                                   n = 50, threshold = 0.1) {
  set.seed(seed)
  pool <- setdiff(universe$genes, universe$module_a)
  candidates <- intersect(c(universe$module_b, universe$module_c, universe$background), pool)
  strength_to_seed <- rowSums(adj[candidates, c(a1_risk, a1_protective), drop = FALSE])
  ord <- candidates[order(strength_to_seed, decreasing = TRUE, na.last = NA)]
  decoy <- sample_or_all(ord, min(n, length(ord)))
  branch <- ifelse(rowSums(adj[decoy, a1_risk, drop = FALSE]) >= rowSums(adj[decoy, a1_protective, drop = FALSE]),
                   "near_risk", "near_protective")
  direction <- ifelse(branch == "near_risk", "protective", "risk")
  names(direction) <- decoy
  list(genes = decoy, branch = setNames(branch, decoy), direction = direction)
}

make_universe <- function(lib, template_key, seed, n_genes = 1000) {
  set.seed(seed)
  template <- choose_module_a(lib$template_payload[[template_key]])
  cell_type <- lib$template_payload[[template_key]]$cell_type
  cell_payload <- lib$cell_payload[[cell_type]]
  a_genes <- template$genes
  available <- setdiff(cell_payload$all_genes, a_genes)
  expr_pool <- cell_payload$expression_pool
  expr_named <- expr_pool
  high_expr_genes <- intersect(names(sort(expr_named, decreasing = TRUE)), available)
  b_genes <- sample_or_all(high_expr_genes, 50)
  remaining <- setdiff(available, b_genes)
  deg_named <- cell_payload$degree
  high_deg_genes <- intersect(names(sort(deg_named, decreasing = TRUE)), remaining)
  c_genes <- sample_or_all(high_deg_genes, 50)
  remaining <- setdiff(remaining, c_genes)
  bg_n <- n_genes - length(a_genes) - length(b_genes) - length(c_genes)
  if (bg_n < 0) stop("Requested network universe is smaller than required A/B/C blocks")
  if (length(remaining) >= bg_n) {
    bg_genes <- sample(remaining, bg_n)
  } else {
    extra_n <- bg_n - length(remaining)
    bg_genes <- c(remaining, sprintf("synthetic_universe_background_%05d", seq_len(extra_n)))
  }
  genes <- c(a_genes, b_genes, c_genes, bg_genes)
  mean_expression <- expr_named[genes]
  miss <- !is.finite(mean_expression)
  if (any(miss)) mean_expression[miss] <- sample(expr_pool, sum(miss), replace = TRUE)
  names(mean_expression) <- genes
  list(genes = genes, mean_expression = mean_expression, module_a = a_genes,
       module_b = b_genes, module_c = c_genes, background = bg_genes,
       module_a_tom = template$tom, cell_type = cell_type)
}

embed_module <- function(adj, module_genes, module_tom) {
  adj[module_genes, module_genes] <- module_tom[module_genes, module_genes]
  diag(adj) <- 1
  adj
}

boost_module_c <- function(adj, c_genes, pool) {
  if (length(c_genes) < 2) return(adj)
  n <- length(c_genes)
  vals <- stats::quantile(pool, probs = c(0.8, 0.99), na.rm = TRUE)
  high_pool <- pool[pool >= vals[1] & pool <= vals[2]]
  block <- matrix(0, n, n, dimnames = list(c_genes, c_genes))
  block[lower.tri(block)] <- sample(high_pool, n * (n - 1) / 2, replace = TRUE)
  block <- block + t(block)
  diag(block) <- 1
  adj[c_genes, c_genes] <- block
  adj
}

make_relevant_network <- function(lib, universe) {
  pool <- lib$cell_payload[[universe$cell_type]]$weight_pool
  adj <- sym_from_pool(universe$genes, pool)
  if (!is.null(universe$background_sparsity) && identical(universe$background_sparsity, "zero_proximal_background")) {
    bg <- universe$background
    bg <- intersect(bg, rownames(adj))
    if (length(bg)) {
      adj[bg, bg] <- pmin(adj[bg, bg] * 0.25, 0.095)
      outside_branch <- setdiff(rownames(adj), universe$module_a)
      adj[bg, outside_branch] <- pmin(adj[bg, outside_branch] * 0.35, 0.095)
      adj[outside_branch, bg] <- t(adj[bg, outside_branch, drop = FALSE])
    }
  }
  adj <- embed_module(adj, universe$module_a, universe$module_a_tom)
  adj <- boost_module_c(adj, universe$module_c, pool)
  adj
}

make_i1_network <- function(lib, universe, seed) {
  set.seed(seed)
  other <- setdiff(names(lib$cell_payload), universe$cell_type)
  ct <- sample(other, 1)
  adj <- sym_from_pool(universe$genes, lib$cell_payload[[ct]]$weight_pool)
  adj <- boost_module_c(adj, universe$module_c, lib$cell_payload[[ct]]$weight_pool)
  adj
}

make_i2_network <- function(lib, universe, relevant, seed) {
  set.seed(seed)
  adj <- relevant
  pool <- lib$cell_payload[[universe$cell_type]]$weight_pool
  a <- universe$module_a
  bg <- sample(universe$background, length(a))
  adj[bg, bg] <- universe$module_a_tom
  low_pool <- pool[pool <= stats::quantile(pool, 0.55, na.rm = TRUE)]
  repl <- matrix(0, length(a), length(a), dimnames = list(a, a))
  repl[lower.tri(repl)] <- sample(low_pool, length(a) * (length(a) - 1) / 2, replace = TRUE)
  repl <- repl + t(repl)
  ring_w <- max(0.11, stats::quantile(pool, 0.65, na.rm = TRUE))
  for (i in seq_along(a)) {
    j <- if (i == length(a)) 1 else i + 1
    repl[i, j] <- ring_w
    repl[j, i] <- ring_w
  }
  diag(repl) <- 1
  adj[a, a] <- repl
  diag(adj) <- 1
  adj
}

make_i3_network <- function(lib, universe, relevant, seed) {
  set.seed(seed)
  adj <- relevant
  a <- universe$module_a
  pool <- lib$cell_payload[[universe$cell_type]]$weight_pool
  low_pool <- pool[pool <= stats::quantile(pool, 0.55, na.rm = TRUE)]
  repl <- matrix(0, length(a), length(a), dimnames = list(a, a))
  repl[lower.tri(repl)] <- sample(low_pool, length(a) * (length(a) - 1) / 2, replace = TRUE)
  repl <- repl + t(repl)
  ring_w <- max(0.11, stats::quantile(pool, 0.65, na.rm = TRUE))
  for (i in seq_along(a)) {
    j <- if (i == length(a)) 1 else i + 1
    repl[i, j] <- ring_w
    repl[j, i] <- ring_w
  }
  diag(repl) <- 1
  adj[a, a] <- repl
  bg <- sample(universe$background, length(a))
  adj[bg, bg] <- universe$module_a_tom
  diag(adj) <- 1
  adj
}

permute_weights <- function(adj, seed) {
  set.seed(seed)
  vals <- adj[lower.tri(adj)]
  m <- matrix(0, nrow(adj), ncol(adj), dimnames = dimnames(adj))
  m[lower.tri(m)] <- sample(vals)
  m <- m + t(m)
  diag(m) <- 1
  m
}

degree_preserving_rewire <- function(adj, seed, threshold = 0.1) {
  set.seed(seed)
  genes <- rownames(adj)
  bin <- adj > threshold
  diag(bin) <- FALSE
  g <- igraph::graph_from_adjacency_matrix(bin * 1, mode = "undirected", diag = FALSE)
  g2 <- igraph::rewire(g, with = igraph::keeping_degseq(niter = max(1000, igraph::ecount(g) * 2)))
  out <- matrix(0, nrow(adj), ncol(adj), dimnames = list(genes, genes))
  el <- igraph::as_edgelist(g2, names = TRUE)
  weights <- adj[lower.tri(adj)]
  weights <- weights[weights > threshold]
  if (nrow(el)) {
    w <- sample(weights, nrow(el), replace = TRUE)
    out[cbind(el[,1], el[,2])] <- w
    out[cbind(el[,2], el[,1])] <- w
  }
  diag(out) <- 1
  out
}

network_metrics <- function(adj, modules, mean_expression, threshold = 0.1) {
  bin <- adj > threshold
  diag(bin) <- FALSE
  degree <- rowSums(bin)
  strength <- rowSums(adj * bin)
  a <- modules$A
  bg <- modules$background
  sub <- adj[a, a, drop = FALSE]
  g <- igraph::graph_from_adjacency_matrix((sub > threshold) * 1, mode = "undirected", diag = FALSE)
  whole_g <- igraph::graph_from_adjacency_matrix(bin * 1, mode = "undirected", diag = FALSE)
  whole_comp <- igraph::components(whole_g)
  a_comp <- igraph::components(g)
  a1_idx <- match(modules$A1, a)
  a2_for_paths <- if (!is.null(modules$A2_all)) modules$A2_all else modules$A2
  a2_idx <- match(a2_for_paths, a)
  cl <- suppressWarnings(mean(igraph::transitivity(g, type = "localundirected", isolates = "zero")))
  bg_cl <- NA_real_
  if (length(bg) > 3) {
    bg_sample <- bg[seq_len(min(length(bg), 300))]
    bg_g <- igraph::graph_from_adjacency_matrix((adj[bg_sample, bg_sample, drop = FALSE] > threshold) * 1,
                                                mode = "undirected", diag = FALSE)
    bg_cl_vals <- suppressWarnings(igraph::transitivity(bg_g, type = "localundirected", isolates = "zero"))
    bg_cl <- stats::median(bg_cl_vals, na.rm = TRUE)
    bg_cl_p40 <- stats::quantile(bg_cl_vals, 0.40, na.rm = TRUE)
  } else {
    bg_cl_p40 <- NA_real_
  }
  paths <- suppressWarnings(igraph::distances(g))
  finite_paths <- paths[is.finite(paths) & paths > 0]
  a1_a2_paths <- as.numeric(paths[a1_idx, a2_idx, drop = FALSE])
  a1_a2_paths <- a1_a2_paths[is.finite(a1_a2_paths)]
  same_comp <- rep(FALSE, length(a2_idx))
  for (i in seq_along(a2_idx)) {
    same_comp[i] <- any(a_comp$membership[a2_idx[i]] == a_comp$membership[a1_idx])
  }
  data.frame(
    density = mean(bin[lower.tri(bin)]),
    isolated_fraction = mean(degree == 0),
    largest_connected_component_fraction = max(whole_comp$csize) / nrow(adj),
    module_a_connected_component_fraction = max(a_comp$csize) / length(a),
    finite_edge_weights = all(is.finite(adj)),
    duplicated_gene_names = anyDuplicated(rownames(adj)) > 0,
    median_degree = stats::median(degree),
    median_strength = stats::median(strength),
    within_a_tom = stats::median(sub[lower.tri(sub)]),
    between_a_tom = stats::median(adj[a, setdiff(rownames(adj), a), drop = FALSE]),
    within_between_a_ratio = stats::median(sub[lower.tri(sub)]) / stats::median(adj[a, setdiff(rownames(adj), a), drop = FALSE]),
    clustering_a = ifelse(is.finite(cl), cl, 0),
    background_clustering_median = ifelse(is.finite(bg_cl), bg_cl, 0),
    background_clustering_p40 = ifelse(is.finite(bg_cl_p40), bg_cl_p40, 0),
    median_a_path = if (length(finite_paths)) stats::median(finite_paths) else NA_real_,
    median_a1_a2_path = if (length(a1_a2_paths)) stats::median(a1_a2_paths) else Inf,
    a1_a2_same_component_fraction = mean(same_comp),
    direct_1hop_fraction = if (length(a1_a2_paths)) mean(apply(paths[a1_idx, a2_idx, drop = FALSE], 2, min, na.rm = TRUE) == 1) else NA_real_,
    two_hop_fraction = if (length(a1_a2_paths)) mean(apply(paths[a1_idx, a2_idx, drop = FALSE], 2, min, na.rm = TRUE) == 2) else NA_real_,
    three_plus_hop_fraction = if (length(a1_a2_paths)) mean(apply(paths[a1_idx, a2_idx, drop = FALSE], 2, min, na.rm = TRUE) >= 3) else NA_real_,
    expr_degree_cor = suppressWarnings(stats::cor(mean_expression, degree, method = "spearman")),
    expr_strength_cor = suppressWarnings(stats::cor(mean_expression, strength, method = "spearman")),
    stringsAsFactors = FALSE
  )
}

matching_qc <- function(relevant, other, mean_expression, modules = NULL) {
  bin_r <- relevant > 0.1
  bin_o <- other > 0.1
  diag(bin_r) <- FALSE
  diag(bin_o) <- FALSE
  strength_r <- rowSums(relevant * bin_r)
  strength_o <- rowSums(other * bin_o)
  out <- data.frame(
    density_diff = abs(mean(bin_r[lower.tri(bin_r)]) - mean(bin_o[lower.tri(bin_o)])),
    median_strength_ratio = stats::median(strength_o) / stats::median(strength_r),
    expression_ks = ks_stat(mean_expression, mean_expression),
    identical_twas = TRUE,
    stringsAsFactors = FALSE
  )
  if (!is.null(modules)) {
    mr <- network_metrics(relevant, modules, mean_expression)
    mo <- network_metrics(other, modules, mean_expression)
    out$within_a_tom_reduction_fraction <- 1 - (mo$within_a_tom / mr$within_a_tom)
    out$a1_a2_path_increase <- mo$median_a1_a2_path - mr$median_a1_a2_path
    out$same_component_fraction_reduction <- 1 - (mo$a1_a2_same_component_fraction / mr$a1_a2_same_component_fraction)
    out$module_disrupted <- out$within_a_tom_reduction_fraction >= 0.25 |
      out$a1_a2_path_increase >= 1 |
      out$same_component_fraction_reduction >= 0.25
  }
  out
}

topology_qc_pass <- function(qc) {
  with(qc, all(density_diff <= 0.025,
               median_strength_ratio >= 0.45,
               median_strength_ratio <= 2.20,
               expression_ks <= 0.22,
               identical_twas,
               module_disrupted))
}

make_replicate <- function(lib, template_key, rep_id, seed, null = FALSE) {
  universe <- make_universe(lib, template_key, seed = seed)
  rel <- make_relevant_network(lib, universe)
  deg_a <- rowSums(universe$module_a_tom > 0.1) - 1
  names(deg_a) <- universe$module_a
  set.seed(seed + 1)
  a1 <- non_extreme_seed_genes(universe$module_a_tom, universe$module_a, 5)
  a1_dir <- assign_a1_directions(a1, rep_id)
  strat <- select_path_stratified_targets(rel, universe$module_a, a1, seed = seed + 11)
  rel <- strat$adj
  degree_tab <- module_degree_table(rel, universe$module_a)
  primary_a2 <- select_degree_capped_primary(strat$A2_intermediate, degree_tab,
                                             seed = seed + 14, target_n = 30)
  c_decoy <- construct_high_degree_decoy(rel, universe, degree_tab, seed = seed + 15)
  degree_audit <- degree_audit_metrics(primary_a2, degree_tab,
                                       setdiff(universe$module_a, c(a1, primary_a2)),
                                       c_decoy)
  a2_dir <- assign_a2_direction(rel, universe$module_a, union(strat$A2_all, primary_a2),
                                a1_dir$risk, a1_dir$protective, seed = seed + 12)
  a2_risk <- names(a2_dir)[a2_dir == "risk"]
  a2_protective <- names(a2_dir)[a2_dir == "protective"]
  decoy <- construct_decoy_module(rel, universe, a1_dir$risk, a1_dir$protective, seed = seed + 13)
  twas <- generate_twas(universe$genes, a1_dir$risk, a1_dir$protective,
                        a2_risk = a2_risk, a2_protective = a2_protective,
                        decoy = decoy$genes, decoy_direction = decoy$direction,
                        seed = seed + 2, null = null)
  i1 <- make_i1_network(lib, universe, seed + 3)
  i2 <- make_i2_network(lib, universe, rel, seed + 4)
  i3 <- make_i3_network(lib, universe, rel, seed + 5)
  bridge <- setdiff(universe$module_a, c(a1, strat$A2_proximal))
  list(rep_id = rep_id, template_key = template_key, universe = universe,
       modules = list(A = universe$module_a, B = universe$module_b, C = universe$module_c,
                      background = universe$background,
                      A1 = a1, A1_risk = a1_dir$risk, A1_protective = a1_dir$protective,
                      A2 = primary_a2, A2_all = strat$A2_all,
                      A2_proximal = strat$A2_proximal,
                      A2_intermediate = strat$A2_intermediate,
                      A2_distal = strat$A2_distal,
                      A2_intermediate_distal = union(strat$A2_intermediate, strat$A2_distal),
                      A2_intermediate_degree_capped = primary_a2,
                      A2_risk = a2_risk, A2_protective = a2_protective,
                      A2_low_degree = intersect(primary_a2, degree_tab$SYMBOL[degree_tab$degree_bin == "low_degree"]),
                      A2_mid_degree = intersect(primary_a2, degree_tab$SYMBOL[degree_tab$degree_bin == "mid_degree"]),
                      A2_high_degree = intersect(primary_a2, degree_tab$SYMBOL[degree_tab$degree_bin == "high_degree"]),
                      A2_extreme_high_degree = intersect(primary_a2, degree_tab$SYMBOL[degree_tab$degree_bin == "extreme_high_degree"]),
                      A2_intermediate_degree_capped_risk = intersect(primary_a2, a2_risk),
                      A2_intermediate_degree_capped_protective = intersect(primary_a2, a2_protective),
                      A_bridge = intersect(bridge, universe$module_a),
                      D_opposite_sign_decoy = decoy$genes,
                      C_high_degree_decoy = c_decoy),
       directions = list(A2 = a2_dir, D = decoy$direction, D_branch = decoy$branch,
                         A1_majority_direction = a1_dir$majority),
       path_qc = strat$metrics,
       degree_qc = degree_audit,
       degree_table = degree_tab,
       path_lengths = strat$paths,
       thinning_log = strat$thinned_edges,
       twas = twas, networks = list(relevant = rel, I1 = i1, I2 = i2, I3 = i3))
}
