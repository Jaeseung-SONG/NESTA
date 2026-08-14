source(file.path("/home/js/NESTA/simulation", "R", "utils.R"))

cell_types <- c("B cells", "Endothelial cells", "Epithelial cells", "Fibroblasts",
                "Myeloid cells", "Proliferating cells", "SMCs or Pericytes",
                "T cells or NK cells")
tom_dir <- "/home/js/Thyroid_disorder/TOM"
coex_dir <- "/home/js/Thyroid_disorder/Coex_Net_Thyr_unsigned"

tom_file_for <- function(cell_type) file.path(tom_dir, paste0(cell_type, "_TOM.rda"))
coex_file_for <- function(cell_type) file.path(coex_dir, paste0(cell_type, ".rds"))

load_tom_dist <- function(cell_type) {
  e <- new.env(parent = emptyenv())
  load(tom_file_for(cell_type), envir = e)
  nm <- ls(e)
  if (!length(nm)) stop("No object in TOM file for ", cell_type)
  obj <- get(nm[1], envir = e)
  if (!inherits(obj, "dist")) stop("TOM object is not dist for ", cell_type)
  obj
}

tom_matrix_for <- function(cell_type, genes) {
  d <- load_tom_dist(cell_type)
  if (attr(d, "Size") != length(genes)) stop("TOM size and module gene count differ for ", cell_type)
  m <- as.matrix(d)
  rownames(m) <- genes
  colnames(m) <- genes
  diag(m) <- 1
  m
}

extract_cell_library <- function(cell_type, threshold = 0.1, max_weight_sample = 200000) {
  suppressPackageStartupMessages(library(Seurat))
  suppressPackageStartupMessages(library(hdWGCNA))
  obj <- readRDS(coex_file_for(cell_type))
  modules <- hdWGCNA::GetModules(obj)
  modules <- modules[!is.na(modules$gene_name), c("gene_name", "module", "color")]
  expr <- Matrix::rowMeans(obj@assays$SCT@data)
  expr <- expr[modules$gene_name]
  genes <- modules$gene_name
  tom <- tom_matrix_for(cell_type, genes)
  off <- tom[lower.tri(tom)]
  set.seed(abs(sum(utf8ToInt(cell_type))) + 17)
  weight_pool <- if (length(off) > max_weight_sample) sample(off, max_weight_sample) else off
  adj_bin <- tom > threshold
  diag(adj_bin) <- FALSE
  degree <- rowSums(adj_bin)
  strength <- rowSums(tom * adj_bin)
  global <- data.frame(
    cell_type = cell_type,
    n_genes = length(genes),
    density = mean(adj_bin[lower.tri(adj_bin)]),
    isolated_fraction = mean(degree == 0),
    median_degree = stats::median(degree),
    median_strength = stats::median(strength),
    expr_degree_cor = suppressWarnings(stats::cor(expr, degree, method = "spearman", use = "complete.obs")),
    expr_strength_cor = suppressWarnings(stats::cor(expr, strength, method = "spearman", use = "complete.obs")),
    stringsAsFactors = FALSE
  )
  module_ids <- setdiff(unique(modules$module), "grey")
  rows <- list()
  module_payload <- list()
  for (mod in module_ids) {
    idx <- which(modules$module == mod)
    size <- length(idx)
    if (size < 2) next
    sub <- tom[idx, idx, drop = FALSE]
    within_vals <- sub[lower.tri(sub)]
    between_vals <- tom[idx, -idx, drop = FALSE]
    g <- igraph::graph_from_adjacency_matrix((sub > threshold) * 1, mode = "undirected", diag = FALSE)
    comps <- igraph::components(g)
    largest_component_size <- if (length(comps$csize)) max(comps$csize) else 0
    clust <- suppressWarnings(mean(igraph::transitivity(g, type = "localundirected", isolates = "zero")))
    paths <- suppressWarnings(igraph::distances(g))
    finite_paths <- paths[is.finite(paths) & paths > 0]
    row <- data.frame(
      cell_type = cell_type,
      module = as.character(mod),
      size = size,
      expression_coverage = mean(is.finite(expr[idx])),
      median_within_tom = stats::median(within_vals, na.rm = TRUE),
      median_between_tom = stats::median(as.numeric(between_vals), na.rm = TRUE),
      mean_degree = mean(degree[idx]),
      median_degree = stats::median(degree[idx]),
      mean_strength = mean(strength[idx]),
      median_strength = stats::median(strength[idx]),
      clustering = ifelse(is.finite(clust), clust, 0),
      largest_component_size = largest_component_size,
      median_path_length = if (length(finite_paths)) stats::median(finite_paths) else NA_real_,
      expr_degree_cor = suppressWarnings(stats::cor(expr[idx], degree[idx], method = "spearman", use = "complete.obs")),
      expr_strength_cor = suppressWarnings(stats::cor(expr[idx], strength[idx], method = "spearman", use = "complete.obs")),
      stringsAsFactors = FALSE
    )
    key <- paste(cell_type, mod, sep = "::")
    rows[[key]] <- row
    module_payload[[key]] <- list(
      cell_type = cell_type,
      module = as.character(mod),
      genes = genes[idx],
      tom = sub,
      expression = expr[idx],
      degree = degree[idx],
      strength = strength[idx]
    )
  }
  metrics <- do.call(rbind, rows)
  if (is.null(metrics)) metrics <- data.frame()
  cell_median <- stats::median(metrics$median_within_tom, na.rm = TRUE)
  metrics$passes_template <- with(metrics,
    size >= 30 & size <= 200 & module != "grey" & median_within_tom > cell_median &
        expression_coverage >= 0.8 & (size - 5) >= 20)
  list(global = global, modules = metrics, payload = module_payload,
       expression_pool = expr[is.finite(expr)], weight_pool = weight_pool,
       all_genes = genes, degree = degree, strength = strength)
}

build_tom_topology_library <- function(out_dir = project_file("results/tom_library"),
                                       threshold = 0.1, max_templates = 40) {
  safe_dir_create(out_dir)
  cache_path <- file.path(out_dir, "tom_topology_library.rds")
  cache_csv <- file.path(out_dir, "empirical_module_templates.csv")
  if (file.exists(cache_path) && file.exists(cache_csv)) {
    return(readRDS(cache_path))
  }
  libs <- list()
  failures <- list()
  for (ct in cell_types) {
    res <- try(extract_cell_library(ct, threshold = threshold), silent = TRUE)
    if (inherits(res, "try-error")) {
      failures[[length(failures) + 1]] <- data.frame(
        cell_type = ct,
        status = "excluded",
        reason = as.character(res),
        stringsAsFactors = FALSE
      )
    } else {
      libs[[ct]] <- res
    }
  }
  if (!length(libs)) stop("No valid TOM libraries could be extracted")
  failure_tab <- if (length(failures)) do.call(rbind, failures) else data.frame(
    cell_type = character(), status = character(), reason = character())
  module_metrics <- do.call(rbind, lapply(libs, `[[`, "modules"))
  global_metrics <- do.call(rbind, lapply(libs, `[[`, "global"))
  # Relaxed binding plan accepts size up to 250 and keeps the prior
  # above-cell-type-median TOM rule; the 170626_1813 32-template library is
  # acceptable and must not be narrowed using NESTA performance.
  module_metrics$cell_type_upper75_within_tom <- ave(module_metrics$median_within_tom, module_metrics$cell_type,
                                                     FUN = function(x) stats::quantile(x, 0.75, na.rm = TRUE))
  module_metrics$passes_template <- with(module_metrics,
    size >= 30 & size <= 250 & module != "grey" &
      expression_coverage >= 0.8 & (size - 5) >= 20 &
      largest_component_size >= 30 &
      median_within_tom > ave(median_within_tom, cell_type,
                              FUN = function(x) stats::median(x, na.rm = TRUE)))
  candidates <- module_metrics[module_metrics$passes_template, , drop = FALSE]
  candidates$size_bin <- cut(candidates$size, breaks = c(29, 60, 100, 200), include.lowest = TRUE)
  candidates$tom_bin <- cut(candidates$median_within_tom,
                            breaks = unique(stats::quantile(candidates$median_within_tom,
                                                             probs = c(0, .33, .66, 1), na.rm = TRUE)),
                            include.lowest = TRUE)
  set.seed(2026061701)
  if (nrow(candidates) > max_templates) {
    candidates$stratum <- paste(candidates$cell_type, candidates$size_bin, candidates$tom_bin, sep = "|")
    keep <- unlist(tapply(seq_len(nrow(candidates)), candidates$stratum, function(i) sample(i, min(length(i), ceiling(max_templates / length(unique(candidates$stratum)))))))
    if (length(keep) > max_templates) keep <- sample(keep, max_templates)
    candidates <- candidates[sort(keep), , drop = FALSE]
  }
  keys <- paste(candidates$cell_type, candidates$module, sep = "::")
  payload <- list()
  for (key in keys) payload[[key]] <- libs[[strsplit(key, "::", fixed = TRUE)[[1]][1]]]$payload[[key]]
  lib <- list(
    threshold = threshold,
    cell_types = cell_types,
    global_metrics = global_metrics,
    module_metrics = module_metrics,
    templates = candidates,
    template_payload = payload,
    cell_payload = lapply(libs, function(x) x[c("expression_pool", "weight_pool", "all_genes", "degree", "strength")])
  )
  atomic_write_csv(global_metrics, file.path(out_dir, "tom_global_metrics.csv"))
  atomic_write_csv(module_metrics, file.path(out_dir, "tom_module_metrics.csv"))
  atomic_write_csv(candidates, file.path(out_dir, "empirical_module_templates.csv"))
  atomic_write_csv(failure_tab, file.path(out_dir, "tom_extraction_failures.csv"))
  atomic_save_rds(lib, file.path(out_dir, "tom_topology_library.rds"))
  lib
}
