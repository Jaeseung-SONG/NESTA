ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, mustWork = TRUE)
}

download_if_missing <- function(url, destfile) {
  if (file.exists(destfile) && file.info(destfile)$size > 0) {
    return(FALSE)
  }
  status <- tryCatch({
    utils::download.file(url, destfile = destfile, mode = "wb", quiet = FALSE)
  }, error = function(e) {
    if (nzchar(Sys.which("curl"))) {
      rc <- system2("curl", c("-L", "--fail", "--retry", "3", "-o", destfile, url))
      if (identical(rc, 0L) && file.exists(destfile) && file.info(destfile)$size > 0) return(0L)
    }
    attr(e, "download_url") <- url
    stop(e)
  })
  if (!identical(status, 0L) || !file.exists(destfile) || file.info(destfile)$size == 0) {
    stop("STRING download failed or produced an empty file: ", url)
  }
  TRUE
}

read_string_info <- function(info_path) {
  info <- data.table::fread(info_path)
  id_col <- intersect(c("#string_protein_id", "protein_external_id", "protein_id", "string_protein_id"), names(info))[1]
  symbol_col <- intersect(c("preferred_name", "preferred.name", "gene_name", "gene"), names(info))[1]
  if (is.na(id_col) || is.na(symbol_col)) {
    stop("Could not identify STRING protein id / gene symbol columns in: ", info_path)
  }
  out <- info[, .(protein_id = as.character(get(id_col)),
                  gene = as.character(get(symbol_col)))]
  out <- out[!is.na(protein_id) & protein_id != "" & !is.na(gene) & gene != ""]
  unique(out)
}

summarize_edge_graph <- function(edges, label, twas_genes = NULL, threshold = NA_real_) {
  if (nrow(edges) == 0) {
    return(data.table::data.table(
      network = label, threshold = threshold, gene_count = 0L, edge_count = 0L,
      density = NA_real_, mean_degree = NA_real_, median_degree = NA_real_,
      mean_weighted_degree = NA_real_, median_weighted_degree = NA_real_,
      twas_gene_overlap = 0L, zero_imputed_nodes = 0L
    ))
  }
  genes <- unique(c(edges$from, edges$to))
  idx <- match(c(edges$from, edges$to), genes)
  deg <- tabulate(idx, nbins = length(genes))
  wdeg <- tapply(rep(edges$weight, 2), idx, sum)
  wdeg <- as.numeric(wdeg[as.character(seq_along(genes))])
  wdeg[is.na(wdeg)] <- 0
  twas_n <- if (is.null(twas_genes)) NA_integer_ else length(intersect(genes, twas_genes))
  data.table::data.table(
    network = label,
    threshold = threshold,
    gene_count = length(genes),
    edge_count = nrow(edges),
    density = ifelse(length(genes) > 1, 2 * nrow(edges) / (length(genes) * (length(genes) - 1)), NA_real_),
    mean_degree = mean(deg),
    median_degree = stats::median(deg),
    q25_degree = as.numeric(stats::quantile(deg, 0.25, names = FALSE)),
    q75_degree = as.numeric(stats::quantile(deg, 0.75, names = FALSE)),
    max_degree = max(deg),
    mean_weighted_degree = mean(wdeg),
    median_weighted_degree = stats::median(wdeg),
    twas_gene_overlap = twas_n,
    zero_imputed_nodes = if (is.null(twas_genes)) NA_integer_ else length(genes) - twas_n
  )
}

build_string_ppi_networks <- function(cache_dir,
                                      output_dir = cache_dir,
                                      version = "12.0",
                                      taxon_id = 9606,
                                      thresholds = c(700),
                                      twas_genes = NULL,
                                      report_path = NULL) {
  cache_dir <- ensure_dir(cache_dir)
  output_dir <- ensure_dir(output_dir)
  links_url <- sprintf("https://stringdb-downloads.org/download/protein.links.v%s/%s.protein.links.v%s.txt.gz",
                       version, taxon_id, version)
  info_url <- sprintf("https://stringdb-downloads.org/download/protein.info.v%s/%s.protein.info.v%s.txt.gz",
                      version, taxon_id, version)
  links_path <- file.path(cache_dir, basename(links_url))
  info_path <- file.path(cache_dir, basename(info_url))
  download_log <- data.table::data.table(
    string_version = version,
    taxon_id = taxon_id,
    url = c(links_url, info_url),
    destfile = c(links_path, info_path),
    downloaded_now = c(download_if_missing(links_url, links_path),
                       download_if_missing(info_url, info_path)),
    size_bytes = file.info(c(links_path, info_path))$size,
    mtime = format(file.info(c(links_path, info_path))$mtime, "%Y-%m-%d %H:%M:%S %Z")
  )
  data.table::fwrite(download_log, file.path(output_dir, "string_download_manifest.tsv"), sep = "\t")

  protein_map <- read_string_info(info_path)
  links <- data.table::fread(links_path)
  raw_edge_count <- nrow(links)
  required <- c("protein1", "protein2", "combined_score")
  missing <- setdiff(required, names(links))
  if (length(missing)) stop("STRING links file is missing columns: ", paste(missing, collapse = ", "))

  links <- links[, .(protein1 = as.character(protein1),
                     protein2 = as.character(protein2),
                     combined_score = as.numeric(combined_score))]
  links <- merge(links, protein_map, by.x = "protein1", by.y = "protein_id", all.x = FALSE)
  data.table::setnames(links, "gene", "gene1")
  links <- merge(links, protein_map, by.x = "protein2", by.y = "protein_id", all.x = FALSE)
  data.table::setnames(links, "gene", "gene2")
  links <- links[!is.na(gene1) & !is.na(gene2) & gene1 != "" & gene2 != "" & gene1 != gene2]
  links[, from := pmin(gene1, gene2)]
  links[, to := pmax(gene1, gene2)]
  links <- links[, .(combined_score = max(combined_score, na.rm = TRUE)), by = .(from, to)]
  links[, weight := combined_score / 1000]

  thresholds <- sort(unique(as.numeric(thresholds)))
  network_files <- list()
  summaries <- list()
  for (thr in thresholds) {
    edges <- links[combined_score >= thr, .(from, to, weight, combined_score)]
    rds_path <- file.path(output_dir, sprintf("STRING_human_v%s_score_ge_%s_edges.rds", version, thr))
    saveRDS(edges, rds_path)
    network_files[[as.character(thr)]] <- rds_path
    summaries[[as.character(thr)]] <- summarize_edge_graph(
      edges,
      sprintf("STRING_score_ge_%s", thr),
      twas_genes = twas_genes,
      threshold = thr
    )[, `:=`(
      string_version = version,
      taxon_id = taxon_id,
      raw_string_edge_count = raw_edge_count,
      mapped_gene_count_before_threshold = length(unique(c(links$from, links$to))),
      mapped_edge_count_before_threshold = nrow(links),
      pre_subset_to_twas_genes = FALSE,
      pre_subset_to_known_gd_markers = FALSE
    )]
  }
  summary_tab <- data.table::rbindlist(summaries, fill = TRUE)
  data.table::fwrite(summary_tab, file.path(output_dir, "string_ppi_all_threshold_full_network_summary.tsv"), sep = "\t")

  if (!is.null(report_path)) {
    lines <- c(
      "# STRING Download And Full Network Build Report",
      "",
      sprintf("STRING version: `%s`.", version),
      sprintf("Taxonomy ID: `%s`.", taxon_id),
      sprintf("Protein links URL: `%s`.", links_url),
      sprintf("Protein info URL: `%s`.", info_url),
      sprintf("Raw STRING protein edge rows: `%s`.", raw_edge_count),
      sprintf("Mapped unique gene count before thresholding: `%s`.", length(unique(c(links$from, links$to)))),
      sprintf("Mapped unique undirected gene-gene edge count before thresholding: `%s`.", nrow(links)),
      "",
      "No TWAS-gene or known-marker subsetting is performed before saving thresholded full STRING graphs.",
      "",
      paste(capture.output(print(summary_tab)), collapse = "\n")
    )
    writeLines(lines, report_path)
  }

  list(
    download_log = download_log,
    network_files = network_files,
    summary = summary_tab,
    all_mapped_edges = links,
    links_path = links_path,
    info_path = info_path
  )
}

read_ppi_graph <- function(path, edge_cutoff = 0) {
  if (!file.exists(path)) stop("PPI network file does not exist: ", path)
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    obj <- readRDS(path)
  } else {
    obj <- data.table::fread(path)
  }
  if (inherits(obj, "igraph")) {
    g <- obj
  } else if (is.matrix(obj) || inherits(obj, "Matrix")) {
    g <- igraph::graph_from_adjacency_matrix(obj, mode = "undirected", weighted = TRUE, diag = FALSE)
  } else if (is.data.frame(obj) || data.table::is.data.table(obj)) {
    dt <- data.table::as.data.table(obj)
    if (all(c("from", "to") %in% names(dt))) {
      weight_col <- intersect(c("weight", "combined_score"), names(dt))[1]
      if (is.na(weight_col)) {
        dt[, weight := 1]
      } else if (weight_col != "weight") {
        data.table::setnames(dt, weight_col, "weight")
      }
      dt[, weight := as.numeric(weight)]
      if (max(dt$weight, na.rm = TRUE) > 1) dt[, weight := weight / 1000]
      dt <- dt[is.finite(weight) & weight > edge_cutoff & from != to]
      g <- igraph::graph_from_data_frame(dt[, .(from, to, weight)], directed = FALSE)
    } else {
      mat <- as.matrix(dt)
      storage.mode(mat) <- "numeric"
      g <- igraph::graph_from_adjacency_matrix(mat, mode = "undirected", weighted = TRUE, diag = FALSE)
    }
  } else {
    stop("Unsupported PPI network object in: ", path)
  }
  igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE,
                   edge.attr.comb = list(weight = "max", "ignore"))
}
