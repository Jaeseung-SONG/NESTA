#!/usr/bin/env Rscript
source(file.path(dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1]))), "lib_config.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: 04_gene_level_overlap_evaluation.R <config.yaml>")
cfg <- read_simple_config(args[1])
rlibs <- cfg_get(cfg, "r_libs", required = FALSE, default = "")
if (nzchar(rlibs)) .libPaths(c(strsplit(rlibs, ":", fixed = TRUE)[[1]], .libPaths()))
suppressPackageStartupMessages(library(data.table))
root <- cfg_get(cfg, "output_dir")
tables <- file.path(root, "tables")

load(cfg_get(cfg, "result2_known_marker_rdata"))
known <- unique(signif.tab[signif.tab$is.known.target == "Known_Marker", "gene", drop = TRUE])
known <- known[!is.na(known) & known != ""]
fwrite(data.table(gene = known), file.path(tables, "result2_known_GD_markers.tsv"), sep = "\t")

gd <- as.data.table(readRDS(cfg_get(cfg, "gd_twas")))
gd <- gd[!is.na(SYMBOL) & SYMBOL != "" & is.finite(TWAS.Z) & is.finite(TWAS.P)]
setorder(gd, SYMBOL, TWAS.P)
gd1 <- gd[!duplicated(SYMBOL), .(Gene = SYMBOL, TWAS.Z = TWAS.Z)]

summary_tab <- fread(file.path(tables, "string_ppi_full_network_summary.tsv"))
coex_topology <- fread(file.path(tables, "coexpression_reference_network_topology.tsv"))
get_string_row <- function(mode) summary_tab[method_threshold_mode == mode][1]
default_row <- get_string_row("string_default")
coex_row <- get_string_row("string_coex_comparable")

cond <- data.table(
  method_id = c("coexpression_reference_mc",
                "string_default_twas_only",
                "string_default_m2_expression_weighted",
                "string_coex_comparable_twas_only",
                "string_coex_comparable_m2_expression_weighted"),
  network_type = c("thyroid_cell_type_coexpression", rep("STRING_PPI", 4)),
  initial_weight_mode = c("m2_expression_weighted", "twas_only", "m2_expression_weighted", "twas_only", "m2_expression_weighted"),
  threshold_mode = c("not_applicable", "string_default", "string_default", "string_coex_comparable", "string_coex_comparable"),
  threshold_value = c(NA_real_, default_row$threshold, default_row$threshold, coex_row$threshold, coex_row$threshold),
  score_dir = c(file.path(root, "nesta_coexpression_reference"),
                file.path(root, "nesta_string/string_default_twas_only"),
                file.path(root, "nesta_string/string_default_m2_expression_weighted"),
                file.path(root, "nesta_string/string_coex_comparable_twas_only"),
                file.path(root, "nesta_string/string_coex_comparable_m2_expression_weighted"))
)
fwrite(cond, file.path(tables, "string_ppi_2x2_condition_manifest_for_evaluation.tsv"), sep = "\t")

read_method <- function(method_id, score_dir) {
  files <- sort(Sys.glob(file.path(score_dir, "*_scores.rds")))
  if (!length(files)) stop("No score files found for ", method_id, " in ", score_dir)
  rows <- rbindlist(lapply(files, function(f) {
    x <- as.data.table(readRDS(f))
    if ("method" %in% names(x)) x <- x[as.character(method) == "mc"]
    if ("Final.Heat" %in% names(x)) {
      gene_col <- if ("SYMBOL" %in% names(x)) "SYMBOL" else if ("Gene" %in% names(x)) "Gene" else "node_id"
      z_col <- if ("TWAS.Z" %in% names(x)) "TWAS.Z" else if ("weight" %in% names(x)) "weight" else NA_character_
      y <- data.table(
        Gene = as.character(x[[gene_col]]),
        Cell_type = as.character(x$Analysis_name),
        Final.Heat = as.numeric(x$Final.Heat)
      )
      if (!is.na(z_col)) y[, TWAS.Z := as.numeric(x[[z_col]])] else y <- merge(y, gd1, by = "Gene", all.x = TRUE)
    } else {
      y <- data.table(
        Gene = as.character(x$node_id),
        Cell_type = as.character(x$Analysis_name),
        Final.Heat = as.numeric(x$F.score)
      )
      y <- merge(y, gd1, by = "Gene", all.x = TRUE)
    }
    y[, source_file := f]
    y
  }), fill = TRUE)
  rows[is.na(TWAS.Z), TWAS.Z := 0]
  rows[!is.na(Gene) & Gene != "" & is.finite(Final.Heat) & is.finite(TWAS.Z)]
}

network_counts <- function(row, scored_gene_count) {
  if (row$network_type == "STRING_PPI") {
    s <- summary_tab[method_threshold_mode == row$threshold_mode][1]
    return(list(nodes = s$gene_count, edges = s$edge_count,
                nonzero = s$twas_gene_overlap, zeros = s$zero_imputed_nodes))
  }
  twas_overlap <- unique(coex_topology$twas_gene_overlap)
  list(nodes = scored_gene_count,
       edges = sum(coex_topology$edge_count, na.rm = TRUE),
       nonzero = length(intersect(unique(gd1$Gene), unique(unlist(lapply(Sys.glob(file.path(row$score_dir, "*_scores.rds")), function(f) {
         x <- as.data.table(readRDS(f))
         if ("node_id" %in% names(x)) x$node_id else if ("SYMBOL" %in% names(x)) x$SYMBOL else character()
       }))))),
       zeros = scored_gene_count - length(intersect(unique(gd1$Gene), unique(unlist(lapply(Sys.glob(file.path(row$score_dir, "*_scores.rds")), function(f) {
         x <- as.data.table(readRDS(f))
         if ("node_id" %in% names(x)) x$node_id else if ("SYMBOL" %in% names(x)) x$SYMBOL else character()
       })))))
  )
}

evaluate <- function(row) {
  rows <- read_method(row$method_id, row$score_dir)
  rows[, delta_NESTA := Final.Heat - TWAS.Z]
  gene <- rows[, .(
    max_abs_Final_Heat = max(abs(Final.Heat)),
    max_abs_delta_NESTA = max(abs(delta_NESTA)),
    n_celltype_rows = .N
  ), by = Gene]
  qfh <- as.numeric(quantile(gene$max_abs_Final_Heat, 0.99, names = FALSE))
  qd <- as.numeric(quantile(gene$max_abs_delta_NESTA, 0.99, names = FALSE))
  gene[, rank_score := pmax(max_abs_Final_Heat / qfh, max_abs_delta_NESTA / qd)]
  gene[, selected := max_abs_Final_Heat >= qfh | max_abs_delta_NESTA >= qd]
  selected <- gene[selected == TRUE, Gene]
  universe <- unique(gene$Gene)
  known_universe <- intersect(known, universe)
  overlap <- intersect(selected, known)
  N <- length(universe)
  K <- length(known_universe)
  n <- length(selected)
  x <- length(overlap)
  counts <- network_counts(row, N)
  fwrite(gene[order(-rank_score)], file.path(tables, paste0(row$method_id, "_gene_level_scores.tsv")), sep = "\t")
  fwrite(data.table(method_id = row$method_id, gene = selected), file.path(tables, paste0(row$method_id, "_selected_genes.tsv")), sep = "\t")
  data.table(
    method_id = row$method_id,
    network_type = row$network_type,
    initial_weight_mode = row$initial_weight_mode,
    threshold_mode = row$threshold_mode,
    threshold_value = row$threshold_value,
    full_network_node_count = counts$nodes,
    full_network_edge_count = counts$edges,
    nonzero_twas_weighted_node_count = counts$nonzero,
    zero_imputed_node_count = counts$zeros,
    scored_gene_count = N,
    q99_Final_Heat_threshold = qfh,
    q99_delta_NESTA_threshold = qd,
    selected_gene_count = n,
    result2_known_marker_count_total = length(known),
    result2_known_marker_count_in_available_gene_universe = K,
    result2_overlap_count = x,
    overlap_fraction_among_selected_genes = ifelse(n > 0, x / n, NA_real_),
    overlap_fraction_among_known_markers = x / length(known),
    enrichment_relative_to_available_gene_universe = ifelse(n > 0 && K > 0, (x / n) / (K / N), NA_real_),
    hypergeometric_p_value = ifelse(n > 0 && K > 0, phyper(x - 1, K, N - K, n, lower.tail = FALSE), NA_real_),
    n_score_files = length(unique(rows$source_file)),
    n_celltype_rows = nrow(rows)
  )
}

primary <- rbindlist(lapply(seq_len(nrow(cond)), function(i) evaluate(cond[i])), fill = TRUE)
fwrite(primary, file.path(tables, "string_vs_coexpression_full_network_gene_level_q99_overlap.tsv"), sep = "\t")

co_n <- primary[method_id == "coexpression_reference_mc", selected_gene_count]
co_sel <- fread(file.path(tables, "coexpression_reference_mc_selected_genes.tsv"))$gene
size_rows <- list()
overlap_rows <- list()
for (method in primary$method_id) {
  gene <- fread(file.path(tables, paste0(method, "_gene_level_scores.tsv")))
  setorder(gene, -rank_score)
  top <- head(gene$Gene, co_n)
  sel <- fread(file.path(tables, paste0(method, "_selected_genes.tsv")))$gene
  size_rows[[method]] <- data.table(
    method_id = method,
    matched_to_method = "coexpression_reference_mc",
    matched_size = co_n,
    result2_overlap_at_matched_size = length(intersect(top, known)),
    overlap_fraction_at_matched_size = length(intersect(top, known)) / co_n
  )
  overlap_rows[[method]] <- data.table(
    method_id = method,
    coexpression_selected_gene_count = length(co_sel),
    method_selected_gene_count = length(sel),
    selected_gene_overlap_with_coexpression = length(intersect(sel, co_sel)),
    overlap_fraction_of_method_selected = ifelse(length(sel) > 0, length(intersect(sel, co_sel)) / length(sel), NA_real_),
    overlap_fraction_of_coexpression_selected = ifelse(length(co_sel) > 0, length(intersect(sel, co_sel)) / length(co_sel), NA_real_)
  )
}
fwrite(rbindlist(size_rows), file.path(tables, "string_vs_coexpression_full_network_size_matched_sensitivity.tsv"), sep = "\t")
fwrite(rbindlist(overlap_rows), file.path(tables, "string_vs_coexpression_full_network_selected_gene_overlap.tsv"), sep = "\t")

print(primary)
