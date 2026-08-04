#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: 05_gene_level_overlap_evaluation.R <config.yaml>")
source(file.path(dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1]))), "lib_config.R"))
suppressPackageStartupMessages(library(data.table))
cfg <- read_simple_config(args[[1]])
pid <- cfg_get(cfg, "phenotype_id")
out <- cfg_get(cfg, "output_dir")
files <- list.files(file.path(out, "nesta"), pattern = paste0("^", pid, "_.*_scores[.]rds$"), full.names = TRUE)
if (!length(files)) stop("No NESTA score files found for phenotype: ", pid)
scores <- rbindlist(lapply(files, function(f) as.data.table(readRDS(f))), fill = TRUE)
z_col <- if ("TWAS.Z" %in% names(scores)) "TWAS.Z" else "weight"
scores[, TWAS.Z.used := as.numeric(get(z_col))]
scores[, Final.Heat := as.numeric(Final.Heat)]
scores <- scores[is.finite(Final.Heat) & is.finite(TWAS.Z.used) & !is.na(SYMBOL) & SYMBOL != ""]
scores[, delta_NESTA := Final.Heat - TWAS.Z.used]
gene <- scores[, .(max_abs_Final_Heat = max(abs(Final.Heat)), max_abs_delta_NESTA = max(abs(delta_NESTA))), by = .(gene = SYMBOL)]
qf <- as.numeric(quantile(gene$max_abs_Final_Heat, 0.99, na.rm = TRUE))
qd <- as.numeric(quantile(gene$max_abs_delta_NESTA, 0.99, na.rm = TRUE))
selected <- gene[max_abs_Final_Heat >= qf | max_abs_delta_NESTA >= qd]
load(cfg_get(cfg, "result2_known_marker_rdata"))
ref <- as.data.table(signif.tab)
known <- sort(unique(ref[is.known.target == "Known_Marker" & !is.na(gene) & gene != "", gene]))
known_in_universe <- intersect(known, gene$gene)
overlap <- intersect(selected$gene, known)
N <- uniqueN(gene$gene); K <- length(known_in_universe); n <- uniqueN(selected$gene); k <- length(overlap)
res <- data.table(
  phenotype_id = pid,
  phenotype_label = cfg_get(cfg, "phenotype_label", required = FALSE, default = pid),
  twas_z_column_used = z_col,
  n_score_files = length(files),
  n_celltype_rows = nrow(scores),
  available_gene_universe = N,
  q99_max_abs_Final_Heat = qf,
  q99_max_abs_delta_NESTA = qd,
  selected_gene_count_q99_union = n,
  result2_known_marker_count_total = length(known),
  result2_known_marker_count_in_universe = K,
  result2_overlap_count = k,
  result2_overlap_fraction_of_selected = k / n,
  result2_overlap_fraction_of_total_known_markers = k / length(known),
  result2_overlap_fraction_of_known_markers_in_universe = k / K,
  enrichment_relative_to_available_gene_universe = (k / n) / (K / N),
  hypergeometric_p_value = phyper(k - 1, K, N - K, n, lower.tail = FALSE)
)
fwrite(res, file.path(out, "tables", paste0("negative_control_", pid, "_gene_level_q99_overlap.tsv")), sep = "\t")
