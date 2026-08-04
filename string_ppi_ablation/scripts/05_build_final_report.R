#!/usr/bin/env Rscript
source(file.path(dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1]))), "lib_config.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: 05_build_final_report.R <config.yaml>")
cfg <- read_simple_config(args[1])
rlibs <- cfg_get(cfg, "r_libs", required = FALSE, default = "")
if (nzchar(rlibs)) .libPaths(c(strsplit(rlibs, ":", fixed = TRUE)[[1]], .libPaths()))
suppressPackageStartupMessages(library(data.table))
root <- cfg_get(cfg, "output_dir")
tables <- file.path(root, "tables")
reports <- file.path(root, "reports")

primary <- fread(file.path(tables, "string_vs_coexpression_full_network_gene_level_q99_overlap.tsv"))
size_match <- fread(file.path(tables, "string_vs_coexpression_full_network_size_matched_sensitivity.tsv"))
sel_overlap <- fread(file.path(tables, "string_vs_coexpression_full_network_selected_gene_overlap.tsv"))
string_summary <- fread(file.path(tables, "string_ppi_full_network_summary.tsv"))
condition_manifest <- fread(file.path(tables, "string_ppi_2x2_condition_manifest.tsv"))
anchor <- fread(file.path(tables, "coexpression_full_network_reference_anchor.tsv"))
run_status <- fread(file.path(tables, "string_ppi_2x2_full_network_nesta_run_status.tsv"))

fmt_tab <- function(dt, cols) {
  paste(capture.output(print(dt[, ..cols])), collapse = "\n")
}

best <- primary[order(-result2_overlap_count, hypergeometric_p_value)][1]
co <- primary[method_id == "coexpression_reference_mc"][1]
string_best <- primary[network_type == "STRING_PPI"][order(-result2_overlap_count, hypergeometric_p_value)][1]
interp <- if (co$result2_overlap_count >= string_best$result2_overlap_count) {
  "The full-network thyroid co-expression reference recovers at least as many Result2 known GD markers as the best full-network STRING/PPI condition under the gene-level q99 rule."
} else {
  "The best full-network STRING/PPI condition recovers more Result2 known GD markers than the co-expression reference under the gene-level q99 rule; this should be reported directly and interpreted with STRING hub/density and marker-set-bias caveats."
}

writeLines(c(
  "# NESTA STRING Full-network Code Audit",
  "",
  "Inspected implementation: `/home/js/NESTA/Analysis/Nesta.R`.",
  "",
  "## `--Is_expression_network` parsing",
  "",
  "`--Is_expression_network` is parsed as a character flag and normalized to logical TRUE/FALSE, accepting TRUE/FALSE and YES/NO. This preserves the default co-expression mode while allowing reviewer-requested `--Is_expression_network NO` calls.",
  "",
  "## Full-network zero-imputation",
  "",
  "The NESTA input score frame is initialized from all vertices in the loaded reference graph. TWAS statistics are left-joined after collapsing duplicate TWAS rows by minimum `TWAS.P`; genes without an eligible TWAS statistic remain in the graph with initial weight 0.",
  "",
  "For `twas_only`, `comb.weight = TWAS.Z` for matched genes and 0 for missing genes.",
  "",
  "For `m2_expression_weighted`, mean expression is left-joined to the full graph node list. Nodes without expression or TWAS statistics are retained and assigned zero expression/TWAS contribution, so they remain zero-initialized rather than being dropped.",
  "",
  "## STRING/PPI mode",
  "",
  "`--Is_expression_network NO` sources `/home/js/NESTA/Analysis/string_ppi_utils.R` and reads full threshold-filtered STRING edge tables through `read_ppi_graph`. The corrected pipeline passes `STRING_human_v12.0_score_ge_<threshold>_edges.rds`, not the archived `*_GD_TWAS_overlap_graph.rds` files.",
  "",
  "Default manuscript behavior is unchanged when no new flags are supplied: expression-network mode defaults to TRUE and `Initial_weight_mode auto` resolves to `m2_expression_weighted` for co-expression networks."
), file.path(reports, "NESTA_STRING_FULL_NETWORK_CODE_AUDIT.md"))

writeLines(c(
  "# STRING/PPI Full-network Ablation Final Report",
  "",
  "## 1. Why the previous STRING/PPI ablation was superseded",
  "",
  "The prior STRING/PPI ablation is not reviewer-facing because its run script generated `*_GD_TWAS_overlap_graph.rds` objects by subsetting STRING to genes with GD TWAS statistics before diffusion. The reviewer-requested comparison requires full-network diffusion with zero-imputed initial weights for genes without TWAS statistics.",
  "",
  "A manifest of the archived prior output is written to `reports/PREVIOUS_STRING_ABLATION_ARCHIVE_MANIFEST.tsv`; the prior output root was moved under `archive_previous_invalid_run/`.",
  "",
  "## 2. What changed in NESTA full-network PPI mode",
  "",
  "`/home/js/NESTA/Analysis/Nesta.R` now accepts `--Is_expression_network NO` and keeps topology-only network vertices while assigning missing TWAS genes zero initial weight. The implementation details are documented in `reports/NESTA_STRING_FULL_NETWORK_CODE_AUDIT.md`.",
  "",
  "## 3. Zero-imputation implementation",
  "",
  "For every node in a loaded diffusion graph, the pipeline creates a full-length initial vector. Matched TWAS genes receive their collapsed TWAS statistic; all other graph nodes receive 0. In M2 mode, expression weighting is applied after a left join, so nodes missing expression metadata are retained with zero contribution.",
  "",
  "## 4. STRING download/build details",
  "",
  "STRING human taxonomy ID 9606 was built from protein links and protein info files. Protein IDs were mapped to preferred gene symbols, duplicate undirected gene-gene edges were collapsed by maximum combined score, and combined scores were converted to 0-1 edge weights. Full processed thresholded edge lists were saved under `string_download/` without TWAS-gene or known-marker subsetting.",
  "",
  "## 5. STRING threshold definitions",
  "",
  fmt_tab(string_summary[, .(method_threshold_mode, threshold, gene_count, edge_count, density, mean_degree, median_degree, twas_gene_overlap, zero_imputed_nodes)], names(string_summary[, .(method_threshold_mode, threshold, gene_count, edge_count, density, mean_degree, median_degree, twas_gene_overlap, zero_imputed_nodes)])),
  "",
  "## 6. 2 x 2 condition manifest",
  "",
  fmt_tab(condition_manifest, names(condition_manifest)),
  "",
  "Run status summary:",
  "",
  fmt_tab(run_status[, .N, by = .(condition, status)], c("condition", "status", "N")),
  "",
  "## 7. Co-expression full-network reference anchor",
  "",
  sprintf("Anchor source: `%s`.", anchor$anchor_source[1]),
  sprintf("Manuscript-style selected genes: `%s`; known GD marker overlap: `%s/%s`.", anchor$manuscript_style_selected_gene_count[1], anchor$manuscript_style_known_marker_overlap[1], primary$result2_known_marker_count_total[1]),
  "The 45/721 check is a manuscript/phenotype-specific threshold check and is intentionally reported separately from the collapse-first q99 reviewer comparison.",
  sprintf("Gene-level q99 co-expression selected genes: `%s`; known GD marker overlap: `%s/%s`.", anchor$gene_level_q99_selected_gene_count[1], anchor$gene_level_q99_known_marker_overlap[1], primary$result2_known_marker_count_total[1]),
  "",
  "## 8. Primary gene-level q99 comparison",
  "",
  fmt_tab(primary[, .(method_id, network_type, initial_weight_mode, threshold_mode, threshold_value, scored_gene_count, selected_gene_count, result2_overlap_count, overlap_fraction_among_selected_genes, enrichment_relative_to_available_gene_universe, hypergeometric_p_value)], c("method_id", "network_type", "initial_weight_mode", "threshold_mode", "threshold_value", "scored_gene_count", "selected_gene_count", "result2_overlap_count", "overlap_fraction_among_selected_genes", "enrichment_relative_to_available_gene_universe", "hypergeometric_p_value")),
  "",
  "## 9. Size-matched sensitivity",
  "",
  fmt_tab(size_match, names(size_match)),
  "",
  "## 10. STRING vs co-expression selected-gene overlap",
  "",
  fmt_tab(sel_overlap, names(sel_overlap)),
  "",
  "## 11. Reviewer-facing interpretation",
  "",
  interp,
  sprintf("Best method by known-marker overlap: `%s` with `%s` overlaps.", best$method_id, best$result2_overlap_count),
  "",
  "## 12. Caveats",
  "",
  "STRING/PPI is a generic literature- and evidence-integrated network and can prioritize immune/high-degree hubs for reasons unrelated to thyroid cell-type co-expression specificity. The Result2 marker set is itself a known-marker benchmark rather than an independent causal gold standard. The primary comparison therefore emphasizes full-network diffusion behavior, collapse-first gene-level selection, enrichment, and size-matched sensitivity rather than forcing a directional conclusion."
), file.path(reports, "STRING_PPI_FULL_NETWORK_ABLATION_FINAL_REPORT.md"))

cat("Wrote final full-network STRING/PPI ablation report.\n")
