#!/usr/bin/env Rscript
source(file.path(dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1]))), "lib_config.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: 03_collect_or_run_coexpression_reference.R <config.yaml>")
cfg <- read_simple_config(args[1])
rlibs <- cfg_get(cfg, "r_libs", required = FALSE, default = "")
if (nzchar(rlibs)) .libPaths(c(strsplit(rlibs, ":", fixed = TRUE)[[1]], .libPaths()))
suppressPackageStartupMessages(library(data.table))
root <- cfg_get(cfg, "output_dir")
out <- file.path(root, "nesta_coexpression_reference")
reports <- file.path(root, "reports")
tables <- file.path(root, "tables")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

reference_dir <- cfg_get(cfg, "coexpression_reference_score_dir")
files <- Sys.glob(file.path(reference_dir, "Graves_res_*_scores.rds"))
anchor_source <- "copied_existing_manuscript_grid_scores"
if (length(files)) {
  ok <- file.copy(files, out, overwrite = TRUE)
  score_files <- file.path(out, basename(files))[ok | file.exists(file.path(out, basename(files)))]
} else {
  anchor_source <- "rerun_full_coexpression_mc"
  score_files <- character()
  refs <- sort(Sys.glob(file.path(cfg_get(cfg, "coexpression_network_dir"), "*.rds")))
  for (ref in refs) {
    cell <- tools::file_path_sans_ext(basename(ref))
    prefix <- paste0("Graves_res_", cell)
    cmd <- c(
      cfg_get(cfg, "nesta_script"),
      "--TWAS_res", cfg_get(cfg, "gd_twas"),
      "--Reference_net", ref,
      "--Is_expression_network", "YES",
      "--Diffuse_grid", "FALSE",
      "--Diffuse_method", cfg_get(cfg, "diffuse_method", required = FALSE, default = "mc"),
      "--Diffuse_nperm", cfg_get(cfg, "diffuse_nperm", required = FALSE, default = "300"),
      "--check_bias", "FALSE",
      "--out_dir", paste0(out, "/"),
      "--prefix", prefix,
      "--Analysis_name", cell,
      "--Initial_weight_mode", "nesta_expression_weighted"
    )
    log <- file.path(out, paste0(safe_slug(cell), ".coexpression_reference.log"))
    rc <- system2("Rscript", cmd, stdout = log, stderr = log)
    if (rc != 0) stop("Co-expression NESTA rerun failed for: ", cell)
    score_files <- c(score_files, file.path(out, paste0(prefix, "_scores.rds")))
  }
}
if (!length(score_files)) stop("No co-expression reference score files available.")

load(cfg_get(cfg, "result2_known_marker_rdata"))
known <- unique(signif.tab[signif.tab$is.known.target == "Known_Marker", "gene", drop = TRUE])
known <- known[!is.na(known) & known != ""]

man <- as.data.table(signif.tab)
man <- man[Phenotype == "Graves" & as.character(method) == "mc"]
man_selected <- unique(man[Diffuse.signif == "Sig" | Delta.signif == "Sig", gene])
man_overlap <- length(intersect(man_selected, known))

gd <- as.data.table(readRDS(cfg_get(cfg, "gd_twas")))
gd <- gd[!is.na(SYMBOL) & SYMBOL != "" & is.finite(TWAS.Z) & is.finite(TWAS.P)]
setorder(gd, SYMBOL, TWAS.P)
gd1 <- gd[!duplicated(SYMBOL), .(Gene = SYMBOL, TWAS.Z = TWAS.Z)]

rows <- rbindlist(lapply(score_files, function(f) {
  x <- as.data.table(readRDS(f))
  if ("method" %in% names(x)) x <- x[as.character(method) == "mc"]
  if ("Final.Heat" %in% names(x)) {
    gene_col <- if ("SYMBOL" %in% names(x)) "SYMBOL" else if ("Gene" %in% names(x)) "Gene" else "node_id"
    z_col <- if ("TWAS.Z" %in% names(x)) "TWAS.Z" else if ("weight" %in% names(x)) "weight" else NA_character_
    y <- data.table(Gene = as.character(x[[gene_col]]), Cell_type = as.character(x$Analysis_name),
                    Final.Heat = as.numeric(x$Final.Heat))
    if (!is.na(z_col)) y[, TWAS.Z := as.numeric(x[[z_col]])] else y <- merge(y, gd1, by = "Gene", all.x = TRUE)
  } else {
    y <- data.table(Gene = as.character(x$node_id), Cell_type = as.character(x$Analysis_name),
                    Final.Heat = as.numeric(x$F.score))
    y <- merge(y, gd1, by = "Gene", all.x = TRUE)
  }
  y
}), fill = TRUE)
rows[is.na(TWAS.Z), TWAS.Z := 0]
rows <- rows[!is.na(Gene) & Gene != "" & is.finite(Final.Heat)]
rows[, delta_NESTA := Final.Heat - TWAS.Z]
gene <- rows[, .(max_abs_Final_Heat = max(abs(Final.Heat)),
                 max_abs_delta_NESTA = max(abs(delta_NESTA))), by = Gene]
qfh <- as.numeric(quantile(gene$max_abs_Final_Heat, 0.99, names = FALSE))
qd <- as.numeric(quantile(gene$max_abs_delta_NESTA, 0.99, names = FALSE))
gene[, selected := max_abs_Final_Heat >= qfh | max_abs_delta_NESTA >= qd]
q99_selected <- gene[selected == TRUE, Gene]
q99_overlap <- length(intersect(q99_selected, known))
N <- uniqueN(gene$Gene)
K <- length(intersect(known, gene$Gene))
n <- length(q99_selected)
x <- q99_overlap
q99_p <- phyper(x - 1, K, N - K, n, lower.tail = FALSE)

fwrite(data.table(
  anchor_source = anchor_source,
  score_file_count = length(score_files),
  manuscript_style_selected_gene_count = length(man_selected),
  manuscript_style_known_marker_overlap = man_overlap,
  expected_manuscript_overlap_note = "Prior reviewer-work expectation was approximately 45/721 for the manuscript/phenotype-specific rule.",
  gene_level_q99_selected_gene_count = n,
  gene_level_q99_known_marker_overlap = x,
  gene_level_q99_enrichment = (x / n) / (K / N),
  gene_level_q99_hypergeometric_p = q99_p
), file.path(tables, "coexpression_full_network_reference_anchor.tsv"), sep = "\t")

warning_line <- if (man_overlap == 45) {
  "The manuscript-style check exactly matches the expected 45 known-marker overlap."
} else {
  sprintf("Warning: the manuscript-style check observed `%s` known-marker overlap, not exactly 45; interpret STRING comparisons using the explicit gene-level q99 rule.", man_overlap)
}

writeLines(c(
  "# Co-expression Full Network Reference Audit",
  "",
  sprintf("Anchor source: `%s`.", anchor_source),
  sprintf("Score files used: `%s`.", length(score_files)),
  sprintf("Reference score directory: `%s`.", reference_dir),
  sprintf("Copied/rerun score destination: `%s`.", out),
  "",
  "## Manuscript-style Co-expression Check",
  "",
  sprintf("Selected gene count under `Phenotype == Graves`, `method == mc`, and `Diffuse.signif == Sig | Delta.signif == Sig`: `%s`.", length(man_selected)),
  sprintf("Result2 known GD marker overlap: `%s/%s`.", man_overlap, length(known)),
  warning_line,
  "",
  "## Gene-level Collapse-first q99 Co-expression Check",
  "",
  sprintf("Selected gene count: `%s`.", n),
  sprintf("Result2 known GD marker overlap: `%s/%s`.", x, length(known)),
  sprintf("Hypergeometric p-value: `%s`.", signif(q99_p, 5)),
  "",
  "This q99 check is the anchor used for the STRING/PPI comparison tables."
), file.path(reports, "COEXPRESSION_FULL_NETWORK_REFERENCE_AUDIT.md"))

cat("Prepared co-expression reference:", length(score_files), "score files.\n")
