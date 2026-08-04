#!/usr/bin/env Rscript
source(file.path(dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1]))), "lib_config.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: 01_download_build_string_ppi.R <config.yaml>")
cfg <- read_simple_config(args[1])
rlibs <- cfg_get(cfg, "r_libs", required = FALSE, default = "")
if (nzchar(rlibs)) .libPaths(c(strsplit(rlibs, ":", fixed = TRUE)[[1]], .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(hdWGCNA)
})
root <- cfg_get(cfg, "output_dir")
ensure_output_tree(root)
tables <- file.path(root, "tables")
reports <- file.path(root, "reports")
cache <- file.path(root, "string_download")
source(cfg_get(cfg, "string_utils"))

gd_twas_path <- cfg_get(cfg, "gd_twas")
twas <- as.data.table(readRDS(gd_twas_path))
stopifnot(all(c("SYMBOL", "TWAS.Z", "TWAS.P") %in% names(twas)))
twas <- twas[!is.na(SYMBOL) & SYMBOL != "" & is.finite(TWAS.Z) & is.finite(TWAS.P)]
twas_genes <- unique(twas$SYMBOL)

grep_script <- function(path, pattern) {
  if (!file.exists(path)) return(character())
  lines <- readLines(path, warn = FALSE)
  hits <- grep(pattern, lines)
  sprintf("%s:%s:%s", path, hits, trimws(lines[hits]))
}
workflow_hits <- c(
  grep_script("/home/js/Thyroid_disorder/Net_anal_Grid.R", "Graves_res\\.rds|Coex_Net_Thyr|Grid_test_0715"),
  grep_script("/home/js/Thyroid_disorder/Basal_genetic_analysis.R", "Graves_res\\.rds")
)
writeLines(c(
  "# GD TWAS Input Audit",
  "",
  sprintf("Verified GD TWAS input: `%s`.", gd_twas_path),
  sprintf("Rows: `%s`.", nrow(twas)),
  sprintf("Unique TWAS gene symbols: `%s`.", length(twas_genes)),
  "",
  "Evidence from manuscript workflow scripts:",
  "",
  paste0("- `", workflow_hits, "`")
), file.path(reports, "GD_TWAS_INPUT_AUDIT.md"))

coex_files <- sort(Sys.glob(file.path(cfg_get(cfg, "coexpression_network_dir"), "*.rds")))
edge_cutoff <- cfg_num(cfg, "coexpression_edge_cutoff", 0.1)
coex <- rbindlist(lapply(coex_files, function(f) {
  cell <- tools::file_path_sans_ext(basename(f))
  obj <- readRDS(f)
  tom <- GetTOM(obj)
  diag(tom) <- 0
  keep <- tom > edge_cutoff
  deg <- rowSums(keep)
  wdeg <- rowSums(tom * keep)
  data.table(
    network = paste0("coexpression_", cell),
    network_family = "coexpression",
    cell_type = cell,
    gene_count = nrow(tom),
    edge_count = sum(upper.tri(tom) & keep),
    density = 2 * sum(upper.tri(tom) & keep) / (nrow(tom) * (nrow(tom) - 1)),
    mean_degree = mean(deg),
    median_degree = median(deg),
    mean_weighted_degree = mean(wdeg),
    median_weighted_degree = median(wdeg),
    edge_cutoff = edge_cutoff,
    twas_gene_overlap = length(intersect(rownames(tom), twas_genes)),
    zero_imputed_nodes = nrow(tom) - length(intersect(rownames(tom), twas_genes)),
    pre_subset_to_twas_genes = FALSE
  )
}), fill = TRUE)
fwrite(coex, file.path(tables, "coexpression_reference_network_topology.tsv"), sep = "\t")

built <- tryCatch({
  build_string_ppi_networks(
    cache_dir = cache,
    output_dir = cache,
    version = cfg_get(cfg, "string_version"),
    taxon_id = cfg_int(cfg, "string_taxon_id", 9606),
    thresholds = cfg_vec_num(cfg, "string_candidate_thresholds"),
    twas_genes = twas_genes,
    report_path = file.path(reports, "STRING_DOWNLOAD_AND_FULL_NETWORK_BUILD_REPORT.md")
  )
}, error = function(e) {
  writeLines(c(
    "# STRING Download Failure",
    "",
    "The full-network STRING/PPI rerun stopped because the STRING data could not be downloaded or read.",
    "",
    sprintf("Error: `%s`.", conditionMessage(e)),
    "",
    "No STRING data were fabricated."
  ), file.path(reports, "STRING_DOWNLOAD_AND_FULL_NETWORK_BUILD_REPORT.md"))
  stop(e)
})

string_all <- copy(built$summary)
string_all[, network_family := "STRING"]
target_mean_degree <- median(coex$mean_degree, na.rm = TRUE)
target_density <- median(coex$density, na.rm = TRUE)
string_all[, topology_distance := abs(log1p(mean_degree) - log1p(target_mean_degree))]
string_all[, density_distance := abs(log1p(density) - log1p(target_density))]
selected_thr <- string_all[order(topology_distance, density_distance, -threshold)][1, threshold]

default_thr <- cfg_num(cfg, "string_default_threshold", 700)
selected_modes <- data.table(
  method_threshold_mode = c("string_default", "string_coex_comparable"),
  threshold = c(default_thr, selected_thr),
  selection_rule = c("conventional STRING high-confidence threshold", "closest candidate to median thyroid co-expression mean degree; density used only as tie-breaker")
)
selected_summary <- merge(selected_modes, string_all, by = "threshold", all.x = TRUE)

fwrite(rbindlist(list(coex, string_all), fill = TRUE),
       file.path(tables, "string_threshold_topology_calibration.tsv"), sep = "\t")
fwrite(selected_summary, file.path(tables, "string_ppi_full_network_summary.tsv"), sep = "\t")

old_report <- readLines(file.path(reports, "STRING_DOWNLOAD_AND_FULL_NETWORK_BUILD_REPORT.md"), warn = FALSE)
writeLines(c(
  old_report,
  "",
  "## Threshold Definitions",
  "",
  sprintf("`string_default`: combined score >= `%s`.", default_thr),
  sprintf("`string_coex_comparable`: combined score >= `%s`.", selected_thr),
  sprintf("Topology calibration target: median co-expression mean degree `%s`; median co-expression density `%s`.",
          signif(target_mean_degree, 5), signif(target_density, 5)),
  "",
  "Candidate threshold calibration table: `tables/string_threshold_topology_calibration.tsv`.",
  "Selected full-network summary table: `tables/string_ppi_full_network_summary.tsv`."
), file.path(reports, "STRING_DOWNLOAD_AND_FULL_NETWORK_BUILD_REPORT.md"))

cat("Selected co-expression-comparable STRING threshold:", selected_thr, "\n")
