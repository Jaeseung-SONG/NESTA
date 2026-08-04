#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: 06_run_rhoge.R <config.yaml>")
source(file.path(dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1]))), "lib_config.R"))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(RHOGE))
cfg <- read_simple_config(args[[1]])
out <- cfg_get(cfg, "output_dir")
rhoge_dir <- ensure_dir(file.path(out, "rhoge"))
report_dir <- ensure_dir(file.path(out, "reports"))
phenotype_id <- cfg_get(cfg, "phenotype_id")
phenotype_label <- cfg_get(cfg, "phenotype_label", required = FALSE, default = phenotype_id)
gd_path <- cfg_get(cfg, "gd_twas")
neg_path <- file.path(out, "twas", paste0(phenotype_id, "_Thyroid_TWAS.rds"))
gd_n <- as.numeric(cfg_get(cfg, "gd_gwas_n"))
neg_n <- as.numeric(cfg_get(cfg, "negative_control_gwas_n"))
gw_p <- as.numeric(cfg_get(cfg, "rhoge_gw_p", required = FALSE, default = "0.05"))
bd_p1 <- as.numeric(cfg_get(cfg, "rhoge_bd_p1", required = FALSE, default = "0.05"))
bd_p2 <- as.numeric(cfg_get(cfg, "rhoge_bd_p2", required = FALSE, default = "0.05"))
min_regions <- as.numeric(cfg_get(cfg, "rhoge_min_regions", required = FALSE, default = "10"))
seed <- as.integer(cfg_get(cfg, "rhoge_seed", required = FALSE, default = "20260729"))

if (!file.exists(gd_path)) stop("GD TWAS file not found: ", gd_path)
if (!file.exists(neg_path)) stop("Negative-control TWAS file not found: ", neg_path)
gd <- as.data.table(readRDS(cfg_get(cfg, "gd_twas")))
neg <- as.data.table(readRDS(neg_path))

prep <- function(x, phenotype) {
  required <- c("FILE", "ID", "CHR", "P0", "P1", "HSQ", "TWAS.Z", "TWAS.P")
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(phenotype, " TWAS is missing required RHOGE columns: ", paste(missing, collapse = ", "))
  if (!("SYMBOL" %in% names(x))) x[, SYMBOL := NA_character_]
  raw_n <- nrow(x)
  y <- copy(x)
  y[, gene_id := as.character(ID)]
  y[, gene_symbol := as.character(SYMBOL)]
  y[, TWAS.Z := as.numeric(TWAS.Z)]
  y[, TWAS.P := as.numeric(TWAS.P)]
  y[, HSQ := as.numeric(HSQ)]
  y[, CHR := as.integer(CHR)]
  y[, P0 := as.numeric(P0)]
  y[, P1 := as.numeric(P1)]
  y <- y[!is.na(gene_id) & gene_id != "" &
           is.finite(TWAS.Z) & is.finite(TWAS.P) &
           is.finite(HSQ) & HSQ > 0 &
           is.finite(CHR) & is.finite(P0) & is.finite(P1)]
  setorder(y, TWAS.P)
  y <- y[!duplicated(gene_id)]
  list(
    rhoge = y[, .(FILE, ID = gene_id, CHR, P0, P1, HSQ, TWAS.Z, TWAS.P)],
    annotated = y[, .(gene_id, gene_symbol, CHR, P0, P1, HSQ, TWAS.Z, TWAS.P, phenotype)],
    raw_n = raw_n,
    kept_n = nrow(y),
    removed_n = raw_n - nrow(y),
    bonf_sig_n = y[TWAS.P < 0.05 / nrow(y), .N],
    nominal_sig_n = y[TWAS.P < 0.05, .N]
  )
}

gd2 <- prep(gd, "GD")
neg2 <- prep(neg, phenotype_id)
overlap <- merge(gd2$annotated, neg2$annotated, by = "gene_id", suffixes = c("_GD", paste0("_", phenotype_id)))

fwrite(gd2$annotated, file.path(rhoge_dir, "GD_thyroid_TWAS_for_RHOGE.tsv"), sep = "\t")
fwrite(neg2$annotated, file.path(rhoge_dir, paste0(phenotype_id, "_thyroid_TWAS_for_RHOGE.tsv")), sep = "\t")
fwrite(overlap, file.path(rhoge_dir, paste0("GD_", phenotype_id, "_harmonized_TWAS_overlap_for_RHOGE.tsv")), sep = "\t")

set.seed(seed)
gw <- as.data.table(RHOGE::rhoge.gw(gd2$rhoge, neg2$rhoge, n1 = gd_n, n2 = neg_n, p = gw_p))
gw[, `:=`(
  analysis = "rhoge.gw",
  direction = "GD_vs_negative_control",
  trait1 = "GD",
  trait2 = phenotype_id,
  p_threshold = gw_p,
  gwas_n_trait1 = gd_n,
  gwas_n_trait2 = neg_n,
  overlapping_genes = length(intersect(gd2$rhoge$ID, neg2$rhoge$ID)),
  seed = seed
)]

set.seed(seed)
gw_rev <- as.data.table(RHOGE::rhoge.gw(neg2$rhoge, gd2$rhoge, n1 = neg_n, n2 = gd_n, p = gw_p))
gw_rev[, `:=`(
  analysis = "rhoge.gw",
  direction = "negative_control_vs_GD",
  trait1 = phenotype_id,
  trait2 = "GD",
  p_threshold = gw_p,
  gwas_n_trait1 = neg_n,
  gwas_n_trait2 = gd_n,
  overlapping_genes = length(intersect(gd2$rhoge$ID, neg2$rhoge$ID)),
  seed = seed
)]
gw_all <- rbindlist(list(gw, gw_rev), fill = TRUE)
fwrite(gw_all, file.path(rhoge_dir, paste0("RHOGE_GW_GD_vs_", phenotype_id, ".tsv")), sep = "\t")

bd_status <- "success"
bd_message <- NA_character_
bd <- tryCatch({
  set.seed(seed)
  as.data.table(RHOGE::rhoge.bd(gd2$rhoge, neg2$rhoge, n1 = gd_n, n2 = neg_n,
                                p1 = bd_p1, p2 = bd_p2, min_regions = min_regions))
}, error = function(e) {
  bd_status <<- "failed"
  bd_message <<- conditionMessage(e)
  data.table()
})
if (nrow(bd)) {
  bd[, `:=`(
    analysis = "rhoge.bd",
    trait1 = "GD",
    trait2 = phenotype_id,
    p1_threshold = bd_p1,
    p2_threshold = bd_p2,
    min_regions = min_regions,
    gwas_n_trait1 = gd_n,
    gwas_n_trait2 = neg_n,
    seed = seed
  )]
}
fwrite(bd, file.path(rhoge_dir, paste0("RHOGE_BD_GD_vs_", phenotype_id, ".tsv")), sep = "\t")

manifest <- data.table(
  field = c(
    "RHOGE_package_version", "GD_TWAS_input_path", "negative_control_TWAS_input_path",
    "gene_identifier", "effect_column", "p_value_column", "GD_GWAS_N",
    "negative_control_GWAS_N", "GD_raw_rows", "GD_kept_rows", "GD_removed_rows",
    "negative_control_raw_rows", "negative_control_kept_rows", "negative_control_removed_rows",
    "overlapping_genes", "GD_bonferroni_significant_TWAS_genes", "GD_nominal_TWAS_genes",
    "negative_control_bonferroni_significant_TWAS_genes", "negative_control_nominal_TWAS_genes",
    "rhoge_gw_p_threshold", "rhoge_bd_p1_threshold", "rhoge_bd_p2_threshold", "rhoge_bd_status",
    "rhoge_bd_message"
  ),
  value = c(
    as.character(utils::packageVersion("RHOGE")), gd_path, neg_path,
    "FUSION ID / Ensembl gene identifier", "TWAS.Z", "TWAS.P", gd_n,
    neg_n, gd2$raw_n, gd2$kept_n, gd2$removed_n,
    neg2$raw_n, neg2$kept_n, neg2$removed_n,
    length(intersect(gd2$rhoge$ID, neg2$rhoge$ID)), gd2$bonf_sig_n, gd2$nominal_sig_n,
    neg2$bonf_sig_n, neg2$nominal_sig_n,
    gw_p, bd_p1, bd_p2, bd_status, bd_message
  )
)
fwrite(manifest, file.path(rhoge_dir, "RHOGE_INPUT_MANIFEST.tsv"), sep = "\t")

rho <- gw$RHOGE[1]
se <- gw$SE[1]
pval <- gw$P[1]
non_sig <- is.finite(pval) && pval >= 0.05
close_zero <- is.finite(rho) && abs(rho) < 0.1
conclusion <- if (non_sig) {
  "The unrelated negative-control GWAS does not show significant expression-mediated genetic correlation with GD at the thyroid TWAS level."
} else {
  "The RHOGE result is statistically significant; do not describe this negative-control comparison as lacking expression-mediated genetic correlation."
}

report <- c(
  "# RHOGE GD vs 46_irnt Report",
  "",
  sprintf("RHOGE package version: `%s`.", as.character(utils::packageVersion("RHOGE"))),
  sprintf("GD TWAS input: `%s`.", gd_path),
  sprintf("Negative-control TWAS input: `%s`.", neg_path),
  "",
  "The manuscript GD thyroid TWAS input was identified from `/home/js/Thyroid_disorder/Net_anal_Grid.R`, which reads `./TWAS_res/Graves_res.rds` for the Graves/NESTA analysis, and `/home/js/Thyroid_disorder/Basal_genetic_analysis.R`, which also reads `./TWAS_res/Graves_res.rds` as `GD_TWAS`.",
  "",
  "Input harmonization used FUSION `ID` as the gene identifier, `TWAS.Z` as the effect/Z-score column, and `TWAS.P` as the p-value column.",
  sprintf("Overlapping finite gene count used for RHOGE input matching: `%s`.", length(intersect(gd2$rhoge$ID, neg2$rhoge$ID))),
  sprintf("Rows removed from GD due to missing, non-finite, non-positive HSQ, or duplicate gene IDs: `%s`.", gd2$removed_n),
  sprintf("Rows removed from 46_irnt due to missing, non-finite, non-positive HSQ, or duplicate gene IDs: `%s`.", neg2$removed_n),
  "",
  "Primary RHOGE genome-wide correlation used `RHOGE::rhoge.gw` at nominal `TWAS.P < 0.05`; the result is symmetric, so the GD->46_irnt and 46_irnt->GD input orders produced the same estimate.",
  sprintf("rho estimate: `%.8f`.", rho),
  sprintf("standard error: `%.8f`.", se),
  sprintf("p-value: `%.8g`.", pval),
  sprintf("Non-significant at alpha 0.05: `%s`.", ifelse(non_sig, "yes", "no")),
  sprintf("Estimate close to zero using |rho| < 0.1: `%s`.", ifelse(close_zero, "yes", "no")),
  "",
  sprintf("Bidirectional RHOGE sensitivity used `RHOGE::rhoge.bd` with `p1 = %.4g`, `p2 = %.4g`, and `min_regions = %s`.", bd_p1, bd_p2, min_regions),
  if (nrow(bd)) paste(capture.output(print(bd[, .(TEST, ESTIMATE, SE, TSTAT, DF, P)])), collapse = "\n") else paste0("Bidirectional sensitivity status: ", bd_status, "; ", bd_message),
  "",
  paste0("Reviewer-facing interpretation: ", conclusion)
)
writeLines(report, file.path(report_dir, "rhoge_GD_vs_46_irnt_report.md"))

message("RHOGE complete: rho=", signif(rho, 6), ", SE=", signif(se, 6), ", P=", signif(pval, 6))
