#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: 01_prepare_neale_sumstats.R <config.yaml>")
source(file.path(dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1]))), "lib_config.R"))
suppressPackageStartupMessages(library(data.table))
cfg <- read_simple_config(args[[1]])
pid <- cfg_get(cfg, "phenotype_id")
out <- cfg_get(cfg, "output_dir")
sumstats <- cfg_get(cfg, "neale_sumstats")
ldsc_dir <- ensure_dir(file.path(out, "twas", "ldsc"))
metadata_dir <- ensure_dir(file.path(out, "input_metadata"))
if (!file.exists(sumstats) || file.info(sumstats)$size == 0) stop("Missing Neale sumstats: ", sumstats)
x <- fread(cmd = paste("zcat", shQuote(sumstats)))
need <- c("variant", "beta", "se", "pval", "n_complete_samples")
if (!all(need %in% names(x))) stop("Missing Neale columns: ", paste(setdiff(need, names(x)), collapse = ", "))
parts <- tstrsplit(x$variant, ":", fixed = TRUE)
if (length(parts) < 4) stop("variant column is not chr:pos:ref:alt")
outdt <- data.table(
  SNP = x$variant,
  CHR = parts[[1]],
  BP = as.integer(parts[[2]]),
  A1 = parts[[4]],
  A2 = parts[[3]],
  Z = as.numeric(x$beta) / as.numeric(x$se),
  P = as.numeric(x$pval),
  N = as.numeric(x$n_complete_samples)
)
outdt <- outdt[is.finite(Z) & is.finite(P) & P > 0 & P <= 1 & is.finite(N)]
out_file <- file.path(ldsc_dir, paste0(pid, ".ldsc_fusion_input.tsv.gz"))
fwrite(outdt, out_file, sep = "\t")
writeLines(c(
  "# Allele QC",
  "",
  "Neale round2 variant parsed as `chr:pos:ref:alt`.",
  "Effect allele A1 is `alt`; non-effect allele A2 is `ref`.",
  "`minor_allele` is not used as the effect allele.",
  "`Z = beta / se`; `P = pval`; `N = n_complete_samples`."
), file.path(metadata_dir, "ALLELE_QC_REPORT.md"))
message("Wrote ", out_file)
