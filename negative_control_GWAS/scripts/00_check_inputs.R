#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: 00_check_inputs.R <config.yaml>")
script_path <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])
source(file.path(dirname(normalizePath(script_path)), "lib_config.R"))
cfg <- read_simple_config(args[[1]])
required_paths <- c(
  "neale_sumstats", "ldsc_munge_script", "fusion_assoc_script",
  "fusion_weights_pos", "fusion_weights_dir", "nesta_script",
  "thyroid_network_dir", "gd_twas", "result2_known_marker_rdata"
)
missing <- character()
for (key in required_paths) {
  path <- cfg_get(cfg, key)
  if (!file.exists(path)) missing <- c(missing, sprintf("%s=%s", key, path))
}
ref_prefix <- cfg_get(cfg, "fusion_ref_ld_chr")
if (!file.exists(paste0(ref_prefix, "1.bim"))) missing <- c(missing, paste0("fusion_ref_ld_chr chromosome files missing for prefix ", ref_prefix))
if (length(missing)) {
  writeLines(c("Missing required inputs:", paste0("- ", missing)))
  quit(status = 1)
}
out <- cfg_get(cfg, "output_dir")
for (d in c("input_metadata", "twas", "nesta", "rhoge", "tables", "reports")) {
  ensure_dir(file.path(out, d))
}
message("All required configured inputs are present.")
