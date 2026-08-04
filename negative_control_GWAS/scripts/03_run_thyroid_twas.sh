#!/usr/bin/env bash
set -euo pipefail
CONFIG="${1:?Usage: 03_run_thyroid_twas.sh <config.yaml>}"
get_cfg() { awk -F: -v k="$1" '$1==k {sub(/^[ \t]+/,"",$2); print $2}' "$CONFIG"; }
PID="$(get_cfg phenotype_id)"
OUT="$(get_cfg output_dir)"
R_LIBS_USER="$(get_cfg r_libs_user)"
FUSION="$(get_cfg fusion_assoc_script)"
WEIGHTS="$(get_cfg fusion_weights_pos)"
WEIGHTS_DIR="$(get_cfg fusion_weights_dir)"
REF="$(get_cfg fusion_ref_ld_chr)"
SUMSTATS="$OUT/twas/ldsc/${PID}.sumstats.gz"
mkdir -p "$OUT/twas/chromosome_outputs"
for chr in $(seq 1 22); do
  R_LIBS_USER="$R_LIBS_USER" Rscript "$FUSION" \
    --sumstats "$SUMSTATS" \
    --weights "$WEIGHTS" \
    --weights_dir "$WEIGHTS_DIR" \
    --ref_ld_chr "$REF" \
    --chr "$chr" \
    --out "$OUT/twas/chromosome_outputs/GTExv8.EUR.Thyroid.${chr}.dat" \
    --perm 100000 \
    > "$OUT/twas/chromosome_outputs/chr${chr}.stdout.log" \
    2> "$OUT/twas/chromosome_outputs/chr${chr}.stderr.log"
done
Rscript -e "library(data.table); files <- Sys.glob(file.path('$OUT/twas/chromosome_outputs','GTExv8.EUR.Thyroid.*.dat')); x <- rbindlist(lapply(files, fread), fill=TRUE); saveRDS(as.data.frame(x), '$OUT/twas/${PID}_Thyroid_TWAS.rds')"
