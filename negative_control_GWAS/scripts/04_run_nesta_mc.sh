#!/usr/bin/env bash
set -euo pipefail
CONFIG="${1:?Usage: 04_run_nesta_mc.sh <config.yaml>}"
get_cfg() { awk -F: -v k="$1" '$1==k {sub(/^[ \t]+/,"",$2); print $2}' "$CONFIG"; }
PID="$(get_cfg phenotype_id)"
OUT="$(get_cfg output_dir)"
R_LIBS_USER="$(get_cfg r_libs_user)"
NESTA="$(get_cfg nesta_script)"
NET_DIR="$(get_cfg thyroid_network_dir)"
TWAS="$OUT/twas/${PID}_Thyroid_TWAS.rds"
mkdir -p "$OUT/nesta"
printf 'cell_type\tstatus\trc\tscore_file\n' > "$OUT/nesta/NESTA_RUN_STATUS.tsv"
for net in "$NET_DIR"/*.rds; do
  base="$(basename "$net" .rds)"
  safe="${base// /_}"
  link="$OUT/nesta/${safe}.rds"
  ln -sf "$net" "$link"
  prefix="${PID}_${safe}"
  score="$OUT/nesta/${prefix}_scores.rds"
  if Rscript -e "x<-readRDS('$score'); stopifnot(all(c('SYMBOL','Final.Heat') %in% names(x)))" >/dev/null 2>&1; then
    printf '%s\tskipped_existing_output\t0\t%s\n' "$safe" "$score" >> "$OUT/nesta/NESTA_RUN_STATUS.tsv"
    continue
  fi
  set +e
  R_LIBS_USER="$R_LIBS_USER" Rscript "$NESTA" \
    --TWAS_res "$TWAS" --Reference_net "$link" \
    --Diffuse_grid FALSE --Diffuse_method mc \
    --out_dir "$OUT/nesta/" --prefix "$prefix" --Analysis_name "$safe" \
    > "$OUT/nesta/${prefix}.stdout.log" 2> "$OUT/nesta/${prefix}.stderr.log"
  rc=$?
  set -e
  if Rscript -e "x<-readRDS('$score'); stopifnot(all(c('SYMBOL','Final.Heat') %in% names(x))); stopifnot('weight' %in% names(x) || 'TWAS.Z' %in% names(x))" >/dev/null 2>&1; then
    if [[ "$rc" -eq 0 ]]; then status=success; else status=score_written_postscore_error; fi
  else
    status=failed
  fi
  printf '%s\t%s\t%s\t%s\n' "$safe" "$status" "$rc" "$score" >> "$OUT/nesta/NESTA_RUN_STATUS.tsv"
done
