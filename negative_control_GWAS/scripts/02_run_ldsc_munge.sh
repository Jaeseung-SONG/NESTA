#!/usr/bin/env bash
set -euo pipefail
CONFIG="${1:?Usage: 02_run_ldsc_munge.sh <config.yaml>}"
DIR="$(cd "$(dirname "$0")" && pwd)"
get_cfg() { awk -F: -v k="$1" '$1==k {sub(/^[ \t]+/,"",$2); print $2}' "$CONFIG"; }
PID="$(get_cfg phenotype_id)"
OUT="$(get_cfg output_dir)"
ENV_PATH="$(get_cfg ldsc_conda_env)"
MUNGE="$(get_cfg ldsc_munge_script)"
IN="$OUT/twas/ldsc/${PID}.ldsc_fusion_input.tsv.gz"
OUT_PREFIX="$OUT/twas/ldsc/${PID}"
cmd=(conda run -p "$ENV_PATH" python "$MUNGE" --sumstats "$IN" --out "$OUT_PREFIX")
if [[ -n "${W_HM3_SNPLIST:-}" ]]; then
  cmd+=(--merge-alleles "$W_HM3_SNPLIST")
fi
"${cmd[@]}" 2> "$OUT/twas/ldsc/munge.stderr.log" | tee "$OUT/twas/ldsc/munge.stdout.log"
