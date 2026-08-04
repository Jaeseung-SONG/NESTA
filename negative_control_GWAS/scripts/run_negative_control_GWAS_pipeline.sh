#!/usr/bin/env bash
set -euo pipefail
CONFIG="${1:-$(cd "$(dirname "$0")/.." && pwd)/config/example_config.yaml}"
DIR="$(cd "$(dirname "$0")" && pwd)"
Rscript "$DIR/00_check_inputs.R" "$CONFIG"
Rscript "$DIR/01_prepare_neale_sumstats.R" "$CONFIG"
"$DIR/02_run_ldsc_munge.sh" "$CONFIG"
"$DIR/03_run_thyroid_twas.sh" "$CONFIG"
"$DIR/04_run_nesta_mc.sh" "$CONFIG"
Rscript "$DIR/05_gene_level_overlap_evaluation.R" "$CONFIG"
Rscript "$DIR/06_run_rhoge.R" "$CONFIG" || true
Rscript "$DIR/07_build_final_reports.R" "$CONFIG"
