#!/usr/bin/env bash
set -euo pipefail

CONFIG="${1:-/home/js/NESTA/string_ppi_ablation/config/example_config.yaml}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

Rscript "$SCRIPT_DIR/00_check_inputs.R" "$CONFIG"
Rscript "$SCRIPT_DIR/01_download_build_string_ppi.R" "$CONFIG"
Rscript "$SCRIPT_DIR/03_collect_or_run_coexpression_reference.R" "$CONFIG"
bash "$SCRIPT_DIR/02_run_nesta_string_ppi_2x2_full_network.sh" "$CONFIG"
Rscript "$SCRIPT_DIR/04_gene_level_overlap_evaluation.R" "$CONFIG"
Rscript "$SCRIPT_DIR/05_build_final_report.R" "$CONFIG"
