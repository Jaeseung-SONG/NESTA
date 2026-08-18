#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
OUTPUT_DIR="${1:-${SCRIPT_DIR}/output}"

mkdir -p "${OUTPUT_DIR}"

for CELL_TYPE in A B; do
  Rscript "${REPO_ROOT}/Analysis/Nesta.R" \
    --TWAS_res "${SCRIPT_DIR}/data/simulated_twas_results.tsv" \
    --Reference_net "${SCRIPT_DIR}/data/cell_type_${CELL_TYPE}_network.tsv" \
    --Is_expression_network NO \
    --Initial_weight_mode nesta_expression_weighted \
    --Mean_expression "${SCRIPT_DIR}/data/cell_type_${CELL_TYPE}_mean_expression.tsv" \
    --Diffuse_method raw \
    --Diffuse_nperm 50 \
    --check_bias FALSE \
    --edge_cutoff 0 \
    --Analysis_name "Cell_type_${CELL_TYPE}" \
    --out_dir "${OUTPUT_DIR}" \
    --prefix "cell_type_${CELL_TYPE}"
done

Rscript "${SCRIPT_DIR}/summarize_results.R" "${OUTPUT_DIR}"
Rscript "${SCRIPT_DIR}/validate_results.R" "${OUTPUT_DIR}"

echo "Tutorial completed: ${OUTPUT_DIR}"
