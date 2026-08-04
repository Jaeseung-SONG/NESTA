#!/usr/bin/env bash
set -euo pipefail

CONFIG="${1:?Usage: 02_run_nesta_string_ppi_2x2_full_network.sh <config.yaml>}"
get_cfg() { awk -v k="$1" 'index($0,k":")==1{sub("^[^:]*:[ \t]*",""); print}' "$CONFIG"; }

ROOT=$(get_cfg output_dir)
NESTA=$(get_cfg nesta_script)
TWAS=$(get_cfg gd_twas)
NETDIR=$(get_cfg coexpression_network_dir)
RLIB=$(get_cfg r_libs)
STRING_VERSION=$(get_cfg string_version)
DEFAULT_THR=$(get_cfg string_default_threshold)
NPERM=$(get_cfg diffuse_nperm)
MAXJOBS=$(get_cfg max_parallel_jobs)

SUMMARY="$ROOT/tables/string_ppi_full_network_summary.tsv"
COEX_THR=$(awk -F'\t' '$2=="string_coex_comparable"{print $1; exit}' "$SUMMARY")
STATUS="$ROOT/tables/string_ppi_2x2_full_network_nesta_run_status.tsv"
MANIFEST="$ROOT/tables/string_ppi_2x2_condition_manifest.tsv"

printf 'condition\tnetwork_type\tthreshold_mode\tthreshold_value\tinitial_weight_mode\tnetwork_file\n' > "$MANIFEST"
printf 'condition\tcell_type\tscore_file\tnetwork_file\tinitial_weight_mode\tfull_network_node_count\tnonzero_twas_weighted_nodes\tzero_imputed_nodes\tedge_count\truntime_seconds\tstatus\tstdout_log\tstderr_log\n' > "$STATUS"

network_count_field() {
  local threshold="$1"; local column="$2"
  awk -F'\t' -v thr="$threshold" -v col="$column" '
    NR==1{for(i=1;i<=NF;i++) h[$i]=i; next}
    $1==thr && $2 ~ /^string_/ {print $(h[col]); exit}
  ' "$SUMMARY"
}

run_one() {
  local condition="$1"; local cell="$2"; local threshold="$3"; local mode="$4"; local expr_ref="$5"; local threshold_mode="$6"
  local network="$ROOT/string_download/STRING_human_v${STRING_VERSION}_score_ge_${threshold}_edges.rds"
  local kernel="$ROOT/string_download/STRING_human_v${STRING_VERSION}_score_ge_${threshold}_full_network_mc_kernel.rds"
  local outdir="$ROOT/nesta_string/$condition"
  local logdir="$outdir/logs"
  mkdir -p "$logdir"
  local prefix="Graves_STRING_${condition}_${cell}"
  local score="$outdir/${prefix}_scores.rds"
  local stdout="$logdir/${cell// /_}.stdout.log"
  local stderr="$logdir/${cell// /_}.stderr.log"
  local nodes edges nonzero zeros start end runtime status
  nodes=$(network_count_field "$threshold" gene_count)
  edges=$(network_count_field "$threshold" edge_count)
  nonzero=$(network_count_field "$threshold" twas_gene_overlap)
  zeros=$(network_count_field "$threshold" zero_imputed_nodes)
  start=$(date +%s)
  local args=("$NESTA"
    --TWAS_res "$TWAS"
    --Reference_net "$network"
    --Is_expression_network NO
    --Diffuse_grid FALSE
    --Diffuse_method mc
    --Diffuse_nperm "$NPERM"
    --check_bias FALSE
    --edge_cutoff 0
    --out_dir "$outdir/"
    --prefix "$prefix"
    --Analysis_name "$cell"
    --Initial_weight_mode "$mode")
  if [[ -s "$kernel" ]]; then
    args+=(--Kernel_rds "$kernel")
  else
    args+=(--Save_kernel_rds "$kernel")
  fi
  if [[ "$mode" == "m2_expression_weighted" ]]; then
    args+=(--Expression_reference_net "$expr_ref")
  fi
  if R_LIBS="$RLIB" Rscript "${args[@]}" > "$stdout" 2> "$stderr"; then
    [[ -s "$score" ]] && status=completed || status=missing_score
  else
    status=failed
  fi
  end=$(date +%s)
  runtime=$((end - start))
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$condition" "$cell" "$score" "$network" "$mode" "$nodes" "$nonzero" "$zeros" "$edges" "$runtime" "$status" "$stdout" "$stderr" >> "$STATUS"
}

register_condition() {
  local condition="$1"; local threshold="$2"; local mode="$3"; local threshold_mode="$4"
  local network="$ROOT/string_download/STRING_human_v${STRING_VERSION}_score_ge_${threshold}_edges.rds"
  printf '%s\tSTRING_PPI\t%s\t%s\t%s\t%s\n' "$condition" "$threshold_mode" "$threshold" "$mode" "$network" >> "$MANIFEST"
  for ref in "$NETDIR"/*.rds; do
    cell=$(basename "$ref" .rds)
    run_one "$condition" "$cell" "$threshold" "$mode" "$ref" "$threshold_mode" &
    while [[ $(jobs -pr | wc -l) -ge "$MAXJOBS" ]]; do sleep 2; done
  done
  wait
}

register_condition string_default_twas_only "$DEFAULT_THR" twas_only string_default
register_condition string_default_m2_expression_weighted "$DEFAULT_THR" m2_expression_weighted string_default
register_condition string_coex_comparable_twas_only "$COEX_THR" twas_only string_coex_comparable
register_condition string_coex_comparable_m2_expression_weighted "$COEX_THR" m2_expression_weighted string_coex_comparable
