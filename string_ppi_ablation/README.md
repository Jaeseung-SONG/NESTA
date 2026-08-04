# STRING/PPI Full-network NESTA Ablation

This package reruns the reviewer-requested comparison of NESTA using submitted thyroid cell-type co-expression networks versus NESTA using a general STRING/PPI network.

The corrected analysis uses full-network diffusion. Networks are not subset to TWAS genes before diffusion. Instead, every node in the threshold-filtered diffusion graph is retained, genes with GD thyroid TWAS statistics receive an initial weight, and genes without TWAS statistics receive initial weight 0.

## Reviewer Motivation

The reviewer asked whether thyroid-specific co-expression networks add useful disease prioritization compared with a generic STRING/PPI network. The intended comparison is:

- Full thyroid co-expression network diffusion.
- Full STRING/PPI network diffusion.

The previous STRING/PPI ablation is superseded because it may have diffused on TWAS-gene-induced STRING subgraphs.

## Design

The STRING/PPI ablation is a 2 x 2 design:

- `string_default_twas_only`: STRING combined score >= 700, TWAS.Z initial weights.
- `string_default_m2_expression_weighted`: STRING combined score >= 700, M2 expression-weighted initial weights.
- `string_coex_comparable_twas_only`: topology-calibrated STRING threshold, TWAS.Z initial weights.
- `string_coex_comparable_m2_expression_weighted`: topology-calibrated STRING threshold, M2 expression-weighted initial weights.

`string_default` uses the conventional high-confidence STRING threshold, combined score >= 700.

`string_coex_comparable` is selected only from network topology. Candidate thresholds are compared against the thyroid co-expression networks using mean degree, with density used as a tie-breaker. Known GD marker recovery is not used to choose this threshold.

## Initial Weights

`twas_only`:

- Genes with TWAS statistics receive `TWAS.Z`.
- Missing genes remain in the graph with initial weight 0.

`m2_expression_weighted`:

- Genes with TWAS statistics receive the submitted M2 expression-weighted initial score.
- Missing TWAS genes remain in the graph with initial weight 0.
- In STRING/PPI mode, the same full STRING graph is reused for each thyroid cell type while expression weighting is computed from that cell type's original thyroid co-expression object.

## Co-expression Reference

The pipeline first looks for original manuscript GD co-expression NESTA score files. If available, it copies them into the corrected output root as the reference anchor. If they are missing, it reruns co-expression NESTA on the original full thyroid co-expression networks under the submitted M2 expression-weighted mode.

The co-expression audit reports both:

- The manuscript-style phenotype-specific check, including the expected approximate 45/721 known-marker overlap check.
- The gene-level collapse-first q99 check used for the STRING comparison.

## Evaluation

The primary benchmark is the Result2 known GD marker set from:

`/home/js/Thyroid_disorder/Sig_gene_per_method/Dat_with_significance.Rdata`

The evaluation is gene-level collapse-first:

1. Load all per-cell-type score files for a method.
2. Use `Final.Heat` as the final score, or `F.score` for grid-format manuscript outputs.
3. Use `TWAS.Z` if present, otherwise use `weight`, otherwise merge the GD TWAS table and zero-fill missing TWAS genes.
4. Compute `delta_NESTA = Final.Heat - TWAS.Z`.
5. Collapse across cell types by gene using maximum absolute `Final.Heat` and maximum absolute `delta_NESTA`.
6. Compute q99 thresholds on each method's gene-level distributions.
7. Select genes passing either q99 threshold.
8. Compare selected genes to Result2 known GD markers.

A size-matched sensitivity ranks genes by the maximum of normalized absolute `Final.Heat` and normalized absolute `delta_NESTA`, then compares each STRING condition at the co-expression selected-gene count.

## Run

From the repository:

```bash
bash /home/js/NESTA/string_ppi_ablation/scripts/run_string_ppi_ablation_full_network_pipeline.sh \
  /home/js/NESTA/string_ppi_ablation/config/example_config.yaml
```

The default output root is:

`/home/js/Thyroid_disorder/revision_zenodo_upload/string_ppi_ablation_full_network`

## Expected Inputs

- NESTA implementation: `/home/js/NESTA/Analysis/Nesta.R`
- STRING helper: `/home/js/NESTA/Analysis/string_ppi_utils.R`
- GD thyroid TWAS: `/home/js/Thyroid_disorder/TWAS_res/Graves_res.rds`
- Thyroid co-expression networks: `/home/js/Thyroid_disorder/Coex_Net_Thyr/`
- Result2 known marker data: `/home/js/Thyroid_disorder/Sig_gene_per_method/Dat_with_significance.Rdata`

## Expected Outputs

Key reports:

- `reports/GD_TWAS_INPUT_AUDIT.md`
- `reports/STRING_DOWNLOAD_AND_FULL_NETWORK_BUILD_REPORT.md`
- `reports/NESTA_STRING_FULL_NETWORK_CODE_AUDIT.md`
- `reports/COEXPRESSION_FULL_NETWORK_REFERENCE_AUDIT.md`
- `reports/STRING_PPI_FULL_NETWORK_ABLATION_FINAL_REPORT.md`

Key tables:

- `tables/string_ppi_full_network_summary.tsv`
- `tables/string_threshold_topology_calibration.tsv`
- `tables/string_vs_coexpression_full_network_gene_level_q99_overlap.tsv`
- `tables/string_vs_coexpression_full_network_size_matched_sensitivity.tsv`
- `tables/string_vs_coexpression_full_network_selected_gene_overlap.tsv`

## Caveats

STRING/PPI is generic and can favor high-degree or literature-rich immune genes. The Result2 marker set is a known-marker benchmark, not a causal gold standard. The analysis should therefore report whichever network performs better under the corrected full-network rule, with hub and marker-set-bias diagnostics considered during interpretation.
