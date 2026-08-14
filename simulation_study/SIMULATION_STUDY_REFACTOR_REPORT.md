# Simulation Study Refactor Report

Date: 2026-08-14

## Summary

The GitHub-facing `simulation_study/` package was refactored so the default public workflow produces the finalized reviewer-facing NESTA simulation outputs and the delta_NESTA-inclusive threshold sensitivity tables for Reviewer #2.

The public framing is now:

- Primary prioritization score: `Final_Heat`
- Complementary re-prioritization score: `delta_NESTA = Final_Heat - TWAS.Z`
- Main sensitivity rules: Final Heat only, delta_NESTA only, Final OR delta_NESTA, Final AND delta_NESTA
- Threshold grid: top 0.5%, 1%, 2%, and 5%
- Top 1% interpretation: conservative operating point, not empirical recovery optimum
- Relaxed-threshold interpretation: recovery/specificity sensitivity analysis

## Files Changed

- Added `simulation_study/README.md`.
- Replaced `simulation_study/run_simulation_study.R` with the default reviewer-facing delta threshold workflow runner.
- Rewrote `simulation_study/R/reproducibility_checks.R` around finalized aggregate checks and compact threshold-table checks.
- Updated `simulation_study/config/final_simulation_study.yaml`.
- Updated `simulation_study/reference_manifest/expected_metrics.yaml`.
- Updated `simulation_study/reference_manifest/reference_result_manifest.csv`.
- Updated lightweight public placeholder modules:
  - `simulation_study/R/generate_network.R`
  - `simulation_study/R/metrics_false_positive.R`
  - `simulation_study/R/metrics_recall_auprc.R`
- Updated `simulation_study/tests/test_reproducibility_smoke.R`.
- Updated `simulation_study/REPRODUCIBILITY_VERIFIED.md`.
- Updated `simulation_study/DELTA_THRESHOLD_SENSITIVITY_STATUS.md`.
- Added public scripts:
  - `simulation_study/scripts/01_run_delta_threshold_sensitivity.R`
  - `simulation_study/scripts/02_compare_reviewer_outputs.R`
  - `simulation_study/scripts/03_make_manuscript_tables.R`
  - `simulation_study/scripts/08_per_gene_export_delta_sensitivity.R`
- Vendored finalized reproducibility code/config under `simulation_study/internal_provenance/finalized_simulation_code/`.

## Archived Or De-Emphasized Material

The older comparator-focused public wrappers/modules were moved out of the default public workflow and into:

`simulation_study/internal_provenance/legacy_public_workflow/`

Moved files:

- `R/run_comparators.R`
- `scripts/01_run_decision_rule_repair.R`
- `scripts/02_run_confirmatory.R`
- `scripts/03_run_strict_top_fraction_audit.R`
- `scripts/04_run_false_positive_audit.R`
- `scripts/05_compare_to_reference.R`
- `scripts/06_make_manuscript_tables.R`

The old blocked post-hoc audit script was removed from the public scripts because it is superseded by the validated per-gene export and threshold sensitivity workflow.

## Terminology Cleanup

Public-facing documentation, default runner outputs, and per-gene export columns now use reviewer-facing terminology:

- `Final_Heat`
- `TWAS.Z`
- `delta_NESTA`
- `expression_weighted_initialization`
- target, signed target, and decoy labels

The old `M2` public output column was removed from the per-gene export. Old internal names remain only inside vendored finalized provenance code where changing function names would risk breaking deterministic reproduction.

## Reproducibility Check

Default workflow run:

`/home/js/NESTA/simulation_study_results/reviewer_delta_threshold_sensitivity_140826_1929`

Final aggregate reproducibility checks passed: 10/10.

All required locked metrics matched exactly within tolerance:

| Scenario | Metric | Observed | Expected |
|---|---:|---:|---:|
| decision_rule_repair | F_top100 | 0.79125 | 0.79125 |
| decision_rule_repair | F_top150 | 0.99750 | 0.99750 |
| decision_rule_repair | F_top200 | 1.00000 | 1.00000 |
| comparator_framed_confirmatory | F_top100 | 0.77950 | 0.77950 |
| comparator_framed_confirmatory | F_top150 | 0.99875 | 0.99875 |
| comparator_framed_confirmatory | F_top200 | 1.00000 | 1.00000 |
| comparator_framed_confirmatory | risk-direction recovery | 0.79625 | 0.79625 |
| comparator_framed_confirmatory | protective-direction recovery | 0.76275 | 0.76275 |
| comparator_framed_confirmatory | opposite-sign decoy selection | 0.00000 | 0.00000 |
| comparator_framed_confirmatory | high-score decoy selection | 0.00000 | 0.00000 |

Compact threshold-table checks passed against the prior validated export.

## Threshold Sensitivity Output Paths

Main local output directory:

`/home/js/NESTA/simulation_study_results/reviewer_delta_threshold_sensitivity_140826_1929`

Key files:

- `DELTA_THRESHOLD_SENSITIVITY_REPORT.md`
- `PER_GENE_EXPORT_REPRODUCIBILITY_REPORT.md`
- `reproduced_final_metric_comparison.tsv`
- `reviewer_final_metric_comparison.tsv`
- `reviewer_compact_threshold_table_comparison.tsv`
- `delta_threshold_sensitivity_summary.tsv`
- `delta_threshold_sensitivity_by_scenario.tsv`
- `delta_threshold_sensitivity_directional.tsv`
- `delta_threshold_sensitivity_decoy_control.tsv`
- `delta_threshold_sensitivity_selected_count_summary.tsv`
- `manuscript_tables/table_delta_threshold_sensitivity.tsv`
- `manuscript_tables/table_delta_threshold_sensitivity_compact.tsv`

Per-gene gzip tables were generated locally under `per_gene_tables/` and are ignored as generated results.

## Compact Table Preview

| Threshold | Rule | Selected | Recall | FPR | Risk Recall | Protective Recall | Opposite Decoy | High-score Decoy |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| top 1% | Final Heat only | 10.0000 | 0.0001875 | 0.0105184 | 0.004375 | 0.004125 | 0 | 0 |
| top 1% | Final OR delta_NESTA | 12.4000 | 0.0005000 | 0.0130316 | 0.006750 | 0.007375 | 0 | 0 |
| top 5% | Final Heat only | 50.0000 | 0.1303750 | 0.0471421 | 0.807375 | 0.782125 | 0 | 0 |
| top 5% | delta_NESTA only | 50.0000 | 0.2190625 | 0.0434079 | 0.925125 | 0.913500 | 0 | 0 |
| top 5% | Final OR delta_NESTA | 57.6275 | 0.2248125 | 0.0511947 | 0.927000 | 0.915875 | 0 | 0 |
| top 5% | Final AND delta_NESTA | 42.3725 | 0.1246250 | 0.0393553 | 0.805500 | 0.779750 | 0 | 0 |

## Verification Commands Run

- `Rscript simulation_study/tests/test_reproducibility_smoke.R`
- `Rscript -e 'invisible(parse(file="simulation_study/run_simulation_study.R")); invisible(parse(file="simulation_study/scripts/08_per_gene_export_delta_sensitivity.R")); invisible(parse(file="simulation_study/R/reproducibility_checks.R")); cat("parse_ok\n")'`
- `Rscript simulation_study/run_simulation_study.R --timestamp 140826_1929 --report-dir /home/js/NESTA/simulation_study_results/reviewer_delta_threshold_sensitivity_140826_1929`

## Caveats

- The finalized code contains old internal function and scenario names because it is vendored for exact deterministic reproduction. Those files are isolated under `internal_provenance/`.
- `comparator_framed_confirmatory` remains as a locked scenario identifier in reproducibility tables because it is the finalized scenario name used by the archived outputs.
- Full per-gene gzip tables are generated locally but should not be committed or copied to Dropbox by default.
- Dropbox syncing is optional via `--sync-dropbox`; it was not requested for the refactor verification run.
