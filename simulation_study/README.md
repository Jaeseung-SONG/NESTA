# NESTA Reviewer-Facing Simulation Study

This package reproduces the finalized NESTA simulation study outputs used for reviewer response. The public workflow is centered on two scores:

- `Final_Heat`: the final propagated NESTA heat score.
- `delta_NESTA`: the complementary re-prioritization score, defined as `Final_Heat - TWAS.Z`.

The default reviewer-facing analysis evaluates whether the top 1% rule is a conservative operating point and how relaxed thresholds change recovery and false-positive behavior. It compares four selection rules across top 0.5%, 1%, 2%, and 5% thresholds:

- Final Heat only
- delta_NESTA only
- Final OR delta_NESTA
- Final AND delta_NESTA

The threshold rule is the manuscript-style two-tailed rule:

```text
abs(Final_Heat) > q_p(abs(Final_Heat))
abs(delta_NESTA) > q_p(abs(delta_NESTA))
```

Direction-aware recovery is also reported: risk targets use the high positive tail, and protective targets use the high negative tail.

## Inputs And Provenance

The finalized deterministic simulation code and frozen configuration are vendored under:

```text
simulation_study/internal_provenance/finalized_simulation_code/
```

This directory is used only to preserve the exact finalized seed schedules, target assignment, decoy assignment, score generation, and aggregate reproducibility logic. The public runner exposes clean reviewer-facing output names and does not make the older comparator-centered framing the default workflow.

Reference reports are listed in:

```text
simulation_study/reference_manifest/
```

## Default Workflow

Run the reviewer-facing simulation workflow from the repository root:

```bash
Rscript simulation_study/run_simulation_study.R \
  --report-dir /home/js/NESTA/simulation_study_results/reviewer_delta_threshold_sensitivity_<DDMMYY_HHMM>
```

The output directory must be new. The workflow:

1. Reuses the finalized code/config/seeds unchanged.
2. Exports compact per-gene per-replicate score/label tables.
3. Recomputes the finalized aggregate reproducibility checks.
4. Stops if the finalized aggregate metrics do not match.
5. Generates the delta_NESTA-inclusive threshold sensitivity tables.
6. Verifies representative compact-table values against the validated rerun/export.

Dropbox syncing is not part of the default public workflow. To request a minimal Dropbox copy, add `--sync-dropbox` or set `NESTA_SYNC_DROPBOX=1`.

## Per-Gene Export

The per-gene export is compact and intentionally excludes large expression matrices and network adjacency matrices. It includes:

- `scenario`
- `replicate_id`
- `base_seed`, `signal_seed`, `nesta_seed`
- `gene_id`
- `Final_Heat`
- `TWAS.Z`
- `delta_NESTA`
- `expression_weighted_initialization`
- true target, risk target, protective target labels
- opposite-sign and high-score decoy labels
- selected topology/path helper labels and degree summaries

Per-gene gzip tables are written under `per_gene_tables/` in the local result directory and are not copied to Dropbox by default.

## Reviewer-Facing Outputs

The main threshold outputs are:

```text
DELTA_THRESHOLD_SENSITIVITY_REPORT.md
delta_threshold_sensitivity_summary.tsv
delta_threshold_sensitivity_by_scenario.tsv
delta_threshold_sensitivity_directional.tsv
delta_threshold_sensitivity_decoy_control.tsv
delta_threshold_sensitivity_selected_count_summary.tsv
manuscript_tables/table_delta_threshold_sensitivity.tsv
manuscript_tables/table_delta_threshold_sensitivity_compact.tsv
```

The compact manuscript table can be reproduced from any completed result directory with:

```bash
Rscript simulation_study/scripts/03_make_manuscript_tables.R <result_dir>
```

## Interpretation

The top 1% threshold is framed as a conservative operating point: it selects few genes, controls false-positive and decoy selection, and is not presented as the empirically optimal recovery threshold. Top 2% and top 5% are sensitivity analyses that show the recall/specificity tradeoff. In the validated rerun/export, relaxed thresholds increased target recovery, especially for delta_NESTA-inclusive OR selection, while opposite-sign and high-score decoy selection remained zero in the compact confirmatory table.

## Comparator-Focused Legacy Material

Older PPR/RWR-centered wrappers and public scripts are no longer part of the default workflow. They were moved to:

```text
simulation_study/internal_provenance/legacy_public_workflow/
```

They are retained only for historical provenance and are not run by the default reviewer-facing simulation command.
