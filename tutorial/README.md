# NESTA Quick-Start Tutorial

This tutorial demonstrates the central NESTA workflow with a compact,
simulation-derived dataset. It is designed to run in under a minute after the
R dependencies are installed; it does **not** reproduce the full manuscript
analysis or the reviewer-facing simulation study.

## What this example shows

The same simulated TWAS results are analyzed with two different cell
type-specific co-expression networks and mean-expression profiles. NESTA
therefore starts from the same genetic evidence but produces different ranked
gene profiles for Cell type A and Cell type B.

The example contains 18 genes and includes positive and negative TWAS signals,
network-connected target genes, two weighted co-expression edge lists, and two
cell type-specific mean-expression tables. The data are synthetic and contain
no Graves' disease or participant-level data. See
[`DATA_PROVENANCE.md`](DATA_PROVENANCE.md) for its relationship to the finalized
simulation study.

## Run the tutorial

From the repository root:

```bash
bash tutorial/run_tutorial.sh
```

To use a different output directory:

```bash
bash tutorial/run_tutorial.sh /path/to/new_or_existing_output_directory
```

The script runs the public `Analysis/Nesta.R` command twice, summarizes the
cell type-specific rankings, and validates the output against the expected
top-five genes.

## Inputs

| File | Required columns | Purpose |
|---|---|---|
| `simulated_twas_results.tsv` | `SYMBOL`, `TWAS.Z`, `TWAS.P` | Shared gene-level genetic evidence |
| `cell_type_*_network.tsv` | `from`, `to`, `weight` | Weighted cell type-specific co-expression network |
| `cell_type_*_mean_expression.tsv` | `SYMBOL`, `Mean_expression` | Cell type-specific expression context |

The extra `truth_label` column in the simulated TWAS table documents how the
toy data were constructed. NESTA does not use it during prioritization.

## Command used for each cell type

```bash
Rscript Analysis/Nesta.R \
  --TWAS_res tutorial/data/simulated_twas_results.tsv \
  --Reference_net tutorial/data/cell_type_A_network.tsv \
  --Is_expression_network NO \
  --Initial_weight_mode nesta_expression_weighted \
  --Mean_expression tutorial/data/cell_type_A_mean_expression.tsv \
  --Diffuse_method raw \
  --edge_cutoff 0 \
  --Analysis_name Cell_type_A \
  --out_dir tutorial/output \
  --prefix cell_type_A
```

`--Is_expression_network NO` here means that the network is supplied as a
lightweight edge list rather than an hdWGCNA/Seurat RDS object. The network and
mean-expression table are nevertheless cell type-specific, and
`--Initial_weight_mode nesta_expression_weighted` applies the same
expression-weighted initialization concept used by NESTA.

## Outputs

The main files are:

```text
tutorial/output/
├── cell_type_A_scores.rds
├── cell_type_A_scores.tsv
├── cell_type_B_scores.rds
├── cell_type_B_scores.tsv
└── cell_type_ranking_comparison.tsv
```

Key output columns are:

| Column | Interpretation |
|---|---|
| `TWAS.Z` | Input gene-level TWAS Z-score; zero for network genes without usable TWAS evidence |
| `Mean_expression` | Cell type-specific mean expression used for initialization |
| `Initial.Heat` | Expression-weighted TWAS signal before diffusion |
| `Final.Heat` | Primary NESTA prioritization score after signed diffusion |
| `delta_NESTA` | Auxiliary score describing the change from the original TWAS Z-score |
| `Analysis_name` | Cell type or analysis label |

Rank genes primarily by the magnitude and direction of `Final.Heat`. Use
`delta_NESTA` to inspect network-driven re-prioritization rather than as a
replacement for the primary score.

## Expected result

The tutorial checks that the two cell types do not produce identical top-five
rankings. Cell type A favors the first simulated module, whereas Cell type B
favors the second. The expected ordering is stored in
`tutorial/expected/top_genes.tsv`.

Small floating-point differences can occur across R and package versions, so
the validation compares gene ordering rather than exact score decimals.

## Full reviewer-facing simulation study

The complete threshold-sensitivity and reproducibility workflow is documented
separately in [`simulation_study/README.md`](../simulation_study/README.md).
That workflow is intended to reproduce revision analyses and is substantially
more computationally demanding than this quick-start tutorial.
