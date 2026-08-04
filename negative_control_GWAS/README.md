# Negative-Control GWAS Analysis

This package reproduces the reviewer-facing negative-control GWAS/TWAS/NESTA analysis for the NESTA revision.

The final selected negative-control phenotype is UK Biobank Neale round2 `46_irnt`, hand grip strength, left. The exploratory four-phenotype screen also tested television time, computer time, and driving time, but those were screening candidates only. They are not included in the final package so the revision-facing materials focus on the selected unrelated phenotype.

## Inputs

Edit `config/example_config.yaml` with local paths to:

- Neale UKBB round2 `46_irnt` summary statistics
- LDSC conda environment and `munge_sumstats.py`
- FUSION `FUSION.assoc_test.R`
- GTEx v8 EUR Thyroid FUSION weights and 1000G EUR LD reference
- Submitted NESTA implementation, normally `Analysis/Nesta.R`
- Thyroid coexpression network RDS files
- Manuscript GD thyroid TWAS result
- Result 2 known-marker RData containing `signif.tab`
- RHOGE R package installation
- Output directory

The example config includes server paths used for the revision run, but every path can be overridden.

## Run

```bash
bash scripts/run_negative_control_GWAS_pipeline.sh config/example_config.yaml
```

The wrapper runs input checks, Neale sumstat preparation, LDSC munging, thyroid TWAS, fixed-method NESTA, gene-level overlap evaluation, RHOGE input preparation, and final report generation.

## Primary Evaluation Rule

The primary negative-control evaluation is gene-level and collapse-first:

1. Load all `46_irnt` NESTA `mc` score files.
2. Use `Final.Heat` as the final NESTA score.
3. Use `TWAS.Z`, or `weight` when NESTA score files store the TWAS Z-score there.
4. Compute `delta_NESTA = Final.Heat - TWAS.Z`.
5. Collapse across cell types by gene using `max(abs(Final.Heat))` and `max(abs(delta_NESTA))`.
6. Compute q99 thresholds on those gene-level distributions.
7. Select genes passing either gene-level q99 threshold.
8. Evaluate overlap against Result 2 known GD markers from `signif.tab`, where `is.known.target == "Known_Marker"`.

This is manuscript-consistent because the manuscript interprets reprioritized results as genes, not repeated gene-cell-type rows.

## Expected Checks

For the completed revision run, the gene-level evaluation should approximately reproduce:

- selected genes: 282
- Result 2 known GD marker reference: 721 genes
- overlap: 12 genes
- enrichment relative to available universe: about 1.006
- hypergeometric p-value: about 0.533

Material disagreement with these values should be investigated before using the results.

## RHOGE

`scripts/06_run_rhoge.R` prepares harmonized GD and `46_irnt` thyroid TWAS inputs using FUSION `ID` as the gene identifier and `TWAS.Z` / `TWAS.P` as the effect and p-value columns, then calls the installed RHOGE R package.

The primary expression-mediated genetic-correlation estimate uses `RHOGE::rhoge.gw` with nominal `TWAS.P < 0.05`. The script runs the input order in both directions; `rhoge.gw` is symmetric and should return the same estimate for GD vs `46_irnt` and `46_irnt` vs GD. It also writes a nominal-threshold `RHOGE::rhoge.bd` bidirectional sensitivity table.

For the completed revision run, RHOGE should approximately reproduce:

- overlapping finite genes: 11288
- genome-wide rho: 0.0203
- standard error: 0.1410
- p-value: 0.886

This supports the reviewer-facing conclusion only when the recomputed p-value remains non-significant.

## Outputs

Expected output folders:

- `input_metadata/`: phenotype metadata, download audit, allele QC
- `twas/`: munged inputs, chromosome TWAS logs, combined thyroid TWAS RDS
- `nesta/`: fixed `mc` NESTA score files and logs
- `tables/`: gene-level overlap and supporting result tables
- `rhoge/`: harmonized RHOGE input files and RHOGE outputs when available
- `reports/`: reviewer-facing reports and run provenance
