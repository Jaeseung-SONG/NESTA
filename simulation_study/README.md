# NESTA Simulation Study

This directory is the compact reproducibility package for the final manuscript-relevant NESTA simulation study.

It preserves four scenarios:

1. Decision-rule repair and topology H robustness.
2. Comparator-framed confirmatory simulation.
3. Strict top-fraction shortlist audit.
4. False-positive control audit.

The submitted NESTA implementation remains binding at `/home/js/NESTA/Analysis/Nesta.R`. The package may add ranking, comparator, and metric wrappers, but must not alter submitted Final Heat values.

## Quick Verification

```bash
Rscript run_simulation_study.R --verify --report-dir /path/to/report
```

## Rerun Scenario Scripts

Scenario rerun wrappers are in `scripts/`. They are intentionally thin and point to the final accepted scenario implementations preserved in `reference_provenance/` and the active local implementation.

## Outputs

Generated verification reports and manuscript tables are written to the requested Dropbox result directory. The `outputs/` directory is kept empty in source control except for `README_KEEP_EMPTY.md`.
