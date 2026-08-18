# Tutorial data provenance

The files under `tutorial/data/` are a compact, deterministic teaching fixture
derived from the design principles of the finalized NESTA simulation study.
They are not copied from the Graves' disease analysis and are not a downsample
of participant-level or expression-matrix data.

The fixture preserves the simulation features needed to explain NESTA:

- signed, gene-level TWAS evidence with positive and negative seed genes;
- weakly associated genes placed near genetically supported genes;
- weighted gene-gene networks with modular structure;
- the same TWAS evidence analyzed in two different cellular contexts; and
- cell type-specific expression profiles that emphasize different modules.

The full finalized simulation code, frozen configuration, seed schedules, and
reviewer-facing validation workflow remain under `simulation_study/`. The
tutorial fixture intentionally omits large expression matrices, adjacency
matrices, repeated simulations, decoy audits, and threshold-sensitivity runs so
that a new user can exercise the public NESTA CLI in seconds rather than hours.

The expected top-five rankings are stored in
`tutorial/expected/top_genes.tsv`. Exact floating-point scores are not frozen
because numerical values may differ slightly across R and package versions.
