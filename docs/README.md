# alarawms/staphit2: Documentation

The alarawms/staphit2 documentation is split into the following pages:

- [Usage](usage.md)
  - An overview of how the pipeline works, how to run it and a description of all of the different command-line flags.
  - Includes a [clonal complex surveillance workflow](usage.md#clonal-complex-surveillance) for building curated CC97 trees from mixed local and public data.
- [Methods and validation](methods.md)
  - Why each sample is dropped or flagged, how the checks (mixed cultures, gene truncation, AMR consensus, cgMLST, clustering thresholds, BEAST2 dating) work, and how they were calibrated on a 684-sample run.
- [Output](output.md)
  - An overview of the different results produced by the pipeline and how to interpret them.
  - Includes documentation for the [iTOL annotation export](output.md#itol) (`bin/export_itol.py`).
