# Baysor C++ Branch

Baysor segments imaging-based spatial transcriptomics data using spatial
position, local gene composition, and optional prior segmentation inputs.

This `cpp` branch is the native C++ implementation of the main Baysor workflow.
It is built with CMake and currently exposes three CLI subcommands:

- `run`
- `preview`
- `segfree`

## What Is Implemented

- `baysor run`
- `baysor preview`
- `baysor segfree`
- CSV / Parquet transcript-table input
- image-mask and polygon prior segmentation
- Xenium-aware input resolution through `experiment.xenium`
- `legacy` and `parquet` output styles
- HTML plotting / diagnostic outputs from `run`

## Start Here

- [Installation](installation.md)
- [How Baysor Segmentation Works](algorithm.md)
- [Running Baysor](run.md)
- [Dataset Preview](preview.md)
- [Segmentation-Free NCVs](segfree.md)
- [Input Data](inputs.md)
- [Outputs](outputs.md)
- [Output Files](output_files.md)
- [Xenium Workflow](xenium.md)
- [Examples](examples.md)
