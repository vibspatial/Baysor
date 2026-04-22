# Changelog

All notable changes to the C++ line of Baysor are documented here.

## [cpp-0.8.1] — 2026-04-22

### Added

- `legacy` and `parquet` output styles for `baysor run`.
- A documented Xenium workflow based on `experiment.xenium` input plus `xeniumranger import-segmentation`.
- User-facing documentation pages for installation, CLI usage, inputs, outputs, preview, segfree, configuration, examples, and Xenium workflows.
- A file-level output reference for the current output bundles.
- Additional molecule clustering modes for `baysor run`:
  - `mrf`
  - `louvain`
  - `leiden`
  - `none`

### Changed

- `run`, `preview`, and `segfree` now resolve `experiment.xenium` automatically.
- `legacy` output now emits Xenium Ranger-friendly fields automatically for Xenium-origin inputs.
- NCV color generation in `run` and `preview` now uses an anchor-first streaming path:
  - learn the NCV basis from anchors
  - fit UMAP on sampled anchors
  - stream exact low-dimensional NCV vectors for all molecules
- Default `n_cells_init` is now prior-aware when a prior segmentation is present.
- `segmentation_counts.loom` writing now uses a row-oriented path instead of a late sparse-matrix duplication step.
- For graph clustering methods, `n_clusters` now controls the final coarse cluster count after anchor communities are merged.
- NCV-based clustering, 2D report UMAPs, and 3D NCV color UMAPs now share a consistent separation between:
  - the spatial neighborhood used to compute NCVs
  - the graph neighborhood used in NCV space

### Fixed

- Reduced clustering memory by removing the per-thread dense gene-by-gene correlation matrix.
- Reduced segmentation memory by switching persistent component gene counts away from dense `double` storage.
- Reduced NCV memory by removing global all-molecule neighborhood materialization and the full all-molecule high-dimensional NCV matrix from the `run` and `preview` paths.
- Reduced Loom write time and memory by avoiding an extra sparse-matrix conversion during output.
- Fixed NCV color / report embedding consistency when `cluster_method=louvain` or `leiden`, so report UMAP geometry no longer depends on an unintended clustering-only NCV neighborhood.
- Fixed streamed NCV projection to use the same logged feature transform as anchor NCV fitting, which restores non-collapsed NCV colors in `run` and `preview`.

## [0.8.0] — 2026-04-17

### Added

- Native C++ root build and test layout driven by CMake.
- 3D boundary estimation and 3D polygon output for volumetric datasets such as STARmap.
- Labeled TIFF prior support in the C++ loader.
- Adaptive NCV anchor selection that no longer fails when no molecules exceed a hard `0.95` confidence cutoff.

### Changed

- Repository structure now reflects the C++ implementation as the primary codebase.
- Example READMEs, Dockerfiles, and CI workflows now target the C++ CLI.

### Fixed

- Multiple segmentation-loop parity issues in the C++ implementation, including prior handling, drop thresholds, clustering, and history refinement behavior.
- STARmap NCV postprocessing now degrades gracefully instead of failing or producing empty anchors.
