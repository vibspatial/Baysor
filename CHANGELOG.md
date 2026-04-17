# Changelog

All notable changes to the C++ line of Baysor are documented here.

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
