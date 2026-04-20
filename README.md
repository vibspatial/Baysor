# Baysor

**Bay**esian **s**egmentation **o**f imaging-based spatial t**r**anscriptomics data

## Overview

Baysor segments imaging-based spatial transcriptomics data using spatial position, local gene composition, and optional prior segmentation masks.

This `cpp` branch contains the first C++ port of Baysor.

The current goal of this branch is to preserve the core segmentation algorithm
of the current Baysor release line on `master` (`v0.7.1`), while improving the
implementation around it:

- native C++17 / CMake build
- substantial performance and memory optimizations
- `legacy` and `parquet` output styles
- Parquet / GeoParquet output support
- direct `experiment.xenium` input resolution
- documented Xenium workflow via `xeniumranger import-segmentation`
- the `run`, `preview`, and `segfree` subcommands in one native binary

Future C++ releases may diverge algorithmically, but this first release is
intended as a faithful C++ implementation of the current Baysor algorithm with
a more efficient runtime and broader modern I/O support.

## Usage

The main CLI entrypoint is:

```bash
./build/baysor run --help
```

Example datasets and runnable commands:

- [Xenium pancreas](examples/Xenium_pancreas_membrane_377/README.md)
- [ISS](examples/iss/README.md)
- [osm-FISH](examples/osm-FISH/README.md)
- [STARmap](examples/STARmap/README.md)

User-facing documentation for this branch:

- [docs/README.md](docs/README.md)

## Highlights

- **Algorithmic continuity**: follows the Baysor `v0.7.1` segmentation
  algorithmic line while reimplementing it in C++.
- **Performance work**: reduces memory pressure in clustering, segmentation,
  NCV computation, and Loom writing, and improves Parquet loading.
- **Modern output support**: keeps the familiar `legacy` bundle and adds a
  `parquet` bundle with Parquet / GeoParquet tables and a 10x-style HDF5 count
  matrix.
- **Xenium workflow**: accepts `experiment.xenium` directly and documents the
  recommended Xenium Explorer handoff through `xeniumranger import-segmentation`.
- **Volumetric support**: includes 3D handling and polygon output for datasets
  such as STARmap.

## Build

### Dependencies

Install a C++17 toolchain plus the libraries required by `find_package()` in [CMakeLists.txt](CMakeLists.txt):

- Eigen3
- OpenMP
- spdlog
- CGAL
- Arrow / Parquet
- HDF5
- nlohmann_json
- libtiff
- GTest for the test target

Several header-only dependencies are fetched automatically by CMake.

### Configure and build

```bash
cmake -S . -B build -DBAYSOR_WITH_TESTS=OFF
cmake --build build --target baysor -j"$(nproc)"
```

### Run tests

```bash
cmake -S . -B build
cmake --build build --target baysor baysor_tests -j"$(nproc)"
ctest --test-dir build --output-on-failure
```

## Citation

If you find Baysor useful for your publication, please cite:

```
Petukhov V, Xu RJ, Soldatov RA, Cadinu P, Khodosevich K, Moffitt JR & Kharchenko PV.
Cell segmentation in imaging-based spatial transcriptomics.
Nat Biotechnol (2021). https://doi.org/10.1038/s41587-021-01044-w
```
