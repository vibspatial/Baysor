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

Install CMake, Ninja, a C++17 toolchain, plus the libraries required by
`find_package()` in [CMakeLists.txt](CMakeLists.txt). Versions are kept
intentionally broad for package-manager builds:

| Dependency | Version note |
| --- | --- |
| CMake | `>= 3.20` |
| C++ compiler | C++17 compiler; GCC 9.4.0 and Visual Studio 2022 are known to work |
| Ninja | Recent Ninja; 1.10.0 is known to work |
| Eigen3 | `>= 3.3` |
| OpenMP | C++ OpenMP target; GCC OpenMP 4.5 is known to work |
| spdlog | Not pinned; 1.5.0 is known to work |
| CGAL | Not pinned; 5.0.2 is known to work |
| Arrow / Parquet | Not pinned; 19.0.1 is known to work; Arrow must include compute, CSV, and Parquet support |
| HDF5 | Not pinned; 1.10.x is known to work |
| nlohmann_json | Not pinned; 3.7.3 is known to work |
| libtiff | Not pinned; 4.1.0 is known to work |
| GTest | Optional tests only; 1.10.0 is known to work |

Several header-only dependencies are fetched automatically by CMake with pinned
tags: `aarand v1.0.2`, `CppKmeans v3.1.1`, `subpar v0.3.1`,
`knncolle v2.3.0`, `CppIrlba v2.0.2`, and `umappp v2.0.1`.

### Configure, build, and install

After dependencies are installed, the standard CMake commands are:

```bash
cmake -S . -B build
cmake --build build
cmake --install build
```

This configures an end-user build: optimized, tests off, and installed to
`./install/bin`.

The same system-package build can be run as one command:

```bash
cmake -P cmake/build_and_install.cmake
```

For a user-space dependency install through vcpkg, set `VCPKG_ROOT` first:

```bash
export VCPKG_ROOT="$HOME/.local/src/vcpkg"
cmake -P cmake/build_and_install.cmake
```

If `VCPKG_ROOT` is set, the helper uses vcpkg. Otherwise it uses system
packages visible to CMake.

### Run tests

Tests are not built by default. To enable them:

```bash
cmake --preset tests
cmake --build --preset tests
ctest --preset tests
```

Detailed platform-specific installation instructions are in [docs/installation.md](docs/installation.md).

## Citation

If you find Baysor useful for your publication, please cite:

```
Petukhov V, Xu RJ, Soldatov RA, Cadinu P, Khodosevich K, Moffitt JR & Kharchenko PV.
Cell segmentation in imaging-based spatial transcriptomics.
Nat Biotechnol (2021). https://doi.org/10.1038/s41587-021-01044-w
```
