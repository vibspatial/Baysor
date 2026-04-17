# Baysor

**Bay**esian **s**egmentation **o**f imaging-based spatial t**r**anscriptomics data

## Overview

Baysor segments imaging-based spatial transcriptomics data using spatial position, local gene composition, and optional prior segmentation masks. This branch contains the native C++ implementation and is built with CMake.

## Usage

The main CLI entrypoint is:

```bash
./build/baysor run --help
```

Example datasets and runnable commands live under:

- `examples/iss`
- `examples/osm-FISH`
- `examples/STARmap`

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
cmake -S . -B build
cmake --build build
```

### Run tests

```bash
ctest --test-dir build --output-on-failure
```

## Citation

If you find Baysor useful for your publication, please cite:

```
Petukhov V, Xu RJ, Soldatov RA, Cadinu P, Khodosevich K, Moffitt JR & Kharchenko PV.
Cell segmentation in imaging-based spatial transcriptomics.
Nat Biotechnol (2021). https://doi.org/10.1038/s41587-021-01044-w
```
