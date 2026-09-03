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

Native integration documentation:

- [Native segmentation API contract](contract.md)

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
| Arrow / Parquet | `>= 19`; 19.0.1 and 20.0.0 are known to build; Arrow must include compute, CSV, and Parquet support |
| HDF5 | Not pinned; 1.10.x is known to work |
| nlohmann_json | Not pinned; 3.7.3 is known to work |
| libtiff | Not pinned; 4.1.0 is known to work |

The supported-version table is intentionally broad. The reproducible macOS
development environment below uses Eigen 5.0.1.

Several header-only dependencies are fetched automatically by CMake with pinned
tags: `aarand v1.0.2`, `CppKmeans v3.1.1`, `subpar v0.3.1`,
`knncolle v2.3.0`, `CppIrlba v2.0.2`, and `umappp v2.0.1`.

### macOS development environment with Micromamba

The following setup is intended for Apple Silicon macOS development when
Homebrew or administrator access is unavailable. It installs Micromamba and all
native dependencies under `~/Applications`; it does not write to system
directories. On an Intel Mac, replace `osx-arm64` with `osx-64` in the download
URL.

First, verify that the Apple command-line developer tools are present:

```bash
xcode-select -p
```

If that command reports that no developer tools are installed, start Apple's
user-level installer and complete the displayed installation:

```bash
xcode-select --install
```

Choose the user-space installation locations. These exports apply only to the
current shell; repeat them in a new shell or add them to the appropriate shell
configuration file.

```bash
export BAYSOR_TOOLS_DIR="${HOME}/Applications"
export BAYSOR_MICROMAMBA="${BAYSOR_TOOLS_DIR}/micromamba/bin/micromamba"
export MAMBA_ROOT_PREFIX="${BAYSOR_TOOLS_DIR}/micromamba-root"
export BAYSOR_DEV_PREFIX="${BAYSOR_TOOLS_DIR}/baysor-dev"
```

Download the official Micromamba binary:

```bash
mkdir -p "${BAYSOR_TOOLS_DIR}/micromamba"
curl -Ls https://micro.mamba.pm/api/micromamba/osx-arm64/latest \
  | tar -xj -C "${BAYSOR_TOOLS_DIR}/micromamba" bin/micromamba
"${BAYSOR_MICROMAMBA}" --version
```

Create the development environment from conda-forge:

```bash
"${BAYSOR_MICROMAMBA}" create -y \
  --prefix "${BAYSOR_DEV_PREFIX}" \
  --channel conda-forge \
  --override-channels \
  cmake \
  ninja \
  pkg-config \
  "eigen=5.0.1" \
  spdlog \
  cgal-cpp \
  "libarrow=20.*" \
  "libparquet=20.*" \
  "hdf5=1.14.*" \
  nlohmann_json \
  libtiff \
  llvm-openmp \
  gtest
```

The explicit Arrow and Parquet packages are important. Current conda-forge C++
packages are named `libarrow` and `libparquet`; do not substitute the legacy
`arrow-cpp` package, which may resolve to Arrow 13 and is incompatible with the
current Baysor source.

Activate the environment in the current Zsh session:

```bash
eval "$("${BAYSOR_MICROMAMBA}" shell hook -s zsh)"
micromamba activate "${BAYSOR_DEV_PREFIX}"
```

To make the paths and Micromamba shell integration available in future Zsh
sessions, add the following block to `~/.zshrc`:

```zsh
export BAYSOR_TOOLS_DIR="${HOME}/Applications"
export BAYSOR_MICROMAMBA="${BAYSOR_TOOLS_DIR}/micromamba/bin/micromamba"
export MAMBA_ROOT_PREFIX="${BAYSOR_TOOLS_DIR}/micromamba-root"
export BAYSOR_DEV_PREFIX="${BAYSOR_TOOLS_DIR}/baysor-dev"

if [[ -x "${BAYSOR_MICROMAMBA}" ]]; then
  eval "$("${BAYSOR_MICROMAMBA}" shell hook -s zsh)"
fi
```

The executable check prevents shell-startup errors if the Micromamba binary is
unavailable. The hook enables `micromamba activate` and `micromamba deactivate`;
it does not activate the Baysor environment automatically. Reload the file and
activate the environment when starting a Baysor development session:

```bash
source ~/.zshrc
micromamba activate "${BAYSOR_DEV_PREFIX}"
```

Do not add the activation command itself to `~/.zshrc` unless every new terminal
should start inside the Baysor development environment.

Verify the relevant tools and libraries before configuring Baysor:

```bash
cmake --version
ninja --version
micromamba list libarrow
micromamba list libparquet
echo "${CONDA_PREFIX}"
```

From the Baysor repository root, configure a clean optimized test build, compile
the CLI and unit-test targets, and run the tests:

```bash
cmake --fresh --preset tests \
  -DCMAKE_PREFIX_PATH="${CONDA_PREFIX}" \
  -DOpenMP_ROOT="${CONDA_PREFIX}"

cmake --build --preset tests --clean-first
ctest --preset tests
```

After this initial configuration and clean build, the regular development loop
is simply:

```bash
cmake --build --preset tests
ctest --preset tests
```

The build is incremental, so only targets affected by source changes are
recompiled. Re-run the `cmake --fresh --preset tests` configuration when native
dependencies, the Micromamba environment, CMake options, or material build
configuration have changed, or when intentionally discarding a stale CMake
cache.

During environment creation, Micromamba warns that conda packages may contain
installation scripts. This is a security notice rather than an installation
failure; continue only if the configured package source (`conda-forge` above) is
trusted.

If configuration reports Parquet 13, or compilation fails around
`GetRecordBatchReader` or `AddKeyValueMetadata`, the environment contains an
incompatible Arrow installation. Upgrade an existing environment explicitly:

```bash
micromamba install \
  --prefix "${BAYSOR_DEV_PREFIX}" \
  --channel conda-forge \
  --override-channels \
  "libarrow=20.*" \
  "libparquet=20.*"
```

Then repeat the `cmake --fresh` command. If CTest says that `baysor_tests`
cannot be found, inspect the preceding build output: that message means
compilation failed before the test executable was created.

### Configure, build, and install

After dependencies are installed, use the same command on Linux, macOS, and
Windows:

```bash
cmake -P cmake/build_and_install.cmake
```

This configures an end-user build: optimized, tests off, and installed to
`./install/bin`. Platform-specific prerequisite commands are in
[docs/installation.md](docs/installation.md). Windows uses vcpkg when
`VCPKG_ROOT` is set; Linux and macOS use system packages by default.

Run the installed binary with:

```bash
./install/bin/baysor --help
```

Detailed installation instructions are in [docs/installation.md](docs/installation.md).

## Citation

If you find Baysor useful for your publication, please cite:

```
Petukhov V, Xu RJ, Soldatov RA, Cadinu P, Khodosevich K, Moffitt JR & Kharchenko PV.
Cell segmentation in imaging-based spatial transcriptomics.
Nat Biotechnol (2021). https://doi.org/10.1038/s41587-021-01044-w
```
