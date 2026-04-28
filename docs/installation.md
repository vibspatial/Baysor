# Installation

The C++ branch is built with CMake. Tests are not built by default; enable them
with `-DBAYSOR_WITH_TESTS=ON` or the `tests` / `vcpkg-tests` presets.

## Required Dependencies

Baysor needs CMake, Ninja, a C++17 toolchain, plus the C++ libraries below.
Only CMake, the C++ standard, and Eigen have explicit minimums in the build.
Other libraries are intentionally not pinned so that system package managers,
Homebrew, and vcpkg can provide compatible versions.

| Dependency | Required / known-working version |
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
| GTest | Optional, only when `BAYSOR_WITH_TESTS=ON`; 1.10.0 is known to work |

Several header-only UMAP dependencies are fetched automatically by CMake with
pinned source tags:

| Header-only dependency | Pinned tag |
| --- | --- |
| `aarand` | `v1.0.2` |
| `CppKmeans` | `v3.1.1` |
| `subpar` | `v0.3.1` |
| `knncolle` | `v2.3.0` |
| `CppIrlba` | `v2.0.2` |
| `umappp` | `v2.0.1` |

## User-Space Build With vcpkg

This is the recommended cross-platform path when you do not want to install most
libraries system-wide. Dependencies are installed under the source tree or build
tree and are ignored by git.

Install vcpkg once:

```bash
git clone https://github.com/microsoft/vcpkg "$HOME/.local/src/vcpkg"
"$HOME/.local/src/vcpkg/bootstrap-vcpkg.sh" -disableMetrics
export VCPKG_ROOT="$HOME/.local/src/vcpkg"
```

On Windows, use PowerShell:

```powershell
git clone https://github.com/microsoft/vcpkg "$env:USERPROFILE\vcpkg"
& "$env:USERPROFILE\vcpkg\bootstrap-vcpkg.bat" -disableMetrics
$env:VCPKG_ROOT = "$env:USERPROFILE\vcpkg"
```

Then build and install with the same command on Linux, macOS, and Windows:

```bash
cmake -P cmake/build_and_install.cmake
```

This installs `baysor` to `./install/bin` by default. The default build is
optimized, leaves tests off, and does not require writing to system directories.

The equivalent manual commands are:

```bash
cmake --preset user-vcpkg
cmake --build --preset user-vcpkg
cmake --install build/user-vcpkg
```

On macOS, install the OpenMP runtime first:

```bash
brew install libomp
```

## System Package Builds

If dependencies are already installed in standard locations, use standard CMake:

```bash
cmake -S . -B build
cmake --build build
cmake --install build
```

Or use the helper script:

```bash
cmake -P cmake/build_and_install.cmake
```

For a custom user-owned prefix, point CMake at it:

```bash
cmake -S . -B build/user \
  -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DBAYSOR_WITH_TESTS=OFF \
  -DCMAKE_PREFIX_PATH="$HOME/.local" \
  -DCMAKE_INSTALL_PREFIX="$HOME/.local"
cmake --build build/user --target baysor
cmake --install build/user
```

## Ubuntu 24.04

The Docker build uses Ubuntu 24.04. Ubuntu's default repositories do not always
provide the Arrow / Parquet development packages in the form expected by this
build, so the commands below use the Apache Arrow apt source.

```bash
sudo apt-get update
sudo apt-get install -y --no-install-recommends \
  ca-certificates \
  lsb-release \
  wget

wget https://packages.apache.org/artifactory/arrow/$(lsb_release --id --short | tr 'A-Z' 'a-z')/apache-arrow-apt-source-latest-$(lsb_release --codename --short).deb
sudo apt-get install -y --no-install-recommends ./apache-arrow-apt-source-latest-$(lsb_release --codename --short).deb
rm ./apache-arrow-apt-source-latest-$(lsb_release --codename --short).deb

sudo apt-get update
sudo apt-get install -y --no-install-recommends \
  build-essential \
  cmake \
  ninja-build \
  pkg-config \
  git \
  libeigen3-dev \
  libomp-dev \
  libspdlog-dev \
  libcgal-dev \
  libarrow-dev \
  libparquet-dev \
  libhdf5-dev \
  nlohmann-json3-dev \
  libtiff-dev

cmake -P cmake/build_and_install.cmake
```

## macOS

```bash
brew install \
  cmake \
  ninja \
  pkg-config \
  eigen \
  libomp \
  spdlog \
  cgal \
  apache-arrow \
  hdf5 \
  nlohmann-json \
  libtiff

cmake -P cmake/build_and_install.cmake
```

## Windows

Use Visual Studio 2022 or newer and the vcpkg manifest:

```powershell
cmake -P cmake/build_and_install.cmake
```

The binary will be under:

```text
install/bin/baysor.exe
```

## Tests

System-package build with tests:

```bash
cmake --preset tests
cmake --build --preset tests
ctest --preset tests
```

User-space vcpkg build with tests:

```bash
cmake --preset vcpkg-tests
cmake --build --preset vcpkg-tests
ctest --preset vcpkg-tests
```

## Troubleshooting Dependencies

The CMake configure step checks each required dependency and prints the package
manager command to install it when it is missing. If a dependency is installed
in a non-standard location, set either:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=/path/to/prefix
```

or the package-specific config directory:

```bash
cmake -S . -B build -DArrow_DIR=/path/to/lib/cmake/arrow
```

## Continuous Integration

The `platforms_build` workflow builds the `baysor` target on Ubuntu, macOS, and
Windows. Ubuntu and macOS use native binary packages so CI does not spend time
building Apache Arrow and Thrift from source; Windows uses the vcpkg manifest.
A separate Ubuntu job enables `BAYSOR_WITH_TESTS=ON`, builds `baysor_tests`, and
runs `ctest`.
