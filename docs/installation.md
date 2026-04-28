# Installation

The C++ branch is built with CMake. Tests are not built by default; enable them
with `-DBAYSOR_WITH_TESTS=ON` or the `tests` / `vcpkg-tests` presets.

## Required Dependencies

Baysor needs CMake, Ninja, a C++17 toolchain, plus:

- Eigen3
- OpenMP
- spdlog
- CGAL
- Arrow / Parquet
- HDF5
- nlohmann_json
- libtiff
- GTest only when `BAYSOR_WITH_TESTS=ON`

Several header-only UMAP dependencies are fetched automatically by CMake.

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

Configure and build:

```bash
cmake --preset vcpkg-release
cmake --build --preset vcpkg-release
```

Install the binary to a user-owned prefix:

```bash
cmake --install build/vcpkg-release --prefix "$HOME/.local"
```

On macOS, install the OpenMP runtime first:

```bash
brew install libomp
```

## System Package Builds

If dependencies are already installed in standard locations:

```bash
cmake --preset release
cmake --build --preset release
```

For a custom user-owned prefix, point CMake at it:

```bash
cmake -S . -B build/release \
  -DCMAKE_BUILD_TYPE=Release \
  -DBAYSOR_WITH_TESTS=OFF \
  -DCMAKE_PREFIX_PATH="$HOME/.local" \
  -DCMAKE_INSTALL_PREFIX="$HOME/.local"
cmake --build build/release --target baysor -j"$(nproc)"
cmake --install build/release
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

cmake --preset release
cmake --build --preset release
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

cmake --preset release
cmake --build --preset release
```

## Windows

Use Visual Studio 2022 or newer and the vcpkg manifest:

```powershell
cmake -S . -B build\windows -G "Visual Studio 17 2022" -A x64 `
  -DBAYSOR_WITH_TESTS=OFF `
  -DCMAKE_TOOLCHAIN_FILE="$env:VCPKG_ROOT\scripts\buildsystems\vcpkg.cmake" `
  -DVCPKG_TARGET_TRIPLET=x64-windows
cmake --build build\windows --config Release --target baysor --parallel
```

The binary will be under:

```text
build/windows/Release/baysor.exe
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
Windows using the vcpkg manifest. A separate Ubuntu job enables
`BAYSOR_WITH_TESTS=ON`, builds `baysor_tests`, and runs `ctest`.
