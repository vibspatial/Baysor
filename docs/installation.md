# Installation

The C++ branch is built with CMake.

## Required Dependencies

Install a C++ toolchain plus the libraries required by
[CMakeLists.txt](../CMakeLists.txt):

- Eigen3
- OpenMP
- spdlog
- CGAL
- Arrow / Parquet
- HDF5
- nlohmann_json
- libtiff
- GTest for the optional test target

Several header-only dependencies are fetched automatically by CMake.

## Ubuntu 24.04

The repo CI and Docker build use Ubuntu 24.04. The commands below match that
working setup.

### 1. Install The Apache Arrow Apt Source

Ubuntu's default repositories do not always provide the Arrow / Parquet
development packages in the form expected by this build, so Baysor uses the
Apache Arrow apt source.

```bash
sudo apt-get update
sudo apt-get install -y --no-install-recommends \
  ca-certificates \
  lsb-release \
  wget

wget https://packages.apache.org/artifactory/arrow/$(lsb_release --id --short | tr 'A-Z' 'a-z')/apache-arrow-apt-source-latest-$(lsb_release --codename --short).deb
sudo apt-get install -y --no-install-recommends ./apache-arrow-apt-source-latest-$(lsb_release --codename --short).deb
rm ./apache-arrow-apt-source-latest-$(lsb_release --codename --short).deb
```

### 2. Install Build Dependencies

```bash
sudo apt-get update
sudo apt-get install -y --no-install-recommends \
  build-essential \
  cmake \
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
```

### 3. Configure And Build The Binary

```bash
cmake -S . -B build 
cmake --build build --target baysor -j"$(nproc)"
```

The main binary will be:

```text
./build/baysor
```

### 4. Optional: Build And Run Tests

```bash
sudo apt-get install -y --no-install-recommends libgtest-dev
cmake -S . -B build
cmake --build build --target baysor baysor_tests -j"$(nproc)"
ctest --test-dir build --output-on-failure
```

## Configure And Build

If the required libraries are already installed on your system:

```bash
cmake -S . -B build -DBAYSOR_WITH_TESTS=OFF
cmake --build build --target baysor -j"$(nproc)"
```

If you want tests too:

```bash
cmake -S . -B build
cmake --build build --target baysor baysor_tests -j"$(nproc)"
ctest --test-dir build --output-on-failure
```

## Common Commands

Show CLI help:

```bash
./build/baysor --help
```

Show run-specific help:

```bash
./build/baysor run --help
```

## Notes

- Xenium support requires the normal build dependencies only; there is no
  separate Xenium-specific build target.
- The `parquet` output style relies on Arrow / Parquet and HDF5 being available
  at build time.
