FROM ubuntu:24.04 AS build

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
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
    libtiff-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /src
COPY . .

RUN cmake -S . -B build -GNinja -DBAYSOR_WITH_TESTS=OFF
RUN cmake --build build --target baysor -j"$(nproc)"

FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    libeigen3-dev \
    libomp-dev \
    libspdlog-dev \
    libcgal-dev \
    libarrow-dev \
    libparquet-dev \
    libhdf5-dev \
    nlohmann-json3-dev \
    libtiff-dev \
    && rm -rf /var/lib/apt/lists/*

COPY --from=build /src/build/baysor /usr/local/bin/baysor

ENTRYPOINT ["/usr/local/bin/baysor"]
