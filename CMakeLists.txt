cmake_minimum_required(VERSION 3.20)
project(baysor VERSION 0.8.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Prefer optimized single-config builds unless the user explicitly chose
# something else. This does not affect multi-config generators.
if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Options
option(BAYSOR_WITH_CUDA "Enable CUDA/GPU support" OFF)
option(BAYSOR_WITH_TESTS "Build tests" ON)
option(BAYSOR_WITH_REPORTING "Build reporting/visualization support" ON)

# Dependencies managed by vcpkg or find_package
# Core (always required)
find_package(Eigen3 3.3 REQUIRED)
find_package(OpenMP REQUIRED)
find_package(spdlog REQUIRED)

# UMAP and dependencies (all header-only, fetched via FetchContent).
# All deps are pre-declared with pinned tags so that FetchContent dedup
# prevents transitive CMakeLists from fetching "master".
include(FetchContent)
FetchContent_Declare(aarand
    GIT_REPOSITORY https://github.com/LTLA/aarand
    GIT_TAG v1.0.2
)
# knncolle v2.3.0 transitively depends on kmeans v3.x and subpar v0.2.1
FetchContent_Declare(kmeans
    GIT_REPOSITORY https://github.com/LTLA/CppKmeans
    GIT_TAG v3.1.1
)
FetchContent_Declare(subpar
    GIT_REPOSITORY https://github.com/LTLA/subpar
    GIT_TAG v0.3.1
)
FetchContent_Declare(knncolle
    GIT_REPOSITORY https://github.com/knncolle/knncolle
    GIT_TAG v2.3.0
)
FetchContent_Declare(irlba
    GIT_REPOSITORY https://github.com/LTLA/CppIrlba
    GIT_TAG v2.0.2
)
FetchContent_Declare(umappp
    GIT_REPOSITORY https://github.com/LTLA/umappp
    GIT_TAG v2.0.1
)
# Populate all in one call; order determines FetchContent dedup priority.
FetchContent_MakeAvailable(aarand kmeans subpar knncolle irlba umappp)

# Patch irlba v2.0.2 bug: Eigen::VectorXd has no begin()/end().
# Replace std::max_element(BS.begin(), BS.end()) with BS.maxCoeff().
set(_irlba_compute "${irlba_SOURCE_DIR}/include/irlba/compute.hpp")
if(EXISTS "${_irlba_compute}")
    file(READ "${_irlba_compute}" _irlba_content)
    string(REPLACE
        "*std::max_element(BS.begin(), BS.end())"
        "BS.maxCoeff()"
        _irlba_content "${_irlba_content}")
    file(WRITE "${_irlba_compute}" "${_irlba_content}")
endif()

# Spatial
find_package(CGAL REQUIRED)
# nanoflann is header-only, vendored in third_party/ or found via find_package
# find_package(nanoflann REQUIRED)

# I/O
find_package(Arrow REQUIRED)
find_package(Parquet REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS C)
find_package(nlohmann_json REQUIRED)
find_package(TIFF REQUIRED)

# CLI
# CLI11 and toml++ are header-only
# find_package(CLI11 REQUIRED)
# find_package(tomlplusplus REQUIRED)

if(BAYSOR_WITH_CUDA)
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
    # find_package(faiss REQUIRED)
    add_compile_definitions(BAYSOR_CUDA_ENABLED)
endif()

# --- Library target ---
add_library(baysor_lib
    # Utils
    src/utils/general.cpp
    src/utils/options.cpp
    src/utils/xenium.cpp

    # Data loading
    src/data_loading/data.cpp
    src/data_loading/prior_segmentation.cpp

    # Processing - distributions
    src/processing/distributions/mv_normal.cpp
    src/processing/distributions/categorical_smoothed.cpp

    # Processing - models
    src/processing/models/adj_list.cpp
    src/processing/models/component.cpp
    src/processing/models/bmm_data.cpp

    # Processing - data_processing
    src/processing/data_processing/triangulation.cpp
    src/processing/data_processing/initialization.cpp
    src/processing/data_processing/neighborhood_composition.cpp
    src/processing/data_processing/noise_estimation.cpp
    src/processing/data_processing/boundary_estimation.cpp
    src/processing/data_processing/umap_wrappers.cpp

    # Processing - bmm_algorithm
    src/processing/bmm_algorithm/bmm_algorithm.cpp
    src/processing/bmm_algorithm/molecule_clustering.cpp
    src/processing/bmm_algorithm/compartment_segmentation.cpp
    src/processing/bmm_algorithm/tracing.cpp
    src/processing/bmm_algorithm/history_analysis.cpp

    # Processing - utils
    src/processing/utils/utils.cpp
    src/processing/utils/convex_hull.cpp

    # Reporting
    src/reporting/color_utils.cpp
    src/reporting/output.cpp
    src/reporting/preview_report.cpp
    src/reporting/run_report.cpp
)

target_include_directories(baysor_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${HDF5_INCLUDE_DIRS}
)

target_link_libraries(baysor_lib PUBLIC
    Eigen3::Eigen
    OpenMP::OpenMP_CXX
    spdlog::spdlog
    CGAL::CGAL
    Arrow::arrow_shared
    Parquet::parquet_shared
    ${HDF5_C_LIBRARIES}
    nlohmann_json::nlohmann_json
    TIFF::TIFF
    umappp
    knncolle::knncolle
)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(baysor_lib PRIVATE -Wno-stringop-overflow)
endif()

if(BAYSOR_WITH_CUDA)
    target_link_libraries(baysor_lib PUBLIC
        CUDA::cudart
        CUDA::cusparse
    )
endif()

# --- Executable ---
add_executable(baysor
    src/cli/main.cpp
)

target_link_libraries(baysor PRIVATE baysor_lib)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(baysor PRIVATE -Wno-stringop-overflow)
endif()

# --- Tests ---
if(BAYSOR_WITH_TESTS)
    enable_testing()
    find_package(GTest REQUIRED)

    add_executable(baysor_tests
        tests/test_main.cpp
    )
    target_link_libraries(baysor_tests PRIVATE baysor_lib GTest::gtest_main)
    target_compile_definitions(baysor_tests PRIVATE
        BAYSOR_TEST_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
    )
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        target_compile_options(baysor_tests PRIVATE -Wno-stringop-overflow)
    endif()
    add_test(NAME baysor_tests COMMAND baysor_tests)
endif()
