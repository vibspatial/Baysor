cmake_minimum_required(VERSION 3.20)

get_filename_component(_source_dir "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)

if(NOT DEFINED BAYSOR_PRESET OR BAYSOR_PRESET STREQUAL "")
    if(DEFINED ENV{VCPKG_ROOT} AND NOT "$ENV{VCPKG_ROOT}" STREQUAL "")
        set(BAYSOR_PRESET "user-vcpkg")
    else()
        set(BAYSOR_PRESET "user")
    endif()
endif()

if(BAYSOR_PRESET STREQUAL "user")
    set(_build_dir "${_source_dir}/build/user")
elseif(BAYSOR_PRESET STREQUAL "user-vcpkg")
    set(_build_dir "${_source_dir}/build/user-vcpkg")
elseif(BAYSOR_PRESET STREQUAL "tests")
    set(_build_dir "${_source_dir}/build/tests")
elseif(BAYSOR_PRESET STREQUAL "vcpkg-tests")
    set(_build_dir "${_source_dir}/build/vcpkg-tests")
else()
    message(FATAL_ERROR
        "Unknown BAYSOR_PRESET='${BAYSOR_PRESET}'. Supported values are: "
        "user, user-vcpkg, tests, vcpkg-tests."
    )
endif()

if(BAYSOR_PRESET MATCHES "^vcpkg" AND (NOT DEFINED ENV{VCPKG_ROOT} OR "$ENV{VCPKG_ROOT}" STREQUAL ""))
    message(FATAL_ERROR
        "BAYSOR_PRESET='${BAYSOR_PRESET}' requires VCPKG_ROOT to point to a "
        "bootstrapped vcpkg checkout."
    )
endif()

function(_baysor_run)
    message(STATUS "Running: ${ARGV}")
    execute_process(
        COMMAND ${ARGV}
        WORKING_DIRECTORY "${_source_dir}"
        RESULT_VARIABLE _result
    )
    if(_result)
        message(FATAL_ERROR "Command failed with exit code ${_result}: ${ARGV}")
    endif()
endfunction()

_baysor_run("${CMAKE_COMMAND}" --preset "${BAYSOR_PRESET}")
_baysor_run("${CMAKE_COMMAND}" --build --preset "${BAYSOR_PRESET}")

if(DEFINED BAYSOR_INSTALL_PREFIX AND NOT BAYSOR_INSTALL_PREFIX STREQUAL "")
    get_filename_component(_install_prefix "${BAYSOR_INSTALL_PREFIX}" ABSOLUTE BASE_DIR "${_source_dir}")
    _baysor_run("${CMAKE_COMMAND}" --install "${_build_dir}" --prefix "${_install_prefix}")
else()
    _baysor_run("${CMAKE_COMMAND}" --install "${_build_dir}")
endif()
