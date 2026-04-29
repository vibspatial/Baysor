cmake_minimum_required(VERSION 3.20)

get_filename_component(_source_dir "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)

if(NOT DEFINED BAYSOR_PRESET OR BAYSOR_PRESET STREQUAL "")
    if(WIN32 AND DEFINED ENV{VCPKG_ROOT} AND NOT "$ENV{VCPKG_ROOT}" STREQUAL "")
        set(BAYSOR_PRESET "user-vcpkg")
    else()
        set(BAYSOR_PRESET "user")
    endif()
endif()

if(BAYSOR_PRESET STREQUAL "user")
    set(_build_dir "${_source_dir}/build/user")
    set(_with_tests OFF)
    set(_with_vcpkg OFF)
elseif(BAYSOR_PRESET STREQUAL "user-vcpkg")
    set(_build_dir "${_source_dir}/build/user-vcpkg")
    set(_with_tests OFF)
    set(_with_vcpkg ON)
elseif(BAYSOR_PRESET STREQUAL "tests")
    set(_build_dir "${_source_dir}/build/tests")
    set(_with_tests ON)
    set(_with_vcpkg OFF)
elseif(BAYSOR_PRESET STREQUAL "vcpkg-tests")
    set(_build_dir "${_source_dir}/build/vcpkg-tests")
    set(_with_tests ON)
    set(_with_vcpkg ON)
else()
    message(FATAL_ERROR
        "Unknown BAYSOR_PRESET='${BAYSOR_PRESET}'. Supported values are: "
        "user, user-vcpkg, tests, vcpkg-tests."
    )
endif()

if(BAYSOR_PRESET MATCHES "vcpkg" AND (NOT DEFINED ENV{VCPKG_ROOT} OR "$ENV{VCPKG_ROOT}" STREQUAL ""))
    message(FATAL_ERROR
        "BAYSOR_PRESET='${BAYSOR_PRESET}' requires VCPKG_ROOT to point to a "
        "bootstrapped vcpkg checkout."
    )
endif()

function(_baysor_run)
    string(REPLACE ";" " " _command_text "${ARGV}")
    message(STATUS "Running: ${_command_text}")
    execute_process(
        COMMAND ${ARGV}
        WORKING_DIRECTORY "${_source_dir}"
        RESULT_VARIABLE _result
    )
    if(_result)
        message(FATAL_ERROR "Command failed with exit code ${_result}: ${ARGV}")
    endif()
endfunction()

find_program(_baysor_ninja NAMES ninja ninja-build)
if(NOT _baysor_ninja)
    if(APPLE)
        set(_baysor_ninja_install_hint [=[
On macOS, install Ninja with Homebrew:
  brew install ninja

Then rerun:
  cmake -P cmake/build_and_install.cmake

If you do not use Homebrew, install Ninja with your package manager and make sure the 'ninja' executable is on PATH.]=])
    elseif(UNIX)
        set(_baysor_ninja_install_hint [=[
Install Ninja with your system package manager, for example:
  sudo apt-get install ninja-build    # Debian/Ubuntu
  sudo dnf install ninja-build        # Fedora
  sudo pacman -S ninja                # Arch

Then rerun:
  cmake -P cmake/build_and_install.cmake]=])
    elseif(WIN32)
        set(_baysor_ninja_install_hint [=[
Install Ninja with your package manager, for example:
  winget install Ninja-build.Ninja
  choco install ninja

Then rerun:
  cmake -P cmake/build_and_install.cmake

Make sure the 'ninja' executable is on PATH.]=])
    else()
        set(_baysor_ninja_install_hint [=[
Install Ninja from https://ninja-build.org/ and make sure the 'ninja'
executable is on PATH. Then rerun:
  cmake -P cmake/build_and_install.cmake]=])
    endif()

    message(FATAL_ERROR
        "Baysor's build helper uses the Ninja CMake generator, but Ninja was "
        "not found on PATH.\n\n"
        "${_baysor_ninja_install_hint}"
    )
endif()

set(_configure_args
    -S "${_source_dir}"
    -B "${_build_dir}"
    -G Ninja
    -DCMAKE_BUILD_TYPE=Release
    -DCMAKE_INSTALL_PREFIX=${_source_dir}/install
    -DBAYSOR_WITH_TESTS=${_with_tests}
)

if(_with_vcpkg)
    list(APPEND _configure_args
        -DCMAKE_TOOLCHAIN_FILE=$ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake
        -DVCPKG_INSTALLED_DIR=${_source_dir}/vcpkg_installed
    )
    if(_with_tests)
        list(APPEND _configure_args -DVCPKG_MANIFEST_FEATURES=tests)
    endif()
endif()

set(_build_targets baysor)
if(_with_tests)
    list(APPEND _build_targets baysor_tests)
endif()

_baysor_run("${CMAKE_COMMAND}" ${_configure_args})
_baysor_run("${CMAKE_COMMAND}" --build "${_build_dir}" --target ${_build_targets})

if(DEFINED BAYSOR_INSTALL_PREFIX AND NOT BAYSOR_INSTALL_PREFIX STREQUAL "")
    get_filename_component(_install_prefix "${BAYSOR_INSTALL_PREFIX}" ABSOLUTE BASE_DIR "${_source_dir}")
    _baysor_run("${CMAKE_COMMAND}" --install "${_build_dir}" --prefix "${_install_prefix}")
else()
    _baysor_run("${CMAKE_COMMAND}" --install "${_build_dir}")
endif()
