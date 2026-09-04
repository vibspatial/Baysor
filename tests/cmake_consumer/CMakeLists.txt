cmake_minimum_required(VERSION 3.20)
project(baysor_add_subdirectory_consumer LANGUAGES CXX)

if(NOT DEFINED BAYSOR_SOURCE_DIR)
    message(FATAL_ERROR "Configure with -DBAYSOR_SOURCE_DIR=/path/to/Baysor")
endif()

add_subdirectory("${BAYSOR_SOURCE_DIR}" baysor EXCLUDE_FROM_ALL)

add_executable(baysor_consumer main.cpp)
target_link_libraries(baysor_consumer PRIVATE baysor::baysor)
