cmake_minimum_required(VERSION 3.11)

project(modular)
include(ExternalProject)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_LIBRARY_OUTPUT_DERICTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DERICTORY ${CMAKE_BINARY_DIR}/bin)

# SET(CMAKE_CXX_FLAGS "-O3 -Wall -frtti -Wextra")
SET(CMAKE_CXX_FLAGS "-O3")
SET(CMAKE_BUILD_TYPE "Release")
# SET(CMAKE_BUILD_TYPE "Debug")
include_directories("/opt/homebrew/include")
link_directories("/opt/homebrew/lib/")

add_subdirectory(src)

include(FetchContent)
find_package(Git REQUIRED)
set(FETCHCONTENT_QUIET ON)
set(FETCHCONTENT_UPDATES_DISCONNECTED ON)
set(BUILD_SHARED_LIBS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)


set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "" FORCE)
set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
FetchContent_Declare(
        benchmark
        GIT_REPOSITORY "https://github.com/google/benchmark.git"
        GIT_TAG aa59d40f8822358d012e854488c2e52857b8a1f0
        PATCH_COMMAND ""
)
FetchContent_MakeAvailable(benchmark)

add_subdirectory(benchmark)

# main loop: add executables for all example files 
file (GLOB EXAMPLES_SRC_FILES CONFIGURE_DEPENDS examples/*.cpp)
foreach (app ${EXAMPLES_SRC_FILES})
    get_filename_component ( exe ${app} NAME_WE )
    add_executable ( ${exe} ${app} )
    set_property(TARGET ${exe} PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DERICTORY}/examples)
    target_link_libraries(${exe} PUBLIC modular)
endforeach()

# test structure organization
include(CTest)
if(BUILD_TESTING)
    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.14.0  
    )
    # Для MSVC: отключим установку GoogleTest
    set(INSTALL_GTEST OFF CACHE BOOL "Disable install" FORCE)
    FetchContent_MakeAvailable(googletest)

    add_subdirectory(test)
    
endif()