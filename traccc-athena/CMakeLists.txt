# Copyright (C) 2023-2025 CERN for the benefit of the ATLAS collaboration

# Set up the project.
cmake_minimum_required(VERSION 3.20)
project(traccc_athena VERSION 0.0.1 LANGUAGES CXX)

# Warning => error
set(CMAKE_COMPILE_WARNING_AS_ERROR ON)
# -O3 for releases
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
# CUDA: no deprecated GPU targets
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} $ENV{CUDAFLAGS} -Wno-deprecated-gpu-targets")

# # Use C++20 with CUDA. Athena by default still uses C++14. But traccc needs
set(CMAKE_CUDA_STANDARD 20 CACHE STRING "CUDA C++ standard")

# Find the base release.
find_package(Athena 25.0 REQUIRED)

# Generate a compile_commands.json file by default.
set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE CACHE BOOL
   "Generate compile_commands.json" FORCE)

# Set up CTest.
atlas_ctest_setup()

# Build options
option(BUILD_CUDA "Build the CUDA version" ON)
option(BUILD_SYCL "Build the SYCL version" OFF)

if(NOT BUILD_CUDA AND NOT BUILD_SYCL)
   message(FATAL_ERROR "no implementation were enabled: set BUILD_CUDA or BUILD_SYCL")
endif()

# Decide where to take traccc from.
option(ATLAS_USE_SYSTEM_TRACCC "Use an existing traccc installation" OFF)
if(ATLAS_USE_SYSTEM_TRACCC)
   find_package(traccc REQUIRED)
else()
   include(cmake/build-traccc.cmake)
endif()

# Set up the analysis project from the repository.
atlas_project(USE Athena ${Athena_VERSION})

# Generate an environment setup script. This call makes sure that whatever
# "external" we use in this analysis project, is set up correctly in the
# runtime environment.
lcg_generate_env(SH_FILE "${CMAKE_BINARY_DIR}/${ATLAS_PLATFORM}/env_setup.sh")
install(FILES "${CMAKE_BINARY_DIR}/${ATLAS_PLATFORM}/env_setup.sh"
   DESTINATION .)

# Set up CPack. This is necessary to be able to build TGZ files correctly from
# the compiled code. Making it possible to use the code on the grid.
atlas_cpack_setup()
