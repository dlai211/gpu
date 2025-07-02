# Copyright (C) 2023 CERN for the benefit of the ATLAS collaboration

# CMake include(s).
cmake_minimum_required( VERSION 3.14 )
include( FetchContent )


# # Silence FetchContent warnings with CMake >=3.24.
if( POLICY CMP0135 )
   cmake_policy( SET CMP0135 NEW )
endif()

# find_package(TBB REQUIRED)

# # Set where to get traccc from.
# # that commit hash points to recent Traccc commit
set( ATLAS_TRACCC_SOURCE
   "GIT_REPOSITORY;https://github.com/acts-project/traccc.git;GIT_TAG;5445580"
   CACHE STRING "Source of the traccc code" )
mark_as_advanced( ATLAS_TRACCC_SOURCE )

# if(NOT EXISTS ${ATLAS_TRACCC_SOURCE})
#    message(FATAL_ERROR "The local traccc source directory does not exist: ${ATLAS_TRACCC_SOURCE}")
# endif()

# Add the local traccc source directory as a subdirectory.
# add_subdirectory(${ATLAS_TRACCC_SOURCE} ${CMAKE_BINARY_DIR}/traccc_binary_dir)


# # Fetch the traccc code.
FetchContent_Declare( traccc ${ATLAS_TRACCC_SOURCE} SYSTEM )

# # Configure the build of traccc.
# Configure the build of traccc.
if(BUILD_CUDA)
  set( TRACCC_BUILD_CUDA TRUE CACHE BOOL "Build traccc with CUDA support" )
else()
  set( TRACCC_BUILD_CUDA FALSE CACHE BOOL "Build traccc with CUDA support" )
  set( VECMEM_BUILD_CUDA_LIBRARY FALSE CACHE BOOL "Build CUDA as part of Vecmem" )
  endif()
if(BUILD_SYCL)
  set( TRACCC_BUILD_SYCL TRUE CACHE BOOL "Build traccc with SYCL support" )
else()
  set( TRACCC_BUILD_SYCL FALSE CACHE BOOL "Build traccc with SYCL support" )
  set( VECMEM_BUILD_SYCL_LIBRARY FALSE CACHE BOOL "Build SYCL as part of Vecmem" )
endif()
set( TRACCC_USE_SYSTEM_THRUST TRUE CACHE BOOL "Build traccc with local Thurst" )
set( TRACCC_BUILD_TESTING FALSE CACHE BOOL "Build traccc tests" )
set( TRACCC_BUILD_BENCHMARKS FALSE CACHE BOOL "Build traccc benchmarks" )
set( TRACCC_BUILD_EXAMPLES FALSE CACHE BOOL "Build traccc examples" )
set( TRACCC_USE_SYSTEM_TBB TRUE CACHE BOOL "Use system TBB" )
set( TRACCC_USE_SYSTEM_EIGEN3 TRUE CACHE BOOL "Use system Eigen" )
set( TRACCC_USE_SYSTEM_ACTS TRUE CACHE BOOL "Use system ACTS" )
set( TRACCC_SETUP_GOOGLETEST FALSE CACHE BOOL "Setup GoogleTest" )
set( TRACCC_SETUP_BENCHMARKS FALSE CACHE BOOL "Setup GoogleBenchmark" )
set( ALGEBRA_PLUGINS_USE_SYSTEM_VC TRUE CACHE BOOL "Use system Vc" )
set( DETRAY_USE_SYSTEM_NLOHMANN TRUE CACHE BOOL "Use system nlohmann_json" )
set( VECMEM_BUILD_HIP_LIBRARY FALSE CACHE BOOL "Build HIP as part of Vecmem" )

atlas_platform_id( ATLAS_PLATFORM )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${ATLAS_PLATFORM}/lib" )

# Configure the build of traccc.
FetchContent_MakeAvailable( traccc )
