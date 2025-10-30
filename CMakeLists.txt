cmake_minimum_required(VERSION 3.18)
project(KokkosProject LANGUAGES CXX)

# -------------------------------------------------------------------
# Backend options
# -------------------------------------------------------------------
option(Kokkos_ENABLE_OPENMP "Enable Kokkos OpenMP backend" ON)
option(Kokkos_ENABLE_CUDA   "Enable Kokkos CUDA backend"   OFF)
option(Kokkos_ENABLE_HIP    "Enable Kokkos HIP backend"    OFF)
option(Kokkos_ENABLE_SERIAL "Enable Kokkos Serial backend" OFF)

# -------------------------------------------------------------------
# Detect / set compiler
# -------------------------------------------------------------------
include(${CMAKE_SOURCE_DIR}/cmake/compilers.cmake)

# -------------------------------------------------------------------
# Configure C++ standard
# -------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# -------------------------------------------------------------------
# Add Kokkos submodule
# -------------------------------------------------------------------
add_subdirectory(${CMAKE_SOURCE_DIR}/externals/kokkos)

# -------------------------------------------------------------------
# Add your code
# -------------------------------------------------------------------
add_subdirectory(src)
add_subdirectory(apps)

