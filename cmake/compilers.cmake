# cmake/compilers.cmake
message(STATUS "=== Detecting compiler and backend for Kokkos ===")

# -------------------------------------------------------------------
# Step 1. Detect available toolchains
# -------------------------------------------------------------------
find_program(NVCC nvcc)
find_program(HIPCC hipcc)
find_program(GPP g++)
find_program(CLANGPP clang++)

set(DETECTED_BACKEND "NONE")

if(NVCC)
  set(DETECTED_BACKEND "CUDA")
elseif(HIPCC)
  set(DETECTED_BACKEND "HIP")
elseif(GPP OR CLANGPP)
  set(DETECTED_BACKEND "OPENMP")
endif()

# -------------------------------------------------------------------
# Step 2. Apply user override (if any)
# -------------------------------------------------------------------
if(Kokkos_ENABLE_CUDA)
  set(DETECTED_BACKEND "CUDA")
elseif(Kokkos_ENABLE_HIP)
  set(DETECTED_BACKEND "HIP")
elseif(Kokkos_ENABLE_OPENMP)
  set(DETECTED_BACKEND "OPENMP")
elseif(Kokkos_ENABLE_SERIAL)
  set(DETECTED_BACKEND "SERIAL")
endif()

message(STATUS "Detected backend: ${DETECTED_BACKEND}")

# -------------------------------------------------------------------
# Step 3. Choose compiler and set flags
# -------------------------------------------------------------------
if(NOT CMAKE_CXX_COMPILER)
  if(DETECTED_BACKEND STREQUAL "CUDA")
    set(Kokkos_ENABLE_CUDA ON CACHE BOOL "" FORCE)
    set(CMAKE_CXX_COMPILER ${CMAKE_SOURCE_DIR}/external/kokkos/bin/nvcc_wrapper CACHE FILEPATH "C++ compiler")
    message(STATUS "Using NVCC wrapper for CUDA backend: ${CMAKE_CXX_COMPILER}")

  elseif(DETECTED_BACKEND STREQUAL "HIP")
    set(Kokkos_ENABLE_HIP ON CACHE BOOL "" FORCE)
    set(CMAKE_CXX_COMPILER ${HIPCC} CACHE FILEPATH "C++ compiler")
    message(STATUS "Using HIP compiler: ${HIPCC}")

  elseif(DETECTED_BACKEND STREQUAL "OPENMP")
    set(Kokkos_ENABLE_OPENMP ON CACHE BOOL "" FORCE)
    if(GPP)
      set(CMAKE_CXX_COMPILER ${GPP} CACHE FILEPATH "C++ compiler")
    elseif(CLANGPP)
      set(CMAKE_CXX_COMPILER ${CLANGPP} CACHE FILEPATH "C++ compiler")
    endif()
    message(STATUS "Using OpenMP compiler: ${CMAKE_CXX_COMPILER}")

  else()
    # No GPU or OpenMP detected — fallback to serial
    set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "" FORCE)
    if(GPP)
      set(CMAKE_CXX_COMPILER ${GPP} CACHE FILEPATH "C++ compiler")
    elseif(CLANGPP)
      set(CMAKE_CXX_COMPILER ${CLANGPP} CACHE FILEPATH "C++ compiler")
    else()
      message(FATAL_ERROR "No suitable compiler found.")
    endif()
    message(STATUS "Falling back to Serial backend.")
  endif()
endif()

message(STATUS "CMAKE_CXX_COMPILER = ${CMAKE_CXX_COMPILER}")

