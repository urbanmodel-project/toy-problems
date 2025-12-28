enable_language(Fortran)
add_library(kokkos_utils
  vector_add.cpp
)

target_include_directories(kokkos_utils
  PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(kokkos_utils
  PUBLIC
    Kokkos::kokkos
)

target_compile_features(kokkos_utils PUBLIC cxx_std_17)


# Static library for Urban physics
add_library(urban STATIC
  UrbanDataAllocator.cpp
  UrbanSolarRad.cpp
  UrbanLongwaveRad.cpp
  UrbanFluxes.cpp
  Urban.cpp
  KokkosWrappers.cpp
  f90-mod/Urban.F90
)

target_include_directories(urban
  PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(urban
  PUBLIC
    Kokkos::kokkos
)

target_compile_features(urban PUBLIC cxx_std_17)

# Place Fortran modules in a known directory for consumers
set_target_properties(urban PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules
)

