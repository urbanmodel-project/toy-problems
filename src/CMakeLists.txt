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
  kokkos_c_wrappers.cpp
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

