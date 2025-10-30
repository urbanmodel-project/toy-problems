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

