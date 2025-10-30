add_executable(kokkos_app main.cpp)

target_link_libraries(kokkos_app
  PRIVATE
    kokkos_utils
)

