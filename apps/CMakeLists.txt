add_executable(driver_vector_add driver_vector_add.cpp)

target_link_libraries(driver_vector_add
  PRIVATE
    kokkos_utils
)

