add_executable(driver_vector_add driver_vector_add.cpp)
add_executable(driver_urban driver_urban.cpp)
add_executable(driver_urban_c_api driver_urban_c_api.cpp)

target_link_libraries(driver_vector_add
  PRIVATE
    kokkos_utils
)

target_link_libraries(driver_urban
  PRIVATE
    urban
)

target_link_libraries(driver_urban_c_api
  PRIVATE
    urban
)

