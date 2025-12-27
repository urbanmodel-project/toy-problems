# Enable Fortran in this subdirectory so CMake picks up the Fortran compiler
enable_language(Fortran)

add_executable(driver_vector_add driver_vector_add.cpp)
add_executable(driver_urban driver_urban.cpp)
add_executable(driver_urban_c_api driver_urban_c_api.cpp)
add_executable(driver_urban_f
  driver_urban_f.f90
  ../src/urban_f.f90
)

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

target_link_libraries(driver_urban_f
  PRIVATE
    urban
)

# Ensure Fortran linker is used and modules are generated in a known dir
set_target_properties(driver_urban_f PROPERTIES
  LINKER_LANGUAGE Fortran
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules
)

