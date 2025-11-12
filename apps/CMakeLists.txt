add_executable(driver_vector_add driver_vector_add.cpp)
add_executable(driver_urban 
  driver_urban.cpp
  ../src/UrbanDataAllocator.cpp
  ../src/UrbanSolarRad.cpp
  ../src/UrbanLongwaveRad.cpp
  )

target_link_libraries(driver_vector_add
  PRIVATE
    kokkos_utils
)

target_include_directories(driver_urban
  PUBLIC
  ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(driver_urban
  PRIVATE
    Kokkos::kokkos
)

