#include <Kokkos_Core.hpp>
#include "kokkos_utils/vector_add.hpp"
#include <iostream>

int main(int argc, char* argv[]) {
  Kokkos::initialize(argc, argv);
  {
#ifdef KOKKOS_ENABLE_CUDA
    std::cout << "Rank: " << (std::getenv("PMI_RANK") ? std::getenv("PMI_RANK") : "0")
              << " using device: " << Kokkos::Cuda().cuda_device_prop().name
              << " (device " << Kokkos::Cuda().cuda_device() << ")"
              << std::endl;
#endif
    kokkos_utils::vector_add(10);
  }
  Kokkos::finalize();
  return 0;
}

