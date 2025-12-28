#include "kokkos_utils/vector_add.hpp"
#include <Kokkos_Core.hpp>
#include <iostream>
#include <mpi.h>

int main(int argc, char *argv[]) {
  // Initialize MPI
  MPI_Init(&argc, &argv);

  int mpi_rank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  Kokkos::initialize(argc, argv);
  {

    // Safest: host-only print guard
    if (mpi_rank == 0) {
// Ensure only master thread prints
#pragma omp master
      {
        std::cout << "=== Kokkos Configuration (rank 0 only) ===" << std::endl;
        Kokkos::print_configuration(std::cout);
#ifdef KOKKOS_ENABLE_CUDA
        std::cout << "Rank: "
                  << (std::getenv("PMI_RANK") ? std::getenv("PMI_RANK") : "0")
                  << " using device: " << Kokkos::Cuda().cuda_device_prop().name
                  << " (device " << Kokkos::Cuda().cuda_device() << ")"
                  << std::endl;
#endif
      }
    }

    kokkos_utils::vector_add(10);
  }
  Kokkos::finalize();
  MPI_Finalize();
  return 0;
}
