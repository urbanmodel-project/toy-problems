#include <Kokkos_Core.hpp>
#include "kokkos_utils/vector_add.hpp"

int main(int argc, char* argv[]) {
  Kokkos::initialize(argc, argv);
  {
    kokkos_utils::vector_add(10);
  }
  Kokkos::finalize();
  return 0;
}

