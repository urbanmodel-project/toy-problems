#include <Kokkos_Core.hpp>
extern "C" {
void urban_kokkos_initialize() {
  int argc = 0;
  char **argv = nullptr;
  Kokkos::initialize(argc, argv);
}
void urban_kokkos_finalize() { Kokkos::finalize(); }
}
