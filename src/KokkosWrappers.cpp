#include <Kokkos_Core.hpp>
extern "C" {
void UrbanKokkosInitialize() {
  int argc = 0;
  char **argv = nullptr;
  Kokkos::initialize(argc, argv);
}
void UrbanKokkosFinalize() { Kokkos::finalize(); }
}
