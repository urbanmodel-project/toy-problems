#include "kokkos_utils/vector_add.hpp"
#include <private/DataTypesImpl.h>
#include <iostream>

namespace kokkos_utils {

void vector_add(int N) {
  Kokkos::View<double *> A("A", N);
  Kokkos::View<double *> B("B", N);
  Kokkos::View<double *> C("C", N);

  Kokkos::parallel_for(
      "InitAB", N, KOKKOS_LAMBDA(const int i) {
        A(i) = i * 1.0;
        B(i) = i * 2.0;
      });

  Kokkos::parallel_for(
      "VecAdd", N, KOKKOS_LAMBDA(const int i) { C(i) = A(i) + B(i); });

  auto hC = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), C);
  std::cout << "C[0] = " << hC(0) << ", C[N-1] = " << hC(N - 1) << std::endl;
}

} // namespace kokkos_utils
