// Common data type aliases for Kokkos views used across Urban module.
#ifndef URBAN_PRIVATE_DATATYPES_HPP
#define URBAN_PRIVATE_DATATYPES_HPP

#include <Kokkos_Core.hpp>

// Scalar types
using R8 = double;
using Real = double;

// Host/device Kokkos view aliases
using Array1DR8 = Kokkos::View<R8 *, Kokkos::DefaultExecutionSpace>;
using HostArray1DR8 = Kokkos::View<R8 *, Kokkos::HostSpace>;

using Array2DR8 = Kokkos::View<R8 **, Kokkos::DefaultExecutionSpace>;
using HostArray2DR8 = Kokkos::View<R8 **, Kokkos::HostSpace>;

using Array3DR8 = Kokkos::View<R8 ***, Kokkos::DefaultExecutionSpace>;
using HostArray3DR8 = Kokkos::View<R8 ***, Kokkos::HostSpace>;

#endif // URBAN_PRIVATE_DATATYPES_HPP
