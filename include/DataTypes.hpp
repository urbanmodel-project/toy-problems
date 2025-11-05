#ifndef URBANXX_DATA_TYPES_HPP
#define URBANXX_DATA_TYPES_HPP
#include "Kokkos_Core.hpp"

namespace URBANXX {

// Standard integer and floating point types
using I4 = std::int32_t; ///< alias for 32-bit integer
using I8 = std::int64_t; ///< alias for 64-bit integer
using R4 = float;        ///< alias for 32-bit (single prec) real
using R8 = double;       ///< alias for 64-bit (double prec) real

/// generic real 64-bit (default) or 32-bit (if -DSINGLE_PRECISION used)
#ifdef URBANXX_SINGLE_PRECISION
using Real = float;
#else
using Real = double;
#endif

// user-defined literal for generic reals
KOKKOS_INLINE_FUNCTION constexpr Real operator""_Real(long double x) {
   return x;
}

// Aliases for Kokkos memory spaces
#ifdef URBANXX_ENABLE_CUDA
using MemSpace = Kokkos::CudaSpace;
#elif URBANXX_ENABLE_HIP
using MemSpace = Kokkos::Experimental::HIPSpace;
#elif URBANXX_ENABLE_SYCL
using MemSpace = Kokkos::Experimental::SYCLDeviceUSMSpace;
#elif URBANXX_ENABLE_OPENMP
using MemSpace = Kokkos::HostSpace;
#elif URBANXX_ENABLE_SERIAL
using MemSpace = Kokkos::HostSpace;
#else
#error "None of URBANXX_ENABLE_CUDA, URBANXX_ENABLE_HIP, URBANXX_ENABLE_OPENMP, \
URBANXX_ENABLE_SERIAL is defined."
#endif

// Aliases for Kokkos memory layouts
#ifdef URBANXX_LAYOUT_RIGHT

using MemLayout    = Kokkos::LayoutRight;
using MemInvLayout = Kokkos::LayoutLeft;

#elif URBANXX_LAYOUT_LEFT

using MemLayout    = Kokkos::LayoutLeft;
using MemInvLayout = Kokkos::LayoutRight;

#else
#error "URBANXX Memory Layout is not defined."
#endif

using HostMemSpace     = Kokkos::HostSpace;
using HostMemLayout    = MemLayout;
using HostMemInvLayout = MemInvLayout;

#define MAKE_URBANXX_VIEW_DIMS(N, V, T, ML, MS)  \
   using N##1D##T = Kokkos::V<T *, ML, MS>;    \
   using N##2D##T = Kokkos::V<T **, ML, MS>;   \
   using N##3D##T = Kokkos::V<T ***, ML, MS>;  \
   using N##4D##T = Kokkos::V<T ****, ML, MS>; \
   using N##5D##T = Kokkos::V<T *****, ML, MS>;

#define MAKE_URBANXX_VIEW_TYPES(N, V, ML, MS) \
   MAKE_URBANXX_VIEW_DIMS(N, V, I4, ML, MS)   \
   MAKE_URBANXX_VIEW_DIMS(N, V, I8, ML, MS)   \
   MAKE_URBANXX_VIEW_DIMS(N, V, R4, ML, MS)   \
   MAKE_URBANXX_VIEW_DIMS(N, V, R8, ML, MS)   \
   MAKE_URBANXX_VIEW_DIMS(N, V, Real, ML, MS)

// Aliases for Kokkos device arrays of various dimensions and types
MAKE_URBANXX_VIEW_TYPES(Array, View, MemLayout, MemSpace)

// Aliases for Kokkos host arrays of various dimensions and types
MAKE_URBANXX_VIEW_TYPES(HostArray, View, HostMemLayout, HostMemSpace)

#undef MAKE_URBANXX_VIEW_TYPES
#undef MAKE_URBANXX_VIEW_DIMS

template <class T>
inline constexpr bool isKokkosArray = Kokkos::is_view<T>::value;

} // end namespace URBANXX

//===----------------------------------------------------------------------===//
#endif
