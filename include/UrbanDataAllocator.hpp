// UranDataAllocator.hpp

#ifndef URBAN_DATA_ALLOCATOR_HPP
#define URBAN_DATA_ALLOCATOR_HPP

#include <UrbanData.hpp> // Required to know the definition of UrbanSharedDataBundle

namespace URBANXX {
/**
 * @class DataAllocator
 * @brief Handles the one-time allocation of all Kokkos Views within the
 * UrbanSharedDataBundle.
 */
class UrbanDataAllocator {
private:
  // Holds a reference to the main data structure to be allocated.
  UrbanSharedDataBundle &data_bundle;

public:
  /**
   * @brief Constructor: Stores a non-const reference to the bundle.
   */
  UrbanDataAllocator(UrbanSharedDataBundle &bundle);

  /**
   * @brief Allocates memory for all Kokkos Views in the bundle on the device.
   */
  void allocate_all_views() const;

  // Optional: Methods for allocating specific sub-groups of views
  // void allocate_geometry() const;
  // void allocate_fluxes() const;
};
} // namespace URBANXX

#endif // DATA_ALLOCATOR_HPP