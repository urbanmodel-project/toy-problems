#include <Kokkos_Core.hpp>
#include <UrbanData.hpp>
#include <UrbanDataAllocator.hpp>

namespace URBANXX {
// Constructor Definition
UrbanDataAllocator::UrbanDataAllocator(UrbanSharedDataBundle &bundle)
    : data_bundle(bundle) {
  // Implementation is simple: just initializes the reference member.
}

// Method Definition: Contains the heavy lifting of allocation
void UrbanDataAllocator::allocate_all_views() const {
  int N_LUN = data_bundle.N_LUN;
  int N_RAD = data_bundle.N_RAD;

  printf("Starting device memory allocation for %d land units and %d bands.\n",
         N_LUN, N_RAD);

  // --- CanyonGeometryData Allocation ---
  data_bundle.geometry.CanyonHwr = Array1DR8("CanyonHwr", N_LUN);
  data_bundle.geometry.ViewFactorSkyForRoad =
      Array1DR8("ViewFactorSkyForRoad", N_LUN);
  // ... all other geometry views ...

  // --- SolarInputData Allocation ---
  data_bundle.input.coszen = Array1DR8("coszen", N_LUN);
  data_bundle.input.SdirHoriz = Array2DR8("SdirHoriz", N_LUN, N_RAD);
  // ... all other input views ...

  // --- AlbedoFluxes Allocation ---
  data_bundle.fluxes.SdirRoad = Array2DR8("SdirRoad", N_LUN, N_RAD);
  // ... all other flux views ...

  printf("All primary Views successfully allocated on device.\n");
}

} // namespace URBANXX