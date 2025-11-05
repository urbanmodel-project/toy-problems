#include <Kokkos_Core.hpp>
#include <UrbanData.hpp>
#include <UrbanDataAllocator.hpp>

namespace URBANXX{
  // Constructor Definition
  UrbanDataAllocator::UrbanDataAllocator(UrbanSharedDataBundle& bundle) 
      : data_bundle(bundle) 
  {
      // Implementation is simple: just initializes the reference member.
  }

  // Method Definition: Contains the heavy lifting of allocation
  void UrbanDataAllocator::allocate_all_views() const {
      int N_LUN = data_bundle.N_LUN;
      int N_RAD = data_bundle.N_RAD;

      printf("Starting device memory allocation for %d land units and %d bands.\n", N_LUN, N_RAD);

      // --- CanyonGeometryData Allocation ---
      data_bundle.geometry.canyon_hwr = Array1DR8("canyon_hwr", N_LUN);
      data_bundle.geometry.vf_sr = Array1DR8("vf_sr", N_LUN);
      // ... all other geometry views ...

      // --- SolarInputData Allocation ---
      data_bundle.input.coszen = Array1DR8("coszen", N_LUN);
      data_bundle.input.sdir_h = Array2DR8("sdir_h", N_LUN, N_RAD);
      // ... all other input views ...
      
      // --- AlbedoFluxes Allocation ---
      data_bundle.fluxes.sdir_road = Array2DR8("sdir_road", N_LUN, N_RAD);
      // ... all other flux views ...

      printf("All primary Views successfully allocated on device.\n");
  }

}