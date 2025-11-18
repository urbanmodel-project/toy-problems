#include <Kokkos_Core.hpp>
#include <UrbanData.hpp>
#include <UrbanFluxes.hpp>
#include <cmath>
#include <iostream>

namespace URBANXX {
UrbanSurfaceFluxes::UrbanSurfaceFluxes(UrbanSharedDataBundle &bundle)
    : data_bundle(bundle) {
      ForcTemp = bundle.input.ForcTemp;
      ForcTh = bundle.input.ForcPotTemp;
      ForcRho = bundle.input.ForcRho;
      ForcQ = bundle.input.ForcSpcHumd;
      ForcPbot = bundle.input.ForcPress;
      ForcU = bundle.input.ForcWindU;
      ForcV = bundle.input.ForcWindV;
    }
} // namespace URBANXX