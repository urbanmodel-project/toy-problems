#ifndef URBAN_FLUXES_HPP
#define URBAN_FLUXES_HPP

#include <DataTypes.hpp>
#include <Kokkos_Core.hpp>
#include <UrbanData.hpp>

namespace URBANXX {

class UrbanSurfaceFluxes {
private:
  Array1DR8 ForcTemp;
  Array1DR8 ForcTh;
  Array1DR8 ForcRho;
  Array1DR8 ForcQ;
  Array1DR8 ForcPbot;
  Array1DR8 ForcU;
  Array1DR8 ForcV;

  Array1DR8 Taf;
  Array1DR8 Qaf;
  UrbanSharedDataBundle data_bundle;

  Array1DR8 Hwr;
  Array1DR8 fPervRoad; // fraction of pervious road w.r.t. total road
  Array1DR8 fRoof;     // fraction of roof w.r.t. grid cell
public:
  UrbanSurfaceFluxes(UrbanSharedDataBundle &bundle);
  void computeSurfaceFluxes();
};

} // namespace URBANXX

#endif