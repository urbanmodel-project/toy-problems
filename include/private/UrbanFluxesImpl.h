#ifndef URBAN_FLUXES_HPP
#define URBAN_FLUXES_HPP

#include <private/DataTypesImpl.h>
#include <Kokkos_Core.hpp>
#include <private/UrbanDataImpl.h>

namespace URBANXX {

class UrbanSurface {

public:
  Array1DR8 Temp;
  Array1DR8 es;
  Array1DR8 esdT;
  Array1DR8 qs;
  Array1DR8 qsdT;

  UrbanSurface(int, Real);
  KOKKOS_INLINE_FUNCTION void ComputeQsat(int, Real);
};

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

  UrbanSurface Roof;
  UrbanSurface SunlitWall;
  UrbanSurface ShadeWall;
  UrbanSurface ImperviousRoad;
  UrbanSurface PerviousRoad;

public:
  UrbanSurfaceFluxes(UrbanSharedDataBundle &bundle);
  void computeSurfaceFluxes();
  void computeNewTafAndQaf(int c, Real, Real, Real, Real, Real &, Real &);
  void computeQsatForSurfaces(int);
};

} // namespace URBANXX

#endif