#ifndef URBAN_DATA_HPP
#define URBAN_DATA_HPP

#include <DataTypes.hpp>
#include <Kokkos_Core.hpp>
#include <cmath>

namespace URBANXX {

// Type Definitions
using ExecutionSpace = Kokkos::DefaultExecutionSpace;

// --- Fortran Type: urbanparams_vars and Canyon Geometry ---
struct CanyonGeometryData {
  HostArray1DR8 CanyonHwrH; // ratio of building height to street width
  Array1DR8 CanyonHwr;      // ratio of building height to street width

  Array1DR8 ViewFactorSkyForRoad;      // view factor of sky for road
  Array1DR8 ViewFactorSkyForWall;      // view factor of sky for one wall
  Array1DR8 ViewFactorRoadToWall;      // view factor from road to wall
  Array1DR8 ViewFactorWallToOtherWall; // view factor from wall to opposing wall
};

// --- Local Variables and Initial Inputs ---
struct SolarInputData {
  HostArray1DR8 CoszenH; // cosine solar zenith angle
  Array1DR8 Coszen;      // cosine solar zenith angle

  HostArray2DR8 SdirHorizH; // direct beam solar radiation on horizontal surface
  Array2DR8 SdirHoriz;      // direct beam solar radiation on horizontal surface

  HostArray2DR8 SdifHorizH; // diffuse solar radiation on horizontal surface
  Array2DR8 SdifHoriz;      // diffuse solar radiation on horizontal surface

  HostArray1DR8 FracSnowH; // fraction of ground covered by snow
  Array1DR8 FracSnow;      // fraction of ground covered by snow
};

// --- Intermediate/Output Fluxes (Shared between subroutines) ---
struct AlbedoFluxes {
  // Incident Fluxes (Output of incident_direct/diffuse)
  Array2DR8 SdirRoad;      // direct incident on road
  Array2DR8 SdifRoad;      // diffuse incident on road
  Array2DR8 SdirSunwall;   // direct incident on sunlit wall
  Array2DR8 SdifSunwall;   // diffuse incident on sunlit wall
  Array2DR8 SdirShadewall; // direct incident on shaded wall
  Array2DR8 SdifShadewall; // diffuse incident on shaded wall

  // Snow-corrected albedos (Output of SnowAlbedo)
  Array2DR8 alb_roof_dir_s;
  Array2DR8 alb_improad_dir_s;
  Array2DR8 alb_perroad_dir_s;
  // ... diffuse versions omitted for brevity ...

  // Absorbed Fluxes (Output of net_solar) (solarabs_vars)
  Array2DR8 sabs_improad_dir; // direct solar absorbed by impervious road
  Array2DR8 sabsSunwall_dir;  // direct solar absorbed by sunwall

  // Reflected Fluxes (Output of net_solar) (Local sref_... variables)
  Array2DR8 sref_roof_dir;   // direct solar reflected by roof (sref_roof_dir)
  Array2DR8 srefSunwall_dir; // direct solar reflected by sunwall
};

// --- The Master Bundle passed between all classes ---
struct UrbanSharedDataBundle {
  CanyonGeometryData geometry;
  SolarInputData input;
  AlbedoFluxes fluxes;
  const int N_LUN;
  const int N_RAD;
};
} // namespace URBANXX
#endif