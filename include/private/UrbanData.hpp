#ifndef URBAN_DATA_HPP
#define URBAN_DATA_HPP

#include <private/DataTypes.hpp>
#include <Kokkos_Core.hpp>
#include <cmath>

namespace URBANXX {

// Type Definitions
using ExecutionSpace = Kokkos::DefaultExecutionSpace;

// --- Fortran Type: urbanparams_vars and Canyon Geometry ---
struct CanyonGeometryData {
  HostArray1DR8 CanyonHwrH; // ratio of building height to street width
  Array1DR8 CanyonHwr;      // ratio of building height to street width

  Array1DR8 ViewFactorSkyFromRoad;  // view factor of sky for road
  Array1DR8 ViewFactorSkyFromWall;  // view factor of sky for one wall
  Array1DR8 ViewFactorRoadFromWall; // view factor from road to wall
  Array1DR8
      ViewFactorOtherWallFromWall;  // view factor from wall to opposing wall
  Array1DR8 ViewFactorWallFromRoad; // view factor of one wall to road
};

#define DECLARE_DUAL_ARRAY(suffix, name)                                       \
  HostArray##suffix name##H;                                                   \
  Array##suffix name;

// --- Local Variables and Initial Inputs ---
struct AtmosphereInputData {
  // cosine solar zenith angle (-)
  DECLARE_DUAL_ARRAY(1DR8, Coszen);

  // direct beam solar radiation on horizontal surface (W/m**2)
  DECLARE_DUAL_ARRAY(2DR8, SdirHoriz);

  // diffuse solar radiation on horizontal surface (W/m**2)
  DECLARE_DUAL_ARRAY(2DR8, SdifHoriz);

  // fraction of ground covered by snow (-)
  DECLARE_DUAL_ARRAY(1DR8, FracSnow);

  // downwelling longwave radiation (W/m**2)
  DECLARE_DUAL_ARRAY(1DR8, DownwellingLongRad);

  // air temperature (K)
  DECLARE_DUAL_ARRAY(1DR8, ForcTemp);

  // potential temperature (Pa)
  DECLARE_DUAL_ARRAY(1DR8, ForcPotTemp);

  // air density (kg/m**3)
  DECLARE_DUAL_ARRAY(1DR8, ForcRho);

  // specific humidity (kg/kg)
  DECLARE_DUAL_ARRAY(1DR8, ForcSpcHumd);

  // atomspheric pressure (Pa)
  DECLARE_DUAL_ARRAY(1DR8, ForcPress);

  // wind speed in east direction (m/s)
  DECLARE_DUAL_ARRAY(1DR8, ForcWindU);

  // wind speed in north direction (m/s)
  DECLARE_DUAL_ARRAY(1DR8, ForcWindV);
};

struct CombinedRoadDataType {
  Array3DR8
      DownwellingShortRad; // downwelling shortwave radiation per unit road area
  Array3DR8 SnowAlbedo;    // snow albedo
  Array3DR8 AlbedoWithSnowEffects; // albedo of road including snow effects
};

struct RoadDataType {
  Array3DR8 SnowAlbedo;            // snow albedo
  Array3DR8 AlbedoWithSnowEffects; // albedo of road including snow effects
  Array3DR8 ReflectedShortRad; // reflected solar radiation per unit ground area
                               // per unit incident flux
  Array3DR8 BaseAlbedo;        // albedo of road without snow
  Array3DR8 AbsorbedShortRad;  // absorbed solar radiation per unit ground area
                               // per unit incident flux
  Array1DR8 Emissivity;        // emissivity
  Array1DR8 Temperature;       // temperature
  Array1DR8 NetLongRad;        // net longwave radiation
  Array1DR8 UpwardLongRad;     // upward longwave radiation
};

struct WallDataType {
  Array3DR8
      DownwellingShortRad; // downwelling shortwave radiation per unit wall area
  Array3DR8
      ReflectedShortRad; // reflected shortave radiation per unit wall area per
  // unit incident flux
  Array3DR8 BaseAlbedo;       // albedo of wall
  Array3DR8 AbsorbedShortRad; // absorbed shortwave radation per unit wall area
                              // per unit incident flux
  Array1DR8 Emissivity;       // emissivity
  Array1DR8 Temperature;      // temperature
  Array1DR8 NetLongRad;       // net longwave radiation
  Array1DR8 UpwardLongRad;    // upward longwave radiation
};

struct RoofDataType {
  Array3DR8 SnowAlbedo;            // snow albedo
  Array3DR8 AlbedoWithSnowEffects; // albedo of roof including snow effects
  Array3DR8 ReflectedShortRad; // reflected solar radiation per unit ground area
                               // per unit incident flux
  Array3DR8 BaseAlbedo;        // albedo of roof without snow
  Array3DR8 AbsorbedShortRad;  // absorbed solar radiation per unit ground area
                               // per unit incident flux
  Array1DR8 Emissivity;        // emissivity
  Array1DR8 Temperature;       // temperature
  Array1DR8 NetLongRad;        // net longwave radiation
  Array1DR8 UpwardLongRad;     // upward longwave radiation
};

// --- The Master Bundle passed between all classes ---
struct UrbanSharedDataBundle {
  CanyonGeometryData geometry;
  AtmosphereInputData input;

  RoofDataType Roof;
  WallDataType SunlitWall;
  WallDataType ShadedWall;
  RoadDataType ImperviousRoad;
  RoadDataType PerviousRoad;
  CombinedRoadDataType CombinedRoad;

  const int N_LUN;
  const int N_RAD_BAND;
};
} // namespace URBANXX
#endif