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

  Array1DR8 ViewFactorSkyFromRoad;  // view factor of sky for road
  Array1DR8 ViewFactorSkyFromWall;  // view factor of sky for one wall
  Array1DR8 ViewFactorRoadFromWall; // view factor from road to wall
  Array1DR8
      ViewFactorOtherWallFromWall;  // view factor from wall to opposing wall
  Array1DR8 ViewFactorWallFromRoad; // view factor of one wall to road
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
};

// --- The Master Bundle passed between all classes ---
struct UrbanSharedDataBundle {
  CanyonGeometryData geometry;
  SolarInputData input;

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