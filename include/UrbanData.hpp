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

template <typename ArrayType> struct RadiationType {
  ArrayType dir;
  ArrayType dif;
};

struct CombinedRoadType {
  RadiationType<Array2DR8>
      DownwellingShortRad; // downwelling shortwave radiation per unit road area
  RadiationType<Array2DR8> SnowAlbedo; // snow albedo
  RadiationType<Array2DR8>
      AlbedoWithSnowEffects; // albedo of road including snow effects
};

struct RoadType {
  RadiationType<Array2DR8> SnowAlbedo; // snow albedo
  RadiationType<Array2DR8>
      AlbedoWithSnowEffects; // albedo of road including snow effects
  RadiationType<Array2DR8>
      ReflectedShortRad; // reflected solar radiation per unit ground area per
                         // unit incident flux
  RadiationType<Array2DR8>
      AbsorbedShortRad; // absorbed solar radiation per unit ground area per
                        // unit incident flux
};

struct WallType {
  RadiationType<Array2DR8>
      DownwellingShortRad; // downwelling shortwave radiation per unit wall area
  RadiationType<Array2DR8>
      ReflectedShortRad; // reflected shortave radiation per unit wall area per
  // unit incident flux
  RadiationType<Array2DR8>
      AbsorbedShortRad; // absorbed shortwave radation per unit wall area per
                        // unit incident flux
};

struct RoofType {
  RadiationType<Array2DR8> SnowAlbedo; // snow albedo
  RadiationType<Array2DR8>
      AlbedoWithSnowEffects; // albedo of roof including snow effects
  RadiationType<Array2DR8>
      ReflectedShortRad; // reflected solar radiation per unit ground area per
                         // unit incident flux
  RadiationType<Array2DR8>
      AbsorbedShortRad; // absorbed solar radiation per unit ground area per
                        // unit incident flux
};

// --- The Master Bundle passed between all classes ---
struct UrbanSharedDataBundle {
  CanyonGeometryData geometry;
  SolarInputData input;

  RoofType Roof;
  WallType SunlitWall;
  WallType ShadedWall;
  RoadType ImperviousRoad;
  RoadType PerviousRoad;
  CombinedRoadType CombinedRoad;

  const int N_LUN;
  const int N_RAD;
};
} // namespace URBANXX
#endif