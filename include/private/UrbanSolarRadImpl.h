#ifndef URBAN_SOLAR_RAD_IMPL_H
#define URBAN_SOLAR_RAD_IMPL_H

#include <Kokkos_Core.hpp>
#include <private/UrbanDataTypesImpl.h>
#include <private/UrbanDataImpl.h>
#include <private/UrbanRadCommonImpl.h>

namespace URBANXX {

class NetSolarRoad {
private:
  Array1DR8 Hwr;
  Array1DR8 ViewFactorSkyFromRoad;
  Array1DR8 ViewFactorWallFromRoad;
  RoadDataType RoadData;
  Real FractionOfTotalRoad;

public:
  NetSolarRoad(CanyonGeometryData *geometry, RoadDataType &roadData, Real);
  void ComputeAbsAndRefRad(RadIndices idx, Real InRad, RadOutput &out) const;
  void ComputeRefRadByComponent(RadIndices idx, Real InRad,
                                RadRefComponents &ref) const;
};

class NetSolarWall {
private:
  Array1DR8 Hwr;
  Array1DR8 ViewFactorSkyFromWall;
  Array1DR8 ViewFactorRoadFromWall;
  Array1DR8 ViewFactorOtherWallFromWall;
  WallDataType WallData;

public:
  NetSolarWall(CanyonGeometryData *geometry, WallDataType &wallData);
  void ComputeAbsAndRefRad(RadIndices idx, Real InRad, RadOutput &out) const;
  void ComputeRefRadByComponent(RadIndices idx, Real InRad,
                                RadRefComponents &ref) const;
};

class UrbanAlbedo {
private:
  UrbanSharedDataBundle &data_bundle;
  NetSolarWall SunlitWallNetSolar;
  NetSolarWall ShadeWallNetSolar;
  NetSolarRoad ImperviousRoadNetSolar;
  NetSolarRoad PerviousRoadNetSolar;

public:
  UrbanAlbedo(UrbanSharedDataBundle &bundle);

  void setSolarInputs();
  void computeIncidentRadiation();
  void computeSnowAlbedo();
  void computeCombinedAlbedo();
  void computeNetSolar();
};
} // namespace URBANXX

#endif // URBAN_SOLAR_RAD_IMPL_H