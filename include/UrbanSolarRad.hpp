#ifndef URBAN_SOLAR_RAD_HPP
#define URBAN_SOLAR_RAD_HPP

#include <DataTypes.hpp>
#include <Kokkos_Core.hpp>
#include <UrbanData.hpp>
#include <UrbanRadCommon.hpp>

namespace URBANXX {

class NetSolarRoad {
private:
  Array1DR8 Hwr;
  Array1DR8 ViewFactorSkyFromRoad;
  Array1DR8 ViewFactorWallFromRoad;
  RoadDataType RoadData;
  Real Weight;

public:
  NetSolarRoad(CanyonGeometryData *geometry, RoadDataType &roadData, Real);
  void ComputeAbsAndRefRad(RadIndices idx, Real InRad, RadOutput &out,
                           bool scale_by_weight) const;
  void ComputeRefRadByComponent(RadIndices idx, Real InRad,
                                RadRefComponents &ref,
                                bool scale_by_weight) const;
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

#endif