#ifndef URBAN_LONGWAVE_RAD_HPP
#define URBAN_LONGWAVE_RAD_HPP

#include <private/DataTypes.hpp>
#include <Kokkos_Core.hpp>
#include <private/UrbanData.hpp>
#include <private/UrbanRadCommon.hpp>

namespace URBANXX {

class NetLongwaveRoad {
private:
  Array1DR8 Hwr;
  Array1DR8 ViewFactorSkyFromRoad;
  Array1DR8 ViewFactorWallFromRoad;
  RoadDataType RoadData;
  Real FractionOfTotalRoad;

public:
  NetLongwaveRoad(CanyonGeometryData *geometry, RoadDataType &roadData, Real);
  void ComputeAbsAndRefRad(RadIndices idx, Real InRad, RadOutput &out) const;
  void ComputeRadByComponent(RadIndices idx, Real InRad,
                             RadRefComponents &ref) const;
  void ComputeRefRadByComponent(RadIndices idx, Real InRad,
                                RadRefComponents &ref) const;
  void ComputeEmiRadByComponent(RadIndices idx, Real InRad,
                                RadRefComponents &ref) const;
};

class NetLongwaveWall {
private:
  Array1DR8 Hwr;
  Array1DR8 ViewFactorSkyFromWall;
  Array1DR8 ViewFactorRoadFromWall;
  Array1DR8 ViewFactorOtherWallFromWall;
  WallDataType WallData;

public:
  NetLongwaveWall(CanyonGeometryData *geometry, WallDataType &wallData);
  void ComputeAbsAndRefRad(RadIndices idx, Real InRad, RadOutput &out) const;
  void ComputeRadByComponent(RadIndices idx, Real InRad,
                             RadRefComponents &ref) const;
  void ComputeRefRadByComponent(RadIndices idx, Real InRad,
                                RadRefComponents &ref) const;
  void ComputeEmiRadByComponent(RadIndices idx, Real InRad,
                                RadRefComponents &ref) const;
};

class UrbanLongwave {
private:
  UrbanSharedDataBundle &data_bundle;
  NetLongwaveWall SunlitWallLwave;
  NetLongwaveWall ShadeWallLwave;
  NetLongwaveRoad ImperviousRoadLwave;
  NetLongwaveRoad PerviousRoadLwave;

public:
  UrbanLongwave(UrbanSharedDataBundle &bundle);
  void setLongwaveInputs();
  void computeNetLongwave();
};
} // namespace URBANXX

#endif