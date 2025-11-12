#ifndef URBAN_LONGWAVE_RAD_HPP
#define URBAN_LONGWAVE_RAD_HPP

#include <DataTypes.hpp>
#include <Kokkos_Core.hpp>
#include <UrbanData.hpp>
#include <UrbanRadCommon.hpp>

namespace URBANXX {

class NetLongwaveRoad {
private:
  Array1DR8 Hwr;
  Array1DR8 ViewFactorSkyFromRoad;
  Array1DR8 ViewFactorWallFromRoad;
  RoadDataType RoadData;
  Real Weight;

public:
  NetLongwaveRoad(CanyonGeometryData *geometry, RoadDataType &roadData, Real);
  void ComputeAbsAndRefRad(RadIndices idx, Real InRad, RadOutput &out,
                           bool scale_by_weight) const;
  void ComputeRadByComponent(RadIndices idx, Real InRad, RadRefComponents &ref,
                             bool scale_by_weight) const;
  void ComputeRefRadByComponent(RadIndices idx, Real InRad,
                                RadRefComponents &ref,
                                bool scale_by_weight) const;
  void ComputeEmiRadByComponent(RadIndices idx, Real InRad,
                                RadRefComponents &ref,
                                bool scale_by_weight) const;
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
};
} // namespace URBANXX

#endif