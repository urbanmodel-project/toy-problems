#ifndef URBAN_ALBEDO_HPP
#define URBAN_ALBEDO_HPP

#include <DataTypes.hpp>
#include <Kokkos_Core.hpp>
#include <UrbanData.hpp>

namespace URBANXX {
struct RadIndices {
  int c;
  int ib;
  int rtype;
};

struct RadOutput {
  Real Abs;
  Real Ref;
};

struct RadRefComponents {
  Real ToSky;
  Real ToRoad;
  Real ToSunwall;
  Real ToShadewall;
  Real ToOtherwall;
};

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

  void set_solar_inputs() const;
  void compute_incident_radiation() const;
  void compute_snow_albedo() const;
  void compute_combined_albedo() const;
  void compute_net_solar() const;
};
} // namespace URBANXX

#endif