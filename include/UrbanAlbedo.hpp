#ifndef URBAN_ALBEDO_HPP
#define URBAN_ALBEDO_HPP

#include <DataTypes.hpp>
#include <Kokkos_Core.hpp>
#include <UrbanData.hpp>

namespace URBANXX {
class NetSolarRoad {
private:
  Array1DR8 Hwr;
  Array1DR8 ViewFactorSkyFromRoad;
  Array1DR8 ViewFactorWallFromRoad;
  RoadDataType RoadData;

public:
  NetSolarRoad(CanyonGeometryData *geometry, RoadDataType &roadData);
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