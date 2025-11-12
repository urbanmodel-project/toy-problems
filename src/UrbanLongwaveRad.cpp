#include <Kokkos_Core.hpp>
#include <UrbanData.hpp>
#include <UrbanLongwaveRad.hpp>
#include <UrbanRadCommon.hpp>
#include <cmath>
#include <iomanip>
#include <iostream>

namespace URBANXX {
NetLongwaveRoad::NetLongwaveRoad(CanyonGeometryData *geometry, RoadDataType &roadData,
                           Real roadWeight)
    : Hwr(geometry->CanyonHwr),
      ViewFactorSkyFromRoad(geometry->ViewFactorSkyFromRoad),
      ViewFactorWallFromRoad(geometry->ViewFactorWallFromRoad),
      RoadData(roadData), Weight(roadWeight) {}

NetLongwaveWall::NetLongwaveWall(CanyonGeometryData *geometry, WallDataType &wallData)
    : Hwr(geometry->CanyonHwr),
      ViewFactorSkyFromWall(geometry->ViewFactorSkyFromWall),
      ViewFactorRoadFromWall(geometry->ViewFactorRoadFromWall),
      ViewFactorOtherWallFromWall(geometry->ViewFactorOtherWallFromWall),
      WallData(wallData) {}

UrbanLongwave::UrbanLongwave(UrbanSharedDataBundle &bundle)
    : data_bundle(bundle),
      SunlitWallLwave(&bundle.geometry, bundle.SunlitWall),
      ShadeWallLwave(&bundle.geometry, bundle.ShadedWall),
      ImperviousRoadLwave(&bundle.geometry, bundle.ImperviousRoad,
                             1.0 - 0.16666667163372040),
      PerviousRoadLwave(&bundle.geometry, bundle.PerviousRoad,
                           0.16666667163372040) {}
}