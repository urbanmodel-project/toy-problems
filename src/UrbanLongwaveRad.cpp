#include <Kokkos_Core.hpp>
#include <UrbanData.hpp>
#include <UrbanLongwaveRad.hpp>
#include <UrbanRadCommon.hpp>
#include <cmath>
#include <iomanip>
#include <iostream>

namespace URBANXX {
NetLongwaveRoad::NetLongwaveRoad(CanyonGeometryData *geometry,
                                 RoadDataType &roadData, Real roadWeight)
    : Hwr(geometry->CanyonHwr),
      ViewFactorSkyFromRoad(geometry->ViewFactorSkyFromRoad),
      ViewFactorWallFromRoad(geometry->ViewFactorWallFromRoad),
      RoadData(roadData), Weight(roadWeight) {}

NetLongwaveWall::NetLongwaveWall(CanyonGeometryData *geometry,
                                 WallDataType &wallData)
    : Hwr(geometry->CanyonHwr),
      ViewFactorSkyFromWall(geometry->ViewFactorSkyFromWall),
      ViewFactorRoadFromWall(geometry->ViewFactorRoadFromWall),
      ViewFactorOtherWallFromWall(geometry->ViewFactorOtherWallFromWall),
      WallData(wallData) {}

UrbanLongwave::UrbanLongwave(UrbanSharedDataBundle &bundle)
    : data_bundle(bundle), SunlitWallLwave(&bundle.geometry, bundle.SunlitWall),
      ShadeWallLwave(&bundle.geometry, bundle.ShadedWall),
      ImperviousRoadLwave(&bundle.geometry, bundle.ImperviousRoad,
                          1.0 - 0.16666667163372040),
      PerviousRoadLwave(&bundle.geometry, bundle.PerviousRoad,
                        0.16666667163372040) {}

KOKKOS_FUNCTION
void NetLongwaveRoad::ComputeAbsAndRefRad(RadIndices idx, Real InRad,
                                          RadOutput &out,
                                          bool scale_by_weight) const {

  Real emiss = RoadData.Emissivity(idx.c);
  Real Temp = RoadData.Temperature(idx.c);

  out.Abs = emiss * InRad;
  out.Ref = (1.0 - emiss) * InRad;
  out.Emi = emiss * std::pow(Temp, 4.0);
  if (scale_by_weight) {
    out.Abs *= Weight;
    out.Ref *= Weight;
    out.Emi *= Weight;
  }
}

KOKKOS_FUNCTION
void NetLongwaveRoad::ComputeRadByComponent(RadIndices idx, Real Data,
                                            RadRefComponents &ref,
                                            bool scale_by_weight) const {

  Real vf_sr = ViewFactorSkyFromRoad(idx.c);
  Real vf_wr = ViewFactorWallFromRoad(idx.c);

  ref.ToSky = Data * vf_sr;
  ref.ToSunwall = Data * vf_wr;
  ref.ToShadewall = Data * vf_wr;
  ref.ToOtherwall = 0.0;
  ref.ToRoad = 0.0;

  if (scale_by_weight) {
    ref.ToSky *= Weight;
    ref.ToSunwall *= Weight;
    ref.ToShadewall *= Weight;
  }
}

KOKKOS_FUNCTION
void NetLongwaveRoad::ComputeRefRadByComponent(RadIndices idx, Real InRad,
                                               RadRefComponents &ref,
                                               bool scale_by_weight) const {

  Real emiss = RoadData.Emissivity(idx.c);
  Real Data = (1.0 - emiss) * InRad;

  ComputeRadByComponent(idx, Data, ref, scale_by_weight);
}

KOKKOS_FUNCTION
void NetLongwaveRoad::ComputeEmiRadByComponent(RadIndices idx, Real InRad,
                                               RadRefComponents &ref,
                                               bool scale_by_weight) const {

  Real emiss = RoadData.Emissivity(idx.c);
  Real Temp = RoadData.Temperature(idx.c);
  Real Data = emiss * std::pow(Temp, 4.0);

  ComputeRadByComponent(idx, Data, ref, scale_by_weight);
}

KOKKOS_FUNCTION
void NetLongwaveWall::ComputeAbsAndRefRad(RadIndices idx, Real InRad,
                                          RadOutput &out) const {

  Real emiss = WallData.Emissivity(idx.c);
  Real Temp = WallData.Temperature(idx.c);

  out.Abs = emiss * InRad;
  out.Ref = (1.0 - emiss) * InRad;
  out.Emi = emiss * std::pow(Temp, 4.0);
}

KOKKOS_FUNCTION
void NetLongwaveWall::ComputeRadByComponent(RadIndices idx, Real Data,
                                            RadRefComponents &ref) const {

  Real vf_sw = ViewFactorSkyFromWall(idx.c);
  Real vf_rw = ViewFactorRoadFromWall(idx.c);
  Real vf_ww = ViewFactorOtherWallFromWall(idx.c);

  ref.ToSky = Data * vf_sw;
  ref.ToRoad = Data * vf_rw;
  ref.ToOtherwall = Data * vf_ww;
  ref.ToSunwall = 0.0;
  ref.ToShadewall = 0.0;
}

KOKKOS_FUNCTION
void NetLongwaveWall::ComputeRefRadByComponent(RadIndices idx, Real InRad,
                                               RadRefComponents &ref) const {

  Real emiss = WallData.Emissivity(idx.c);
  Real Data = (1.0 - emiss) * InRad;

  ComputeRadByComponent(idx, Data, ref);
}

KOKKOS_FUNCTION
void NetLongwaveWall::ComputeEmiRadByComponent(RadIndices idx, Real InRad,
                                               RadRefComponents &ref) const {

  Real emiss = WallData.Emissivity(idx.c);
  Real Temp = WallData.Temperature(idx.c);
  Real Data = emiss * std::pow(Temp, 4.0);

  ComputeRadByComponent(idx, Data, ref);
}

} // namespace URBANXX
