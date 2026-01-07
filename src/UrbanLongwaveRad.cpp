#include <Kokkos_Core.hpp>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <private/UrbanDataImpl.h>
#include <private/UrbanLongwaveRadImpl.h>
#include <private/UrbanRadCommonImpl.h>
#include <stdexcept>
#include <string>

#define STEBOL 5.670374419e-8
#define PERVIOUS_ROAD_FRACTION 0.16666667163372040

namespace URBANXX {
NetLongwaveRoad::NetLongwaveRoad(CanyonGeometryData *geometry,
                                 RoadDataType &roadData, Real roadWeight)
    : Hwr(geometry->CanyonHwr),
      ViewFactorSkyFromRoad(geometry->ViewFactorSkyFromRoad),
      ViewFactorWallFromRoad(geometry->ViewFactorWallFromRoad),
      RoadData(roadData), FractionOfTotalRoad(roadWeight) {}

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
                          1.0 - PERVIOUS_ROAD_FRACTION),
      PerviousRoadLwave(&bundle.geometry, bundle.PerviousRoad,
                        PERVIOUS_ROAD_FRACTION) {}

KOKKOS_FUNCTION
void NetLongwaveRoad::ComputeAbsAndRefRad(RadIndices idx, Real InRad,
                                          RadOutput &out) const {

  Real emiss = RoadData.Emissivity(idx.c);
  Real Temp = RoadData.Temperature(idx.c);

  out.Abs = emiss * InRad;
  out.Ref = (1.0 - emiss) * InRad;
  out.Emi = emiss * STEBOL * std::pow(Temp, 4.0);

  out.AbsByWt = out.Abs * FractionOfTotalRoad;
  out.RefByWt = out.Ref * FractionOfTotalRoad;
  out.EmiByWt = out.Emi * FractionOfTotalRoad;
}

KOKKOS_FUNCTION
void NetLongwaveRoad::ComputeRadByComponent(RadIndices idx, Real Data,
                                            RadRefComponents &ref) const {

  Real vf_sr = ViewFactorSkyFromRoad(idx.c);
  Real vf_wr = ViewFactorWallFromRoad(idx.c);

  ref.ToSky = Data * vf_sr;
  ref.ToSunwall = Data * vf_wr;
  ref.ToShadewall = Data * vf_wr;
  ref.ToOtherwall = 0.0;
  ref.ToRoad = 0.0;

  ref.ToSkyByWt = ref.ToSky * FractionOfTotalRoad;
  ref.ToSunwallByWt = ref.ToSunwall * FractionOfTotalRoad;
  ref.ToShadewallByWt = ref.ToShadewall * FractionOfTotalRoad;
  ref.ToOtherwallByWt = ref.ToOtherwall * FractionOfTotalRoad;
  ref.ToRoadByWt = ref.ToRoad * FractionOfTotalRoad;
}

KOKKOS_FUNCTION
void NetLongwaveRoad::ComputeRefRadByComponent(RadIndices idx, Real InRad,
                                               RadRefComponents &ref) const {

  Real emiss = RoadData.Emissivity(idx.c);
  Real Data = (1.0 - emiss) * InRad;

  ComputeRadByComponent(idx, Data, ref);
}

KOKKOS_FUNCTION
void NetLongwaveRoad::ComputeEmiRadByComponent(RadIndices idx, Real InRad,
                                               RadRefComponents &ref) const {

  Real emiss = RoadData.Emissivity(idx.c);
  Real Temp = RoadData.Temperature(idx.c);
  Real Data = emiss * STEBOL * std::pow(Temp, 4.0);

  ComputeRadByComponent(idx, Data, ref);
}

KOKKOS_FUNCTION
void NetLongwaveWall::ComputeAbsAndRefRad(RadIndices idx, Real InRad,
                                          RadOutput &out) const {

  Real emiss = WallData.Emissivity(idx.c);
  Real Temp = WallData.Temperature(idx.c);

  out.Abs = emiss * InRad;
  out.Ref = (1.0 - emiss) * InRad;
  out.Emi = emiss * STEBOL * std::pow(Temp, 4.0);

  out.AbsByWt = out.Abs;
  out.RefByWt = out.Ref;
  out.EmiByWt = out.Emi;
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

  ref.ToSkyByWt = ref.ToSky;
  ref.ToSunwallByWt = ref.ToSunwall;
  ref.ToShadewallByWt = ref.ToShadewall;
  ref.ToOtherwallByWt = ref.ToOtherwall;
  ref.ToRoadByWt = ref.ToRoad;
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
  Real Data = emiss * STEBOL * std::pow(Temp, 4.0);

  ComputeRadByComponent(idx, Data, ref);
}

void UrbanLongwave::setLongwaveInputs() {

  std::cout << "Setting longwave inputs\n";

  int N_LUN = data_bundle.N_LUN;
  HostArray1DR8 DwLongH = data_bundle.input.DownwellingLongRadH;
  Array1DR8 DwLong = data_bundle.input.DownwellingLongRad;

  const Real lwdown = 432.79580327766450;

  for (int i = 0; i < N_LUN; i++) {
    DwLongH(i) = lwdown;
  }

  Kokkos::deep_copy(DwLong, DwLongH);
}

void UrbanLongwave::computeNetLongwave() {
  int N_LUN = data_bundle.N_LUN;
  auto & vf_sr_view = data_bundle.geometry.ViewFactorSkyFromRoad;
  auto & vf_sw_view = data_bundle.geometry.ViewFactorSkyFromWall;
  auto & hwr_view = data_bundle.geometry.CanyonHwr;
  auto & dwlong_view = data_bundle.input.DownwellingLongRad;

  if (vf_sr_view.extent(0) == 0 || vf_sw_view.extent(0) == 0 ||
      hwr_view.extent(0) == 0 || dwlong_view.extent(0) == 0) {
    throw std::runtime_error(
        std::string("UrbanLongwave::computeNetLongwave: one or more views have "
                    "zero extent:") +
        " vf_sr=" + std::to_string(vf_sr_view.extent(0)) +
        " vf_sw=" + std::to_string(vf_sw_view.extent(0)) +
        " hwr=" + std::to_string(hwr_view.extent(0)) +
        " dwlong=" + std::to_string(dwlong_view.extent(0)));
  }
  Kokkos::parallel_for(
      "ComputeNetLongwave", N_LUN, KOKKOS_LAMBDA(const int c) {
        Array1DR8 vf_sr = data_bundle.geometry.ViewFactorSkyFromRoad;
        Array1DR8 vf_sw = data_bundle.geometry.ViewFactorSkyFromWall;
        Array1DR8 hwr = data_bundle.geometry.CanyonHwr;

        Array1DR8 DwLong = data_bundle.input.DownwellingLongRad;

        Array1DR8 NetImproad = data_bundle.ImperviousRoad.NetLongRad;
        Array1DR8 NetPerroad = data_bundle.PerviousRoad.NetLongRad;
        Array1DR8 NetSunwall = data_bundle.SunlitWall.NetLongRad;
        Array1DR8 NetShadewall = data_bundle.ShadedWall.NetLongRad;

        Array1DR8 UpImproad = data_bundle.ImperviousRoad.UpwardLongRad;
        Array1DR8 UpPerroad = data_bundle.PerviousRoad.UpwardLongRad;
        Array1DR8 UpSunwall = data_bundle.SunlitWall.UpwardLongRad;
        Array1DR8 UpShadewall = data_bundle.ShadedWall.UpwardLongRad;

        RadIndices idx{c, 0, 0};

        Real LtotForRoad = DwLong(c) * vf_sr(c);

        // +++++++++++++++++++++++++++++++++++++++++++++++++
        // Impervious road
        // +++++++++++++++++++++++++++++++++++++++++++++++++
        RadOutput ImpRoadOut;
        RadRefComponents ImpRoadRef, ImpRoadEmi;

        // Impervious road: absorbed and reflected
        ImperviousRoadLwave.ComputeAbsAndRefRad(idx, LtotForRoad, ImpRoadOut);
        // Impervious road: reflected for each components
        ImperviousRoadLwave.ComputeRefRadByComponent(idx, LtotForRoad,
                                                     ImpRoadRef);

        // Impervious road: emitted for each components
        ImperviousRoadLwave.ComputeEmiRadByComponent(idx, LtotForRoad,
                                                     ImpRoadEmi);

        // +++++++++++++++++++++++++++++++++++++++++++++++++
        // Pervious road
        // +++++++++++++++++++++++++++++++++++++++++++++++++
        RadOutput PerRoadOut, PerRoadOutByWt;
        RadRefComponents PerRoadRef, PerRoadRefByWt;
        RadRefComponents PerRoadEmi, PerRoadEmiByWt;

        // Impervious road: absorbed and reflected
        PerviousRoadLwave.ComputeAbsAndRefRad(idx, LtotForRoad, PerRoadOut);

        // Impervious road: reflected for each components
        PerviousRoadLwave.ComputeRefRadByComponent(idx, LtotForRoad,
                                                   PerRoadRef);

        // Impervious road: emitted for each components
        PerviousRoadLwave.ComputeEmiRadByComponent(idx, LtotForRoad,
                                                   PerRoadEmi);

        // +++++++++++++++++++++++++++++++++++++++++++++++++
        // Combining data from impervious and pervious road
        // +++++++++++++++++++++++++++++++++++++++++++++++++

        // absorbed, reflected, and emitted radation
        Real RoadAbs, RoadRef, RoadEmi;
        RoadAbs = ImpRoadOut.AbsByWt + PerRoadOut.AbsByWt;
        RoadRef = ImpRoadOut.RefByWt + PerRoadOut.RefByWt;
        RoadEmi = ImpRoadOut.EmiByWt + PerRoadOut.EmiByWt;

        // reflected radation to sky, sunwall, and shadewall
        Real RoadRefToSky, RoadRefToSunwall, RoadRefToShadewall;
        RoadRefToSky = ImpRoadRef.ToSkyByWt + PerRoadRef.ToSkyByWt;
        RoadRefToSunwall = ImpRoadRef.ToSunwallByWt + PerRoadRef.ToSunwallByWt;
        RoadRefToShadewall =
            ImpRoadRef.ToShadewallByWt + PerRoadRef.ToShadewallByWt;

        // emitted radiation to sky, sunwall, and shadewall
        Real RoadEmiToSky, RoadEmiToSunwall, RoadEmiToShadewall;
        RoadEmiToSky = ImpRoadEmi.ToSkyByWt + PerRoadEmi.ToSkyByWt;
        RoadEmiToSunwall = ImpRoadEmi.ToSunwallByWt + PerRoadEmi.ToSunwallByWt;
        RoadEmiToShadewall =
            ImpRoadEmi.ToShadewallByWt + PerRoadEmi.ToShadewallByWt;

        // +++++++++++++++++++++++++++++++++++++++++++++++++
        // Sunlit wall
        // +++++++++++++++++++++++++++++++++++++++++++++++++
        Real LtotForSunwall = DwLong(c) * vf_sw(c);
        RadOutput SunwallOut;
        RadRefComponents SunwallRef, SunwallEmi;

        SunlitWallLwave.ComputeAbsAndRefRad(idx, LtotForSunwall, SunwallOut);
        SunlitWallLwave.ComputeRefRadByComponent(idx, LtotForSunwall,
                                                 SunwallRef);
        SunlitWallLwave.ComputeEmiRadByComponent(idx, LtotForSunwall,
                                                 SunwallEmi);

        // +++++++++++++++++++++++++++++++++++++++++++++++++
        // Shaded wall
        // +++++++++++++++++++++++++++++++++++++++++++++++++
        Real LtotForShadewall = DwLong(c) * vf_sw(c);

        RadOutput ShadewallOut;
        RadRefComponents ShadewallRef, ShadewallEmi;

        ShadeWallLwave.ComputeAbsAndRefRad(idx, LtotForSunwall, ShadewallOut);
        ShadeWallLwave.ComputeRefRadByComponent(idx, LtotForShadewall,
                                                ShadewallRef);
        ShadeWallLwave.ComputeEmiRadByComponent(idx, LtotForSunwall,
                                                ShadewallEmi);

        // Cummulative absorbed and reflected radiation for the four
        // surfaces

        NetImproad(c) = ImpRoadOut.Emi - ImpRoadOut.Abs;
        NetPerroad(c) = PerRoadOut.Emi - PerRoadOut.Abs;
        NetSunwall(c) = SunwallOut.Emi - SunwallOut.Abs;
        NetShadewall(c) = ShadewallOut.Emi - ShadewallOut.Abs;

        UpImproad(c) = ImpRoadRef.ToSky + ImpRoadEmi.ToSky;
        UpPerroad(c) = PerRoadRef.ToSky + PerRoadEmi.ToSky;
        UpSunwall(c) = SunwallRef.ToSky + SunwallEmi.ToSky;
        UpShadewall(c) = ShadewallRef.ToSky + ShadewallEmi.ToSky;

        const int max_iter = 50;
        for (int iter = 0; iter < max_iter; iter++) {

          // step(1)
          //  For roads
          LtotForRoad = (SunwallRef.ToRoad + SunwallEmi.ToRoad +
                         ShadewallRef.ToRoad + ShadewallEmi.ToRoad) *
                        hwr(c);

          ImperviousRoadLwave.ComputeAbsAndRefRad(idx, LtotForRoad, ImpRoadOut);
          PerviousRoadLwave.ComputeAbsAndRefRad(idx, LtotForRoad, PerRoadOut);

          RoadAbs = ImpRoadOut.AbsByWt + PerRoadOut.AbsByWt;
          RoadRef = ImpRoadOut.RefByWt + PerRoadOut.RefByWt;

          //  For sunwall
          LtotForSunwall = (RoadRefToSunwall + RoadEmiToSunwall) / hwr(c) +
                           ShadewallRef.ToOtherwall + ShadewallEmi.ToOtherwall;
          SunlitWallLwave.ComputeAbsAndRefRad(idx, LtotForSunwall, SunwallOut);

          //  For shadwall
          LtotForShadewall =
              (RoadRefToShadewall + RoadEmiToShadewall) / hwr(c) +
              SunwallRef.ToOtherwall + SunwallEmi.ToOtherwall;
          ShadeWallLwave.ComputeAbsAndRefRad(idx, LtotForSunwall, ShadewallOut);

          //  set emitted values to zero, so they are not counted multiple times
          //  within the iteration loop
          SunwallEmi.ToRoad = 0.0;
          SunwallEmi.ToOtherwall = 0.0;
          ShadewallEmi.ToRoad = 0.0;
          ShadewallEmi.ToOtherwall = 0.0;
          RoadEmiToSunwall = 0.0;
          RoadEmiToShadewall = 0.0;

          // step (2)
          NetImproad(c) -= ImpRoadOut.Abs;
          NetPerroad(c) -= PerRoadOut.Abs;
          NetSunwall(c) -= SunwallOut.Abs;
          NetShadewall(c) -= ShadewallOut.Abs;

          // step (3)
          ImperviousRoadLwave.ComputeRefRadByComponent(idx, LtotForRoad,
                                                       ImpRoadRef);
          PerviousRoadLwave.ComputeRefRadByComponent(idx, LtotForRoad,
                                                     PerRoadRef);

          RoadRefToSky = ImpRoadRef.ToSkyByWt + PerRoadRef.ToSkyByWt;
          RoadRefToSunwall =
              ImpRoadRef.ToSunwallByWt + PerRoadRef.ToSunwallByWt;
          RoadRefToShadewall =
              ImpRoadRef.ToShadewallByWt + PerRoadRef.ToShadewallByWt;

          SunlitWallLwave.ComputeRefRadByComponent(idx, LtotForSunwall,
                                                   SunwallRef);
          ShadeWallLwave.ComputeRefRadByComponent(idx, LtotForShadewall,
                                                  ShadewallRef);

          // step (4)
          UpImproad(c) += ImpRoadRef.ToSky;
          UpPerroad(c) += PerRoadRef.ToSky;
          UpSunwall(c) += SunwallRef.ToSky;
          UpShadewall(c) += ShadewallRef.ToSky;

          Real crit = std::max({RoadAbs, SunwallOut.Abs, ShadewallOut.Abs});
          const Real errcrit = 0.001;
          if (crit < errcrit)
            break;
        }
      });
}

} // namespace URBANXX
