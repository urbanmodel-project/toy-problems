#include <Kokkos_Core.hpp>
#include <UrbanData.hpp>
#include <UrbanRadCommon.hpp>
#include <UrbanSolarRad.hpp>
#include <cmath>
#include <iomanip>
#include <iostream>

namespace URBANXX {

NetSolarRoad::NetSolarRoad(CanyonGeometryData *geometry, RoadDataType &roadData,
                           Real roadWeight)
    : Hwr(geometry->CanyonHwr),
      ViewFactorSkyFromRoad(geometry->ViewFactorSkyFromRoad),
      ViewFactorWallFromRoad(geometry->ViewFactorWallFromRoad),
      RoadData(roadData), Weight(roadWeight) {}

NetSolarWall::NetSolarWall(CanyonGeometryData *geometry, WallDataType &wallData)
    : Hwr(geometry->CanyonHwr),
      ViewFactorSkyFromWall(geometry->ViewFactorSkyFromWall),
      ViewFactorRoadFromWall(geometry->ViewFactorRoadFromWall),
      ViewFactorOtherWallFromWall(geometry->ViewFactorOtherWallFromWall),
      WallData(wallData) {}

UrbanAlbedo::UrbanAlbedo(UrbanSharedDataBundle &bundle)
    : data_bundle(bundle),
      SunlitWallNetSolar(&bundle.geometry, bundle.SunlitWall),
      ShadeWallNetSolar(&bundle.geometry, bundle.ShadedWall),
      ImperviousRoadNetSolar(&bundle.geometry, bundle.ImperviousRoad,
                             1.0 - 0.16666667163372040),
      PerviousRoadNetSolar(&bundle.geometry, bundle.PerviousRoad,
                           0.16666667163372040) {}

KOKKOS_FUNCTION
void NetSolarRoad::ComputeAbsAndRefRad(RadIndices idx, Real InRad,
                                       RadOutput &out) const {

  Real albedo = RoadData.BaseAlbedo(idx.c, idx.ib, idx.rtype);

  out.Abs = (1.0 - albedo) * InRad;
  out.Ref = albedo * InRad;
  out.Emi = 0.0;

  out.AbsByWt = out.Abs * Weight;
  out.RefByWt = out.Ref * Weight;
  out.EmiByWt = out.Emi * Weight;
}

KOKKOS_FUNCTION
void NetSolarRoad::ComputeRefRadByComponent(RadIndices idx, Real InRad,
                                            RadRefComponents &ref) const {

  Real albedo = RoadData.BaseAlbedo(idx.c, idx.ib, idx.rtype);
  Real vf_sr = ViewFactorSkyFromRoad(idx.c);
  Real vf_wr = ViewFactorWallFromRoad(idx.c);

  Real RefRad = albedo * InRad;

  ref.ToSky = RefRad * vf_sr;
  ref.ToSunwall = RefRad * vf_wr;
  ref.ToShadewall = RefRad * vf_wr;
  ref.ToOtherwall = 0.0;
  ref.ToRoad = 0.0;

  ref.ToSkyByWt = ref.ToSky * Weight;
  ref.ToSunwallByWt = ref.ToSunwall * Weight;
  ref.ToShadewallByWt = ref.ToShadewall * Weight;
  ref.ToOtherwallByWt = ref.ToOtherwall * Weight;
  ref.ToRoadByWt = ref.ToRoad * Weight;
}

KOKKOS_FUNCTION
void NetSolarWall::ComputeAbsAndRefRad(RadIndices idx, Real InRad,
                                       RadOutput &out) const {

  Real albedo = WallData.BaseAlbedo(idx.c, idx.ib, idx.rtype);

  out.Abs = (1.0 - albedo) * InRad;
  out.Ref = albedo * InRad;
  out.Emi = 0.0;

  out.AbsByWt = out.Abs;
  out.RefByWt = out.Ref;
  out.EmiByWt = out.Emi;
}

KOKKOS_FUNCTION
void NetSolarWall::ComputeRefRadByComponent(RadIndices idx, Real InRad,
                                            RadRefComponents &ref) const {

  Real albedo = WallData.BaseAlbedo(idx.c, idx.ib, idx.rtype);
  Real vf_sw = ViewFactorSkyFromWall(idx.c);
  Real vf_rw = ViewFactorRoadFromWall(idx.c);
  Real vf_ww = ViewFactorOtherWallFromWall(idx.c);

  Real RefRad = albedo * InRad;

  ref.ToSky = RefRad * vf_sw;
  ref.ToRoad = RefRad * vf_rw;
  ref.ToOtherwall = RefRad * vf_ww;
  ref.ToSunwall = 0.0;
  ref.ToShadewall = 0.0;

  ref.ToSkyByWt = ref.ToSky;
  ref.ToSunwallByWt = ref.ToSunwall;
  ref.ToShadewallByWt = ref.ToShadewall;
  ref.ToOtherwallByWt = ref.ToOtherwall;
  ref.ToRoadByWt = ref.ToRoad;
}

template <typename ViewType>
void print_view_1d(const ViewType &view, const std::string &name = "") {
  auto h_view = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), view);
  std::cout << name << " = [";
  for (std::size_t i = 0; i < h_view.extent(0); ++i) {
    std::cout << std::setprecision(15) << h_view(i);
    if (i + 1 < h_view.extent(0))
      std::cout << ", ";
  }
  std::cout << "]\n";
}

template <typename ViewType>
void print_view_2d(const ViewType &view, const std::string &name = "") {
  auto h_view = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), view);
  std::cout << name << " (" << h_view.extent(0) << "x" << h_view.extent(1)
            << "):\n";
  for (std::size_t i = 0; i < h_view.extent(0); ++i) {
    for (std::size_t j = 0; j < h_view.extent(1); ++j)
      std::cout << std::setprecision(15) << h_view(i, j) << " ";
    std::cout << "\n";
  }
}

template <typename ViewType>
void print_view_3d(const ViewType &view, const std::string &name = "") {
  auto h_view = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), view);
  std::cout << name << " (" << h_view.extent(0) << "x" << h_view.extent(1)
            << "x" << h_view.extent(2) << "):\n";
  for (std::size_t i = 0; i < h_view.extent(0); ++i) {
    for (std::size_t j = 0; j < h_view.extent(1); ++j) {
      for (std::size_t k = 0; k < h_view.extent(2); ++k)
        std::cout << std::setprecision(15) << h_view(i, j, k) << " ";
      std::cout << "\n";
    }
    std::cout << "\n";
  }
}

void UrbanAlbedo::setSolarInputs() {

  std::cout << "Setting solar inputs\n";

  int N_LUN = data_bundle.N_LUN;
  HostArray1DR8 CoszenH = data_bundle.input.CoszenH;
  Array1DR8 Coszen = data_bundle.input.Coszen;

  const Real coszen = 7.9054122593736065E-003;

  for (int i = 0; i < N_LUN; i++) {
    CoszenH(i) = coszen;
  }

  Kokkos::deep_copy(Coszen, CoszenH);

  if (0)
    print_view_1d(Coszen);
}

void UrbanAlbedo::computeIncidentRadiation() {
  std::cout << "In compute_incident_radiation\n";

  const Real rpi = M_PI;

  int N_LUN = data_bundle.N_LUN;
  int N_RAD_BAND = data_bundle.N_RAD_BAND;

  Kokkos::parallel_for(
      "ComputeIncidentRadiation", N_LUN, KOKKOS_LAMBDA(const int c) {
        const Real coszen = data_bundle.input.Coszen(c);
        const Real hwr = data_bundle.geometry.CanyonHwr(c);

        Array3DR8 s_sunlitwall = data_bundle.SunlitWall.DownwellingShortRad;
        Array3DR8 s_shadewall = data_bundle.ShadedWall.DownwellingShortRad;
        Array3DR8 s_road = data_bundle.CombinedRoad.DownwellingShortRad;

        Array1DR8 vf_sr = data_bundle.geometry.ViewFactorSkyFromRoad;
        Array1DR8 vf_sw = data_bundle.geometry.ViewFactorSkyFromWall;

        // the incident direct and diffuse radiation for VIS and NIR bands is
        // assumed to be unity
        std::vector<Real> sdir(N_RAD_BAND, 1.0);
        std::vector<Real> sdif(N_RAD_BAND, 1.0);

        if (coszen > 0) {
          const Real tiny = 1.0e-6;
          const Real zen = std::acos(coszen);
          const Real z = std::max(zen, tiny);
          const Real val = std::min(1.0 / (hwr * std::tan(z)), 1.0);
          const Real theta0 = std::asin(val);
          const Real tanzen = std::tan(zen);
          const Real costheta0 = std::cos(theta0);
          const Real theta0OverPi = theta0 / rpi;

          for (int ib = 0; ib < N_RAD_BAND; ib++) {
            // director radiation
            s_shadewall(c, ib, 0) = 0.0; // eqn. 2.15
            s_road(c, ib, 0) = sdir[ib] * (2.0 * theta0OverPi -
                                           2.0 / rpi * hwr * tanzen *
                                               (1.0 - costheta0)); // eqn 2.17
            s_sunlitwall(c, ib, 0) =
                2.0 * sdir[ib] *
                ((1.0 / hwr) * (0.5 - theta0OverPi) +
                 (1.0 / rpi) * tanzen * (1.0 - costheta0)); // eqn. 2.16

            // diffuse radiation
            s_road(c, ib, 1) = sdif[ib] * vf_sr(c);       // eqn 2.30
            s_sunlitwall(c, ib, 1) = sdif[ib] * vf_sw(c); // eqn 2.32
            s_shadewall(c, ib, 1) = sdif[ib] * vf_sw(c);  // eqn 2.31
          }
        }
      });

  if (0) {
    printf("data_bundle.ShadeWall.DownwellingShortRad\n");
    print_view_3d(data_bundle.ShadedWall.DownwellingShortRad);
    printf("data_bundle.SunlitWall.DownwellingShortRad\n");
    print_view_3d(data_bundle.SunlitWall.DownwellingShortRad);
    printf("data_bundle.CombinedRoad.DownwellingShortRad\n");
    print_view_3d(data_bundle.CombinedRoad.DownwellingShortRad);
  }
}

void UrbanAlbedo::computeSnowAlbedo() {
  std::cout << "In compute_snow_albedo\n";

  const Real rpi = M_PI;
  const Real SNOW_ALBEDO_VIS = 0.66;
  const Real SNOW_ALBEDO_NIR = 0.56;

  int N_LUN = data_bundle.N_LUN;
  int N_RAD_BAND = data_bundle.N_RAD_BAND;
  Kokkos::parallel_for(
      "ComputeIncidentRadiation", N_LUN, KOKKOS_LAMBDA(const int c) {
        const Real coszen = data_bundle.input.Coszen(c);

        Array3DR8 albsn_roof = data_bundle.Roof.SnowAlbedo;
        // Array2DR8 albsni_roof = data_bundle.Roof.SnowAlbedo.dif;

        Array3DR8 albsn_improad = data_bundle.ImperviousRoad.SnowAlbedo;
        // Array2DR8 albsni_improad = data_bundle.ImperviousRoad.SnowAlbedo.dif;

        Array3DR8 albsn_perroad = data_bundle.PerviousRoad.SnowAlbedo;
        // Array2DR8 albsni_perroad = data_bundle.PerviousRoad.SnowAlbedo.dif;

        if (coszen > 0) {
          albsn_roof(0, c, 0) = SNOW_ALBEDO_VIS;
          albsn_roof(1, c, 0) = SNOW_ALBEDO_NIR;
          albsn_roof(0, c, 1) = SNOW_ALBEDO_VIS;
          albsn_roof(1, c, 1) = SNOW_ALBEDO_NIR;

          albsn_improad(0, c, 0) = SNOW_ALBEDO_VIS;
          albsn_improad(1, c, 0) = SNOW_ALBEDO_NIR;
          albsn_improad(0, c, 1) = SNOW_ALBEDO_VIS;
          albsn_improad(1, c, 1) = SNOW_ALBEDO_NIR;

          albsn_perroad(0, c, 0) = SNOW_ALBEDO_VIS;
          albsn_perroad(1, c, 0) = SNOW_ALBEDO_NIR;
          albsn_perroad(0, c, 1) = SNOW_ALBEDO_VIS;
          albsn_perroad(1, c, 1) = SNOW_ALBEDO_NIR;
        }
      });
}

void UrbanAlbedo::computeCombinedAlbedo() {
  std::cout << "In compute_combined_albedo\n";

  const Real rpi = M_PI;
  const Real SNOW_ALBEDO_VIS = 0.66;
  const Real SNOW_ALBEDO_NIR = 0.56;

  int N_LUN = data_bundle.N_LUN;
  int N_RAD_BAND = data_bundle.N_RAD_BAND;
  Kokkos::parallel_for(
      "ComputeIncidentRadiation", N_LUN, KOKKOS_LAMBDA(const int c) {
        const Real coszen = data_bundle.input.Coszen(c);

        Array3DR8 albsn_roof = data_bundle.Roof.SnowAlbedo;
        // Array2DR8 albsni_roof = data_bundle.Roof.SnowAlbedo.dif;
        Array3DR8 alb_roof = data_bundle.Roof.BaseAlbedo;
        // Array2DR8 alb_roof_dif = data_bundle.Roof.BaseAlbedo.dif;
        Array3DR8 alb_roof_s = data_bundle.Roof.AlbedoWithSnowEffects;
        // Array2DR8 alb_roof_dif_s =
        // data_bundle.Roof.AlbedoWithSnowEffects.dif;

        Array3DR8 albsn_improad = data_bundle.ImperviousRoad.SnowAlbedo;
        // Array2DR8 albsni_improad = data_bundle.ImperviousRoad.SnowAlbedo.dif;
        Array3DR8 alb_improad = data_bundle.ImperviousRoad.BaseAlbedo;
        // Array2DR8 alb_improad_dif =
        // data_bundle.ImperviousRoad.BaseAlbedo.dif;
        Array3DR8 alb_improad_s =
            data_bundle.ImperviousRoad.AlbedoWithSnowEffects;
        // Array2DR8 alb_improad_dif_s =
        //     data_bundle.ImperviousRoad.AlbedoWithSnowEffects.dif;

        Array3DR8 albsn_perroad = data_bundle.PerviousRoad.SnowAlbedo;
        // Array2DR8 albsni_perroad = data_bundle.PerviousRoad.SnowAlbedo.dif;
        Array3DR8 alb_perroad = data_bundle.PerviousRoad.BaseAlbedo;
        // Array2DR8 alb_perroad_dif = data_bundle.PerviousRoad.BaseAlbedo.dif;
        Array3DR8 alb_perroad_s =
            data_bundle.PerviousRoad.AlbedoWithSnowEffects;
        // Array2DR8 alb_perroad_dif_s =
        //     data_bundle.PerviousRoad.AlbedoWithSnowEffects.dif;

        const Real frac_sno = 0.0;

        for (int ib = 0; ib < N_RAD_BAND; ib++) {
          alb_roof_s(c, ib, 0) = alb_roof(c, ib, 0) * (1.0 - frac_sno) +
                                 albsn_roof(c, ib, 0) * frac_sno;
          alb_roof_s(c, ib, 1) = alb_roof(c, ib, 1) * (1.0 - frac_sno) +
                                 albsn_roof(c, ib, 0) * frac_sno;

          alb_improad_s(c, ib, 0) = alb_improad(c, ib, 0) * (1.0 - frac_sno) +
                                    albsn_improad(c, ib, 0) * frac_sno;
          alb_improad_s(c, ib, 1) = alb_improad(c, ib, 1) * (1.0 - frac_sno) +
                                    albsn_improad(c, ib, 1) * frac_sno;

          alb_perroad_s(c, ib, 0) = alb_perroad(c, ib, 0) * (1.0 - frac_sno) +
                                    albsn_perroad(c, ib, 1) * frac_sno;
          alb_perroad_s(c, ib, 1) = alb_perroad(c, ib, 0) * (1.0 - frac_sno) +
                                    albsn_perroad(c, ib, 1) * frac_sno;
        }
      });

  if (0) {
    print_view_3d(data_bundle.Roof.AlbedoWithSnowEffects);
    print_view_3d(data_bundle.ImperviousRoad.AlbedoWithSnowEffects);
    print_view_3d(data_bundle.PerviousRoad.AlbedoWithSnowEffects);
  }
}

void UrbanAlbedo::computeNetSolar() {
  std::cout << "In compute_net_solar\n";

  int N_LUN = data_bundle.N_LUN;
  int N_RAD_BAND = data_bundle.N_RAD_BAND;

  Kokkos::parallel_for(
      "ComputeNetSolar", N_LUN, KOKKOS_LAMBDA(const int c) {
        const Real coszen = data_bundle.input.Coszen(c);
        if (coszen > 0) {

          Array3DR8 SradRoad = data_bundle.CombinedRoad.DownwellingShortRad;
          Array3DR8 SradSunwall = data_bundle.SunlitWall.DownwellingShortRad;
          Array3DR8 SradShadewall = data_bundle.ShadedWall.DownwellingShortRad;

          Array3DR8 AbsImproad = data_bundle.ImperviousRoad.AbsorbedShortRad;
          Array3DR8 AbsPerroad = data_bundle.PerviousRoad.AbsorbedShortRad;
          Array3DR8 AbsSunwall = data_bundle.SunlitWall.AbsorbedShortRad;
          Array3DR8 AbsShadewall = data_bundle.ShadedWall.AbsorbedShortRad;

          Array3DR8 RefImproad = data_bundle.ImperviousRoad.ReflectedShortRad;
          Array3DR8 RefPerroad = data_bundle.PerviousRoad.ReflectedShortRad;
          Array3DR8 RefSunwall = data_bundle.SunlitWall.ReflectedShortRad;
          Array3DR8 RefShadewall = data_bundle.ShadedWall.ReflectedShortRad;

          Array1DR8 hwr = data_bundle.geometry.CanyonHwr;

          for (int rtype = 0; rtype < 2; rtype++) {

            for (int ib = 0; ib < N_RAD_BAND; ib++) {

              RadIndices idx{c, ib, rtype};

              // Impervious road
              RadOutput ImpRoadOut;
              RadRefComponents ImpRoadRef;

              Real StotForRoad = SradRoad(c, ib, rtype);

              ImperviousRoadNetSolar.ComputeAbsAndRefRad(idx, StotForRoad,
                                                         ImpRoadOut);

              ImperviousRoadNetSolar.ComputeRefRadByComponent(idx, StotForRoad,
                                                              ImpRoadRef);

              // Pervious road
              RadOutput PerRoadOut;
              RadRefComponents PerRoadRef;

              PerviousRoadNetSolar.ComputeAbsAndRefRad(idx, StotForRoad,
                                                       PerRoadOut);

              PerviousRoadNetSolar.ComputeRefRadByComponent(idx, StotForRoad,
                                                            PerRoadRef);

              // Combining data from impervious and pervious road
              Real RoadAbs, RoadRef;
              Real RoadRefToSky, RoadRefToSunwall, RoadRefToShadewall;

              RoadAbs = ImpRoadOut.Abs + PerRoadOut.AbsByWt;
              RoadRef = ImpRoadOut.RefByWt + PerRoadOut.RefByWt;
              RoadRefToSky = ImpRoadRef.ToSkyByWt + PerRoadRef.ToSkyByWt;
              RoadRefToSunwall =
                  ImpRoadRef.ToSunwallByWt + PerRoadRef.ToSunwallByWt;
              RoadRefToShadewall =
                  ImpRoadRef.ToShadewallByWt + PerRoadRef.ToShadewallByWt;

              // Sunlit wall
              Real StotForSunwall = SradSunwall(c, ib, rtype);

              RadOutput SunwallOut;
              RadRefComponents SunwallRef;

              SunlitWallNetSolar.ComputeAbsAndRefRad(idx, StotForSunwall,
                                                     SunwallOut);
              SunlitWallNetSolar.ComputeRefRadByComponent(
                  idx, SradSunwall(c, ib, rtype), SunwallRef);

              // Shadedwall
              Real shadewall_a, shadewall_r;
              Real shadewall_r_sky, shadewall_r_road, shadewall_r_sunwall;
              Real StotForShadewall = SradShadewall(c, ib, rtype);

              RadOutput ShadewallOut;
              RadRefComponents ShadewallRef;

              ShadeWallNetSolar.ComputeAbsAndRefRad(idx, StotForSunwall,
                                                    ShadewallOut);
              ShadeWallNetSolar.ComputeRefRadByComponent(idx, StotForShadewall,
                                                         ShadewallRef);

              // Cummulative absorbed and reflected radiation for the four
              // surfaces

              AbsImproad(c, ib, rtype) = ImpRoadOut.Abs;
              AbsPerroad(c, ib, rtype) = PerRoadOut.Abs;
              AbsSunwall(c, ib, rtype) = SunwallOut.Abs;
              AbsShadewall(c, ib, rtype) = ShadewallOut.Abs;

              RefImproad(c, ib, rtype) = ImpRoadRef.ToSky;
              RefPerroad(c, ib, rtype) = PerRoadRef.ToSky;
              RefSunwall(c, ib, rtype) = SunwallRef.ToSky;
              RefShadewall(c, ib, rtype) = ShadewallRef.ToSky;

              const int max_iter = 50;
              for (int iter = 0; iter < max_iter; iter++) {
                // step(1)
                StotForRoad =
                    (SunwallRef.ToRoad + ShadewallRef.ToRoad) * hwr(c);

                ImperviousRoadNetSolar.ComputeAbsAndRefRad(idx, StotForRoad,
                                                           ImpRoadOut);

                PerviousRoadNetSolar.ComputeAbsAndRefRad(idx, StotForRoad,
                                                         PerRoadOut);

                RoadAbs = ImpRoadOut.AbsByWt + PerRoadOut.AbsByWt;
                RoadRef = ImpRoadOut.RefByWt + PerRoadOut.RefByWt;

                StotForSunwall =
                    RoadRefToSunwall / hwr(c) + ShadewallRef.ToOtherwall;
                SunlitWallNetSolar.ComputeAbsAndRefRad(idx, StotForSunwall,
                                                       SunwallOut);

                StotForShadewall =
                    RoadRefToShadewall / hwr(c) + SunwallRef.ToOtherwall;
                ShadeWallNetSolar.ComputeAbsAndRefRad(idx, StotForShadewall,
                                                      ShadewallOut);

                // step (2)
                AbsImproad(c, ib, rtype) += ImpRoadOut.Abs;
                AbsPerroad(c, ib, rtype) += PerRoadOut.Abs;
                AbsSunwall(c, ib, rtype) += SunwallOut.Abs;
                AbsShadewall(c, ib, rtype) += ShadewallOut.Abs;

                // step (3)
                ImperviousRoadNetSolar.ComputeRefRadByComponent(
                    idx, StotForRoad, ImpRoadRef);

                PerviousRoadNetSolar.ComputeRefRadByComponent(idx, StotForRoad,
                                                              PerRoadRef);

                RoadRefToSky = ImpRoadRef.ToSkyByWt + PerRoadRef.ToSkyByWt;
                RoadRefToSunwall =
                    ImpRoadRef.ToSunwallByWt + PerRoadRef.ToSunwallByWt;
                RoadRefToShadewall =
                    ImpRoadRef.ToShadewallByWt + PerRoadRef.ToShadewallByWt;

                SunlitWallNetSolar.ComputeRefRadByComponent(idx, StotForSunwall,
                                                            SunwallRef);

                ShadeWallNetSolar.ComputeRefRadByComponent(
                    idx, StotForShadewall, ShadewallRef);

                // step (4)
                RefImproad(c, ib, rtype) += ImpRoadRef.ToSky;
                RefPerroad(c, ib, rtype) += PerRoadRef.ToSky;
                RefSunwall(c, ib, rtype) += SunwallRef.ToSky;
                RefShadewall(c, ib, rtype) += ShadewallRef.ToSky;

                Real crit =
                    std::max({RoadAbs, SunwallOut.Abs, ShadewallOut.Abs});
                if (c == 0)
                  printf("(%d, %d, %d): %d %e %e %e\n", c, ib, rtype, iter,
                         RoadAbs, SunwallOut.Abs, ShadewallOut.Abs);
                const Real errcrit = 0.00001;
                if (crit < errcrit)
                  break;
              }
            }
          }

        } else {
        }
      });

  if (0) {
    // print_view_2d(data_bundle.Roof.AlbedoWithSnowEffects.dir);
  }
}

} // namespace URBANXX