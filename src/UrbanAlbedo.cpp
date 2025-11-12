#include <Kokkos_Core.hpp>
#include <UrbanAlbedo.hpp>
#include <UrbanData.hpp>
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
                                       RadOutput &out,
                                       bool scale_by_weight) const {

  Real albedo = RoadData.BaseAlbedo(idx.c, idx.ib, idx.rtype);

  out.Abs = (1.0 - albedo) * InRad;
  out.Ref = albedo * InRad;
  if (scale_by_weight) {
    out.Abs *= Weight;
    out.Ref *= Weight;
  }
}

KOKKOS_FUNCTION
void NetSolarRoad::ComputeRefRadByComponent(RadIndices idx, Real InRad,
                                            RadRefComponents &ref,
                                            bool scale_by_weight) const {

  Real albedo = RoadData.BaseAlbedo(idx.c, idx.ib, idx.rtype);
  Real vf_sr = ViewFactorSkyFromRoad(idx.c);
  Real vf_wr = ViewFactorWallFromRoad(idx.c);

  Real RefRad = albedo * InRad;

  ref.ToSky = RefRad * vf_sr;
  ref.ToSunwall = RefRad * vf_wr;
  ref.ToShadewall = RefRad * vf_wr;
  ref.ToOtherwall = 0.0;
  ref.ToRoad = 0.0;

  if (scale_by_weight) {
    ref.ToSky *= Weight;
    ref.ToSunwall *= Weight;
    ref.ToShadewall *= Weight;
  }
}

KOKKOS_FUNCTION
void NetSolarWall::ComputeAbsAndRefRad(RadIndices idx, Real InRad,
                                       RadOutput &out) const {

  Real albedo = WallData.BaseAlbedo(idx.c, idx.ib, idx.rtype);

  out.Abs = (1.0 - albedo) * InRad;
  out.Ref = albedo * InRad;
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

void UrbanAlbedo::set_solar_inputs() const {

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

void UrbanAlbedo::compute_incident_radiation() const {
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

void UrbanAlbedo::compute_snow_albedo() const {
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

void UrbanAlbedo::compute_combined_albedo() const {
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

void UrbanAlbedo::compute_net_solar() const {
  std::cout << "In compute_net_solar\n";

  int N_LUN = data_bundle.N_LUN;
  int N_RAD_BAND = data_bundle.N_RAD_BAND;

  Kokkos::parallel_for(
      "ComputeNetSolar", N_LUN, KOKKOS_LAMBDA(const int c) {
        const Real coszen = data_bundle.input.Coszen(c);
        if (coszen > 0) {

          Array3DR8 srad_road = data_bundle.CombinedRoad.DownwellingShortRad;
          Array3DR8 srad_sunlitwall =
              data_bundle.SunlitWall.DownwellingShortRad;
          Array3DR8 srad_shadewall = data_bundle.ShadedWall.DownwellingShortRad;

          Array3DR8 sabs_improad = data_bundle.ImperviousRoad.AbsorbedShortRad;
          Array3DR8 sabs_perroad = data_bundle.PerviousRoad.AbsorbedShortRad;
          Array3DR8 sabs_sunwall = data_bundle.SunlitWall.AbsorbedShortRad;
          Array3DR8 sabs_shadewall = data_bundle.ShadedWall.AbsorbedShortRad;

          Array3DR8 sref_improad = data_bundle.ImperviousRoad.ReflectedShortRad;
          Array3DR8 sref_perroad = data_bundle.PerviousRoad.ReflectedShortRad;
          Array3DR8 sref_sunwall = data_bundle.SunlitWall.ReflectedShortRad;
          Array3DR8 sref_shadewall = data_bundle.ShadedWall.ReflectedShortRad;

          Array1DR8 hwr = data_bundle.geometry.CanyonHwr;

          bool scale_by_weight = true;

          for (int rtype = 0; rtype < 2; rtype++) {

            for (int ib = 0; ib < N_RAD_BAND; ib++) {

              RadIndices idx{c, ib, rtype};

              // Impervious road
              RadOutput improad_out, improad_out_by_wt;
              RadRefComponents improad_ref, improad_ref_by_wt;

              Real stot_for_road = srad_road(c, ib, rtype);

              ImperviousRoadNetSolar.ComputeAbsAndRefRad(
                  idx, stot_for_road, improad_out, !scale_by_weight);

              ImperviousRoadNetSolar.ComputeAbsAndRefRad(
                  idx, stot_for_road, improad_out_by_wt, scale_by_weight);

              ImperviousRoadNetSolar.ComputeRefRadByComponent(
                  idx, stot_for_road, improad_ref, !scale_by_weight);

              ImperviousRoadNetSolar.ComputeRefRadByComponent(
                  idx, stot_for_road, improad_ref_by_wt, scale_by_weight);

              // Pervious road
              RadOutput perroad_out, perroad_out_by_wt;
              RadRefComponents perroad_ref, perroad_ref_by_wt;

              PerviousRoadNetSolar.ComputeAbsAndRefRad(
                  idx, stot_for_road, perroad_out, !scale_by_weight);

              PerviousRoadNetSolar.ComputeAbsAndRefRad(
                  idx, stot_for_road, perroad_out_by_wt, scale_by_weight);

              PerviousRoadNetSolar.ComputeRefRadByComponent(
                  idx, stot_for_road, perroad_ref, !scale_by_weight);

              PerviousRoadNetSolar.ComputeRefRadByComponent(
                  idx, stot_for_road, perroad_ref_by_wt, scale_by_weight);

              // Combining data from impervious and pervious road
              Real road_a, road_r;
              Real road_r_sky, road_r_sunwall, road_r_shadewall;

              road_a = improad_out_by_wt.Abs + perroad_out_by_wt.Abs;
              road_r = improad_out_by_wt.Ref + perroad_out_by_wt.Ref;
              road_r_sky = improad_ref_by_wt.ToSky + perroad_ref_by_wt.ToSky;
              road_r_sunwall =
                  improad_ref_by_wt.ToSunwall + perroad_ref_by_wt.ToSunwall;
              road_r_shadewall =
                  improad_ref_by_wt.ToShadewall + perroad_ref_by_wt.ToShadewall;

              // Sunlit wall
              Real stot_for_sunwall = srad_sunlitwall(c, ib, rtype);

              RadOutput sunwall_out;
              RadRefComponents sunwall_ref;

              SunlitWallNetSolar.ComputeAbsAndRefRad(idx, stot_for_sunwall,
                                                     sunwall_out);
              SunlitWallNetSolar.ComputeRefRadByComponent(
                  idx, srad_sunlitwall(c, ib, rtype), sunwall_ref);

              // Shadedwall
              Real shadewall_a, shadewall_r;
              Real shadewall_r_sky, shadewall_r_road, shadewall_r_sunwall;
              Real stot_for_shadewall = srad_shadewall(c, ib, rtype);

              RadOutput shadewall_out;
              RadRefComponents shadewall_ref;

              ShadeWallNetSolar.ComputeAbsAndRefRad(idx, stot_for_sunwall,
                                                    shadewall_out);
              ShadeWallNetSolar.ComputeRefRadByComponent(
                  idx, stot_for_shadewall, shadewall_ref);

              // Cummulative absorbed and reflected radiation for the four
              // surfaces

              sabs_improad(c, ib, rtype) = improad_out.Abs;
              sabs_perroad(c, ib, rtype) = perroad_out.Abs;
              sabs_sunwall(c, ib, rtype) = sunwall_out.Abs;
              sabs_shadewall(c, ib, rtype) = shadewall_out.Abs;

              sref_improad(c, ib, rtype) = improad_out.Ref;
              sref_perroad(c, ib, rtype) = perroad_out.Ref;
              sref_sunwall(c, ib, rtype) = sunwall_out.Ref;
              sref_shadewall(c, ib, rtype) = shadewall_out.Ref;

              const int max_iter = 50;
              for (int iter = 0; iter < max_iter; iter++) {
                // step(1)
                stot_for_road =
                    (sunwall_ref.ToRoad + shadewall_ref.ToRoad) * hwr(c);

                ImperviousRoadNetSolar.ComputeAbsAndRefRad(
                    idx, stot_for_road, improad_out, !scale_by_weight);

                ImperviousRoadNetSolar.ComputeAbsAndRefRad(
                    idx, stot_for_road, improad_out_by_wt, scale_by_weight);

                PerviousRoadNetSolar.ComputeAbsAndRefRad(
                    idx, stot_for_road, perroad_out, !scale_by_weight);

                PerviousRoadNetSolar.ComputeAbsAndRefRad(
                    idx, stot_for_road, perroad_out_by_wt, scale_by_weight);

                road_a = improad_out_by_wt.Abs + perroad_out_by_wt.Abs;
                road_r = improad_out_by_wt.Ref + perroad_out_by_wt.Ref;

                stot_for_sunwall =
                    road_r_sunwall / hwr(c) + shadewall_ref.ToSunwall;
                SunlitWallNetSolar.ComputeAbsAndRefRad(idx, stot_for_sunwall,
                                                       sunwall_out);

                stot_for_shadewall =
                    road_r_shadewall / hwr(c) + sunwall_ref.ToShadewall;
                ShadeWallNetSolar.ComputeAbsAndRefRad(idx, stot_for_shadewall,
                                                      shadewall_out);

                // step (2)
                sabs_improad(c, ib, rtype) += improad_out.Abs;
                sabs_perroad(c, ib, rtype) += perroad_out.Abs;
                sabs_sunwall(c, ib, rtype) += sunwall_out.Abs;
                sabs_shadewall(c, ib, rtype) += shadewall_out.Abs;

                // step (3)
                ImperviousRoadNetSolar.ComputeRefRadByComponent(
                    idx, stot_for_road, improad_ref_by_wt, scale_by_weight);

                PerviousRoadNetSolar.ComputeRefRadByComponent(
                    idx, stot_for_road, perroad_ref_by_wt, scale_by_weight);

                road_r_sky = improad_ref_by_wt.ToSky + perroad_ref_by_wt.ToSky;
                road_r_sunwall =
                    improad_ref_by_wt.ToSunwall + perroad_ref_by_wt.ToSunwall;
                road_r_shadewall = improad_ref_by_wt.ToShadewall +
                                   perroad_ref_by_wt.ToShadewall;

                SunlitWallNetSolar.ComputeRefRadByComponent(
                    idx, stot_for_sunwall, sunwall_ref);

                ShadeWallNetSolar.ComputeRefRadByComponent(
                    idx, stot_for_shadewall, shadewall_ref);

                // step (4)
                sref_improad(c, ib, rtype) += improad_out.Ref;
                sref_perroad(c, ib, rtype) += perroad_out.Ref;
                sref_sunwall(c, ib, rtype) += sunwall_out.Ref;
                sref_shadewall(c, ib, rtype) += shadewall_out.Ref;

                Real crit =
                    std::max({road_a, sunwall_out.Abs, shadewall_out.Abs});
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