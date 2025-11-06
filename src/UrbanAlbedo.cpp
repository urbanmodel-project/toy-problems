#include <Kokkos_Core.hpp>
#include <UrbanAlbedo.hpp>
#include <UrbanData.hpp>
#include <cmath>
#include <iomanip>
#include <iostream>

namespace URBANXX {

UrbanAlbedo::UrbanAlbedo(UrbanSharedDataBundle &bundle) : data_bundle(bundle) {}

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
  int N_RAD = data_bundle.N_RAD;

  Kokkos::parallel_for(
      "ComputeIncidentRadiation", N_LUN, KOKKOS_LAMBDA(const int c) {
        const Real coszen = data_bundle.input.Coszen(c);
        const Real hwr = data_bundle.geometry.CanyonHwr(c);

        Array2DR8 sdir_sunlitwall =
            data_bundle.SunlitWall.DownwellingShortRad.dir;
        Array2DR8 sdif_sunlitwall =
            data_bundle.SunlitWall.DownwellingShortRad.dif;

        Array2DR8 sdir_shadewall =
            data_bundle.ShadedWall.DownwellingShortRad.dir;
        Array2DR8 sdif_shadewall =
            data_bundle.ShadedWall.DownwellingShortRad.dif;

        Array2DR8 sdir_road = data_bundle.CombinedRoad.DownwellingShortRad.dir;
        Array2DR8 sdif_road = data_bundle.CombinedRoad.DownwellingShortRad.dif;

        Array1DR8 vf_sr = data_bundle.geometry.ViewFactorSkyForRoad;
        Array1DR8 vf_sw = data_bundle.geometry.ViewFactorSkyForWall;

        // the incident direct and diffuse radiation for VIS and NIR bands is
        // assumed to be unity
        std::vector<Real> sdir(N_RAD, 1.0);
        std::vector<Real> sdif(N_RAD, 1.0);

        if (coszen > 0) {
          const Real tiny = 1.0e-6;
          const Real zen = std::acos(coszen);
          const Real z = std::max(zen, tiny);
          const Real val = std::min(1.0 / (hwr * std::tan(z)), 1.0);
          const Real theta0 = std::asin(val);
          const Real tanzen = std::tan(zen);
          const Real costheta0 = std::cos(theta0);
          const Real theta0OverPi = theta0 / rpi;

          for (int ib = 0; ib < N_RAD; ib++) {
            // director radiation
            sdir_shadewall(ib, c) = 0.0; // eqn. 2.15
            sdir_road(ib, c) = sdir[ib] * (2.0 * theta0OverPi -
                                           2.0 / rpi * hwr * tanzen *
                                               (1.0 - costheta0)); // eqn 2.17
            sdir_sunlitwall(ib, c) =
                2.0 * sdir[ib] *
                ((1.0 / hwr) * (0.5 - theta0OverPi) +
                 (1.0 / rpi) * tanzen * (1.0 - costheta0)); // eqn. 2.16

            // diffuse radiation
            sdif_road(ib, c) = sdif[ib] * vf_sr(c);       // eqn 2.30
            sdif_sunlitwall(ib, c) = sdif[ib] * vf_sw(c); // eqn 2.32
            sdif_shadewall(ib, c) = sdif[ib] * vf_sw(c);  // eqn 2.31
          }
        }
      });

  if (0) {
    print_view_2d(data_bundle.ShadedWall.DownwellingShortRad.dir);
    print_view_2d(data_bundle.SunlitWall.DownwellingShortRad.dir);
    print_view_2d(data_bundle.CombinedRoad.DownwellingShortRad.dir);

    print_view_2d(data_bundle.ShadedWall.DownwellingShortRad.dif);
    print_view_2d(data_bundle.SunlitWall.DownwellingShortRad.dif);
    print_view_2d(data_bundle.CombinedRoad.DownwellingShortRad.dif);
  }
}

void UrbanAlbedo::compute_snow_albedo() const {
  std::cout << "In compute_snow_albedo\n";

  const Real rpi = M_PI;
  const Real SNOW_ALBEDO_VIS = 0.66;
  const Real SNOW_ALBEDO_NIR = 0.56;

  int N_LUN = data_bundle.N_LUN;
  int N_RAD = data_bundle.N_RAD;
  Kokkos::parallel_for(
      "ComputeIncidentRadiation", N_LUN, KOKKOS_LAMBDA(const int c) {
        const Real coszen = data_bundle.input.Coszen(c);

        Array2DR8 albsnd_roof = data_bundle.Roof.SnowAlbedo.dir;
        Array2DR8 albsni_roof = data_bundle.Roof.SnowAlbedo.dif;

        Array2DR8 albsnd_improad = data_bundle.ImperviousRoad.SnowAlbedo.dir;
        Array2DR8 albsni_improad = data_bundle.ImperviousRoad.SnowAlbedo.dif;

        Array2DR8 albsnd_perroad = data_bundle.PerviousRoad.SnowAlbedo.dir;
        Array2DR8 albsni_perroad = data_bundle.PerviousRoad.SnowAlbedo.dif;

        if (coszen > 0) {
          albsnd_roof(0, c) = SNOW_ALBEDO_VIS;
          albsnd_roof(1, c) = SNOW_ALBEDO_NIR;
          albsni_roof(0, c) = SNOW_ALBEDO_VIS;
          albsni_roof(1, c) = SNOW_ALBEDO_NIR;

          albsnd_improad(0, c) = SNOW_ALBEDO_VIS;
          albsnd_improad(1, c) = SNOW_ALBEDO_NIR;
          albsni_improad(0, c) = SNOW_ALBEDO_VIS;
          albsni_improad(1, c) = SNOW_ALBEDO_NIR;

          albsnd_perroad(0, c) = SNOW_ALBEDO_VIS;
          albsnd_perroad(1, c) = SNOW_ALBEDO_NIR;
          albsni_perroad(0, c) = SNOW_ALBEDO_VIS;
          albsni_perroad(1, c) = SNOW_ALBEDO_NIR;
        }
      });
}

void UrbanAlbedo::compute_combined_albedo() const {
  std::cout << "In compute_combined_albedo\n";

  const Real rpi = M_PI;
  const Real SNOW_ALBEDO_VIS = 0.66;
  const Real SNOW_ALBEDO_NIR = 0.56;

  int N_LUN = data_bundle.N_LUN;
  int N_RAD = data_bundle.N_RAD;
  Kokkos::parallel_for(
      "ComputeIncidentRadiation", N_LUN, KOKKOS_LAMBDA(const int c) {
        const Real coszen = data_bundle.input.Coszen(c);

        Array2DR8 albsnd_roof = data_bundle.Roof.SnowAlbedo.dir;
        Array2DR8 albsni_roof = data_bundle.Roof.SnowAlbedo.dif;
        Array2DR8 alb_roof_dir = data_bundle.Roof.BaseAlbedo.dir;
        Array2DR8 alb_roof_dif = data_bundle.Roof.BaseAlbedo.dif;
        Array2DR8 alb_roof_dir_s = data_bundle.Roof.AlbedoWithSnowEffects.dir;
        Array2DR8 alb_roof_dif_s = data_bundle.Roof.AlbedoWithSnowEffects.dif;

        Array2DR8 albsnd_improad = data_bundle.ImperviousRoad.SnowAlbedo.dir;
        Array2DR8 albsni_improad = data_bundle.ImperviousRoad.SnowAlbedo.dif;
        Array2DR8 alb_improad_dir = data_bundle.ImperviousRoad.BaseAlbedo.dir;
        Array2DR8 alb_improad_dif = data_bundle.ImperviousRoad.BaseAlbedo.dif;
        Array2DR8 alb_improad_dir_s =
            data_bundle.ImperviousRoad.AlbedoWithSnowEffects.dir;
        Array2DR8 alb_improad_dif_s =
            data_bundle.ImperviousRoad.AlbedoWithSnowEffects.dif;

        Array2DR8 albsnd_perroad = data_bundle.PerviousRoad.SnowAlbedo.dir;
        Array2DR8 albsni_perroad = data_bundle.PerviousRoad.SnowAlbedo.dif;
        Array2DR8 alb_perroad_dir = data_bundle.PerviousRoad.BaseAlbedo.dir;
        Array2DR8 alb_perroad_dif = data_bundle.PerviousRoad.BaseAlbedo.dif;
        Array2DR8 alb_perroad_dir_s =
            data_bundle.PerviousRoad.AlbedoWithSnowEffects.dir;
        Array2DR8 alb_perroad_dif_s =
            data_bundle.PerviousRoad.AlbedoWithSnowEffects.dif;

        const Real frac_sno = 0.0;

        for (int ib = 0; ib < N_RAD; ib++) {
          alb_roof_dir_s(ib, c) = alb_roof_dir(ib, c) * (1.0 - frac_sno) +
                                  albsnd_roof(ib, c) * frac_sno;
          alb_roof_dif_s(ib, c) = alb_roof_dif(ib, c) * (1.0 - frac_sno) +
                                  albsni_roof(ib, c) * frac_sno;

          alb_improad_dir_s(ib, c) = alb_improad_dir(ib, c) * (1.0 - frac_sno) +
                                     albsnd_improad(ib, c) * frac_sno;
          alb_improad_dif_s(ib, c) = alb_improad_dir(ib, c) * (1.0 - frac_sno) +
                                     albsni_improad(ib, c) * frac_sno;

          alb_perroad_dir_s(ib, c) = alb_perroad_dir(ib, c) * (1.0 - frac_sno) +
                                     albsnd_perroad(ib, c) * frac_sno;
          alb_perroad_dif_s(ib, c) = alb_perroad_dir(ib, c) * (1.0 - frac_sno) +
                                     albsni_perroad(ib, c) * frac_sno;
        }
      });

  if (0) {
    print_view_2d(data_bundle.Roof.AlbedoWithSnowEffects.dir);
    print_view_2d(data_bundle.ImperviousRoad.AlbedoWithSnowEffects.dir);
    print_view_2d(data_bundle.PerviousRoad.AlbedoWithSnowEffects.dir);

    print_view_2d(data_bundle.Roof.AlbedoWithSnowEffects.dif);
    print_view_2d(data_bundle.ImperviousRoad.AlbedoWithSnowEffects.dif);
    print_view_2d(data_bundle.PerviousRoad.AlbedoWithSnowEffects.dif);
  }
}

} // namespace URBANXX