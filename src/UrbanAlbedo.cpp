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
        Array2DR8 sdir_shadewall =
            data_bundle.ShadedWall.DownwellingShortRad.dir;
        Array2DR8 sdir_road = data_bundle.CombinedRoad.DownwellingShortRad.dir;

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
            sdir_shadewall(ib, c) = 0.0; // eqn. 2.15
            sdir_road(ib, c) = sdir[ib] * (2.0 * theta0OverPi -
                                           2.0 / rpi * hwr * tanzen *
                                               (1.0 - costheta0)); // eqn 2.17
            sdir_sunlitwall(ib, c) =
                2.0 * sdir[ib] *
                ((1.0 / hwr) * (0.5 - theta0OverPi) +
                 (1.0 / rpi) * tanzen * (1.0 - costheta0)); // eqn. 2.16
          }
        }
      });

  if (0) {
    print_view_2d(data_bundle.ShadedWall.DownwellingShortRad.dir);
    print_view_2d(data_bundle.SunlitWall.DownwellingShortRad.dir);
    print_view_2d(data_bundle.CombinedRoad.DownwellingShortRad.dir);
  }
}

} // namespace URBANXX