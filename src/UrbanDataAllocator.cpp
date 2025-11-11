#include <Kokkos_Core.hpp>
#include <UrbanData.hpp>
#include <UrbanDataAllocator.hpp>
#include <iomanip>
#include <iostream>

namespace URBANXX {
// Constructor Definition
UrbanDataAllocator::UrbanDataAllocator(UrbanSharedDataBundle &bundle)
    : data_bundle(bundle) {
  // Implementation is simple: just initializes the reference member.
}

#define ALLOCATE_VIEW(viewname, type, ...)                                     \
  viewname = type(#viewname, __VA_ARGS__);

#define RADIATION_TYPE_ALLOCATE_VIEW(base, type, ...)                          \
  ALLOCATE_VIEW(base.dir, type, __VA_ARGS__);                                  \
  ALLOCATE_VIEW(base.dif, type, __VA_ARGS__)

void CanyonGeometryAllocateViews(int N_LUN, CanyonGeometryData &geometry) {

  ALLOCATE_VIEW(geometry.CanyonHwrH, HostArray1DR8, N_LUN);
  ALLOCATE_VIEW(geometry.CanyonHwr, Array1DR8, N_LUN);

  ALLOCATE_VIEW(geometry.ViewFactorSkyFromRoad, Array1DR8, N_LUN);
  ALLOCATE_VIEW(geometry.ViewFactorSkyFromWall, Array1DR8, N_LUN);
  ALLOCATE_VIEW(geometry.ViewFactorRoadFromWall, Array1DR8, N_LUN);
  ALLOCATE_VIEW(geometry.ViewFactorOtherWallFromWall, Array1DR8, N_LUN);
  ALLOCATE_VIEW(geometry.ViewFactorWallFromRoad, Array1DR8, N_LUN);

  printf("All CanopyGeomtry Views successfully allocated on device.\n");
}

void SolarInputDataAllocateViews(int N_LUN, int N_RAD, SolarInputData &solar) {

  ALLOCATE_VIEW(solar.CoszenH, HostArray1DR8, N_LUN);
  ALLOCATE_VIEW(solar.Coszen, Array1DR8, N_LUN);

  ALLOCATE_VIEW(solar.SdirHorizH, HostArray2DR8, N_RAD, N_LUN);
  ALLOCATE_VIEW(solar.SdirHoriz, Array2DR8, N_RAD, N_LUN);

  ALLOCATE_VIEW(solar.SdifHorizH, HostArray2DR8, N_RAD, N_LUN);
  ALLOCATE_VIEW(solar.SdifHoriz, Array2DR8, N_RAD, N_LUN);

  ALLOCATE_VIEW(solar.FracSnowH, HostArray1DR8, N_LUN);
  ALLOCATE_VIEW(solar.FracSnow, Array1DR8, N_LUN);

  printf("All SolarInputData Views successfully allocated on device.\n");
}

void CombinedRoadDataTypeAllocateViews(int N_LUN, int N_RAD,
                                       CombinedRoadDataType &combineRoad) {
  RADIATION_TYPE_ALLOCATE_VIEW(combineRoad.DownwellingShortRad, Array2DR8,
                               N_RAD, N_LUN);
  RADIATION_TYPE_ALLOCATE_VIEW(combineRoad.SnowAlbedo, Array2DR8, N_RAD, N_LUN);
  RADIATION_TYPE_ALLOCATE_VIEW(combineRoad.AlbedoWithSnowEffects, Array2DR8,
                               N_RAD, N_LUN);
}

void RoadDataTypeAllocateViews(int N_LUN, int N_RAD, RoadDataType &road) {
  RADIATION_TYPE_ALLOCATE_VIEW(road.SnowAlbedo, Array2DR8, N_RAD, N_LUN);
  RADIATION_TYPE_ALLOCATE_VIEW(road.AlbedoWithSnowEffects, Array2DR8, N_RAD,
                               N_LUN);
  RADIATION_TYPE_ALLOCATE_VIEW(road.BaseAlbedo, Array2DR8, N_RAD, N_LUN);
  RADIATION_TYPE_ALLOCATE_VIEW(road.ReflectedShortRad, Array2DR8, N_RAD, N_LUN);
  RADIATION_TYPE_ALLOCATE_VIEW(road.AbsorbedShortRad, Array2DR8, N_RAD, N_LUN);
}

void WallDataTypeAllocateViews(int N_LUN, int N_RAD, WallDataType &wall) {
  RADIATION_TYPE_ALLOCATE_VIEW(wall.DownwellingShortRad, Array2DR8, N_RAD,
                               N_LUN);
  RADIATION_TYPE_ALLOCATE_VIEW(wall.BaseAlbedo, Array2DR8, N_RAD, N_LUN);
  RADIATION_TYPE_ALLOCATE_VIEW(wall.ReflectedShortRad, Array2DR8, N_RAD, N_LUN);
  RADIATION_TYPE_ALLOCATE_VIEW(wall.AbsorbedShortRad, Array2DR8, N_RAD, N_LUN);
}

void RoofDataTypeAllocateViews(int N_LUN, int N_RAD, RoofDataType &roof) {
  RADIATION_TYPE_ALLOCATE_VIEW(roof.SnowAlbedo, Array2DR8, N_RAD, N_LUN);
  RADIATION_TYPE_ALLOCATE_VIEW(roof.AlbedoWithSnowEffects, Array2DR8, N_RAD,
                               N_LUN);
  RADIATION_TYPE_ALLOCATE_VIEW(roof.BaseAlbedo, Array2DR8, N_RAD, N_LUN);
  RADIATION_TYPE_ALLOCATE_VIEW(roof.ReflectedShortRad, Array2DR8, N_RAD, N_LUN);
  RADIATION_TYPE_ALLOCATE_VIEW(roof.AbsorbedShortRad, Array2DR8, N_RAD, N_LUN);
}

// Method Definition: Contains the heavy lifting of allocation
void UrbanDataAllocator::allocate_all_views() const {
  int N_LUN = data_bundle.N_LUN;
  int N_RAD = data_bundle.N_RAD;

  printf("Starting device memory allocation for %d land units and %d bands.\n",
         N_LUN, N_RAD);

  CanyonGeometryAllocateViews(N_LUN, data_bundle.geometry);
  SolarInputDataAllocateViews(N_LUN, N_RAD, data_bundle.input);
  CombinedRoadDataTypeAllocateViews(N_LUN, N_RAD, data_bundle.CombinedRoad);
  WallDataTypeAllocateViews(N_LUN, N_RAD, data_bundle.SunlitWall);
  WallDataTypeAllocateViews(N_LUN, N_RAD, data_bundle.ShadedWall);
  CombinedRoadDataTypeAllocateViews(N_LUN, N_RAD, data_bundle.CombinedRoad);
  RoadDataTypeAllocateViews(N_LUN, N_RAD, data_bundle.ImperviousRoad);
  RoadDataTypeAllocateViews(N_LUN, N_RAD, data_bundle.PerviousRoad);
  RoofDataTypeAllocateViews(N_LUN, N_RAD, data_bundle.Roof);

  printf("All primary Views successfully allocated on device.\n");
}

template <typename ViewType>
void print_view_1d(const ViewType &view, const std::string &name = "") {
  auto h_view = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), view);
  std::cout << name << ":\n";
  for (std::size_t i = 0; i < h_view.extent(0); ++i) {
    std::cout << i << " " << std::setprecision(15) << h_view(i) << "\n";
  }
  std::cout << "\n";
}

void UrbanDataAllocator::initialize_canyon_geometry() const {

  int N_LUN = data_bundle.N_LUN;
  int N_RAD = data_bundle.N_RAD;

  const Real HWR_TBD = 4.80000019073486;
  // const Real HWR_HD  = 1.60000002384186;
  // const Real HWR_MD  = 0.600000023841858;

  HostArray1DR8 CanyonHwrH = data_bundle.geometry.CanyonHwrH;
  Array1DR8 CanyonHwr = data_bundle.geometry.CanyonHwr;
  for (int i = 0; i < N_LUN; i++) {
    CanyonHwrH(i) = HWR_TBD;
  }

  Kokkos::deep_copy(CanyonHwr, CanyonHwrH);

  if (0)
    print_view_1d(CanyonHwr, "CanyonHwr");

  Kokkos::parallel_for(
      "ComputingViewFactors", N_LUN, KOKKOS_LAMBDA(const int c) {
        Array1DR8 sr = data_bundle.geometry.ViewFactorSkyFromRoad;
        Array1DR8 sw = data_bundle.geometry.ViewFactorSkyFromWall;
        Array1DR8 rw = data_bundle.geometry.ViewFactorRoadFromWall;
        Array1DR8 ww = data_bundle.geometry.ViewFactorOtherWallFromWall;
        Array1DR8 wr = data_bundle.geometry.ViewFactorWallFromRoad;

        const Real hwr = data_bundle.geometry.CanyonHwr(c);

        const Real sqrt_term = std::sqrt(hwr * hwr + 1.0);

        sr(c) = sqrt_term - hwr;                     // eqn 2.25
        wr(c) = 0.5 * (1.0 - sr(c));                 // eqn 2.27
        sw(c) = 0.5 * (hwr + 1.0 - sqrt_term) / hwr; // eqn 2.24
        rw(c) = sw(c);                               // eqn 2.27
        ww(c) = 1.0 - sw(c) - rw(c);                 // eqn 2.28
      });
}

void UrbanDataAllocator::initialize_properties() const {

  int N_LUN = data_bundle.N_LUN;
  int N_RAD = data_bundle.N_RAD;

  const Real ALB_IMPROAD_DIR = 0.230000004172325;
  const Real ALB_IMPROAD_DIF = ALB_IMPROAD_DIR;
  const Real ALB_PERROAD_DIR = 0.0799999982118607;
  const Real ALB_PERROAD_DIF = ALB_PERROAD_DIR;
  const Real ALB_ROOF_DIR = 0.254999995231628;
  const Real ALB_ROOF_DIF = ALB_ROOF_DIR;
  const Real ALB_WALL_DIR = 0.200000002980232;
  const Real ALB_WALL_DIF = ALB_WALL_DIR;

  Kokkos::parallel_for(
      "SetParameters", N_LUN, KOKKOS_LAMBDA(const int c) {
        Array2DR8 alb_roof_dir = data_bundle.Roof.BaseAlbedo.dir;
        Array2DR8 alb_roof_dif = data_bundle.Roof.BaseAlbedo.dif;
        Array2DR8 alb_improad_dir = data_bundle.ImperviousRoad.BaseAlbedo.dir;
        Array2DR8 alb_improad_dif = data_bundle.ImperviousRoad.BaseAlbedo.dif;
        Array2DR8 alb_perroad_dir = data_bundle.PerviousRoad.BaseAlbedo.dir;
        Array2DR8 alb_perroad_dif = data_bundle.PerviousRoad.BaseAlbedo.dif;
        Array2DR8 alb_sunwall_dir = data_bundle.SunlitWall.BaseAlbedo.dir;
        Array2DR8 alb_sunwall_dif = data_bundle.SunlitWall.BaseAlbedo.dif;
        Array2DR8 alb_shdwall_dir = data_bundle.ShadedWall.BaseAlbedo.dir;
        Array2DR8 alb_shdwall_dif = data_bundle.ShadedWall.BaseAlbedo.dif;

        for (int ib = 0; ib < N_RAD; ib++) {
          alb_roof_dir(ib, c) = ALB_ROOF_DIR;
          alb_roof_dif(ib, c) = ALB_ROOF_DIF;

          alb_improad_dir(ib, c) = ALB_IMPROAD_DIR;
          alb_improad_dif(ib, c) = ALB_IMPROAD_DIF;

          alb_perroad_dir(ib, c) = ALB_PERROAD_DIR;
          alb_perroad_dif(ib, c) = ALB_PERROAD_DIF;

          alb_sunwall_dir(ib, c) = ALB_WALL_DIR;
          alb_sunwall_dif(ib, c) = ALB_WALL_DIF;

          alb_shdwall_dir(ib, c) = ALB_WALL_DIR;
          alb_shdwall_dif(ib, c) = ALB_WALL_DIF;
        }
      });
}

} // namespace URBANXX
