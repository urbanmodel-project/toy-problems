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

  ALLOCATE_VIEW(geometry.ViewFactorSkyForRoad, Array1DR8, N_LUN);
  ALLOCATE_VIEW(geometry.ViewFactorSkyForWall, Array1DR8, N_LUN);
  ALLOCATE_VIEW(geometry.ViewFactorRoadToWall, Array1DR8, N_LUN);
  ALLOCATE_VIEW(geometry.ViewFactorWallToOtherWall, Array1DR8, N_LUN);

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

void CombinedRoadTypeAllocateViews(int N_LUN, int N_RAD,
                                   CombinedRoadType &combineRoad) {
  RADIATION_TYPE_ALLOCATE_VIEW(combineRoad.DownwellingShortRad, Array2DR8,
                               N_RAD, N_LUN);
  RADIATION_TYPE_ALLOCATE_VIEW(combineRoad.SnowAlbedo, Array2DR8, N_RAD, N_LUN);
  RADIATION_TYPE_ALLOCATE_VIEW(combineRoad.AlbedoWithSnowEffects, Array2DR8,
                               N_RAD, N_LUN);
}

void RoadTypeAllocateViews(int N_LUN, int N_RAD, RoadType &road) {
  RADIATION_TYPE_ALLOCATE_VIEW(road.SnowAlbedo, Array2DR8, N_RAD, N_LUN);
  RADIATION_TYPE_ALLOCATE_VIEW(road.AlbedoWithSnowEffects, Array2DR8, N_RAD,
                               N_LUN);
  RADIATION_TYPE_ALLOCATE_VIEW(road.ReflectedShortRad, Array2DR8, N_RAD, N_LUN);
  RADIATION_TYPE_ALLOCATE_VIEW(road.AbsorbedShortRad, Array2DR8, N_RAD, N_LUN);
}

void WallTypeAllocateViews(int N_LUN, int N_RAD, WallType &wall) {
  RADIATION_TYPE_ALLOCATE_VIEW(wall.DownwellingShortRad, Array2DR8, N_RAD,
                               N_LUN);
  RADIATION_TYPE_ALLOCATE_VIEW(wall.ReflectedShortRad, Array2DR8, N_RAD, N_LUN);
  RADIATION_TYPE_ALLOCATE_VIEW(wall.AbsorbedShortRad, Array2DR8, N_RAD, N_LUN);
}

void RoofTypeAllocateViews(int N_LUN, int N_RAD, RoofType &roof) {
  RADIATION_TYPE_ALLOCATE_VIEW(roof.SnowAlbedo, Array2DR8, N_RAD, N_LUN);
  RADIATION_TYPE_ALLOCATE_VIEW(roof.AlbedoWithSnowEffects, Array2DR8, N_RAD,
                               N_LUN);
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
  CombinedRoadTypeAllocateViews(N_LUN, N_RAD, data_bundle.CombinedRoad);
  WallTypeAllocateViews(N_LUN, N_RAD, data_bundle.SunlitWall);
  WallTypeAllocateViews(N_LUN, N_RAD, data_bundle.ShadedWall);
  CombinedRoadTypeAllocateViews(N_LUN, N_RAD, data_bundle.CombinedRoad);
  RoadTypeAllocateViews(N_LUN, N_RAD, data_bundle.ImperviousRoad);
  RoadTypeAllocateViews(N_LUN, N_RAD, data_bundle.PerviousRoad);

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
}

} // namespace URBANXX
