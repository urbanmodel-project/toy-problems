#include <Kokkos_Core.hpp>
#include <Urban.h>
#include <mpi.h>

#include <cstdio>
#include <vector>

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  Kokkos::initialize(argc, argv);
  int rank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  const int N_LUN = 3;
  const int N_RAD_BAND = 2;

  UrbanConfig cfg{};
  cfg.N_LUN = N_LUN;
  cfg.N_RAD_BAND = N_RAD_BAND;

  UrbanType sim = nullptr;
  if (UrbanCreate(&cfg, &sim) != URBAN_SUCCESS) {
    if (rank == 0)
      std::fprintf(stderr, "UrbanCreate failed\n");
    Kokkos::finalize();
    MPI_Finalize();
    return 1;
  }

  // Set options via API rather than config struct fields
  {
    UrbanErrorCode rc = UrbanInitialize(sim);
    if (rc != URBAN_SUCCESS) {
      if (rank == 0)
        std::fprintf(stderr, "UrbanInitialize failed: %s\n",
                     UrbanGetErrorString(rc));
      Kokkos::finalize();
      MPI_Finalize();
      return 1;
    }
  }

  std::vector<double> solar(N_LUN, 400.0);
  std::vector<double> longwave(N_LUN, 300.0);
  std::vector<double> air_temp(N_LUN, 290.0);
  std::vector<double> wind(N_LUN, 2.0);

  UrbanInputs in{};
  in.solar_down = {solar.data(), solar.size()};
  in.longwave_down = {longwave.data(), longwave.size()};
  in.air_temp = {air_temp.data(), air_temp.size()};
  in.wind_speed = {wind.data(), wind.size()};
  {
    UrbanErrorCode rc = UrbanSetInputs(sim, &in);
    if (rc != URBAN_SUCCESS) {
      if (rank == 0)
        std::fprintf(stderr, "UrbanSetInputs failed: %s\n",
                     UrbanGetErrorString(rc));
      Kokkos::finalize();
      MPI_Finalize();
      return 1;
    }
  }

  {
    UrbanErrorCode rc = UrbanStep(sim);
    if (rc != URBAN_SUCCESS) {
      if (rank == 0)
        std::fprintf(stderr, "UrbanStep failed: %s\n", UrbanGetErrorString(rc));
      Kokkos::finalize();
      MPI_Finalize();
      return 1;
    }
  }

  std::vector<double> sw(N_LUN), lw(N_LUN), flux(N_LUN);
  UrbanOutputs out{};
  out.net_shortwave = {sw.data(), sw.size()};
  out.net_longwave = {lw.data(), lw.size()};
  out.surface_flux = {flux.data(), flux.size()};
  {
    UrbanErrorCode rc = UrbanGetOutputs(sim, &out);
    if (rc != URBAN_SUCCESS) {
      if (rank == 0)
        std::fprintf(stderr, "UrbanGetOutputs failed: %s\n",
                     UrbanGetErrorString(rc));
      Kokkos::finalize();
      MPI_Finalize();
      return 1;
    }
  }

  if (rank == 0) {
    std::printf("net_sw: ");
    for (int i = 0; i < N_LUN; ++i)
      std::printf("%g ", sw[i]);
    std::printf("\nnet_lw: ");
    for (int i = 0; i < N_LUN; ++i)
      std::printf("%g ", lw[i]);
    std::printf("\nflux: ");
    for (int i = 0; i < N_LUN; ++i)
      std::printf("%g ", flux[i]);
    std::printf("\n");
  }

  {
    UrbanErrorCode rc = UrbanDestroy(&sim);
    if (rc != URBAN_SUCCESS) {
      if (rank == 0)
        std::fprintf(stderr, "UrbanDestroy failed: %s\n",
                     UrbanGetErrorString(rc));
      Kokkos::finalize();
      MPI_Finalize();
      return 1;
    }
  }
  Kokkos::finalize();
  MPI_Finalize();
  return 0;
}
