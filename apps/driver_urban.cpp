#include <mpi.h>
#include <Kokkos_Core.hpp>
#include <UrbanData.hpp>
#include <UrbanDataAllocator.hpp>
#include <UrbanSolarRad.hpp>
#include <UrbanLongwaveRad.hpp>
#include <iostream>

using namespace URBANXX;

int main(int argc, char* argv[]) {

    const int N_LUN = 1000; // Number of Land Units
    const int N_RAD_BAND = 2;    // Number of Radiation Bands (Visible/NIR) [cite: 1591]
    const double DT = 60.0; // Time step (not used, but common context)

    // Initialize MPI
    MPI_Init(&argc, &argv);

    int mpi_rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    Kokkos::initialize(argc, argv);
    {
        // Safest: host-only print guard
        if (mpi_rank == 0) {
        // Ensure only master thread prints
        #pragma omp master
        {
            std::cout << "=== Kokkos Configuration (rank 0 only) ===" << std::endl;
            Kokkos::print_configuration(std::cout);
        }

        const int N_LUN = 3;
        const int N_RAD_BAND = 2;
        UrbanSharedDataBundle simulation_bundle = {
                .N_LUN = N_LUN,
                .N_RAD_BAND = N_RAD_BAND
            };

        UrbanDataAllocator allocator(simulation_bundle);
        allocator.allocate_all_views();
        allocator.initialize_canyon_geometry();
        allocator.initialize_properties();
        allocator.initialize_states();

        UrbanAlbedo albedo_physics(simulation_bundle);
        UrbanLongwave longwave_physics(simulation_bundle);

        if (1) {
            albedo_physics.setSolarInputs();
            albedo_physics.computeIncidentRadiation();
            albedo_physics.computeSnowAlbedo();
            albedo_physics.computeCombinedAlbedo();
            albedo_physics.computeNetSolar();
        }

    }
        
    }

    Kokkos::finalize();
    MPI_Finalize();
    return 0;
}
