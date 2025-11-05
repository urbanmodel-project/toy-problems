#ifndef URBAN_DATA_HPP
#define URBAN_DATA_HPP

#include <Kokkos_Core.hpp>
#include <DataTypes.hpp>
#include <cmath>

namespace URBANXX {

  // Type Definitions
  using ExecutionSpace = Kokkos::DefaultExecutionSpace;

  // --- Fortran Type: urbanparams_vars and Canyon Geometry ---
  struct CanyonGeometryData {
      Array1DR8 canyon_hwr;      // ratio of building height to street width (lun_pp%canyon_hwr) [cite: 1478]
      Array1DR8 vf_sr;           // view factor of sky for road (urbanparams_vars%vf_sr) [cite: 1136, 1340]
      Array1DR8 vf_sw;           // view factor of sky for one wall (urbanparams_vars%vf_sw) [cite: 1136]
      Array1DR8 vf_rw;           // view factor from road to wall (urbanparams_vars%vf_rw) [cite: 1340]
      Array1DR8 vf_ww;           // view factor from wall to opposing wall (urbanparams_vars%vf_ww) [cite: 1340]
  };

  // --- Local Variables and Initial Inputs ---
  struct SolarInputData {
      Array1DR8 coszen;          // cosine solar zenith angle (coszen) [cite: 1444]
      Array2DR8 sdir_h;          // direct beam solar radiation on horizontal surface (sdir) [cite: 1446]
      Array2DR8 sdif_h;          // diffuse solar radiation on horizontal surface (sdif) [cite: 1447]
      Array1DR8 frac_sno;        // fraction of ground covered by snow (col_ws%frac_sno) [cite: 1480]
  };

  // --- Intermediate/Output Fluxes (Shared between subroutines) ---
  struct AlbedoFluxes {
      // Incident Fluxes (Output of incident_direct/diffuse)
      Array2DR8 sdir_road;       // direct incident on road [cite: 1449]
      Array2DR8 sdif_road;       // diffuse incident on road [cite: 1450]
      Array2DR8 sdir_sunwall;    // direct incident on sunlit wall [cite: 1451]
      Array2DR8 sdif_sunwall;    // diffuse incident on sunlit wall [cite: 1452]
      Array2DR8 sdir_shadewall;  // direct incident on shaded wall [cite: 1453]
      Array2DR8 sdif_shadewall;  // diffuse incident on shaded wall [cite: 1454]

      // Snow-corrected albedos (Output of SnowAlbedo)
      Array2DR8 alb_roof_dir_s;
      Array2DR8 alb_improad_dir_s;
      Array2DR8 alb_perroad_dir_s;
      // ... diffuse versions omitted for brevity ...

      // Absorbed Fluxes (Output of net_solar) (solarabs_vars)
      Array2DR8 sabs_improad_dir; // direct solar absorbed by impervious road (solarabs_vars%sabs_improad_dir_lun) [cite: 1496]
      Array2DR8 sabs_sunwall_dir; // direct solar absorbed by sunwall (solarabs_vars%sabs_sunwall_dir_lun) [cite: 1492]
      // ... diffuse and other absorbed fluxes omitted for brevity ...
      
      // Reflected Fluxes (Output of net_solar) (Local sref_... variables)
      Array2DR8 sref_roof_dir;    // direct solar reflected by roof (sref_roof_dir) [cite: 1466]
      Array2DR8 sref_sunwall_dir; // direct solar reflected by sunwall (sref_sunwall_dir) [cite: 1468]
      // ... diffuse and other reflected fluxes omitted for brevity ...
  };

  // --- The Master Bundle passed between all classes ---
  struct UrbanSharedDataBundle {
      CanyonGeometryData geometry;
      SolarInputData input;
      AlbedoFluxes fluxes;
      const int N_LUN;
      const int N_RAD;
  };
}
#endif