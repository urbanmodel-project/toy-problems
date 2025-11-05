#ifndef URBAN_DATA_HPP
#define URBAN_DATA_HPP

#include <DataTypes.hpp>
#include <Kokkos_Core.hpp>
#include <cmath>

namespace URBANXX {

// Type Definitions
using ExecutionSpace = Kokkos::DefaultExecutionSpace;

// --- Fortran Type: urbanparams_vars and Canyon Geometry ---
struct CanyonGeometryData {
  Array1DR8 canyon_hwr; // ratio of building height to street width
  Array1DR8 vf_sr;      // view factor of sky for road
  Array1DR8 vf_sw;      // view factor of sky for one wall
  Array1DR8 vf_rw;      // view factor from road to wall
  Array1DR8 vf_ww;      // view factor from wall to opposing wall
};

// --- Local Variables and Initial Inputs ---
struct SolarInputData {
  Array1DR8 coszen; // cosine solar zenith angle (coszen)
  Array2DR8 sdir_h; // direct beam solar radiation on horizontal surface (sdir)
  Array2DR8 sdif_h; // diffuse solar radiation on horizontal surface (sdif)
  Array1DR8 frac_sno; // fraction of ground covered by snow (col_ws%frac_sno)
};

// --- Intermediate/Output Fluxes (Shared between subroutines) ---
struct AlbedoFluxes {
  // Incident Fluxes (Output of incident_direct/diffuse)
  Array2DR8 sdir_road;      // direct incident on road
  Array2DR8 sdif_road;      // diffuse incident on road
  Array2DR8 sdir_sunwall;   // direct incident on sunlit wall
  Array2DR8 sdif_sunwall;   // diffuse incident on sunlit wall
  Array2DR8 sdir_shadewall; // direct incident on shaded wall
  Array2DR8 sdif_shadewall; // diffuse incident on shaded wall

  // Snow-corrected albedos (Output of SnowAlbedo)
  Array2DR8 alb_roof_dir_s;
  Array2DR8 alb_improad_dir_s;
  Array2DR8 alb_perroad_dir_s;
  // ... diffuse versions omitted for brevity ...

  // Absorbed Fluxes (Output of net_solar) (solarabs_vars)
  Array2DR8 sabs_improad_dir; // direct solar absorbed by impervious road
  Array2DR8 sabs_sunwall_dir; // direct solar absorbed by sunwall

  // Reflected Fluxes (Output of net_solar) (Local sref_... variables)
  Array2DR8 sref_roof_dir;    // direct solar reflected by roof (sref_roof_dir)
  Array2DR8 sref_sunwall_dir; // direct solar reflected by sunwall
};

// --- The Master Bundle passed between all classes ---
struct UrbanSharedDataBundle {
  CanyonGeometryData geometry;
  SolarInputData input;
  AlbedoFluxes fluxes;
  const int N_LUN;
  const int N_RAD;
};
} // namespace URBANXX
#endif