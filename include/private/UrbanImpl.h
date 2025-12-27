// Internal implementation header for Urban C API handle
#ifndef URBAN_PRIVATE_URBANIMPL_H
#define URBAN_PRIVATE_URBANIMPL_H

#include <urban.h>

#include <private/UrbanData.hpp>
#include <private/UrbanDataAllocator.hpp>
#include <private/UrbanSolarRad.hpp>
#include <private/UrbanLongwaveRad.hpp>
#include <private/UrbanFluxes.hpp>

#include <memory>

// Opaque handle backing struct used by the C API
struct UrbanOpaque_ {
  URBANXX::UrbanSharedDataBundle bundle;
  std::unique_ptr<URBANXX::UrbanDataAllocator> allocator;
  std::unique_ptr<URBANXX::UrbanAlbedo> albedo;
  std::unique_ptr<URBANXX::UrbanLongwave> longwave;
  std::unique_ptr<URBANXX::UrbanSurfaceFluxes> fluxes;

  UrbanLogFn logger = nullptr;
  void* logger_ud = nullptr;

  explicit UrbanOpaque_(int n_lun, int n_rad)
      : bundle{.geometry{}, .input{}, .Roof{}, .SunlitWall{}, .ShadedWall{},
               .ImperviousRoad{}, .PerviousRoad{}, .CombinedRoad{},
               .N_LUN = n_lun, .N_RAD_BAND = n_rad} {}
};

#endif // URBAN_PRIVATE_URBANIMPL_H