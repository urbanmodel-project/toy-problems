// Internal implementation header for Urban C API handle
#ifndef URBAN_PRIVATE_URBANIMPL_H
#define URBAN_PRIVATE_URBANIMPL_H

#include <Urban.h>

#include <private/UrbanDataImpl.h>
#include <private/UrbanDataAllocatorImpl.h>
#include <private/UrbanSolarRadImpl.h>
#include <private/UrbanLongwaveRadImpl.h>
#include <private/UrbanFluxesImpl.h>

#include <memory>

// Opaque handle backing struct used by the C API
struct _p_UrbanType {
  URBANXX::UrbanSharedDataBundle bundle;
  std::unique_ptr<URBANXX::UrbanDataAllocator> allocator;
  std::unique_ptr<URBANXX::UrbanAlbedo> albedo;
  std::unique_ptr<URBANXX::UrbanLongwave> longwave;
  std::unique_ptr<URBANXX::UrbanSurfaceFluxes> fluxes;

  UrbanLogFn logger = nullptr;
  void* logger_ud = nullptr;

  explicit _p_UrbanType(int n_lun, int n_rad)
      : bundle{URBANXX::CanyonGeometryData{},
               URBANXX::AtmosphereInputData{},
               URBANXX::RoofDataType{},
               URBANXX::WallDataType{},
               URBANXX::WallDataType{},
               URBANXX::RoadDataType{},
               URBANXX::RoadDataType{},
               URBANXX::CombinedRoadDataType{},
               n_lun,
               n_rad} {}
};

#endif // URBAN_PRIVATE_URBANIMPL_H