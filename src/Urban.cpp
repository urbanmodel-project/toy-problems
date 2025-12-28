#include <Urban.h>

#include <private/UrbanDataAllocatorImpl.h>
#include <private/UrbanDataImpl.h>
#include <private/UrbanFluxesImpl.h>
#include <private/UrbanLongwaveRadImpl.h>
#include <private/UrbanSolarRadImpl.h>

#include <Kokkos_Core.hpp>
#include <memory>
#include <stdexcept>
#include <vector>

#include <private/UrbanImpl.h>

using namespace URBANXX;

static void urban_log(_p_UrbanType *h, int level, const char *msg) {
  if (h && h->logger)
    h->logger(level, msg, h->logger_ud);
}

extern "C" {

const char *UrbanGetErrorString(UrbanErrorCode err) {
  switch (err) {
  case URBAN_SUCCESS:
    return "success";
  case URBAN_ERR_INVALID_ARGUMENT:
    return "invalid argument";
  case URBAN_ERR_NOT_INITIALIZED:
    return "not initialized";
  case URBAN_ERR_INTERNAL:
    return "internal error";
  default:
    return "unknown error";
  }
}

UrbanErrorCode UrbanSetLogger(UrbanType handle, UrbanLogFn fn,
                              void *user_data) {
  if (!handle)
    return URBAN_ERR_INVALID_ARGUMENT;
  handle->logger = fn;
  handle->logger_ud = user_data;
  return URBAN_SUCCESS;
}

UrbanErrorCode UrbanCreate(const UrbanConfig *cfg, UrbanType *out) {
  if (!cfg || !out)
    return URBAN_ERR_INVALID_ARGUMENT;
  if (cfg->N_LUN <= 0 || cfg->N_RAD_BAND <= 0)
    return URBAN_ERR_INVALID_ARGUMENT;
  try {
    _p_UrbanType *h = new _p_UrbanType(cfg->N_LUN, cfg->N_RAD_BAND);
    h->allocator = std::make_unique<UrbanDataAllocator>(h->bundle);
    // Defer physics construction until after allocate_all_views
    h->albedo = nullptr;
    h->longwave = nullptr;
    h->fluxes = nullptr;
    *out = h;
    return URBAN_SUCCESS;
  } catch (...) {
    return URBAN_ERR_INTERNAL;
  }
}

UrbanErrorCode UrbanDestroy(UrbanType *handle) {
  if (!handle || !*handle)
    return URBAN_ERR_INVALID_ARGUMENT;
  try {
    delete *handle;
    *handle = nullptr;
    return URBAN_SUCCESS;
  } catch (...) {
    return URBAN_ERR_INTERNAL;
  }
}

UrbanErrorCode UrbanInitialize(UrbanType handle) {
  if (!handle || !handle->allocator)
    return URBAN_ERR_INVALID_ARGUMENT;
  try {
    handle->allocator->allocate_all_views();
    handle->allocator->initialize_canyon_geometry();
    handle->allocator->initialize_properties();
    handle->allocator->initialize_states();
    // Now that views are allocated, construct physics modules
    handle->albedo = std::make_unique<UrbanAlbedo>(handle->bundle);
    handle->longwave = std::make_unique<UrbanLongwave>(handle->bundle);
    handle->fluxes = std::make_unique<UrbanSurfaceFluxes>(handle->bundle);
    urban_log(handle, 1, "Urban initialized");
    return URBAN_SUCCESS;
  } catch (...) {
    return URBAN_ERR_INTERNAL;
  }
}

static void copy_1d_to_dual(const UrbanArrayD &arr, HostArray1DR8 &host,
                            Array1DR8 &dev) {
  if (!arr.data)
    return;
  const size_t n = arr.size;
  if (host.extent(0) < (int)n)
    return; // simplistic bounds guard
  for (size_t i = 0; i < n; ++i)
    host(i) = static_cast<R8>(arr.data[i]);
  Kokkos::deep_copy(dev, host);
}

UrbanErrorCode UrbanSetInputs(UrbanType handle, const UrbanInputs *in) {
  if (!handle || !in)
    return URBAN_ERR_INVALID_ARGUMENT;
  try {
    // Map a minimal subset of inputs to available fields.
    // Shortwave: distribute to direct beam band 0; set diffuse to zero
    {
      auto &H = handle->bundle.input.SdirHorizH;
      auto &D = handle->bundle.input.SdirHoriz;
      // If rank-2, fill [i][0] with solar_down and [i][1]=0
      const size_t n = in->solar_down.size;
      // Zero the host diffuse array to avoid copying uninitialized data
      Kokkos::deep_copy(handle->bundle.input.SdifHorizH, 0.0);
      Kokkos::deep_copy(H, 0.0);
      for (size_t i = 0; i < n; ++i) {
        H(i, 0) = static_cast<R8>(in->solar_down.data[i]);
        H(i, 1) = 0.0;
      }
      Kokkos::deep_copy(D, H);
      // Set diffuse to zero
      Kokkos::deep_copy(handle->bundle.input.SdifHoriz,
                        handle->bundle.input.SdifHorizH);
    }

    // Longwave downwelling
    copy_1d_to_dual(in->longwave_down, handle->bundle.input.DownwellingLongRadH,
                    handle->bundle.input.DownwellingLongRad);
    // Air temperature
    copy_1d_to_dual(in->air_temp, handle->bundle.input.ForcTempH,
                    handle->bundle.input.ForcTemp);
    // Wind speed (map to U, set V=0)
    copy_1d_to_dual(in->wind_speed, handle->bundle.input.ForcWindUH,
                    handle->bundle.input.ForcWindU);
    Kokkos::deep_copy(handle->bundle.input.ForcWindV,
                      handle->bundle.input.ForcWindVH);

    urban_log(handle, 1, "Inputs set");
    return URBAN_SUCCESS;
  } catch (...) {
    return URBAN_ERR_INTERNAL;
  }
}

UrbanErrorCode UrbanStep(UrbanType handle) {
  if (!handle)
    return URBAN_ERR_INVALID_ARGUMENT;
  if (!handle->albedo || !handle->longwave || !handle->fluxes)
    return URBAN_ERR_INVALID_ARGUMENT;
  try {
    handle->longwave->setLongwaveInputs();
    handle->longwave->computeNetLongwave();

    handle->albedo->setSolarInputs();
    handle->albedo->computeIncidentRadiation();
    handle->albedo->computeSnowAlbedo();
    handle->albedo->computeCombinedAlbedo();
    handle->albedo->computeNetSolar();

    handle->fluxes->computeSurfaceFluxes();

    urban_log(handle, 1, "Step completed");
    return URBAN_SUCCESS;
  } catch (const std::runtime_error &e) {
    urban_log(handle, 0, e.what());
    return URBAN_ERR_INVALID_ARGUMENT;
  } catch (...) {
    return URBAN_ERR_INTERNAL;
  }
}

UrbanErrorCode UrbanGetOutputs(UrbanType handle, UrbanOutputs *out) {
  if (!handle || !out)
    return URBAN_ERR_INVALID_ARGUMENT;
  try {
    const int n = handle->bundle.N_LUN;
    // Ensure buffers provided are large enough
    if (!out->net_shortwave.data || out->net_shortwave.size < (size_t)n)
      return URBAN_ERR_INVALID_ARGUMENT;
    if (!out->net_longwave.data || out->net_longwave.size < (size_t)n)
      return URBAN_ERR_INVALID_ARGUMENT;
    if (!out->surface_flux.data || out->surface_flux.size < (size_t)n)
      return URBAN_ERR_INVALID_ARGUMENT;

    // Copy representative outputs using proper host mirrors
    {
      auto dev = handle->bundle.ImperviousRoad.NetLongRad;
      auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), dev);
      for (int i = 0; i < n; ++i)
        out->net_longwave.data[i] = static_cast<double>(host(i));
    }

    {
      auto dev = handle->bundle.CombinedRoad.DownwellingShortRad;
      auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), dev);
      for (int i = 0; i < n; ++i)
        out->net_shortwave.data[i] = static_cast<double>(host(i, 0, 0));
    }

    {
      auto dev = handle->bundle.input.ForcTemp;
      auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), dev);
      for (int i = 0; i < n; ++i)
        out->surface_flux.data[i] = static_cast<double>(host(i));
    }

    urban_log(handle, 1, "Outputs retrieved");
    return URBAN_SUCCESS;
  } catch (...) {
    return URBAN_ERR_INTERNAL;
  }
}

UrbanErrorCode UrbanSetOptionInt(UrbanType, const char *, int) {
  return URBAN_SUCCESS;
}
UrbanErrorCode UrbanSetOptionDouble(UrbanType, const char *, double) {
  return URBAN_SUCCESS;
}
UrbanErrorCode UrbanSetOptionBool(UrbanType, const char *, bool) {
  return URBAN_SUCCESS;
}

UrbanErrorCode UrbanCopyOutputs(UrbanType handle, double *net_sw,
                                size_t net_sw_size, double *net_lw,
                                size_t net_lw_size, double *flux,
                                size_t flux_size) {
  UrbanOutputs out{
      {net_sw, net_sw_size}, {net_lw, net_lw_size}, {flux, flux_size}};
  return UrbanGetOutputs(handle, &out);
}

} // extern "C"
