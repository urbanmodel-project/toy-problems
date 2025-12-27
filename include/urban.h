#ifndef URBAN_H
#define URBAN_H

/* C-only includes */
#include <stddef.h> /* size_t */
#include <stdint.h> /* uint32_t, int32_t */
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Versioning */
#define URBAN_VERSION_MAJOR 0
#define URBAN_VERSION_MINOR 1
#define URBAN_VERSION_PATCH 0

/* Visibility (ok for static libs; future-proof for shared) */
#if defined(_WIN32)
  #if defined(URBAN_BUILD_SHARED)
    #define URBAN_EXTERN __declspec(dllexport)
  #elif defined(URBAN_USE_SHARED)
    #define URBAN_EXTERN __declspec(dllimport)
  #else
    #define URBAN_EXTERN
  #endif
#else
  #if defined(URBAN_BUILD_SHARED)
    #define URBAN_EXTERN __attribute__((visibility("default")))
  #else
    #define URBAN_EXTERN
  #endif
#endif

/* Error/status codes */
typedef enum UrbanErrorCode {
  URBAN_SUCCESS = 0,
  URBAN_ERR_INVALID_ARGUMENT = 1,
  URBAN_ERR_NOT_INITIALIZED = 2,
  URBAN_ERR_INTERNAL = 3
} UrbanErrorCode;

/* Opaque handle to a simulation instance */
typedef struct _p_UrbanType* UrbanType;

/* Optional logging callback */
typedef void (*UrbanLogFn)(int level, const char* message, void* user_data);

/* Basic numeric array views for C ABI stability */
typedef struct UrbanArrayD {
  double* data;
  size_t  size;
} UrbanArrayD;

typedef struct UrbanArrayI32 {
  int32_t* data;
  size_t   size;
} UrbanArrayI32;

/* Public configuration (minimal, extend as needed) */
typedef struct UrbanConfig {
  int32_t N_LUN;         /* land units */
  int32_t N_RAD_BAND;    /* radiation bands */
  /* Optional tuning flags; keep ABI simple */
  bool    enable_openmp;
  int32_t omp_num_threads;
} UrbanConfig;

/* Inputs provided per step or per initialization */
typedef struct UrbanInputs {
  /* Example environmental forcings; expand as needed */
  UrbanArrayD solar_down;     /* incoming shortwave */
  UrbanArrayD longwave_down;  /* incoming longwave */
  UrbanArrayD air_temp;       /* ambient air temperature */
  UrbanArrayD wind_speed;     /* wind speed */
  /* Add additional input arrays or scalars here */
} UrbanInputs;

/* Outputs retrieved after step() */
typedef struct UrbanOutputs {
  UrbanArrayD net_shortwave;  /* per unit */
  UrbanArrayD net_longwave;   /* per unit */
  UrbanArrayD surface_flux;   /* per unit aggregated flux */
  /* Add more outputs as needed */
} UrbanOutputs;

/* Lifecycle */
URBAN_EXTERN UrbanErrorCode UrbanCreate(const UrbanConfig* cfg, UrbanType* out);
URBAN_EXTERN UrbanErrorCode UrbanDestroy(UrbanType* handle);

/* Optional utilities */
URBAN_EXTERN UrbanErrorCode UrbanSetLogger(UrbanType handle, UrbanLogFn fn, void* user_data);
URBAN_EXTERN const char* UrbanGetErrorString(UrbanErrorCode err);

/* Initialization and data binding */
URBAN_EXTERN UrbanErrorCode UrbanInitialize(UrbanType handle);
/* Bind or update inputs; arrays are borrowed for duration of call unless documented otherwise */
URBAN_EXTERN UrbanErrorCode UrbanSetInputs(UrbanType handle, const UrbanInputs* in);

/* Advance physics; internal orchestrates albedo, longwave, fluxes */
URBAN_EXTERN UrbanErrorCode UrbanStep(UrbanType handle);

/* Retrieve outputs; caller supplies buffers or receives borrowed views */
URBAN_EXTERN UrbanErrorCode UrbanGetOutputs(UrbanType handle, UrbanOutputs* out);

/* Optional: simple option setters to avoid ABI churn */
URBAN_EXTERN UrbanErrorCode UrbanSetOptionInt(UrbanType handle, const char* name, int value);
URBAN_EXTERN UrbanErrorCode UrbanSetOptionDouble(UrbanType handle, const char* name, double value);
URBAN_EXTERN UrbanErrorCode UrbanSetOptionBool(UrbanType handle, const char* name, bool value);

/* Optional: copy-out helpers if you prefer owned outputs */
URBAN_EXTERN UrbanErrorCode UrbanCopyOutputs(UrbanType handle,
                                       double* net_sw, size_t net_sw_size,
                                       double* net_lw, size_t net_lw_size,
                                       double* flux,  size_t flux_size);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* URBAN_H */
