// Public C API for the Urban library
#ifndef URBAN_H
#define URBAN_H

#include <stddef.h> // size_t

#ifdef __cplusplus
extern "C" {
#endif

// Visibility / linkage macro for C API functions
#ifdef __cplusplus
#define URBAN_EXTERN extern "C"
#else
#define URBAN_EXTERN extern
#endif

#ifndef __cplusplus
#include <stdbool.h>
#endif

// Opaque handle type (pointer to private struct)
typedef struct _p_UrbanType* UrbanType;

// Error codes for API calls
typedef enum {
	URBAN_SUCCESS = 0,
	URBAN_ERR_INVALID_ARGUMENT = 1,
	URBAN_ERR_NOT_INITIALIZED = 2,
	URBAN_ERR_INTERNAL = 3
} UrbanErrorCode;

// Simple array view used for inputs/outputs
typedef struct {
	double* data;
	size_t  size;
} UrbanArrayD;

// Configuration for creating an Urban handle
typedef struct {
	int N_LUN;      // number of land units
	int N_RAD_BAND; // number of shortwave radiation bands
} UrbanConfig;

// Input fields provided to the model
typedef struct {
	UrbanArrayD solar_down;    // downwelling shortwave on horizontal surface
	UrbanArrayD longwave_down; // downwelling longwave radiation
	UrbanArrayD air_temp;      // air temperature (K)
	UrbanArrayD wind_speed;    // wind speed magnitude (m/s)
} UrbanInputs;

// Output fields retrieved from the model
typedef struct {
	UrbanArrayD net_shortwave; // representative shortwave metric
	UrbanArrayD net_longwave;  // representative longwave metric
	UrbanArrayD surface_flux;  // representative surface flux metric
} UrbanOutputs;

// Logger callback function type
typedef void (*UrbanLogFn)(int level, const char* message, void* user_data);

// Utility
URBAN_EXTERN const char* UrbanGetErrorString(UrbanErrorCode err);

// Lifecycle
URBAN_EXTERN UrbanErrorCode UrbanCreate(const UrbanConfig* cfg, UrbanType* out);
URBAN_EXTERN UrbanErrorCode UrbanDestroy(UrbanType* handle);
URBAN_EXTERN UrbanErrorCode UrbanInitialize(UrbanType handle);

// I/O
URBAN_EXTERN UrbanErrorCode UrbanSetInputs(UrbanType handle, const UrbanInputs* in);
URBAN_EXTERN UrbanErrorCode UrbanGetOutputs(UrbanType handle, UrbanOutputs* out);

// Stepping
URBAN_EXTERN UrbanErrorCode UrbanStep(UrbanType handle);

// Options
URBAN_EXTERN UrbanErrorCode UrbanSetOptionInt(UrbanType handle, const char* name, int value);
URBAN_EXTERN UrbanErrorCode UrbanSetOptionDouble(UrbanType handle, const char* name, double value);
URBAN_EXTERN UrbanErrorCode UrbanSetOptionBool(UrbanType handle, const char* name, bool value);

// Logging
URBAN_EXTERN UrbanErrorCode UrbanSetLogger(UrbanType handle, UrbanLogFn fn, void* user_data);

// Convenience: copy outputs directly into buffers
URBAN_EXTERN UrbanErrorCode UrbanCopyOutputs(UrbanType handle,
																						 double* net_sw, size_t net_sw_size,
																						 double* net_lw, size_t net_lw_size,
																						 double* flux,  size_t flux_size);

#ifdef __cplusplus
} // extern "C"
#endif

#endif // URBAN_H
