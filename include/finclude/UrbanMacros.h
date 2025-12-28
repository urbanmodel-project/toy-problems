#ifndef URBAN_MACROS_F90
#define URBAN_MACROS_F90

#define URBAN_SUCCESS 0
#define CallA(stmt) \
  call stmt; \
  if (status /= 0) stop 1

#endif /* URBAN_MACROS_F90 */
