#ifndef URBAN_MACROS_F90
#define URBAN_MACROS_F90

#define URBAN_SUCCESS 0
#define URBAN_CALL_OR_STOP(stmt) \
  call stmt
#define CallA(stmt) \
  URBAN_CALL_OR_STOP(stmt); \
  if (status /= URBAN_SUCCESS) then; \
    print *, 'API call failed: ', status; \
    stop 1; \
  end if

#endif /* URBAN_MACROS_F90 */
