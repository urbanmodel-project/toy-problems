#ifndef URBAN_MACROS_F90
#define URBAN_MACROS_F90

#define CallA(stmt) call stmt; if (status /= 0) then; print *, 'CallA failed: ', status; stop 1; end if

#endif /* URBAN_MACROS_F90 */
