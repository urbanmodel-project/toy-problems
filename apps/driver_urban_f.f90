program driver_urban_f
#include <finclude/UrbanMacros.h>
  use, intrinsic :: iso_c_binding
  use urban
  implicit none

  integer(kind=c_int), parameter :: N_LUN = 3
  integer(kind=c_int), parameter :: N_RAD_BAND = 2

  type(UrbanConfig) :: cfg
  type(UrbanType) :: sim
  integer(kind=c_int) :: status
  real(kind=c_double) :: solar(N_LUN), longwave(N_LUN), air_temp(N_LUN), wind(N_LUN)
  integer :: i
  type(UrbanInputs) :: in
  real(kind=c_double) :: sw(N_LUN), lw(N_LUN), flux(N_LUN)
  type(UrbanOutputs) :: out

  call UrbanKokkosInitialize()

  CallA(UrbanConfigSetNLun(cfg, N_LUN, status))
  CallA(UrbanConfigSetNRadBand(cfg, N_RAD_BAND, status))

  CallA(UrbanCreate(cfg, sim, status))

  CallA(UrbanInitialize(sim, status))

  do i = 1, N_LUN
    solar(i) = 400.0d0
    longwave(i) = 300.0d0
    air_temp(i) = 290.0d0
    wind(i) = 2.0d0
  end do

  CallA(UrbanInputsSetSolarDown(in, solar, status))
  CallA(UrbanInputsSetLongwaveDown(in, longwave, status))
  CallA(UrbanInputsSetAirTemp(in, air_temp, status))
  CallA(UrbanInputsSetWindSpeed(in, wind, status))

  CallA(UrbanSetInputs(sim, in, status))

  CallA(UrbanStep(sim, status))

  CallA(UrbanOutputsSetNetShortwave(out, sw, status))
  CallA(UrbanOutputsSetNetLongwave(out, lw, status))
  CallA(UrbanOutputsSetSurfaceFlux(out, flux, status))

  CallA(UrbanGetOutputs(sim, out, status))

  print *, 'net_sw:', sw
  print *, 'net_lw:', lw
  print *, 'flux  :', flux

  CallA(UrbanDestroy(sim, status))

  call UrbanKokkosFinalize()

end program driver_urban_f
