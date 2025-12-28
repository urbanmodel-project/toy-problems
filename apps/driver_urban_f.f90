program driver_urban_f
#include <finclude/UrbanMacros.h>
  use, intrinsic :: iso_c_binding
  use urban
  implicit none

  integer(c_int), parameter :: N_LUN = 3
  integer(c_int), parameter :: N_RAD_BAND = 2

  type(UrbanConfig_c) :: cfg
  type(UrbanType) :: sim
  integer(c_int) :: status
  real(c_double) :: solar(N_LUN), longwave(N_LUN), air_temp(N_LUN), wind(N_LUN)
  integer :: i
  type(UrbanInputs_c) :: in
  real(c_double) :: sw(N_LUN), lw(N_LUN), flux(N_LUN)
  type(UrbanOutputs_c) :: out

  call UrbanKokkosInitialize()

  cfg%N_LUN = N_LUN
  cfg%N_RAD_BAND = N_RAD_BAND
  cfg%enable_openmp = .true.
  cfg%omp_num_threads = 1

  CallA(UrbanCreate(cfg, sim, status))

  CallA(UrbanInitialize(sim, status))

  do i = 1, N_LUN
    solar(i) = 400.0d0
    longwave(i) = 300.0d0
    air_temp(i) = 290.0d0
    wind(i) = 2.0d0
  end do

  in%solar_down = make_array_d(solar)
  in%longwave_down = make_array_d(longwave)
  in%air_temp = make_array_d(air_temp)
  in%wind_speed = make_array_d(wind)

  CallA(UrbanSetInputs(sim, in, status))

  CallA(UrbanStep(sim, status))

  out%net_shortwave = make_array_d(sw)
  out%net_longwave  = make_array_d(lw)
  out%surface_flux  = make_array_d(flux)

  CallA(UrbanGetOutputs(sim, out, status))

  print *, 'net_sw:', sw
  print *, 'net_lw:', lw
  print *, 'flux  :', flux

  CallA(UrbanDestroy(sim, status))

  call UrbanKokkosFinalize()

end program driver_urban_f
