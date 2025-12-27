program driver_urban_f
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

  call urban_kokkos_initialize()

  cfg%N_LUN = N_LUN
  cfg%N_RAD_BAND = N_RAD_BAND
  cfg%enable_openmp = .true.
  cfg%omp_num_threads = 1

  status = UrbanCreate(cfg, sim)
  if (status /= 0) then
    print *, 'UrbanCreate failed, status=', status
    call urban_kokkos_finalize()
    stop 1
  end if

  status = UrbanInitialize(sim)
  if (status /= 0) then
    print *, 'UrbanInitialize failed, status=', status
    call urban_kokkos_finalize()
    stop 1
  end if

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

  status = UrbanSetInputs(sim, in)
  if (status /= 0) print *, 'UrbanSetInputs status=', status

  status = UrbanStep(sim)
  if (status /= 0) print *, 'UrbanStep status=', status

  out%net_shortwave = make_array_d(sw)
  out%net_longwave  = make_array_d(lw)
  out%surface_flux  = make_array_d(flux)

  status = UrbanGetOutputs(sim, out)
  if (status /= 0) print *, 'UrbanGetOutputs status=', status

  print *, 'net_sw:', sw
  print *, 'net_lw:', lw
  print *, 'flux  :', flux

  status = UrbanDestroy(sim)
  if (status /= 0) print *, 'UrbanDestroy status=', status

  call urban_kokkos_finalize()

end program driver_urban_f
