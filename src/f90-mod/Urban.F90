module urban
  use, intrinsic :: iso_c_binding
  implicit none

  ! Status codes
  integer(kind=c_int), parameter :: URBAN_SUCCESS = 0

  ! Fortran-friendly opaque wrapper around C UrbanType
  type :: UrbanType
    type(c_ptr) :: c_ptr
  end type UrbanType

  ! Opaque handle type for C pointer
  type, bind(C) :: UrbanArrayD_c
    type(c_ptr) :: data
    integer(kind=c_size_t) :: size
  end type UrbanArrayD_c

  ! UrbanConfig: C-interoperable configuration structure
  ! WARNING: Do not directly access or modify member variables (N_LUN,
  ! N_RAD_BAND). Use the provided setter functions:
  !   - UrbanConfigSetNLun(cfg, n_lun, status)
  !   - UrbanConfigSetNRadBand(cfg, n_rad_band, status)
  ! Direct field access bypasses proper validation and may cause errors.
  type, bind(C) :: UrbanConfig
    integer(kind=c_int) :: N_LUN
    integer(kind=c_int) :: N_RAD_BAND
  end type UrbanConfig

  ! UrbanInputs: C-interoperable input data structure
  ! WARNING: Do not directly access or modify member variables (solar_down,
  ! longwave_down, air_temp, wind_speed). Use the provided setter functions:
  !   - UrbanInputsSetSolarDown(in, array, status)
  !   - UrbanInputsSetLongwaveDown(in, array, status)
  !   - UrbanInputsSetAirTemp(in, array, status)
  !   - UrbanInputsSetWindSpeed(in, array, status)
  ! Direct field access bypasses proper array handling and may cause errors.
  type, bind(C) :: UrbanInputs
    type(UrbanArrayD_c) :: solar_down
    type(UrbanArrayD_c) :: longwave_down
    type(UrbanArrayD_c) :: air_temp
    type(UrbanArrayD_c) :: wind_speed
  end type UrbanInputs

  ! UrbanOutputs: C-interoperable output data structure
  ! WARNING: Do not directly access or modify member variables (net_shortwave,
  ! net_longwave, surface_flux). Use the provided getter functions before
  ! calling UrbanGetOutputs:
  !   - UrbanOutputsGetNetShortwave(out, array, status)
  !   - UrbanOutputsGetNetLongwave(out, array, status)
  !   - UrbanOutputsGetSurfaceFlux(out, array, status)
  ! Direct field access bypasses proper array handling and may cause errors.
  type, bind(C) :: UrbanOutputs
    type(UrbanArrayD_c) :: net_shortwave
    type(UrbanArrayD_c) :: net_longwave
    type(UrbanArrayD_c) :: surface_flux
  end type UrbanOutputs

  interface
    ! Kokkos C wrapper interfaces
    subroutine UrbanKokkosInitialize() bind(C, name="UrbanKokkosInitialize")
    end subroutine UrbanKokkosInitialize

    subroutine UrbanKokkosFinalize() bind(C, name="UrbanKokkosFinalize")
    end subroutine UrbanKokkosFinalize

    function UrbanGetErrorString(err) bind(C, name="UrbanGetErrorString") result(msg)
      import :: c_int, c_ptr
      integer(kind=c_int), value :: err
      type(c_ptr) :: msg
    end function UrbanGetErrorString

    function UrbanSetLogger_c(handle, fn, user_data) bind(C, name="UrbanSetLogger") result(status)
      import :: c_int, c_ptr, c_funptr
      type(c_ptr), value :: handle
      type(c_funptr), value :: fn
      type(c_ptr), value :: user_data
      integer(kind=c_int) :: status
    end function UrbanSetLogger_c

    function UrbanCreate_c(cfg, out) bind(C, name="UrbanCreate") result(status)
      import :: UrbanConfig, c_int, c_ptr
      type(UrbanConfig), intent(in) :: cfg
      type(c_ptr), intent(out) :: out
      integer(kind=c_int) :: status
    end function UrbanCreate_c

    function UrbanDestroy_c(handle) bind(C, name="UrbanDestroy") result(status)
      import :: c_ptr, c_int
      ! Pointer-to-pointer: C signature is UrbanDestroy(UrbanType* handle)
      ! Pass-by-reference so C can null out the pointer.
      type(c_ptr), intent(inout) :: handle
      integer(kind=c_int) :: status
    end function UrbanDestroy_c

    function UrbanInitialize_c(handle) bind(C, name="UrbanInitialize") result(status)
      import :: c_ptr, c_int
      type(c_ptr), value :: handle
      integer(kind=c_int) :: status
    end function UrbanInitialize_c

    function UrbanSetInputs_c(handle, in) bind(C, name="UrbanSetInputs") result(status)
      import :: c_ptr, c_int, UrbanInputs
      type(c_ptr), value :: handle
      type(UrbanInputs), intent(in) :: in
      integer(kind=c_int) :: status
    end function UrbanSetInputs_c

    function UrbanStep_c(handle) bind(C, name="UrbanStep") result(status)
      import :: c_ptr, c_int
      type(c_ptr), value :: handle
      integer(kind=c_int) :: status
    end function UrbanStep_c

    function UrbanGetOutputs_c(handle, out) bind(C, name="UrbanGetOutputs") result(status)
      import :: c_ptr, c_int, UrbanOutputs
      type(c_ptr), value :: handle
      type(UrbanOutputs), intent(inout) :: out
      integer(kind=c_int) :: status
    end function UrbanGetOutputs_c
  end interface

contains

  ! Helper to build UrbanArrayD_c from a Fortran real array
  ! Note: Do not pass temporaries or array expressions here.
  ! The actual argument must be a persistent, contiguous array whose
  ! lifetime spans all C calls that use the returned pointer.
  function make_array_d(a) result(x)
    real(kind=c_double), intent(in), target :: a(:)
    type(UrbanArrayD_c) :: x
    x%data = c_loc(a(1))
    x%size = size(a, kind=c_size_t)
  end function make_array_d

  ! Setter functions for UrbanConfig
  subroutine UrbanConfigSetNLun(cfg, n_lun, status)
    type(UrbanConfig), intent(inout) :: cfg
    integer(kind=c_int), intent(in) :: n_lun
    integer(kind=c_int), intent(out) :: status
    cfg%N_LUN = n_lun
    status = URBAN_SUCCESS
  end subroutine UrbanConfigSetNLun

  subroutine UrbanConfigSetNRadBand(cfg, n_rad_band, status)
    type(UrbanConfig), intent(inout) :: cfg
    integer(kind=c_int), intent(in) :: n_rad_band
    integer(kind=c_int), intent(out) :: status
    cfg%N_RAD_BAND = n_rad_band
    status = URBAN_SUCCESS
  end subroutine UrbanConfigSetNRadBand

  ! Setter functions for UrbanInputs
  subroutine UrbanInputsSetSolarDown(in, a, status)
    type(UrbanInputs), intent(inout) :: in
    real(kind=c_double), intent(in), target :: a(:)
    integer(kind=c_int), intent(out) :: status
    in%solar_down = make_array_d(a)
    status = URBAN_SUCCESS
  end subroutine UrbanInputsSetSolarDown

  subroutine UrbanInputsSetLongwaveDown(in, a, status)
    type(UrbanInputs), intent(inout) :: in
    real(kind=c_double), intent(in), target :: a(:)
    integer(kind=c_int), intent(out) :: status
    in%longwave_down = make_array_d(a)
    status = URBAN_SUCCESS
  end subroutine UrbanInputsSetLongwaveDown

  subroutine UrbanInputsSetAirTemp(in, a, status)
    type(UrbanInputs), intent(inout) :: in
    real(kind=c_double), intent(in), target :: a(:)
    integer(kind=c_int), intent(out) :: status
    in%air_temp = make_array_d(a)
    status = URBAN_SUCCESS
  end subroutine UrbanInputsSetAirTemp

  subroutine UrbanInputsSetWindSpeed(in, a, status)
    type(UrbanInputs), intent(inout) :: in
    real(kind=c_double), intent(in), target :: a(:)
    integer(kind=c_int), intent(out) :: status
    in%wind_speed = make_array_d(a)
    status = URBAN_SUCCESS
  end subroutine UrbanInputsSetWindSpeed

  ! Getter functions for UrbanOutputs
  subroutine UrbanOutputsGetNetShortwave(out, a, status)
    type(UrbanOutputs), intent(inout) :: out
    real(kind=c_double), intent(in), target :: a(:)
    integer(kind=c_int), intent(out) :: status
    out%net_shortwave = make_array_d(a)
    status = URBAN_SUCCESS
  end subroutine UrbanOutputsGetNetShortwave

  subroutine UrbanOutputsGetNetLongwave(out, a, status)
    type(UrbanOutputs), intent(inout) :: out
    real(kind=c_double), intent(in), target :: a(:)
    integer(kind=c_int), intent(out) :: status
    out%net_longwave = make_array_d(a)
    status = URBAN_SUCCESS
  end subroutine UrbanOutputsGetNetLongwave

  subroutine UrbanOutputsGetSurfaceFlux(out, a, status)
    type(UrbanOutputs), intent(inout) :: out
    real(kind=c_double), intent(in), target :: a(:)
    integer(kind=c_int), intent(out) :: status
    out%surface_flux = make_array_d(a)
    status = URBAN_SUCCESS
  end subroutine UrbanOutputsGetSurfaceFlux

  ! Fortran-friendly wrappers that accept UrbanType and forward to C bindings
  subroutine UrbanCreate(cfg, sim, status)
    type(UrbanConfig), intent(in) :: cfg
    type(UrbanType), intent(inout) :: sim
    integer(kind=c_int), intent(out) :: status
    status = UrbanCreate_c(cfg, sim%c_ptr)
  end subroutine UrbanCreate

  ! Set a logging callback. `fn` must be a C-interoperable procedure with bind(C)
  ! signature: subroutine log(level, message, user_data) bind(C)
  !   integer(c_int), value :: level
  !   type(c_ptr), value :: message
  !   type(c_ptr), value :: user_data
  subroutine UrbanSetLogger(sim, fn, user_data, status)
    type(UrbanType), intent(inout) :: sim
    type(c_funptr), intent(in) :: fn
    type(c_ptr), intent(in) :: user_data
    integer(kind=c_int), intent(out) :: status
    status = UrbanSetLogger_c(sim%c_ptr, fn, user_data)
  end subroutine UrbanSetLogger

  subroutine UrbanDestroy(sim, status)
    type(UrbanType), intent(inout) :: sim
    integer(kind=c_int), intent(out) :: status
    status = UrbanDestroy_c(sim%c_ptr)
    if (status == 0_c_int) sim%c_ptr = c_null_ptr
  end subroutine UrbanDestroy

  subroutine UrbanInitialize(sim, status)
    type(UrbanType), intent(in) :: sim
    integer(kind=c_int), intent(out) :: status
    status = UrbanInitialize_c(sim%c_ptr)
  end subroutine UrbanInitialize

  subroutine UrbanSetInputs(sim, in, status)
    type(UrbanType), intent(in) :: sim
    type(UrbanInputs), intent(in) :: in
    integer(kind=c_int), intent(out) :: status
    status = UrbanSetInputs_c(sim%c_ptr, in)
  end subroutine UrbanSetInputs

  subroutine UrbanStep(sim, status)
    type(UrbanType), intent(in) :: sim
    integer(kind=c_int), intent(out) :: status
    status = UrbanStep_c(sim%c_ptr)
  end subroutine UrbanStep

  subroutine UrbanGetOutputs(sim, out, status)
    type(UrbanType), intent(in) :: sim
    type(UrbanOutputs), intent(inout) :: out
    integer(kind=c_int), intent(out) :: status
    status = UrbanGetOutputs_c(sim%c_ptr, out)
  end subroutine UrbanGetOutputs

end module urban
