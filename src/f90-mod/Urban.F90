module urban_f
  use, intrinsic :: iso_c_binding
  implicit none

  ! Opaque handle type for C pointer
  type, bind(C) :: UrbanArrayD_c
    type(c_ptr)        :: data
    integer(c_size_t)  :: size
  end type UrbanArrayD_c

  type, bind(C) :: UrbanConfig_c
    integer(c_int) :: N_LUN
    integer(c_int) :: N_RAD_BAND
    logical(c_bool) :: enable_openmp
    integer(c_int) :: omp_num_threads
  end type UrbanConfig_c

  type, bind(C) :: UrbanInputs_c
    type(UrbanArrayD_c) :: solar_down
    type(UrbanArrayD_c) :: longwave_down
    type(UrbanArrayD_c) :: air_temp
    type(UrbanArrayD_c) :: wind_speed
  end type UrbanInputs_c

  type, bind(C) :: UrbanOutputs_c
    type(UrbanArrayD_c) :: net_shortwave
    type(UrbanArrayD_c) :: net_longwave
    type(UrbanArrayD_c) :: surface_flux
  end type UrbanOutputs_c

  interface
    ! Kokkos C wrapper interfaces
    subroutine urban_kokkos_initialize() bind(C, name="urban_kokkos_initialize")
    end subroutine urban_kokkos_initialize

    subroutine urban_kokkos_finalize() bind(C, name="urban_kokkos_finalize")
    end subroutine urban_kokkos_finalize

    function UrbanGetErrorString(err) bind(C, name="UrbanGetErrorString") result(msg)
      import :: c_int, c_ptr
      integer(c_int), value :: err
      type(c_ptr) :: msg
    end function UrbanGetErrorString

    function UrbanCreate(cfg, out) bind(C, name="UrbanCreate") result(status)
      import :: UrbanConfig_c, c_int, c_ptr
      type(UrbanConfig_c), intent(in) :: cfg
      type(c_ptr) :: out
      integer(c_int) :: status
    end function UrbanCreate

    function UrbanDestroy(handle) bind(C, name="UrbanDestroy") result(status)
      import :: c_ptr, c_int
      type(c_ptr) :: handle
      integer(c_int) :: status
    end function UrbanDestroy

    function UrbanInitialize(handle) bind(C, name="UrbanInitialize") result(status)
      import :: c_ptr, c_int
      type(c_ptr), value :: handle
      integer(c_int) :: status
    end function UrbanInitialize

    function UrbanSetInputs(handle, in) bind(C, name="UrbanSetInputs") result(status)
      import :: c_ptr, c_int, UrbanInputs_c
      type(c_ptr), value :: handle
      type(UrbanInputs_c), intent(in) :: in
      integer(c_int) :: status
    end function UrbanSetInputs

    function UrbanStep(handle) bind(C, name="UrbanStep") result(status)
      import :: c_ptr, c_int
      type(c_ptr), value :: handle
      integer(c_int) :: status
    end function UrbanStep

    function UrbanGetOutputs(handle, out) bind(C, name="UrbanGetOutputs") result(status)
      import :: c_ptr, c_int, UrbanOutputs_c
      type(c_ptr), value :: handle
      type(UrbanOutputs_c), intent(inout) :: out
      integer(c_int) :: status
    end function UrbanGetOutputs
  end interface

contains

  ! Helper to build UrbanArrayD_c from a Fortran real array
  pure function make_array_d(a) result(x)
    real(c_double), intent(in), target :: a(:)
    type(UrbanArrayD_c) :: x
    x%data = c_loc(a(1))
    x%size = size(a, kind=c_size_t)
  end function make_array_d

end module urban_f
