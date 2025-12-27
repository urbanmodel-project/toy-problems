module urban
  use, intrinsic :: iso_c_binding
  implicit none

  ! Fortran-friendly opaque wrapper around C UrbanType
  type :: UrbanType
    type(c_ptr) :: c_ptr
  end type UrbanType

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

    function UrbanSetLogger_c(handle, fn, user_data) bind(C, name="UrbanSetLogger") result(status)
      import :: c_int, c_ptr, c_funptr
      type(c_ptr),  value :: handle
      type(c_funptr), value :: fn
      type(c_ptr),  value :: user_data
      integer(c_int) :: status
    end function UrbanSetLogger_c

    function UrbanCreate_c(cfg, out) bind(C, name="UrbanCreate") result(status)
      import :: UrbanConfig_c, c_int, c_ptr
      type(UrbanConfig_c), intent(in) :: cfg
      type(c_ptr), intent(out) :: out
      integer(c_int) :: status
    end function UrbanCreate_c

    function UrbanDestroy_c(handle) bind(C, name="UrbanDestroy") result(status)
      import :: c_ptr, c_int
      ! Pointer-to-pointer: C signature is UrbanDestroy(UrbanType* handle)
      ! Pass-by-reference so C can null out the pointer.
      type(c_ptr), intent(inout) :: handle
      integer(c_int) :: status
    end function UrbanDestroy_c

    function UrbanInitialize_c(handle) bind(C, name="UrbanInitialize") result(status)
      import :: c_ptr, c_int
      type(c_ptr), value :: handle
      integer(c_int) :: status
    end function UrbanInitialize_c

    function UrbanSetInputs_c(handle, in) bind(C, name="UrbanSetInputs") result(status)
      import :: c_ptr, c_int, UrbanInputs_c
      type(c_ptr), value :: handle
      type(UrbanInputs_c), intent(in) :: in
      integer(c_int) :: status
    end function UrbanSetInputs_c

    function UrbanStep_c(handle) bind(C, name="UrbanStep") result(status)
      import :: c_ptr, c_int
      type(c_ptr), value :: handle
      integer(c_int) :: status
    end function UrbanStep_c

    function UrbanGetOutputs_c(handle, out) bind(C, name="UrbanGetOutputs") result(status)
      import :: c_ptr, c_int, UrbanOutputs_c
      type(c_ptr), value :: handle
      type(UrbanOutputs_c), intent(inout) :: out
      integer(c_int) :: status
    end function UrbanGetOutputs_c
  end interface

contains

  ! Helper to build UrbanArrayD_c from a Fortran real array
  ! Note: Do not pass temporaries or array expressions here.
  ! The actual argument must be a persistent, contiguous array whose
  ! lifetime spans all C calls that use the returned pointer.
  function make_array_d(a) result(x)
    real(c_double), intent(in), target :: a(:)
    type(UrbanArrayD_c) :: x
    x%data = c_loc(a(1))
    x%size = size(a, kind=c_size_t)
  end function make_array_d

  ! Fortran-friendly wrappers that accept UrbanType and forward to C bindings
  subroutine UrbanCreate(cfg, sim, status)
    type(UrbanConfig_c), intent(in) :: cfg
    type(UrbanType), intent(inout) :: sim
    integer(c_int), intent(out) :: status
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
    type(c_ptr),    intent(in) :: user_data
    integer(c_int), intent(out) :: status
    status = UrbanSetLogger_c(sim%c_ptr, fn, user_data)
  end subroutine UrbanSetLogger

  subroutine UrbanDestroy(sim, status)
    type(UrbanType), intent(inout) :: sim
    integer(c_int), intent(out) :: status
    status = UrbanDestroy_c(sim%c_ptr)
    if (status == 0_c_int) sim%c_ptr = c_null_ptr
  end subroutine UrbanDestroy

  subroutine UrbanInitialize(sim, status)
    type(UrbanType), intent(in) :: sim
    integer(c_int), intent(out) :: status
    status = UrbanInitialize_c(sim%c_ptr)
  end subroutine UrbanInitialize

  subroutine UrbanSetInputs(sim, in, status)
    type(UrbanType), intent(in) :: sim
    type(UrbanInputs_c), intent(in) :: in
    integer(c_int), intent(out) :: status
    status = UrbanSetInputs_c(sim%c_ptr, in)
  end subroutine UrbanSetInputs

  subroutine UrbanStep(sim, status)
    type(UrbanType), intent(in) :: sim
    integer(c_int), intent(out) :: status
    status = UrbanStep_c(sim%c_ptr)
  end subroutine UrbanStep

  subroutine UrbanGetOutputs(sim, out, status)
    type(UrbanType), intent(in) :: sim
    type(UrbanOutputs_c), intent(inout) :: out
    integer(c_int), intent(out) :: status
    status = UrbanGetOutputs_c(sim%c_ptr, out)
  end subroutine UrbanGetOutputs

end module urban
