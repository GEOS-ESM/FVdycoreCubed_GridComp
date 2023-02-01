#define ASSERT_(a) if(.not. a) then; print *, __FILE__, __LINE__; call MPI_Abort(MPI_COMM_WORLD, 1, mpierr); endif

program fv3_driver

  use fms_mod, only: fms_init, fms_end
  use fv_control_mod, only: fv_init1, fv_init2, fv_end
  use fv_arrays_mod , only: fv_atmos_type
  use fv_state_utilities_mod, only: set_resolution_dependent_geos_defaults

  implicit none

  include "mpif.h"

  integer, parameter :: nx = 1, ny = 6
  logical, parameter :: fix_mass = .true.
  real, parameter :: dt = 450.0

  type(fv_atmos_type), allocatable :: FV_Atm(:)
  logical, allocatable :: grids_on_this_pe(:)
  integer :: mpierr, irank, nranks, p_split

  call MPI_Init(mpierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, irank, mpierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nranks, mpierr)
  ASSERT_(nx*ny == nranks)

  call fms_init(MPI_COMM_WORLD)
  p_split = 1
  call fv_init1(FV_Atm, DT, grids_on_this_pe, p_split)
  FV_Atm(1)%flagstruct%ntiles = 6
  FV_Atm(1)%flagstruct%npx = nx + 1
  FV_Atm(1)%flagstruct%npy = ny/6 + 1
  call set_resolution_dependent_geos_defaults(fix_mass, dt, FV_Atm)
  call fv_init2(FV_Atm, DT, grids_on_this_pe, p_split)

  call fv_end(FV_Atm, grids_on_this_pe, .false.)
  ! call fms_end
  call MPI_Finalize(mpierr)

end program fv3_driver
