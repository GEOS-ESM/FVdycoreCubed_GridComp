#define ASSERT_(a) if(.not. a) then; print *, __FILE__, __LINE__; call MPI_Abort(MPI_COMM_WORLD, 1, mpierr); endif

program main

  use iso_c_binding, only: c_ptr, c_null_ptr

  use MAPL_ConstantsMod, only: MAPL_KAPPA, MAPL_CP, MAPL_RVAP, MAPL_RGAS

  use fms_mod, only: fms_init, fms_end
  use fv_control_mod, only: fv_init1, fv_init2, fv_end
  use fv_arrays_mod , only: fv_atmos_type
  use fv_state_utilities_mod, only: set_resolution_dependent_geos_defaults
  use fv_dynamics_mod, only: fv_dynamics

  use grid_bounds_mod, only: GridBounds_T
  use domain_dim_mod, only: DomainDim_T
  use input_scalars_mod, only: InputScalars_T
  use input_arrays_mod, only: InputArrays_T
  use geos_gtfv3_interface_mod, only: geos_gtfv3_interface_f

  implicit none

  include 'mpif.h'

  integer, parameter :: NITER = 1, NTILES = 6, NUM_TRACERS = 7
  logical, parameter :: fix_mass = .true.
  real, parameter :: time_total = -1

  integer :: irank, nranks, mpierr
  character(len=256) :: input_file

  type(GridBounds_T) :: bd
  type(DomainDim_T) :: dim
  type(InputScalars_T) :: scalars
  type(InputArrays_T) :: arr

  integer :: ctr, p_split
  character(len=23) :: dt_iso
  real :: start, finish, zvir
  type(fv_atmos_type), allocatable :: FV_Atm(:)
  logical, allocatable :: grids_on_this_pe(:)

  ! Start
  call MPI_Init(mpierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nranks, mpierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, irank, mpierr)

  ! Read input data - scalars
  write(input_file, '(a23, i0.2, a4)') 'input-data/scalar_data.', irank, '.bin'
  scalars = InputScalars_T(input_file)
  call scalars%wr1te()
  bd = scalars%get_grid_bounds()
  dim = scalars%get_domain_dimensions()
  ASSERT_(scalars%layout_1 * scalars%layout_2 * NTILES == nranks)

  ! Initialize FV3
  call fms_init(MPI_COMM_WORLD)
  p_split = 1
  call fv_init1(FV_Atm, scalars%dt, grids_on_this_pe, p_split)
  FV_Atm(1)%flagstruct%ntiles = NTILES
  FV_Atm(1)%flagstruct%npx = dim%npx + 1
  FV_Atm(1)%flagstruct%npy = dim%npy + 1
  call set_resolution_dependent_geos_defaults(fix_mass, scalars%dt, FV_Atm)
  call fv_init2(FV_Atm, scalars%dt, grids_on_this_pe, p_split)
  ASSERT_(FV_Atm(1)%flagstruct%hydrostatic .eqv. .false.)
  
  ! Read input data - arrays
  write(input_file, '(a22, i0.2, a4)') 'input-data/array_data.', irank, '.bin'
  arr = InputArrays_T(input_file, bd, dim, scalars%nq_tot)
  ! call arr%wr1te()

 !  if (irank == 0) then
 !     call get_date_time_isoformat(dt_iso)
 !     write(*, '(a2,1x,a23,a18)') 'F:', dt_iso, ' --calling fortran interface'
 ! end if

 if (irank == 0) print*, 'irank, sum(u), sum(v), sum(w), sum(delz)'
 call MPI_Barrier(MPI_COMM_WORLD, mpierr)
 print *, irank, sum(arr%u), sum (arr%v), sum(arr%w), sum(arr%delz)
 call MPI_Barrier(MPI_COMM_WORLD, mpierr)

 ! Run FV3
 zvir = MAPL_RVAP/MAPL_RGAS - 1
 do ctr = 1, NITER
    call cpu_time(start)
    call fv_dynamics( &
         ! Input
         FV_Atm(1)%npx, FV_Atm(1)%npy, FV_Atm(1)%npz, NUM_TRACERS, FV_Atm(1)%ng, scalars%dt, &
         FV_Atm(1)%flagstruct%consv_te, FV_Atm(1)%flagstruct%fill, FV_Atm(1)%flagstruct%reproduce_sum, &
         MAPL_KAPPA, MAPL_CP, zvir, &
         scalars%ptop, scalars%ks, NUM_TRACERS, &
         FV_Atm(1)%flagstruct%k_split, FV_Atm(1)%flagstruct%n_split, FV_Atm(1)%flagstruct%q_split, &
         ! Input/Output
         arr%u, arr%v, arr%w, arr%delz, &
         ! Input
         FV_Atm(1)%flagstruct%hydrostatic, &
         ! Input/Output
         arr%pt, arr%delp, arr%q(:, :, :, 1:NUM_TRACERS), &
         arr%ps, arr%pe, arr%pk, arr%peln, arr%pkz, &
         arr%phis, arr%varflt, arr%q_con, arr%omga, &
         arr%ua, arr%va, arr%uc, arr%vc, &
         ! Input
         arr%ak, arr%bk, &
         ! Input/Output
         arr%mfx, arr%mfy, arr%cx, arr%cy, &
         ! Input
         FV_Atm(1)%ze0, FV_Atm(1)%flagstruct%hybrid_z, &
         FV_Atm(1)%gridstruct, FV_Atm(1)%flagstruct, &
         FV_Atm(1)%neststruct, FV_Atm(1)%idiag, FV_Atm(1)%bd, &
         FV_Atm(1)%parent_grid, FV_Atm(1)%domain, &
         ! Input/Output
         arr%diss_est, arr%u_dt, arr%v_dt, arr%w_dt, arr%t_dt, &
         ! Input
         time_total)
    call cpu_time(finish)
    print *, irank, ', fv_dynamics: time taken = ', finish - start, 's'
 end do
 
 if (irank == 0) print*, 'irank, sum(u), sum(v), sum(w), sum(delz)'
 call MPI_Barrier(MPI_COMM_WORLD, mpierr)
 print *, irank, sum(arr%u), sum (arr%v), sum(arr%w), sum(arr%delz)
 call MPI_Barrier(MPI_COMM_WORLD, mpierr)

 call fv_end(FV_Atm, grids_on_this_pe, .false.)
 call MPI_Finalize(mpierr)

contains

  subroutine get_date_time_isoformat(dt_iso)
    ! Arguments
    character(len=23), intent(out) :: dt_iso
    ! Locals
    character(len=8) :: date
    character(len=10) :: time
    character(len=*), parameter :: fmt = '(a4,a1,a2,a1,a2,a1,a2,a1,a2,a1,a6)'

    call date_and_time(date=date, time=time)
    write(dt_iso, fmt) date(1:4), '-', date(5:6), '-', date(7:8), 'T', time(1:2), ':', time(3:4), ':', time(5:10)

  end subroutine get_date_time_isoformat

end program main
