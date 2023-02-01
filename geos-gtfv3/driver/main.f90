  program main

    use iso_c_binding, only: c_ptr, c_null_ptr

    use grid_bounds_mod, only: GridBounds_T
    use domain_dim_mod, only: DomainDim_T
    use input_scalars_mod, only: InputScalars_T
    use input_arrays_mod, only: InputArrays_T
    use geos_gtfv3_interface_mod, only: geos_gtfv3_interface_f

    use MAPL_ConstantsMod, only: MAPL_CP, MAPL_RGAS, MAPL_RVAP, MAPL_GRAV, MAPL_RADIUS, &
                                 MAPL_KAPPA, MAPL_PI_R8, MAPL_ALHL, MAPL_PSDRY
    use fv_control_mod, only: fv_init1, fv_init2
    use fv_arrays_mod , only: fv_atmos_type, FVPRC
    use fv_dynamics_mod, only: fv_dynamics
    use fms_mod, only: fms_init
    use fv_state_utilities_mod, only: set_resolution_dependent_geos_defaults
    use ieee_exceptions, only: ieee_get_halting_mode, ieee_set_halting_mode, ieee_all
    
    implicit none

    include 'mpif.h'

    integer, parameter :: NITER = 2, NTILES = 6

    integer :: rank, size, mpierr
    character(len=256) :: input_file

    type(GridBounds_T) :: bd
    type(DomainDim_T) :: dim
    type(InputScalars_T) :: scalars
    type(InputArrays_T) :: gt_arr
    type(InputArrays_T) :: fn_arr

    real(FVPRC) :: kappa    ! kappa
    real(FVPRC) :: cp       ! heat capacity of air at constant pressure
    real(FVPRC) :: zvir
    logical :: fix_mass = .true.
    integer :: p_split=1
    type(fv_atmos_type), allocatable, save :: FV_Atm(:)
    logical, allocatable, save             :: grids_on_this_pe(:)

    logical :: halting_mode(5)
    integer :: ctr
    character(len=23) :: dt_iso
    real :: start, finish
    logical :: run_fn = .false.
    logical :: run_gt = .false
    
    ! Start
    call MPI_Init(mpierr)
    call MPI_Comm_size(MPI_COMM_WORLD, size, mpierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, mpierr)
    run_gt = size.eq.1.or.size.eq.6
    run_fn = size.eq.6.or.size.eq.72
    
    ! Read input data
    write(input_file, '(a23, i1, a4)') 'input-data/scalar_data.', rank, '.bin'
    scalars = InputScalars_T(input_file)
    call scalars%wr1te()
    bd = scalars%get_grid_bounds()
    dim = scalars%get_domain_dimensions()
    write(input_file, '(a22, i1, a4)') 'input-data/array_data.', rank, '.bin'
    gt_arr = InputArrays_T(input_file, bd, dim, scalars%nq_tot)
    fn_arr = InputArrays_T(input_file, bd, dim, scalars%nq_tot)
    
    ! Init Fn FV3
    if (run_fn) then
      call fms_init(MPI_COMM_WORLD)
      call fv_init1(FV_Atm, scalars%dt, grids_on_this_pe, p_split) !ksplit
      FV_Atm(1)%flagstruct%npx    = dim%npx
      FV_Atm(1)%flagstruct%npy    = dim%npy
      FV_Atm(1)%flagstruct%npz    = dim%npz
      call set_resolution_dependent_geos_defaults(fix_mass, scalars%dt, FV_Atm)
      call fv_init2(FV_Atm, scalars%dt, grids_on_this_pe, p_split)
    endif
      
    kappa  = MAPL_KAPPA
    cp     = MAPL_CP
    zvir   = MAPL_RVAP/MAPL_RGAS - 1.   ! RWV/RAIR-1
    ! call arr%wr1te()

    if (rank == 0) then
      call get_date_time_isoformat(dt_iso)
      write(*, '(a2,1x,a23,a18)') 'F:', dt_iso, ' --calling fortran interface'
    end if

    ! if (rank == 0) print*, 'rank, sum(u), sum(v), sum(w), sum(delz)'
    ! call MPI_Barrier(MPI_COMM_WORLD, mpierr)
    ! print *, rank, sum(arr%u), sum (arr%v), sum(arr%w), sum(arr%delz)
    ! call MPI_Barrier(MPI_COMM_WORLD, mpierr)

    do ctr = 1, NITER
        ! Run GT
        if (run_gt) then
 
          ! A workaround to the issue of SIGFPE abort during importing of numpy, is to
          ! disable trapping of floating point exceptions temporarily, call the interface
          ! to the Python function and resume trapping
          call ieee_get_halting_mode(ieee_all, halting_mode)
          call ieee_set_halting_mode(ieee_all, .false.)

          call cpu_time(start)
          call geos_gtfv3_interface_f( &
              MPI_COMM_WORLD, &
              dim%npx, dim%npy, dim%npz, NTILES, &
              bd%is, bd%ie, bd%js, bd%je, &
              bd%isd, bd%ied, bd%jsd, bd%jed, &
              scalars%dt, scalars%nq_tot, scalars%ng, scalars%ptop, scalars%ks, &
              scalars%layout_1, scalars%layout_2, &
              scalars%adiabatic, &

              ! Input/Output
              gt_arr%u, gt_arr%v, gt_arr%w, gt_arr%delz, &
              gt_arr%pt, gt_arr%delp, gt_arr%q, &
              gt_arr%ps, gt_arr%pe, gt_arr%pk, gt_arr%peln, gt_arr%pkz, &
              gt_arr%phis, gt_arr%q_con, gt_arr%omga, &
              gt_arr%ua, gt_arr%va, gt_arr%uc, gt_arr%vc, &

              ! Input
              gt_arr%ak, gt_arr%bk, &

              ! Input/Output
              gt_arr%mfx, gt_arr%mfy, gt_arr%cx, gt_arr%cy, gt_arr%diss_est)
          call cpu_time(finish)
          print *, ctr, '- GT time taken on rank ', rank, ': ', finish-start

          call ieee_set_halting_mode(ieee_all, halting_mode)

        endif
    
        if (run_fn) then
          call cpu_time(start)
          call fv_dynamics( &
              dim%npx, dim%npy, dim%npz, FV_Atm(1)%ncnst, FV_Atm(1)%ng,   &
              scalars%dt,                                                 &
              FV_Atm(1)%flagstruct%consv_te, FV_Atm(1)%flagstruct%fill,   &
              FV_Atm(1)%flagstruct%reproduce_sum,                         &
              kappa, cp, zvir, scalars%ptop, scalars%ks, FV_Atm(1)%ncnst, &
              FV_Atm(1)%flagstruct%n_split, FV_Atm(1)%flagstruct%q_split, & !n_split, q_split
              fn_arr%u, fn_arr%v, fn_arr%w, fn_arr%delz,           &
              FV_Atm(1)%flagstruct%hydrostatic,                    &
              fn_arr%pt, fn_arr%delp, fn_arr%q, fn_arr%ps,         &
              fn_arr%pe, fn_arr%pk, fn_arr%peln, fn_arr%pkz,       &
              fn_arr%phis, fn_arr%q_con, fn_arr%omga,              &
              fn_arr%ua, fn_arr%va, fn_arr%uc, fn_arr%vc,          &
              fn_arr%ak, fn_arr%bk,                                &
              fn_arr%mfx, fn_arr%mfy, fn_arr%cx, fn_arr%cy,        &
              FV_Atm(1)%ze0,                                       &
              FV_Atm(1)%flagstruct%hybrid_z,                       &
              FV_Atm(1)%gridstruct, FV_Atm(1)%flagstruct,          &
              FV_Atm(1)%neststruct, FV_Atm(1)%idiag, FV_Atm(1)%bd, FV_Atm(1)%parent_grid, FV_Atm(1)%domain, &
              fn_arr%diss_est)
          call cpu_time(finish)
          print *, ctr, '- FN time taken on rank ', rank, ': ', finish-start
        endif

    end do

    ! if (rank == 0) print*, 'rank, sum(u), sum(v), sum(w), sum(delz)'
    ! call MPI_Barrier(MPI_COMM_WORLD, mpierr)
    ! print *, rank, sum(arr%u), sum (arr%v), sum(arr%w), sum(arr%delz)
    ! call MPI_Barrier(MPI_COMM_WORLD, mpierr)

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
