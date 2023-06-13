program gtfv3_driver

  use iso_c_binding, only: c_ptr, c_null_ptr

  use grid_bounds_mod, only: GridBounds_T
  use domain_dim_mod, only: DomainDim_T
  use input_scalars_mod, only: InputScalars_T
  use input_arrays_mod, only: InputArrays_T
  use geos_gtfv3_interface_mod, only: geos_gtfv3_interface_f, geos_gtfv3_interface_finalize_f

  implicit none

  include 'mpif.h'

  integer, parameter :: NITER = 2, NTILES = 6, NUM_TRACERS = 7

  integer :: irank, nranks, mpierr
  character(len=256) :: input_file

  type(GridBounds_T) :: bd
  type(DomainDim_T) :: dim
  type(InputScalars_T) :: scalars
  type(InputArrays_T) :: arr

  integer :: ctr
  character(len=23) :: dt_iso
  real :: start, finish
  real, dimension(NITER) :: timings

  ! Start
  call MPI_Init(mpierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nranks, mpierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, irank, mpierr)

  ! Read input data - scalars
  write(input_file, '(a23, i1, a4)') 'input-data/scalar_data.', irank, '.bin'
  scalars = InputScalars_T(input_file)
  call scalars%wr1te()
  bd = scalars%get_grid_bounds()
  dim = scalars%get_domain_dimensions()

  ! Read input data - arrays
  write(input_file, '(a22, i1, a4)') 'input-data/array_data.', irank, '.bin'
  arr = InputArrays_T(input_file, bd, dim, scalars%nq_tot)
  ! call arr%wr1te()

  if (irank == 0) then
     call get_date_time_isoformat(dt_iso)
     write(*, '(a2,1x,a23,a18)') 'F:', dt_iso, ' --calling fortran interface'
 end if

 if (irank == 0) print*, 'irank, sum(u), sum(v), sum(w), sum(delz)'
 call MPI_Barrier(MPI_COMM_WORLD, mpierr)
 print *, irank, sum(arr%u), sum (arr%v), sum(arr%w), sum(arr%delz)
 call MPI_Barrier(MPI_COMM_WORLD, mpierr)

 ! Run gtFV3
 do ctr = 1, NITER
    call cpu_time(start)
    call geos_gtfv3_interface_f( &
         MPI_COMM_WORLD, &
         dim%npx, dim%npy, dim%npz, NTILES, &
         bd%is, bd%ie, bd%js, bd%je, &
         bd%isd, bd%ied, bd%jsd, bd%jed, &
         scalars%dt, NUM_TRACERS, scalars%ng, scalars%ptop, scalars%ks, &
         scalars%layout_1, scalars%layout_2, &
         scalars%adiabatic, &

         ! Input/Output
         arr%u, arr%v, arr%w, arr%delz, &
         arr%pt, arr%delp, arr%q(:,:,:,1:NUM_TRACERS), &
         arr%ps, arr%pe, arr%pk, arr%peln, arr%pkz, &
         arr%phis, arr%q_con, arr%omga, &
         arr%ua, arr%va, arr%uc, arr%vc, &

         ! Input
         arr%ak, arr%bk, &

         ! Input/Output
         arr%mfx, arr%mfy, arr%cx, arr%cy, arr%diss_est)
    call cpu_time(finish)
    timings(ctr) = finish-start
 end do

 call geos_gtfv3_interface_finalize_f()

 print *, ctr, '- timings on rank', irank, ': ', timings

 if (irank == 0) print*, 'rank, sum(u), sum(v), sum(w), sum(delz)'
 call MPI_Barrier(MPI_COMM_WORLD, mpierr)
 print *, irank, sum(arr%u), sum (arr%v), sum(arr%w), sum(arr%delz)
 call MPI_Barrier(MPI_COMM_WORLD, mpierr)

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

end program gtfv3_driver
