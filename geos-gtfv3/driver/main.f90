program main

  use iso_c_binding, only: c_ptr, c_null_ptr

  use grid_bounds_mod, only: GridBounds_T
  use domain_dim_mod, only: DomainDim_T
  use input_scalars_mod, only: InputScalars_T
  use input_arrays_mod, only: InputArrays_T
  use geos_gtfv3_interface_mod, only: geos_gtfv3_interface_f

  implicit none

  include 'mpif.h'

  integer :: rank, size, mpierr, bin_file_handle
  character(len=256) :: input_file, input_arrays_file

  integer :: npx, npy, npz, nq_tot, ng
  type(GridBounds_T) :: bd
  type(DomainDim_T) :: dim
  type(InputScalars_T) :: scalars
  type(InputArrays_T) :: arr
  type(c_ptr) :: grid_data = c_null_ptr

  integer :: dt(8), ctr
  character(len=23) :: date_time_s
  character(len=*), parameter :: iso8601 = '(i4, 2("-", i2.2), "T", 2(i0.2, ":"), i0.2, ".", i0.3)'
  real :: start, finish

  ! Start
  call MPI_Init(mpierr)
  call MPI_Comm_size(MPI_COMM_WORLD, size, mpierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, mpierr)

  ! Read input data
  write(input_file, '(a23, i1, a4)') 'input-data/scalar_data.', rank, '.bin'
  scalars = InputScalars_T(input_file)
  ! call scalars%wr1te()
  bd = scalars%get_grid_bounds()
  dim = scalars%get_domain_dimensions()
  write(input_file, '(a22, i1, a4)') 'input-data/array_data.', rank, '.bin'
  arr = InputArrays_T(input_file, bd, dim, scalars%ncnst)
  ! call arr%wr1te()

  if (rank == 0) then
     call date_and_time(values=dt)
     write(date_time_s, iso8601) dt(1:3), dt(5:8)
     write(*, '(a2,1x,a23,1x,a27)') 'F:', date_time_s, '--calling fortran interface'
  end if

  scalars%kord_tm = -9
  scalars%tau = 0.0

  ! if (rank == 0) print*, 'rank, sum(u), sum(v), sum(w), sum(delz)'
  ! call MPI_Barrier(MPI_COMM_WORLD, mpierr)
  ! print *, rank, sum(arr%u), sum (arr%v), sum(arr%w), sum(arr%delz)
  ! call MPI_Barrier(MPI_COMM_WORLD, mpierr)

  do ctr = 1, 10
     call cpu_time(start)
     call geos_gtfv3_interface_f( &
          MPI_COMM_WORLD, &
          dim%npx, dim%npy, dim%npz, &
          bd%is, bd%ie, bd%js, bd%je, &
          bd%isd, bd%ied, bd%jsd, bd%jed, &
          scalars%dt, scalars%nq_tot, scalars%ng, scalars%ptop, scalars%ks, &
          scalars%layout_1, scalars%layout_2, &
          scalars%adiabatic, &
          scalars%hydrostatic, scalars%z_tracer, scalars%make_nh, scalars%fv_debug, &
          scalars%reproduce_sum, scalars%do_sat_adj, &
          scalars%do_vort_damp, scalars%rf_fast, scalars%fill, &
          scalars%ncnst, scalars%n_split, scalars%k_split, &
          scalars%fv_sg_adj, scalars%n_sponge, scalars%n_zfilter, scalars%nwat, &
          scalars%hord_tr, scalars%hord_tm, scalars%hord_dp, scalars%hord_mt, scalars%hord_vt, &
          scalars%nord, scalars%kord_tm, scalars%kord_tr, scalars%kord_wz, scalars%kord_mt, &
          scalars%d_ext, scalars%beta, scalars%vtdm4, scalars%ke_bg, &
          scalars%d_con, scalars%d2_bg, scalars%d2_bg_k1, scalars%d2_bg_k2, &
          scalars%p_fac, scalars%a_imp, scalars%dddmp, scalars%d4_bg, &
          scalars%rf_cutoff, scalars%tau, scalars%consv_te, &
          
          arr%u, arr%v, arr%w, arr%delz, &
          arr%pt, arr%delp, arr%q, &
          arr%ps, arr%pe, arr%pk, arr%peln, arr%pkz, &
          arr%phis, arr%q_con, arr%omga, &
          arr%ua, arr%va, arr%uc, arr%vc, &
          
          arr%ak, arr%bk, &
          
          arr%mfx, arr%mfy, arr%cx, arr%cy, arr%diss_est, &
          
          arr%dx, arr%dy, arr%dxa, arr%dya, arr%dxc, arr%dyc, &
          arr%rdx, arr%rdy, arr%rdxa, arr%rdya, arr%rdxc, arr%rdyc, &
          
          arr%cosa, arr%cosa_s, arr%sina_u, arr%sina_v, arr%cosa_u, arr%cosa_v, &
          arr%rsin2, arr%rsina, arr%rsin_u, arr%rsin_v, &
          arr%sin_sg, arr%cos_sg, &
          arr%area, arr%rarea, arr%rarea_c, arr%f0, arr%fC, &
          arr%del6_u, arr%del6_v, arr%divg_u, arr%divg_v, &
          arr%agrid, arr%bgrid, arr%a11, arr%a12, arr%a21, arr%a22, &
          arr%edge_e, arr%edge_w, arr%edge_n, arr%edge_s, &
          arr%nested, arr%stretched_grid, arr%da_min, arr%da_min_c)
     call cpu_time(finish)
     print *, ctr, '- time taken on rank ', rank, ': ', finish-start
  end do

  if (rank == 0) print*, 'rank, sum(u), sum(v), sum(w), sum(delz)'
  call MPI_Barrier(MPI_COMM_WORLD, mpierr)
  print *, rank, sum(arr%u), sum (arr%v), sum(arr%w), sum(arr%delz)
  call MPI_Barrier(MPI_COMM_WORLD, mpierr)

  call MPI_Finalize(mpierr)

end program main
