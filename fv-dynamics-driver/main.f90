program main

  use grid_bounds_mod, only: GridBounds_T
  use domain_dim_mod, only: DomainDim_T
  use input_scalars_mod, only: InputScalars_T
  use input_arrays_mod, only: InputArrays_T
  use fv_dynamics_interface_mod, only: fv_dynamics_interface
  
  implicit none

  include 'mpif.h'

  integer :: rank, size, mpierr, bin_file_handle
  character(len=256) :: input_file, input_arrays_file
  
  integer :: npx, npy, npz, nq_tot, ng
  type(GridBounds_T) :: bd
  type(DomainDim_T) :: dim
  type(InputScalars_T) :: scalars
  type(InputArrays_T) :: arr
  
  ! Start
  call MPI_Init(mpierr)
  call MPI_Comm_size(MPI_COMM_WORLD, size, mpierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, mpierr)

  ! Read input data
  write(input_file, '(a23, i1, a4)') 'input-data/scalar_data.', rank, '.bin'
  scalars = InputScalars_T(input_file)
  call scalars%wr1te()
  bd = scalars%get_grid_bounds()
  dim = scalars%get_domain_dimensions()
  write(input_file, '(a22, i1, a4)') 'input-data/array_data.', rank, '.bin'
  arr = InputArrays_T(input_file, bd, dim)
  
  call fv_dynamics_interface( &
       MPI_COMM_WORLD, &
       dim%npx, dim%npy, dim%npz, scalars%nq_tot, scalars%ng, &
       bd%isd, bd%ied, bd%jsd, bd%jed, scalars%dt, &
       scalars%consv_te, scalars%fill, scalars%reproduce_sum, scalars%kappa, &
       scalars%cp_air, scalars%zvir, scalars%ptop, scalars%ks, &
       scalars%ncnst, scalars%n_split, scalars%q_split, &
       arr%u, arr%v, arr%w)
  
  call MPI_Finalize(mpierr)
  
end program main
