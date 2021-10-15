module input_scalars_mod

  use grid_bounds_mod, only: GridBounds_T
  use domain_dim_mod, only: DomainDim_T

  implicit none
  
  private

  public InputScalars_T

  type InputScalars_T
     type(GridBounds_T), private :: bd
     type(DomainDim_T), private :: dim
     integer :: nq_tot, ng
     real :: dt, consv_te, kappa, cp_air, zvir, ptop
     integer :: ks, ncnst, n_split, q_split
     logical :: fill, reproduce_sum, hydrostatic, hybrid_z
   contains
     procedure :: wr1te
     procedure :: get_grid_bounds
     procedure :: get_domain_dimensions
  end type InputScalars_T

  interface InputScalars_T
     procedure :: read_data_from_file
  end interface InputScalars_T

contains

  function read_data_from_file(file_name) result(scalars)

    ! Arguments
    character(len=*), intent(in) :: file_name
    type(InputScalars_T) :: scalars ! output

    ! Locals
    integer :: file_handle

    ! Start
    open(newunit = file_handle, file = file_name, form = 'unformatted', status = 'old')
    read(file_handle) scalars%dim%npx, scalars%dim%npy, scalars%dim%npz
    read(file_handle) scalars%bd%isd, scalars%bd%ied, scalars%bd%jsd, scalars%bd%jed
    read(file_handle) &
         scalars%nq_tot, scalars%ng, scalars%dt, &
         scalars%consv_te, scalars%fill, scalars%reproduce_sum, &
         scalars%kappa, scalars%cp_air, scalars%zvir, scalars%ptop, scalars%ks, &
         scalars%ncnst, scalars%n_split, scalars%q_split, &
         scalars%hydrostatic, scalars%hybrid_z
    close(file_handle)
    
  end function read_data_from_file

  subroutine wr1te(self)

    ! Arguments
    class(InputScalars_T), intent(in) :: self

    call self%bd%wr1te()
    call self%dim%wr1te()
    print *, 'nq_tot: ', self%nq_tot
    print *, 'ng: ', self%ng
    print *, 'dt: ', self%dt
    print *, 'consv_te: ', self%consv_te
    print *, 'fill: ', self%fill
    print *, 'reproduce_sum: ', self%reproduce_sum
    print *, 'kappa: ', self%kappa
    print *, 'cp_air: ', self%cp_air
    print *, 'zvir: ', self%zvir
    print *, 'ptop: ', self%ptop
    print *, 'ks: ', self%ks
    print *, 'ncnst: ', self%ncnst
    print *, 'n_split: ', self%n_split
    print *, 'q_split: ', self%q_split
    print *, 'hydrostatic: ', self%hydrostatic
    print *, 'hybrid_z: ', self%hybrid_z
    
  end subroutine wr1te

  function get_grid_bounds(self) result(bd)

    ! Arguments
    class(InputScalars_T), intent(in) :: self
    type(GridBounds_T) :: bd ! output

    bd = self%bd

  end function get_grid_bounds

  function get_domain_dimensions(self) result(dim)

    ! Arguments
    class(InputScalars_T), intent(in) :: self
    type(DomainDim_T) :: dim ! output

    dim = self%dim

  end function get_domain_dimensions
  
end module input_scalars_mod
