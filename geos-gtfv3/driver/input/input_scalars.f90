module input_scalars_mod

  use grid_bounds_mod, only: GridBounds_T
  use domain_dim_mod, only: DomainDim_T

  implicit none

  private

  public InputScalars_T

  type InputScalars_T
     type(GridBounds_T), private :: bd
     type(DomainDim_T), private :: dim
     real :: dt
     integer :: nq_tot, ng
     real :: ptop
     integer :: ks, layout_1, layout_2
     logical :: adiabatic
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
    integer :: layout(2)

    ! Start
    open(newunit = file_handle, file = file_name, form = 'unformatted', status = 'old')
    read(file_handle) scalars%dim%npx, scalars%dim%npy, scalars%dim%npz
    read(file_handle) scalars%bd%is, scalars%bd%ie, scalars%bd%js, scalars%bd%je
    read(file_handle) scalars%bd%isd, scalars%bd%ied, scalars%bd%jsd, scalars%bd%jed
    read(file_handle) &
         scalars%dt, scalars%nq_tot, scalars%ng, scalars%ptop, scalars%ks, &
         layout, scalars%adiabatic
    close(file_handle)
    scalars%layout_1 = layout(1)
    scalars%layout_2 = layout(2)

  end function read_data_from_file

  subroutine wr1te(self)

    ! Arguments
    class(InputScalars_T), intent(in) :: self

    call self%bd%wr1te()
    call self%dim%wr1te()
    print *, 'dt: ', self%dt
    print *, 'nq_tot/ng: ', self%nq_tot, self%ng
    print *, 'ptop/ks: ', self%ptop, self%ks
    print *, 'layout: ', '[', self%layout_1, ',', self%layout_2, ']'
    print *, 'adiabatic: ', self%adiabatic
    
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
