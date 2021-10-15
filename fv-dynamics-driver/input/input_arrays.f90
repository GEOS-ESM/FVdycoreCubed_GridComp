module input_arrays_mod

  use grid_bounds_mod, only: GridBounds_T
  use domain_dim_mod, only: DomainDim_T
  
  implicit none

  public InputArrays_T

  type InputArrays_T
     real, allocatable, dimension(:,:,:) :: u, v, w
  end type InputArrays_T

  interface InputArrays_T
     procedure :: initialize_and_read_data_
  end interface InputArrays_T

contains

  function initialize_and_read_data_(file_name, bd, dim) result(arrays)

    ! Arguments
    character(len=*), intent(in) :: file_name
    type(GridBounds_T), intent(in) :: bd
    type(DomainDim_T), intent(in) :: dim
    type(InputArrays_T) :: arrays ! output

    ! Locals
    integer :: file_handle
    real, parameter :: undef = 1.0e15
    
    ! Start

    ! Allocate memory
    allocate(arrays%u(bd%isd:bd%ied, bd%jsd:bd%jed+1, dim%npz), source = undef)
    allocate(arrays%v(bd%isd:bd%ied+1, bd%jsd:bd%jed, dim%npz), source = undef)
    allocate(arrays%w(bd%isd:bd%ied,bd%jsd:bd%jed,1:dim%npz), source = undef)

    ! Now read data
    open(newunit = file_handle, file = file_name, form = 'unformatted', status = 'old')
    read(file_handle) arrays%u, arrays%v, arrays%w
    
  end function initialize_and_read_data_

end module input_arrays_mod
