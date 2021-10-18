module input_arrays_mod

  use grid_bounds_mod, only: GridBounds_T
  use domain_dim_mod, only: DomainDim_T

  implicit none

  public InputArrays_T

  type InputArrays_T
     real, allocatable, dimension(:,:,:) :: u, v, w, delz
     real, allocatable :: pt(:,:,:), delp(:,:,:), q(:,:,:,:)
     real, allocatable :: ps(:,:), pe(:,:,:), pk(:,:,:), peln(:,:,:), pkz(:,:,:)
     real, allocatable :: phis(:,:), q_con(:,:,:), omga(:,:,:)
     real, allocatable :: ua(:,:,:), va(:,:,:), uc(:,:,:), vc(:,:,:)
     real, allocatable :: ak(:), bk(:)
     real, allocatable :: mfx(:,:,:), mfy(:,:,:), cx(:,:,:), cy(:,:,:)
  end type InputArrays_T

  interface InputArrays_T
     procedure :: initialize_and_read_data_
  end interface InputArrays_T

contains

  function initialize_and_read_data_(file_name, bd, dim, ncnst) result(arr)

    ! Arguments
    character(len=*), intent(in) :: file_name
    type(GridBounds_T), intent(in) :: bd
    type(DomainDim_T), intent(in) :: dim
    integer, intent(in) :: ncnst
    type(InputArrays_T) :: arr ! output

    ! Locals
    integer :: file_handle
    real, parameter :: undef = 1.0e15

    ! Start

    ! Allocate memory
    !# u, v, w, delz
    allocate(arr%u(bd%isd:bd%ied, bd%jsd:bd%jed+1, 1:dim%npz), source = undef)
    allocate(arr%v(bd%isd:bd%ied+1, bd%jsd:bd%jed, 1:dim%npz), source = undef)
    allocate(arr%w(1, 1, 1), source = undef)
    allocate(arr%delz(1, 1, 1), source = undef)
    !# pt, delp, q
    allocate(arr%pt(bd%isd:bd%ied, bd%jsd:bd%jed, dim%npz), source = undef)
    allocate(arr%delp(bd%isd:bd%ied, bd%jsd:bd%jed, dim%npz), source = undef)
    allocate(arr%q(bd%isd:bd%ied, bd%jsd:bd%jed, dim%npz, ncnst), source = undef)
    !# ps, pe, pk, peln, pkz
    allocate(arr%ps(bd%isd:bd%ied, bd%jsd:bd%jed), source = undef)
    allocate(arr%pe(bd%is-1:bd%ie+1, dim%npz+1, bd%js-1:bd%je+1), source = undef)
    allocate(arr%pk(bd%is:bd%ie, bd%js:bd%je, dim%npz+1), source = undef)
    allocate(arr%peln(bd%is:bd%ie, dim%npz+1, bd%js:bd%je), source = undef)
    allocate(arr%pkz(bd%is:bd%ie, bd%js:bd%je, dim%npz), source = undef)
    !# phis, q_con, omga
    allocate(arr%phis(bd%isd:bd%ied, bd%jsd:bd%jed), source = undef)
    allocate(arr%q_con(1, 1, 1), source = undef)
    allocate(arr%omga(bd%isd:bd%ied, bd%jsd:bd%jed, dim%npz), source = undef)
    !# ua, va, uc, vc
    allocate(arr%ua(bd%isd:bd%ied, bd%jsd:bd%jed, dim%npz), source = undef)
    allocate(arr%va(bd%isd:bd%ied, bd%jsd:bd%jed, dim%npz), source = undef)
    allocate(arr%uc(bd%isd:bd%ied+1, bd%jsd:bd%jed, dim%npz), source = undef)
    allocate(arr%vc(bd%isd:bd%ied, bd%jsd:bd%jed+1, dim%npz), source = undef)
    !# ak, bk
    allocate(arr%ak(dim%npz+1), source = undef)
    allocate(arr%bk(dim%npz+1), source = undef)
    !# mfx, mfy, cx, cy
    allocate(arr%mfx(bd%is:bd%ie+1, bd%js:bd%je, dim%npz), source = undef)
    allocate(arr%mfy(bd%is:bd%ie, bd%js:bd%je+1, dim%npz), source = undef)
    allocate(arr%cx(bd%is:bd%ie+1, bd%jsd:bd%jed, dim%npz), source = undef)
    allocate(arr%cy(bd%isd:bd%ied, bd%js:bd%je+1, dim%npz), source = undef)

    ! Now read data
    open(newunit = file_handle, file = file_name, form = 'unformatted', status = 'old')
    read(file_handle) &
         arr%u, arr%v, arr%w, arr%delz, &
         arr%pt, arr%delp, arr%q, &
         arr%ps, arr%pe, arr%pk, arr%peln, arr%pkz, &
         arr%phis, arr%q_con, arr%omga, &
         arr%ua, arr%va, arr%uc, arr%vc, &
         arr%ak, arr%bk, &
         arr%mfx, arr%mfy, arr%cx, arr%cy
    close(file_handle)

    print *, 'u, v, w, delz: ', sum(arr%u), sum(arr%v), sum(arr%w), sum(arr%delz)
    print *, 'pt, delp, q: ', sum(arr%pt), sum(arr%delp), sum(arr%q)
    print *, 'ps, pe, pk, peln, pkz: ', sum(arr%ps), sum(arr%pe), sum(arr%pk), sum(arr%peln), sum(arr%pkz)
    print *, 'phis, q_con, omga: ', sum(arr%phis), sum(arr%q_con), sum(arr%omga)
    print *, 'ua, va, uc, vc: ', sum(arr%ua), sum(arr%va), sum(arr%uc), sum(arr%vc)
    print *, 'ak, bk: ', sum(arr%ak), sum(arr%bk)
    print *, 'mfx, mfy, cx, cy: ', sum(arr%mfx), sum(arr%mfy), sum(arr%cx), sum(arr%cy)

  end function initialize_and_read_data_

end module input_arrays_mod
