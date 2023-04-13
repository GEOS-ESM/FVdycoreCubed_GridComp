module input_arrays_mod

  use grid_bounds_mod, only: GridBounds_T
  use domain_dim_mod, only: DomainDim_T
  use iso_fortran_env, only: real64

  implicit none

  public InputArrays_T

  type InputArrays_T
     real, allocatable, dimension(:,:,:) :: u, v, w, delz
     real, allocatable :: pt(:,:,:), delp(:,:,:), q(:,:,:,:)
     real, allocatable :: ps(:,:), pe(:,:,:), pk(:,:,:), peln(:,:,:), pkz(:,:,:)
     real, allocatable :: phis(:,:), varflt(:,:), q_con(:,:,:), omga(:,:,:)
     real, allocatable :: ua(:,:,:), va(:,:,:), uc(:,:,:), vc(:,:,:)
     real, allocatable :: ak(:), bk(:)
     real, allocatable :: mfx(:,:,:), mfy(:,:,:), cx(:,:,:), cy(:,:,:)
     real, allocatable :: diss_est(:,:,:), u_dt(:,:,:), v_dt(:,:,:), w_dt(:,:,:), t_dt(:,:,:)
   contains
     procedure :: wr1te
  end type InputArrays_T

  interface InputArrays_T
     procedure :: initialize_and_read_data_
  end interface InputArrays_T

contains

  function initialize_and_read_data_(file_name, bd, dim, num_tracers) result(arr)

    ! Arguments
    character(len=*), intent(in) :: file_name
    type(GridBounds_T), intent(in) :: bd
    type(DomainDim_T), intent(in) :: dim
    integer, intent(in) :: num_tracers
    type(InputArrays_T) :: arr ! output

    ! Locals
    real, parameter :: undef = 1.0e15
    real(kind=real64), parameter :: undef_real64 = 1.0e15
    integer :: file_handle, npz

    ! Start
    npz = dim%npz

    ! Allocate memory
    !# u, v, w, delz
    allocate(arr%u(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz), source = undef)
    allocate(arr%v(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz), source = undef)
    allocate(arr%w(bd%isd:bd%ied, bd%jsd:bd%jed, 1:npz), source = undef) ! ??
    allocate(arr%delz(bd%isd:bd%ied, bd%jsd:bd%jed, 1:npz), source = undef) ! ??
    !# pt, delp, q
    allocate(arr%pt(bd%isd:bd%ied, bd%jsd:bd%jed, npz), source = undef)
    allocate(arr%delp(bd%isd:bd%ied, bd%jsd:bd%jed, npz), source = undef)
    allocate(arr%q(bd%isd:bd%ied, bd%jsd:bd%jed, npz, num_tracers), source = undef)
    !# ps, pe, pk, peln, pkz
    allocate(arr%ps(bd%isd:bd%ied, bd%jsd:bd%jed), source = undef)
    allocate(arr%pe(bd%is-1:bd%ie+1, npz+1, bd%js-1:bd%je+1), source = undef)
    allocate(arr%pk(bd%is:bd%ie, bd%js:bd%je, npz+1), source = undef)
    allocate(arr%peln(bd%is:bd%ie, npz+1, bd%js:bd%je), source = undef)
    allocate(arr%pkz(bd%is:bd%ie, bd%js:bd%je, npz), source = undef)
    !# phis, varflt, q_con, omga
    allocate(arr%phis(bd%isd:bd%ied, bd%jsd:bd%jed), source = undef)
    allocate(arr%varflt(bd%is:bd%ie, bd%js:bd%je), source = undef)
    allocate(arr%q_con(bd%isd:bd%ied, bd%jsd:bd%jed, 1:npz), source = undef) ! ??
    allocate(arr%omga(bd%isd:bd%ied, bd%jsd:bd%jed, npz), source = undef)
    !# ua, va, uc, vc
    allocate(arr%ua(bd%isd:bd%ied, bd%jsd:bd%jed, npz), source = undef)
    allocate(arr%va(bd%isd:bd%ied, bd%jsd:bd%jed, npz), source = undef)
    allocate(arr%uc(bd%isd:bd%ied+1, bd%jsd:bd%jed, npz), source = undef)
    allocate(arr%vc(bd%isd:bd%ied, bd%jsd:bd%jed+1, npz), source = undef)
    !# ak, bk
    allocate(arr%ak(npz+1), source = undef)
    allocate(arr%bk(npz+1), source = undef)
    !# mfx, mfy, cx, cy, diss_est
    allocate(arr%mfx(bd%is:bd%ie+1, bd%js:bd%je, npz), source = undef)
    allocate(arr%mfy(bd%is:bd%ie, bd%js:bd%je+1, npz), source = undef)
    allocate(arr%cx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz), source = undef)
    allocate(arr%cy(bd%isd:bd%ied, bd%js:bd%je+1, npz), source = undef)
    !# diss_est, u_dt, v_dt, w_dt, t_dt
    allocate(arr%diss_est(bd%isd:bd%ied, bd%jsd:bd%jed, npz), source = undef)
    allocate(arr%u_dt(bd%is:bd%ie, bd%js:bd%je, npz), source = undef)
    allocate(arr%v_dt(bd%is:bd%ie, bd%js:bd%je, npz), source = undef)
    allocate(arr%w_dt(bd%is:bd%ie, bd%js:bd%je, npz), source = undef)
    allocate(arr%t_dt(bd%is:bd%ie, bd%js:bd%je, npz), source = undef)

    ! Now read data
    open(newunit = file_handle, file = file_name, form = 'unformatted', status = 'old')
    read(file_handle) &
         arr%u, arr%v, arr%w, arr%delz, &
         arr%pt, arr%delp, arr%q, &
         arr%ps, arr%pe, arr%pk, arr%peln, arr%pkz, &
         arr%phis, arr%q_con, arr%omga, & ! no varflt yet
         arr%ua, arr%va, arr%uc, arr%vc, &
         arr%ak, arr%bk, &
         arr%mfx, arr%mfy, arr%cx, arr%cy, arr%diss_est ! no u/v/w/t_dt
    close(file_handle)
    ! print *, 'phis:', shape(arr%phis), sum(arr%phis), &
    !      arr%phis(3,1), arr%phis(4,1), arr%phis(9,1), arr%phis(14,1), arr%phis(15,1), &
    !      arr%phis(3,7), arr%phis(4,7), arr%phis(9,7), arr%phis(14,7), arr%phis(15,7)


  end function initialize_and_read_data_

  subroutine wr1te(self)

    ! Arguments
    class(InputArrays_T), intent(in) :: self

    ! Start
    print *, 'u, v, w, delz: ', sum(self%u), sum(self%v), sum(self%w), sum(self%delz)
    print *, 'pt, delp, q: ', sum(self%pt), sum(self%delp), sum(self%q)
    print *, 'ps, pe, pk, peln, pkz: ', &
         sum(self%ps), sum(self%pe), sum(self%pk), sum(self%peln), sum(self%pkz), &
         shape(self%pe), shape(self%peln)
    print *, 'phis, varflt, q_con, omga: ', sum(self%phis), sum(self%varflt), sum(self%q_con), sum(self%omga)
    print *, 'q_con shape: ', shape(self%q_con)
    print *, 'ua, va, uc, vc: ', sum(self%ua), sum(self%va), sum(self%uc), sum(self%vc)
    print *, 'ak, bk: ', sum(self%ak), sum(self%bk)
    print *, 'mfx, mfy, cx, cy: ', sum(self%mfx), sum(self%mfy), sum(self%cx), sum(self%cy)
    print *, 'diss_est, u_dt, v_dt, w_dt: ', sum(self%diss_est), sum(self%u_dt), sum(self%v_dt), sum(self%w_dt), sum(self%t_dt)

  end subroutine wr1te

end module input_arrays_mod
