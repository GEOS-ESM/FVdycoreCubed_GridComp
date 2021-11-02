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
     real, allocatable :: phis(:,:), q_con(:,:,:), omga(:,:,:)
     real, allocatable :: ua(:,:,:), va(:,:,:), uc(:,:,:), vc(:,:,:)
     real, allocatable :: ak(:), bk(:)
     real, allocatable :: mfx(:,:,:), mfy(:,:,:), cx(:,:,:), cy(:,:,:)
     real, allocatable :: diss_est(:,:,:)
     ! gridstruct variables
     real, allocatable :: dx(:,:), dy(:,:), dxa(:,:), dya(:,:), dxc(:,:), dyc(:,:)
     real, allocatable :: rdx(:,:), rdy(:,:), rdxa(:,:), rdya(:,:), rdxc(:,:), rdyc(:,:)
     real, allocatable :: cosa(:,:), cosa_s(:,:)
     real, allocatable :: sina_u(:,:), sina_v(:,:), cosa_u(:,:), cosa_v(:,:)
     real, allocatable :: rsin2(:,:), rsina(:,:), rsin_u(:,:), rsin_v(:,:)
     real, allocatable :: sin_sg(:,:,:), cos_sg(:,:,:)
     real, allocatable :: area(:,:), rarea(:,:), rarea_c(:,:), f0(:,:), fC(:,:)
     real, allocatable :: del6_u(:,:), del6_v(:,:), divg_u(:,:), divg_v(:,:)
     real, allocatable :: agrid(:,:,:), bgrid(:,:,:)
     real(kind=real64), allocatable :: edge_e(:), edge_w(:), edge_n(:), edge_s(:)
     logical :: nested, stretched_grid
     real(kind=real64) :: da_min, da_min_c
   contains
     procedure :: wr1te
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
    integer :: file_handle, npz
    real, parameter :: undef = 1.0e15
    real(kind=real64) :: undef_real64 = 1.0e15

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
    allocate(arr%q(bd%isd:bd%ied, bd%jsd:bd%jed, npz, ncnst), source = undef)
    !# ps, pe, pk, peln, pkz
    allocate(arr%ps(bd%isd:bd%ied, bd%jsd:bd%jed), source = undef)
    allocate(arr%pe(bd%is-1:bd%ie+1, npz+1, bd%js-1:bd%je+1), source = undef)
    allocate(arr%pk(bd%is:bd%ie, bd%js:bd%je, npz+1), source = undef)
    allocate(arr%peln(bd%is:bd%ie, npz+1, bd%js:bd%je), source = undef)
    allocate(arr%pkz(bd%is:bd%ie, bd%js:bd%je, npz), source = undef)
    !# phis, q_con, omga
    allocate(arr%phis(bd%isd:bd%ied, bd%jsd:bd%jed), source = undef)
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
    !# mfx, mfy, cx, cy
    allocate(arr%mfx(bd%is:bd%ie+1, bd%js:bd%je, npz), source = undef)
    allocate(arr%mfy(bd%is:bd%ie, bd%js:bd%je+1, npz), source = undef)
    allocate(arr%cx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz), source = undef)
    allocate(arr%cy(bd%isd:bd%ied, bd%js:bd%je+1, npz), source = undef)
    !# diss_est
    allocate(arr%diss_est(bd%isd:bd%ied, bd%jsd:bd%jed, npz), source = undef)
    !# dx, dy, dxa, dya, dxc, dyc
    allocate(arr%dx(bd%isd:bd%ied, bd%jsd:bd%jed+1), source = undef)
    allocate(arr%dy(bd%isd:bd%ied+1, bd%jsd:bd%jed), source = undef)
    allocate(arr%dxa(bd%isd:bd%ied, bd%jsd:bd%jed), source = undef)
    allocate(arr%dya(bd%isd:bd%ied, bd%jsd:bd%jed), source = undef)
    allocate(arr%dxc(bd%isd:bd%ied+1, bd%jsd:bd%jed), source = undef)
    allocate(arr%dyc(bd%isd:bd%ied, bd%jsd:bd%jed+1), source = undef)
    !# rdx, rdy, rdxa, rdya, rdxc, rdyc
    allocate(arr%rdx(bd%isd:bd%ied, bd%jsd:bd%jed+1), source = undef)
    allocate(arr%rdy(bd%isd:bd%ied+1, bd%jsd:bd%jed), source = undef)
    allocate(arr%rdxa(bd%isd:bd%ied, bd%jsd:bd%jed), source = undef)
    allocate(arr%rdya(bd%isd:bd%ied, bd%jsd:bd%jed), source = undef)
    allocate(arr%rdxc(bd%isd:bd%ied+1, bd%jsd:bd%jed), source = undef)
    allocate(arr%rdyc(bd%isd:bd%ied, bd%jsd:bd%jed+1), source = undef)
    !# cosa, cosa_s
    allocate(arr%cosa(bd%isd:bd%ied+1, bd%jsd:bd%jed+1), source = undef)
    allocate(arr%cosa_s(bd%isd:bd%ied, bd%jsd:bd%jed), source = undef)
    !# sina_u, sina_v, cosa_u, cosa_v
    allocate(arr%sina_u(bd%isd:bd%ied+1, bd%jsd:bd%jed), source = undef)
    allocate(arr%sina_v(bd%isd:bd%ied, bd%jsd:bd%jed+1), source = undef)
    allocate(arr%cosa_u(bd%isd:bd%ied+1, bd%jsd:bd%jed), source = undef)
    allocate(arr%cosa_v(bd%isd:bd%ied, bd%jsd:bd%jed+1), source = undef)
    !# rsin2, rsina, rsin_u, rsin_v
    allocate(arr%rsin2(bd%isd:bd%ied, bd%jsd:bd%jed), source = undef)
    allocate(arr%rsina(bd%is:bd%ie+1, bd%js:bd%je+1), source = undef)
    allocate(arr%rsin_u(bd%isd:bd%ied+1, bd%jsd:bd%jed), source = undef)
    allocate(arr%rsin_v(bd%isd:bd%ied, bd%jsd:bd%jed+1), source = undef)
    !# sin_sg, cos_sg
    allocate(arr%sin_sg(bd%isd:bd%ied, bd%jsd:bd%jed, 9), source = undef)
    allocate(arr%cos_sg(bd%isd:bd%ied, bd%jsd:bd%jed, 9), source = undef)
    !# area, rarea, rarea_c, f0, fC
    allocate(arr%area(bd%isd:bd%ied, bd%jsd:bd%jed), source = undef)
    allocate(arr%rarea(bd%isd:bd%ied, bd%jsd:bd%jed), source = undef)
    allocate(arr%rarea_c(bd%isd:bd%ied+1, bd%jsd:bd%jed+1), source = undef)
    allocate(arr%f0(bd%isd:bd%ied, bd%jsd:bd%jed), source = undef)
    allocate(arr%fC(bd%isd:bd%ied+1, bd%jsd:bd%jed+1), source = undef)
    !# del6_u/v, divg_u/v
    allocate(arr%del6_u(bd%isd:bd%ied, bd%jsd:bd%jed+1), source = undef)
    allocate(arr%del6_v(bd%isd:bd%ied+1, bd%jsd:bd%jed), source = undef)
    allocate(arr%divg_u(bd%isd:bd%ied, bd%jsd:bd%jed+1), source = undef)
    allocate(arr%divg_v(bd%isd:bd%ied+1, bd%jsd:bd%jed), source = undef)
    !# agrid
    allocate(arr%agrid(bd%isd:bd%ied, bd%jsd:bd%jed, 2), source = undef)
    allocate(arr%bgrid(bd%isd:bd%ied+1, bd%jsd:bd%jed+1, 2), source = undef)
    !# edge_e/w/n/s
    allocate(arr%edge_e(dim%npy), source = undef_real64)
    allocate(arr%edge_w(dim%npy), source = undef_real64)
    allocate(arr%edge_n(dim%npx), source = undef_real64)
    allocate(arr%edge_s(dim%npx), source = undef_real64)

    ! Now read data
    open(newunit = file_handle, file = file_name, form = 'unformatted', status = 'old')
    read(file_handle) &
         arr%u, arr%v, arr%w, arr%delz, &
         arr%pt, arr%delp, arr%q, &
         arr%ps, arr%pe, arr%pk, arr%peln, arr%pkz, &
         arr%phis, arr%q_con, arr%omga, &
         arr%ua, arr%va, arr%uc, arr%vc, &
         arr%ak, arr%bk, &
         arr%mfx, arr%mfy, arr%cx, arr%cy, &
         arr%diss_est, &
         arr%dx, arr%dy, arr%dxa, arr%dya, arr%dxc, arr%dyc, &
         arr%rdx, arr%rdy, arr%rdxa, arr%rdya, arr%rdxc, arr%rdyc, &
         arr%cosa, arr%cosa_s, arr%sina_u, arr%sina_v, arr%cosa_u, arr%cosa_v, &
         arr%rsin2, arr%rsina, arr%rsin_u, arr%rsin_v, &
         arr%sin_sg, arr%cos_sg, arr%area, arr%rarea, arr%rarea_c, arr%f0, arr%fC, &
         arr%del6_u, arr%del6_v, arr%divg_u, arr%divg_v, &
         arr%agrid, arr%bgrid, &
         arr%edge_e, arr%edge_w, arr%edge_n, arr%edge_s, &
         arr%nested, arr%stretched_grid, arr%da_min, arr%da_min_c
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
    print *, 'phis, q_con, omga: ', sum(self%phis), sum(self%q_con), sum(self%omga)
    print *, 'q_con shape: ', shape(self%q_con)
    print *, 'ua, va, uc, vc: ', sum(self%ua), sum(self%va), sum(self%uc), sum(self%vc)
    print *, 'ak, bk: ', sum(self%ak), sum(self%bk)
    print *, 'mfx, mfy, cx, cy: ', sum(self%mfx), sum(self%mfy), sum(self%cx), sum(self%cy)
    print *, 'diss_est: ', sum(self%diss_est)
    print *, 'dx: ', &
         sum(self%dx), sum(self%dy), sum(self%dxa), sum(self%dya), sum(self%dxc), sum(self%dyc)
    print *, 'rdx: ', &
         sum(self%rdx), sum(self%rdy), sum(self%rdxa), &
         sum(self%rdya), sum(self%rdxc), sum(self%rdyc)
    print *, 'cosa: ', sum(self%cosa), sum(self%cosa_s), &
         sum(self%sina_u), sum(self%sina_v), sum(self%cosa_u), sum(self%cosa_v)
    print *, 'rsin2: ', sum(self%rsin2), sum(self%rsina), sum(self%rsin_u), sum(self%rsin_v)
    print *, 'sin_sg: ', sum(self%sin_sg), sum(self%cos_sg)
    print *, 'area: ', sum(self%area), sum(self%rarea), sum(self%rarea_c), &
         sum(self%f0), sum(self%fC)
    print *, 'del6_u: ', sum(self%del6_u), sum(self%del6_v), sum(self%divg_u), sum(self%divg_v)
    print *, 'agrid, bgrid: ', sum(self%agrid), sum(self%bgrid)
    print *, 'edge_e/w/n/s: ', &
         shape(self%edge_e), shape(self%edge_w), shape(self%edge_n), shape(self%edge_s), &
         sum(self%edge_e), sum(self%edge_w), sum(self%edge_n), sum(self%edge_s)
    print *, 'nested: ', self%nested, self%stretched_grid, self%da_min, self%da_min_c
    print *, 'delp:', self%delp(4,10,14), self%delp(5,10,14), self%delp(4,11,14)
    print *, 'delz:', self%delz(-2,10,4), self%delz(-1,10,4), self%delz(-2,11,4)
    print *, 'delz:', self%delz(4,10,14), self%delz(5,10,14), self%delz(4,11,14)

  end subroutine wr1te

end module input_arrays_mod
