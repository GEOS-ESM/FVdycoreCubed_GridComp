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
     logical :: hydrostatic, z_tracer, make_nh, fv_debug, reproduce_sum
     logical :: do_sat_adj, do_vort_damp, rf_fast, fill
     integer :: ncnst, n_split, k_split, fv_sg_adj, n_sponge, n_zfilter, nwat
     integer :: hord_tr, hord_tm, hord_dp, hord_mt, hord_vt, nord
     integer :: kord_tm, kord_tr, kord_wz, kord_mt
     real :: d_ext, beta, vtdm4, ke_bg, d_con, d2_bg, d2_bg_k1, d2_bg_k2
     real :: p_fac, a_imp, dddmp, d4_bg, rf_cutoff, tau, consv_te
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
         layout, scalars%adiabatic, &
         !-logical
         scalars%hydrostatic, scalars%z_tracer, scalars%make_nh, &
         scalars%fv_debug, scalars%reproduce_sum, scalars%do_sat_adj, &
         scalars%do_vort_damp, scalars%rf_fast, scalars%fill, &
         !-integer
         scalars%ncnst, &
         scalars%n_split, scalars%k_split, scalars%fv_sg_adj, &
         scalars%n_sponge, scalars%n_zfilter, scalars%nwat, &
         scalars%hord_tr, scalars%hord_tm, scalars%hord_dp, scalars%hord_mt, scalars%hord_vt, &
         scalars%nord, scalars%kord_tm, scalars%kord_tr, scalars%kord_wz, scalars%kord_mt, &
         !-real
         scalars%d_ext, scalars%beta, scalars%vtdm4, scalars%ke_bg, &
         scalars%d_con, scalars%d2_bg, scalars%d2_bg_k1, scalars%d2_bg_k2, &
         scalars%p_fac, scalars%a_imp, scalars%dddmp, scalars%d4_bg, &
         scalars%rf_cutoff, scalars%tau, scalars%consv_te
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
    print *, 'nq_tot: ', self%nq_tot
    print *, 'ng: ', self%ng
    print *, 'ptop: ', self%ptop
    print *, 'ks: ', self%ks
    print *, 'layout: ', '[', self%layout_1, ',', self%layout_2, ']'
    print *, 'adiabatic: ', self%adiabatic
    print *, 'hydrostatic: ', self%hydrostatic
    print *, 'z_tracer: ', self%z_tracer
    print *, 'make_nh: ', self%make_nh
    print *, 'fv_debug: ', self%fv_debug
    print *, 'reproduce_sum: ', self%reproduce_sum
    print *, 'do_sat_adj: ', self%do_sat_adj
    print *, 'do_vort_damp: ', self%do_vort_damp
    print *, 'rf_fast: ', self%rf_fast
    print *, 'fill: ', self%fill
    print *, 'ncnst: ', self%ncnst
    print *, 'n_split: ', self%n_split
    print *, 'k_split: ', self%k_split
    print *, 'fv_sg_adj: ', self%fv_sg_adj
    print *, 'n_sponge: ', self%n_sponge
    print *, 'n_zfilter: ', self%n_zfilter
    print *, 'nwat: ', self%nwat
    print *, 'hord_tr: ', self%hord_tr
    print *, 'hord_tm: ', self%hord_tm
    print *, 'hord_dp: ', self%hord_dp
    print *, 'hord_mt: ', self%hord_mt
    print *, 'hord_vt: ', self%hord_vt
    print *, 'nord: ', self%nord
    print *, 'kord_tm: ', self%kord_tm
    print *, 'kord_tr: ', self%kord_tr
    print *, 'kord_wz: ', self%kord_wz
    print *, 'kord_mt: ', self%kord_mt
    print *, 'd_ext: ', self%d_ext
    print *, 'beta: ', self%beta
    print *, 'vtdm4: ', self%vtdm4
    print *, 'ke_bg: ', self%ke_bg
    print *, 'd_con: ', self%d_con
    print *, 'd2_bg: ', self%d2_bg
    print *, 'd2_bg_k1: ', self%d2_bg_k1
    print *, 'd2_bg_k2: ', self%d2_bg_k2
    print *, 'p_fac: ', self%p_fac
    print *, 'a_imp: ', self%a_imp
    print *, 'dddmp: ', self%dddmp
    print *, 'd4_bg: ', self%d4_bg
    print *, 'rf_cutoff: ', self%rf_cutoff
    print *, 'tau: ', self%tau
    print *, 'consv_te: ', self%consv_te

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
