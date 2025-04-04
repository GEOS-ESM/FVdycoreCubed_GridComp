module pyfv3_interface_mod

   use iso_c_binding, only: c_int, c_float, c_double, c_bool
   use fv_arrays_mod, only: fv_flags_type

   implicit none

   private
   public :: pyfv3_interface_f_init
   public :: pyfv3_interface_f_run
   public :: pyfv3_interface_f_finalize
   public :: fv_flags_interface_type
   public :: make_fv_flags_C_interop

   !-----------------------------------------------------------------------
   ! Shadow C interoperable config struct for FV. See `fv_arrays.f90` for
   ! the original structure, docs and default values
   !-----------------------------------------------------------------------
   type, bind(c) :: fv_flags_interface_type

      ! Fortran flagstruct
      ! Skipped: character(len=80) :: grid_name = 'Gnomonic'
      ! Skipped: character(len=120):: grid_file = 'Inline'
      integer(kind=c_int) :: grid_type
      integer(kind=c_int) :: hord_mt
      integer(kind=c_int) :: kord_mt
      integer(kind=c_int) :: kord_wz
      integer(kind=c_int) :: hord_vt
      integer(kind=c_int) :: hord_tm
      integer(kind=c_int) :: hord_dp
      integer(kind=c_int) :: kord_tm
      integer(kind=c_int) :: hord_tr
      integer(kind=c_int) :: kord_tr
      real(kind=c_float) :: scale_z
      real(kind=c_float) :: w_max
      real(kind=c_float) :: z_min
      real(kind=c_float) :: lim_fac
      integer(kind=c_int) :: nord
      integer(kind=c_int) :: nord_tr
      real(kind=c_float) :: dddmp
      real(kind=c_float) :: d2_bg
      real(kind=c_float) :: d4_bg
      real(kind=c_float) :: vtdm4
      real(kind=c_float) :: trdm2
      real(kind=c_float) :: d2_bg_k1
      real(kind=c_float) :: d2_bg_k2
      real(kind=c_float) :: d2_divg_max_k1
      real(kind=c_float) :: d2_divg_max_k2
      real(kind=c_float) :: damp_k_k1
      real(kind=c_float) :: damp_k_k2
      integer(kind=c_int) :: n_zs_filter
      integer(kind=c_int) :: nord_zs_filter
      logical(kind=c_bool) :: full_zs_filter
      logical(kind=c_bool) :: rf_fast
      logical(kind=c_bool) :: Beljaars_TOFD
      logical(kind=c_bool) :: consv_am
      logical(kind=c_bool) :: do_sat_adj
      logical(kind=c_bool) :: do_f3d
      logical(kind=c_bool) :: no_dycore
      logical(kind=c_bool) :: convert_ke
      logical(kind=c_bool) :: do_vort_damp
      logical(kind=c_bool) :: use_old_omega
      real(kind=c_float) :: beta
      integer(kind=c_int) :: n_zfilter
      integer(kind=c_int) :: n_sponge
      real(kind=c_float) :: d_ext
      integer(kind=c_int) :: nwat
      logical(kind=c_bool) :: warm_start
      logical(kind=c_bool) :: inline_q
      logical(kind=c_bool) :: adiabatic
      real(kind=c_float) :: shift_fac
      logical(kind=c_bool) :: do_schmidt
      real(kind=c_float) :: stretch_fac
      real(kind=c_float) :: target_lat
      real(kind=c_float) :: target_lon
      logical(kind=c_bool) :: reset_eta
      real(kind=c_float) :: p_fac
      real(kind=c_float) :: a_imp
      real(kind=c_float) :: dz_min
      integer(kind=c_int) :: n_split
      integer(kind=c_int) :: m_split
      integer(kind=c_int) :: k_split
      logical(kind=c_bool) :: use_logp
      integer(kind=c_int) :: q_split
      integer(kind=c_int) :: print_freq
      logical(kind=c_bool) :: write_3d_diags
      integer(kind=c_int) :: npx
      integer(kind=c_int) :: npy
      integer(kind=c_int) :: npz
      integer(kind=c_int) :: npz_rst
      integer(kind=c_int) :: ncnst
      integer(kind=c_int) :: pnats
      integer(kind=c_int) :: dnats
      integer(kind=c_int) :: ntiles
      integer(kind=c_int) :: ndims
      integer(kind=c_int) :: nf_omega
      integer(kind=c_int) :: fv_sg_adj
      integer(kind=c_int) :: na_init
      logical(kind=c_bool) :: nudge_dz
      real(kind=c_float) :: p_ref
      real(kind=c_float) :: dry_mass
      integer(kind=c_int) :: nt_prog
      integer(kind=c_int) :: nt_phys
      real(kind=c_float) :: tau_h2o
      real(kind=c_float) :: delt_max
      real(kind=c_float) :: d_con
      real(kind=c_float) :: ke_bg
      real(kind=c_float) :: consv_te
      real(kind=c_float) :: tau
      real(kind=c_float) :: rf_cutoff
      logical(kind=c_bool) :: filter_phys
      logical(kind=c_bool) :: dwind_2d
      logical(kind=c_bool) :: breed_vortex_inline
      logical(kind=c_bool) :: range_warn
      logical(kind=c_bool) :: fill
      logical(kind=c_bool) :: fill_dp
      logical(kind=c_bool) :: fill_wz
      logical(kind=c_bool) :: check_negative
      logical(kind=c_bool) :: non_ortho
      logical(kind=c_bool) :: moist_phys
      logical(kind=c_bool) :: do_Held_Suarez
      logical(kind=c_bool) :: do_reed_physics
      logical(kind=c_bool) :: reed_cond_only
      logical(kind=c_bool) :: reproduce_sum
      logical(kind=c_bool) :: adjust_dry_mass
      logical(kind=c_bool) :: fv_debug
      logical(kind=c_bool) :: srf_init
      logical(kind=c_bool) :: mountain
      logical(kind=c_bool) :: old_divg_damp
      integer(kind=c_int) :: remap_option
      integer(kind=c_int) :: gmao_remap
      logical(kind=c_bool) :: z_tracer
      logical(kind=c_bool) :: fv_land
      logical(kind=c_bool) :: nudge
      logical(kind=c_bool) :: nudge_ic
      logical(kind=c_bool) :: ncep_ic
      logical(kind=c_bool) :: nggps_ic
      logical(kind=c_bool) :: ecmwf_ic
      logical(kind=c_bool) :: gfs_phil
      logical(kind=c_bool) :: agrid_vel_rst
      logical(kind=c_bool) :: use_new_ncep
      logical(kind=c_bool) :: use_ncep_phy
      logical(kind=c_bool) :: fv_diag_ic
      logical(kind=c_bool) :: external_ic
      logical(kind=c_bool) :: external_eta
      logical(kind=c_bool) :: read_increment
      logical(kind=c_bool) :: do_skeb
      integer(kind=c_int) :: skeb_npass
      ! Skipped: character(len=128) :: res_latlon_dynamics
      ! Skipped: character(len=128) :: res_latlon_tracers
      logical(kind=c_bool) :: hydrostatic
      logical(kind=c_bool) :: phys_hydrostatic
      logical(kind=c_bool) :: use_hydro_pressure
      logical(kind=c_bool) :: do_uni_zfull
      logical(kind=c_bool) :: hybrid_z
      logical(kind=c_bool) :: Make_NH
      logical(kind=c_bool) :: make_hybrid_z
      logical(kind=c_bool) :: nudge_qv
      real(kind=c_float) :: add_noise
      integer(kind=c_int) :: a2b_ord
      integer(kind=c_int) :: c2l_ord
      real(kind=c_float) :: dx_const
      real(kind=c_float) :: dy_const
      real(kind=c_float) :: deglat
      real(kind=c_double) :: deglon_start
      logical(kind=c_bool)  :: adj_mass_vmr
      logical(kind=c_bool) :: compute_coords_locally
      ! Grid information
      integer(kind=c_int) :: layout_x
      integer(kind=c_int) :: layout_y
      ! Magic number
      integer(kind=c_int) :: make_fv_flags_C_interop = 123456789
   end type


   interface

      subroutine pyfv3_interface_f_init( &
         fv_flags, &
         comm, &
         npx, npy, npz, ntiles, &
         is, ie, js, je, isd, ied, jsd, jed, &
         bdt, nq_tot, ak, bk &
         ) bind(c, name='pyfv3_interface_c_init')

         import c_int, c_float, c_double, fv_flags_interface_type

         implicit none
         type(fv_flags_interface_type), intent(in) :: fv_flags
         integer(kind=c_int), value, intent(in) :: comm
         integer(kind=c_int), value, intent(in) :: npx, npy, npz, ntiles
         integer(kind=c_int), value, intent(in) :: is, ie, js, je
         integer(kind=c_int), value, intent(in) :: isd, ied, jsd, jed
         real(kind=c_float), value, intent(in) :: bdt ! large time step
         integer(kind=c_int), value, intent(in) :: nq_tot ! transported tracers
         real(kind=c_float), dimension(*), intent(in) :: ak, bk

      end subroutine pyfv3_interface_f_init

      subroutine pyfv3_interface_f_run( &
      ! Input
         comm, &
         npx, npy, npz, ntiles, &
         is, ie, js, je, &
         isd, ied, jsd, jed, &
         bdt, nq_tot, ng, ptop, ks, layout_1, layout_2, &
         adiabatic, &

      ! Input/Output
         u, v, w, delz, &
         pt, delp, q, &
         ps, pe, pk, peln, pkz, & ! these 5 can be re-computed from delp and ptop
         phis, q_con, omga, ua, va, uc, vc, &

      ! Input/Output
         mfx, mfy, cx, cy, diss_est &
         ) bind(c, name='pyfv3_interface_c_run')

         import c_int, c_float, c_double

         implicit none

         ! Input
         integer(kind=c_int), value, intent(in) :: comm
         integer(kind=c_int), value, intent(in) :: npx, npy, npz, ntiles
         integer(kind=c_int), value, intent(in) :: is, ie, js, je
         integer(kind=c_int), value, intent(in) :: isd, ied, jsd, jed
         real(kind=c_float), value, intent(in) :: bdt ! large time step
         integer(kind=c_int), value, intent(in) :: nq_tot ! transported tracers
         integer(kind=c_int), value, intent(in) :: ng
         real(kind=c_float), value, intent(in) :: ptop
         integer(kind=c_int), value, intent(in) :: ks
         integer(kind=c_int), value, intent(in) :: layout_1, layout_2
         logical, value, intent(in) :: adiabatic

         ! Input/Output
         real(kind=c_float), dimension(*), intent(inout) :: u, v, w, delz
         real(kind=c_float), dimension(*), intent(inout) :: pt, delp, q
         ! These 5 can be re-computed from delp and ptop
         real(kind=c_float), dimension(*), intent(inout) :: ps, pe, pk, peln, pkz
         real(kind=c_float), dimension(*), intent(inout) :: phis, q_con, omga, ua, va, uc, vc

         ! Input/Output
         real(kind=c_float), dimension(*), intent(inout) :: mfx, mfy, cx, cy, diss_est

      end subroutine pyfv3_interface_f_run

      subroutine pyfv3_interface_f_finalize() bind(c, name='pyfv3_interface_c_finalize')
      end subroutine pyfv3_interface_f_finalize

   end interface

contains

   subroutine make_fv_flags_C_interop(fv_flags, layout, c_fv_flags)
      type(fv_flags_type), intent(in) :: fv_flags
      integer, intent(in) :: layout(2)
      type(fv_flags_interface_type), intent(out) :: c_fv_flags

      c_fv_flags%grid_type = fv_flags%grid_type
      c_fv_flags%hord_mt = fv_flags%hord_mt
      c_fv_flags%kord_mt = fv_flags%kord_mt
      c_fv_flags%kord_wz = fv_flags%kord_wz
      c_fv_flags%hord_vt = fv_flags%hord_vt
      c_fv_flags%hord_tm = fv_flags%hord_tm
      c_fv_flags%hord_dp = fv_flags%hord_dp
      c_fv_flags%kord_tm = fv_flags%kord_tm
      c_fv_flags%hord_tr = fv_flags%hord_tr
      c_fv_flags%kord_tr = fv_flags%kord_tr
      c_fv_flags%scale_z = fv_flags%scale_z
      c_fv_flags%w_max = fv_flags%w_max
      c_fv_flags%z_min = fv_flags%z_min
      c_fv_flags%lim_fac = fv_flags%lim_fac
      c_fv_flags%nord = fv_flags%nord
      c_fv_flags%nord_tr = fv_flags%nord_tr
      c_fv_flags%dddmp = fv_flags%dddmp
      c_fv_flags%d2_bg = fv_flags%d2_bg
      c_fv_flags%d4_bg = fv_flags%d4_bg
      c_fv_flags%vtdm4 = fv_flags%vtdm4
      c_fv_flags%trdm2 = fv_flags%trdm2
      c_fv_flags%d2_bg_k1 = fv_flags%d2_bg_k1
      c_fv_flags%d2_bg_k2 = fv_flags%d2_bg_k2
      c_fv_flags%d2_divg_max_k1 = fv_flags%d2_divg_max_k1
      c_fv_flags%d2_divg_max_k2 = fv_flags%d2_divg_max_k2
      c_fv_flags%damp_k_k1 = fv_flags%damp_k_k1
      c_fv_flags%damp_k_k2 = fv_flags%damp_k_k2
      c_fv_flags%n_zs_filter = fv_flags%n_zs_filter
      c_fv_flags%nord_zs_filter = fv_flags%nord_zs_filter
      c_fv_flags%full_zs_filter = fv_flags%full_zs_filter
      c_fv_flags%rf_fast = fv_flags%RF_fast
      c_fv_flags%Beljaars_TOFD = fv_flags%Beljaars_TOFD
      c_fv_flags%consv_am = fv_flags%consv_am
      c_fv_flags%do_sat_adj = fv_flags%do_sat_adj
      c_fv_flags%do_f3d = fv_flags%do_f3d
      c_fv_flags%no_dycore = fv_flags%no_dycore
      c_fv_flags%convert_ke = fv_flags%convert_ke
      c_fv_flags%do_vort_damp = fv_flags%do_vort_damp
      c_fv_flags%use_old_omega = fv_flags%use_old_omega
      c_fv_flags%beta = fv_flags%beta
      c_fv_flags%n_zfilter = fv_flags%n_zfilter
      c_fv_flags%n_sponge = fv_flags%n_sponge
      c_fv_flags%d_ext = fv_flags%d_ext
      c_fv_flags%nwat = fv_flags%nwat
      c_fv_flags%warm_start = fv_flags%warm_start
      c_fv_flags%inline_q = fv_flags%inline_q
      c_fv_flags%adiabatic = fv_flags%adiabatic
      c_fv_flags%shift_fac = fv_flags%shift_fac
      c_fv_flags%do_schmidt = fv_flags%do_schmidt
      c_fv_flags%stretch_fac = fv_flags%stretch_fac
      c_fv_flags%target_lat = fv_flags%target_lat
      c_fv_flags%target_lon = fv_flags%target_lon
      c_fv_flags%reset_eta = fv_flags%reset_eta
      c_fv_flags%p_fac = fv_flags%p_fac
      c_fv_flags%a_imp = fv_flags%a_imp
      c_fv_flags%dz_min = fv_flags%dz_min
      c_fv_flags%n_split = fv_flags%n_split
      c_fv_flags%m_split = fv_flags%m_split
      c_fv_flags%k_split = fv_flags%k_split
      c_fv_flags%use_logp = fv_flags%use_logp
      c_fv_flags%q_split = fv_flags%q_split
      c_fv_flags%print_freq = fv_flags%print_freq
      c_fv_flags%write_3d_diags = fv_flags%write_3d_diags
      c_fv_flags%npx = fv_flags%npx
      c_fv_flags%npy = fv_flags%npy
      c_fv_flags%npz = fv_flags%npz
      c_fv_flags%npz_rst = fv_flags%npz_rst
      c_fv_flags%ncnst = fv_flags%ncnst
      c_fv_flags%pnats = fv_flags%pnats
      c_fv_flags%dnats = fv_flags%dnats
      c_fv_flags%ntiles = fv_flags%ntiles
      c_fv_flags%ndims = fv_flags%ndims
      c_fv_flags%nf_omega = fv_flags%nf_omega
      c_fv_flags%fv_sg_adj = fv_flags%fv_sg_adj
      c_fv_flags%na_init = fv_flags%na_init
      c_fv_flags%nudge_dz = fv_flags%nudge_dz
      c_fv_flags%p_ref = fv_flags%p_ref
      c_fv_flags%dry_mass = fv_flags%dry_mass
      c_fv_flags%nt_prog = fv_flags%nt_prog
      c_fv_flags%nt_phys = fv_flags%nt_phys
      c_fv_flags%tau_h2o = fv_flags%tau_h2o
      c_fv_flags%delt_max = fv_flags%delt_max
      c_fv_flags%d_con = fv_flags%d_con
      c_fv_flags%ke_bg = fv_flags%ke_bg
      c_fv_flags%consv_te = fv_flags%consv_te
      c_fv_flags%tau = fv_flags%tau
      c_fv_flags%rf_cutoff = fv_flags%rf_cutoff
      c_fv_flags%filter_phys = fv_flags%filter_phys
      c_fv_flags%dwind_2d = fv_flags%dwind_2d
      c_fv_flags%breed_vortex_inline = fv_flags%breed_vortex_inline
      c_fv_flags%range_warn = fv_flags%range_warn
      c_fv_flags%fill = fv_flags%fill
      c_fv_flags%fill_dp = fv_flags%fill_dp
      c_fv_flags%fill_wz = fv_flags%fill_wz
      c_fv_flags%check_negative = fv_flags%check_negative
      c_fv_flags%non_ortho = fv_flags%non_ortho
      c_fv_flags%moist_phys = fv_flags%moist_phys
      c_fv_flags%do_Held_Suarez = fv_flags%do_Held_Suarez
      c_fv_flags%do_reed_physics = fv_flags%do_reed_physics
      c_fv_flags%reed_cond_only = fv_flags%reed_cond_only
      c_fv_flags%reproduce_sum = fv_flags%reproduce_sum
      c_fv_flags%adjust_dry_mass = fv_flags%adjust_dry_mass
      c_fv_flags%fv_debug = fv_flags%fv_debug
      c_fv_flags%srf_init = fv_flags%srf_init
      c_fv_flags%mountain = fv_flags%mountain
      c_fv_flags%old_divg_damp = fv_flags%old_divg_damp
      c_fv_flags%remap_option = fv_flags%remap_option
      c_fv_flags%gmao_remap = fv_flags%gmao_remap
      c_fv_flags%z_tracer = fv_flags%z_tracer
      c_fv_flags%fv_land = fv_flags%fv_land
      c_fv_flags%nudge = fv_flags%nudge
      c_fv_flags%nudge_ic = fv_flags%nudge_ic
      c_fv_flags%ncep_ic = fv_flags%ncep_ic
      c_fv_flags%nggps_ic = fv_flags%nggps_ic
      c_fv_flags%ecmwf_ic = fv_flags%ecmwf_ic
      c_fv_flags%gfs_phil = fv_flags%gfs_phil
      c_fv_flags%agrid_vel_rst = fv_flags%agrid_vel_rst
      c_fv_flags%use_new_ncep = fv_flags%use_new_ncep
      c_fv_flags%use_ncep_phy = fv_flags%use_ncep_phy
      c_fv_flags%fv_diag_ic = fv_flags%fv_diag_ic
      c_fv_flags%external_ic = fv_flags%external_ic
      c_fv_flags%external_eta = fv_flags%external_eta
      c_fv_flags%read_increment = fv_flags%read_increment
      c_fv_flags%do_skeb = fv_flags%do_skeb
      c_fv_flags%skeb_npass = fv_flags%skeb_npass
      c_fv_flags%hydrostatic = fv_flags%hydrostatic
      c_fv_flags%phys_hydrostatic = fv_flags%phys_hydrostatic
      c_fv_flags%use_hydro_pressure = fv_flags%use_hydro_pressure
      c_fv_flags%do_uni_zfull = fv_flags%do_uni_zfull
      c_fv_flags%hybrid_z = fv_flags%hybrid_z
      c_fv_flags%Make_NH = fv_flags%Make_NH
      c_fv_flags%make_hybrid_z = fv_flags%make_hybrid_z
      c_fv_flags%nudge_qv = fv_flags%nudge_qv
      c_fv_flags%add_noise = fv_flags%add_noise
      c_fv_flags%a2b_ord = fv_flags%a2b_ord
      c_fv_flags%c2l_ord = fv_flags%c2l_ord
      c_fv_flags%dx_const = fv_flags%dx_const
      c_fv_flags%dy_const = fv_flags%dy_const
      c_fv_flags%deglat = fv_flags%deglat
      c_fv_flags%deglon_start = fv_flags%deglon_start
      c_fv_flags%adj_mass_vmr = fv_flags%adj_mass_vmr
      c_fv_flags%compute_coords_locally = fv_flags%compute_coords_locally
      c_fv_flags%layout_x = layout(1)
      c_fv_flags%layout_y = layout(2)
   end subroutine make_fv_flags_C_interop

end module pyfv3_interface_mod
