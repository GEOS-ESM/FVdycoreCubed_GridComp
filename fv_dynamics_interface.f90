module fv_dynamics_interface_mod

  use iso_c_binding, only: c_int, c_float, c_double

  implicit none

  private
  public :: fv_dynamics_interface_f

  interface

     subroutine fv_dynamics_interface_f( &
          ! Input
          comm, &
          npx, npy, npz, &
          is, ie, js, je, &
          isd, ied, jsd, jed, &
          bdt, nq_tot, ng, ptop, ks, layout_1, layout_2, &
          adiabatic, &
          ! flagstruct data
          hydrostatic, z_tracer, make_nh, fv_debug, &
          reproduce_sum, do_sat_adj, do_vort_damp, rf_fast, fill, &
          ncnst, n_split, k_split, fv_sg_adj, n_sponge, n_zfilter, nwat, &
          hord_tr, hord_tm, hord_dp, hord_mt, hord_vt, &
          nord, kord_tm, kord_tr, kord_wz, kord_mt, &
          d_ext, beta, vtdm4, ke_bg, d_con, d2_bg, d2_bg_k1, d2_bg_k2, &
          p_fac, a_imp, dddmp, d4_bg, rf_cutoff, tau, consv_te, &

          ! Input/Output
          u, v, w, delz, &
          pt, delp, q, &
          ps, pe, pk, peln, pkz, & ! these 5 can be re-computed from delp and ptop
          phis, q_con, omga, ua, va, uc, vc, &

          ! Input
          ak, bk, &

          ! Input/Output
          mfx, mfy, cx, cy, diss_est, &

          ! Input (gridstruct data)
          dx, dy, dxa, dya, dxc, dyc, rdx, rdy, rdxa, rdya, rdxc, rdyc, &
          cosa, cosa_s, sina_u, sina_v, cosa_u, cosa_v, rsin2, rsina, rsin_u, rsin_v, &
          sin_sg, cos_sg, &
          area, rarea, rarea_c, f0, fC, del6_u, del6_v, divg_u, divg_v, &
          agrid, bgrid, a11, a12, a21, a22, &
          edge_e, edge_w, edge_n, edge_s, nested, stretched_grid, da_min, da_min_c &
          ) bind(c, name='fv_dynamics_interface_c')

       import c_int, c_float, c_double

       implicit none

       ! Input
       integer(kind=c_int), value, intent(in) :: comm
       integer(kind=c_int), value, intent(in) :: npx, npy, npz
       integer(kind=c_int), value, intent(in) :: is, ie, js, je
       integer(kind=c_int), value, intent(in) :: isd, ied, jsd, jed
       real(kind=c_float), value, intent(in) :: bdt ! large time step
       integer(kind=c_int), value, intent(in) :: nq_tot ! transported tracers
       integer(kind=c_int), value, intent(in) :: ng
       real(kind=c_float), value, intent(in) :: ptop
       integer(kind=c_int), value, intent(in) :: ks
       integer(kind=c_int), value, intent(in) :: layout_1, layout_2
       logical, value, intent(in) :: adiabatic

       logical, value, intent(in) :: hydrostatic
       logical, value, intent(in) :: z_tracer
       logical, value, intent(in) :: make_nh
       logical, value, intent(in) :: fv_debug
       logical, value, intent(in) :: reproduce_sum
       logical, value, intent(in) :: do_sat_adj
       logical, value, intent(in) :: do_vort_damp
       logical, value, intent(in) :: rf_fast
       logical, value, intent(in) :: fill

       integer(kind=c_int), value, intent(in) :: ncnst
       integer(kind=c_int), value, intent(in) :: n_split
       integer(kind=c_int), value, intent(in) :: k_split
       integer(kind=c_int), value, intent(in) :: fv_sg_adj
       integer(kind=c_int), value, intent(in) :: n_sponge
       integer(kind=c_int), value, intent(in) :: n_zfilter
       integer(kind=c_int), value, intent(in) :: nwat
       integer(kind=c_int), value, intent(in) :: hord_tr
       integer(kind=c_int), value, intent(in) :: hord_tm
       integer(kind=c_int), value, intent(in) :: hord_dp
       integer(kind=c_int), value, intent(in) :: hord_mt
       integer(kind=c_int), value, intent(in) :: hord_vt
       integer(kind=c_int), value, intent(in) :: nord
       integer(kind=c_int), value, intent(in) :: kord_tm
       integer(kind=c_int), value, intent(in) :: kord_tr
       integer(kind=c_int), value, intent(in) :: kord_wz
       integer(kind=c_int), value, intent(in) :: kord_mt
       real(kind=c_float), value, intent(in) :: d_ext
       real(kind=c_float), value, intent(in) :: beta
       real(kind=c_float), value, intent(in) :: vtdm4
       real(kind=c_float), value, intent(in) :: ke_bg
       real(kind=c_float), value, intent(in) :: d_con
       real(kind=c_float), value, intent(in) :: d2_bg
       real(kind=c_float), value, intent(in) :: d2_bg_k1
       real(kind=c_float), value, intent(in) :: d2_bg_k2
       real(kind=c_float), value, intent(in) :: p_fac
       real(kind=c_float), value, intent(in) :: a_imp
       real(kind=c_float), value, intent(in) :: dddmp
       real(kind=c_float), value, intent(in) :: d4_bg
       real(kind=c_float), value, intent(in) :: rf_cutoff
       real(kind=c_float), value, intent(in) :: tau
       real(kind=c_float), value, intent(in) :: consv_te

       ! Input/Output
       real(kind=c_float), dimension(*), intent(inout) :: u, v, w, delz
       real(kind=c_float), dimension(*), intent(inout) :: pt, delp, q
       ! These 5 can be re-computed from delp and ptop
       real(kind=c_float), dimension(*), intent(inout) :: ps, pe, pk, peln, pkz
       real(kind=c_float), dimension(*), intent(inout) :: phis, q_con, omga, ua, va, uc, vc

       ! Input
       real(kind=c_float), dimension(*), intent(in) :: ak, bk

       ! Input/Output
       real(kind=c_float), dimension(*), intent(inout) :: mfx, mfy, cx, cy, diss_est

       ! Input
       real(kind=c_float), dimension(*), intent(in) :: dx, dy, dxa, dya, dxc, dyc
       real(kind=c_float), dimension(*), intent(in) :: rdx, rdy, rdxa, rdya, rdxc, rdyc
       real(kind=c_float), dimension(*), intent(in) :: cosa, cosa_s, sina_u, sina_v, cosa_u, cosa_v
       real(kind=c_float), dimension(*), intent(in) :: rsin2, rsina, rsin_u, rsin_v, sin_sg, cos_sg
       real(kind=c_float), dimension(*), intent(in) :: area, rarea, rarea_c, f0, fC
       real(kind=c_float), dimension(*), intent(in) :: del6_u, del6_v, divg_u, divg_v
       real(kind=c_float), dimension(*), intent(in) :: agrid, bgrid, a11, a12, a21, a22
       real(kind=c_double), dimension(*), intent(in) :: edge_e, edge_w, edge_n, edge_s
       logical, value, intent(in) :: nested, stretched_grid
       real(kind=c_double), value, intent(in) :: da_min, da_min_c

     end subroutine fv_dynamics_interface_f

  end interface

end module fv_dynamics_interface_mod
