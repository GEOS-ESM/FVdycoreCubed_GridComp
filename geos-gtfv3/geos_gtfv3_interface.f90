module geos_gtfv3_interface_mod

  use iso_c_binding, only: c_int, c_float, c_double

  implicit none

  private
  public :: geos_gtfv3_interface_init_f
  public :: geos_gtfv3_interface_f
  public :: geos_gtfv3_interface_finalize_f

  interface

    subroutine geos_gtfv3_interface_init_f(comm, &
      npx, npy, npz, ntiles, &
      is, ie, js, je, &
      isd, ied, jsd, jed, &
      bdt, nq_tot) bind(c, name='geos_gtfv3_interface_init_c')
      import c_int, c_float, c_double

      implicit none
      integer(kind=c_int), value, intent(in) :: comm
      integer(kind=c_int), value, intent(in) :: npx, npy, npz, ntiles
      integer(kind=c_int), value, intent(in) :: is, ie, js, je
      integer(kind=c_int), value, intent(in) :: isd, ied, jsd, jed
      real(kind=c_float), value, intent(in) :: bdt ! large time step
      integer(kind=c_int), value, intent(in) :: nq_tot ! transported tracers
    end subroutine geos_gtfv3_interface_init_f

     subroutine geos_gtfv3_interface_f( &
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

          ! Input
          ak, bk, &

          ! Input/Output
          mfx, mfy, cx, cy, diss_est &
          ) bind(c, name='geos_gtfv3_interface_c')

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

       ! Input
       real(kind=c_float), dimension(*), intent(in) :: ak, bk

       ! Input/Output
       real(kind=c_float), dimension(*), intent(inout) :: mfx, mfy, cx, cy, diss_est

     end subroutine geos_gtfv3_interface_f

     subroutine geos_gtfv3_interface_finalize_f() bind(c, name='geos_gtfv3_interface_finalize_c')
     end subroutine geos_gtfv3_interface_finalize_f

  end interface

end module geos_gtfv3_interface_mod
