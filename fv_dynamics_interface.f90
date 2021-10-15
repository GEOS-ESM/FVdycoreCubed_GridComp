module fv_dynamics_interface_mod

  use iso_c_binding, only: c_int, c_float

  implicit none

  private
  public :: fv_dynamics_interface
  
  interface

     subroutine fv_dynamics_interface( &
          ! Input
          comm, &
          npx, npy, npz, nq_tot, ng, &
          isd, ied, jsd, jed, bdt, &
          consv_te, fill, reproduce_sum, &
          kappa, cp_air, zvir, ptop, ks, ncnst, n_split, q_split, &
          ! Input/Output
          u, v, w &
          ) bind(c, name = 'fv_dynamics_interface')

     ! subroutine fv_dynamics_interface( &
     !      ! Input
     !      comm, &
     !      npx, npy, npz, nq_tot, ng, &
     !      isd, ied, jsd, jed, bdt, &
     !      consv_te, fill, reproduce_sum, &
     !      kappa, cp_air, zvir, ptop, ks, ncnst, n_split, q_split, &
     !      ! Input/Output
     !      u, v, w, delz, &
     !      ! Input
     !      hydrostatic, &
     !      ! Input/Output
     !      pt, delp, q, &
     !      ps, pe, pk, peln, pkz, & ! these 5 can be re-computed from delp and ptop
     !      phis, q_con, omga, ua, va, uc, vc, &
     !      ! Input
     !      ak, bk, &
     !      ! Input/Output
     !      mfx, mfy, cx, cy, &
     !      ze0, &
     !      ! Input
     !      hybrid_z &
     !      ) bind(c, name='fv_dynamics_interface')

       import c_int, c_float

       implicit none

       ! Input
       integer(kind=c_int), value, intent(in) :: comm
       integer(kind=c_int), value, intent(in) :: npx, npy, npz
       integer(kind=c_int), value, intent(in) :: nq_tot ! transported tracers
       integer(kind=c_int), value, intent(in) :: ng
       integer(kind=c_int), value, intent(in) :: isd, ied, jsd, jed
       real(kind=c_float), value, intent(in) :: bdt ! large time step
       real(kind=c_float), value, intent(in) :: consv_te
       logical, value, intent(in) :: fill
       logical, value, intent(in) :: reproduce_sum
       real(kind=c_float), value, intent(in) :: kappa
       real(kind=c_float), value, intent(in) :: cp_air
       real(kind=c_float), value, intent(in) :: zvir
       real(kind=c_float), value, intent(in) :: ptop
       integer(kind=c_int), value, intent(in) :: ks
       integer(kind=c_int), value, intent(in) :: ncnst
       integer(kind=c_int), value, intent(in) :: n_split
       integer(kind=c_int), value, intent(in) :: q_split

       ! Input/Output
       real(kind=c_float), dimension(*), intent(inout) :: u, v, w
       ! real(kind=c_float), dimension(*), intent(inout) :: delz

       ! ! Input
       ! logical, value, intent(in) :: hydrostatic

       ! ! Input/Output
       ! real(kind=c_float), dimension(*), intent(inout) :: pt ! temperature (K)
       ! real(kind=c_float), dimension(*), intent(inout) :: delp ! pressure thickness (pascal)
       ! real(kind=c_float), dimension(*), intent(inout) :: q ! specific humidity and constituents
       ! ! The next 5 can be re-computed from delp and ptop
       ! real(kind=c_float), dimension(*), intent(inout) :: ps ! surface pressure (pascal)
       ! real(kind=c_float), dimension(*), intent(inout) :: pe ! edge pressure (pascal)
       ! real(kind=c_float), dimension(*), intent(inout) :: pk ! pe**kappa
       ! real(kind=c_float), dimension(*), intent(inout) :: peln ! ln(pe)
       ! real(kind=c_float), dimension(*), intent(inout) :: pkz ! finite-volume mean pk
       ! ! The above 5
       ! real(kind=c_float), dimension(*), intent(inout) :: phis ! surface geopotential (g*Z_surf)
       ! real(kind=c_float), dimension(*), intent(inout) :: q_con
       ! real(kind=c_float), dimension(*), intent(inout) :: omga
       ! real(kind=c_float), dimension(*), intent(inout) :: ua, va
       ! real(kind=c_float), dimension(*), intent(inout) :: uc, vc ! c-grid winds
       
       ! ! Input
       ! real(kind=c_float), dimension(*), intent(in) :: ak, bk

       ! ! Input/Output
       ! real(kind=c_float), dimension(*), intent(inout) :: mfx, mfy ! accumulated mass flux
       ! real(kind=c_float), dimension(*), intent(inout) :: cx, cy ! accumulated courant number
       ! ! TODO: Is this used anywhere?
       ! real(kind=c_float), dimension(*), intent(inout) :: ze0 ! height at edges (m); non-hydrostatic

       ! ! Input
       ! logical, intent(in) :: hybrid_z ! used for remapping

     end subroutine fv_dynamics_interface
       
  end interface

end module fv_dynamics_interface_mod
