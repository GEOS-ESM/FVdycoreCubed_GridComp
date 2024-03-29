#include "MAPL_Generic.h"

module fv_regridding_utils

   use ESMF 
   use fv_arrays_mod,     only: fv_atmos_type, fv_grid_type, fv_grid_bounds_type, FVPRC, REAL4, REAL8
   use fv_diagnostics_mod,only: prt_maxmin
   use fv_mp_mod,         only: is_master, ng
   use fv_mapz_mod,       only: mappm
   use mpp_mod,            only: mpp_error, FATAL, NOTE, mpp_broadcast,mpp_npes
   use MAPL

   implicit none

   private

   public remap_scalar
   public fv_rst
   public copy_fv_rst

   real(FVPRC), parameter :: PI           = MAPL_PI_R8
   real(FVPRC), parameter :: OMEGA        = MAPL_OMEGA
   real(FVPRC), parameter :: GRAV         = MAPL_GRAV
   real(FVPRC), parameter :: KAPPA        = MAPL_KAPPA
   real(FVPRC), parameter :: RDGAS        = MAPL_RGAS
   real(FVPRC), parameter :: RVGAS        = MAPL_RVAP
   real(FVPRC), parameter :: CP_AIR       = MAPL_CP
   real(FVPRC), parameter:: zvir = rvgas/rdgas - 1.

   type fv_var
      character(len=128)   :: name
      integer              :: nlev
      integer              :: n_ungrid
      integer              :: rank
      real(FVPRC), allocatable :: ptr2d(:,:) 
      real(FVPRC), allocatable :: ptr3d(:,:,:) 
      real(FVPRC), allocatable :: ptr4d(:,:,:,:)
      contains
         procedure :: alloc_var 
         procedure :: dealloc_var
   end type fv_var

   type fv_rst
      integer               :: ungrid_size
      logical               :: has_edge
      logical               :: has_center
      character(len=1024)   :: file_name
      logical               :: have_descriptor
      type(fv_var), allocatable :: vars(:)
   end type fv_rst
      

contains

 subroutine dealloc_var(this)
    class(fv_var), intent(inout) :: this
    if (this%rank==2) then
       deallocate(this%ptr2d)
    else if (this%rank==3) then
       deallocate(this%ptr3d)
    else if (this%rank==4) then
       deallocate(this%ptr4d)
    end if
 end subroutine dealloc_var

 subroutine alloc_var(this,is,ie,js,je,km,rc)
    class(fv_var), intent(inout) :: this
    integer, intent(in) :: is,ie,js,je
    integer, intent(in), optional :: km
    integer, intent(out), optional :: rc
 
    integer :: status
    integer :: km_use

    if (this%rank==2) then
       allocate(this%ptr2d(is:ie,js:je),source=0.0_FVPRC)
    else if (this%rank==3) then
       if (this%n_ungrid > 0) then
          allocate(this%ptr3d(is:ie,js:je,this%n_ungrid),source=0.0_FVPRC)
       else if (this%n_ungrid == 0) then
          if (present(km)) then 
             km_use = km
          else
             km_use = this%nlev
          end if
          allocate(this%ptr3d(is:ie,js:je,km_use),source=0.0_FVPRC)
       end if
    else if (this%rank == 4) then
       if (present(km)) then 
          km_use = km
       else
          km_use = this%nlev
       end if
       allocate(this%ptr4d(is:ie,js:je,km_use,this%n_ungrid),source=0.0_FVPRC)
    end if
    _RETURN(_SUCCESS)

 end subroutine alloc_var

 subroutine copy_fv_rst(in_rst,out_rst)
  type(fv_rst), pointer, intent(inout) :: in_rst(:)
  type(fv_rst), pointer, intent(inout) :: out_rst(:)
  
  integer :: ifile,ivar
  allocate(out_rst(size(in_rst)) )
  do ifile=1,size(in_rst)
     allocate( out_rst(ifile)%vars(size(in_rst(ifile)%vars) ) )
     out_rst(ifile)%file_name=in_rst(ifile)%file_name
     out_rst(ifile)%have_descriptor=in_rst(ifile)%have_descriptor
     out_rst(ifile)%ungrid_size=in_rst(ifile)%ungrid_size
     out_rst(ifile)%has_center=in_rst(ifile)%has_center
     out_rst(ifile)%has_edge=in_rst(ifile)%has_edge
     do ivar=1,size(in_rst(ifile)%vars)
        out_rst(ifile)%vars(ivar)%name=in_rst(ifile)%vars(ivar)%name
        out_rst(ifile)%vars(ivar)%nlev=in_rst(ifile)%vars(ivar)%nlev
        out_rst(ifile)%vars(ivar)%n_ungrid=in_rst(ifile)%vars(ivar)%n_ungrid
        out_rst(ifile)%vars(ivar)%rank=in_rst(ifile)%vars(ivar)%rank
     enddo
  enddo
   
 end subroutine copy_fv_rst

 subroutine remap_scalar(im, jm, km, npz, nq, ncnst, ak0, bk0, psc, gzc, ta, qa, Atm, in_fv_rst,out_fv_rst)
  type(fv_atmos_type), intent(inout) :: Atm
  integer, intent(in):: im, jm, km, npz, nq, ncnst
  real(FVPRC),    intent(in):: ak0(km+1), bk0(km+1)
  real(FVPRC),    intent(in), dimension(Atm%bd%is:Atm%bd%ie,Atm%bd%js:Atm%bd%je):: psc, gzc
  real(FVPRC),    intent(in), dimension(Atm%bd%is:Atm%bd%ie,Atm%bd%js:Atm%bd%je,km):: ta
  real(FVPRC),    intent(in), dimension(Atm%bd%is:Atm%bd%ie,Atm%bd%js:Atm%bd%je,km,ncnst):: qa
  type(fv_rst), pointer,   intent(inout) :: in_fv_rst(:)
  type(fv_rst), pointer,   intent(inout) :: out_fv_rst(:)
! local:
  real(FVPRC), dimension(Atm%bd%is:Atm%bd%ie,km):: tp
  real(FVPRC), dimension(Atm%bd%is:Atm%bd%ie,km+1):: pe0, pn0
  real(FVPRC), dimension(Atm%bd%is:Atm%bd%ie,npz):: qn1
  real(FVPRC), dimension(Atm%bd%is:Atm%bd%ie,npz+1):: pe1, pn1
  real(FVPRC), dimension(Atm%bd%is:Atm%bd%ie,npz+1):: q_edge_old,q_edge_new
  real(FVPRC) pt0(km), gz(km+1), pk0(km+1)
  real(FVPRC) qp( Atm%bd%is:Atm%bd%ie,km,ncnst)
  real(FVPRC) qp1(Atm%bd%is:Atm%bd%ie,km)
  real(FVPRC) pst, p1, p2, alpha, rdg
  integer i,j,k, iq
  integer  sphum,ifile,ivar,n_ungrid
  integer :: is,  ie,  js,  je
  integer :: isd, ied, jsd, jed
  logical :: doVert

  is  = Atm%bd%is
  ie  = Atm%bd%ie
  js  = Atm%bd%js
  je  = Atm%bd%je
  isd = Atm%bd%isd
  ied = Atm%bd%ied
  jsd = Atm%bd%jsd
  jed = Atm%bd%jed
  sphum   = 1
  if ( sphum/=1 ) then
       call mpp_error(FATAL,'SPHUM must be 1st tracer')
  endif

  do j=js,je
     do i=is,ie

       do iq=1,ncnst
          do k=1,km
             qp(i,k,iq) = qa(i,j,k,iq)
          enddo
       enddo
       do k=1,km
          tp(i,k) = ta(i,j,k)*(1.+zvir*qp(i,k,sphum))
       enddo

! Tracers:

       do k=1,km+1
          pe0(i,k) = ak0(k) + bk0(k)*psc(i,j)
          pn0(i,k) = log(pe0(i,k))
            pk0(k) = pe0(i,k)**kappa
       enddo

! * Adjust interpolated ps to model terrain
       gz(km+1) = gzc(i,j)
       do k=km,1,-1
           gz(k) = gz(k+1) + rdgas*tp(i,k)*(pn0(i,k+1)-pn0(i,k))
       enddo
! Only lowest layer potential temp is needed
          pt0(km) = tp(i,km)/(pk0(km+1)-pk0(km))*(kappa*(pn0(i,km+1)-pn0(i,km)))
       if( Atm%phis(i,j)>gzc(i,j) ) then
           do k=km,1,-1
              if( Atm%phis(i,j) <  gz(k)  .and.    &
                  Atm%phis(i,j) >= gz(k+1) ) then
                  pst = pk0(k) + (pk0(k+1)-pk0(k))*(gz(k)-Atm%phis(i,j))/(gz(k)-gz(k+1))
                  go to 123
              endif
           enddo
       else
! Extrapolation into the ground
           pst = pk0(km+1) + (gzc(i,j)-Atm%phis(i,j))/(cp_air*pt0(km))
       endif

123    Atm%ps(i,j) = pst**(1./kappa)

     enddo   !i-loop

     do i=is,ie
        pe1(i,1) = Atm%ak(1)
        pn1(i,1) = log(pe1(i,1))
     enddo
     do k=2,npz+1
       do i=is,ie
          pe1(i,k) = Atm%ak(k) + Atm%bk(k)*Atm%ps(i,j)
          pn1(i,k) = log(pe1(i,k))
       enddo
     enddo

! * Compute delp
     do k=1,npz
        do i=is,ie
           Atm%delp(i,j,k) = pe1(i,k+1) - pe1(i,k)
        enddo
     enddo

!---------------
! map tracers
!----------------
     do iq=1,ncnst
        call mappm(km, pe0, qp(is,1,iq), npz, pe1,  qn1, is,ie, 0, 11, Atm%ptop)
        do k=1,npz
           do i=is,ie
              Atm%q(i,j,k,iq) = qn1(i,k)
           enddo
        enddo
     enddo

!---------------
! map extra 3d variables
!----------------
     
     do ifile=1,size(out_fv_rst)

        if (out_fv_rst(ifile)%have_descriptor) then
           do ivar=1,size(out_fv_rst(ifile)%vars)
              if (out_fv_rst(ifile)%vars(ivar)%rank==2) then
                 out_fv_rst(ifile)%vars(ivar)%ptr2d(is:ie,j)=in_fv_rst(ifile)%vars(ivar)%ptr2d(is:ie,j)
              else if (out_fv_rst(ifile)%vars(ivar)%rank==3) then
                 if (out_fv_rst(ifile)%vars(ivar)%nLev==npz) then
                    do k=1,in_fv_rst(ifile)%vars(ivar)%nLev
                       qp1(is:ie,k)=in_fv_rst(ifile)%vars(ivar)%ptr3d(is:ie,j,k)
                    enddo
                    call mappm(km, pe0, qp1, npz, pe1,  qn1, is,ie, 0, 11, Atm%ptop)
                    do k=1,npz
                       do i=is,ie
                          out_fv_rst(ifile)%vars(ivar)%ptr3d(i,j,k) = qn1(i,k)
                       enddo
                    enddo
                 else if (out_fv_rst(ifile)%vars(ivar)%nLev==npz+1) then
                    do k=1,in_fv_rst(ifile)%vars(ivar)%nLev
                       q_edge_old(is:ie,k)=in_fv_rst(ifile)%vars(ivar)%ptr3d(is:ie,j,k)
                    enddo
                    call remap_edge(q_edge_old,q_edge_new,is,ie,km,npz,Atm%ak,Atm%bk)
                    do k=1,npz+1
                       do i=is,ie
                          out_fv_rst(ifile)%vars(ivar)%ptr3d(i,j,k) = q_edge_new(i,k)
                       enddo
                    enddo
                 end if
              else if (out_fv_rst(ifile)%vars(ivar)%rank==4) then
                 do n_ungrid=1,out_fv_rst(ifile)%vars(ivar)%n_ungrid
                    if (out_fv_rst(ifile)%vars(ivar)%nLev==npz) then
                       do k=1,in_fv_rst(ifile)%vars(ivar)%nLev
                          qp1(is:ie,k)=in_fv_rst(ifile)%vars(ivar)%ptr4d(is:ie,j,k,n_ungrid)
                       enddo
                       call mappm(km, pe0, qp1, npz, pe1,  qn1, is,ie, 0, 11, Atm%ptop)
                       do k=1,npz
                          do i=is,ie
                             out_fv_rst(ifile)%vars(ivar)%ptr4d(i,j,k,n_ungrid) = qn1(i,k)
                          enddo
                       enddo
                    else if (out_fv_rst(ifile)%vars(ivar)%nLev==npz+1) then
                       do k=1,in_fv_rst(ifile)%vars(ivar)%nLev
                          q_edge_old(is:ie,k)=in_fv_rst(ifile)%vars(ivar)%ptr4d(is:ie,j,k,n_ungrid)
                       enddo
                       call remap_edge(q_edge_old,q_edge_new,is,ie,km,npz,Atm%ak,Atm%bk)
                       do k=1,npz+1
                          do i=is,ie
                             out_fv_rst(ifile)%vars(ivar)%ptr4d(i,j,k,n_ungrid) = q_edge_new(i,k)
                          enddo
                       enddo
                    end if
                 end do
              end if
           enddo
        else
           do ivar=1,size(out_fv_rst(ifile)%vars)
              out_fv_rst(ifile)%vars(ivar)%ptr3d(is:ie,j,:)=in_fv_rst(ifile)%vars(ivar)%ptr3d(is:ie,j,:)
           end do
        end if
     enddo

!-------------------------------------------------------------
! map virtual temperature using geopotential conserving scheme.
!-------------------------------------------------------------
     call mappm(km, pn0, tp, npz, pn1, qn1, is,ie, 1, 9, Atm%ptop)
     do k=1,npz
        do i=is,ie
           Atm%pt(i,j,k) = qn1(i,k)/(1.+zvir*Atm%q(i,j,k,sphum))
        enddo
     enddo

  enddo

  call prt_maxmin('PS_model', Atm%ps, is, ie, js, je, ng, 1, 0.01_FVPRC)

  if (is_master()) write(*,*) 'done remap_scalar'

end subroutine remap_scalar

subroutine remap_edge(q1,q2,is,ie,km,kn,ak,bk)

!  q1,km - old levels
!  q2,kn - new levels
   integer, intent(in) :: is,ie,km,kn
   real(FVPRC),intent(in) :: ak(kn+1), bk(kn+1)
   real(FVPRC),intent(in) :: q1(is:ie,km+1)
   real(FVPRC),intent(out) :: q2(is:ie,kn+1)

   integer i,k
   do i=is,ie
      do k=1,kn+1
         if (bk(k)==0.0) then
            q2(i,k)=0.0
         else
            q2(i,k)=bk(k)*q1(i,km+1)
         end if
      enddo
   enddo

end subroutine remap_edge

end module fv_regridding_utils

