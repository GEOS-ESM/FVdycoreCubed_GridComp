module GEOS_FV3_UtilitiesMod

  implicit none
  private

  public :: a2d2c


  contains
   
   subroutine A2D2C(U,V,npz,getC)

! Move A-Grid winds/tendencies oriented on lat/lon to the D-grid, or Optionally C-grid cubed-sphere orientation
! will return d-grid winds unless user asks for c-grid winds

      use mpp_domains_mod,   only: mpp_update_domains, mpp_get_boundary, DGRID_NE
      use mpp_parameter_mod, only: AGRID
      use FV_StateMod,       only: fv_atm
      use fv_arrays_mod,     only: REAL4, REAL8, FVPRC
      use fv_mp_mod,         only: is,js,ie,je,isd,jsd,ied,jed,ng
      use sw_core_mod,       only: d2a2c_vect
      implicit none

      real,    intent(INOUT)           :: U(:,:,:)
      real,    intent(INOUT)           :: V(:,:,:)
      integer, intent(   IN)           :: npz
      logical, intent(   IN)           :: getC

!local variables
      integer :: i,j,k, im2,jm2

      real(REAL8) :: ud(is:ie,js:je+1,npz)
      real(REAL8) :: vd(is:ie+1,js:je,npz)

      real(FVPRC) :: ut(isd:ied, jsd:jed)
      real(FVPRC) :: vt(isd:ied, jsd:jed)

      real(REAL8) :: v3(is-1:ie+1,js-1:je+1,3)
      real(REAL8) :: ue(is-1:ie+1,js  :je+1,3)    ! 3D winds at edges
      real(REAL8) :: ve(is  :ie+1,js-1:je+1,3)    ! 3D winds at edges
      real(REAL8), dimension(is:ie):: ut1, ut2, ut3
      real(REAL8), dimension(js:je):: vt1, vt2, vt3

      real(FVPRC) :: uctemp(isd:ied+1,jsd:jed  ,npz)
      real(FVPRC) :: vctemp(isd:ied  ,jsd:jed+1,npz)
      real(FVPRC) ::  utemp(isd:ied  ,jsd:jed+1,npz)
      real(FVPRC) ::  vtemp(isd:ied+1,jsd:jed  ,npz)
      real(FVPRC) :: uatemp(isd:ied,jsd:jed,npz)
      real(FVPRC) :: vatemp(isd:ied,jsd:jed,npz)
!#else
!      real*8 :: uctemp(isd:ied+1,jsd:jed  ,npz)
!!      real*8 :: vctemp(isd:ied  ,jsd:jed+1,npz)
!#endif

      real(FVPRC) :: wbuffer(js:je,npz)
      real(FVPRC) :: sbuffer(is:ie,npz)
      real(FVPRC) :: ebuffer(js:je,npz)
      real(FVPRC) :: nbuffer(is:ie,npz)
      integer     :: npx, npy
      integer :: STATUS

      npx = FV_Atm(1)%npx
      npy = FV_Atm(1)%npy

      uatemp = 0.0d0
      vatemp = 0.0d0

      uatemp(is:ie,js:je,:) = U
      vatemp(is:ie,js:je,:) = V

      im2 = (npx-1)/2
      jm2 = (npy-1)/2

! Cubed-Sphere
      call mpp_update_domains(uatemp, FV_Atm(1)%domain, complete=.false.)
      call mpp_update_domains(vatemp, FV_Atm(1)%domain, complete=.true.)
      do k=1, npz
! Compute 3D wind tendency on A grid
         do j=js-1,je+1
            do i=is-1,ie+1
               v3(i,j,1) = uatemp(i,j,k)*fv_atm(1)%gridstruct%vlon(i,j,1) + vatemp(i,j,k)*fv_atm(1)%gridstruct%vlat(i,j,1)
               v3(i,j,2) = uatemp(i,j,k)*fv_atm(1)%gridstruct%vlon(i,j,2) + vatemp(i,j,k)*fv_atm(1)%gridstruct%vlat(i,j,2)
               v3(i,j,3) = uatemp(i,j,k)*fv_atm(1)%gridstruct%vlon(i,j,3) + vatemp(i,j,k)*fv_atm(1)%gridstruct%vlat(i,j,3)
            enddo
         enddo
! A --> D
! Interpolate to cell edges
         do j=js,je+1
            do i=is-1,ie+1
               ue(i,j,1) = v3(i,j-1,1) + v3(i,j,1)
               ue(i,j,2) = v3(i,j-1,2) + v3(i,j,2)
               ue(i,j,3) = v3(i,j-1,3) + v3(i,j,3)
            enddo
         enddo

         do j=js-1,je+1
            do i=is,ie+1
               ve(i,j,1) = v3(i-1,j,1) + v3(i,j,1)
               ve(i,j,2) = v3(i-1,j,2) + v3(i,j,2)
               ve(i,j,3) = v3(i-1,j,3) + v3(i,j,3)
            enddo
         enddo
! --- E_W edges (for v-wind):
         if ( is==1 ) then
            i = 1
            do j=js,je
               if ( j>jm2 ) then
                  vt1(j) = fv_atm(1)%gridstruct%edge_vect_w(j)*ve(i,j-1,1)+(1.-fv_atm(1)%gridstruct%edge_vect_w(j))*ve(i,j,1)
                  vt2(j) = fv_atm(1)%gridstruct%edge_vect_w(j)*ve(i,j-1,2)+(1.-fv_atm(1)%gridstruct%edge_vect_w(j))*ve(i,j,2)
                  vt3(j) = fv_atm(1)%gridstruct%edge_vect_w(j)*ve(i,j-1,3)+(1.-fv_atm(1)%gridstruct%edge_vect_w(j))*ve(i,j,3)
               else
                  vt1(j) = fv_atm(1)%gridstruct%edge_vect_w(j)*ve(i,j+1,1)+(1.-fv_atm(1)%gridstruct%edge_vect_w(j))*ve(i,j,1)
                  vt2(j) = fv_atm(1)%gridstruct%edge_vect_w(j)*ve(i,j+1,2)+(1.-fv_atm(1)%gridstruct%edge_vect_w(j))*ve(i,j,2)
                  vt3(j) = fv_atm(1)%gridstruct%edge_vect_w(j)*ve(i,j+1,3)+(1.-fv_atm(1)%gridstruct%edge_vect_w(j))*ve(i,j,3)
               endif
            enddo
            do j=js,je
               ve(i,j,1) = vt1(j)
               ve(i,j,2) = vt2(j)
               ve(i,j,3) = vt3(j)
            enddo
         endif
         if ( (ie+1)==npx ) then
            i = npx
            do j=js,je
               if ( j>jm2 ) then
                  vt1(j) = fv_atm(1)%gridstruct%edge_vect_e(j)*ve(i,j-1,1)+(1.-fv_atm(1)%gridstruct%edge_vect_e(j))*ve(i,j,1)
                  vt2(j) = fv_atm(1)%gridstruct%edge_vect_e(j)*ve(i,j-1,2)+(1.-fv_atm(1)%gridstruct%edge_vect_e(j))*ve(i,j,2)
                  vt3(j) = fv_atm(1)%gridstruct%edge_vect_e(j)*ve(i,j-1,3)+(1.-fv_atm(1)%gridstruct%edge_vect_e(j))*ve(i,j,3)
               else
                  vt1(j) = fv_atm(1)%gridstruct%edge_vect_e(j)*ve(i,j+1,1)+(1.-fv_atm(1)%gridstruct%edge_vect_e(j))*ve(i,j,1)
                  vt2(j) = fv_atm(1)%gridstruct%edge_vect_e(j)*ve(i,j+1,2)+(1.-fv_atm(1)%gridstruct%edge_vect_e(j))*ve(i,j,2)
                  vt3(j) = fv_atm(1)%gridstruct%edge_vect_e(j)*ve(i,j+1,3)+(1.-fv_atm(1)%gridstruct%edge_vect_e(j))*ve(i,j,3)
               endif
            enddo
            do j=js,je
               ve(i,j,1) = vt1(j)
               ve(i,j,2) = vt2(j)
               ve(i,j,3) = vt3(j)
            enddo
         endif
! N-S edges (for u-wind):
         if ( js==1 ) then
            j = 1
            do i=is,ie
               if ( i>im2 ) then
                  ut1(i) = fv_atm(1)%gridstruct%edge_vect_s(i)*ue(i-1,j,1)+(1.-fv_atm(1)%gridstruct%edge_vect_s(i))*ue(i,j,1)
                  ut2(i) = fv_atm(1)%gridstruct%edge_vect_s(i)*ue(i-1,j,2)+(1.-fv_atm(1)%gridstruct%edge_vect_s(i))*ue(i,j,2)
                  ut3(i) = fv_atm(1)%gridstruct%edge_vect_s(i)*ue(i-1,j,3)+(1.-fv_atm(1)%gridstruct%edge_vect_s(i))*ue(i,j,3)
               else
                  ut1(i) = fv_atm(1)%gridstruct%edge_vect_s(i)*ue(i+1,j,1)+(1.-fv_atm(1)%gridstruct%edge_vect_s(i))*ue(i,j,1)
                  ut2(i) = fv_atm(1)%gridstruct%edge_vect_s(i)*ue(i+1,j,2)+(1.-fv_atm(1)%gridstruct%edge_vect_s(i))*ue(i,j,2)
                  ut3(i) = fv_atm(1)%gridstruct%edge_vect_s(i)*ue(i+1,j,3)+(1.-fv_atm(1)%gridstruct%edge_vect_s(i))*ue(i,j,3)
               endif
            enddo
            do i=is,ie
               ue(i,j,1) = ut1(i)
               ue(i,j,2) = ut2(i)
               ue(i,j,3) = ut3(i)
            enddo
         endif

         if ( (je+1)==npy ) then
            j = npy
            do i=is,ie
               if ( i>im2 ) then
                  ut1(i) = fv_atm(1)%gridstruct%edge_vect_n(i)*ue(i-1,j,1)+(1.-fv_atm(1)%gridstruct%edge_vect_n(i))*ue(i,j,1)
                  ut2(i) = fv_atm(1)%gridstruct%edge_vect_n(i)*ue(i-1,j,2)+(1.-fv_atm(1)%gridstruct%edge_vect_n(i))*ue(i,j,2)
                  ut3(i) = fv_atm(1)%gridstruct%edge_vect_n(i)*ue(i-1,j,3)+(1.-fv_atm(1)%gridstruct%edge_vect_n(i))*ue(i,j,3)
               else
                  ut1(i) = fv_atm(1)%gridstruct%edge_vect_n(i)*ue(i+1,j,1)+(1.-fv_atm(1)%gridstruct%edge_vect_n(i))*ue(i,j,1)
                  ut2(i) = fv_atm(1)%gridstruct%edge_vect_n(i)*ue(i+1,j,2)+(1.-fv_atm(1)%gridstruct%edge_vect_n(i))*ue(i,j,2)
                  ut3(i) = fv_atm(1)%gridstruct%edge_vect_n(i)*ue(i+1,j,3)+(1.-fv_atm(1)%gridstruct%edge_vect_n(i))*ue(i,j,3)
               endif
            enddo
            do i=is,ie
               ue(i,j,1) = ut1(i)
               ue(i,j,2) = ut2(i)
               ue(i,j,3) = ut3(i)
            enddo
         endif
! Update:
         do j=js,je+1
            do i=is,ie
               ud(i,j,k) = 0.5*( ue(i,j,1)*fv_atm(1)%gridstruct%es(1,i,j,1) +  &
                     ue(i,j,2)*fv_atm(1)%gridstruct%es(2,i,j,1) +  &
                     ue(i,j,3)*fv_atm(1)%gridstruct%es(3,i,j,1) )
            enddo
         enddo
         do j=js,je
            do i=is,ie+1
               vd(i,j,k) = 0.5*( ve(i,j,1)*fv_atm(1)%gridstruct%ew(1,i,j,2) +  &
                     ve(i,j,2)*fv_atm(1)%gridstruct%ew(2,i,j,2) +  &
                     ve(i,j,3)*fv_atm(1)%gridstruct%ew(3,i,j,2) )
            enddo
         enddo

      enddo         ! k-loop

      if (getC) then

         ! Now we have D-Grid winds, need to make call to d2a2c_vect

         utemp = 0.0d0
         vtemp = 0.0d0
         utemp(is:ie,js:je,:) = ud(is:ie,js:je,:)
         vtemp(is:ie,js:je,:) = vd(is:ie,js:je,:)

         ! update shared edges
         call mpp_get_boundary(utemp, vtemp, FV_Atm(1)%domain, &
                               wbuffery=wbuffer, ebuffery=ebuffer, &
                               sbufferx=sbuffer, nbufferx=nbuffer, &
                               gridtype=DGRID_NE, complete=.true. )
         do k=1,npz
            do i=is,ie
               utemp(i,je+1,k) = nbuffer(i,k)
            enddo
            do j=js,je
               vtemp(ie+1,j,k) = ebuffer(j,k)
            enddo
         enddo

         call mpp_update_domains(utemp, vtemp, FV_Atm(1)%domain, gridtype=DGRID_NE, complete=.true.)
         do k=1,npz
            call d2a2c_vect(utemp(:,:,k),  vtemp(:,:,k), &
                      uatemp(:,:,k), vatemp(:,:,k), &
                      uctemp(:,:,k), vctemp(:,:,k), ut, vt, .true., &
                      fv_atm(1)%gridstruct,fv_atm(1)%bd,npx,npy,.false.,0)
         enddo

         U(:,:,:) = uctemp(is:ie,js:je,:)
         V(:,:,:) = vctemp(is:ie,js:je,:)

      else

         ! return d-grid winds
         u = ud(is:ie,js:je,:)
         v = vd(is:ie,js:je,:)

      end if

   end subroutine A2D2C

end module GEOS_FV3_UtilitiesMod
