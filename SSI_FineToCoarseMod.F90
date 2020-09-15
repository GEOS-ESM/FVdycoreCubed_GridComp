#include "MAPL_ErrLog.h"
! Copy arrays from fine decomp to coarse decomp

module SSI_FineToCoarse
   use ESMF
   use MAPL

   use SSI_TypeMod, only : SSI_Type

   private 

   public :: SSI_CopyFineToCoarse, SSI_BundleCopyFineToCoarse
   public :: SSI_copy_ptr_f2c
!   public :: SSI_StateSync

   interface SSI_CopyFineToCoarse
      module procedure SSI_CopyFineToCoarse_R4_2 
      module procedure SSI_CopyFineToCoarse_R8_2 
      module procedure SSI_CopyFineToCoarse_R4_3 
      module procedure SSI_CopyFineToCoarse_R8_3 
   end interface

   interface SSI_BundleCopyFineToCoarse
      module procedure SSI_BundleCopyFineToCoarse_Index_R4_3 
      module procedure SSI_BundleCopyFineToCoarse_Index_R8_3 
      module procedure SSI_BundleCopyFineToCoarse_Name_R4_3 
      module procedure SSI_BundleCopyFineToCoarse_Name_R8_3 
   end interface

   interface SSI_copy_ptr_f2c
      module procedure SSI_copy_ptr_f2c_R4_2 
      module procedure SSI_copy_ptr_f2c_R8_2 
      module procedure SSI_copy_ptr_f2c_R4_3 
      module procedure SSI_copy_ptr_f2c_R8_3 
   end interface

   contains

#define IDENTITY(x)x
#define    SUB__(N,A)   SUB___(N,A)
#define SUB___(N,A) IDENTITY(N)IDENTITY(_)IDENTITY(A)
#define NAME_ SSI_CopyFineToCoarse
#define NAME_BUNDLE_INDEX_ SSI_BundleCopyFineToCoarse_Index
#define NAME_BUNDLE_NAME_ SSI_BundleCopyFineToCoarse_Name
#define NAME_COPY_ SSI_copy_ptr_f2c
#define COPY____(C,R,P) IDENTITY(C)IDENTITY(R),IDENTITY(P)

!--------------------------
#undef TKR_
#undef DIMENSIONS_
#undef TYPEKIND_
#undef SUB_
#undef RANGE_
#define TKR_ R4_2
#define DIMENSIONS_ (:,:)
#define RANGE_ (is:ie,js:je)
#define TYPEKIND_ ESMF_KIND_R4
#define SUB_ SUB__(NAME_,TKR_)
#define COPY_ COPY____(coarse_Array,RANGE_,farrayPtr)
#include "SSI_CopyFineToCoarse.H"
#undef SUB_
#define SUB_ SUB__(NAME_COPY_,TKR_)
#include "SSI_copy_ptr_f2c.H"

!--------------------------
#undef TKR_
#undef DIMENSIONS_
#undef TYPEKIND_
#undef SUB_
#undef RANGE_
#define TKR_ R8_2
#define DIMENSIONS_ (:,:)
#define RANGE_ (is:ie,js:je)
#define TYPEKIND_ ESMF_KIND_R8
#define SUB_ SUB__(NAME_,TKR_)
#define COPY_ COPY____(coarse_Array,RANGE_,farrayPtr)
#include "SSI_CopyFineToCoarse.H"
#undef SUB_
#define SUB_ SUB__(NAME_COPY_,TKR_)
#include "SSI_copy_ptr_f2c.H"

!--------------------------
#undef TKR_
#undef DIMENSIONS_
#undef TYPEKIND_
#undef SUB_
#undef RANGE_
#define TKR_ R4_3
#define DIMENSIONS_ (:,:,:)
#define RANGE_ (is:ie,js:je,1:km)
#define TYPEKIND_ ESMF_KIND_R4
#define SUB_ SUB__(NAME_,TKR_)
#define COPY_ COPY____(coarse_Array,RANGE_,farrayPtr)
#include "SSI_CopyFineToCoarse.H"
#undef SUB_
#define SUB_ SUB__(NAME_COPY_,TKR_)
#include "SSI_copy_ptr_f2c.H"

!--------------------------
#undef TKR_
#undef DIMENSIONS_
#undef TYPEKIND_
#undef SUB_
#undef RANGE_
#define TKR_ R8_3
#define DIMENSIONS_ (:,:,:)
#define RANGE_ (is:ie,js:je,1:km)
#define TYPEKIND_ ESMF_KIND_R8
#define SUB_ SUB__(NAME_,TKR_)
#define COPY_ COPY____(coarse_Array,RANGE_,farrayPtr)
#include "SSI_CopyFineToCoarse.H"
#undef SUB_
#define SUB_ SUB__(NAME_COPY_,TKR_)
#include "SSI_copy_ptr_f2c.H"

!-----BUNDLE---------------------
!--------------------------
#undef TKR_
#undef DIMENSIONS_
#undef TYPEKIND_
#undef SUB_INDEX_
#undef SUB_NAME_
#undef RANGE_
#define TKR_ R4_3
#define DIMENSIONS_ (:,:,:)
#define RANGE_ (is:ie,js:je,1:km)
#define TYPEKIND_ ESMF_KIND_R4
#define    SUB_INDEX_           SUB__(NAME_BUNDLE_INDEX_,TKR_)
#define    SUB_NAME_           SUB__(NAME_BUNDLE_NAME_,TKR_)
#define COPY_ COPY____(coarse_Array,RANGE_,farrayPtr)
#include "SSI_BundleCopyFineToCoarse.H"

!--------------------------
#undef TKR_
#undef DIMENSIONS_
#undef TYPEKIND_
#undef SUB_INDEX_
#undef SUB_NAME_
#undef RANGE_
#define TKR_ R8_3
#define DIMENSIONS_ (:,:,:)
#define RANGE_ (is:ie,js:je,1:km)
#define TYPEKIND_ ESMF_KIND_R8
#define    SUB_INDEX_           SUB__(NAME_BUNDLE_INDEX_,TKR_)
#define    SUB_NAME_           SUB__(NAME_BUNDLE_NAME_,TKR_)
#define COPY_ COPY____(coarse_Array,RANGE_,farrayPtr)
#include "SSI_BundleCopyFineToCoarse.H"

!subroutine SSI_StateSync(state, rc)
!
!   !use ESMF
!   !use MAPL
!   implicit none
!
!   type(ESMF_State), intent(inout) :: state
!   integer, optional, intent(out) :: rc
!
!!local
!   type(ESMF_Field) :: field
!   character(len=ESMF_MAXSTR)      :: IAm='SSI_StateSync'
!   integer :: status
!   integer :: itemCount
!   character(len=ESMF_MAXSTR), allocatable :: itemNameList(:)
!   type(ESMF_FieldBundle) :: tradv
!   integer :: i, ii, numTracers, dimCount
!   type(ESMF_StateIntent_Flag) :: stateintent
!   real, pointer :: temp2d(:,:)
!   real, pointer :: temp3d(:,:,:)
!   type(ESMF_FieldStatus_Flag) :: field_status
!
!   call ESMF_StateGet(state, itemCount=itemCount, rc=status)
!   VERIFY_(STATUS)
!   allocate(itemNameList(itemCount), stat=status)
!   VERIFY_(STATUS)
!   call ESMF_StateGet(state, itemNameList=itemNameList, rc=status)
!   VERIFY_(STATUS)
!   do i = 1, itemCount
!      if(trim(itemNameList(i))=='AK' .or. trim(itemNameList(i))=='BK' .or. trim(itemNameList(i))=='PREF') cycle
!      !print *, __FILE__, i, itemCount,trim(itemNameList(i))
!      if(trim(itemNameList(i))=='TRADV') then
!         call ESMF_StateGet(state, trim(itemNameList(i)), tradv, rc=status)
!         VERIFY_(STATUS)
!         call ESMF_FieldBundleGet(tradv,fieldCount=numTracers,rc=status)
!         VERIFY_(STATUS)
!         do ii=1,numTracers
!            call ESMF_FieldBundleGet(tradv,fieldIndex=ii,field=field,rc=status)
!            VERIFY_(status)
!            call SSI_FieldSync(field, rc=status)
!            VERIFY_(status)
!         enddo
!      else
!         call ESMF_StateGet(state, trim(itemNameList(i)), field, rc=status)
!         VERIFY_(STATUS)
!         call ESMF_StateGet(state, stateintent=stateintent, rc=status)
!         VERIFY_(STATUS)
!         if(stateintent == ESMF_STATEINTENT_EXPORT) then
!            call ESMF_FieldGet(field, status=field_status, rc=status)
!            VERIFY_(status)
!            if(field_status == ESMF_FIELDSTATUS_COMPLETE) then
!               call ESMF_FieldGet(field, dimCount=dimCount, rc=status)
!               VERIFY_(status)
!               if(dimCount == 2) then
!                  call ESMF_FieldGet(field, farrayPtr=temp2d, rc=status)
!                  VERIFY_(STATUS)
!                  if(associated(temp2d)) then
!                     call SSI_FieldSync(field, rc=status)
!                  endif
!               else
!                  call ESMF_FieldGet(field, farrayPtr=temp3d, rc=status)
!                  VERIFY_(STATUS)
!                  if(associated(temp3d)) then
!                     call SSI_FieldSync(field, rc=status)
!                  endif
!               endif  !dimCount == 2
!            endif   !field_status == ESMF_FIELDSTATUS_COMPLETE
!         else
!            call SSI_FieldSync(field, rc=status)
!            VERIFY_(status)
!         endif  !stateintent == ESMF_STATEINTENT_EXPORT 
!      endif !trim(itemNameList(i))=='TRADV'
!   enddo
!   !print *, '===================================================='
!
!   RETURN_(ESMF_SUCCESS)
!
!   contains
!
!   subroutine SSI_FieldSync(field, rc)
!      type(ESMF_Field), intent(inout) :: field
!      integer, optional, intent(out) :: rc
!
!      integer, allocatable :: arrayImg(:)
!      integer              :: ssiLocalDeCount
!      type(ESMF_Array) :: array
!      character(len=ESMF_MAXSTR)      :: IAm='SSI_FieldSync'
!      integer :: status
!
!      call ESMF_AttributeGet(field, name='SSI_ARRAY_SIZE', &
!           value=ssiLocalDeCount, rc=status)
!      VERIFY_(STATUS)
!      allocate(arrayImg(ssiLocalDeCount), stat=status)
!      VERIFY_(STATUS)
!      call ESMF_AttributeGet(field, name='SSI_ARRAY_SAVED', &
!           valueList=arrayImg, rc=status)
!      VERIFY_(STATUS)
!      array = transfer(arrayimg,array)
!      call ESMF_ArraySync(array,rc=status)
!      VERIFY_(STATUS)
!      RETURN_(ESMF_SUCCESS)
!   end subroutine SSI_FieldSync
!
!end subroutine SSI_StateSync

end module SSI_FineToCoarse
