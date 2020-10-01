#include "MAPL_ErrLog.h"
! Copy arrays from fine decomp to coarse decomp

module SSI_CoarseToFine
   use ESMF
   use MAPL

   use SSI_TypeMod, only : SSI_Type
   use fv_timing_mod,       only: timing_on, timing_off

   interface SSI_CopyCoarseToFine
      module procedure SSI_CopyCoarseToFine_R4_2 
      module procedure SSI_CopyCoarseToFine_R8_2 
      module procedure SSI_CopyCoarseToFine_R4_3 
      module procedure SSI_CopyCoarseToFine_R8_3 
   end interface

   interface SSI_BundleCopyCoarseToFine
      module procedure SSI_BundleCopyCoarseToFine_Index_R4_3 
      module procedure SSI_BundleCopyCoarseToFine_Index_R8_3 
      module procedure SSI_BundleCopyCoarseToFine_Name_R4_3 
      module procedure SSI_BundleCopyCoarseToFine_Name_R8_3 
   end interface

   interface SSI_copy_ptr_c2f
      module procedure SSI_copy_ptr_c2f_R4_2
      module procedure SSI_copy_ptr_c2f_R8_2
      module procedure SSI_copy_ptr_c2f_R4_3
      module procedure SSI_copy_ptr_c2f_R8_3
   end interface

   contains

#define IDENTITY(x)x
#define    SUB__(N,A)   SUB___(N,A)
#define SUB___(N,A) IDENTITY(N)IDENTITY(_)IDENTITY(A)
#define NAME_ SSI_CopyCoarseToFine
#define NAME_BUNDLE_INDEX_ SSI_BundleCopyCoarseToFine_Index
#define NAME_BUNDLE_NAME_ SSI_BundleCopyCoarseToFine_Name
#define NAME_COPY_ SSI_copy_ptr_c2f
#define COPY____(P,C,R) IDENTITY(P),IDENTITY(C)IDENTITY(R)

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
#define    SUB_           SUB__(NAME_,TKR_)
#define COPY_ COPY____(farrayPtr,coarse_Array,RANGE_)
#include "SSI_CopyCoarseToFine.H"
#undef SUB_
#define SUB_ SUB__(NAME_COPY_,TKR_)
#include "SSI_copy_ptr_c2f.H"

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
#define    SUB_           SUB__(NAME_,TKR_)
#define COPY_ COPY____(farrayPtr,coarse_Array,RANGE_)
#include "SSI_CopyCoarseToFine.H"
#undef SUB_
#define SUB_ SUB__(NAME_COPY_,TKR_)
#include "SSI_copy_ptr_c2f.H"

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
#define    SUB_           SUB__(NAME_,TKR_)
#define COPY_ COPY____(farrayPtr,coarse_Array,RANGE_)
#include "SSI_CopyCoarseToFine.H"
#undef SUB_
#define SUB_ SUB__(NAME_COPY_,TKR_)
#include "SSI_copy_ptr_c2f.H"

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
#define    SUB_           SUB__(NAME_,TKR_)
#define COPY_ COPY____(farrayPtr,coarse_Array,RANGE_)
#include "SSI_CopyCoarseToFine.H"
#undef SUB_
#define SUB_ SUB__(NAME_COPY_,TKR_)
#include "SSI_copy_ptr_c2f.H"

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
#define COPY_ COPY____(farrayPtr,coarse_Array,RANGE_)
#include "SSI_BundleCopyCoarseToFine.H"

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
#define COPY_ COPY____(farrayPtr,coarse_Array,RANGE_)
#include "SSI_BundleCopyCoarseToFine.H"

end module SSI_CoarseToFine
