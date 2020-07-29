#include "MAPL_ErrLog.h"
! Copy arrays from fine decomp to coarse decomp

module SSI_FineToCoarse
   use ESMF
   use MAPL

   private 

   public :: SSI_CopyFineToCoarse, SSI_BundleCopyFineToCoarse

   interface SSI_CopyFineToCoarse
      module procedure SSI_CopyFineToCoarse_R4_2 
      module procedure SSI_CopyFineToCoarse_R8_2 
      module procedure SSI_CopyFineToCoarse_R4_3 
      module procedure SSI_CopyFineToCoarse_R8_3 
   end interface

   interface SSI_BundleCopyFineToCoarse
      module procedure SSI_BundleCopyFineToCoarse_R4_3 
      module procedure SSI_BundleCopyFineToCoarse_R8_3 
   end interface

   contains

#define IDENTITY(x)x
#define    SUB__(N,A)   SUB___(N,A)
#define SUB___(N,A) IDENTITY(N)IDENTITY(_)IDENTITY(A)
#define NAME_ SSI_CopyFineToCoarse
#define NAME_BUNDLE_ SSI_BundleCopyFineToCoarse
#define COPY____(C,R,P) IDENTITY(C)IDENTITY(R)IDENTITY(=)IDENTITY(P)

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

!-----BUNDLE---------------------
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
#define    SUB_           SUB__(NAME_BUNDLE_,TKR_)
#define COPY_ COPY____(coarse_Array,RANGE_,farrayPtr)
#include "SSI_BundleCopyFineToCoarse.H"

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
#define    SUB_           SUB__(NAME_BUNDLE_,TKR_)
#define COPY_ COPY____(coarse_Array,RANGE_,farrayPtr)
#include "SSI_BundleCopyFineToCoarse.H"

end module SSI_FineToCoarse
