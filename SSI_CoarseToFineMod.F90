#include "MAPL_ErrLog.h"
! Copy arrays from fine decomp to coarse decomp

module SSI_CoarseToFine
   use ESMF
   use MAPL

   interface SSI_CopyCoarseToFine
      module procedure SSI_CopyCoarseToFine_R4_2 
      module procedure SSI_CopyCoarseToFine_R8_2 
      module procedure SSI_CopyCoarseToFine_R4_3 
      module procedure SSI_CopyCoarseToFine_R8_3 
   end interface

   interface SSI_BundleCopyCoarseToFine
      module procedure SSI_BundleCopyCoarseToFine_R4_3 
      module procedure SSI_BundleCopyCoarseToFine_R8_3 
   end interface

   contains

#define IDENTITY(x)x
#define    SUB__(N,A)   SUB___(N,A)
#define SUB___(N,A) IDENTITY(N)IDENTITY(_)IDENTITY(A)
#define NAME_ SSI_CopyCoarseToFine
#define NAME_BUNDLE_ SSI_BundleCopyCoarseToFine
#define COPY____(P,C,R) IDENTITY(P)IDENTITY(=)IDENTITY(C)IDENTITY(R)

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
#define COPY_ COPY____(farrayPtr,coarse_Array,RANGE_)
#include "SSI_BundleCopyCoarseToFine.H"

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
#define COPY_ COPY____(farrayPtr,coarse_Array,RANGE_)
#include "SSI_BundleCopyCoarseToFine.H"

end module SSI_CoarseToFine
