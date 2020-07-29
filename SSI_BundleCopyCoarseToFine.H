
subroutine SUB_(bundle, coarse_Array, index, rc)

   use ESMF
   implicit none

   type(ESMF_FieldBundle) :: bundle
   integer, intent(in) :: index 
   real(TYPEKIND_), intent(inout) :: coarse_Array DIMENSIONS_
   integer, optional, intent(out) :: rc
    
!local
   character(len=ESMF_MAXSTR)      :: IAm='SSI_BundleCopyFineToCoarse'
#include "copy_header.h"
   call ESMF_FieldBundleGet(bundle, fieldIndex=index, field=field, rc=status)
   VERIFY_(STATUS)
   
#include "range_def.h"

    RETURN_(ESMF_SUCCESS)

end subroutine SUB_
