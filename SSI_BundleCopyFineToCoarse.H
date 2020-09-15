
subroutine SUB_INDEX_(bundle, coarse_Array, index, f2c_SSI_arr_map, rc)

   use ESMF
   implicit none

   type(ESMF_FieldBundle) :: bundle
   integer, intent(in) :: index 
   real(TYPEKIND_), intent(inout) :: coarse_Array DIMENSIONS_
   type(SSI_Type), intent(in) :: f2c_SSI_arr_map
   integer, optional, intent(out) :: rc
    
!local
   character(len=ESMF_MAXSTR)      :: IAm='SSI_BundleCopyFineToCoarse_Index'
#include "SSI_Copy_Header.h"

   call ESMF_FieldBundleGet(bundle, fieldIndex=index, field=field, rc=status)
   VERIFY_(STATUS)
   
#include "SSI_Range_Def.h"

    RETURN_(ESMF_SUCCESS)

end subroutine SUB_INDEX_

subroutine SUB_NAME_(bundle, coarse_Array, name, f2c_SSI_arr_map, rc)

   use ESMF
   implicit none

   type(ESMF_FieldBundle) :: bundle
   character(len=ESMF_MAXSTR), intent(in) :: name
   real(TYPEKIND_), intent(inout) :: coarse_Array DIMENSIONS_
   type(SSI_Type), intent(in) :: f2c_SSI_arr_map
   integer, optional, intent(out) :: rc
    
!local
   character(len=ESMF_MAXSTR)      :: IAm='SSI_BundleCopyFineToCoarse_Name'
#include "SSI_Copy_Header.h"

   call ESMF_FieldBundleGet(bundle, fieldName=trim(name), field=field, rc=status)
   VERIFY_(STATUS)
   
#include "SSI_Range_Def.h"

    RETURN_(ESMF_SUCCESS)

end subroutine SUB_NAME_
