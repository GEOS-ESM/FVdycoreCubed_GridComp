
subroutine SUB_(state, coarse_Array, name, f2c_SSI_arr_map, rc)

   use ESMF
   implicit none

   type(ESMF_State), intent(in) :: state
   real(TYPEKIND_), intent(in) :: coarse_Array DIMENSIONS_
   character(len=*), intent(in) :: name 
   type(SSI_Type), intent(in) :: f2c_SSI_arr_map
   integer, optional, intent(out) :: rc
    
!local
   character(len=ESMF_MAXSTR)      :: IAm='SSI_CopyFineToCoarse'

#include "SSI_Copy_Header.h"

   call ESMF_StateGet(state, trim(name), field, rc=status)
   VERIFY_(STATUS)
   
#include "SSI_Range_Def.h"

    RETURN_(ESMF_SUCCESS)

end subroutine SUB_
