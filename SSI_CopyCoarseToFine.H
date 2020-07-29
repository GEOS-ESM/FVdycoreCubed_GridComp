
subroutine SUB_(state, coarse_Array, name, rc)

   use ESMF
   implicit none

   type(ESMF_State), intent(in) :: state
   real(TYPEKIND_), intent(in) :: coarse_Array DIMENSIONS_
   character(len=*), intent(in) :: name 
   integer, optional, intent(out) :: rc
    
!local
   character(len=ESMF_MAXSTR)      :: IAm='SSI_CopyFineToCoarse'
#include "copy_header.h"

   call ESMF_StateGet(state, trim(name), field, rc=status)
   VERIFY_(STATUS)
   
#include "range_def.h"

    RETURN_(ESMF_SUCCESS)

end subroutine SUB_
