! copy array sub-section to ptr
subroutine SUB_(farrayPtr, coarse_Array, rc)

   implicit none

   real(TYPEKIND_), intent(in) :: coarse_Array DIMENSIONS_
   real(TYPEKIND_), pointer :: farrayPtr DIMENSIONS_
   integer, optional, intent(out) :: rc

   character(len=ESMF_MAXSTR) :: Iam = 'SSI_copy_ptr_c2'
   integer :: status


   farrayPtr = coarse_Array

   RETURN_(ESMF_SUCCESS)

end subroutine SUB_
