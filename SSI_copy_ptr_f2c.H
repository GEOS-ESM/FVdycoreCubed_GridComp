! copy ptr to array sub-section
subroutine SUB_(coarse_Array, farrayPtr, rc)

   implicit none

   real(TYPEKIND_) :: coarse_Array DIMENSIONS_
   real(TYPEKIND_) :: farrayPtr DIMENSIONS_
   integer, optional, intent(out) :: rc

   character(len=ESMF_MAXSTR) :: Iam = 'SSI_copy_ptr_f2c'
   integer :: status

   coarse_Array = farrayPtr

   RETURN_(ESMF_SUCCESS)

end subroutine SUB_
