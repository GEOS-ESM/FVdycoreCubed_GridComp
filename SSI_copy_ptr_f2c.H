! copy ptr to array sub-section
subroutine SUB_(coarse_Array, farrayPtr, rc)

   implicit none

   real(TYPEKIND_) :: coarse_Array DIMENSIONS_
   real(TYPEKIND_), pointer :: farrayPtr DIMENSIONS_
   integer, optional, intent(out) :: rc

   character(len=ESMF_MAXSTR) :: Iam = 'SSI_copy_ptr_f2c'
   integer :: status
   !integer :: ndim, ndim1
   !integer, allocatable :: arrsize(:), arrsize1(:)
   !integer, save :: i=1


   coarse_Array = farrayPtr

         !ndim = size(shape(farrayPtr))
         !ndim1 = size(shape(coarse_Array))
         !allocate(arrsize(ndim))
         !allocate(arrsize1(ndim1))
         !arrsize = shape(farrayPtr)
         !arrsize1 = shape(coarse_Array)
         !print *, __FILE__, arrsize, 'i= ', i
         !print *, __FILE__, arrsize1, 'i= ', i
         !i = i+1
   RETURN_(ESMF_SUCCESS)

end subroutine SUB_
