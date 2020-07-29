   type(ESMF_LocalArray), allocatable :: localArrayList(:)
   integer, allocatable :: arrayImg(:), localDeToDeMap(:)
   integer              :: ssiLocalDeCount
   type(ESMF_Field) :: field
   type(ESMF_Array) :: array
   real(TYPEKIND_), pointer :: farrayPtr DIMENSIONS_
   integer :: status
   integer :: ndim, is, ie, js, je, km
   integer, allocatable :: arrsize(:)
   type(ESMF_VM) :: vm
   integer :: localPet, nthreads
   integer :: start_array_loc, end_array_loc, array_loc
