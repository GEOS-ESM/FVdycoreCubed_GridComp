   type(ESMF_LocalArray), allocatable :: localArrayList(:)
   integer              :: ssiLocalDeCount
   type(ESMF_Field) :: field
   real(TYPEKIND_), pointer :: farrayPtr DIMENSIONS_
   integer :: status
   integer, allocatable :: arrsize(:)
   type(ESMF_VM) :: vm
   integer :: localPet, nthreads
   integer :: arr_loc
   integer :: nth_x, nth_y, nnx, nny, npet_x, npet_y, pet_id_x, pet_id_y
   integer :: ndim, is, ie, js, je, km, ith, jth
   integer, allocatable :: gcImg(:)

   integer :: nx, npx
   character(ESMF_MAXSTR) :: local_name
