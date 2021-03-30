    call timing_on('DATA_COPY')

   call ESMF_FieldGet(field, ssiLocalDeCount=ssiLocalDeCount, rc=status)
   VERIFY_(STATUS)
   allocate(localArrayList(ssiLocalDeCount), stat=status)
   VERIFY_(STATUS)
   call ESMF_FieldGet(field, localarrayList=localArrayList, rc=status)
   VERIFY_(STATUS)

   call ESMF_FieldGet(field,name=local_name, rc=status)
   VERIFY_(STATUS)
   nth_x = f2c_SSI_arr_map%nth_x
   nth_y = f2c_SSI_arr_map%nth_y
   nnx = f2c_SSI_arr_map%nnx
   nny = f2c_SSI_arr_map%nny
   npet_x = f2c_SSI_arr_map%npet_x
   npet_y = f2c_SSI_arr_map%npet_y
   pet_id_x = f2c_SSI_arr_map%pet_id_x
   pet_id_y = f2c_SSI_arr_map%pet_id_y

   npx = f2c_SSI_arr_map%npx
   nx  = f2c_SSI_arr_map%nx

   do jth = 1, nth_y
      if (jth == 1) then
         !js = f2c_SSI_arr_map%js
         js = 1
      else
         js = je + 1
      end if
!$omp parallel do  &
!$omp private(is, ie, je, arr_loc, rc, farrayPtr, status, ndim, arrsize, km) &
!$omp default(shared)
      do ith = 1, nth_x
         if (ith == 1 .and. jth == 1) then
            ! first fine PET whose DE is the first in coarse will
            ! always reference first local array in localArrayList
            arr_loc = 1
         else
            arr_loc = ith + pet_id_x*nth_x + (pet_id_y*nth_y+jth-1)*nnx
         end if
         call ESMF_LocalArrayGet(localArrayList(arr_Loc), farrayPtr=farrayPtr, &
            rc=status) 
         !VERIFY_(STATUS)
         ndim = size(shape(farrayPtr))
         allocate(arrsize(ndim))
         arrsize = shape(farrayPtr)

         is = (npx/nx) * (ith-1) + 1
         ie = (npx/nx) * ith     

         je = js + arrsize(2) - 1
         if (ndim == 3) km = arrsize(3)
         call NAME_COPY_(COPY_, rc=status)
         !VERIFY_(STATUS)
         !COPY_
         deallocate(arrsize)
      end do
!$omp end parallel do
   end do

   deallocate(localArrayList)

   call timing_off('DATA_COPY')

