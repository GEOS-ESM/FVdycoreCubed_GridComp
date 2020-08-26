   call ESMF_AttributeGet(field, name='SSI_ARRAY_SIZE', &
        value=ssiLocalDeCount, rc=status)
   VERIFY_(STATUS)
   allocate(arrayImg(ssiLocalDeCount), stat=status)
   VERIFY_(STATUS)
   call ESMF_AttributeGet(field, name='SSI_ARRAY_SAVED', &
        valueList=arrayImg, rc=status)
   VERIFY_(STATUS)
   array = transfer(arrayimg,array)
   allocate(localDeToDeMap(ssiLocalDeCount), stat=status)
   VERIFY_(STATUS)
   allocate(localArrayList(ssiLocalDeCount), stat=status)
   VERIFY_(STATUS)
   call ESMF_ArrayGet(array, localDeToDeMap=localDeToDeMap, &
        localarrayList=localArrayList, rc=status)
   VERIFY_(STATUS)
   call ESMF_VMGetCurrent(vm, rc=status)
   VERIFY_(STATUS)
   call ESMF_VMGet(vm, localPet=localPet, rc=status)
   VERIFY_(STATUS)
   call ESMF_VMGet(vm, pet=localPet, peCount=nthreads, rc=status)
   VERIFY_(STATUS)

   nth_x = f2c_SSI_arr_map%nth_x
   nth_y = f2c_SSI_arr_map%nth_y
   nnx = f2c_SSI_arr_map%nnx
   nny = f2c_SSI_arr_map%nny
   npet_x = f2c_SSI_arr_map%npet_x
   npet_y = f2c_SSI_arr_map%npet_y
   pet_id_x = f2c_SSI_arr_map%pet_id_x
   pet_id_y = f2c_SSI_arr_map%pet_id_y
   
   !print *, __FILE__, __LINE__, localPet, nthreads, ssiLocalDeCount

   do jth = 1, nth_y
      if (jth == 1) then
         !js = f2c_SSI_arr_map%js
         js = 1
      else
         js = je + 1
      end if
      do ith = 1, nth_x
         if (ith == 1 .and. jth == 1) then
            !is = f2c_SSI_arr_map%is
            is = 1
            ! first fine PET whose DE is the first in coarse will
            ! always reference first local array in localArrayList
            arr_loc = 1
            !print '(9i3,i4)', localPet, ith, jth,pet_id_x,pet_id_y,nth_x,nth_y,nnx,nny, arr_loc
         else
            arr_loc = ith + pet_id_x*nth_x + (pet_id_y*nth_y+jth-1)*nnx
            is = ie + 1
            !print '(9i3,i4)', localPet, ith, jth,pet_id_x,pet_id_y,nth_x,nth_y,nnx,nny, arr_loc
         end if
         call ESMF_LocalArrayGet(localArrayList(arr_Loc), farrayPtr=farrayPtr, &
            rc=status) 
         VERIFY_(STATUS)
         ndim = size(shape(farrayPtr))
         allocate(arrsize(ndim))
         arrsize = shape(farrayPtr)
         ie = is + arrsize(1) - 1
         je = js + arrsize(2) - 1
         if (ndim == 3) km = arrsize(3)
         !print *, __FILE__, localPet, arr_loc, 'is ', is,ie
         !print *, __FILE__, localPet, arr_loc, 'js ', js,je
         call NAME_COPY_(COPY_, rc=status)
         VERIFY_(STATUS)
         !COPY_
       !if(localPet == 0) then
       !   write(113,*) shape(coarse_Array), is, ie, js, je
       !   write(113,*) shape(coarse_Array), arr_loc, ith, jth
       !   write(113,*) coarse_Array RANGE_
       !   write(113,*) '==============='
       !   write(113,*) shape(farrayPtr)
       !   write(113,*) farrayPtr
       !   write(113,*) '==============='
       !   write(113,*) '==============='
       !endif
       !if(localPet == 1) then
       !   write(114,*) shape(coarse_Array), is, ie, js, je
       !   write(114,*) coarse_Array RANGE_
       !   write(114,*) '==============='
       !   write(114,*) shape(farrayPtr)
       !   write(114,*) farrayPtr
       !   write(114,*) '==============='
       !   write(114,*) '==============='
       !endif
       !if(localPet == 2) then
       !   write(115,*) shape(coarse_Array), is, ie, js, je
       !   write(115,*) coarse_Array RANGE_
       !   write(115,*) '==============='
       !   write(115,*) shape(farrayPtr)
       !   write(115,*) farrayPtr
       !   write(115,*) '==============='
       !   write(115,*) '==============='
       !endif
       !if(localPet == 3) then
       !   write(116,*) shape(coarse_Array), is, ie, js, je
       !   write(116,*) coarse_Array RANGE_
       !   write(116,*) '==============='
       !   write(116,*) shape(farrayPtr)
       !   write(116,*) farrayPtr
       !   write(116,*) '==============='
       !   write(116,*) '==============='
       !endif
         !call ESMF_VMBarrier(vm, rc=status)
         !VERIFY_(STATUS)
         deallocate(arrsize)
      end do
   end do
