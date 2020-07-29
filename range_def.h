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


   start_array_loc = mod(localPet*nthreads,ssiLocalDeCount)+1
   end_array_loc   = start_array_loc + nthreads - 1
   print *, __FILE__, __LINE__, localPet, start_array_loc, end_array_loc, nthreads,ssiLocalDeCount

   do array_loc = start_array_loc, end_array_loc
      if (array_loc == start_array_loc) then
         call ESMF_LocalArrayGet(localArrayList(1), farrayPtr=farrayPtr, &
            rc=status)
         VERIFY_(STATUS)
         ndim = size(shape(farrayPtr))
         allocate(arrsize(ndim))
         arrsize = shape(farrayPtr)
         is = 1
         ie = arrsize(1)
         js = 1
         je = arrsize(2)
         if (ndim == 3) km = arrsize(3)
      else
         call ESMF_LocalArrayGet(localArrayList(array_loc), farrayPtr=farrayPtr, &
            rc=status)
         VERIFY_(STATUS)
         ndim = size(shape(farrayPtr))
         allocate(arrsize(ndim))
         arrsize = shape(farrayPtr)
         is = ie + 1
         ie = is + arrsize(1) - 1
         js = je + 1
         je = js + arrsize(2) - 1
         if (ndim == 3) km = arrsize(3)
      end if
      COPY_
      deallocate(arrsize)
   end do
