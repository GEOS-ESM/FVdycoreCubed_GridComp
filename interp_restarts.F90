#define VERIFY_(A)   if((A)/=0) then; PRINT *, 'Interp_restarts.x', __LINE__; call MPI_Abort(A); endif

program interp_restarts

!--------------------------------------------------------------------!
! purpose: driver for interpolation of GEOS FV and Moist restarts    !
!          to the cubed-sphere grid with optional vertical levels    !
!--------------------------------------------------------------------!
   use ESMF
   use mpp_mod,        only: mpp_error, FATAL, NOTE, mpp_root_pe, mpp_broadcast
   use fms_mod,        only: print_memory_usage, fms_init, fms_end, file_exist
   use fv_control_mod, only: fv_init1, fv_init2, fv_end
   use fv_arrays_mod,  only: fv_atmos_type, REAL4, REAL8, FVPRC
   use fv_mp_mod,      only: is_master, ng, mp_gather, tile
   use fv_regrid_c2c, only: get_geos_ic
   use fv_regridding_utils
   use fv_grid_utils_mod, only: ptop_min
   use init_hydro_mod, only: p_var
   use constants_mod,  only: pi, omega, grav, kappa, rdgas, rvgas, cp_air
   use fv_diagnostics_mod,only: prt_mxm
! use fv_eta_mod,     only: set_eta
   use m_set_eta,     only: set_eta
   use memutils_mod, only: print_memuse_stats
   use MAPL
   use pflogger, only: pfl_initialize => initialize
   use gFTL_StringVector
   use gFTL_StringIntegerMap
   use rs_scaleMod

   implicit none

#include "mpif.h"

   type(fv_atmos_type), allocatable, save :: FV_Atm(:)
   logical, allocatable, save             :: grids_on_this_pe(:)

   real, parameter:: zvir = rvgas/rdgas - 1.

   character(ESMF_MAXSTR) :: fname1, str, astr
#ifndef __GFORTRAN__
   external :: getarg, iargc
   integer iargc
#endif
   real(FVPRC) :: dt

   real(ESMF_KIND_R8), allocatable :: r8_ak(:)
   real(ESMF_KIND_R8), allocatable :: r8_bk(:)
   real(ESMF_KIND_R8), allocatable :: r8_akbk(:)

   real(ESMF_KIND_R4), pointer :: r4_local(:,:,:)
   real(ESMF_KIND_R8), pointer :: r8_local(:,:,:), pt_local(:,:,:)
   real(ESMF_KIND_R4), pointer :: r4_local2D(:,:)

   integer i,j,k,n,iq, ihydro
   integer im,jm,km,nq
   real(ESMF_KIND_R8) :: ptop
   real(ESMF_KIND_R8) :: pint

   integer :: is,ie, js,je
   integer :: ks
   integer :: status

   integer :: nmoist
   type(Netcdf4_Fileformatter) :: InFmt,OutFmt
   type(FileMetadata), allocatable :: InCfg(:),OutCfg(:)
   integer :: nVars,imc,jmc
   character(62) :: vname
   type(StringVector) :: moist_variables,all_moist_vars
   type(StringVectorIterator) :: siter
   type(StringIntegerMap) :: tracer_names
   ! bma added
   integer :: p_split, npx, npy, npz, ivar, lcnt_var, iq0
   integer :: n_args,n_files,ifile,nlev,n_output
   character(len=ESMF_MAXPATHLEN), allocatable :: extra_files(:),extra_output(:)
   type(fv_rst), pointer :: rst_files(:) => null()
   type(ArrDescr) :: ArrDes
   integer        :: info
   logical        :: amWriter
   integer :: isl,iel,jsl,jel,npes_x,npes_y,n_writers,n_readers
   type(ESMF_Grid) :: grid
   logical :: scale_rst
   character(len=:), pointer :: var_name
   type(StringVariableMap), pointer :: variables
   type(Variable), pointer :: myVariable
   type(StringVariableMapIterator) :: var_iter
   type(StringVector), pointer :: var_dimensions
   character(len=:), pointer :: dname
   integer :: dim1,ndims
   type(CubedSphereGridFactory) :: csfactory
   real, allocatable :: schmidt_parameters(:)

! Start up FMS/MPP
   print_memory_usage = .true.
   call fms_init()
   call ESMF_Initialize(logKindFlag=ESMF_LOGKIND_NONE,mpiCommunicator=MPI_COMM_WORLD)
   p_split = 1
   call fv_init1(FV_Atm, dt, grids_on_this_pe, p_split)
   call print_memuse_stats('interp_restarts: fms_init')

   n_args = command_argument_count()
   n_files = 0
   n_output = 0
   n_writers=1
   n_readers=1
   ihydro = 1
   scale_rst = .true.
   do i=1,n_args
     call get_command_argument(i,str)
     select case(trim(str))
     case ('-im')
        call get_command_argument(i+1,astr)
        read(astr,*)npx
     case('-lm')
        call get_command_argument(i+1,astr)
        read(astr,*)npz
     case('-do_hydro')
        call get_command_argument(i+1,astr)
        read(astr,*)ihydro
     case('-input_files')
        do j=i+1,n_args
           call get_command_argument(j,astr)
           if ( index(astr,'-') .ne. 0) then
              exit
           end if
           n_files=n_files+1
        enddo
        allocate(extra_files(n_files))
        nq = 0
        do j=i+1,i+n_files
           nq=nq+1
           call get_command_argument(j,extra_files(nq))
        enddo
     case('-output_files')
        do j=i+1,n_args
           call get_command_argument(j,astr)
           if ( index(astr,'-') .ne. 0) then
              exit
           end if
           n_output=n_output+1
        enddo
        allocate(extra_output(n_output))
        nq = 0
        do j=i+1,i+n_output
           nq=nq+1
           call get_command_argument(j,extra_output(nq))
        enddo
     case('-nreader')
        call get_command_argument(i+1,astr)
        read(astr,*)n_readers
     case('-nwriter')
        call get_command_argument(i+1,astr)
        read(astr,*)n_writers
     case('-scalers')
        call get_command_argument(i+1,astr)
        if (trim(astr) == "T") then
           scale_rst=.true.
        else if (trim(astr) == "F") then
           scale_rst=.false.
        else
           write(*,*)'bad argument to scalers, will scale by default'
        end if
     case('-stretched_grid')
        allocate(schmidt_parameters(3))
        call get_command_argument(i+1,astr)
        read(astr,*)schmidt_parameters(1)
        call get_command_argument(i+2,astr)
        read(astr,*)schmidt_parameters(2)
        call get_command_argument(i+3,astr)
        read(astr,*)schmidt_parameters(3)
     end select
   end do


   npx = npx+1
   FV_Atm(1)%flagstruct%npx=npx
   npy = npx
   FV_Atm(1)%flagstruct%npy=npy

   FV_Atm(1)%flagstruct%npz=npz
   FV_Atm(1)%flagstruct%ntiles = 6

   FV_Atm(1)%flagstruct%hydrostatic = .true.
   if (ihydro == 0) FV_Atm(1)%flagstruct%hydrostatic = .false.
   FV_Atm(1)%flagstruct%Make_NH = .false.
   if (.not. FV_Atm(1)%flagstruct%hydrostatic) FV_Atm(1)%flagstruct%Make_NH = .true.
   if (allocated(schmidt_parameters)) then
      FV_Atm(1)%flagstruct%do_schmidt = .true.
      FV_Atm(1)%flagstruct%target_lon=schmidt_parameters(1)
      FV_Atm(1)%flagstruct%target_lat=schmidt_parameters(2)
      FV_Atm(1)%flagstruct%stretch_fac=schmidt_parameters(3)
   end if

   if (n_files > 0) allocate(rst_files(n_files))

! Initialize SHMEM in MAPL
   call pfl_initialize()
   call MAPL_GetNodeInfo (comm=MPI_COMM_WORLD, rc=status)
   call MAPL_InitializeShmem (rc=status)

                       write(fv_atm(1)%flagstruct%grid_file, "('c',i2.2,'_mosaic.nc')") npx-1
   if (npx-1 >=   100) write(fv_atm(1)%flagstruct%grid_file, "('c',i3.3,'_mosaic.nc')") npx-1
   if (npx-1 >=  1000) write(fv_atm(1)%flagstruct%grid_file, "('c',i4.4,'_mosaic.nc')") npx-1
   if (npx-1 >= 10000) write(fv_atm(1)%flagstruct%grid_file, "('c',i5.5,'_mosaic.nc')") npx-1
   dt = 1800
   call fv_init2(FV_Atm, dt, grids_on_this_pe, p_split)

   if (size(extra_files) > 0) then
      if (size(extra_files) /= size(extra_output)) call mpp_error(FATAL, 'the number of extra input and output file names must be same size')
   end if
   call print_memuse_stats('interp_restarts: fv_init')

! Determine Total Number of Tracers (MOIST, GOCART, PCHEM, ANA)
! -------------------------------------------------------------
   nmoist  = 0

   call print_memuse_stats('interp_restarts: rs_count')
   call mpp_broadcast(nmoist, mpp_root_pe())

   if (is_master()) print*, 'HYDROSTATIC : ', FV_Atm(1)%flagstruct%hydrostatic
   if (is_master()) print*, 'Make_NH     : ', FV_Atm(1)%flagstruct%Make_NH
   if (is_master()) print*, 'Tracers     : ', FV_Atm(1)%ncnst

! Need to get ak/bk
   if( file_exist("fvcore_internal_restart_in") ) then
      call InFmt%open("fvcore_internal_restart_in",pFIO_READ,rc=status)
      allocate(InCfg(1))
      InCfg(1) = InFmt%read()
      im = InCfg(1)%get_dimension('lon')
      jm = InCfg(1)%get_dimension('lat')
      km = InCfg(1)%get_dimension('lev')
      call InFmt%close()
      deallocate(InCfg)
   else
      call mpp_error(FATAL, 'ABORT: fvcore_internal_restart_in does not exist')
   endif

   if( file_exist("moist_internal_restart_in") ) then
      call InFmt%open("moist_internal_restart_in",pFIO_READ,rc=status)
      allocate(InCfg(1))
      InCfg(1) = InFmt%read()
      call MAPL_IOCountLevels(InCfg(1),nmoist)
      all_moist_vars = MAPL_IOGetNonDimVars(InCfg(1),rc=status)
      siter = all_moist_vars%begin()
      variables => InCfg(1)%get_variables()
      lcnt_var=2
      do while (siter /= all_moist_vars%end())
         var_name => siter%get()
         myVariable => variables%at(var_name)
         var_dimensions => myVariable%get_dimensions()
         ndims = var_dimensions%size()
         if (ndims==2) nmoist=nmoist-1
         if (ndims==3) then
            call moist_variables%push_back(trim(var_name))
            if (trim(var_name)=='Q') then
               iq0=1
               call tracer_names%insert(trim(var_name),iq0)
            else
               iq0=lcnt_var
               call tracer_names%insert(trim(var_name),iq0)
               lcnt_var=lcnt_var+1
            end if
         end if
         call siter%next()
      end do
      call InFmt%close()
      deallocate(InCfg)
   else
      call mpp_error(FATAL, 'ABORT: moist_internal_restart_in does not exist')
   endif

   call print_memuse_stats('interp_restarts: rs_count')
   call mpp_broadcast(nmoist, mpp_root_pe())

   allocate ( r8_ak(npz+1) )
   allocate ( r8_bk(npz+1) )
   call set_eta(npz,ks,ptop,pint,r8_ak,r8_bk)
   FV_Atm(1)%ak = r8_ak
   FV_Atm(1)%bk = r8_bk
   deallocate ( r8_ak,r8_bk )
   nq = nmoist
   FV_Atm(1)%ncnst = nq/km
   if( is_master() ) then
      print *
      write(6,100)
100      format(2x,' k ','      A(k)    ',2x,' B(k)   ',2x,'  Pref    ',2x,'  DelP',/, &
            1x,'----',3x,'----------',2x,'--------',2x,'----------',2x,'---------' )
      k=1
      write(6,101) k,FV_Atm(1)%ak(k)*0.01, FV_Atm(1)%bk(k), FV_Atm(1)%ak(k)*0.01 + 1000.0*FV_Atm(1)%bk(k)
      do k=2,ubound(FV_Atm(1)%ak,1)
         write(6,102) k,FV_Atm(1)%ak(k)*0.01, FV_Atm(1)%bk(k), FV_Atm(1)%ak(k)*0.01 + 1000.0*FV_Atm(1)%bk(k), &
               (FV_Atm(1)%ak(k)-FV_Atm(1)%ak(k-1))*0.01 + 1000.0*(FV_Atm(1)%bk(k)-FV_Atm(1)%bk(k-1))
      enddo
      print *
101      format(2x,i3,2x,f10.6,2x,f8.4,2x,f10.4)
102      format(2x,i3,2x,f10.6,2x,f8.4,2x,f10.4,3x,f8.4)
103      format(2x,a,i6,3x,a,f7.2,a)
      write(6,103) 'Total Number of Tracers in  MOIST: ',nmoist ,'(/KM = ',float(nmoist) /float(km),')'
      print *
   endif

   do i=1,n_files

      if (file_exist(trim(extra_files(i)))) then
         rst_files(i)%file_name=trim(extra_files(i))

         call InFmt%open(trim(extra_files(i)),pFIO_READ,rc=status)
         allocate(InCfg(1))
         InCfg(1) = InFmt%read()
         call MAPL_IOCountNonDimVars(InCfg(1),nVars)

         if (InCfg(1)%has_dimension("unknown_dim1")) then
            rst_files(i)%ungrid_size = InCfg(1)%get_dimension("unknown_dim1")
         else
            rst_files(i)%ungrid_size = -1
         end if
         rst_files(i)%has_center = InCfg(1)%has_dimension("lev")
         rst_files(i)%has_edge = InCfg(1)%has_dimension("edge")
         allocate(rst_files(i)%vars(nVars))

         variables => InCfg(1)%get_variables()
         var_iter = variables%begin()
         n=0
         do while (var_iter /= variables%end())

            var_name => var_iter%key()
            myVariable => var_iter%value()
            if (.not.InCfg(1)%is_coordinate_variable(var_name)) then
               n=n+1
               var_dimensions => myVariable%get_dimensions()
               ndims = var_dimensions%size()
               rst_files(i)%vars(n)%name=trim(var_name)
               rst_files(i)%vars(n)%n_ungrid=0
               if (ndims ==2) then
                  rst_files(i)%vars(n)%nlev=1
                  rst_files(i)%vars(n)%rank=2
               else if (ndims==3) then
                  rst_files(i)%vars(n)%rank=3
                  dname => myVariable%get_ith_dimension(3)
                  dim1=InCfg(1)%get_dimension(dname)
                  if (dim1 == km) then
                      rst_files(i)%vars(n)%nlev=npz
                  else if (dim1 == km+1) then
                      rst_files(i)%vars(n)%nlev=npz+1
                  else
                      stop
                  end if
               else if (ndims==4) then
                  rst_files(i)%vars(n)%rank=4
                  dname => myVariable%get_ith_dimension(3)
                  dim1=InCfg(1)%get_dimension(dname)
                  if (dim1 == km) then
                      rst_files(i)%vars(n)%nlev=npz
                  else if (dim1 == km+1) then
                      rst_files(i)%vars(n)%nlev=npz+1
                  else
                      stop
                  end if

                  dname => myVariable%get_ith_dimension(4)
                  rst_files(i)%vars(n)%n_ungrid = InCfg(1)%get_dimension(dname)
               end if
            end if
            call var_iter%next()
         enddo

         rst_files(i)%have_descriptor=.true.
         call InFmt%close()
         deallocate(InCfg)
      else
         call mpp_error(FATAL, 'ABORT: '//trim(extra_files(i))//' does not exist')
      end if

   end do

   call print_memuse_stats('interp_restarts: begining get_external_ic')

   npes_x=fv_atm(1)%layout(1)
   npes_y=fv_atm(1)%layout(2)
   is = FV_Atm(1)%bd%isc
   ie = FV_Atm(1)%bd%iec
   js = FV_Atm(1)%bd%jsc
   je = FV_Atm(1)%bd%jec
   isl=is
   iel=ie
   jsl=(npx-1)*(tile-1)+js
   jel=(npx-1)*(tile-1)+je

   call ArrDescrInit(Arrdes,MPI_COMM_WORLD,npx-1,(npx-1)*6,npz,npes_x,npes_y*6,n_readers,n_writers,isl,iel,jsl,jel,rc=status)
   call ArrDescrSet(arrdes,offset=0_MPI_OFFSET_KIND)
   if (allocated(schmidt_parameters)) then
      csfactory = CubedSphereGridFactory(im_world=npx-1,lm=npz,nx=npes_x,ny=npes_y,stretch_factor=schmidt_parameters(3), &
                  target_lon=schmidt_parameters(1),target_lat=schmidt_parameters(2))
   else
      csfactory = CubedSphereGridFactory(im_world=npx-1,lm=npz,nx=npes_x,ny=npes_y)
   end if
   grid = grid_manager%make_grid(csfactory,rc=status)

   FV_Atm(1)%flagstruct%Make_NH = .false. ! Do this after rescaling
   if (jm == 6*im) then
      call get_geos_ic( FV_Atm, rst_files, .true., grid)
   else
      call get_geos_ic( FV_Atm, rst_files, .false., grid)
   endif
   FV_Atm(1)%flagstruct%Make_NH = .true. ! Reset this for later

   if (scale_rst) then
      call scale_drymass(fv_atm,tracer_names,rc=status)
      VERIFY_(status)
   end if

   if (FV_Atm(1)%flagstruct%Make_NH) then
      if (is_master()) print*, 'Updating FV3 NonHydrostatic State'
      do k=1,npz
      do j=js,je
      do i=is,ie
      FV_Atm(1)%w(i,j,k) = 0.0
      FV_Atm(1)%delz(i,j,k) = (-MAPL_RGAS/MAPL_GRAV)*FV_Atm(1)%pt(i,j,k)*(log(FV_Atm(1)%pe(i,k+1,j))-log(FV_Atm(1)%pe(i,k,j)))
      FV_Atm(1)%pkz(i,j,k)  = exp( MAPL_KAPPA*log((-MAPL_RGAS/MAPL_GRAV)*FV_Atm(1)%delp(i,j,k)*FV_Atm(1)%pt(i,j,k)*    &
                                       (1.0+(MAPL_RVAP/MAPL_RGAS - 1.)*FV_Atm(1)%q(i,j,k,1))/FV_Atm(1)%delz(i,j,k)) )
      enddo
      enddo
      enddo
   endif

   allocate(pt_local(is:ie,js:je,npz))
   pt_local=0.0d0
   do k=1,npz
! Convert to Potential Temperature
      pt_local(is:ie,js:je,k) = FV_Atm(1)%pt(is:ie,js:je,k)/FV_Atm(1)%pkz(is:ie,js:je,k)
   enddo

   amWriter = arrdes%writers_comm/=MPI_COMM_NULL
   call MPI_Info_create(info,status)

   call print_memuse_stats('interp_restarts: going to write restarts')

! write fvcore_internal_rst
   if( file_exist("fvcore_internal_restart_in") ) then

      write(fname1, "('fvcore_internal_rst_c',i4.4,'_',i3.3,'L')") npx-1,npz
      if (is_master()) print*, 'Writing : ', TRIM(fname1)

      call InFmt%open("fvcore_internal_restart_in",pFIO_READ,rc=status)
      allocate(InCfg(1),OutCfg(1))
      InCfg(1)=InFmt%read(rc=status)
      call MAPL_IOCountNonDimVars(InCfg(1),nVars)

      if (AmWriter) then
         imc = npx-1
         jmc = imc*6
         call MAPL_IOChangeRes(InCfg(1),OutCfg(1),(/'lon ','lat ','lev ','edge'/),(/imc,jmc,npz,npz+1/),rc=status)
         if (allocated(schmidt_parameters)) then
             call add_stretch_params(OutCfg(1),schmidt_parameters)
         end if

         ! if dz and w were not in the original file add them
         ! they need to be in there for the restart

         if (.not.fv_atm(1)%flagstruct%hydrostatic) then
            ! fix thic
            !call MAPL_NCIOAddVar(ncioOut,"DZ",(/lonid,latid,levid/),6,units="m",long_name="height_thickness",rc=status)
            !call MAPL_NCIOAddVar(ncioOut,"W",(/lonid,latid,levid/),6,units="m s-1",long_name="vertical_velocity",rc=status)
         endif
         call OutFmt%create_par(fname1,comm=arrdes%writers_comm,info=info,rc=status)
         call OutFmt%write(OutCfg(1),rc=status)
      end if


! AK and BK
      allocate ( r8_akbk(npz+1) )
      r8_akbk = FV_Atm(1)%ak
      if (AmWriter) call MAPL_VarWrite(OutFmt,"AK",r8_akbk)
      r8_akbk = FV_Atm(1)%bk
      if (AmWriter) call MAPL_VarWrite(OutFmt,"BK",r8_akbk)
      deallocate ( r8_akbk )

      allocate(r8_local(is:ie,js:je,npz+1))

! U
      if (is_master()) print*, 'Writing : ', TRIM(fname1), ' U'
      call prt_mxm('U', FV_Atm(1)%u, is, ie, js, je, FV_Atm(1)%ng, npz, 1.0_FVPRC, FV_Atm(1)%gridstruct%area_64, FV_Atm(1)%domain)
      r8_local(is:ie,js:je,1:npz) = FV_Atm(1)%u(is:ie,js:je,1:npz)
      call MAPL_VarWrite(OutFmt,"U",r8_local(is:ie,js:je,1:npz),arrdes=arrdes,rc=status)
      VERIFY_(status)
! V
      if (is_master()) print*, 'Writing : ', TRIM(fname1), ' V'
      call prt_mxm('V', FV_Atm(1)%v, is, ie, js, je, FV_Atm(1)%ng, npz, 1.0_FVPRC, FV_Atm(1)%gridstruct%area_64, FV_Atm(1)%domain)
      r8_local(is:ie,js:je,1:npz) = FV_Atm(1)%v(is:ie,js:je,1:npz)
      call MAPL_VarWrite(OutFmt,"V",r8_local(is:ie,js:je,1:npz),arrdes=arrdes,rc=status)
      VERIFY_(status)
! PT
      if (is_master()) print*, 'Writing : ', TRIM(fname1), ' PT'
      call prt_mxm('T', FV_Atm(1)%pt, is, ie, js, je, FV_Atm(1)%ng, npz, 1.0_FVPRC, FV_Atm(1)%gridstruct%area_64, FV_Atm(1)%domain)

      call MAPL_VarWrite(OutFmt,"PT",pt_local(is:ie,js:je,1:npz),arrdes=arrdes,rc=status)
       VERIFY_(status)

! PE
      if (is_master()) print*, 'Writing : ', TRIM(fname1), ' PE'
 !!!  call prt_mxm('PE', FV_Atm(1)%pe, is, ie, js, je, FV_Atm(1)%ng, npz+1, 1.0, FV_Atm(1)%gridstruct%area_64, FV_Atm(1)%domain)
      do k=1,npz+1
         r8_local(is:ie,js:je,k) = FV_Atm(1)%pe(is:ie,k,js:je)
      enddo
      call MAPL_VarWrite(OutFmt,"PE",r8_local(is:ie,js:je,1:npz+1),arrdes=arrdes,rc=status)
      VERIFY_(status)
! PKZ
      if (is_master()) print*, 'Writing : ', TRIM(fname1), ' PKZ'
 !!!  call prt_mxm('PKZ', FV_Atm(1)%pkz, is, ie, js, je, FV_Atm(1)%ng, npz, 1.0, FV_Atm(1)%gridstruct%area_64, FV_Atm(1)%domain)
      r8_local(is:ie,js:je,1:npz) = FV_Atm(1)%pkz(is:ie,js:je,1:npz)
      call MAPL_VarWrite(OutFmt,"PKZ",r8_local(is:ie,js:je,1:npz),arrdes=arrdes,rc=status)
      VERIFY_(status)

      if (.not. fv_atm(1)%flagstruct%hydrostatic) then
! DZ
         if (is_master()) print*, 'Writing : ', TRIM(fname1), ' DZ'
 !!!     call prt_mxm('DZ', FV_Atm(1)%delz, is, ie, js, je, FV_Atm(1)%ng, npz, 1.0, FV_Atm(1)%gridstruct%area_64, FV_Atm(1)%domain)
         r8_local(is:ie,js:je,1:npz) = FV_Atm(1)%delz(is:ie,js:je,1:npz)
         call MAPL_VarWrite(OutFmt,"DZ",r8_local(is:ie,js:je,1:npz),arrdes=arrdes,rc=status)
         VERIFY_(status)

! W
         if (is_master()) print*, 'Writing : ', TRIM(fname1), ' W'
 !!!     call prt_mxm('W', FV_Atm(1)%w, is, ie, js, je, FV_Atm(1)%ng, npz, 1.0, FV_Atm(1)%gridstruct%area_64, FV_Atm(1)%domain)
         r8_local(is:ie,js:je,1:npz) = FV_Atm(1)%w(is:ie,js:je,1:npz)
         call MAPL_VarWrite(OutFmt,"W",r8_local(is:ie,js:je,1:npz),arrdes=arrdes,rc=status)
         VERIFY_(status)
      endif

      if (AmWriter) call OutFmt%close()
      deallocate(InCfg,OutCfg)

      deallocate (r8_local)
      deallocate (pt_local)

   endif

! MOIST
!
      allocate(r4_local(is:ie,js:je,npz+1))
      allocate(r4_local2D(is:ie,js:je))

      if( file_exist("moist_internal_restart_in") ) then
         write(fname1, "('moist_internal_rst_c',i4.4,'_',i3.3,'L')") npx-1,npz
         if (is_master()) print*, 'Writing : ', TRIM(fname1)
         imc = npx-1
         jmc = imc*6
         call InFmt%open("moist_internal_restart_in",pFIO_READ,rc=status)
         allocate(InCfg(1),OutCfg(1))
         InCfg(1)=InFmt%read(rc=status)
         call MAPL_IOChangeRes(InCfg(1),OutCfg(1),(/'lon','lat','lev'/),(/imc,jmc,npz/),rc=status)
         if (allocated(schmidt_parameters)) then
             call add_stretch_params(OutCfg(1),schmidt_parameters)
         end if
         if (AmWriter) then
            call OutFmt%create_par(fname1,comm=arrdes%writers_comm,info=info,rc=status)
            call OutFmt%write(OutCfg(1),rc=status)
         end if
         deallocate(InCfg)
         call InFmt%close()
      end if
      siter = all_moist_Vars%begin()
      Variables => OutCfg(1)%get_variables()
      lcnt_var=2
      ivar=0
      do while (siter /= all_moist_vars%end())
         ivar=ivar+1
         var_name => siter%get()
         myVariable => variables%at(var_name)
         var_dimensions => myVariable%get_dimensions()
         ndims = var_dimensions%size()
         if (is_master()) print*, 'Writing : ', TRIM(fname1), ' ', ivar
         if (ndims==2) then
            r4_local2d(is:ie,js:je)=0.0
            call MAPL_VarWrite(OutFmt,trim(var_name),r4_local2d(is:ie,js:je),arrdes=arrdes,rc=status)
            VERIFY_(status)
         else if (ndims==3) then
            if (trim(var_name)=='Q') then
               iq0=1
            else
               iq0=lcnt_var
               lcnt_var=lcnt_var+1
            end if
            r4_local(is:ie,js:je,1:npz) = FV_Atm(1)%q(is:ie,js:je,:,iq0)
            call MAPL_VarWrite(OutFmt,triM(var_name),r4_local(is:ie,js:je,1:npz),arrdes=arrdes,rc=status)
            VERIFY_(status)
         end if
         call siter%next()
      end do
      if (AmWriter) call OutFmt%close()
      deallocate(outCfg)
      deallocate(r4_local)

! extra restarts
!
      do ifile=1,size(rst_files)

         if (is_master()) write(*,*)'Writing results of ',trim(rst_files(ifile)%file_name)
         fname1=extra_output(ifile)
         if (is_master()) print*, 'Writing : ', TRIM(fname1)
         call ArrDescrSet(arrdes,offset=0_MPI_OFFSET_KIND)
         if (AmWriter) then
            imc = npx-1
            jmc = imc*6
            call InFmt%open(trim(rst_files(ifile)%file_name),pFIO_READ,rc=status)
            allocate(InCfg(1),OutCfg(1))
            InCfg(1)=InFmt%read(rc=status)

            if ((rst_files(ifile)%has_edge .eqv. .false.) &
            .and. (rst_files(ifile)%has_center .eqv. .false.) &
            .and. (rst_files(ifile)%ungrid_size == -1)) then
                call MAPL_IOChangeRes(InCfg(1),OutCfg(1),['lon','lat'],[imc,jmc],rc=status)
            else if ((rst_files(ifile)%has_edge .eqv. .false.) &
            .and. (rst_files(ifile)%has_center .eqv. .true.) &
            .and. (rst_files(ifile)%ungrid_size == -1)) then
               call MAPL_IOChangeRes(InCfg(1),OutCfg(1),['lon','lat','lev'],[imc,jmc,npz],rc=status)
            else if ((rst_files(ifile)%has_edge .eqv. .true.) &
            .and. (rst_files(ifile)%has_center .eqv. .true.) &
            .and. (rst_files(ifile)%ungrid_size == -1)) then
               call MAPL_IOChangeRes(InCfg(1),OutCfg(1),['lon ','lat ','lev ','edge'],[imc,jmc,npz,npz+1],rc=status)
            else if ((rst_files(ifile)%has_edge .eqv. .false.) &
            .and. (rst_files(ifile)%has_center .eqv. .true.) &
            .and. (rst_files(ifile)%ungrid_size > 0)) then
               call MAPL_IOChangeRes(InCfg(1),OutCfg(1),[character(len=12):: 'lon ','lat ','lev ','unknown_dim1'],[imc,jmc,npz,rst_files(ifile)%ungrid_size],rc=status)
            end if
            if (allocated(schmidt_parameters)) then
                call add_stretch_params(OutCfg(1),schmidt_parameters)
            end if

            call OutFmt%create_par(fname1,comm=arrdes%writers_comm,info=info,rc=status)
            call OutFmt%write(OutCfg(1),rc=status)
            deallocate(InCfg,OutCfg)
            call InFmt%close()
         end if
         do iq=1,size(rst_files(ifile)%vars)
            vname = trim(rst_files(ifile)%vars(iq)%name)
            if (is_master()) print*, 'Writing : ', TRIM(vname), ' ', iq
            nlev = rst_files(ifile)%vars(iq)%nlev
            allocate(r4_local(is:ie,js:je,nlev))
            if (rst_files(ifile)%vars(iq)%rank ==2) then
               r4_local2d(is:ie,js:je)=rst_files(ifile)%vars(iq)%ptr2d(is:ie,js:je)
               call MAPL_VarWrite(OutFmt,vname,r4_local2d(is:ie,js:je),arrdes=arrdes)
            else if (rst_files(ifile)%vars(iq)%rank ==3) then
               r4_local(is:ie,js:je,1:nlev)=rst_files(ifile)%vars(iq)%ptr3d(is:ie,js:je,1:nlev)
               call MAPL_VarWrite(OutFmt,vname,r4_local(is:ie,js:je,1:nlev),arrdes)
            else if (rst_files(ifile)%vars(iq)%rank ==4) then
               do n=1,size(rst_files(ifile)%vars(iq)%ptr4d,4)
                  do k=1,size(rst_files(ifile)%vars(iq)%ptr4d,3)
                     r4_local2d(is:ie,js:je)=rst_files(ifile)%vars(iq)%ptr4d(is:ie,js:je,k,n)
                     call MAPL_VarWrite(OutFmt,vname,r4_local2d(is:ie,js:je),arrdes=arrdes,lev=k,offset2=n)
                  enddo
               enddo
            end if
         end do
         deallocate(r4_local)
         if (AmWriter) then
            call OutFmt%close()
         end if

      end do

      deallocate(r4_local2D)

! Finalize SHMEM in MAPL
   call MAPL_FinalizeShmem (rc=status)

   call fv_end(fv_atm, grids_on_this_pe, .false.)
   call fms_end()

contains

   subroutine rs_count( filename,nrecs )
      implicit none
      integer        nrecs,rc
      character(*)   filename

      open  (55,file=trim(filename),form='unformatted',access='sequential')
      nrecs =  0
      rc =  0
      do while (rc.eq.0)
         read (55,iostat=rc)
         if( rc.eq.0 ) nrecs = nrecs + 1
      enddo
      close (55)

      return
   end subroutine rs_count

   subroutine add_stretch_params(meta,stretch_parameters)
      type(FileMetadata), intent(inout) :: meta
      real, intent(in) :: stretch_parameters(3)

      call meta%add_attribute('TARGET_LON',stretch_parameters(1))
      call meta%add_attribute('TARGET_LAT',stretch_parameters(2))
      call meta%add_attribute('STRETCH_FACTOR',stretch_parameters(3))

   end subroutine add_stretch_params

end program interp_restarts
