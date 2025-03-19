#include "MAPL_ErrLog.h"
module fv_regrid_c2c

#ifdef MAPL_MODE
#define DEALLOCGLOB_(A) if(associated(A)) then;A=0;if(MAPL_ShmInitialized) then; call MAPL_DeAllocNodeArray(A,rc=status);else; deallocate(A);endif;NULLIFY(A);endif
#endif

   use fms_mod,            only: file_exist, read_data, field_exist
   use mpp_mod,            only: mpp_error, FATAL
   use mpp_domains_mod,    only: domain2d, mpp_update_domains, mpp_get_boundary, DGRID_NE
   use tracer_manager_mod, only: get_tracer_names, get_number_tracers, get_tracer_index
   use field_manager_mod,  only: MODEL_ATMOS

   use MAPL
   use gFTL_StringVector
   use gFTL_StringIntegerMap
   use, intrinsic :: iso_fortran_env, only: REAL64, REAL32

   use fv_arrays_mod,     only: fv_atmos_type, fv_grid_type, R_GRID, FVPRC
   use fv_diagnostics_mod,only: prt_mxm, prt_maxmin
   use fv_mp_mod,         only: is_master
   use fv_grid_utils_mod, only: ptop_min
   use fv_mapz_mod,       only: map_scalar
   use init_hydro_mod,    only: p_var
   use mpp_mod,           only: mpp_pe
   use memutils_mod, only: print_memuse_stats
   use fv_regridding_utils
   use ESMF

   implicit none

#include "mpif.h"

   private

   real(FVPRC), parameter :: PI           = MAPL_PI_R8
   real(FVPRC), parameter :: OMEGA        = MAPL_OMEGA
   real(FVPRC), parameter :: GRAV         = MAPL_GRAV
   real(FVPRC), parameter :: KAPPA        = MAPL_KAPPA
   real(FVPRC), parameter :: RDGAS        = MAPL_RGAS
   real(FVPRC), parameter :: RVGAS        = MAPL_RVAP
   real(FVPRC), parameter :: CP_AIR       = MAPL_CP

   real(FVPRC), parameter:: zvir = rvgas/rdgas - 1.

   public get_geos_ic

   integer, SAVE :: tile, npes_x, npes_y

   integer :: status
   integer :: IUNIT=15
   integer :: OUNIT=16

   interface read_topo_file
      module procedure read_topo_file_r4
      module procedure read_topo_file_r8
   end interface
                        
contains
                        
   subroutine read_topo_file_r4(fname,output,grid,rc)
      character(len=*), intent(in) :: fname
      type(ESMF_Grid), intent(in) :: grid
      real(REAL32), intent(inout) :: output(:,:)
      integer, intent(out), optional :: rc
      
      integer :: status,dims(3),funit
      integer :: rank
      type(ESMF_VM) :: vm
      real(REAL32), allocatable :: input(:,:)

      call ESMF_VMGetCurrent(vm,_RC)
      call ESMF_VMGet(vm,localPet=rank,_RC)
      call MAPL_GridGet(grid,globalCellCountPerDim=dims,_RC)
      if (rank ==0) then
         allocate(input(dims(1),dims(2)))
         open(newunit=funit,file=trim(fname),form='unformatted',iostat=status)
         _VERIFY(status)
         read(funit)input
      else
         allocate(input(0,0))
      end if
      call ArrayScatter(local_array=output,global_array=input,grid=grid,_RC)
      _RETURN(_SUCCESS)
   end subroutine read_topo_file_r4

   subroutine read_topo_file_r8(fname,output,grid,rc)
      character(len=*), intent(in) :: fname
      type(ESMF_Grid), intent(in) :: grid
      real(REAL64), intent(inout) :: output(:,:)
      integer, intent(out), optional :: rc
      
      integer :: status,dims(3),funit
      real, allocatable :: input(:,:)
      integer :: rank
      type(ESMF_VM) :: vm
      real(REAL64), allocatable :: input_r8(:,:)

      call ESMF_VMGetCurrent(vm,_RC)
      call ESMF_VMGet(vm,localPet=rank,_RC)
      call MAPL_GridGet(grid,globalCellCountPerDim=dims,_RC)
      if (rank ==0) then
         allocate(input(dims(1),dims(2)),input_r8(dims(2),dims(2)))
         open(newunit=funit,file=trim(fname),form='unformatted',iostat=status)
         _VERIFY(status)
         read(funit)input
         input_r8 = input      
      else
         allocate(input(0,0),input_r8(0,0))
      end if
      call ArrayScatter(local_array=output,global_array=input_r8,grid=grid,_RC)
      _RETURN(_SUCCESS)
   end subroutine read_topo_file_r8

   subroutine get_geos_ic( Atm_i, Atm, grid_i, grid, Arrdes_i, extra_rst )

      type(fv_atmos_type), intent(inout) :: Atm_i(:), Atm(:)
      type(ESMF_Grid), intent(inout) :: grid_i, grid
      type(ArrDescr), intent(inout) :: ArrDes_i
      type(fv_rst), pointer, intent(inout) :: extra_rst(:)
      integer i,j

      do i=1,size(extra_rst)
         do j=1,size(extra_rst(i)%vars)
            call extra_rst(i)%vars(j)%alloc_var(Atm(1)%bd%isd,Atm(1)%bd%ied,Atm(1)%bd%jsd,Atm(1)%bd%jed)
         enddo
      enddo

      call get_geos_cubed_ic( Atm_i, Atm, grid_i, grid, Arrdes_i, extra_rst)

      Atm(1)%flagstruct%dry_mass = MAPL_PSDRY
      Atm(1)%flagstruct%adjust_dry_mass = .true.

     ! Atm(1)%pt must be Dry T
     ! nh_pkz=.false. since FV3 restarts right now expect the hydrostatic pkz ALWAYS !
      call p_var(Atm(1)%npz, Atm(1)%bd%is, Atm(1)%bd%ie, Atm(1)%bd%js, Atm(1)%bd%je, Atm(1)%ak(1),  ptop_min,    &
            Atm(1)%delp, Atm(1)%delz, Atm(1)%pt, Atm(1)%ps,               &
            Atm(1)%pe,   Atm(1)%peln, Atm(1)%pk, Atm(1)%pkz,              &
            kappa, Atm(1)%q, Atm(1)%ng, Atm(1)%ncnst, dble(Atm(1)%gridstruct%area),Atm(1)%flagstruct%dry_mass,    &
            Atm(1)%flagstruct%adjust_dry_mass, Atm(1)%flagstruct%mountain, Atm(1)%flagstruct%moist_phys,   &
            Atm(1)%flagstruct%hydrostatic, Atm(1)%flagstruct%nwat, Atm(1)%domain, Atm(1)%flagstruct%make_nh, nh_pkz=.false.)
       if (Atm(1)%flagstruct%make_nh) Atm(1)%w = tiny(1.0_FVPRC)

   end subroutine get_geos_ic

   subroutine get_geos_cubed_ic( Atm_i, Atm, grid_i, grid, Arrdes_i, extra_rst )
      type(fv_atmos_type), intent(inout) :: Atm_i(:), Atm(:)
      type(ESMF_Grid), intent(inout) :: grid_i, grid
      type(ArrDescr), intent(inout) :: ArrDes_i
      type(fv_rst), pointer, intent(inout) :: extra_rst(:)

      character(len=128) :: fname, fname1
      real(FVPRC), allocatable:: pkz0(:,:)
      real(FVPRC), allocatable:: ps0(:,:), gz0(:,:), t0(:,:,:), q0(:,:,:), qlev(:,:)
      real(FVPRC), allocatable:: pe0(:,:,:), u0(:,:,:), v0(:,:,:), w0(:,:,:), dz0(:,:,:)
      real(FVPRC), allocatable:: ak0(:), bk0(:)
      integer :: i, j, k, l, iq, im, jm, km, npx, npy, npz
      integer :: ntiles=6
      character (len=8) :: imc, jmc

      real(FVPRC), allocatable:: psc(:,:), gzc(:,:)
      real(FVPRC), allocatable:: tp(:,:,:), qp(:,:,:,:)
      real(FVPRC), allocatable:: ud(:,:,:), vd(:,:,:)
      real(FVPRC), allocatable:: wp(:,:,:), dzp(:,:,:), dpp(:,:,:)

      real(REAL64), allocatable :: akbk_r8(:)

      integer :: is_i,ie_i, js_i,je_i
      integer :: isd_i,ied_i, jsd_i,jed_i

      integer :: is,ie, js,je
      integer :: isd,ied, jsd,jed
      integer :: ng

      real(FVPRC), allocatable :: ebuffer(:,:)
      real(FVPRC), allocatable :: nbuffer(:,:)
      real(FVPRC), allocatable :: wbuffer(:,:)
      real(FVPRC), allocatable :: sbuffer(:,:)

      integer            :: filetype
      type(Netcdf4_Fileformatter) :: formatter
      type(FileMetadata), allocatable :: cfg(:)
      integer            :: nDims, nVars, ivar, n_ungrid
      character(len=128) :: vname
      integer            :: lvar_cnt,ifile,nlev
      type(fv_rst), pointer :: tracer_bundles(:) => null()

!bma added
      character(len=:), pointer :: var_name
      type(StringVariableMap), pointer :: variables
      type(Variable), pointer :: myVariable
      type(StringVector) :: all_moist_vars 
      type(StringVector), pointer :: var_dimensions
      type(StringVectorIterator) :: siter
      type(StringVector) :: moist_variables
      type(StringIntegerMap) :: moist_tracers
      integer, pointer :: iptr
      character(len=:), pointer :: cptr
      type(StringIntegerMapIterator) :: iter
      type(ESMF_Grid) :: gridIn
      class(AbstractRegridder), pointer :: regridder=>null()

!--------------------------------------------------------------------!
! create ESMF regridder
!--------------------------------------------------------------------!
      regridder => new_regridder_manager%make_regridder(grid_i,grid,REGRID_METHOD_BILINEAR,rc=status)

!--------------------------------------------------------------------!
! initialize cubed sphere grid: in                                   !
!--------------------------------------------------------------------!
      isd_i = Atm_i(1)%bd%isd
      ied_i = Atm_i(1)%bd%ied
      jsd_i = Atm_i(1)%bd%jsd
      jed_i = Atm_i(1)%bd%jed      
      is_i  = Atm_i(1)%bd%is
      ie_i  = Atm_i(1)%bd%ie
      js_i  = Atm_i(1)%bd%js
      je_i  = Atm_i(1)%bd%je

!--------------------------------------------------------------------!
! initialize cubed sphere grid: out                                  !
!--------------------------------------------------------------------!
      npx = Atm(1)%npx
      npy = Atm(1)%npy
      npz = Atm(1)%npz
      isd = Atm(1)%bd%isd
      ied = Atm(1)%bd%ied
      jsd = Atm(1)%bd%jsd
      jed = Atm(1)%bd%jed
      is  = Atm(1)%bd%is
      ie  = Atm(1)%bd%ie 
      js  = Atm(1)%bd%js 
      je  = Atm(1)%bd%je 
      ng  = Atm(1)%ng

! Zero out all initial tracer fields:
      if (allocated(Atm(1)%q)) deallocate( Atm(1)%q )
      allocate  ( Atm(1)%q(Atm(1)%bd%isd:Atm(1)%bd%ied,Atm(1)%bd%jsd:Atm(1)%bd%jed,Atm(1)%npz,Atm(1)%ncnst) )
      Atm(1)%q = 0.
! Read input FV core restart file
      fname = "fvcore_internal_restart_in"
      if( file_exist(fname) ) then

         allocate(cfg(1))
         call formatter%open(fname,pFIO_READ,rc=status)
         cfg(1) = formatter%read(rc=status)
         im =cfg(1)%get_dimension('lon',rc=status)
         jm =cfg(1)%get_dimension('lat',rc=status)
         km =cfg(1)%get_dimension('lev',rc=status)

         if(is_master()) write(*,*) 'Using GEOS restart:', fname
         if(is_master()) write(*,*) 'External IC dimensions:', im   , jm       , km
         if(is_master()) write(*,*) 'Interpolating to      :', npx-1, (npy-1)*6, npz

         allocate ( ak0(km+1) )
         allocate ( bk0(km+1) )
         allocate ( akbk_r8(km+1) )
         call MAPL_VarRead(formatter,"AK",akbk_r8)
         ak0 = akbk_r8
         call MAPL_VarRead(formatter,"BK",akbk_r8)
         bk0 = akbk_r8
         deallocate ( akbk_r8 )
         call print_memuse_stats('get_geos_cubed_ic: read ak/bk')
         close (IUNIT)

         if( is_master() ) then
           print *
           write(*,*) 'Input Vertical Grid'     
           write(*,*) '--------------------'
           write(6,100)
100           format(2x,' k ','      A(k)    ',2x,' B(k)   ',2x,'  Pref    ',2x,'  DelP',/, &
                 1x,'----',3x,'----------',2x,'--------',2x,'----------',2x,'---------' )
           k=1
           write(6,101) k,ak0(k)*0.01, bk0(k), ak0(k)*0.01 + 1000.0*bk0(k)
           do k=2,ubound(ak0,1)
              write(6,102) k,ak0(k)*0.01, bk0(k), ak0(k)*0.01 + 1000.0*bk0(k), &
                    (ak0(k)-ak0(k-1))*0.01 + 1000.0*(bk0(k)-bk0(k-1))
           enddo
           print *
101        format(2x,i3,2x,f10.6,2x,f8.4,2x,f10.4)
102        format(2x,i3,2x,f10.6,2x,f8.4,2x,f10.4,3x,f8.4)
         endif

! Read U
         allocate (  u0(isd_i:ied_i,jsd_i:jed_i+1,km) )
         u0(:,:,:) = 0.0 
         do k=1,km
            call MAPL_VarRead(formatter,"U",u0(is_i:ie_i,js_i:je_i,k),arrdes=Arrdes_i,lev=k)
         enddo
         call print_memuse_stats('get_geos_cubed_ic: read U')
! Read V
         allocate (  v0(isd_i:ied_i+1,jsd_i:jed_i,km) )
         v0(:,:,:) = 0.0
         do k=1,km
            call MAPL_VarRead(formatter,"V",v0(is_i:ie_i,js_i:je_i,k),arrdes=Arrdes_i,lev=k)
         enddo
         call print_memuse_stats('get_geos_cubed_ic: read V')
         allocate ( sbuffer(is_i:ie_i,km) )
         allocate ( wbuffer(js_i:je_i,km) )
         allocate ( nbuffer(is_i:ie_i,km) )
         allocate ( ebuffer(js_i:je_i,km) )
         call mpp_get_boundary(u0, v0, Atm_i(1)%domain, &
               wbuffery=wbuffer, ebuffery=ebuffer, &
               sbufferx=sbuffer, nbufferx=nbuffer, &
               gridtype=DGRID_NE )
         do k=1,km
            do i=is_i,ie_i    
               u0(i,je_i+1,k) = nbuffer(i,k)
            enddo
            do j=js_i,je_i
               v0(ie_i+1,j,k) = ebuffer(j,k)
            enddo
         enddo
         deallocate ( sbuffer )
         deallocate ( wbuffer )
         deallocate ( nbuffer )
         deallocate ( ebuffer )
         call mpp_update_domains( u0, v0, Atm_i(1)%domain, gridtype=DGRID_NE, complete=.true. )
         call prt_maxmin(' U_geos', u0, is_i, ie_i  , js_i, je_i+1, Atm_i(1)%ng, km, 1.0_FVPRC)
         call prt_maxmin(' V_geos', v0, is_i, ie_i+1, js_i, je_i  , Atm_i(1)%ng, km, 1.0_FVPRC)
         allocate ( ud(isd:ied  ,jsd:jed+1,km) )
         allocate ( vd(isd:ied+1,jsd:jed  ,km) )
!------------------------------------------------------------------!
! D->A : regrid : A-> D interpolation for U and V components
!------------------------------------------------------------------!
         do k=1,km
            call d2a2d(u0(:,:,k), v0(:,:,k), ud(:,:,k), vd(:,:,k), &
                       Atm_i(1), Atm(1), regridder)
         enddo ! npz
! ------------------------------------------------------------------------------------------------------------------------------------------!
         allocate ( sbuffer(is:ie,km) )
         allocate ( wbuffer(js:je,km) )
         allocate ( nbuffer(is:ie,km) )
         allocate ( ebuffer(js:je,km) )
         call mpp_get_boundary(ud, vd, Atm(1)%domain, &
               wbuffery=wbuffer, ebuffery=ebuffer, &
               sbufferx=sbuffer, nbufferx=nbuffer, &
               gridtype=DGRID_NE )
         do k=1,km
            do i=is,ie    
               ud(i,je+1,k) = nbuffer(i,k)
            enddo 
            do j=js,je
               vd(ie+1,j,k) = ebuffer(j,k)
            enddo
         enddo
         deallocate ( sbuffer )
         deallocate ( wbuffer )
         deallocate ( nbuffer )
         deallocate ( ebuffer )
         call mpp_update_domains( ud, vd, Atm(1)%domain, gridtype=DGRID_NE, complete=.true. )
         call prt_maxmin(' U_rgrd', ud, is, ie  , js, je+1, Atm(1)%ng, km, 1.0_FVPRC)
         call prt_maxmin(' V_rgrd', vd, is, ie+1, js, je  , Atm(1)%ng, km, 1.0_FVPRC)
      !  ! ud and vd and new winds on D-Grid: Cubed-Sphere oriented
         deallocate ( v0 )
         deallocate ( u0 )

! Read W
         if (.not. Atm(1)%flagstruct%hydrostatic) then
         allocate (  w0(is_i:ie_i,js_i:je_i,km) )
         allocate (  wp(is  :ie  ,js  :je  ,km) )
         do k=1,km
            call MAPL_VarRead(formatter,"W",w0(is_i:ie_i,js_i:je_i,k),arrdes=Arrdes_i,lev=k)
            call regridder%regrid(w0(is_i:ie_i,js_i:je_i,k),wp(:,:,k),rc=status)
         enddo
         call prt_maxmin(' W_geos', w0, is_i, ie_i, js_i, je_i, 0, km, 1.0_FVPRC)
         call prt_maxmin(' W_rgrd', wp, is, ie, js, je, 0, km, 1.0_FVPRC)
         deallocate ( w0 )
         call print_memuse_stats('get_geos_cubed_ic: read W')
         endif
! Read DZ
         if (.not. Atm(1)%flagstruct%hydrostatic) then
         allocate (  dz0(is_i:ie_i,js_i:je_i,km) )
         allocate (  dzp(is  :ie  ,js  :je  ,km) )
         do k=1,km
            call MAPL_VarRead(formatter,"DZ",dz0(is_i:ie_i,js_i:je_i,k),arrdes=Arrdes_i,lev=k)
            call regridder%regrid(dz0(is_i:ie_i,js_i:je_i,k),dzp(:,:,k),rc=status)
         enddo
         call prt_maxmin('DZ_geos', dz0, is_i, ie_i, js_i, je_i, 0, km, 1.0_FVPRC)
         call prt_maxmin('DZ_rgrd', dzp, is, ie, js, je, 0, km, 1.0_FVPRC)
         deallocate ( dz0 )
         call print_memuse_stats('get_geos_cubed_ic: read DZ')
         endif
! Read PT
         allocate (  t0(is_i:ie_i,js_i:je_i,km) )
         do k=1,km
            call MAPL_VarRead(formatter,"PT",t0(is_i:ie_i,js_i:je_i,k),arrdes=Arrdes_i,lev=k)
         enddo
         call prt_maxmin('PT_geos', t0, is_i, ie_i, js_i, je_i, 0, km, 1.0_FVPRC)
         call print_memuse_stats('get_geos_cubed_ic: read T')
! Read PE 
         allocate ( pe0(is_i:ie_i,js_i:je_i,km+1) )
         do k=1,km+1
           call MAPL_VarRead(formatter,"PE",pe0(is_i:ie_i,js_i:je_i,k),arrdes=Arrdes_i,lev=k)
         enddo
! Get PS
         allocate ( ps0(is_i:ie_i,js_i:je_i) )
         ps0(:,:) = pe0(:,:,km+1)
! Read PKZ
         allocate ( pkz0(is_i:ie_i,js_i:je_i) )
         do k=1,km
            call MAPL_VarRead(formatter,"PKZ",pkz0(is_i:ie_i,js_i:je_i),arrdes=Arrdes_i,lev=k)
            t0(:,:,k) = t0(:,:,k)*pkz0(:,:)
         enddo
         call prt_maxmin( ' T_geos', t0, is_i, ie_i, js_i, je_i, 0, km, 1.0_FVPRC)
         deallocate ( pkz0 )
         call print_memuse_stats('get_geos_cubed_ic: converted T')

! Horiz Interp for T
         allocate (  tp(is:ie,js:je,km) )
         do k=1,km
            call regridder%regrid(t0(is_i:ie_i,js_i:je_i,k),tp(:,:,k),rc=status)
         enddo
         call prt_maxmin( ' T_rgrd', tp, is, ie, js, je, 0, km, 1.0_FVPRC)
         deallocate ( t0 )
         call print_memuse_stats('get_geos_cubed_ic: converted T')

         call formatter%close()
         deallocate(cfg)

         write(imc, "(i8)") im
         write(jmc, "(i8)") jm
         imc = adjustl(imc)
         jmc = adjustl(jmc)

! Read input topography
         write(fname1, "('topo_DYN_ave_',a,'x',a,'.data')") trim(imc), trim(jmc)
         if (.not. file_exist(fname1)) then
            call mpp_error(FATAL,'get_geos_cubed_ic: cannot find topo_DYN_ave file')
         endif
         allocate ( gz0(is_i:ie_i,js_i:je_i) )
         call print_memuse_stats('get_geos_cubed_ic: '//TRIM(fname1)//' being read')
         call read_topo_file(fname1,gz0(is_i:ie_i,js_i:je_i),grid_i)
         gz0 = gz0*grav

! Horiz Interp for surface pressure 
         allocate( psc(is:ie,js:je) )
         call prt_maxmin('PS_geos', ps0, is_i, ie_i, js_i, je_i, 0, 1, 1.0_FVPRC)
         call regridder%regrid(ps0(is_i:ie_i,js_i:je_i),psc(is:ie,js:je),rc=status)
         deallocate ( ps0 )
! Horiz Interp for surface height
         allocate( gzc(is:ie,js:je) )
         call prt_maxmin('GZ_geos', gz0, is_i, ie_i, js_i, je_i, 0, 1, 1.0/grav)
         call regridder%regrid(gz0(is_i:ie_i,js_i:je_i),gzc(is:ie,js:je),rc=status)
         deallocate ( gz0 )

         call prt_maxmin('PS_rgrd', psc, is, ie, js, je, 0, 1, 1.0_FVPRC)
         call prt_maxmin('GZ_rgrd', gzc, is, ie, js, je, 0, 1, 1.0/grav)

! Horiz Interp for Q
         allocate ( q0(is_i:ie_i,js_i:je_i,km+1) )
         allocate ( qp(is:ie,js:je,km,Atm(1)%ncnst) )
         q0(:,:,:) = 0.0
         qp(:,:,:,:) = 0.0

! Horiz Interp for moist tracers
! is there a moist restart file to interpolate?
! Read in tracers: only sphum at this point
         if( file_exist("moist_internal_restart_in") ) then
            if (is_master()) print*, ''
            if (is_master()) print*, 'Regridding moist_internal_restart_in'

            call MAPL_NCIOGetFileType("moist_internal_restart_in",filetype)

            lvar_cnt = 0
            allocate(cfg(1))
            call formatter%open("moist_internal_restart_in",pFIO_READ,rc=status)
            cfg(1) = formatter%read(rc=status)
            all_moist_vars = MAPL_IOGetNonDimVars(cfg(1))
            siter = all_moist_vars%begin()
            Variables => cfg(1)%get_variables()
            do while(siter /= all_moist_vars%end())
               var_name => siter%get()
               myVariable => variables%at(var_name)
               var_dimensions => myVariable%get_dimensions()
               ndims = var_dimensions%size()
               if (ndims == 3) call moist_variables%push_back(trim(var_name))
               call siter%next()
            enddo
            if (moist_variables%size() /= atm(1)%ncnst) call mpp_error(FATAL,'Wrong number of variables in moist file') 

            lvar_cnt=0
            do ivar=1,Atm(1)%ncnst
               vname = moist_variables%at(ivar)
               if (trim(vname)=='Q') then
                  lvar_cnt=lvar_cnt+1
               elseif (trim(vname)=='QLLS') then
                  lvar_cnt=lvar_cnt+1
               elseif (trim(vname)=='QLCN') then
                  lvar_cnt=lvar_cnt+1
               elseif (trim(vname)=='QILS') then
                  lvar_cnt=lvar_cnt+1
               elseif (trim(vname)=='QICN') then
                  lvar_cnt=lvar_cnt+1
               elseif (trim(vname)=='QRAIN') then
                  lvar_cnt=lvar_cnt+1
               elseif (trim(vname)=='QSNOW') then
                  lvar_cnt=lvar_cnt+1
               elseif (trim(vname)=='QGRAUPEL') then
                  lvar_cnt=lvar_cnt+1
               endif
            enddo
            Atm(1)%flagstruct%nwat = lvar_cnt
            lvar_cnt=lvar_cnt+1

            do ivar=1,Atm(1)%ncnst
               vname = moist_variables%at(ivar)
               if (trim(vname)=='Q') then
                  iq=1
               elseif (trim(vname)=='QLLS') then
                  iq=2
               elseif (trim(vname)=='QLCN') then
                  iq=3
               elseif (trim(vname)=='QILS') then
                  iq=4
               elseif (trim(vname)=='QICN') then
                  iq=5
               elseif (trim(vname)=='QRAIN') then
                  iq=6
               elseif (trim(vname)=='QSNOW') then
                  iq=7
               elseif (trim(vname)=='QGRAUPEL') then
                  iq=8
               else
                  iq=lvar_cnt
                  lvar_cnt=lvar_cnt+1
               end if
               call moist_tracers%insert(trim(vname),iq)
               do k=1,km
                  call MAPL_VarRead(formatter,vname,q0(is_i:ie_i,js_i:je_i,k),arrdes=Arrdes_i,lev=k)
                  call regridder%regrid(q0(is_i:ie_i,js_i:je_i,k),qp(:,:,k,iq),rc=status)
               enddo
               call prt_maxmin( trim(vname)//'_geos_moist', q0, is_i, ie_i, js_i, je_i, 0, km, 1._FVPRC)
            enddo

            call formatter%close()
            deallocate(cfg)

         end if

! Horiz Interp for extra tracers
! make copy of input on input levs

        call copy_fv_rst(extra_rst,tracer_bundles)
        do i=1,size(extra_rst)
           do j=1,size(extra_rst(i)%vars)
              if (extra_rst(i)%vars(j)%nLev/=1) then
                 if (extra_rst(i)%vars(j)%nLev == npz) then 
                    tracer_bundles(i)%vars(j)%nLev=km
                    call tracer_bundles(i)%vars(j)%alloc_var(is,ie,js,je,km)
                 else if (extra_rst(i)%vars(j)%nLev == npz+1) then
                    tracer_bundles(i)%vars(j)%nLev=km+1
                    call tracer_bundles(i)%vars(j)%alloc_var(is,ie,js,je,km+1)
                 end if    
              else
                 call tracer_bundles(i)%vars(j)%alloc_Var(is,ie,js,je) 
              end if
           enddo
        enddo

        do ifile=1,size(tracer_bundles)
            if (is_master()) print*, ''
            if (is_master()) print*, 'Regridding: ',trim(tracer_bundles(ifile)%file_name)

            allocate(cfg(1))
            call formatter%open(trim(tracer_bundles(ifile)%file_name),pFIO_READ,rc=status)
            cfg(1) = formatter%read(rc=status)
            call MAPL_IOCountNonDimVars(cfg(1),nvars,rc=status)

            allocate ( qlev(is_i:ie_i,js_i:je_i) )
            qlev(:,:) = 0.0

            do ivar=1,size(tracer_bundles(ifile)%vars)
               nlev=tracer_bundles(ifile)%vars(ivar)%nLev
               vname = trim(tracer_bundles(ifile)%vars(ivar)%name)
               if (tracer_bundles(ifile)%vars(ivar)%rank ==2) then
                  call MAPL_VarRead(formatter,vname,qlev(is_i:ie_i,js_i:je_i),arrdes=Arrdes_i)
                  q0(is_i:ie_i,js_i:je_i,1) = qlev(is_i:ie_i,js_i:je_i)
                  call regridder%regrid(qlev(is_i:ie_i,js_i:je_i),tracer_bundles(ifile)%vars(ivar)%ptr2d(is:ie,js:je),rc=status)
                  call prt_maxmin( trim(vname)//'_geos_'//trim(tracer_bundles(ifile)%file_name), q0, is_i, ie_i, js_i, je_i, 0, 1, 1._FVPRC)
               else if (tracer_bundles(ifile)%vars(ivar)%rank ==3) then
                  do k=1,nlev
                     call MAPL_VarRead(formatter,vname,qlev(is_i:ie_i,js_i:je_i),arrdes=Arrdes_i,lev=k)
                     q0(is_i:ie_i,js_i:je_i,k) = qlev(is_i:ie_i,js_i:je_i)
                     call regridder%regrid(qlev(is_i:ie_i,js_i:je_i),tracer_bundles(ifile)%vars(ivar)%ptr3d(is:ie,js:je,k),rc=status)
                  enddo
                  call prt_maxmin( trim(vname)//'_geos_'//trim(tracer_bundles(ifile)%file_name), q0, is_i, ie_i, js_i, je_i, 0, nlev, 1._FVPRC)
               else if (tracer_bundles(ifile)%vars(ivar)%rank ==4) then
                  do n_ungrid=1,tracer_bundles(ifile)%vars(ivar)%n_ungrid
                     do k=1,nlev
                        call MAPL_VarRead(formatter,vname,qlev(is_i:ie_i,js_i:je_i),arrdes=Arrdes_i,lev=k,offset2=n_ungrid)
                        q0(is_i:ie_i,js_i:je_i,k) = qlev(is_i:ie_i,js_i:je_i)       
                        call regridder%regrid(qlev(is_i:ie_i,js_i:je_i),tracer_bundles(ifile)%vars(ivar)%ptr4d(is:ie,js:je,k,n_ungrid),rc=status)
                     enddo
                     call prt_maxmin( trim(vname)//'_geos_'//trim(tracer_bundles(ifile)%file_name), q0, is_i, ie_i, js_i, je_i, 0, nlev, 1._FVPRC)
                  enddo
               end if
            enddo

            call formatter%close()
            deallocate(cfg)
            deallocate(qlev)

         enddo
         deallocate ( q0 )
                   
         if (is_master()) print*, ''
         if (is_master()) print*, 'Vertical Remapping: '
! Vert remap for scalars
         call read_topo_file('topo_dynave.data',Atm(1)%phis(is:ie,js:je),grid)
         call mpp_update_domains(Atm(1)%phis, Atm(1)%domain)
         Atm(1)%phis = Atm(1)%phis*grav
         call remap_scalar(im, jm, km, npz, Atm(1)%ncnst, Atm(1)%ncnst, ak0, bk0, psc, gzc, &
                           tp, wp, dzp, qp, Atm(1), tracer_bundles, extra_rst )

         if (.not. Atm(1)%flagstruct%hydrostatic) then
            deallocate ( wp )
            deallocate ( dzp )
         endif
         deallocate ( tp )
         deallocate ( qp )
         call print_memuse_stats('get_geos_cubed_ic: remap_scalar')
! Horz/Vert remap for U/V
         call remap_winds(is,ie, js, je, isd,ied, jsd,jed, km, npz, ak0, bk0, psc, ud, vd, Atm(1))
         deallocate ( ud )
         deallocate ( vd )
         call print_memuse_stats('get_geos_cubed_ic: remap_winds')

      else
         call mpp_error(FATAL,'==> Error from get_geos_ic:        &
               & Expected file '//trim(fname)//' does not exist')
      endif

      deallocate( psc, gzc )

      if (allocated(bk0)) deallocate ( bk0 )
      if (allocated(ak0)) deallocate ( ak0 )

      call prt_maxmin('GZ_model', Atm(1)%phis, is, ie, js, je, ng, 1, 1.0/grav)
      call prt_maxmin('PS_model', Atm(1)%ps, is, ie, js, je, ng, 1, 0.01_FVPRC)
      call prt_maxmin('DP_model', Atm(1)%delp, is, ie, js, je, ng, npz, 1.0_FVPRC)
      call prt_maxmin(' U_model', Atm(1)%u, is, ie, js, je, ng, npz, 1.0_FVPRC)
      call prt_maxmin(' V_model', Atm(1)%v, is, ie, js, je, ng, npz, 1.0_FVPRC)
      call prt_maxmin(' T_model', Atm(1)%pt, is, ie, js, je, ng, npz, 1.0_FVPRC)
      if (.not. Atm(1)%flagstruct%hydrostatic) then
        call prt_maxmin(' W_model', Atm(1)%w, is, ie, js, je, ng, npz, 1.0_FVPRC)
        call prt_maxmin('DZ_model', Atm(1)%delz, is, ie, js, je, ng, npz, 1.0_FVPRC)
      endif

! Range check the MOIST tracers
! Iterate over tracer names
  
      iter = moist_tracers%begin()
      do while (iter /= moist_tracers%end())
         iptr => iter%value()
         cptr => iter%key()
         if (.not.match(cptr)) then 
            do k=1,npz
               do j=js,je
                  do i=is,ie
                     Atm(1)%q(i,j,k,iptr) = MIN(Atm(1)%q(i,j,k,iptr),1.d0)
                     Atm(1)%q(i,j,k,iptr) = MAX(Atm(1)%q(i,j,k,iptr),0.d0)
                  enddo
               enddo
            enddo
         endif
         call iter%next()
      enddo
      do iq=1,Atm(1)%ncnst
         call prt_maxmin('QP_model', Atm(1)%q(is:ie,js:je,1:npz,iq), is, ie, js, je, 0, npz, 1._FVPRC)
      enddo

      do i=1,size(tracer_bundles)
         do j=1,size(tracer_bundles(i)%vars)
            call tracer_bundles(i)%vars(j)%dealloc_var()
         enddo
      enddo
      deallocate(tracer_bundles)

      contains
         function match(var_name) result(inList)
            character(len=*) :: var_name
            logical :: inList
            integer :: i
            character(len=10) :: exclude_vars(3)
            exclude_vars(1)="Q"
            exclude_vars(2)="NCPL"
            exclude_vars(3)="NCPI"
            inList = .false.
            do i=1,size(exclude_vars)
               if (trim(exclude_vars(i))==trim(var_name)) inList = .true.
            enddo
         end function match

   end subroutine get_geos_cubed_ic

            subroutine remap_winds(is,ie, js, je, isd,ied, jsd,jed, km, npz, ak0, bk0, psc, ud, vd, Atm)
               integer, parameter :: kord_mt = 9
               integer, parameter :: ikord_mt = 1
               type(fv_atmos_type), intent(inout) :: Atm
               integer, intent(in):: is,ie, js, je, isd,ied, jsd,jed, km, npz
               real(FVPRC),    intent(in):: ak0(km+1), bk0(km+1)
               real(FVPRC),    intent(in):: psc(is:ie,js:je)
               real(FVPRC),    intent(in), dimension(isd:ied  ,jsd:jed+1,km):: ud
               real(FVPRC),    intent(in), dimension(isd:ied+1,jsd:jed  ,km):: vd
! local:
               real(FVPRC), dimension(isd:ied,jsd:jed,npz):: ut, vt   ! winds
               real(FVPRC), dimension(isd:ied,jsd:jed, km+1):: pe0
               real(FVPRC), dimension(isd:ied,jsd:jed,npz+1):: pe1
               real(FVPRC), dimension(is:ie, km+1) :: pe0d
               real(FVPRC), dimension(is:ie,npz+1) :: pe1d
               real(FVPRC), dimension(is:ie, km  ) :: dpe0
               real(FVPRC), dimension(is:ie,npz  ) :: dpe1
               real(FVPRC), dimension(is:ie,npz):: qn1
               integer :: kord(1)
               integer i,j,k

               kord(ikord_mt) = kord_mt

               ut = 0.0
               vt = 0.0

               call prt_mxm('REMAP_WINDS: UD', ud, is,ie, js,je, Atm%ng, km, 1.0_FVPRC, Atm%gridstruct%area_64, Atm%domain)
               call prt_mxm('REMAP_WINDS: VD', vd, is,ie, js,je, Atm%ng, km, 1.0_FVPRC, Atm%gridstruct%area_64, Atm%domain)

               do k=1,km+1
                  do j=js,je
                     do i=is,ie
                        pe0(i,j,k) = ak0(k) + bk0(k)*psc(i,j)
                     enddo
                  enddo
               enddo
               call mpp_update_domains(pe0, Atm%domain) 

               do k=1,npz+1
                  do j=js,je
                     do i=is,ie
                        pe1(i,j,k) = Atm%ak(k) + Atm%bk(k)*Atm%ps(i,j)
                     enddo
                  enddo
               enddo
               call mpp_update_domains(pe1, Atm%domain)

               do 5000 j=js,je

!------
! map u
!------
                  do k=1,km+1
                     do i=is,ie
                        pe0d(i,k) = 0.5*(pe0(i,j-1,k)+pe0(i,j,k))
                     enddo
                  enddo
                  do k=1,km
                     do i=is,ie
                        dpe0(i,k) = pe0d(i,k+1)-pe0d(i,k)         
                     enddo
                  enddo
                  do k=1,npz+1
                     do i=is,ie
                        pe1d(i,k) = 0.5*(pe1(i,j-1,k)+pe1(i,j,k))                    
                     enddo
                  enddo
                  do k=1,npz
                     do i=is,ie
                        dpe1(i,k) = pe1d(i,k+1)-pe1d(i,k)             
                     enddo
                  enddo
                  call map_scalar( km, pe0d, ud(is:ie,j,1:km),    &
                                  npz, pe1d, qn1,                 &
                                  dpe0, dpe1, is, ie,             &
                                  j,  is, ie, j, j, -1, kord(ikord_mt), &
                        optional_top=.true., optional_bot=.true.)
                  do k=1,npz
                     do i=is,ie
                        ut(i,j,k) = qn1(i,k)
                     enddo
                  enddo
!------
! map v
!------
                  do k=1,km+1
                     do i=is,ie
                        pe0d(i,k) = 0.5*(pe0(i-1,j,k)+pe0(i,j,k))
                     enddo
                  enddo
                  do k=1,km
                     do i=is,ie
                        dpe0(i,k) = pe0d(i,k+1)-pe0d(i,k)
                     enddo
                  enddo
                  do k=1,npz+1
                     do i=is,ie
                        pe1d(i,k) = 0.5*(pe1(i-1,j,k)+pe1(i,j,k))                    
                     enddo
                  enddo
                  do k=1,npz
                     do i=is,ie
                        dpe1(i,k) = pe1d(i,k+1)-pe1d(i,k)         
                     enddo
                  enddo
                  call map_scalar( km, pe0d, vd(is:ie,j,1:km),    &
                                  npz, pe1d, qn1,                 &
                                  dpe0, dpe1, is, ie,             &
                                  j,  is, ie, j, j, -1, kord(ikord_mt), &
                        optional_top=.true., optional_bot=.true.)
                  do k=1,npz
                     do i=is,ie
                        vt(i,j,k) = qn1(i,k)
                     enddo
                  enddo

5000           continue

               call prt_mxm('REMAP_WINDS: U', ut, is, ie, js, je, Atm%ng, npz, 1.0_FVPRC, Atm%gridstruct%area_64, Atm%domain)
               call prt_mxm('REMAP_WINDS: V', vt, is, ie, js, je, Atm%ng, npz, 1.0_FVPRC, Atm%gridstruct%area_64, Atm%domain)

               Atm%u(isd:ied,jsd:jed,1:npz) = ut
               Atm%v(isd:ied,jsd:jed,1:npz) = vt

               end subroutine remap_winds

! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

subroutine xyz_to_dgrid(v3, ud, vd, npx, npy, is, ie, js, je, isd, ied, jsd, jed, gridstruct)

! Move A-Grid xyz winds to the D-grid cubed-sphere orientation
    
! !INPUT/OUTPUT PARAMETERS:
      integer, intent(in)        :: npx, npy, is, ie, js, je, isd, ied, jsd, jed
      real(REAL64)               :: v3(3, isd:ied  ,jsd:jed )
      real(FVPRC), intent(inout) :: ud(isd:ied,jsd:jed+1) ! U-Wind
      real(FVPRC), intent(inout) :: vd(isd:ied+1,jsd:jed) ! V-Wind
      type(fv_grid_type), intent(IN), target :: gridstruct
! !Local Variables 
      integer :: i,j, im2,jm2

      real(REAL64) :: ue(is-1:ie+1,js  :je+1,3)    ! 3D winds at edges
      real(REAL64) :: ve(is  :ie+1,js-1:je+1,3)    ! 3D winds at edges
      real(REAL64), dimension(is:ie):: ut1, ut2, ut3
      real(REAL64), dimension(js:je):: vt1, vt2, vt3

      im2 = (npx-1)/2
      jm2 = (npy-1)/2

! A --> D
! Interpolate to cell edges
       do j=js,je+1
          do i=is-1,ie+1
             ue(i,j,1) = 0.5*(v3(1,i,j-1) + v3(1,i,j))
             ue(i,j,2) = 0.5*(v3(2,i,j-1) + v3(2,i,j))
             ue(i,j,3) = 0.5*(v3(3,i,j-1) + v3(3,i,j))
          enddo
       enddo
       do j=js-1,je+1
          do i=is,ie+1
             ve(i,j,1) = 0.5*(v3(1,i-1,j) + v3(1,i,j))
             ve(i,j,2) = 0.5*(v3(2,i-1,j) + v3(2,i,j))
             ve(i,j,3) = 0.5*(v3(3,i-1,j) + v3(3,i,j))
          enddo
       enddo
! --- E_W edges (for v-wind):
     if ( is==1 ) then
       i = 1
       do j=js,je
        if ( j>jm2 ) then
             vt1(j) = gridstruct%edge_vect_w(j)*ve(i,j-1,1)+(1.-gridstruct%edge_vect_w(j))*ve(i,j,1)
             vt2(j) = gridstruct%edge_vect_w(j)*ve(i,j-1,2)+(1.-gridstruct%edge_vect_w(j))*ve(i,j,2)
             vt3(j) = gridstruct%edge_vect_w(j)*ve(i,j-1,3)+(1.-gridstruct%edge_vect_w(j))*ve(i,j,3)
        else
             vt1(j) = gridstruct%edge_vect_w(j)*ve(i,j+1,1)+(1.-gridstruct%edge_vect_w(j))*ve(i,j,1)
             vt2(j) = gridstruct%edge_vect_w(j)*ve(i,j+1,2)+(1.-gridstruct%edge_vect_w(j))*ve(i,j,2)
             vt3(j) = gridstruct%edge_vect_w(j)*ve(i,j+1,3)+(1.-gridstruct%edge_vect_w(j))*ve(i,j,3)
        endif
       enddo
       do j=js,je
          ve(i,j,1) = vt1(j)
          ve(i,j,2) = vt2(j)
          ve(i,j,3) = vt3(j)
       enddo
     endif
     if ( (ie+1)==npx ) then
       i = npx
       do j=js,je
        if ( j>jm2 ) then
             vt1(j) = gridstruct%edge_vect_e(j)*ve(i,j-1,1)+(1.-gridstruct%edge_vect_e(j))*ve(i,j,1)
             vt2(j) = gridstruct%edge_vect_e(j)*ve(i,j-1,2)+(1.-gridstruct%edge_vect_e(j))*ve(i,j,2)
             vt3(j) = gridstruct%edge_vect_e(j)*ve(i,j-1,3)+(1.-gridstruct%edge_vect_e(j))*ve(i,j,3)
        else
             vt1(j) = gridstruct%edge_vect_e(j)*ve(i,j+1,1)+(1.-gridstruct%edge_vect_e(j))*ve(i,j,1)
             vt2(j) = gridstruct%edge_vect_e(j)*ve(i,j+1,2)+(1.-gridstruct%edge_vect_e(j))*ve(i,j,2)
             vt3(j) = gridstruct%edge_vect_e(j)*ve(i,j+1,3)+(1.-gridstruct%edge_vect_e(j))*ve(i,j,3)
        endif
       enddo
       do j=js,je
          ve(i,j,1) = vt1(j)
          ve(i,j,2) = vt2(j)
          ve(i,j,3) = vt3(j)
       enddo
     endif
! N-S edges (for u-wind):
     if ( js==1 ) then
       j = 1
       do i=is,ie
        if ( i>im2 ) then
             ut1(i) = gridstruct%edge_vect_s(i)*ue(i-1,j,1)+(1.-gridstruct%edge_vect_s(i))*ue(i,j,1)
             ut2(i) = gridstruct%edge_vect_s(i)*ue(i-1,j,2)+(1.-gridstruct%edge_vect_s(i))*ue(i,j,2)
             ut3(i) = gridstruct%edge_vect_s(i)*ue(i-1,j,3)+(1.-gridstruct%edge_vect_s(i))*ue(i,j,3)
        else
             ut1(i) = gridstruct%edge_vect_s(i)*ue(i+1,j,1)+(1.-gridstruct%edge_vect_s(i))*ue(i,j,1)
             ut2(i) = gridstruct%edge_vect_s(i)*ue(i+1,j,2)+(1.-gridstruct%edge_vect_s(i))*ue(i,j,2)
             ut3(i) = gridstruct%edge_vect_s(i)*ue(i+1,j,3)+(1.-gridstruct%edge_vect_s(i))*ue(i,j,3)
        endif
       enddo
       do i=is,ie
          ue(i,j,1) = ut1(i)
          ue(i,j,2) = ut2(i)
          ue(i,j,3) = ut3(i)
       enddo
     endif
     if ( (je+1)==npy ) then
       j = npy
       do i=is,ie
        if ( i>im2 ) then
             ut1(i) = gridstruct%edge_vect_n(i)*ue(i-1,j,1)+(1.-gridstruct%edge_vect_n(i))*ue(i,j,1)
             ut2(i) = gridstruct%edge_vect_n(i)*ue(i-1,j,2)+(1.-gridstruct%edge_vect_n(i))*ue(i,j,2)
             ut3(i) = gridstruct%edge_vect_n(i)*ue(i-1,j,3)+(1.-gridstruct%edge_vect_n(i))*ue(i,j,3)
        else
             ut1(i) = gridstruct%edge_vect_n(i)*ue(i+1,j,1)+(1.-gridstruct%edge_vect_n(i))*ue(i,j,1)
             ut2(i) = gridstruct%edge_vect_n(i)*ue(i+1,j,2)+(1.-gridstruct%edge_vect_n(i))*ue(i,j,2)
             ut3(i) = gridstruct%edge_vect_n(i)*ue(i+1,j,3)+(1.-gridstruct%edge_vect_n(i))*ue(i,j,3)
        endif
       enddo
       do i=is,ie
          ue(i,j,1) = ut1(i)
          ue(i,j,2) = ut2(i)
          ue(i,j,3) = ut3(i)
       enddo
     endif
! Update:
       do j=js,je+1
          do i=is,ie
             ud(i,j) = ue(i,j,1)*gridstruct%es(1,i,j,1) +  &
                       ue(i,j,2)*gridstruct%es(2,i,j,1) +  &
                       ue(i,j,3)*gridstruct%es(3,i,j,1) 
          enddo
       enddo
       do j=js,je
          do i=is,ie+1
             vd(i,j) = ve(i,j,1)*gridstruct%ew(1,i,j,2) +  &
                       ve(i,j,2)*gridstruct%ew(2,i,j,2) +  &
                       ve(i,j,3)*gridstruct%ew(3,i,j,2) 
          enddo
       enddo

end subroutine xyz_to_dgrid

  subroutine d2a2d(ui, vi, uo, vo, Atm_i, Atm, regridder)
    !------------------------------------------------------------------!

    type(fv_atmos_type), intent(inout) :: Atm_i, Atm

    class(AbstractRegridder), pointer :: regridder

    real(FVPRC), dimension(Atm_i%bd%isd:Atm_i%bd%ied  ,Atm_i%bd%jsd:Atm_i%bd%jed+1), intent(in) :: ui
    real(FVPRC), dimension(Atm_i%bd%isd:Atm_i%bd%ied+1,Atm_i%bd%jsd:Atm_i%bd%jed  ), intent(in) :: vi

    real(FVPRC), dimension(Atm%bd%isd:Atm%bd%ied  ,Atm%bd%jsd:Atm%bd%jed+1), intent(inout) :: uo
    real(FVPRC), dimension(Atm%bd%isd:Atm%bd%ied+1,Atm%bd%jsd:Atm%bd%jed  ), intent(inout) :: vo

    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    real(REAL64), dimension(Atm_i%bd%isc:Atm_i%bd%iec  ,Atm_i%bd%jsc:Atm_i%bd%jec+1) :: wu
    real(REAL64), dimension(Atm_i%bd%isc:Atm_i%bd%iec+1,Atm_i%bd%jsc:Atm_i%bd%jec  ) :: wv

    real(REAL64), dimension(3,Atm_i%bd%isd:Atm_i%bd%ied,Atm_i%bd%jsd:Atm_i%bd%jed) :: va_xyz_i
    real(REAL32), dimension(  Atm_i%bd%isd:Atm_i%bd%ied,Atm_i%bd%jsd:Atm_i%bd%jed) :: tmp_i

    real(REAL64), dimension(3,Atm%bd%isd:Atm%bd%ied,Atm%bd%jsd:Atm%bd%jed) :: va_xyz_o
    real(REAL32), dimension(  Atm%bd%isd:Atm%bd%ied,Atm%bd%jsd:Atm%bd%jed) :: tmp_o

    real(REAL64) :: utmp, vtmp, ua, va

    integer :: i,j,n
    integer :: is, ie, js, je

    is = Atm_i%bd%is
    ie = Atm_i%bd%ie
    js = Atm_i%bd%js
    je = Atm_i%bd%je

!------------------------------------------------------------------!
! D -> A: co-variant u,d to contra-variant ua, va                  !
! vorticity preserving                                             !
!------------------------------------------------------------------!
       do j=js,je+1
          do i=is,ie
             wu(i,j) = ui(i,j)*Atm_i%gridstruct%dx(i,j)
          enddo
       enddo
       do j=js,je
          do i=is,ie+1
             wv(i,j) = vi(i,j)*Atm_i%gridstruct%dy(i,j)
          enddo
       enddo
       va_xyz_i = 0.0
       do j=js,je
          do i=is,ie
             !------------------------------------------------------------!
             ! D-Grid Cubed winds to cell center co-variant winds
             !------------------------------------------------------------!
             utmp = (wu(i,j) + wu(i,j+1)) * Atm_i%gridstruct%rdxa(i,j)
             vtmp = (wv(i,j) + wv(i+1,j)) * Atm_i%gridstruct%rdya(i,j)
             !------------------------------------------------------------!
             ! Cubed (cell center co-variant winds) to lat-lon:
             !------------------------------------------------------------!
             ua = Atm_i%gridstruct%a11(i,j)*utmp + Atm_i%gridstruct%a12(i,j)*vtmp
             va = Atm_i%gridstruct%a21(i,j)*utmp + Atm_i%gridstruct%a22(i,j)*vtmp
             !------------------------------------------------------------!
             ! latlon oriented vector winds:                              !
             !------------------------------------------------------------!
             va_xyz_i(1,i,j) = ua * Atm_i%gridstruct%vlon(i,j,1) + &
                               va * Atm_i%gridstruct%vlat(i,j,1)
             va_xyz_i(2,i,j) = ua * Atm_i%gridstruct%vlon(i,j,2) + &
                               va * Atm_i%gridstruct%vlat(i,j,2)
             va_xyz_i(3,i,j) = ua * Atm_i%gridstruct%vlon(i,j,3) + &
                               va * Atm_i%gridstruct%vlat(i,j,3)
          enddo
       enddo
!------------------------------------------------------------!
! Regrid from i->0 on A-grid vector winds
!------------------------------------------------------------!
       do n=1,3
          tmp_i = va_xyz_i(n,:,:)
          call regridder%regrid(tmp_i(       is:       ie,       js:       je), &
                                tmp_o(Atm%bd%is:Atm%bd%ie,Atm%bd%js:Atm%bd%je), rc=status)
          call mpp_update_domains(tmp_o, Atm%domain) 
          va_xyz_o(n,:,:) = tmp_o
       enddo
!------------------------------------------------------------!
! Move new winds from A-grid vector to D-grid co-variant
!------------------------------------------------------------!
       call xyz_to_dgrid(va_xyz_o, uo, vo, Atm%npx, Atm%npy, &
                         Atm%bd%is , Atm%bd%ie , Atm%bd%js , Atm%bd%je , &
                         Atm%bd%isd, Atm%bd%ied, Atm%bd%jsd, Atm%bd%jed, &
                         Atm%gridstruct)

  end subroutine d2a2d
                        
                  end module fv_regrid_c2c

