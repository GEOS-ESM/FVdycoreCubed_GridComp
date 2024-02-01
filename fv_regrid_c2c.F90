#include "MAPL_ErrLog.h"
module fv_regrid_c2c

#ifdef MAPL_MODE
#define DEALLOCGLOB_(A) if(associated(A)) then;A=0;if(MAPL_ShmInitialized) then; call MAPL_DeAllocNodeArray(A,rc=status);else; deallocate(A);endif;NULLIFY(A);endif
#endif

   use fv_arrays_mod,  only: REAL4, REAL8, FVPRC, R_GRID
   use fms_mod,            only: file_exist, read_data, field_exist
   use fms_io_mod,         only: get_tile_string, field_size
   use mpp_mod,            only: mpp_error, FATAL, NOTE, mpp_broadcast,mpp_npes
   use mpp_parameter_mod,  only: AGRID_PARAM=>AGRID
   use mpp_domains_mod,    only: mpp_get_tile_id, domain2d, mpp_update_domains, mpp_get_boundary, DGRID_NE
   use tracer_manager_mod, only: get_tracer_names, get_number_tracers, get_tracer_index
   use field_manager_mod,  only: MODEL_ATMOS

   use MAPL
   use gFTL_StringVector
   use gFTL_StringIntegerMap
   use, intrinsic :: iso_fortran_env, only: REAL64, REAL32

   use fv_arrays_mod,     only: fv_atmos_type, fv_grid_type, fv_grid_bounds_type, FVPRC, REAL4, REAL8
   use fv_diagnostics_mod,only: prt_maxmin
   use fv_mp_mod,         only: is_master, ng, mp_barrier, mp_gather, mp_bcst, &
         is,js,ie,je, isd,jsd,ied,jed, fill_corners, YDir
   use fv_grid_utils_mod, only: ptop_min
   use fv_grid_utils_mod, only: normalize_vect
   use fv_grid_utils_mod, only: direct_transform
   use fv_mapz_mod,       only: mappm
   use fv_surf_map_mod,   only: surfdrv
   use fv_timing_mod,     only: timing_on, timing_off
   use init_hydro_mod,    only: p_var
   use fv_timing_mod,      only: timing_on, timing_off
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
      real(real4), intent(inout) :: output(:,:)
      integer, intent(out), optional :: rc
      
      integer :: status,dims(3),funit
      integer :: rank
      type(ESMF_VM) :: vm
      real(real4), allocatable :: input(:,:)

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
      real(real8), intent(inout) :: output(:,:)
      integer, intent(out), optional :: rc
      
      integer :: status,dims(3),funit
      real, allocatable :: input(:,:)
      integer :: rank
      type(ESMF_VM) :: vm
      real(real8), allocatable :: input_r8(:,:)

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

   subroutine get_geos_ic( Atm, extra_rst, rstcube, gridOut, do_schmidt_in, schmidt_parameters_in)

      type(fv_atmos_type), intent(inout) :: Atm(:)
      type(fv_rst), pointer, intent(inout) :: extra_rst(:)
      logical :: rstcube
      type(ESMF_Grid), intent(inout) :: gridOut
      logical, intent(in) :: do_schmidt_in
      real, intent(in) :: schmidt_parameters_in(3) 
      real(FVPRC):: alpha
      integer i,j
      logical :: cubed_sphere,fv_diag_ic

! * Initialize coriolis param:

      do j=jsd,jed+1
         do i=isd,ied+1
            Atm(1)%gridstruct%fc(i,j) = 2.*omega*( -1.*cos(Atm(1)%gridstruct%grid(i,j,1))*cos(Atm(1)%gridstruct%grid(i,j,2))*sin(alpha) + &
                  sin(Atm(1)%gridstruct%grid(i,j,2))*cos(alpha) )
         enddo
      enddo

      do j=jsd,jed
         do i=isd,ied
            Atm(1)%gridstruct%f0(i,j) = 2.*omega*( -1.*cos(Atm(1)%gridstruct%agrid(i,j,1))*cos(Atm(1)%gridstruct%agrid(i,j,2))*sin(alpha) + &
                  sin(Atm(1)%gridstruct%agrid(i,j,2))*cos(alpha) )
         enddo
      enddo

      cubed_sphere=.true.
      fv_diag_ic=.false.
      call mpp_update_domains( Atm(1)%gridstruct%f0, Atm(1)%domain )
      if ( cubed_sphere ) call fill_corners(Atm(1)%gridstruct%f0, Atm(1)%npx, Atm(1)%npy, YDir)

      Atm(1)%phis = 0.

      do i=1,size(extra_rst)
         do j=1,size(extra_rst(i)%vars)
            call extra_rst(i)%vars(j)%alloc_var(isd,ied,jsd,jed)
         enddo
      enddo

      if (allocated(Atm(1)%q)) deallocate( Atm(1)%q )
      allocate  ( Atm(1)%q(isd:ied,jsd:jed,Atm(1)%npz,Atm(1)%ncnst) )
      call get_geos_cubed_ic( Atm, extra_rst, gridOut, do_schmidt_in, schmidt_parameters_in )

      call prt_maxmin('T', Atm(1)%pt, is, ie, js, je, ng, Atm(1)%npz, 1.0_FVPRC)

      Atm(1)%flagstruct%moist_phys=.false.
      call p_var(Atm(1)%npz,  is, ie, js, je, Atm(1)%ak(1),  ptop_min,         &
            Atm(1)%delp, Atm(1)%delz, Atm(1)%pt, Atm(1)%ps,               &
            Atm(1)%pe,   Atm(1)%peln, Atm(1)%pk, Atm(1)%pkz,              &
            kappa, Atm(1)%q, ng, Atm(1)%ncnst, dble(Atm(1)%gridstruct%area),Atm(1)%flagstruct%dry_mass,           &
            Atm(1)%flagstruct%adjust_dry_mass, Atm(1)%flagstruct%mountain, Atm(1)%flagstruct%moist_phys,   &
            Atm(1)%flagstruct%hydrostatic, Atm(1)%flagstruct%nwat, Atm(1)%domain,Atm(1)%flagstruct%make_nh, do_pkz=.false.)

   end subroutine get_geos_ic

   subroutine get_geos_cubed_ic( Atm, extra_rst, gridOut, do_schmidt_in, schmidt_parameters_in )
      use GHOST_CUBSPH_mod,  only : A_grid, ghost_cubsph_update
      use CUB2CUB_mod,    only: get_c2c_weight,                 &
            interpolate_c2c
      type(fv_atmos_type), intent(inout) :: Atm(:)
      type(fv_rst), pointer, intent(inout) :: extra_rst(:)
      type(ESMF_Grid), intent(inout) :: gridOut
      logical, intent(in) :: do_schmidt_in
      real, intent(in) :: schmidt_parameters_in(3)

      character(len=128) :: fname, fname1
      real(FVPRC), allocatable:: pkz0(:,:)
      real(FVPRC), allocatable:: ps0(:,:), gz0(:,:), t0(:,:,:), q0(:,:,:),qlev(:,:)
      real(FVPRC), allocatable:: u0(:,:,:), v0(:,:,:)
      real(FVPRC), allocatable:: ak0(:), bk0(:)
      integer :: i, j, k, l, iq, im, jm, km, npts, npx, npy, npz
      integer :: ntiles=6
      character (len=8) :: imc, jmc

      real(FVPRC) psc(is:ie,js:je)
      real(FVPRC) gzc(is:ie,js:je)
      real(FVPRC), allocatable:: tp(:,:,:), qp(:,:,:,:)
      real(FVPRC), allocatable:: ua(:,:,:), va(:,:,:)

      real(REAL64), allocatable :: akbk_r8(:)

      real(REAL64), dimension(:,:,:,:), allocatable :: corner_in, corner_out

      integer :: is_i,ie_i, js_i,je_i
      integer :: isd_i,ied_i, jsd_i,jed_i
      integer :: ng_i = 5
      type(domain2d) :: domain_i

      real(FVPRC), allocatable :: ebuffer(:,:)
      real(FVPRC), allocatable :: nbuffer(:,:)
      real(FVPRC), allocatable :: wbuffer(:,:)
      real(FVPRC), allocatable :: sbuffer(:,:)

      integer (kind=MPI_OFFSET_KIND) :: offset

      integer            :: filetype,nqmap
      type(Netcdf4_Fileformatter) :: formatter
      type(FileMetadata), allocatable :: cfg(:)
      integer            :: nDims, nVars, ivar, n_ungrid
      character(len=128) :: vname
      integer            :: tileoff,lvar_cnt,ifile,nlev
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
      type(CubedSphereGridFactory) :: cs_factory
      type(ESMF_Grid) :: gridIn
      class(AbstractRegridder), pointer :: regridder=>null()
      type(ArrDescr) :: input_arrdescr
      integer :: n_writers,n_readers,input_i1,input_in,input_j1,input_jn

      npx = Atm(1)%npx
      npy = Atm(1)%npy
      npz = Atm(1)%npz

! Zero out all initial tracer fields:
      Atm(1)%q = 0.

! Read input FV core restart file
      fname = "fvcore_internal_restart_in"

      call ESMF_AttributeGet(gridOut,name="num_reader",value=n_readers)
      if( file_exist(fname) ) then

         allocate(cfg(1))
         call formatter%open(fname,pFIO_READ,rc=status)
         cfg(1) = formatter%read(rc=status)
         im =cfg(1)%get_dimension('lon',rc=status)
         jm =cfg(1)%get_dimension('lat',rc=status)
         km =cfg(1)%get_dimension('lev',rc=status)

         if(is_master()) write(*,*) 'Using GEOS restart:', fname

         if ( file_exist(fname) ) then
            if(is_master())  write(*,*) 'External IC dimensions:', im   , jm       , km
            if(is_master())  write(*,*) 'Interpolating to      :', npx-1, (npy-1)*6, npz
         else
            call mpp_error(FATAL,'==> Error from get_geos_ic:        &
                  & field not found')
         endif

         if (do_schmidt_in) then
            cs_factory = CubedSphereGridFactory(nx=Atm(1)%layout(1),ny=Atm(1)%layout(2),im_world=im,lm=km, &
                         stretch_factor=schmidt_parameters_in(3), &
                             target_lon=schmidt_parameters_in(1),&
                             target_lat=schmidt_parameters_in(2))
         else
            cs_factory = CubedSphereGridFactory(nx=Atm(1)%layout(1),ny=Atm(1)%layout(2),im_world=im,lm=km)
         end if
         gridIn = grid_manager%make_grid(cs_factory,rc=status)
         call MAPL_Grid_Interior(gridIn,input_i1,input_in,input_j1,input_jn)
         call ArrDescrInit(input_arrdescr,MPI_COMM_WORLD,im,im*6,km, Atm(1)%layout(1),Atm(1)%layout(2)*6,n_readers, &
              n_readers,input_i1,input_in,input_j1,input_jn,rc=status)
         regridder => new_regridder_manager%make_regridder(gridIn,gridOut,REGRID_METHOD_BILINEAR,rc=status)
!--------------------------------------------------------------------!
! setup input cubed-sphere domain                                    !
!--------------------------------------------------------------------!
         call mpp_domain_decomp(domain_i,im+1,(jm/ntiles)+1,ntiles,ng_i,ntiles, &
               is_i,ie_i,js_i,je_i,isd_i,ied_i,jsd_i,jed_i,tile)
!--------------------------------------------------------------------!
! initialize cubed sphere grid: in                                   !
!--------------------------------------------------------------------!
         allocate(corner_in(2,is_i:ie_i+1,js_i:je_i+1,tile:tile))
         if (do_schmidt_in) then
            call init_cubsph_grid(im+1, is_i,ie_i, js_i,je_i, ntiles, corner_in, do_schmidt_in, real(schmidt_parameters_in,kind=R_GRID))
         else
            call init_cubsph_grid(im+1, is_i,ie_i, js_i,je_i, ntiles, corner_in, do_schmidt_in)
         end if
         call print_memuse_stats('get_geos_cubed_ic: init corner_in')
!--------------------------------------------------------------------!
! initialize cubed sphere grid: out                                  !
!--------------------------------------------------------------------!
         allocate(corner_out(2,is:ie+1,js:je+1,tile:tile))
         do l=tile,tile
            do j=js,je+1
               do i=is,ie+1
                  corner_out(1,i,j,l) = Atm(1)%gridstruct%grid(i,j,1)
                  corner_out(2,i,j,l) = Atm(1)%gridstruct%grid(i,j,2)
               enddo
            enddo
         enddo
         call print_memuse_stats('get_geos_cubed_ic: init corner_out')
         npts = im
         call print_memuse_stats('get_geos_cubed_ic: get_c2c_weight')

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

! Read U
         allocate (  u0(isd_i:ied_i,jsd_i:jed_i+1,km) )
         u0(:,:,:) = 0.0 
         tileoff = (tile-1)*(jm/ntiles)
         do k=1,km
            call MAPL_VarRead(formatter,"U",u0(is_i:ie_i,js_i:je_i,k),arrdes=input_arrdescr,lev=k)
         enddo
         call print_memuse_stats('get_geos_cubed_ic: read U')
! Read V
         allocate (  v0(isd_i:ied_i+1,jsd_i:jed_i,km) )
         v0(:,:,:) = 0.0
         tileoff = (tile-1)*(jm/ntiles)
         do k=1,km
            call MAPL_VarRead(formatter,"V",v0(is_i:ie_i,js_i:je_i,k),arrdes=input_arrdescr,lev=k)
         enddo
         call print_memuse_stats('get_geos_cubed_ic: read V')
         allocate ( sbuffer(is_i:ie_i,km) )
         allocate ( wbuffer(js_i:je_i,km) )
         allocate ( nbuffer(is_i:ie_i,km) )
         allocate ( ebuffer(js_i:je_i,km) )
         call mpp_get_boundary(u0, v0, domain_i, &
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
         call mpp_update_domains( u0, v0, domain_i, gridtype=DGRID_NE, complete=.true. )
         call prt_maxmin(' U_geos', u0, is_i, ie_i, js_i, je_i, ng_i, km, 1.0_FVPRC)
         call prt_maxmin(' V_geos', v0, is_i, ie_i, js_i, je_i, ng_i, km, 1.0_FVPRC)
         allocate ( ua(is:ie,js:je,km) )
         allocate ( va(is:ie,js:je,km) )
         call interp_c2c_vect(npts, npts, npx-1, npy-1, km, ntiles, domain_i, &
               is,ie, js,je, isd_i,ied_i, jsd_i,jed_i, is_i,ie_i, js_i,je_i, &
               u0, v0, ua, va, regridder, corner_in(:,is_i:ie_i+1,js_i:je_i+1,tile), corner_out, Atm(1)%gridstruct)
         deallocate ( v0 )
         deallocate ( u0 )
         deallocate ( corner_in )
         deallocate ( corner_out )
! Read T
         allocate (  t0(isd_i:ied_i,jsd_i:jed_i,km) )
         t0(:,:,:) = 0.0
         tileoff = (tile-1)*(jm/ntiles)
         do k=1,km
            call MAPL_VarRead(formatter,"PT",t0(is_i:ie_i,js_i:je_i,k),arrdes=input_arrdescr,lev=k)
         enddo
         call print_memuse_stats('get_geos_cubed_ic: read T')
! Read PE at Surface only
         allocate ( ps0(isd_i:ied_i,jsd_i:jed_i) )
         ps0(:,:) = 0.0
         tileoff = (tile-1)*(jm/ntiles)
         call MAPL_VarRead(formatter,"PE",ps0(is_i:ie_i,js_i:je_i),arrdes=input_arrdescr,lev=km+1)
         call mpp_update_domains(ps0, domain_i)
! Read PKZ
         allocate ( pkz0(isd_i:ied_i,jsd_i:jed_i) )
         pkz0(:,:) = 0.0
         tileoff = (tile-1)*(jm/ntiles)
         do k=1,km
            call MAPL_VarRead(formatter,"PKZ",pkz0(is_i:ie_i,js_i:je_i),arrdes=input_arrdescr,lev=k)
            t0(is_i:ie_i,js_i:je_i,k) = t0(is_i:ie_i,js_i:je_i,k)*pkz0(is_i:ie_i,js_i:je_i)
         enddo
         call print_memuse_stats('get_geos_cubed_ic: converted T')
         deallocate ( pkz0 )

         call formatter%close()
         deallocate(cfg)

         allocate ( gz0(isd_i:ied_i,jsd_i:jed_i) )
         gz0(:,:) = 0.0

         write(imc, "(i8)") im
         write(jmc, "(i8)") jm
         imc = adjustl(imc)
         jmc = adjustl(jmc)

         write(fname1, "('topo_DYN_ave_',a,'x',a,'.data')") trim(imc), trim(jmc)
         if (.not. file_exist(fname1)) then
            call mpp_error(FATAL,'get_geos_cubed_ic: cannot find topo_DYN_ave file')
         endif
         call print_memuse_stats('get_geos_cubed_ic: '//TRIM(fname1)//' being read')
         offset = 4
         call read_topo_file(fname1,gz0(is_i:ie_i,js_i:je_i),gridIn)
         !call parallel_read_file_r4(fname1, npts, is_i,ie_i, js_i,je_i, 1, offset, gz0(is_i:ie_i,js_i:je_i))
         call mpp_update_domains(gz0, domain_i)
         gz0 = gz0*grav

! Read cubed-sphere phis from file since IMPORT is not ready yet
         offset = 4
         !call parallel_read_file_r4('topo_dynave.data', Atm(1)%npx-1, is,ie, js,je, 1, offset, Atm(1)%phis(is:ie,js:je))
         call read_topo_file('topo_dynave.data',atm(1)%phis(is:ie,js:je),gridOut)
         call mpp_update_domains(Atm(1)%phis, Atm(1)%domain)
         Atm(1)%phis = Atm(1)%phis*grav
         call print_memuse_stats('get_geos_cubed_ic: phis')

! Horiz Interp for surface pressure 
         call prt_maxmin('PS_geos', ps0, is_i, ie_i, js_i, je_i, ng_i, 1, 1.0_FVPRC)
         call regridder%regrid(ps0(is_i:ie_i,js_i:je_i),psc(is:ie,js:je),rc=status)
         deallocate ( ps0 )
! Horiz Interp for surface height
         call prt_maxmin('GZ_geos', gz0, is_i, ie_i, js_i, je_i, ng_i, 1, 1.0/grav)
         call regridder%regrid(gz0(is_i:ie_i,js_i:je_i),gzc(is:ie,js:je),rc=status)
         deallocate ( gz0 )

! Horiz Interp for Q
         allocate ( q0(isd_i:ied_i,jsd_i:jed_i,km) )
         allocate ( qp(is:ie,js:je,km,Atm(1)%ncnst) )
         q0(:,:,:) = 0.0
         qp(:,:,:,:) = 0.0

! Horiz Interp for moist tracers
! is there a moist restart file to interpolate?
! Read in tracers: only sphum at this point
         if( file_exist("moist_internal_restart_in") ) then
            if (is_master()) print*, 'Trying to interpolate moist_internal_restart_in'

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
            tileoff = (tile-1)*(jm/ntiles)
            lvar_cnt=2

            do ivar=1,Atm(1)%ncnst
               vname = moist_variables%at(ivar)
               if (trim(vname)=='Q') then
                  iq=1
               else
                  iq=lvar_cnt
                  lvar_cnt=lvar_cnt+1
               end if
               call moist_tracers%insert(trim(vname),iq)
               do k=1,km
                  call MAPL_VarRead(formatter,vname,q0(is_i:ie_i,js_i:je_i,k),arrdes=input_arrdescr,lev=k)
                  call regridder%regrid(q0(is_i:ie_i,js_i:je_i,k),qp(:,:,k,iq),rc=status)
               enddo
               call prt_maxmin( 'Q_geos_moist', q0, is_i, ie_i, js_i, je_i, ng_i, km, 1._FVPRC)
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
            if (is_master()) print*, 'Trying to interpolate: ',trim(tracer_bundles(ifile)%file_name)

            allocate(cfg(1))
            call formatter%open(trim(tracer_bundles(ifile)%file_name),pFIO_READ,rc=status)
            cfg(1) = formatter%read(rc=status)
            call MAPL_IOCountNonDimVars(cfg(1),nvars,rc=status)
            tileoff = (tile-1)*(jm/ntiles)

            allocate ( qlev(isd_i:ied_i,jsd_i:jed_i) )
            qlev(:,:) = 0.0

            do ivar=1,size(tracer_bundles(ifile)%vars)
               nlev=tracer_bundles(ifile)%vars(ivar)%nLev
               vname = trim(tracer_bundles(ifile)%vars(ivar)%name)
               if (tracer_bundles(ifile)%vars(ivar)%rank ==2) then
                  call MAPL_VarRead(formatter,vname,qlev(is_i:ie_i,js_i:je_i),arrdes=input_arrdescr)
                  call regridder%regrid(qlev(is_i:ie_i,js_i:je_i),tracer_bundles(ifile)%vars(ivar)%ptr2d(is:ie,js:je),rc=status)
               else if (tracer_bundles(ifile)%vars(ivar)%rank ==3) then
                  do k=1,nlev
                     call MAPL_VarRead(formatter,vname,qlev(is_i:ie_i,js_i:je_i),arrdes=input_arrdescr,lev=k)
                     call regridder%regrid(qlev(is_i:ie_i,js_i:je_i),tracer_bundles(ifile)%vars(ivar)%ptr3d(is:ie,js:je,k),rc=status)
                  enddo
               else if (tracer_bundles(ifile)%vars(ivar)%rank ==4) then
                  do n_ungrid=1,tracer_bundles(ifile)%vars(ivar)%n_ungrid
                     do k=1,nlev
                        call MAPL_VarRead(formatter,vname,qlev(is_i:ie_i,js_i:je_i),arrdes=input_arrdescr,lev=k,offset2=n_ungrid)
                        call regridder%regrid(qlev(is_i:ie_i,js_i:je_i),tracer_bundles(ifile)%vars(ivar)%ptr4d(is:ie,js:je,k,n_ungrid),rc=status)
                     enddo
                  enddo
               end if
               !call prt_maxmin( 'Q_geos_gocart', q0, is_i, ie_i, js_i, je_i, ng_i, km, 1._FVPRC)
            enddo

            call formatter%close()
            deallocate(cfg)
            deallocate(qlev)

         enddo
                   
! Horiz Interp for T
         deallocate ( q0 )
         call mpp_update_domains(t0, domain_i)
         call prt_maxmin( 'T_geos', t0, is_i, ie_i, js_i, je_i, ng_i, km, 1.0_FVPRC)
         allocate (  tp(is:ie,js:je,km) )
         do k=1,km
            call regridder%regrid(t0(is_i:ie_i,js_i:je_i,k),tp(:,:,k),rc=status)
         enddo
         deallocate ( t0 )

! Horz/Vert remap for scalars
         nqmap =  Atm(1)%ncnst

         call remap_scalar(im, jm, km, npz, nqmap, nqmap, ak0, bk0, psc, gzc, tp, qp, Atm(1),tracer_bundles,extra_rst)

         deallocate ( tp )
         deallocate ( qp )
         call print_memuse_stats('get_geos_cubed_ic: remap_scalar')
! Horz/Vert remap for U/V
         call remap_winds(im, jm, km, npz, ak0, bk0, psc, ua, va, Atm)
         deallocate ( ua )
         deallocate ( va )
         call print_memuse_stats('get_geos_cubed_ic: remap_winds')

      else
         call mpp_error(FATAL,'==> Error from get_geos_ic:        &
               & Expected file '//trim(fname)//' does not exist')
      endif

      if (allocated(bk0)) deallocate ( bk0 )
      if (allocated(ak0)) deallocate ( ak0 )

! Finished, let's check the results !

      call prt_maxmin('GZ_model', Atm(1)%phis, is, ie, js, je, ng, 1, 1.0/grav)
      call prt_maxmin('PS_model', Atm(1)%ps, is, ie, js, je, ng, 1, 0.01_FVPRC)
      call prt_maxmin('DP_model', Atm(1)%delp, is, ie, js, je, ng, npz, 1.0_FVPRC)
      call prt_maxmin(' U_model', Atm(1)%u, is, ie, js, je, ng, npz, 1.0_FVPRC)
      call prt_maxmin(' V_model', Atm(1)%v, is, ie, js, je, ng, npz, 1.0_FVPRC)
      call prt_maxmin('PT_model', Atm(1)%pt, is, ie, js, je, ng, npz, 1.0_FVPRC)
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

 subroutine remap_coef( im, jm, lon, lat, id1, id2, jdc, s2c, agrid, bd )

   type(fv_grid_bounds_type), intent(IN) :: bd
  integer, intent(in):: im, jm
  real(FVPRC),    intent(in):: lon(im), lat(jm)
  real(FVPRC),    intent(out):: s2c(bd%is:bd%ie,bd%js:bd%je,4)
  integer, intent(out), dimension(bd%is:bd%ie,bd%js:bd%je):: id1, id2, jdc
  real(FVPRC),    intent(in):: agrid(bd%isd:bd%ied,bd%jsd:bd%jed,2)
! local:
  real(FVPRC) :: rdlon(im)
  real(FVPRC) :: rdlat(jm)
  real(FVPRC):: a1, b1
  integer i,j, i1, i2, jc, i0, j0

  integer :: is,  ie,  js,  je
  integer :: isd, ied, jsd, jed

  is  = bd%is
  ie  = bd%ie
  js  = bd%js
  je  = bd%je
  isd = bd%isd
  ied = bd%ied
  jsd = bd%jsd
  jed = bd%jed
  do i=1,im-1
     rdlon(i) = 1. / (lon(i+1) - lon(i))
  enddo
     rdlon(im) = 1. / (lon(1) + 2.*pi - lon(im))

  do j=1,jm-1
     rdlat(j) = 1. / (lat(j+1) - lat(j))
  enddo

! * Interpolate to cubed sphere cell center
  do 5000 j=js,je

     do i=is,ie

       if ( agrid(i,j,1)>lon(im) ) then
            i1 = im;     i2 = 1
            a1 = (agrid(i,j,1)-lon(im)) * rdlon(im)
       elseif ( agrid(i,j,1)<lon(1) ) then
            i1 = im;     i2 = 1
            a1 = (agrid(i,j,1)+2.*pi-lon(im)) * rdlon(im)
       else
            do i0=1,im-1
            if ( agrid(i,j,1)>=lon(i0) .and. agrid(i,j,1)<=lon(i0+1) ) then
               i1 = i0;  i2 = i0+1
               a1 = (agrid(i,j,1)-lon(i1)) * rdlon(i0)
               go to 111
            endif
            enddo
       endif
111    continue       
       if ( agrid(i,j,2)<lat(1) ) then
            jc = 1
            b1 = 0.
       elseif ( agrid(i,j,2)>lat(jm) ) then
            jc = jm-1
            b1 = 1.
       else
          do j0=1,jm-1
          if ( agrid(i,j,2)>=lat(j0) .and. agrid(i,j,2)<=lat(j0+1) ) then
               jc = j0
               b1 = (agrid(i,j,2)-lat(jc)) * rdlat(jc)
               go to 222
          endif
          enddo
       endif
222    continue

       if ( a1<0.0 .or. a1>1.0 .or.  b1<0.0 .or. b1>1.0 ) then
            write(*,*) 'gid=', mpp_pe(), i,j,a1, b1
       endif

       s2c(i,j,1) = (1.-a1) * (1.-b1)
       s2c(i,j,2) =     a1  * (1.-b1)
       s2c(i,j,3) =     a1  *     b1
       s2c(i,j,4) = (1.-a1) *     b1
       id1(i,j) = i1
       id2(i,j) = i2
       jdc(i,j) = jc
     enddo   !i-loop
5000 continue   ! j-loop

 end subroutine remap_coef

            subroutine remap_winds(im, jm, km, npz, ak0, bk0, psc, ua, va, Atm)
               type(fv_atmos_type), intent(inout) :: Atm(:)
               integer, intent(in):: im, jm, km, npz
               real(FVPRC),    intent(in):: ak0(km+1), bk0(km+1)
               real(FVPRC),    intent(in):: psc(is:ie,js:je)
               real(FVPRC),    intent(in), dimension(is:ie,js:je,km):: ua, va
! local:
               real(FVPRC), dimension(isd:ied,jsd:jed,npz):: ut, vt   ! winds
               real(FVPRC), dimension(is:ie, km+1):: pe0
               real(FVPRC), dimension(is:ie,npz+1):: pe1
               real(FVPRC), dimension(is:ie,npz):: qn1
               integer i,j,k

               ut = 0.0
               vt = 0.0
               do 5000 j=js,je

                  do k=1,km+1
                     do i=is,ie
                        pe0(i,k) = ak0(k) + bk0(k)*psc(i,j)
                     enddo
                  enddo

                  do k=1,npz+1
                     do i=is,ie
                        pe1(i,k) = Atm(1)%ak(k) + Atm(1)%bk(k)*Atm(1)%ps(i,j)
                     enddo
                  enddo

!------
! map u
!------
                  call mappm(km, pe0, ua(is:ie,j,1:km), npz, pe1, qn1, is,ie, -1, 4, Atm(1)%ptop)
                  do k=1,npz
                     do i=is,ie
                        ut(i,j,k) = qn1(i,k)
                     enddo
                  enddo
!------
! map v
!------
                  call mappm(km, pe0, va(is:ie,j,1:km), npz, pe1, qn1, is,ie, -1, 4, Atm(1)%ptop)
                  do k=1,npz
                     do i=is,ie
                        vt(i,j,k) = qn1(i,k)
                     enddo
                  enddo

5000              continue

                  call prt_maxmin('UT', ut, is, ie, js, je, ng, npz, 1.0_FVPRC)
                  call prt_maxmin('VT', vt, is, ie, js, je, ng, npz, 1.0_FVPRC)

!----------------------------------------------
! winds: lat-lon ON A to Cubed-D transformation:
!----------------------------------------------
                  call cubed_a2d(Atm(1)%npx, Atm(1)%npy, npz, ut, vt, Atm(1)%u, Atm(1)%v, Atm(1)%gridstruct, &
                  Atm(1)%domain, Atm(1)%bd )

                  if (is_master()) write(*,*) 'done remap_winds'

               end subroutine remap_winds


               subroutine remap_wz(im, jm, km, npz, mg, ak0, bk0, psc, wa, wz, Atm)
                  type(fv_atmos_type), intent(inout) :: Atm(:)
                  integer, intent(in):: im, jm, km, npz
                  integer, intent(in):: mg     ! mg = 0 for delz; mg=3 for w
                  real(FVPRC),    intent(in):: ak0(km+1), bk0(km+1)
                  real(FVPRC),    intent(in):: psc(is:ie,js:je)
                  real(FVPRC),    intent(in), dimension(is:ie,js:je,km):: wa
                  real(FVPRC),   intent(out):: wz(is-mg:ie+mg,js-mg:je+mg,npz)
! local:
                  real(FVPRC), dimension(is:ie, km+1):: pe0
                  real(FVPRC), dimension(is:ie,npz+1):: pe1
                  real(FVPRC), dimension(is:ie,npz):: qn1
                  integer i,j,k

                  do 5000 j=js,je

                     do k=1,km+1
                        do i=is,ie
                           pe0(i,k) = ak0(k) + bk0(k)*psc(i,j)
                        enddo
                     enddo

                     do k=1,npz+1
                        do i=is,ie
                           pe1(i,k) = Atm(1)%ak(k) + Atm(1)%bk(k)*Atm(1)%ps(i,j)
                        enddo
                     enddo

!------
! map w
!------
                     call mappm(km, pe0, wa(is:ie,j,1:km), npz, pe1, qn1, is,ie, -1, 4, Atm(1)%ptop)
                     do k=1,npz
                        do i=is,ie
                           wz(i,j,k) = qn1(i,k)
                        enddo
                     enddo

5000                 continue

! call prt_maxmin('WZ', wz, is, ie, js, je, mg, npz, 1._FVPRC, is_master())
! if (is_master()) write(*,*) 'done remap_wz'

                  end subroutine remap_wz

  subroutine remap_xyz( im, jbeg, jend, jm, km, npz, nq, ncnst, lon, lat, ak0, bk0, ps0, gz0,   &
                        ua, va, ta, qa, Atm )

  type(fv_atmos_type), intent(inout), target :: Atm
  integer, intent(in):: im, jm, km, npz, nq, ncnst
  integer, intent(in):: jbeg, jend
  real(FVPRC),    intent(in):: lon(im), lat(jm), ak0(km+1), bk0(km+1)
  real(FVPRC),    intent(in):: gz0(im,jbeg:jend), ps0(im,jbeg:jend)
  real(FVPRC),    intent(in), dimension(im,jbeg:jend,km):: ua, va, ta
  real(FVPRC),    intent(in), dimension(im,jbeg:jend,km,ncnst):: qa

  real(FVPRC), pointer, dimension(:,:,:) :: agrid

! local:
  real(FVPRC), dimension(Atm%bd%isd:Atm%bd%ied,Atm%bd%jsd:Atm%bd%jed,npz):: ut, vt   ! winds 
  real(FVPRC), dimension(Atm%bd%is:Atm%bd%ie,km):: up, vp, tp
  real(FVPRC), dimension(Atm%bd%is:Atm%bd%ie,km+1):: pe0, pn0
  real(FVPRC) pt0(km), gz(km+1), pk0(km+1)
  real(FVPRC) qp(Atm%bd%is:Atm%bd%ie,km,ncnst)
  real(FVPRC), dimension(Atm%bd%is:Atm%bd%ie,npz):: qn1
  real(FVPRC), dimension(Atm%bd%is:Atm%bd%ie,npz+1):: pe1, pn1
  real(FVPRC) :: rdlon(im)
  real(FVPRC) :: rdlat(jm)
  real(FVPRC):: a1, b1, c1, c2, c3, c4
  real(FVPRC):: gzc, psc, pst
  integer i,j,k, i1, i2, jc, i0, j0, iq
! integer  sphum, liq_wat, ice_wat, cld_amt
  integer  sphum
  integer :: is,  ie,  js,  je
  integer :: isd, ied, jsd, jed

  is  = Atm%bd%is
  ie  = Atm%bd%ie
  js  = Atm%bd%js
  je  = Atm%bd%je
  isd = Atm%bd%isd
  ied = Atm%bd%ied
  jsd = Atm%bd%jsd
  jed = Atm%bd%jed

  !!NOTE: Only Atm is used in this routine.
  agrid => Atm%gridstruct%agrid

  sphum   = get_tracer_index(MODEL_ATMOS, 'sphum')
! liq_wat = get_tracer_index(MODEL_ATMOS, 'liq_wat')
! ice_wat = get_tracer_index(MODEL_ATMOS, 'ice_wat')
! cld_amt = get_tracer_index(MODEL_ATMOS, 'cld_amt')

   if ( sphum/=1 ) then
        call mpp_error(FATAL,'SPHUM must be 1st tracer')
   endif
  pk0(1) = ak0(1)**kappa

  do i=1,im-1
     rdlon(i) = 1. / (lon(i+1) - lon(i))
  enddo
     rdlon(im) = 1. / (lon(1) + 2.*pi - lon(im))

  do j=1,jm-1
     rdlat(j) = 1. / (lat(j+1) - lat(j))
  enddo

! * Interpolate to cubed sphere cell center
  do 5000 j=js,je

     do i=is,ie
        pe0(i,1) = ak0(1)
        pn0(i,1) = log(ak0(1))
     enddo

     do i=is,ie

       if ( agrid(i,j,1)>lon(im) ) then
            i1 = im;     i2 = 1
            a1 = (agrid(i,j,1)-lon(im)) * rdlon(im)
       elseif ( agrid(i,j,1)<lon(1) ) then
            i1 = im;     i2 = 1
            a1 = (agrid(i,j,1)+2.*pi-lon(im)) * rdlon(im)
       else
            do i0=1,im-1
            if ( agrid(i,j,1)>=lon(i0) .and. agrid(i,j,1)<=lon(i0+1) ) then
               i1 = i0;  i2 = i0+1
               a1 = (agrid(i,j,1)-lon(i1)) * rdlon(i0)
               go to 111
            endif
            enddo
       endif

111    continue
       if ( agrid(i,j,2)<lat(1) ) then
            jc = 1
            b1 = 0.
       elseif ( agrid(i,j,2)>lat(jm) ) then
            jc = jm-1
            b1 = 1.
       else
          do j0=1,jm-1
          if ( agrid(i,j,2)>=lat(j0) .and. agrid(i,j,2)<=lat(j0+1) ) then
               jc = j0
               b1 = (agrid(i,j,2)-lat(jc)) * rdlat(jc)
               go to 222
          endif
          enddo
       endif
222    continue

#ifndef DEBUG_REMAP
       if ( a1<0.0 .or. a1>1.0 .or.  b1<0.0 .or. b1>1.0 ) then
            write(*,*) i,j,a1, b1
       endif
#endif
       c1 = (1.-a1) * (1.-b1)
       c2 =     a1  * (1.-b1)
       c3 =     a1  *     b1
       c4 = (1.-a1) *     b1

! Interpolated surface pressure
       psc = c1*ps0(i1,jc  ) + c2*ps0(i2,jc  ) +    &
             c3*ps0(i2,jc+1) + c4*ps0(i1,jc+1)

! Interpolated surface geopotential
       gzc = c1*gz0(i1,jc  ) + c2*gz0(i2,jc  ) +    &
             c3*gz0(i2,jc+1) + c4*gz0(i1,jc+1)

! 3D fields:
       do iq=1,ncnst
!          if ( iq==sphum .or. iq==liq_wat .or. iq==ice_wat .or. iq==cld_amt ) then
          do k=1,km
             qp(i,k,iq) = c1*qa(i1,jc,  k,iq) + c2*qa(i2,jc,  k,iq) +  &
                          c3*qa(i2,jc+1,k,iq) + c4*qa(i1,jc+1,k,iq)
          enddo
!          endif
       enddo
       do k=1,km
          up(i,k) = c1*ua(i1,jc,  k) + c2*ua(i2,jc,  k) +  &
                    c3*ua(i2,jc+1,k) + c4*ua(i1,jc+1,k)
          vp(i,k) = c1*va(i1,jc,  k) + c2*va(i2,jc,  k) +  &
                    c3*va(i2,jc+1,k) + c4*va(i1,jc+1,k)
          tp(i,k) = c1*ta(i1,jc,  k) + c2*ta(i2,jc,  k) +  &
                    c3*ta(i2,jc+1,k) + c4*ta(i1,jc+1,k)
! Virtual effect:
          tp(i,k) = tp(i,k)*(1.+zvir*qp(i,k,sphum))
       enddo
! Tracers:

       do k=2,km+1
          pe0(i,k) = ak0(k) + bk0(k)*psc
          pn0(i,k) = log(pe0(i,k))
          pk0(k) = pe0(i,k)**kappa
       enddo

#ifdef USE_DATA_ZS
       Atm%  ps(i,j) = psc
       Atm%phis(i,j) = gzc
#else

! * Adjust interpolated ps to model terrain
       gz(km+1) = gzc
       do k=km,1,-1
           gz(k) = gz(k+1) + rdgas*tp(i,k)*(pn0(i,k+1)-pn0(i,k))
       enddo
! Only lowest layer potential temp is needed
          pt0(km) = tp(i,km)/(pk0(km+1)-pk0(km))*(kappa*(pn0(i,km+1)-pn0(i,km)))
       if( Atm%phis(i,j)>gzc ) then
           do k=km,1,-1
              if( Atm%phis(i,j) <  gz(k)  .and.    &
                  Atm%phis(i,j) >= gz(k+1) ) then
                  pst = pk0(k) + (pk0(k+1)-pk0(k))*(gz(k)-Atm%phis(i,j))/(gz(k)-gz(k+1))
                  go to 123
              endif
           enddo
       else
! Extrapolation into the ground
           pst = pk0(km+1) + (gzc-Atm%phis(i,j))/(cp_air*pt0(km))
       endif

123    Atm%ps(i,j) = pst**(1./kappa)
#endif
     enddo   !i-loop


! * Compute delp from ps
     do i=is,ie
        pe1(i,1) = Atm%ak(1)
        pn1(i,1) = log(pe1(i,1))
     enddo
     do k=2,npz+1
       do i=is,ie
          pe1(i,k) = Atm%ak(k) + Atm%bk(k)*Atm%ps(i,j)
          pn1(i,k) = log(pe1(i,k))
       enddo
     enddo

     do k=1,npz
        do i=is,ie
           Atm%delp(i,j,k) = pe1(i,k+1) - pe1(i,k)
        enddo
     enddo

! Use kord=9 for winds; kord=11 for tracers
!------
! map u
!------
      call mappm(km, pe0, up, npz, pe1, qn1, is,ie, -1, 9, Atm%ptop)
      do k=1,npz
         do i=is,ie
            ut(i,j,k) = qn1(i,k)
         enddo
      enddo
!------
! map v
!------
      call mappm(km, pe0, vp, npz, pe1, qn1, is,ie, -1, 9, Atm%ptop)
      do k=1,npz
         do i=is,ie
            vt(i,j,k) = qn1(i,k)
         enddo
      enddo

!---------------
! map tracers
!----------------
      do iq=1,ncnst
! Note: AM2 physics tracers only
!         if ( iq==sphum .or. iq==liq_wat .or. iq==ice_wat .or. iq==cld_amt ) then
         call mappm(km, pe0, qp(is,1,iq), npz, pe1,  qn1, is,ie, 0, 11, Atm%ptop)
         do k=1,npz
            do i=is,ie
               Atm%q(i,j,k,iq) = qn1(i,k)
            enddo
         enddo
!         endif
      enddo

!-------------------------------------------------------------
! map virtual temperature using geopotential conserving scheme.
!-------------------------------------------------------------
      call mappm(km, pn0, tp, npz, pn1, qn1, is,ie, 1, 9, Atm%ptop)
      do k=1,npz
         do i=is,ie
            Atm%pt(i,j,k) = qn1(i,k)/(1.+zvir*Atm%q(i,j,k,sphum))
         enddo
      enddo

5000 continue

  call prt_maxmin('PS_model', Atm%ps, is, ie, js, je, ng, 1, 0.01_FVPRC)
  call prt_maxmin('UT', ut, is, ie, js, je, ng, npz, 1._FVPRC)
  call prt_maxmin('VT', vt, is, ie, js, je, ng, npz, 1._FVPRC)

!----------------------------------------------
! winds: lat-lon ON A to Cubed-D transformation:
!----------------------------------------------
  call cubed_a2d(Atm%npx, Atm%npy, npz, ut, vt, Atm%u, Atm%v, Atm%gridstruct, Atm%domain, Atm%bd )

  if (is_master()) write(*,*) 'done remap_xyz'

 end subroutine remap_xyz



                     subroutine init_cubsph_grid(npts, is,ie, js,je, ntiles, sph_corner, do_schmidt, schmidt_parameters)  
                        use mpp_mod, only: mpp_root_pe
!------------------------------------------------------------------!
! read/generate cubed sphere grid                                  !
! calculate cell center from cell corners                          !
!                                                                  !
! input:                                                           !
! npts, is,ie, js,je, ntiles       number of grid points and tiles !
!                                                                  !
! output:                                                          !
! sph_corner             cell corners in spherical coor            !
!------------------------------------------------------------------!
                        use GHOST_CUBSPH_mod, only: B_grid, ghost_cubsph_update             
                        use fv_grid_utils_mod, only : gnomonic_grids
                        use fv_grid_tools_mod, only : mirror_grid

                        integer, intent(in) :: npts, is,ie, js,je, ntiles
                        real*8, dimension(2,is:ie+1,js:je+1), intent(out) :: sph_corner
                        logical, intent(in) :: do_schmidt
                        real(kind=R_GRID), optional, intent(in) :: schmidt_parameters(3)
!------------------------------------------------------------------!
! local variables                                                  !
!------------------------------------------------------------------!
                        integer :: i, j, n, masterproc
                        real*8, pointer :: xs(:,:), ys(:,:)
                        real*8, pointer :: grid_in(:,:,:,:)
                        integer :: grid_type = 0
!------------------------------------------------------------------!
! create sph_corner                                                !
!------------------------------------------------------------------!
#ifdef FVREGRID_MAPL_MODE
! allocate global arrays (preferable in shared memory)
                        masterproc = mpp_root_pe()
                        if(MAPL_ShmInitialized) then
                           if (is_master()) write(*,*) 'Using MAPL_Shmem in external_ic: init_cubsph_grid' 
                           call MAPL_AllocNodeArray(grid_in,Shp=(/npts,npts,2,ntiles/),rc=STATUS)
                        else
                           if (is_master()) write(*,*) 'WARNING... in external_ic: Global grid allocate'
                           allocate( grid_in(npts,npts,2,ntiles) )
                        endif
                        if (is_master()) then
                          allocate( xs(npts,npts) )
                          allocate( ys(npts,npts) )
                          call gnomonic_grids(grid_type, npts-1, xs, ys)
                          do j=1,npts
                             do i=1,npts
                                grid_in(i,j,1,1) = xs(i,j)
                                grid_in(i,j,2,1) = ys(i,j)
                             enddo
                          enddo
                          deallocate ( xs )
                          deallocate ( ys )
                        endif
                        if(MAPL_ShmInitialized) then
                          call MAPL_SyncSharedMemory(rc=STATUS)
                          call MAPL_BroadcastToNodes( grid_in, N=size(grid_in), ROOT=masterproc, RC=status)
                          call MAPL_SyncSharedMemory(rc=STATUS)
                        else
                          call mpp_broadcast(grid_in, size(grid_in), masterproc)
                        endif
#else
                        allocate( xs(npts,npts) )
                        allocate( ys(npts,npts) )
                        allocate( grid_in(npts,npts,2,ntiles) )
                        call gnomonic_grids(grid_type, npts-1, xs, ys)
                        do j=1,npts
                           do i=1,npts
                              grid_in(i,j,1,1) = xs(i,j)
                              grid_in(i,j,2,1) = ys(i,j)
                           enddo
                        enddo
                        deallocate ( xs )
                        deallocate ( ys )
#endif

! mirror_grid assumes that the tile=1 is centered on equator and greenwich meridian Lon[-pi,pi]
                        call mirror_grid(grid_in, 0, npts, npts, 2, ntiles)

                        sph_corner(1,is:ie+1,js:je+1) = grid_in(is:ie+1,js:je+1,1,tile)
                        sph_corner(2,is:ie+1,js:je+1) = grid_in(is:ie+1,js:je+1,2,tile)
                        do j=js,je+1
                           do i=is,ie+1
!---------------------------------
! Shift the corner away from Japan
!---------------------------------
! This will result in the corner close to east coast of China
                              if ( .not.do_schmidt ) &
                                 sph_corner(1,i,j) = sph_corner(1,i,j) - pi/18.
                              if ( sph_corner(1,i,j) < 0. )              &
                                   sph_corner(1,i,j) = sph_corner(1,i,j) + 2.*pi
                              if (ABS(sph_corner(1,i,j)) < 1.e-10) sph_corner(1,i,j) = 0.0
                              if (ABS(sph_corner(2,i,j)) < 1.e-10) sph_corner(2,i,j) = 0.0
                           enddo
                        enddo
!------------------------
! Schmidt transformation:
!------------------------
                        if ( do_schmidt ) then
                           call direct_transform(schmidt_parameters(3), is, ie+1, js, je+1, &
                                schmidt_parameters(1), schmidt_parameters(2), &
                                tile, sph_corner(1,is:ie+1,js:je+1), sph_corner(2,is:ie+1,js:je+1))
                        endif

#ifdef FVREGRID_MAPL_MODE
                        call MAPL_SyncSharedMemory(rc=STATUS)
                        DEALLOCGLOB_(grid_in)
                        call MAPL_SyncSharedMemory(rc=STATUS)
#else
                        deallocate ( grid_in )
#endif
!------------------------------------------------------------------!
! do halo update                                                   !
!------------------------------------------------------------------!

                     end subroutine init_cubsph_grid

                     subroutine interp_c2c_vect(npx_in, npy_in, npx_out, npy_out, npz, ntiles, domain_i, &
                           is,ie, js,je, isd_i,ied_i, jsd_i,jed_i, is_i,ie_i, js_i,je_i, &
                           u_in, v_in, u_out, v_out, regridder, corner_in, corner_out, gridstruct)
                        use GRID_UTILS_mod, only: latlon2xyz
                        use GRID_UTILS_mod,   only: get_dx, get_dxa, get_dy, get_dya,     &
                              get_center_vect, get_west_vect,       &
                              get_south_vect, get_cosa_center
                        use FLOW_PROJ_mod,    only: d2a, d2a_vect, a2d_vect
                        use GHOST_CUBSPH_mod,  only : A_grid, ghost_cubsph_update
                        use CUB2CUB_mod,    only: do_c2c_interpolation
                        integer, intent(IN) :: npx_in, npy_in, npx_out, npy_out, npz, ntiles
                        type(domain2d), intent(INOUT) :: domain_i
                        integer, intent(IN) :: is,ie, js,je, isd_i,ied_i, jsd_i,jed_i, is_i,ie_i, js_i,je_i
                        class(AbstractRegridder), pointer :: regridder
                        real(REAL64),  intent(IN) :: corner_in(2,is_i:ie_i+1,js_i:je_i+1)
                        real(REAL64),  intent(IN) :: corner_out(2,is:ie+1,js:je+1)
                        real(FVPRC),  intent(IN) ::  u_in(isd_i:ied_i,jsd_i:jed_i+1,npz)
                        real(FVPRC),  intent(IN) ::  v_in(isd_i:ied_i+1,jsd_i:jed_i,npz)
                        real(FVPRC),  intent(OUT):: u_out(is:ie,js:je,npz)
                        real(FVPRC),  intent(OUT):: v_out(is:ie,js:je,npz)
                        type(fv_grid_type), intent(IN), target :: gridstruct

                        real(FVPRC) :: tmp(isd_i:ied_i,jsd_i:jed_i)

                        integer :: i,j,k,n
                        real(REAL64), dimension(:,:,:), allocatable :: vxyz_in, vxyz_out
                        real(REAL32), dimension(:,:), allocatable :: tmp_in, tmp_out
                        real(REAL64), dimension(:,:,:), allocatable :: ec1, ec2, ew1, ew2, es1, es2
                        real(REAL64), dimension(:,:)  , allocatable :: dx, dy, dxa, dya, rdxa, rdya, cosa_s, sina_s
                        real(FVPRC) :: u1, v1, vx, vy, vz
                        real(REAL64) :: pc1(3), pc2(3)

                        real(REAL64), dimension(:,:,:), allocatable :: xyz_corner_in, xyz_corner_out

                        character(len=64) :: strTxt

                        tmp = 0
!------------------------------------------------------------------!
! calculate xyz cell corners and cell centers                      !
!------------------------------------------------------------------!
                        allocate(xyz_corner_in (3, isd_i:ied_i+1, jsd_i:jed_i+1), &
                              xyz_corner_out(3, is   :ie   +1, js   :je   +1))
                        do j=js_i,je_i+1
                           do i=is_i,ie_i+1
                              call latlon2xyz(corner_in(:,i,j), xyz_corner_in(:,i,j))
                           enddo
                        enddo
                        do j=js,je+1
                           do i=is,ie+1
                              call latlon2xyz(corner_out(:,i,j), xyz_corner_out(:,i,j))
                           enddo
                        enddo
                        call print_memuse_stats('interp_c2c_vect: GRIDS')

!----------------------------------------------------------!
! allocate horizontal flow variables                       !
!----------------------------------------------------------!
                        allocate(vxyz_in(3,isd_i:ied_i,jsd_i:jed_i))
                        allocate(dx(isd_i:ied_i,jsd_i:jed_i+1), dxa(isd_i:ied_i,jsd_i:jed_i), rdxa(isd_i:ied_i,jsd_i:jed_i))
                        allocate(dy(isd_i:ied_i+1,jsd_i:jed_i), dya(isd_i:ied_i,jsd_i:jed_i), rdya(isd_i:ied_i,jsd_i:jed_i))
                        allocate(ec1(3,isd_i:ied_i,jsd_i:jed_i), ec2(3,isd_i:ied_i,jsd_i:jed_i))
                        allocate(cosa_s(isd_i:ied_i,jsd_i:jed_i), sina_s(isd_i:ied_i,jsd_i:jed_i))
                        allocate(vxyz_out(3,is:ie,js:je))
                        allocate(tmp_in(is_i:ie_i,js_i:je_i))
                        allocate(tmp_out(is:ie,js:je))

!-------------------------------------------------------!
! geometrical properties of input grid                  !
!-------------------------------------------------------!
                        call get_dx (xyz_corner_in(:,isd_i:ied_i+1,jsd_i:jed_i+1), &
                              isd_i,ied_i  ,jsd_i,jed_i, &
                              is_i ,ie_i   ,js_i ,je_i , dx)
                        call get_dxa(xyz_corner_in(:,isd_i:ied_i+1,jsd_i:jed_i+1), &
                              isd_i,ied_i  ,jsd_i,jed_i, &
                              is_i ,ie_i   ,js_i ,je_i , dxa, rdxa=rdxa)
                        call get_dy (xyz_corner_in(:,isd_i:ied_i+1,jsd_i:jed_i+1), &
                              isd_i,ied_i  ,jsd_i,jed_i, &
                              is_i ,ie_i   ,js_i ,je_i , dy)
                        call get_dya(xyz_corner_in(:,isd_i:ied_i+1,jsd_i:jed_i+1), &
                              isd_i,ied_i  ,jsd_i,jed_i, &
                              is_i ,ie_i   ,js_i ,je_i , dya, rdya=rdya)
                        call get_center_vect(xyz_corner_in(:,isd_i:ied_i+1,jsd_i:jed_i+1), &
                              isd_i,ied_i  ,jsd_i,jed_i, &
                              is_i ,ie_i   ,js_i ,je_i , ec1, ec2)
                        call get_cosa_center(ec1, ec2, isd_i,ied_i  ,jsd_i,jed_i, &
                              is_i ,ie_i   ,js_i ,je_i , cosa_s, sina_s)

! Flow interpolation for U and V components
                        do k=1,npz
!-------------------------------------------------------!
! calculate flow vector for a-grid                      !
!-------------------------------------------------------!
                           call d2a_vect(DBLE(u_in(:,:,k)), DBLE(v_in(:,:,k)), DBLE(dx), DBLE(dy), DBLE(rdxa), DBLE(rdya), DBLE(cosa_s), DBLE(ec1), DBLE(ec2), &
                                 isd_i, ied_i, jsd_i, jed_i, 1, 1,           &
                                 is_i , ie_i , js_i , je_i , 1, 1,           &
                                 vxyz_in(:,isd_i:ied_i,jsd_i:jed_i))
                           write(strTxt,'(A,i3.3)') 'interp_c2c_vect: INPUT D2A level:', k
                           if (k==npz) call print_memuse_stats(strTxt)
!----------------------------------------------------------!
! ghost cell update of vxyz_in                               !
!----------------------------------------------------------!
                           do n=1,3
                              tmp(is_i:ie_i,js_i:je_i) = vxyz_in(n,is_i:ie_i,js_i:je_i)
                              call mpp_update_domains(tmp, domain_i)
                              vxyz_in(n,:,:) = tmp
                           enddo
                           do n=1,3
                              tmp_in=vxyz_in(n,is_i:ie_i,js_i:je_i)
                              call regridder%regrid(tmp_in,tmp_out,rc=status)
                              vxyz_out(n,:,:)=tmp_out
                           enddo
                           do j=js,je
                              do i=is,ie
                                 vx = vxyz_out(1,i,j)
                                 vy = vxyz_out(2,i,j)
                                 vz = vxyz_out(3,i,j)
!----------------------------------------------------------!
! convert flow vector to wind vectors on A-Grid            !
!----------------------------------------------------------!
                                 pc1(:)=xyz_corner_out(:,i+1,j)+xyz_corner_out(:,i+1,j+1)          &
                                       -xyz_corner_out(:,i  ,j)-xyz_corner_out(:,i  ,j+1)
                                 pc2(:)=xyz_corner_out(:,i,j+1)+xyz_corner_out(:,i+1,j+1)          &
                                       -xyz_corner_out(:,i,j  )-xyz_corner_out(:,i+1,j  )
                                 call normalize_vect(pc1(:))
                                 call normalize_vect(pc2(:))
                                 u1 = vx*pc1(1) + vy*pc1(2) + vz*pc1(3)
                                 v1 = vx*pc2(1) + vy*pc2(2) + vz*pc2(3)
!----------------------------------------------------------!
! rotate wind vectors from cubed to latlon orientation     !
!----------------------------------------------------------!
                                 u_out(i,j,k) = 2.0*(gridstruct%a11(i,j)*u1 + gridstruct%a12(i,j)*v1)
                                 v_out(i,j,k) = 2.0*(gridstruct%a21(i,j)*u1 + gridstruct%a22(i,j)*v1)
                              enddo
                           enddo
                           write(strTxt,'(A,i3.3)') 'interp_c2c_vect: OUTPUT A-grid C2L level:', k
                           if (k==npz) call print_memuse_stats(strTxt)
                        enddo ! npz

                        deallocate(dx, dy, dxa, dya, rdxa, rdya, ec1, ec2, cosa_s, sina_s)
                        deallocate(vxyz_in)
                        deallocate(vxyz_out)
                        deallocate ( xyz_corner_in, xyz_corner_out )

                     end subroutine interp_c2c_vect

                     Function REVERSE(A) Result(B)
                        real(FVPRC), Intent(In) :: A(:,:)
                        real(FVPRC) :: B(Size(A,1),Size(A,2))

                        Integer :: i, n

                        n = Size(A, 1)

                        Do i = 1, n
                           B(i,:) = A(1+n-i,:)
                        End Do

                     End Function REVERSE

 subroutine cubed_a2d( npx, npy, npz, ua, va, u, v, gridstruct, fv_domain, bd )

! Purpose; Transform wind on A grid to D grid

  use mpp_domains_mod,    only: mpp_update_domains

  type(fv_grid_bounds_type), intent(IN) :: bd
  integer, intent(in):: npx, npy, npz
  real(FVPRC), intent(inout), dimension(bd%isd:bd%ied,bd%jsd:bd%jed,npz):: ua, va
  real(FVPRC), intent(out):: u(bd%isd:bd%ied,  bd%jsd:bd%jed+1,npz)
  real(FVPRC), intent(out):: v(bd%isd:bd%ied+1,bd%jsd:bd%jed  ,npz)
  type(fv_grid_type), intent(IN), target :: gridstruct
  type(domain2d), intent(INOUT) :: fv_domain
! local:
  real(FVPRC) v3(3,bd%is-1:bd%ie+1,bd%js-1:bd%je+1)
  real(FVPRC) ue(3,bd%is-1:bd%ie+1,bd%js:bd%je+1)    ! 3D winds at edges
  real(FVPRC) ve(3,bd%is:bd%ie+1,bd%js-1:bd%je+1)    ! 3D winds at edges
  real(FVPRC), dimension(bd%is:bd%ie):: ut1, ut2, ut3
  real(FVPRC), dimension(bd%js:bd%je):: vt1, vt2, vt3
  integer i, j, k, im2, jm2

  real(REAL64), pointer, dimension(:,:,:)   :: vlon, vlat
  real(REAL64), pointer, dimension(:)       :: edge_vect_w, edge_vect_e, edge_vect_s, edge_vect_n
  real(REAL64), pointer, dimension(:,:,:,:) :: ew, es

  integer :: is,  ie,  js,  je
  integer :: isd, ied, jsd, jed

  is  = bd%is
  ie  = bd%ie
  js  = bd%js
  je  = bd%je
  isd = bd%isd
  ied = bd%ied
  jsd = bd%jsd
  jed = bd%jed
  vlon => gridstruct%vlon
  vlat => gridstruct%vlat

  edge_vect_w => gridstruct%edge_vect_w
  edge_vect_e => gridstruct%edge_vect_e
  edge_vect_s => gridstruct%edge_vect_s
  edge_vect_n => gridstruct%edge_vect_n

  ew => gridstruct%ew
  es => gridstruct%es

  call mpp_update_domains(ua, fv_domain, complete=.false.)
  call mpp_update_domains(va, fv_domain, complete=.true.)

    im2 = (npx-1)/2
    jm2 = (npy-1)/2
    do k=1, npz
! Compute 3D wind on A grid
       do j=js-1,je+1
          do i=is-1,ie+1
             !v3(1,i,j) = ua(i,j,k)*vlon(1,i,j) + va(i,j,k)*vlat(1,i,j)
             !v3(2,i,j) = ua(i,j,k)*vlon(2,i,j) + va(i,j,k)*vlat(2,i,j)
             !v3(3,i,j) = ua(i,j,k)*vlon(3,i,j) + va(i,j,k)*vlat(3,i,j)
             v3(1,i,j) = ua(i,j,k)*vlon(i,j,1) + va(i,j,k)*vlat(i,j,1)
             v3(2,i,j) = ua(i,j,k)*vlon(i,j,2) + va(i,j,k)*vlat(i,j,2)
             v3(3,i,j) = ua(i,j,k)*vlon(i,j,3) + va(i,j,k)*vlat(i,j,3)
          enddo
       enddo

! A --> D
! Interpolate to cell edges
       do j=js,je+1
          do i=is-1,ie+1
             ue(1,i,j) = 0.5*(v3(1,i,j-1) + v3(1,i,j))
             ue(2,i,j) = 0.5*(v3(2,i,j-1) + v3(2,i,j))
             ue(3,i,j) = 0.5*(v3(3,i,j-1) + v3(3,i,j))
          enddo
       enddo

       do j=js-1,je+1
          do i=is,ie+1
             ve(1,i,j) = 0.5*(v3(1,i-1,j) + v3(1,i,j))
             ve(2,i,j) = 0.5*(v3(2,i-1,j) + v3(2,i,j))
             ve(3,i,j) = 0.5*(v3(3,i-1,j) + v3(3,i,j))
          enddo
       enddo

! --- E_W edges (for v-wind):
     if (.not. gridstruct%nested) then
     if ( is==1) then
       i = 1
       do j=js,je
        if ( j>jm2 ) then
             vt1(j) = edge_vect_w(j)*ve(1,i,j-1)+(1.-edge_vect_w(j))*ve(1,i,j)
             vt2(j) = edge_vect_w(j)*ve(2,i,j-1)+(1.-edge_vect_w(j))*ve(2,i,j)
             vt3(j) = edge_vect_w(j)*ve(3,i,j-1)+(1.-edge_vect_w(j))*ve(3,i,j)
        else
             vt1(j) = edge_vect_w(j)*ve(1,i,j+1)+(1.-edge_vect_w(j))*ve(1,i,j)
             vt2(j) = edge_vect_w(j)*ve(2,i,j+1)+(1.-edge_vect_w(j))*ve(2,i,j)
             vt3(j) = edge_vect_w(j)*ve(3,i,j+1)+(1.-edge_vect_w(j))*ve(3,i,j)
        endif
       enddo
       do j=js,je
          ve(1,i,j) = vt1(j)
          ve(2,i,j) = vt2(j)
          ve(3,i,j) = vt3(j)
       enddo
     endif

     if ( (ie+1)==npx ) then
       i = npx
       do j=js,je
        if ( j>jm2 ) then
             vt1(j) = edge_vect_e(j)*ve(1,i,j-1)+(1.-edge_vect_e(j))*ve(1,i,j)
             vt2(j) = edge_vect_e(j)*ve(2,i,j-1)+(1.-edge_vect_e(j))*ve(2,i,j)
             vt3(j) = edge_vect_e(j)*ve(3,i,j-1)+(1.-edge_vect_e(j))*ve(3,i,j)
        else
             vt1(j) = edge_vect_e(j)*ve(1,i,j+1)+(1.-edge_vect_e(j))*ve(1,i,j)
             vt2(j) = edge_vect_e(j)*ve(2,i,j+1)+(1.-edge_vect_e(j))*ve(2,i,j)
             vt3(j) = edge_vect_e(j)*ve(3,i,j+1)+(1.-edge_vect_e(j))*ve(3,i,j)
        endif
       enddo
       do j=js,je
          ve(1,i,j) = vt1(j)
          ve(2,i,j) = vt2(j)
          ve(3,i,j) = vt3(j)
       enddo
     endif

! N-S edges (for u-wind):
     if ( js==1 ) then
       j = 1
       do i=is,ie
        if ( i>im2 ) then
             ut1(i) = edge_vect_s(i)*ue(1,i-1,j)+(1.-edge_vect_s(i))*ue(1,i,j)
             ut2(i) = edge_vect_s(i)*ue(2,i-1,j)+(1.-edge_vect_s(i))*ue(2,i,j)
             ut3(i) = edge_vect_s(i)*ue(3,i-1,j)+(1.-edge_vect_s(i))*ue(3,i,j)
        else
             ut1(i) = edge_vect_s(i)*ue(1,i+1,j)+(1.-edge_vect_s(i))*ue(1,i,j)
             ut2(i) = edge_vect_s(i)*ue(2,i+1,j)+(1.-edge_vect_s(i))*ue(2,i,j)
             ut3(i) = edge_vect_s(i)*ue(3,i+1,j)+(1.-edge_vect_s(i))*ue(3,i,j)
        endif
       enddo
       do i=is,ie
          ue(1,i,j) = ut1(i)
          ue(2,i,j) = ut2(i)
          ue(3,i,j) = ut3(i)
       enddo
     endif

     if ( (je+1)==npy ) then
       j = npy
       do i=is,ie
        if ( i>im2 ) then
             ut1(i) = edge_vect_n(i)*ue(1,i-1,j)+(1.-edge_vect_n(i))*ue(1,i,j)
             ut2(i) = edge_vect_n(i)*ue(2,i-1,j)+(1.-edge_vect_n(i))*ue(2,i,j)
             ut3(i) = edge_vect_n(i)*ue(3,i-1,j)+(1.-edge_vect_n(i))*ue(3,i,j)
        else
             ut1(i) = edge_vect_n(i)*ue(1,i+1,j)+(1.-edge_vect_n(i))*ue(1,i,j)
             ut2(i) = edge_vect_n(i)*ue(2,i+1,j)+(1.-edge_vect_n(i))*ue(2,i,j)
             ut3(i) = edge_vect_n(i)*ue(3,i+1,j)+(1.-edge_vect_n(i))*ue(3,i,j)
        endif
       enddo
       do i=is,ie
          ue(1,i,j) = ut1(i)
          ue(2,i,j) = ut2(i)
          ue(3,i,j) = ut3(i)
       enddo
     endif

     endif ! .not. nested
     do j=js,je+1
        do i=is,ie
           u(i,j,k) =  ue(1,i,j)*es(1,i,j,1) +  &
                       ue(2,i,j)*es(2,i,j,1) +  &
                       ue(3,i,j)*es(3,i,j,1)
        enddo
     enddo
     do j=js,je
        do i=is,ie+1
           v(i,j,k) = ve(1,i,j)*ew(1,i,j,2) +  &
                      ve(2,i,j)*ew(2,i,j,2) +  &
                      ve(3,i,j)*ew(3,i,j,2)
        enddo
     enddo

   enddo         ! k-loop

 end subroutine cubed_a2d


                     subroutine d2a3d(u, v,  ua,   va,  im,  jm, km, lon)
                        integer, intent(in):: im, jm, km           ! Dimensions
                        real(FVPRC), intent(in ) :: lon(im)
                        real(FVPRC), intent(in ), dimension(im,jm,km):: u, v
                        real(FVPRC), intent(out), dimension(im,jm,km):: ua, va
! local
                        real(FVPRC) :: coslon(im),sinlon(im)    ! Sine and cosine in longitude
                        integer i, j, k
                        integer imh
                        real(FVPRC) un, vn, us, vs

                        integer :: ks, ke

                        imh = im/2

                        do i=1,im
                           sinlon(i) = sin(lon(i))
                           coslon(i) = cos(lon(i))
                        enddo

                        do k=1,km
                           do j=2,jm-1
                              do i=1,im
                                 ua(i,j,k) = 0.5*(u(i,j,k) + u(i,j+1,k))
                              enddo
                           enddo

                           do j=2,jm-1
                              do i=1,im-1
                                 va(i,j,k) = 0.5*(v(i,j,k) + v(i+1,j,k))
                              enddo
                              va(im,j,k) = 0.5*(v(im,j,k) + v(1,j,k))
                           enddo

! Projection at SP
                           us = 0.
                           vs = 0.
                           do i=1,imh
                              us = us + (ua(i+imh,2,k)-ua(i,2,k))*sinlon(i)      &
                                    + (va(i,2,k)-va(i+imh,2,k))*coslon(i)
                              vs = vs + (ua(i+imh,2,k)-ua(i,2,k))*coslon(i)      &
                                    + (va(i+imh,2,k)-va(i,2,k))*sinlon(i)
                           enddo
                           us = us/im
                           vs = vs/im
                           do i=1,imh
                              ua(i,1,k)   = -us*sinlon(i) - vs*coslon(i)
                              va(i,1,k)   =  us*coslon(i) - vs*sinlon(i)
                              ua(i+imh,1,k)   = -ua(i,1,k)
                              va(i+imh,1,k)   = -va(i,1,k)
                           enddo

! Projection at NP
                           un = 0.
                           vn = 0.
                           do i=1,imh
                              un = un + (ua(i+imh,jm-1,k)-ua(i,jm-1,k))*sinlon(i)    &
                                    + (va(i+imh,jm-1,k)-va(i,jm-1,k))*coslon(i)
                              vn = vn + (ua(i,jm-1,k)-ua(i+imh,jm-1,k))*coslon(i)    &
                                    + (va(i+imh,jm-1,k)-va(i,jm-1,k))*sinlon(i)
                           enddo

                           un = un/im
                           vn = vn/im
                           do i=1,imh
                              ua(i,jm,k) = -un*sinlon(i) + vn*coslon(i)
                              va(i,jm,k) = -un*coslon(i) - vn*sinlon(i)
                              ua(i+imh,jm,k) = -ua(i,jm,k)
                              va(i+imh,jm,k) = -va(i,jm,k)
                           enddo
                        enddo

                     end subroutine d2a3d



                     subroutine pmaxmin( qname, a, im, jm, fac )

                        integer, intent(in):: im, jm
                        character(len=*) :: qname
                        integer i, j
                        real(FVPRC) a(im,jm)

                        real(FVPRC) qmin(jm), qmax(jm)
                        real(FVPRC) pmax, pmin
                        real(FVPRC) fac                     ! multiplication factor

                        do j=1,jm
                           pmax = a(1,j)
                           pmin = a(1,j)
                           do i=2,im
                              pmax = max(pmax, a(i,j))
                              pmin = min(pmin, a(i,j))
                           enddo
                           qmax(j) = pmax
                           qmin(j) = pmin
                        enddo
!
! Now find max/min of amax/amin
!
                        pmax = qmax(1)
                        pmin = qmin(1)
                        do j=2,jm
                           pmax = max(pmax, qmax(j))
                           pmin = min(pmin, qmin(j))
                        enddo

                        write(*,*) qname, ' max = ', pmax*fac, ' min = ', pmin*fac

                     end subroutine pmaxmin

                     subroutine pmaxmin4d( qname, a, im, jm, km, lm, fac )

                        character*(*)  qname 
                        integer, intent(in):: im, jm, km, lm
                        integer i, j, k, l
                        real(FVPRC) a(im,jm,km,lm)

                        real(FVPRC) qmin(jm), qmax(jm)
                        real(FVPRC) pmax, pmin
                        real(FVPRC) fac                     ! multiplication factor

                        pmax = a(1,1,1,1)
                        pmin = a(1,1,1,1)
                        do l=1,lm
                           do k=1,km
                              do j=1,jm
                                 do i=1,im
                                    pmax = max(pmax, a(i,j,k,l))
                                    pmin = min(pmin, a(i,j,k,l))
                                 enddo
                              enddo
                           enddo
                        enddo

                        write(*,*) qname, ' max = ', pmax*fac, ' min = ', pmin*fac

                     end subroutine pmaxmin4d


                     subroutine mpp_domain_decomp(domain,npx,npy,nregions,ng,grid_type, &
                           is,ie,js,je,isd,ied,jsd,jed,tile)
                        use mpp_domains_mod, only: mpp_domains_init, MPP_DOMAIN_TIME, &
                              mpp_define_mosaic, mpp_get_compute_domain, mpp_get_data_domain, &
                              mpp_domains_set_stack_size, mpp_define_layout
                        use mpp_mod,         only : mpp_pe
                        type(domain2D), intent(OUT) :: domain
                        integer, intent(IN)  :: npx,npy,nregions,ng,grid_type
                        integer, intent(OUT) :: is,ie,js,je,isd,ied,jsd,jed,tile

                        integer :: layout(2)
                        integer, allocatable :: pe_start(:), pe_end(:)

                        integer :: num_contact, ntiles, npes_per_tile
                        integer, allocatable, dimension(:)       :: npes_tile, tile1, tile2
                        integer, allocatable, dimension(:)       :: istart1, iend1, jstart1, jend1
                        integer, allocatable, dimension(:)       :: istart2, iend2, jstart2, jend2
                        integer, allocatable, dimension(:,:)     :: layout2D, global_indices

                        character*80 :: evalue
                        integer :: ios,nx,ny,n,num_alloc
                        character(len=32) :: type = "unknown"
                        logical :: is_symmetry
                        integer :: npes

                        npes = mpp_npes()

                        nx = npx-1
                        ny = npy-1

                        call print_memuse_stats('external_ic:mpp_domain_decomp: top')

                        call mpp_domains_init(MPP_DOMAIN_TIME)

                        call print_memuse_stats('external_ic:mpp_domain_decomp: mpp_domains_init')

! call mpp_domains_set_stack_size(10000)
! call mpp_domains_set_stack_size(900000)
! call mpp_domains_set_stack_size(1500000)
                        call mpp_domains_set_stack_size(3000000)

                        select case(nregions)
                        case ( 1 )  ! Lat-Lon "cyclic"

                           select case (grid_type)
                           case (3)   ! Lat-Lon "cyclic"
                              type="Lat-Lon: cyclic"
                              ntiles = 4
                              num_contact = 8
                              if( mod(npes,ntiles) .NE. 0 ) then
                                 call mpp_error(NOTE,'TEST_MPP_DOMAINS: for Cyclic mosaic, npes should be multiple of ntiles. ' // &
                                       'No test is done for Cyclic mosaic. ' )
                                 return
                              end if
                              npes_per_tile = npes/ntiles
                              call mpp_define_layout( (/1,npx-1,1,npy-1/), npes_per_tile, layout )
                              layout = (/1,npes_per_tile/) ! force decomp only in Lat-Direction
                           case (4)   ! Cartesian, double periodic
                              type="Cartesian: double periodic"
                              ntiles = 1
                              num_contact = 2
                              npes_per_tile = npes/ntiles
                              call mpp_define_layout( (/1,npx-1,1,npy-1/), npes_per_tile, layout )
                           case (5)   ! latlon patch
                              type="Lat-Lon: patch"
                              ntiles = 1
                              num_contact = 0
                              npes_per_tile = npes/ntiles
                              call mpp_define_layout( (/1,npx-1,1,npy-1/), npes_per_tile, layout )
                           case (6)   ! latlon strip
                              type="Lat-Lon: strip"
                              ntiles = 1
                              num_contact = 1
                              npes_per_tile = npes/ntiles
                              call mpp_define_layout( (/1,npx-1,1,npy-1/), npes_per_tile, layout )
                           case (7)   ! Cartesian, channel
                              type="Cartesian: channel"
                              ntiles = 1
                              num_contact = 1
                              npes_per_tile = npes/ntiles
                              call mpp_define_layout( (/1,npx-1,1,npy-1/), npes_per_tile, layout )
                           end select

                        case ( 6 )  ! Cubed-Sphere
                           type="Cubic: cubed-sphere"
                           ntiles = 6
                           num_contact = 12
!--- cubic grid always have six tiles, so npes should be multiple of 6
                           if( mod(npes,ntiles) .NE. 0 .OR. npx-1 .NE. npy-1) then
                              call mpp_error(NOTE,'mpp_domain_decomp: for Cubic_grid mosaic, npes should be multiple of ntiles(6) ' // &
                                    'and npx-1 should equal npy-1, mpp_domain_decomp is NOT done for Cubic-grid mosaic. ' )
                              return
                           end if
                           npes_per_tile = npes/ntiles
                           call  mpp_define_layout( (/1,npx-1,1,npy-1/), npes_per_tile, layout )
                           npes_x = layout(1)
                           npes_y = layout(2)
                           if ( (npx/npes_x < ng) .or. (npy/npes_y < ng) ) then
                              write(*,310) npes_x, npes_y, npx/npes_x, npy/npes_y
310                           format('Invalid layout, NPES_X:',i4.4,'NPES_Y:',i4.4,'ncells_X:',i4.4,'ncells_Y:',i4.4)
                              call mpp_error(FATAL, 'mpp_domain_decomp: bad decomp')
                           endif
                        case default
                           call mpp_error(FATAL, 'mpp_domain_decomp: no such test: '//type)
                        end select

                        call print_memuse_stats('external_ic:mpp_domain_decomp: mpp_define_layout')

                        allocate(layout2D(2,ntiles), global_indices(4,ntiles), npes_tile(ntiles) )
                        allocate(pe_start(ntiles),pe_end(ntiles))
                        npes_tile = npes_per_tile
                        do n = 1, ntiles
                           global_indices(:,n) = (/1,npx-1,1,npy-1/)
                           layout2D(:,n)         = layout
                           pe_start(n) = (n-1)*layout(1)*layout(2)
                           pe_end(n)   = pe_start(n) + layout(1)*layout(2) -1
                        end do
                        num_alloc=max(1,num_contact)
                        allocate(tile1(num_alloc), tile2(num_alloc) )
                        allocate(istart1(num_alloc), iend1(num_alloc), jstart1(num_alloc), jend1(num_alloc) )
                        allocate(istart2(num_alloc), iend2(num_alloc), jstart2(num_alloc), jend2(num_alloc) )

                        is_symmetry = .true.

                        call print_memuse_stats('external_ic:mpp_domain_decomp: allocates 1')

                        select case(nregions)
                        case ( 1 )

                           select case (grid_type)
                           case (3)   ! Lat-Lon "cyclic"
!--- Contact line 1, between tile 1 (EAST) and tile 2 (WEST)
                              tile1(1) = 1; tile2(1) = 2
                              istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
                              istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
!--- Contact line 2, between tile 1 (SOUTH) and tile 3 (NORTH)  --- cyclic
                              tile1(2) = 1; tile2(2) = 3
                              istart1(2) = 1;  iend1(2) = nx; jstart1(2) = 1;   jend1(2) = 1
                              istart2(2) = 1;  iend2(2) = nx; jstart2(2) = ny;  jend2(2) = ny
!--- Contact line 3, between tile 1 (WEST) and tile 2 (EAST) --- cyclic
                              tile1(3) = 1; tile2(3) = 2
                              istart1(3) = 1;  iend1(3) = 1;  jstart1(3) = 1;  jend1(3) = ny
                              istart2(3) = nx; iend2(3) = nx; jstart2(3) = 1;  jend2(3) = ny
!--- Contact line 4, between tile 1 (NORTH) and tile 3 (SOUTH)
                              tile1(4) = 1; tile2(4) = 3
                              istart1(4) = 1;  iend1(4) = nx; jstart1(4) = ny;  jend1(4) = ny
                              istart2(4) = 1;  iend2(4) = nx; jstart2(4) = 1;   jend2(4) = 1
!--- Contact line 5, between tile 2 (SOUTH) and tile 4 (NORTH) --- cyclic
                              tile1(5) = 2; tile2(5) = 4
                              istart1(5) = 1;  iend1(5) = nx; jstart1(5) = 1;  jend1(5) = 1
                              istart2(5) = 1;  iend2(5) = nx; jstart2(5) = ny; jend2(5) = ny
!--- Contact line 6, between tile 2 (NORTH) and tile 4 (SOUTH)
                              tile1(6) = 2; tile2(6) = 4
                              istart1(6) = 1;  iend1(6) = nx; jstart1(6) = ny;  jend1(6) = ny
                              istart2(6) = 1;  iend2(6) = nx; jstart2(6) = 1;   jend2(6) = 1
!--- Contact line 7, between tile 3 (EAST) and tile 4 (WEST)
                              tile1(7) = 3; tile2(7) = 4
                              istart1(7) = nx; iend1(7) = nx; jstart1(7) = 1;  jend1(7) = ny
                              istart2(7) = 1;  iend2(7) = 1;  jstart2(7) = 1;  jend2(7) = ny
!--- Contact line 8, between tile 3 (WEST) and tile 4 (EAST) --- cyclic
                              tile1(8) = 3; tile2(8) = 4
                              istart1(8) = 1;  iend1(8) = 1;  jstart1(8) = 1;  jend1(8) = ny
                              istart2(8) = nx; iend2(8) = nx; jstart2(8) = 1;  jend2(8) = ny
                              is_symmetry = .false.
                           case (4)   ! Cartesian, double periodic
!--- Contact line 1, between tile 1 (EAST) and tile 1 (WEST)
                              tile1(1) = 1; tile2(1) = 1
                              istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
                              istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
!--- Contact line 2, between tile 1 (SOUTH) and tile 1 (NORTH)  --- cyclic
                              tile1(2) = 1; tile2(2) = 1
                              istart1(2) = 1;  iend1(2) = nx; jstart1(2) = 1;   jend1(2) = 1
                              istart2(2) = 1;  iend2(2) = nx; jstart2(2) = ny;  jend2(2) = ny
                           case (5)   ! latlon patch

                           case (6)   !latlon strip
!--- Contact line 1, between tile 1 (EAST) and tile 1 (WEST)
                              tile1(1) = 1; tile2(1) = 1
                              istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
                              istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
                           case (7)   ! Cartesian, channel
!--- Contact line 1, between tile 1 (EAST) and tile 1 (WEST)
                              tile1(1) = 1; tile2(1) = 1
                              istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
                              istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
                           end select

                        case ( 6 )  ! Cubed-Sphere
!--- Contact line 1, between tile 1 (EAST) and tile 2 (WEST)
                           tile1(1) = 1; tile2(1) = 2
                           istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
                           istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
!--- Contact line 2, between tile 1 (NORTH) and tile 3 (WEST)
                           tile1(2) = 1; tile2(2) = 3
                           istart1(2) = 1;  iend1(2) = nx; jstart1(2) = ny; jend1(2) = ny
                           istart2(2) = 1;  iend2(2) = 1;  jstart2(2) = ny; jend2(2) = 1
!--- Contact line 3, between tile 1 (WEST) and tile 5 (NORTH)
                           tile1(3) = 1; tile2(3) = 5
                           istart1(3) = 1;  iend1(3) = 1;  jstart1(3) = 1;  jend1(3) = ny
                           istart2(3) = nx; iend2(3) = 1;  jstart2(3) = ny; jend2(3) = ny
!--- Contact line 4, between tile 1 (SOUTH) and tile 6 (NORTH)
                           tile1(4) = 1; tile2(4) = 6
                           istart1(4) = 1;  iend1(4) = nx; jstart1(4) = 1;  jend1(4) = 1
                           istart2(4) = 1;  iend2(4) = nx; jstart2(4) = ny; jend2(4) = ny
!--- Contact line 5, between tile 2 (NORTH) and tile 3 (SOUTH)
                           tile1(5) = 2; tile2(5) = 3
                           istart1(5) = 1;  iend1(5) = nx; jstart1(5) = ny; jend1(5) = ny
                           istart2(5) = 1;  iend2(5) = nx; jstart2(5) = 1;  jend2(5) = 1
!--- Contact line 6, between tile 2 (EAST) and tile 4 (SOUTH)
                           tile1(6) = 2; tile2(6) = 4
                           istart1(6) = nx; iend1(6) = nx; jstart1(6) = 1;  jend1(6) = ny
                           istart2(6) = nx; iend2(6) = 1;  jstart2(6) = 1;  jend2(6) = 1
!--- Contact line 7, between tile 2 (SOUTH) and tile 6 (EAST)
                           tile1(7) = 2; tile2(7) = 6
                           istart1(7) = 1;  iend1(7) = nx; jstart1(7) = 1;  jend1(7) = 1
                           istart2(7) = nx; iend2(7) = nx; jstart2(7) = ny; jend2(7) = 1
!--- Contact line 8, between tile 3 (EAST) and tile 4 (WEST)
                           tile1(8) = 3; tile2(8) = 4
                           istart1(8) = nx; iend1(8) = nx; jstart1(8) = 1;  jend1(8) = ny
                           istart2(8) = 1;  iend2(8) = 1;  jstart2(8) = 1;  jend2(8) = ny
!--- Contact line 9, between tile 3 (NORTH) and tile 5 (WEST)
                           tile1(9) = 3; tile2(9) = 5
                           istart1(9) = 1;  iend1(9) = nx; jstart1(9) = ny; jend1(9) = ny
                           istart2(9) = 1;  iend2(9) = 1;  jstart2(9) = ny; jend2(9) = 1
!--- Contact line 10, between tile 4 (NORTH) and tile 5 (SOUTH)
                           tile1(10) = 4; tile2(10) = 5
                           istart1(10) = 1;  iend1(10) = nx; jstart1(10) = ny; jend1(10) = ny
                           istart2(10) = 1;  iend2(10) = nx; jstart2(10) = 1;  jend2(10) = 1
!--- Contact line 11, between tile 4 (EAST) and tile 6 (SOUTH)
                           tile1(11) = 4; tile2(11) = 6
                           istart1(11) = nx; iend1(11) = nx; jstart1(11) = 1;  jend1(11) = ny
                           istart2(11) = nx; iend2(11) = 1;  jstart2(11) = 1;  jend2(11) = 1
!--- Contact line 12, between tile 5 (EAST) and tile 6 (WEST)
                           tile1(12) = 5; tile2(12) = 6
                           istart1(12) = nx; iend1(12) = nx; jstart1(12) = 1;  jend1(12) = ny
                           istart2(12) = 1;  iend2(12) = 1;  jstart2(12) = 1;  jend2(12) = ny
                        end select

                        call mpp_define_mosaic(global_indices, layout2D, domain, ntiles, num_contact, tile1, tile2, &
                              istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,      &
                              pe_start=pe_start, pe_end=pe_end, symmetry=is_symmetry,              &
                              shalo = ng, nhalo = ng, whalo = ng, ehalo = ng, name = type)
                        call print_memuse_stats('external_ic:mpp_domain_decomp: mpp_define_mosaic')

                        deallocate(pe_start,pe_end)

!--- find the tile number
                        tile = mpp_pe()/npes_per_tile+1
                        call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
                        call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )

                        call print_memuse_stats('external_ic:mpp_domain_decomp: mpp_get domains')

                     end subroutine mpp_domain_decomp
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------
                        
                  end module fv_regrid_c2c

