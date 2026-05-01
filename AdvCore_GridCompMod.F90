#include "MAPL.h"

!BOP
!MODULE: AdvCore_GridCompMod
!DESCRIPTION:
!    This a MAPL component that can be used in
!    either with offline or online applications to advect an arbitrary set
!    of constituents.
!
! \paragraph{Scientific Description:}
!
!   The advection scheme used is that from the FVdycore grid-point
!   dynamical core.  It runs on a sphere and uses finite-volume
!   discretization techniques. The advection is time split into a
!   horizontal phase that is assumed to be vertically Lagrangian and a
!   vertical remap phase. A complete description of the core from
!   which this component is taken may be found in:
!
!   \begin{quote}
!   Lin, S.-J. 2004, A vertically Lagrangian Finite-Volume Dynamical
!   Core for Global Models. {\em Mon. Wea. Rev.}, {\bf 132}, 2293-2307.
!   \end{quote}
!
!  \paragraph{Code Implementation:}
!
!    It code uses the MAPL (http://MAPLCode.org/maplwiki/) to
!    encapsulate the FV advection scheme as an ESMF gridded component
!    using the ESMF paradigm of initialize, run and finalize methods,
!    and their SetServices routine. As in all ESMF codes, only
!    SetServices is public and the interface consists of of a Clock
!    and Import and Export states.  The import state includes a
!    specialized description of the motion field in terms of C-grid
!    winds and mass fluxes. These are assumed to have been accumulated
!    over the time interval specified in the resource file. The
!    default of this interval is 1800 seconds. The layer pressure
!    thicknesses in the import state are assumed to be the
!    instantaneous values valid at the beginning of this interval.  If
!    these thicknesses are friendly they will be updated to values
!    valid at the end of the interval, consistent with the given
!    motion field.  Mixing ratios of the constituents to be advected
!    are placed ESMF Fields within an ESMF Bundle in the Import
!    state. Each Field in the Bundle is tested for ``Friendliness'' to
!    advection; if friendly it is advected and its values updated.
!
!    Currently no Export capability is implemented.

!INTERFACE:
module AdvCore_GridCompMod

   !USES:
   use ESMF
   use MAPL
   use mapl3g_GridGet, only: GridGet
   use m_set_eta,       only: set_eta
   use mpp_mod,         only: mpp_pe, mpp_root_pe
   use fv_arrays_mod,   only: fv_atmos_type, FVPRC, REAL4, REAL8
   use fms_mod,         only: fms_init, set_domain, nullify_domain
   use fv_control_mod,  only: fv_init1, fv_init2, fv_end
   use fv_tracer2d_mod, only: offline_tracer_advection
   use fv_mp_mod,       only: is,ie, js,je, is_master, tile
   use fv_grid_utils_mod, only: g_sum_r8

   use fv_diagnostics_mod, only: prt_maxmin, prt_minmax

   USE FV_StateMod,     only: AdvCoreTracers => T_TRACERS
   USE FV_StateMod,     only: FV_Atm

   implicit none
   private

   integer :: QSPLIT
   integer :: nx, ny
   integer :: npes_x, npes_y
   logical :: FV3_DynCoreIsRunning=.false.
   integer :: AdvCore_Advection=1
   logical :: rpt_mass=.false.
   logical :: DEBUG_ADV = .false.
   real(FVPRC) :: dt

   integer, parameter :: ntiles_per_pe = 1

   ! Tracer I/O History stuff
   integer, parameter :: ntracers=38
   integer :: ntracer
   character(len=ESMF_MAXSTR) :: myTracer
   character(len=ESMF_MAXSTR) :: tMassStr
   logical, save :: firstRun=.true.

   !PUBLIC MEMBER FUNCTIONS:

   public SetServices
   logical, allocatable, save :: grids_on_my_pe(:)
   !EOP

contains

   !BOP
   !IROUTINE: SetServices - Externally visible registration routine
   !INTERFACE:
   subroutine SetServices(gc, rc)
      !ARGUMENTS:
      type(ESMF_GridComp), intent(inout) :: gc
      integer, optional,   intent(  out) :: rc

      !DESCRIPTION:
      ! User-supplied setservices routine.
      ! The register routine sets the subroutines to be called
      ! as the init, run, and finalize routines.  Note that those are
      ! private to the module.
      !EOP

      character(len=ESMF_MAXSTR) :: IAm
      integer :: status
      character(len=ESMF_MAXSTR) :: comp_name
      character(len=ESMF_MAXSTR) :: dycore
      type(ESMF_VM) :: VM
      integer :: comm, ndt
      integer :: p_split=1

      ! Get my name and set-up traceback handle
      call ESMF_GridCompGet(gc, name=comp_name, vm=vm, _RC)
      Iam = trim(comp_name) // 'SetServices'

#include "AdvCore_Import___.h"
      call MAPL_AddImportSpec(gc, &
           short_name='TRADV', &
           long_name='advected_quantities', &
           units='unknown', &
           datatype=MAPL_BundleItem, _RC)

#include "AdvCore_Export___.h"
      ! 3D Tracers
      do ntracer=1,ntracers
         write(myTracer, "('TEST_TRACER',i5.5)") ntracer-1
         call MAPL_AddExportSpec(gc, &
              short_name=trim(myTracer), &
              long_name=trim(myTracer), &
              units='1', &
              dims=MAPL_DimsHorzVert, &
              vlocation=MAPL_VLocationCenter, _RC)
      enddo

      ! Set the Profiling timers
      ! Register methods with MAPL
      call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_INITIALIZE, Initialize, _RC)
      call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_RUN, Run, _RC)
      call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_FINALIZE, Finalize, _RC)

      ! Check if AdvCore is running without FV3_DynCoreIsRunning, if yes then setup the MAPL Grid
      call MAPL_GridCompGetResource(gc, 'DYCORE:', value=dycore, default="", _RC)

      if(adjustl(DYCORE)=="FV3") then
         FV3_DynCoreIsRunning = .true.
         AdvCore_Advection = 0
      endif
      if(adjustl(DYCORE)=="FV3+ADV") then
         FV3_DynCoreIsRunning = .true.
         AdvCore_Advection = 1
      endif

      call MAPL_GridCompGetResource(gc, 'AdvCore_Advection:', value=AdvCore_Advection, default=AdvCore_Advection, _RC)
      call MAPL_GridCompGetResource(gc, 'DEBUG_ADV:', value=DEBUG_ADV, default=.FALSE., _RC)

      ! Start up FMS/MPP
      !-------------------------------------------
      call ESMF_VMGet(vm, mpiCommunicator=comm, _RC)
      call fms_init(comm)

      if (.NOT. FV3_DynCoreIsRunning) then
          ! Make sure FV3 is setup
          call fv_init1(FV_Atm, dt, grids_on_my_pe, p_split)
          ! Get Domain decomposition
         call MAPL_GridCompGetResource(gc, 'NX:', value=nx, default=0, _RC)
         FV_Atm(1)%layout(1) = nx
         call MAPL_GridCompGetResource(gc, 'NY:', value=ny, default=0, _RC)
         if (FV_Atm(1)%flagstruct%grid_type == 4) then
            FV_Atm(1)%layout(2) = ny
         else
            FV_Atm(1)%layout(2) = ny / 6
         end if
         ! Get Resolution Information
         ! FV grid dimensions setup from MAPL
         call MAPL_GridCompGetResource(gc, 'IM:', value=FV_Atm(1)%flagstruct%npx, default=32, _RC)
         call MAPL_GridCompGetResource(gc, 'JM:', value=FV_Atm(1)%flagstruct%npy, default=192, _RC)
         call MAPL_GridCompGetResource(gc, 'LM:', value=FV_Atm(1)%flagstruct%npz, default=72, _RC)

         ! FV likes npx;npy in terms of cell vertices
         if (FV_Atm(1)%flagstruct%npy == 6*FV_Atm(1)%flagstruct%npx) then
            FV_Atm(1)%flagstruct%ntiles = 6
            FV_Atm(1)%flagstruct%npy    = FV_Atm(1)%flagstruct%npx+1
            FV_Atm(1)%flagstruct%npx    = FV_Atm(1)%flagstruct%npx+1
         else
            FV_Atm(1)%flagstruct%ntiles = 1
            FV_Atm(1)%flagstruct%npy    = FV_Atm(1)%flagstruct%npy+1
            FV_Atm(1)%flagstruct%npx    = FV_Atm(1)%flagstruct%npx+1
         endif
      endif

      call MAPL_GridCompGetResource(gc, 'RUN_DT:', value=ndt, default=0, _RC)
      DT = ndt

      call MAPL_GridCompGetResource(gc, 'ADV_CORE_REPORT_TRACER_MASS:', value=rpt_mass, default=rpt_mass, _RC)
      call MAPL_GridCompGetResource(gc, 'ADV_QSPLIT:', value=QSPLIT, default=0, _RC)

      ! Start up FV if AdvCore is running without FV3_DynCoreIsRunning
      if (.NOT. FV3_DynCoreIsRunning) then
         call fv_init2(FV_Atm, dt, grids_on_my_pe, p_split)
      end if

      ! Ending with a Generic SetServices call is a MAPL requirement
      _RETURN(_SUCCESS)
   end subroutine SetServices

   !BOP
   !IROUTINE: Initialize - initialization routine
   !INTERFACE:
   subroutine Initialize(gc, import, export, clock, rc)
      !INPUT/OUTPUT PARAMETERS:
      type(ESMF_GridComp), intent(inout) :: gc  ! Gridded component
      type(ESMF_State), intent(inout) :: import ! Import state
      type(ESMF_State), intent(inout) :: export ! Export state
      type(ESMF_Clock), intent(inout) :: clock  ! The clock
      !OUTPUT PARAMETERS:
      integer, optional, intent(out) :: rc      ! Error code

      !DESCRIPTION:
      ! This initialization routine creates the import and export states,
      ! as well as the internal state, which is attached to the component.
      ! It also determines the distribution (and therefore the grid)
      ! and performs allocations of persistent data,
      !EOP

      !BOC
      character(len=ESMF_MAXSTR) :: IAm, comp_name
      type(ESMF_Config) :: cf
      type(ESMF_VM) :: vm
      type(ESMF_Grid) :: grid
      real, pointer :: temp2d(:,:)
      logical :: gridCreated
      integer :: is, ie, js, je, status

      ! Get the target components name and set-up traceback handle.
      Iam = "Initialize"
      call ESMF_GridCompGet(gc, name=comp_name, config=cf, vm=vm, _RC)
      Iam = trim(comp_name) // trim(Iam)

      ! Retrieve the pointer to the state
      gridCreated=.false.
      call ESMF_GridCompGet(gc, grid=grid, rc=status)
      if (status == ESMF_SUCCESS) then
         call ESMF_GridValidate(grid, rc=status)
         if (status==ESMF_SUCCESS) gridCreated = .true.
      end if

      if (.not. gridCreated) call MAPL_GridCreate(gc, _RC)

      ! Compute Grid-Cell Area
      if (.NOT. FV3_DynCoreIsRunning) then
         is = FV_Atm(1)%bd%isc
         ie = FV_Atm(1)%bd%iec
         js = FV_Atm(1)%bd%jsc
         je = FV_Atm(1)%bd%jec
         call MAPL_GetPointer(export, temp2d, 'AREA', ALLOC=.TRUE., _RC)
         temp2d = FV_Atm(1)%gridstruct%area(is:ie, js:je)
      endif

      _RETURN(_SUCCESS)
   end subroutine Initialize


!BOP
!
! !IROUTINE: Run - run routine
!
! !INTERFACE:
!
      subroutine Run(gc, import, export, clock, RC)
!
! !INPUT/OUTPUT PARAMETERS:
      type(ESMF_GridComp), intent(inout) :: gc     ! Gridded component
      type(ESMF_State),    intent(inout) :: import ! Import state
      type(ESMF_State),    intent(inout) :: export ! Export state
      type(ESMF_Clock),    intent(inout) :: clock  ! The clock
!
! !OUTPUT PARAMETERS:
      integer, optional,   intent(  out) :: RC     ! Error code
!
! !DESCRIPTION:
!
! The Run method advanced the advection one long time step, as
! specified in the configuration.  This may be broken down int a
! number of internal, small steps, also configurable.
!
!EOP
!=============================================================================
!BOC
! !LOCAL VARIABLES:
      character(len=ESMF_MAXSTR)    :: IAm
      integer                       :: STATUS
      character(len=ESMF_MAXSTR)    :: COMP_NAME
      type (ESMF_Grid)              :: ESMFGRID

! Imports
      REAL(REAL8), POINTER, DIMENSION(:,:,:)   :: iCX
      REAL(REAL8), POINTER, DIMENSION(:,:,:)   :: iCY
      REAL(REAL8), POINTER, DIMENSION(:,:,:)   :: iMFX
      REAL(REAL8), POINTER, DIMENSION(:,:,:)   :: iMFY
      REAL(REAL8), POINTER, DIMENSION(:,:,:)   :: iPLE0
      REAL(REAL8), POINTER, DIMENSION(:,:,:)   :: iPLE1

! Locals
      REAL(FVPRC), POINTER, DIMENSION(:,:,:)   :: CX
      REAL(FVPRC), POINTER, DIMENSION(:,:,:)   :: CY
      REAL(FVPRC), POINTER, DIMENSION(:,:,:)   :: MFX
      REAL(FVPRC), POINTER, DIMENSION(:,:,:)   :: MFY
      REAL(FVPRC), POINTER, DIMENSION(:,:,:)   :: PLE0
      REAL(FVPRC), POINTER, DIMENSION(:,:,:)   :: PLE1
      REAL(FVPRC), POINTER, DIMENSION(:,:,:,:) :: TRACERS
      REAL(REAL8), allocatable :: TMASS0(:)
      REAL(REAL8), allocatable :: TMASS1(:)
      TYPE(AdvCoreTracers), POINTER :: advTracers(:)
      type(ESMF_FieldBundle) :: TRADV
      type(ESMF_Field)       :: field
      type(ESMF_Array)       :: array
      INTEGER :: IM, JM, LM, N, NQ, LS
! Temporaries for exports/tracers
      REAL, POINTER :: temp3D(:,:,:)
      real(REAL4),        pointer     :: tracer_r4 (:,:,:)
      real(REAL8),        pointer     :: tracer_r8 (:,:,:)
      character(len=ESMF_MAXSTR)    :: fieldName
      type(ESMF_TypeKind_Flag)      :: kind
      character(len=ESMF_MAXSTR)    :: STRING
! Excluding tracers
      type(ESMF_FieldBundle), save        :: bundleAdv
      type (ESMF_Config)                  :: CF
      logical                             :: adjustTracers
      type(ESMF_Alarm)                    :: predictorAlarm
      type(ESMF_Grid)                     :: bgrid
      integer, save                       :: nq_saved = 0
      integer                             :: i,j,k
      integer                             :: nqt
      logical                             :: tend
      logical                             :: exclude
      character(len=ESMF_MAXSTR)          :: tmpstring
      character(len=ESMF_MAXSTR)          :: adjustTracerMode
      character(len=ESMF_MAXSTR), allocatable :: xlist(:)
      character(len=ESMF_MAXSTR), allocatable :: biggerlist(:)
      integer, parameter                  :: XLIST_MAX = 60

      real(FVPRC), allocatable           :: DEBUG_ARRAY(:,:,:)
      real(FVPRC) :: fac1    = 1.0

! Get my name and set-up traceback handle
! ---------------------------------------

      Iam = 'Run'
      call ESMF_GridCompGet( gc, name=COMP_NAME, CONFIG=CF, grid=ESMFGRID, RC=STATUS )
      VERIFY_(STATUS)
      Iam = trim(COMP_NAME) // Iam

      if (AdvCore_Advection>0) then

! Get parameters from generic state.
!-----------------------------------
      call MAPL_GridCompGet(gc, grid=ESMFGRID, num_levels=LM, _RC)
      call GridGet(ESMFGRID, im=IM, jm=JM, _RC)

      CALL MAPL_GetPointer(import, iPLE0, 'PLE0', ALLOC = .TRUE., RC=STATUS)
      VERIFY_(STATUS)
      CALL MAPL_GetPointer(import, iPLE1, 'PLE1', ALLOC = .TRUE., RC=STATUS)
      VERIFY_(STATUS)
      CALL MAPL_GetPointer(import, iMFX,   'MFX', ALLOC = .TRUE., RC=STATUS)
      VERIFY_(STATUS)
      CALL MAPL_GetPointer(import, iMFY,   'MFY', ALLOC = .TRUE., RC=STATUS)
      VERIFY_(STATUS)
      CALL MAPL_GetPointer(import, iCX,     'CX', ALLOC = .TRUE., RC=STATUS)
      VERIFY_(STATUS)
      CALL MAPL_GetPointer(import, iCY,     'CY', ALLOC = .TRUE., RC=STATUS)
      VERIFY_(STATUS)

      ALLOCATE( PLE0(IM,JM,LM+1) )
      ALLOCATE( PLE1(IM,JM,LM+1) )
      ALLOCATE(  MFX(IM,JM,LM  ) )
      ALLOCATE(  MFY(IM,JM,LM  ) )
      ALLOCATE(   CX(IM,JM,LM  ) )
      ALLOCATE(   CY(IM,JM,LM  ) )

      PLE0 = iPLE0
      PLE1 = iPLE1
       MFX = iMFX
       MFY = iMFY
        CX = iCX
        CY = iCY

      ! The quantities to be advected come as friendlies in a bundle
      !  in the import state.
      !--------------------------------------------------------------

      call ESMF_StateGet(import, "TRADV", TRADV, rc=STATUS)
      VERIFY_(STATUS)

      !-------------------------------------------------------------------
      ! ALT: this section attempts to limit the amount of advected tracers
      !-------------------------------------------------------------------
      adjustTracers = .false.
      call MAPL_GridCompGetResource ( gc, 'EXCLUDE_ADVECTION_TRACERS:', &
           value=adjustTracerMode, &
           default='ALWAYS', rc=status )
      VERIFY_(STATUS)
      if (adjustTracerMode == 'ALWAYS') then
         adjustTracers = .true.
      else if (adjustTracerMode == 'PREDICTOR') then
         !get PredictorAlarm from clock
         call ESMF_ClockGetAlarm(clock, alarmName='PredictorAlarm', &
              alarm=PredictorAlarm, rc=status)
         if (status == ESMF_SUCCESS) then
            !check if ringing
            if (ESMF_AlarmIsRinging(predictorAlarm)) then
               adjustTracers = .true.
            end if
         end if
      else if (adjustTracerMode == 'NO') then
         ! Proceed without warning
         adjustTracers = .false.
      else
         call WRITE_PARALLEL('Invalid option, ignored')
         adjustTracers = .false.
      end if
      if (adjustTracers) then
         if (firstRun) then
            ! get the list of excluded tracers from resource
            allocate(xlist(XLIST_MAX), stat=status)
            VERIFY_(STATUS)
            n = 0
            call ESMF_ConfigFindLabel ( CF,'EXCLUDE_ADVECTION_TRACERS_LIST:',rc=STATUS )
            if(STATUS==ESMF_SUCCESS) then

               tend  = .false.
               do while (.not.tend)
                  call ESMF_ConfigGetAttribute (CF,value=tmpstring,default='',rc=STATUS) !ALT: we don't check return status!!!
                  if (tmpstring /= '')  then
                     n = n + 1
                     if (n > size(xlist)) then
                        allocate( biggerlist(2*n), stat=status )
                        VERIFY_(STATUS)
                        biggerlist(1:n-1)=xlist
                        call move_alloc(from=biggerlist, to=xlist)
                     end if
                     xlist(n) = tmpstring
                  end if
                  call ESMF_ConfigNextLine(CF,tableEnd=tend,rc=STATUS )
                  VERIFY_(STATUS)
               enddo
            end if

            ! Count the number of tracers
            !---------------------
            call ESMF_FieldBundleGet(TRADV, grid=bgrid,fieldCount=nqt,  RC=STATUS)
            VERIFY_(STATUS)
            BundleAdv = ESMF_FieldBundleCreate ( name='xTRADV', rc=STATUS )
            VERIFY_(STATUS)
            call ESMF_FieldBundleSet ( BundleAdv, grid=bgrid, rc=STATUS )
            VERIFY_(STATUS)
            !loop over NQ in TRADV
            do i = 1, nqt
               !get field from TRADV and its name
               call ESMF_FieldBundleGet(TRADV, fieldIndex=i, field=field, rc=status)
               VERIFY_(STATUS)
               call ESMF_FieldGet(FIELD, name=fieldname, RC=STATUS)
               VERIFY_(STATUS)
               !exclude everything that is not cloud/water species
               if ( (FV3_DynCoreIsRunning      ) .and. &
                   ( (TRIM(fieldname) == 'Q'       ) .or. &
                     (TRIM(fieldname) == 'QLCN'    ) .or. &
                     (TRIM(fieldname) == 'QLLS'    ) .or. &
                     (TRIM(fieldname) == 'QICN'    ) .or. &
                     (TRIM(fieldname) == 'QILS'    ) .or. &
                     (TRIM(fieldname) == 'CLCN'    ) .or. &
                     (TRIM(fieldname) == 'CLLS'    ) .or. &
                     (TRIM(fieldname) == 'NCPL'    ) .or. &
                     (TRIM(fieldname) == 'NCPI'    ) .or. &
                     (TRIM(fieldname) == 'QRAIN'   ) .or. &
                     (TRIM(fieldname) == 'QSNOW'   ) .or. &
                     (TRIM(fieldname) == 'QGRAUPEL') ) ) then
                   ! write(STRING,'(A,A)') "ADV is excluding ", TRIM(fieldname)
                   ! call WRITE_PARALLEL( trim(STRING)   )
                     n = n + 1
                     if (n > size(xlist)) then
                        allocate( biggerlist(2*n), stat=status )
                        VERIFY_(STATUS)
                        biggerlist(1:n-1)=xlist
                        call move_alloc(from=biggerlist, to=xlist)
                     end if
                     xlist(n) = TRIM(fieldname)
               end if
               !loop over exclude_list
               exclude = .false.
               do j = 1, n
                  if (fieldname == xlist(j)) then
                     exclude = .true.
                     exit
                  end if
               end do
               if (.not. exclude) then
                  call MAPL_FieldBundleAdd(BundleAdv, FIELD, RC=STATUS)
                  VERIFY_(STATUS)
               end if
            end do

            if (allocated(xlist)) then
               deallocate(xlist)
            end if

            if (allocated(biggerlist)) then
               deallocate(biggerlist)
            end if

            firstRun=.false.
         end if ! firstRun
         TRADV = bundleAdv
      end if ! adjustTracers

      call ESMF_FieldBundleGet(TRADV, fieldCount=NQ,    rc=STATUS)
      VERIFY_(STATUS)

      if (NQ > 0) then
         ! We allocate a list of tracers big enough to hold all items in the bundle
         !-------------------------------------------------------------------------
         ALLOCATE( TRACERS(IM,JM,LM,NQ),stat=STATUS )
         VERIFY_(STATUS)
         ALLOCATE( advTracers(NQ),stat=STATUS )
         VERIFY_(STATUS)

         if (NQ /= NQ_SAVED) then
            write(STRING,'(A,I5,A)') "AdvCore is Advecting the following ", nq, " tracers:"
            call WRITE_PARALLEL( trim(STRING)   )
         end if

         ! Go through the bundle copying the friendlies into the tracer list.
         !-------------------------------------------------------------------------
         do N=1,NQ
            call ESMF_FieldBundleGet (TRADV, fieldIndex=N, field=FIELD, RC=STATUS)
            VERIFY_(STATUS)
            call ESMF_FieldGet  (field, array=array, name=fieldName, RC=STATUS)
            VERIFY_(STATUS)
            call ESMF_ArrayGet(array,typekind=kind, rc=status )
            VERIFY_(STATUS)
            advTracers(N)%is_r4 = (kind == ESMF_TYPEKIND_R4)   ! Is real*4?
            advTracers(N)%tName = fieldName

            if (NQ /= NQ_SAVED) call WRITE_PARALLEL( trim('--'//fieldName) )

            if (advTracers(N)%is_r4) then
               call ESMF_ArrayGet(array,farrayptr=tracer_r4, rc=status )
               VERIFY_(STATUS)
               advTracers(N)%content_r4 => tracer_r4
               TRACERS(:,:,:,N) = advTracers(N)%content_r4
            else
               call ESMF_ArrayGet(array,farrayptr=tracer_r8, rc=status )
               VERIFY_(STATUS)
               advTracers(N)%content => tracer_r8
               TRACERS(:,:,:,N) = advTracers(N)%content
            end if
         end do

         if (NQ /= NQ_SAVED) then
            NQ_SAVED = NQ
         end if

         ! Get Tracer Mass before advection
         !---------------------------------
         if (rpt_mass) then
         allocate( TMASS0(NQ) )
         call global_integral(TMASS0, TRACERS, PLE0, IM,JM,LM, NQ)
         endif

         ! Run FV3 advection
         !------------------
         call offline_tracer_advection(TRACERS, PLE0, PLE1, MFX, MFY, CX, CY, &
                                       FV_Atm(1)%gridstruct, FV_Atm(1)%flagstruct, FV_Atm(1)%bd, &
                                       FV_Atm(1)%domain, FV_Atm(1)%npx, FV_Atm(1)%npy, FV_Atm(1)%npz,   &
                                       NQ, dt, QSPLIT)

         ! Get Tracer Mass after advection
         !--------------------------------
         if (rpt_mass) then
         allocate( TMASS1(NQ) )
         call global_integral(TMASS1, TRACERS, PLE1, IM,JM,LM, NQ)
         endif

         ! Conserve Specific Mass of Constituents Keeping Mixing_Ratio Constant WRT_Dry_Air
         ! --------------------------------------------------------------------------------
         if (rpt_mass) then
         do N=1,NQ
            if (TMASS1(N) > 0.0) then
            if (ABS((TMASS0(N)-TMASS1(N))/TMASS1(N)) >= epsilon(1.0_REAL4)) then
              if (is_master()) write(6,125) trim(advTracers(N)%tName), (TMASS1(N)-TMASS0(N))/TMASS0(N)
            !!TRACERS(:,:,:,N) = TRACERS(:,:,:,N) * TMASS0(N)/TMASS1(N)
            end if
            125 format('Mass Conservation Adjustment in AdvCore:'2x,A,2x,g21.14)
            end if
         end do
         deallocate( TMASS0 )
         deallocate( TMASS1 )
         endif

         ! Go through the bundle copying tracers back to the bundle.
         !-------------------------------------------------------------------------
         do N=1,NQ
            if (advTracers(N)%is_r4) then
               advTracers(N)%content_r4 = TRACERS(:,:,:,N)
            else
               advTracers(N)%content    = TRACERS(:,:,:,N)
            end if

            !-----------------------------------------------
            !--> Fill Export States
            !--> This section is used for diagnostics only.
            !--> It has no effect on CTM experiments.
            !-----------------------------------------------
            if (N<=min(ntracers,NQ)) then
               write(myTracer, "('TEST_TRACER',i5.5)") N-1
               call MAPL_GetPointer(export, temp3D, TRIM(myTracer), rc=status)
               VERIFY_(STATUS)
               if (associated(temp3D)) temp3D = TRACERS(:,:,:,N)
            endif
         enddo



         ! Clean negative tracers and check
         !-------------------------------------------------------------------------
         if (DEBUG_ADV) then
           prt_minmax     = DEBUG_ADV
           if (mpp_pe()==0) print*,''
           if (mpp_pe()==0) print*,'-------------- FV3 Tracer Debug After ADV --------------'
           allocate( DEBUG_ARRAY(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,FV_Atm(1)%npz) )
         endif
         do n=1,NQ
            if (advTracers(n)%is_r4) then
               where (advTracers(n)%content_r4 < tiny(0.0))
                      advTracers(n)%content_r4 = 0.0
               end where
               if (DEBUG_ADV) DEBUG_ARRAY = advTracers(n)%content_r4
            else
               where (advTracers(n)%content < tiny(0.0))
                      advTracers(n)%content = 0.0
               end where
               if (DEBUG_ADV) DEBUG_ARRAY = advTracers(n)%content
            endif
            if (DEBUG_ADV) then
               call prt_maxmin(TRIM(advTracers(n)%tname), DEBUG_ARRAY, FV_Atm(1)%bd%isc, FV_Atm(1)%bd%iec, FV_Atm(1)%bd%jsc, FV_Atm(1)%bd%jec, 0, FV_Atm(1)%npz, fac1)
            endif
         enddo
         if (DEBUG_ADV) then
           deallocate ( DEBUG_ARRAY )
           if (mpp_pe()==0) print*,'-------------- FV3 Tracer Debug After ADV --------------'
           if (mpp_pe()==0) print*,''
           prt_minmax     = .false.
         endif

         ! Deallocate the list of tracers
         !-------------------------------------------------------------------------
         deallocate( TRACERS, stat=STATUS )
         VERIFY_(STATUS)
         deallocate( advTracers, stat=STATUS )
         VERIFY_(STATUS)

      end if ! NQ > 0

      DEALLOCATE( PLE0 )
      DEALLOCATE( PLE1 )
      DEALLOCATE(  MFX )
      DEALLOCATE(  MFY )
      DEALLOCATE(   CX )
      DEALLOCATE(   CY )

      end if ! AdvCore_Advection

      RETURN_(ESMF_SUCCESS)

      end subroutine Run
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Finalize - user supplied finalize routine
!
! !INTERFACE:
!
  subroutine Finalize(gc, import, export, clock, RC)
!
! !INPUT/OUTPUT PARAMETERS:
      type(ESMF_GridComp), intent(inout) :: gc     ! Gridded component
      type(ESMF_State),    intent(inout) :: import ! Import state
      type(ESMF_State),    intent(inout) :: export ! Export state
      type(ESMF_Clock),    intent(inout) :: clock  ! The clock
!
! !OUTPUT PARAMETERS:
      integer, optional,   intent(  out) :: RC     ! Error code
!
! !DESCRIPTION:
!    Finalize merely destroys the FVadv object that was created in Initialize
!    and releases the space for the persistent data .
!
!EOP
!=============================================================================
!BOC
! !LOCAL VARIABLES:

      character(len=ESMF_MAXSTR)    :: IAm
      integer                       :: STATUS
      character(len=ESMF_MAXSTR)    :: COMP_NAME

! Get my name and set-up traceback handle
! ---------------------------------------

      Iam = 'Finalize'
      call ESMF_GridCompGet( gc, NAME=COMP_NAME, RC=STATUS )
      VERIFY_(STATUS)
      Iam = trim(COMP_NAME) // TRIM(Iam)

      ! Clean up FV if AdvCore is running without FV3_DynCoreIsRunning
      !--------------------------------------------------
      if (.NOT. FV3_DynCoreIsRunning) then
         call fv_end(FV_Atm, grids_on_my_pe, .false.)
      endif

      RETURN_(ESMF_SUCCESS)
      end subroutine Finalize


subroutine global_integral (QG,Q,PLE,IM,JM,KM,NQ)

      real(REAL8), intent(OUT)   :: QG(NQ)
      real(FVPRC), intent(IN)    :: Q(IM,JM,KM,NQ)
      real(FVPRC), intent(IN)    :: PLE(IM,JM,KM+1)
      integer,     intent(IN)    :: IM,JM,KM,NQ
! Locals
      integer   :: k,n
      real(REAL8), allocatable ::    dp(:,:,:)
      real(REAL8), allocatable :: qsum1(:,:)
      real(REAL8) :: mass

      allocate(    dp(im,jm,km) )
      allocate( qsum1(im,jm)    )

! Compute Pressure Thickness
! --------------------------
      do k=1,KM
         dp(:,:,k) = PLE(:,:,k+1)-PLE(:,:,k)
      enddo

! Compute Global Mass
! -------------------
      qsum1(:,:) = 0.d0
      do k=1,KM
         qsum1(:,:) = qsum1(:,:) + dp(:,:,k)
      enddo
      mass = g_sum_r8(FV_Atm(1)%domain, qsum1, is,ie, js,je, FV_Atm(1)%ng, FV_Atm(1)%gridstruct%area_64, 1)

! Loop over Tracers
! -----------------
     do n=1,NQ
        qsum1(:,:) = 0.d0
        do k=1,KM
           qsum1(:,:) = qsum1(:,:) + Q(:,:,k,n)*dp(:,:,k)
        enddo
        qg(n) = g_sum_r8(FV_Atm(1)%domain, qsum1, is,ie, js,je, FV_Atm(1)%ng, FV_Atm(1)%gridstruct%area_64, 1)
        if (mass > 0.0) qg(n) = qg(n)/mass
     enddo

     deallocate( dp )
     deallocate( qsum1 )

end subroutine global_integral

end module AdvCore_GridCompMod
