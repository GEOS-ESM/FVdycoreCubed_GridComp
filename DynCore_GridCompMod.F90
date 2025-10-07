!  $id: DynCore_GridCompMod.F90,v 1.1.1.1 2007/05/29 12:26:20 atrayanov Exp $

#include "MAPL_Generic.h"

!#define SCALAR_WINDS
!#define INC_WINDS

!-----------------------------------------------------------------------
!              ESMA - Earth System Modeling Applications
!-----------------------------------------------------------------------
module FVdycoreCubed_GridComp

   !BOP
   !MODULE: FVdycoreCubed_GridComp --- Dynamical Core Grid Component

   !USES:
   use ESMF
   use mapl_ErrorHandlingMod, only: MAPL_Verify, MAPL_Assert, MAPL_Return, MAPL_VRFY

   use ESMFL_Mod, only: ESMFL_StateGetPointerToData, ESMFL_BundleGetPointerToData, MAPL_AreaMean
   use MAPL_Constants, only: MAPL_RADIUS, MAPL_CP, MAPL_PI, MAPL_PI_R8, MAPL_OMEGA, MAPL_KAPPA
   use MAPL_Constants, only: MAPL_P00, MAPL_GRAV, MAPL_RGAS, MAPL_RVAP, MAPL_CPVAP, MAPL_O3MW, MAPL_AIRMW
   use MAPL_Constants, only: MAPL_VectorField, MAPL_BundleItem
   use MAPL_Constants, only: MAPL_RestartSkip, MAPL_RestartRequired, MAPL_InitialRestart
   ! use MAPL_Constants, only: MAPL_VLocationCenter, MAPL_VLocationEdge, MAPL_VLocationNone
   use MAPL_Constants, only: MAPL_DimsHorzVert, MAPL_DimsHorzOnly, MAPL_DimsVertOnly

   use MAPL_GenericMod, only: MAPL_MetaComp, MAPL_TimerAdd, MAPL_TimerOn
   use MAPL_GenericMod, only: MAPL_TimerAdd, MAPL_TimerOn, MAPL_TimerOff
   use MAPL_GenericMod, only: MAPL_GenericFinalize
   use MAPL_GenericMod, only: MAPL_GetObjectFromGC, MAPL_GetResource, MAPL_GridCreate, MAPL_Get
   ! use MAPL_GenericMod, only: MAPL_AddImportSpec, MAPL_AddExportSpec, MAPL_AddInternalSpec
   use MAPL_AbstractRegridderMod, only: AbstractRegridder
   use MAPL_SunMod, only: MAPL_SunOrbit, MAPL_SunGetInsolation
   use MAPL_BaseMod, only: MAPL_AttributeSet, MAPL_FieldBundleAdd
   use MAPL_BaseMod, only: MAPL_UNDEF, MAPL_RemapBounds
   use MAPL_GridManagerMod, only: grid_manager
   use MAPL_RegridderManagerMod, only: regridder_manager
   use MAPL_RegridMethods, only: REGRID_METHOD_BILINEAR
   use MAPL_CFIOMod, only: MAPL_CFIORead
   use MAPL_MemUtilsMod, only: MAPL_MemUtilsWrite
   use MAPL_FieldPointerUtilities, only: MAPL_FieldDestroy
   use MAPL_MaxMinMod, only: MAPL_MaxMin
   use MAPL_CommsMod, only: MAPL_AM_I_ROOT, ArrayGather

   use FileIOSharedMod, only: WRITE_PARALLEL

   use mapl3g_generic, only: MAPL_GridCompGet, MAPL_GridCompGetResource
   use mapl3g_generic, only: MAPL_GridCompSetEntryPoint, MAPL_GridCompGetInternalState
   use mapl3g_generic, only: MAPL_GridCompAddSpec, MAPL_STATEITEM_FIELDBUNDLE
   use mapl3g_VerticalStaggerLoc, only: VERTICAL_STAGGER_NONE, VERTICAL_STAGGER_CENTER, VERTICAL_STAGGER_EDGE
   use mapl3g_Geom_API, only: MAPL_GridGet
   use mapl3g_State_API, only: MAPL_StateGetPointer
   use pflogger, only: logger_t => logger

   use m_set_eta, only: set_eta

   ! FV Specific Module
   use fv_arrays_mod,  only: REAL4, REAL8, FVPRC
   !use fv_grid_tools_mod, only: grid_type
   use FV_StateMod, only : &
        FV_Atm,                                   &
        FV_To_State, State_To_FV, DEBUG_FV_STATE, &
        DynTracers      => T_TRACERS,             &
        DynVars         => T_FVDYCORE_VARS,       &
        DynGrid         => T_FVDYCORE_GRID,       &
        DynState        => T_FVDYCORE_STATE,      &
        DynSetup        => FV_Setup,              &
        DynInit         => FV_InitState,          &
        DynRun          => FV_Run,                &
        DynFinalize     => FV_Finalize,           &
        getAllWinds     => fv_getAllWinds,        &
        getVorticity    => fv_getVorticity,       &
        getDivergence   => fv_getDivergence,      &
        fillMassFluxes  => fv_fillMassFluxes,     &
        computeMassFluxes => fv_computeMassFluxes,&
        getVerticalMassFlux => fv_getVerticalMassFlux,&
        getOmega        => fv_getOmega,           &
        getEPV          => fv_getEPV,             &
        getPKZ          => fv_getPKZ,             &
        getDELZ         => fv_getDELZ,            &
        getQ            => fv_getQ,               &
        Agrid_To_Native => INTERP_AGRID_TO_DGRID, &
        DYN_COLDSTART   => COLDSTART,             &
        DYN_CASE        => CASE_ID,               &
        DYN_DEBUG       => DEBUG,                 &
        HYDROSTATIC     => FV_HYDROSTATIC,        &
        fv_getUpdraftHelicity, DEBUG_DYN,  DEBUG_ADV, &
        ADIABATIC, SW_DYNAMICS, AdvCore_Advection
   use m_topo_remap, only: dyn_topo_remap
   use CubeGridPrototype, only: register_grid_and_regridders

   !PUBLIC MEMBER FUNCTIONS:
   implicit none
   private

   ! Include the MPI library definitons:
   include 'mpif.h'

   type(ESMF_FieldBundle), save :: BundleAdv
   integer :: NXQ = 0
   logical :: overwrite_Q = .true.
   logical :: DEBUG_TQ_ERRORS

   public :: SetServices ! Register component methods

   !DESCRIPTION: This module implements the Dynamical Core as
   !               an ESMF gridded component.

   ! \paragraph*{Overview}
   !
   !   This module contains an ESMF wrapper for a generic
   !   Dynamical Core.
   !
   ! \paragraph*{Internal State}
   !
   !  FVdycore maintains an internal state consisting of the
   !  following fields:  control variables
   !
   !   \begin{itemize}
   !     \item {\tt U}:    U winds on the native grid  (m/s)
   !     \item {\tt V}:    V winds on the native grid (m/s)
   !     \item {\tt PT}:   Dry Potential Temperature (T/PKZ)
   !     \item {\tt PE}:   Edge pressures
   !     \item {\tt Q}:    Tracers
   !     \item {\tt PKZ}:  Consistent mean for p$^\kappa$
   !     \item {\tt DZ}:   Height thickness (Non-Hydrostatic)
   !   \end{itemize}
   !
   !  as well as a GRID (to be mentioned later)
   !  and same additional run-specific variables
   !
   ! Note: {\tt PT} is not updated if the flag {\tt CONVT} is true.
   !
   ! The internal state is updated each time FVdycore is called.
   !
   ! \paragraph*{Import State}
   !
   ! The import state consists of the tendencies of the
   ! control variables plus the surface geopotential heights:
   !
   !   \begin{itemize}
   !     \item {\tt DUDT}:    U wind tendency on a A-grid (m/s)
   !     \item {\tt DVDT}:    V wind tendency on a A-grid (m/s)
   !     \item {\tt DTDT}:    Delta-pressure-weighted temperature tendency
   !     \item {\tt DPEDT}:   Edge pressure tendency
   !     \item {\tt PHIS}:    Surface Geopotential Heights
   !     \item {\tt DWDT}:    V wind tendency on a A-grid (m/s)
   !   \end{itemize}
   !
   ! These are by definition on an A-grid and have an XY
   ! domain decomposition.
   !
   ! \paragraph*{Export State}
   !
   !   The export state can provide the following variables:
   !
   !   \begin{itemize}
   !     \item {\tt U}:          U winds on a A-grid (m/s) [Lat-Lon Oriented Flow]
   !     \item {\tt V}:          V winds on a A-grid (m/s) [Lat-Lon Oriented Flow]
   !     \item {\tt U\_AGRID}:   U winds on a A-grid (m/s)
   !     \item {\tt V\_AGRID}:   V winds on a A-grid (m/s)
   !     \item {\tt U\_CGRID}:   U winds on a C-grid (m/s)
   !     \item {\tt V\_CGRID}:   V winds on a C-grid (m/s)
   !     \item {\tt U\_DGRID}:   U winds on a D-grid (m/s)
   !     \item {\tt V\_DGRID}:   V winds on a D-grid (m/s)
   !     \item {\tt T}:         Temperature (K)
   !     \item {\tt Q}:         Tracers
   !     \item {\tt TH}:        Potential Temperature (K)
   !     \item {\tt ZL}:        Mid-Layer Heights (m)
   !     \item {\tt ZLE}:       Edge Heights (m)
   !     \item {\tt PLE}:       Edge pressures (Pa)
   !     \item {\tt PLK}:       P$^\kappa$ at Mid-Layers
   !     \item {\tt PKE}:       P$^\kappa$ at Edges
   !     \item {\tt OMEGA}:     Vertical pressure velocity (pa/s)
   !     \item {\tt PV}:        Ertel's Potential Vorticity (m$^2$ / kg*s)
   !     \item {\tt DUDT}:      U-wind Tendency (m/s/s)
   !     \item {\tt DVDT}:      V-wind Tendency (m/s/s)
   !     \item {\tt DTDT}:      Mass-Weighted Temperature Tendency (Pa K/s)
   !   \end{itemize}
   !
   !   All variables are on an A-grid with points at the poles, and have an XY decomposition.
   !
   ! \paragraph*{Grids and Decompositions}
   !
   !   The current version supports only a 1D latitude-based
   !   decomposition of the domain (with OMP task-parallelism
   !   in the vertical, resulting in reasonable scalability
   !   on large PE configurations).  In the near future it will
   !   support a 2D domain decomposition, in which import and
   !   export state are decomposed in longitude and latitude,
   !   while the internal state (for the most part) is
   !   decomposed in latitude and level.  When needed,
   !   the data is redistributed (``transposed'') internally.
   !
   !   There are two fundamental ESMF grids in use;
   !   \begin{itemize}
   !     \item {GRIDXY}: longitude-latitude ESMF grid (public)
   !     \item {GRIDYZ}: A latitude-level cross-sectional
   !                     decomposition (private to this module)
   !   \end{itemize}
   !
   !   PILGRIM will be used for communication until ESMF has
   !   sufficient functionality and performance to take over
   !   the task.  The use of pilgrim requires a call to
   !   {\tt INIT\_SPMD} to set SPMD parameters, decompositions,
   !   etc.
   !
   ! \paragraph*{Required Files}
   !
   !  The following files are needed for a standard restart run:
   !
   !  \begin{itemize}
   !    \item Layout file
   !      \begin{itemize}
   !        \item {\tt nprxy\_x, nprxy\_y, npryz\_y, npryz\_z}:
   !          process dimensions in XY and YZ.
   !        \item {\tt imxy, jmxy, jmyz, kmyz}: distributions for XY and YZ
   !        \item {\tt iord, jord}: the order of the lon. and lat. algorithms
   !        \item {\tt dtime}:  The large (advection) time step
   !        \item {\tt nsplit}: the ratio between the large and small time step
   !          (possibly zero for automatic determination),
   !      \end{itemize}
   !    \item Restart file
   !      \begin{itemize}
   !        \item date in standard format yy, mm, dd, hh, mm, ss
   !        \item dimensions im, jm, km, nq
   !        \item control variables {\tt U, V, PT, PE, Q}
   !      \end{itemize}
   !    \item Topography file
   !
   !  \end{itemize}
   !
   ! \paragraph*{Future Additions}
   !
   !  \begin{itemize}
   !    \item  Conservation of energy (CONSV  == .TRUE. )
   !    \item  2D decomposition (requires transposes in the coupler)
   !    \item  Use r8 instead of r4 (currently supported in StopGap)
   !  \end{itemize}
   !
   !EOP

   ! !REVISION HISTORY:
   ! 11Jul2003  Sawyer    From Trayanov/da Silva EVAC
   ! 23Jul2003  Sawyer    First informal tiptoe-through
   ! 29Jul2003  Sawyer    Modifications based on comments from 23Jul2003
   ! 28Aug2003  Sawyer    First check-in; Internal state to D-grid
   ! 15Sep2003  Sawyer    Extensive bug fixes, revisions
   ! 24Sep2003  Sawyer    Modified names; corrected weighting of T, Q
   ! 22Oct2003  Sawyer    pmgrid removed (data now in spmd\_dyn)
   ! 25Nov2003  Sawyer    Optimization for 1D decomposition
   ! 03Dec2003  Sawyer    Switched over to specified decompositions
   ! 04Dec2003  Sawyer    Moved T_FVDYCORE_GRID to dynamics_vars
   ! 21Jan2004  Takacs    Modified Import/Export, Added Generic State, Added TOPO utility
   ! 20Sep2004  Sawyer    Revised cd_core, trac2d interfaces, refactoring
   ! 06Oct2004  Sawyer    More refactoring, removed spmd_dyn
   ! 17Feb2005  Sawyer    Added Ertel's potential vorticity to diagnostics
   ! 20Mar2005  Sawyer    Tracers are now pointers into import state
   ! 12Apr2005  Sawyer    Extensive changes to minimize tracer memory
   ! 18May2005  Sawyer    Put FVdycore_wrapper in separate file; CAM/GEOS5 merge
   ! 16Nov2005  Takacs    Added option for DCADJ, Merge with Daedalus_p5
   ! 18Jan2006  Putman    Added mass fluxes to export state
   ! 24Jul2012  Todling   Revisit intermittent replay (corrections for cubed)

   integer,  parameter :: r8           = REAL8
   integer,  parameter :: r4           = REAL4

   real(r4), parameter :: RADIUS       = MAPL_RADIUS
   real(r4), parameter :: CP           = MAPL_CP
   real(r4), parameter :: PI           = MAPL_PI_R8
   real(r4), parameter :: OMEGA        = MAPL_OMEGA
   real(r4), parameter :: KAPPA        = MAPL_KAPPA
   real(r4), parameter :: P00          = MAPL_P00
   real(r4), parameter :: GRAV         = MAPL_GRAV
   real(r4), parameter :: RGAS         = MAPL_RGAS
   real(r4), parameter :: RVAP         = MAPL_RVAP
   real(r4), parameter :: EPS          = RVAP/RGAS-1.0

   integer,  parameter :: TIME_TO_RUN  = 1
   integer,  parameter :: CHECK_MAXMIN = 2

   integer :: I, J, K  !  Default declaration for loops.

   ! Tracer I/O History stuff
   integer, parameter :: ntracers=11
   integer, parameter :: plevs(5) = [850, 700, 600, 500, 300]

   ! Wrapper for extracting internal state

   type dyn_wrap
      type(DynState), pointer :: dyn_state
   end type dyn_wrap

   interface addTracer
      module procedure addTracer_r4
      module procedure addTracer_r8
   end interface addTracer

   interface Write_Profile
      module procedure Write_Profile_R4
      module procedure Write_Profile_R8
      module procedure Write_Profile_2d_R4
      module procedure Write_Profile_2d_R8
   end interface Write_Profile

   logical :: DO_ADD_INCS = .true.

contains

   !BOP

   !IROUTINE: SetServices

   !DESCRIPTION:  SetServices registers Initialize, Run, and Finalize
   !   methods for FV. Two stages of the FV run method are registered. The
   !   first one does the dynamics calculations, and the second adds
   !   increments from external sources that appear in the Import state.
   !   SetServices also creates a private internal state in which FV
   !   keeps invariant or auxilliary state variables, as well as pointers to
   !   the true state variables. The MAPL internal state contains the
   !   true state variables and is managed by MAPL.
   !
   !  The component uses all three states (Import, Export
   !  and Internal), in addition to a Private (non-ESMF) Internal state. All
   !  three are managed by MAPL.
   !
   !  The Private Internal state contains invariant
   !  quantities defined by an FV specific routine, as well as pointers
   !  to the true state variables, kept in the MAPL Internal state.
   !  The MAPL Internal is kept at FV's real*8 precision.
   !
   !  The Import State conatins tendencies to be added in the second
   !  run stage, the geopotential at the lower boundary, and a bundle
   !  of Friendly tracers to be advected. The Import and Export states
   !  are both at the default precision.

   !INTERFACE:

   Subroutine SetServices(gc, rc)

      !ARGUMENTS:
      type(ESMF_GridComp), intent(inout) :: gc     ! gridded component
      integer, intent(out), optional :: rc     ! return code

      !DESCRIPTION: Set services (register) for the FVCAM Dynamical Core GridComp
      !EOP

      type(DynState), pointer :: dyn_internal_state
      type(DYN_wrap) :: wrap
      character(len=:), allocatable :: layout_file
      character(len=ESMF_MAXSTR) :: myTracer
      class(logger_t), pointer :: logger
      integer :: FV3_STANDALONE, ilev, itracer, status

      call MAPL_GridCompGet(gc, logger=logger, _RC)
      call logger%info("SetServices:: starting...")

      ! Wrap gridcomp's private state and store it in gc
      allocate(dyn_internal_state, _STAT)
      wrap%dyn_state => dyn_internal_state
      call ESMF_UserCompSetInternalState(gc, 'DYNstate', wrap, _RC)

#include "DynCore_Import___.h"
#include "DynCore_Export___.h"
#include "DynCore_Internal___.h"

      call MAPL_GridCompAddSpec(gc, &
           state_intent=ESMF_STATEINTENT_IMPORT, &
           short_name='TRADV', &
           standard_name='advected_quantities', &
           ! pchakrab: TODO - we shouldn't need dims and vstagger for a bundle
           dims="xyz", &
           vstagger=VERTICAL_STAGGER_NONE, &
           units='unknown', &
           itemtype=MAPL_STATEITEM_FIELDBUNDLE, _RC)

#ifdef SKIP_TRACERS
      do itracer = 1, ntracers
         do ilev = 1, size(plevs)
            write(myTracer, "('Q',i5.5,'_',i3.3)") itracer-1, plevs(ilev)
            call MAPL_AddExportSpec(gc, &
                 short_name=TRIM(myTracer), &
                 long_name =TRIM(myTracer), &
                 units='1', &
                 dims=MAPL_DimsHorzOnly, &
                 vlocation=MAPL_VLocationNone, _RC)
         enddo
         write(myTracer, "('Q',i5.5)") itracer-1
         call MAPL_AddExportSpec(gc, &
              short_name=TRIM(myTracer), &
              long_name=TRIM(myTracer), &
              units='1', &
              dims=MAPL_DimsHorzVert, &
              vlocation=MAPL_VLocationCenter, _RC)
      enddo
#endif

      ! pchakrab: TODO: DO WE STILL NEED THIS COMMENT?
      !ALT: technically the first 2 records of "old" style FV restart have
      !     6 ints: YYYY MM DD H M S
      !     5 ints: I,J,K, KS (num true pressure levels), NQ (num tracers) headers

      ! Set the Profiling timers
      call MAPL_TimerAdd(gc, name="INITIALIZE", _RC)
      call MAPL_TimerAdd(gc, name="RUN", _RC)
      call MAPL_TimerAdd(gc, name="RUN2", _RC)
      call MAPL_TimerAdd(gc, name="-DYN_INIT", _RC)
      call MAPL_TimerAdd(gc, name="--FMS_INIT", _RC)
      call MAPL_TimerAdd(gc, name="--FV_INIT", _RC)
      call MAPL_TimerAdd(gc, name="-DYN_ANA", _RC)
      call MAPL_TimerAdd(gc, name="-DYN_PROLOGUE", _RC)
      call MAPL_TimerAdd(gc, name="-DYN_CORE", _RC)
      call MAPL_TimerAdd(gc, name="-DYN_EPILOGUE", _RC)
      call MAPL_TimerAdd(gc, name="--FV_DYNAMICS", _RC)
      call MAPL_TimerAdd(gc, name="--MASS_FIX", _RC)
      call MAPL_TimerAdd(gc, name="FINALIZE", _RC)

      ! Register services for this component
      call MAPL_GridCompSetEntryPoint(gc, ESMF_Method_Initialize,  Initialize, _RC)
      call MAPL_GridCompSetEntryPoint(gc, ESMF_Method_Run, Run, phase_name="Run", _RC)
      call MAPL_GridCompSetEntryPoint(gc, ESMF_Method_Run, RunAddIncs, phase_name="RunAddIncs", _RC)
      call MAPL_GridCompSetEntryPoint(gc, ESMF_Method_Finalize, Finalize, _RC)
      !  call MAPL_GridCompSetEntryPoint(gc, ESMF_SETREADRESTART, Coldstart, _RC)

      ! Setup FMS/FV3
      call MAPL_GridCompGetResource(gc, "LAYOUT", layout_file, default="fvcore_layout.rc", _RC)
      call DynSetup(gc, layout_file, _RC)

      ! Register prototype of cubed sphere grid and associated regridders
      call register_grid_and_regridders()

      ! At this point check if FV is standalone and init the grid
      call MAPL_GridCompGetResource(gc, "FV3_STANDALONE", FV3_STANDALONE, default=0, _RC)
      if (FV3_STANDALONE /= 0) then
         ! call MAPL_GridCreate(gc, _RC)
         call MAPL_GridCompAddSpec(gc, &
              state_intent=ESMF_STATEINTENT_EXPORT, &
              short_name='TRADVEX', &
              standard_name='advected_quantities', &
              dims="xyz", &
              vstagger=VERTICAL_STAGGER_NONE, &
              units='unknown', &
              itemtype=MAPL_STATEITEM_FIELDBUNDLE, _RC)
      endif

      call MAPL_GridCompGetResource(gc, "DEBUG_DYN", DEBUG_DYN, default=.false., _RC)
      call MAPL_GridCompGetResource(gc, "DEBUG_ADV", DEBUG_ADV, default=.false., _RC)
      call MAPL_GridCompGetResource(gc, "DEBUG_TQ_ERRORS", DEBUG_TQ_ERRORS, default=.false., _RC)

      call logger%info("SetServices:: ...complete")
      _RETURN(_SUCCESS)

   end subroutine SetServices

   subroutine Initialize(gc, import, export, clock, rc)
      !ARGUMENTS:
      type(ESMF_GridComp):: gc   ! composite gridded component
      type(ESMF_State) :: import ! import state
      type(ESMF_State) :: export ! export state
      type(ESMF_Clock) :: clock  ! the clock
      integer, intent(out) :: rc ! Error code, 0 all is well, error otherwise

      type(DYN_wrap) :: wrap
      type(DynState), pointer :: state
      type(DynGrid), pointer :: DycoreGrid

      type(ESMF_Field) :: field
      type(ESMF_State) :: internal
      type(ESMF_TimeInterval) :: intv
      type(ESMF_Alarm) :: alarm
      type(ESMF_FieldBundle) :: tradv, tradvex

      character(len=:), allocatable :: layout_file, ReplayMode

      real(r4), pointer :: pref(:)
      real(r4), pointer :: ple(:,:,:)
      real(r4), pointer :: u(:,:,:), v(:,:,:), t(:,:,:)
      real(r4), pointer :: temp2d(:,:)

      real(r8), pointer ::  ak(:), bk(:)
      real(r4), pointer ::  ak4(:), bk4(:)
      real(r8), pointer ::  ud(:,:,:), vd(:,:,:)
      real(r8), pointer ::  pe(:,:,:), pt(:,:,:), pk(:,:,:)
      real(r8), allocatable ::  ur(:,:,:), vr(:,:,:)

      real :: DNS_INTERVAL
      integer :: ColdRestart=0
      integer :: ifirst, ilast, jfirst, jlast, km
      integer :: i, numTracers, fv3_standalone, status

      class(logger_t), pointer :: logger

      call MAPL_GridCompGet(gc, logger=logger, _RC)
      call logger%info("Initialize:: starting...")

      ! ! Start the timers
      ! call MAPL_TimerOn(MAPL, "TOTAL")
      ! call MAPL_TimerOn(MAPL, "INITIALIZE")

      ! Get the private state
      call ESMF_UserCompGetInternalState(gc, 'DYNstate', wrap, _RC)
      state => wrap%dyn_state

      DycoreGrid  => state%grid ! direct handle to grid

      call MAPL_GridCompGetResource(gc, "LAYOUT", layout_file, default="fvcore_layout.rc", _RC)
      call MAPL_GridCompGetResource(gc, "DO_ADD_INCS", DO_ADD_INCS, default=DO_ADD_INCS, _RC)

      ! Check for ColdStart from the configuration
      call MAPL_GridCompGetResource(gc, "COLDSTART", ColdRestart, default=0, _RC)
      if (ColdRestart /= 0 ) then
         call Coldstart(gc, import, export, clock, _RC)
      endif

      ! Set Private Internal State from Restart File
      call MAPL_GridCompGetInternalState(gc, internal, _RC)

      ! call MAPL_TimerOn(MAPL, "-DYN_INIT")
      call DynInit(state, clock, import, gc, _RC)
      ! call MAPL_TimerOff(MAPL, "-DYN_INIT")

      ! Create PLE and PREF EXPORT Coupling (Needs to be done only once per run)
      call MAPL_GetPointer(export, ak4, "AK", _RC)
      call MAPL_GetPointer(export, bk4, "BK", _RC)
      call MAPL_GetPointer(internal, ak, "AK", _RC)
      call MAPL_GetPointer(internal, bk, "BK", _RC)
      call MAPL_GetPointer(export, pref, "PREF", ALLOC=.true., _RC)
      ak4 = ak
      bk4 = bk
      pref = ak + bk * P00

      call MAPL_GetPointer(internal, ud, "U", _RC)
      call MAPL_GetPointer(internal, vd, "V", _RC)
      call MAPL_GetPointer(internal, pe, "PE", _RC)
      call MAPL_GetPointer(internal, pt, "PT", _RC)
      call MAPL_GetPointer(internal, pk, "PKZ", _RC)

      call MAPL_GetPointer(export, ple, "PLE", ALLOC=.true., _RC)
      call MAPL_GetPointer(export, u, "U", ALLOC=.true., _RC)
      call MAPL_GetPointer(export, v, "V", ALLOC=.true., _RC)
      call MAPL_GetPointer(export, t, "T", ALLOC=.true., _RC)

      ! Create A-Grid Winds
      ifirst = state%grid%is
      ilast  = state%grid%ie
      jfirst = state%grid%js
      jlast  = state%grid%je
      km     = state%grid%npz

      allocate(ur(ifirst:ilast,jfirst:jlast,km))
      allocate(vr(ifirst:ilast,jfirst:jlast,km))
      call getAllWinds(ud, vd, ur=ur, vr=vr)
      u = ur
      v = vr
      t = pt*pk
      ple = pe
      deallocate(ur, vr)

      ! Fill Grid-Cell Area Delta-X/Y
      call MAPL_GetPointer(export, temp2d, "DXC", ALLOC=.true., _RC)
      temp2d = DycoreGrid%dxc

      call MAPL_GetPointer(export, temp2d, "DYC", ALLOC=.true., _RC)
      temp2d = DycoreGrid%dyc

      call MAPL_GetPointer(export, temp2d, "AREA", ALLOC=.true., _RC)
      temp2d = DycoreGrid%area

      ! ======================================================================
      !ALT: the next section addresses the problem when export variables have been
      !     assigned values during Initialize. To prevent "connected" exports
      !     being overwritten by DEFAULT in the Import spec in the other component
      !     we label them as being "initailized by restart". A better solution
      !     would be to move the computation to phase 2 of Initialize and
      !     eliminate this section alltogether
      ! ======================================================================
      call ESMF_StateGet(export, "PREF", field, _RC)
      call MAPL_AttributeSet(field, NAME="MAPL_InitStatus", VALUE=MAPL_InitialRestart, _RC)

      call ESMF_StateGet(export, "PLE", field, _RC)
      call MAPL_AttributeSet(field, NAME="MAPL_InitStatus", VALUE=MAPL_InitialRestart, _RC)

      call ESMF_StateGet(export, "U", field, _RC)
      call MAPL_AttributeSet(field, NAME="MAPL_InitStatus", VALUE=MAPL_InitialRestart, _RC)

      call ESMF_StateGet(export, "V", field, _RC)
      call MAPL_AttributeSet(field, NAME="MAPL_InitStatus", VALUE=MAPL_InitialRestart, _RC)

      call ESMF_StateGet(export, "T", field, _RC)
      call MAPL_AttributeSet(field, NAME="MAPL_InitStatus", VALUE=MAPL_InitialRestart, _RC)

      call MAPL_GridCompGetResource(gc, "FV3_STANDALONE", fv3_standalone, default=0, _RC)
      if (fv3_standalone /= 0) then
         call ESMF_StateGet(import, "TRADV", tradv, _RC)
         call ESMF_StateGet(export, "TRADVEX", tradvex, _RC)
         call ESMF_FieldBundleGet(tradv, fieldCount=numTracers, _RC)
         do i=1,numTracers
            call ESMF_FieldBundleGet(tradv, fieldIndex=i, field=field, _RC)
            call MAPL_FieldBundleAdd(tradvex, field, _RC)
         enddo
      end if

      !=====Begin intemittent replay=======================

      ! Set the intermittent replay alarm, if needed.
      ! Note that it is a non-sticky alarm
      ! and is set to ringing on first step. So it will
      ! work whether the clock is backed-up and ticked
      ! or not.

      call MAPL_GridCompGetResource(gc, "REPLAY_MODE", ReplayMode, default="NoReplay", _RC)
      if (adjustl(ReplayMode) == "Intermittent") then
         call MAPL_GridCompGetResource(gc, "REPLAY_INTERVAL", DNS_INTERVAL, default=21600., _RC)
         call ESMF_TimeIntervalSet(intv, s=nint(DNS_INTERVAL), _RC)
         alarm = ESMF_AlarmCreate(name="INTERMITTENT", clock=clock, ringInterval=intv, sticky=.false., _RC)
         call ESMF_AlarmRingerOn(alarm, _RC)
      end if

      !========End intermittent replay========================

      ! call MAPL_TimerOff(MAPL,"INITIALIZE")
      ! call MAPL_TimerOff(MAPL,"TOTAL")

      _RETURN(_SUCCESS)
   end subroutine Initialize

   !BOP
   !IROUTINE: Run
   !DESCRIPTION: This is the first Run stage of FV. It is the container
   !    for the dycore calculations. Subroutines from the core are
   !    invoked to do most of the work. A second run method, descibed below,
   !    adds the import tendencies from external sources to the FV
   !    variables.
   !
   !    In addition to computing and adding all dynamical contributions
   !    to the FV variables (i.e., winds, pressures, and temperatures),
   !    this method advects an arbitrary number of  tracers. These appear
   !    in a ``Friendly'' bundle in the IMPORT state and are updated with
   !    the advective tendency.

   !INTERFACE:

   subroutine Run(gc, import, export, clock, rc)

      !ARGUMENTS:
      type(ESMF_GridComp) :: gc
      type(ESMF_State) :: import
      type(ESMF_State) :: export
      type(ESMF_Clock) :: clock
      integer, intent(out) :: rc
      !EOP

      integer :: status
      type(ESMF_FieldBundle) :: bundle, ana_bundle
      type(ESMF_Field) :: field, ana_field
      type(ESMF_Config) :: cf
      type(ESMF_Alarm) :: alarm
      type(ESMF_Grid) :: esmfgrid, ana_grid
      type(ESMF_Time) :: current_time, RefTime
      class(AbstractRegridder), pointer :: L2C, C2L

      type(MAPL_MetaComp), pointer :: MAPL

      type(DYN_wrap) :: wrap
      type(DynState), pointer :: STATE
      type(DynGrid), pointer :: grid
      type(DynVars), pointer :: vars

      integer  :: NQ
      integer  :: IM, JM, KM
      integer  :: NKE, NPHI
      integer  :: NUMVARS
      integer  :: ifirstxy, ilastxy, jfirstxy, jlastxy
      integer  :: kend, i, j, K, L, n
      integer  :: im_replay,jm_replay
      logical, parameter :: convt = .false. ! Until this is run with full physics
      logical  :: is_shutoff, is_ringing

      real(r8),     pointer :: phisxy(:,:)
      real(kind=4), pointer ::   phis(:,:)

      real(r8), allocatable ::    plk(:,:,:) ! pl**kappa
      real(r8), allocatable ::   pkxy(:,:,:) ! pe**kappa
      real(r8), allocatable ::    pe0(:,:,:) ! edge-level pressure before dynamics
      real(r8), allocatable ::    pe1(:,:,:) ! edge-level pressure after dynamics
      real(r8), allocatable ::     pl(:,:,:) ! mid-level pressure
      real(r8), allocatable :: tempxy(:,:,:) ! mid-level temperature

      real(r8), allocatable ::     ua(:,:,:) ! temporary array
      real(r8), allocatable ::     va(:,:,:) ! temporary array
      real(r8), allocatable ::     uc(:,:,:) ! temporary array
      real(r8), allocatable ::     vc(:,:,:) ! temporary array
      real(r8), allocatable ::    uc0(:,:,:) ! temporary array
      real(r8), allocatable ::    vc0(:,:,:) ! temporary array
      real(r8), allocatable ::     ur(:,:,:) ! temporary array
      real(r8), allocatable ::     vr(:,:,:) ! temporary array
      real(r8), allocatable ::     qv(:,:,:) ! temporary array
      real(r8), allocatable ::     ql(:,:,:) ! temporary array
      real(r8), allocatable ::     qi(:,:,:) ! temporary array
      real(r8), allocatable ::     qr(:,:,:) ! temporary array
      real(r8), allocatable ::     qs(:,:,:) ! temporary array
      real(r8), allocatable ::     qg(:,:,:) ! temporary array
      real(r8), allocatable ::  qdnew(:,:,:) ! temporary array
      real(r8), allocatable ::  qdold(:,:,:) ! temporary array
      real(r8), allocatable ::  qvold(:,:,:) ! temporary array
      real(r8), allocatable ::  qlold(:,:,:) ! temporary array
      real(r8), allocatable ::  qiold(:,:,:) ! temporary array
      real(r8), allocatable ::  qrold(:,:,:) ! temporary array
      real(r8), allocatable ::  qsold(:,:,:) ! temporary array
      real(r8), allocatable ::  qgold(:,:,:) ! temporary array
      real(r8), allocatable ::delpold(:,:,:) ! temporary array
      real(r8), allocatable ::     ox(:,:,:) ! temporary array
      real(r8), allocatable ::     zl(:,:,:) ! temporary array
      real(r8), allocatable ::    zle(:,:,:) ! temporary array
      real(r8), allocatable ::  logpe(:,:,:) ! temporary array
      real(r8), allocatable ::   delp(:,:,:) ! temporary array
      real(r8), allocatable ::   dudt(:,:,:) ! temporary array
      real(r8), allocatable ::   dvdt(:,:,:) ! temporary array
      real(r8), allocatable ::   dtdt(:,:,:) ! temporary array
      real(r8), allocatable ::   dqdt(:,:,:) ! temporary array
      real(r8), allocatable ::  dthdt(:,:,:) ! temporary array
      real(r8), allocatable ::  ddpdt(:,:,:) ! temporary array
      real(r8), allocatable ::  dpedt(:,:,:) ! temporary array
      real(FVPRC), allocatable :: tmp3d (:,:,:) ! temporary array
      real(FVPRC), allocatable ::  vort (:,:,:) ! temporary array
      real(FVPRC), allocatable ::  divg (:,:,:) ! temporary array
      real(r8), allocatable ::     dmdt(:,:) ! temporary array

      real(r8), allocatable :: qsum1 (:,:)   ! Vertically Integrated Variable
      real(r4), allocatable :: qsum2 (:,:)   ! Vertically Integrated Variable

      real(r8), allocatable :: penrg (:,:)   ! Vertically Integrated Cp*T
      real(r8), allocatable :: kenrg (:,:)   ! Vertically Integrated K
      real(r8), allocatable :: tenrg (:,:)   ! PHIS*(Psurf-Ptop)
      real(r8), allocatable :: penrg0(:,:)   ! Vertically Integrated Cp*T
      real(r8), allocatable :: kenrg0(:,:)   ! Vertically Integrated K
      real(r8), allocatable :: tenrg0(:,:)   ! PHIS*(Psurf-Ptop)
      real(r8), allocatable :: kedyn (:,:)
      real(r8), allocatable :: pedyn (:,:)
      real(r8), allocatable :: tedyn (:,:)

      real(kind=4), allocatable :: dqvdtanaint1(:,:)
      real(kind=4), allocatable :: dqvdtanaint2(:,:)
      real(kind=4), allocatable :: dqldtanaint1(:,:)
      real(kind=4), allocatable :: dqldtanaint2(:,:)
      real(kind=4), allocatable :: dqidtanaint1(:,:)
      real(kind=4), allocatable :: dqidtanaint2(:,:)
      real(kind=4), allocatable :: doxdtanaint1(:,:)
      real(kind=4), allocatable :: doxdtanaint2(:,:)
      real(kind=4), allocatable :: dthdtanaint1(:,:)
      real(kind=4), allocatable :: dthdtanaint2(:,:)

      real(kind=4), allocatable :: tropp1(:,:)   ! Tropopause Pressure
      real(kind=4), allocatable :: tropp2(:,:)   ! Tropopause Pressure
      real(kind=4), allocatable :: tropp3(:,:)   ! Tropopause Pressure
      real(kind=4), allocatable :: tropt (:,:)   ! Tropopause Temperature
      real(kind=4), allocatable :: tropq (:,:)   ! Tropopause Specific Humidity

      real(r8), allocatable :: omaxyz(:,:,:) ! vertical pressure velocity (pa/sec)
      real(r8), allocatable :: epvxyz(:,:,:) ! ertel's potential vorticity

      real(r8), allocatable :: cxxyz(:,:,:)  ! Accumulated eastward courant numbers
      real(r8), allocatable :: cyxyz(:,:,:)  ! Accumulated northward courant numbers
      real(r8), allocatable :: mfxxyz(:,:,:) ! Accumulated eastward mass flux
      real(r8), allocatable :: mfyxyz(:,:,:) ! Accumulated northward mass flux
      real(r8), allocatable :: mfzxyz(:,:,:) ! Accumulated vertical mass flux

      real(FVPRC)              :: dt            ! Dynamics time step
      real(r8), allocatable :: trsum1(:)     ! Global Sum of Tracers before Add_Incs
      real(r8), allocatable :: trsum2(:)     ! Global Sum of Tracers after  Add_Incs

      real(kind=4), pointer ::      dudtana(:,:,:)
      real(kind=4), pointer ::      dvdtana(:,:,:)
      real(kind=4), pointer ::      dtdtana(:,:,:)
      real(kind=4), pointer ::     ddpdtana(:,:,:)
      real(kind=4), pointer ::       qctmp (:,:,:)
      real(kind=4), pointer ::       dqldt (:,:,:)
      real(kind=4), pointer ::       dqidt (:,:,:)
      real(kind=4), pointer ::       doxdt (:,:,:)
      real(kind=4), pointer ::      dqvana (:,:,:)
      real(kind=4), pointer ::      dqlana (:,:,:)
      real(kind=4), pointer ::      dqiana (:,:,:)
      real(kind=4), pointer ::      dqrana (:,:,:)
      real(kind=4), pointer ::      dqsana (:,:,:)
      real(kind=4), pointer ::      dqgana (:,:,:)
      real(kind=4), pointer ::      doxana (:,:,:)
      real(kind=4), pointer ::       temp3d(:,:,:)
      real(kind=4), pointer ::       vtmp3d(:,:,:)
      real(kind=4), pointer ::         area(:,:)
      real(kind=4), pointer ::       temp2d(:,:)
      real(kind=4), pointer ::       tempu (:,:)
      real(kind=4), pointer ::       tempv (:,:)
      real(kind=4), allocatable ::   cubetemp3d(:,:,:)
      real(kind=4), allocatable ::   cubevtmp3d(:,:,:)

      real(kind=4), pointer :: uh25(:,:)
      real(kind=4), pointer :: uh03(:,:)
      real(kind=4), pointer :: srh01(:,:)
      real(kind=4), pointer :: srh03(:,:)
      real(kind=4), pointer :: srh25(:,:)

      real(r8),     allocatable ::   uatmp(:,:,:)
      real(r8),     allocatable ::   vatmp(:,:,:)
      real(r8),     allocatable ::   udtmp(:,:,:)
      real(r8),     allocatable ::   vdtmp(:,:,:)

      character(len=ESMF_MAXSTR), ALLOCATABLE       :: NAMES (:)
      character(len=ESMF_MAXSTR), ALLOCATABLE       :: NAMES0(:)
      character(len=ESMF_MAXSTR) :: IAm
      character(len=ESMF_MAXSTR) :: comp_name
      character(len=ESMF_MAXSTR) :: STRING
      character(len=ESMF_MAXSTR) :: ReplayFile
      character(len=ESMF_MAXSTR) :: ReplayType
      character(len=ESMF_MAXSTR) :: ReplayMode
      character(len=ESMF_MAXSTR) :: cremap,tremap
      character(len=ESMF_MAXSTR) :: uname,vname,tname,qname,psname,dpname,o3name,rgrid,tvar

      type(MAPL_SunOrbit)        :: ORBIT
      real(kind=4), pointer      :: LATS(:,:)
      real(kind=4), pointer      :: LONS(:,:)
      real(kind=4), allocatable  ::  ZTH(:,:)
      real(kind=4), allocatable  ::  SLR(:,:)

      real                  :: rc_blend_p_above
      real                  :: rc_blend_p_below
      real                  :: sclinc
      integer               :: rc_blend

      real                  :: HGT_SURFACE

      character(len=ESMF_MAXSTR) :: ANA_IS_WEIGHTED
      logical                    ::     IS_WEIGHTED

      type(DynTracers)            :: qqq       ! Specific Humidity
      type(DynTracers)            :: ooo       ! ox
      logical LCONSV, LFILL
      integer  CONSV,  FILL
      integer nx_ana, ny_ana

      logical, save                       :: firstime=.true.
      logical                             :: adjustTracers
      type(ESMF_Alarm)                    :: predictorAlarm
      type(ESMF_Grid)                     :: bgrid
      integer                             :: pos
      integer                             :: nqt
      logical                             :: tend
      logical                             :: exclude
      character(len=ESMF_MAXSTR)          :: tmpstring
      character(len=ESMF_MAXSTR)          :: fieldname
      character(len=ESMF_MAXSTR)          :: adjustTracerMode
      character(len=ESMF_MAXSTR), allocatable :: xlist(:)
      character(len=ESMF_MAXSTR), allocatable :: biggerlist(:)
      integer, parameter                  :: XLIST_MAX = 60
      logical                             :: isPresent

      logical                             :: doEnergetics
      logical                             :: doTropvars

      integer :: FV3_STANDALONE

      character(len=ESMF_MAXSTR) :: myTracer
      integer :: itracer
      real(kind=r8) :: t1, t2, dyn_run_timer

      Iam = "Run"
      call ESMF_GridCompGet(gc, name=comp_name, CONFIG=cf, grid=esmfgrid, _RC)
      Iam = trim(comp_name) // trim(Iam)

      call ESMF_GridValidate(esmfgrid, _RC)

      ! Retrieve the pointer to the generic state
      call MAPL_GetObjectFromGC(gc, MAPL, _RC)

      call MAPL_TimerOn(MAPL, "TOTAL")
      call MAPL_TimerOn(MAPL, "RUN")

      call MAPL_Get(MAPL, LONS=LONS, LATS=LATS, _RC)

      call MAPL_GetPointer(export, temp2d, 'LONS', _RC)
      if( associated(temp2D) ) temp2d = LONS
      call MAPL_GetPointer(export, temp2d, 'LATS', _RC)
      if( associated(temp2D) ) temp2d = LATS

      ! Retrieve the pointer to the internal state
      call ESMF_UserCompGetInternalState(gc, 'DYNstate', wrap, status)
      VERIFY_(STATUS)
      state => wrap%dyn_state

      vars  => state%vars   ! direct handle to control variables
      grid  => state%grid   ! direct handle to grid
      dt    =  state%dt     ! dynamics time step (large)

      ifirstxy = grid%is
      ilastxy  = grid%ie
      jfirstxy = grid%js
      jlastxy  = grid%je

      im       = grid%npx
      jm       = grid%npy
      km       = grid%npz

      is_ringing = ESMF_AlarmIsRinging(STATE%ALARMS(TIME_TO_RUN), _RC)
      if (.not. is_ringing) return

      ! Allocate Arrays
      ! ---------------
      ALLOCATE(   delp(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(   dudt(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(   dvdt(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(   dtdt(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(   dqdt(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(  dthdt(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(  ddpdt(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(  dpedt(ifirstxy:ilastxy,jfirstxy:jlastxy,km+1) )
      ALLOCATE( tempxy(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(    pe0(ifirstxy:ilastxy,jfirstxy:jlastxy,km+1) )
      ALLOCATE(    pe1(ifirstxy:ilastxy,jfirstxy:jlastxy,km+1) )
      ALLOCATE(     pl(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(     ua(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(     va(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(     uc(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(     vc(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(    uc0(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(    vc0(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(     ur(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(     vr(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(     qv(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(     ql(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(     qi(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(     qr(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(     qs(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(     qg(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(     ox(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )

      ALLOCATE(  qsum1(ifirstxy:ilastxy,jfirstxy:jlastxy)    )
      ALLOCATE(  qsum2(ifirstxy:ilastxy,jfirstxy:jlastxy)    )

      ALLOCATE(   dmdt(ifirstxy:ilastxy,jfirstxy:jlastxy)    )

      doEnergetics=.false.
      call MAPL_GetPointer(export, temp2D, 'KEANA', _RC)
      if(associated(temp2D)) doEnergetics=.true.
      call MAPL_GetPointer(export, temp2D, 'PEANA', _RC)
      if(associated(temp2D)) doEnergetics=.true.
      call MAPL_GetPointer(export, temp2D, 'TEANA', _RC)
      if(associated(temp2D)) doEnergetics=.true.
      call MAPL_GetPointer(export, temp2D, 'KEDYN', _RC)
      if(associated(temp2D)) doEnergetics=.true.
      call MAPL_GetPointer(export, temp2D, 'PEDYN', _RC)
      if(associated(temp2D)) doEnergetics=.true.
      call MAPL_GetPointer(export, temp2D, 'TEDYN', _RC)
      if(associated(temp2D)) doEnergetics=.true.
      if (doEnergetics) then
         ALLOCATE(  kedyn(ifirstxy:ilastxy,jfirstxy:jlastxy)    )
         ALLOCATE(  pedyn(ifirstxy:ilastxy,jfirstxy:jlastxy)    )
         ALLOCATE(  tedyn(ifirstxy:ilastxy,jfirstxy:jlastxy)    )
         ALLOCATE(  kenrg(ifirstxy:ilastxy,jfirstxy:jlastxy)    )
         ALLOCATE(  penrg(ifirstxy:ilastxy,jfirstxy:jlastxy)    )
         ALLOCATE(  tenrg(ifirstxy:ilastxy,jfirstxy:jlastxy)    )
         ALLOCATE( kenrg0(ifirstxy:ilastxy,jfirstxy:jlastxy)    )
         ALLOCATE( penrg0(ifirstxy:ilastxy,jfirstxy:jlastxy)    )
         ALLOCATE( tenrg0(ifirstxy:ilastxy,jfirstxy:jlastxy)    )
      endif

      ALLOCATE(   vort(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(   divg(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )

      ALLOCATE(  tmp3d(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )

      ALLOCATE( phisxy   (ifirstxy:ilastxy,jfirstxy:jlastxy     ) )
      ALLOCATE(    plk   (ifirstxy:ilastxy,jfirstxy:jlastxy,km  ) )
      ALLOCATE(   pkxy   (ifirstxy:ilastxy,jfirstxy:jlastxy,km+1) )
      ALLOCATE(     zl   (ifirstxy:ilastxy,jfirstxy:jlastxy,km  ) )
      ALLOCATE(    zle   (ifirstxy:ilastxy,jfirstxy:jlastxy,km+1) )
      ALLOCATE(  logpe   (ifirstxy:ilastxy,jfirstxy:jlastxy,km+1) )
      ALLOCATE( omaxyz   (ifirstxy:ilastxy,jfirstxy:jlastxy,km  ) )
      ALLOCATE( epvxyz   (ifirstxy:ilastxy,jfirstxy:jlastxy,km  ) )
      ALLOCATE(  cxxyz   (ifirstxy:ilastxy,jfirstxy:jlastxy,km  ) )
      ALLOCATE(  cyxyz   (ifirstxy:ilastxy,jfirstxy:jlastxy,km  ) )
      ALLOCATE( mfxxyz   (ifirstxy:ilastxy,jfirstxy:jlastxy,km  ) )
      ALLOCATE( mfyxyz   (ifirstxy:ilastxy,jfirstxy:jlastxy,km  ) )
      ALLOCATE( mfzxyz   (ifirstxy:ilastxy,jfirstxy:jlastxy,km+1) )

      ! Report advected friendlies
      call ESMF_StateGet(import, 'TRADV', bundle, _RC)

      !-------------------------------------------------------------------
      ! ALT: this section attempts to limit the amount of advected tracers
      !-------------------------------------------------------------------
      adjustTracers = .false.
      call MAPL_GetResource(MAPL, adjustTracerMode, 'EXCLUDE_ADVECTION_TRACERS:', default='ALWAYS', _RC)
      if (adjustTracerMode == 'ALWAYS') then
         adjustTracers = .true.
      else if (adjustTracerMode == 'PREDICTOR') then
         !get PredictorAlarm from clock
         call ESMF_ClockGetAlarm(clock, alarmName='PredictorAlarm', alarm=PredictorAlarm, rc=status)
         if (status == ESMF_SUCCESS) then
            !check if ringing
            if (ESMF_AlarmIsRinging(predictorAlarm)) then
               adjustTracers = .true.
            end if
         end if
      else
         call WRITE_PARALLEL('Invalid option, ignored')
         adjustTracers = .false.
      end if
      if (adjustTracers) then
         if (firstime) then
            firstime=.false.
            ! get the list of excluded tracers from resource
            n = 0
            call ESMF_ConfigFindLabel(cf, 'EXCLUDE_ADVECTION_TRACERS_LIST:', isPresent=isPresent, _RC)
            if(isPresent .or. (AdvCore_Advection >= 1)) then
               tend  = .false.
               allocate(xlist(XLIST_MAX), stat=status)
               VERIFY_(STATUS)
               if (isPresent) then
                  do while (.not.tend)
                     call ESMF_ConfigGetAttribute(cf, value=tmpstring, default='', rc=status) !ALT: we don't check return status
                     if (tmpstring /= '')  then
                        n = n + 1
                        if (n > size(xlist)) then
                           allocate(biggerlist(2*n), _STAT)
                           biggerlist(1:n-1)=xlist
                           call move_alloc(from=biggerlist, to=xlist)
                        end if
                        xlist(n) = tmpstring
                     end if
                     call ESMF_ConfigNextLine(cf, tableEnd=tend, _RC)
                  enddo
               endif
            end if

            ! Count the number of tracers
            call ESMF_FieldBundleGet(bundle, grid=bgrid, fieldCount=nqt, _RC)
            BundleAdv = ESMF_FieldBundleCreate(name='xTRADV', _RC)
            call ESMF_FieldBundleSet(BundleAdv, grid=bgrid, _RC)
            !loop over NQ in TRADV
            do i = 1, nqt
               !get field from TRADV and its name
               call ESMF_FieldBundleGet(bundle, fieldIndex=i, field=field, _RC)
               call ESMF_FieldGet(FIELD, name=fieldname, _RC)
               !exclude everything that is not cloud/water species
               if ( (AdvCore_Advection >= 1    ) .and. &
                    (trim(fieldname) /= 'Q'       ) .and. &
                    (trim(fieldname) /= 'QLCN'    ) .and. &
                    (trim(fieldname) /= 'QLLS'    ) .and. &
                    (trim(fieldname) /= 'QICN'    ) .and. &
                    (trim(fieldname) /= 'QILS'    ) .and. &
                    (trim(fieldname) /= 'CLCN'    ) .and. &
                    (trim(fieldname) /= 'CLLS'    ) .and. &
                    (trim(fieldname) /= 'NCPL'    ) .and. &
                    (trim(fieldname) /= 'NCPI'    ) .and. &
                    (trim(fieldname) /= 'QRAIN'   ) .and. &
                    (trim(fieldname) /= 'QSNOW'   ) .and. &
                    (trim(fieldname) /= 'QGRAUPEL') ) then
                  write(STRING,'(A,A)') "FV3+ADV is excluding ", trim(fieldname)
                  call WRITE_PARALLEL(trim(STRING))
                  n = n + 1
                  if (n > size(xlist)) then
                     allocate(biggerlist(2*n), _STAT)
                     biggerlist(1:n-1)=xlist
                     call move_alloc(from=biggerlist, to=xlist)
                  end if
                  xlist(n) = trim(fieldname)
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
                  call MAPL_FieldBundleAdd(BundleAdv, field, _RC)
               end if
            end do

            if (allocated(xlist)) then
               !   ! Just in case xlist was allocated, but nothing was in it, could have garbage
               !   if (n > 0) then
               !      call ESMF_FieldBundleRemove(BUNDLE, fieldNameList=xlist, &
               !         relaxedFlag=.true., rc=status)
               !      VERIFY_(STATUS)
               !   end if
               deallocate(xlist)
            end if

         end if ! firstime
         bundle = BundleAdv ! replace TRADV
      else
         BundleAdv = bundle ! replace with TRADV
      end if ! adjustTracers

      call ESMF_FieldBundleGet(bundle, fieldCount=NQ, _RC)

      if (NQ > 0) then
         allocate(NAMES(NQ), _STAT)
         call ESMF_FieldBundleGet(bundle, itemorderflag=ESMF_ITEMORDER_ADDORDER, fieldNameList=NAMES, _RC)
         if( .not.allocated(names0) ) then
            allocate(NAMES0(NQ), _STAT)
            NAMES0 = NAMES
         endif
      endif

      ! Surface Geopotential from import state
      call MAPL_GetPointer(import, PHIS, 'PHIS', _RC)

      phisxy = real(phis,kind=r8)

      ! Get tracers from import State (Note: Contains Updates from Analysis)
      call PULL_Q (STATE, import, qqq, NXQ, RC=rc)

      !-----------------------------
      ! end of fewer_tracers-section
      !-----------------------------

      do k=1,size(names)
         pos = index(names(k), '::')
         if(pos > 0) then
            if( (names(k)(pos+2:))=='OX' ) ooo = vars%tracer(k)
         elseif(names(k)=='Q') then
            qqq = vars%tracer(k)
         end if
      end do

      ! WMP Begin REPLAY/ANA section
      call ESMF_ConfigGetAttribute(cf, FV3_STANDALONE, label="FV3_STANDALONE:", default=0, _RC)
      if (FV3_STANDALONE == 0) then
         call MAPL_TimerOn(MAPL, "-DYN_ANA")
         call ESMF_ClockGetAlarm(clock, 'ReplayShutOff', alarm, _RC)
         is_shutoff = ESMF_AlarmIsRinging(alarm, _RC)
      else
         is_shutoff = .true.
      end if

      if (.not. is_shutoff) then
         ! If requested, do Intermittent Replay

         call MAPL_GetResource(MAPL, ReplayMode, 'REPLAY_MODE:', default="NoReplay", _RC)

         REPLAYING: if(adjustl(ReplayMode)=="Intermittent") then

            ! If replay alarm is ringing, we need to reset state
            call ESMF_ClockGetAlarm(clock, 'INTERMITTENT', alarm, _RC)
            call ESMF_ClockGet(clock, CurrTime=current_time, _RC)

            is_ringing = ESMF_AlarmIsRinging(alarm, _RC)

            RefTime = current_time

            call check_replay_time_(is_ringing)
            TIME_TO_REPLAY: if(is_ringing) then

               call ESMF_AlarmRingerOff(alarm, _RC)

               ! Read in file name of field to replay to and all other relavant resources
               call MAPL_GetResource(MAPL, ReplayFile, 'REPLAY_FILE:', _RC)
               call MAPL_GetResource(MAPL, ReplayType, 'REPLAY_TYPE:', default="FULL", _RC)
               call MAPL_GetResource(MAPL, im_replay, label="REPLAY_IM:", _RC)
               call MAPL_GetResource(MAPL, jm_replay, label="REPLAY_JM:", _RC)

               call MAPL_GetResource(MAPL, psname, label="REPLAY_PSNAME:", default="NULL",  _RC)
               call MAPL_GetResource(MAPL, dpname, label="REPLAY_DPNAME:", default="delp",  _RC)
               call MAPL_GetResource(MAPL, uname, label="REPLAY_UNAME:", default="uwnd",   _RC)
               call MAPL_GetResource(MAPL, vname, label="REPLAY_VNAME:", default="vwnd",   _RC)
               call MAPL_GetResource(MAPL, tname, label="REPLAY_TNAME:", default="theta",  _RC)
               call MAPL_GetResource(MAPL, qname, label="REPLAY_QNAME:", default="sphu",   _RC)
               call MAPL_GetResource(MAPL, o3name, label="REPLAY_O3NAME:", default="ozone", _RC)

               call MAPL_GetResource(MAPL, rgrid, label="REPLAY_GRID:", default="D-GRID", _RC)
               call MAPL_GetResource(MAPL, tvar, label="REPLAY_TVAR:", default="THETAV", _RC)

               call MAPL_GetResource(MAPL, CREMAP, label="REPLAY_REMAP:", default="no", _RC)
               call MAPL_GetResource(MAPL, TREMAP, label="REPLAY_REMAP_ALL_TRACERS:", default="yes", _RC)

               call MAPL_GetResource(MAPL, rc_blend, 'REPLAY_BLEND:', default=0, _RC)
               call MAPL_GetResource(MAPL, rc_blend_p_above, 'REPLAY_BLEND_P_ABOVE:', default=10.0, _RC)
               call MAPL_GetResource(MAPL, rc_blend_p_below, 'REPLAY_BLEND_P_BELOW:', default=100.0, _RC)

               call MAPL_GetResource(MAPL, sclinc, label='SCLINC:', default=1.0, _RC)

               ! Read the fields to be reset into a bundle
               call ESMF_ConfigGetAttribute(cf, nx_ana, label ='NX:', _RC)
               call ESMF_ConfigGetAttribute(cf, ny_ana, label ='NY:', _RC)

               block
                  use MAPL_LatLonGridFactoryMod
                  ana_grid = grid_manager%make_grid( &
                       & LatLonGridFactory(im_world=IM_REPLAY, jm_world=JM_REPLAY, lm=km, &
                       & nx=nx_ana, ny=ny_ana, rc=status))
                  VERIFY_(status)
               end block

               ana_bundle = ESMF_FieldBundleCreate(_RC)
               call ESMF_FieldBundleSet(ana_bundle, grid=ana_grid, _RC)

               call MAPL_CFIORead(ReplayFile, RefTime, ana_bundle, _RC)

               ! Create transform from lat-lon to cubed
               l2c => regridder_manager%make_regridder(ana_grid, esmfgrid, REGRID_METHOD_BILINEAR, _RC)

               ! Fill the state variables from the bundle only if
               ! the corresponding fields are there

               ! soon dump_n_splash will go; we'll have instead:
               !    call get_inc_on_anagrid_ - this will convert the internal state to
               !      ana-grid, diff with what's in file and produce what incremental_
               !      normally works from - a knob will tell incremental_ where fields
               !      are in memory or need reading from file.
               !    call incremental_
               !    call state_remap_
               if (trim(ReplayType)=='FULL') then
                  call dump_n_splash_
               else
                  call incremental_
               endif
               call state_remap_

               ! Done with replay; clean-up
               call ESMF_FieldBundleGet(ana_bundle, FieldCount=NUMVARS, _RC)

               do k=1,NUMVARS
                  call ESMF_FieldBundleGet(ana_bundle, k, ana_field, _RC)
                  call MAPL_FieldDestroy(ana_field, _RC)
               end do

               call ESMF_FieldBundleDestroy(ana_bundle, _RC)
            end if TIME_TO_REPLAY
         end if REPLAYING

         ! Create Local Copy of QV and OX (Contains Updates from Analysis)
         ox = 0.0
         qv = 0.0

         if (.not. ADIABATIC) then
            do k=1,size(names)

               pos = index(names(k),'::')
               if(pos > 0) then
                  if( (names(k)(pos+2:))=='OX' ) then
                     if ( (ooo%is_r4) .and. associated(ooo%content_r4) ) then
                        if (size(ox)==size(ooo%content_r4)) then
                           ox = ooo%content_r4
                        endif
                     elseif (associated(ooo%content)) then
                        if (size(ox)==size(ooo%content)) then
                           ox = ooo%content
                        endif
                     endif
                  endif
               endif

               if( trim(names(k))=='Q'  ) then
                  if ( (qqq%is_r4) .and. associated(qqq%content_r4) ) then
                     if (size(qv)==size(qqq%content_r4)) then
                        qv = qqq%content_r4
                     endif
                  elseif (associated(qqq%content)) then
                     if (size(qv)==size(qqq%content)) then
                        qv = qqq%content
                     endif
                  endif
                  _ASSERT(all(qv >= 0.0),'Before AnaAddIncs: negative or nan water vapor detected')
               endif

            enddo
         endif

         ! Diagnostics Before Analysis Increments are Added
         ALLOCATE(delpold(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
         ALLOCATE(  qdnew(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
         ALLOCATE(  qdold(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
         ALLOCATE(  qvold(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
         ALLOCATE(  qlold(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
         ALLOCATE(  qiold(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
         ALLOCATE(  qrold(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
         ALLOCATE(  qsold(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
         ALLOCATE(  qgold(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )

         call MAPL_GetPointer(import, dqvana, 'DQVANA', _RC)   ! Get QV Increment from Analysis
         call MAPL_GetPointer(import, dqlana, 'DQLANA', _RC)   ! Get QL Increment from Analysis
         call MAPL_GetPointer(import, dqiana, 'DQIANA', _RC)   ! Get QI Increment from Analysis
         call MAPL_GetPointer(import, dqrana, 'DQRANA', _RC)   ! Get QR Increment from Analysis
         call MAPL_GetPointer(import, dqsana, 'DQSANA', _RC)   ! Get QS Increment from Analysis
         call MAPL_GetPointer(import, dqgana, 'DQGANA', _RC)   ! Get QG Increment from Analysis
         call MAPL_GetPointer(import, doxana, 'DOXANA', _RC)   ! Get OX Increment from Analysis

         QL = 0.0
         QI = 0.0
         QR = 0.0
         QS = 0.0
         QG = 0.0
         do N = 1,size(names)
            if( trim(names(N)).eq.'QLCN' .or. &
                 trim(names(N)).eq.'QLLS' ) then
               if( state%vars%tracer(N)%is_r4 ) then
                  QL = QL + state%vars%tracer(N)%content_r4
               else
                  QL = QL + state%vars%tracer(N)%content
               endif
            endif
            if( trim(names(N)).eq.'QICN' .or. &
                 trim(names(N)).eq.'QILS' ) then
               if( state%vars%tracer(N)%is_r4 ) then
                  QI = QI + state%vars%tracer(N)%content_r4
               else
                  QI = QI + state%vars%tracer(N)%content
               endif
            endif
            if( trim(names(N)).eq.'QRAIN' ) then
               if( state%vars%tracer(N)%is_r4 ) then
                  QR = state%vars%tracer(N)%content_r4
               else
                  QR = state%vars%tracer(N)%content
               endif
            endif
            if( trim(names(N)).eq.'QSNOW' ) then
               if( state%vars%tracer(N)%is_r4 ) then
                  QS = state%vars%tracer(N)%content_r4
               else
                  QS = state%vars%tracer(N)%content
               endif
            endif
            if( trim(names(N)).eq.'QGRAUPEL' ) then
               if( state%vars%tracer(N)%is_r4 ) then
                  QG = state%vars%tracer(N)%content_r4
               else
                  QG = state%vars%tracer(N)%content
               endif
            endif
         enddo
         QVOLD = QV-DQVANA
         QLOLD = QL-DQLANA
         QIOLD = QI-DQIANA
         QROLD = QR-DQRANA
         QSOLD = QS-DQSANA
         QGOLD = QG-DQGANA

         ! Get A-grid winds
         call getAllWinds(vars%u, vars%v, UR=ur, VR=vr)

         delp   = vars%pe(:,:,2:)  -vars%pe(:,:,:km)   ! Pressure Thickness
         dmdt   = vars%pe(:,:,km+1)-vars%pe(:,:,1)     ! Psurf-Ptop
         tempxy = vars%pt * (1.0+eps*(qv-dqvana))       ! Compute THV Before Analysis Update

         if (doEnergetics) then
            call Energetics(ur, vr, tempxy, vars%pe, delp, vars%pkz, phisxy, kenrg, penrg, tenrg)
         end if

         ! DUDTANA
         call MAPL_GetPointer(export, dudtana, 'DUDTANA', _RC)
         if( associated(dudtana) ) dudtana = ur

         ! DVDTANA
         call MAPL_GetPointer(export, dvdtana, 'DVDTANA', _RC)
         if( associated(dvdtana) ) dvdtana = vr

         ! DTDTANA
         call MAPL_GetPointer(export, dtdtana, 'DTDTANA', _RC)
         if( associated(dtdtana) ) dtdtana = vars%pt * vars%pkz

         ! DDELPDTANA
         call MAPL_GetPointer(export, ddpdtana, 'DDELPDTANA', _RC)
         if( associated(ddpdtana) ) ddpdtana = delp

         ! DTHVDTANAINT
         call MAPL_GetPointer(export, temp2D, 'DTHVDTANAINT', _RC)
         if( associated(temp2D) ) then
            tempxy       = vars%pt*(1+eps*(qv-dqvana))   ! Set tempxy = TH*QVold (Before Analysis Update)
            ALLOCATE( dthdtanaint1(ifirstxy:ilastxy,jfirstxy:jlastxy) )
            ALLOCATE( dthdtanaint2(ifirstxy:ilastxy,jfirstxy:jlastxy) )
            dthdtanaint1 = 0.0
            do k=1,km
               dthdtanaint1 = dthdtanaint1 + tempxy(:,:,k)*delp(:,:,k)
            enddo
         endif

         ! DQVDTANAINT
         call MAPL_GetPointer(export, temp2D, 'DQVDTANAINT', _RC)
         if( associated(temp2D) ) then
            ALLOCATE( dqvdtanaint1(ifirstxy:ilastxy,jfirstxy:jlastxy) )
            ALLOCATE( dqvdtanaint2(ifirstxy:ilastxy,jfirstxy:jlastxy) )
            tempxy       = qv-dqvana   ! Set tempxy = QVold (Before Analysis Update)
            dqvdtanaint1 = 0.0
            do k=1,km
               dqvdtanaint1 = dqvdtanaint1 + tempxy(:,:,k)*delp(:,:,k)
            enddo
         endif

         ! DQLDTANAINT
         call MAPL_GetPointer(export, temp2D, 'DQLDTANAINT', _RC)
         if( associated(temp2D) ) then
            ALLOCATE( dqldtanaint1(ifirstxy:ilastxy,jfirstxy:jlastxy) )
            ALLOCATE( dqldtanaint2(ifirstxy:ilastxy,jfirstxy:jlastxy) )
            dqldtanaint1 = 0.0
            do N = 1,size(names)
               if( trim(names(N)).eq.'QLCN' .or. &
                    trim(names(N)).eq.'QLLS' ) then
                  do k=1,km
                     if( state%vars%tracer(N)%is_r4 ) then
                        dqldtanaint1 = dqldtanaint1 + state%vars%tracer(N)%content_r4(:,:,k)*delp(:,:,k)
                     else
                        dqldtanaint1 = dqldtanaint1 + state%vars%tracer(N)%content   (:,:,k)*delp(:,:,k)
                     endif
                  enddo
               endif
            enddo
            do k=1,km
               dqldtanaint1 = dqldtanaint1 - dqlana(:,:,k)*delp(:,:,k)
            enddo
         endif

         ! DQIDTANAINT
         call MAPL_GetPointer(export, temp2D, 'DQIDTANAINT', _RC)
         if( associated(temp2D) ) then
            ALLOCATE( dqidtanaint1(ifirstxy:ilastxy,jfirstxy:jlastxy) )
            ALLOCATE( dqidtanaint2(ifirstxy:ilastxy,jfirstxy:jlastxy) )
            dqidtanaint1 = 0.0
            do N = 1,size(names)
               if( trim(names(N)).eq.'QICN' .or. &
                    trim(names(N)).eq.'QILS' ) then
                  do k=1,km
                     if( state%vars%tracer(N)%is_r4 ) then
                        dqidtanaint1 = dqidtanaint1 + state%vars%tracer(N)%content_r4(:,:,k)*delp(:,:,k)
                     else
                        dqidtanaint1 = dqidtanaint1 + state%vars%tracer(N)%content   (:,:,k)*delp(:,:,k)
                     endif
                  enddo
               endif
            enddo
            do k=1,km
               dqidtanaint1 = dqidtanaint1 - dqiana(:,:,k)*delp(:,:,k)
            enddo
         endif

         ! DOXDTANAINT
         call MAPL_GetPointer(export, temp2D, 'DOXDTANAINT', _RC)
         if( associated(temp2D) ) then
            ALLOCATE( doxdtanaint1(ifirstxy:ilastxy,jfirstxy:jlastxy) )
            ALLOCATE( doxdtanaint2(ifirstxy:ilastxy,jfirstxy:jlastxy) )
            tempxy       = OX-doxana   ! Set tempxy = OXold (Before Analysis Update)
            doxdtanaint1 = 0.0
            do k=1,km
               doxdtanaint1 = doxdtanaint1 + tempxy(:,:,k)*delp(:,:,k)
            enddo
         endif

         ! Add Diabatic Forcing from Analysis to State Variables
         if (vars%nwat >= 6) then
            QDOLD = 1.0 - (QVOLD+QLOLD+QIOLD+QROLD+QSOLD+QGOLD)
            QDNEW = 1.0 - (QV   +QL   +QI   +QR   +QS   +QG   )
         else
            QDOLD = 1.0 - (QVOLD+QLOLD+QIOLD)
            QDNEW = 1.0 - (QV   +QL   +QI   )
         endif
         call MAPL_GetPointer(export, area, 'AREA', _RC)

         allocate( trsum1(nq) )
         allocate( trsum2(nq) )

         call MAPL_GetResource(MAPL, ANA_IS_WEIGHTED, label="ANA_IS_WEIGHTED:", default='NO', _RC)
         ANA_IS_WEIGHTED = ESMF_UtilStringUpperCase(ANA_IS_WEIGHTED)
         IS_WEIGHTED =   adjustl(ANA_IS_WEIGHTED)=="YES" .or. adjustl(ANA_IS_WEIGHTED)=="NO"
         _ASSERT(IS_WEIGHTED ,'needs informative message')
         IS_WEIGHTED = adjustl(ANA_IS_WEIGHTED)=="YES"

         ! Add Analysis Tendencies
         delpold = delp                            ! Old Pressure Thickness

         call ADD_INCS(MAPL, STATE, import, DT, IS_WEIGHTED=IS_WEIGHTED)

         if (DYN_DEBUG) call DEBUG_FV_STATE('ANA ADD_INCS', STATE)

         delp = vars%pe(:,:,2:)-vars%pe(:,:,:km)   ! Updated Pressure Thickness

         ! Compute Old Global Sums of Tracers over Locations where Mass has changed
         if ((.not. ADIABATIC)) then
            do n=1,NQ
               qsum1(:,:) = 0.0_r8
               if( STATE%VARS%TRACER(N)%IS_R4 ) then
                  do k=1,km
                     where( delp(:,:,k).ne.delpold(:,:,k) )
                        qsum1(:,:) = qsum1(:,:) + state%vars%tracer(n)%content_r4(:,:,k)*delpold(:,:,k)
                     end where
                  enddo
               else
                  do k=1,km
                     where( delp(:,:,k).ne.delpold(:,:,k) )
                        qsum1(:,:) = qsum1(:,:) + state%vars%tracer(n)%content   (:,:,k)*delpold(:,:,k)
                     end where
                  enddo
               endif
               where( qsum1.ne.0.0_r8 )
                  qsum2 = qsum1
               elsewhere
                  qsum2 = MAPL_UNDEF
               end where
               call MAPL_AreaMean(TRSUM1(n), qsum2, area, esmfgrid, _RC)
            enddo
         endif

         ! Update Specific Mass of Aerosol Constituents Keeping Mixing_Ratio Constant WRT_Dry_Air After ANA Updates
         if ((.not. ADIABATIC)) then
            do n=1,NQ
               if( (trim(names(n)).ne.'Q'   ) .and. &
                    (trim(names(n)).ne.'QLLS') .and. &
                    (trim(names(n)).ne.'QLCN') .and. &
                    (trim(names(n)).ne.'QILS') .and. &
                    (trim(names(n)).ne.'QICN') .and. &
                    (trim(names(n)).ne.'CLLS') .and. &
                    (trim(names(n)).ne.'CLCN') .and. &
                    (trim(names(n)).ne.'QRAIN') .and. &
                    (trim(names(n)).ne.'QSNOW') .and. &
                    (trim(names(n)).ne.'QGRAUPEL') ) then
                  if( STATE%VARS%TRACER(N)%IS_R4 ) then
                     state%vars%tracer(n)%content_r4 = state%vars%tracer(n)%content_r4 * ( QDNEW/QDOLD )
                  else
                     state%vars%tracer(n)%content    = state%vars%tracer(n)%content    * ( QDNEW/QDOLD )
                  endif
               endif
            enddo
         endif

         ! Compute New Global Sums of Tracers over Locations where Mass has changed
         if ((.not. ADIABATIC)) then
            do n=1,NQ
               qsum1(:,:) = 0.0_r8
               if( STATE%VARS%TRACER(N)%IS_R4 ) then
                  do k=1,km
                     where( delp(:,:,k).ne.delpold(:,:,k) )
                        qsum1(:,:) = qsum1(:,:) + state%vars%tracer(n)%content_r4(:,:,k)*delp(:,:,k)
                     end where
                  enddo
               else
                  do k=1,km
                     where( delp(:,:,k).ne.delpold(:,:,k) )
                        qsum1(:,:) = qsum1(:,:) + state%vars%tracer(n)%content   (:,:,k)*delp(:,:,k)
                     end where
                  enddo
               endif
               where( qsum1.ne.0.0_r8 )
                  qsum2 = qsum1
               elsewhere
                  qsum2 = MAPL_UNDEF
               end where
               call MAPL_AreaMean(TRSUM2(n), qsum2, area, esmfgrid, _RC)
            enddo
         endif

         ! Ensure Conservation of Global Mass of Aerosol Constituents After ANA Updates
         if ((.not. ADIABATIC)) then
            do n=1,NQ
               if( (trim(names(n)).ne.'Q'   ) .and. &
                    (trim(names(n)).ne.'QLLS') .and. &
                    (trim(names(n)).ne.'QLCN') .and. &
                    (trim(names(n)).ne.'QILS') .and. &
                    (trim(names(n)).ne.'QICN') .and. &
                    (trim(names(n)).ne.'CLLS') .and. &
                    (trim(names(n)).ne.'CLCN') .and. &
                    (trim(names(n)).ne.'QRAIN') .and. &
                    (trim(names(n)).ne.'QSNOW') .and. &
                    (trim(names(n)).ne.'QGRAUPEL')       ) then

                  if( real(trsum1(n),kind=4).ne.MAPL_UNDEF .and. &
                       real(trsum2(n),kind=4).ne.MAPL_UNDEF       ) then
                     trsum2(n) = real( trsum1(n)/trsum2(n),kind=4)
                  else
                     trsum2(n) = 1.0d0
                  endif
                  ! IF (MAPL_AM_I_ROOT()) print *, trim(names(n)),' ratio is: ',trsum2(n)

                  if( STATE%VARS%TRACER(N)%IS_R4 ) then
                     do k=1,km
                        where( delp(:,:,k).ne.delpold(:,:,k) )
                           state%vars%tracer(n)%content_r4(:,:,k) = state%vars%tracer(n)%content_r4(:,:,k) * trsum2(n)
                        end where
                     enddo
                  else
                     do k=1,km
                        where( delp(:,:,k).ne.delpold(:,:,k) )
                           state%vars%tracer(n)%content   (:,:,k) = state%vars%tracer(n)%content   (:,:,k) * trsum2(n)
                        end where
                     enddo
                  endif
               endif
            enddo
         endif

         deallocate( trsum1 )
         deallocate( trsum2 )

         ! Update Local Copy of QV and OX to account for Global Sum Adjustment
         do k=1,size(names)
            pos = index(names(k),'::')
            if(pos > 0) then
               if( (names(k)(pos+2:))=='OX' ) then
                  if ( ooo%is_r4 ) then
                     ox = ooo%content_r4
                  else
                     ox = ooo%content
                  endif
               endif
            endif
            if( trim(names(k))=='Q'  ) then
               if ( qqq%is_r4 ) then
                  qv = qqq%content_r4
               else
                  qv = qqq%content
               endif
               _ASSERT(all(qv >= 0.0),'After AnaAddIncs: negative or nan water vapor detected')
            endif
         enddo

         ! Diagnostics After Analysis Increments are Added
         call MAPL_GetPointer(export, temp2D, 'DMDTANA', _RC)
         if( associated(temp2D) ) temp2D = ( (vars%pe(:,:,km+1)-vars%pe(:,:,1)) - dmdt )/(grav*dt)

         call getAllWinds(vars%u, vars%v, UC=uc0, VC=vc0, UR=ur, VR=vr)

         dmdt = vars%pe(:,:,km+1)-vars%pe(:,:,1)     ! Psurf-Ptop

         ! DUDTANA
         call MAPL_GetPointer(export, dudtana, 'DUDTANA', _RC)
         if( associated(dudtana) ) then
            dudtana = (ur-dudtana)/dt
         endif

         ! DVDTANA
         call MAPL_GetPointer(export, dvdtana, 'DVDTANA', _RC)
         if( associated(dvdtana) ) then
            dvdtana = (vr-dvdtana)/dt
         endif

         ! DTDTANA
         call MAPL_GetPointer(export, dtdtana, 'DTDTANA', _RC)
         if( associated(dtdtana) ) then
            dtdtana = ((vars%pt*vars%pkz)-dtdtana)/dt
         endif

         ! DDELPDTANA
         call MAPL_GetPointer(export, ddpdtana, 'DDELPDTANA', _RC)
         if( associated(ddpdtana) ) then
            ddpdtana = (delp-ddpdtana)/dt
         endif

         ! DTHVDTANAINT
         ! ------------
         call MAPL_GetPointer(export, temp2D, 'DTHVDTANAINT', _RC)
         if( associated(temp2D) ) then
            tempxy       = vars%pt*(1+eps*qv)   ! Set tempxy = TH*QVnew (After Analysis Update)
            dthdtanaint2 = 0.0
            do k=1,km
               dthdtanaint2 = dthdtanaint2 + tempxy(:,:,k)*delp(:,:,k)
            enddo
            temp2D       = (dthdtanaint2-dthdtanaint1) * MAPL_P00**MAPL_KAPPA / (MAPL_GRAV*DT)
            DEALLOCATE( dthdtanaint1 )
            DEALLOCATE( dthdtanaint2 )
         endif

         ! DQVDTANAINT
         call MAPL_GetPointer(export, temp2D, 'DQVDTANAINT', _RC)
         if( associated(temp2D) ) then
            tempxy       = qv         ! Set tempxy = QNEW (After Analysis Update)
            dqvdtanaint2 = 0.0
            do k=1,km
               dqvdtanaint2 = dqvdtanaint2 + tempxy(:,:,k)*delp(:,:,k)
            enddo
            temp2D       = (dqvdtanaint2-dqvdtanaint1) / (MAPL_GRAV*DT)
            DEALLOCATE( dqvdtanaint1 )
            DEALLOCATE( dqvdtanaint2 )
         endif

         ! DQLDTANAINT
         call MAPL_GetPointer(export, temp2D, 'DQLDTANAINT', _RC)
         if( associated(temp2D) ) then
            dqldtanaint2 = 0.0
            do N = 1,size(names)
               if( trim(names(N)).eq.'QLCN' .or. &
                    trim(names(N)).eq.'QLLS' ) then
                  do k=1,km
                     if( state%vars%tracer(N)%is_r4 ) then
                        dqldtanaint2 = dqldtanaint2 + state%vars%tracer(N)%content_r4(:,:,k)*delp(:,:,k)
                     else
                        dqldtanaint2 = dqldtanaint2 + state%vars%tracer(N)%content   (:,:,k)*delp(:,:,k)
                     endif
                  enddo
               endif
            enddo
            temp2D = (dqldtanaint2-dqldtanaint1) / (MAPL_GRAV*DT)
            DEALLOCATE( dqldtanaint1 )
            DEALLOCATE( dqldtanaint2 )
         endif

         ! DQIDTANAINT
         call MAPL_GetPointer(export, temp2D, 'DQIDTANAINT', _RC)
         if( associated(temp2D) ) then
            dqidtanaint2 = 0.0
            do N = 1,size(names)
               if( trim(names(N)).eq.'QICN' .or. &
                    trim(names(N)).eq.'QILS' ) then
                  do k=1,km
                     if( state%vars%tracer(N)%is_r4 ) then
                        dqidtanaint2 = dqidtanaint2 + state%vars%tracer(N)%content_r4(:,:,k)*delp(:,:,k)
                     else
                        dqidtanaint2 = dqidtanaint2 + state%vars%tracer(N)%content   (:,:,k)*delp(:,:,k)
                     endif
                  enddo
               endif
            enddo
            temp2D = (dqidtanaint2-dqidtanaint1) / (MAPL_GRAV*DT)
            DEALLOCATE( dqidtanaint1 )
            DEALLOCATE( dqidtanaint2 )
         endif

         ! DOXDTANAINT
         call MAPL_GetPointer(export, temp2D, 'DOXDTANAINT', _RC)
         if( associated(temp2D) ) then
            tempxy       = ox         ! Set tempxy = OXnew (After Analysis Update)
            doxdtanaint2 = 0.0
            do k=1,km
               doxdtanaint2 = doxdtanaint2 + tempxy(:,:,k)*delp(:,:,k)
            enddo
            temp2D = (doxdtanaint2-doxdtanaint1) * (MAPL_O3MW/MAPL_AIRMW) / (MAPL_GRAV*DT)
            DEALLOCATE( doxdtanaint1 )
            DEALLOCATE( doxdtanaint2 )
         endif

         DEALLOCATE( delpold)
         DEALLOCATE( qdnew  )
         DEALLOCATE( qdold  )
         DEALLOCATE( qvold  )
         DEALLOCATE( qlold  )
         DEALLOCATE( qiold  )
         DEALLOCATE( qrold  )
         DEALLOCATE( qsold  )
         DEALLOCATE( qgold  )

         ! WMP End ANA section
      else ! REPLAY/ANA is_shutoff

         ox = 0.0
         qv = 0.0
         if (.not. ADIABATIC) then
            do k=1,size(names)
               pos = index(names(k),'::')
               if(pos > 0) then
                  if( (names(k)(pos+2:))=='OX' ) then
                     if ( (ooo%is_r4) .and. associated(ooo%content_r4) ) then
                        if (size(ox)==size(ooo%content_r4)) then
                           ox = ooo%content_r4
                        endif
                     elseif (associated(ooo%content)) then
                        if (size(ox)==size(ooo%content)) then
                           ox = ooo%content
                        endif
                     endif
                  endif
               endif
               if( trim(names(k))=='Q'  ) then
                  if ( (qqq%is_r4) .and. associated(qqq%content_r4) ) then
                     if (size(qv)==size(qqq%content_r4)) then
                        qv = qqq%content_r4
                     endif
                  elseif (associated(qqq%content)) then
                     if (size(qv)==size(qqq%content)) then
                        qv = qqq%content
                     endif
                  endif
                  _ASSERT(all(qv >= 0.0),'DYN_ANA: negative or nan water vapor detected')
               endif
            enddo
         endif
         call getAllWinds(vars%u, vars%v, UC=uc0, VC=vc0, UR=ur, VR=vr)
         delp   = vars%pe(:,:,2:)  -vars%pe(:,:,:km)   ! Pressure Thickness
         dmdt   = vars%pe(:,:,km+1)-vars%pe(:,:,1)     ! Psurf-Ptop
         tempxy = vars%pt * (1.0+eps*qv)
         if (doEnergetics) &
              call Energetics (ur,vr,tempxy,vars%pe,delp,vars%pkz,phisxy,kenrg,penrg,tenrg)

      endif
      if (FV3_STANDALONE == 0) then
         call MAPL_TimerOff(MAPL,"-DYN_ANA")
      endif


      call MAPL_TimerOn(MAPL,"-DYN_PROLOGUE")
      ! Create FV Thermodynamic Variables
      tempxy = vars%pt * vars%pkz      ! Compute Dry Temperature

      ! Initialize Diagnostic Dynamics Tendencies
      dpedt  = vars%pe      ! Edge Pressure      Tendency
      ddpdt  =    delp      ! Pressure Thickness Tendency
      dudt   =     ur       ! U-Wind on A-Grid   Tendency
      dvdt   =     vr       ! V-Wind on A-Grid   Tendency
      dtdt   = tempxy       ! Dry Temperature    Tendency
      dqdt   =     qv       ! Specific Humidity  Tendency
      dthdt  = vars%pt*(1.0+eps*qv)*delp

      call FILLOUT3(export, 'QV_DYN_IN',       qv, _RC)
      call FILLOUT3(export, 'T_DYN_IN',    tempxy, _RC)
      call FILLOUT3(export, 'U_DYN_IN',        ur, _RC)
      call FILLOUT3(export, 'V_DYN_IN',        vr, _RC)
      call FILLOUT3(export, 'PLE_DYN_IN', vars%pe, _RC)

      ! Initialize 3-D Tracer Dynamics Tendencies
      call MAPL_GetPointer( export,dqldt,'DQLDTDYN', _RC)
      call MAPL_GetPointer( export,dqidt,'DQIDTDYN', _RC)
      call MAPL_GetPointer( export,doxdt,'DOXDTDYN', _RC)

      if (allocated(names)) then

         if( associated(dqldt) ) then
            dqldt = 0.0
            do k = 1,size(names)
               if( trim(names(k)).eq.'QLCN' .or. &
                    trim(names(k)).eq.'QLLS' ) then
                  if( state%vars%tracer(k)%is_r4 ) then
                     if (size(dqldt)==size(state%vars%tracer(k)%content_r4)) &
                          dqldt = dqldt - state%vars%tracer(k)%content_r4
                  else
                     if (size(dqldt)==size(state%vars%tracer(k)%content)) &
                          dqldt = dqldt - state%vars%tracer(k)%content
                  endif
               endif
            enddo
         endif

         if( associated(dqidt) ) then
            dqidt = 0.0
            do k = 1,size(names)
               if( trim(names(k)).eq.'QICN' .or. &
                    trim(names(k)).eq.'QILS' ) then
                  if( state%vars%tracer(k)%is_r4 ) then
                     if (size(dqidt)==size(state%vars%tracer(k)%content_r4)) &
                          dqidt = dqidt - state%vars%tracer(k)%content_r4
                  else
                     if (size(dqidt)==size(state%vars%tracer(k)%content)) &
                          dqidt = dqidt - state%vars%tracer(k)%content
                  endif
               endif
            enddo
         endif

         if( associated(doxdt) ) then
            doxdt = 0.0
            do k = 1,size(names)
               pos = index(names(k),'::')
               if(pos > 0) then
                  if( (names(k)(pos+2:))=='OX' ) then
                     if( state%vars%tracer(k)%is_r4 ) then
                        if (size(doxdt)==size(state%vars%tracer(k)%content_r4)) &
                             doxdt = doxdt - state%vars%tracer(k)%content_r4
                     else
                        if (size(doxdt)==size(state%vars%tracer(k)%content)) &
                             doxdt = doxdt - state%vars%tracer(k)%content
                     endif
                  endif
               endif
            enddo
         endif
      endif

      ! Initialize 2-D Vertically Integrated Tracer Dynamics Tendencies
      call MAPL_GetPointer(export, temp2D, 'DQVDTDYNINT', _RC)
      if( associated(temp2D) ) then
         temp2d = 0.0
         do k=1,km
            temp2d = temp2d - qv(:,:,k)*delp(:,:,k)
         enddo
      endif

      call MAPL_GetPointer(export, temp2D, 'DQLDTDYNINT', _RC)
      if( associated(temp2D) ) then
         temp2d = 0.0
         do N = 1,size(names)
            if( trim(names(N)).eq.'QLCN' .or. &
                 trim(names(N)).eq.'QLLS' ) then
               if( state%vars%tracer(N)%is_r4 ) then
                  do k=1,km
                     temp2d = temp2d - state%vars%tracer(N)%content_r4(:,:,k)*delp(:,:,k)
                  enddo
               else
                  do k=1,km
                     temp2d = temp2d - state%vars%tracer(N)%content(:,:,k)*delp(:,:,k)
                  enddo
               endif
            endif
         enddo
      endif

      call MAPL_GetPointer(export, temp2D, 'DQIDTDYNINT', _RC)
      if( associated(temp2D) ) then
         temp2d = 0.0
         do N = 1,size(names)
            if( trim(names(N)).eq.'QICN' .or. &
                 trim(names(N)).eq.'QILS' ) then
               if( state%vars%tracer(N)%is_r4 ) then
                  do k=1,km
                     temp2d = temp2d - state%vars%tracer(N)%content_r4(:,:,k)*delp(:,:,k)
                  enddo
               else
                  do k=1,km
                     temp2d = temp2d - state%vars%tracer(N)%content(:,:,k)*delp(:,:,k)
                  enddo
               endif
            endif
         enddo
      endif

      call MAPL_GetPointer(export, temp2D, 'DOXDTDYNINT', _RC)
      if( associated(temp2D) ) then
         temp2d = 0.0
         do N = 1,size(names)
            pos = index(names(N),'::')
            if(pos > 0) then
               if( (names(N)(pos+2:))=='OX' ) then
                  if( state%vars%tracer(N)%is_r4 ) then
                     do k=1,km
                        temp2d = temp2d - state%vars%tracer(N)%content_r4(:,:,k)*delp(:,:,k)
                     enddo
                  else
                     do k=1,km
                        temp2d = temp2d - state%vars%tracer(N)%content(:,:,k)*delp(:,:,k)
                     enddo
                  endif
               endif
            endif
         enddo
      endif

      ! Compute Energetics After Analysis (and Before Dycore)
      tempxy = vars%pt * (1.0+eps*qv)       ! Compute THV After Analysis Update

      if (doEnergetics) then
         call Energetics(ur, vr, tempxy, vars%pe, delp, vars%pkz, phisxy, kenrg0, penrg0, tenrg0)
         kenrg = (kenrg0-kenrg)/DT
         penrg = (penrg0-penrg)/DT
         tenrg = (tenrg0-tenrg)/DT
         call FILLOUT2(export, 'KEANA', kenrg, _RC)
         call FILLOUT2(export, 'PEANA', penrg, _RC)
         call FILLOUT2(export, 'TEANA', tenrg, _RC)
      endif

      ! Call Wrapper (DynRun) for FVDycore
      call MAPL_GetResource( MAPL, CONSV, 'CONSV:', default=1, _RC)
      call MAPL_GetResource( MAPL,  FILL,  'FILL:', default=0, _RC)

      LCONSV = CONSV.eq.1
      LFILL  =  FILL.eq.1

      ! Get pressures before dynamics
      pe0=vars%pe

      call MAPL_TimerOff(MAPL, "-DYN_PROLOGUE")

      !-------------------------------------------------------

      call MAPL_TimerOn(MAPL, "-DYN_CORE")
      t1 = MPI_Wtime(status)
      call DynRun(STATE, export, clock, gc, PLE0=pe0, _RC)
      t2 = MPI_Wtime(status)
      dyn_run_timer = t2-t1
      call MAPL_TimerOff(MAPL, "-DYN_CORE")

      call MAPL_TimerOn(MAPL, "-DYN_EPILOGUE")
      ! Computational diagnostics
      call MAPL_GetPointer(export, temp2d, 'TIME_IN_DYN', _RC)
      if(associated(temp2d)) temp2d = dyn_run_timer
      call MAPL_GetPointer(export, temp2d, 'PID', _RC)
      if(associated(temp2d)) temp2d = 0 !WMP need to get from MAPL gid

      !#define DEBUG_WINDS
#if defined(DEBUG_WINDS)
      call Write_Profile(grid, vars%u, 'U-after-DynRun')
      call Write_Profile(grid, vars%v, 'V-after-DynRun')
#endif
      plk  = exp( kappa * log( 0.5*(vars%pe(:,:,1:km)+vars%pe(:,:,2:km+1)) ) )
      pkxy = exp( kappa * log( vars%pe ) )

      !----------------------------------------------------------------------------

      if (SW_DYNAMICS) then

         call MAPL_GetPointer(export, temp2d, 'PHIS', _RC)
         if(associated(temp2d)) temp2d = phisxy

         call MAPL_GetPointer(export, temp2d, 'PS', _RC)
         if(associated(temp2d)) temp2d =  vars%pe(:,:,km+1)/GRAV

         call getAllWinds(vars%u, vars%v, UA=ua, VA=va, UC=uc, VC=vc, UR=ur, VR=vr)
         call FILLOUT3(export, 'U_DGRID', vars%u  , _RC)
         call FILLOUT3(export, 'V_DGRID', vars%v  , _RC)
         call FILLOUT3(export, 'U_CGRID', uc      , _RC)
         call FILLOUT3(export, 'V_CGRID', vc      , _RC)
         call FILLOUT3(export, 'U_AGRID', ua      , _RC)
         call FILLOUT3(export, 'V_AGRID', va      , _RC)

         call FILLOUT3(export, 'U'      , ur      , _RC)
         call FILLOUT3(export, 'V'      , vr      , _RC)

      else               ! .not. SW_DYNAMICS

         ! Load Local Variable with Vapor Specific Humidity
         if ((.not. ADIABATIC) .and. (STATE%GRID%NQ > 0)) then
            if ( qqq%is_r4 ) then
               if (size(qv)==size(qqq%content_r4)) qv = qqq%content_r4
            else
               if (size(qv)==size(qqq%content)   ) qv = qqq%content
            endif
            _ASSERT(all(qv >= 0.0),'After DynRun: negative or nan water vapor detected')
         else
            qv = 0.0
         endif

         ! Vertically Integrated THV Tendency Diagnostic
         delp  = ( vars%pe(:,:,2:) - vars%pe(:,:,:km) )
         dthdt = ( vars%pt*(1.0+eps*qv)*delp-dthdt )/dt

         call MAPL_GetPointer(export, temp2d, 'DTHVDTDYNINT', _RC)
         if(associated(temp2d)) then
            qsum1 = 0.0
            do k=1,km
               qsum1 = qsum1 + dthdt(:,:,k)
            enddo
            temp2d = qsum1 * (MAPL_P00**MAPL_KAPPA) / grav
         end if

         ! Compute Dry Theta and T with Unified Poles
         tempxy  = vars%pt * vars%pkz

         ! Compute Mid-Layer Pressure and Pressure Thickness
         delp = ( vars%pe(:,:,2:) - vars%pe(:,:,:km) )
         pl   = ( vars%pe(:,:,2:) + vars%pe(:,:,:km) ) * 0.5

         ! Get all wind derivatives
         call getAllWinds(vars%u, vars%v, UA=ua, VA=va, UC=uc, VC=vc, UR=ur, VR=vr)
         call getVorticity(vars%u, vars%v, vort)
         call getDivergence(uc, vc, divg)

         ! Compute absolute vorticity on the D grid
         call getEPV(vars%pt, vort, ua, va, epvxyz)
         call MAPL_GetPointer(export, temp3D, 'EPV', _RC)
         if(associated(temp3d)) temp3d = epvxyz*(p00**kappa)

         ! Compute Tropopause Pressure, Temperature, and Moisture
         doTropvars=.false.
         call MAPL_GetPointer(export, temp2D, 'TROPP_THERMAL', _RC)
         if(associated(temp2D)) doTropvars=.true.
         call MAPL_GetPointer(export, temp2D, 'TROPP_EPV', _RC)
         if(associated(temp2D)) doTropvars=.true.
         call MAPL_GetPointer(export, temp2D, 'TROPP_BLENDED', _RC)
         if(associated(temp2D)) doTropvars=.true.
         call MAPL_GetPointer(export, temp2D, 'TROPK_BLENDED', _RC)
         if(associated(temp2D)) doTropvars=.true.
         call MAPL_GetPointer(export, temp2D, 'TROPT', _RC)
         if(associated(temp2D)) doTropvars=.true.
         call MAPL_GetPointer(export, temp2D, 'TROPQ', _RC)
         if(associated(temp2D)) doTropvars=.true.

         if (doTropvars) then
            allocate( tropp1 (ifirstxy:ilastxy,jfirstxy:jlastxy) )
            allocate( tropp2 (ifirstxy:ilastxy,jfirstxy:jlastxy) )
            allocate( tropp3 (ifirstxy:ilastxy,jfirstxy:jlastxy) )
            allocate( tropt  (ifirstxy:ilastxy,jfirstxy:jlastxy) )
            allocate( tropq  (ifirstxy:ilastxy,jfirstxy:jlastxy) )
            call tropovars( &
                 ilastxy-ifirstxy+1,jlastxy-jfirstxy+1,km, &
                 real(vars%pe            ,kind=4),         &
                 real(pl                 ,kind=4),         &
                 real(tempxy             ,kind=4),         &
                 real(qv                 ,kind=4),         &
                 real(epvxyz*(p00**kappa),kind=4),         &
                 tropp1,tropp2,tropp3,tropt,tropq          )

            ! get blended index
            call MAPL_GetPointer(export, temp2D, 'TROPK_BLENDED', _RC)
            if( associated(temp2D) ) then
               kend = km
               do j=jfirstxy,jlastxy
                  do i=ifirstxy,ilastxy
                     if (tropp3(i,j) .NE. MAPL_UNDEF) then
                        kend = 1
                        do while (vars%pe(i,j,kend).LE.tropp3(i,j))
                           kend = kend+1
                        enddo
                     else
                        kend = 1
                        do while (vars%pe(i,j,kend).LE.40000.0)
                           kend = kend+1
                        enddo
                     endif
                     temp2D(i-ifirstxy+1,j-jfirstxy+1) = kend
                  enddo
               enddo
            endif

            call MAPL_GetPointer(export, temp2D, 'TROPP_THERMAL', _RC)
            if(associated(temp2D)) temp2D = tropp1

            call MAPL_GetPointer(export, temp2D, 'TROPP_EPV', _RC)
            if(associated(temp2D)) temp2D = tropp2

            call MAPL_GetPointer(export, temp2D, 'TROPP_BLENDED', _RC)
            if(associated(temp2D)) temp2D = tropp3

            call MAPL_GetPointer(export, temp2D, 'TROPT', _RC)
            if(associated(temp2D)) temp2D = tropt

            call MAPL_GetPointer(export, temp2D, 'TROPQ', _RC)
            if(associated(temp2D)) temp2D = tropq

            deallocate( tropp1 )
            deallocate( tropp2 )
            deallocate( tropp3 )
            deallocate( tropt  )
            deallocate( tropq  )
         endif

         ! Get Cubed-Sphere Wind Exports
         call FILLOUT3(export, 'U_DGRID', vars%u, _RC)
         call FILLOUT3(export, 'V_DGRID', vars%v, _RC)
         call FILLOUT3(export, 'U_CGRID', uc    , _RC)
         call FILLOUT3(export, 'V_CGRID', vc    , _RC)
         call FILLOUT3(export, 'U_AGRID', ua    , _RC)
         call FILLOUT3(export, 'V_AGRID', va    , _RC)

         if (DEBUG_DYN) then
            call MAPL_MaxMin('DYN: Q_AF_DYN ', qv)
            call MAPL_MaxMin('DYN: T_AF_DYN ', tempxy)
            call MAPL_MaxMin('DYN: U_AF_DYN ', ua)
            call MAPL_MaxMin('DYN: V_AF_DYN ', va)
         endif

         ! Compute Diagnostic Dynamics Tendencies
         !  (Note: initial values of d(m,u,v,T,q)/dt are progs m,u,v,T,q)
         dmdt = ( vars%pe(:,:,km+1)-vars%pe(:,:,1) - dmdt )/(grav*dt)

         dudt = (    ur-dudt )/dt
         dvdt = (    vr-dvdt )/dt
         dtdt = (  tempxy-dtdt )/dt
         dqdt = (      qv-dqdt )/dt

         dpedt = ( vars%pe - dpedt )/dt
         ddpdt = ( delp - ddpdt )/dt ! Pressure Thickness Tendency


         call FILLOUT3(export, 'DELP'      , delp , _RC)
         call FILLOUT3(export, 'DUDTDYN'   , dudt , _RC)
         call FILLOUT3(export, 'DVDTDYN'   , dvdt , _RC)
         call FILLOUT3(export, 'DTDTDYN'   , dtdt , _RC)
         call FILLOUT3(export, 'DQVDTDYN'  , dqdt , _RC)
         call FILLOUT3(export, 'DDELPDTDYN', ddpdt, _RC)
         call FILLOUT3(export, 'DPLEDTDYN' , dpedt, _RC)

         ! fill pressure exports (PLE0: Before) & (PLE1: After) from FV3
         call FILLOUT3r8(export, 'PLE0', pe0, _RC)
         pe1=vars%pe
         call FILLOUT3r8(export, 'PLE1', pe1, _RC)

         if (AdvCore_Advection==2) then
            ! Compute time-centered C-Grid Courant Numbers and Mass Fluxes on Cubed Orientation
            uc0 = 0.5*(uc +uc0)
            vc0 = 0.5*(vc +vc0)
            pe0 = 0.5*(pe1+pe0)
            call computeMassFluxes(uc0, vc0, pe0, mfxxyz, mfyxyz, cxxyz, cyxyz, dt)
            call FILLOUT3r8(export, 'CX'  , cxxyz  , _RC)
            call FILLOUT3r8(export, 'CY'  , cyxyz  , _RC)
            call FILLOUT3r8(export, 'MFX' , mfxxyz , _RC)
            call FILLOUT3r8(export, 'MFY' , mfyxyz , _RC)
         else
            ! Fill Advection C-Grid Courant Numbers and Mass Fluxes on Cubed Orientation from FV3 DynCore
            call fillMassFluxes(mfxxyz, mfyxyz, cxxyz, cyxyz)
            call FILLOUT3r8(export, 'CX'  , cxxyz  , _RC)
            call FILLOUT3r8(export, 'CY'  , cyxyz  , _RC)
            call FILLOUT3r8(export, 'MFX' , mfxxyz , _RC)
            call FILLOUT3r8(export, 'MFY' , mfyxyz , _RC)
         endif

         call FILLOUT3(export, 'CU' ,  cxxyz , _RC)
         call FILLOUT3(export, 'CV' ,  cyxyz , _RC)
         call FILLOUT3(export, 'MX' , mfxxyz , _RC)
         call FILLOUT3(export, 'MY' , mfyxyz , _RC)

         ! Compute and return the vertical mass flux
         call getVerticalMassFlux(mfxxyz, mfyxyz, mfzxyz, dt)
         call FILLOUT3r8(export, 'MFZ' , mfzxyz , _RC)
         call FILLOUT3(export, 'MZ' , mfzxyz , _RC)

         call FILLOUT3(export, 'U'      , ur      , _RC)
         call FILLOUT3(export, 'V'      , vr      , _RC)
         call FILLOUT3(export, 'T'      , tempxy  , _RC)
         call FILLOUT3(export, 'Q'      , qv      , _RC)
         call FILLOUT3(export, 'PL'     , pl      , _RC)
         call FILLOUT3(export, 'PLE'    , vars%pe , _RC)
         call FILLOUT3(export, 'PLK'    , plk     , _RC)
         call FILLOUT3(export, 'PKE'    , pkxy    , _RC)
         call FILLOUT3(export, 'PT'     , vars%pt , _RC)
         call FILLOUT3(export, 'PE'     , vars%pe , _RC)

#ifdef SKIP_TRACERS
         do itracer = 1, ntracers
            write(myTracer, "('Q',i5.5)") itracer-1
            call MAPL_GetPointer(export, temp3D, TRIM(myTracer), _RC)
            if((associated(temp3d)) .and. (NQ>=itracer)) then
               if (state%vars%tracer(itracer)%is_r4) then
                  temp3d = state%vars%tracer(itracer)%content_r4
               else
                  temp3d = state%vars%tracer(itracer)%content
               endif
            endif
         enddo
#endif

         call MAPL_GetPointer(export, temp3D, 'PV', _RC)
         if(associated(temp3d)) temp3d = epvxyz/vars%pt

         call MAPL_GetPointer(export, temp3D, 'S', _RC)
         if(associated(temp3d)) temp3d = tempxy*cp

         call MAPL_GetPointer(export, temp3d, 'TH', _RC)
         !   if(associated(temp3d)) temp3d = vars%pt*(p00**kappa)
         if(associated(temp3d)) then
            temp3d = (tempxy)*(p00/(0.5*(vars%pe(:,:,1:km)+vars%pe(:,:,2:km+1))))**kappa
         end if

         call MAPL_GetPointer(export, temp2d, 'DMDTDYN', _RC)
         if(associated(temp2d)) temp2d = dmdt

         ! Compute 3-D Tracer Dynamics Tendencies
         call MAPL_GetPointer(export, qctmp, 'QC', _RC)

         if( associated(qctmp) ) then
            qctmp = 0.0
            do k = 1,size(names)
               if( trim(names(k)).eq.'QLCN' .or. &
                    trim(names(k)).eq.'QILS' .or. &
                    trim(names(k)).eq.'QICN' .or. &
                    trim(names(k)).eq.'QLLS' ) then
                  if( state%vars%tracer(k)%is_r4 ) then
                     if (size(dqldt)==size(state%vars%tracer(k)%content_r4)) &
                          qctmp = qctmp + state%vars%tracer(k)%content_r4
                  else
                     if (size(dqldt)==size(state%vars%tracer(k)%content)) &
                          qctmp = qctmp + state%vars%tracer(k)%content
                  endif
               endif
            enddo
         endif

         if( associated(dqldt) ) then
            do N = 1,size(names)
               if( trim(names(N)).eq.'QLCN' .or. &
                    trim(names(N)).eq.'QLLS' ) then
                  if( state%vars%tracer(N)%is_r4 ) then
                     dqldt = dqldt + state%vars%tracer(N)%content_r4
                  else
                     dqldt = dqldt + state%vars%tracer(N)%content
                  endif
               endif
            enddo
            dqldt = dqldt/dt
         endif

         if( associated(dqidt) ) then
            do N = 1,size(names)
               if( trim(names(N)).eq.'QICN' .or. &
                    trim(names(N)).eq.'QILS' ) then
                  if( state%vars%tracer(N)%is_r4 ) then
                     dqidt = dqidt + state%vars%tracer(N)%content_r4
                  else
                     dqidt = dqidt + state%vars%tracer(N)%content
                  endif
               endif
            enddo
            dqidt = dqidt/dt
         endif

         if( associated(doxdt) ) then
            do N = 1,size(names)
               pos = index(names(N),'::')
               if(pos > 0) then
                  if( (names(N)(pos+2:))=='OX' ) then
                     if( state%vars%tracer(N)%is_r4 ) then
                        doxdt = doxdt + state%vars%tracer(N)%content_r4
                     else
                        doxdt = doxdt + state%vars%tracer(N)%content
                     endif
                  endif
               endif
            enddo
            doxdt = doxdt/dt
         endif

         ! Compute 2-D Vertically Integrated Tracer Dynamics Tendencies
         call MAPL_GetPointer(export, temp2D, 'DQVDTDYNINT', _RC)
         if( associated(temp2D) ) then
            do k=1,km
               temp2d = temp2d + qv(:,:,k)*delp(:,:,k)
            enddo
            temp2d = temp2d/(grav*dt)
         endif

         call MAPL_GetPointer(export, temp2D, 'DQLDTDYNINT', _RC)
         if( associated(temp2D) ) then
            do N = 1,size(names)
               if( trim(names(N)).eq.'QLCN' .or. &
                    trim(names(N)).eq.'QLLS' ) then
                  if( state%vars%tracer(N)%is_r4 ) then
                     do k=1,km
                        temp2d = temp2d + state%vars%tracer(N)%content_r4(:,:,k)*delp(:,:,k)
                     enddo
                  else
                     do k=1,km
                        temp2d = temp2d + state%vars%tracer(N)%content(:,:,k)*delp(:,:,k)
                     enddo
                  endif
               endif
            enddo
            temp2d = temp2d/(grav*dt)
         endif

         call MAPL_GetPointer(export, temp2D, 'DQIDTDYNINT', _RC)
         if( associated(temp2D) ) then
            do N = 1,size(names)
               if( trim(names(N)).eq.'QICN' .or. &
                    trim(names(N)).eq.'QILS' ) then
                  if( state%vars%tracer(N)%is_r4 ) then
                     do k=1,km
                        temp2d = temp2d + state%vars%tracer(N)%content_r4(:,:,k)*delp(:,:,k)
                     enddo
                  else
                     do k=1,km
                        temp2d = temp2d + state%vars%tracer(N)%content(:,:,k)*delp(:,:,k)
                     enddo
                  endif
               endif
            enddo
            temp2d = temp2d/(grav*dt)
         endif

         call MAPL_GetPointer(export, temp2D, 'DOXDTDYNINT', _RC)
         if( associated(temp2D) ) then
            do N = 1,size(names)
               pos = index(names(N),'::')
               if(pos > 0) then
                  if( (names(N)(pos+2:))=='OX' ) then
                     if( state%vars%tracer(N)%is_r4 ) then
                        do k=1,km
                           temp2d = temp2d + state%vars%tracer(N)%content_r4(:,:,k)*delp(:,:,k)
                        enddo
                     else
                        do k=1,km
                           temp2d = temp2d + state%vars%tracer(N)%content(:,:,k)*delp(:,:,k)
                        enddo
                     endif
                  endif
               endif
            enddo
            temp2d = temp2d * (MAPL_O3MW/MAPL_AIRMW) / (MAPL_GRAV*DT)
         endif

         ! Virtual temperature
         tempxy =  tempxy*(1.0+eps*qv)

         call MAPL_GetPointer(export, temp3D, 'TV', _RC)
         if(associated(temp3D)) temp3D = tempxy

         ! Fluxes: UCPT & VCPT
         !--------------------
         call MAPL_GetPointer(export, temp2d, 'UCPT', _RC)
         if(associated(temp2d)) then
            temp2d = 0.0
            do k=1,km
               temp2d = temp2d + ur(:,:,k)*tempxy(:,:,k)*delp(:,:,k)
            enddo
            temp2d = temp2d*(cp/grav)
         end if

         call MAPL_GetPointer(export, temp2d, 'VCPT', _RC)
         if(associated(temp2d)) then
            temp2d = 0.0
            do k=1,km
               temp2d = temp2d + vr(:,:,k)*tempxy(:,:,k)*delp(:,:,k)
            enddo
            temp2d = temp2d*(cp/grav)
         end if

         ! Compute Energetics After Dycore
         tempxy = vars%pt*(1.0+eps*qv)  ! Convert TH to THV

         call MAPL_GetPointer(export, temp3d, 'THV', _RC)
         if(associated(temp3d)) temp3d = tempxy

         if (doEnergetics) then
            call Energetics(ur, vr, tempxy, vars%pe, delp, vars%pkz, phisxy, kenrg, penrg, tenrg)
            kedyn   = (kenrg -kenrg0)/DT
            pedyn   = (penrg -penrg0)/DT
            tedyn   = (tenrg -tenrg0)/DT
            call MAPL_GetPointer(export, temp2d, 'KEDYN', _RC)
            if(associated(temp2d)) temp2d = kedyn
            call MAPL_GetPointer(export, temp2d, 'PEDYN', _RC)
            if(associated(temp2d)) temp2d = pedyn
            call MAPL_GetPointer(export, temp2d, 'TEDYN', _RC)
            if(associated(temp2d)) temp2d = tedyn
         endif

         ! Compute/Get Omega
         call getOmega(omaxyz)

         ! Fluxes: UKE & VKE
         call MAPL_GetPointer(export, tempu, 'UKE', _RC)
         call MAPL_GetPointer(export, tempv, 'VKE', _RC)

         if(associated(tempu) .or. associated(tempv)) then
            tmp3d = 0.5*(ur**2 + vr**2)
         end if

         if(associated(tempu)) then
            tempu = 0.0
            do k=1,km
               tempu = tempu + ur(:,:,k)*tmp3d(:,:,k)*delp(:,:,k)
            enddo
            tempu = tempu / grav
         end if

         if(associated(tempv)) then
            tempv = 0.0
            do k=1,km
               tempv = tempv + vr(:,:,k)*tmp3d(:,:,k)*delp(:,:,k)
            enddo
            tempv = tempv / grav
         end if

         ! Fluxes: UQV & VQV
         call MAPL_GetPointer(export, temp2d, 'UQV', _RC)
         if(associated(temp2d)) then
            temp2d = 0.0
            do k=1,km
               temp2d = temp2d + ur(:,:,k)*QV(:,:,k)*delp(:,:,k)
            enddo
            temp2d = temp2d / grav
         end if

         call MAPL_GetPointer(export, temp2d, 'VQV', _RC)
         if(associated(temp2d)) then
            temp2d = 0.0
            do k=1,km
               temp2d = temp2d + vr(:,:,k)*QV(:,:,k)*delp(:,:,k)
            enddo
            temp2d = temp2d / grav
         end if

         ! Fluxes: UQL & VQL
         call MAPL_GetPointer(export, temp2d, 'UQL', _RC)
         if(associated(temp2d)) then
            temp2d = 0.0
            do N = 1,size(names)
               if( trim(names(n)).eq.'QLCN' .or. &
                    trim(names(n)).eq.'QLLS' ) then
                  do k=1,km
                     if( state%vars%tracer(n)%is_r4 ) then
                        temp2d = temp2d + ur(:,:,k)*state%vars%tracer(n)%content_r4(:,:,k)*delp(:,:,k)
                     else
                        temp2d = temp2d + ur(:,:,k)*state%vars%tracer(n)%content   (:,:,k)*delp(:,:,k)
                     endif
                  enddo
               endif
            enddo
            temp2d = temp2d / grav
         end if

         call MAPL_GetPointer(export, temp2d, 'VQL', _RC)
         if(associated(temp2d)) then
            temp2d = 0.0
            do N = 1,size(names)
               if( trim(names(n)).eq.'QLCN' .or. &
                    trim(names(n)).eq.'QLLS' ) then
                  do k=1,km
                     if( state%vars%tracer(n)%is_r4 ) then
                        temp2d = temp2d + vr(:,:,k)*state%vars%tracer(n)%content_r4(:,:,k)*delp(:,:,k)
                     else
                        temp2d = temp2d + vr(:,:,k)*state%vars%tracer(n)%content   (:,:,k)*delp(:,:,k)
                     endif
                  enddo
               endif
            enddo
            temp2d = temp2d / grav
         end if

         ! Fluxes: UQI & VQI
         call MAPL_GetPointer(export, temp2d, 'UQI', _RC)
         if(associated(temp2d)) then
            temp2d = 0.0
            do N = 1,size(names)
               if( trim(names(n)).eq.'QICN' .or. &
                    trim(names(n)).eq.'QILS' ) then
                  do k=1,km
                     if( state%vars%tracer(n)%is_r4 ) then
                        temp2d = temp2d + ur(:,:,k)*state%vars%tracer(n)%content_r4(:,:,k)*delp(:,:,k)
                     else
                        temp2d = temp2d + ur(:,:,k)*state%vars%tracer(n)%content   (:,:,k)*delp(:,:,k)
                     endif
                  enddo
               endif
            enddo
            temp2d = temp2d / grav
         end if

         call MAPL_GetPointer(export, temp2d, 'VQI', _RC)
         if(associated(temp2d)) then
            temp2d = 0.0
            do N = 1,size(names)
               if( trim(names(n)).eq.'QICN' .or. &
                    trim(names(n)).eq.'QILS' ) then
                  do k=1,km
                     if( state%vars%tracer(n)%is_r4 ) then
                        temp2d = temp2d + vr(:,:,k)*state%vars%tracer(n)%content_r4(:,:,k)*delp(:,:,k)
                     else
                        temp2d = temp2d + vr(:,:,k)*state%vars%tracer(n)%content   (:,:,k)*delp(:,:,k)
                     endif
                  enddo
               endif
            enddo
            temp2d = temp2d / grav
         end if

         ! Height related diagnostics
         zle(:,:,km+1) = phisxy(:,:)
         do k=km,1,-1
            zle(:,:,k) = zle(:,:,k+1) + cp*tempxy(:,:,k)*( pkxy(:,:,k+1)-pkxy(:,:,k) )
         enddo
         zle = zle/grav

         call MAPL_GetPointer(export, temp3d, 'ZLE', _RC)
         if(associated(temp3d)) temp3d = zle

         call MAPL_GetPointer(export, temp3d, 'ZL', _RC)
         if(associated(temp3d)) temp3d = 0.5*( zle(:,:,:km)+zle(:,:,2:) )

         call MAPL_GetPointer(export, temp3d, 'S', _RC)
         if(associated(temp3d)) temp3d = temp3d + grav*(0.5*( zle(:,:,:km)+zle(:,:,2:) ))

         ! Fluxes: UPHI & VPHI
         call MAPL_GetPointer(export, tempu, 'UPHI', _RC)
         call MAPL_GetPointer(export, tempv, 'VPHI', _RC)

         if( associated(tempu).or.associated(tempv) ) zl = 0.5*( zle(:,:,:km)+zle(:,:,2:) )

         if(associated(tempu)) then
            tempu = 0.0
            do k=1,km
               tempu = tempu + ur(:,:,k)*zl(:,:,k)*delp(:,:,k)
            enddo
         end if

         if(associated(tempv)) then
            tempv = 0.0
            do k=1,km
               tempv = tempv + vr(:,:,k)*zl(:,:,k)*delp(:,:,k)
            enddo
         end if

         ! Fill Surface and Near-Surface Variables
         HGT_SURFACE = 50.0
         if (km .eq. 72) HGT_SURFACE =  0.0
         call MAPL_GetResource(MAPL, HGT_SURFACE, label="HGT_SURFACE:", default=HGT_SURFACE, _RC)
         if ( HGT_SURFACE .gt. 0.0 ) then
            ! Near surface height for surface
            call MAPL_GetPointer(export, temp2d, 'DZ', _RC)
            if(associated(temp2d)) temp2d = HGT_SURFACE

            ! Get the height above the surface
            do k=1,km+1
               zle(:,:,k) = zle(:,:,k) - zle(:,:,km+1)
            enddo

            call MAPL_GetPointer(export, temp2d, 'PS', _RC)
            if(associated(temp2d)) temp2d =  vars%pe(:,:,km+1)

            call MAPL_GetPointer(export, temp2d, 'US', _RC)
            if(associated(temp2d)) then
               call VertInterp(temp2d, ur, -zle, -HGT_SURFACE, _RC)
            end if

            call MAPL_GetPointer(export, temp2d, 'VS', _RC)
            if(associated(temp2d)) then
               call VertInterp(temp2d, vr, -zle, -HGT_SURFACE, _RC)
            end if

            call MAPL_GetPointer(export, temp2d, 'TA', _RC)
            if(associated(temp2d)) then
               tempxy  = vars%pt * vars%pkz
               call VertInterp(temp2d, tempxy, -zle, -HGT_SURFACE, _RC)
            end if

            call MAPL_GetPointer(export, temp2d, 'QA', _RC)
            if(associated(temp2d)) then
               call VertInterp(temp2d, qv, -zle, -HGT_SURFACE, positive_definite=.true., _RC)
            end if

            call MAPL_GetPointer(export, temp2d, 'SPEED', _RC)
            if(associated(temp2d)) then
               call VertInterp(temp2d, sqrt(ur**2 + vr**2), -zle, -HGT_SURFACE, _RC)
            end if
         else
            ! Fill Surface with Lowest Model Level Variables
            call MAPL_GetPointer(export, temp2d, 'DZ', _RC)
            if(associated(temp2d)) temp2d = 0.5*( zle(:,:,km)-zle(:,:,km+1) )

            call MAPL_GetPointer(export, temp2d, 'PS', _RC)
            if(associated(temp2d)) temp2d = vars%pe(:,:,km+1)

            call MAPL_GetPointer(export, temp2d, 'US', _RC)
            if(associated(temp2d)) temp2d = ur(:,:,km)

            call MAPL_GetPointer(export, temp2d, 'VS', _RC)
            if(associated(temp2d)) temp2d = vr(:,:,km)

            call MAPL_GetPointer(export, temp2d, 'TA', _RC)
            if(associated(temp2d)) then
               tempxy = vars%pt * vars%pkz
               temp2d = tempxy(:,:,km)
            endif

            call MAPL_GetPointer(export, temp2d, 'QA', _RC)
            if(associated(temp2d)) temp2d = qv(:,:,km)

            call MAPL_GetPointer(export, temp2d, 'SPEED', _RC)
            if(associated(temp2d)) temp2d = sqrt( ur(:,:,km)**2 + vr(:,:,km)**2 )
         endif


         call MAPL_GetPointer(export, temp2d, 'WSPD_10M', _RC)
         if(associated(temp2d)) then
            call VertInterp(temp2d, sqrt(ur**2 + vr**2), -zle, -10.0, _RC)
         end if

         if (.not. HYDROSTATIC) then
            call MAPL_GetPointer(export, temp2d, 'VVEL_UP_100_1000', _RC)
            if(associated(temp2d)) then
               temp2d = vars%w(ifirstxy:ilastxy,jfirstxy:jlastxy,km)
               do k=km-1,1,-1
                  do j=jfirstxy,jlastxy
                     do i=ifirstxy,ilastxy
                        if ( (vars%w(i,j,k) > temp2d(i-ifirstxy+1,j-jfirstxy+1)) .and. &
                             (vars%pe(i,j,k) >= 10000.0) ) then
                           temp2d(i-ifirstxy+1,j-jfirstxy+1) = vars%w(i,j,k)
                        endif
                     enddo
                  enddo
               enddo
            end if
            call MAPL_GetPointer(export, temp2d, 'VVEL_DN_100_1000', _RC)
            if(associated(temp2d)) then
               temp2d = vars%w(ifirstxy:ilastxy,jfirstxy:jlastxy,km)
               do k=km-1,1,-1
                  do j=jfirstxy,jlastxy
                     do i=ifirstxy,ilastxy
                        if ( (vars%w(i,j,k) < temp2d(i-ifirstxy+1,j-jfirstxy+1)) .and. &
                             (vars%pe(i,j,k) >= 10000.0) ) then
                           temp2d(i-ifirstxy+1,j-jfirstxy+1) = vars%w(i,j,k)
                        endif
                     enddo
                  enddo
               enddo
            end if
         end if

         ! Updraft Helicty Exports
         call MAPL_GetPointer(export,  uh25, 'UH25', ALLOC=.TRUE., _RC)
         call MAPL_GetPointer(export,  uh03, 'UH03', ALLOC=.TRUE., _RC)
         call MAPL_GetPointer(export, srh01,'SRH01', ALLOC=.TRUE., _RC)
         call MAPL_GetPointer(export, srh03,'SRH03', ALLOC=.TRUE., _RC)
         call MAPL_GetPointer(export, srh25,'SRH25', ALLOC=.TRUE., _RC)
         ! Per WMP, this calculation is not useful if running hydrostatic
         if (.not. HYDROSTATIC) then
            if( associated( uh25) .or. associated( uh03) .or. &
                 associated(srh01) .or. associated(srh03) .or. associated(srh25) ) then
               call fv_getUpdraftHelicity(uh25, uh03, srh01, srh03, srh25)
            endif
         endif

         ! Divergence Exports
         logpe = log(vars%pe)

         call MAPL_GetPointer(export, temp3d, 'DIVG', _RC)
         if(associated(temp3d)) temp3d = divg

         call MAPL_GetPointer(export, temp2d, 'DIVG200', _RC)
         if(associated(temp2d)) then
            call VertInterp(temp2d, dble(divg), logpe, log(20000.), _RC)
         end if

         call MAPL_GetPointer(export, temp2d, 'DIVG500', _RC)
         if(associated(temp2d)) then
            call VertInterp(temp2d, dble(divg), logpe, log(50000.), _RC)
         end if

         call MAPL_GetPointer(export, temp2d, 'DIVG700', _RC)
         if(associated(temp2d)) then
            call VertInterp(temp2d, dble(divg), logpe, log(70000.), _RC)
         end if

         call MAPL_GetPointer(export, temp2d, 'DIVG850', _RC)
         if(associated(temp2d)) then
            call VertInterp(temp2d, dble(divg), logpe, log(85000.), _RC)
         end if

         ! Vorticity Exports
         call MAPL_GetPointer(export, temp3d, 'VORT', _RC)
         if(associated(temp3d)) temp3d = vort

         call MAPL_GetPointer(export, temp2d, 'VORT200', _RC)
         if(associated(temp2d)) then
            call VertInterp(temp2d, dble(vort), logpe, log(20000.), _RC)
         end if

         call MAPL_GetPointer(export, temp2d, 'VORT500', _RC)
         if(associated(temp2d)) then
            call VertInterp(temp2d, dble(vort), logpe, log(50000.), _RC)
         end if

         call MAPL_GetPointer(export, temp2d, 'VORT700', _RC)
         if(associated(temp2d)) then
            call VertInterp(temp2d, dble(vort), logpe, log(70000.), _RC)
         end if

         call MAPL_GetPointer(export, temp2d, 'VORT850', _RC)
         if(associated(temp2d)) then
            call VertInterp(temp2d, dble(vort), logpe, log(85000.), _RC)
         end if

         ! Vertical Velocity Exports
         call FILLOUT3(export, 'OMEGA', omaxyz, _RC)

         call MAPL_GetPointer(export, temp2d, 'OMEGA850', _RC)
         if(associated(temp2d)) then
            call VertInterp(temp2d, omaxyz, logpe, log(85000.), _RC)
         end if

         call MAPL_GetPointer(export, temp2d, 'OMEGA700', _RC)
         if(associated(temp2d)) then
            call VertInterp(temp2d, omaxyz, logpe, log(70000.), _RC)
         end if

         call MAPL_GetPointer(export, temp2d, 'OMEGA500', _RC)
         if(associated(temp2d)) then
            call VertInterp(temp2d, omaxyz, logpe, log(50000.), _RC)
         end if

         call MAPL_GetPointer(export, temp2d, 'OMEGA200', _RC)
         if(associated(temp2d)) then
            call VertInterp(temp2d, omaxyz, logpe, log(20000.), _RC)
         end if

         call MAPL_GetPointer(export, temp2d, 'OMEGA10', _RC)
         if(associated(temp2d)) then
            call VertInterp(temp2d, omaxyz, logpe, log(1000.), _RC)
         end if

         if (.not. HYDROSTATIC) then

            call FILLOUT3(export, 'DELZ', vars%dz(ifirstxy:ilastxy,jfirstxy:jlastxy,:), _RC)
            call FILLOUT3(export, 'W', vars%w(ifirstxy:ilastxy,jfirstxy:jlastxy,:), _RC)

            call MAPL_GetPointer(export, temp2d, 'W850', _RC)
            if(associated(temp2d)) then
               call VertInterp(temp2d, vars%w(ifirstxy:ilastxy,jfirstxy:jlastxy,:), logpe, log(85000.), _RC)
            end if

            call MAPL_GetPointer(export, temp2d, 'W500', _RC)
            if(associated(temp2d)) then
               call VertInterp(temp2d, vars%w(ifirstxy:ilastxy,jfirstxy:jlastxy,:), logpe, log(50000.), _RC)
            end if

            call MAPL_GetPointer(export, temp2d, 'W200', _RC)
            if(associated(temp2d)) then
               call VertInterp(temp2d, vars%w(ifirstxy:ilastxy,jfirstxy:jlastxy,:), logpe, log(20000.), _RC)
            end if

            call MAPL_GetPointer(export, temp2d, 'W10', _RC)
            if(associated(temp2d)) then
               call VertInterp(temp2d, vars%w(ifirstxy:ilastxy,jfirstxy:jlastxy,:), logpe, log(1000.), _RC)
            end if
         endif

      end if   ! SW_DYNAMICS

      call MAPL_TimerOff(MAPL, "-DYN_EPILOGUE")

      ! De-Allocate Arrays

      if (doEnergetics) then
         deallocate( KEDYN  )
         deallocate( PEDYN  )
         deallocate( TEDYN  )
         deallocate( KENRG  )
         deallocate( PENRG  )
         deallocate( TENRG  )
         deallocate( KENRG0 )
         deallocate( PENRG0 )
         deallocate( TENRG0 )
      endif

      deallocate( qsum1 )
      deallocate( qsum2 )

      deallocate( zl     )
      deallocate( zle    )
      deallocate( logpe  )
      deallocate( plk    )
      deallocate( pkxy   )
      deallocate( vort   )
      deallocate( divg   )
      deallocate( tmp3d  )
      deallocate( omaxyz )
      deallocate( epvxyz )
      deallocate(  cxxyz )
      deallocate(  cyxyz )
      deallocate( mfxxyz )
      deallocate( mfyxyz )
      deallocate( mfzxyz )
      deallocate( tempxy )
      deallocate( pe0    )
      deallocate( pe1    )
      deallocate( pl     )
      deallocate( ua     )
      deallocate( va     )
      deallocate( uc     )
      deallocate( vc     )
      deallocate( uc0    )
      deallocate( vc0    )
      deallocate( ur     )
      deallocate( vr     )
      deallocate( qv     )
      deallocate( ql     )
      deallocate( qi     )
      deallocate( qr     )
      deallocate( qs     )
      deallocate( qg     )
      deallocate( ox     )
      deallocate( delp   )
      deallocate( dmdt   )
      deallocate( dudt   )
      deallocate( dvdt   )
      deallocate( dtdt   )
      deallocate( dqdt   )
      deallocate( dthdt  )
      deallocate( dpedt  )
      deallocate( ddpdt  )
      deallocate( phisxy )
      if (allocated(names)) deallocate( names  )
      if (allocated(names0)) deallocate( names0  )

      call freeTracers(state)

      call MAPL_TimerOff(MAPL, "RUN")
      call MAPL_TimerOff(MAPL, "TOTAL")

      !if (ADIABATIC) then
      !  ! Fill Exports
      !   call RunAddIncs(gc, import, export, clock, rc)
      !endif

      _RETURN(_SUCCESS)

   contains

      subroutine check_replay_time_(lring)

         logical :: lring

         integer :: REPLAY_REF_DATE, REPLAY_REF_TIME, REPLAY_REF_TGAP
         integer :: REF_TIME(6), REF_TGAP(6)
         type (ESMF_TimeInterval)  :: RefTGap

         call MAPL_GetResource(MAPL, ReplayType, 'REPLAY_TYPE:', default="FULL", _RC)
         !  if (trim(ReplayType) == "FULL") return

         call MAPL_GetResource(MAPL, REPLAY_REF_DATE, label = 'REPLAY_REF_DATE:', default=-1, _RC)
         call MAPL_GetResource(MAPL, REPLAY_REF_TIME, label = 'REPLAY_REF_TIME:', default=-1, _RC)
         call MAPL_GetResource(MAPL, REPLAY_REF_TGAP, label = 'REPLAY_REF_TGAP:', default=-1, _RC)

         if(REPLAY_REF_DATE==-1.or.REPLAY_REF_TIME==-1) return

         REF_TIME(1) =     REPLAY_REF_DATE/10000
         REF_TIME(2) = mod(REPLAY_REF_DATE,10000)/100
         REF_TIME(3) = mod(REPLAY_REF_DATE,100)
         REF_TIME(4) =     REPLAY_REF_TIME/10000
         REF_TIME(5) = mod(REPLAY_REF_TIME,10000)/100
         REF_TIME(6) = mod(REPLAY_REF_TIME,100)

         ! set replay time
         call ESMF_TimeSet(RefTime, &
              YY =  REF_TIME(1), &
              MM =  REF_TIME(2), &
              DD =  REF_TIME(3), &
              H  =  REF_TIME(4), &
              M  =  REF_TIME(5), &
              S  =  REF_TIME(6), _RC)
         if (REPLAY_REF_TGAP>0) then
            REF_TGAP    = 0
            REF_TGAP(4) =     REPLAY_REF_TGAP/10000
            REF_TGAP(5) = mod(REPLAY_REF_TGAP,10000)/100
            REF_TGAP(6) = mod(REPLAY_REF_TGAP,100)
            call ESMF_TimeIntervalSet(RefTGap, &
                 YY = REF_TGAP(1), &
                 MM = REF_TGAP(2), &
                 D = REF_TGAP(3), &
                 H = REF_TGAP(4), &
                 M = REF_TGAP(5), &
                 S = REF_TGAP(6), &
                 startTime = current_time, _RC)
            RefTime = RefTime - RefTGap
         endif

         ! check if it's time to replay
         if(RefTime==current_time) then
            lring=.true.
         else
            lring=.false.
         endif

         ! In this case, increment RefTime to proper time
         if (REPLAY_REF_TGAP>0) then
            RefTime = current_time + RefTGap
         endif

      end subroutine check_replay_time_

      subroutine dump_n_splash_

         real(kind=4), pointer :: XTMP2d (:,:) =>NULL()
         real(kind=4), pointer :: XTMP3d(:,:,:)=>NULL()
         real(kind=4), pointer :: YTMP3d(:,:,:)=>NULL()
         real(r8), allocatable :: ana_thv (:,:,:)
         real(r8), allocatable :: ana_phis  (:,:)
         real(r8), allocatable :: ana_pkxy  (:,:,:)
         real(r8), allocatable :: ana_pkz   (:,:,:)
         real(r8), allocatable :: ana_dp    (:,:,:)
         real(r8), allocatable :: ana_pe    (:,:,:)
         real(r8), allocatable :: ana_qq    (:,:,:,:)
         real(r8), allocatable :: ana_pt    (:,:,:)
         real(r8), allocatable :: ana_u     (:,:,:)
         real(r8), allocatable :: ana_v     (:,:,:)
         real(r4), allocatable :: aux3d     (:,:,:)
         real(r4), allocatable :: UAtmpR4   (:,:,:)
         real(r4), allocatable :: VAtmpR4   (:,:,:)

         character(len=ESMF_MAXSTR) :: NAME
         real(r4), pointer :: ptr3dr4   (:,:,:)
         real(r8), pointer :: ptr3dr8   (:,:,:)
         integer :: iwind,rank,icnt
         integer :: iib,iie,jjb,jje,nq3d
         integer, parameter :: iapproach=2 ! handle pressure more carefully
         logical :: do_remap, remap_all_tracers

         do_remap = (cremap=="yes" .or. cremap=="YES")
         remap_all_tracers = (tremap=="yes" .or. tremap=="YES")
         nq3d=2 ! this routine only updates QV and OX
         iib = lbound(vars%pe,1)
         iie = ubound(vars%pe,1)
         jjb = lbound(vars%pe,2)
         jje = ubound(vars%pe,2)
         allocate(   ana_thv (iib:iie,jjb:jje,km  ) )
         allocate(   ana_pkxy(iib:iie,jjb:jje,km+1) )
         allocate(   ana_pkz (iib:iie,jjb:jje,km  ) )
         allocate(    ana_dp (iib:iie,jjb:jje,km  ) )
         allocate(    ana_pe (iib:iie,jjb:jje,km+1) )
         allocate(    ana_qq (iib:iie,jjb:jje,km  ,nq3d) )
         allocate(    ana_pt (iib:iie,jjb:jje,km  ) )
         allocate(     ana_u (grid%is:grid%ie  ,grid%js:grid%je+1,km) )
         allocate(     ana_v (grid%is:grid%ie+1,grid%js:grid%je  ,km) )
         ! U
         iwind=0
         if( trim(uname).ne.'NULL' ) then
            call ESMFL_BundleGetPointerToData(ana_bundle, trim(uname), XTMP3d, _RC)
            iwind=iwind+1
         endif
         ! V
         if( trim(vname).ne.'NULL' ) then
            call ESMFL_BundleGetPointerToData(ana_bundle, trim(vname), YTMP3D, _RC)
            iwind=iwind+1
         endif

         ! calculate d-grid winds
         if(iwind==0) then
            ana_u = vars%u(grid%is:grid%ie,grid%js:grid%je,1:km)
            ana_v = vars%v(grid%is:grid%ie,grid%js:grid%je,1:km)
         else if(iwind==1) then
            status=1
            call WRITE_PARALLEL('cannot handle single wind component')
            VERIFY_(STATUS)
         else if (iwind==2) then
#ifdef INC_WINDS
            if (iapproach==1) then
#endif /* INC_WINDS */
               allocate(cubeTEMP3D(grid%is:grid%ie,grid%js:grid%je,km) )
               allocate(cubeVTMP3D(grid%is:grid%ie,grid%js:grid%je,km) )
#ifdef SCALAR_WINDS
               call WRITE_PARALLEL('Replaying winds as scalars')
               call l2c%regrid(XTMP3d, cubeTEMP3D, _RC)
               call l2c%regrid(YTMP3d, cubeVTMP3D, _RC)
#else
               call WRITE_PARALLEL('Replaying winds')
               call l2c%regrid(XTMP3d, YTMP3d, cubeTEMP3d, cubeVTMP3d, rc=status)
#endif /* SCALAR_WINDS */
               allocate( UAtmp(grid%is:grid%ie  ,grid%js:grid%je  ,km) )
               allocate( VAtmp(grid%is:grid%ie  ,grid%js:grid%je  ,km) )
               UAtmp = cubetemp3d ! A-grid winds on cube
               VAtmp = cubevtmp3d ! A-grid winds on cube
               deallocate(cubeTEMP3D)
               deallocate(cubeVTMP3D)
               allocate( UDtmp(grid%is:grid%ie  ,grid%js:grid%je+1,km) )
               allocate( VDtmp(grid%is:grid%ie+1,grid%js:grid%je  ,km) )
               call Agrid_To_Native( UAtmp, VAtmp, UDtmp, VDtmp ) ! Calculate D-grid winds from rotated A-grid winds
               ana_u = UDtmp(grid%is:grid%ie,grid%js:grid%je,1:km)
               ana_v = VDtmp(grid%is:grid%ie,grid%js:grid%je,1:km)
               deallocate(udtmp,vdtmp)
               deallocate(uatmp,vatmp)
#ifdef INC_WINDS
            else ! approach 2: operate on increments
               allocate(cubeTEMP3D(grid%is:grid%ie,grid%js:grid%je,km) )
               allocate(cubeVTMP3D(grid%is:grid%ie,grid%js:grid%je,km) )
               allocate( UAtmpR4(grid%is:grid%ie  ,grid%js:grid%je  ,km) )
               allocate( VAtmpR4(grid%is:grid%ie  ,grid%js:grid%je  ,km) )
               ! get background A-grid winds
               call getAllWinds(vars%u, vars%v, UR=ana_u, VR=ana_v)
               ! transform background A-grid winds to lat-lon
               call regridder_manager%make_regridder(esmfgrid, ana_grid, REGRID_METHOD_BILINEAR, _RC)
               cubeTEMP3d = ana_u(grid%is:grid%ie,grid%js:grid%je,1:km) ! copy to satisfy interface below
               cubeVTMP3d = ana_v(grid%is:grid%ie,grid%js:grid%je,1:km) ! copy to satisfy interface below
               call c2l%regrid(cubeTEMP3d, cubeVTMP3d, UAtmpR4, VAtmpR4, _RC)
               ! calculate unrotated analysis increments of lat-lon U/V-A-grid winds
               UAtmpR4 = XTMP3d-UAtmpR4
               UAtmpR4 = VTMP3d-VAtmpR4
               ! convert the lat-lon A-grid wind increment back to the cubed
               call WRITE_PARALLEL('Replaying winds')
               call l2c%regrid(UAtmpR4, VAtmpR4, cubeTEMP3d, cubeVTMP3d, _RC)
               ! convert cubed wind increment to D-grid
               allocate( UDtmp(grid%is:grid%ie  ,grid%js:grid%je+1,km) )
               allocate( VDtmp(grid%is:grid%ie+1,grid%js:grid%je  ,km) )
               deallocate(ana_u,ana_v)
               allocate( ana_u(grid%is:grid%ie  ,grid%js:grid%je  ,km) )
               allocate( ana_v(grid%is:grid%ie  ,grid%js:grid%je  ,km) )
               ana_u = cubeTEMP3d ! need this to satisfy interface below
               ana_v = cubeVTMP3d ! need this to satisfy interface below
               call Agrid_To_Native(ana_u, ana_v, UDtmp, VDtmp) ! Calculate D-grid winds from rotated A-grid winds
               ! update winds: rotate, cubed, D-grid analyzed winds
               deallocate(ana_u,ana_v)
               allocate( ana_u(grid%is:grid%ie  ,grid%js:grid%je+1,km) )
               allocate( ana_v(grid%is:grid%ie+1,grid%js:grid%je  ,km) )
               ana_u = vars%u + UDtmp
               ana_v = vars%v + VDtmp
               ! clean up
               deallocate(VDtmp)
               deallocate(UDtmp)
               deallocate(UAtmpR4)
               deallocate(VAtmpR4)
               deallocate(cubeVTMP3D)
               deallocate(cubeTEMP3D)
            endif
#endif /* INC_WINDS */
         endif

         ! PE or PS
         if( trim(dpname).ne.'NULL' ) then
            call ESMFL_BundleGetPointerToData(ana_bundle, trim(dpname), XTMP3d, _RC)
            call WRITE_PARALLEL('Replaying '//trim(dpname))
            if ( iapproach == 1 ) then ! convert lat-lon delp to cubed and proceed
               allocate(cubeTEMP3D(size(vars%pe,1),size(vars%pe,2),km))
               call l2c%regrid(XTMP3d, cubeTEMP3D, _RC)
               ana_dp=cubeTEMP3D
               deallocate(cubeTEMP3D)
            else
               ! just because pressure is such delicate beast: convert cubed delp
               ! to lat-lon, calculate an increment in lat-lon, convert increment
               ! on delp to cubed, and create cubed version of analyzed delp
               allocate(aux3d (size(XTMP3d,1),size(XTMP3d,2),km))
               allocate(cubeTEMP3D(size(vars%pe,1),size(vars%pe,2),km))
               ! delp on the cube
               cubeTEMP3D(:,:,:) = vars%pe(:,:,2:)-vars%pe(:,:,:km)
               ! transform cubed delp
               c2l => regridder_manager%make_regridder(esmfgrid, ana_grid, REGRID_METHOD_BILINEAR, _RC)
               call c2l%regrid(cubeTEMP3D, aux3d, _RC)
               ! calculate delp increment on lat-lon and transform it to cubed
               aux3d = XTMP3d - aux3d
               call l2c%regrid(aux3d, cubeTEMP3D, _RC)
               ! delp analysis on the cube (careful since want to preserve
               ! precision in delp to the best extent possible)
               ana_dp = vars%pe(:,:,2:)-vars%pe(:,:,:km) + cubeTEMP3D
               deallocate(aux3d)
               deallocate(cubeTEMP3D)
            endif
            ana_pe(:,:,1) = grid%ak(1)
            do k=2,km+1
               ana_pe(:,:,k) = ana_pe(:,:,k-1) + ana_dp(:,:,k-1)
            enddo
            pkxy = ana_pe**kappa
            do k=1,km
               ana_pkz(:,:,k) = ( pkxy(:,:,k+1)-pkxy(:,:,k) ) &
                    / ( kappa*( log(ana_pe(:,:,k+1))-log(ana_pe(:,:,k))) )
            enddo
         else
            if( trim(psname).ne.'NULL' ) then
               call ESMFL_BundleGetPointerToData(ana_bundle, trim(psname), XTMP2D, _RC)
               call WRITE_PARALLEL('Replaying '//trim(psname))
               allocate(cubeTEMP3D(size(vars%pe,1),size(vars%pe,2),1))
               allocate(     aux3D(size(XTMP2d ,1),size(XTMP2d ,2),1))
               if ( iapproach == 1 ) then ! convert lat-lon delp to cubed and proceed
                  aux3d(:,:,1)=XTMP2D ! rank-2 interface to HorzT does not work
                  call l2c%regrid(aux3d, cubeTEMP3D, _RC)
               else
                  ! operate on increment to ps
                  ! transform cubed delp
                  cubeTEMP3D(:,:,1) = vars%pe(:,:,km+1) ! cubed ps
                  c2l => regridder_manager%make_regridder(esmfgrid, ana_grid, REGRID_METHOD_BILINEAR, _RC)
                  call c2l%regrid(cubeTEMP3D, aux3d, _RC)
                  ! increment to ps on the lat-lon
                  aux3d(:,:,1) = XTMP2D - aux3d(:,:,1)
                  ! lat-lon increment to ps converted to the cube
                  call l2c%regrid(aux3d, cubeTEMP3D, _RC)
                  ! ps update on the cube
                  cubeTEMP3d(:,:,1) = vars%pe(:,:,km+1) + cubeTEMP3D(:,:,1)
               endif
               do k=1,km+1
                  ana_pe(:,:,k) = grid%ak(k) + cubeTEMP3d(:,:,1)*grid%bk(k)
               enddo
               deallocate(aux3D)
               deallocate(cubeTEMP3D)
               do k=2,km+1
                  ana_dp(:,:,k-1) = ana_pe(:,:,k) - ana_pe(:,:,k-1)
               enddo
               pkxy = ana_pe**kappa
               do k=1,km
                  ana_pkz(:,:,k) = ( pkxy(:,:,k+1)-pkxy(:,:,k) ) &
                       / ( kappa*( log(ana_pe(:,:,k+1))-log(ana_pe(:,:,k))) )
               enddo
            else
               ana_pe  = vars%pe
               ana_pkz = vars%pkz
            endif
         endif

         ! O3
         if( trim(o3name).ne.'NULL' ) then
            call ESMFL_BundleGetPointerToData(ana_bundle, trim(o3name), XTMP3d, _RC)
            allocate(cubeTEMP3D(size(vars%pe,1),size(vars%pe,2),km))
            call l2c%regrid(XTMP3d, cubeTEMP3D, _RC)

            ! Ozone needs to be adjusted to OX
            call WRITE_PARALLEL('Replaying '//trim(o3name))

            call MAPL_Get(MAPL, LONS=LONS, LATS=LATS, ORBIT=ORBIT, _RC)

            allocate( ZTH( size(LONS,1),size(LONS,2) ) )
            allocate( SLR( size(LONS,1),size(LONS,2) ) )

            call MAPL_SunGetInsolation(LONS, LATS, ORBIT, ZTH, SLR, clock=clock, _RC)

            pl = ( vars%pe(:,:,2:) + vars%pe(:,:,:km) ) * 0.5

            do L=1,km
               if( ooo%is_r4 ) then
                  where(PL(:,:,L) >= 100.0 .or. ZTH <= 0.0) &
                       ooo%content_r4(:,:,L) = max(0.,cubeTEMP3D(:,:,L)*(MAPL_AIRMW/MAPL_O3MW)*1.0E-6)
               else
                  where(PL(:,:,L) >= 100.0 .or. ZTH <= 0.0) &
                       ooo%content   (:,:,L) = max(0.,cubeTEMP3D(:,:,L)*(MAPL_AIRMW/MAPL_O3MW)*1.0E-6)
               endif
            enddo

            deallocate( ZTH, SLR )
            deallocate(cubeTEMP3D)
         endif
         if( ooo%is_r4 ) then ! ana_qq(2) used as aux var to hold ox
            ana_qq(:,:,:,2) = ooo%content_r4
         else
            ana_qq(:,:,:,2) = ooo%content
         endif

         ! QV
         if( trim(qname).ne.'NULL' ) then
            call ESMFL_BundleGetPointerToData(ana_bundle, trim(qname), XTMP3d, _RC)
            allocate(cubeTEMP3D(size(vars%pe,1),size(vars%pe,2),km))
            call l2c%regrid(XTMP3d, cubeTEMP3D, _RC)
            call WRITE_PARALLEL('Replaying '//trim(qname))
            if( qqq%is_r4 ) then
               qqq%content_r4 = max(0.,cubeTEMP3D)
            else
               qqq%content    = max(0.,cubeTEMP3D)
            endif
            deallocate(cubeTEMP3D)
         endif
         if( qqq%is_r4 ) then ! ana_qq(1) used as aux var to calculate pt/pthv
            ana_qq(:,:,:,1) = qqq%content_r4
         else
            ana_qq(:,:,:,1) = qqq%content
         endif

         ! PT
         if( trim(tname).ne.'NULL' ) then
            call ESMFL_BundleGetPointerToData(ana_bundle, trim(tname), XTMP3d, _RC)
            allocate(cubeTEMP3D(size(ana_thv,1),size(ana_thv,2),km))
            call l2c%regrid(XTMP3d, cubeTEMP3D, _RC)
            call WRITE_PARALLEL('Replaying '//trim(tname)// '; treated as '//trim(tvar))
            if( trim(tvar).eq.'THETAV' ) ana_thv = cubeTEMP3D
            if( trim(tvar).eq.'TV'     ) ana_thv = cubeTEMP3D/ana_pkz
            if( trim(tvar).eq.'THETA' .or. &
                 trim(tvar).eq.'T'      ) then
               if( trim(tvar).eq.'THETA' ) ana_thv = cubeTEMP3D*(1.0+eps*ana_qq(:,:,:,1))
               if( trim(tvar).eq.'T'     ) ana_thv = cubeTEMP3D*(1.0+eps*ana_qq(:,:,:,1))/ana_pkz
            endif
            deallocate(cubeTEMP3D)
            ana_pt  = ana_thv/(1.0+eps*ana_qq(:,:,:,1))
         else
            ana_thv = vars%pt*(1.0+eps*ana_qq(:,:,:,1))
            ana_pt  = vars%pt
         endif

         ! Refresh vars ("update" them)
         vars%u   = ana_u(grid%is:grid%ie,grid%js:grid%je,:)
         vars%v   = ana_v(grid%is:grid%ie,grid%js:grid%je,:)
         vars%pe  = ana_pe
         vars%pkz = ana_pkz
         vars%pt  = ana_pt

         ! clean up
         deallocate( ana_v       )
         deallocate( ana_u       )
         deallocate( ana_pt      )
         deallocate( ana_qq      )
         deallocate( ana_dp      )
         deallocate( ana_pe      )
         deallocate( ana_pkz     )
         deallocate( ana_pkxy    )
         deallocate( ana_thv     )

         call WRITE_PARALLEL('Dump_n_Splash Replay Done')
      end subroutine dump_n_splash_

      subroutine incremental_

         real(r8), allocatable :: dpkxy  (:,:,:)
         real(r8), allocatable :: dpkz   (:,:,:)
         real(r8), allocatable :: dpe    (:,:,:)
         real(r8), allocatable :: dqqv   (:,:,:)
         real(r8), allocatable :: dqox   (:,:,:)
         real(r8), allocatable :: dth    (:,:,:)
         real(r8), allocatable :: du     (:,:,:)
         real(r8), allocatable :: dv     (:,:,:)
         real(r4), allocatable :: aux3d  (:,:,:)
         integer :: iib,iie,jjb,jje
         integer :: iwind
         logical :: allhere,iamr4

         iib = lbound(vars%pe,1)
         iie = ubound(vars%pe,1)
         jjb = lbound(vars%pe,2)
         jje = ubound(vars%pe,2)
         allocate( dpkxy(iib:iie,jjb:jje,km+1) )
         allocate( dpkz (iib:iie,jjb:jje,km  ) )
         allocate(  dpe (iib:iie,jjb:jje,km+1) )
         allocate( dqqv (iib:iie,jjb:jje,km  ) )
         allocate( dqox (iib:iie,jjb:jje,km  ) )
         allocate(  dth (iib:iie,jjb:jje,km  ) )
         allocate(   du (grid%is:grid%ie  ,grid%js:grid%je+1,km) )
         allocate(   dv (grid%is:grid%ie+1,grid%js:grid%je  ,km) )
         dpkxy=0.0d0
         dpkz =0.0d0
         dpe  =0.0d0
         dqqv =0.0d0
         dqox =0.0d0
         dth  =0.0d0
         du   =0.0d0
         dv   =0.0d0

         allhere = trim(uname ).ne.'NULL'.and.trim(vname ).ne.'NULL'.and. &
              trim(o3name).ne.'NULL'.and. &
              trim(tname ).ne.'NULL'.and.trim(qname ).ne.'NULL'
         if(.not.allhere) then
            call WRITE_PARALLEL('Not all varibles needed for replay are available')
            status = 999
            VERIFY_(status)
         endif
         call WRITE_PARALLEL('Starting incremental replay')

         ! U
         iwind=0
         if( trim(uname).ne.'NULL' ) then
            call ESMFL_BundleGetPointerToData(ana_bundle, trim(uname), TEMP3D, _RC)
            iwind=iwind+1
         endif
         ! V
         if( trim(vname).ne.'NULL' ) then
            call ESMFL_BundleGetPointerToData(ana_bundle, trim(vname), VTMP3D, _RC)
            iwind=iwind+1
         endif

         ! calculate d-grid winds
         if(iwind==1) then
            status=1
            print *, 'cannot handle single wind component'
            VERIFY_(STATUS)
         else if (iwind==2) then
            allocate(cubeTEMP3D(grid%is:grid%ie,grid%js:grid%je,km) )
            allocate(cubeVTMP3D(grid%is:grid%ie,grid%js:grid%je,km) )
#ifdef SCALAR_WINDS
            call WRITE_PARALLEL('Replaying increment of winds as scalars')
            call l2c%regrid(TEMP3D, cubeTEMP3D, _RC)
            call l2c%regrid(VTMP3D, cubeVTMP3D, _RC)
#else
            call WRITE_PARALLEL('Replaying increment of winds')
            call l2c%regrid(TEMP3d, VTMP3d, cubeTEMP3d, cubeVTMP3d, _RC)
#endif /* SCALAR_WINDS */
            allocate( UAtmp(grid%is:grid%ie  ,grid%js:grid%je  ,km) )
            allocate( VAtmp(grid%is:grid%ie  ,grid%js:grid%je  ,km) )
            UAtmp = cubetemp3d ! A-grid winds on cube
            VAtmp = cubevtmp3d ! A-grid winds on cube
            call Agrid_To_Native(UAtmp, VAtmp, du, dv) ! Calculate D-grid winds from rotated A-grid winds
            deallocate(uatmp,vatmp)
            deallocate(cubeTEMP3D)
            deallocate(cubeVTMP3D)
         endif

         ! DELP
         if( trim(psname)=='NULL' .and. trim(dpname).ne.'NULL' ) then
            call ESMFL_BundleGetPointerToData(ana_bundle, trim(dpname), TEMP3D, _RC)
            call WRITE_PARALLEL('Replaying increment of '//trim(dpname))
            allocate(cubeTEMP3D(size(vars%pe,1),size(vars%pe,2),km))
            call l2c%regrid(TEMP3D, cubeTEMP3D, _RC)
            dpe(:,:,1) = 0.0
            do k=2,km+1
               dpe(:,:,k) = dpe(:,:,k-1) + cubeTEMP3D(:,:,k-1)
            enddo
            deallocate(cubeTEMP3D)

            pkxy =            (vars%pe)** kappa
            dpkxy = kappa*(pkxy/vars%pe)*dpe
            do k=1,km
               dpkz(:,:,k) = (  (    dpkxy (:,:,k+1) -   dpkxy(:,:,k) )* &
                    log((vars%pe (:,:,k+1))/(vars%pe(:,:,k) )) &
                    -  (     pkxy (:,:,k+1) -    pkxy(:,:,k) )* &
                    (     dpe  (:,:,k+1) * vars%pe(:,:,k) &
                    -     dpe  (:,:,k)   * vars%pe(:,:,k+1) ) &
                    / (vars%pe(:,:,k+1)*vars%pe(:,:,k)) &
                    )  / (kappa*( log(vars%pe(:,:,k+1)/vars%pe(:,:,k)) )**2)
            enddo
         endif

         ! PS
         if( trim(psname)/='NULL' .and. trim(dpname)=='NULL' ) then
            call ESMFL_BundleGetPointerToData(ana_bundle, trim(psname), TEMP2D, _RC)
            call WRITE_PARALLEL('Replaying increment of '//trim(psname))
            allocate(cubeTEMP3D(size(vars%pe,1),size(vars%pe,2),1))
            allocate(     aux3D(size( TEMP2D,1),size( TEMP2D,2),1))
            aux3d(:,:,1) = TEMP2D ! same trick of putting in rank-3 array for transforms
            call l2c%regrid(aux3d, cubeTEMP3D, _RC)
            do k=2,km+1
               dpe(:,:,k-1) =  grid%ak(k) - grid%ak(k-1) + cubeTEMP3d(:,:,1)*(grid%bk(k)-grid%bk(k-1))
            enddo
            deallocate(     aux3d)
            deallocate(cubeTEMP3D)

            pkxy =            (vars%pe)** kappa
            dpkxy = kappa*(pkxy/vars%pe)*dpe
            do k=1,km
               dpkz(:,:,k) = (  (    dpkxy (:,:,k+1) -   dpkxy(:,:,k) )* &
                    log((vars%pe (:,:,k+1))/(vars%pe(:,:,k) )) &
                    -  (     pkxy (:,:,k+1) -    pkxy(:,:,k) )* &
                    (     dpe  (:,:,k+1) * vars%pe(:,:,k) &
                    -     dpe  (:,:,k)   * vars%pe(:,:,k+1) ) &
                    / (vars%pe(:,:,k+1)*vars%pe(:,:,k)) &
                    )  / (kappa*( log(vars%pe(:,:,k+1)/vars%pe(:,:,k)) )**2)
            enddo
         endif

         ! O3
         if( trim(o3name).ne.'NULL' ) then
            call ESMFL_BundleGetPointerToData(ana_bundle, trim(o3name), TEMP3D, _RC)
            allocate(cubeTEMP3D(size(vars%pe,1),size(vars%pe,2),km))
            call l2c%regrid(TEMP3D, cubeTEMP3D, _RC)

            ! Ozone needs to be adjusted to OX
            call WRITE_PARALLEL('Replaying increment of '//trim(o3name))

            call MAPL_Get(MAPL, LONS=LONS, LATS=LATS, ORBIT=ORBIT, _RC)

            allocate( ZTH( size(LONS,1),size(LONS,2) ) )
            allocate( SLR( size(LONS,1),size(LONS,2) ) )

            call MAPL_SunGetInsolation(LONS, LATS, ORBIT, ZTH, SLR, clock=clock, _RC)

            pl = ( vars%pe(:,:,2:) + vars%pe(:,:,:km) ) * 0.5

            do L=1,km
               where(PL(:,:,L) >= 100.0 .or. ZTH <= 0.0) &
                    dqox(:,:,L) = cubeTEMP3D(:,:,L)*(MAPL_AIRMW/MAPL_O3MW)*1.0E-6
            enddo

            deallocate( ZTH, SLR )
            deallocate(cubeTEMP3D)
         endif

         ! QV
         if( trim(qname).ne.'NULL' ) then
            call ESMFL_BundleGetPointerToData(ana_bundle, trim(qname), TEMP3D, _RC)
            allocate(cubeTEMP3D(size(vars%pe,1),size(vars%pe,2),km))
            call l2c%regrid(TEMP3D, cubeTEMP3D, _RC)
            call WRITE_PARALLEL('Replaying increment of '//trim(qname))
            dqqv = cubeTEMP3D
            deallocate(cubeTEMP3D)
         endif

         ! PT
         if( trim(tname).ne.'NULL' ) then
            if(trim(tvar).ne.'TV') then
               call WRITE_PARALLEL('Error: Cannot Replay TVAR '//trim(tvar))
               STATUS=99
               VERIFY_(STATUS)
            endif
            if(trim(tname).ne.'tv') then
               call WRITE_PARALLEL('Error: Cannot Replay TNAME '//trim(tname))
               STATUS=99
               VERIFY_(STATUS)
            endif
            call ESMFL_BundleGetPointerToData(ana_bundle, trim(tname), TEMP3D, _RC)
            allocate(cubeTEMP3D(size(vars%pe,1),size(vars%pe,2),km))
            call l2c%regrid(TEMP3D, cubeTEMP3D, _RC)
            call WRITE_PARALLEL('Replaying increment of '//trim(tname))
            ! have an incremental change to virtual temperature;
            ! want an incremental change to dry potential temperature
            ! calculate first incremental change to t-dry (save in dth for now)
            if( qqq%is_r4 ) then
               dth = (cubeTEMP3D - eps*vars%pt*vars%pkz*dqqv)/(1.0+eps*qqq%content_r4)
            else
               dth = (cubeTEMP3D - eps*vars%pt*vars%pkz*dqqv)/(1.0+eps*qqq%content   )
            endif
            ! finally calculate increment to dry theta
            dth = (dth - vars%pt*dpkz)/vars%pkz
            deallocate(cubeTEMP3D)
         endif

         ! Only at the end, apply incremental correction to pressure,
         ! potential temperature and water vapor
         vars%u   = vars%u   + sclinc * du(grid%is:grid%ie,grid%js:grid%je,1:km)
         vars%v   = vars%v   + sclinc * dv(grid%is:grid%ie,grid%js:grid%je,1:km)
         pkxy =     pkxy + sclinc * dpkxy
         vars%pkz = vars%pkz + sclinc * dpkz
         vars%pe  = vars%pe  + sclinc * dpe
         vars%pt  = vars%pt  + sclinc * dth
         if( qqq%is_r4 ) then  ! protection for negative qv is slightly inconsistent w/ update of temperature
            qqq%content_r4 = max(0.0_r4,qqq%content_r4 + sclinc*dqqv)
         else
            qqq%content    = max(0.0_r8,qqq%content    + sclinc*dqqv)
         endif
         if( ooo%is_r4 ) then  ! brute-force protection against non-zero values
            ooo%content_r4 = max(0.0_r4,ooo%content_r4 + sclinc*dqox)
         else
            ooo%content    = max(0.0_r8,ooo%content    + sclinc*dqox)
         end if

         ! clean up
         deallocate( du,dv   )
         deallocate( dth     )
         deallocate( dqox    )
         deallocate( dqqv    )
         deallocate( dpe     )
         deallocate( dpkz    )
         deallocate( dpkxy   )

         call WRITE_PARALLEL('Incremental replay complete')
      end subroutine incremental_

      subroutine state_remap_

         real(kind=4), pointer :: XTMP2d (:,:) =>NULL()
         real(kind=4), pointer :: XTMP3d(:,:,:)=>NULL()
         real(kind=4), pointer :: YTMP3d(:,:,:)=>NULL()
         real(r8), allocatable :: ana_thv (:,:,:)
         real(r8), allocatable :: ana_phis  (:,:)
         real(r8), allocatable :: ana_qq    (:,:,:,:)
         real(r8), allocatable :: ana_u     (:,:,:)
         real(r8), allocatable :: ana_v     (:,:,:)
         real(r4), allocatable :: aux3d     (:,:,:)
         !
         character(len=ESMF_MAXSTR) :: NAME
         real(r4), pointer :: ptr3dr4   (:,:,:)
         real(r8), pointer :: ptr3dr8   (:,:,:)
         integer :: iwind,icnt,nq3d,rank
         integer :: iib,iie,jjb,jje
         logical :: do_remap,remap_all_tracers

         do_remap = (cremap=="yes" .or. cremap=="YES")
         if (.not. do_remap) return

         remap_all_tracers = (tremap=="yes" .or. tremap=="YES")
         nq3d=2 ! at a minimum it will remap QV and OX
         if(do_remap.and.remap_all_tracers) then
            nq3d=0
            do N=1,NQ
               call ESMF_FieldBundleGet(BUNDLE, N, Field, _RC)
               call ESMF_FieldGet(Field, dimCount = rank, _RC)
               if (rank==2) cycle
               if (rank==3) nq3d=nq3d+1
            enddo
            write(STRING,'(A,I5,A)') "Found  ", nq3d, " 3d-tracers to remap"
            call WRITE_PARALLEL( trim(STRING)   )
         endif
         if (nq3d<2) then
            call WRITE_PARALLEL('state_remap: invalid number of tracers')
            status=999
            VERIFY_(STATUS)
         endif

         iib = lbound(vars%pe,1)
         iie = ubound(vars%pe,1)
         jjb = lbound(vars%pe,2)
         jje = ubound(vars%pe,2)

         allocate( ana_thv(iib:iie,jjb:jje,km  ) )
         allocate( ana_qq (iib:iie,jjb:jje,km  ,nq3d) )
         allocate(ana_phis(size(vars%pe,1),size(vars%pe,2)))

         if( qqq%is_r4 ) then
            ana_thv = vars%pt*(1.0+eps*qqq%content_r4(:,:,:))
         else
            ana_thv = vars%pt*(1.0+eps*qqq%content   (:,:,:))
         endif

         call WRITE_PARALLEL('Replay start remapping')
         !
         call ESMFL_BundleGetPointerToData(ana_bundle, 'phis', XTMP2D, _RC)
         allocate(cubeTEMP3D(size(vars%pe,1),size(vars%pe,2),1))
         allocate(     aux3D(size(XTMP2D ,1),size(XTMP2D ,2),1))
         aux3d(:,:,1)=XTMP2D ! this is a trick since the 2d interface to the transform has not worked for me (RT)
         call l2c%regrid(aux3D, cubeTEMP3D, _RC)
         ana_phis=cubeTEMP3D(:,:,1)
         deallocate(     aux3D)
         deallocate(cubeTEMP3D)
         !
         if (remap_all_tracers) then
            icnt=0
            do N=1,NQ
               call ESMF_FieldBundleGet(BUNDLE, N, Field, _RC)
               call ESMF_FieldGet(Field, NAME=NAME, dimCount=rank, _RC)
               if (rank==2) cycle
               if (rank==3) then
                  icnt=icnt+1
                  if (icnt>nq3d) then
                     call WRITE_PARALLEL('state_remap: number of tracers exceeds known value')
                     status=999
                     VERIFY_(STATUS)
                  endif
                  call ESMFL_BundleGetPointerToData(BUNDLE, NAME, ptr3dr4, _RC)
                  ana_qq(:,:,:,icnt) = ptr3dr4
               endif
            enddo
            if (icnt/=nq3d) then
               call WRITE_PARALLEL('state_remap: inconsitent number of tracers')
               status=999
               VERIFY_(STATUS)
            endif
         else
            if( qqq%is_r4 ) then
               ana_qq(:,:,:,1) = qqq%content_r4(:,:,:)
            else
               ana_qq(:,:,:,1) = qqq%content   (:,:,:)
            endif
            if( ooo%is_r4 ) then
               ana_qq(:,:,:,2) = ooo%content_r4(:,:,:)
            else
               ana_qq(:,:,:,2) = ooo%content   (:,:,:)
            endif
         endif ! remap_all_tracers

         call dyn_topo_remap ( vars%pe, vars%u, vars%v, ana_thv, ana_qq, ana_phis, phisxy, &
              grid%ak, grid%bk, size(ana_thv,1), size(ana_thv,2), km, nq3d )

         if (remap_all_tracers) then
            icnt=0
            do N=1,NQ
               call ESMF_FieldBundleGet(BUNDLE, N, Field, _RC)
               call ESMF_FieldGet(Field, NAME=NAME, dimCount=rank, _RC)
               if (rank==2) cycle
               if (rank==3) then
                  icnt=icnt+1
                  call ESMFL_BundleGetPointerToData(BUNDLE, NAME, ptr3dr4, _RC)
                  ptr3dr4 = ana_qq(:,:,:,icnt)
                  if(trim(NAME)=="Q") then
                     if( qqq%is_r4 ) then
                        qqq%content_r4(:,:,:) = ana_qq(:,:,:,icnt)
                     else
                        qqq%content   (:,:,:) = ana_qq(:,:,:,icnt)
                     endif
                  endif
                  if(trim(NAME)=="OX") then
                     if( ooo%is_r4 ) then
                        ooo%content_r4(:,:,:) = ana_qq(:,:,:,icnt)
                     else
                        ooo%content   (:,:,:) = ana_qq(:,:,:,icnt)
                     endif
                  endif
               endif
            enddo
         else
            if( qqq%is_r4 ) then
               qqq%content_r4(:,:,:) = ana_qq(:,:,:,1)
            else
               qqq%content   (:,:,:) = ana_qq(:,:,:,1)
            endif
            if( ooo%is_r4 ) then
               ooo%content_r4(:,:,:) = ana_qq(:,:,:,2)
            else
               ooo%content   (:,:,:) = ana_qq(:,:,:,2)
            endif
         endif ! remap_all_tracers

         if( qqq%is_r4 ) then
            vars%pt=ana_thv(:,:,:)/(1.0+eps*qqq%content_r4(:,:,:))
         else
            vars%pt=ana_thv(:,:,:)/(1.0+eps*qqq%content   (:,:,:))
         endif

         pkxy = vars%pe**kappa
         do k=1,km
            vars%pkz(:,:,k) = ( pkxy(:,:,k+1)-pkxy(:,:,k) ) &
                 / ( kappa*( log(vars%pe(:,:,k+1))-log(vars%pe(:,:,k)) ) )
         enddo

         call WRITE_PARALLEL('Replay done remapping')

         deallocate(ana_qq)
         deallocate(ana_thv)
         deallocate(ana_phis)
      end subroutine state_remap_

   end subroutine Run

   subroutine PULL_Q(STATE, import, QQQ, iNXQ, InFieldName, RC)

      type (DynState)        :: STATE
      type (ESMF_State)              :: import
      type (DynTracers)               :: QQQ       ! Specific Humidity
      integer,           intent(IN)  :: iNXQ
      character(len=*), optional, intent(IN) :: InFieldName
      integer, optional, intent(OUT) :: RC

      integer                          :: STATUS
      character(len=ESMF_MAXSTR)       :: IAm="Pull_Q"
      character(len=ESMF_MAXSTR)       :: FIELDNAME, QFieldName
      type (ESMF_FieldBundle)          :: BUNDLE
      type (ESMF_Field)                :: field
      type (ESMF_Array)                :: array
      type (ESMF_TypeKind_Flag)        :: kind
      real(r4),              pointer   :: ptr_r4(:,:,:)
      real(r8),              pointer   :: ptr_r8(:,:,:)
      integer                          :: N,NQ
      integer                          :: i1,in,j1,jn,im,jm,km


      QFieldName = "Q"
      if (present(InFieldName)) QFieldName=InFieldName

      i1 = state%grid%is
      in = state%grid%ie
      j1 = state%grid%js
      jn = state%grid%je
      im = state%grid%npx
      jm = state%grid%npy
      km = state%grid%npz

      BUNDLE = bundleAdv

      ! Count the friendlies
      call ESMF_FieldBundleGet(BUNDLE, fieldCount=NQ, RC=STATUS)
      VERIFY_(STATUS)

      NQ = NQ + iNXQ
      STATE%GRID%NQ = NQ       ! GRID%NQ is now the "official" NQ

      ! Tracer pointer array
      IF( ASSOCIATED( STATE%VARS%tracer ) ) then
         call freeTracers(state)
      ENDIF

      ALLOCATE(STATE%VARS%tracer(nq), STAT=STATUS)
      VERIFY_(STATUS)

      DO n = 1, NQ-iNXQ
         call ESMF_FieldBundleGet(bundle, fieldIndex=n, field=field, rc=status)
         VERIFY_(STATUS)
         call ESMF_FieldGet(FIELD, Array=Array, name=fieldname, RC=STATUS)
         VERIFY_(STATUS)
         call ESMF_ArrayGet(array,typekind=kind,rc=status)
         VERIFY_(STATUS)

         STATE%VARS%TRACER(N)%IS_R4  = (kind == ESMF_TYPEKIND_R4)   ! Is real*4?

         STATE%VARS%TRACER(N)%TNAME = fieldname

         if ( STATE%VARS%TRACER(N)%IS_R4 ) then
            call ESMF_ArrayGet(array, localDE=0, farrayptr=ptr_r4, rc=status)
            VERIFY_(STATUS)
            state%vars%tracer(n)%content_r4 => MAPL_RemapBounds(PTR_R4, i1,in,j1,jn, &
                 1, km)
            if (fieldname == QFieldName) then
               qqq%is_r4 = .true.
               qqq%content_r4 => state%vars%tracer(n)%content_r4
            end if

         else

            call ESMF_ArrayGet(array, localDE=0, farrayptr=ptr_r8, rc=status)
            VERIFY_(STATUS)

            state%vars%tracer(n)%content => PTR_R8
            if (fieldname == QFieldName) then
               qqq%is_r4   = .false.
               qqq%content => state%vars%tracer(n)%content
            end if

         endif
      END DO

   end subroutine PULL_Q

   !BOP
   !IROUTINE: RunAddIncs

   !DESCRIPTION: This is the second registered stage of FV.
   !    It calls an Fv supplied routine to add external contributions
   !    to FV's state variables. It does not touch the Friendly tracers.
   !    It also computes additional diagnostics and updates the
   !    FV internal state to reflect the added tendencies.

   !INTERFACE:
   subroutine RunAddIncs(gc, import, export, clock, rc)

      !ARGUMENTS:
      type(ESMF_GridComp) :: gc
      type(ESMF_State) :: import
      type(ESMF_State) :: export
      type(ESMF_Clock) :: clock
      integer, intent(out) :: rc
      !EOP

      integer                                          :: status
      character(len=ESMF_MAXSTR) :: IAm

      type (MAPL_MetaComp), pointer :: genstate

      type (DYN_wrap) :: wrap
      type (DynState), pointer :: STATE
      type (DynGrid),  pointer :: GRID
      type (DynVars),  pointer :: VARS
      type (DynTracers)                 :: qqq     ! Specific Humidity

      real(r8), allocatable :: penrg (:,:)   ! Vertically Integrated Cp*T
      real(r8), allocatable :: kenrg (:,:)   ! Vertically Integrated K
      real(r8), allocatable :: tenrg (:,:)   ! PHIS*(Psurf-Ptop)
      real(r8), allocatable :: penrg0(:,:)   ! Vertically Integrated Cp*T
      real(r8), allocatable :: kenrg0(:,:)   ! Vertically Integrated K
      real(r8), allocatable :: tenrg0(:,:)   ! PHIS*(Psurf-Ptop)

      real(r8),     pointer :: phisxy(:,:)
      real(r4),     pointer ::   phis(:,:)
      real(r8), allocatable ::    slp(:,:)
      real(r8), allocatable ::  H1000(:,:)
      real(r8), allocatable ::  H850 (:,:)
      real(r8), allocatable ::  H500 (:,:)
      real(r8), allocatable ::  tmp3d(:,:,:)
      real(r8), allocatable ::    plk(:,:,:)
      real(r8), allocatable ::    pke(:,:,:)
      real(r8), allocatable ::     pl(:,:,:)
      real(r8), allocatable ::     ua(:,:,:)
      real(r8), allocatable ::     va(:,:,:)
      real(r8), allocatable ::     uc(:,:,:)
      real(r8), allocatable ::     vc(:,:,:)
      real(r8), allocatable ::     ur(:,:,:)
      real(r8), allocatable ::     vr(:,:,:)
      real(r8), allocatable ::     qv(:,:,:)
      real(r8), allocatable ::     dp(:,:,:)
      real(r8), allocatable ::    thv(:,:,:)
      real(r8), allocatable ::    zle(:,:,:)
      real(r8), allocatable :: tempxy(:,:,:)
      real(r8)              :: TMAX, TMIN

      real(r8), allocatable ::  logpl(:,:,:)
      real(r8), allocatable ::  logpe(:,:,:)
      real(r8), allocatable ::  logps(:,:)

      real(FVPRC)              :: dt

      real(r4), pointer     :: QOLD(:,:,:)
      real(r4), pointer     :: temp3d(:,:,:)
      real(r4), pointer     :: temp2d(:,:  )

      integer ifirstxy, ilastxy
      integer jfirstxy, jlastxy
      integer im,jm,km, iNXQ
      real(r4), pointer     :: ztemp1(:,:  )
      real(r4), pointer     :: ztemp2(:,:  )
      real(r4), pointer     :: ztemp3(:,:  )

      real(kind=4), allocatable :: dthdtphyint1(:,:)
      real(kind=4), allocatable :: dthdtphyint2(:,:)

      logical :: doEnergetics

      integer i,j,k

      character(len=ESMF_MAXSTR) :: comp_name

      Iam = "RunAddIncs"
      call ESMF_GridCompGet( gc, name=comp_name, RC=STATUS )
      VERIFY_(STATUS)
      Iam = trim(comp_name) // trim(Iam)

      ! Retrieve the pointer to the generic state
      call MAPL_GetObjectFromGC (gc, GENSTATE,  RC=STATUS )
      VERIFY_(STATUS)

      call MAPL_TimerOn(GENSTATE,"TOTAL")
      call MAPL_TimerOn(GENSTATE,"RUN2")

      ! Retrieve the pointer to the internal state
      call ESMF_UserCompGetInternalState(gc, 'DYNstate', wrap, status)
      VERIFY_(STATUS)
      state => wrap%dyn_state

      vars  => state%vars   ! direct handle to control variables
      grid  => state%grid   ! direct handle to grid
      dt    =  state%dt     ! dynamics time step (large)

      ifirstxy = grid%is
      ilastxy  = grid%ie
      jfirstxy = grid%js
      jlastxy  = grid%je

      im  = grid%npx
      jm  = grid%npy
      km  = grid%npz
      iNXQ = 0

      if (.not. SW_DYNAMICS) then

         ALLOCATE( dthdtphyint1(ifirstxy:ilastxy,jfirstxy:jlastxy) )
         ALLOCATE( dthdtphyint2(ifirstxy:ilastxy,jfirstxy:jlastxy) )

         doEnergetics=.false.
         call MAPL_GetPointer(export,temp2D,'KE'   ,rc=status); VERIFY_(STATUS)
         if(associated(temp2D)) doEnergetics=.true.
         call MAPL_GetPointer(export,temp2D,'KEPHY',rc=status); VERIFY_(STATUS)
         if(associated(temp2D)) doEnergetics=.true.
         call MAPL_GetPointer(export,temp2D,'PEPHY',rc=status); VERIFY_(STATUS)
         if(associated(temp2D)) doEnergetics=.true.
         call MAPL_GetPointer(export,temp2D,'TEPHY',rc=status); VERIFY_(STATUS)
         if(associated(temp2D)) doEnergetics=.true.
         if (doEnergetics) then
            ALLOCATE(  kenrg(ifirstxy:ilastxy,jfirstxy:jlastxy) )
            ALLOCATE(  penrg(ifirstxy:ilastxy,jfirstxy:jlastxy) )
            ALLOCATE(  tenrg(ifirstxy:ilastxy,jfirstxy:jlastxy) )
            ALLOCATE( kenrg0(ifirstxy:ilastxy,jfirstxy:jlastxy) )
            ALLOCATE( penrg0(ifirstxy:ilastxy,jfirstxy:jlastxy) )
            ALLOCATE( tenrg0(ifirstxy:ilastxy,jfirstxy:jlastxy) )
         endif

         ALLOCATE(  tmp3d(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
         ALLOCATE( phisxy(ifirstxy:ilastxy,jfirstxy:jlastxy) )
         ALLOCATE(  logps(ifirstxy:ilastxy,jfirstxy:jlastxy) )

         ALLOCATE(     ua(ifirstxy:ilastxy,jfirstxy:jlastxy,km)   )
         ALLOCATE(     va(ifirstxy:ilastxy,jfirstxy:jlastxy,km)   )
         ALLOCATE(     uc(ifirstxy:ilastxy,jfirstxy:jlastxy,km)   )
         ALLOCATE(     vc(ifirstxy:ilastxy,jfirstxy:jlastxy,km)   )
         ALLOCATE(     ur(ifirstxy:ilastxy,jfirstxy:jlastxy,km)   )
         ALLOCATE(     vr(ifirstxy:ilastxy,jfirstxy:jlastxy,km)   )
         ALLOCATE(     qv(ifirstxy:ilastxy,jfirstxy:jlastxy,km)   )
         ALLOCATE(     pl(ifirstxy:ilastxy,jfirstxy:jlastxy,km)   )
         ALLOCATE(  logpl(ifirstxy:ilastxy,jfirstxy:jlastxy,km)   )
         ALLOCATE(     dp(ifirstxy:ilastxy,jfirstxy:jlastxy,km)   )
         ALLOCATE(    thv(ifirstxy:ilastxy,jfirstxy:jlastxy,km)   )
         ALLOCATE( tempxy(ifirstxy:ilastxy,jfirstxy:jlastxy,km)   )

         ALLOCATE(    plk(ifirstxy:ilastxy,jfirstxy:jlastxy,km)   )
         ALLOCATE(    pke(ifirstxy:ilastxy,jfirstxy:jlastxy,km+1) )
         ALLOCATE(  logpe(ifirstxy:ilastxy,jfirstxy:jlastxy,km+1) )
         ALLOCATE(    zle(ifirstxy:ilastxy,jfirstxy:jlastxy,km+1) )


         call MAPL_GetPointer ( import, PHIS, 'PHIS', RC=STATUS )
         VERIFY_(STATUS)

         phisxy = real(phis,kind=r8)

         ! Compute Pressure Thickness
         dp = ( vars%pe(:,:,2:) - vars%pe (:,:,:km) )

         ! Load Specific Humidity
         call MAPL_GetPointer(export,QOLD,'Q',  rc=status)

         call PULL_Q ( STATE, import, qqq, iNXQ, RC=rc )
         if ((.not. ADIABATIC) .and. (STATE%GRID%NQ > 0)) then
            if ( (qqq%is_r4) .and. (associated(qqq%content_r4)) ) then
               if (size(qv)==size(qqq%content_r4)) qv = qqq%content_r4
            elseif (associated(qqq%content)) then
               if (size(qv)==size(qqq%content)) qv = qqq%content
            endif
            _ASSERT(all(qv >= 0.0),'RunAddIncs: negative or nan water vapor detected')
         else
            qv = 0.0
         endif

         ! Compute Energetics Before Diabatic Forcing
         if (associated(QOLD)) then
            thv = vars%pt*(1.0+eps*QOLD)
         else
            thv = vars%pt
         endif

         if (doEnergetics) then
            call getAllWinds(vars%u, vars%v, UA=ua, VA=va, UC=uc, VC=vc, UR=ur, VR=vr)
            call Energetics (ur,vr,thv,vars%pe,dp,vars%pkz,phisxy,kenrg0,penrg0,tenrg0)
         endif

         ! DTHVDTPHYINT
         call MAPL_GetPointer ( export, temp2D, 'DTHVDTPHYINT', rc=status )
         VERIFY_(STATUS)
         if( associated(temp2D) ) then
            dthdtphyint1 = 0.0
            do k=1,km
               dthdtphyint1 = dthdtphyint1 + thv(:,:,k)*dp(:,:,k)
            enddo
         endif

         ! Add Diabatic Forcing to State Variables
         call MAPL_TimerOn (GENSTATE,"PHYS_ADD_INCS")
         call ADD_INCS ( GENSTATE,STATE,import,DT )
         call MAPL_TimerOff(GENSTATE,"PHYS_ADD_INCS")


         if (DYN_DEBUG) call DEBUG_FV_STATE('PHYSICS ADD_INCS',STATE)

         ! Update Mid-Layer Pressure and Pressure Thickness
         dp = ( vars%pe(:,:,2:) - vars%pe (:,:,:km) )
         pl = ( vars%pe(:,:,2:) + vars%pe (:,:,:km) )*0.5

         logpl = log(pl)
         logpe = log(vars%pe)
         logps = log(vars%pe(:,:,km+1))

         ! Get Cubed-Sphere Wind Exports
         call getAllWinds(vars%u, vars%v, UA=ua, VA=va, UC=uc, VC=vc, UR=ur, VR=vr)
         call FILLOUT3 (export, 'U_DGRID', vars%u  , rc=status); VERIFY_(STATUS)
         call FILLOUT3 (export, 'V_DGRID', vars%v  , rc=status); VERIFY_(STATUS)
         call FILLOUT3 (export, 'U_CGRID', uc      , rc=status); VERIFY_(STATUS)
         call FILLOUT3 (export, 'V_CGRID', vc      , rc=status); VERIFY_(STATUS)
         call FILLOUT3 (export, 'U_AGRID', ua      , rc=status); VERIFY_(STATUS)
         call FILLOUT3 (export, 'V_AGRID', va      , rc=status); VERIFY_(STATUS)

         ! Compute Energetics After Diabatic Forcing
         thv = vars%pt*(1.0+eps*qv)

#if defined(DEBUG_VPT)
         call Write_Profile(grid, thv, 'VPT')
#endif

         if (doEnergetics) then
            call Energetics (ur,vr,thv,vars%pe,dp,vars%pkz,phisxy,kenrg,penrg,tenrg)
            call MAPL_GetPointer(export,temp2d,'KE',  rc=status)
            VERIFY_(STATUS)
            if(associated(temp2d)) temp2d = kenrg
            kenrg = (kenrg-kenrg0)/DT
            penrg = (penrg-penrg0)/DT
            tenrg = (tenrg-tenrg0)/DT
            call FILLOUT2 (export, 'KEPHY', kenrg, rc=status); VERIFY_(STATUS)
            call FILLOUT2 (export, 'PEPHY', penrg, rc=status); VERIFY_(STATUS)
            call FILLOUT2 (export, 'TEPHY', tenrg, rc=status); VERIFY_(STATUS)
         endif

         ! DTHVDTPHYINT
         call MAPL_GetPointer ( export, temp2D, 'DTHVDTPHYINT', rc=status )
         VERIFY_(STATUS)
         if( associated(temp2D) ) then
            dthdtphyint2 = 0.0
            do k=1,km
               dthdtphyint2 = dthdtphyint2 + thv(:,:,k)*dp(:,:,k)
            enddo
            temp2D       = (dthdtphyint2-dthdtphyint1) * MAPL_P00**MAPL_KAPPA / (MAPL_GRAV*DT)
         endif

         plk = exp( kappa * log( 0.5*(vars%pe(:,:,1:km)+vars%pe(:,:,2:km+1)) ) )
         pke = exp( kappa * log( vars%pe ) )

         tempxy = vars%pt * vars%pkz   ! Dry Temperature

         !#if defined(DEBUG_T)
         !  call Write_Profile(grid, tempxy, 'T')
         !#endif

         if (DEBUG_DYN) then
            call MAPL_MaxMin('DYN: Q_AF_INC ', qv)
            call MAPL_MaxMin('DYN: T_AF_INC ', tempxy, pmax=TMAX, pmin=TMIN)
            call MAPL_MaxMin('DYN: U_AF_INC ', ua)
            call MAPL_MaxMin('DYN: V_AF_INC ', va)
            if (TMIN <= 130.0_r8) call Write_Profile(grid, tempxy, 'TAFINC')
            if (TMAX >= 333.0_r8) call Write_Profile(grid, tempxy, 'TAFINC')
         endif

         call FILLOUT3 (export, 'DELP'   , dp      , rc=status); VERIFY_(STATUS)
         call FILLOUT3 (export, 'U'      , ur      , rc=status); VERIFY_(STATUS)
         call FILLOUT3 (export, 'V'      , vr      , rc=status); VERIFY_(STATUS)
         call FILLOUT3 (export, 'T'      , tempxy  , rc=status); VERIFY_(STATUS)
         call FILLOUT3 (export, 'Q'      , qv      , rc=status); VERIFY_(STATUS)
         call FILLOUT3 (export, 'PL'     , pl      , rc=status); VERIFY_(STATUS)
         call FILLOUT3 (export, 'PLE'    , vars%pe , rc=status); VERIFY_(STATUS)
         call FILLOUT3 (export, 'PLK'    , plk     , rc=status); VERIFY_(STATUS)
         call FILLOUT3 (export, 'PKE'    , pke     , rc=status); VERIFY_(STATUS)
         call FILLOUT3 (export, 'THV'    , thv     , rc=status); VERIFY_(STATUS)
         call FILLOUT3 (export, 'PT'     , vars%pt , rc=status); VERIFY_(STATUS)
         call FILLOUT3 (export, 'PE'     , vars%pe , rc=status); VERIFY_(STATUS)

         call MAPL_GetPointer(export,temp3d,'TH',rc=status)
         VERIFY_(STATUS)
         if(associated(temp3d)) temp3d = (tempxy)*(p00/(0.5*(vars%pe(:,:,1:km)+vars%pe(:,:,2:km+1))))**kappa

#ifdef SKIP_TRACERS
         do itracer=1,ntracers
            write(myTracer, "('Q',i5.5)") itracer-1
            call MAPL_GetPointer(export, temp3D, TRIM(myTracer), rc=status)
            VERIFY_(STATUS)
            if((associated(temp3d)) .and. (STATE%GRID%NQ>=itracer)) then
               if (state%vars%tracer(itracer)%is_r4) then
                  temp3d = state%vars%tracer(itracer)%content_r4
               else
                  temp3d = state%vars%tracer(itracer)%content
               endif
            endif
         enddo
#endif

         ! Compute Edge Heights
         zle(:,:,km+1) = phisxy(:,:)
         do k=km,1,-1
            zle(:,:,k) = zle(:,:,k+1) + cp*thv(:,:,k)*( pke(:,:,k+1)-pke(:,:,k) )
         enddo
         zle(:,:,:) = zle(:,:,:)/grav

         call FILLOUT3 (export, 'ZLE', zle, rc=status); VERIFY_(STATUS)

         ! Compute Mid-Layer Heights
         call MAPL_GetPointer(export,temp3d,'ZL',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp3d)) temp3d = 0.5*( zle(:,:,2:) + zle(:,:,:km) )


         ! Fill Single Level Variables
         call MAPL_GetPointer(export,temp2d,'Z700',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) then
            call VertInterp(temp2d,zle*grav,logpe,log(70000.)  , rc=status)
            VERIFY_(STATUS)
         end if

         call MAPL_GetPointer(export,temp2d,'Z500',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) then
            call VertInterp(temp2d,zle*grav,logpe,log(50000.)  , rc=status)
            VERIFY_(STATUS)
         end if

         call MAPL_GetPointer(export,temp2d,'Z300',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) then
            call VertInterp(temp2d,zle*grav,logpe,log(30000.)  , rc=status)
            VERIFY_(STATUS)
         end if

         call MAPL_GetPointer(export,temp2d,'H100',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) then
            call VertInterp(temp2d,zle,logpe,log(10000.)  , rc=status)
            VERIFY_(STATUS)
         end if

         call MAPL_GetPointer(export,temp2d,'H200',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) then
            call VertInterp(temp2d,zle,logpe,log(20000.)  , rc=status)
            VERIFY_(STATUS)
         end if

         call MAPL_GetPointer(export,temp2d,'H250',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) then
            call VertInterp(temp2d,zle,logpe,log(25000.)  , rc=status)
            VERIFY_(STATUS)
         end if

         call MAPL_GetPointer(export,temp2d,'H300',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) then
            call VertInterp(temp2d,zle,logpe,log(30000.)  , rc=status)
            VERIFY_(STATUS)
         end if

         call MAPL_GetPointer(export,temp2d,'H500',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) then
            call VertInterp(temp2d,zle,logpe,log(50000.)  , rc=status)
            VERIFY_(STATUS)
         end if

         call MAPL_GetPointer(export,temp2d,'H700',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) then
            call VertInterp(temp2d,zle,logpe,log(70000.)  , rc=status)
            VERIFY_(STATUS)
         end if

         call MAPL_GetPointer(export,temp2d,'H850',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) then
            call VertInterp(temp2d,zle,logpe,log(85000.)  , rc=status)
            VERIFY_(STATUS)
         end if

         call MAPL_GetPointer(export,temp2d,'H1000',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) then
            call VertInterp(temp2d,zle,logpe,log(100000.)  , rc=status)
            VERIFY_(STATUS)
         end if

         call MAPL_GetPointer(export,temp2d,'U50M',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) then
            call VertInterp(temp2d,ur,-zle,-50., rc=status)
            VERIFY_(STATUS)
         end if

         call MAPL_GetPointer(export,temp2d,'V50M',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) then
            call VertInterp(temp2d,vr,-zle,-50., rc=status)
            VERIFY_(STATUS)
         end if

         call MAPL_GetPointer(export,temp2d,'U100',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) then
            call VertInterp(temp2d,ur,logpe,log(10000.)  ,  rc=status)
            VERIFY_(STATUS)
         end if

         call MAPL_GetPointer(export,temp2d,'U200',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) then
            call VertInterp(temp2d,ur,logpe,log(20000.)  ,  rc=status)
            VERIFY_(STATUS)
         end if

         call MAPL_GetPointer(export,temp2d,'U250',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) then
            call VertInterp(temp2d,ur,logpe,log(25000.)  ,  rc=status)
            VERIFY_(STATUS)
         end if

         call MAPL_GetPointer(export,temp2d,'U300',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) then
            call VertInterp(temp2d,ur,logpe,log(30000.)  ,  rc=status)
            VERIFY_(STATUS)
         end if

         call MAPL_GetPointer(export,temp2d,'U500',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) then
            call VertInterp(temp2d,ur,logpe,log(50000.)  ,  rc=status)
            VERIFY_(STATUS)
         end if

         call MAPL_GetPointer(export,temp2d,'U700',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) then
            call VertInterp(temp2d,ur,logpe,log(70000.)  ,  rc=status)
            VERIFY_(STATUS)
         end if

         call MAPL_GetPointer(export,temp2d,'U850',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) then
            call VertInterp(temp2d,ur,logpe,log(85000.)  ,  rc=status)
            VERIFY_(STATUS)
         end if

         call MAPL_GetPointer(export,temp2d,'V100',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) then
            call VertInterp(temp2d,vr,logpe,log(10000.)  ,  rc=status)
            VERIFY_(STATUS)
         end if

         call MAPL_GetPointer(export,temp2d,'V200',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) then
            call VertInterp(temp2d,vr,logpe,log(20000.)  ,  rc=status)
            VERIFY_(STATUS)
         end if

         call MAPL_GetPointer(export,temp2d,'V250',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) then
            call VertInterp(temp2d,vr,logpe,log(25000.)  ,  rc=status)
            VERIFY_(STATUS)
         end if

         call MAPL_GetPointer(export,temp2d,'V300',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) then
            call VertInterp(temp2d,vr,logpe,log(30000.)  ,  rc=status)
            VERIFY_(STATUS)
         end if

         call MAPL_GetPointer(export,temp2d,'V500',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) then
            call VertInterp(temp2d,vr,logpe,log(50000.)  ,  rc=status)
            VERIFY_(STATUS)
         end if

         call MAPL_GetPointer(export,temp2d,'V700',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) then
            call VertInterp(temp2d,vr,logpe,log(70000.)  ,  rc=status)
            VERIFY_(STATUS)
         end if

         call MAPL_GetPointer(export,temp2d,'V850',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) then
            call VertInterp(temp2d,vr,logpe,log(85000.)  ,  rc=status)
            VERIFY_(STATUS)
         end if

         call MAPL_GetPointer(export,temp2d,'T100',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) then
            call VertInterp(temp2d,tempxy,logpe,log(10000.)  , rc=status)
            VERIFY_(STATUS)
         end if

         call MAPL_GetPointer(export,temp2d,'T200',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) then
            call VertInterp(temp2d,tempxy,logpe,log(20000.)  , rc=status)
            VERIFY_(STATUS)
         end if

         call MAPL_GetPointer(export,temp2d,'T250',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) then
            call VertInterp(temp2d,tempxy,logpe,log(25000.)  , rc=status)
            VERIFY_(STATUS)
         end if

         call MAPL_GetPointer(export,temp2d,'T300',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) then
            call VertInterp(temp2d,tempxy,logpe,log(30000.)  , rc=status)
            VERIFY_(STATUS)
         end if

         call MAPL_GetPointer(export,temp2d,'T500',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) then
            call VertInterp(temp2d,tempxy,logpe,log(50000.)  ,  rc=status)
            VERIFY_(STATUS)
         end if

         call MAPL_GetPointer(export,temp2d,'T700',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) then
            call VertInterp(temp2d,tempxy,logpe,log(70000.)  ,  rc=status)
            VERIFY_(STATUS)
         end if

         call MAPL_GetPointer(export,temp2d,'T850',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) then
            call VertInterp(temp2d,tempxy,logpe,log(85000.)  ,  rc=status)
            VERIFY_(STATUS)
         end if

         call MAPL_GetPointer(export,temp2d,'Q100',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) then
            call VertInterp(temp2d,qv,logpe,log(10000.)  ,  positive_definite=.true., rc=status)
            VERIFY_(STATUS)
         end if

         call MAPL_GetPointer(export,temp2d,'Q200',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) then
            call VertInterp(temp2d,qv,logpe,log(20000.)  ,  positive_definite=.true., rc=status)
            VERIFY_(STATUS)
         end if

         call MAPL_GetPointer(export,temp2d,'Q250',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) then
            call VertInterp(temp2d,qv,logpe,log(25000.)  ,  positive_definite=.true., rc=status)
            VERIFY_(STATUS)
         end if

         call MAPL_GetPointer(export,temp2d,'Q300',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) then
            call VertInterp(temp2d,qv,logpe,log(30000.)  ,  positive_definite=.true., rc=status)
            VERIFY_(STATUS)
         end if

         call MAPL_GetPointer(export,temp2d,'Q500',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) then
            call VertInterp(temp2d,qv,logpe,log(50000.)  ,  positive_definite=.true., rc=status)
            VERIFY_(STATUS)
         end if

         call MAPL_GetPointer(export,temp2d,'Q700',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) then
            call VertInterp(temp2d,qv,logpe,log(70000.)  ,  positive_definite=.true., rc=status)
            VERIFY_(STATUS)
         end if

         call MAPL_GetPointer(export,temp2d,'Q850',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) then
            call VertInterp(temp2d,qv,logpe,log(85000.)  ,  positive_definite=.true., rc=status)
            VERIFY_(STATUS)
         end if

         ! Fill Model Top Level Variables
         call MAPL_GetPointer(export,temp2d,'UTOP', rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) temp2d = ur(:,:,1)

         call MAPL_GetPointer(export,temp2d,'VTOP', rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) temp2d = vr(:,:,1)

         call MAPL_GetPointer(export,temp2d,'TTOP', rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) temp2d = tempxy(:,:,1)

         call MAPL_GetPointer(export,temp2d,'DELPTOP', rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) temp2d = dp(:,:,1)

         ! Compute Surface Pressure
         call MAPL_GetPointer(export,temp2d,'PS',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) temp2d=vars%pe(:,:,km+1)

         ! Get the height above the surface
         do k=1,km+1
            zle(:,:,k) = zle(:,:,k) - zle(:,:,km+1)
         enddo

         call MAPL_GetPointer(export,temp3d,'ZLE0',rc=status)
         VERIFY_(STATUS)
         if(associated(temp3d)) temp3d = zle

         call MAPL_GetPointer(export,temp3d,'ZL0' ,rc=status)
         VERIFY_(STATUS)
         if(associated(temp3d)) temp3d = 0.5*( zle(:,:,:km)+zle(:,:,2:) )

         ! Compute Vertically Averaged T,U
         call MAPL_GetPointer(export,temp2d,'TAVE',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) then
            temp2d = 0.0
            do k=1,km
               temp2d = temp2d + tempxy(:,:,k)*dp(:,:,k)
            enddo
            temp2d = temp2d / (vars%pe(:,:,km+1)-vars%pe(:,:,1))
         endif

         call MAPL_GetPointer(export,temp2d,'UAVE',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp2d)) then
            temp2d = 0.0
            do k=1,km
               temp2d = temp2d + ur(:,:,k)*dp(:,:,k)
            enddo
            temp2d = temp2d / (vars%pe(:,:,km+1)-vars%pe(:,:,1))
         endif

         ! Convert T to Tv
         tempxy = tempxy*(1.0+eps*qv)

         call MAPL_GetPointer(export,temp3d,'TV',  rc=status)
         VERIFY_(STATUS)
         if(associated(temp3d)) temp3d=tempxy

         ! Compute Sea-Level Pressure
         call MAPL_GetPointer(export,temp2d,'SLP'  ,rc=status)
         VERIFY_(STATUS)
         call MAPL_GetPointer(export,Ztemp1,'H1000',rc=status)
         VERIFY_(STATUS)
         call MAPL_GetPointer(export,Ztemp2,'H850' ,rc=status)
         VERIFY_(STATUS)
         call MAPL_GetPointer(export,Ztemp3,'H500' ,rc=status)
         VERIFY_(STATUS)

         if(associated(temp2d) .or. associated(ztemp1) &
              .or. associated(ztemp2) &
              .or. associated(ztemp3) ) then
            ALLOCATE(  slp(ifirstxy:ilastxy,jfirstxy:jlastxy) )
            ALLOCATE(H1000(ifirstxy:ilastxy,jfirstxy:jlastxy) )
            ALLOCATE(H850 (ifirstxy:ilastxy,jfirstxy:jlastxy) )
            ALLOCATE(H500 (ifirstxy:ilastxy,jfirstxy:jlastxy) )
            do j=jfirstxy,jlastxy
               do i=ifirstxy,ilastxy
                  call get_slp ( km,vars%pe (i,j,  km+1),phisxy(i,j),  slp(i,j), &
                       vars%pe (i,j,1:km+1),                        &
                       vars%pkz(i,j,1:km  ),                        &
                       tempxy(i,j,1:km  ),                        &
                       H1000(i,j), H850(i,j), H500(i,j)          )
               enddo
            enddo

            !#define DEBUG_SLP
#if defined(DEBUG_SLP)
            call Write_Profile(grid, slp/100.0, 'SLP')
#endif

            if(associated(temp2d)) temp2d = slp
            if(associated(ztemp1)) where( ztemp1.eq.MAPL_UNDEF ) ztemp1 = H1000
            if(associated(ztemp2)) where( ztemp2.eq.MAPL_UNDEF ) ztemp2 = H850
            if(associated(ztemp3)) where( ztemp3.eq.MAPL_UNDEF ) ztemp3 = H500
            DEALLOCATE(slp,H1000,H850,H500)
         end if

         ! Deallocate Memory
         if (doEnergetics) then
            DEALLOCATE(  kenrg )
            DEALLOCATE(  penrg )
            DEALLOCATE(  tenrg )
            DEALLOCATE( kenrg0 )
            DEALLOCATE( penrg0 )
            DEALLOCATE( tenrg0 )
         endif

         DEALLOCATE(  tmp3d )

         DEALLOCATE( phisxy )

         DEALLOCATE(     ua )
         DEALLOCATE(     va )
         DEALLOCATE(     uc )
         DEALLOCATE(     vc )
         DEALLOCATE(     ur )
         DEALLOCATE(     vr )
         DEALLOCATE(     qv )
         DEALLOCATE(     pl )
         DEALLOCATE(     dp )
         DEALLOCATE( tempxy )

         DEALLOCATE(    thv )
         DEALLOCATE(    plk )
         DEALLOCATE(    pke )
         DEALLOCATE(  logpl )
         DEALLOCATE(  logpe )
         DEALLOCATE(  logps )
         DEALLOCATE(    zle )
         DEALLOCATE( dthdtphyint1 )
         DEALLOCATE( dthdtphyint2 )

         call freeTracers(state)

      end if ! .not. SW_DYNAMICS

      call MAPL_TimerOff(GENSTATE,"RUN2")
      call MAPL_TimerOff(GENSTATE,"TOTAL")

      _RETURN(_SUCCESS)
   end subroutine RunAddIncs

   subroutine ADD_INCS ( MAPL,STATE,import,DT,IS_WEIGHTED,RC )

      use fms_mod, only: set_domain, nullify_domain
      use fv_diagnostics_mod, only: prt_maxmin
      use time_manager_mod,   only: time_type
      use fv_update_phys_mod, only: fv_update_phys

      !INPUT PARAMETERS:
      type (MAPL_MetaComp)                   :: MAPL
      type(DynState), pointer                :: STATE
      type(ESMF_State),       intent(INOUT)  :: import
      real(FVPRC),            intent(IN   )  :: DT
      integer,  optional,     intent(OUT  )  :: RC
      logical,  optional,     intent(IN   )  :: is_weighted

      !DESCRIPTION:  This routine adds the tendencies to the state,
      !              weighted appropriately by the time step.  Temperature
      !              tendencies are pressure weighted (ie., DELP*DT/Dt).
      !              All tendencies are on the A-grid, and have an XY decomposition.

      integer               :: status
      logical               :: is_weighted_

      integer               :: II, JJ, I, J, L
      integer               :: is ,ie , js ,je , km
      integer               :: isd,ied, jsd,jed
      real(r4), allocatable :: fvQOLD(:,:,:), QTEND(:,:,:)
      real(r8), allocatable :: DPNEW(:,:,:),DPOLD(:,:,:)

      real(REAL8), allocatable :: tend_ua(:,:,:), tend_va(:,:,:)
      real(REAL8), allocatable :: tend_un(:,:,:), tend_vn(:,:,:)

      real(FVPRC), allocatable :: u_dt(:,:,:), v_dt(:,:,:), t_dt(:,:,:)

      real(r4), pointer :: tend(:,:,:)

      real(r4), pointer, dimension(:,:)   :: LONS
      real(r4), pointer, dimension(:,:)   :: LATS

      type(DynTracers)      :: qqq       ! Specific Humidity
      real(FVPRC), allocatable :: Q(:,:,:,:), CVM(:,:,:)
      integer :: n, nwat_tracers, nwat, sphum, liq_wat, ice_wat, rainwat, snowwat, graupel
      real, parameter:: c_ice = 1972.            !< heat capacity of ice at -15.C
      real, parameter:: c_liq = 4.1855e+3        !< GFS: heat capacity of water at 0C
      real, parameter:: c_vap = MAPL_CPVAP       !< 1846.
      real, parameter:: c_air = MAPL_CP

      character(len=ESMF_MAXSTR)         :: IAm="ADD_INCS"
      real(FVPRC) :: fac

      type (time_type) :: Time_Nudge

      if(present(is_weighted)) then
         is_weighted_ = is_weighted
      else
         is_weighted_ = .true.
      endif

      is = state%grid%is
      ie = state%grid%ie
      js = state%grid%js
      je = state%grid%je
      km = state%grid%npz

      isd = state%grid%isd
      ied = state%grid%ied
      jsd = state%grid%jsd
      jed = state%grid%jed

      call MAPL_Get( MAPL, LONS=LONS, LATS=LATS, RC=STATUS )
      VERIFY_(STATUS)

      ! **********************************************************************
      ! ****  Use QV from FV3 init when coldstarting idealized cases      ****
      ! **********************************************************************

      ! Determine how many water species we have
      nwat = state%vars%nwat
      nwat_tracers = 0
      if ((nwat==0) .AND. (.not. ADIABATIC)) then
         do n=1,STATE%GRID%NQ
            if (TRIM(state%vars%tracer(n)%tname) == 'Q'       ) nwat_tracers = nwat_tracers + 1
            if (TRIM(state%vars%tracer(n)%tname) == 'QLCN'    ) nwat_tracers = nwat_tracers + 1
            if (TRIM(state%vars%tracer(n)%tname) == 'QLLS'    ) nwat_tracers = nwat_tracers + 1
            if (TRIM(state%vars%tracer(n)%tname) == 'QICN'    ) nwat_tracers = nwat_tracers + 1
            if (TRIM(state%vars%tracer(n)%tname) == 'QILS'    ) nwat_tracers = nwat_tracers + 1
         enddo
         ! We must have these first 5 at a minimum
         _ASSERT(nwat_tracers == 5, 'expecting 5 water species: Q QLCN QLLS QICN QILS')
         ! Check for QRAIN, QSNOW, QGRAUPEL
         do n=1,STATE%GRID%NQ
            if (TRIM(state%vars%tracer(n)%tname) == 'QRAIN'   ) nwat_tracers = nwat_tracers + 1
            if (TRIM(state%vars%tracer(n)%tname) == 'QSNOW'   ) nwat_tracers = nwat_tracers + 1
            if (TRIM(state%vars%tracer(n)%tname) == 'QGRAUPEL') nwat_tracers = nwat_tracers + 1
         enddo
         if (nwat_tracers >= 5) nwat = 3 ! STATE has QV, QLIQ, QICE
         if (nwat_tracers == 8) nwat = 6 ! STATE has QV, QLIQ, QICE, QRAIN, QSNOW, QGRAUPEL
      endif
      if (.not. ADIABATIC) then
         _ASSERT(nwat >= 1, 'expecting water species (nwat) to match')
      endif

      select case(nwat)
      case(1)
         sphum   = 1
         liq_wat = -1
         ice_wat = -1
         rainwat = -1
         snowwat = -1
         graupel = -1
      case(3)
         sphum   = 1
         liq_wat = 2
         ice_wat = 3
         rainwat = -1
         snowwat = -1
         graupel = -1
      case(6:7)
         sphum   = 1
         liq_wat = 2
         ice_wat = 3
         rainwat = 4
         snowwat = 5
         graupel = 6
      end select

      if (nwat >= 1) then
         ALLOCATE(   Q(is:ie,js:je,1:km,nwat) )
         ALLOCATE( CVM(is:ie,js:je,1:km) )
         Q(:,:,:,:) = 0.0
         call PULL_Q ( STATE, import, qqq, NXQ, InFieldName='Q', RC=rc )
         if (DYN_COLDSTART .and. overwrite_Q .and. (.not. ADIABATIC)) then
            ! USE Q computed by FV3
            call getQ(Q(:,:,:,sphum), 'Q')
            overwrite_Q=.false.
            call WRITE_PARALLEL("Using QV from FV3 Initial Conditions")
            fac = 1.0
            call prt_maxmin('AI Q', Q(:,:,:,sphum),  is, ie, js, je, 0, km, fac)
            if ( (qqq%is_r4) .and. (associated(qqq%content_r4)) ) then
               if (size(Q(:,:,:,sphum))==size(qqq%content_r4)) qqq%content_r4 = Q(:,:,:,sphum)
            elseif (associated(qqq%content)) then
               if (size(Q(:,:,:,sphum))==size(qqq%content)) qqq%content = Q(:,:,:,sphum)
            endif
         else
            ! Grab QV from imports
            if ( (qqq%is_r4) .and. (associated(qqq%content_r4)) ) then
               if (size(Q(:,:,:,sphum))==size(qqq%content_r4)) Q(:,:,:,sphum) = qqq%content_r4
            elseif (associated(qqq%content)) then
               if (size(Q(:,:,:,sphum))==size(qqq%content)) Q(:,:,:,sphum) = qqq%content
            endif
         endif
      endif
      if (nwat >= 3) then
         ! Grab QLIQ from imports
         call PULL_Q ( STATE, import, qqq, NXQ, InFieldName='QLLS', RC=rc )
         if ( (qqq%is_r4) .and. (associated(qqq%content_r4)) ) then
            if (size(Q(:,:,:,liq_wat))==size(qqq%content_r4)) Q(:,:,:,liq_wat) = Q(:,:,:,liq_wat) + qqq%content_r4
         elseif (associated(qqq%content)) then
            if (size(Q(:,:,:,liq_wat))==size(qqq%content)) Q(:,:,:,liq_wat) = Q(:,:,:,liq_wat) + qqq%content
         endif
         call PULL_Q ( STATE, import, qqq, NXQ, InFieldName='QLCN', RC=rc )
         if ( (qqq%is_r4) .and. (associated(qqq%content_r4)) ) then
            if (size(Q(:,:,:,liq_wat))==size(qqq%content_r4)) Q(:,:,:,liq_wat) = Q(:,:,:,liq_wat) + qqq%content_r4
         elseif (associated(qqq%content)) then
            if (size(Q(:,:,:,liq_wat))==size(qqq%content)) Q(:,:,:,liq_wat) = Q(:,:,:,liq_wat) + qqq%content
         endif
         ! Grab QICE from imports
         call PULL_Q ( STATE, import, qqq, NXQ, InFieldName='QILS', RC=rc )
         if ( (qqq%is_r4) .and. (associated(qqq%content_r4)) ) then
            if (size(Q(:,:,:,ice_wat))==size(qqq%content_r4)) Q(:,:,:,ice_wat) = Q(:,:,:,ice_wat) + qqq%content_r4
         elseif (associated(qqq%content)) then
            if (size(Q(:,:,:,ice_wat))==size(qqq%content)) Q(:,:,:,ice_wat) = Q(:,:,:,ice_wat) + qqq%content
         endif
         call PULL_Q ( STATE, import, qqq, NXQ, InFieldName='QICN', RC=rc )
         if ( (qqq%is_r4) .and. (associated(qqq%content_r4)) ) then
            if (size(Q(:,:,:,ice_wat))==size(qqq%content_r4)) Q(:,:,:,ice_wat) = Q(:,:,:,ice_wat) + qqq%content_r4
         elseif (associated(qqq%content)) then
            if (size(Q(:,:,:,ice_wat))==size(qqq%content)) Q(:,:,:,ice_wat) = Q(:,:,:,ice_wat) + qqq%content
         endif
      endif
      if (nwat >= 6) then
         ! Grab RAIN from imports
         call PULL_Q ( STATE, import, qqq, NXQ, InFieldName='QRAIN', RC=rc )
         if ( (qqq%is_r4) .and. (associated(qqq%content_r4)) ) then
            if (size(Q(:,:,:,rainwat))==size(qqq%content_r4)) Q(:,:,:,rainwat) = qqq%content_r4
         elseif (associated(qqq%content)) then
            if (size(Q(:,:,:,rainwat))==size(qqq%content)) Q(:,:,:,rainwat) = qqq%content
         endif
         ! Grab SNOW from imports
         call PULL_Q ( STATE, import, qqq, NXQ, InFieldName='QSNOW', RC=rc )
         if ( (qqq%is_r4) .and. (associated(qqq%content_r4)) ) then
            if (size(Q(:,:,:,snowwat))==size(qqq%content_r4)) Q(:,:,:,snowwat) = qqq%content_r4
         elseif (associated(qqq%content)) then
            if (size(Q(:,:,:,snowwat))==size(qqq%content)) Q(:,:,:,snowwat) = qqq%content
         endif
         ! Grab GRAUPEL from imports
         call PULL_Q ( STATE, import, qqq, NXQ, InFieldName='QGRAUPEL', RC=rc )
         if ( (qqq%is_r4) .and. (associated(qqq%content_r4)) ) then
            if (size(Q(:,:,:,graupel))==size(qqq%content_r4)) Q(:,:,:,graupel) = qqq%content_r4
         elseif (associated(qqq%content)) then
            if (size(Q(:,:,:,graupel))==size(qqq%content)) Q(:,:,:,graupel) = qqq%content
         endif
      endif

      if ( (.not. ADIABATIC) .and. (DO_ADD_INCS) ) then

         ! **********************************************************************
         ! ****                      Wind Tendencies                         ****
         ! ****         Note: State Variables are on the D-Grid,             ****
         ! ****        while import Tendencies are on the A-Grid             ****
         ! **********************************************************************

         ALLOCATE( tend_ua(is:ie  ,js:je  ,km) )
         ALLOCATE( tend_va(is:ie  ,js:je  ,km) )
         ALLOCATE( tend_un(is:ie  ,js:je+1,km) )
         ALLOCATE( tend_vn(is:ie+1,js:je  ,km) )

         call ESMFL_StateGetPointerToData ( import,TEND,'DUDT',RC=STATUS )
         VERIFY_(STATUS)

         tend_ua(is:ie,js:je,1:km) = tend

         call ESMFL_StateGetPointerToData ( import,TEND,'DVDT',RC=STATUS )
         VERIFY_(STATUS)

         tend_va(is:ie,js:je,1:km) = tend

         !if (.not. HYDROSTATIC ) then
         !  call ESMFL_StateGetPointerToData ( import,TEND,'DWDT',RC=STATUS )
         !  VERIFY_(STATUS)
         !  STATE%VARS%W = STATE%VARS%W + DT*TEND(is:ie,js:je,1:km)
         !endif

         ! Put the wind tendencies on the Native Dynamics grid
         call Agrid_To_Native( tend_ua, tend_va, tend_un, tend_vn, &
              wind_increment_limiter = 800.d0/86400.d0 )


         ! Add the wind tendencies to the control variables
         STATE%VARS%U = STATE%VARS%U + DT*TEND_UN(is:ie,js:je,1:km)
         STATE%VARS%V = STATE%VARS%V + DT*TEND_VN(is:ie,js:je,1:km)

         DEALLOCATE( tend_ua )
         DEALLOCATE( tend_va )

         ! **********************************************************************
         ! ****           Compute Old Pressure Thickness                     ****
         ! **********************************************************************

         ALLOCATE( DPOLD(is:ie,js:je,km) )

         if(is_weighted_) then
            do k=1,km
               DPOLD(:,:,k) = ( state%vars%pe(:,:,k+1)-state%vars%pe(:,:,k) )
            enddo
         else
            DPOLD = 1.0
         end if

         ! **********************************************************************
         ! ****                     Update Edge Pressures                    ****
         ! **********************************************************************

         call ESMFL_StateGetPointerToData ( import,TEND,'DPEDT',RC=STATUS )
         VERIFY_(STATUS)

         STATE%VARS%PE = STATE%VARS%PE + DT*TEND

         ! **********************************************************************
         ! ****           Compute New Pressure Thickness                     ****
         ! **********************************************************************

         ALLOCATE( DPNEW(is:ie,js:je,km) )

         if(is_weighted_) then
            do k=1,km
               DPNEW(:,:,k) = ( state%vars%pe(:,:,k+1)-state%vars%pe(:,:,k) )
            enddo
         else
            DPNEW = 1.0
         end if

         ! *********************************************************************
         ! ****                  Dry Temperature Tendency                   ****
         ! ****                  ------------------------                   ****
         ! ****  Note: State  Variable is Potential Temperature T/P**kappa  ****
         ! ****        import Variable is a) D/Dt (T)     , IS_WEIGHTED=.F. ****
         ! ****                           b) D/Dt (T*DELP), IS_WEIGHTED=.T. ****
         ! *********************************************************************

         call ESMFL_StateGetPointerToData ( import,TEND,'DTDT',RC=STATUS )
         VERIFY_(STATUS)

         !if (DYN_DEBUG) then
         !   call prt_maxmin('AI PT1', STATE%VARS%PT ,  is, ie, js, je, 0, km, 1.d00, MAPL_AM_I_ROOT())
         !endif

         select case (nwat)
         case (6:7)
            CVM = (1.-( Q(:,:,:,  sphum)+Q(:,:,:,liq_wat)+Q(:,:,:,rainwat)+Q(:,:,:,ice_wat)+&
                 Q(:,:,:,snowwat)+Q(:,:,:,graupel) )               )*c_air + &
                 (Q(:,:,:,  sphum)                                  )*c_vap + &
                 (Q(:,:,:,liq_wat)+Q(:,:,:,rainwat)                 )*c_liq + &
                 (Q(:,:,:,ice_wat)+Q(:,:,:,snowwat)+Q(:,:,:,graupel))*c_ice
         case (3)
            CVM = (1.-( Q(:,:,:,  sphum)+Q(:,:,:,liq_wat)+Q(:,:,:,ice_wat) ) )*c_air + &
                 (Q(:,:,:,  sphum)                                     )*c_vap + &
                 (Q(:,:,:,liq_wat)                                     )*c_liq + &
                 (Q(:,:,:,ice_wat)                                     )*c_ice
         case default
            CVM = MAPL_CP
         end select

         ! Make previous PT into just T
         STATE%VARS%PT = STATE%VARS%PT*STATE%VARS%PKZ

         if (.not. HYDROSTATIC ) then
            ! remove old T from DZ
            STATE%VARS%DZ = STATE%VARS%DZ / STATE%VARS%PT

            ! Update T
            STATE%VARS%PT =  STATE%VARS%PT                         *DPOLD
            STATE%VARS%PT = (STATE%VARS%PT + DT*TEND*(MAPL_CP/CVM))/DPNEW

            ! update DZ with new T
            STATE%VARS%DZ = STATE%VARS%DZ * STATE%VARS%PT
         else
            ! Update T
            STATE%VARS%PT =  STATE%VARS%PT                         *DPOLD
            STATE%VARS%PT = (STATE%VARS%PT + DT*TEND*(MAPL_CP/CVM))/DPNEW
         endif

         if (DEBUG_TQ_ERRORS) then
            do L=1,KM
               do J=js,je
                  do I=is,ie
                     if ( (STATE%VARS%PT(I,J,L) > 333.0) .OR. (STATE%VARS%PT(I,J,L)/=STATE%VARS%PT(I,J,L)) .OR. &
                          (Q(I,J,L,sphum  ) < 0.0) .OR. (Q(I,J,L,sphum  )/=Q(I,J,L,sphum  )) .OR. &
                          (Q(I,J,L,liq_wat) < 0.0) .OR. (Q(I,J,L,liq_wat)/=Q(I,J,L,liq_wat)) .OR. &
                          (Q(I,J,L,ice_wat) < 0.0) .OR. (Q(I,J,L,ice_wat)/=Q(I,J,L,ice_wat)) .OR. &
                          (Q(I,J,L,rainwat) < 0.0) .OR. (Q(I,J,L,rainwat)/=Q(I,J,L,rainwat)) .OR. &
                          (Q(I,J,L,snowwat) < 0.0) .OR. (Q(I,J,L,snowwat)/=Q(I,J,L,snowwat)) .OR. &
                          (Q(I,J,L,graupel) < 0.0) .OR. (Q(I,J,L,graupel)/=Q(I,J,L,graupel)) ) then
                        print *, "T or Q  spike detected : ", STATE%VARS%PT(I,J,L)
                        print *, "  Temp  ANA|PHY  Increment : ", (DT*TEND(I,J,L)*(MAPL_CP/CVM(I,J,L)))/DPNEW(I,J,L)
                        print *, "    IN ADD_INCS inside DYN   "
                        II=I-is+1
                        JJ=J-js+1
                        print *, "  Latitude       =", LATS(II,JJ)*180.0/MAPL_PI
                        print *, "  Longitude      =", LONS(II,JJ)*180.0/MAPL_PI
                        print *, "  Pressure (mb)  =", 0.5*(STATE%VARS%PE(I,J,L+1)+STATE%VARS%PE(I,J,L))/100.0

                        print *, "  UWND =", STATE%VARS%U(I,J,L), " UINC =", DT*TEND_UN(I,J,L)
                        print *, "  VWND =", STATE%VARS%V(I,J,L), " VINC =", DT*TEND_VN(I,J,L)
                        if (nwat >= 6) then
                           print *, "  QV=", Q(I,J,L,sphum  ), "  QL=", Q(I,J,L,liq_wat), "  QI=", Q(I,J,L,ice_wat)
                           print *, "  QR=", Q(I,J,L,rainwat), "  QS=", Q(I,J,L,snowwat), "  QG=", Q(I,J,L,graupel)
                        end if
                     endif
                  end do ! IM loop
               end do ! JM loop
            end do ! LM loop
         endif

         DEALLOCATE( tend_un )
         DEALLOCATE( tend_vn )


         ! Update PKZ from hydrostatic pressures
         !  This isn't entirely necessary, FV3 overwrites this in fv_dynamics
         !  but we have to get back to PT here
         !!   call getPKZ(STATE%VARS%PKZ,STATE%VARS%PT,Q,STATE%VARS%PE,STATE%VARS%DZ,HYDROSTATIC)
         call getPKZ(STATE%VARS%PKZ,STATE%VARS%PE)

         ! Make T back into PT
         STATE%VARS%PT = STATE%VARS%PT/STATE%VARS%PKZ

         !if (DYN_DEBUG) then
         !call prt_maxmin('AI PT2', STATE%VARS%PT ,  is, ie, js, je, 0, km, 1.d00, MAPL_AM_I_ROOT())
         !endif

         DEALLOCATE (DPNEW)
         DEALLOCATE (DPOLD)

      endif ! .not. Adiabatic

      if (ALLOCATED(Q  )) DEALLOCATE( Q   )
      if (ALLOCATED(CVM)) DEALLOCATE( CVM )

      return

   end subroutine ADD_INCS

   subroutine FILLOUT3r8(export, name, V, RC)
      type (ESMF_State),  intent(inout) :: export
      character(len=*),   intent(IN   ) :: name
      real(r8),           intent(IN   ) :: V(:,:,:)
      integer, optional,  intent(  out) :: rc

      real(r8), pointer          :: CPL(:,:,:)
      integer                    :: status
      character(len=ESMF_MAXSTR) :: IAm="Fillout3r8"

      call MAPL_GetPointer(export, cpl, name, RC=STATUS)
      VERIFY_(STATUS)
      if(associated(cpl)) cpl=v
   end subroutine FILLOUT3r8

   subroutine FILLOUT3(export, name, V, RC)
      type (ESMF_State),  intent(inout) :: export
      character(len=*),   intent(IN   ) :: name
      real(r8),           intent(IN   ) :: V(:,:,:)
      integer, optional,  intent(  out) :: rc

      real(r4), pointer          :: CPL(:,:,:)
      integer                    :: status
      character(len=ESMF_MAXSTR) :: IAm="Fillout3"

      call MAPL_GetPointer(export, cpl, name, RC=STATUS)
      VERIFY_(STATUS)
      if(associated(cpl)) cpl=v
   end subroutine FILLOUT3

   subroutine FILLOUT2(export, name, V, rc)
     type (ESMF_State),  intent(inout) :: export
     character(len=*),   intent(IN   ) :: name
     real(r8),           intent(IN   ) :: V(:,:)
     integer, optional,  intent(  out) :: rc

     real(kind=4), pointer      :: CPL(:,:)
     integer                    :: status
     character(len=ESMF_MAXSTR) :: IAm="Fillout2"

     call MAPL_GetPointer(export, cpl, name, RC=STATUS)
     VERIFY_(STATUS)
     if(associated(cpl)) cpl=v

     return
   end subroutine FILLOUT2

   subroutine Energetics (ua,va,thv,ple,delp,pk,phiS,keint,peint,teint,ke,cpt,gze)

      real(8), optional, intent(out) ::   ke(:,:,:)
      real(8), optional, intent(out) ::  cpt(:,:,:)
      real(8), optional, intent(out) ::  gze(:,:,:)
      real(8)   ua(:,:,:)
      real(8)   va(:,:,:)
      real(8)  thv(:,:,:)
      real(8)  ple(:,:,:)
      real(8) delp(:,:,:)
      real(8)   pk(:,:,:)
      real(8)   keint(:,:)
      real(8)   peint(:,:)
      real(8)   teint(:,:)
      real(8) phiS(:,:)

      real(8) kinetic, potential
      integer i,ifirst,ilast
      integer j,jfirst,jlast
      integer km,k

      real(8), allocatable ::   pke(:,:,:)
      real(8), allocatable ::  phiT(:,:)

      ifirst = lbound( ua,1 )
      ilast  = ubound( ua,1 )
      jfirst = lbound( ua,2 )
      jlast  = ubound( ua,2 )
      km     = ubound( ua,3 )

      allocate( pke  ( ifirst:ilast, jfirst:jlast , 1:km+1 ) )
      allocate( phiT ( ifirst:ilast, jfirst:jlast ) )

      ! Compute Model Edge Heights
      pke  = ple**kappa
      phiT = phiS
      if( present(gze) ) gze(:,:,km+1) = phiS
      do k=km,1,-1
         phiT = phiT + cp*thv(:,:,k)*( pke(:,:,k+1)-pke(:,:,k) )
         if( present(gze) ) gze(:,:,k) = phiT
      enddo

      ! Compute Energetics:  Cp*Tv + K + PHI
      keint = 0.0
      peint = 0.0
      do k=1,km
         do j=jfirst,jlast
            do i=ifirst,ilast
               kinetic      = 0.5_r8*( ua(i,j,k)**2 + va(i,j,k)**2 )
               potential    =  cp*thv(i,j,k)*pk(i,j,k)
               keint(i,j)   =   keint(i,j) +   kinetic  *delp(i,j,k)
               peint(i,j)   =   peint(i,j) +   potential*delp(i,j,k)
               if( present(ke)  )  ke(i,j,k) = kinetic
               if( present(cpt) ) cpt(i,j,k) = potential
            enddo
         enddo
      enddo
      keint(:,:) =    keint(:,:)/grav
      peint(:,:) =    peint(:,:)/grav
      teint(:,:) = (phiS(:,:)*ple(:,:,km+1)-phiT(:,:)*ple(:,:,1))/grav

      deallocate ( pke  )
      deallocate ( phiT )

      return
   end subroutine Energetics

   !BOP
   !IROUTINE: Finalize

   !DESCRIPTION: Writes restarts and cleans-up through MAPL\_GenericFinalize and
   !   deallocates memory from the Private Internal state.

   !INTERFACE:
   subroutine Finalize(gc, import, export, clock, rc)

      !ARGUMENTS:
      type (ESMF_GridComp) :: gc
      type (ESMF_State) :: import
      type (ESMF_State) :: export
      type (ESMF_Clock) :: clock
      integer, intent(out) :: rc
      !EOP

      type (DYN_wrap) :: wrap
      type (DynState), pointer  :: STATE

      character(len=ESMF_MAXSTR)        :: IAm
      character(len=ESMF_MAXSTR)        :: comp_name
      integer                           :: status

      type (MAPL_MetaComp),     pointer :: MAPL
      type (ESMF_Config)                :: cf

      Iam = "Finalize"
      call ESMF_GridCompGet( gc, name=comp_name, config=cf, RC=STATUS )
      VERIFY_(STATUS)
      Iam = trim(comp_name) // Iam

      ! Retrieve the pointer to the state
      call MAPL_GetObjectFromGC (gc, MAPL,  RC=STATUS )
      VERIFY_(STATUS)

      call MAPL_TimerOn(MAPL,"TOTAL")
      call MAPL_TimerOn(MAPL,"FINALIZE")

      ! Retrieve the pointer to the state
      call ESMF_UserCompGetInternalState(gc, 'DYNstate', wrap, status)
      VERIFY_(STATUS)

      state => wrap%dyn_state

      call DynFinalize( STATE )

      ! Call Generic Finalize
      call MAPL_TimerOff(MAPL,"FINALIZE")
      call MAPL_TimerOff(MAPL,"TOTAL")

      call MAPL_GenericFinalize ( gc, import, export, clock,  RC=STATUS)
      VERIFY_(STATUS)

      _RETURN(_SUCCESS)

   contains

      subroutine PRINT_TIMES(TIMES,DAYS)
         integer(kind=8), intent(INOUT) :: TIMES(:,:)
         real(r8),        intent(IN   ) :: DAYS
         TIMES = 0
         return
      end subroutine PRINT_TIMES

   end subroutine FINALIZE

   subroutine get_slp ( km,ps,phis,slp,pe,pk,tv,H1000,H850,H500)
      implicit   none

      integer  km
      real(r8)   pk(km)    ! layer-mean P**kappa
      real(r8)   tv(km)    ! layer-mean virtual Temperature
      real(r8)   pe(km+1)  ! press at layer edges (Pa)
      real(r8)   ps        ! surface pressure (Pa)
      real(r8) phis        ! surface geopotential
      real(r8)  slp        ! sea-level pressure (hPa)
      real(r8)  H1000      ! 1000mb height
      real(r8)  H850       !  850mb height
      real(r8)  H500       !  500mb height
      real(r8)  tstar                 ! extrapolated temperature (K)
      real(r8) p_bot
      real(r8) tref                   ! Reference virtual temperature (K)
      real(r8) pref                   ! Reference pressure level (Pa)
      real(r8) pkref                  ! Reference pressure level (Pa) ** kappa
      real(r8) dp1, dp2

      real(r8), parameter :: gamma    = 6.5e-3
      real(r8), parameter :: p_offset = 15000.
      real(r8), parameter :: gg       = gamma/MAPL_GRAV

      real(r8), parameter :: factor   = MAPL_grav / ( MAPL_Rgas * gamma )
      real(r8), parameter :: yfactor  = MAPL_Rgas * gg

      integer k_bot, k, k1, k2

      p_bot = ps - p_offset
      k_bot = -1

      do k = km, 2, -1
         if ( pe(k+1) .lt. p_bot ) then
            k_bot = k
            exit
         endif
      enddo

      k1    = k_bot - 1
      k2    = k_bot
      dp1   = pe(k_bot)   - pe(k_bot-1)
      dp2   = pe(k_bot+1) - pe(k_bot)
      pkref = ( pk(k1)*dp1 + pk(k2)*dp2 ) / (dp1+dp2)
      tref  = ( tv(k1)*dp1 + tv(k2)*dp2 ) / (dp1+dp2)
      pref  = 0.5 * ( pe(k_bot+1) + pe(k_bot-1) )
      tstar = tref*( ps/pref )**yfactor

      slp   = ps*( 1.0+gg*phis/tstar )**factor
      H1000 = (phis/MAPL_grav) - (tstar/gamma)*((100000.0/ps)**(1./factor)-1.0)
      H850  = (phis/MAPL_grav) - (tstar/gamma)*(( 85000.0/ps)**(1./factor)-1.0)
      H500  = (phis/MAPL_grav) - (tstar/gamma)*(( 50000.0/ps)**(1./factor)-1.0)

      return
   end subroutine get_slp

   subroutine VertInterp(v2,v3,ple,pp,positive_definite,rc)
      real(r4), intent(OUT) :: v2(:,:)
      real(r8), intent(IN ) :: v3(:,:,:)
      real(r8), intent(IN ) :: ple(:,:,:)
      real    , intent(IN ) :: pp
      logical, optional, intent(IN ) :: positive_definite
      integer, optional, intent(OUT) :: rc

      real, dimension(size(v2,1),size(v2,2)) :: al,PT,PB
      integer km
      logical edge

      character*(10) :: Iam='VertInterp'

      km   = size(ple,3)-1
      edge = size(v3,3)==km+1

      _ASSERT(edge .or. size(v3,3)==km,'needs informative message')

      v2   = MAPL_UNDEF

      if(EDGE) then
         pb   = ple(:,:,km+1)
         do k=km,1,-1
            pt = ple(:,:,k)
            if(all(pb<pp)) exit
            where(pp>pt .and. pp<=pb)
               al = (pb-pp)/(pb-pt)
               v2 = v3(:,:,k)*al + v3(:,:,k+1)*(1.0-al)
            end where
            pb = pt
         end do
      else
         pb = 0.5*(ple(:,:,km)+ple(:,:,km+1))
         do k=km,2,-1
            pt = 0.5*(ple(:,:,k-1)+ple(:,:,k))
            if(all(pb<pp)) exit
            where( (pp>pt.and.pp<=pb) )
               al = (pb-pp)/(pb-pt)
               v2 = v3(:,:,k-1)*al + v3(:,:,k)*(1.0-al)
            end where
            pb = pt
         end do
         pt = 0.5*(ple(:,:,km)+ple(:,:,km-1))
         pb = 0.5*(ple(:,:,km)+ple(:,:,km+1))
         where( (pp>pb.and.pp<=ple(:,:,km+1)) )
            v2 = v3(:,:,km)
         end where
      end if

      if (present(positive_definite)) then
         if (positive_definite) then
            where (v2 < tiny(0.0))
               v2 = 0.0
            endwhere
         endif
      endif

      _RETURN(_SUCCESS)
   end subroutine VertInterp

   !BOP
   !IROUTINE: Coldstart

   !DESCRIPTION:
   !   Routine to coldstart from an isothermal state of rest.
   !   The temperature can be specified in the config, otherwise
   !   it is 300K. The surface pressure is assumed to be 1000 hPa.

   !INTERFACE:
   subroutine Coldstart(gc, import, export, clock, rc)

      USE sw, only : sw_phis=>surface_geopotential
      USE sw, only : sw_hght=>height
      USE sw, only : sw_uwnd=>u_wind
      USE sw, only : sw_vwnd=>v_wind
      USE jw, only : temperature, u_wind, v_wind, surface_geopotential
      USE jw, only : tracer_q, tracer_q1_q2, tracer_q3
      USE testcases_3_4_5_6, only : advection, Rossby_Haurwitz, mountain_Rossby, gravity_wave

      !ARGUMENTS:
      type(ESMF_GridComp), intent(inout) :: gc
      type(ESMF_State), intent(inout) :: import
      type(ESMF_State), intent(inout) :: export
      type (ESMF_Clock), intent(inout) :: clock
      integer, intent(out), optional :: rc
      !EOP

      integer                           :: status

      type(ESMF_State) :: internal

      real(REAL8), pointer :: AK1(:), BK1(:), AK(:), BK(:)
      real(REAL8), pointer :: U(:,:,:), V(:,:,:), PT(:,:,:)
      real(REAL8), pointer :: PE1(:,:,:), PE(:,:,:), PKZ(:,:,:)
      real(REAL4), pointer :: phis(:,:)
      real(REAL4), pointer :: lons(:,:), lats(:,:)

      integer :: i, j, k, n, L
      ! integer :: IS, IE, JS, JE, KS, KE, IM, JM, KM, LS
      integer :: is, ie, js, je, ks, ke, im, jm, km, ls
      integer :: case_id, case_rotation, case_tracers

      real :: T0
      real(REAL8) :: dummy_1, dummy_2, dummy_3, dummy_4, dummy_5, dummy_6
      real(REAL8) :: dz, ztop, height, pressure
      real(REAL8) :: LONc,LATc
      real(REAL8) :: eta, eta_top, rot_ang
      real(REAL8) :: ptop, pint
      real(REAL8), allocatable :: PS(:,:)

      type(DYN_wrap) :: wrap
      type(DynState), pointer :: state
      type(DynGrid),  pointer :: grid

      logical :: perturb
      logical :: ak_is_missing = .false.
      logical :: bk_is_missing = .false.
      integer :: FV3_STANDALONE
      logical :: isPresent

      ! Tracer Stuff
      real(REAL4), pointer :: tracer(:,:,:)
      real(REAL8), allocatable :: Q5(:,:,:)
      real(REAL8), allocatable :: Q6(:,:,:)
      type(ESMF_Grid) :: esmfgrid
      type(ESMF_FieldBundle) :: tradv_bundle
      character(len=ESMF_MAXSTR) :: fieldname
      character(len=ESMF_MAXSTR) :: string
      real(REAL8), parameter :: r0_6=0.6
      real(REAL8), parameter :: r1_0=1.0

      call ESMF_UserCompGetInternalState(gc, 'DYNstate', wrap, _RC)
      state => wrap%dyn_state
      grid  => state%grid ! direct handle to grid

      call MAPL_GridCompGetResource(gc, "T0", T0, default=273., _RC)
      call MAPL_GridCompGetInternalState(gc, internal, _RC)
      call MAPL_GridCompGet(gc, grid=esmfgrid, _RC)
      call MAPL_GridGet(esmfgrid, latitudes=lats, longitudes=lons, _RC)

      if (FV_Atm(1)%flagstruct%grid_type == 4) then
         ! Doubly-Period setup based on first LAT/LON coordinate
         lons(:,:) =  0.0
         lats(:,:) = 15.0*PI/180.0
      endif

      call MAPL_GetPointer(internal, U, "U", _RC) ! A-Grid U Wind
      is = lbound(U,1); ie = ubound(U,1)
      js = lbound(U,2); je = ubound(U,2)
      ks = lbound(U,3); ke = ubound(U,3)
      km = ke-ks+1
      call MAPL_GetPointer(internal, V, "V", _RC) ! A-Grid V Wind
      call MAPL_GetPointer(internal, PT, "PT", _RC) ! potential temperature
      call MAPL_GetPointer(internal, PE1, "PE", _RC) ! edge pressures - 1 based
      PE(is:ie, js:je, 0:km) => PE1(is:ie, js:je, 1:km+1)
      call MAPL_GetPointer(internal, PKZ , "PKZ", _RC) ! presssure ^ kappa at mid-layers
      call MAPL_GetPointer(internal, ak1, "AK" , _RC) ! AK for vertical coordinate - 1 based
      ak(0:km) => ak1(1:km+1)
      call MAPL_GetPointer(internal, bk1, "BK" , _RC) ! BK for vertical coordinate - 1 based
      bk(0:km) => bk1(1:km+1)
      call MAPL_GetPointer(import, phis, "PHIS", _RC) ! surface geopotential

      U = 0.0

      ALLOCATE( PS(is:ie,js:je) )

      call MAPL_GridCompGetResource(gc, "IM", IM, default=0, _RC)
      call MAPL_GridCompGetResource(gc, "JM", JM, default=0, _RC)

      if (KM<=2) then   ! Shallow Water

         call MAPL_GridCompGetResource(gc, "CASE_ID", case_id, default=1, _RC)
         DYN_CASE = case_id

         do j=js,je
            do i=is,ie
               LONc = lons(i,j)
               LATc = lats(i,j)
               U(i,j,1)  = sw_uwnd(LONc,LATc,case_id)
               V(i,j,1)  = sw_vwnd(LONc,LATc,case_id)
               PE(i,j,0) = sw_phis(LONc,LATc,case_id)
               PE(i,j,1) = sw_hght(LONc,LATc,case_id)
               phis(i,j) = PE(i,j,0)
            enddo
         enddo

      else              ! 3-D Baroclinic

         U(is:ie,js:je,KE) = .001*abs(lats(:,:))
         V = 0.0

         ! pchakrab - TODO: read ak from resource file, if present
         ! call ESMF_ConfigFindLabel( cf, 'AK:', isPresent=isPresent, rc = status )
         ! VERIFY_(STATUS)
         ! if (isPresent) then
         !    do L = 0, SIZE(AK)-1
         !       call ESMF_ConfigNextLine  ( CF, rc=STATUS )
         !       call ESMF_ConfigGetAttribute( cf, AK(L), rc = status )
         !       VERIFY_(STATUS)
         !    enddo
         ! else
         ak_is_missing = .true.
         ! endif

         ! pchakrab - TODO: read bk from resource file, if present
         ! call ESMF_ConfigFindLabel( cf, 'BK:', isPresent=isPresent, rc = status )
         ! VERIFY_(STATUS)
         ! if (isPresent) then
         !    do L = 0, SIZE(bk)-1
         !       call ESMF_ConfigNextLine  ( CF, rc=STATUS )
         !       call ESMF_ConfigGetAttribute( cf, BK(L), rc = status )
         !       VERIFY_(STATUS)
         !    enddo
         ! else
         bk_is_missing = .true.
         ! endif

         if (ak_is_missing .or. bk_is_missing) call set_eta(km, ls, ptop, pint, AK, BK)
         _ASSERT(any(AK /= 0.0) .or. any(BK /= 0.0), "needs informative message")

         do L = lbound(PE,3), ubound(PE,3)
            PE(:,:,L) = AK(L) + BK(L)*MAPL_P00
         enddo
         PKZ = 0.5*( PE(:,:,lbound(PE,3):ubound(PE,3)-1) + PE(:,:,lbound(PE,3)+1:ubound(PE,3)) )
         PKZ = PKZ**MAPL_KAPPA
         PT = T0/PKZ

         ! Check if running standalone model
         call MAPL_GridCompGetResource(gc, "FV3_STANDALONE", FV3_STANDALONE, default=0, _RC)

         ! 3D Baroclinic Test Cases
         call MAPL_GridCompGetResource(gc, "CASE_ID", case_id, default=0, _RC)
         call MAPL_GridCompGetResource(gc, "CASE_ROTATION", case_rotation, default=0 , _RC)
         call MAPL_GridCompGetResource(gc, "CASE_TRACERS" , case_tracers, default=1234, _RC)
         DYN_CASE = case_id

         write(string,'(A,I5,A)') "Initializing CASE_ID ", case_id, " in FVcubed:"
         call WRITE_PARALLEL( trim(string) )

         ! Parse case_rotation
         if (case_rotation == -1) rot_ang =  0
         if (case_rotation ==  0) rot_ang =  0
         if (case_rotation ==  1) rot_ang = 15
         if (case_rotation ==  2) rot_ang = 30
         if (case_rotation ==  3) rot_ang = 45
         if (case_rotation ==  4) rot_ang = 60
         if (case_rotation ==  5) rot_ang = 75
         if (case_rotation ==  6) rot_ang = 90
         if (case_rotation == -1) then
            grid%f_coriolis_angle = -999
         else
            grid%f_coriolis_angle = rot_ang*PI/180.0
         endif

         if (case_id == 1) then ! Steady State

            perturb = .false.
            do k=ks,ke
               eta = 0.5*( (ak(k-1)+ak(k))/1.e5 + bk(k-1)+bk(k) )
               do j=js,je
                  do i=is,ie
                     LONc = lons(i,j)
                     LATc = lats(i,j)
                     U(i,j,k) = u_wind(LONc,LATc,eta,perturb,rot_ang)
                     V(i,j,k) = v_wind(LONc,LATc,eta,perturb,rot_ang)
                     if (k==KS) phis(i,j) = surface_geopotential(LONc,LATc,rot_ang)
                     PT(i,j,k) = temperature(LONc,LATc,eta,rot_ang)
                  enddo
               enddo
            enddo
            PT = PT/PKZ

         elseif (case_id == 2) then ! Baroclinic Wave

            perturb = .true.
            do k=ks,ke
               eta = 0.5*( (ak(k-1)+ak(k))/1.e5 + bk(k-1)+bk(k) )
               do j=js,je
                  do i=is,ie
                     LONc = lons(i,j)
                     LATc = lats(i,j)
                     U(i,j,k) = u_wind(LONc,LATc,eta,perturb,rot_ang)
                     V(i,j,k) = v_wind(LONc,LATc,eta,perturb,rot_ang)
                     if (k==KS) phis(i,j) = surface_geopotential(LONc,LATc,rot_ang)
                     PT(i,j,k) = temperature(LONc,LATc,eta,rot_ang)
                     !if (grid_type==4) then
                     !  if (k==KS) then
                     !     T_PERTURB = (SIN(PI*FLOAT(i-1)/FLOAT(IE-IS))**4.0) * &
                     !                 (SIN(PI*FLOAT(j-1)/FLOAT(JE-JS))**4.0)
                     !     print*, i, j, T_PERTURB
                     !     PT(i,j,k) = PT(i,j,k) + T_PERTURB
                     !  endif
                     !endif
                  enddo
               enddo
            enddo
            PT = PT/PKZ

         elseif (case_id == 3) then ! Advection

            !PURE_ADVECTION = .true.

            allocate(Q5(is:ie, js:je, 0:KM-1), _STAT)
            allocate(Q6(is:ie, js:je, 0:KM-1), _STAT)

            ztop = 12000.0
            dz   = ztop/KM
            do k=ks,ke
               height = (ztop - 0.5*dz) - (k)*dz  ! Layer middle height
               do j=js,je
                  do i=is,ie
                     LONc = lons(i,j)
                     LATc = lats(i,j)
                     call  advection('56', LONc, LATc, height, rot_ang,  &
                          dummy_1, dummy_2, dummy_3, dummy_4, &
                          PS(i,j), Q5(i,j,k), Q6(i,j,k))
                     U(i,j,k)  = dummy_1
                     V(i,j,k)  = dummy_2
                     PT(i,j,k) = dummy_3
                     phis(i,j) = dummy_4
                  enddo
               enddo
            enddo
            do L=lbound(PE,3),ubound(PE,3)
               PE(:,:,L) = AK(L) + BK(L)*PS(:,:)
            enddo

            do k=ks,ke
               do j=js,je
                  do i=is,ie
                     PKZ(i,j,k) = ( (PE(i,j,k+1)**kappa) - (PE(i,j,k)**kappa) ) /  &
                            ( kappa*(log(PE(i,j,k+1))-log(PE(i,j,k))) )
                  enddo
               enddo
            enddo

            PT = PT/PKZ

         elseif (case_id == 4) then ! 3D Rossby-Haurwitz

            do j=js,je
               do i=is,ie
                  LONc = lons(i,j)
                  LATc = lats(i,j)
                  pressure = 500.
                  call Rossby_Haurwitz(LONc,LATc, pressure, dummy_1, dummy_2, dummy_3, dummy_4, PS(i,j))
                  U(i,j,1)  = dummy_1
                  V(i,j,1)  = dummy_2
                  PT(i,j,1) = dummy_3
                  phis(i,j) = dummy_4
               enddo
            enddo
            do L=lbound(PE,3),ubound(PE,3)
               PE(:,:,L) = AK(L) + BK(L)*PS(:,:)
            enddo
            do k=ks,ke
               do j=js,je
                  do i=is,ie
                     LONc = lons(i,j)
                     LATc = lats(i,j)
                     pressure = 0.5*(PE(i,j,k)+PE(i,j,k+1))
                     call Rossby_Haurwitz(LONc,LATc, pressure, dummy_1, dummy_2, dummy_3, dummy_4, PS(i,j))
                     U(i,j,k)  = dummy_1
                     V(i,j,k)  = dummy_2
                     PT(i,j,k) = dummy_3
                     phis(i,j) = dummy_4
                  enddo
               enddo
            enddo

            do k=ks,ke
               do j=js,je
                  do i=is,ie
                     PKZ(i,j,k) = ( (PE(i,j,k+1)**kappa) - (PE(i,j,k)**kappa) ) /  &
                          ( kappa*(log(PE(i,j,k+1))-log(PE(i,j,k))) )
                  enddo
               enddo
            enddo
            PT = PT/PKZ

         elseif (case_id == 5) then ! Mountain-Induced Rossby Wave

            do k=ks,ke
               do j=js,je
                  do i=is,ie
                     LONc = lons(i,j)
                     LATc = lats(i,j)
                     pressure = 0.5*(PE(i,j,k)+PE(i,j,k+1))
                     call mountain_Rossby(case_rotation,LONc,LATc, pressure, dummy_1, dummy_2, dummy_3, dummy_4, PS(i,j))
                     U(i,j,k)  = dummy_1
                     V(i,j,k)  = dummy_2
                     PT(i,j,k) = dummy_3
                     phis(i,j) = dummy_4
                  enddo
               enddo
            enddo
            do L=lbound(PE,3),ubound(PE,3)
               PE(:,:,L) = AK(L) + BK(L)*PS(:,:)
            enddo

            do k=ks,ke
               do j=js,je
                  do i=is,ie
                     PKZ(i,j,k) = ( (PE(i,j,k+1)**kappa) - (PE(i,j,k)**kappa) ) /  &
                          ( kappa*(log(PE(i,j,k+1))-log(PE(i,j,k))) )
                  enddo
               enddo
            enddo

            PT = PT/PKZ

         elseif (case_id == 6) then ! Gravity Waves

            ! case_rotation index has different meaning for this test
            if (case_rotation < 3) then
               grid%f_coriolis_angle = -999
            else
               grid%f_coriolis_angle = 0.0
            endif
            ! Get ICs
            ztop = 10000.d0
            dz   = ztop/KM
            do k=ks,ke
               height = (ztop - 0.5d0*dz) - (k)*dz  ! Layer middle height
               do j=js,je
                  do i=is,ie
                     LONc = lons(i,j)
                     LATc = lats(i,j)
                     call gravity_wave(case_rotation, LONc,LATc, height, dummy_1, dummy_2, dummy_3, dummy_4, PS(i,j))
                     U(i,j,k)  = dummy_1
                     V(i,j,k)  = dummy_2
                     PT(i,j,k) = dummy_3
                     phis(i,j) = dummy_4
                  enddo
               enddo
            enddo
            ! Reconstruct Edge Pressures and AK BK arrays for rotation=0, otherwise use values from set_eta which are OK
            if (case_rotation == 0) then
               PTOP = 27381.905d0
               do k=lbound(PE,3),ubound(PE,3)
                  height = ztop - k*dz  ! Layer edge height
                  do j=js,je
                     do i=is,ie
                        LONc = lons(i,j)
                        LATc = lats(i,j)
                        call gravity_wave(case_rotation, LONc,LATc, height, dummy_1, dummy_2, dummy_3, dummy_4, dummy_5, pressure=dummy_6)
                        PE(i,j,k) = dummy_6
                        eta     = PE(i,j,k)/PS(i,j)
                        eta_top = PTOP/PS(i,j)
                        BK(k) = (eta - eta_top)/(1.d0 - eta_top)
                        AK(k) = 100000.d0 * (eta - BK(k))
                     enddo
                  enddo
               enddo
            endif
            ! Update PE, PKZ and PT
            do L=lbound(PE,3),ubound(PE,3)
               PE(:,:,L) = AK(L) + BK(L)*PS(:,:)
            enddo

            do k=ks,ke
               do j=js,je
                  do i=is,ie
                     PKZ(i,j,k) = ( (PE(i,j,k+1)**kappa) - (PE(i,j,k)**kappa) ) /  &
                          ( kappa*(log(PE(i,j,k+1))-log(PE(i,j,k))) )
                  enddo
               enddo
            enddo

            PT = PT/PKZ

         endif ! case_id

         ! Parse Tracers
         if (FV3_STANDALONE /= 0) then
            call ESMF_StateGet(import, "TRADV", tradv_bundle, _RC)
            call MAPL_GridCompGet(gc, grid=esmfgrid, _RC)

            allocate(tracer(is:ie, js:je, 1:KM), _STAT)
            tracer(:,:,:)  = 0.0
            fieldname = 'Q'
            call addTracer(state, tradv_bundle, tracer, esmfgrid, fieldname)

            if (case_tracers /= 1234) then

               do n=1,case_tracers
                  tracer(:,:,:) = 0.0
                  write(fieldname, "('Q',i3.3)") n
                  call addTracer(state, tradv_bundle, tracer, esmfgrid, fieldname)
               enddo

            else

               !-------------------
               !     tracer q1
               !-------------------
               tracer(:,:,:) = 0.0
               do k=ks,ke
                  eta = 0.5*( (ak(k-1)+ak(k))/1.e5 + bk(k-1)+bk(k) )
                  do j=js,je
                     do i=is,ie
                        LONc = lons(i,j)
                        LATc = lats(i,j)
                        dummy_1 = tracer_q1_q2(LONc,LATc,eta,rot_ang,r0_6)
                        tracer(i,j,k) = dummy_1
                     enddo
                  enddo
               enddo
               fieldname = 'Q1'
               call addTracer(state, tradv_bundle, tracer, esmfgrid, fieldname)

               !-------------------
               !     tracer q2
               !-------------------
               do k=ks,ke
                  eta = 0.5*( (ak(k-1)+ak(k))/1.e5 + bk(k-1)+bk(k) )
                  do j=js,je
                     do i=is,ie
                        LONc = lons(i,j)
                        LATc = lats(i,j)
                        dummy_1 = tracer_q1_q2(LONc,LATc,eta,rot_ang,r1_0)
                        tracer(i,j,k) = dummy_1
                     enddo
                  enddo
               enddo
               fieldname = 'Q2'
               call addTracer(state, tradv_bundle, tracer, esmfgrid, fieldname)

               !-------------------
               !     tracer q3
               !-------------------
               do k=ks,ke
                  eta = 0.5*( (ak(k-1)+ak(k))/1.e5 + bk(k-1)+bk(k) )
                  do j=js,je
                     do i=is,ie
                        LONc = lons(i,j)
                        LATc = lats(i,j)
                        dummy_1 = tracer_q3(LONc,LATc,eta,rot_ang)
                        tracer(i,j,k) = dummy_1
                     enddo
                  enddo
               enddo
               fieldname = 'Q3'
               call addTracer(state, tradv_bundle, tracer, esmfgrid, fieldname)

               !-------------------
               !     tracer q4
               !-------------------
               tracer(:,:,:)  = 1.0_r4
               fieldname = 'Q4'
               call addTracer(state, tradv_bundle, tracer, esmfgrid, fieldname)

               !-------------------
               !     tracer q5
               !-------------------
               if (allocated(Q5)) then
                  tracer(:,:,:)  = Q5(:,:,:)
                  fieldname = 'Q5'
                  call addTracer(state, tradv_bundle, tracer, esmfgrid, fieldname)
                  deallocate(Q5, _STAT)
               endif

               !-------------------
               !     tracer q6
               !-------------------
               if (allocated(Q6)) then
                  tracer(:,:,:)  = Q6(:,:,:)
                  fieldname = 'Q6'
                  call addTracer(state, tradv_bundle, tracer, esmfgrid, fieldname)
                  deallocate(Q6, _STAT)
               endif

            endif

            deallocate(tracer, _STAT)

         endif ! parse tracers

      endif

      deallocate(PS)

      DYN_COLDSTART=.true.

      _RETURN(_SUCCESS)

   end subroutine Coldstart

#ifdef MY_SET_ETA
   subroutine set_eta(km, ptop, ak, bk)

      integer,  intent(in   )::  km          ! vertical dimension
      real(REAL8),   intent(  out):: ptop         ! model top (Pa)
      real(REAL8),   intent(inout):: ak(km+1)
      real(REAL8),   intent(inout):: bk(km+1)

      ! local
      real(REAL8) a20_01(21),b20_01(21)      ! NCAR Colloquium 20-levels N=0.01
      real(REAL8) a20_0178(21),b20_0178(21)  ! NCAR Colloquium 20-levels N=0.0178
      real(REAL8) a26(27),b26(27)            ! NCAR Colloquium 26-levels
      real(REAL8) a72(73), b72(73)           ! GEOS-5 72-levels
      real(REAL8) a137(138), b137(138)       ! GEOS-5 137-levels

      real(REAL8) :: p0=1000.E2
      real(REAL8) :: pc=200.E2
      real(REAL8) pt, pint, lnpe, dlnp
      real(REAL8) press(km+1)
      integer  k, ks

      data a20_01 / &
           0.27381905404907E+05,  0.26590539035976E+05,  0.25752394878279E+05,  0.24865429808716E+05, &
           0.23927536347865E+05,  0.22936541085572E+05,  0.21890203071294E+05,  0.20786212168493E+05, &
           0.19622187372385E+05,  0.18395675090318E+05,  0.17104147384052E+05,  0.15745000173179E+05, &
           0.14315551398919E+05,  0.12813039147516E+05,  0.11234619732416E+05,  0.95773657344247E+04, &
           0.78382639990006E+04,  0.60142135898353E+04,  0.41020236978492E+04,  0.20984115047143E+04, &
           0.00000000000000E+00 /

      data b20_01 / &
           0.00000000000000E+00,  0.28901070149364E-01,  0.59510487036309E-01,  0.91902866472543E-01, &
           0.12615517459290E+00,  0.16234678535331E+00,  0.20055953931639E+00,  0.24087780374962E+00, &
           0.28338853406205E+00,  0.32818133660555E+00,  0.37534853286773E+00,  0.42498522508382E+00, &
           0.47718936329560E+00,  0.53206181388604E+00,  0.58970642961892E+00,  0.65023012121324E+00, &
           0.71374293048299E+00,  0.78035810507338E+00,  0.85019217482527E+00,  0.92336502980036E+00, &
           0.10000000000000E+01 /

      data a20_0178 / &
           0.32021324453921E+05,  0.31137565415634E+05,  0.30202026400316E+05,  0.29211673587770E+05, &
           0.28163295404433E+05,  0.27053492108706E+05,  0.25878664766072E+05,  0.24635003578258E+05, &
           0.23318475528610E+05,  0.21924811303582E+05,  0.20449491447964E+05,  0.18887731708932E+05, &
           0.17234467521390E+05,  0.15484337584307E+05,  0.13631666474783E+05,  0.11670446243450E+05, &
           0.95943169315531E+04,  0.73965459465018E+04,  0.50700062290314E+04,  0.26071531411601E+04, &
           0.00000000000000E+00 /

      data b20_0178 / &
           0.00000000000000E+00,  0.27599078219223E-01,  0.56815203138214E-01,  0.87743118501982E-01, &
           0.12048311914891E+00,  0.15514137625266E+00,  0.19183028162025E+00,  0.23066881216269E+00, &
           0.27178291572025E+00,  0.31530591949337E+00,  0.36137896240390E+00,  0.41015145278854E+00, &
           0.46178155290889E+00,  0.51643669184922E+00,  0.57429410846515E+00,  0.63554142614418E+00, &
           0.70037726124166E+00,  0.76901186716541E+00,  0.84166781619770E+00,  0.91858072126555E+00, &
           0.10000000000000E+01 /


      data a26 / &
           219.4067,   489.5209,   988.2418,  1805.2010,  2983.7240,  4462.3340,   &
           6160.5870,  7851.2430,  7731.2710,  7590.1310,  7424.0860,   &
           7228.7440,  6998.9330,  6728.5740,  6410.5090,  6036.3220,   &
           5596.1110,  5078.2250,  4468.9600,  3752.1910,  2908.9490,   &
           2084.739,   1334.443,    708.499,   252.1360,  0.0, 0.0     /

      data b26 / &
           0.0, 0.0, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,&
           0.0000000, 0.01505309, 0.03276228, 0.05359622, 0.07810627,      &
           0.1069411, 0.1408637, 0.1807720, 0.2277220, 0.2829562,       &
           0.3479364, 0.4243822, 0.5143168, 0.6201202, 0.7235355,       &
           0.8176768, 0.8962153, 0.9534761, 0.9851122, 1.0000000        /

      data a72 / &
           1.0000000,       2.0000002,       3.2700005,       4.7585009,       6.6000011, &
           8.9345014,       11.970302,       15.949503,       21.134903,       27.852606, &
           36.504108,       47.580610,       61.677911,       79.513413,       101.94402, &
           130.05102,       165.07903,       208.49704,       262.02105,       327.64307, &
           407.65710,       504.68010,       621.68012,       761.98417,       929.29420, &
           1127.6902,       1364.3402,       1645.7103,       1979.1604,       2373.0405, &
           2836.7806,       3381.0007,       4017.5409,       4764.3911,       5638.7912, &
           6660.3412,       7851.2316,       9236.5722,       10866.302,       12783.703, &
           15039.303,       17693.003,       20119.201,       21686.501,       22436.301, &
           22389.800,       21877.598,       21214.998,       20325.898,       19309.696, &
           18161.897,       16960.896,       15625.996,       14290.995,       12869.594, &
           11895.862,       10918.171,       9936.5219,       8909.9925,       7883.4220, &
           7062.1982,       6436.2637,       5805.3211,       5169.6110,       4533.9010, &
           3898.2009,       3257.0809,       2609.2006,       1961.3106,       1313.4804, &
           659.37527,       4.8048257,       0.0000000 /

      data b72 / &
           0.0000000,       0.0000000,       0.0000000,       0.0000000,       0.0000000, &
           0.0000000,       0.0000000,       0.0000000,       0.0000000,       0.0000000, &
           0.0000000,       0.0000000,       0.0000000,       0.0000000,       0.0000000, &
           0.0000000,       0.0000000,       0.0000000,       0.0000000,       0.0000000, &
           0.0000000,       0.0000000,       0.0000000,       0.0000000,       0.0000000, &
           0.0000000,       0.0000000,       0.0000000,       0.0000000,       0.0000000, &
           0.0000000,       0.0000000,       0.0000000,       0.0000000,       0.0000000, &
           0.0000000,       0.0000000,       0.0000000,       0.0000000,       0.0000000, &
           0.0000000,   8.1754130e-09,    0.0069600246,     0.028010041,     0.063720063, &
           0.11360208,      0.15622409,      0.20035011,      0.24674112,      0.29440312, &
           0.34338113,      0.39289115,      0.44374018,      0.49459020,      0.54630418, &
           0.58104151,      0.61581843,      0.65063492,      0.68589990,      0.72116594, &
           0.74937819,      0.77063753,      0.79194696,      0.81330397,      0.83466097, &
           0.85601798,      0.87742898,      0.89890800,      0.92038701,      0.94186501, &
           0.96340602,      0.98495195,       1.0000000 /

      data a137 / &
           1.000000, 2.000365, 3.102241, 4.666084, 6.827977, 9.746966, 13.605424, 18.608931, 24.985718, 32.985710,  &
           42.879242, 54.955463, 69.520576, 86.895882, 107.415741, 131.425507, 159.279404, 191.338562, 227.968948, 269.539581,  &
           316.420746, 368.982361, 427.592499, 492.616028, 564.413452, 643.339905, 729.744141, 823.967834, 926.344910, 1037.20117,  &
           1156.853638, 1285.610352, 1423.770142, 1571.622925, 1729.448975, 1897.519287, 2076.095947, 2265.431641, 2465.770508, 2677.348145,  &
           2900.391357, 3135.119385, 3381.743652, 3640.468262, 3911.490479, 4194.930664, 4490.817383, 4799.149414, 5119.895020, 5452.990723,  &
           5798.344727, 6156.074219, 6526.946777, 6911.870605, 7311.869141, 7727.412109, 8159.354004, 8608.525391, 9076.400391, 9562.682617,  &
           10065.978516, 10584.631836, 11116.662109, 11660.067383, 12211.547852, 12766.873047, 13324.668945, 13881.331055, 14432.139648, 14975.615234,  &
           15508.256836, 16026.115234, 16527.322266, 17008.789062, 17467.613281, 17901.621094, 18308.433594, 18685.718750, 19031.289062, 19343.511719,  &
           19620.042969, 19859.390625, 20059.931641, 20219.664062, 20337.863281, 20412.308594, 20442.078125, 20425.718750, 20361.816406, 20249.511719,  &
           20087.085938, 19874.025391, 19608.572266, 19290.226562, 18917.460938, 18489.707031, 18006.925781, 17471.839844, 16888.687500, 16262.046875,  &
           15596.695312, 14898.453125, 14173.324219, 13427.769531, 12668.257812, 11901.339844, 11133.304688, 10370.175781, 9617.515625, 8880.453125,  &
           8163.375000, 7470.343750, 6804.421875, 6168.531250, 5564.382812, 4993.796875, 4457.375000, 3955.960938, 3489.234375, 3057.265625,  &
           2659.140625, 2294.242188, 1961.500000, 1659.476562, 1387.546875, 1143.250000, 926.507812, 734.992188, 568.062500, 424.414062,  &
           302.476562, 202.484375, 122.101562, 62.781250, 22.835938, 3.757813, 0.000000, 0.000000/

      data b137 / &
           0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,  &
           0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,  &
           0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,  &
           0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,  &
           0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,  &
           0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000007, 0.000024, 0.000059, 0.000112, 0.000199,  &
           0.000340, 0.000562, 0.000890, 0.001353, 0.001992, 0.002857, 0.003971, 0.005378, 0.007133, 0.009261,  &
           0.011806, 0.014816, 0.018318, 0.022355, 0.026964, 0.032176, 0.038026, 0.044548, 0.051773, 0.059728,  &
           0.068448, 0.077958, 0.088286, 0.099462, 0.111505, 0.124448, 0.138313, 0.153125, 0.168910, 0.185689,  &
           0.203491, 0.222333, 0.242244, 0.263242, 0.285354, 0.308598, 0.332939, 0.358254, 0.384363, 0.411125,  &
           0.438391, 0.466003, 0.493800, 0.521619, 0.549301, 0.576692, 0.603648, 0.630036, 0.655736, 0.680643,  &
           0.704669, 0.727739, 0.749797, 0.770798, 0.790717, 0.809536, 0.827256, 0.843881, 0.859432, 0.873929,  &
           0.887408, 0.899900, 0.911448, 0.922096, 0.931881, 0.940860, 0.949064, 0.956550, 0.963352, 0.969513,  &
           0.975078, 0.980072, 0.984542, 0.988500, 0.991984, 0.995003, 0.997630, 1.000000/

      SELECT CASE(km)

      CASE(20)

         do k=1,km+1
            ak(k) = a20_0178(k)
            bk(k) = b20_0178(k)
         enddo
         ! Search KS
         ks = 0
         do k=1,km
            if(bk(k) > 0) then
               ks = k-1
               goto 120
            endif
         enddo
120      continue

      CASE(26)

         do k=1,km+1
            ak(k) = a26(k)
            bk(k) = b26(k)
         enddo
         ! Search KS
         ks = 0
         do k=1,km
            if(bk(k) > 0) then
               ks = k-1
               goto 126
            endif
         enddo
126      continue

      CASE(40)
         !--------------------------------------------------
         ! Pure sigma-coordinate with uniform spacing in "z"
         !--------------------------------------------------
         ptop = 27381.905404907        ! model top pressure (pascal)
         press(1) = ptop
         press(km+1) = p0
         dlnp = (log(p0) - log(ptop)) / real(km)

         lnpe = log(press(km+1))
         do k=km,2,-1
            lnpe = lnpe - dlnp
            press(k) = exp(lnpe)
         enddo

         ! Search KS
         ks = 0
         do k=1,km
            if(press(k) >= pc) then
               ks = k-1
               goto 140
            endif
         enddo
140      continue

         if(ks /= 0) then
            do k=1,ks
               ak(k) = press(k)
               bk(k) = 0.
            enddo
         endif

         pint = press(ks+1)
         do k=ks+1,km
            ak(k) =  pint*(press(km)-press(k))/(press(km)-pint)
            bk(k) = (press(k) - ak(k)) / press(km+1)
         enddo
         ak(km+1) = 0.
         bk(km+1) = 1.

      CASE(60)
         !--------------------------------------------------
         ! Pure sigma-coordinate with uniform spacing in "z"
         !--------------------------------------------------
         ptop = 25499.234876157        ! model top pressure (pascal)
         press(1) = ptop
         press(km+1) = p0
         dlnp = (log(p0) - log(ptop)) / real(km)

         lnpe = log(press(km+1))
         do k=km,2,-1
            lnpe = lnpe - dlnp
            press(k) = exp(lnpe)
         enddo

         ! Search KS
         ks = 0
         do k=1,km
            if(press(k) >= pc) then
               ks = k-1
               goto 160
            endif
         enddo
160      continue

         if(ks /= 0) then
            do k=1,ks
               ak(k) = press(k)
               bk(k) = 0.
            enddo
         endif

         pint = press(ks+1)
         do k=ks+1,km
            ak(k) =  pint*(press(km)-press(k))/(press(km)-pint)
            bk(k) = (press(k) - ak(k)) / press(km+1)
         enddo
         ak(km+1) = 0.
         bk(km+1) = 1.

      CASE(72)

         do k=1,km+1
            ak(k) = a72(k)
            bk(k) = b72(k)
         enddo
         ! Search KS
         ks = 0
         do k=1,km
            if(bk(k) > 0) then
               ks = k-1
               goto 172
            endif
         enddo
172      continue

      CASE(137)

         do k=1,km+1
            ak(k) = a137(k)
            bk(k) = b137(k)
         enddo
         ! Search KS
         ks = 0
         do k=1,km
            if(bk(k) > 0) then
               ks = k-1
               goto 137
            endif
         enddo
137      continue

      CASE DEFAULT

         print*, 'Bad KM in FVdycoreCubed_GridComp:set_eta', km

      END SELECT

   end subroutine set_eta
#endif

   subroutine addTracer_r8(state, bundle, var, grid, fieldname)
      type(DynState), pointer :: state
      type(ESMF_FieldBundle) :: bundle
      real(r8), pointer :: var(:, :, :)
      type(ESMF_Grid) :: grid
      type(ESMF_DistGrid) :: dist_grid
      character(len=ESMF_MAXSTR) :: fieldname

      integer :: nq,rc,status
      type(DynTracers), pointer :: t(:)
      type(ESMF_Field) :: field
      real(r8), pointer :: ptr(:, :, :)

      call ESMF_GridGet(grid, distGrid=dist_grid, _RC)
      call ESMF_FieldBundleGet(bundle, fieldCount=nq, _RC)

      nq = nq + 1

      field = ESMF_FieldCreate(grid, var, datacopyflag=ESMF_DATACOPY_VALUE, name=fieldname, _RC)
      ! pchakrab - TODO: Need to replace MAPL_VLocationCenter with VERTICAL_STAGGER_CENTER
      ! pchakrab - TODO: need a replacement for MAPL_DimsHorzVert
      ! call ESMF_AttributeSet(field, name='VLOCATION', value=MAPL_VLocationCenter, _RC)
      ! call ESMF_AttributeSet(field, name='DIMS', value=MAPL_DimsHorzVert, rc=status)
      call MAPL_FieldBundleAdd(bundle, field, _RC)

      if (nq == 1) then
         allocate(state%VARS%tracer(nq), _STAT)
         call ESMF_FieldGet(field, localDE=0, farrayptr=ptr, _RC)
         state%vars%tracer(nq)%content => ptr
         state%vars%tracer(nq)%is_r4 = .false.
      else
         allocate(t(nq))
         t(1:nq-1) = state%vars%tracer
         deallocate(state%vars%tracer)
         state%vars%tracer => t
         call ESMF_FieldGet(field, localDE=0, farrayptr=ptr, _RC)
         state%vars%tracer(nq)%content => ptr
         state%vars%tracer(nq  )%is_r4 = .false.
      endif

      state%grid%nq = nq

      _RETURN(_SUCCESS)
   end subroutine addTracer_r8

   subroutine addTracer_r4(state, bundle, var, grid, fieldname)
      type(DynState), pointer :: state
      type(ESMF_FieldBundle) :: bundle
      real(r4), pointer :: var(:, :, :)
      type(ESMF_Grid) :: grid
      type(ESMF_DistGrid) :: dist_grid
      character(len=ESMF_MAXSTR) :: fieldname

      integer :: nq,rc,status
      type(DynTracers), pointer :: t(:)
      type(ESMF_Field) :: field
      real(r4), pointer :: ptr(:, :, :)

      call ESMF_GridGet(grid, distGrid=dist_grid, _RC)
      call ESMF_FieldBundleGet(bundle, fieldCount=NQ, _RC)

      nq = nq + 1

      field = ESMF_FieldCreate(grid, var, datacopyflag=ESMF_DATACOPY_VALUE, name=fieldname, _RC)
      ! pchakrab - TODO: Need to replace MAPL_VLocationCenter with VERTICAL_STAGGER_CENTER
      ! pchakrab - TODO: need a replacement for MAPL_DimsHorzVert
      ! call ESMF_AttributeSet(field, name='VLOCATION', value=MAPL_VLocationCenter, _RC)
      ! call ESMF_AttributeSet(field, name='DIMS', value=MAPL_DimsHorzVert, _RC)
      call MAPL_FieldBundleAdd(bundle, field, _RC)

      if (NQ == 1) then
         allocate(state%VARS%tracer(nq), _STAT)
         call ESMF_FieldGet(field, localDE=0, farrayptr=ptr, _RC)
         state%vars%tracer(nq)%content_r4 => ptr
         state%vars%tracer(nq)%is_r4 = .true.
      else
         allocate(t(nq))
         t(1:nq-1) = state%vars%tracer
         deallocate(state%vars%tracer)
         state%vars%tracer => t
         call ESMF_FieldGet(field, localDE=0, farrayptr=ptr, _RC)
         state%vars%tracer(nq)%content_r4 => ptr
         state%vars%tracer(nq  )%is_r4 = .true.
      endif

      state%grid%nq = nq

      _RETURN(_SUCCESS)
   end subroutine addTracer_r4

   subroutine freeTracers(state)
      type (DynState) :: STATE

      if (associated(STATE%VARS%tracer)) then
         DEALLOCATE( STATE%VARS%tracer)   ! Comment out to output tracer to checkpoint file
         NULLIFY( STATE%VARS%tracer)
      end if

      return
   end subroutine freeTracers

   Subroutine Write_Profile_2d_R8(grid, arr, name)
      type (DynGrid),   intent(IN) :: grid
      real(r8),         intent(IN) :: arr(grid%is:grid%ie,grid%js:grid%je)
      character(len=*), intent(IN) :: name

      integer  :: istrt,iend, jstrt,jend
      integer  :: im, jm
      real(r8) :: arr_global(grid%npx,grid%ntiles*grid%npy)
      real(r8) :: rng(3)
      real(r8) :: GSUM

      real(kind=ESMF_KIND_R8)     :: locArr(grid%is:grid%ie,grid%js:grid%je)
      real(kind=ESMF_KIND_R8)     :: glbArr(grid%npx,grid%ntiles*grid%npy)

      istrt = grid%is
      iend  = grid%ie
      jstrt = grid%js
      jend  = grid%je
      im    = grid%npx
      jm    = grid%npy*grid%ntiles

      !call write_parallel('GlobalSUm')
      locArr(:,:) = arr(:,:)
      call ArrayGather(locArr, glbArr, grid%grid)
      arr_global(:,:) = glbArr

      IF (MAPL_AM_I_ROOT()) Then
         rng(1) = MINVAL(MINVAL(arr_global,DIM=1),DIM=1)
         rng(2) = MAXVAL(MAXVAL(arr_global,DIM=1),DIM=1)
         rng(3) = SUM(SUM(arr_global,DIM=1),DIM=1)/(IM*JM)
         GSUM     = SUM(SUM(arr_global,DIM=1),DIM=1)

         print*,'***********'
         print*,'stats for ',trim(name)

         Write(*,'(3(f21.9,1x))')rng(:)
         !   Write(*,"('GlobalSum: ',f21.9)") GSUM
         print*,'***********'
         print*,' '
      End IF

   End Subroutine Write_Profile_2d_R8

   Subroutine Write_Profile_2d_R4(grid, arr, name)
      type (DynGrid),   intent(IN) :: grid
      real(r4),         intent(IN) :: arr(grid%is:grid%ie,grid%js:grid%je)
      character(len=*), intent(IN) :: name

      integer  :: istrt,iend, jstrt,jend
      integer  :: im, jm
      real(r4) :: arr_global(grid%npx,grid%ntiles*grid%npy)
      real(r4) :: rng(3)
      real(r4) :: GSUM

      real(kind=ESMF_KIND_R4)     :: locArr(grid%is:grid%ie,grid%js:grid%je)
      real(kind=ESMF_KIND_R4)     :: glbArr(grid%npx,grid%ntiles*grid%npy)

      istrt = grid%is
      iend  = grid%ie
      jstrt = grid%js
      jend  = grid%je
      im    = grid%npx
      jm    = grid%npy*grid%ntiles

      ! call write_parallel('GlobalSUm')
      locArr(:,:) = arr(:,:)
      call ArrayGather(locArr, glbArr, grid%grid)
      arr_global(:,:) = glbArr

      IF (MAPL_AM_I_ROOT()) Then
         rng(1) = MINVAL(MINVAL(arr_global,DIM=1),DIM=1)
         rng(2) = MAXVAL(MAXVAL(arr_global,DIM=1),DIM=1)
         rng(3) = SUM(SUM(arr_global,DIM=1),DIM=1)/(IM*JM)
         GSUM     = SUM(SUM(arr_global,DIM=1),DIM=1)

         print*,'***********'
         print*,'stats for ',trim(name)

         Write(*,'(3(f21.9,1x))')rng(:)
         !    Write(*,"('GlobalSum: ',f21.9)") GSUM
         print*,'***********'
         print*,' '
      End IF

   End Subroutine Write_Profile_2d_R4

   Subroutine Write_Profile_R8(grid, arr, name)
      type (DynGrid),   intent(IN) :: grid
      real(r8),         intent(IN) :: arr(grid%is:grid%ie,grid%js:grid%je,1:grid%npz)
      character(len=*), intent(IN) :: name

      integer  :: istrt,iend, jstrt,jend, kstrt,kend
      integer  :: im, jm, km, k
      real(r8), allocatable :: arr_global(:,:,:)
      real(r8) :: rng(3,grid%npz)
      real(r8) :: GSUM
      logical :: amIRoot

      real(kind=ESMF_KIND_R8)     :: locArr(grid%is:grid%ie,grid%js:grid%je)
      istrt = grid%is
      iend  = grid%ie
      jstrt = grid%js
      jend  = grid%je
      kstrt = 1
      kend  = grid%npz
      im    = grid%npx
      jm    = grid%npy*grid%ntiles
      km    = grid%npz

      amIRoot = MAPL_AM_I_ROOT()
      if (amIRoot) then
         allocate(arr_global(grid%npx,grid%ntiles*grid%npy,km))
      else
         allocate(arr_global(1,1,km))
      end if

      ! call write_parallel('GlobalSUm')
      do k=kstrt,kend
         locArr(:,:) = arr(:,:,k)
         call ArrayGather(locArr, arr_global(:,:,k), grid%grid)
      enddo

      IF (amIRoot) Then
         rng(1,:) = MINVAL(MINVAL(arr_global,DIM=1),DIM=1)
         rng(2,:) = MAXVAL(MAXVAL(arr_global,DIM=1),DIM=1)
         rng(3,:) = SUM(SUM(arr_global,DIM=1),DIM=1)/(IM*JM)
         GSUM     = SUM(SUM(SUM(arr_global,DIM=1),DIM=1),DIM=1)

         print*,'***********'
         print*,'stats for ',trim(name)

         Do k = 1, km
            Write(*,'(a,i4.0,3(f21.9,1x))')'k:',k,rng(:,k)
         End Do
         !    Write(*,"('GlobalSum: ',f21.9)") GSUM
         print*,'***********'
         print*,' '
      End IF

      deallocate(arr_global)

   End Subroutine Write_Profile_R8

   Subroutine Write_Profile_R4(grid, arr, name, delp)
      type (DynGrid),   intent(IN) :: grid
      real(r4),         intent(IN) :: arr(grid%is:grid%ie,grid%js:grid%je,1:grid%npz)
      character(len=*), intent(IN) :: name
      real(r8), optional, intent(IN) :: delp(grid%is:grid%ie,grid%js:grid%je,1:grid%npz)

      integer  :: istrt,iend, jstrt,jend, kstrt,kend
      integer  :: im, jm, km, k
      real(r4), allocatable :: arr_global(:,:,:)
      real(r4) :: rng(3,grid%npz)
      real(r8) :: gsum_p
      real(r4) :: GSUM
      logical :: amIRoot

      real(kind=ESMF_KIND_R8)     :: locArr(grid%is:grid%ie,grid%js:grid%je)
      real(kind=ESMF_KIND_R8), allocatable :: glbArr(:,:)

      istrt = grid%is
      iend  = grid%ie
      jstrt = grid%js
      jend  = grid%je
      kstrt = 1
      kend  = grid%npz
      im    = grid%npx
      jm    = grid%npy*grid%ntiles
      km    = grid%npz

      amIRoot = MAPL_AM_I_ROOT()
      if (amIRoot) then
         allocate(arr_global(grid%npx,grid%ntiles*grid%npy,km))
         allocate(glbArr(grid%npx,grid%ntiles*grid%npy))
      else
         allocate(arr_global(1,1,km))
         allocate(glbArr(1,1))
      end if

      do k=kstrt,kend
         locArr(:,:) = arr(:,:,k)
         call ArrayGather(locArr, glbArr, grid%grid)
         if (amIRoot) then
            arr_global(:,:,k) = glbArr
         end if
      enddo
      IF (amIRoot) Then
         rng(1,:) = MINVAL(MINVAL(arr_global,DIM=1),DIM=1)
         rng(2,:) = MAXVAL(MAXVAL(arr_global,DIM=1),DIM=1)
         rng(3,:) = SUM(SUM(arr_global,DIM=1),DIM=1)/(IM*JM)
         print*,'***********'
         print*,'stats for ',trim(name)
         Do k = 1, km
            Write(*,'(a,i4.0,3(f21.9,1x))')'k:',k,rng(:,k)
         End Do
         print*,'***********'
         print*,' '
      End IF

      if (present(delp)) then
         gsum_p = 0
         do k=kstrt,kend
            locArr(:,:) = arr(:,:,k)*grid%area(:,:)*delp(:,:,k)
            call ArrayGather(locArr, glbArr, grid%grid)
            if (amIRoot) then
               arr_global(:,:,k) = glbArr
            end if
            locArr(:,:) = delp(:,:,k)
            call ArrayGather(locArr, glbArr, grid%grid)
            if (amIRoot) then
               gsum_p = gsum_p + SUM(SUM(glbArr,DIM=1),DIM=1)
            end if
         enddo
         IF (amIRoot) Then
            GSUM     = SUM(SUM(SUM(arr_global,DIM=1),DIM=1),DIM=1)
            print*,'***********'
            Write(*,"('GlobalSum: ',e21.9)") GSUM/(grid%globalarea*gsum_p)
            print*,'***********'
            print*,' '
         End IF
      endif

      deallocate(arr_global, glbArr)

   End Subroutine Write_Profile_R4

   function R8_TO_R4(dbl_var)
      real(REAL8), intent(IN) :: dbl_var(:,:)
      real(REAL4) :: R8_TO_R4(LBOUND(dbl_var,1):UBOUND(dbl_var,1), LBOUND(dbl_var,2):UBOUND(dbl_var,2))
      integer :: i, j

      real(REAL8), parameter :: eps = 1.e-15_REAL8
      real(REAL8), parameter :: big = 1.e15_REAL8

      do j=LBOUND(dbl_var,2),UBOUND(dbl_var,2)
         do i=LBOUND(dbl_var,1),UBOUND(dbl_var,1)
            R8_TO_R4(i,j) = SIGN(MIN(big,MAX(eps,ABS(dbl_var(i,j)))),dbl_var(i,j))
         enddo
      enddo
   end function R8_TO_R4

   function R4_TO_R8(sngl_var)
      real(REAL4), intent(IN) :: sngl_var(:,:)
      real(REAL8)  :: R4_TO_R8(LBOUND(sngl_var,1):UBOUND(sngl_var,1), LBOUND(sngl_var,2):UBOUND(sngl_var,2))
      integer :: i, j

      real(REAL4), parameter :: eps = 1.e-15_REAL4
      real(REAL4), parameter :: big = 1.e15_REAL4

      do j=LBOUND(sngl_var,2),UBOUND(sngl_var,2)
         do i=LBOUND(sngl_var,1),UBOUND(sngl_var,1)
            R4_TO_R8(i,j) = SIGN(MIN(big,MAX(eps,ABS(sngl_var(i,j)))),sngl_var(i,j))
         enddo
      enddo
   end function R4_TO_R8

end module FVdycoreCubed_GridComp

subroutine SetServices(gc, rc)
   use ESMF
   use FVdycoreCubed_GridComp, only : mySetservices=>SetServices
   type(ESMF_GridComp) :: gc
   integer, intent(out) :: rc
   call mySetServices(gc, rc=rc)
end subroutine SetServices
