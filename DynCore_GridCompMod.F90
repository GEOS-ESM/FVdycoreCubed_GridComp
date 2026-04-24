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
   use mapl_ErrorHandlingMod, only: MAPL_Verify, MAPL_Assert, MAPL_Return

   use MAPL_Constants, only: MAPL_RADIUS, MAPL_CP, MAPL_PI, MAPL_PI_R8, MAPL_OMEGA, MAPL_KAPPA
   use MAPL_Constants, only: MAPL_P00, MAPL_GRAV, MAPL_RGAS, MAPL_RVAP, MAPL_CPVAP, MAPL_O3MW, MAPL_AIRMW
   use MAPL_Constants, only: MAPL_VectorField ! pchakrab: TODO - need MAPL3 equivalent
   use MAPL_Constants, only: MAPL_UNDEFINED_REAL

   use FileIOSharedMod, only: WRITE_PARALLEL

   use mapl3g_Utilities, only: MAPL_MaxMin, MAPL_AreaMean
   use mapl3g_Utilities_Comms_API, only: MAPL_Am_I_Root, MAPL_ArrayGather
   use mapl3g_generic, only: MAPL_GridCompSetGeometry
   use mapl3g_generic, only: MAPL_GridCompGet, MAPL_GridCompGetResource
   use mapl3g_generic, only: MAPL_GridCompSetEntryPoint, MAPL_GridCompGetInternalState
   use mapl3g_generic, only: MAPL_GridCompAddSpec, MAPL_STATEITEM_SERVICE
   use mapl3g_generic, only: MAPL_UserCompSetInternalState, MAPL_UserCompGetInternalState
   use mapl3g_generic, only: MAPL_GridCompTimerStart, MAPL_GridCompTimerStop
   use mapl3g_VerticalStaggerLoc, only: VERTICAL_STAGGER_NONE, VERTICAL_STAGGER_CENTER, VERTICAL_STAGGER_EDGE
   use mapl3g_Geom_API, only: MAPL_GridGetCoordinates
   use mapl3g_State_API, only: MAPL_StateGetPointer
   use mapl3g_Field_API, only: MAPL_FieldCreate
   use mapl3g_FieldBundle_API, only: MAPL_FieldBundleAdd
   use mapl3g_RestartModes, only: MAPL_RESTART_SKIP, MAPL_RESTART_REQUIRED

   use pflogger, only: logger_t => logger
   use gftl2_StringVector, only: StringVector

   use m_set_eta, only: set_eta

   ! FV Specific Module
   use fv_arrays_mod, only: REAL4, REAL8, FVPRC
   !use fv_grid_tools_mod, only: grid_type
   use FV_StateMod, only : &
        FV_Atm, &
        FV_To_State, State_To_FV, DEBUG_FV_STATE, &
        DynTracers => T_TRACERS, &
        DynVars => T_FVDYCORE_VARS, &
        DynGrid => T_FVDYCORE_GRID, &
        DynState => T_FVDYCORE_STATE, &
        DynSetup => FV_Setup, &
        DynInit => FV_InitState, &
        DynRun => FV_Run, &
        DynFinalize => FV_Finalize, &
        getAllWinds => fv_getAllWinds, &
        getVorticity => fv_getVorticity, &
        getDivergence => fv_getDivergence, &
        fillMassFluxes => fv_fillMassFluxes, &
        computeMassFluxes => fv_computeMassFluxes,&
        getVerticalMassFlux => fv_getVerticalMassFlux,&
        getOmega => fv_getOmega, &
        getEPV => fv_getEPV, &
        getPKZ => fv_getPKZ, &
        getDELZ => fv_getDELZ, &
        getQ => fv_getQ, &
        Agrid_To_Native => INTERP_AGRID_TO_DGRID, &
        DYN_COLDSTART => COLDSTART, &
        DYN_CASE => CASE_ID, &
        DYN_DEBUG => DEBUG, &
        HYDROSTATIC => FV_HYDROSTATIC, &
        fv_getUpdraftHelicity, DEBUG_DYN, DEBUG_ADV, &
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

   integer, parameter :: r8 = REAL8
   integer, parameter :: r4 = REAL4

   real(kind=r4), parameter :: RADIUS = MAPL_RADIUS
   real(kind=r4), parameter :: CP = MAPL_CP
   real(kind=r4), parameter :: PI = MAPL_PI_R8
   real(kind=r4), parameter :: OMEGA = MAPL_OMEGA
   real(kind=r4), parameter :: KAPPA = MAPL_KAPPA
   real(kind=r4), parameter :: P00 = MAPL_P00
   real(kind=r4), parameter :: GRAV = MAPL_GRAV
   real(kind=r4), parameter :: RGAS = MAPL_RGAS
   real(kind=r4), parameter :: RVAP = MAPL_RVAP
   real(kind=r4), parameter :: EPS = RVAP / RGAS - 1.0

   integer, parameter :: TIME_TO_RUN = 1
   integer, parameter :: CHECK_MAXMIN = 2

   ! Tracer I/O History stuff
   integer, parameter :: ntracers = 11
   integer, parameter :: plevs(5) = [850, 700, 600, 500, 300]

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

   character(*), parameter :: PRIVATE_STATE = "DYN_STATE"

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

   subroutine SetServices(gc, rc)

      !ARGUMENTS:
      type(ESMF_GridComp) :: gc
      integer, intent(out) :: rc

      !DESCRIPTION: Set services (register) for the FVCAM Dynamical Core GridComp
      !EOP

#ifdef SKIP_TRACERS
      character(len=ESMF_MAXSTR) :: myTracer
      integer :: ilev, itracer
#endif
      logical :: FV3_STANDALONE
      integer :: status

      ! Wrap gridcomp's private state and store it in gc
      _SET_NAMED_PRIVATE_STATE(gc, DynState, PRIVATE_STATE)

#include "DynCore_Import___.h"
#include "DynCore_Export___.h"
#include "DynCore_Internal___.h"

      ! pchakrab - DynCore is the provider here, so the service items are not needed
      ! NOTE: SERVICE, irrespective of whether you are a provider or subscriber, adds the bundle
      ! to BOTH the export and import states
      call MAPL_GridCompAddSpec(gc, &
           state_intent=ESMF_STATEINTENT_EXPORT, &
           short_name='TRADV', &
           standard_name='advected_quantities', &
      ! pchakrab: TODO - we shouldn't need dims and vstagger for a bundle
           dims="xyz", &
           vstagger=VERTICAL_STAGGER_NONE, &
           units='unknown', &
           itemtype=MAPL_STATEITEM_SERVICE, _RC)

#ifdef SKIP_TRACERS
      do itracer = 1, ntracers
         do ilev = 1, size(plevs)
            write(myTracer, "('Q',i5.5,'_',i3.3)") itracer - 1, plevs(ilev)
            call MAPL_AddExportSpec(gc, &
                 short_name=TRIM(myTracer), &
                 long_name=TRIM(myTracer), &
                 units='1', &
                 dims=MAPL_DimsHorzOnly, &
                 vlocation=MAPL_VLocationNone, _RC)
         end do
         write(myTracer, "('Q',i5.5)") itracer - 1
         call MAPL_AddExportSpec(gc, &
              short_name=TRIM(myTracer), &
              long_name=TRIM(myTracer), &
              units='1', &
              dims=MAPL_DimsHorzVert, &
              vlocation=MAPL_VLocationCenter, _RC)
      end do
#endif

      ! pchakrab: TODO: DO WE STILL NEED THIS COMMENT?
      !ALT: technically the first 2 records of "old" style FV restart have
      !     6 ints: YYYY MM DD H M S
      !     5 ints: I,J,K, KS (num true pressure levels), NQ (num tracers) headers

      ! Register services for this component
      call MAPL_GridCompSetEntryPoint(gc, ESMF_Method_Initialize, Initialize, _RC)
      call MAPL_GridCompSetEntryPoint(gc, ESMF_Method_Run, Run, phase_name="Run", _RC)
      call MAPL_GridCompSetEntryPoint(gc, ESMF_Method_Run, RunAddIncs, phase_name="RunAddIncs", _RC)
      call MAPL_GridCompSetEntryPoint(gc, ESMF_Method_Finalize, Finalize, _RC)
      !  call MAPL_GridCompSetEntryPoint(gc, ESMF_SETREADRESTART, Coldstart, _RC)

      ! Setup geometry
      call MAPL_GridCompSetGeometry(gc, _RC)

      ! Register prototype of cubed sphere grid and associated regridders
      call register_grid_and_regridders()

      ! At this point check if FV is standalone and init the grid
      call MAPL_GridCompGetResource(gc, "FV3_STANDALONE", FV3_STANDALONE, default=.false., _RC)
      if (FV3_STANDALONE) then
         call MAPL_GridCompAddSpec(gc, &
              state_intent=ESMF_STATEINTENT_EXPORT, &
              short_name='TRADVEX', &
              standard_name='advected_quantities', &
         ! pchakrab - TODO: we shouldn't need dims and vstagger for a bundle
              dims="xyz", &
              vstagger=VERTICAL_STAGGER_NONE, &
              units='unknown', &
              itemtype=MAPL_STATEITEM_SERVICE, _RC)
      end if

      call MAPL_GridCompGetResource(gc, "DEBUG_DYN", DEBUG_DYN, default=.false., _RC)
      call MAPL_GridCompGetResource(gc, "DEBUG_ADV", DEBUG_ADV, default=.false., _RC)
      call MAPL_GridCompGetResource(gc, "DEBUG_TQ_ERRORS", DEBUG_TQ_ERRORS, default=.false., _RC)

      _RETURN(_SUCCESS)

   end subroutine SetServices

   subroutine Initialize(gc, import, export, clock, rc)
      !ARGUMENTS:
      type(ESMF_GridComp) :: gc ! composite gridded component
      type(ESMF_State) :: import ! import state
      type(ESMF_State) :: export ! export state
      type(ESMF_Clock) :: clock ! the clock
      integer, intent(out) :: rc ! Error code, 0 all is well, error otherwise

      type(DynState), pointer :: self

      type(ESMF_Field) :: field
      type(ESMF_State) :: internal
      type(ESMF_FieldBundle) :: tradv, tradvex

      real(kind=r4), pointer :: pref(:)
      real(kind=r4), pointer :: u(:, :, :), v(:, :, :), t(:, :, :)
      real(kind=r4), pointer :: temp2d(:, :)

      real(kind=r8), pointer :: ak(:), bk(:)
      real(kind=r8), pointer :: ud(:, :, :), vd(:, :, :)
      real(kind=r8), pointer :: pt(:, :, :), pk(:, :, :)
      real(kind=r8), allocatable :: ur(:, :, :), vr(:, :, :) ! rotated winds

      logical :: ColdRestart, FV3_STANDALONE
      integer :: ifirst, ilast, jfirst, jlast, km
      integer :: i, numTracers, status

      ! Setup FMS/FV3
      call MAPL_GridCompTimerStart(gc, "DYN_SETUP", _RC)
      call DynSetup(gc, _RC)
      call MAPL_GridCompTimerStop(gc, "DYN_SETUP", _RC)

      ! Get the private state
      _GET_NAMED_PRIVATE_STATE(gc, DynState, PRIVATE_STATE, self)

      call MAPL_GridCompGetResource(gc, "DO_ADD_INCS", DO_ADD_INCS, default=DO_ADD_INCS, _RC)

      ! Check for ColdStart from the configuration
      call MAPL_GridCompGetResource(gc, "COLDSTART", ColdRestart, default=.false., _RC)
      if (ColdRestart) then
         call COLDSTART(gc, import, export, clock, _RC)
      end if

      ! Set Private Internal State from Restart File
      call MAPL_GridCompGetInternalState(gc, internal, _RC)

      call MAPL_GridCompTimerStart(gc, "DYN_INIT", _RC)
      call DynInit(self, clock, internal, import, gc, _RC)
      call MAPL_GridCompTimerStop(gc, "DYN_INIT", _RC)

      ! Create PREF EXPORT Coupling (Needs to be done only once per run)
      call MAPL_StateGetPointer(internal, ak, "AK", _RC)
      call MAPL_StateGetPointer(internal, bk, "BK", _RC)
      call MAPL_StateGetPointer(export, pref, "PREF", _RC)
      if (associated(pref)) pref = ak + bk * P00

      ! Create A-Grid Winds
      call MAPL_StateGetPointer(export, u, "U", _RC)
      call MAPL_StateGetPointer(export, v, "V", _RC)
      if (associated(u) .and. associated(v)) then
         ifirst = self%grid%is
         ilast = self%grid%ie
         jfirst = self%grid%js
         jlast = self%grid%je
         km = self%grid%npz
         allocate(ur(ifirst:ilast, jfirst:jlast, km))
         allocate(vr(ifirst:ilast, jfirst:jlast, km))
         call MAPL_StateGetPointer(internal, ud, "U", _RC)
         call MAPL_StateGetPointer(internal, vd, "V", _RC)
         call getAllWinds(ud, vd, ur=ur, vr=vr)
         u = ur
         v = vr
         deallocate(ur, vr)
      end if

      ! Temperature
      call MAPL_StateGetPointer(export, t, "T", _RC)
      if (associated(t)) then
         call MAPL_StateGetPointer(internal, pt, "PT", _RC)
         call MAPL_StateGetPointer(internal, pk, "PKZ", _RC)
         t = pt * pk
      end if

      ! Fill Grid-Cell Area Delta-X/Y
      call MAPL_StateGetPointer(export, temp2d, "DXC", _RC)
      if (associated(temp2d)) temp2d = self%grid%dxc

      call MAPL_StateGetPointer(export, temp2d, "DYC", _RC)
      if (associated(temp2d)) temp2d = self%grid%dyc

      call MAPL_StateGetPointer(export, temp2d, "AREA", _RC)
      if (associated(temp2d)) temp2d = self%grid%area

      ! ======================================================================
      !ALT: the next section addresses the problem when export variables have been
      !     assigned values during Initialize. To prevent "connected" exports
      !     being overwritten by DEFAULT in the Import spec in the other component
      !     we label them as being "initailized by restart". A better solution
      !     would be to move the computation to phase 2 of Initialize and
      !     eliminate this section alltogether
      ! ======================================================================
      ! ! pchakrab: TODO - do we need to port the following to MAPL3
      ! call ESMF_StateGet(export, "PREF", field, _RC)
      ! call MAPL_AttributeSet(field, NAME="MAPL_InitStatus", VALUE=MAPL_InitialRestart, _RC)

      ! call ESMF_StateGet(export, "PLE", field, _RC)
      ! call MAPL_AttributeSet(field, NAME="MAPL_InitStatus", VALUE=MAPL_InitialRestart, _RC)

      ! call ESMF_StateGet(export, "U", field, _RC)
      ! call MAPL_AttributeSet(field, NAME="MAPL_InitStatus", VALUE=MAPL_InitialRestart, _RC)

      ! call ESMF_StateGet(export, "V", field, _RC)
      ! call MAPL_AttributeSet(field, NAME="MAPL_InitStatus", VALUE=MAPL_InitialRestart, _RC)

      ! call ESMF_StateGet(export, "T", field, _RC)
      ! call MAPL_AttributeSet(field, NAME="MAPL_InitStatus", VALUE=MAPL_InitialRestart, _RC)

      call MAPL_GridCompGetResource(gc, "FV3_STANDALONE", FV3_STANDALONE, default=.false., _RC)
      if (FV3_STANDALONE) then
         call ESMF_StateGet(import, "TRADV", tradv, _RC)
         call ESMF_StateGet(export, "TRADVEX", tradvex, _RC)
         call ESMF_FieldBundleGet(tradv, fieldCount=numTracers, _RC)
         do i = 1, numTracers
            call ESMF_FieldBundleGet(tradv, fieldIndex=i, field=field, _RC)
            call MAPL_FieldBundleAdd(tradvex, [field], _RC)
         end do
      end if

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

      integer :: status, comm
      type(ESMF_VM) :: vm
      type(ESMF_FieldBundle) :: bundle
      type(ESMF_Field) :: field
      type(ESMF_Alarm) :: alarm
      type(ESMF_Grid) :: esmfgrid
      type(ESMF_HConfig) :: hconfig

      type(DynState), pointer :: self
      type(DynGrid), pointer :: grid
      type(DynVars), pointer :: vars

      integer :: nq
      integer :: im, jm, km
      integer :: numvars
      integer :: ifirstxy, ilastxy, jfirstxy, jlastxy
      integer :: kend, i, j, k, n
      ! pchakrab - convt is not used anywhere
      logical, parameter :: convt = .false. ! Until this is run with full physics
      logical :: is_shutoff, is_ringing

      real(kind=r8), pointer :: phisxy(:, :)
      real(kind=4), pointer :: phis(:, :)

      real(kind=r8), allocatable :: plk(:, :, :) ! pl**kappa
      real(kind=r8), allocatable :: pkxy(:, :, :) ! pe**kappa
      real(kind=r8), allocatable :: pe0(:, :, :) ! edge-level pressure before dynamics
      real(kind=r8), allocatable :: pe1(:, :, :) ! edge-level pressure after dynamics
      real(kind=r8), allocatable :: pl(:, :, :) ! mid-level pressure
      real(kind=r8), allocatable :: tempxy(:, :, :) ! mid-level temperature

      real(kind=r8), allocatable :: ua(:, :, :) ! temporary array
      real(kind=r8), allocatable :: va(:, :, :)
      real(kind=r8), allocatable :: uc(:, :, :)
      real(kind=r8), allocatable :: vc(:, :, :)
      real(kind=r8), allocatable :: uc0(:, :, :)
      real(kind=r8), allocatable :: vc0(:, :, :)
      real(kind=r8), allocatable :: ur(:, :, :)
      real(kind=r8), allocatable :: vr(:, :, :)
      real(kind=r8), allocatable :: qv(:, :, :)
      real(kind=r8), allocatable :: ql(:, :, :)
      real(kind=r8), allocatable :: qi(:, :, :)
      real(kind=r8), allocatable :: qr(:, :, :)
      real(kind=r8), allocatable :: qs(:, :, :)
      real(kind=r8), allocatable :: qg(:, :, :)
      real(kind=r8), allocatable :: qdnew(:, :, :)
      real(kind=r8), allocatable :: qdold(:, :, :)
      real(kind=r8), allocatable :: qvold(:, :, :)
      real(kind=r8), allocatable :: qlold(:, :, :)
      real(kind=r8), allocatable :: qiold(:, :, :)
      real(kind=r8), allocatable :: qrold(:, :, :)
      real(kind=r8), allocatable :: qsold(:, :, :)
      real(kind=r8), allocatable :: qgold(:, :, :)
      real(kind=r8), allocatable :: delpold(:, :, :)
      real(kind=r8), allocatable :: ox(:, :, :)
      real(kind=r8), allocatable :: zl(:, :, :)
      real(kind=r8), allocatable :: zle(:, :, :)
      real(kind=r8), allocatable :: logpe(:, :, :)
      real(kind=r8), allocatable :: delp(:, :, :)
      real(kind=r8), allocatable :: dudt(:, :, :)
      real(kind=r8), allocatable :: dvdt(:, :, :)
      real(kind=r8), allocatable :: dtdt(:, :, :)
      real(kind=r8), allocatable :: dqdt(:, :, :)
      real(kind=r8), allocatable :: dthdt(:, :, :)
      real(kind=r8), allocatable :: ddpdt(:, :, :)
      real(kind=r8), allocatable :: dpedt(:, :, :)
      real(kind=FVPRC), allocatable :: tmp3d(:, :, :)
      real(kind=FVPRC), allocatable :: vort(:, :, :)
      real(kind=FVPRC), allocatable :: divg(:, :, :)
      real(kind=r8), allocatable :: dmdt(:, :)

      real(kind=r8), allocatable :: qsum1(:, :) ! Vertically Integrated Variable
      real(kind=r4), allocatable :: qsum2(:, :) ! Vertically Integrated Variable

      real(kind=r8), allocatable :: penrg(:, :) ! Vertically Integrated Cp*T
      real(kind=r8), allocatable :: kenrg(:, :) ! Vertically Integrated K
      real(kind=r8), allocatable :: tenrg(:, :) ! PHIS*(Psurf-Ptop)
      real(kind=r8), allocatable :: penrg0(:, :) ! Vertically Integrated Cp*T
      real(kind=r8), allocatable :: kenrg0(:, :) ! Vertically Integrated K
      real(kind=r8), allocatable :: tenrg0(:, :) ! PHIS*(Psurf-Ptop)
      real(kind=r8), allocatable :: kedyn(:, :)
      real(kind=r8), allocatable :: pedyn(:, :)
      real(kind=r8), allocatable :: tedyn(:, :)

      real(kind=r4), allocatable :: dqvdtanaint1(:, :)
      real(kind=r4), allocatable :: dqvdtanaint2(:, :)
      real(kind=r4), allocatable :: dqldtanaint1(:, :)
      real(kind=r4), allocatable :: dqldtanaint2(:, :)
      real(kind=r4), allocatable :: dqidtanaint1(:, :)
      real(kind=r4), allocatable :: dqidtanaint2(:, :)
      real(kind=r4), allocatable :: doxdtanaint1(:, :)
      real(kind=r4), allocatable :: doxdtanaint2(:, :)
      real(kind=r4), allocatable :: dthdtanaint1(:, :)
      real(kind=r4), allocatable :: dthdtanaint2(:, :)

      real(kind=r4), allocatable :: tropp1(:, :) ! Tropopause Pressure
      real(kind=r4), allocatable :: tropp2(:, :) ! Tropopause Pressure
      real(kind=r4), allocatable :: tropp3(:, :) ! Tropopause Pressure
      real(kind=r4), allocatable :: tropt(:, :) ! Tropopause Temperature
      real(kind=r4), allocatable :: tropq(:, :) ! Tropopause Specific Humidity

      real(kind=r8), allocatable :: omaxyz(:, :, :) ! vertical pressure velocity (pa/sec)
      real(kind=r8), allocatable :: epvxyz(:, :, :) ! ertel's potential vorticity

      real(kind=r8), allocatable :: cxxyz(:, :, :) ! Accumulated eastward courant numbers
      real(kind=r8), allocatable :: cyxyz(:, :, :) ! Accumulated northward courant numbers
      real(kind=r8), allocatable :: mfxxyz(:, :, :) ! Accumulated eastward mass flux
      real(kind=r8), allocatable :: mfyxyz(:, :, :) ! Accumulated northward mass flux
      real(kind=r8), allocatable :: mfzxyz(:, :, :) ! Accumulated vertical mass flux

      real(kind=FVPRC) :: dt ! Dynamics time step
      real(kind=r8), allocatable :: trsum1(:) ! Global Sum of Tracers before Add_Incs
      real(kind=r8), allocatable :: trsum2(:) ! Global Sum of Tracers after  Add_Incs

      real(kind=r4), pointer :: dudtana(:, :, :)
      real(kind=r4), pointer :: dvdtana(:, :, :)
      real(kind=r4), pointer :: dtdtana(:, :, :)
      real(kind=r4), pointer :: ddpdtana(:, :, :)
      real(kind=r4), pointer :: qctmp(:, :, :)
      real(kind=r4), pointer :: dqldt(:, :, :)
      real(kind=r4), pointer :: dqidt(:, :, :)
      real(kind=r4), pointer :: doxdt(:, :, :)
      real(kind=r4), pointer :: dqvana(:, :, :)
      real(kind=r4), pointer :: dqlana(:, :, :)
      real(kind=r4), pointer :: dqiana(:, :, :)
      real(kind=r4), pointer :: dqrana(:, :, :)
      real(kind=r4), pointer :: dqsana(:, :, :)
      real(kind=r4), pointer :: dqgana(:, :, :)
      real(kind=r4), pointer :: doxana(:, :, :)
      real(kind=r4), pointer :: temp3d(:, :, :)
      real(kind=r4), pointer :: vtmp3d(:, :, :)
      real(kind=r4), pointer :: area(:, :)
      real(kind=r4), pointer :: temp2d(:, :)
      real(kind=r4), pointer :: tempu(:, :)
      real(kind=r4), pointer :: tempv(:, :)

      real(kind=r4), pointer :: uh25(:, :), uh03(:, :)
      real(kind=r4), pointer :: srh01(:, :), srh03(:, :), srh25(:, :), shr06(:, :)
      real(kind=r4), allocatable :: uh25tmp(:, :), uh03tmp(:, :)
      real(kind=r4), allocatable :: srh01tmp(:, :), srh03tmp(:, :), srh25tmp(:, :), shr06tmp(:, :)

      real(kind=r8), allocatable :: uatmp(:, :, :)
      real(kind=r8), allocatable :: vatmp(:, :, :)
      real(kind=r8), allocatable :: udtmp(:, :, :)
      real(kind=r8), allocatable :: vdtmp(:, :, :)

      character(len=ESMF_MAXSTR), allocatable :: names(:), names0(:)

      real(kind=r4), allocatable :: lats(:, :), lons(:, :)
      ! real(r4), allocatable :: ZTH(:,:), SLR(:,:) - used in some commented code

      real :: HGT_SURFACE

      character(len=:), allocatable :: ANA_IS_WEIGHTED
      logical :: is_weighted

      type(DynTracers) :: qqq ! Specific Humidity
      type(DynTracers) :: ooo ! ox
      logical :: LCONSV, LFILL
      integer :: CONSV, FILL

      logical, save :: firstime = .true.
      logical :: adjust_tracers, exclude, do_energetics, do_tropvars
      logical :: FV3_STANDALONE
      integer :: pos, nqt
#ifdef SKIP_TRACERS
      integer :: itracer
#endif
      type(ESMF_Alarm) :: predictor_alarm
      type(ESMF_Grid) :: bgrid
      character(len=ESMF_MAXSTR) :: fieldname
      character(len=ESMF_MAXSTR), allocatable :: biggerlist(:)
      character(len=:), allocatable :: adjust_tracer_mode, xlist(:)
      real(kind=r8) :: t1, t2, dyn_run_timer
      class(logger_t), pointer :: logger

      call ESMF_VMGetCurrent(vm, _RC)
      call ESMF_VMGet(vm, mpiCommunicator=comm, _RC)

      call MAPL_GridCompGet(gc, grid=esmfgrid, hconfig=hconfig, logger=logger, _RC)
      call ESMF_GridValidate(esmfgrid, _RC)

      call MAPL_GridGetCoordinates(esmfgrid, longitudes=lons, latitudes=lats, _RC)
      call MAPL_StateGetPointer(export, temp2d, "LONS", _RC)
      if (associated(temp2d)) temp2d = lons
      call MAPL_StateGetPointer(export, temp2d, "LATS", _RC)
      if (associated(temp2d)) temp2d = lats

      ! Retrieve the pointer to the internal state
      _GET_NAMED_PRIVATE_STATE(gc, DynState, PRIVATE_STATE, self)

      vars => self%vars ! direct handle to control variables
      grid => self%grid ! direct handle to grid
      dt = self%dt ! dynamics time step (large)

      ifirstxy = grid%is
      ilastxy = grid%ie
      jfirstxy = grid%js
      jlastxy = grid%je
      im = grid%npx
      jm = grid%npy
      km = grid%npz

      is_ringing = ESMF_AlarmIsRinging(self%ALARMS(TIME_TO_RUN), _RC)
      if (.not.is_ringing) return

      ! Allocate Arrays
      allocate(delp(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
      allocate(dudt(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
      allocate(dvdt(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
      allocate(dtdt(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
      allocate(dqdt(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
      allocate(dthdt(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
      allocate(ddpdt(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
      allocate(dpedt(ifirstxy:ilastxy, jfirstxy:jlastxy, km + 1))
      allocate(tempxy(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
      allocate(pe0(ifirstxy:ilastxy, jfirstxy:jlastxy, km + 1))
      allocate(pe1(ifirstxy:ilastxy, jfirstxy:jlastxy, km + 1))
      allocate(pl(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
      allocate(ua(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
      allocate(va(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
      allocate(uc(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
      allocate(vc(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
      allocate(uc0(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
      allocate(vc0(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
      allocate(ur(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
      allocate(vr(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
      allocate(qv(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
      allocate(ql(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
      allocate(qi(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
      allocate(qr(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
      allocate(qs(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
      allocate(qg(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
      allocate(ox(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
      allocate(qsum1(ifirstxy:ilastxy, jfirstxy:jlastxy))
      allocate(qsum2(ifirstxy:ilastxy, jfirstxy:jlastxy))
      allocate(dmdt(ifirstxy:ilastxy, jfirstxy:jlastxy))

      do_energetics = .false.
      call MAPL_StateGetPointer(export, temp2d, "KEANA", _RC)
      if (associated(temp2d)) do_energetics = .true.
      call MAPL_StateGetPointer(export, temp2d, "PEANA", _RC)
      if (associated(temp2d)) do_energetics = .true.
      call MAPL_StateGetPointer(export, temp2d, "TEANA", _RC)
      if (associated(temp2d)) do_energetics = .true.
      call MAPL_StateGetPointer(export, temp2d, "KEDYN", _RC)
      if (associated(temp2d)) do_energetics = .true.
      call MAPL_StateGetPointer(export, temp2d, "PEDYN", _RC)
      if (associated(temp2d)) do_energetics = .true.
      call MAPL_StateGetPointer(export, temp2d, "TEDYN", _RC)
      if (associated(temp2d)) do_energetics = .true.
      if (do_energetics) then
         allocate(kedyn(ifirstxy:ilastxy, jfirstxy:jlastxy))
         allocate(pedyn(ifirstxy:ilastxy, jfirstxy:jlastxy))
         allocate(tedyn(ifirstxy:ilastxy, jfirstxy:jlastxy))
         allocate(kenrg(ifirstxy:ilastxy, jfirstxy:jlastxy))
         allocate(penrg(ifirstxy:ilastxy, jfirstxy:jlastxy))
         allocate(tenrg(ifirstxy:ilastxy, jfirstxy:jlastxy))
         allocate(kenrg0(ifirstxy:ilastxy, jfirstxy:jlastxy))
         allocate(penrg0(ifirstxy:ilastxy, jfirstxy:jlastxy))
         allocate(tenrg0(ifirstxy:ilastxy, jfirstxy:jlastxy))
      end if

      allocate(vort(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
      allocate(divg(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
      allocate(tmp3d(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
      allocate(phisxy(ifirstxy:ilastxy, jfirstxy:jlastxy))
      allocate(plk(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
      allocate(pkxy(ifirstxy:ilastxy, jfirstxy:jlastxy, km + 1))
      allocate(zl(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
      allocate(zle(ifirstxy:ilastxy, jfirstxy:jlastxy, km + 1))
      allocate(logpe(ifirstxy:ilastxy, jfirstxy:jlastxy, km + 1))
      allocate(omaxyz(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
      allocate(epvxyz(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
      allocate(cxxyz(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
      allocate(cyxyz(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
      allocate(mfxxyz(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
      allocate(mfyxyz(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
      allocate(mfzxyz(ifirstxy:ilastxy, jfirstxy:jlastxy, km + 1))

      ! Report advected friendlies
      ! call ESMF_StateGet(import, "TRADV", bundle, _RC)
      ! TODO: pchakrab - we don't want to modify import data, so at this point
      ! we should copy the tracer field data to the other bundle, from export
      call ESMF_StateGet(export, "TRADV", bundle, _RC)

      ! ALT: this section attempts to limit the amount of advected tracers
      adjust_tracers = .false.
      call MAPL_GridCompGetResource(gc, "EXCLUDE_ADVECTION_TRACERS", adjust_tracer_mode, default="ALWAYS", _RC)
      if (adjust_tracer_mode == "ALWAYS") then
         adjust_tracers = .true.
      else if (adjust_tracer_mode == "PREDICTOR") then
         !get predictor_alarm from clock
         call ESMF_ClockGetAlarm(clock, alarmName="PredictorAlarm", alarm=predictor_alarm, rc=status)
         if (status == ESMF_SUCCESS) then
            !check if ringing
            if (ESMF_AlarmIsRinging(predictor_alarm)) then
               adjust_tracers = .true.
            end if
         end if
      else
         call logger%info("Run:: Invalid value specified for EXCLUDE_ADVECTION_TRACERS, ignored")
         adjust_tracers = .false.
      end if
      if (adjust_tracers) then
         if (firstime) then
            firstime = .false.
            xlist = ESMF_HConfigAsStringSeq(hconfig, keyString="EXCLUDE_ADVECTION_TRACERS_LIST", stringLen=ESMF_MAXSTR, _RC)
            if (allocated(xlist)) n = size(xlist)

            ! Count the number of tracers
            call ESMF_FieldBundleGet(bundle, grid=bgrid, fieldCount=nqt, _RC)
            BundleAdv = ESMF_FieldBundleCreate(name="xTRADV", _RC)
            call ESMF_FieldBundleSet(BundleAdv, grid=bgrid, _RC)
            ! loop over NQ in TRADV
            do i = 1, nqt
               !get field from TRADV and its name
               call ESMF_FieldBundleGet(bundle, fieldIndex=i, field=field, _RC)
               call ESMF_FieldGet(field, name=fieldname, _RC)
               !exclude everything that is not cloud/water species
               if ((AdvCore_Advection >= 1) .and. &
                    (TRIM(fieldname) /= "Q") .and. &
                    (TRIM(fieldname) /= "QLCN") .and. &
                    (TRIM(fieldname) /= "QLLS") .and. &
                    (TRIM(fieldname) /= "QICN") .and. &
                    (TRIM(fieldname) /= "QILS") .and. &
                    (TRIM(fieldname) /= "CLCN") .and. &
                    (TRIM(fieldname) /= "CLLS") .and. &
                    (TRIM(fieldname) /= "NCPL") .and. &
                    (TRIM(fieldname) /= "NCPI") .and. &
                    (TRIM(fieldname) /= "QRAIN") .and. &
                    (TRIM(fieldname) /= "QSNOW") .and. &
                    (TRIM(fieldname) /= "QGRAUPEL")) then
                  call logger%info("Run:: FV3+ADV is excluding " // TRIM(fieldname))

                  n = n + 1
                  if (n > size(xlist)) then
                     allocate(biggerlist(2 * n), _STAT)
                     biggerlist(1:n - 1) = xlist
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
               if (.not.exclude) then
                  call MAPL_FieldBundleAdd(BundleAdv, [field], _RC)
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
      end if ! adjust_tracers

      call ESMF_FieldBundleGet(bundle, fieldCount=nq, _RC)

      if (nq > 0) then
         allocate(names(nq), _STAT)
         call ESMF_FieldBundleGet(bundle, itemorderflag=ESMF_ITEMORDER_ADDORDER, fieldNameList=names, _RC)
         if (.not.allocated(names0)) then
            allocate(names0(nq), _STAT)
            names0 = names
         end if
      end if

      ! Surface Geopotential from import state
      call MAPL_StateGetPointer(import, phis, "PHIS", _RC)

      phisxy = real(phis, kind=r8)

      ! Get tracers from import State (Note: Contains Updates from Analysis)
      call PULL_Q(self, import, qqq, NXQ, _RC)

      !-----------------------------
      ! end of fewer_tracers-section
      !-----------------------------

      do k = 1, nq ! size(names)
         pos = index(names(k), "::")
         if (pos > 0) then
            if ((names(k)(pos + 2:)) == "OX") ooo = vars%tracer(k)
         elseif (names(k) == "Q") then
            qqq = vars%tracer(k)
         end if
      end do

      ! WMP Begin REPLAY/ANA section
      call MAPL_GridCompGetResource(gc, "FV3_STANDALONE", FV3_STANDALONE, default=.false., _RC)
      is_shutoff = .true.
      ! if (.not. FV3_STANDALONE) then
      !    call MAPL_GridCompTimerStart(gc, "DYN_ANA", _RC)
      !    block
      !       integer :: iter, alarm_count
      !       type(ESMF_Alarm), allocatable :: alarm_list(:)
      !       character(len=ESMF_MAXSTR) :: alarm_name
      !       call ESMF_ClockGetAlarmList(clock, ESMF_ALARMLIST_ALL, alarmCount=alarm_count, _RC)
      !       _HERE, "alarm count: ", alarm_count
      !       allocate(alarm_list(alarm_count), _STAT)
      !       call ESMF_ClockGetAlarmList(clock, ESMF_ALARMLIST_ALL, alarmList=alarm_list, _RC)
      !       do iter = 1, alarm_count
      !          call ESMF_AlarmGet(alarm_list(iter), name=alarm_name, _RC)
      !          _HERE, "alarm ", iter, ": ", trim(alarm_name)
      !       end do
      !    end block
      !    call ESMF_ClockGetAlarm(clock, "REPLAYShutOff", alarm, _RC)
      !    is_shutoff = ESMF_AlarmIsRinging(alarm, _RC)
      ! end if

      is_shutoff = .true. ! TODO: pchakrab - Until we figure out the issue with alarm creation
      if (.not. is_shutoff) then
         ! If requested, do Intermittent Replay
         ! pchakrab - we are not doing this anymore, right?

         ! Create Local Copy of QV and OX (Contains Updates from Analysis)
         ox = 0.0
         qv = 0.0

         if (.not. ADIABATIC) then
            do k = 1, size(names)

               pos = index(names(k), "::")
               if (pos > 0) then
                  if ((names(k)(pos + 2:)) == "OX") then
                     if ((ooo%is_r4) .and. associated(ooo%content_r4)) then
                        if (size(ox) == size(ooo%content_r4)) then
                           ox = ooo%content_r4
                        end if
                     elseif (associated(ooo%content)) then
                        if (size(ox) == size(ooo%content)) then
                           ox = ooo%content
                        end if
                     end if
                  end if
               end if

               if (TRIM(names(k)) == "Q") then
                  if ((qqq%is_r4) .and. associated(qqq%content_r4)) then
                     if (size(qv) == size(qqq%content_r4)) then
                        qv = qqq%content_r4
                     end if
                  elseif (associated(qqq%content)) then
                     if (size(qv) == size(qqq%content)) then
                        qv = qqq%content
                     end if
                  end if
                  _ASSERT(all(qv >= 0.0), "Before AnaAddIncs: negative or nan water vapor detected")
               end if

            end do
         end if

         ! Diagnostics Before Analysis Increments are Added
         allocate(delpold(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
         allocate(qdnew(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
         allocate(qdold(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
         allocate(qvold(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
         allocate(qlold(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
         allocate(qiold(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
         allocate(qrold(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
         allocate(qsold(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
         allocate(qgold(ifirstxy:ilastxy, jfirstxy:jlastxy, km))

         call MAPL_StateGetPointer(import, dqvana, "DQVANA", _RC) ! Get QV Increment from Analysis
         call MAPL_StateGetPointer(import, dqlana, "DQLANA", _RC) ! Get QL Increment from Analysis
         call MAPL_StateGetPointer(import, dqiana, "DQIANA", _RC) ! Get QI Increment from Analysis
         call MAPL_StateGetPointer(import, dqrana, "DQRANA", _RC) ! Get QR Increment from Analysis
         call MAPL_StateGetPointer(import, dqsana, "DQSANA", _RC) ! Get QS Increment from Analysis
         call MAPL_StateGetPointer(import, dqgana, "DQGANA", _RC) ! Get QG Increment from Analysis
         call MAPL_StateGetPointer(import, doxana, "DOXANA", _RC) ! Get OX Increment from Analysis

         ql = 0.0
         qi = 0.0
         qr = 0.0
         qs = 0.0
         qg = 0.0
         do n = 1, size(names)
            if (TRIM(names(n)) == "QLCN" .or. TRIM(names(n)) == "QLLS") then
               if (self%vars%tracer(n)%is_r4) then
                  ql = ql + self%vars%tracer(n)%content_r4
               else
                  ql = ql + self%vars%tracer(n)%content
               end if
            end if
            if (TRIM(names(n)) == "QICN" .or. TRIM(names(n)) == "QILS") then
               if (self%vars%tracer(n)%is_r4) then
                  qi = qi + self%vars%tracer(n)%content_r4
               else
                  qi = qi + self%vars%tracer(n)%content
               end if
            end if
            if (TRIM(names(n)) == "QRAIN") then
               if (self%vars%tracer(n)%is_r4) then
                  qr = self%vars%tracer(n)%content_r4
               else
                  qr = self%vars%tracer(n)%content
               end if
            end if
            if (TRIM(names(n)) == "QSNOW") then
               if (self%vars%tracer(n)%is_r4) then
                  qs = self%vars%tracer(n)%content_r4
               else
                  qs = self%vars%tracer(n)%content
               end if
            end if
            if (TRIM(names(n)) == "QGRAUPEL") then
               if (self%vars%tracer(n)%is_r4) then
                  qg = self%vars%tracer(n)%content_r4
               else
                  qg = self%vars%tracer(n)%content
               end if
            end if
         end do
         qvold = qv - dqvana
         qlold = ql - dqlana
         qiold = qi - dqiana
         qrold = qr - dqrana
         qsold = qs - dqsana
         qgold = qg - dqgana

         ! Get A-grid winds
         call getAllWinds(vars%u, vars%v, ur=ur, vr=vr)

         delp = vars%pe(:, :, 2:) - vars%pe(:, :, :km) ! Pressure Thickness
         dmdt = vars%pe(:, :, km + 1) - vars%pe(:, :, 1) ! Psurf-Ptop
         tempxy = vars%pt * (1.0 + EPS * (qv - dqvana)) ! Compute THV Before Analysis Update

         if (do_energetics) then
            call Energetics(ur, vr, tempxy, vars%pe, delp, vars%pkz, phisxy, kenrg, penrg, tenrg)
         end if

         ! DUDTANA
         call MAPL_StateGetPointer(export, dudtana, "DUDTANA", _RC)
         if (associated(dudtana)) dudtana = ur

         ! DVDTANA
         call MAPL_StateGetPointer(export, dvdtana, "DVDTANA", _RC)
         if (associated(dvdtana)) dvdtana = vr

         ! DTDTANA
         call MAPL_StateGetPointer(export, dtdtana, "DTDTANA", _RC)
         if (associated(dtdtana)) dtdtana = vars%pt * vars%pkz

         ! DDELPDTANA
         call MAPL_StateGetPointer(export, ddpdtana, "DDELPDTANA", _RC)
         if (associated(ddpdtana)) ddpdtana = delp

         ! DTHVDTANAINT
         call MAPL_StateGetPointer(export, temp2d, "DTHVDTANAINT", _RC)
         if (associated(temp2d)) then
            tempxy = vars%pt * (1 + EPS * (qv - dqvana)) ! Set tempxy = TH*QVold (Before Analysis Update)
            allocate(dthdtanaint1(ifirstxy:ilastxy, jfirstxy:jlastxy))
            allocate(dthdtanaint2(ifirstxy:ilastxy, jfirstxy:jlastxy))
            dthdtanaint1 = 0.0
            do k = 1, km
               dthdtanaint1 = dthdtanaint1 + tempxy(:, :, k) * delp(:, :, k)
            end do
         end if

         ! DQVDTANAINT
         call MAPL_StateGetPointer(export, temp2d, "DQVDTANAINT", _RC)
         if (associated(temp2d)) then
            allocate(dqvdtanaint1(ifirstxy:ilastxy, jfirstxy:jlastxy))
            allocate(dqvdtanaint2(ifirstxy:ilastxy, jfirstxy:jlastxy))
            tempxy = qv - dqvana ! Set tempxy = QVold (Before Analysis Update)
            dqvdtanaint1 = 0.0
            do k = 1, km
               dqvdtanaint1 = dqvdtanaint1 + tempxy(:, :, k) * delp(:, :, k)
            end do
         end if

         ! DQLDTANAINT
         call MAPL_StateGetPointer(export, temp2d, "DQLDTANAINT", _RC)
         if (associated(temp2d)) then
            allocate(dqldtanaint1(ifirstxy:ilastxy, jfirstxy:jlastxy))
            allocate(dqldtanaint2(ifirstxy:ilastxy, jfirstxy:jlastxy))
            dqldtanaint1 = 0.0
            do n = 1, size(names)
               if (TRIM(names(n)) == "QLCN" .or. TRIM(names(n)) == "QLLS") then
                  do k = 1, km
                     if (self%vars%tracer(n)%is_r4) then
                        dqldtanaint1 = dqldtanaint1 + self%vars%tracer(n)%content_r4(:, :, k) * delp(:, :, k)
                     else
                        dqldtanaint1 = dqldtanaint1 + self%vars%tracer(n)%content(:, :, k) * delp(:, :, k)
                     end if
                  end do
               end if
            end do
            do k = 1, km
               dqldtanaint1 = dqldtanaint1 - dqlana(:, :, k) * delp(:, :, k)
            end do
         end if

         ! DQIDTANAINT
         call MAPL_StateGetPointer(export, temp2d, "DQIDTANAINT", _RC)
         if (associated(temp2d)) then
            allocate(dqidtanaint1(ifirstxy:ilastxy, jfirstxy:jlastxy))
            allocate(dqidtanaint2(ifirstxy:ilastxy, jfirstxy:jlastxy))
            dqidtanaint1 = 0.0
            do n = 1, size(names)
               if (TRIM(names(n)) == "QICN" .or. TRIM(names(n)) == "QILS") then
                  do k = 1, km
                     if (self%vars%tracer(n)%is_r4) then
                        dqidtanaint1 = dqidtanaint1 + self%vars%tracer(n)%content_r4(:, :, k) * delp(:, :, k)
                     else
                        dqidtanaint1 = dqidtanaint1 + self%vars%tracer(n)%content(:, :, k) * delp(:, :, k)
                     end if
                  end do
               end if
            end do
            do k = 1, km
               dqidtanaint1 = dqidtanaint1 - dqiana(:, :, k) * delp(:, :, k)
            end do
         end if

         ! DOXDTANAINT
         call MAPL_StateGetPointer(export, temp2d, "DOXDTANAINT", _RC)
         if (associated(temp2d)) then
            allocate(doxdtanaint1(ifirstxy:ilastxy, jfirstxy:jlastxy))
            allocate(doxdtanaint2(ifirstxy:ilastxy, jfirstxy:jlastxy))
            tempxy = ox - doxana ! Set tempxy = OXold (Before Analysis Update)
            doxdtanaint1 = 0.0
            do k = 1, km
               doxdtanaint1 = doxdtanaint1 + tempxy(:, :, k) * delp(:, :, k)
            end do
         end if

         ! Add Diabatic Forcing from Analysis to State Variables
         if (vars%nwat >= 6) then
            qdold = 1.0 - (qvold + qlold + qiold + qrold + qsold + qgold)
            qdnew = 1.0 - (qv + ql + qi + qr + qs + qg)
         else
            qdold = 1.0 - (qvold + qlold + qiold)
            qdnew = 1.0 - (qv + ql + qi)
         end if
         call MAPL_StateGetPointer(export, area, "AREA", _RC)

         allocate(trsum1(nq))
         allocate(trsum2(nq))

         call MAPL_GridCompGetResource(gc, "ANA_IS_WEIGHTED", ANA_IS_WEIGHTED, default="NO", _RC)
         ANA_IS_WEIGHTED = ESMF_UtilStringUpperCase(ANA_IS_WEIGHTED)
         is_weighted = adjustl(ANA_IS_WEIGHTED) == "YES" .or. adjustl(ANA_IS_WEIGHTED) == "NO"
         _ASSERT(is_weighted, "needs informative message")
         is_weighted = adjustl(ANA_IS_WEIGHTED) == "YES"

         ! Add Analysis Tendencies
         delpold = delp ! Old Pressure Thickness

         call ADD_INCS(esmfgrid, self, import, dt, is_weighted=is_weighted)

         if (DYN_DEBUG) call DEBUG_FV_STATE("ANA ADD_INCS", self)

         delp = vars%pe(:, :, 2:) - vars%pe(:, :, :km) ! Updated Pressure Thickness

         ! Compute Old Global Sums of Tracers over Locations where Mass has changed
         if ((.not.ADIABATIC)) then
            do n = 1, nq
               qsum1(:, :) = 0.0_r8
               if (self%vars%tracer(n)%is_r4) then
                  do k = 1, km
                     where (delp(:, :, k) /= delpold(:, :, k))
                        qsum1(:, :) = qsum1(:, :) + self%vars%tracer(n)%content_r4(:, :, k) * delpold(:, :, k)
                     end where
                  end do
               else
                  do k = 1, km
                     where (delp(:, :, k) /= delpold(:, :, k))
                        qsum1(:, :) = qsum1(:, :) + self%vars%tracer(n)%content(:, :, k) * delpold(:, :, k)
                     end where
                  end do
               end if
               where (qsum1 /= 0.0_r8)
                  qsum2 = qsum1
                  elsewhere
                  qsum2 = MAPL_UNDEFINED_REAL
               end where
               trsum1(n) = MAPL_AreaMean(qsum2, area, comm, _RC)
            end do
         end if

         ! Update Specific Mass of Aerosol Constituents Keeping Mixing_Ratio Constant WRT_Dry_Air After ANA Updates
         if ((.not. ADIABATIC)) then
            do n = 1, nq
               if ((TRIM(names(n)) /= "Q") .and. &
                    (TRIM(names(n)) /= "QLLS") .and. &
                    (TRIM(names(n)) /= "QLCN") .and. &
                    (TRIM(names(n)) /= "QILS") .and. &
                    (TRIM(names(n)) /= "QICN") .and. &
                    (TRIM(names(n)) /= "CLLS") .and. &
                    (TRIM(names(n)) /= "CLCN") .and. &
                    (TRIM(names(n)) /= "QRAIN") .and. &
                    (TRIM(names(n)) /= "QSNOW") .and. &
                    (TRIM(names(n)) /= "QGRAUPEL")) then
                  if (self%vars%tracer(n)%is_r4) then
                     self%vars%tracer(n)%content_r4 = self%vars%tracer(n)%content_r4 * (qdnew / qdold)
                  else
                     self%vars%tracer(n)%content = self%vars%tracer(n)%content * (qdnew / qdold)
                  end if
               end if
            end do
         end if

         ! Compute New Global Sums of Tracers over Locations where Mass has changed
         if ((.not. ADIABATIC)) then
            do n = 1, nq
               qsum1(:, :) = 0.0_r8
               if (self%vars%tracer(n)%is_r4) then
                  do k = 1, km
                     where (delp(:, :, k) /= delpold(:, :, k))
                        qsum1(:, :) = qsum1(:, :) + self%vars%tracer(n)%content_r4(:, :, k) * delp(:, :, k)
                     end where
                  end do
               else
                  do k = 1, km
                     where (delp(:, :, k) /= delpold(:, :, k))
                        qsum1(:, :) = qsum1(:, :) + self%vars%tracer(n)%content(:, :, k) * delp(:, :, k)
                     end where
                  end do
               end if
               where (qsum1 /= 0.0_r8)
                  qsum2 = qsum1
                  elsewhere
                  qsum2 = MAPL_UNDEFINED_REAL
               end where
               trsum2(n) = MAPL_AreaMean(qsum2, area, comm, _RC)
            end do
         end if

         ! Ensure Conservation of Global Mass of Aerosol Constituents After ANA Updates
         if ((.not. ADIABATIC)) then
            do n = 1, nq
               if ((TRIM(names(n)) /= "Q") .and. &
                    (TRIM(names(n)) /= "QLLS") .and. &
                    (TRIM(names(n)) /= "QLCN") .and. &
                    (TRIM(names(n)) /= "QILS") .and. &
                    (TRIM(names(n)) /= "QICN") .and. &
                    (TRIM(names(n)) /= "CLLS") .and. &
                    (TRIM(names(n)) /= "CLCN") .and. &
                    (TRIM(names(n)) /= "QRAIN") .and. &
                    (TRIM(names(n)) /= "QSNOW") .and. &
                    (TRIM(names(n)) /= "QGRAUPEL")) then

                  if (real(trsum1(n), kind=4) /= MAPL_UNDEFINED_REAL .and. &
                       real(trsum2(n), kind=4) /= MAPL_UNDEFINED_REAL) then
                     trsum2(n) = real(trsum1(n) / trsum2(n), kind=4)
                  else
                     trsum2(n) = 1.0d0
                  end if
                  ! IF (MAPL_AM_I_ROOT()) print *, trim(names(n))," ratio is: ",trsum2(n)

                  if (self%vars%tracer(n)%is_r4) then
                     do k = 1, km
                        where (delp(:, :, k) /= delpold(:, :, k))
                           self%vars%tracer(n)%content_r4(:, :, k) = self%vars%tracer(n)%content_r4(:, :, k) * trsum2(n)
                        end where
                     end do
                  else
                     do k = 1, km
                        where (delp(:, :, k) /= delpold(:, :, k))
                           self%vars%tracer(n)%content(:, :, k) = self%vars%tracer(n)%content(:, :, k) * trsum2(n)
                        end where
                     end do
                  end if
               end if
            end do
         end if

         deallocate(trsum1)
         deallocate(trsum2)

         ! Update Local Copy of QV and OX to account for Global Sum Adjustment
         do k = 1, size(names)
            pos = index(names(k), "::")
            if (pos > 0) then
               if ((names(k)(pos + 2:)) == "OX") then
                  if (ooo%is_r4) then
                     ox = ooo%content_r4
                  else
                     ox = ooo%content
                  end if
               end if
            end if
            if (TRIM(names(k)) == "Q") then
               if (qqq%is_r4) then
                  qv = qqq%content_r4
               else
                  qv = qqq%content
               end if
               _ASSERT(all(qv >= 0.0), "After AnaAddIncs: negative or nan water vapor detected")
            end if
         end do

         ! Diagnostics After Analysis Increments are Added
         call MAPL_StateGetPointer(export, temp2d, "DMDTANA", _RC)
         if (associated(temp2d)) temp2d = ((vars%pe(:, :, km + 1) - vars%pe(:, :, 1)) - dmdt) / (GRAV * dt)

         call getAllWinds(vars%u, vars%v, uc=uc0, vc=vc0, ur=ur, vr=vr)

         dmdt = vars%pe(:, :, km + 1) - vars%pe(:, :, 1) ! Psurf-Ptop

         ! DUDTANA
         call MAPL_StateGetPointer(export, dudtana, "DUDTANA", _RC)
         if (associated(dudtana)) then
            dudtana = (ur - dudtana) / dt
         end if

         ! DVDTANA
         call MAPL_StateGetPointer(export, dvdtana, "DVDTANA", _RC)
         if (associated(dvdtana)) then
            dvdtana = (vr - dvdtana) / dt
         end if

         ! DTDTANA
         call MAPL_StateGetPointer(export, dtdtana, "DTDTANA", _RC)
         if (associated(dtdtana)) then
            dtdtana = ((vars%pt * vars%pkz) - dtdtana) / dt
         end if

         ! DDELPDTANA
         call MAPL_StateGetPointer(export, ddpdtana, "DDELPDTANA", _RC)
         if (associated(ddpdtana)) then
            ddpdtana = (delp - ddpdtana) / dt
         end if

         ! DTHVDTANAINT
         ! ------------
         call MAPL_StateGetPointer(export, temp2d, "DTHVDTANAINT", _RC)
         if (associated(temp2d)) then
            tempxy = vars%pt * (1 + EPS * qv) ! Set tempxy = TH*QVnew (After Analysis Update)
            dthdtanaint2 = 0.0
            do k = 1, km
               dthdtanaint2 = dthdtanaint2 + tempxy(:, :, k) * delp(:, :, k)
            end do
            temp2d = (dthdtanaint2 - dthdtanaint1) * MAPL_P00**MAPL_KAPPA / (MAPL_GRAV * dt)
            deallocate(dthdtanaint1)
            deallocate(dthdtanaint2)
         end if

         ! DQVDTANAINT
         call MAPL_StateGetPointer(export, temp2d, "DQVDTANAINT", _RC)
         if (associated(temp2d)) then
            tempxy = qv ! Set tempxy = QNEW (After Analysis Update)
            dqvdtanaint2 = 0.0
            do k = 1, km
               dqvdtanaint2 = dqvdtanaint2 + tempxy(:, :, k) * delp(:, :, k)
            end do
            temp2d = (dqvdtanaint2 - dqvdtanaint1) / (MAPL_GRAV * dt)
            deallocate(dqvdtanaint1)
            deallocate(dqvdtanaint2)
         end if

         ! DQLDTANAINT
         call MAPL_StateGetPointer(export, temp2d, "DQLDTANAINT", _RC)
         if (associated(temp2d)) then
            dqldtanaint2 = 0.0
            do n = 1, size(names)
               if (TRIM(names(n)) == "QLCN" .or. &
                    TRIM(names(n)) == "QLLS") then
                  do k = 1, km
                     if (self%vars%tracer(n)%is_r4) then
                        dqldtanaint2 = dqldtanaint2 + self%vars%tracer(n)%content_r4(:, :, k) * delp(:, :, k)
                     else
                        dqldtanaint2 = dqldtanaint2 + self%vars%tracer(n)%content(:, :, k) * delp(:, :, k)
                     end if
                  end do
               end if
            end do
            temp2d = (dqldtanaint2 - dqldtanaint1) / (MAPL_GRAV * dt)
            deallocate(dqldtanaint1)
            deallocate(dqldtanaint2)
         end if

         ! DQIDTANAINT
         call MAPL_StateGetPointer(export, temp2d, "DQIDTANAINT", _RC)
         if (associated(temp2d)) then
            dqidtanaint2 = 0.0
            do n = 1, size(names)
               if (TRIM(names(n)) == "QICN" .or. &
                    TRIM(names(n)) == "QILS") then
                  do k = 1, km
                     if (self%vars%tracer(n)%is_r4) then
                        dqidtanaint2 = dqidtanaint2 + self%vars%tracer(n)%content_r4(:, :, k) * delp(:, :, k)
                     else
                        dqidtanaint2 = dqidtanaint2 + self%vars%tracer(n)%content(:, :, k) * delp(:, :, k)
                     end if
                  end do
               end if
            end do
            temp2d = (dqidtanaint2 - dqidtanaint1) / (MAPL_GRAV * dt)
            deallocate(dqidtanaint1)
            deallocate(dqidtanaint2)
         end if

         ! DOXDTANAINT
         call MAPL_StateGetPointer(export, temp2d, "DOXDTANAINT", _RC)
         if (associated(temp2d)) then
            tempxy = ox ! Set tempxy = OXnew (After Analysis Update)
            doxdtanaint2 = 0.0
            do k = 1, km
               doxdtanaint2 = doxdtanaint2 + tempxy(:, :, k) * delp(:, :, k)
            end do
            temp2d = (doxdtanaint2 - doxdtanaint1) * (MAPL_O3MW / MAPL_AIRMW) / (MAPL_GRAV * dt)
            deallocate(doxdtanaint1)
            deallocate(doxdtanaint2)
         end if

         deallocate(delpold)
         deallocate(qdnew)
         deallocate(qdold)
         deallocate(qvold)
         deallocate(qlold)
         deallocate(qiold)
         deallocate(qrold)
         deallocate(qsold)
         deallocate(qgold)

         ! WMP End ANA section
      else ! REPLAY/ANA is_shutoff

         ox = 0.0
         qv = 0.0
         if (.not. ADIABATIC) then
            do k = 1, size(names)
               pos = index(names(k), "::")
               if (pos > 0) then
                  if ((names(k)(pos + 2:)) == "OX") then
                     if ((ooo%is_r4) .and. associated(ooo%content_r4)) then
                        if (size(ox) == size(ooo%content_r4)) then
                           ox = ooo%content_r4
                        end if
                     elseif (associated(ooo%content)) then
                        if (size(ox) == size(ooo%content)) then
                           ox = ooo%content
                        end if
                     end if
                  end if
               end if
               if (TRIM(names(k)) == "Q") then
                  if ((qqq%is_r4) .and. associated(qqq%content_r4)) then
                     if (size(qv) == size(qqq%content_r4)) then
                        qv = qqq%content_r4
                     end if
                  elseif (associated(qqq%content)) then
                     if (size(qv) == size(qqq%content)) then
                        qv = qqq%content
                     end if
                  end if
                  _ASSERT(all(qv >= 0.0), "DYN_ANA: negative or nan water vapor detected")
               end if
            end do
         end if
         call getAllWinds(vars%u, vars%v, uc=uc0, vc=vc0, ur=ur, vr=vr)
         delp = vars%pe(:, :, 2:) - vars%pe(:, :, :km) ! Pressure Thickness
         dmdt = vars%pe(:, :, km + 1) - vars%pe(:, :, 1) ! Psurf-Ptop
         tempxy = vars%pt * (1.0 + EPS * qv)
         if (do_energetics) then
            call Energetics(ur, vr, tempxy, vars%pe, delp, vars%pkz, phisxy, kenrg, penrg, tenrg)
         end if

      end if
      if (.not. FV3_STANDALONE) then
         ! call MAPL_GridCompTimerStop(gc, "DYN_ANA", _RC)
      end if

      call MAPL_GridCompTimerStart(gc, "DYN_PROLOGUE", _RC)
      ! Create FV Thermodynamic Variables
      tempxy = vars%pt * vars%pkz ! Compute Dry Temperature

      ! Initialize Diagnostic Dynamics Tendencies
      dpedt = vars%pe ! Edge Pressure      Tendency
      ddpdt = delp ! Pressure Thickness Tendency
      dudt = ur ! U-Wind on A-Grid   Tendency
      dvdt = vr ! V-Wind on A-Grid   Tendency
      dtdt = tempxy ! Dry Temperature    Tendency
      dqdt = qv ! Specific Humidity  Tendency
      dthdt = vars%pt * (1.0 + EPS * qv) * delp

      call FILLOUT3(export, "QV_DYN_IN", qv, _RC)
      call FILLOUT3(export, "T_DYN_IN", tempxy, _RC)
      call FILLOUT3(export, "U_DYN_IN", ur, _RC)
      call FILLOUT3(export, "V_DYN_IN", vr, _RC)
      call FILLOUT3(export, "PLE_DYN_IN", vars%pe, _RC)
      call FILLOUT3(export, "PLE4", vars%pe, _RC)

      ! Initialize 3-D Tracer Dynamics Tendencies
      call MAPL_StateGetPointer(export, dqldt, "DQLDTDYN", _RC)
      call MAPL_StateGetPointer(export, dqidt, "DQIDTDYN", _RC)
      call MAPL_StateGetPointer(export, doxdt, "DOXDTDYN", _RC)

      if (allocated(names)) then

         if (associated(dqldt)) then
            dqldt = 0.0
            do k = 1, size(names)
               if (TRIM(names(k)) == "QLCN" .or. &
                    TRIM(names(k)) == "QLLS") then
                  if (self%vars%tracer(k)%is_r4) then
                     if (size(dqldt) == size(self%vars%tracer(k)%content_r4)) &
                          dqldt = dqldt - self%vars%tracer(k)%content_r4
                  else
                     if (size(dqldt) == size(self%vars%tracer(k)%content)) &
                          dqldt = dqldt - self%vars%tracer(k)%content
                  end if
               end if
            end do
         end if

         if (associated(dqidt)) then
            dqidt = 0.0
            do k = 1, size(names)
               if (TRIM(names(k)) == "QICN" .or. &
                    TRIM(names(k)) == "QILS") then
                  if (self%vars%tracer(k)%is_r4) then
                     if (size(dqidt) == size(self%vars%tracer(k)%content_r4)) &
                          dqidt = dqidt - self%vars%tracer(k)%content_r4
                  else
                     if (size(dqidt) == size(self%vars%tracer(k)%content)) &
                          dqidt = dqidt - self%vars%tracer(k)%content
                  end if
               end if
            end do
         end if

         if (associated(doxdt)) then
            doxdt = 0.0
            do k = 1, size(names)
               pos = index(names(k), "::")
               if (pos > 0) then
                  if ((names(k)(pos + 2:)) == "OX") then
                     if (self%vars%tracer(k)%is_r4) then
                        if (size(doxdt) == size(self%vars%tracer(k)%content_r4)) &
                             doxdt = doxdt - self%vars%tracer(k)%content_r4
                     else
                        if (size(doxdt) == size(self%vars%tracer(k)%content)) &
                             doxdt = doxdt - self%vars%tracer(k)%content
                     end if
                  end if
               end if
            end do
         end if
      end if

      ! Initialize 2-D Vertically Integrated Tracer Dynamics Tendencies
      call MAPL_StateGetPointer(export, temp2d, "DQVDTDYNINT", _RC)
      if (associated(temp2d)) then
         temp2d = 0.0
         do k = 1, km
            temp2d = temp2d - qv(:, :, k) * delp(:, :, k)
         end do
      end if

      call MAPL_StateGetPointer(export, temp2d, "DQLDTDYNINT", _RC)
      if (associated(temp2d)) then
         temp2d = 0.0
         do n = 1, size(names)
            if (TRIM(names(n)) == "QLCN" .or. &
                 TRIM(names(n)) == "QLLS") then
               if (self%vars%tracer(n)%is_r4) then
                  do k = 1, km
                     temp2d = temp2d - self%vars%tracer(n)%content_r4(:, :, k) * delp(:, :, k)
                  end do
               else
                  do k = 1, km
                     temp2d = temp2d - self%vars%tracer(n)%content(:, :, k) * delp(:, :, k)
                  end do
               end if
            end if
         end do
      end if

      call MAPL_StateGetPointer(export, temp2d, "DQIDTDYNINT", _RC)
      if (associated(temp2d)) then
         temp2d = 0.0
         do n = 1, size(names)
            if (TRIM(names(n)) == "QICN" .or. &
                 TRIM(names(n)) == "QILS") then
               if (self%vars%tracer(n)%is_r4) then
                  do k = 1, km
                     temp2d = temp2d - self%vars%tracer(n)%content_r4(:, :, k) * delp(:, :, k)
                  end do
               else
                  do k = 1, km
                     temp2d = temp2d - self%vars%tracer(n)%content(:, :, k) * delp(:, :, k)
                  end do
               end if
            end if
         end do
      end if

      call MAPL_StateGetPointer(export, temp2d, "DOXDTDYNINT", _RC)
      if (associated(temp2d)) then
         temp2d = 0.0
         do n = 1, size(names)
            pos = index(names(n), "::")
            if (pos > 0) then
               if ((names(n)(pos + 2:)) == "OX") then
                  if (self%vars%tracer(n)%is_r4) then
                     do k = 1, km
                        temp2d = temp2d - self%vars%tracer(n)%content_r4(:, :, k) * delp(:, :, k)
                     end do
                  else
                     do k = 1, km
                        temp2d = temp2d - self%vars%tracer(n)%content(:, :, k) * delp(:, :, k)
                     end do
                  end if
               end if
            end if
         end do
      end if

      ! Compute Energetics After Analysis (and Before Dycore)
      tempxy = vars%pt * (1.0 + EPS * qv) ! Compute THV After Analysis Update

      if (do_energetics) then
         call Energetics(ur, vr, tempxy, vars%pe, delp, vars%pkz, phisxy, kenrg0, penrg0, tenrg0)
         kenrg = (kenrg0 - kenrg) / dt
         penrg = (penrg0 - penrg) / dt
         tenrg = (tenrg0 - tenrg) / dt
         call FILLOUT2(export, "KEANA", kenrg, _RC)
         call FILLOUT2(export, "PEANA", penrg, _RC)
         call FILLOUT2(export, "TEANA", tenrg, _RC)
      end if

      ! Call Wrapper (DynRun) for FVDycore
      call MAPL_GridCompGetResource(gc, "CONSV", CONSV, default=1, _RC)
      call MAPL_GridCompGetResource(gc, "FILL", FILL, default=0, _RC)

      LCONSV = CONSV == 1
      LFILL = FILL == 1

      ! Get pressures before dynamics
      pe0 = vars%pe

      call MAPL_GridCompTimerStop(gc, "DYN_PROLOGUE", _RC)

      !-------------------------------------------------------

      call MAPL_GridCompTimerStart(gc, "DYN_CORE", _RC)
      t1 = MPI_Wtime(status)
      call DynRun(self, export, clock, gc, PLE0=pe0, _RC)
      t2 = MPI_Wtime(status)
      dyn_run_timer = t2 - t1
      call MAPL_GridCompTimerStop(gc, "DYN_CORE", _RC)

      call MAPL_GridCompTimerStart(gc, "DYN_EPILOGUE", _RC)
      ! Computational diagnostics
      call MAPL_StateGetPointer(export, temp2d, "TIME_IN_DYN", _RC)
      if (associated(temp2d)) temp2d = dyn_run_timer
      call MAPL_StateGetPointer(export, temp2d, "PID", _RC)
      if (associated(temp2d)) temp2d = 0 !WMP need to get from MAPL gid

      !#define DEBUG_WINDS
#if defined(DEBUG_WINDS)
      call Write_Profile(grid, vars%u, "U-after-DynRun")
      call Write_Profile(grid, vars%v, "V-after-DynRun")
#endif
      plk = exp(KAPPA * log(0.5 * (vars%pe(:, :, 1:km) + vars%pe(:, :, 2:km + 1))))
      pkxy = exp(KAPPA * log(vars%pe))

      !----------------------------------------------------------------------------

      if (SW_DYNAMICS) then

         call MAPL_StateGetPointer(export, temp2d, "PHIS", _RC)
         if (associated(temp2d)) temp2d = phisxy

         call MAPL_StateGetPointer(export, temp2d, "PS", _RC)
         if (associated(temp2d)) temp2d = vars%pe(:, :, km + 1) / GRAV

         call getAllWinds(vars%u, vars%v, ua=ua, va=va, uc=uc, vc=vc, ur=ur, vr=vr)
         call FILLOUT3(export, "U_DGRID", vars%u, _RC)
         call FILLOUT3(export, "V_DGRID", vars%v, _RC)
         call FILLOUT3(export, "U_CGRID", uc, _RC)
         call FILLOUT3(export, "V_CGRID", vc, _RC)
         call FILLOUT3(export, "U_AGRID", ua, _RC)
         call FILLOUT3(export, "V_AGRID", va, _RC)

         call FILLOUT3(export, "U", ur, _RC)
         call FILLOUT3(export, "V", vr, _RC)

      else ! .not. SW_DYNAMICS

         ! Load Local Variable with Vapor Specific Humidity
         if ((.not. ADIABATIC) .and. (self%grid%nq > 0)) then
            if (qqq%is_r4) then
               if (size(qv) == size(qqq%content_r4)) qv = qqq%content_r4
            else
               if (size(qv) == size(qqq%content)) qv = qqq%content
            end if
            _ASSERT(all(qv >= 0.0), "After DynRun: negative or nan water vapor detected")
         else
            qv = 0.0
         end if

         ! Vertically Integrated THV Tendency Diagnostic
         delp = (vars%pe(:, :, 2:) - vars%pe(:, :, :km))
         dthdt = (vars%pt * (1.0 + EPS * qv) * delp - dthdt) / dt

         call MAPL_StateGetPointer(export, temp2d, "DTHVDTDYNINT", _RC)
         if (associated(temp2d)) then
            qsum1 = 0.0
            do k = 1, km
               qsum1 = qsum1 + dthdt(:, :, k)
            end do
            temp2d = qsum1 * (MAPL_P00**MAPL_KAPPA) / GRAV
         end if

         ! Compute Dry Theta and T with Unified Poles
         tempxy = vars%pt * vars%pkz

         ! Compute Mid-Layer Pressure and Pressure Thickness
         delp = (vars%pe(:, :, 2:) - vars%pe(:, :, :km))
         pl = (vars%pe(:, :, 2:) + vars%pe(:, :, :km)) * 0.5

         ! Get all wind derivatives
         call getAllWinds(vars%u, vars%v, ua=ua, va=va, uc=uc, vc=vc, ur=ur, vr=vr)
         call getVorticity(vars%u, vars%v, vort)
         call getDivergence(uc, vc, divg)

         ! Compute absolute vorticity on the D grid
         call getEPV(vars%pt, vort, ua, va, epvxyz)
         call MAPL_StateGetPointer(export, temp3d, "EPV", _RC)
         if (associated(temp3d)) temp3d = epvxyz * (P00**KAPPA)

         ! Compute Tropopause Pressure, Temperature, and Moisture
         do_tropvars = .false.
         call MAPL_StateGetPointer(export, temp2d, "TROPP_THERMAL", _RC)
         if (associated(temp2d)) do_tropvars = .true.
         call MAPL_StateGetPointer(export, temp2d, "TROPP_EPV", _RC)
         if (associated(temp2d)) do_tropvars = .true.
         call MAPL_StateGetPointer(export, temp2d, "TROPP_BLENDED", _RC)
         if (associated(temp2d)) do_tropvars = .true.
         call MAPL_StateGetPointer(export, temp2d, "TROPK_BLENDED", _RC)
         if (associated(temp2d)) do_tropvars = .true.
         call MAPL_StateGetPointer(export, temp2d, "TROPT", _RC)
         if (associated(temp2d)) do_tropvars = .true.
         call MAPL_StateGetPointer(export, temp2d, "TROPQ", _RC)
         if (associated(temp2d)) do_tropvars = .true.

         if (do_tropvars) then
            allocate(tropp1(ifirstxy:ilastxy, jfirstxy:jlastxy))
            allocate(tropp2(ifirstxy:ilastxy, jfirstxy:jlastxy))
            allocate(tropp3(ifirstxy:ilastxy, jfirstxy:jlastxy))
            allocate(tropt(ifirstxy:ilastxy, jfirstxy:jlastxy))
            allocate(tropq(ifirstxy:ilastxy, jfirstxy:jlastxy))
            call tropovars( &
                 ilastxy - ifirstxy + 1, jlastxy - jfirstxy + 1, km, &
                 real(vars%pe, kind=4), &
                 real(pl, kind=4), &
                 real(tempxy, kind=4), &
                 real(qv, kind=4), &
                 real(epvxyz * (P00**KAPPA), kind=4), &
                 tropp1, tropp2, tropp3, tropt, tropq)

            ! get blended index
            call MAPL_StateGetPointer(export, temp2d, 'TROPK_BLENDED', _RC)
            if (associated(temp2d)) then
               kend = km
               do j = jfirstxy, jlastxy
                  do i = ifirstxy, ilastxy
                     if (tropp3(i, j) /= MAPL_UNDEFINED_REAL) then
                        kend = 1
                        do while (vars%pe(i, j, kend) <= tropp3(i, j))
                           kend = kend + 1
                        end do
                     else
                        kend = 1
                        do while (vars%pe(i, j, kend) <= 40000.0)
                           kend = kend + 1
                        end do
                     end if
                     temp2d(i - ifirstxy + 1, j - jfirstxy + 1) = kend
                  end do
               end do
            end if

            call MAPL_StateGetPointer(export, temp2d, 'TROPP_THERMAL', _RC)
            if (associated(temp2d)) temp2d = tropp1

            call MAPL_StateGetPointer(export, temp2d, 'TROPP_EPV', _RC)
            if (associated(temp2d)) temp2d = tropp2

            call MAPL_StateGetPointer(export, temp2d, 'TROPP_BLENDED', _RC)
            if (associated(temp2d)) temp2d = tropp3

            call MAPL_StateGetPointer(export, temp2d, 'TROPT', _RC)
            if (associated(temp2d)) temp2d = tropt

            call MAPL_StateGetPointer(export, temp2d, 'TROPQ', _RC)
            if (associated(temp2d)) temp2d = tropq

            deallocate(tropp1)
            deallocate(tropp2)
            deallocate(tropp3)
            deallocate(tropt)
            deallocate(tropq)
         end if

         ! Get Cubed-Sphere Wind Exports
         call FILLOUT3(export, 'U_DGRID', vars%u, _RC)
         call FILLOUT3(export, 'V_DGRID', vars%v, _RC)
         call FILLOUT3(export, 'U_CGRID', uc, _RC)
         call FILLOUT3(export, 'V_CGRID', vc, _RC)
         call FILLOUT3(export, 'U_AGRID', ua, _RC)
         call FILLOUT3(export, 'V_AGRID', va, _RC)

         ! if (DEBUG_DYN) then
         !    block
         !       real :: maxmin(2)
         !       maxmin = MAPL_MaxMin(qv, comm, _RC)
         !       call logger%info("max/min(Q_AF_DYN): %f/%f", maxmin(1), maxmin(2))
         !       maxmin = MAPL_MaxMin(tempxy, comm, _RC)
         !       call logger%info("max/min(T_AF_DYN): %f/%f", maxmin(1), maxmin(2))
         !       maxmin = MAPL_MaxMin(ua, comm, _RC)
         !       call logger%info("max/min(U_AF_DYN): %f/%f", maxmin(1), maxmin(2))
         !       maxmin = MAPL_MaxMin(va, comm, _RC)
         !       call logger%info("max/min(V_AF_DYN): %f/%f", maxmin(1), maxmin(2))
         !    end block
         ! end if

         ! Compute Diagnostic Dynamics Tendencies
         !  (Note: initial values of d(m,u,v,T,q)/dt are progs m,u,v,T,q)
         dmdt = (vars%pe(:, :, km + 1) - vars%pe(:, :, 1) - dmdt) / (GRAV * dt)

         dudt = (ur - dudt) / dt
         dvdt = (vr - dvdt) / dt
         dtdt = (tempxy - dtdt) / dt
         dqdt = (qv - dqdt) / dt

         dpedt = (vars%pe - dpedt) / dt
         ddpdt = (delp - ddpdt) / dt ! Pressure Thickness Tendency

         call FILLOUT3(export, 'DELP', delp, _RC)
         call FILLOUT3(export, 'DUDTDYN', dudt, _RC)
         call FILLOUT3(export, 'DVDTDYN', dvdt, _RC)
         call FILLOUT3(export, 'DTDTDYN', dtdt, _RC)
         call FILLOUT3(export, 'DQVDTDYN', dqdt, _RC)
         call FILLOUT3(export, 'DDELPDTDYN', ddpdt, _RC)
         ! pchakrab - TODO: figure out the issue with DPLEDTDYN
         ! call FILLOUT3(export, 'DPLEDTDYN' , dpedt, _RC)

         ! fill pressure exports (PLE0: Before) & (PLE1: After) from FV3
         call FILLOUT3r8(export, 'PLE0', pe0, _RC)
         pe1 = vars%pe
         call FILLOUT3r8(export, 'PLE1', pe1, _RC)

         if (AdvCore_Advection == 2) then
            ! Compute time-centered C-Grid Courant Numbers and Mass Fluxes on Cubed Orientation
            uc0 = 0.5 * (uc + uc0)
            vc0 = 0.5 * (vc + vc0)
            pe0 = 0.5 * (pe1 + pe0)
            call computeMassFluxes(uc0, vc0, pe0, mfxxyz, mfyxyz, cxxyz, cyxyz, dt)
            call FILLOUT3r8(export, 'CX', cxxyz, _RC)
            call FILLOUT3r8(export, 'CY', cyxyz, _RC)
            call FILLOUT3r8(export, 'MFX', mfxxyz, _RC)
            call FILLOUT3r8(export, 'MFY', mfyxyz, _RC)
         else
            ! Fill Advection C-Grid Courant Numbers and Mass Fluxes on Cubed Orientation from FV3 DynCore
            call fillMassFluxes(mfxxyz, mfyxyz, cxxyz, cyxyz)
            call FILLOUT3r8(export, 'CX', cxxyz, _RC)
            call FILLOUT3r8(export, 'CY', cyxyz, _RC)
            call FILLOUT3r8(export, 'MFX', mfxxyz, _RC)
            call FILLOUT3r8(export, 'MFY', mfyxyz, _RC)
         end if

         call FILLOUT3(export, 'CU', cxxyz, _RC)
         call FILLOUT3(export, 'CV', cyxyz, _RC)
         call FILLOUT3(export, 'MX', mfxxyz, _RC)
         call FILLOUT3(export, 'MY', mfyxyz, _RC)

         ! Compute and return the vertical mass flux
         call getVerticalMassFlux(mfxxyz, mfyxyz, mfzxyz, dt)
         call FILLOUT3r8(export, 'MFZ', mfzxyz, _RC)
         call FILLOUT3(export, 'MZ', mfzxyz, _RC)

         call FILLOUT3(export, 'U', ur, _RC)
         call FILLOUT3(export, 'V', vr, _RC)
         call FILLOUT3(export, 'T', tempxy, _RC)
         call FILLOUT3(export, 'Q', qv, _RC)
         call FILLOUT3(export, 'PL', pl, _RC)
         call FILLOUT3(export, 'PLK', plk, _RC)
         call FILLOUT3(export, 'PKE', pkxy, _RC)
         call FILLOUT3(export, 'PT', vars%pt, _RC)
         call FILLOUT3(export, 'PE', vars%pe, _RC)

#ifdef SKIP_TRACERS
         do itracer = 1, ntracers
            write(myTracer, "('Q',i5.5)") itracer - 1
            call MAPL_StateGetPointer(export, temp3d, TRIM(myTracer), _RC)
            if ((associated(temp3d)) .and. (nq >= itracer)) then
               if (self%vars%tracer(itracer)%is_r4) then
                  temp3d = self%vars%tracer(itracer)%content_r4
               else
                  temp3d = self%vars%tracer(itracer)%content
               end if
            end if
         end do
#endif

         call MAPL_StateGetPointer(export, temp3d, 'PV', _RC)
         if (associated(temp3d)) temp3d = epvxyz / vars%pt

         call MAPL_StateGetPointer(export, temp3d, 'S', _RC)
         if (associated(temp3d)) temp3d = tempxy * CP

         call MAPL_StateGetPointer(export, temp3d, 'TH', _RC)
         !   if(associated(temp3d)) temp3d = vars%pt*(p00**kappa)
         if (associated(temp3d)) then
            temp3d = (tempxy) * (P00 / (0.5 * (vars%pe(:, :, 1:km) + vars%pe(:, :, 2:km + 1))))**KAPPA
         end if

         call MAPL_StateGetPointer(export, temp2d, 'DMDTDYN', _RC)
         if (associated(temp2d)) temp2d = dmdt

         ! Compute 3-D Tracer Dynamics Tendencies
         call MAPL_StateGetPointer(export, qctmp, 'QC', _RC)

         if (associated(qctmp)) then
            qctmp = 0.0
            do k = 1, size(names)
               if (TRIM(names(k)) == 'QLCN' .or. &
                    TRIM(names(k)) == 'QILS' .or. &
                    TRIM(names(k)) == 'QICN' .or. &
                    TRIM(names(k)) == 'QLLS') then
                  if (self%vars%tracer(k)%is_r4) then
                     if (size(dqldt) == size(self%vars%tracer(k)%content_r4)) &
                          qctmp = qctmp + self%vars%tracer(k)%content_r4
                  else
                     if (size(dqldt) == size(self%vars%tracer(k)%content)) &
                          qctmp = qctmp + self%vars%tracer(k)%content
                  end if
               end if
            end do
         end if

         if (associated(dqldt)) then
            do n = 1, size(names)
               if (TRIM(names(n)) == 'QLCN' .or. &
                    TRIM(names(n)) == 'QLLS') then
                  if (self%vars%tracer(n)%is_r4) then
                     dqldt = dqldt + self%vars%tracer(n)%content_r4
                  else
                     dqldt = dqldt + self%vars%tracer(n)%content
                  end if
               end if
            end do
            dqldt = dqldt / dt
         end if

         if (associated(dqidt)) then
            do n = 1, size(names)
               if (TRIM(names(n)) == 'QICN' .or. &
                    TRIM(names(n)) == 'QILS') then
                  if (self%vars%tracer(n)%is_r4) then
                     dqidt = dqidt + self%vars%tracer(n)%content_r4
                  else
                     dqidt = dqidt + self%vars%tracer(n)%content
                  end if
               end if
            end do
            dqidt = dqidt / dt
         end if

         if (associated(doxdt)) then
            do n = 1, size(names)
               pos = index(names(n), '::')
               if (pos > 0) then
                  if ((names(n)(pos + 2:)) == 'OX') then
                     if (self%vars%tracer(n)%is_r4) then
                        doxdt = doxdt + self%vars%tracer(n)%content_r4
                     else
                        doxdt = doxdt + self%vars%tracer(n)%content
                     end if
                  end if
               end if
            end do
            doxdt = doxdt / dt
         end if

         ! Compute 2-D Vertically Integrated Tracer Dynamics Tendencies
         call MAPL_StateGetPointer(export, temp2d, 'DQVDTDYNINT', _RC)
         if (associated(temp2d)) then
            do k = 1, km
               temp2d = temp2d + qv(:, :, k) * delp(:, :, k)
            end do
            temp2d = temp2d / (GRAV * dt)
         end if

         call MAPL_StateGetPointer(export, temp2d, 'DQLDTDYNINT', _RC)
         if (associated(temp2d)) then
            do n = 1, size(names)
               if (TRIM(names(n)) == 'QLCN' .or. &
                    TRIM(names(n)) == 'QLLS') then
                  if (self%vars%tracer(n)%is_r4) then
                     do k = 1, km
                        temp2d = temp2d + self%vars%tracer(n)%content_r4(:, :, k) * delp(:, :, k)
                     end do
                  else
                     do k = 1, km
                        temp2d = temp2d + self%vars%tracer(n)%content(:, :, k) * delp(:, :, k)
                     end do
                  end if
               end if
            end do
            temp2d = temp2d / (GRAV * dt)
         end if

         call MAPL_StateGetPointer(export, temp2d, 'DQIDTDYNINT', _RC)
         if (associated(temp2d)) then
            do n = 1, size(names)
               if (TRIM(names(n)) == 'QICN' .or. &
                    TRIM(names(n)) == 'QILS') then
                  if (self%vars%tracer(n)%is_r4) then
                     do k = 1, km
                        temp2d = temp2d + self%vars%tracer(n)%content_r4(:, :, k) * delp(:, :, k)
                     end do
                  else
                     do k = 1, km
                        temp2d = temp2d + self%vars%tracer(n)%content(:, :, k) * delp(:, :, k)
                     end do
                  end if
               end if
            end do
            temp2d = temp2d / (GRAV * dt)
         end if

         call MAPL_StateGetPointer(export, temp2d, 'DOXDTDYNINT', _RC)
         if (associated(temp2d)) then
            do n = 1, size(names)
               pos = index(names(n), '::')
               if (pos > 0) then
                  if ((names(n)(pos + 2:)) == 'OX') then
                     if (self%vars%tracer(n)%is_r4) then
                        do k = 1, km
                           temp2d = temp2d + self%vars%tracer(n)%content_r4(:, :, k) * delp(:, :, k)
                        end do
                     else
                        do k = 1, km
                           temp2d = temp2d + self%vars%tracer(n)%content(:, :, k) * delp(:, :, k)
                        end do
                     end if
                  end if
               end if
            end do
            temp2d = temp2d * (MAPL_O3MW / MAPL_AIRMW) / (MAPL_GRAV * dt)
         end if

         ! Virtual temperature
         tempxy = tempxy * (1.0 + EPS * qv)

         call MAPL_StateGetPointer(export, temp3d, 'TV', _RC)
         if (associated(temp3d)) temp3d = tempxy

         ! Fluxes: UCPT & VCPT
         !--------------------
         call MAPL_StateGetPointer(export, temp2d, 'UCPT', _RC)
         if (associated(temp2d)) then
            temp2d = 0.0
            do k = 1, km
               temp2d = temp2d + ur(:, :, k) * tempxy(:, :, k) * delp(:, :, k)
            end do
            temp2d = temp2d * (CP / GRAV)
         end if

         call MAPL_StateGetPointer(export, temp2d, 'VCPT', _RC)
         if (associated(temp2d)) then
            temp2d = 0.0
            do k = 1, km
               temp2d = temp2d + vr(:, :, k) * tempxy(:, :, k) * delp(:, :, k)
            end do
            temp2d = temp2d * (CP / GRAV)
         end if

         ! Compute Energetics After Dycore
         tempxy = vars%pt * (1.0 + EPS * qv) ! Convert TH to THV

         call MAPL_StateGetPointer(export, temp3d, 'THV', _RC)
         if (associated(temp3d)) temp3d = tempxy

         if (do_energetics) then
            call Energetics(ur, vr, tempxy, vars%pe, delp, vars%pkz, phisxy, kenrg, penrg, tenrg)
            kedyn = (kenrg - kenrg0) / dt
            pedyn = (penrg - penrg0) / dt
            tedyn = (tenrg - tenrg0) / dt
            call MAPL_StateGetPointer(export, temp2d, 'KEDYN', _RC)
            if (associated(temp2d)) temp2d = kedyn
            call MAPL_StateGetPointer(export, temp2d, 'PEDYN', _RC)
            if (associated(temp2d)) temp2d = pedyn
            call MAPL_StateGetPointer(export, temp2d, 'TEDYN', _RC)
            if (associated(temp2d)) temp2d = tedyn
         end if

         ! Compute/Get Omega
         call getOmega(omaxyz)

         ! Fluxes: UKE & VKE
         call MAPL_StateGetPointer(export, tempu, 'UKE', _RC)
         call MAPL_StateGetPointer(export, tempv, 'VKE', _RC)

         if (associated(tempu) .or. associated(tempv)) then
            tmp3d = 0.5 * (ur**2 + vr**2)
         end if

         if (associated(tempu)) then
            tempu = 0.0
            do k = 1, km
               tempu = tempu + ur(:, :, k) * tmp3d(:, :, k) * delp(:, :, k)
            end do
            tempu = tempu / GRAV
         end if

         if (associated(tempv)) then
            tempv = 0.0
            do k = 1, km
               tempv = tempv + vr(:, :, k) * tmp3d(:, :, k) * delp(:, :, k)
            end do
            tempv = tempv / GRAV
         end if

         ! Fluxes: UQV & VQV
         call MAPL_StateGetPointer(export, temp2d, 'UQV', _RC)
         if (associated(temp2d)) then
            temp2d = 0.0
            do k = 1, km
               temp2d = temp2d + ur(:, :, k) * qv(:, :, k) * delp(:, :, k)
            end do
            temp2d = temp2d / GRAV
         end if

         call MAPL_StateGetPointer(export, temp2d, 'VQV', _RC)
         if (associated(temp2d)) then
            temp2d = 0.0
            do k = 1, km
               temp2d = temp2d + vr(:, :, k) * qv(:, :, k) * delp(:, :, k)
            end do
            temp2d = temp2d / GRAV
         end if

         ! Fluxes: UQL & VQL
         call MAPL_StateGetPointer(export, temp2d, 'UQL', _RC)
         if (associated(temp2d)) then
            temp2d = 0.0
            do n = 1, size(names)
               if (TRIM(names(n)) == 'QLCN' .or. &
                    TRIM(names(n)) == 'QLLS') then
                  do k = 1, km
                     if (self%vars%tracer(n)%is_r4) then
                        temp2d = temp2d + ur(:, :, k) * self%vars%tracer(n)%content_r4(:, :, k) * delp(:, :, k)
                     else
                        temp2d = temp2d + ur(:, :, k) * self%vars%tracer(n)%content(:, :, k) * delp(:, :, k)
                     end if
                  end do
               end if
            end do
            temp2d = temp2d / GRAV
         end if

         call MAPL_StateGetPointer(export, temp2d, 'VQL', _RC)
         if (associated(temp2d)) then
            temp2d = 0.0
            do n = 1, size(names)
               if (TRIM(names(n)) == 'QLCN' .or. &
                    TRIM(names(n)) == 'QLLS') then
                  do k = 1, km
                     if (self%vars%tracer(n)%is_r4) then
                        temp2d = temp2d + vr(:, :, k) * self%vars%tracer(n)%content_r4(:, :, k) * delp(:, :, k)
                     else
                        temp2d = temp2d + vr(:, :, k) * self%vars%tracer(n)%content(:, :, k) * delp(:, :, k)
                     end if
                  end do
               end if
            end do
            temp2d = temp2d / GRAV
         end if

         ! Fluxes: UQI & VQI
         call MAPL_StateGetPointer(export, temp2d, 'UQI', _RC)
         if (associated(temp2d)) then
            temp2d = 0.0
            do n = 1, size(names)
               if (TRIM(names(n)) == 'QICN' .or. &
                    TRIM(names(n)) == 'QILS') then
                  do k = 1, km
                     if (self%vars%tracer(n)%is_r4) then
                        temp2d = temp2d + ur(:, :, k) * self%vars%tracer(n)%content_r4(:, :, k) * delp(:, :, k)
                     else
                        temp2d = temp2d + ur(:, :, k) * self%vars%tracer(n)%content(:, :, k) * delp(:, :, k)
                     end if
                  end do
               end if
            end do
            temp2d = temp2d / GRAV
         end if

         call MAPL_StateGetPointer(export, temp2d, 'VQI', _RC)
         if (associated(temp2d)) then
            temp2d = 0.0
            do n = 1, size(names)
               if (TRIM(names(n)) == 'QICN' .or. &
                    TRIM(names(n)) == 'QILS') then
                  do k = 1, km
                     if (self%vars%tracer(n)%is_r4) then
                        temp2d = temp2d + vr(:, :, k) * self%vars%tracer(n)%content_r4(:, :, k) * delp(:, :, k)
                     else
                        temp2d = temp2d + vr(:, :, k) * self%vars%tracer(n)%content(:, :, k) * delp(:, :, k)
                     end if
                  end do
               end if
            end do
            temp2d = temp2d / GRAV
         end if

         ! Height related diagnostics
         zle(:, :, km + 1) = phisxy(:, :)
         do k = km, 1, -1
            zle(:, :, k) = zle(:, :, k + 1) + CP * tempxy(:, :, k) * (pkxy(:, :, k + 1) - pkxy(:, :, k))
         end do
         zle = zle / GRAV

         call MAPL_StateGetPointer(export, temp3d, 'ZLE', _RC)
         if (associated(temp3d)) temp3d = zle

         call MAPL_StateGetPointer(export, temp3d, 'ZL', _RC)
         if (associated(temp3d)) temp3d = 0.5 * (zle(:, :, :km) + zle(:, :, 2:))

         call MAPL_StateGetPointer(export, temp3d, 'S', _RC)
         if (associated(temp3d)) temp3d = temp3d + GRAV * (0.5 * (zle(:, :, :km) + zle(:, :, 2:)))

         ! Fluxes: UPHI & VPHI
         call MAPL_StateGetPointer(export, tempu, 'UPHI', _RC)
         call MAPL_StateGetPointer(export, tempv, 'VPHI', _RC)

         if (associated(tempu) .or. associated(tempv)) zl = 0.5 * (zle(:, :, :km) + zle(:, :, 2:))

         if (associated(tempu)) then
            tempu = 0.0
            do k = 1, km
               tempu = tempu + ur(:, :, k) * zl(:, :, k) * delp(:, :, k)
            end do
         end if

         if (associated(tempv)) then
            tempv = 0.0
            do k = 1, km
               tempv = tempv + vr(:, :, k) * zl(:, :, k) * delp(:, :, k)
            end do
         end if

         ! Fill Surface and Near-Surface Variables
         HGT_SURFACE = 50.0
         if (km == 72) HGT_SURFACE = 0.0
         call MAPL_GridCompGetResource(gc, "HGT_SURFACE", HGT_SURFACE, default=HGT_SURFACE, _RC)
         if (HGT_SURFACE > 0.0) then
            ! Near surface height for surface
            call MAPL_StateGetPointer(export, temp2d, 'DZ', _RC)
            if (associated(temp2d)) temp2d = HGT_SURFACE

            ! Get the height above the surface
            do k = 1, km + 1
               zle(:, :, k) = zle(:, :, k) - zle(:, :, km + 1)
            end do

            call MAPL_StateGetPointer(export, temp2d, 'PS', _RC)
            if (associated(temp2d)) temp2d = vars%pe(:, :, km + 1)

            call MAPL_StateGetPointer(export, temp2d, 'US', _RC)
            if (associated(temp2d)) then
               call VertInterp(temp2d, ur, -zle, -HGT_SURFACE, _RC)
            end if

            call MAPL_StateGetPointer(export, temp2d, 'VS', _RC)
            if (associated(temp2d)) then
               call VertInterp(temp2d, vr, -zle, -HGT_SURFACE, _RC)
            end if

            call MAPL_StateGetPointer(export, temp2d, 'TA', _RC)
            if (associated(temp2d)) then
               tempxy = vars%pt * vars%pkz
               call VertInterp(temp2d, tempxy, -zle, -HGT_SURFACE, _RC)
            end if

            call MAPL_StateGetPointer(export, temp2d, 'QA', _RC)
            if (associated(temp2d)) then
               call VertInterp(temp2d, qv, -zle, -HGT_SURFACE, positive_definite=.true., _RC)
            end if

            call MAPL_StateGetPointer(export, temp2d, 'SPEED', _RC)
            if (associated(temp2d)) then
               call VertInterp(temp2d, sqrt(ur**2 + vr**2), -zle, -HGT_SURFACE, _RC)
            end if
         else
            ! Fill Surface with Lowest Model Level Variables
            call MAPL_StateGetPointer(export, temp2d, 'DZ', _RC)
            if (associated(temp2d)) temp2d = 0.5 * (zle(:, :, km) - zle(:, :, km + 1))

            call MAPL_StateGetPointer(export, temp2d, 'PS', _RC)
            if (associated(temp2d)) temp2d = vars%pe(:, :, km + 1)

            call MAPL_StateGetPointer(export, temp2d, 'US', _RC)
            if (associated(temp2d)) temp2d = ur(:, :, km)

            call MAPL_StateGetPointer(export, temp2d, 'VS', _RC)
            if (associated(temp2d)) temp2d = vr(:, :, km)

            call MAPL_StateGetPointer(export, temp2d, 'TA', _RC)
            if (associated(temp2d)) then
               tempxy = vars%pt * vars%pkz
               temp2d = tempxy(:, :, km)
            end if

            call MAPL_StateGetPointer(export, temp2d, 'QA', _RC)
            if (associated(temp2d)) temp2d = qv(:, :, km)

            call MAPL_StateGetPointer(export, temp2d, 'SPEED', _RC)
            if (associated(temp2d)) temp2d = sqrt(ur(:, :, km)**2 + vr(:, :, km)**2)
         end if

         call MAPL_StateGetPointer(export, temp2d, 'WSPD_10M', _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, sqrt(ur**2 + vr**2), -zle, -10.0, _RC)
         end if

         if (.not. HYDROSTATIC) then
            call MAPL_StateGetPointer(export, temp2d, 'VVEL_UP_100_1000', _RC)
            if (associated(temp2d)) then
               temp2d = vars%w(ifirstxy:ilastxy, jfirstxy:jlastxy, km)
               do k = km - 1, 1, -1
                  do j = jfirstxy, jlastxy
                     do i = ifirstxy, ilastxy
                        if ((vars%w(i, j, k) > temp2d(i - ifirstxy + 1, j - jfirstxy + 1)) .and. &
                             (vars%pe(i, j, k) >= 10000.0)) then
                           temp2d(i - ifirstxy + 1, j - jfirstxy + 1) = vars%w(i, j, k)
                        end if
                     end do
                  end do
               end do
            end if
            call MAPL_StateGetPointer(export, temp2d, 'VVEL_DN_100_1000', _RC)
            if (associated(temp2d)) then
               temp2d = vars%w(ifirstxy:ilastxy, jfirstxy:jlastxy, km)
               do k = km - 1, 1, -1
                  do j = jfirstxy, jlastxy
                     do i = ifirstxy, ilastxy
                        if ((vars%w(i, j, k) < temp2d(i - ifirstxy + 1, j - jfirstxy + 1)) .and. &
                             (vars%pe(i, j, k) >= 10000.0)) then
                           temp2d(i - ifirstxy + 1, j - jfirstxy + 1) = vars%w(i, j, k)
                        end if
                     end do
                  end do
               end do
            end if
         end if

         ! Updraft Helicty Exports
         call MAPL_StateGetPointer(export, uh25, 'UH25', _RC)
         call MAPL_StateGetPointer(export, uh03, 'UH03', _RC)
         call MAPL_StateGetPointer(export, srh01, 'SRH01', _RC)
         call MAPL_StateGetPointer(export, srh03, 'SRH03', _RC)
         call MAPL_StateGetPointer(export, srh25, 'SRH25', _RC)
         ! Per WMP, this calculation is not useful if running hydrostatic
         if (associated(uh25) .or. associated(uh03) .or. &
              associated(srh01) .or. associated(srh03) .or. associated(srh25) .or. &
              associated(shr06)) then
            if (.not. HYDROSTATIC) then
               call fv_getUpdraftHelicity(uh25tmp, uh03tmp, srh01tmp, srh03tmp, srh25tmp, shr06tmp)
               if (associated(uh25)) uh25 = uh25tmp
               if (associated(uh03)) uh03 = uh03tmp
               if (associated(srh01)) srh01 = srh01tmp
               if (associated(srh03)) srh03 = srh03tmp
               if (associated(srh25)) srh25 = srh25tmp
               if (associated(shr06)) shr06 = shr06tmp
            end if
         end if

         ! Divergence Exports
         logpe = log(vars%pe)

         call MAPL_StateGetPointer(export, temp3d, 'DIVG', _RC)
         if (associated(temp3d)) temp3d = divg

         call MAPL_StateGetPointer(export, temp2d, 'DIVG200', _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, dble(divg), logpe, log(20000.), _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, 'DIVG500', _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, dble(divg), logpe, log(50000.), _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, 'DIVG700', _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, dble(divg), logpe, log(70000.), _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, 'DIVG850', _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, dble(divg), logpe, log(85000.), _RC)
         end if

         ! Vorticity Exports
         call MAPL_StateGetPointer(export, temp3d, 'VORT', _RC)
         if (associated(temp3d)) temp3d = vort

         call MAPL_StateGetPointer(export, temp2d, 'VORT200', _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, dble(vort), logpe, log(20000.), _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, 'VORT500', _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, dble(vort), logpe, log(50000.), _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, 'VORT700', _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, dble(vort), logpe, log(70000.), _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, 'VORT850', _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, dble(vort), logpe, log(85000.), _RC)
         end if

         ! Vertical Velocity Exports
         call FILLOUT3(export, 'OMEGA', omaxyz, _RC)

         call MAPL_StateGetPointer(export, temp2d, 'OMEGA850', _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, omaxyz, logpe, log(85000.), _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, 'OMEGA700', _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, omaxyz, logpe, log(70000.), _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, 'OMEGA500', _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, omaxyz, logpe, log(50000.), _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, 'OMEGA200', _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, omaxyz, logpe, log(20000.), _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, 'OMEGA10', _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, omaxyz, logpe, log(1000.), _RC)
         end if

         if (.not. HYDROSTATIC) then

            call FILLOUT3(export, 'DELZ', vars%dz(ifirstxy:ilastxy, jfirstxy:jlastxy, :), _RC)
            call FILLOUT3(export, 'W', vars%w(ifirstxy:ilastxy, jfirstxy:jlastxy, :), _RC)

            call MAPL_StateGetPointer(export, temp2d, 'W850', _RC)
            if (associated(temp2d)) then
               call VertInterp(temp2d, vars%w(ifirstxy:ilastxy, jfirstxy:jlastxy, :), logpe, log(85000.), _RC)
            end if

            call MAPL_StateGetPointer(export, temp2d, 'W500', _RC)
            if (associated(temp2d)) then
               call VertInterp(temp2d, vars%w(ifirstxy:ilastxy, jfirstxy:jlastxy, :), logpe, log(50000.), _RC)
            end if

            call MAPL_StateGetPointer(export, temp2d, 'W200', _RC)
            if (associated(temp2d)) then
               call VertInterp(temp2d, vars%w(ifirstxy:ilastxy, jfirstxy:jlastxy, :), logpe, log(20000.), _RC)
            end if

            call MAPL_StateGetPointer(export, temp2d, 'W10', _RC)
            if (associated(temp2d)) then
               call VertInterp(temp2d, vars%w(ifirstxy:ilastxy, jfirstxy:jlastxy, :), logpe, log(1000.), _RC)
            end if
         end if

      end if ! SW_DYNAMICS

      call MAPL_GridCompTimerStop(gc, "DYN_EPILOGUE", _RC)

      ! De-Allocate Arrays

      if (do_energetics) then
         deallocate(kedyn)
         deallocate(pedyn)
         deallocate(tedyn)
         deallocate(kenrg)
         deallocate(penrg)
         deallocate(tenrg)
         deallocate(kenrg0)
         deallocate(penrg0)
         deallocate(tenrg0)
      end if

      deallocate(qsum1)
      deallocate(qsum2)

      deallocate(zl)
      deallocate(zle)
      deallocate(logpe)
      deallocate(plk)
      deallocate(pkxy)
      deallocate(vort)
      deallocate(divg)
      deallocate(tmp3d)
      deallocate(omaxyz)
      deallocate(epvxyz)
      deallocate(cxxyz)
      deallocate(cyxyz)
      deallocate(mfxxyz)
      deallocate(mfyxyz)
      deallocate(mfzxyz)
      deallocate(tempxy)
      deallocate(pe0)
      deallocate(pe1)
      deallocate(pl)
      deallocate(ua)
      deallocate(va)
      deallocate(uc)
      deallocate(vc)
      deallocate(uc0)
      deallocate(vc0)
      deallocate(ur)
      deallocate(vr)
      deallocate(qv)
      deallocate(ql)
      deallocate(qi)
      deallocate(qr)
      deallocate(qs)
      deallocate(qg)
      deallocate(ox)
      deallocate(delp)
      deallocate(dmdt)
      deallocate(dudt)
      deallocate(dvdt)
      deallocate(dtdt)
      deallocate(dqdt)
      deallocate(dthdt)
      deallocate(dpedt)
      deallocate(ddpdt)
      deallocate(phisxy)
      if (allocated(names)) deallocate(names)
      if (allocated(names0)) deallocate(names0)

      call freeTracers(self)

      !if (ADIABATIC) then
      !  ! Fill Exports
      !   call RunAddIncs(gc, import, export, clock, rc)
      !endif

      _RETURN(_SUCCESS)

   end subroutine Run

   subroutine PULL_Q(self, import, qqq, iNXQ, InFieldName, rc)
      type(DynState) :: self
      type(ESMF_State) :: import
      type(DynTracers) :: qqq ! Specific Humidity
      integer, intent(in) :: iNXQ
      character(len=*), optional, intent(in) :: InFieldName
      integer, optional, intent(out) :: rc

      character(len=ESMF_MAXSTR) :: fieldname, QFieldName
      type(ESMF_FieldBundle) :: bundle
      type(ESMF_Field) :: field
      type(ESMF_Array) :: array
      type(ESMF_TypeKind_Flag) :: kind
      real(kind=r4), pointer, contiguous :: ptr_r4(:, :, :)
      real(kind=r8), pointer :: ptr_r8(:, :, :)
      integer :: n, nq
      integer :: i1, in, j1, jn, im, jm, km
      integer :: status

      QFieldName = "Q"
      if (present(InFieldName)) QFieldName = InFieldName

      i1 = self%grid%is
      in = self%grid%ie
      j1 = self%grid%js
      jn = self%grid%je
      im = self%grid%npx
      jm = self%grid%npy
      km = self%grid%npz

      bundle = BundleAdv

      ! Count the friendlies
      call ESMF_FieldBundleGet(bundle, fieldCount=nq, _RC)

      nq = nq + iNXQ
      self%grid%nq = nq ! grid%NQ is now the "official" NQ

      ! Tracer pointer array
      if (associated(self%vars%tracer)) then
         call freeTracers(self)
      end if

      allocate(self%vars%tracer(nq), _STAT)

      do n = 1, nq - iNXQ
         call ESMF_FieldBundleGet(bundle, fieldIndex=n, field=field, _RC)
         call ESMF_FieldGet(field, array=array, name=fieldname, _RC)
         call ESMF_ArrayGet(array, typekind=kind, _RC)
         self%vars%tracer(n)%is_r4 = (kind == ESMF_TYPEKIND_R4) ! Is real*4?
         self%vars%tracer(n)%TNAME = fieldname
         if (self%vars%tracer(n)%is_r4) then
            call ESMF_ArrayGet(array, localDE=0, farrayptr=ptr_r4, _RC)
            self%vars%tracer(n)%content_r4(i1:in, j1:jn, 1:km) => ptr_r4
            if (fieldname == QFieldName) then
               qqq%is_r4 = .true.
               qqq%content_r4 => self%vars%tracer(n)%content_r4
            end if
         else
            call ESMF_ArrayGet(array, localDE=0, farrayptr=ptr_r8, _RC)
            self%vars%tracer(n)%content => ptr_r8
            if (fieldname == QFieldName) then
               qqq%is_r4 = .false.
               qqq%content => self%vars%tracer(n)%content
            end if
         end if
      end do

      _RETURN(_SUCCESS)
      _UNUSED_DUMMY(import)
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

      type(DynState), pointer :: self
      type(DynGrid), pointer :: grid
      type(DynVars), pointer :: vars
      type(DynTracers) :: qqq ! Specific Humidity
      type(ESMF_Grid) :: esmfgrid

      real(kind=r8), allocatable :: penrg(:, :) ! Vertically Integrated Cp*T
      real(kind=r8), allocatable :: kenrg(:, :) ! Vertically Integrated K
      real(kind=r8), allocatable :: tenrg(:, :) ! PHIS*(Psurf-Ptop)
      real(kind=r8), allocatable :: penrg0(:, :) ! Vertically Integrated Cp*T
      real(kind=r8), allocatable :: kenrg0(:, :) ! Vertically Integrated K
      real(kind=r8), allocatable :: tenrg0(:, :) ! PHIS*(Psurf-Ptop)

      real(kind=r8), pointer :: phisxy(:, :)
      real(kind=r4), pointer :: phis(:, :)
      real(kind=r8), allocatable :: slp(:, :)
      real(kind=r8), allocatable :: H1000(:, :)
      real(kind=r8), allocatable :: H850(:, :)
      real(kind=r8), allocatable :: H500(:, :)
      real(kind=r8), allocatable :: tmp3d(:, :, :)
      real(kind=r8), allocatable :: plk(:, :, :)
      real(kind=r8), allocatable :: pke(:, :, :)
      real(kind=r8), allocatable :: pl(:, :, :)
      real(kind=r8), allocatable :: ua(:, :, :)
      real(kind=r8), allocatable :: va(:, :, :)
      real(kind=r8), allocatable :: uc(:, :, :)
      real(kind=r8), allocatable :: vc(:, :, :)
      real(kind=r8), allocatable :: ur(:, :, :)
      real(kind=r8), allocatable :: vr(:, :, :)
      real(kind=r8), allocatable :: qv(:, :, :)
      real(kind=r8), allocatable :: dp(:, :, :)
      real(kind=r8), allocatable :: thv(:, :, :)
      real(kind=r8), allocatable :: zle(:, :, :)
      real(kind=r8), allocatable :: tempxy(:, :, :)
      real(kind=r8) :: tmax, tmin

      real(kind=r8), allocatable :: logpl(:, :, :)
      real(kind=r8), allocatable :: logpe(:, :, :)
      real(kind=r8), allocatable :: logps(:, :)

      real(kind=r4), pointer :: qold(:, :, :)
      real(kind=r4), pointer :: temp3d(:, :, :)
      real(kind=r4), pointer :: temp2d(:, :)

      real(kind=r4), pointer :: ztemp1(:, :)
      real(kind=r4), pointer :: ztemp2(:, :)
      real(kind=r4), pointer :: ztemp3(:, :)

      real(kind=r4), allocatable :: dthdtphyint1(:, :)
      real(kind=r4), allocatable :: dthdtphyint2(:, :)

      real(kind=FVPRC) :: dt
      logical :: do_energetics
      integer :: ifirstxy, ilastxy, jfirstxy, jlastxy
      integer :: im, jm, km, iNXQ
      integer :: i, j, k, status
      class(logger_t), pointer :: logger

      call MAPL_GridCompGet(gc, grid=esmfgrid, logger=logger, _RC)

      ! Retrieve the pointer to the internal state
      _GET_NAMED_PRIVATE_STATE(gc, DynState, PRIVATE_STATE, self)

      vars => self%vars ! direct handle to control variables
      grid => self%grid ! direct handle to grid
      dt = self%dt ! dynamics time step (large)

      ifirstxy = grid%is
      ilastxy = grid%ie
      jfirstxy = grid%js
      jlastxy = grid%je

      im = grid%npx
      jm = grid%npy
      km = grid%npz
      iNXQ = 0

      if (.not. SW_DYNAMICS) then

         allocate(dthdtphyint1(ifirstxy:ilastxy, jfirstxy:jlastxy))
         allocate(dthdtphyint2(ifirstxy:ilastxy, jfirstxy:jlastxy))

         do_energetics = .false.
         call MAPL_StateGetPointer(export, temp2d, "KE", _RC)
         if (associated(temp2d)) do_energetics = .true.
         call MAPL_StateGetPointer(export, temp2d, "KEPHY", _RC)
         if (associated(temp2d)) do_energetics = .true.
         call MAPL_StateGetPointer(export, temp2d, "PEPHY", _RC)
         if (associated(temp2d)) do_energetics = .true.
         call MAPL_StateGetPointer(export, temp2d, "TEPHY", _RC)
         if (associated(temp2d)) do_energetics = .true.
         if (do_energetics) then
            allocate(kenrg(ifirstxy:ilastxy, jfirstxy:jlastxy))
            allocate(penrg(ifirstxy:ilastxy, jfirstxy:jlastxy))
            allocate(tenrg(ifirstxy:ilastxy, jfirstxy:jlastxy))
            allocate(kenrg0(ifirstxy:ilastxy, jfirstxy:jlastxy))
            allocate(penrg0(ifirstxy:ilastxy, jfirstxy:jlastxy))
            allocate(tenrg0(ifirstxy:ilastxy, jfirstxy:jlastxy))
         end if

         allocate(tmp3d(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
         allocate(phisxy(ifirstxy:ilastxy, jfirstxy:jlastxy))
         allocate(logps(ifirstxy:ilastxy, jfirstxy:jlastxy))

         allocate(ua(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
         allocate(va(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
         allocate(uc(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
         allocate(vc(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
         allocate(ur(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
         allocate(vr(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
         allocate(qv(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
         allocate(pl(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
         allocate(logpl(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
         allocate(dp(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
         allocate(thv(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
         allocate(tempxy(ifirstxy:ilastxy, jfirstxy:jlastxy, km))

         allocate(plk(ifirstxy:ilastxy, jfirstxy:jlastxy, km))
         allocate(pke(ifirstxy:ilastxy, jfirstxy:jlastxy, km + 1))
         allocate(logpe(ifirstxy:ilastxy, jfirstxy:jlastxy, km + 1))
         allocate(zle(ifirstxy:ilastxy, jfirstxy:jlastxy, km + 1))

         call MAPL_StateGetPointer(import, phis, "PHIS", _RC)

         phisxy = real(phis, kind=r8)

         ! Compute Pressure Thickness
         dp = (vars%pe(:, :, 2:) - vars%pe(:, :, :km))

         ! Load Specific Humidity
         call MAPL_StateGetPointer(export, qold, 'Q', _RC)

         call PULL_Q(self, import, qqq, iNXQ, _RC)
         if ((.not. ADIABATIC) .and. (self%grid%nq > 0)) then
            if ((qqq%is_r4) .and. (associated(qqq%content_r4))) then
               if (size(qv) == size(qqq%content_r4)) qv = qqq%content_r4
            elseif (associated(qqq%content)) then
               if (size(qv) == size(qqq%content)) qv = qqq%content
            end if
            _ASSERT(all(qv >= 0.0), "RunAddIncs: negative or nan water vapor detected")
         else
            qv = 0.0
         end if

         ! Compute Energetics Before Diabatic Forcing
         if (associated(qold)) then
            thv = vars%pt * (1.0 + EPS * qold)
         else
            thv = vars%pt
         end if

         if (do_energetics) then
            call getAllWinds(vars%u, vars%v, ua=ua, va=va, uc=uc, vc=vc, ur=ur, vr=vr)
            call Energetics(ur, vr, thv, vars%pe, dp, vars%pkz, phisxy, kenrg0, penrg0, tenrg0)
         end if

         ! DTHVDTPHYINT
         call MAPL_StateGetPointer(export, temp2d, "DTHVDTPHYINT", _RC)
         if (associated(temp2d)) then
            dthdtphyint1 = 0.0
            do k = 1, km
               dthdtphyint1 = dthdtphyint1 + thv(:, :, k) * dp(:, :, k)
            end do
         end if

         ! Add Diabatic Forcing to State Variables
         call MAPL_GridCompTimerStart(gc, "PHYS_ADD_INCS", _RC)
         call ADD_INCS(esmfgrid, self, import, dt)
         call MAPL_GridCompTimerStop(gc, "PHYS_ADD_INCS", _RC)

         if (DYN_DEBUG) call DEBUG_FV_STATE('PHYSICS ADD_INCS', self)

         ! Update Mid-Layer Pressure and Pressure Thickness
         dp = (vars%pe(:, :, 2:) - vars%pe(:, :, :km))
         pl = (vars%pe(:, :, 2:) + vars%pe(:, :, :km)) * 0.5

         logpl = log(pl)
         logpe = log(vars%pe)
         logps = log(vars%pe(:, :, km + 1))

         ! Get Cubed-Sphere Wind Exports
         call getAllWinds(vars%u, vars%v, ua=ua, va=va, uc=uc, vc=vc, ur=ur, vr=vr)
         call FILLOUT3(export, 'U_DGRID', vars%u, _RC)
         call FILLOUT3(export, 'V_DGRID', vars%v, _RC)
         call FILLOUT3(export, 'U_CGRID', uc, _RC)
         call FILLOUT3(export, 'V_CGRID', vc, _RC)
         call FILLOUT3(export, 'U_AGRID', ua, _RC)
         call FILLOUT3(export, 'V_AGRID', va, _RC)

         ! Compute Energetics After Diabatic Forcing
         thv = vars%pt * (1.0 + EPS * qv)

#if defined(DEBUG_VPT)
         call Write_Profile(grid, thv, 'VPT')
#endif

         if (do_energetics) then
            call Energetics(ur, vr, thv, vars%pe, dp, vars%pkz, phisxy, kenrg, penrg, tenrg)
            call MAPL_StateGetPointer(export, temp2d, "KE", _RC)
            if (associated(temp2d)) temp2d = kenrg
            kenrg = (kenrg - kenrg0) / dt
            penrg = (penrg - penrg0) / dt
            tenrg = (tenrg - tenrg0) / dt
            call FILLOUT2(export, "KEPHY", kenrg, _RC)
            call FILLOUT2(export, "PEPHY", penrg, _RC)
            call FILLOUT2(export, "TEPHY", tenrg, _RC)
         end if

         ! DTHVDTPHYINT
         call MAPL_StateGetPointer(export, temp2d, "DTHVDTPHYINT", _RC)
         if (associated(temp2d)) then
            dthdtphyint2 = 0.0
            do k = 1, km
               dthdtphyint2 = dthdtphyint2 + thv(:, :, k) * dp(:, :, k)
            end do
            temp2d = (dthdtphyint2 - dthdtphyint1) * MAPL_P00**MAPL_KAPPA / (MAPL_GRAV * dt)
         end if

         plk = exp(KAPPA * log(0.5 * (vars%pe(:, :, 1:km) + vars%pe(:, :, 2:km + 1))))
         pke = exp(KAPPA * log(vars%pe))

         tempxy = vars%pt * vars%pkz ! Dry Temperature

#if defined(DEBUG_T)
         call Write_Profile(grid, tempxy, 'T')
#endif

         ! if (DEBUG_DYN) then
         !    block
         !       type(ESMF_VM) :: vm
         !       integer :: comm
         !       real :: maxmin(2)

         !       call ESMF_VMGetCurrent(vm, _RC)
         !       call ESMF_VMGet(vm, mpiCommunicator=comm, _RC)
         !       maxmin = MAPL_MaxMin(qv, comm, _RC)
         !       call logger%info("max/min(Q_AF_INC): %f/%f", maxmin(1), maxmin(2))
         !       maxmin = MAPL_MaxMin(tempxy, comm, _RC)
         !       call logger%info("max/min(T_AF_INC): %f/%f", maxmin(1), maxmin(2))
         !       maxmin = MAPL_MaxMin(ua, comm, _RC)
         !       call logger%info("max/min(U_AF_INC): %f/%f", maxmin(1), maxmin(2))
         !       maxmin = MAPL_MaxMin(va, comm, _RC)
         !       call logger%info("max/min(V_AF_INC): %f/%f", maxmin(1), maxmin(2))
         !    end block
         ! end if

         call FILLOUT3(export, "DELP", dp, _RC)
         call FILLOUT3(export, "U", ur, _RC)
         call FILLOUT3(export, "V", vr, _RC)
         call FILLOUT3(export, "T", tempxy, _RC)
         call FILLOUT3(export, "Q", qv, _RC)
         call FILLOUT3(export, "PL", pl, _RC)
         call FILLOUT3(export, "PLK", plk, _RC)
         call FILLOUT3(export, "PKE", pke, _RC)
         call FILLOUT3(export, "THV", thv, _RC)
         call FILLOUT3(export, "PT", vars%pt, _RC)
         call FILLOUT3(export, "PE", vars%pe, _RC)

         call MAPL_StateGetPointer(export, temp3d, "TH", _RC)
         if (associated(temp3d)) temp3d = (tempxy) &
              * (P00 / (0.5 * (vars%pe(:, :, 1:km) + vars%pe(:, :, 2:km + 1))))**KAPPA

#ifdef SKIP_TRACERS
         do itracer = 1, ntracers
            write(myTracer, "('Q',i5.5)") itracer - 1
            call MAPL_StateGetPointer(export, temp3d, TRIM(myTracer), _RC)
            if ((associated(temp3d)) .and. (self%grid%nq >= itracer)) then
               if (self%vars%tracer(itracer)%is_r4) then
                  temp3d = self%vars%tracer(itracer)%content_r4
               else
                  temp3d = self%vars%tracer(itracer)%content
               end if
            end if
         end do
#endif

         ! Compute Edge Heights
         zle(:, :, km + 1) = phisxy(:, :)
         do k = km, 1, -1
            zle(:, :, k) = zle(:, :, k + 1) + CP * thv(:, :, k) * (pke(:, :, k + 1) - pke(:, :, k))
         end do
         zle(:, :, :) = zle(:, :, :) / GRAV

         call FILLOUT3(export, "ZLE", zle, _RC)

         ! Compute Mid-Layer Heights
         call MAPL_StateGetPointer(export, temp3d, "ZL", _RC)
         if (associated(temp3d)) temp3d = 0.5 * (zle(:, :, 2:) + zle(:, :, :km))

         ! Fill Single Level Variables
         call MAPL_StateGetPointer(export, temp2d, "Z700", _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, zle * GRAV, logpe, log(70000.), _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, 'Z500', _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, zle * GRAV, logpe, log(50000.), _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, 'Z300', _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, zle * GRAV, logpe, log(30000.), _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, "H100", _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, zle, logpe, log(10000.), _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, "H200", _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, zle, logpe, log(20000.), _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, "H250", _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, zle, logpe, log(25000.), _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, "H300", _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, zle, logpe, log(30000.), _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, "H500", _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, zle, logpe, log(50000.), _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, "H700", _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, zle, logpe, log(70000.), _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, "H850", _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, zle, logpe, log(85000.), _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, "H1000", _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, zle, logpe, log(100000.), _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, "U50M", _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, ur, -zle, -50., _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, "V50M", _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, vr, -zle, -50., _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, "U100", _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, ur, logpe, log(10000.), _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, "U200", _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, ur, logpe, log(20000.), _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, "U250", _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, ur, logpe, log(25000.), _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, "U300", _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, ur, logpe, log(30000.), _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, "U500", _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, ur, logpe, log(50000.), _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, "U700", _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, ur, logpe, log(70000.), _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, "U850", _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, ur, logpe, log(85000.), _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, "V100", _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, vr, logpe, log(10000.), _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, "V200", _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, vr, logpe, log(20000.), _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, "V250", _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, vr, logpe, log(25000.), _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, "V300", _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, vr, logpe, log(30000.), _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, "V500", _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, vr, logpe, log(50000.), _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, "V700", _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, vr, logpe, log(70000.), _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, "V850", _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, vr, logpe, log(85000.), _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, "T100", _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, tempxy, logpe, log(10000.), _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, "T200", _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, tempxy, logpe, log(20000.), _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, "T250", _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, tempxy, logpe, log(25000.), _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, "T300", _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, tempxy, logpe, log(30000.), _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, "T500", _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, tempxy, logpe, log(50000.), _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, "T700", _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, tempxy, logpe, log(70000.), _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, "T850", _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, tempxy, logpe, log(85000.), _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, "Q100", _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, qv, logpe, log(10000.), positive_definite=.true., _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, "Q200", _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, qv, logpe, log(20000.), positive_definite=.true., _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, "Q250", _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, qv, logpe, log(25000.), positive_definite=.true., _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, "Q300", _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, qv, logpe, log(30000.), positive_definite=.true., _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, "Q500", _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, qv, logpe, log(50000.), positive_definite=.true., _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, "Q700", _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, qv, logpe, log(70000.), positive_definite=.true., _RC)
         end if

         call MAPL_StateGetPointer(export, temp2d, "Q850", _RC)
         if (associated(temp2d)) then
            call VertInterp(temp2d, qv, logpe, log(85000.), positive_definite=.true., _RC)
         end if

         ! Fill Model Top Level Variables
         call MAPL_StateGetPointer(export, temp2d, "UTOP", _RC)
         if (associated(temp2d)) temp2d = ur(:, :, 1)

         call MAPL_StateGetPointer(export, temp2d, "VTOP", _RC)
         if (associated(temp2d)) temp2d = vr(:, :, 1)

         call MAPL_StateGetPointer(export, temp2d, "TTOP", _RC)
         if (associated(temp2d)) temp2d = tempxy(:, :, 1)

         call MAPL_StateGetPointer(export, temp2d, "DELPTOP", _RC)
         if (associated(temp2d)) temp2d = dp(:, :, 1)

         ! Compute Surface Pressure
         call MAPL_StateGetPointer(export, temp2d, "PS", _RC)
         if (associated(temp2d)) temp2d = vars%pe(:, :, km + 1)

         ! Get the height above the surface
         do k = 1, km + 1
            zle(:, :, k) = zle(:, :, k) - zle(:, :, km + 1)
         end do

         call MAPL_StateGetPointer(export, temp3d, "ZLE0", _RC)
         if (associated(temp3d)) temp3d = zle

         call MAPL_StateGetPointer(export, temp3d, "ZL0", _RC)
         if (associated(temp3d)) temp3d = 0.5 * (zle(:, :, :km) + zle(:, :, 2:))

         ! Compute Vertically Averaged T,U
         call MAPL_StateGetPointer(export, temp2d, "TAVE", _RC)
         if (associated(temp2d)) then
            temp2d = 0.0
            do k = 1, km
               temp2d = temp2d + tempxy(:, :, k) * dp(:, :, k)
            end do
            temp2d = temp2d / (vars%pe(:, :, km + 1) - vars%pe(:, :, 1))
         end if

         call MAPL_StateGetPointer(export, temp2d, "UAVE", _RC)
         if (associated(temp2d)) then
            temp2d = 0.0
            do k = 1, km
               temp2d = temp2d + ur(:, :, k) * dp(:, :, k)
            end do
            temp2d = temp2d / (vars%pe(:, :, km + 1) - vars%pe(:, :, 1))
         end if

         ! Convert T to Tv
         tempxy = tempxy * (1.0 + EPS * qv)

         call MAPL_StateGetPointer(export, temp3d, "TV", _RC)
         if (associated(temp3d)) temp3d = tempxy

         ! Compute Sea-Level Pressure
         call MAPL_StateGetPointer(export, temp2d, "SLP", _RC)
         call MAPL_StateGetPointer(export, ztemp1, "H1000", _RC)
         call MAPL_StateGetPointer(export, ztemp2, "H850", _RC)
         call MAPL_StateGetPointer(export, ztemp3, "H500", _RC)

         if (associated(temp2d) &
              .or. associated(ztemp1) &
              .or. associated(ztemp2) &
              .or. associated(ztemp3)) then
            allocate(slp(ifirstxy:ilastxy, jfirstxy:jlastxy))
            allocate(H1000(ifirstxy:ilastxy, jfirstxy:jlastxy))
            allocate(H850(ifirstxy:ilastxy, jfirstxy:jlastxy))
            allocate(H500(ifirstxy:ilastxy, jfirstxy:jlastxy))
            do j = jfirstxy, jlastxy
               do i = ifirstxy, ilastxy
                  call get_slp(km, vars%pe(i, j, km + 1), phisxy(i, j), slp(i, j), &
                       vars%pe(i, j, 1:km + 1), &
                       vars%pkz(i, j, 1:km), &
                       tempxy(i, j, 1:km), &
                       H1000(i, j), H850(i, j), H500(i, j))
               end do
            end do

            !#define DEBUG_SLP
#if defined(DEBUG_SLP)
            call Write_Profile(grid, slp / 100.0, 'SLP')
#endif

            if (associated(temp2d)) temp2d = slp
            if (associated(ztemp1)) where (ztemp1 == MAPL_UNDEFINED_REAL) ztemp1 = H1000
            if (associated(ztemp2)) where (ztemp2 == MAPL_UNDEFINED_REAL) ztemp2 = H850
            if (associated(ztemp3)) where (ztemp3 == MAPL_UNDEFINED_REAL) ztemp3 = H500
            deallocate(slp, H1000, H850, H500)
         end if

         ! Deallocate Memory
         if (do_energetics) then
            deallocate(kenrg)
            deallocate(penrg)
            deallocate(tenrg)
            deallocate(kenrg0)
            deallocate(penrg0)
            deallocate(tenrg0)
         end if

         deallocate(tmp3d)

         deallocate(phisxy)

         deallocate(ua)
         deallocate(va)
         deallocate(uc)
         deallocate(vc)
         deallocate(ur)
         deallocate(vr)
         deallocate(qv)
         deallocate(pl)
         deallocate(dp)
         deallocate(tempxy)

         deallocate(thv)
         deallocate(plk)
         deallocate(pke)
         deallocate(logpl)
         deallocate(logpe)
         deallocate(logps)
         deallocate(zle)
         deallocate(dthdtphyint1)
         deallocate(dthdtphyint2)

         call freeTracers(self)

      end if ! .not. SW_DYNAMICS

      _RETURN(_SUCCESS)
      _UNUSED_DUMMY(clock)

   end subroutine RunAddIncs

   subroutine ADD_INCS(esmfgrid, self, import, dt, is_weighted, rc)

      use fms_mod, only: set_domain, nullify_domain
      use fv_diagnostics_mod, only: prt_maxmin
      use time_manager_mod, only: time_type
      use fv_update_phys_mod, only: fv_update_phys

      !INPUT PARAMETERS:
      type(ESMF_Grid), intent(in) :: esmfgrid
      type(DynState), pointer :: self
      type(ESMF_State), intent(inout) :: import
      real(kind=FVPRC), intent(in) :: dt
      logical, optional, intent(in) :: is_weighted
      integer, optional, intent(out) :: rc

      !DESCRIPTION:  This routine adds the tendencies to the state,
      !              weighted appropriately by the time step.  Temperature
      !              tendencies are pressure weighted (ie., DELP*DT/Dt).
      !              All tendencies are on the A-grid, and have an XY decomposition.

      logical :: is_weighted_

      integer :: II, JJ, i, j, k, L
      integer :: is, ie, js, je, km
      integer :: isd, ied, jsd, jed
      real(kind=r4), pointer :: tend(:, :, :)
      real(kind=r4), allocatable, dimension(:, :) :: lons, lats
      real(kind=r8), allocatable :: DPNEW(:, :, :), DPOLD(:, :, :)
      real(kind=r8), allocatable :: tend_ua(:, :, :), tend_va(:, :, :)
      real(kind=r8), allocatable :: tend_un(:, :, :), tend_vn(:, :, :)
      real(kind=FVPRC), allocatable :: Q(:, :, :, :), CVM(:, :, :)

      type(DynTracers) :: qqq ! Specific Humidity
      integer :: n, nwat_tracers, nwat, sphum, liq_wat, ice_wat, rainwat, snowwat, graupel
      real, parameter :: c_ice = 1972. !< heat capacity of ice at -15.C
      real, parameter :: c_liq = 4.1855e+3 !< GFS: heat capacity of water at 0C
      real, parameter :: c_vap = MAPL_CPVAP !< 1846.
      real, parameter :: c_air = MAPL_CP
      real(kind=FVPRC) :: fac
      integer :: status

      is_weighted_ = .true.
      if (present(is_weighted)) then
         is_weighted_ = is_weighted
      end if

      is = self%grid%is
      ie = self%grid%ie
      js = self%grid%js
      je = self%grid%je
      km = self%grid%npz

      isd = self%grid%isd
      ied = self%grid%ied
      jsd = self%grid%jsd
      jed = self%grid%jed

      ! call MAPL_Get( MAPL, LONS=LONS, LATS=LATS, RC=STATUS )
      ! VERIFY_(STATUS)
      call MAPL_GridGetCoordinates(esmfgrid, latitudes=lats, longitudes=lons, _RC)

      ! **********************************************************************
      ! ****  Use QV from FV3 init when coldstarting idealized cases      ****
      ! **********************************************************************

      ! Determine how many water species we have
      nwat = self%vars%nwat
      nwat_tracers = 0
      if ((nwat == 0) .and. (.not. ADIABATIC)) then
         do n = 1, self%grid%nq
            if (TRIM(self%vars%tracer(n)%TNAME) == 'Q') nwat_tracers = nwat_tracers + 1
            if (TRIM(self%vars%tracer(n)%TNAME) == 'QLCN') nwat_tracers = nwat_tracers + 1
            if (TRIM(self%vars%tracer(n)%TNAME) == 'QLLS') nwat_tracers = nwat_tracers + 1
            if (TRIM(self%vars%tracer(n)%TNAME) == 'QICN') nwat_tracers = nwat_tracers + 1
            if (TRIM(self%vars%tracer(n)%TNAME) == 'QILS') nwat_tracers = nwat_tracers + 1
         end do
         ! We must have these first 5 at a minimum
         _ASSERT(nwat_tracers == 5, 'expecting 5 water species: Q QLCN QLLS QICN QILS')
         ! Check for QRAIN, QSNOW, QGRAUPEL
         do n = 1, self%grid%nq
            if (TRIM(self%vars%tracer(n)%TNAME) == 'QRAIN') nwat_tracers = nwat_tracers + 1
            if (TRIM(self%vars%tracer(n)%TNAME) == 'QSNOW') nwat_tracers = nwat_tracers + 1
            if (TRIM(self%vars%tracer(n)%TNAME) == 'QGRAUPEL') nwat_tracers = nwat_tracers + 1
         end do
         if (nwat_tracers >= 5) nwat = 3 ! state has QV, QLIQ, QICE
         if (nwat_tracers == 8) nwat = 6 ! state has QV, QLIQ, QICE, QRAIN, QSNOW, QGRAUPEL
      end if
      if (.not. ADIABATIC) then
         _ASSERT(nwat >= 1, 'expecting water species (nwat) to match')
      end if

      select case (nwat)
      case (1)
         sphum = 1
         liq_wat = -1
         ice_wat = -1
         rainwat = -1
         snowwat = -1
         graupel = -1
      case (3)
         sphum = 1
         liq_wat = 2
         ice_wat = 3
         rainwat = -1
         snowwat = -1
         graupel = -1
      case (6:7)
         sphum = 1
         liq_wat = 2
         ice_wat = 3
         rainwat = 4
         snowwat = 5
         graupel = 6
      end select

      if (nwat >= 1) then
         allocate(Q(is:ie, js:je, 1:km, nwat))
         allocate(CVM(is:ie, js:je, 1:km))
         Q(:, :, :, :) = 0.0
         call PULL_Q(self, import, qqq, NXQ, InFieldName='Q', rc=rc)
         if (DYN_COLDSTART .and. overwrite_Q .and. (.not. ADIABATIC)) then
            ! USE Q computed by FV3
            call getQ(Q(:, :, :, sphum), 'Q')
            overwrite_Q = .false.
            call WRITE_PARALLEL("Using QV from FV3 Initial Conditions")
            fac = 1.0
            call prt_maxmin('AI Q', Q(:, :, :, sphum), is, ie, js, je, 0, km, fac)
            if ((qqq%is_r4) .and. (associated(qqq%content_r4))) then
               if (size(Q(:, :, :, sphum)) == size(qqq%content_r4)) qqq%content_r4 = Q(:, :, :, sphum)
            elseif (associated(qqq%content)) then
               if (size(Q(:, :, :, sphum)) == size(qqq%content)) qqq%content = Q(:, :, :, sphum)
            end if
         else
            ! Grab QV from imports
            if ((qqq%is_r4) .and. (associated(qqq%content_r4))) then
               if (size(Q(:, :, :, sphum)) == size(qqq%content_r4)) Q(:, :, :, sphum) = qqq%content_r4
            elseif (associated(qqq%content)) then
               if (size(Q(:, :, :, sphum)) == size(qqq%content)) Q(:, :, :, sphum) = qqq%content
            end if
         end if
      end if
      if (nwat >= 3) then
         ! Grab QLIQ from imports
         call PULL_Q(self, import, qqq, NXQ, InFieldName='QLLS', rc=rc)
         if ((qqq%is_r4) .and. (associated(qqq%content_r4))) then
            if (size(Q(:, :, :, liq_wat)) == size(qqq%content_r4)) Q(:, :, :, liq_wat) = Q(:, :, :, liq_wat) &
                 + qqq%content_r4
         elseif (associated(qqq%content)) then
            if (size(Q(:, :, :, liq_wat)) == size(qqq%content)) Q(:, :, :, liq_wat) = Q(:, :, :, liq_wat) + qqq%content
         end if
         call PULL_Q(self, import, qqq, NXQ, InFieldName='QLCN', rc=rc)
         if ((qqq%is_r4) .and. (associated(qqq%content_r4))) then
            if (size(Q(:, :, :, liq_wat)) == size(qqq%content_r4)) Q(:, :, :, liq_wat) = Q(:, :, :, liq_wat) &
                 + qqq%content_r4
         elseif (associated(qqq%content)) then
            if (size(Q(:, :, :, liq_wat)) == size(qqq%content)) Q(:, :, :, liq_wat) = Q(:, :, :, liq_wat) + qqq%content
         end if
         ! Grab QICE from imports
         call PULL_Q(self, import, qqq, NXQ, InFieldName='QILS', rc=rc)
         if ((qqq%is_r4) .and. (associated(qqq%content_r4))) then
            if (size(Q(:, :, :, ice_wat)) == size(qqq%content_r4)) Q(:, :, :, ice_wat) = Q(:, :, :, ice_wat) &
                 + qqq%content_r4
         elseif (associated(qqq%content)) then
            if (size(Q(:, :, :, ice_wat)) == size(qqq%content)) Q(:, :, :, ice_wat) = Q(:, :, :, ice_wat) + qqq%content
         end if
         call PULL_Q(self, import, qqq, NXQ, InFieldName='QICN', rc=rc)
         if ((qqq%is_r4) .and. (associated(qqq%content_r4))) then
            if (size(Q(:, :, :, ice_wat)) == size(qqq%content_r4)) Q(:, :, :, ice_wat) = Q(:, :, :, ice_wat) &
                 + qqq%content_r4
         elseif (associated(qqq%content)) then
            if (size(Q(:, :, :, ice_wat)) == size(qqq%content)) Q(:, :, :, ice_wat) = Q(:, :, :, ice_wat) + qqq%content
         end if
      end if
      if (nwat >= 6) then
         ! Grab RAIN from imports
         call PULL_Q(self, import, qqq, NXQ, InFieldName='QRAIN', rc=rc)
         if ((qqq%is_r4) .and. (associated(qqq%content_r4))) then
            if (size(Q(:, :, :, rainwat)) == size(qqq%content_r4)) Q(:, :, :, rainwat) = qqq%content_r4
         elseif (associated(qqq%content)) then
            if (size(Q(:, :, :, rainwat)) == size(qqq%content)) Q(:, :, :, rainwat) = qqq%content
         end if
         ! Grab SNOW from imports
         call PULL_Q(self, import, qqq, NXQ, InFieldName='QSNOW', rc=rc)
         if ((qqq%is_r4) .and. (associated(qqq%content_r4))) then
            if (size(Q(:, :, :, snowwat)) == size(qqq%content_r4)) Q(:, :, :, snowwat) = qqq%content_r4
         elseif (associated(qqq%content)) then
            if (size(Q(:, :, :, snowwat)) == size(qqq%content)) Q(:, :, :, snowwat) = qqq%content
         end if
         ! Grab GRAUPEL from imports
         call PULL_Q(self, import, qqq, NXQ, InFieldName='QGRAUPEL', rc=rc)
         if ((qqq%is_r4) .and. (associated(qqq%content_r4))) then
            if (size(Q(:, :, :, graupel)) == size(qqq%content_r4)) Q(:, :, :, graupel) = qqq%content_r4
         elseif (associated(qqq%content)) then
            if (size(Q(:, :, :, graupel)) == size(qqq%content)) Q(:, :, :, graupel) = qqq%content
         end if
      end if

      if ((.not. ADIABATIC) .and. (DO_ADD_INCS)) then

         ! **********************************************************************
         ! ****                      Wind Tendencies                         ****
         ! ****         Note: State Variables are on the D-Grid,             ****
         ! ****        while import Tendencies are on the A-Grid             ****
         ! **********************************************************************

         allocate(tend_ua(is:ie, js:je, km))
         allocate(tend_va(is:ie, js:je, km))
         allocate(tend_un(is:ie, js:je + 1, km))
         allocate(tend_vn(is:ie + 1, js:je, km))

         call MAPL_StateGetPointer(import, tend, "DUDT", _RC)
         tend_ua(is:ie, js:je, 1:km) = tend

         call MAPL_StateGetPointer(import, tend, "DVDT", _RC)
         tend_va(is:ie, js:je, 1:km) = tend

         !if (.not. HYDROSTATIC ) then
         !  call MAPL_StateGetPointer(import, TEND, 'DWDT', _RC)
         !  self%vars%W = self%vars%W + DT*TEND(is:ie,js:je,1:km)
         !endif

         ! Put the wind tendencies on the Native Dynamics grid
         call Agrid_To_Native(tend_ua, tend_va, tend_un, tend_vn, wind_increment_limiter=800.d0 / 86400.d0)

         ! Add the wind tendencies to the control variables
         self%vars%u = self%vars%u + dt * tend_un(is:ie, js:je, 1:km)
         self%vars%v = self%vars%v + dt * tend_vn(is:ie, js:je, 1:km)

         deallocate(tend_ua)
         deallocate(tend_va)

         ! **********************************************************************
         ! ****           Compute Old Pressure Thickness                     ****
         ! **********************************************************************

         allocate(DPOLD(is:ie, js:je, km))

         if (is_weighted_) then
            do k = 1, km
               DPOLD(:, :, k) = (self%vars%pe(:, :, k + 1) - self%vars%pe(:, :, k))
            end do
         else
            DPOLD = 1.0
         end if

         ! **********************************************************************
         ! ****                     Update Edge Pressures                    ****
         ! **********************************************************************

         call MAPL_StateGetPointer(import, tend, "DPEDT", _RC)
         self%vars%pe = self%vars%pe + dt * tend

         ! **********************************************************************
         ! ****           Compute New Pressure Thickness                     ****
         ! **********************************************************************

         allocate(DPNEW(is:ie, js:je, km))

         if (is_weighted_) then
            do k = 1, km
               DPNEW(:, :, k) = (self%vars%pe(:, :, k + 1) - self%vars%pe(:, :, k))
            end do
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

         call MAPL_StateGetPointer(import, tend, "DTDT", _RC)

         !if (DYN_DEBUG) then
         !   call prt_maxmin('AI PT1', self%vars%PT ,  is, ie, js, je, 0, km, 1.d00, MAPL_AM_I_ROOT())
         !endif

         select case (nwat)
         case (6:7)
            CVM = (1. - (Q(:, :, :, sphum) + Q(:, :, :, liq_wat) + Q(:, :, :, rainwat) + Q(:, :, :, ice_wat) +&
                 Q(:, :, :, snowwat) + Q(:, :, :, graupel))) * c_air + &
                 (Q(:, :, :, sphum)) * c_vap + &
                 (Q(:, :, :, liq_wat) + Q(:, :, :, rainwat)) * c_liq + &
                 (Q(:, :, :, ice_wat) + Q(:, :, :, snowwat) + Q(:, :, :, graupel)) * c_ice
         case (3)
            CVM = (1. - (Q(:, :, :, sphum) + Q(:, :, :, liq_wat) + Q(:, :, :, ice_wat))) * c_air + &
                 (Q(:, :, :, sphum)) * c_vap + &
                 (Q(:, :, :, liq_wat)) * c_liq + &
                 (Q(:, :, :, ice_wat)) * c_ice
         case default
            CVM = MAPL_CP
         end select

         ! Make previous PT into just T
         self%vars%pt = self%vars%pt * self%vars%pkz

         if (.not. HYDROSTATIC) then
            ! remove old T from DZ
            self%vars%dz = self%vars%dz / self%vars%pt

            ! Update T
            self%vars%pt = self%vars%pt * DPOLD
            self%vars%pt = (self%vars%pt + dt * tend * (MAPL_CP / CVM)) / DPNEW

            ! update DZ with new T
            self%vars%dz = self%vars%dz * self%vars%pt
         else
            ! Update T
            self%vars%pt = self%vars%pt * DPOLD
            self%vars%pt = (self%vars%pt + dt * tend * (MAPL_CP / CVM)) / DPNEW
         end if

         if (DEBUG_TQ_ERRORS) then
            do L = 1, km
               do j = js, je
                  do i = is, ie
                     if ((self%vars%pt(i, j, L) > 333.0) .or. (self%vars%pt(i, j, L)/=self%vars%pt(i, j, L)) .or. &
                          (Q(i, j, L, sphum) < 0.0) .or. (Q(i, j, L, sphum)/=Q(i, j, L, sphum)) .or. &
                          (Q(i, j, L, liq_wat) < 0.0) .or. (Q(i, j, L, liq_wat)/=Q(i, j, L, liq_wat)) .or. &
                          (Q(i, j, L, ice_wat) < 0.0) .or. (Q(i, j, L, ice_wat)/=Q(i, j, L, ice_wat)) .or. &
                          (Q(i, j, L, rainwat) < 0.0) .or. (Q(i, j, L, rainwat)/=Q(i, j, L, rainwat)) .or. &
                          (Q(i, j, L, snowwat) < 0.0) .or. (Q(i, j, L, snowwat)/=Q(i, j, L, snowwat)) .or. &
                          (Q(i, j, L, graupel) < 0.0) .or. (Q(i, j, L, graupel)/=Q(i, j, L, graupel))) then
                        print *, "T or Q  spike detected : ", self%vars%pt(i, j, L)
                        print *, "  Temp  ANA|PHY  Increment : ", (dt * tend(i, j, L) * (MAPL_CP / CVM(i, j, L))) / &
                             DPNEW(i, j, L)
                        print *, "    IN ADD_INCS inside DYN   "
                        II = i - is + 1
                        JJ = j - js + 1
                        print *, "  Latitude       =", lats(II, JJ) * 180.0 / MAPL_PI
                        print *, "  Longitude      =", lons(II, JJ) * 180.0 / MAPL_PI
                        print *, "  Pressure (mb)  =", 0.5 * (self%vars%pe(i, j, L + 1) + self%vars%pe(i, j, L)) / 100.0

                        print *, "  UWND =", self%vars%u(i, j, L), " UINC =", dt * tend_un(i, j, L)
                        print *, "  VWND =", self%vars%v(i, j, L), " VINC =", dt * tend_vn(i, j, L)
                        if (nwat >= 6) then
                           print *, "  QV=", Q(i, j, L, sphum), "  QL=", Q(i, j, L, liq_wat), "  QI=", Q(i, j, L, &
                                ice_wat)
                           print *, "  QR=", Q(i, j, L, rainwat), "  QS=", Q(i, j, L, snowwat), "  QG=", Q(i, j, L, &
                                graupel)
                        end if
                     end if
                  end do ! IM loop
               end do ! JM loop
            end do ! LM loop
         end if

         deallocate(tend_un)
         deallocate(tend_vn)

         ! Update PKZ from hydrostatic pressures
         !  This isn't entirely necessary, FV3 overwrites this in fv_dynamics
         !  but we have to get back to PT here
         !!   call getPKZ(self%vars%PKZ,self%vars%PT,Q,self%vars%PE,self%vars%DZ,HYDROSTATIC)
         call getPKZ(self%vars%pkz, self%vars%pe)

         ! Make T back into PT
         self%vars%pt = self%vars%pt / self%vars%pkz

         !if (DYN_DEBUG) then
         !call prt_maxmin('AI PT2', self%vars%PT ,  is, ie, js, je, 0, km, 1.d00, MAPL_AM_I_ROOT())
         !endif

         deallocate(DPNEW)
         deallocate(DPOLD)

      end if ! .not. Adiabatic

      if (allocated(Q)) deallocate(Q)
      if (allocated(CVM)) deallocate(CVM)

      _RETURN(_SUCCESS)

   end subroutine ADD_INCS

   subroutine FILLOUT3r8(export, name, v, rc)
      type(ESMF_State), intent(inout) :: export
      character(len=*), intent(in) :: name
      real(kind=r8), intent(in) :: v(:, :, :)
      integer, optional, intent(out) :: rc

      real(kind=r8), pointer :: cpl(:, :, :)
      integer :: status

      call MAPL_StateGetPointer(export, cpl, name, _RC)
      if (associated(cpl)) cpl = v

      _RETURN(_SUCCESS)
   end subroutine FILLOUT3r8

   subroutine FILLOUT3(export, name, v, rc)
      type(ESMF_State), intent(inout) :: export
      character(len=*), intent(in) :: name
      real(kind=r8), intent(in) :: v(:, :, :)
      integer, optional, intent(out) :: rc

      real(kind=r4), pointer :: cpl(:, :, :)
      integer :: status

      call MAPL_StateGetPointer(export, cpl, name, _RC)
      if (associated(cpl)) cpl = v

      _RETURN(_SUCCESS)
   end subroutine FILLOUT3

   subroutine FILLOUT2(export, name, v, rc)
      type(ESMF_State), intent(inout) :: export
      character(len=*), intent(in) :: name
      real(kind=r8), intent(in) :: v(:, :)
      integer, optional, intent(out) :: rc

      real(kind=r4), pointer :: cpl(:, :)
      integer :: status

      call MAPL_StateGetPointer(export, cpl, name, _RC)
      if (associated(cpl)) cpl = v

      _RETURN(_SUCCESS)
   end subroutine FILLOUT2

   subroutine Energetics(ua, va, thv, ple, delp, pk, phis, keint, peint, teint, ke, cpt, gze)

      real(kind=8), optional, intent(out) :: ke(:, :, :)
      real(kind=8), optional, intent(out) :: cpt(:, :, :)
      real(kind=8), optional, intent(out) :: gze(:, :, :)
      real(kind=8) :: ua(:, :, :)
      real(kind=8) :: va(:, :, :)
      real(kind=8) :: thv(:, :, :)
      real(kind=8) :: ple(:, :, :)
      real(kind=8) :: delp(:, :, :)
      real(kind=8) :: pk(:, :, :)
      real(kind=8) :: keint(:, :)
      real(kind=8) :: peint(:, :)
      real(kind=8) :: teint(:, :)
      real(kind=8) :: phis(:, :)

      real(kind=8) :: kinetic, potential
      integer :: i, ifirst, ilast
      integer :: j, jfirst, jlast
      integer :: km, k

      real(kind=8), allocatable :: pke(:, :, :)
      real(kind=8), allocatable :: phiT(:, :)

      ifirst = lbound(ua, 1)
      ilast = ubound(ua, 1)
      jfirst = lbound(ua, 2)
      jlast = ubound(ua, 2)
      km = ubound(ua, 3)

      allocate(pke(ifirst:ilast, jfirst:jlast, 1:km + 1))
      allocate(phiT(ifirst:ilast, jfirst:jlast))

      ! Compute Model Edge Heights
      pke = ple**KAPPA
      phiT = phis
      if (present(gze)) gze(:, :, km + 1) = phis
      do k = km, 1, -1
         phiT = phiT + CP * thv(:, :, k) * (pke(:, :, k + 1) - pke(:, :, k))
         if (present(gze)) gze(:, :, k) = phiT
      end do

      ! Compute Energetics:  Cp*Tv + K + PHI
      keint = 0.0
      peint = 0.0
      do k = 1, km
         do j = jfirst, jlast
            do i = ifirst, ilast
               kinetic = 0.5_r8 * (ua(i, j, k)**2 + va(i, j, k)**2)
               potential = CP * thv(i, j, k) * pk(i, j, k)
               keint(i, j) = keint(i, j) + kinetic * delp(i, j, k)
               peint(i, j) = peint(i, j) + potential * delp(i, j, k)
               if (present(ke)) ke(i, j, k) = kinetic
               if (present(cpt)) cpt(i, j, k) = potential
            end do
         end do
      end do
      keint(:, :) = keint(:, :) / GRAV
      peint(:, :) = peint(:, :) / GRAV
      teint(:, :) = (phis(:, :) * ple(:, :, km + 1) - phiT(:, :) * ple(:, :, 1)) / GRAV

      deallocate(pke)
      deallocate(phiT)

      return
   end subroutine Energetics

   !BOP
   !IROUTINE: Finalize

   !DESCRIPTION: Writes restarts and cleans-up through MAPL\_GenericFinalize and
   !   deallocates memory from the Private Internal state.

   !INTERFACE:
   subroutine Finalize(gc, import, export, clock, rc)
      !ARGUMENTS:
      type(ESMF_GridComp) :: gc
      type(ESMF_State) :: import
      type(ESMF_State) :: export
      type(ESMF_Clock) :: clock
      integer, intent(out) :: rc
      !EOP

      type(DynState), pointer :: self
      integer :: status

      ! Retrieve the pointer to the state
      _GET_NAMED_PRIVATE_STATE(gc, DynState, PRIVATE_STATE, self)

      call DynFinalize(self)

      _RETURN(_SUCCESS)
      _UNUSED_DUMMY(import)
      _UNUSED_DUMMY(export)
      _UNUSED_DUMMY(clock)
   end subroutine Finalize

   subroutine get_slp(km, ps, phis, slp, pe, pk, tv, H1000, H850, H500)
      integer :: km
      real(kind=r8) :: pk(km) ! layer-mean P**kappa
      real(kind=r8) :: tv(km) ! layer-mean virtual Temperature
      real(kind=r8) :: pe(km + 1) ! press at layer edges (Pa)
      real(kind=r8) :: ps ! surface pressure (Pa)
      real(kind=r8) :: phis ! surface geopotential
      real(kind=r8) :: slp ! sea-level pressure (hPa)
      real(kind=r8) :: H1000 ! 1000mb height
      real(kind=r8) :: H850 !  850mb height
      real(kind=r8) :: H500 !  500mb height
      real(kind=r8) :: tstar ! extrapolated temperature (K)
      real(kind=r8) :: p_bot
      real(kind=r8) :: tref ! Reference virtual temperature (K)
      real(kind=r8) :: pref ! Reference pressure level (Pa)
      real(kind=r8) :: pkref ! Reference pressure level (Pa) ** kappa
      real(kind=r8) :: dp1, dp2

      real(kind=r8), parameter :: gamma = 6.5e-3
      real(kind=r8), parameter :: p_offset = 15000.
      real(kind=r8), parameter :: gg = gamma / MAPL_GRAV

      real(kind=r8), parameter :: factor = MAPL_GRAV / (MAPL_RGAS * gamma)
      real(kind=r8), parameter :: yfactor = MAPL_RGAS * gg

      integer :: k_bot, k, k1, k2

      p_bot = ps - p_offset
      k_bot = -1

      do k = km, 2, -1
         if (pe(k + 1) < p_bot) then
            k_bot = k
            exit
         end if
      end do

      k1 = k_bot - 1
      k2 = k_bot
      dp1 = pe(k_bot) - pe(k_bot - 1)
      dp2 = pe(k_bot + 1) - pe(k_bot)
      pkref = (pk(k1) * dp1 + pk(k2) * dp2) / (dp1 + dp2)
      tref = (tv(k1) * dp1 + tv(k2) * dp2) / (dp1 + dp2)
      pref = 0.5 * (pe(k_bot + 1) + pe(k_bot - 1))
      tstar = tref * (ps / pref)**yfactor

      slp = ps * (1.0 + gg * phis / tstar)**factor
      H1000 = (phis / MAPL_GRAV) - (tstar / gamma) * ((100000.0 / ps)**(1. / factor) - 1.0)
      H850 = (phis / MAPL_GRAV) - (tstar / gamma) * ((85000.0 / ps)**(1. / factor) - 1.0)
      H500 = (phis / MAPL_GRAV) - (tstar / gamma) * ((50000.0 / ps)**(1. / factor) - 1.0)

      return
   end subroutine get_slp

   subroutine VertInterp(v2, v3, ple, pp, positive_definite, rc)
      real(kind=r4), intent(out) :: v2(:, :)
      real(kind=r8), intent(in) :: v3(:, :, :)
      real(kind=r8), intent(in) :: ple(:, :, :)
      real, intent(in) :: pp
      logical, optional, intent(in) :: positive_definite
      integer, optional, intent(out) :: rc

      real, dimension(size(v2, 1), size(v2, 2)) :: al, pt, PB
      integer :: k, km
      logical :: edge

      km = size(ple, 3) - 1
      edge = size(v3, 3) == km + 1

      _ASSERT(edge .or. size(v3, 3) == km, 'needs informative message')

      v2 = MAPL_UNDEFINED_REAL

      if (edge) then
         PB = ple(:, :, km + 1)
         do k = km, 1, -1
            pt = ple(:, :, k)
            if (all(PB < pp)) exit
            where (pp > pt .and. pp <= PB)
               al = (PB - pp) / (PB - pt)
               v2 = v3(:, :, k) * al + v3(:, :, k + 1) * (1.0 - al)
            end where
            PB = pt
         end do
      else
         PB = 0.5 * (ple(:, :, km) + ple(:, :, km + 1))
         do k = km, 2, -1
            pt = 0.5 * (ple(:, :, k - 1) + ple(:, :, k))
            if (all(PB < pp)) exit
            where ((pp > pt .and. pp <= PB))
               al = (PB - pp) / (PB - pt)
               v2 = v3(:, :, k - 1) * al + v3(:, :, k) * (1.0 - al)
            end where
            PB = pt
         end do
         pt = 0.5 * (ple(:, :, km) + ple(:, :, km - 1))
         PB = 0.5 * (ple(:, :, km) + ple(:, :, km + 1))
         where ((pp > PB .and. pp <= ple(:, :, km + 1)))
            v2 = v3(:, :, km)
         end where
      end if

      if (present(positive_definite)) then
         if (positive_definite) then
            where (v2 < tiny(0.0))
               v2 = 0.0
            end where
         end if
      end if

      _RETURN(_SUCCESS)
   end subroutine VertInterp

   !BOP
   !IROUTINE: Coldstart

   !DESCRIPTION:
   !   Routine to coldstart from an isothermal state of rest.
   !   The temperature can be specified in the config, otherwise
   !   it is 300K. The surface pressure is assumed to be 1000 hPa.

   !INTERFACE:
   subroutine COLDSTART(gc, import, export, clock, rc)

      use sw, only : sw_phis => surface_geopotential
      use sw, only : sw_hght => height
      use sw, only : sw_uwnd => u_wind
      use sw, only : sw_vwnd => v_wind
      use jw, only : temperature, u_wind, v_wind, surface_geopotential
      use jw, only : tracer_q, tracer_q1_q2, tracer_q3
      use testcases_3_4_5_6, only : advection, Rossby_Haurwitz, mountain_Rossby, gravity_wave

      !ARGUMENTS:
      type(ESMF_GridComp), intent(inout) :: gc
      type(ESMF_State), intent(inout) :: import
      type(ESMF_State), intent(inout) :: export
      type(ESMF_Clock), intent(inout) :: clock
      integer, intent(out), optional :: rc
      !EOP

      integer :: status

      type(ESMF_State) :: internal

      real(kind=REAL8), pointer :: ak1(:), bk1(:), ak(:), bk(:)
      real(kind=REAL8), pointer :: u(:, :, :), v(:, :, :), pt(:, :, :)
      real(kind=REAL8), pointer :: pe1(:, :, :), pkz(:, :, :)
      real(kind=REAL8), allocatable :: pe(:, :, :)
      real(kind=REAL4), pointer :: phis(:, :)
      real(kind=REAL4), allocatable :: lons(:, :), lats(:, :)

      integer :: i, j, k, n, L
      integer :: is, ie, js, je, ks, ke, im, jm, km, ls
      integer :: CASE_ID, case_rotation, case_tracers
      integer :: num_levels

      real :: T0
      real(kind=REAL8) :: dummy_1, dummy_2, dummy_3, dummy_4, dummy_5, dummy_6
      real(kind=REAL8) :: dz, ztop, height, pressure
      real(kind=REAL8) :: LONc, LATc
      real(kind=REAL8) :: eta, eta_top, rot_ang
      real(kind=REAL8) :: ptop, pint
      real(kind=REAL8), allocatable :: ps(:, :)

      type(DynState), pointer :: self
      type(DynGrid), pointer :: grid

      logical :: perturb
      logical :: ak_is_missing = .false.
      logical :: bk_is_missing = .false.
      logical :: FV3_STANDALONE
      ! logical :: isPresent - used in some code that has been commented out

      ! Tracer Stuff
      real(kind=REAL4), pointer :: tracer(:, :, :)
      real(kind=REAL8), allocatable :: Q5(:, :, :)
      real(kind=REAL8), allocatable :: Q6(:, :, :)
      type(ESMF_Geom) :: geom
      type(ESMF_Grid) :: esmfgrid
      type(ESMF_FieldBundle) :: tradv_bundle
      character(len=ESMF_MAXSTR) :: fieldname, string
      real(kind=REAL8), parameter :: r0_6 = 0.6
      real(kind=REAL8), parameter :: r1_0 = 1.0

      class(logger_t), pointer :: logger

      _GET_NAMED_PRIVATE_STATE(gc, DynState, PRIVATE_STATE, self)
      grid => self%grid ! direct handle to grid

      call MAPL_GridCompGet(gc, logger=logger, _RC)
      call MAPL_GridCompGetResource(gc, "T0", T0, default=273., _RC)
      call MAPL_GridCompGetInternalState(gc, internal, _RC)

      call MAPL_GridCompGet(gc, grid=esmfgrid, _RC)
      call MAPL_GridGetCoordinates(esmfgrid, latitudes=lats, longitudes=lons, _RC)
      if (FV_Atm(1)%flagstruct%grid_type == 4) then
         ! Doubly-Period setup based on first LAT/LON coordinate
         lons(:, :) = 0.0
         lats(:, :) = 15.0 * PI / 180.0
      end if

      call MAPL_StateGetPointer(internal, u, "U", _RC) ! A-Grid U Wind
      is = lbound(u, 1)
      ie = ubound(u, 1)
      js = lbound(u, 2)
      je = ubound(u, 2)
      ks = lbound(u, 3)
      ke = ubound(u, 3)
      km = ke - ks + 1
      call MAPL_StateGetPointer(internal, v, "V", _RC) ! A-Grid V Wind
      call MAPL_StateGetPointer(internal, pt, "PT", _RC) ! potential temperature
      ! call MAPL_StateGetPointer(internal, PE1, "PE", _RC) ! edge pressures - 1 based
      ! pchakrab - gfortran has issues with the following rank
      ! remapping (ifort doesn't) , so we allocate a new array
      ! PE(is:ie, js:je, 0:km) => PE1(is:ie, js:je, 1:km+1)
      allocate(pe(is:ie, js:je, 0:km), _STAT)
      call MAPL_StateGetPointer(internal, pkz, "PKZ", _RC) ! presssure ^ kappa at mid-layers
      call MAPL_StateGetPointer(internal, ak1, "AK", _RC) ! AK for vertical coordinate - 1 based
      ak(0:km) => ak1(1:km + 1)
      call MAPL_StateGetPointer(internal, bk1, "BK", _RC) ! BK for vertical coordinate - 1 based
      bk(0:km) => bk1(1:km + 1)
      call MAPL_StateGetPointer(import, phis, "PHIS", _RC) ! surface geopotential

      u = 0.0

      allocate(ps(is:ie, js:je))

      call MAPL_GridCompGetResource(gc, "IM", im, default=0, _RC)
      call MAPL_GridCompGetResource(gc, "JM", jm, default=0, _RC)

      if (km <= 2) then ! Shallow Water

         call MAPL_GridCompGetResource(gc, "CASE_ID", CASE_ID, default=1, _RC)
         DYN_CASE = CASE_ID

         do j = js, je
            do i = is, ie
               LONc = lons(i, j)
               LATc = lats(i, j)
               u(i, j, 1) = sw_uwnd(LONc, LATc, CASE_ID)
               v(i, j, 1) = sw_vwnd(LONc, LATc, CASE_ID)
               pe(i, j, 0) = sw_phis(LONc, LATc, CASE_ID)
               pe(i, j, 1) = sw_hght(LONc, LATc, CASE_ID)
               phis(i, j) = pe(i, j, 0)
            end do
         end do

      else ! 3-D Baroclinic

         u(is:ie, js:je, ke) = .001 * abs(lats(:, :))
         v = 0.0

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

         if (ak_is_missing .or. bk_is_missing) call set_eta(km, ls, ptop, pint, ak, bk)
         _ASSERT(any(ak /= 0.0) .or. any(bk /= 0.0), "needs informative message")

         do L = lbound(pe, 3), ubound(pe, 3)
            pe(:, :, L) = ak(L) + bk(L) * MAPL_P00
         end do
         pkz = 0.5 * (pe(:, :, lbound(pe, 3):ubound(pe, 3) - 1) + pe(:, :, lbound(pe, 3) + 1:ubound(pe, 3)))
         pkz = pkz**MAPL_KAPPA
         pt = T0 / pkz

         ! Check if running standalone model
         call MAPL_GridCompGetResource(gc, "FV3_STANDALONE", FV3_STANDALONE, default=.false., _RC)

         ! 3D Baroclinic Test Cases
         call MAPL_GridCompGetResource(gc, "CASE_ID", CASE_ID, default=0, _RC)
         call MAPL_GridCompGetResource(gc, "CASE_ROTATION", case_rotation, default=0, _RC)
         call MAPL_GridCompGetResource(gc, "CASE_TRACERS", case_tracers, default=1234, _RC)
         DYN_CASE = CASE_ID

         write(string, '(A,I5,A)') "Initializing CASE_ID ", CASE_ID, " in FVcubed:"
         call logger%info(TRIM(string))

         ! Parse case_rotation
         if (case_rotation == -1) rot_ang = 0
         if (case_rotation == 0) rot_ang = 0
         if (case_rotation == 1) rot_ang = 15
         if (case_rotation == 2) rot_ang = 30
         if (case_rotation == 3) rot_ang = 45
         if (case_rotation == 4) rot_ang = 60
         if (case_rotation == 5) rot_ang = 75
         if (case_rotation == 6) rot_ang = 90
         if (case_rotation == -1) then
            grid%f_coriolis_angle = -999
         else
            grid%f_coriolis_angle = rot_ang * PI / 180.0
         end if

         if (CASE_ID == 1) then ! Steady State

            perturb = .false.
            do k = ks, ke
               eta = 0.5 * ((ak(k - 1) + ak(k)) / 1.e5 + bk(k - 1) + bk(k))
               do j = js, je
                  do i = is, ie
                     LONc = lons(i, j)
                     LATc = lats(i, j)
                     u(i, j, k) = u_wind(LONc, LATc, eta, perturb, rot_ang)
                     v(i, j, k) = v_wind(LONc, LATc, eta, perturb, rot_ang)
                     if (k == ks) phis(i, j) = surface_geopotential(LONc, LATc, rot_ang)
                     pt(i, j, k) = temperature(LONc, LATc, eta, rot_ang)
                  end do
               end do
            end do
            pt = pt / pkz

         elseif (CASE_ID == 2) then ! Baroclinic Wave

            perturb = .true.
            do k = ks, ke
               eta = 0.5 * ((ak(k - 1) + ak(k)) / 1.e5 + bk(k - 1) + bk(k))
               do j = js, je
                  do i = is, ie
                     LONc = lons(i, j)
                     LATc = lats(i, j)
                     u(i, j, k) = u_wind(LONc, LATc, eta, perturb, rot_ang)
                     v(i, j, k) = v_wind(LONc, LATc, eta, perturb, rot_ang)
                     if (k == ks) phis(i, j) = surface_geopotential(LONc, LATc, rot_ang)
                     pt(i, j, k) = temperature(LONc, LATc, eta, rot_ang)
                     !if (grid_type==4) then
                     !  if (k==KS) then
                     !     T_PERTURB = (SIN(PI*FLOAT(i-1)/FLOAT(IE-IS))**4.0) * &
                     !                 (SIN(PI*FLOAT(j-1)/FLOAT(JE-JS))**4.0)
                     !     print*, i, j, T_PERTURB
                     !     PT(i,j,k) = PT(i,j,k) + T_PERTURB
                     !  endif
                     !endif
                  end do
               end do
            end do
            pt = pt / pkz

         elseif (CASE_ID == 3) then ! Advection

            !PURE_ADVECTION = .true.

            allocate(Q5(is:ie, js:je, 0:km - 1), _STAT)
            allocate(Q6(is:ie, js:je, 0:km - 1), _STAT)

            ztop = 12000.0
            dz = ztop / km
            do k = ks, ke
               height = (ztop - 0.5 * dz) - (k) * dz ! Layer middle height
               do j = js, je
                  do i = is, ie
                     LONc = lons(i, j)
                     LATc = lats(i, j)
                     call advection('56', LONc, LATc, height, rot_ang, &
                          dummy_1, dummy_2, dummy_3, dummy_4, &
                          ps(i, j), Q5(i, j, k), Q6(i, j, k))
                     u(i, j, k) = dummy_1
                     v(i, j, k) = dummy_2
                     pt(i, j, k) = dummy_3
                     phis(i, j) = dummy_4
                  end do
               end do
            end do
            do L = lbound(pe, 3), ubound(pe, 3)
               pe(:, :, L) = ak(L) + bk(L) * ps(:, :)
            end do

            do k = ks, ke
               do j = js, je
                  do i = is, ie
                     pkz(i, j, k) = ((pe(i, j, k + 1)**KAPPA) - (pe(i, j, k)**KAPPA)) / &
                          (KAPPA * (log(pe(i, j, k + 1)) - log(pe(i, j, k))))
                  end do
               end do
            end do

            pt = pt / pkz

         elseif (CASE_ID == 4) then ! 3D Rossby-Haurwitz

            do j = js, je
               do i = is, ie
                  LONc = lons(i, j)
                  LATc = lats(i, j)
                  pressure = 500.
                  call Rossby_Haurwitz(LONc, LATc, pressure, dummy_1, dummy_2, dummy_3, dummy_4, ps(i, j))
                  u(i, j, 1) = dummy_1
                  v(i, j, 1) = dummy_2
                  pt(i, j, 1) = dummy_3
                  phis(i, j) = dummy_4
               end do
            end do
            do L = lbound(pe, 3), ubound(pe, 3)
               pe(:, :, L) = ak(L) + bk(L) * ps(:, :)
            end do
            do k = ks, ke
               do j = js, je
                  do i = is, ie
                     LONc = lons(i, j)
                     LATc = lats(i, j)
                     pressure = 0.5 * (pe(i, j, k) + pe(i, j, k + 1))
                     call Rossby_Haurwitz(LONc, LATc, pressure, dummy_1, dummy_2, dummy_3, dummy_4, ps(i, j))
                     u(i, j, k) = dummy_1
                     v(i, j, k) = dummy_2
                     pt(i, j, k) = dummy_3
                     phis(i, j) = dummy_4
                  end do
               end do
            end do

            do k = ks, ke
               do j = js, je
                  do i = is, ie
                     pkz(i, j, k) = ((pe(i, j, k + 1)**KAPPA) - (pe(i, j, k)**KAPPA)) / &
                          (KAPPA * (log(pe(i, j, k + 1)) - log(pe(i, j, k))))
                  end do
               end do
            end do
            pt = pt / pkz

         elseif (CASE_ID == 5) then ! Mountain-Induced Rossby Wave

            do k = ks, ke
               do j = js, je
                  do i = is, ie
                     LONc = lons(i, j)
                     LATc = lats(i, j)
                     pressure = 0.5 * (pe(i, j, k) + pe(i, j, k + 1))
                     call mountain_Rossby(case_rotation, LONc, LATc, pressure, dummy_1, dummy_2, dummy_3, dummy_4, ps(i&
                          , j))
                     u(i, j, k) = dummy_1
                     v(i, j, k) = dummy_2
                     pt(i, j, k) = dummy_3
                     phis(i, j) = dummy_4
                  end do
               end do
            end do
            do L = lbound(pe, 3), ubound(pe, 3)
               pe(:, :, L) = ak(L) + bk(L) * ps(:, :)
            end do

            do k = ks, ke
               do j = js, je
                  do i = is, ie
                     pkz(i, j, k) = ((pe(i, j, k + 1)**KAPPA) - (pe(i, j, k)**KAPPA)) / &
                          (KAPPA * (log(pe(i, j, k + 1)) - log(pe(i, j, k))))
                  end do
               end do
            end do

            pt = pt / pkz

         elseif (CASE_ID == 6) then ! Gravity Waves

            ! case_rotation index has different meaning for this test
            if (case_rotation < 3) then
               grid%f_coriolis_angle = -999
            else
               grid%f_coriolis_angle = 0.0
            end if
            ! Get ICs
            ztop = 10000.d0
            dz = ztop / km
            do k = ks, ke
               height = (ztop - 0.5d0 * dz) - (k) * dz ! Layer middle height
               do j = js, je
                  do i = is, ie
                     LONc = lons(i, j)
                     LATc = lats(i, j)
                     call gravity_wave(case_rotation, LONc, LATc, height, dummy_1, dummy_2, dummy_3, dummy_4, ps(i, j))
                     u(i, j, k) = dummy_1
                     v(i, j, k) = dummy_2
                     pt(i, j, k) = dummy_3
                     phis(i, j) = dummy_4
                  end do
               end do
            end do
            ! Reconstruct Edge Pressures and AK BK arrays for rotation=0, otherwise use values from set_eta which are OK
            if (case_rotation == 0) then
               ptop = 27381.905d0
               do k = lbound(pe, 3), ubound(pe, 3)
                  height = ztop - k * dz ! Layer edge height
                  do j = js, je
                     do i = is, ie
                        LONc = lons(i, j)
                        LATc = lats(i, j)
                        call gravity_wave(case_rotation, LONc, LATc, height, dummy_1, dummy_2, dummy_3, dummy_4, &
                             dummy_5, pressure=dummy_6)
                        pe(i, j, k) = dummy_6
                        eta = pe(i, j, k) / ps(i, j)
                        eta_top = ptop / ps(i, j)
                        bk(k) = (eta - eta_top) / (1.d0 - eta_top)
                        ak(k) = 100000.d0 * (eta - bk(k))
                     end do
                  end do
               end do
            end if
            ! Update PE, PKZ and PT
            do L = lbound(pe, 3), ubound(pe, 3)
               pe(:, :, L) = ak(L) + bk(L) * ps(:, :)
            end do

            do k = ks, ke
               do j = js, je
                  do i = is, ie
                     pkz(i, j, k) = ((pe(i, j, k + 1)**KAPPA) - (pe(i, j, k)**KAPPA)) / &
                          (KAPPA * (log(pe(i, j, k + 1)) - log(pe(i, j, k))))
                  end do
               end do
            end do

            pt = pt / pkz

         end if ! case_id

         ! Parse Tracers
         if (FV3_STANDALONE) then
            call ESMF_StateGet(import, "TRADV", tradv_bundle, _RC)
            call MAPL_GridCompGet(gc, geom=geom, num_levels=num_levels, _RC)

            allocate(tracer(is:ie, js:je, 1:km), _STAT)
            tracer(:, :, :) = 0.0
            fieldname = 'Q'
            call addTracer(self, tradv_bundle, tracer, geom, num_levels, fieldname, _RC)

            if (case_tracers /= 1234) then

               do n = 1, case_tracers
                  tracer(:, :, :) = 0.0
                  write(fieldname, "('Q',i3.3)") n
                  call addTracer(self, tradv_bundle, tracer, geom, num_levels, fieldname, _RC)
               end do

            else

               !-------------------
               !     tracer q1
               !-------------------
               tracer(:, :, :) = 0.0
               do k = ks, ke
                  eta = 0.5 * ((ak(k - 1) + ak(k)) / 1.e5 + bk(k - 1) + bk(k))
                  do j = js, je
                     do i = is, ie
                        LONc = lons(i, j)
                        LATc = lats(i, j)
                        dummy_1 = tracer_q1_q2(LONc, LATc, eta, rot_ang, r0_6)
                        tracer(i, j, k) = dummy_1
                     end do
                  end do
               end do
               fieldname = 'Q1'
               call addTracer(self, tradv_bundle, tracer, geom, num_levels, fieldname)

               !-------------------
               !     tracer q2
               !-------------------
               do k = ks, ke
                  eta = 0.5 * ((ak(k - 1) + ak(k)) / 1.e5 + bk(k - 1) + bk(k))
                  do j = js, je
                     do i = is, ie
                        LONc = lons(i, j)
                        LATc = lats(i, j)
                        dummy_1 = tracer_q1_q2(LONc, LATc, eta, rot_ang, r1_0)
                        tracer(i, j, k) = dummy_1
                     end do
                  end do
               end do
               fieldname = 'Q2'
               call addTracer(self, tradv_bundle, tracer, geom, num_levels, fieldname, _RC)

               !-------------------
               !     tracer q3
               !-------------------
               do k = ks, ke
                  eta = 0.5 * ((ak(k - 1) + ak(k)) / 1.e5 + bk(k - 1) + bk(k))
                  do j = js, je
                     do i = is, ie
                        LONc = lons(i, j)
                        LATc = lats(i, j)
                        dummy_1 = tracer_q3(LONc, LATc, eta, rot_ang)
                        tracer(i, j, k) = dummy_1
                     end do
                  end do
               end do
               fieldname = 'Q3'
               call addTracer(self, tradv_bundle, tracer, geom, num_levels, fieldname, _RC)

               !-------------------
               !     tracer q4
               !-------------------
               tracer(:, :, :) = 1.0_r4
               fieldname = 'Q4'
               call addTracer(self, tradv_bundle, tracer, geom, num_levels, fieldname, _RC)

               !-------------------
               !     tracer q5
               !-------------------
               if (allocated(Q5)) then
                  tracer(:, :, :) = Q5(:, :, :)
                  fieldname = 'Q5'
                  call addTracer(self, tradv_bundle, tracer, geom, num_levels, fieldname, _RC)
                  deallocate(Q5, _STAT)
               end if

               !-------------------
               !     tracer q6
               !-------------------
               if (allocated(Q6)) then
                  tracer(:, :, :) = Q6(:, :, :)
                  fieldname = 'Q6'
                  call addTracer(self, tradv_bundle, tracer, geom, num_levels, fieldname, _RC)
                  deallocate(Q6, _STAT)
               end if

            end if

            deallocate(tracer, _STAT)

         end if ! parse tracers

      end if

      call MAPL_StateGetPointer(internal, pe1, "PE", _RC) ! edge pressures - 1 based
      pe1(is:ie, js:je, 1:km + 1) = pe(is:ie, js:je, 0:km)
      deallocate(pe, ps)
      DYN_COLDSTART = .true.

      _RETURN(_SUCCESS)
      _UNUSED_DUMMY(export)
      _UNUSED_DUMMY(clock)
   end subroutine COLDSTART

#ifdef MY_SET_ETA
   subroutine set_eta(km, ptop, ak, bk)

      integer, intent(in) :: km ! vertical dimension
      real(kind=REAL8), intent(out) :: ptop ! model top (Pa)
      real(kind=REAL8), intent(inout) :: ak(km + 1)
      real(kind=REAL8), intent(inout) :: bk(km + 1)

      ! local
      real(kind=REAL8) :: a20_01(21), b20_01(21) ! NCAR Colloquium 20-levels N=0.01
      real(kind=REAL8) :: a20_0178(21), b20_0178(21) ! NCAR Colloquium 20-levels N=0.0178
      real(kind=REAL8) :: a26(27), b26(27) ! NCAR Colloquium 26-levels
      real(kind=REAL8) :: a72(73), b72(73) ! GEOS-5 72-levels
      real(kind=REAL8) :: a137(138), b137(138) ! GEOS-5 137-levels

      real(kind=REAL8) :: p0 = 1000.E2
      real(kind=REAL8) :: pc = 200.E2
      real(kind=REAL8) :: pt, pint, lnpe, dlnp
      real(kind=REAL8) :: press(km + 1)
      integer :: k, ks

      data a20_01 / &
           0.27381905404907E+05, 0.26590539035976E+05, 0.25752394878279E+05, 0.24865429808716E+05, &
           0.23927536347865E+05, 0.22936541085572E+05, 0.21890203071294E+05, 0.20786212168493E+05, &
           0.19622187372385E+05, 0.18395675090318E+05, 0.17104147384052E+05, 0.15745000173179E+05, &
           0.14315551398919E+05, 0.12813039147516E+05, 0.11234619732416E+05, 0.95773657344247E+04, &
           0.78382639990006E+04, 0.60142135898353E+04, 0.41020236978492E+04, 0.20984115047143E+04, &
           0.00000000000000E+00 /

      data b20_01 / &
           0.00000000000000E+00, 0.28901070149364E-01, 0.59510487036309E-01, 0.91902866472543E-01, &
           0.12615517459290E+00, 0.16234678535331E+00, 0.20055953931639E+00, 0.24087780374962E+00, &
           0.28338853406205E+00, 0.32818133660555E+00, 0.37534853286773E+00, 0.42498522508382E+00, &
           0.47718936329560E+00, 0.53206181388604E+00, 0.58970642961892E+00, 0.65023012121324E+00, &
           0.71374293048299E+00, 0.78035810507338E+00, 0.85019217482527E+00, 0.92336502980036E+00, &
           0.10000000000000E+01 /

      data a20_0178 / &
           0.32021324453921E+05, 0.31137565415634E+05, 0.30202026400316E+05, 0.29211673587770E+05, &
           0.28163295404433E+05, 0.27053492108706E+05, 0.25878664766072E+05, 0.24635003578258E+05, &
           0.23318475528610E+05, 0.21924811303582E+05, 0.20449491447964E+05, 0.18887731708932E+05, &
           0.17234467521390E+05, 0.15484337584307E+05, 0.13631666474783E+05, 0.11670446243450E+05, &
           0.95943169315531E+04, 0.73965459465018E+04, 0.50700062290314E+04, 0.26071531411601E+04, &
           0.00000000000000E+00 /

      data b20_0178 / &
           0.00000000000000E+00, 0.27599078219223E-01, 0.56815203138214E-01, 0.87743118501982E-01, &
           0.12048311914891E+00, 0.15514137625266E+00, 0.19183028162025E+00, 0.23066881216269E+00, &
           0.27178291572025E+00, 0.31530591949337E+00, 0.36137896240390E+00, 0.41015145278854E+00, &
           0.46178155290889E+00, 0.51643669184922E+00, 0.57429410846515E+00, 0.63554142614418E+00, &
           0.70037726124166E+00, 0.76901186716541E+00, 0.84166781619770E+00, 0.91858072126555E+00, &
           0.10000000000000E+01 /

      data a26 / &
           219.4067, 489.5209, 988.2418, 1805.2010, 2983.7240, 4462.3340, &
           6160.5870, 7851.2430, 7731.2710, 7590.1310, 7424.0860, &
           7228.7440, 6998.9330, 6728.5740, 6410.5090, 6036.3220, &
           5596.1110, 5078.2250, 4468.9600, 3752.1910, 2908.9490, &
           2084.739, 1334.443, 708.499, 252.1360, 0.0, 0.0 /

      data b26 / &
           0.0, 0.0, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,&
           0.0000000, 0.01505309, 0.03276228, 0.05359622, 0.07810627, &
           0.1069411, 0.1408637, 0.1807720, 0.2277220, 0.2829562, &
           0.3479364, 0.4243822, 0.5143168, 0.6201202, 0.7235355, &
           0.8176768, 0.8962153, 0.9534761, 0.9851122, 1.0000000 /

      data a72 / &
           1.0000000, 2.0000002, 3.2700005, 4.7585009, 6.6000011, &
           8.9345014, 11.970302, 15.949503, 21.134903, 27.852606, &
           36.504108, 47.580610, 61.677911, 79.513413, 101.94402, &
           130.05102, 165.07903, 208.49704, 262.02105, 327.64307, &
           407.65710, 504.68010, 621.68012, 761.98417, 929.29420, &
           1127.6902, 1364.3402, 1645.7103, 1979.1604, 2373.0405, &
           2836.7806, 3381.0007, 4017.5409, 4764.3911, 5638.7912, &
           6660.3412, 7851.2316, 9236.5722, 10866.302, 12783.703, &
           15039.303, 17693.003, 20119.201, 21686.501, 22436.301, &
           22389.800, 21877.598, 21214.998, 20325.898, 19309.696, &
           18161.897, 16960.896, 15625.996, 14290.995, 12869.594, &
           11895.862, 10918.171, 9936.5219, 8909.9925, 7883.4220, &
           7062.1982, 6436.2637, 5805.3211, 5169.6110, 4533.9010, &
           3898.2009, 3257.0809, 2609.2006, 1961.3106, 1313.4804, &
           659.37527, 4.8048257, 0.0000000 /

      data b72 / &
           0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, &
           0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, &
           0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, &
           0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, &
           0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, &
           0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, &
           0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, &
           0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, &
           0.0000000, 8.1754130e-09, 0.0069600246, 0.028010041, 0.063720063, &
           0.11360208, 0.15622409, 0.20035011, 0.24674112, 0.29440312, &
           0.34338113, 0.39289115, 0.44374018, 0.49459020, 0.54630418, &
           0.58104151, 0.61581843, 0.65063492, 0.68589990, 0.72116594, &
           0.74937819, 0.77063753, 0.79194696, 0.81330397, 0.83466097, &
           0.85601798, 0.87742898, 0.89890800, 0.92038701, 0.94186501, &
           0.96340602, 0.98495195, 1.0000000 /

      data a137 / &
           1.000000, 2.000365, 3.102241, 4.666084, 6.827977, 9.746966, 13.605424, 18.608931, 24.985718, 32.985710, &
           42.879242, 54.955463, 69.520576, 86.895882, 107.415741, 131.425507, 159.279404, 191.338562, 227.968948, &
           269.539581, &
           316.420746, 368.982361, 427.592499, 492.616028, 564.413452, 643.339905, 729.744141, 823.967834, 926.344910, &
           1037.20117, &
           1156.853638, 1285.610352, 1423.770142, 1571.622925, 1729.448975, 1897.519287, 2076.095947, 2265.431641, &
           2465.770508, 2677.348145, &
           2900.391357, 3135.119385, 3381.743652, 3640.468262, 3911.490479, 4194.930664, 4490.817383, 4799.149414, &
           5119.895020, 5452.990723, &
           5798.344727, 6156.074219, 6526.946777, 6911.870605, 7311.869141, 7727.412109, 8159.354004, 8608.525391, &
           9076.400391, 9562.682617, &
           10065.978516, 10584.631836, 11116.662109, 11660.067383, 12211.547852, 12766.873047, 13324.668945, &
           13881.331055, 14432.139648, 14975.615234, &
           15508.256836, 16026.115234, 16527.322266, 17008.789062, 17467.613281, 17901.621094, 18308.433594, &
           18685.718750, 19031.289062, 19343.511719, &
           19620.042969, 19859.390625, 20059.931641, 20219.664062, 20337.863281, 20412.308594, 20442.078125, &
           20425.718750, 20361.816406, 20249.511719, &
           20087.085938, 19874.025391, 19608.572266, 19290.226562, 18917.460938, 18489.707031, 18006.925781, &
           17471.839844, 16888.687500, 16262.046875, &
           15596.695312, 14898.453125, 14173.324219, 13427.769531, 12668.257812, 11901.339844, 11133.304688, &
           10370.175781, 9617.515625, 8880.453125, &
           8163.375000, 7470.343750, 6804.421875, 6168.531250, 5564.382812, 4993.796875, 4457.375000, 3955.960938, &
           3489.234375, 3057.265625, &
           2659.140625, 2294.242188, 1961.500000, 1659.476562, 1387.546875, 1143.250000, 926.507812, 734.992188, &
           568.062500, 424.414062, &
           302.476562, 202.484375, 122.101562, 62.781250, 22.835938, 3.757813, 0.000000, 0.000000/

      data b137 / &
           0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, &
           0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, &
           0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, &
           0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, &
           0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, &
           0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000007, 0.000024, 0.000059, 0.000112, 0.000199, &
           0.000340, 0.000562, 0.000890, 0.001353, 0.001992, 0.002857, 0.003971, 0.005378, 0.007133, 0.009261, &
           0.011806, 0.014816, 0.018318, 0.022355, 0.026964, 0.032176, 0.038026, 0.044548, 0.051773, 0.059728, &
           0.068448, 0.077958, 0.088286, 0.099462, 0.111505, 0.124448, 0.138313, 0.153125, 0.168910, 0.185689, &
           0.203491, 0.222333, 0.242244, 0.263242, 0.285354, 0.308598, 0.332939, 0.358254, 0.384363, 0.411125, &
           0.438391, 0.466003, 0.493800, 0.521619, 0.549301, 0.576692, 0.603648, 0.630036, 0.655736, 0.680643, &
           0.704669, 0.727739, 0.749797, 0.770798, 0.790717, 0.809536, 0.827256, 0.843881, 0.859432, 0.873929, &
           0.887408, 0.899900, 0.911448, 0.922096, 0.931881, 0.940860, 0.949064, 0.956550, 0.963352, 0.969513, &
           0.975078, 0.980072, 0.984542, 0.988500, 0.991984, 0.995003, 0.997630, 1.000000/

      select case (km)

      case (20)

         do k = 1, km + 1
            ak(k) = a20_0178(k)
            bk(k) = b20_0178(k)
         end do
         ! Search KS
         ks = 0
         do k = 1, km
            if (bk(k) > 0) then
               ks = k - 1
               goto 120
            end if
         end do
         120 continue

      case (26)

         do k = 1, km + 1
            ak(k) = a26(k)
            bk(k) = b26(k)
         end do
         ! Search KS
         ks = 0
         do k = 1, km
            if (bk(k) > 0) then
               ks = k - 1
               goto 126
            end if
         end do
         126 continue

      case (40)
         !--------------------------------------------------
         ! Pure sigma-coordinate with uniform spacing in "z"
         !--------------------------------------------------
         ptop = 27381.905404907 ! model top pressure (pascal)
         press(1) = ptop
         press(km + 1) = p0
         dlnp = (log(p0) - log(ptop)) / real(km)

         lnpe = log(press(km + 1))
         do k = km, 2, -1
            lnpe = lnpe - dlnp
            press(k) = exp(lnpe)
         end do

         ! Search KS
         ks = 0
         do k = 1, km
            if (press(k) >= pc) then
               ks = k - 1
               goto 140
            end if
         end do
         140 continue

         if (ks /= 0) then
            do k = 1, ks
               ak(k) = press(k)
               bk(k) = 0.
            end do
         end if

         pint = press(ks + 1)
         do k = ks + 1, km
            ak(k) = pint * (press(km) - press(k)) / (press(km) - pint)
            bk(k) = (press(k) - ak(k)) / press(km + 1)
         end do
         ak(km + 1) = 0.
         bk(km + 1) = 1.

      case (60)
         !--------------------------------------------------
         ! Pure sigma-coordinate with uniform spacing in "z"
         !--------------------------------------------------
         ptop = 25499.234876157 ! model top pressure (pascal)
         press(1) = ptop
         press(km + 1) = p0
         dlnp = (log(p0) - log(ptop)) / real(km)

         lnpe = log(press(km + 1))
         do k = km, 2, -1
            lnpe = lnpe - dlnp
            press(k) = exp(lnpe)
         end do

         ! Search KS
         ks = 0
         do k = 1, km
            if (press(k) >= pc) then
               ks = k - 1
               goto 160
            end if
         end do
         160 continue

         if (ks /= 0) then
            do k = 1, ks
               ak(k) = press(k)
               bk(k) = 0.
            end do
         end if

         pint = press(ks + 1)
         do k = ks + 1, km
            ak(k) = pint * (press(km) - press(k)) / (press(km) - pint)
            bk(k) = (press(k) - ak(k)) / press(km + 1)
         end do
         ak(km + 1) = 0.
         bk(km + 1) = 1.

      case (72)

         do k = 1, km + 1
            ak(k) = a72(k)
            bk(k) = b72(k)
         end do
         ! Search KS
         ks = 0
         do k = 1, km
            if (bk(k) > 0) then
               ks = k - 1
               goto 172
            end if
         end do
         172 continue

      case (137)

         do k = 1, km + 1
            ak(k) = a137(k)
            bk(k) = b137(k)
         end do
         ! Search KS
         ks = 0
         do k = 1, km
            if (bk(k) > 0) then
               ks = k - 1
               goto 137
            end if
         end do
         137 continue

      case default

         print*, 'Bad KM in FVdycoreCubed_GridComp:set_eta', km

      end select

   end subroutine set_eta
#endif

   subroutine addTracer_r8(self, bundle, var, geom, num_levels, fieldname, rc)
      type(DynState), pointer :: self
      type(ESMF_FieldBundle) :: bundle
      real(kind=r8), pointer :: var(:, :, :)
      type(ESMF_Geom), intent(in) :: geom
      integer, intent(in) :: num_levels
      character(len=*), intent(in) :: fieldname
      integer, optional, intent(out) :: rc

      integer :: nq, status
      type(DynTracers), pointer :: t(:)
      type(ESMF_Field) :: field
      real(kind=r8), pointer :: ptr(:, :, :)

      call ESMF_FieldBundleGet(bundle, fieldCount=nq, _RC)
      nq = nq + 1

      field = MAPL_FieldCreate( &
           geom, typekind=ESMF_TYPEKIND_R8, &
           name=fieldname, &
           num_levels=num_levels, &
           vert_staggerloc=VERTICAL_STAGGER_CENTER, &
           units="unknown_units", &
           standard_name="unknown_standard_name", &
           long_name="unknown_long_name", _RC)
      call ESMF_FieldGet(field, localDE=0, farrayptr=ptr, _RC)
      ptr = var
      call MAPL_FieldBundleAdd(bundle, [field], _RC)

      if (nq == 1) then
         allocate(self%vars%tracer(nq), _STAT)
         self%vars%tracer(nq)%content => ptr
         self%vars%tracer(nq)%is_r4 = .false.
      else
         allocate(t(nq))
         t(1:nq - 1) = self%vars%tracer
         deallocate(self%vars%tracer)
         self%vars%tracer => t
         self%vars%tracer(nq)%content => ptr
         self%vars%tracer(nq)%is_r4 = .false.
      end if

      self%grid%nq = nq

      _RETURN(_SUCCESS)
   end subroutine addTracer_r8

   subroutine addTracer_r4(self, bundle, var, geom, num_levels, fieldname, rc)
      type(DynState), pointer :: self
      type(ESMF_FieldBundle) :: bundle
      real(kind=r4), pointer :: var(:, :, :)
      type(ESMF_Geom), intent(in) :: geom
      integer, intent(in) :: num_levels
      character(len=*), intent(in) :: fieldname
      integer, optional, intent(out) :: rc

      integer :: nq, status
      type(DynTracers), pointer :: t(:)
      type(ESMF_Field) :: field
      real(kind=r4), pointer :: ptr(:, :, :)

      call ESMF_FieldBundleGet(bundle, fieldCount=nq, _RC)
      nq = nq + 1

      field = MAPL_FieldCreate( &
           geom, typekind=ESMF_TYPEKIND_R4, &
           name=fieldname, &
           num_levels=num_levels, &
           vert_staggerloc=VERTICAL_STAGGER_CENTER, &
           units="unknown_units", &
           standard_name="unknown_standard_name", &
           long_name="unknown_long_name", _RC)
      call ESMF_FieldGet(field, localDE=0, farrayptr=ptr, _RC)
      ptr = var
      call MAPL_FieldBundleAdd(bundle, [field], _RC)

      if (nq == 1) then
         allocate(self%vars%tracer(nq), _STAT)
         self%vars%tracer(nq)%content_r4 => ptr
         self%vars%tracer(nq)%is_r4 = .true.
      else
         allocate(t(nq))
         t(1:nq - 1) = self%vars%tracer
         deallocate(self%vars%tracer)
         self%vars%tracer => t
         self%vars%tracer(nq)%content_r4 => ptr
         self%vars%tracer(nq)%is_r4 = .true.
      end if

      self%grid%nq = nq

      _RETURN(_SUCCESS)
   end subroutine addTracer_r4

   subroutine freeTracers(self)
      type(DynState) :: self

      if (associated(self%vars%tracer)) then
         deallocate(self%vars%tracer) ! Comment out to output tracer to checkpoint file
         nullify(self%vars%tracer)
      end if

      return
   end subroutine freeTracers

   subroutine Write_Profile_2d_R8(grid, arr, name)
      type(DynGrid), intent(in) :: grid
      real(kind=r8), intent(in) :: arr(grid%is:grid%ie, grid%js:grid%je)
      character(len=*), intent(in) :: name

      integer :: istrt, iend, jstrt, jend
      integer :: im, jm
      real(kind=r8) :: arr_global(grid%npx, grid%ntiles * grid%npy)
      real(kind=r8) :: rng(3)
      real(kind=r8) :: GSUM

      real(kind=ESMF_KIND_R8) :: locArr(grid%is:grid%ie, grid%js:grid%je)
      real(kind=ESMF_KIND_R8) :: glbArr(grid%npx, grid%ntiles * grid%npy)

      istrt = grid%is
      iend = grid%ie
      jstrt = grid%js
      jend = grid%je
      im = grid%npx
      jm = grid%npy * grid%ntiles

      !call write_parallel('GlobalSUm')
      locArr(:, :) = arr(:, :)
      call MAPL_ArrayGather(locArr, glbArr, grid%grid)
      arr_global(:, :) = glbArr

      if (MAPL_Am_I_Root()) then
         rng(1) = MINVAL(MINVAL(arr_global, DIM=1), DIM=1)
         rng(2) = MAXVAL(MAXVAL(arr_global, DIM=1), DIM=1)
         rng(3) = SUM(SUM(arr_global, DIM=1), DIM=1) / (im * jm)
         GSUM = SUM(SUM(arr_global, DIM=1), DIM=1)

         print*, '***********'
         print*, 'stats for ', TRIM(name)

         write(*, '(3(f21.9,1x))')rng(:)
         !   Write(*,"('GlobalSum: ',f21.9)") GSUM
         print*, '***********'
         print*, ' '
      end if

   end subroutine Write_Profile_2d_R8

   subroutine Write_Profile_2d_R4(grid, arr, name)
      type(DynGrid), intent(in) :: grid
      real(kind=r4), intent(in) :: arr(grid%is:grid%ie, grid%js:grid%je)
      character(len=*), intent(in) :: name

      integer :: istrt, iend, jstrt, jend
      integer :: im, jm
      real(kind=r4) :: arr_global(grid%npx, grid%ntiles * grid%npy)
      real(kind=r4) :: rng(3)
      real(kind=r4) :: GSUM

      real(kind=ESMF_KIND_R4) :: locArr(grid%is:grid%ie, grid%js:grid%je)
      real(kind=ESMF_KIND_R4) :: glbArr(grid%npx, grid%ntiles * grid%npy)

      istrt = grid%is
      iend = grid%ie
      jstrt = grid%js
      jend = grid%je
      im = grid%npx
      jm = grid%npy * grid%ntiles

      ! call write_parallel('GlobalSUm')
      locArr(:, :) = arr(:, :)
      call MAPL_ArrayGather(locArr, glbArr, grid%grid)
      arr_global(:, :) = glbArr

      if (MAPL_Am_I_Root()) then
         rng(1) = MINVAL(MINVAL(arr_global, DIM=1), DIM=1)
         rng(2) = MAXVAL(MAXVAL(arr_global, DIM=1), DIM=1)
         rng(3) = SUM(SUM(arr_global, DIM=1), DIM=1) / (im * jm)
         GSUM = SUM(SUM(arr_global, DIM=1), DIM=1)

         print*, '***********'
         print*, 'stats for ', TRIM(name)

         write(*, '(3(f21.9,1x))')rng(:)
         !    Write(*,"('GlobalSum: ',f21.9)") GSUM
         print*, '***********'
         print*, ' '
      end if

   end subroutine Write_Profile_2d_R4

   subroutine Write_Profile_R8(grid, arr, name)
      type(DynGrid), intent(in) :: grid
      real(kind=r8), intent(in) :: arr(grid%is:grid%ie, grid%js:grid%je, 1:grid%npz)
      character(len=*), intent(in) :: name

      integer :: istrt, iend, jstrt, jend, kstrt, kend
      integer :: im, jm, km, k
      real(kind=r8), allocatable :: arr_global(:, :, :)
      real(kind=r8) :: rng(3, grid%npz)
      real(kind=r8) :: GSUM
      logical :: amIRoot

      real(kind=ESMF_KIND_R8) :: locArr(grid%is:grid%ie, grid%js:grid%je)
      istrt = grid%is
      iend = grid%ie
      jstrt = grid%js
      jend = grid%je
      kstrt = 1
      kend = grid%npz
      im = grid%npx
      jm = grid%npy * grid%ntiles
      km = grid%npz

      amIRoot = MAPL_Am_I_Root()
      if (amIRoot) then
         allocate(arr_global(grid%npx, grid%ntiles * grid%npy, km))
      else
         allocate(arr_global(1, 1, km))
      end if

      ! call write_parallel('GlobalSUm')
      do k = kstrt, kend
         locArr(:, :) = arr(:, :, k)
         call MAPL_ArrayGather(locArr, arr_global(:, :, k), grid%grid)
      end do

      if (amIRoot) then
         rng(1, :) = MINVAL(MINVAL(arr_global, DIM=1), DIM=1)
         rng(2, :) = MAXVAL(MAXVAL(arr_global, DIM=1), DIM=1)
         rng(3, :) = SUM(SUM(arr_global, DIM=1), DIM=1) / (im * jm)
         GSUM = SUM(SUM(SUM(arr_global, DIM=1), DIM=1), DIM=1)

         print*, '***********'
         print*, 'stats for ', TRIM(name)

         do k = 1, km
            write(*, '(a,i4.0,3(f21.9,1x))')'k:', k, rng(:, k)
         end do
         !    Write(*,"('GlobalSum: ',f21.9)") GSUM
         print*, '***********'
         print*, ' '
      end if

      deallocate(arr_global)

   end subroutine Write_Profile_R8

   subroutine Write_Profile_R4(grid, arr, name, delp)
      type(DynGrid), intent(in) :: grid
      real(kind=r4), intent(in) :: arr(grid%is:grid%ie, grid%js:grid%je, 1:grid%npz)
      character(len=*), intent(in) :: name
      real(kind=r8), optional, intent(in) :: delp(grid%is:grid%ie, grid%js:grid%je, 1:grid%npz)

      integer :: istrt, iend, jstrt, jend, kstrt, kend
      integer :: im, jm, km, k
      real(kind=r4), allocatable :: arr_global(:, :, :)
      real(kind=r4) :: rng(3, grid%npz)
      real(kind=r8) :: gsum_p
      real(kind=r4) :: GSUM
      logical :: amIRoot

      real(kind=ESMF_KIND_R8) :: locArr(grid%is:grid%ie, grid%js:grid%je)
      real(kind=ESMF_KIND_R8), allocatable :: glbArr(:, :)

      istrt = grid%is
      iend = grid%ie
      jstrt = grid%js
      jend = grid%je
      kstrt = 1
      kend = grid%npz
      im = grid%npx
      jm = grid%npy * grid%ntiles
      km = grid%npz

      amIRoot = MAPL_Am_I_Root()
      if (amIRoot) then
         allocate(arr_global(grid%npx, grid%ntiles * grid%npy, km))
         allocate(glbArr(grid%npx, grid%ntiles * grid%npy))
      else
         allocate(arr_global(1, 1, km))
         allocate(glbArr(1, 1))
      end if

      do k = kstrt, kend
         locArr(:, :) = arr(:, :, k)
         call MAPL_ArrayGather(locArr, glbArr, grid%grid)
         if (amIRoot) then
            arr_global(:, :, k) = glbArr
         end if
      end do
      if (amIRoot) then
         rng(1, :) = MINVAL(MINVAL(arr_global, DIM=1), DIM=1)
         rng(2, :) = MAXVAL(MAXVAL(arr_global, DIM=1), DIM=1)
         rng(3, :) = SUM(SUM(arr_global, DIM=1), DIM=1) / (im * jm)
         print*, '***********'
         print*, 'stats for ', TRIM(name)
         do k = 1, km
            write(*, '(a,i4.0,3(f21.9,1x))')'k:', k, rng(:, k)
         end do
         print*, '***********'
         print*, ' '
      end if

      if (present(delp)) then
         gsum_p = 0
         do k = kstrt, kend
            locArr(:, :) = arr(:, :, k) * grid%area(:, :) * delp(:, :, k)
            call MAPL_ArrayGather(locArr, glbArr, grid%grid)
            if (amIRoot) then
               arr_global(:, :, k) = glbArr
            end if
            locArr(:, :) = delp(:, :, k)
            call MAPL_ArrayGather(locArr, glbArr, grid%grid)
            if (amIRoot) then
               gsum_p = gsum_p + SUM(SUM(glbArr, DIM=1), DIM=1)
            end if
         end do
         if (amIRoot) then
            GSUM = SUM(SUM(SUM(arr_global, DIM=1), DIM=1), DIM=1)
            print*, '***********'
            write(*, "('GlobalSum: ',e21.9)") GSUM / (grid%globalarea * gsum_p)
            print*, '***********'
            print*, ' '
         end if
      end if

      deallocate(arr_global, glbArr)

   end subroutine Write_Profile_R4

   function R8_TO_R4(dbl_var)
      real(kind=REAL8), intent(in) :: dbl_var(:, :)
      real(kind=REAL4) :: R8_TO_R4(lbound(dbl_var, 1):ubound(dbl_var, 1), lbound(dbl_var, 2):ubound(dbl_var, 2))
      integer :: i, j

      real(kind=REAL8), parameter :: EPS = 1.e-15_REAL8
      real(kind=REAL8), parameter :: big = 1.e15_REAL8

      do j = lbound(dbl_var, 2), ubound(dbl_var, 2)
         do i = lbound(dbl_var, 1), ubound(dbl_var, 1)
            R8_TO_R4(i, j) = SIGN(MIN(big, MAX(EPS, abs(dbl_var(i, j)))), dbl_var(i, j))
         end do
      end do
   end function R8_TO_R4

   function R4_TO_R8(sngl_var)
      real(kind=REAL4), intent(in) :: sngl_var(:, :)
      real(kind=REAL8) :: R4_TO_R8(lbound(sngl_var, 1):ubound(sngl_var, 1), lbound(sngl_var, 2):ubound(sngl_var, 2))
      integer :: i, j

      real(kind=REAL4), parameter :: EPS = 1.e-15_REAL4
      real(kind=REAL4), parameter :: big = 1.e15_REAL4

      do j = lbound(sngl_var, 2), ubound(sngl_var, 2)
         do i = lbound(sngl_var, 1), ubound(sngl_var, 1)
            R4_TO_R8(i, j) = SIGN(MIN(big, MAX(EPS, abs(sngl_var(i, j)))), sngl_var(i, j))
         end do
      end do
   end function R4_TO_R8

end module FVdycoreCubed_GridComp

subroutine SetServices(gc, rc)
   use ESMF
   use FVdycoreCubed_GridComp, only : mySetservices => SetServices
   type(ESMF_GridComp) :: gc
   integer, intent(out) :: rc
   call mySetservices(gc, rc=rc)
end subroutine SetServices
