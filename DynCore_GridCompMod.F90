!  $id: DynCore_GridCompMod.F90,v 1.1.1.1 2007/05/29 12:26:20 atrayanov Exp $

#include "MAPL_Generic.h"

!#define SCALAR_WINDS
!#define INC_WINDS

!-----------------------------------------------------------------------
!              ESMA - Earth System Modeling Applications
!-----------------------------------------------------------------------
   Module CoarseFVdycoreCubed_GridComp

!BOP
!
! !MODULE: FVdycoreCubed_GridComp --- Dynamical Core Grid Component
!

! !USES:

   use ESMF                ! ESMF base class
   use MAPL                ! GEOS base class
   use m_set_eta,       only: set_eta

   use SSI_FineToCoarse, only: SSI_CopyFineToCoarse, SSI_BundleCopyFineToCoarse
   use SSI_CoarseToFine, only: SSI_CopyCoarseToFine, SSI_BundleCopyCoarseToFine
! FV Specific Module
   use fv_arrays_mod,  only: REAL4, REAL8, FVPRC
   !use fv_grid_tools_mod, only: grid_type
   use FV_StateMod, only : FV_Atm,                                   &
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
                           fv_getUpdraftHelicity,                    &
                           ADIABATIC, SW_DYNAMICS, AdvCore_Advection, &
                           INTERNAL_FineToCoarse, INTERNAL_CoarseToFine
   use m_topo_remap, only: dyn_topo_remap
   use CubeGridPrototype, only: register_grid_and_regridders

! !PUBLIC MEMBER FUNCTIONS:

  implicit none
  private

  ! Include the MPI library definitons:
  include 'mpif.h'

  type(ESMF_FieldBundle), save :: bundleAdv
  integer :: NXQ = 0
  logical :: overwrite_Q = .true.

  public  SetServices      ! Register component methods
  public coarse_setvm
  public DYN_wrap

! !DESCRIPTION: This module implements the Dynamical Core as
!               an ESMF gridded component.
!
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
!
! !REVISION HISTORY:
!
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
!
!----------------------------------------------------------------------

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
! -------------------------------------
    integer, parameter         :: nlevs=5
    integer, parameter         :: ntracers=11
    integer                    :: nlev, ntracer
    integer                    :: plevs(nlevs)
    character(len=ESMF_MAXSTR) :: myTracer
    data plevs /850,700,600,500,300/

! Wrapper for extracting internal state
! -------------------------------------

  type DYN_wrap
     type (DynState), pointer :: DYN_STATE
  end type DYN_wrap

  interface addTracer
     module procedure addTracer_r4
     module procedure addTracer_r8
  end interface

  interface Write_Profile
     module procedure Write_Profile_R4
     module procedure Write_Profile_R8
     module procedure Write_Profile_2d_R4
     module procedure Write_Profile_2d_R8
  end interface

  real(kind=8) :: t1, t2
  real(kind=8) :: dyn_run_timer

  logical :: DO_ADD_INCS = .true.

  type(ESMF_State) :: internal

contains

   subroutine coarse_setvm(gc, rc)
      type(ESMF_GridComp) :: gc
      integer, intent(out) :: rc

      type(ESMF_VM) :: vm
      logical       :: pthreadsEnabled
      logical       :: ssiSharedMemoryEnabled
      !integer       :: ssiMaxPetCount
      character(len=160)    :: msg
      integer       :: nthreads
      integer       :: nth_x, nth_y
      integer :: status
      character(len=ESMF_MAXSTR) :: Iam = "coarse_setvm"
      type (MAPL_MetaComp), pointer :: MAPL  => NULL()
      integer, allocatable :: gcImg(:)
      integer :: itemCount, esmf_stacksize
      type(ESMF_GridComp) :: fineGC
      !type(ESMF_Config)  :: cf

      ! Initialize return code
      rc = ESMF_SUCCESS

! Retrieve fine GC
! ---------------------------------
      call ESMF_AttributeGet(GC, name='GC_IMAGE', itemCount=itemCount, rc=status)
      VERIFY_(STATUS)
      allocate(gcImg(itemCount), stat=status)
      VERIFY_(STATUS)
      call ESMF_AttributeGet(GC, name='GC_IMAGE', valueList=gcImg, rc=status)
      VERIFY_(STATUS)
      fineGC = transfer(gcImg, fineGC)
      deallocate(gcImg, stat=status)
      VERIFY_(STATUS)

! Retrieve the pointer to the state
! ---------------------------------
      call MAPL_GetObjectFromGC (fineGC, MAPL,  RC=STATUS )
      VERIFY_(STATUS)

      ! The following call will give each PET as many PEs as nthreads.
      ! This will reduce the number of PETs that are
      ! executing the component, but each PET will have multipe PEs available,
      ! e.g. to do user-level OpenMP threading.

      ! First test whether ESMF-threading is supported on this machine
      call ESMF_VMGetCurrent(vm, rc=status)
      VERIFY_(STATUS)
      call ESMF_VMGet(vm, pthreadsEnabledFlag=pthreadsEnabled, &
        ssiSharedMemoryEnabledFlag=ssiSharedMemoryEnabled, rc=status)
      VERIFY_(STATUS)
      _ASSERT(ssiSharedMemoryEnabled, 'ESMF built with Shared Memory Required')
      _ASSERT(pthreadsEnabled, 'ESMF built with Pthreads Enabled Required')
! nth_x = coarsening factor in X-direction
        call MAPL_GetResource( MAPL, nth_x, 'NTH_X:', default=1, RC=STATUS )
        VERIFY_(STATUS)
! nth_y = coarsening factor in Y-direction
        call MAPL_GetResource( MAPL, nth_y, 'NTH_Y:', default=1, RC=STATUS )
        VERIFY_(STATUS)
! esmf_stacksize for the main Pthread
        call MAPL_GetResource( MAPL, esmf_stacksize, 'ESMF_STACKSIZE:', &
                 default=20971520, RC=STATUS )
        VERIFY_(STATUS)
! nthreads = num threads to use in dyncore
        nthreads = nth_x*nth_y
        call ESMF_GridCompSetVMMaxPEs(gc, maxPeCountPerPet=nthreads, &
                 pthreadMinStackSize=esmf_stacksize, rc=status)
        VERIFY_(STATUS)

      RETURN_(ESMF_SUCCESS)

   end subroutine

!----------------------------------------------------------------------
!BOP
!
! !IROUTINE: SetServices

! !DESCRIPTION:  SetServices registers Initialize, Run, and Finalize
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
!
!
!
! !INTERFACE:

   Subroutine SetServices ( gc, rc )

! !ARGUMENTS:

   type(ESMF_GridComp) :: gc     ! gridded component
   integer, intent(out)     :: rc     ! return code


! !DESCRIPTION: Set services (register) for the FVCAM Dynamical Core
!               Grid Component.
!
!EOP
!----------------------------------------------------------------------

    type (DynState), pointer :: state
    type (DYN_wrap)                  :: wrap

    integer                          :: status
   character(len=ESMF_MAXSTR) :: Iam = "CoarseSetServices"
   integer, allocatable :: gcImg(:)
   integer :: itemCount
   type(ESMF_GridComp) :: fineGC
   type(MAPL_MetaComp), pointer :: MAPL
   integer :: nx, ny, nnx, nny, nth_x, nth_y, agcm_im

! Retrieve fine GC
! ---------------------------------
    call ESMF_AttributeGet(GC, name='GC_IMAGE', itemCount=itemCount, rc=status)
    VERIFY_(STATUS)
    allocate(gcImg(itemCount), stat=status)
    VERIFY_(STATUS)
    call ESMF_AttributeGet(GC, name='GC_IMAGE', valueList=gcImg, rc=status)
    VERIFY_(STATUS)
    fineGC = transfer(gcImg, fineGC)
    deallocate(gcImg, stat=status)
    VERIFY_(STATUS)

! Retrieve the pointer to the state
! ---------------------------------

    call MAPL_GetObjectFromGC (fineGC, MAPL,  RC=STATUS )
    VERIFY_(STATUS)

! Checks aligments of number of threads, procs/node, etc.

    call MAPL_GetResource( MAPL, nx, 'NX:', default=0, RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, nth_x, 'NTH_X:', default=1, RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, ny, 'NY:', default=0, RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, nth_y, 'NTH_Y:', default=1, RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, nnx, 'NNX:', default=1, RC=STATUS          )
    VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, nny, 'NNY:', default=1, RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, agcm_im, 'AGCM_IM:', RC=STATUS )
    VERIFY_(STATUS)

    _ASSERT(mod(nx*ny, nnx*nny) == 0, 'num_procs/node must evenly divide total num_procs')
    _ASSERT(mod(nnx, nth_x) == 0, 'coarsening factor in X-direction must evenly divide num_procs/node in X-direction')
    _ASSERT(mod(nx, nth_x) == 0, 'coarsening factor in X-direction must evenly divide num_procs in X-direction')
    _ASSERT(mod(nny, nth_y) == 0, 'coarsening factor in Y-direction must evenly divide num_procs/node in Y-direction')
    _ASSERT(mod(ny/6, nth_y) == 0, 'coarsening factor in Y-direction must evenly divide num_procs in Y-direction')
    _ASSERT(mod(agcm_im, nx) == 0, 'subdomain size in X-direction must be equal')


! Allocate this instance of the internal state and put it in wrapper.
! -------------------------------------------------------------------

    allocate( state, stat=status )
    VERIFY_(STATUS)
    wrap%dyn_state => state

! Save pointer to the wrapped internal state in the GC
! ----------------------------------------------------

    call ESMF_UserCompSetInternalState ( GC,'DYNstate',wrap,status )
    VERIFY_(STATUS)

! Register services for this component
! ------------------------------------
    call ESMF_GridCompSetEntryPoint (gc, ESMF_METHOD_INITIALIZE, &
       userRoutine=set_esmf_internal_state, PHASE=1, rc=status)
    VERIFY_(STATUS)

    call ESMF_GridCompSetEntryPoint (gc, ESMF_METHOD_INITIALIZE, &
       userRoutine=Initialize, PHASE=2, rc=status)
    VERIFY_(STATUS)

    call ESMF_GridCompSetEntryPoint (gc, ESMF_METHOD_INITIALIZE, &
       userRoutine=get_esmf_internal_state, PHASE=3, rc=status)
    VERIFY_(STATUS)

    call ESMF_GridCompSetEntryPoint (gc, ESMF_METHOD_RUN, &
       userRoutine=Run, PHASE=1, rc=status)
    VERIFY_(STATUS)

    call ESMF_GridCompSetEntryPoint (gc, ESMF_METHOD_RUN, &
       userRoutine=RunAddIncs, PHASE=2, rc=status)
    VERIFY_(STATUS)

    call ESMF_GridCompSetEntryPoint (gc, ESMF_METHOD_FINALIZE, &
       userRoutine=Finalize, rc=status)
    VERIFY_(STATUS)

    call ESMF_UserCompGetInternalState(gc, 'DYNstate', wrap, status)
    VERIFY_(STATUS)
    state => wrap%dyn_state

    call DynSetup(GC, state, rc=status)
    VERIFY_(STATUS)
    RETURN_(ESMF_SUCCESS)

  end subroutine SetServices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine set_esmf_internal_state(gc, import, export, clock, rc)
     type(ESMF_GridComp) :: gc       !  gridded component
     type(ESMF_State)    :: import   ! import state
     type(ESMF_State)    :: export   ! export state
     type(ESMF_Clock)    :: clock    ! the clock

     integer, intent(out)     :: rc       ! Error code:
     character(len=ESMF_MAXSTR)         :: IAm

     IAm = 'set_esmf_internal_state'

     internal = import

     RETURN_(ESMF_SUCCESS)

  end subroutine set_esmf_internal_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine get_esmf_internal_state(gc, import, export, clock, rc)
     type(ESMF_GridComp) :: gc       !  gridded component
     type(ESMF_State)    :: import   ! import state
     type(ESMF_State)    :: export   ! export state
     type(ESMF_Clock)    :: clock    ! the clock

     integer, intent(out)     :: rc       ! Error code:
     character(len=ESMF_MAXSTR)         :: IAm

     IAm = 'get_esmf_internal_state'

     export = internal

     RETURN_(ESMF_SUCCESS)

  end subroutine get_esmf_internal_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine Initialize ( gc, import, export, clock, rc )

! !ARGUMENTS:

  type(ESMF_GridComp) :: gc       ! composite gridded component
  type(ESMF_State) :: import   ! import state
  type(ESMF_State) :: export   ! export state
  type(ESMF_Clock) :: clock    ! the clock

  integer, intent(out)  :: rc       ! Error code:
                                                 ! = 0 all is well
                                                 ! otherwise, error
  type (ESMF_Config)                 :: cf

  type (DYN_wrap)                    :: wrap
  type (DynState),  pointer  :: STATE

  type (MAPL_MetaComp),      pointer :: mapl

  character (len=ESMF_MAXSTR)        :: layout_file

  type (ESMF_Field)                  :: field

  real(r8), pointer                  ::  ak(:)
  real(r8), pointer                  ::  bk(:)
  real(r8), pointer                  ::  ud(:,:,:)
  real(r8), pointer                  ::  vd(:,:,:)
  real(r8), pointer                  ::  pe(:,:,:)
  real(r8), pointer                  ::  pt(:,:,:)
  real(r8), pointer                  ::  pk(:,:,:)

  real(r8), allocatable              ::  ur(:,:,:)
  real(r8), allocatable              ::  vr(:,:,:)

  real(r4), pointer                  :: ple(:,:,:)
  real(r4), pointer                  ::   u(:,:,:)
  real(r4), pointer                  ::   v(:,:,:)
  real(r4), pointer                  ::   t(:,:,:)

  character(len=ESMF_MAXSTR)         :: ReplayMode
  real                               :: DNS_INTERVAL
  type (ESMF_TimeInterval)           :: Intv
  type (ESMF_Alarm)                  :: Alarm
  integer                            :: ColdRestart=0

  integer                            :: status
  character(len=ESMF_MAXSTR)         :: IAm
  character(len=ESMF_MAXSTR)         :: COMP_NAME

  type (DynGrid),  pointer           :: DycoreGrid

  real(r4), pointer                      :: temp2d(:,:)
  real(r4), pointer                      :: temp3d(:,:,:) => Null()

  integer                            :: ifirst
  integer                            :: ilast
  integer                            :: jfirst
  integer                            :: jlast
  integer                            :: km
  type(ESMF_FieldBundle)             :: tradv, tradvex
  integer                            :: i,numTracers,fv3_standalone

! Begin coarse GC
  integer, allocatable :: gcImg(:)
  integer :: itemCount
  type(ESMF_GridComp) :: fineGC
! End coarse GC

! Begin
!------

    Iam = "CoarseInitialize"
    call ESMF_GridCompGet( GC, name=COMP_NAME, CONFIG=CF, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

! Begin coarse GC
! Retrieve fine GC
! ---------------------------------
    call ESMF_AttributeGet(GC, name='GC_IMAGE', itemCount=itemCount, rc=status)
    VERIFY_(STATUS)
    allocate(gcImg(itemCount), stat=status)
    VERIFY_(STATUS)
    call ESMF_AttributeGet(GC, name='GC_IMAGE', valueList=gcImg, rc=status)
    VERIFY_(STATUS)
    fineGC = transfer(gcImg, fineGC)
    deallocate(gcImg, stat=status)
    VERIFY_(STATUS)
! End coarse GC

! Retrieve the pointer to the state
! ---------------------------------

    call MAPL_GetObjectFromGC (fineGC, MAPL,  RC=STATUS )
    VERIFY_(STATUS)

! Get the private internal state
!-------------------------------

    call ESMF_UserCompGetInternalState(gc, 'DYNstate', wrap, status)
    VERIFY_(STATUS)
    state => wrap%dyn_state

    DycoreGrid  => state%grid   ! direct handle to grid

! move DynSetup here to be able to use the coarse GC
    !call DynSetup(GC, state, rc=status)
    !call DynSetup(GC, rc=status)
    !VERIFY_(STATUS)

! Get file names from the configuration
!--------------------------------------

!BOR
! !RESOURCE_ITEM: none :: name of layout file
    call MAPL_GetResource ( MAPL, layout_file, 'LAYOUT:', default='fvcore_layout.rc', rc=status )
!EOR
    VERIFY_(STATUS)

    call MAPL_GetResource ( MAPL, DO_ADD_INCS, 'DO_ADD_INCS:', default=DO_ADD_INCS, rc=status )
    VERIFY_(STATUS)

! Check for ColdStart from the configuration
!--------------------------------------
    call MAPL_GetResource ( MAPL, ColdRestart, 'COLDSTART:', default=0, rc=status )
    VERIFY_(STATUS)
    if (ColdRestart /=0 ) then
      call Coldstart( gc, import, export, clock, rc=STATUS )
      VERIFY_(STATUS)
    endif

    !call MAPL_TimerOn(MAPL,"-DYN_INIT")
    call DynInit ( STATE, CLOCK, INTERNAL, IMPORT, GC, status)
    VERIFY_(STATUS)
    !call MAPL_TimerOff(MAPL,"-DYN_INIT")

! Create A-Grid Winds
! -------------------
    ifirst = state%grid%is
    ilast  = state%grid%ie
    jfirst = state%grid%js
    jlast  = state%grid%je
    km     = state%grid%npz

    allocate( UR(ifirst:ilast,jfirst:jlast,km) )
    allocate( VR(ifirst:ilast,jfirst:jlast,km) )

    call getAllWinds( state%vars%u, state%vars%v, UR=UR, VR=VR)

      !U = UR
      !V = VR
      !T = PT*PK
    !PLE = PE

    if(.not.associated(temp3d)) then
       allocate(temp3d(ifirst:ilast,jfirst:jlast,1:km), stat=status)
       VERIFY_(STATUS)
    endif
    temp3d = UR
    call SSI_CopyCoarseToFine(export, temp3d, 'U', STATE%f2c_SSI_arr_map, rc=status)
    VERIFY_(STATUS)
    temp3d = VR
    call SSI_CopyCoarseToFine(export, temp3d, 'V', STATE%f2c_SSI_arr_map, rc=status)
    VERIFY_(STATUS)
    temp3d = state%vars%pt * state%vars%pkz
    call SSI_CopyCoarseToFine(export, temp3d, 'T', STATE%f2c_SSI_arr_map, rc=status)
    VERIFY_(STATUS)
    if (associated(temp3d)) then
       deallocate(temp3d)
    endif
    if(.not.associated(temp3d)) then
       allocate(temp3d(ifirst:ilast,jfirst:jlast,1:km+1), stat=status)
       VERIFY_(STATUS)
    endif
    temp3d = state%vars%pe
    call SSI_CopyCoarseToFine(export, temp3d, 'PLE', STATE%f2c_SSI_arr_map, rc=status)
    VERIFY_(STATUS)

    deallocate( UR )
    deallocate( VR )

! Fill Grid-Cell Area Delta-X/Y
! -----------------------------

    if(.not.associated(temp2d)) then
       allocate(temp2d(ifirst:ilast,jfirst:jlast), stat=status)
       VERIFY_(STATUS)
    endif

    temp2d = DycoreGrid%dxc
    call SSI_CopyCoarseToFine(export, temp2d, 'DXC', STATE%f2c_SSI_arr_map, rc=status)
    VERIFY_(STATUS)

    temp2d = DycoreGrid%dyc
    call SSI_CopyCoarseToFine(export, temp2d, 'DYC', STATE%f2c_SSI_arr_map, rc=status)
    VERIFY_(STATUS)

    temp2d = DycoreGrid%area
    call SSI_CopyCoarseToFine(export, temp2d, 'AREA', STATE%f2c_SSI_arr_map, rc=status)
    VERIFY_(STATUS)

!=====Begin intemittent replay=======================

! Set the intermittent replay alarm, if needed.
! Note that it is a non-sticky alarm
! and is set to ringing on first step. So it will
! work whether the clock is backed-up and ticked
! or not.

    call MAPL_GetResource(MAPL, ReplayMode, 'REPLAY_MODE:', default="NoReplay", RC=STATUS )
    VERIFY_(STATUS)

    if(adjustl(ReplayMode)=="Intermittent") then
       call MAPL_GetResource(MAPL, DNS_INTERVAL,'REPLAY_INTERVAL:', default=21600., RC=STATUS )
       VERIFY_(STATUS)
       call ESMF_TimeIntervalSet(Intv, S=nint(DNS_INTERVAL), RC=STATUS)
       VERIFY_(STATUS)

       ALARM = ESMF_AlarmCreate(name='INTERMITTENT', clock=CLOCK,      &
                                ringInterval=Intv, sticky=.false.,    &
                                                            RC=STATUS )
       VERIFY_(STATUS)
       call ESMF_AlarmRingerOn(ALARM, rc=status)
       VERIFY_(STATUS)
    end if

!========End intermittent replay========================

    !call MAPL_TimerOff(MAPL,"INITIALIZE")
    !call MAPL_TimerOff(MAPL,"TOTAL")

    RETURN_(ESMF_SUCCESS)
  end subroutine Initialize



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!BOP

! !IROUTINE: Run

! !DESCRIPTION: This is the first Run stage of FV. It is the container
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
!
!
! !INTERFACE:

subroutine Run(gc, import, export, clock, rc)

! !ARGUMENTS:

  type(ESMF_GridComp) :: gc
  type (ESMF_State) :: import
  type (ESMF_State) :: export
  type (ESMF_Clock) :: clock
  integer, intent(out)   :: rc

!EOP

! !Local Variables:

    integer                                          :: status
    type (ESMF_FieldBundle)                          :: bundle
    type (ESMF_FieldBundle)                          :: ANA_Bundle
    type (ESMF_Field)                                :: field
    type (ESMF_Field)                                :: ANA_field
    type (ESMF_Config)                               :: cf
    type (ESMF_Alarm)                                :: Alarm
    type (ESMF_Grid)                                 :: ESMFGRID
    type (ESMF_Grid)                                 :: ANAgrid
    type (ESMF_Time)                                 :: currentTime
    type (ESMF_Time)                                 :: RefTime
    class (AbstractRegridder), pointer :: L2C
    class (AbstractRegridder), pointer :: C2L

    type (MAPL_MetaComp), pointer :: mapl

    type (DYN_wrap) :: wrap
    type (DynState), pointer :: STATE
    type (DynGrid),  pointer :: GRID
    type (DynVars),  pointer :: VARS

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
    real(kind=4), pointer ::   phis(:,:) => Null()

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

    real(r8), allocatable, target :: ke    (:,:,:) ! Kinetic    Energy

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

    real(kind=4), allocatable :: dummy (:,:,:) ! Dummy 3-D  Variable
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

    real(kind=4), allocatable ::      dudtana(:,:,:)
    real(kind=4), allocatable ::      dvdtana(:,:,:)
    real(kind=4), allocatable ::      dtdtana(:,:,:)
    real(kind=4), allocatable ::     ddpdtana(:,:,:)
    real(kind=4), pointer ::       qctmp (:,:,:)
    real(kind=4), allocatable ::       dqldt (:,:,:)
    real(kind=4), allocatable ::       dqidt (:,:,:)
    real(kind=4), allocatable ::       doxdt (:,:,:)
    real(kind=4), pointer ::      dqvana (:,:,:) => Null()
    real(kind=4), pointer ::      dqlana (:,:,:) => Null()
    real(kind=4), pointer ::      dqiana (:,:,:) => Null()
    real(kind=4), pointer ::      dqrana (:,:,:) => Null()
    real(kind=4), pointer ::      dqsana (:,:,:) => Null()
    real(kind=4), pointer ::      dqgana (:,:,:) => Null()
    real(kind=4), pointer ::      doxana (:,:,:) => Null()
    real(kind=4), pointer ::       temp3d(:,:,:)
    real(kind=4), pointer ::       vtmp3d(:,:,:)
    real(kind=4), pointer ::         area(:,:)
    real(kind=4), pointer ::       temp2d(:,:)
    real(kind=4), pointer ::       dummy2d(:,:) => Null()
    real(kind=4), pointer ::       dummy3d(:,:,:) => Null()
    real(kind=4), pointer ::       dummy3d_kmplus1(:,:,:) => Null()
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
    character(len=ESMF_MAXSTR) :: COMP_NAME
    character(len=ESMF_MAXSTR) :: STRING
    character(len=ESMF_MAXSTR) :: ReplayFile
    character(len=ESMF_MAXSTR) :: ReplayType
    character(len=ESMF_MAXSTR) :: ReplayMode
    character(len=ESMF_MAXSTR) :: cremap,tremap
    character(len=ESMF_MAXSTR) :: uname,vname,tname,qname,psname,dpname,o3name,rgrid,tvar

    type (MAPL_SunOrbit)       :: ORBIT
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
    logical, save                       :: firstime_tracer_alloc=.true.
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

    logical  :: uphi_associated=.false., vphi_associated=.false.
    logical  :: uke_associated=.false. , vke_associated=.false.

    integer, allocatable :: gcImg(:)
    integer :: itemCount
    type(ESMF_GridComp) :: fineGC

! Retrieve the pointer to the internal state
! ------------------------------------------

  call ESMF_UserCompGetInternalState(gc, 'DYNstate', wrap, status)
  VERIFY_(STATUS)
  state => wrap%dyn_state

  Iam = "CoarseRun"
  call ESMF_GridCompGet( GC, name=COMP_NAME, CONFIG=CF, grid=ESMFGRID, RC=STATUS )
  VERIFY_(STATUS)
  Iam = trim(COMP_NAME) // trim(Iam)

  call ESMF_GridValidate(ESMFGRID,RC=STATUS)
  VERIFY_(STATUS)

! Retrieve fine GC
! ---------------------------------
  call ESMF_AttributeGet(GC, name='GC_IMAGE', itemCount=itemCount, rc=status)
  VERIFY_(STATUS)
  allocate(gcImg(itemCount), stat=status)
  VERIFY_(STATUS)
  call ESMF_AttributeGet(GC, name='GC_IMAGE', valueList=gcImg, rc=status)
  VERIFY_(STATUS)
  fineGC = transfer(gcImg, fineGC)
  deallocate(gcImg, stat=status)
  VERIFY_(STATUS)

! Retrieve the pointer to the generic state
! -----------------------------------------

  call MAPL_GetObjectFromGC (fineGC, MAPL,  RC=STATUS )
  VERIFY_(STATUS)

  call INTERNAL_FineToCoarse(STATE, internal, rc=status)
  VERIFY_(status)

  !call MAPL_TimerOn(MAPL,"TOTAL")
  !call MAPL_TimerOn(MAPL,"RUN")

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

  is_ringing = ESMF_AlarmIsRinging( STATE%ALARMS(TIME_TO_RUN),rc=status); VERIFY_(status)
  if (.not. is_ringing) return

! Allocate Arrays
! ---------------
      ALLOCATE(  dudtana(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(  dvdtana(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(  dtdtana(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(  ddpdtana(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(  dqldt(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(  dqidt(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(  doxdt(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(  dummy(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(delpold(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(  qdnew(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(  qdold(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(  qvold(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(  qlold(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(  qiold(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(  qrold(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(  qsold(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(  qgold(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(     ke(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
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
      call MAPL_GetPointer(export,temp2D,'KEANA',rc=status); VERIFY_(STATUS)
      if(associated(temp2D)) doEnergetics=.true.
      call MAPL_GetPointer(export,temp2D,'PEANA',rc=status); VERIFY_(STATUS)
      if(associated(temp2D)) doEnergetics=.true.
      call MAPL_GetPointer(export,temp2D,'TEANA',rc=status); VERIFY_(STATUS)
      if(associated(temp2D)) doEnergetics=.true.
      call MAPL_GetPointer(export,temp2D,'KEDYN',rc=status); VERIFY_(STATUS)
      if(associated(temp2D)) doEnergetics=.true.
      call MAPL_GetPointer(export,temp2D,'PEDYN',rc=status); VERIFY_(STATUS)
      if(associated(temp2D)) doEnergetics=.true.
      call MAPL_GetPointer(export,temp2D,'TEDYN',rc=status); VERIFY_(STATUS)
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
      ALLOCATE( omaxyz   (ifirstxy:ilastxy,jfirstxy:jlastxy,km  ) )
      ALLOCATE( epvxyz   (ifirstxy:ilastxy,jfirstxy:jlastxy,km  ) )
      ALLOCATE(  cxxyz   (ifirstxy:ilastxy,jfirstxy:jlastxy,km  ) )
      ALLOCATE(  cyxyz   (ifirstxy:ilastxy,jfirstxy:jlastxy,km  ) )
      ALLOCATE( mfxxyz   (ifirstxy:ilastxy,jfirstxy:jlastxy,km  ) )
      ALLOCATE( mfyxyz   (ifirstxy:ilastxy,jfirstxy:jlastxy,km  ) )
      ALLOCATE( mfzxyz   (ifirstxy:ilastxy,jfirstxy:jlastxy,km+1) )

! Report advected friendlies
!---------------------------

      call ESMF_StateGet ( IMPORT, 'TRADV' , BUNDLE,   RC=STATUS )
      VERIFY_(STATUS)

      !-------------------------------------------------------------------
      ! ALT: this section attempts to limit the amount of advected tracers
      !-------------------------------------------------------------------
      adjustTracers = .false.
      call MAPL_GetResource ( MAPL, adjustTracerMode, &
           'EXCLUDE_ADVECTION_TRACERS:', &
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
      else
         call WRITE_PARALLEL('Invalid option, ignored')
         adjustTracers = .false.
      end if
      if (adjustTracers) then
         if (firstime) then
            firstime=.false.
            ! get the list of excluded tracers from resource
            n = 0
            call ESMF_ConfigFindLabel ( CF,'EXCLUDE_ADVECTION_TRACERS_LIST:',isPresent=isPresent,rc=STATUS )
            VERIFY_(STATUS)
            if(isPresent) then

               tend  = .false.
               allocate(xlist(XLIST_MAX), stat=status)
               VERIFY_(STATUS)
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
            call ESMF_FieldBundleGet(BUNDLE, grid=bgrid,fieldCount=nqt,  RC=STATUS)
            VERIFY_(STATUS)
            BundleAdv = ESMF_FieldBundleCreate ( name='xTRADV', rc=STATUS )
            VERIFY_(STATUS)
            call ESMF_FieldBundleSet ( BundleAdv, grid=bgrid, rc=STATUS )
            VERIFY_(STATUS)
            !loop over NQ in TRADV
            do i = 1, nqt
               !get field from TRADV and its name
               call ESMF_FieldBundleGet(bundle, fieldIndex=i, field=field, rc=status)
               VERIFY_(STATUS)
               call ESMF_FieldGet(FIELD, name=fieldname, RC=STATUS)
               VERIFY_(STATUS)
               !exclude everything that is not cloud/water species
               if ( (AdvCore_Advection >= 1    ) .and. &
                 (TRIM(fieldname) /= 'Q'       ) .and. &
                 (TRIM(fieldname) /= 'QLCN'    ) .and. &
                 (TRIM(fieldname) /= 'QLLS'    ) .and. &
                 (TRIM(fieldname) /= 'QICN'    ) .and. &
                 (TRIM(fieldname) /= 'QILS'    ) .and. &
                 (TRIM(fieldname) /= 'CLCN'    ) .and. &
                 (TRIM(fieldname) /= 'CLLS'    ) .and. &
                 (TRIM(fieldname) /= 'NCPL'    ) .and. &
                 (TRIM(fieldname) /= 'NCPI'    ) .and. &
                 (TRIM(fieldname) /= 'QRAIN'   ) .and. &
                 (TRIM(fieldname) /= 'QSNOW'   ) .and. &
                 (TRIM(fieldname) /= 'QGRAUPEL') ) then
                    write(STRING,'(A,A)') "FV3+ADV is excluding ", TRIM(fieldname)
                    call WRITE_PARALLEL( trim(STRING)   )
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
           !   ! Just in case xlist was allocated, but nothing was in it, could have garbage
           !   if (n > 0) then
           !      call ESMF_FieldBundleRemove(BUNDLE, fieldNameList=xlist, &
           !         relaxedFlag=.true., rc=status)
           !      VERIFY_(STATUS)
           !   end if
               deallocate(xlist)
            end if

         end if ! firstime
         BUNDLE = bundleAdv ! replace TRADV
      else
         bundleAdv = BUNDLE ! replace with TRADV
      end if ! adjustTracers

      call ESMF_FieldBundleGet ( BUNDLE, fieldCount=NQ, RC=STATUS )
      VERIFY_(STATUS)

!AOO move tracer allocation to here
      if (firstime_tracer_alloc) then
         firstime_tracer_alloc = .false.
         call allocateTracers(state, import, rc=status)
         VERIFY_(STATUS)
      endif

      if (NQ > 0) then
        allocate( NAMES(NQ),STAT=STATUS )
        VERIFY_(STATUS)
        call ESMF_FieldBundleGet ( BUNDLE, itemorderflag=ESMF_ITEMORDER_ADDORDER, fieldNameList=NAMES, rc=STATUS )
        VERIFY_(STATUS)
        if( .not.allocated( names0 ) ) then
             allocate( NAMES0(NQ),STAT=STATUS )
             VERIFY_(STATUS)
             NAMES0 = NAMES
        endif
      endif

! Surface Geopotential from IMPORT state
!---------------------------------------

      !call MAPL_GetPointer ( IMPORT, PHIS, 'PHIS', RC=STATUS )
      !VERIFY_(STATUS)
      if(.not.associated(phis))  then
         allocate(phis(ifirstxy:ilastxy,jfirstxy:jlastxy), stat=status)
         VERIFY_(STATUS)
      endif
      call SSI_CopyFineToCoarse(import, phis, 'PHIS', STATE%f2c_SSI_arr_map, _RC)

      phisxy = real(phis,kind=r8)

! Get tracers from IMPORT State (Note: Contains Updates from Analysis)
!---------------------------------------------------------------------
      call PULL_Q ( STATE, IMPORT, qqq, NXQ, RC=rc )

!-----------------------------
! end of fewer_tracers-section
!-----------------------------

      do k=1,size(names)
         pos = index(names(k),'::')
         if(pos > 0) then
           if( (names(k)(pos+2:))=='OX' ) ooo = vars%tracer(k)
         elseif(names(k)=='Q') then
            qqq = vars%tracer(k)
         end if
      end do

! WMP Begin REPLAY/ANA section
!-----------------------------

      call ESMF_ConfigGetAttribute ( CF, FV3_STANDALONE, Label="FV3_STANDALONE:", default=0, RC=STATUS)
      VERIFY_(STATUS)
      if (FV3_STANDALONE == 0) then
         !call MAPL_TimerOn(MAPL,"-DYN_ANA")
         call ESMF_ClockGetAlarm(Clock,'ReplayShutOff',Alarm,rc=Status)
         VERIFY_(status)
         is_shutoff = ESMF_AlarmIsRinging( Alarm,rc=Status)
         VERIFY_(status)
      else
         is_shutoff = .true.
      end if

      if (.not. is_shutoff) then
! If requested, do Intermittent Replay
!-------------------------------------

      call MAPL_GetResource(MAPL, ReplayMode, 'REPLAY_MODE:', default="NoReplay", RC=STATUS )
      VERIFY_(STATUS)

      REPLAYING: if(adjustl(ReplayMode)=="Intermittent") then

! If replay alarm is ringing, we need to reset state
!---------------------------------------------------

         call ESMF_ClockGetAlarm(Clock,'INTERMITTENT',Alarm,rc=Status)
         VERIFY_(status)
         call ESMF_ClockGet(Clock, CurrTime=currentTIME, rc=status)
         VERIFY_(status)

         is_ringing = ESMF_AlarmIsRinging( Alarm,rc=status )
         VERIFY_(status)

         RefTime = currentTime

         call check_replay_time_(is_ringing)
         TIME_TO_REPLAY: if(is_ringing) then

            call ESMF_AlarmRingerOff(Alarm, __RC__)

!           Read in file name of field to replay to and all other relavant resources
!           ------------------------------------------------------------------------
            call MAPL_GetResource ( MAPL,ReplayFile,'REPLAY_FILE:', RC=STATUS )
            VERIFY_(status)
            call MAPL_GetResource ( MAPL,ReplayType,'REPLAY_TYPE:', Default="FULL", RC=STATUS )
            VERIFY_(status)

            call MAPL_GetResource ( MAPL, im_replay, Label="REPLAY_IM:", RC=status )
            VERIFY_(STATUS)
            call MAPL_GetResource ( MAPL, jm_replay, Label="REPLAY_JM:", RC=status )
            VERIFY_(STATUS)

            call MAPL_GetResource ( MAPL, psname, Label="REPLAY_PSNAME:", Default="NULL",  RC=status )
            VERIFY_(STATUS)
            call MAPL_GetResource ( MAPL, dpname, Label="REPLAY_DPNAME:", Default="delp",  RC=status )
            VERIFY_(STATUS)
            call MAPL_GetResource ( MAPL,  uname, Label="REPLAY_UNAME:", Default="uwnd",   RC=status )
            VERIFY_(STATUS)
            call MAPL_GetResource ( MAPL,  vname, Label="REPLAY_VNAME:", Default="vwnd",   RC=status )
            VERIFY_(STATUS)
            call MAPL_GetResource ( MAPL,  tname, Label="REPLAY_TNAME:", Default="theta",  RC=status )
            VERIFY_(STATUS)
            call MAPL_GetResource ( MAPL,  qname, Label="REPLAY_QNAME:", Default="sphu",   RC=status )
            VERIFY_(STATUS)
            call MAPL_GetResource ( MAPL, o3name, Label="REPLAY_O3NAME:", Default="ozone", RC=status )
            VERIFY_(STATUS)

            call MAPL_GetResource ( MAPL,  rgrid, Label="REPLAY_GRID:", Default="D-GRID", RC=status )
            VERIFY_(STATUS)
            call MAPL_GetResource ( MAPL,   tvar, Label="REPLAY_TVAR:", Default="THETAV", RC=status )
            VERIFY_(STATUS)

            call MAPL_GetResource ( MAPL, CREMAP, LABEL="REPLAY_REMAP:", default="no", RC=status )
            VERIFY_(STATUS)
            call MAPL_GetResource ( MAPL, TREMAP, LABEL="REPLAY_REMAP_ALL_TRACERS:", default="yes", RC=status )
            VERIFY_(STATUS)

            call MAPL_GetResource ( MAPL, rc_blend,         'REPLAY_BLEND:', default=  0  , RC=STATUS )
            VERIFY_(STATUS)
            call MAPL_GetResource ( MAPL, rc_blend_p_above, 'REPLAY_BLEND_P_ABOVE:', default= 10.0, RC=STATUS )
            VERIFY_(STATUS)
            call MAPL_GetResource ( MAPL, rc_blend_p_below, 'REPLAY_BLEND_P_BELOW:', default=100.0, RC=STATUS )
            VERIFY_(STATUS)

            call MAPL_GetResource ( MAPL, sclinc, label ='SCLINC:', default=1.0, RC=STATUS )


           ! Read the fields to be reset into a bundle
           !------------------------------------------

            call ESMF_ConfigGetAttribute( CF, nx_ana, label ='NX:', rc = STATUS )
            VERIFY_(STATUS)
            call ESMF_ConfigGetAttribute( CF, ny_ana, label ='NY:', rc = STATUS )
            VERIFY_(STATUS)

            block
              use MAPL_LatLonGridFactoryMod

              ANAgrid = grid_manager%make_grid( &
                   & LatLonGridFactory(im_world=IM_REPLAY, jm_world=JM_REPLAY, lm=km, &
                   & nx=nx_ana, ny=ny_ana, rc=status))
              VERIFY_(STATUS)
            end block

            ANA_Bundle = ESMF_FieldBundleCreate( RC=STATUS)
            VERIFY_(STATUS)
            call ESMF_FieldBundleSet(ANA_Bundle, grid=ANAGRID, RC=STATUS)
            VERIFY_(STATUS)

            call MAPL_CFIORead(ReplayFile, RefTime, ANA_Bundle, RC=STATUS)
            VERIFY_(STATUS)

!           Create transform from lat-lon to cubed
!           --------------------------------------
            l2c => regridder_manager%make_regridder(ANAGrid, ESMFGRID, REGRID_METHOD_BILINEAR, RC=STATUS)
            VERIFY_(STATUS)

!           Fill the state variables from the bundle only if
!           the corresponding fields are there
!           -------------------------------------------------

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

!           Done with replay; clean-up
!           --------------------------

            call ESMF_FieldBundleGet(ANA_Bundle , FieldCount=NUMVARS,      RC=STATUS)
            VERIFY_(STATUS)

            do k=1,NUMVARS
               call ESMF_FieldBundleGet (ANA_Bundle, k, ANA_FIELD,    RC=STATUS)
               VERIFY_(STATUS)
               call MAPL_FieldDestroy   (ANA_Field,                   RC=STATUS)
               VERIFY_(STATUS)
            end do

            call ESMF_FieldBundleDestroy(ANA_Bundle,                       RC=STATUS)
            VERIFY_(STATUS)


            end if TIME_TO_REPLAY
      end if REPLAYING

! Create Local Copy of QV and OX (Contains Updates from Analysis)
!----------------------------------------------------------------

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
                   _ASSERT(all(qv >= 0.0),'negative water vapor detected')
                endif
             elseif (associated(qqq%content)) then
                if (size(qv)==size(qqq%content)) then
                   qv = qqq%content
                   _ASSERT(all(qv >= 0.0),'negative water vapor detected')
                endif
             endif
         endif

       enddo
    endif

! Diagnostics Before Analysis Increments are Added
!-------------------------------------------------

      !call MAPL_GetPointer ( IMPORT, dqvana, 'DQVANA', RC=STATUS )   ! Get QV Increment from Analysis
      !VERIFY_(STATUS)
      !call MAPL_GetPointer ( IMPORT, dqlana, 'DQLANA', RC=STATUS )   ! Get QL Increment from Analysis
      !VERIFY_(STATUS)
      !call MAPL_GetPointer ( IMPORT, dqiana, 'DQIANA', RC=STATUS )   ! Get QI Increment from Analysis
      !VERIFY_(STATUS)
      !call MAPL_GetPointer ( IMPORT, dqrana, 'DQRANA', RC=STATUS )   ! Get QR Increment from Analysis
      !VERIFY_(STATUS)
      !call MAPL_GetPointer ( IMPORT, dqsana, 'DQSANA', RC=STATUS )   ! Get QS Increment from Analysis
      !VERIFY_(STATUS)
      !call MAPL_GetPointer ( IMPORT, dqgana, 'DQGANA', RC=STATUS )   ! Get QG Increment from Analysis
      !VERIFY_(STATUS)
      !call MAPL_GetPointer ( IMPORT, doxana, 'DOXANA', RC=STATUS )   ! Get OX Increment from Analysis
      !VERIFY_(STATUS)
      if(.not.associated(dqvana))  then
         allocate(dqvana(ifirstxy:ilastxy,jfirstxy:jlastxy,km), stat=status)
         VERIFY_(STATUS)
      endif
      call SSI_CopyFineToCoarse(import, dqvana, 'DQVANA', STATE%f2c_SSI_arr_map, rc=status)
      VERIFY_(STATUS)

      if(.not.associated(dqlana))  then
         allocate(dqlana(ifirstxy:ilastxy,jfirstxy:jlastxy,km), stat=status)
         VERIFY_(STATUS)
      endif
      call SSI_CopyFineToCoarse(import, dqlana, 'DQLANA', STATE%f2c_SSI_arr_map, rc=status)
      VERIFY_(STATUS)

      if(.not.associated(dqiana))  then
         allocate(dqiana(ifirstxy:ilastxy,jfirstxy:jlastxy,km), stat=status)
         VERIFY_(STATUS)
      endif
      call SSI_CopyFineToCoarse(import, dqiana, 'DQIANA', STATE%f2c_SSI_arr_map, rc=status)
      VERIFY_(STATUS)

      if(.not.associated(dqrana))  then
         allocate(dqrana(ifirstxy:ilastxy,jfirstxy:jlastxy,km), stat=status)
         VERIFY_(STATUS)
      endif
      call SSI_CopyFineToCoarse(import, dqrana, 'DQRANA', STATE%f2c_SSI_arr_map, rc=status)
      VERIFY_(STATUS)

      if(.not.associated(dqsana))  then
         allocate(dqsana(ifirstxy:ilastxy,jfirstxy:jlastxy,km), stat=status)
         VERIFY_(STATUS)
      endif
      call SSI_CopyFineToCoarse(import, dqsana, 'DQSANA', STATE%f2c_SSI_arr_map, rc=status)
      VERIFY_(STATUS)

      if(.not.associated(dqgana))  then
         allocate(dqgana(ifirstxy:ilastxy,jfirstxy:jlastxy,km), stat=status)
         VERIFY_(STATUS)
      endif
      call SSI_CopyFineToCoarse(import, dqgana, 'DQGANA', STATE%f2c_SSI_arr_map, rc=status)
      VERIFY_(STATUS)

      if(.not.associated(doxana))  then
         allocate(doxana(ifirstxy:ilastxy,jfirstxy:jlastxy,km), stat=status)
         VERIFY_(STATUS)
      endif
      call SSI_CopyFineToCoarse(import, doxana, 'DOXANA', STATE%f2c_SSI_arr_map, rc=status)
      VERIFY_(STATUS)

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

!! Get A-grid winds
!! ----------------
      call getAllWinds(vars%u, vars%v, UR=ur, VR=vr)

      delp   = vars%pe(:,:,2:)  -vars%pe(:,:,:km)   ! Pressure Thickness
      dmdt   = vars%pe(:,:,km+1)-vars%pe(:,:,1)     ! Psurf-Ptop
      tempxy = vars%pt * (1.0+eps*(qv-dqvana))       ! Compute THV Before Analysis Update

      if (doEnergetics) &
      call Energetics (ur,vr,tempxy,vars%pe,delp,vars%pkz,phisxy,kenrg,penrg,tenrg)

! Pointers to copy back from coarse to fine as needed
      if(.not.associated(dummy3d)) then
         allocate(dummy3d(ifirstxy:ilastxy,jfirstxy:jlastxy,km), stat=status)
         VERIFY_(STATUS)
      endif
      if(.not.associated(dummy3d_kmplus1)) then
         allocate(dummy3d_kmplus1(ifirstxy:ilastxy,jfirstxy:jlastxy,km+1), stat=status)
         VERIFY_(STATUS)
      endif
      if(.not.associated(dummy2d)) then
         allocate(dummy2d(ifirstxy:ilastxy,jfirstxy:jlastxy), stat=status)
         VERIFY_(STATUS)
      endif
! DUDTANA
! -------
      call MAPL_GetPointer ( export, temp3d, 'DUDTANA', rc=status )
      VERIFY_(STATUS)
      if( associated(temp3d) ) then !dudtana = ur
         dudtana = ur
         dummy3d = dudtana
         call SSI_CopyCoarseToFine(export, dummy3d, 'DUDTANA', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      endif

! DVDTANA
! -------
      call MAPL_GetPointer ( export, temp3d, 'DVDTANA', rc=status )
      VERIFY_(STATUS)
      if( associated(temp3d) ) then !dvdtana = vr
         dvdtana = vr
         dummy3d = dvdtana
         call SSI_CopyCoarseToFine(export, dummy3d, 'DVDTANA', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      endif

! DTDTANA
! -------
      call MAPL_GetPointer ( export, temp3d, 'DTDTANA', rc=status )
      VERIFY_(STATUS)
      if( associated(temp3d) ) then !dtdtana = vars%pt * vars%pkz
         dtdtana = vars%pt * vars%pkz
         dummy3d = dtdtana
         call SSI_CopyCoarseToFine(export, dummy3d, 'DTDTANA', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      endif

! DDELPDTANA
! ----------
      call MAPL_GetPointer ( export, temp3d, 'DDELPDTANA', rc=status )
      VERIFY_(STATUS)
      if( associated(temp3d) ) then !ddpdtana = delp
         ddpdtana = delp
         dummy3d = ddpdtana
         call SSI_CopyCoarseToFine(export, dummy3d, 'DDELPDTANA', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      endif

! DTHVDTANAINT
! ------------
      call MAPL_GetPointer ( export, temp2D, 'DTHVDTANAINT', rc=status )
      VERIFY_(STATUS)
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
! -----------
      call MAPL_GetPointer ( export, temp2D, 'DQVDTANAINT', rc=status )
      VERIFY_(STATUS)
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
! -----------
      call MAPL_GetPointer ( export, temp2D, 'DQLDTANAINT', rc=status )
      VERIFY_(STATUS)
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
! -----------
      call MAPL_GetPointer ( export, temp2D, 'DQIDTANAINT', rc=status )
      VERIFY_(STATUS)
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
! -----------
      call MAPL_GetPointer ( export, temp2D, 'DOXDTANAINT', rc=status )
      VERIFY_(STATUS)
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
! -----------------------------------------------------

      if (vars%nwat >= 6) then
        QDOLD = 1.0 - (QVOLD+QLOLD+QIOLD+QROLD+QSOLD+QGOLD)
        QDNEW = 1.0 - (QV   +QL   +QI   +QR   +QS   +QG   )
      else
        QDOLD = 1.0 - (QVOLD+QLOLD+QIOLD)
        QDNEW = 1.0 - (QV   +QL   +QI   )
      endif
      call MAPL_GetPointer(export, area, 'AREA', rc=status)
      VERIFY_(STATUS)

      allocate( trsum1(nq) )
      allocate( trsum2(nq) )

      call MAPL_GetResource(MAPL, ANA_IS_WEIGHTED, Label="ANA_IS_WEIGHTED:", default='NO', RC=STATUS)
      VERIFY_(STATUS)
           ANA_IS_WEIGHTED = ESMF_UtilStringUpperCase(ANA_IS_WEIGHTED)
               IS_WEIGHTED =   adjustl(ANA_IS_WEIGHTED)=="YES" .or. adjustl(ANA_IS_WEIGHTED)=="NO"
      _ASSERT( IS_WEIGHTED ,'needs informative message')
               IS_WEIGHTED =   adjustl(ANA_IS_WEIGHTED)=="YES"

      ! Add Analysis Tendencies
      ! -----------------------
      delpold = delp                            ! Old Pressure Thickness

      call ADD_INCS ( STATE,IMPORT,internal,DT,IS_WEIGHTED=IS_WEIGHTED )

      if (DYN_DEBUG) call DEBUG_FV_STATE('ANA ADD_INCS',STATE)

      delp = vars%pe(:,:,2:)-vars%pe(:,:,:km)   ! Updated Pressure Thickness

      ! Compute Old Global Sums of Tracers over Locations where Mass has changed
      ! ------------------------------------------------------------------------
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
         call MAPL_AreaMean( TRSUM1(n), qsum2, real(grid%area), esmfgrid, rc=STATUS )
         VERIFY_(STATUS)
      enddo
      endif

      ! Update Specific Mass of Aerosol Constituents Keeping Mixing_Ratio Constant WRT_Dry_Air After ANA Updates
      ! --------------------------------------------------------------------------------------------------------
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
      ! ------------------------------------------------------------------------
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
         call MAPL_AreaMean( TRSUM2(n), qsum2, real(grid%area), esmfgrid, rc=STATUS )
         VERIFY_(STATUS)
      enddo
      endif

      ! Ensure Conservation of Global Mass of Aerosol Constituents After ANA Updates
      ! ----------------------------------------------------------------------------
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
!--------------------------------------------------------------------

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
         endif
      enddo

! Diagnostics After Analysis Increments are Added
!------------------------------------------------

      call MAPL_GetPointer ( export, temp2D, 'DMDTANA', rc=status )
      VERIFY_(STATUS)
      if( associated(temp2D) ) then
         dummy2d = ( (vars%pe(:,:,km+1)-vars%pe(:,:,1)) - dmdt )/(grav*dt)
         call SSI_CopyCoarseToFine(export, temp2D, 'DMDTANA', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      endif

      call getAllWinds(vars%u, vars%v, UC=uc0, VC=vc0, UR=ur, VR=vr)

      dmdt = vars%pe(:,:,km+1)-vars%pe(:,:,1)     ! Psurf-Ptop

! DUDTANA
! -------
      call MAPL_GetPointer ( export, temp3d, 'DUDTANA', rc=status )
      VERIFY_(STATUS)
      if( associated(temp3d) ) then
         dudtana = (ur-dudtana)/dt
         dummy3d = dudtana
         call SSI_CopyCoarseToFine(export, dummy3d, 'DUDTANA', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      endif

! DVDTANA
! -------
      call MAPL_GetPointer ( export, temp3d, 'DVDTANA', rc=status )
      VERIFY_(STATUS)
      if( associated(temp3d) ) then
         dvdtana = (vr-dvdtana)/dt
         dummy3d = dvdtana
         call SSI_CopyCoarseToFine(export, dummy3d, 'DVDTANA', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      endif

! DTDTANA
! -------
      call MAPL_GetPointer ( export, temp3d, 'DTDTANA', rc=status )
      VERIFY_(STATUS)
      if( associated(temp3d) ) then
         dummy   =  vars%pt*vars%pkz
         dtdtana = (dummy-dtdtana)/dt
         dummy3d = dtdtana
         call SSI_CopyCoarseToFine(export, dummy3d, 'DTDTANA', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      endif

! DDELPDTANA
! ----------
      call MAPL_GetPointer ( export, temp3d, 'DDELPDTANA', rc=status )
      VERIFY_(STATUS)
      if( associated(temp3d) ) then
         !dummy   =  delp
         ddpdtana = (delp-ddpdtana)/dt
         dummy3d = ddpdtana
         call SSI_CopyCoarseToFine(export, dummy3d, 'DDELPDTANA', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      endif

! DTHVDTANAINT
! ------------
      call MAPL_GetPointer ( export, temp2D, 'DTHVDTANAINT', rc=status )
      VERIFY_(STATUS)
      if( associated(temp2D) ) then
          tempxy       = vars%pt*(1+eps*qv)   ! Set tempxy = TH*QVnew (After Analysis Update)
          dthdtanaint2 = 0.0
          do k=1,km
          dthdtanaint2 = dthdtanaint2 + tempxy(:,:,k)*delp(:,:,k)
          enddo
          dummy2d       = (dthdtanaint2-dthdtanaint1) * MAPL_P00**MAPL_KAPPA / (MAPL_GRAV*DT)
          call SSI_CopyCoarseToFine(export, dummy2d, 'DTHVDTANAINT', STATE%f2c_SSI_arr_map, _RC)
          DEALLOCATE( dthdtanaint1 )
          DEALLOCATE( dthdtanaint2 )
      endif

! DQVDTANAINT
! -----------
      call MAPL_GetPointer ( export, temp2D, 'DQVDTANAINT', rc=status )
      VERIFY_(STATUS)
      if( associated(temp2D) ) then
          tempxy       = qv         ! Set tempxy = QNEW (After Analysis Update)
          dqvdtanaint2 = 0.0
          do k=1,km
          dqvdtanaint2 = dqvdtanaint2 + tempxy(:,:,k)*delp(:,:,k)
          enddo
          dummy2d       = (dqvdtanaint2-dqvdtanaint1) / (MAPL_GRAV*DT)
          call SSI_CopyCoarseToFine(export, dummy2d, 'DQVDTANAINT', STATE%f2c_SSI_arr_map, _RC)
          DEALLOCATE( dqvdtanaint1 )
          DEALLOCATE( dqvdtanaint2 )
      endif

! DQLDTANAINT
! -----------
      call MAPL_GetPointer ( export, temp2D, 'DQLDTANAINT', rc=status )
      VERIFY_(STATUS)
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
          dummy2d = (dqldtanaint2-dqldtanaint1) / (MAPL_GRAV*DT)
          call SSI_CopyCoarseToFine(export, dummy2d, 'DQLDTANAINT', STATE%f2c_SSI_arr_map, _RC)
          DEALLOCATE( dqldtanaint1 )
          DEALLOCATE( dqldtanaint2 )
      endif

! DQIDTANAINT
! -----------
      call MAPL_GetPointer ( export, temp2D, 'DQIDTANAINT', rc=status )
      VERIFY_(STATUS)
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
          dummy2d = (dqidtanaint2-dqidtanaint1) / (MAPL_GRAV*DT)
          call SSI_CopyCoarseToFine(export, dummy2d, 'DQIDTANAINT', STATE%f2c_SSI_arr_map, _RC)
          DEALLOCATE( dqidtanaint1 )
          DEALLOCATE( dqidtanaint2 )
      endif

! DOXDTANAINT
! -----------
      call MAPL_GetPointer ( export, temp2D, 'DOXDTANAINT', rc=status )
      VERIFY_(STATUS)
      if( associated(temp2D) ) then
          tempxy       = ox         ! Set tempxy = OXnew (After Analysis Update)
          doxdtanaint2 = 0.0
          do k=1,km
          doxdtanaint2 = doxdtanaint2 + tempxy(:,:,k)*delp(:,:,k)
          enddo
          dummy2d = (doxdtanaint2-doxdtanaint1) * (MAPL_O3MW/MAPL_AIRMW) / (MAPL_GRAV*DT)
          call SSI_CopyCoarseToFine(export, dummy2d, 'DOXDTANAINT', STATE%f2c_SSI_arr_map, _RC)
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
                   _ASSERT(all(qv >= 0.0),'negative water vapor detected')
                endif
             elseif (associated(qqq%content)) then
                if (size(qv)==size(qqq%content)) then
                   qv = qqq%content
                   _ASSERT(all(qv >= 0.0),'negative water vapor detected')
                endif
             endif
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
         !call MAPL_TimerOff(MAPL,"-DYN_ANA")
      endif


      !call MAPL_TimerOn(MAPL,"-DYN_PROLOGUE")
! Create FV Thermodynamic Variables
!----------------------------------

      tempxy = vars%pt * vars%pkz      ! Compute Dry Temperature

! Initialize Diagnostic Dynamics Tendencies
! -----------------------------------------

      dpedt  = vars%pe      ! Edge Pressure      Tendency
      ddpdt  =    delp      ! Pressure Thickness Tendency
      dudt   =     ur       ! U-Wind on A-Grid   Tendency
      dvdt   =     vr       ! V-Wind on A-Grid   Tendency
      dtdt   = tempxy       ! Dry Temperature    Tendency
      dqdt   =     qv       ! Specific Humidity  Tendency
      dthdt  = vars%pt*(1.0+eps*qv)*delp

      call FILLOUT3 (export,  'QV_DYN_IN',      qv, STATE, _RC)
      call FILLOUT3 (export,   'T_DYN_IN',  tempxy, STATE, _RC)
      call FILLOUT3 (export,   'U_DYN_IN',      ur, STATE, _RC)
      call FILLOUT3 (export,   'V_DYN_IN',      vr, STATE, _RC)
      call FILLOUT3 (export, 'PLE_DYN_IN', vars%pe, STATE, _RC)

! Initialize 3-D Tracer Dynamics Tendencies
! -----------------------------------------

     if (allocated(names)) then

      call MAPL_GetPointer( export,temp3d,'DQLDTDYN', rc=status )
      VERIFY_(STATUS)
      if( associated(temp3d) ) then
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
          dummy3d = dqldt
          call SSI_CopyCoarseToFine(export, dummy3d, 'DQLDTDYN', STATE%f2c_SSI_arr_map, rc=status)
          VERIFY_(STATUS)
      endif

      call MAPL_GetPointer( export,temp3d,'DQIDTDYN', rc=status )
      VERIFY_(STATUS)
      if( associated(temp3d) ) then
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
          dummy3d = dqidt
          call SSI_CopyCoarseToFine(export, dummy3d, 'DQIDTDYN', STATE%f2c_SSI_arr_map, rc=status)
          VERIFY_(STATUS)
      endif

      call MAPL_GetPointer( export,temp3d,'DOXDTDYN', rc=status )
      VERIFY_(STATUS)
      if( associated(temp3d) ) then
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
          dummy3d = doxdt
          call SSI_CopyCoarseToFine(export, dummy3d, 'DOXDTDYN', STATE%f2c_SSI_arr_map, rc=status)
          VERIFY_(STATUS)
      endif
    endif

! Initialize 2-D Vertically Integrated Tracer Dynamics Tendencies
! ---------------------------------------------------------------

      call MAPL_GetPointer ( export, temp2D, 'DQVDTDYNINT', rc=status )
      VERIFY_(STATUS)
      if( associated(temp2D) ) then
          dummy2d = 0.0
          do k=1,km
          dummy2d = dummy2d - qv(:,:,k)*delp(:,:,k)
          enddo
          call SSI_CopyCoarseToFine(export, dummy2d, 'DQVDTDYNINT', STATE%f2c_SSI_arr_map, rc=status)
          VERIFY_(STATUS)
      endif

      call MAPL_GetPointer ( export, temp2D, 'DQLDTDYNINT', rc=status )
      VERIFY_(STATUS)
      if( associated(temp2D) ) then
          dummy2d = 0.0
          do N = 1,size(names)
             if( trim(names(N)).eq.'QLCN' .or. &
                 trim(names(N)).eq.'QLLS' ) then
                 if( state%vars%tracer(N)%is_r4 ) then
                     do k=1,km
                     dummy2d = dummy2d - state%vars%tracer(N)%content_r4(:,:,k)*delp(:,:,k)
                     enddo
                 else
                     do k=1,km
                     dummy2d = dummy2d - state%vars%tracer(N)%content(:,:,k)*delp(:,:,k)
                     enddo
                 endif
             endif
          enddo
          call SSI_CopyCoarseToFine(export, dummy2d, 'DQLDTDYNINT', STATE%f2c_SSI_arr_map, rc=status)
          VERIFY_(STATUS)
      endif

      call MAPL_GetPointer ( export, temp2D, 'DQIDTDYNINT', rc=status )
      VERIFY_(STATUS)
      if( associated(temp2D) ) then
          dummy2d = 0.0
          do N = 1,size(names)
             if( trim(names(N)).eq.'QICN' .or. &
                 trim(names(N)).eq.'QILS' ) then
                 if( state%vars%tracer(N)%is_r4 ) then
                     do k=1,km
                     dummy2d = dummy2d - state%vars%tracer(N)%content_r4(:,:,k)*delp(:,:,k)
                     enddo
                 else
                     do k=1,km
                     dummy2d = dummy2d - state%vars%tracer(N)%content(:,:,k)*delp(:,:,k)
                     enddo
                 endif
             endif
          enddo
          call SSI_CopyCoarseToFine(export, dummy2d, 'DQIDTDYNINT', STATE%f2c_SSI_arr_map, rc=status)
          VERIFY_(STATUS)
      endif

      call MAPL_GetPointer ( export, temp2D, 'DOXDTDYNINT', rc=status )
      VERIFY_(STATUS)
      if( associated(temp2D) ) then
          dummy2d = 0.0
          do N = 1,size(names)
             pos = index(names(N),'::')
             if(pos > 0) then
               if( (names(N)(pos+2:))=='OX' ) then
                 if( state%vars%tracer(N)%is_r4 ) then
                     do k=1,km
                     dummy2d = dummy2d - state%vars%tracer(N)%content_r4(:,:,k)*delp(:,:,k)
                     enddo
                 else
                     do k=1,km
                     dummy2d = dummy2d - state%vars%tracer(N)%content(:,:,k)*delp(:,:,k)
                     enddo
                 endif
               endif
             endif
          enddo
          call SSI_CopyCoarseToFine(export, dummy2d, 'DOXDTDYNINT', STATE%f2c_SSI_arr_map, rc=status)
          VERIFY_(STATUS)
      endif

! Compute Energetics After Analysis (and Before Dycore)
! -----------------------------------------------------

    tempxy = vars%pt * (1.0+eps*qv)       ! Compute THV After Analysis Update

    if (doEnergetics) then
    call Energetics (ur,vr,tempxy,vars%pe,delp,vars%pkz,phisxy, kenrg0,penrg0,tenrg0)
    kenrg = (kenrg0-kenrg)/DT
    penrg = (penrg0-penrg)/DT
    tenrg = (tenrg0-tenrg)/DT
    call FILLOUT2 (export, 'KEANA', kenrg, STATE, _RC)
    call FILLOUT2 (export, 'PEANA', penrg, STATE, _RC)
    call FILLOUT2 (export, 'TEANA', tenrg, STATE, _RC)
    endif

! Call Wrapper (DynRun) for FVDycore
! ----------------------------------
      call MAPL_GetResource( MAPL, CONSV, 'CONSV:', default=1, RC=STATUS )
      VERIFY_(STATUS)
      call MAPL_GetResource( MAPL,  FILL,  'FILL:', default=0, RC=STATUS )
      VERIFY_(STATUS)

      LCONSV = CONSV.eq.1
      LFILL  =  FILL.eq.1

! Fill pressures before dynamics export
!-------------------------------------------------------
      pe0=vars%pe
      call FILLOUT3r8 (export, 'PLE0', pe0, STATE, _RC)

      !call MAPL_TimerOff(MAPL,"-DYN_PROLOGUE")

!-------------------------------------------------------

      !call MAPL_TimerOn(MAPL,"-DYN_CORE")
      t1 = MPI_Wtime(status)
      call DynRun (GC, STATE, EXPORT, CLOCK, internal, import, RC=STATUS)
      !call DynRun (GC, STATE, EXPORT, CLOCK, internal, RC=STATUS)
      VERIFY_(STATUS)
      t2 = MPI_Wtime(status)
      dyn_run_timer = t2-t1
      !call MAPL_TimerOff(MAPL,"-DYN_CORE")

      !call MAPL_TimerOn(MAPL,"-DYN_EPILOGUE")
! Computational diagnostics
! --------------------------
    call MAPL_GetPointer(export,temp2d,'DYNTIMER',rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       dummy2d = dyn_run_timer
       call SSI_CopyCoarseToFine(export, dummy2d, 'DYNTIMER', STATE%f2c_SSI_arr_map, rc=status)
       VERIFY_(STATUS)
    endif
    call MAPL_GetPointer(export,temp2d,'PID',rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       dummy2d = 0 !WMP need to get from MAPL gid
       call SSI_CopyCoarseToFine(export, dummy2d, 'PID', STATE%f2c_SSI_arr_map, rc=status)
       VERIFY_(STATUS)
    endif

!#define DEBUG_WINDS
#if defined(DEBUG_WINDS)
  call Write_Profile(grid, vars%u, 'U-after-DynRun')
  call Write_Profile(grid, vars%v, 'V-after-DynRun')
#endif
      plk  = exp( kappa * log( 0.5*(vars%pe(:,:,1:km)+vars%pe(:,:,2:km+1)) ) )
      pkxy = exp( kappa * log( vars%pe ) )

!----------------------------------------------------------------------------

    if (SW_DYNAMICS) then

      call MAPL_GetPointer(export,temp2d,'PHIS', rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
          dummy2d = phisxy
          call SSI_CopyCoarseToFine(export, dummy2d, 'PHIS', STATE%f2c_SSI_arr_map, rc=status)
          VERIFY_(STATUS)
      endif

      call MAPL_GetPointer(export,temp2d,'PS',  rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
          dummy2d =  vars%pe(:,:,km+1)/GRAV
          call SSI_CopyCoarseToFine(export, dummy2d, 'PS', STATE%f2c_SSI_arr_map, rc=status)
          VERIFY_(STATUS)
      endif

      call getAllWinds(vars%u, vars%v, UA=ua, VA=va, UC=uc, VC=vc, UR=ur, VR=vr)
      call FILLOUT3 (export, 'U_DGRID', vars%u  , STATE, rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'V_DGRID', vars%v  , STATE, rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'U_CGRID', uc      , STATE, rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'V_CGRID', vc      , STATE, rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'U_AGRID', ua      , STATE, rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'V_AGRID', va      , STATE, rc=status); VERIFY_(STATUS)

      call FILLOUT3 (export, 'U'      , ur      , STATE, rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'V'      , vr      , STATE, rc=status); VERIFY_(STATUS)


    else               ! .not. SW_DYNAMICS

! Load Local Variable with Vapor Specific Humidity
! ------------------------------------------------

    if ((.not. ADIABATIC) .and. (STATE%GRID%NQ > 0)) then
      if ( qqq%is_r4 ) then
         if (size(qv)==size(qqq%content_r4)) qv = qqq%content_r4
      else
         if (size(qv)==size(qqq%content)   ) qv = qqq%content
      endif
    else
      qv = 0.0
    endif

! Vertically Integrated THV Tendency Diagnostic
! ---------------------------------------------
      delp  = ( vars%pe(:,:,2:) - vars%pe(:,:,:km) )
      dthdt = ( vars%pt*(1.0+eps*qv)*delp-dthdt )/dt

      call MAPL_GetPointer(export,temp2d,'DTHVDTDYNINT', rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         qsum1 = 0.0
      do k=1,km
         qsum1 = qsum1 + dthdt(:,:,k)
      enddo
         dummy2d = qsum1 * (MAPL_P00**MAPL_KAPPA) / grav
         call SSI_CopyCoarseToFine(export, dummy2d, 'DTHVDTDYNINT', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      end if

! Compute Dry Theta and T with Unified Poles

      tempxy  = vars%pt * vars%pkz

! Compute Mid-Layer Pressure and Pressure Thickness
! -------------------------------------------------

      delp = ( vars%pe(:,:,2:) - vars%pe(:,:,:km) )
      pl   = ( vars%pe(:,:,2:) + vars%pe(:,:,:km) ) * 0.5

! Get all wind derivatives
! ------------------------
      call getAllWinds(vars%u, vars%v, UA=ua, VA=va, UC=uc, VC=vc, UR=ur, VR=vr)
      call getVorticity(vars%u, vars%v, vort)
      call getDivergence(uc, vc, divg)

! Compute absolute vorticity on the D grid
! -------------------------------------------------
      call getEPV(vars%pt,vort,ua,va,epvxyz)
      call MAPL_GetPointer(export, temp3D, 'EPV', rc=status)
      VERIFY_(STATUS)
      if(associated(temp3d)) then
         dummy3d = epvxyz*(p00**kappa)
         call SSI_CopyCoarseToFine(export, dummy3d, 'EPV', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      end if

! Compute Tropopause Pressure, Temperature, and Moisture
! ------------------------------------------------------

      doTropvars=.false.
      call MAPL_GetPointer(export,temp2D,'TROPP_THERMAL',rc=status); VERIFY_(STATUS)
      if(associated(temp2D)) doTropvars=.true.
      call MAPL_GetPointer(export,temp2D,'TROPP_EPV',rc=status); VERIFY_(STATUS)
      if(associated(temp2D)) doTropvars=.true.
      call MAPL_GetPointer(export,temp2D,'TROPP_BLENDED',rc=status); VERIFY_(STATUS)
      if(associated(temp2D)) doTropvars=.true.
      call MAPL_GetPointer(export,temp2D,'TROPK_BLENDED',rc=status); VERIFY_(STATUS)
      if(associated(temp2D)) doTropvars=.true.
      call MAPL_GetPointer(export,temp2D,'TROPT',rc=status); VERIFY_(STATUS)
      if(associated(temp2D)) doTropvars=.true.
      call MAPL_GetPointer(export,temp2D,'TROPQ',rc=status); VERIFY_(STATUS)
      if(associated(temp2D)) doTropvars=.true.

      if (doTropvars) then
         ALLOCATE( tropp1 (ifirstxy:ilastxy,jfirstxy:jlastxy) )
         ALLOCATE( tropp2 (ifirstxy:ilastxy,jfirstxy:jlastxy) )
         ALLOCATE( tropp3 (ifirstxy:ilastxy,jfirstxy:jlastxy) )
         ALLOCATE( tropt  (ifirstxy:ilastxy,jfirstxy:jlastxy) )
         ALLOCATE( tropq  (ifirstxy:ilastxy,jfirstxy:jlastxy) )
         call tropovars ( ilastxy-ifirstxy+1,jlastxy-jfirstxy+1,km, &
                          real(vars%pe            ,kind=4),         &
                          real(pl                 ,kind=4),         &
                          real(tempxy             ,kind=4),         &
                          real(qv                 ,kind=4),         &
                          real(epvxyz*(p00**kappa),kind=4),         &
                          tropp1,tropp2,tropp3,tropt,tropq          )

         ! get blended index
         call MAPL_GetPointer(export,temp2D,'TROPK_BLENDED',rc=status); VERIFY_(STATUS)
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
                  !dummy2d(i-ifirstxy+1,j-jfirstxy+1) = kend
                  dummy2d(i,j) = kend
               enddo
            enddo
            call SSI_CopyCoarseToFine(export, dummy2d, 'TROPK_BLENDED', STATE%f2c_SSI_arr_map, _RC)
         endif

         call MAPL_GetPointer(export,temp2D,'TROPP_THERMAL',rc=status)
         VERIFY_(STATUS)
         if(associated(temp2D)) then
            dummy2d = tropp1
            call SSI_CopyCoarseToFine(export, dummy2d, 'TROPP_THERMAL', STATE%f2c_SSI_arr_map, rc=status)
            VERIFY_(STATUS)
         end if

         call MAPL_GetPointer(export,temp2D,'TROPP_EPV',rc=status)
         VERIFY_(STATUS)
         if(associated(temp2D)) then
            dummy2d = tropp2
            call SSI_CopyCoarseToFine(export, dummy2d, 'TROPP_EPV', STATE%f2c_SSI_arr_map, rc=status)
            VERIFY_(STATUS)
         end if

         call MAPL_GetPointer(export,temp2D,'TROPP_BLENDED',rc=status)
         VERIFY_(STATUS)
         if(associated(temp2D)) then
            dummy2d = tropp3
            call SSI_CopyCoarseToFine(export, dummy2d, 'TROPP_BLENDED', STATE%f2c_SSI_arr_map, rc=status)
            VERIFY_(STATUS)
         end if

         call MAPL_GetPointer(export,temp2D,'TROPT',rc=status)
         VERIFY_(STATUS)
         if(associated(temp2D)) then
            dummy2d = tropt
            call SSI_CopyCoarseToFine(export, dummy2d, 'TROPT', STATE%f2c_SSI_arr_map, rc=status)
            VERIFY_(STATUS)
         end if

         call MAPL_GetPointer(export,temp2D,'TROPQ',rc=status)
         VERIFY_(STATUS)
         if(associated(temp2D)) then
            dummy2d = tropq
            call SSI_CopyCoarseToFine(export, dummy2d, 'TROPQ', STATE%f2c_SSI_arr_map, rc=status)
            VERIFY_(STATUS)
         end if

         DEALLOCATE( tropp1 )
         DEALLOCATE( tropp2 )
         DEALLOCATE( tropp3 )
         DEALLOCATE( tropt  )
         DEALLOCATE( tropq  )
      endif

! Get Cubed-Sphere Wind Exports
! -----------------------------
      call FILLOUT3 (export, 'U_DGRID', vars%u  , STATE, _RC)
      call FILLOUT3 (export, 'V_DGRID', vars%v  , STATE, _RC)
      call FILLOUT3 (export, 'U_CGRID', uc      , STATE, _RC)
      call FILLOUT3 (export, 'V_CGRID', vc      , STATE, _RC)
      call FILLOUT3 (export, 'U_AGRID', ua      , STATE, _RC)
      call FILLOUT3 (export, 'V_AGRID', va      , STATE, _RC)

! Compute Diagnostic Dynamics Tendencies
!  (Note: initial values of d(m,u,v,T,q)/dt are progs m,u,v,T,q)
! --------------------------------------------------------------

      dmdt = ( vars%pe(:,:,km+1)-vars%pe(:,:,1) - dmdt )/(grav*dt)

      dudt = (    ur-dudt )/dt
      dvdt = (    vr-dvdt )/dt
      dtdt = (  tempxy-dtdt )/dt
      dqdt = (      qv-dqdt )/dt

      dpedt = ( vars%pe - dpedt )/dt
      ddpdt = ( delp - ddpdt )/dt ! Pressure Thickness Tendency


      call FILLOUT3 (export, 'DELP'      ,delp , STATE, _RC)
      call FILLOUT3 (export, 'DUDTDYN'   ,dudt , STATE, _RC)
      call FILLOUT3 (export, 'DVDTDYN'   ,dvdt , STATE, _RC)
      call FILLOUT3 (export, 'DTDTDYN'   ,dtdt , STATE, _RC)
      call FILLOUT3 (export, 'DQVDTDYN'  ,dqdt , STATE, _RC)
      call FILLOUT3 (export, 'DDELPDTDYN',ddpdt, STATE, _RC)
      call FILLOUT3 (export, 'DPLEDTDYN' ,dpedt, STATE, _RC)

      pe1=vars%pe
      call FILLOUT3r8 (export, 'PLE1', pe1    , STATE, _RC)

      if (AdvCore_Advection==2) then
      ! Compute time-centered C-Grid Courant Numbers and Mass Fluxes on Cubed Orientation
        uc0 = 0.5*(uc +uc0)
        vc0 = 0.5*(vc +vc0)
        pe0 = 0.5*(pe1+pe0)
        call computeMassFluxes(uc0, vc0, pe0, mfxxyz, mfyxyz, cxxyz, cyxyz, dt)
        call FILLOUT3r8 (export, 'CX'  , cxxyz  , STATE, _RC)
        call FILLOUT3r8 (export, 'CY'  , cyxyz  , STATE, _RC)
        call FILLOUT3r8 (export, 'MFX' , mfxxyz , STATE, _RC)
        call FILLOUT3r8 (export, 'MFY' , mfyxyz , STATE, _RC)
      else
      ! Fill Advection C-Grid Courant Numbers and Mass Fluxes on Cubed Orientation from FV3 DynCore
        call fillMassFluxes(mfxxyz, mfyxyz, cxxyz, cyxyz)
        call FILLOUT3r8 (export, 'CX'  , cxxyz  , STATE, _RC)
        call FILLOUT3r8 (export, 'CY'  , cyxyz  , STATE, _RC)
        call FILLOUT3r8 (export, 'MFX' , mfxxyz , STATE, _RC)
        call FILLOUT3r8 (export, 'MFY' , mfyxyz , STATE, _RC)
      endif

      call FILLOUT3 (export, 'CU' ,  cxxyz , STATE, _RC)
      call FILLOUT3 (export, 'CV' ,  cyxyz , STATE, _RC)
      call FILLOUT3 (export, 'MX' , mfxxyz , STATE, _RC)
      call FILLOUT3 (export, 'MY' , mfyxyz , STATE, _RC)

! Compute and return the vertical mass flux
      call getVerticalMassFlux(mfxxyz, mfyxyz, mfzxyz, dt)
      call FILLOUT3r8 (export, 'MFZ' , mfzxyz , STATE, _RC)

      call FILLOUT3 (export, 'U'      , ur     , STATE, _RC)
      call FILLOUT3 (export, 'V'      , vr    , STATE, _RC)
      call FILLOUT3 (export, 'T'      , tempxy, STATE, _RC)
      call FILLOUT3 (export, 'Q'      , qv   , STATE, _RC)
      call FILLOUT3 (export, 'PL'     , pl  , STATE, _RC)
      call FILLOUT3 (export, 'PLE'    , vars%pe , STATE, _RC)
      call FILLOUT3 (export, 'PLK'    , plk    , STATE, _RC)
      call FILLOUT3 (export, 'PKE'    , pkxy  , STATE, _RC)
      call FILLOUT3 (export, 'PT'     , vars%pt , STATE, _RC)
      call FILLOUT3 (export, 'PE'     , vars%pe, STATE, _RC)


#ifdef SKIP_TRACERS
      do ntracer=1,ntracers
         write(myTracer, "('Q',i5.5)") ntracer-1
         call MAPL_GetPointer(export, temp3D, TRIM(myTracer), rc=status)
         VERIFY_(STATUS)
         if((associated(temp3d)) .and. (NQ>=ntracer)) then
            if (state%vars%tracer(ntracer)%is_r4) then
               dummy3d = state%vars%tracer(ntracer)%content_r4
            else
               dummy3d = state%vars%tracer(ntracer)%content
            endif
            call SSI_CopyCoarseToFine(export, dummy3d, TRIM(myTracer), STATE%f2c_SSI_arr_map, _RC)
         endif
      enddo
#endif

      call MAPL_GetPointer(export, temp3D, 'PV', rc=status)
      VERIFY_(STATUS)
      if(associated(temp3d)) then
         dummy3d = epvxyz/vars%pt
         call SSI_CopyCoarseToFine(export, dummy3d, 'PV', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      end if

      call MAPL_GetPointer(export, temp3D, 'S', rc=status)
      VERIFY_(STATUS)
      if(associated(temp3d)) then
         dummy3d = tempxy*cp
         call SSI_CopyCoarseToFine(export, dummy3d, 'S', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      end if

      call MAPL_GetPointer(export, temp3d, 'TH',rc=status)
      VERIFY_(STATUS)
  !   if(associated(temp3d)) temp3d = vars%pt*(p00**kappa)
      if(associated(temp3d)) then
         dummy3d = (tempxy)*(p00/(0.5*(vars%pe(:,:,1:km)+vars%pe(:,:,2:km+1))))**kappa
         call SSI_CopyCoarseToFine(export, dummy3d, 'TH', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      end if

      call MAPL_GetPointer(export, temp2d, 'DMDTDYN',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         dummy2d = dmdt
         call SSI_CopyCoarseToFine(export, dummy2d, 'DMDTDYN', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      end if

! Compute 3-D Tracer Dynamics Tendencies
! --------------------------------------


      call MAPL_GetPointer(export,temp3d,'QC'      , rc=status )
      VERIFY_(STATUS)

      if( associated(temp3d) ) then
          dummy3d = 0.0
          do k = 1,size(names)
             if( trim(names(k)).eq.'QLCN' .or. &
                 trim(names(k)).eq.'QILS' .or. &
                 trim(names(k)).eq.'QICN' .or. &
                 trim(names(k)).eq.'QLLS' ) then
                 if( state%vars%tracer(k)%is_r4 ) then
                     if (size(dqldt)==size(state%vars%tracer(k)%content_r4)) &
                              dummy3d = dummy3d + state%vars%tracer(k)%content_r4
                 else
                     if (size(dqldt)==size(state%vars%tracer(k)%content)) &
                              dummy3d = dummy3d + state%vars%tracer(k)%content
                 endif
             endif
          enddo
          call SSI_CopyCoarseToFine(export, dummy3d, 'QC', STATE%f2c_SSI_arr_map, rc=status)
          VERIFY_(STATUS)
      endif


      call MAPL_GetPointer( export,temp3d,'DQLDTDYN', rc=status )
      VERIFY_(STATUS)
      if( associated(temp3d) ) then
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
          dummy3d = dqldt
          call SSI_CopyCoarseToFine(export, dummy3d, 'DQLDTDYN', STATE%f2c_SSI_arr_map, rc=status)
          VERIFY_(STATUS)
      endif

      call MAPL_GetPointer( export,temp3d,'DQIDTDYN', rc=status )
      VERIFY_(STATUS)
      if( associated(temp3d) ) then
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
          dummy3d = dqidt
          call SSI_CopyCoarseToFine(export, dummy3d, 'DQIDTDYN', STATE%f2c_SSI_arr_map, rc=status)
          VERIFY_(STATUS)
      endif

      call MAPL_GetPointer( export,temp3d,'DOXDTDYN', rc=status )
      VERIFY_(STATUS)
      if( associated(temp3d) ) then
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
          dummy3d = doxdt
          call SSI_CopyCoarseToFine(export, dummy3d, 'DOXDTDYN', STATE%f2c_SSI_arr_map, rc=status)
          VERIFY_(STATUS)
      endif

! Compute 2-D Vertically Integrated Tracer Dynamics Tendencies
! ------------------------------------------------------------

      call MAPL_GetPointer ( export, temp2D, 'DQVDTDYNINT', rc=status )
      VERIFY_(STATUS)
      if( associated(temp2D) ) then
          dummy2d = 0.
          do k=1,km
          dummy2d = dummy2d + qv(:,:,k)*delp(:,:,k)
          enddo
          dummy2d = dummy2d/(grav*dt)
          call SSI_CopyCoarseToFine(export, dummy2d, 'DQVDTDYNINT', STATE%f2c_SSI_arr_map, rc=status)
          VERIFY_(STATUS)
      endif

      call MAPL_GetPointer ( export, temp2D, 'DQLDTDYNINT', rc=status )
      VERIFY_(STATUS)
      if( associated(temp2D) ) then
          dummy2d = 0.
          do N = 1,size(names)
             if( trim(names(N)).eq.'QLCN' .or. &
                 trim(names(N)).eq.'QLLS' ) then
                 if( state%vars%tracer(N)%is_r4 ) then
                     do k=1,km
                     dummy2d = dummy2d + state%vars%tracer(N)%content_r4(:,:,k)*delp(:,:,k)
                     enddo
                 else
                     do k=1,km
                     dummy2d = dummy2d + state%vars%tracer(N)%content(:,:,k)*delp(:,:,k)
                     enddo
                 endif
             endif
          enddo
          dummy2d = dummy2d/(grav*dt)
          call SSI_CopyCoarseToFine(export, dummy2d, 'DQLDTDYNINT', STATE%f2c_SSI_arr_map, rc=status)
          VERIFY_(STATUS)
      endif

      call MAPL_GetPointer ( export, temp2D, 'DQIDTDYNINT', rc=status )
      VERIFY_(STATUS)
      if( associated(temp2D) ) then
          dummy2d = 0.
          do N = 1,size(names)
             if( trim(names(N)).eq.'QICN' .or. &
                 trim(names(N)).eq.'QILS' ) then
                 if( state%vars%tracer(N)%is_r4 ) then
                     do k=1,km
                     dummy2d = dummy2d + state%vars%tracer(N)%content_r4(:,:,k)*delp(:,:,k)
                     enddo
                 else
                     do k=1,km
                     dummy2d = dummy2d + state%vars%tracer(N)%content(:,:,k)*delp(:,:,k)
                     enddo
                 endif
             endif
          enddo
          dummy2d = dummy2d/(grav*dt)
          call SSI_CopyCoarseToFine(export, dummy2d, 'DQIDTDYNINT', STATE%f2c_SSI_arr_map, rc=status)
          VERIFY_(STATUS)
      endif
      call MAPL_GetPointer ( export, temp2D, 'DOXDTDYNINT', rc=status )
      VERIFY_(STATUS)
      if( associated(temp2D) ) then
          dummy2d = 0.
          do N = 1,size(names)
             pos = index(names(N),'::')
             if(pos > 0) then
               if( (names(N)(pos+2:))=='OX' ) then
                 if( state%vars%tracer(N)%is_r4 ) then
                     do k=1,km
                     dummy2d = dummy2d + state%vars%tracer(N)%content_r4(:,:,k)*delp(:,:,k)
                     enddo
                 else
                     do k=1,km
                     dummy2d = dummy2d + state%vars%tracer(N)%content(:,:,k)*delp(:,:,k)
                     enddo
                 endif
               endif
             endif
          enddo
          dummy2d = dummy2d * (MAPL_O3MW/MAPL_AIRMW) / (MAPL_GRAV*DT)
          call SSI_CopyCoarseToFine(export, dummy2d, 'DOXDTDYNINT', STATE%f2c_SSI_arr_map, rc=status)
          VERIFY_(STATUS)
      endif

! Virtual temperature
! -------------------

      tempxy =  tempxy*(1.0+eps*qv)

      call MAPL_GetPointer(export,temp3D,'TV'   ,rc=status)
      VERIFY_(STATUS)
      if(associated(temp3D)) then
          dummy3d = tempxy
          call SSI_CopyCoarseToFine(export, dummy3d, 'TV', STATE%f2c_SSI_arr_map, rc=status)
          VERIFY_(STATUS)
      endif

! Fluxes: UCPT & VCPT
!--------------------
      call MAPL_GetPointer(export,temp2d,'UCPT',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         dummy2d = 0.0
         do k=1,km
            dummy2d = dummy2d + ur(:,:,k)*tempxy(:,:,k)*delp(:,:,k)
         enddo
         dummy2d = dummy2d*(cp/grav)
          call SSI_CopyCoarseToFine(export, dummy2d, 'UCPT', STATE%f2c_SSI_arr_map, rc=status)
          VERIFY_(STATUS)
      end if

      call MAPL_GetPointer(export,temp2d,'VCPT',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         dummy2d = 0.0
         do k=1,km
            dummy2d = dummy2d + vr(:,:,k)*tempxy(:,:,k)*delp(:,:,k)
         enddo
         dummy2d = dummy2d*(cp/grav)
          call SSI_CopyCoarseToFine(export, dummy2d, 'VCPT', STATE%f2c_SSI_arr_map, rc=status)
          VERIFY_(STATUS)
      end if

! Compute Energetics After Dycore
! -------------------------------

      tempxy = vars%pt*(1.0+eps*qv)  ! Convert TH to THV

      call MAPL_GetPointer(export,temp3d,'THV',rc=status)
      VERIFY_(STATUS)
      if(associated(temp3d)) then
          dummy3d = tempxy
          call SSI_CopyCoarseToFine(export, dummy3d, 'THV', STATE%f2c_SSI_arr_map, rc=status)
          VERIFY_(STATUS)
      endif

      if (doEnergetics) then
      call Energetics (ur,vr,tempxy,vars%pe,delp,vars%pkz,phisxy,kenrg,penrg,tenrg)
      kedyn   = (kenrg -kenrg0)/DT
      pedyn   = (penrg -penrg0)/DT
      tedyn   = (tenrg -tenrg0)/DT

      call MAPL_GetPointer(export,temp2d,'KEDYN',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
          dummy2d = kedyn
          call SSI_CopyCoarseToFine(export, dummy2d, 'KEDYN', STATE%f2c_SSI_arr_map, rc=status)
          VERIFY_(STATUS)
      endif

      call MAPL_GetPointer(export,temp2d,'PEDYN',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
          dummy2d = pedyn
          call SSI_CopyCoarseToFine(export, dummy2d, 'PEDYN', STATE%f2c_SSI_arr_map, rc=status)
          VERIFY_(STATUS)
      endif

      call MAPL_GetPointer(export,temp2d,'TEDYN',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
          dummy2d = tedyn
          call SSI_CopyCoarseToFine(export, dummy2d, 'TEDYN', STATE%f2c_SSI_arr_map, rc=status)
          VERIFY_(STATUS)
      endif
      endif

! Compute/Get Omega
! --------------------------
      zle(:,:,km+1) = phisxy(:,:)
      do k=km,1,-1
        zle(:,:,k) = zle(:,:,k+1) + cp*tempxy(:,:,k)*( pkxy(:,:,k+1)-pkxy(:,:,k) )
      enddo
      zle = zle/grav
      call getOmega ( omaxyz )

! Fluxes: UKE & VKE
! -----------------
      call MAPL_GetPointer(export,temp2d,'UKE',rc=status); VERIFY_(STATUS)
      if(associated(temp2d) ) uke_associated = .true.
      call MAPL_GetPointer(export,temp2d,'VKE',rc=status); VERIFY_(STATUS)
      if(associated(temp2d) ) vke_associated = .true.

      if(uke_associated .or. vke_associated) then
         ke = 0.5*(ur**2 + vr**2)
      end if

      if(uke_associated) then
         dummy2d = 0.0
         do k=1,km
            dummy2d = dummy2d + ur(:,:,k)*ke(:,:,k)*delp(:,:,k)
         enddo
         dummy2d = dummy2d / grav
         call SSI_CopyCoarseToFine(export, dummy2d, 'UKE', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      end if

      if(vke_associated) then
         dummy2d = 0.0
         do k=1,km
            dummy2d = dummy2d + vr(:,:,k)*ke(:,:,k)*delp(:,:,k)
         enddo
         dummy2d = dummy2d / grav
         call SSI_CopyCoarseToFine(export, dummy2d, 'VKE', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      end if

! Fluxes: UQV & VQV
! -----------------
      call MAPL_GetPointer(export,temp2d,'UQV',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         dummy2d = 0.0
         do k=1,km
            dummy2d = dummy2d + ur(:,:,k)*QV(:,:,k)*delp(:,:,k)
         enddo
         dummy2d = dummy2d / grav
         call SSI_CopyCoarseToFine(export, dummy2d, 'UQV', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      end if

      call MAPL_GetPointer(export,temp2d,'VQV',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         dummy2d = 0.0
         do k=1,km
            dummy2d = dummy2d + vr(:,:,k)*QV(:,:,k)*delp(:,:,k)
         enddo
         dummy2d = dummy2d / grav
         call SSI_CopyCoarseToFine(export, dummy2d, 'VQV', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      end if

! Fluxes: UQL & VQL
! -----------------
      call MAPL_GetPointer(export,temp2d,'UQL',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         dummy2d = 0.0
          do N = 1,size(names)
             if( trim(names(n)).eq.'QLCN' .or. &
                 trim(names(n)).eq.'QLLS' ) then
                 do k=1,km
                 if( state%vars%tracer(n)%is_r4 ) then
                      dummy2d = dummy2d + ur(:,:,k)*state%vars%tracer(n)%content_r4(:,:,k)*delp(:,:,k)
                 else
                      dummy2d = dummy2d + ur(:,:,k)*state%vars%tracer(n)%content   (:,:,k)*delp(:,:,k)
                 endif
                 enddo
             endif
          enddo
         dummy2d = dummy2d / grav
         call SSI_CopyCoarseToFine(export, dummy2d, 'UQL', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      end if

      call MAPL_GetPointer(export,temp2d,'VQL',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         dummy2d = 0.0
          do N = 1,size(names)
             if( trim(names(n)).eq.'QLCN' .or. &
                 trim(names(n)).eq.'QLLS' ) then
                 do k=1,km
                 if( state%vars%tracer(n)%is_r4 ) then
                      dummy2d = dummy2d + vr(:,:,k)*state%vars%tracer(n)%content_r4(:,:,k)*delp(:,:,k)
                 else
                      dummy2d = dummy2d + vr(:,:,k)*state%vars%tracer(n)%content   (:,:,k)*delp(:,:,k)
                 endif
                 enddo
             endif
          enddo
         dummy2d = dummy2d / grav
         call SSI_CopyCoarseToFine(export, dummy2d, 'VQL', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      end if

! Fluxes: UQI & VQI
! -----------------
      call MAPL_GetPointer(export,temp2d,'UQI',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         dummy2d = 0.0
          do N = 1,size(names)
             if( trim(names(n)).eq.'QICN' .or. &
                 trim(names(n)).eq.'QILS' ) then
                 do k=1,km
                 if( state%vars%tracer(n)%is_r4 ) then
                      dummy2d = dummy2d + ur(:,:,k)*state%vars%tracer(n)%content_r4(:,:,k)*delp(:,:,k)
                 else
                      dummy2d = dummy2d + ur(:,:,k)*state%vars%tracer(n)%content   (:,:,k)*delp(:,:,k)
                 endif
                 enddo
             endif
          enddo
         dummy2d = dummy2d / grav
         call SSI_CopyCoarseToFine(export, dummy2d, 'UQI', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      end if

      call MAPL_GetPointer(export,temp2d,'VQI',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         dummy2d = 0.0
          do N = 1,size(names)
             if( trim(names(n)).eq.'QICN' .or. &
                 trim(names(n)).eq.'QILS' ) then
                 do k=1,km
                 if( state%vars%tracer(n)%is_r4 ) then
                      dummy2d = dummy2d + vr(:,:,k)*state%vars%tracer(n)%content_r4(:,:,k)*delp(:,:,k)
                 else
                      dummy2d = dummy2d + vr(:,:,k)*state%vars%tracer(n)%content   (:,:,k)*delp(:,:,k)
                 endif
                 enddo
             endif
          enddo
         dummy2d = dummy2d / grav
         call SSI_CopyCoarseToFine(export, dummy2d, 'VQI', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      end if

! Height related diagnostics
! --------------------------
      zle(:,:,km+1) = phisxy(:,:)
      do k=km,1,-1
        zle(:,:,k) = zle(:,:,k+1) + cp*tempxy(:,:,k)*( pkxy(:,:,k+1)-pkxy(:,:,k) )
      enddo
      zle = zle/grav

      call MAPL_GetPointer(export,temp3d,'ZLE',rc=status)
      VERIFY_(STATUS)
      if(associated(temp3d)) then
         dummy3d_kmplus1 = zle
         call SSI_CopyCoarseToFine(export, dummy3d_kmplus1, 'ZLE', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      end if

      call MAPL_GetPointer(export,temp3d,'ZL' ,rc=status)
      VERIFY_(STATUS)
      if(associated(temp3d)) then
         dummy3d = 0.5*( zle(:,:,:km)+zle(:,:,2:) )
         call SSI_CopyCoarseToFine(export, dummy3d, 'ZL', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      end if

      call MAPL_GetPointer(export,temp3d,'S'  ,rc=status)
      VERIFY_(STATUS)
      if(associated(temp3d)) then
         call SSI_CopyFineToCoarse(export, dummy3d, 'S', STATE%f2c_SSI_arr_map, rc=status)
         dummy3d = dummy3d + grav*(0.5*( zle(:,:,:km)+zle(:,:,2:) ))
         call SSI_CopyCoarseToFine(export, dummy3d, 'S', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      end if

! Fluxes: UPHI & VPHI
! -------------------
      call MAPL_GetPointer(export,temp2d,'UPHI',rc=status); VERIFY_(STATUS)
      if(associated(temp2d) ) uphi_associated = .true.
      call MAPL_GetPointer(export,temp2d,'VPHI',rc=status); VERIFY_(STATUS)
      if(associated(temp2d) ) vphi_associated = .true.

      if( uphi_associated .or. vphi_associated) zl = 0.5*( zle(:,:,:km)+zle(:,:,2:) )

      if(uphi_associated) then
         dummy2d = 0.0
         do k=1,km
            dummy2d = dummy2d + ur(:,:,k)*zl(:,:,k)*delp(:,:,k)
         enddo
         call SSI_CopyCoarseToFine(export, dummy2d, 'UPHI', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      end if

      if(vphi_associated) then
         dummy2d = 0.0
         do k=1,km
            dummy2d = dummy2d + vr(:,:,k)*zl(:,:,k)*delp(:,:,k)
         enddo
         call SSI_CopyCoarseToFine(export, dummy2d, 'VPHI', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      end if

      call MAPL_GetResource ( MAPL, HGT_SURFACE, Label="HGT_SURFACE:", DEFAULT= 50.0, RC=STATUS)
      VERIFY_(STATUS)

! Fill Surface and Near-Surface Variables
! ----------------------------------------------
   if ( (KM .ne. 72) .and. (HGT_SURFACE .gt. 0.0) ) then
     ! Near surface height for surface
     ! -------------------------------
      call MAPL_GetPointer(export,temp2d,'DZ', rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         dummy2d = HGT_SURFACE
         call SSI_CopyCoarseToFine(export, dummy2d, 'DZ', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      endif

    ! Get the height above the surface
      do k=1,km+1
         zle(:,:,k) = zle(:,:,k) - zle(:,:,km+1)
      enddo

      call MAPL_GetPointer(export,temp2d,'PS',  rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         dummy2d =  vars%pe(:,:,km+1)
         call SSI_CopyCoarseToFine(export, dummy2d, 'PS', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      endif

      call MAPL_GetPointer(export,temp2d,'US',  rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         call VertInterp(dummy2d,ur,-zle,-HGT_SURFACE, status)
         VERIFY_(STATUS)
         call SSI_CopyCoarseToFine(export, dummy2d, 'US', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      end if

      call MAPL_GetPointer(export,temp2d,'VS'   ,rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         call VertInterp(dummy2d,vr,-zle,-HGT_SURFACE, status)
         VERIFY_(STATUS)
         call SSI_CopyCoarseToFine(export, dummy2d, 'VS', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      end if

      call MAPL_GetPointer(export,temp2d,'TA'   ,rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         tempxy  = vars%pt * vars%pkz
         call VertInterp(dummy2d,tempxy,-zle,-HGT_SURFACE, status)
         VERIFY_(STATUS)
         call SSI_CopyCoarseToFine(export, dummy2d, 'TA', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      end if

      call MAPL_GetPointer(export,temp2d,'QA'   ,rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         call VertInterp(dummy2d,qv,-zle,-HGT_SURFACE, status)
         VERIFY_(STATUS)
         call SSI_CopyCoarseToFine(export, dummy2d, 'QA', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      end if

      call MAPL_GetPointer(export,temp2d,'SPEED',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         call VertInterp(dummy2d,sqrt(ur**2 + vr**2),-zle,-HGT_SURFACE, status)
         VERIFY_(STATUS)
         call SSI_CopyCoarseToFine(export, dummy2d, 'SPEED', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      end if
    else
! Fill Surface with Lowest Model Level Variables
! ----------------------------------------------
      call MAPL_GetPointer(export,temp2d,'DZ', rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         dummy2d = 0.5*( zle(:,:,km)-zle(:,:,km+1) )
         call SSI_CopyCoarseToFine(export, dummy2d, 'DZ', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      end if

      call MAPL_GetPointer(export,temp2d,'PS',  rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
          dummy2d =  vars%pe(:,:,km+1)
          call SSI_CopyCoarseToFine(export, dummy2d, 'PS', STATE%f2c_SSI_arr_map, rc=status)
          VERIFY_(STATUS)
      endif

      call MAPL_GetPointer(export,temp2d,'US',  rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
          dummy2d =       ur(:,:,km)
          call SSI_CopyCoarseToFine(export, dummy2d, 'US', STATE%f2c_SSI_arr_map, rc=status)
          VERIFY_(STATUS)
      endif

      call MAPL_GetPointer(export,temp2d,'VS'   ,rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
          dummy2d =       vr(:,:,km)
          call SSI_CopyCoarseToFine(export, dummy2d, 'VS', STATE%f2c_SSI_arr_map, rc=status)
          VERIFY_(STATUS)
      endif

      call MAPL_GetPointer(export,temp2d,'TA'   ,rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
          tempxy  = vars%pt * vars%pkz
          dummy2d =   tempxy(:,:,km)
          call SSI_CopyCoarseToFine(export, dummy2d, 'TA', STATE%f2c_SSI_arr_map, rc=status)
          VERIFY_(STATUS)
      endif

      call MAPL_GetPointer(export,temp2d,'QA'   ,rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
          dummy2d =       qv(:,:,km)
          call SSI_CopyCoarseToFine(export, dummy2d, 'QA', STATE%f2c_SSI_arr_map, rc=status)
          VERIFY_(STATUS)
      endif

      call MAPL_GetPointer(export,temp2d,'SPEED',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
          dummy2d = sqrt( ur(:,:,km)**2 + vr(:,:,km)**2 )
          call SSI_CopyCoarseToFine(export, dummy2d, 'SPEED', STATE%f2c_SSI_arr_map, rc=status)
          VERIFY_(STATUS)
      endif
   endif

   call MAPL_GetPointer(export,temp2d,'WSPD_10M',rc=status)        
   VERIFY_(STATUS)
   if(associated(temp2d)) then
       call VertInterp(temp2d,sqrt(ur**2 + vr**2),-zle,-10.0, status)
       VERIFY_(STATUS)
   end if

   if (.not. HYDROSTATIC) then
   call MAPL_GetPointer(export,temp2d,'VVEL_UP_100_1000',rc=status)
   VERIFY_(STATUS)
   if(associated(temp2d)) then
       dummy2d = vars%w(ifirstxy:ilastxy,jfirstxy:jlastxy,km)
       do k=km-1,1,-1
          do j=jfirstxy,jlastxy
             do i=ifirstxy,ilastxy
                if ( (vars%w(i,j,k) > dummy2d(i-ifirstxy+1,j-jfirstxy+1)) .and. &
                     (vars%pe(i,j,k) >= 10000.0) ) then
                   dummy2d(i-ifirstxy+1,j-jfirstxy+1) = vars%w(i,j,k)
                endif
             enddo
          enddo
       enddo
       call SSI_CopyCoarseToFine(export, dummy2d, 'VVEL_UP_100_1000', STATE%f2c_SSI_arr_map, _RC)
   end if
   call MAPL_GetPointer(export,temp2d,'VVEL_DN_100_1000',rc=status)
   VERIFY_(STATUS)
   if(associated(temp2d)) then
       dummy2d = vars%w(ifirstxy:ilastxy,jfirstxy:jlastxy,km)
       do k=km-1,1,-1
          do j=jfirstxy,jlastxy
             do i=ifirstxy,ilastxy
                if ( (vars%w(i,j,k) < dummy2d(i-ifirstxy+1,j-jfirstxy+1)) .and. &
                     (vars%pe(i,j,k) >= 10000.0) ) then
                   dummy2d(i-ifirstxy+1,j-jfirstxy+1) = vars%w(i,j,k)
                endif
             enddo
          enddo
       enddo 
       call SSI_CopyCoarseToFine(export, dummy2d, 'VVEL_DN_100_1000', STATE%f2c_SSI_arr_map, _RC)
   end if
   end if

! Updraft Helicty Exports

      call MAPL_GetPointer(export,  uh25, 'UH25', ALLOC=.TRUE., rc=status); VERIFY_(STATUS)
      call MAPL_GetPointer(export,  uh03, 'UH03', ALLOC=.TRUE., rc=status); VERIFY_(STATUS)
      call MAPL_GetPointer(export, srh01,'SRH01', ALLOC=.TRUE., rc=status); VERIFY_(STATUS)
      call MAPL_GetPointer(export, srh03,'SRH03', ALLOC=.TRUE., rc=status); VERIFY_(STATUS)
      call MAPL_GetPointer(export, srh25,'SRH25', ALLOC=.TRUE., rc=status); VERIFY_(STATUS)
      ! Per WMP, this calculation is not useful if running hydrostatic
      if (.not. HYDROSTATIC) then
         if( associated( uh25) .or. associated( uh03) .or. &
            associated(srh01) .or. associated(srh03) .or. associated(srh25) ) then
            call fv_getUpdraftHelicity(uh25, uh03, srh01, srh03, srh25)
            dummy2d = uh25
            call SSI_CopyCoarseToFine(export, dummy2d, 'UH25', STATE%f2c_SSI_arr_map, _RC)
            dummy2d = uh03
            call SSI_CopyCoarseToFine(export, dummy2d, 'UH03', STATE%f2c_SSI_arr_map, _RC)
            dummy2d = srh01
            call SSI_CopyCoarseToFine(export, dummy2d, 'SRH01', STATE%f2c_SSI_arr_map, _RC)
            dummy2d = srh03
            call SSI_CopyCoarseToFine(export, dummy2d, 'SRH03', STATE%f2c_SSI_arr_map, _RC)
            dummy2d = srh25
            call SSI_CopyCoarseToFine(export, dummy2d, 'SRH25', STATE%f2c_SSI_arr_map, _RC)
         endif
      endif

! Divergence Exports

      zle = log(vars%pe)

      call MAPL_GetPointer(export,temp3d,'DIVG',  rc=status)
      VERIFY_(STATUS)
      if(associated(temp3d)) then
         dummy3d = divg
         call SSI_CopyCoarseToFine(export, dummy3d, 'DIVG', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      endif

      call MAPL_GetPointer(export,temp2d,'DIVG200',  rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         call VertInterp(dummy2d,dble(divg),zle,log(20000.)  ,  status)
         VERIFY_(STATUS)
         call SSI_CopyCoarseToFine(export, dummy2d, 'DIVG200', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      end if

      call MAPL_GetPointer(export,temp2d,'DIVG500',  rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         call VertInterp(dummy2d,dble(divg),zle,log(50000.)  ,  status)
         VERIFY_(STATUS)
         call SSI_CopyCoarseToFine(export, dummy2d, 'DIVG500', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      end if

      call MAPL_GetPointer(export,temp2d,'DIVG700',  rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         call VertInterp(dummy2d,dble(divg),zle,log(70000.)  ,  status)
         VERIFY_(STATUS)
         call SSI_CopyCoarseToFine(export, dummy2d, 'DIVG700', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      end if

      call MAPL_GetPointer(export,temp2d,'DIVG850',  rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         call VertInterp(dummy2d,dble(divg),zle,log(85000.)  ,  status)
         VERIFY_(STATUS)
         call SSI_CopyCoarseToFine(export, dummy2d, 'DIVG850', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
       end if

! Vorticity Exports

      call MAPL_GetPointer(export,temp3d,'VORT',  rc=status)
      VERIFY_(STATUS)
      if(associated(temp3d)) then
         dummy3d = vort
         call SSI_CopyCoarseToFine(export, dummy3d, 'VORT', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
       end if

      call MAPL_GetPointer(export,temp2d,'VORT200',  rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         call VertInterp(dummy2d,dble(vort),zle,log(20000.)  ,  status)
         VERIFY_(STATUS)
         call SSI_CopyCoarseToFine(export, dummy2d, 'VORT200', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      end if

      call MAPL_GetPointer(export,temp2d,'VORT500',  rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         call VertInterp(dummy2d,dble(vort),zle,log(50000.)  ,  status)
         VERIFY_(STATUS)
         call SSI_CopyCoarseToFine(export, dummy2d, 'VORT500', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      end if

      call MAPL_GetPointer(export,temp2d,'VORT700',  rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         call VertInterp(dummy2d,dble(vort),zle,log(70000.)  ,  status)
         VERIFY_(STATUS)
         call SSI_CopyCoarseToFine(export, dummy2d, 'VORT700', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      end if

      call MAPL_GetPointer(export,temp2d,'VORT850',  rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         call VertInterp(dummy2d,dble(vort),zle,log(85000.)  ,  status)
         VERIFY_(STATUS)
         call SSI_CopyCoarseToFine(export, dummy2d, 'VORT850', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
       end if

! Vertical Velocity Exports

      call FILLOUT3 (export, 'OMEGA'  , omaxyz     , STATE, rc=status)
      VERIFY_(STATUS)

      call MAPL_GetPointer(export,temp2d,'OMEGA850', rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         call VertInterp(dummy2d,omaxyz,zle,log(85000.)  , status)
         VERIFY_(STATUS)
         call SSI_CopyCoarseToFine(export, dummy2d, 'OMEGA850', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      end if

      call MAPL_GetPointer(export,temp2d,'OMEGA700', rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         call VertInterp(dummy2d,omaxyz,zle,log(70000.)  , status)
         VERIFY_(STATUS)
         call SSI_CopyCoarseToFine(export, dummy2d, 'OMEGA700', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      end if

      call MAPL_GetPointer(export,temp2d,'OMEGA500', rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         call VertInterp(dummy2d,omaxyz,zle,log(50000.)  , status)
         VERIFY_(STATUS)
         call SSI_CopyCoarseToFine(export, dummy2d, 'OMEGA500', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      end if

      call MAPL_GetPointer(export,temp2d,'OMEGA200', rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         call VertInterp(dummy2d,omaxyz,zle,log(20000.)  , status)
         VERIFY_(STATUS)
         call SSI_CopyCoarseToFine(export, dummy2d, 'OMEGA200', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      end if

      call MAPL_GetPointer(export,temp2d,'OMEGA10', rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         call VertInterp(dummy2d,omaxyz,zle,log(1000.)  , status)
         VERIFY_(STATUS)
         call SSI_CopyCoarseToFine(export, dummy2d, 'OMEGA10', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      end if

      if (.not. HYDROSTATIC) then
      call FILLOUT3 (export, 'W'  , vars%w(ifirstxy:ilastxy,jfirstxy:jlastxy,:)     , STATE, rc=status)
      VERIFY_(STATUS)

      call MAPL_GetPointer(export,temp2d,'W850', rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         call VertInterp(temp2d,vars%w(ifirstxy:ilastxy,jfirstxy:jlastxy,:),zle,log(85000.)  , status)
         VERIFY_(STATUS)
      end if

      call MAPL_GetPointer(export,temp2d,'W850', rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         call VertInterp(dummy2d,vars%w(ifirstxy:ilastxy,jfirstxy:jlastxy,:),zle,log(85000.)  , status)
         VERIFY_(STATUS)
         call SSI_CopyCoarseToFine(export, dummy2d, 'W850', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      end if

      call MAPL_GetPointer(export,temp2d,'W500', rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         call VertInterp(dummy2d,vars%w(ifirstxy:ilastxy,jfirstxy:jlastxy,:),zle,log(50000.)  , status)
         VERIFY_(STATUS)
         call SSI_CopyCoarseToFine(export, dummy2d, 'W500', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      end if

      call MAPL_GetPointer(export,temp2d,'W200', rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         call VertInterp(dummy2d,vars%w(ifirstxy:ilastxy,jfirstxy:jlastxy,:),zle,log(20000.)  , status)
         VERIFY_(STATUS)
         call SSI_CopyCoarseToFine(export, dummy2d, 'W200', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      end if

      call MAPL_GetPointer(export,temp2d,'W10', rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         call VertInterp(dummy2d,vars%w(ifirstxy:ilastxy,jfirstxy:jlastxy,:),zle,log(1000.)  , status)
         VERIFY_(STATUS)
         call SSI_CopyCoarseToFine(export, dummy2d, 'W10', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      end if
      endif

     end if   ! SW_DYNAMICS

     call PUSH_Q(STATE, import, _RC)
      !call MAPL_TimerOff(MAPL,"-DYN_EPILOGUE")

! De-Allocate Arrays
! ------------------

      DEALLOCATE( dudtana      )
      DEALLOCATE( dvdtana      )
      DEALLOCATE( dtdtana      )
      DEALLOCATE( ddpdtana     )
      DEALLOCATE( dqldt        )
      DEALLOCATE( dqidt        )
      DEALLOCATE( doxdt        )
      DEALLOCATE( dummy        )

      if (doEnergetics) then
      DEALLOCATE( KEDYN  )
      DEALLOCATE( PEDYN  )
      DEALLOCATE( TEDYN  )
      DEALLOCATE( KENRG  )
      DEALLOCATE( PENRG  )
      DEALLOCATE( TENRG  )
      DEALLOCATE( KENRG0 )
      DEALLOCATE( PENRG0 )
      DEALLOCATE( TENRG0 )
      endif

      DEALLOCATE( qsum1 )
      DEALLOCATE( qsum2 )

      DEALLOCATE( zl     )
      DEALLOCATE( zle    )
      DEALLOCATE( plk    )
      DEALLOCATE( pkxy   )
      DEALLOCATE( vort   )
      DEALLOCATE( divg   )
      DEALLOCATE( tmp3d  )
      DEALLOCATE( omaxyz )
      DEALLOCATE( epvxyz )
      DEALLOCATE(  cxxyz )
      DEALLOCATE(  cyxyz )
      DEALLOCATE( mfxxyz )
      DEALLOCATE( mfyxyz )
      DEALLOCATE( mfzxyz )
      DEALLOCATE( tempxy )
      DEALLOCATE( pe0    )
      DEALLOCATE( pe1    )
      DEALLOCATE( pl     )
      DEALLOCATE( ua     )
      DEALLOCATE( va     )
      DEALLOCATE( uc     )
      DEALLOCATE( vc     )
      DEALLOCATE( uc0    )
      DEALLOCATE( vc0    )
      DEALLOCATE( ur     )
      DEALLOCATE( vr     )
      DEALLOCATE( qv     )
      DEALLOCATE( ql     )
      DEALLOCATE( qi     )
      DEALLOCATE( qr     )
      DEALLOCATE( qs     )
      DEALLOCATE( qg     )
      DEALLOCATE( ox     )
      DEALLOCATE( delp   )
      DEALLOCATE( dmdt   )
      DEALLOCATE( dudt   )
      DEALLOCATE( dvdt   )
      DEALLOCATE( dtdt   )
      DEALLOCATE( dqdt   )
      DEALLOCATE( dthdt  )
      DEALLOCATE( dpedt  )
      DEALLOCATE( ddpdt  )
      DEALLOCATE( phisxy )
      if (allocated(names)) DEALLOCATE( names  )
      if (allocated(names0)) DEALLOCATE( names0  )

      !call freeTracers(state)

  !call MAPL_TimerOff(MAPL,"RUN")
  !call MAPL_TimerOff(MAPL,"TOTAL")

 !if (ADIABATIC) then
 !  ! Fill Exports
 !   call RunAddIncs(gc, import, export, clock, rc)
 !endif

  RETURN_(ESMF_SUCCESS)

contains

subroutine check_replay_time_(lring)

   logical :: lring
   integer :: REPLAY_REF_DATE, REPLAY_REF_TIME, REPLAY_REF_TGAP
   integer :: REF_TIME(6), REF_TGAP(6)
   type (ESMF_TimeInterval)  :: RefTGap

   call MAPL_GetResource(MAPL, ReplayType, 'REPLAY_TYPE:', default="FULL", rc=status )
!  if (trim(ReplayType) == "FULL") return

   CALL MAPL_GetResource( MAPL, REPLAY_REF_DATE, label = 'REPLAY_REF_DATE:', Default=-1, rc=status )
   CALL MAPL_GetResource( MAPL, REPLAY_REF_TIME, label = 'REPLAY_REF_TIME:', Default=-1, rc=status )
   CALL MAPL_GetResource( MAPL, REPLAY_REF_TGAP, label = 'REPLAY_REF_TGAP:', Default=-1, rc=status )

   if(REPLAY_REF_DATE==-1.or.REPLAY_REF_TIME==-1) return

   REF_TIME(1) =     REPLAY_REF_DATE/10000
   REF_TIME(2) = mod(REPLAY_REF_DATE,10000)/100
   REF_TIME(3) = mod(REPLAY_REF_DATE,100)
   REF_TIME(4) =     REPLAY_REF_TIME/10000
   REF_TIME(5) = mod(REPLAY_REF_TIME,10000)/100
   REF_TIME(6) = mod(REPLAY_REF_TIME,100)

! set replay time
! ---------------
   call ESMF_TimeSet(  RefTime, YY =  REF_TIME(1), &
                                MM =  REF_TIME(2), &
                                DD =  REF_TIME(3), &
                                H  =  REF_TIME(4), &
                                M  =  REF_TIME(5), &
                                S  =  REF_TIME(6), rc=status ); VERIFY_(STATUS)
  if (REPLAY_REF_TGAP>0) then
      REF_TGAP    = 0
      REF_TGAP(4) =     REPLAY_REF_TGAP/10000
      REF_TGAP(5) = mod(REPLAY_REF_TGAP,10000)/100
      REF_TGAP(6) = mod(REPLAY_REF_TGAP,100)
      call ESMF_TimeIntervalSet(  RefTGap, YY = REF_TGAP(1), &
                                           MM = REF_TGAP(2), &
                                            D = REF_TGAP(3), &
                                            H = REF_TGAP(4), &
                                            M = REF_TGAP(5), &
                                            S = REF_TGAP(6), &
                                    startTime = currentTime, &
                                                rc = STATUS  ); VERIFY_(STATUS)

      RefTime = RefTime - RefTGap
  endif

! check if it's time to replay
! ----------------------------
  if(RefTime==currentTime) then
     lring=.true.
  else
     lring=.false.
  endif

! In this case, increment RefTime to proper time
! ----------------------------------------------
  if (REPLAY_REF_TGAP>0) then
      RefTime = currentTime + RefTGap
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
!
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
      call ESMFL_BundleGetPointertoData(ANA_Bundle,trim(uname),XTMP3d, RC=STATUS)
      VERIFY_(STATUS)
      iwind=iwind+1
    endif
! V
    if( trim(vname).ne.'NULL' ) then
       call ESMFL_BundleGetPointertoData(ANA_Bundle,trim(vname),YTMP3D, RC=STATUS)
       VERIFY_(STATUS)
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
          call l2c%regrid(XTMP3d, cubeTEMP3D, RC=STATUS )
          VERIFY_(STATUS)
          call l2c%regrid(YTMP3d, cubeVTMP3D, RC=STATUS )
          VERIFY_(STATUS)
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
          call getAllWinds (vars%u,vars%v,UR=ana_u,VR=ana_v)
          ! transform background A-grid winds to lat-lon
          call regridder_manager%make_regridder(ESMFGRID, ANAGrid, REGRID_METHOD_BILINEAR, RC=STATUS)
          VERIFY_(STATUS)
          cubeTEMP3d = ana_u(grid%is:grid%ie,grid%js:grid%je,1:km) ! copy to satisfy interface below
          cubeVTMP3d = ana_v(grid%is:grid%ie,grid%js:grid%je,1:km) ! copy to satisfy interface below
          call c2l%regrid(cubeTEMP3d, cubeVTMP3d, UAtmpR4,    VAtmpR4, RC=STATUS)
          VERIFY_(STATUS)
          ! calculate unrotated analysis increments of lat-lon U/V-A-grid winds
          UAtmpR4 = XTMP3d-UAtmpR4
          UAtmpR4 = VTMP3d-VAtmpR4
          ! convert the lat-lon A-grid wind increment back to the cubed
          call WRITE_PARALLEL('Replaying winds')
          call l2c%regrid(UAtmpR4,    VAtmpR4, cubeTEMP3d, cubeVTMP3d, RC=STATUS)
          ! convert cubed wind increment to D-grid
          allocate( UDtmp(grid%is:grid%ie  ,grid%js:grid%je+1,km) )
          allocate( VDtmp(grid%is:grid%ie+1,grid%js:grid%je  ,km) )
          deallocate(ana_u,ana_v)
          allocate( ana_u(grid%is:grid%ie  ,grid%js:grid%je  ,km) )
          allocate( ana_v(grid%is:grid%ie  ,grid%js:grid%je  ,km) )
          ana_u = cubeTEMP3d ! need this to satisfy interface below
          ana_v = cubeVTMP3d ! need this to satisfy interface below
          call Agrid_To_Native( ana_u, ana_v, UDtmp, VDtmp ) ! Calculate D-grid winds from rotated A-grid winds
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
       call ESMFL_BundleGetPointertoData(ANA_Bundle,trim(dpname),XTMP3d, RC=STATUS)
       VERIFY_(STATUS)
       call WRITE_PARALLEL('Replaying '//trim(dpname))
       if ( iapproach == 1 ) then ! convert lat-lon delp to cubed and proceed
          allocate(cubeTEMP3D(size(vars%pe,1),size(vars%pe,2),km))
          call l2c%regrid(XTMP3d, cubeTEMP3D, RC=STATUS )
          VERIFY_(STATUS)
          ana_dp=cubeTEMP3D
          deallocate(cubeTEMP3D)
       else ! just because pressure is such delicate beast: convert cubed delp
            ! to lat-lon, calculate an increment in lat-lon, convert increment
            ! on delp to cubed, and create cubed version of analyzed delp
            allocate(aux3d (size(XTMP3d,1),size(XTMP3d,2),km))
            allocate(cubeTEMP3D(size(vars%pe,1),size(vars%pe,2),km))
            ! delp on the cube
            cubeTEMP3D(:,:,:) = vars%pe(:,:,2:)-vars%pe(:,:,:km)
            ! transform cubed delp
            c2l => regridder_manager%make_regridder(ESMFGRID, ANAGrid, REGRID_METHOD_BILINEAR, RC=STATUS )
            VERIFY_(STATUS)
            call c2l%regrid(cubeTEMP3D, aux3d, RC=STATUS )
            VERIFY_(STATUS)
            ! calculate delp increment on lat-lon and transform it to cubed
            aux3d = XTMP3d - aux3d
            call l2c%regrid(aux3d, cubeTEMP3D, RC=STATUS )
            VERIFY_(STATUS)
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
          call ESMFL_BundleGetPointertoData(ANA_Bundle,trim(psname),XTMP2D, RC=STATUS)
          VERIFY_(STATUS)
          call WRITE_PARALLEL('Replaying '//trim(psname))
          allocate(cubeTEMP3D(size(vars%pe,1),size(vars%pe,2),1))
          allocate(     aux3D(size(XTMP2d ,1),size(XTMP2d ,2),1))
          if ( iapproach == 1 ) then ! convert lat-lon delp to cubed and proceed
             aux3d(:,:,1)=XTMP2D ! rank-2 interface to HorzT does not work
             call l2c%regrid(aux3d, cubeTEMP3D, RC=STATUS )
             VERIFY_(STATUS)
          else ! operate on increment to ps
             ! transform cubed delp
             cubeTEMP3D(:,:,1) = vars%pe(:,:,km+1) ! cubed ps
             c2l => regridder_manager%make_regridder(ESMFGRID, ANAGrid, REGRID_METHOD_BILINEAR, RC=STATUS )
             VERIFY_(STATUS)
             call c2l%regrid(cubeTEMP3D, aux3d, RC=STATUS )
             VERIFY_(STATUS)
             ! increment to ps on the lat-lon
             aux3d(:,:,1) = XTMP2D - aux3d(:,:,1)
             ! lat-lon increment to ps converted to the cube
             call l2c%regrid(aux3d, cubeTEMP3D, RC=STATUS )
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
       call ESMFL_BundleGetPointertoData(ANA_Bundle,trim(o3name),XTMP3d, RC=STATUS)
       VERIFY_(STATUS)
       allocate(cubeTEMP3D(size(vars%pe,1),size(vars%pe,2),km))
       call l2c%regrid(XTMP3d, cubeTEMP3D, RC=STATUS )
       VERIFY_(STATUS)

!      Ozone needs to be adjusted to OX
!      --------------------------------
       call WRITE_PARALLEL('Replaying '//trim(o3name))

       call MAPL_Get(MAPL, ORBIT=ORBIT, RC=STATUS )
       VERIFY_(STATUS)

       if(.not.associated(LATS)) ALLOCATE(LATS(grid%is:grid%ie, grid%js:grid%je),stat=status)
       VERIFY_(STATUS)
    call SSI_CopyFineToCoarse(internal, LATS, 'LATS', STATE%f2c_SSI_arr_map, rc=status)
    VERIFY_(STATUS)

       if(.not.associated(LONS)) ALLOCATE(LONS(grid%is:grid%ie, grid%js:grid%je),stat=status)
       VERIFY_(STATUS)
    call SSI_CopyFineToCoarse(internal, LONS, 'LONS', STATE%f2c_SSI_arr_map, rc=status)
    VERIFY_(STATUS)

       allocate( ZTH( size(LONS,1),size(LONS,2) ) )
       allocate( SLR( size(LONS,1),size(LONS,2) ) )

       call MAPL_SunGetInsolation( LONS,LATS,ORBIT,ZTH,SLR, CLOCK=CLOCK,RC=STATUS  )
       VERIFY_(STATUS)

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
       call ESMFL_BundleGetPointertoData(ANA_Bundle,trim(qname),XTMP3d, RC=STATUS)
       VERIFY_(STATUS)
       allocate(cubeTEMP3D(size(vars%pe,1),size(vars%pe,2),km))
       call l2c%regrid(XTMP3d, cubeTEMP3D, RC=STATUS )
       VERIFY_(STATUS)
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
       call ESMFL_BundleGetPointertoData(ANA_Bundle,trim(tname),XTMP3d, RC=STATUS)
       VERIFY_(STATUS)
       allocate(cubeTEMP3D(size(ana_thv,1),size(ana_thv,2),km))
       call l2c%regrid(XTMP3d, cubeTEMP3D, RC=STATUS )
       VERIFY_(STATUS)
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

!   Refresh vars ("update" them)
!   -------------
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
       call ESMFL_BundleGetPointertoData(ANA_Bundle,trim(uname),TEMP3D, RC=STATUS)
       VERIFY_(STATUS)
       iwind=iwind+1
    endif
! V
    if( trim(vname).ne.'NULL' ) then
       call ESMFL_BundleGetPointertoData(ANA_Bundle,trim(vname),VTMP3D, RC=STATUS)
       VERIFY_(STATUS)
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
       call l2c%regrid(TEMP3D, cubeTEMP3D, RC=STATUS )
       VERIFY_(STATUS)
       call l2c%regrid(VTMP3D, cubeVTMP3D, RC=STATUS )
       VERIFY_(STATUS)
#else
       call WRITE_PARALLEL('Replaying increment of winds')
       call l2c%regrid(TEMP3d,     VTMP3d, cubeTEMP3d, cubeVTMP3d, RC=STATUS)
#endif /* SCALAR_WINDS */
       allocate( UAtmp(grid%is:grid%ie  ,grid%js:grid%je  ,km) )
       allocate( VAtmp(grid%is:grid%ie  ,grid%js:grid%je  ,km) )
       UAtmp = cubetemp3d ! A-grid winds on cube
       VAtmp = cubevtmp3d ! A-grid winds on cube
       call Agrid_To_Native( UAtmp, VAtmp, du, dv )       ! Calculate D-grid winds from rotated A-grid winds
       deallocate(uatmp,vatmp)
       deallocate(cubeTEMP3D)
       deallocate(cubeVTMP3D)
    endif

! DELP
    if( trim(psname)=='NULL' .and. trim(dpname).ne.'NULL' ) then
       call ESMFL_BundleGetPointertoData(ANA_Bundle,trim(dpname),TEMP3D, RC=STATUS)
       VERIFY_(STATUS)
       call WRITE_PARALLEL('Replaying increment of '//trim(dpname))
       allocate(cubeTEMP3D(size(vars%pe,1),size(vars%pe,2),km))
       call l2c%regrid(TEMP3D, cubeTEMP3D, RC=STATUS )
       VERIFY_(STATUS)
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
       call ESMFL_BundleGetPointertoData(ANA_Bundle,trim(psname),TEMP2D, RC=STATUS)
       VERIFY_(STATUS)
       call WRITE_PARALLEL('Replaying increment of '//trim(psname))
       allocate(cubeTEMP3D(size(vars%pe,1),size(vars%pe,2),1))
       allocate(     aux3D(size( TEMP2D,1),size( TEMP2D,2),1))
       aux3d(:,:,1) = TEMP2D ! same trick of putting in rank-3 array for transforms
       call l2c%regrid(aux3d, cubeTEMP3D, RC=STATUS )
       VERIFY_(STATUS)
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
       call ESMFL_BundleGetPointertoData(ANA_Bundle,trim(o3name),TEMP3D, RC=STATUS)
       VERIFY_(STATUS)
       allocate(cubeTEMP3D(size(vars%pe,1),size(vars%pe,2),km))
       call l2c%regrid(TEMP3D, cubeTEMP3D, RC=STATUS )
       VERIFY_(STATUS)

!      Ozone needs to be adjusted to OX
!      --------------------------------
       call WRITE_PARALLEL('Replaying increment of '//trim(o3name))

       call MAPL_Get(MAPL, ORBIT=ORBIT, RC=STATUS )
       VERIFY_(STATUS)

       if(.not.associated(LATS)) ALLOCATE(LATS(grid%is:grid%ie, grid%js:grid%je),stat=status)
       VERIFY_(STATUS)
    call SSI_CopyFineToCoarse(internal, LATS, 'LATS', STATE%f2c_SSI_arr_map, rc=status)
    VERIFY_(STATUS)

       if(.not.associated(LONS)) ALLOCATE(LONS(grid%is:grid%ie, grid%js:grid%je),stat=status)
       VERIFY_(STATUS)
    call SSI_CopyFineToCoarse(internal, LONS, 'LONS', STATE%f2c_SSI_arr_map, rc=status)
    VERIFY_(STATUS)

       allocate( ZTH( size(LONS,1),size(LONS,2) ) )
       allocate( SLR( size(LONS,1),size(LONS,2) ) )

       call MAPL_SunGetInsolation( LONS,LATS,ORBIT,ZTH,SLR, CLOCK=CLOCK,RC=STATUS  )
       VERIFY_(STATUS)

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
       call ESMFL_BundleGetPointertoData(ANA_Bundle,trim(qname),TEMP3D, RC=STATUS)
       VERIFY_(STATUS)
       allocate(cubeTEMP3D(size(vars%pe,1),size(vars%pe,2),km))
       call l2c%regrid(TEMP3D, cubeTEMP3D, RC=STATUS )
       VERIFY_(STATUS)
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
       call ESMFL_BundleGetPointertoData(ANA_Bundle,trim(tname),TEMP3D, RC=STATUS)
       VERIFY_(STATUS)
       allocate(cubeTEMP3D(size(vars%pe,1),size(vars%pe,2),km))
       call l2c%regrid(TEMP3D, cubeTEMP3D, RC=STATUS )
       VERIFY_(STATUS)
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
! ----------------------------------------------------------
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
          call ESMF_FieldBundleGet(BUNDLE, N, Field, RC=STATUS )
          call ESMF_FieldGet(Field, dimCount = rank, RC=STATUS )
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
    call ESMFL_BundleGetPointertoData(ANA_Bundle,'phis',XTMP2D, RC=STATUS)
    VERIFY_(STATUS)
    allocate(cubeTEMP3D(size(vars%pe,1),size(vars%pe,2),1))
    allocate(     aux3D(size(XTMP2D ,1),size(XTMP2D ,2),1))
    aux3d(:,:,1)=XTMP2D ! this is a trick since the 2d interface to the transform has not worked for me (RT)
    call l2c%regrid(aux3D, cubeTEMP3D, RC=STATUS )
    VERIFY_(STATUS)
    ana_phis=cubeTEMP3D(:,:,1)
    deallocate(     aux3D)
    deallocate(cubeTEMP3D)
!
    if (remap_all_tracers) then
       icnt=0
       do N=1,NQ
          call ESMF_FieldBundleGet(BUNDLE, N, Field, RC=STATUS )
          call ESMF_FieldGet(Field, NAME=NAME, dimCount=rank, RC=STATUS )
          if (rank==2) cycle
          if (rank==3) then
             icnt=icnt+1
             if (icnt>nq3d) then
                 call WRITE_PARALLEL('state_remap: number of tracers exceeds known value')
                 status=999
                 VERIFY_(STATUS)
             endif
             call ESMFL_BundleGetPointerToData(BUNDLE, NAME, ptr3dr4, RC=STATUS )
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
          call ESMF_FieldBundleGet(BUNDLE, N, Field, RC=STATUS )
          call ESMF_FieldGet(Field, NAME=NAME, dimCount=rank, RC=STATUS )
          if (rank==2) cycle
          if (rank==3) then
             icnt=icnt+1
             call ESMFL_BundleGetPointerToData(BUNDLE, NAME, ptr3dr4, RC=STATUS )
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

end subroutine RUN

!-----------------------------------------------------------------------

  subroutine PULL_Q(STATE, IMPORT, QQQ, iNXQ, InFieldName, RC)

    type (DynState)        :: STATE
    type (ESMF_State)              :: IMPORT
    type (DynTracers)               :: QQQ       ! Specific Humidity
    integer,           intent(IN)  :: iNXQ
    character(len=*), optional, intent(IN) :: InFieldName
    integer, optional, intent(OUT) :: RC

    integer                          :: STATUS
    character(len=ESMF_MAXSTR)       :: IAm="CoarsePull_Q"
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
!---------------------

    call ESMF_FieldBundleGet(BUNDLE, fieldCount=NQ, RC=STATUS)
    VERIFY_(STATUS)

               NQ = NQ + iNXQ
    STATE%GRID%NQ = NQ       ! GRID%NQ is now the "official" NQ

!
! Tracer pointer array
!
    !IF( ASSOCIATED( STATE%VARS%tracer ) ) then
    !    call freeTracers(state)
    !ENDIF

    !ALLOCATE(STATE%VARS%tracer(nq), STAT=STATUS)
    !VERIFY_(STATUS)

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
          !call ESMF_ArrayGet(array, localDE=0, farrayptr=ptr_r4, rc=status)
          !VERIFY_(STATUS)
          !state%vars%tracer(n)%content_r4 => MAPL_RemapBounds(PTR_R4, &
          !      i1,in,j1,jn, 1, km)

          call SSI_BundleCopyFineToCoarse(bundle, &
                    state%vars%tracer(n)%content_r4, n, STATE%f2c_SSI_arr_map, rc=status)
          VERIFY_(STATUS)

          if (fieldname == QFieldName) then
             qqq%is_r4 = .true.
             qqq%content_r4 => state%vars%tracer(n)%content_r4
          end if

       else

          !call ESMF_ArrayGet(array, localDE=0, farrayptr=ptr_r8, rc=status)
          !VERIFY_(STATUS)
          !state%vars%tracer(n)%content => PTR_R8

          call SSI_BundleCopyFineToCoarse(bundle, &
                    state%vars%tracer(n)%content, n, STATE%f2c_SSI_arr_map, rc=status)
          VERIFY_(STATUS)

          if (fieldname == QFieldName) then
             qqq%is_r4   = .false.
             qqq%content => state%vars%tracer(n)%content
          end if

       endif
     END DO

     RETURN_(ESMF_SUCCESS)

  end subroutine PULL_Q

!-----------------------------------------------------------------------

  subroutine PUSH_Q(STATE, IMPORT, RC)

    type (DynState)        :: STATE
    type (ESMF_State)              :: IMPORT
    integer, optional, intent(OUT) :: RC

    integer                          :: STATUS
    character(len=ESMF_MAXSTR)       :: IAm="CoarsePush_Q"
    type (ESMF_FieldBundle)          :: BUNDLE
    type (ESMF_Field)                :: field
    integer                          :: N,NQ

    BUNDLE = bundleAdv


! Count the friendlies
!---------------------

    call ESMF_FieldBundleGet(BUNDLE, fieldCount=NQ, RC=STATUS)
    VERIFY_(STATUS)

!
    DO n = 1, NQ

       if ( STATE%VARS%TRACER(N)%IS_R4 ) then

          call SSI_BundleCopyCoarseToFine(bundle, &
                    state%vars%tracer(n)%content_r4, n, STATE%f2c_SSI_arr_map, rc=status)
          VERIFY_(STATUS)

       else

          call SSI_BundleCopyCoarseToFine(bundle, &
                    state%vars%tracer(n)%content, n, STATE%f2c_SSI_arr_map, rc=status)
          VERIFY_(STATUS)

       endif
     END DO

     RETURN_(ESMF_SUCCESS)

  end subroutine PUSH_Q

!-----------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!BOP

! !IROUTINE: RunAddIncs

! !DESCRIPTION: This is the second registered stage of FV.
!    It calls an Fv supplied routine to add external contributions
!    to FV's state variables. It does not touch the Friendly tracers.
!    It also computes additional diagnostics and updates the
!    FV internal state to reflect the added tendencies.
!
!
! !INTERFACE:

  subroutine RunAddIncs(gc, import, export, clock, rc)

! !ARGUMENTS:

    type(ESMF_GridComp) :: gc
    type (ESMF_State) :: import
    type (ESMF_State) :: export
    type (ESMF_Clock) :: clock
    integer, intent(out) :: rc

!EOP

! !Local Variables:

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
    real(r4),     pointer ::   phis(:,:) => Null()
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

    real(r8), allocatable ::  logpl(:,:,:)
    real(r8), allocatable ::  logpe(:,:,:)
    real(r8), allocatable ::  logps(:,:)

    real(FVPRC)              :: dt

    real(r4), pointer     :: QOLD(:,:,:)
    real(r4), pointer     :: temp3d(:,:,:)
    real(r4), pointer     :: dummy3d(:,:,:) => Null()
    real(r4), pointer     :: temp2d(:,:  )
    real(r4), pointer     :: dummy2d(:,:  ) => Null()

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

    character(len=ESMF_MAXSTR) :: COMP_NAME
    integer, allocatable :: gcImg(:)
    integer :: itemCount
    type(ESMF_GridComp) :: fineGC

    Iam = "RunAddIncs"
    call ESMF_GridCompGet( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // trim(Iam)

! Retrieve fine GC
! ---------------------------------
    call ESMF_AttributeGet(GC, name='GC_IMAGE', itemCount=itemCount, rc=status)
    VERIFY_(STATUS)
    allocate(gcImg(itemCount), stat=status)
    VERIFY_(STATUS)
    call ESMF_AttributeGet(GC, name='GC_IMAGE', valueList=gcImg, rc=status)
    VERIFY_(STATUS)
    fineGC = transfer(gcImg, fineGC)
    deallocate(gcImg, stat=status)
    VERIFY_(STATUS)

! Retrieve the pointer to the generic state
! -----------------------------------------

    call MAPL_GetObjectFromGC (fineGC, GENSTATE,  RC=STATUS )
    VERIFY_(STATUS)

    !call MAPL_TimerOn(GENSTATE,"TOTAL")
    !call MAPL_TimerOn(GENSTATE,"RUN2")

! Retrieve the pointer to the internal state
! ------------------------------------------

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

    if(.not.associated(dummy3d)) then
       ALLOCATE(dummy3d(ifirstxy:ilastxy,jfirstxy:jlastxy,km),stat=status)
       VERIFY_(STATUS)
    endif
    if(.not.associated(dummy2d)) then
       ALLOCATE(dummy2d(ifirstxy:ilastxy,jfirstxy:jlastxy),stat=status)
       VERIFY_(STATUS)
    endif

    !call MAPL_GetPointer ( IMPORT, PHIS, 'PHIS', RC=STATUS )
    !VERIFY_(STATUS)
    if(.not.associated(phis))  then
       allocate(phis(ifirstxy:ilastxy,jfirstxy:jlastxy), stat=status)
       VERIFY_(STATUS)
    endif
    call SSI_CopyFineToCoarse(import, phis, 'PHIS', STATE%f2c_SSI_arr_map, rc=status)
    VERIFY_(STATUS)

    phisxy = real(phis,kind=r8)

! Compute Pressure Thickness
! --------------------------

    dp = ( vars%pe(:,:,2:) - vars%pe (:,:,:km) )

! Load Specific Humidity
! ----------------------


    call MAPL_GetPointer(export,QOLD,'Q',  rc=status)

    call PULL_Q ( STATE, IMPORT, qqq, iNXQ, RC=rc )
    if ((.not. ADIABATIC) .and. (STATE%GRID%NQ > 0)) then
      if ( (qqq%is_r4) .and. (associated(qqq%content_r4)) ) then
       if (size(qv)==size(qqq%content_r4)) qv = qqq%content_r4
      elseif (associated(qqq%content)) then
       if (size(qv)==size(qqq%content)) qv = qqq%content
      endif
    else
      qv = 0.0
    endif

! Compute Energetics Before Diabatic Forcing
! ------------------------------------------
    if (associated(QOLD)) then
        call SSI_CopyFineToCoarse(export, dummy3d, 'Q', STATE%f2c_SSI_arr_map, _RC)
       thv = vars%pt*(1.0+eps*dummy3d)
    else
       thv = vars%pt
    endif

    if (doEnergetics) then
      call getAllWinds(vars%u, vars%v, UA=ua, VA=va, UC=uc, VC=vc, UR=ur, VR=vr)
      call Energetics (ur,vr,thv,vars%pe,dp,vars%pkz,phisxy,kenrg0,penrg0,tenrg0)
    endif

! DTHVDTPHYINT
! ------------
      call MAPL_GetPointer ( export, temp2D, 'DTHVDTPHYINT', rc=status )
      VERIFY_(STATUS)
      if( associated(temp2D) ) then
          dthdtphyint1 = 0.0
          do k=1,km
          dthdtphyint1 = dthdtphyint1 + thv(:,:,k)*dp(:,:,k)
          enddo
      endif

! Add Diabatic Forcing to State Variables
! ---------------------------------------

    call ADD_INCS ( STATE,IMPORT, internal, DT  )

    if (DYN_DEBUG) call DEBUG_FV_STATE('PHYSICS ADD_INCS',STATE)

! Update Mid-Layer Pressure and Pressure Thickness
! ------------------------------------------------

    dp = ( vars%pe(:,:,2:) - vars%pe (:,:,:km) )
    pl = ( vars%pe(:,:,2:) + vars%pe (:,:,:km) )*0.5

    logpl = log(pl)
    logpe = log(vars%pe)
    logps = log(vars%pe(:,:,km+1))

! Get Cubed-Sphere Wind Exports
! -----------------------------
    call getAllWinds(vars%u, vars%v, UA=ua, VA=va, UC=uc, VC=vc, UR=ur, VR=vr)
    call FILLOUT3 (export, 'U_DGRID', vars%u  , STATE, rc=status); VERIFY_(STATUS)
    call FILLOUT3 (export, 'V_DGRID', vars%v  , STATE, rc=status); VERIFY_(STATUS)
    call FILLOUT3 (export, 'U_CGRID', uc      , STATE, rc=status); VERIFY_(STATUS)
    call FILLOUT3 (export, 'V_CGRID', vc      , STATE, rc=status); VERIFY_(STATUS)
    call FILLOUT3 (export, 'U_AGRID', ua      , STATE, rc=status); VERIFY_(STATUS)
    call FILLOUT3 (export, 'V_AGRID', va      , STATE, rc=status); VERIFY_(STATUS)

! Compute Energetics After Diabatic Forcing
! -----------------------------------------

    thv = vars%pt*(1.0+eps*qv)

#if defined(DEBUG_VPT)
  call Write_Profile(grid, thv, 'VPT')
#endif

    if (doEnergetics) then
      call Energetics (ur,vr,thv,vars%pe,dp,vars%pkz,phisxy,kenrg,penrg,tenrg)
      call MAPL_GetPointer(export,temp2d,'KE',  rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then ! temp2d = kenrg
          call SSI_CopyCoarseToFine(export, kenrg, 'KE', STATE%f2c_SSI_arr_map, rc=status)
         VERIFY_(STATUS)
      endif
      kenrg = (kenrg-kenrg0)/DT
      penrg = (penrg-penrg0)/DT
      tenrg = (tenrg-tenrg0)/DT
      call FILLOUT2 (export, 'KEPHY', kenrg, STATE, rc=status); VERIFY_(STATUS)
      call FILLOUT2 (export, 'PEPHY', penrg, STATE, rc=status); VERIFY_(STATUS)
      call FILLOUT2 (export, 'TEPHY', tenrg, STATE, rc=status); VERIFY_(STATUS)
    endif

! DTHVDTPHYINT
! ------------
      call MAPL_GetPointer ( export, temp2D, 'DTHVDTPHYINT', rc=status )
      VERIFY_(STATUS)
      if( associated(temp2D) ) then
          dthdtphyint2 = 0.0
          do k=1,km
          dthdtphyint2 = dthdtphyint2 + thv(:,:,k)*dp(:,:,k)
          enddo
          dummy2d       = (dthdtphyint2-dthdtphyint1) * MAPL_P00**MAPL_KAPPA / (MAPL_GRAV*DT)
          call SSI_CopyCoarseToFine(export, dummy2d, 'DTHVDTPHYINT', STATE%f2c_SSI_arr_map, _RC)
      endif

    plk = exp( kappa * log( 0.5*(vars%pe(:,:,1:km)+vars%pe(:,:,2:km+1)) ) )
    pke = exp( kappa * log( vars%pe ) )

    tempxy = vars%pt * vars%pkz   ! Dry Temperature

#if defined(DEBUG_T)
  call Write_Profile(grid, tempxy, 'T')
#endif

    call FILLOUT3 (export, 'DELP'   , dp      , STATE, rc=status); VERIFY_(STATUS)
    call FILLOUT3 (export, 'U'      , ur      , STATE, rc=status); VERIFY_(STATUS)
    call FILLOUT3 (export, 'V'      , vr      , STATE, rc=status); VERIFY_(STATUS)
    call FILLOUT3 (export, 'T'      , tempxy  , STATE, rc=status); VERIFY_(STATUS)
    call FILLOUT3 (export, 'Q'      , qv      , STATE, rc=status); VERIFY_(STATUS)
    call FILLOUT3 (export, 'PL'     , pl      , STATE, rc=status); VERIFY_(STATUS)
    call FILLOUT3 (export, 'PLE'    , vars%pe , STATE, rc=status); VERIFY_(STATUS)
    call FILLOUT3 (export, 'PLK'    , plk     , STATE, rc=status); VERIFY_(STATUS)
    call FILLOUT3 (export, 'PKE'    , pke     , STATE, rc=status); VERIFY_(STATUS)
    call FILLOUT3 (export, 'THV'    , thv     , STATE, rc=status); VERIFY_(STATUS)
    call FILLOUT3 (export, 'PT'     , vars%pt , STATE, rc=status); VERIFY_(STATUS)
    call FILLOUT3 (export, 'PE'     , vars%pe , STATE, rc=status); VERIFY_(STATUS)

    call MAPL_GetPointer(export,temp3d,'TH',rc=status)
    VERIFY_(STATUS)
    if(associated(temp3d)) then
       dummy3d = (tempxy)*(p00/(0.5*(vars%pe(:,:,1:km)+vars%pe(:,:,2:km+1))))**kappa
       call SSI_CopyCoarseToFine(export, dummy3d, 'TH', STATE%f2c_SSI_arr_map, _RC)
    endif

#ifdef SKIP_TRACERS
      do ntracer=1,ntracers
         write(myTracer, "('Q',i5.5)") ntracer-1
         call MAPL_GetPointer(export, temp3D, TRIM(myTracer), rc=status)
         VERIFY_(STATUS)
         if((associated(temp3d)) .and. (STATE%GRID%NQ>=ntracer)) then
            if (state%vars%tracer(ntracer)%is_r4) then
               dummy3d = state%vars%tracer(ntracer)%content_r4
            else
               dummy3d = state%vars%tracer(ntracer)%content
            endif
            call SSI_CopyCoarseToFine(export, dummy3d, TRIM(myTracer), STATE%f2c_SSI_arr_map, _RC)
         endif
      enddo
#endif

! Compute Edge Heights
! --------------------

    zle(:,:,km+1) = phisxy(:,:)
    do k=km,1,-1
       zle(:,:,k) = zle(:,:,k+1) + cp*thv(:,:,k)*( pke(:,:,k+1)-pke(:,:,k) )
    enddo
       zle(:,:,:) = zle(:,:,:)/grav

    call FILLOUT3 (export, 'ZLE', zle, STATE, rc=status); VERIFY_(STATUS)

! Compute Mid-Layer Heights
! -------------------------

    call MAPL_GetPointer(export,temp3d,'ZL',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp3d)) then
       dummy3d = 0.5*( zle(:,:,2:) + zle(:,:,:km) )
       call SSI_CopyCoarseToFine(export, dummy3d, 'ZL', STATE%f2c_SSI_arr_map, _RC)
    endif

    pke = log(vars%pe)


! Fill Single Level Variables
! ---------------------------
 
    call MAPL_GetPointer(export,temp2d,'U200',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(dummy2d,ur,pke,log(20000.)  ,  status)
       VERIFY_(STATUS)
       call SSI_CopyCoarseToFine(export, dummy2d, 'U200', STATE%f2c_SSI_arr_map, rc=status)
       VERIFY_(STATUS)
    end if
 
    call MAPL_GetPointer(export,temp2d,'U250',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(dummy2d,ur,pke,log(25000.)  ,  status)
       VERIFY_(STATUS)
       call SSI_CopyCoarseToFine(export, dummy2d, 'U250', STATE%f2c_SSI_arr_map, rc=status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'U500',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(dummy2d,ur,pke,log(50000.)  ,  status)
       VERIFY_(STATUS)
       call SSI_CopyCoarseToFine(export, dummy2d, 'U500', STATE%f2c_SSI_arr_map, rc=status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'U700',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(dummy2d,ur,pke,log(70000.)  ,  status)
       VERIFY_(STATUS)
       call SSI_CopyCoarseToFine(export, dummy2d, 'U700', STATE%f2c_SSI_arr_map, rc=status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'U850',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(dummy2d,ur,pke,log(85000.)  ,  status)
       VERIFY_(STATUS)
       call SSI_CopyCoarseToFine(export, dummy2d, 'U850', STATE%f2c_SSI_arr_map, rc=status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'V200',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(dummy2d,vr,pke,log(20000.)  ,  status)
       VERIFY_(STATUS)
       call SSI_CopyCoarseToFine(export, dummy2d, 'V200', STATE%f2c_SSI_arr_map, rc=status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'V250',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(dummy2d,vr,pke,log(25000.)  ,  status)
       VERIFY_(STATUS)
       call SSI_CopyCoarseToFine(export, dummy2d, 'V250', STATE%f2c_SSI_arr_map, rc=status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'V500',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(dummy2d,vr,pke,log(50000.)  ,  status)
       VERIFY_(STATUS)
       call SSI_CopyCoarseToFine(export, dummy2d, 'V500', STATE%f2c_SSI_arr_map, rc=status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'V700',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(dummy2d,vr,pke,log(70000.)  ,  status)
       VERIFY_(STATUS)
       call SSI_CopyCoarseToFine(export, dummy2d, 'V700', STATE%f2c_SSI_arr_map, rc=status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'V850',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(dummy2d,vr,pke,log(85000.)  ,  status)
       VERIFY_(STATUS)
       call SSI_CopyCoarseToFine(export, dummy2d, 'V850', STATE%f2c_SSI_arr_map, rc=status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'T250',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(dummy2d,tempxy,pke,log(25000.)  , status)
       VERIFY_(STATUS)
       call SSI_CopyCoarseToFine(export, dummy2d, 'T250', STATE%f2c_SSI_arr_map, rc=status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'T300',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(dummy2d,tempxy,pke,log(30000.)  , status)
       VERIFY_(STATUS)
       call SSI_CopyCoarseToFine(export, dummy2d, 'T300', STATE%f2c_SSI_arr_map, rc=status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'T500',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(dummy2d,tempxy,pke,log(50000.)  ,  status)
       VERIFY_(STATUS)
       call SSI_CopyCoarseToFine(export, dummy2d, 'T500', STATE%f2c_SSI_arr_map, rc=status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'T700',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(dummy2d,tempxy,pke,log(70000.)  ,  status)
       VERIFY_(STATUS)
       call SSI_CopyCoarseToFine(export, dummy2d, 'T700', STATE%f2c_SSI_arr_map, rc=status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'T850',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(dummy2d,tempxy,pke,log(85000.)  ,  status)
       VERIFY_(STATUS)
       call SSI_CopyCoarseToFine(export, dummy2d, 'T850', STATE%f2c_SSI_arr_map, rc=status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'Q250',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(dummy2d,qv,pke,log(25000.)  ,  status)
       VERIFY_(STATUS)
       call SSI_CopyCoarseToFine(export, dummy2d, 'Q250', STATE%f2c_SSI_arr_map, rc=status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'Q500',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(dummy2d,qv,pke,log(50000.)  ,  status)
       VERIFY_(STATUS)
       call SSI_CopyCoarseToFine(export, dummy2d, 'Q500', STATE%f2c_SSI_arr_map, rc=status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'Q850',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(dummy2d,qv,pke,log(85000.)  ,  status)
       VERIFY_(STATUS)
       call SSI_CopyCoarseToFine(export, dummy2d, 'Q850', STATE%f2c_SSI_arr_map, rc=status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'Z700',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(dummy2d,zle*grav,pke,log(70000.)  , status)
       VERIFY_(STATUS)
       call SSI_CopyCoarseToFine(export, dummy2d, 'Z700', STATE%f2c_SSI_arr_map, rc=status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'Z500',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(dummy2d,zle*grav,pke,log(50000.)  , status)
       VERIFY_(STATUS)
       call SSI_CopyCoarseToFine(export, dummy2d, 'Z500', STATE%f2c_SSI_arr_map, rc=status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'Z300',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(dummy2d,zle*grav,pke,log(30000.)  , status)
       VERIFY_(STATUS)
       call SSI_CopyCoarseToFine(export, dummy2d, 'Z300', STATE%f2c_SSI_arr_map, rc=status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'H250',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(dummy2d,zle,pke,log(25000.)  , status)
       VERIFY_(STATUS)
       call SSI_CopyCoarseToFine(export, dummy2d, 'H250', STATE%f2c_SSI_arr_map, rc=status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'H300',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(dummy2d,zle,pke,log(30000.)  , status)
       VERIFY_(STATUS)
       call SSI_CopyCoarseToFine(export, dummy2d, 'H300', STATE%f2c_SSI_arr_map, rc=status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'H500',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(dummy2d,zle,pke,log(50000.)  , status)
       VERIFY_(STATUS)
       call SSI_CopyCoarseToFine(export, dummy2d, 'H500', STATE%f2c_SSI_arr_map, rc=status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'H700',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(dummy2d,zle,pke,log(70000.)  , status)
       VERIFY_(STATUS)
       call SSI_CopyCoarseToFine(export, dummy2d, 'H700', STATE%f2c_SSI_arr_map, rc=status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'H850',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(dummy2d,zle,pke,log(85000.)  , status)
       VERIFY_(STATUS)
       call SSI_CopyCoarseToFine(export, dummy2d, 'H850', STATE%f2c_SSI_arr_map, rc=status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'H1000',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(dummy2d,zle,pke,log(100000.)  , status)
       VERIFY_(STATUS)
       call SSI_CopyCoarseToFine(export, dummy2d, 'H1000', STATE%f2c_SSI_arr_map, rc=status)
       VERIFY_(STATUS)
    end if

! Fill Model Top Level Variables
! ---------------------------------------
    call MAPL_GetPointer(export,temp2d,'UTOP', rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       dummy2d = ur(:,:,1)
       call SSI_CopyCoarseToFine(export, dummy2d, 'UTOP', STATE%f2c_SSI_arr_map, rc=status)
       VERIFY_(STATUS)
    endif

    call MAPL_GetPointer(export,temp2d,'VTOP', rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       dummy2d = vr(:,:,1)
       call SSI_CopyCoarseToFine(export, dummy2d, 'VTOP', STATE%f2c_SSI_arr_map, rc=status)
       VERIFY_(STATUS)
    endif

    call MAPL_GetPointer(export,temp2d,'TTOP', rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       dummy2d = tempxy(:,:,1)
       call SSI_CopyCoarseToFine(export, dummy2d, 'TTOP', STATE%f2c_SSI_arr_map, rc=status)
       VERIFY_(STATUS)
    endif

    call MAPL_GetPointer(export,temp2d,'DELPTOP', rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       dummy2d = dp(:,:,1)
       call SSI_CopyCoarseToFine(export, dummy2d, 'DELPTOP', STATE%f2c_SSI_arr_map, rc=status)
       VERIFY_(STATUS)
    endif

! Compute Heights Above Surface
! -----------------------------
    do k=1,km+1
    zle(:,:,k) = zle(:,:,k) - zle(:,:,km+1)
    enddo
    
    call MAPL_GetPointer(export,temp2d,'U50M',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(dummy2d,ur,-zle,-50., status)
       VERIFY_(STATUS)
       call SSI_CopyCoarseToFine(export, dummy2d, 'U50M', STATE%f2c_SSI_arr_map, rc=status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'V50M',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(dummy2d,vr,-zle,-50., status)
       VERIFY_(STATUS)
       call SSI_CopyCoarseToFine(export, dummy2d, 'V50M', STATE%f2c_SSI_arr_map, rc=status)
       VERIFY_(STATUS)
    end if

! Compute Surface Pressure
! ------------------------

    call MAPL_GetPointer(export,temp2d,'PS',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       dummy2d = vars%pe(:,:,km+1)
       call SSI_CopyCoarseToFine(export, dummy2d, 'PS', STATE%f2c_SSI_arr_map, rc=status)
       VERIFY_(STATUS)
    endif

! Compute Vertically Averaged T,U
! -------------------------------
    call MAPL_GetPointer(export,temp2d,'TAVE',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       dummy2d = 0.0
       do k=1,km
       dummy2d = dummy2d + tempxy(:,:,k)*dp(:,:,k)
       enddo
       dummy2d = dummy2d / (vars%pe(:,:,km+1)-vars%pe(:,:,1))
       call SSI_CopyCoarseToFine(export, dummy2d, 'TAVE', STATE%f2c_SSI_arr_map, rc=status)
       VERIFY_(STATUS)
    endif

    call MAPL_GetPointer(export,temp2d,'UAVE',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       dummy2d = 0.0
       do k=1,km
       dummy2d = dummy2d + ur(:,:,k)*dp(:,:,k)
       enddo
       dummy2d = dummy2d / (vars%pe(:,:,km+1)-vars%pe(:,:,1))
       call SSI_CopyCoarseToFine(export, dummy2d, 'UAVE', STATE%f2c_SSI_arr_map, rc=status)
       VERIFY_(STATUS)
    endif

! Convert T to Tv
! ---------------

    tempxy = tempxy*(1.0+eps*qv)

    call MAPL_GetPointer(export,temp3d,'TV',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp3d)) then
       dummy3d = tempxy
       call SSI_CopyCoarseToFine(export, dummy3d, 'TV', STATE%f2c_SSI_arr_map, rc=status)
       VERIFY_(STATUS)
    endif

! Compute Sea-Level Pressure
! --------------------------
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

       if(associated(temp2d)) then
          dummy2d = slp
          call SSI_CopyCoarseToFine(export, dummy2d, 'SLP', STATE%f2c_SSI_arr_map, rc=status)
          VERIFY_(STATUS)
       endif
       ! first copy current ztemp1,2,3 from fine to coarse so that
       ! the 'where' statement can be properly executed and then copy
       ! back from corase to fine.
       !if(associated(ztemp1)) where( ztemp1.eq.MAPL_UNDEF ) ztemp1 = H1000
       if(associated(ztemp1)) then
          call SSI_CopyFineToCoarse(export, dummy2d, 'H1000', STATE%f2c_SSI_arr_map, rc=status)
          VERIFY_(STATUS)
          where( dummy2d.eq.MAPL_UNDEF ) dummy2d = H1000
          call SSI_CopyCoarseToFine(export, dummy2d, 'H1000', STATE%f2c_SSI_arr_map, rc=status)
          VERIFY_(STATUS)
       endif
       !if(associated(ztemp2)) where( ztemp2.eq.MAPL_UNDEF ) ztemp2 = H850
       if(associated(ztemp2)) then
          call SSI_CopyFineToCoarse(export, dummy2d, 'H850', STATE%f2c_SSI_arr_map, rc=status)
          VERIFY_(STATUS)
          where( dummy2d.eq.MAPL_UNDEF ) dummy2d = H850
          call SSI_CopyCoarseToFine(export, dummy2d, 'H850', STATE%f2c_SSI_arr_map, rc=status)
          VERIFY_(STATUS)
       endif
       !if(associated(ztemp3)) where( ztemp3.eq.MAPL_UNDEF ) ztemp3 = H500
       if(associated(ztemp3)) then
          call SSI_CopyFineToCoarse(export, dummy2d, 'H500', STATE%f2c_SSI_arr_map, rc=status)
          VERIFY_(STATUS)
          where( dummy2d.eq.MAPL_UNDEF ) dummy2d = H500
          call SSI_CopyCoarseToFine(export, dummy2d, 'H500', STATE%f2c_SSI_arr_map, rc=status)
          VERIFY_(STATUS)
       endif
       DEALLOCATE(slp,H1000,H850,H500)
    end if

! Deallocate Memory
! -----------------

    DEALLOCATE(  kenrg )
    DEALLOCATE(  penrg )
    DEALLOCATE(  tenrg )
    DEALLOCATE( kenrg0 )
    DEALLOCATE( penrg0 )
    DEALLOCATE( tenrg0 )
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
    DEALLOCATE(    pke )
    DEALLOCATE(  logpl )
    DEALLOCATE(  logpe )
    DEALLOCATE(  logps )
    DEALLOCATE(    zle )
    DEALLOCATE( dthdtphyint1 )
    DEALLOCATE( dthdtphyint2 )

    !call freeTracers(state)

    end if ! .not. SW_DYNAMICS

    !call MAPL_TimerOff(GENSTATE,"RUN2")
    !call MAPL_TimerOff(GENSTATE,"TOTAL")

    RETURN_(ESMF_SUCCESS)
end subroutine RunAddIncs

!-----------------------------------------------------------------------
  subroutine ADD_INCS ( STATE,IMPORT,internal,DT,IS_WEIGHTED,RC )

   use fms_mod, only: set_domain, nullify_domain
   use fv_diagnostics_mod, only: prt_maxmin
   use time_manager_mod,   only: time_type
   use fv_update_phys_mod, only: fv_update_phys
!
! !INPUT PARAMETERS:

   type(DynState), pointer                :: STATE
   type(ESMF_State),       intent(INOUT)  :: IMPORT
   type(ESMF_State),       intent(INOUT)  :: internal
   real(FVPRC),            intent(IN   )  :: DT
   integer,  optional,     intent(OUT  )  :: RC
   logical,  optional,     intent(IN   )  :: is_weighted

!
! !DESCRIPTION:  This routine adds the tendencies to the state,
!                weighted appropriately by the time step.  Temperature
!                tendencies are pressure weighted (ie., DELP*DT/Dt).
!                All tendencies are on the A-grid, and have an XY decomposition.
!

    integer               :: status
    logical               :: is_weighted_

    integer               :: is ,ie , js ,je , km
    integer               :: isd,ied, jsd,jed
    real(r4), allocatable :: fvQOLD(:,:,:), QTEND(:,:,:)
    real(r8), allocatable :: DPNEW(:,:,:),DPOLD(:,:,:)

    real(REAL8), allocatable :: tend_ua(:,:,:), tend_va(:,:,:)
    real(REAL8), allocatable :: tend_un(:,:,:), tend_vn(:,:,:)

    real(FVPRC), allocatable :: u_dt(:,:,:), v_dt(:,:,:), t_dt(:,:,:)

    real(kind=4), pointer :: tend(:,:,:)
    real(kind=4), pointer :: tend_kp1(:,:,:) => Null()

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
       if (nwat_tracers >= 5) nwat = 1 ! STATE has QV only
       if (.not. HYDROSTATIC) then
          if (nwat_tracers >= 5) nwat = 3 ! STATE has QV, QLIQ, QICE
          if (nwat_tracers == 8) nwat = 6 ! STATE has QV, QLIQ, QICE, QRAIN, QSNOW, QGRAUPEL
       endif
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
    call PULL_Q ( STATE, IMPORT, qqq, NXQ, InFieldName='Q', RC=rc )
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
    call PULL_Q ( STATE, IMPORT, qqq, NXQ, InFieldName='QLLS', RC=rc )
    if ( (qqq%is_r4) .and. (associated(qqq%content_r4)) ) then
       if (size(Q(:,:,:,liq_wat))==size(qqq%content_r4)) Q(:,:,:,liq_wat) = Q(:,:,:,liq_wat) + qqq%content_r4
    elseif (associated(qqq%content)) then
       if (size(Q(:,:,:,liq_wat))==size(qqq%content)) Q(:,:,:,liq_wat) = Q(:,:,:,liq_wat) + qqq%content
    endif
    call PULL_Q ( STATE, IMPORT, qqq, NXQ, InFieldName='QLCN', RC=rc )
    if ( (qqq%is_r4) .and. (associated(qqq%content_r4)) ) then
       if (size(Q(:,:,:,liq_wat))==size(qqq%content_r4)) Q(:,:,:,liq_wat) = Q(:,:,:,liq_wat) + qqq%content_r4
    elseif (associated(qqq%content)) then
       if (size(Q(:,:,:,liq_wat))==size(qqq%content)) Q(:,:,:,liq_wat) = Q(:,:,:,liq_wat) + qqq%content
    endif
    ! Grab QICE from imports
    call PULL_Q ( STATE, IMPORT, qqq, NXQ, InFieldName='QILS', RC=rc )
    if ( (qqq%is_r4) .and. (associated(qqq%content_r4)) ) then
       if (size(Q(:,:,:,ice_wat))==size(qqq%content_r4)) Q(:,:,:,ice_wat) = Q(:,:,:,ice_wat) + qqq%content_r4
    elseif (associated(qqq%content)) then
       if (size(Q(:,:,:,ice_wat))==size(qqq%content)) Q(:,:,:,ice_wat) = Q(:,:,:,ice_wat) + qqq%content
    endif
    call PULL_Q ( STATE, IMPORT, qqq, NXQ, InFieldName='QICN', RC=rc )
    if ( (qqq%is_r4) .and. (associated(qqq%content_r4)) ) then
       if (size(Q(:,:,:,ice_wat))==size(qqq%content_r4)) Q(:,:,:,ice_wat) = Q(:,:,:,ice_wat) + qqq%content_r4
    elseif (associated(qqq%content)) then
       if (size(Q(:,:,:,ice_wat))==size(qqq%content)) Q(:,:,:,ice_wat) = Q(:,:,:,ice_wat) + qqq%content
    endif
    endif
    if (nwat >= 6) then
    ! Grab RAIN from imports
    call PULL_Q ( STATE, IMPORT, qqq, NXQ, InFieldName='QRAIN', RC=rc )
    if ( (qqq%is_r4) .and. (associated(qqq%content_r4)) ) then
       if (size(Q(:,:,:,rainwat))==size(qqq%content_r4)) Q(:,:,:,rainwat) = qqq%content_r4
    elseif (associated(qqq%content)) then
       if (size(Q(:,:,:,rainwat))==size(qqq%content)) Q(:,:,:,rainwat) = qqq%content
    endif
    ! Grab SNOW from imports
    call PULL_Q ( STATE, IMPORT, qqq, NXQ, InFieldName='QSNOW', RC=rc )
    if ( (qqq%is_r4) .and. (associated(qqq%content_r4)) ) then
       if (size(Q(:,:,:,snowwat))==size(qqq%content_r4)) Q(:,:,:,snowwat) = qqq%content_r4
    elseif (associated(qqq%content)) then
       if (size(Q(:,:,:,snowwat))==size(qqq%content)) Q(:,:,:,snowwat) = qqq%content
    endif
    ! Grab GRAUPEL from imports
    call PULL_Q ( STATE, IMPORT, qqq, NXQ, InFieldName='QGRAUPEL', RC=rc )
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
       ! ****        while IMPORT Tendencies are on the A-Grid             ****
       ! **********************************************************************

       ALLOCATE( tend_ua(is:ie  ,js:je  ,km) )
       ALLOCATE( tend_va(is:ie  ,js:je  ,km) )
       ALLOCATE( tend_un(is:ie  ,js:je+1,km) )
       ALLOCATE( tend_vn(is:ie+1,js:je  ,km) )
       !call ESMFL_StateGetPointerToData ( IMPORT,TEND,'DUDT',RC=STATUS )
       !VERIFY_(STATUS)
       if(.not.associated(tend)) then
          allocate(tend(is:ie,js:je,km), stat=status)
          VERIFY_(STATUS)
       endif
       call SSI_CopyFineToCoarse(import, tend, 'DUDT', STATE%f2c_SSI_arr_map, _RC)

       tend_ua(is:ie,js:je,1:km) = tend

       !call ESMFL_StateGetPointerToData ( IMPORT,TEND,'DVDT',RC=STATUS )
       !VERIFY_(STATUS)
       call SSI_CopyFineToCoarse(import, tend, 'DVDT', STATE%f2c_SSI_arr_map, _RC)

       tend_va(is:ie,js:je,1:km) = tend

      !if (.not. HYDROSTATIC ) then
      !  call ESMFL_StateGetPointerToData ( IMPORT,TEND,'DWDT',RC=STATUS )
      !  VERIFY_(STATUS)
      !  STATE%VARS%W = STATE%VARS%W + DT*TEND(is:ie,js:je,1:km)
      !endif

       ! Put the wind tendencies on the Native Dynamics grid
       ! ---------------------------------------------------
       call Agrid_To_Native( tend_ua, tend_va, tend_un, tend_vn, &
                             wind_increment_limiter = 800.d0/86400.d0 )


       ! Add the wind tendencies to the control variables
       ! ------------------------------------------------
       STATE%VARS%U = STATE%VARS%U + DT*TEND_UN(is:ie,js:je,1:km)
       STATE%VARS%V = STATE%VARS%V + DT*TEND_VN(is:ie,js:je,1:km)

       DEALLOCATE( tend_ua )
       DEALLOCATE( tend_va )
       DEALLOCATE( tend_un )
       DEALLOCATE( tend_vn )

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

       !call ESMFL_StateGetPointerToData ( IMPORT,TEND,'DPEDT',RC=STATUS )
       !VERIFY_(STATUS)
       if(.not.associated(tend_kp1)) then
          allocate(tend_kp1(is:ie,js:je,km+1), stat=status)
          VERIFY_(STATUS)
       endif
       call SSI_CopyFineToCoarse(import, tend_kp1, 'DPEDT', STATE%f2c_SSI_arr_map, _RC)

       STATE%VARS%PE = STATE%VARS%PE + DT*TEND_kp1
       !STATE%VARS%PE = STATE%VARS%PE + DT*TEND

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
       ! ****        IMPORT Variable is a) D/Dt (T)     , IS_WEIGHTED=.F. ****
       ! ****                           b) D/Dt (T*DELP), IS_WEIGHTED=.T. ****
       ! *********************************************************************

       !call ESMFL_StateGetPointerToData ( IMPORT,TEND,'DTDT',RC=STATUS )
       !VERIFY_(STATUS)
       call SSI_CopyFineToCoarse(import, tend, 'DTDT', STATE%f2c_SSI_arr_map, _RC)

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


    call INTERNAL_CoarseToFine(STATE, internal, _RC)

    if (ALLOCATED(Q  )) DEALLOCATE( Q   )
    if (ALLOCATED(CVM)) DEALLOCATE( CVM )

   return

 end subroutine ADD_INCS



   subroutine FILLOUT3r8(export, name, V, STATE, RC)
     type (ESMF_State),  intent(inout) :: export
     character(len=*),   intent(IN   ) :: name
     real(r8),           intent(IN   ) :: V(:,:,:)
     integer, optional,  intent(  out) :: rc
     type(DynState), pointer :: STATE

     real(r8), pointer          :: CPL(:,:,:)
     integer                    :: status
     character(len=ESMF_MAXSTR) :: IAm="CoarseFillout3r8"

     call MAPL_GetPointer(export, cpl, name, RC=STATUS)
     VERIFY_(STATUS)
     if(associated(cpl)) then
        call SSI_CopyCoarseToFine(export, v, name, STATE%f2c_SSI_arr_map, rc=status)
        VERIFY_(STATUS)
     endif

   end subroutine FILLOUT3r8

   subroutine FILLOUT3(export, name, V, STATE, RC)
     type (ESMF_State),  intent(inout) :: export
     character(len=*),   intent(IN   ) :: name
     real(r8),           intent(IN   ) :: V(:,:,:)
     integer, optional,  intent(  out) :: rc
     type(DynState), pointer :: STATE

     real(r4), pointer          :: CPL(:,:,:)
     real(r4), pointer          :: dummy(:,:,:)
     integer                    :: status
     character(len=ESMF_MAXSTR) :: IAm="CoarseFillout3"
     integer :: dimen(3)

     call MAPL_GetPointer(export, cpl, name, RC=STATUS)
     VERIFY_(STATUS)
     if(associated(cpl)) then
        dimen = shape(V)
        allocate(dummy(dimen(1),dimen(2),dimen(3)), stat=status)
        VERIFY_(STATUS)
        dummy = V
        call SSI_CopyCoarseToFine(export, dummy, name, STATE%f2c_SSI_arr_map, rc=status)
        VERIFY_(STATUS)
        deallocate(dummy, stat=status)
        VERIFY_(STATUS)
     endif

   end subroutine FILLOUT3

!-----------------------------------------------------------------------

   subroutine FILLOUT2(export, name, V, STATE, rc)
     type (ESMF_State),  intent(inout) :: export
     character(len=*),   intent(IN   ) :: name
     real(r8),           intent(IN   ) :: V(:,:)
     integer, optional,  intent(  out) :: rc
     type(DynState), pointer :: STATE

     real(kind=4), pointer      :: CPL(:,:)
     real(r4), pointer          :: dummy(:,:)
     integer                    :: status
     character(len=ESMF_MAXSTR) :: IAm="CoarseFillout2"
     integer :: dimen(2)

     call MAPL_GetPointer(export, cpl, name, RC=STATUS)
     VERIFY_(STATUS)
     if(associated(cpl)) then
        dimen = shape(V)
        allocate(dummy(dimen(1),dimen(2)), stat=status)
        VERIFY_(STATUS)
        dummy = V
        call SSI_CopyCoarseToFine(export, dummy, name, STATE%f2c_SSI_arr_map, rc=status)
        VERIFY_(STATUS)
        deallocate(dummy, stat=status)
        VERIFY_(STATUS)
     endif

     return
   end subroutine FILLOUT2

!-----------------------------------------------------------------------

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
! --------------------------
    pke  = ple**kappa
    phiT = phiS
    if( present(gze) ) gze(:,:,km+1) = phiS
    do k=km,1,-1
                       phiT = phiT + cp*thv(:,:,k)*( pke(:,:,k+1)-pke(:,:,k) )
    if( present(gze) ) gze(:,:,k) = phiT
    enddo

! Compute Energetics:  Cp*Tv + K + PHI
! ------------------------------------
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP

! !IROUTINE: Finalize

! !DESCRIPTION: Writes restarts and cleans-up through MAPL\_GenericFinalize and
!   deallocates memory from the Private Internal state.
!
! !INTERFACE:

subroutine Finalize(gc, import, export, clock, rc)

! !ARGUMENTS:

    type(ESMF_GridComp) :: gc
    type (ESMF_State) :: import
    type (ESMF_State) :: export
    type (ESMF_Clock) :: clock
    integer, intent(out) :: rc

!EOP

! Local variables
    type (DYN_wrap) :: wrap
    type (DynState), pointer  :: STATE

    character(len=ESMF_MAXSTR)        :: IAm
    character(len=ESMF_MAXSTR)        :: COMP_NAME
    integer                           :: status

    !type (MAPL_MetaComp),     pointer :: MAPL
    type (ESMF_Config)                :: cf


! BEGIN

    Iam = "CoarseFinalize"
    call ESMF_GridCompGet( GC, name=COMP_NAME, config=cf, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

! Retrieve the pointer to the state
! ---------------------------------

    !call MAPL_GetObjectFromGC (GC, MAPL,  RC=STATUS )
    !VERIFY_(STATUS)

    !call MAPL_TimerOn(MAPL,"TOTAL")
    !call MAPL_TimerOn(MAPL,"FINALIZE")

! Retrieve the pointer to the state
!----------------------------------

    call ESMF_UserCompGetInternalState(gc, 'DYNstate', wrap, status)
    VERIFY_(STATUS)

    state => wrap%dyn_state

    call DynFinalize( STATE )
! Call Generic Finalize
!----------------------

    !call MAPL_TimerOff(MAPL,"FINALIZE")
    !call MAPL_TimerOff(MAPL,"TOTAL")

    !call MAPL_GenericFinalize ( GC, IMPORT, EXPORT, CLOCK,  RC=STATUS)
    !VERIFY_(STATUS)

  RETURN_(ESMF_SUCCESS)

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine PRINT_TIMES(TIMES,DAYS)
  integer(kind=8), intent(INOUT) :: TIMES(:,:)
  real(r8),        intent(IN   ) :: DAYS
  TIMES = 0

  return
 end subroutine PRINT_TIMES

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
       tref = ( tv(k1)*dp1 + tv(k2)*dp2 ) / (dp1+dp2)
       pref = 0.5 * ( pe(k_bot+1) + pe(k_bot-1) )
      tstar = tref*( ps/pref )**yfactor

      slp   = ps*( 1.0+gg*phis/tstar )**factor
      H1000 = (phis/MAPL_grav) - (tstar/gamma)*((100000.0/ps)**(1./factor)-1.0)
      H850  = (phis/MAPL_grav) - (tstar/gamma)*(( 85000.0/ps)**(1./factor)-1.0)
      H500  = (phis/MAPL_grav) - (tstar/gamma)*(( 50000.0/ps)**(1./factor)-1.0)
      return
  end subroutine get_slp

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine VertInterp(v2,v3,ple,pp,rc)

    real(r4), intent(OUT) :: v2(:,:)
    real(r8), intent(IN ) :: v3(:,:,:)
    real(r8), intent(IN ) :: ple(:,:,:)
    real    , intent(IN ) :: pp
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

    RETURN_(ESMF_SUCCESS)
  end subroutine VertInterp


!BOP

! !IROUTINE: Coldstart

! !DESCRIPTION:
!   Routine to coldstart from an isothermal state of rest.
!   The temperature can be specified in the config, otherwise
!   it is 300K. The surface pressure is assumed to be 1000 hPa.
!
! !INTERFACE:

subroutine Coldstart(gc, import, export, clock, rc)

    USE sw, only : sw_phis=>surface_geopotential
    USE sw, only : sw_hght=>height
    USE sw, only : sw_uwnd=>u_wind
    USE sw, only : sw_vwnd=>v_wind
    USE jw, only : temperature, u_wind, v_wind, surface_geopotential
    USE jw, only : tracer_q, tracer_q1_q2, tracer_q3
    USE testcases_3_4_5_6, only : advection, Rossby_Haurwitz, mountain_Rossby, gravity_wave

! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: gc
    type(ESMF_State),    intent(inout) :: import
    type(ESMF_State),    intent(inout) :: export
    type (ESMF_Clock),   intent(inout) :: clock
    integer, intent(out), optional     :: rc

!EOP

    character(len=ESMF_MAXSTR)        :: IAm="Coldstart"
    character(len=ESMF_MAXSTR)        :: COMP_NAME
    integer                           :: status

    type (MAPL_MetaComp),     pointer :: MAPL
    type (ESMF_State)                 :: INTERNAL

    real(REAL8), pointer                 :: AK(:), BK(:)
    real(REAL8), pointer                 :: U      (:,:,:)
    real(REAL8), pointer                 :: V      (:,:,:)
    real(REAL8), pointer                 :: PT     (:,:,:)
    real(REAL8), pointer                 :: PE     (:,:,:)
    real(REAL8), pointer                 :: PKZ    (:,:,:)
    real(kind=4), pointer             :: phis   (:,:)
    real(REAL4), pointer                     :: LONS   (:,:)
    real(REAL4), pointer                     :: LATS   (:,:)
    real                              :: T0
    integer                           :: L
    type(ESMF_Config)                 :: CF
    integer                           :: i,j,k,n
    integer                           :: IS,IE, JS,JE, KS,KE, IM,JM,KM, LS

    integer                     :: case_id
    integer                     :: case_rotation
    integer                     :: case_tracers

    real(REAL8) :: dummy_1, dummy_2, dummy_3, dummy_4, dummy_5, dummy_6
    real(REAL8) :: dz, ztop, height, pressure
    real(REAL8) :: LONc,LATc
    real(REAL8) :: eta, eta_top, rot_ang
    real(REAL8) :: ptop, pint
    real(REAL8), allocatable :: PS(:,:)
    logical :: perturb
    logical :: ak_is_missing = .false.
    logical :: bk_is_missing = .false.
    integer :: FV3_STANDALONE

    type (DYN_wrap) :: wrap
    type (DynState), pointer :: STATE
    type (DynGrid),  pointer :: GRID

    logical :: isPresent

! Tracer Stuff
    real(r4), pointer                :: TRACER(:,:,:)
    real(REAL8), allocatable            :: Q5(:,:,:)
    real(REAL8), allocatable            :: Q6(:,:,:)
    type (ESMF_Grid)                 :: esmfGRID
    type (ESMF_FieldBundle)          :: TRADV_BUNDLE
    character(len=ESMF_MAXSTR)       :: FIELDNAME
    character(len=ESMF_MAXSTR)       :: STRING
    real(REAL8), parameter    :: r0_6=0.6
    real(REAL8), parameter    :: r1_0=1.0
    integer :: NQ
    type (ESMF_GridComp) :: fineGC
    integer, allocatable :: gcImg(:)
    integer :: itemCount

! Begin

    call ESMF_GridCompGet( GC, name=COMP_NAME, CONFIG=CF, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // trim(Iam)

! Retrieve fine GC
! ---------------------------------
      call ESMF_AttributeGet(GC, name='GC_IMAGE', itemCount=itemCount, rc=status)
      VERIFY_(STATUS)
      allocate(gcImg(itemCount), stat=status)
      VERIFY_(STATUS)
      call ESMF_AttributeGet(GC, name='GC_IMAGE', valueList=gcImg, rc=status)
      VERIFY_(STATUS)
      fineGC = transfer(gcImg, fineGC)
      deallocate(gcImg, stat=status)
      VERIFY_(STATUS)

! Retrieve the pointer to the state
! ---------------------------------

    call MAPL_GetObjectFromGC (fineGC, MAPL,  RC=STATUS )
    VERIFY_(STATUS)

    call ESMF_UserCompGetInternalState(GC, 'DYNstate', wrap, status)
    VERIFY_(STATUS)
    state => wrap%dyn_state
    grid  => state%grid   ! direct handle to grid

    IS = FV_Atm(1)%bd%isc
    IE = FV_Atm(1)%bd%iec
    JS = FV_Atm(1)%bd%jsc
    JE = FV_Atm(1)%bd%jec
    KS = 1
    KE = FV_Atm(1)%npz
    KM = KE-KS+1

!BOR
! !RESOURCE_ITEM: K :: Value of isothermal temperature on coldstart
    call MAPL_GetResource ( MAPL, T0, 'T0:', default=273., RC=STATUS )
    VERIFY_(STATUS)
!EOR
    !call MAPL_Get ( MAPL,                &
    !       INTERNAL_ESMF_STATE=INTERNAL, &
    !       lats = LATS,                  &
    !       lons = LONS,                  &
    !                           RC=STATUS )
    !VERIFY_(STATUS)

   allocate(LONS(is:ie,js:je), stat=status)
   VERIFY_(STATUS)
   call SSI_CopyFineToCoarse(INTERNAL, LONS, 'LONS', STATE%f2c_SSI_arr_map, rc=status)
   VERIFY_(STATUS)
   allocate(LATS(is:ie,js:je), stat=status)
   VERIFY_(STATUS)
   call SSI_CopyFineToCoarse(INTERNAL, LATS, 'LATS', STATE%f2c_SSI_arr_map, rc=status)
   VERIFY_(STATUS)

   if (FV_Atm(1)%flagstruct%grid_type == 4) then
    ! Doubly-Period setup based on first LAT/LON coordinate
     LONS(:,:) =  0.0
     LATS(:,:) = 15.0*PI/180.0
   endif

!! A-Grid U Wind
!        call MAPL_GetPointer(Internal,U,'U'  ,rc=STATUS)
!        VERIFY_(STATUS)
!! A-Grid V Wind
!        call MAPL_GetPointer(Internal,V,'V'  ,rc=STATUS)
!! Surface Geopotential
!        call MAPL_GetPointer ( IMPORT, phis, 'PHIS', RC=STATUS )
!        VERIFY_(STATUS)
!! Potential-Temperature
!        call MAPL_GetPointer(Internal,PT,'PT',rc=STATUS)
!        VERIFY_(STATUS)
!! Edge Pressures
!        call MAPL_GetPointer(Internal,PE  ,'PE',rc=STATUS)
!        VERIFY_(STATUS)
!! Presssure ^ kappa at mid-layers
!        call MAPL_GetPointer(Internal,PKZ ,'PKZ',rc=STATUS)
!        VERIFY_(STATUS)
! AK and BK for vertical coordinate
        call MAPL_GetPointer(Internal,ak  ,'AK' ,rc=STATUS)
        VERIFY_(STATUS)
        call MAPL_GetPointer(Internal,bk  ,'BK' ,rc=STATUS)
        VERIFY_(STATUS)

   allocate(U(is:ie,js:je,1:km), stat=status)
   VERIFY_(STATUS)
   allocate(V(is:ie,js:je,1:km), stat=status)
   VERIFY_(STATUS)
   allocate(PT(is:ie,js:je,1:km), stat=status)
   VERIFY_(STATUS)
   allocate(PE(is:ie,js:je,0:km), stat=status)
   VERIFY_(STATUS)
   allocate(PKZ(is:ie,js:je,1:km), stat=status)
   VERIFY_(STATUS)
   allocate(phis(is:ie,js:je), stat=status)
   VERIFY_(STATUS)

    U = 0.0

    !IS = lbound(U,1)
    !IE = ubound(U,1)
    !JS = lbound(U,2)
    !JE = ubound(U,2)
    !KS = lbound(U,3)
    !KE = ubound(U,3)
    !KM = KE-KS+1

    ALLOCATE( PS(IS:IE,JS:JE) )

    call ESMF_ConfigGetAttribute( cf, IM, label='IM:', default=0 , rc = rc )
    call ESMF_ConfigGetAttribute( cf, JM, label='JM:', default=0 , rc = rc )

  if (KM<=2) then   ! Shallow Water

      call ESMF_ConfigGetAttribute( cf, case_id, label='CASE_ID:', default=1 , rc = rc )
      DYN_CASE = case_id

       do j=JS,JE
          do i=IS,IE
             LONc = LONS(i,j)
             LATc = LATS(i,j)
             U(i,j,1)  = sw_uwnd(LONc,LATc,case_id)
             V(i,j,1)  = sw_vwnd(LONc,LATc,case_id)
             PE(i,j,0) = sw_phis(LONc,LATc,case_id)
             PE(i,j,1) = sw_hght(LONc,LATc,case_id)
             phis(i,j) = PE(i,j,0)
          enddo
       enddo

  else              ! 3-D Baroclinic

    U(IS:IE,JS:JE,KE) = .001*abs(lats(:,:))
    V = 0.0

    call ESMF_ConfigFindLabel( cf, 'AK:', isPresent=isPresent, rc = status )
    VERIFY_(STATUS)
    if (isPresent) then
       do L = 0, SIZE(AK)-1
          call ESMF_ConfigNextLine  ( CF, rc=STATUS )
          call ESMF_ConfigGetAttribute( cf, AK(L), rc = status )
          VERIFY_(STATUS)
       enddo
    else
       ak_is_missing = .true.
    endif

    call ESMF_ConfigFindLabel( cf, 'BK:', isPresent=isPresent, rc = status )
    VERIFY_(STATUS)
    if (isPresent) then
       do L = 0, SIZE(bk)-1
          call ESMF_ConfigNextLine  ( CF, rc=STATUS )
          call ESMF_ConfigGetAttribute( cf, BK(L), rc = status )
          VERIFY_(STATUS)
       enddo
    else
       bk_is_missing = .true.
    endif

    if (ak_is_missing .or. bk_is_missing) call set_eta(km, ls, ptop, pint, AK, BK)

    _ASSERT(ANY(AK /= 0.0) .or. ANY(BK /= 0.0),'needs informative message')
    do L=lbound(PE,3),ubound(PE,3)
       PE(:,:,L) = AK(L) + BK(L)*MAPL_P00
    enddo

    PKZ = 0.5*(PE(:,:,lbound(PE,3)  :ubound(PE,3)-1) + &
               PE(:,:,lbound(PE,3)+1:ubound(PE,3)  ) )
    PKZ = PKZ**MAPL_KAPPA

    PT = T0/PKZ

! Check if running standalone model
    call ESMF_ConfigGetAttribute ( CF, FV3_STANDALONE, Label="FV3_STANDALONE:", default=0, RC=STATUS)
    VERIFY_(STATUS)

! 3D Baroclinic Test Cases

    call ESMF_ConfigGetAttribute( cf, case_id      , label='CASE_ID:'      , default=0 , rc = rc )
    call ESMF_ConfigGetAttribute( cf, case_rotation, label='CASE_ROTATION:', default=0 , rc = rc )
    call ESMF_ConfigGetAttribute( cf, case_tracers , label='CASE_TRACERS:' , default=1234, rc = rc )
    DYN_CASE = case_id

    write(STRING,'(A,I5,A)') "Initializing CASE_ID ", case_id, " in FVcubed:"
    call WRITE_PARALLEL( trim(STRING) )


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
      do k=KS,KE
         eta = 0.5*( (ak(k-1)+ak(k))/1.e5 + bk(k-1)+bk(k) )
         do j=JS,JE
            do i=IS,IE
               LONc = LONS(i,j)
               LATc = LATS(i,j)
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
      do k=KS,KE
         eta = 0.5*( (ak(k-1)+ak(k))/1.e5 + bk(k-1)+bk(k) )
         do j=JS,JE
            do i=IS,IE
               LONc = LONS(i,j)
               LATc = LATS(i,j)
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

      allocate( Q5(IS:IE, JS:JE, 0:KM-1), STAT=STATUS)
      VERIFY_(STATUS)
      allocate( Q6(IS:IE, JS:JE, 0:KM-1), STAT=STATUS)
      VERIFY_(STATUS)

      ztop = 12000.0
      dz   = ztop/KM
      do k=KS,KE
         height = (ztop - 0.5*dz) - (k)*dz  ! Layer middle height
         do j=JS,JE
            do i=IS,IE
               LONc = LONS(i,j)
               LATc = LATS(i,j)
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


      do k=KS,KE
         do j=JS,JE
            do i=IS,IE
               PKZ(i,j,k) = ( (PE(i,j,k+1)**kappa) - (PE(i,j,k)**kappa) ) /  &
                            ( kappa*(log(PE(i,j,k+1))-log(PE(i,j,k))) )

            enddo
         enddo
      enddo

      PT = PT/PKZ

    elseif (case_id == 4) then ! 3D Rossby-Haurwitz

      do j=JS,JE
         do i=IS,IE
            LONc = LONS(i,j)
            LATc = LATS(i,j)
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
      do k=KS,KE
         do j=JS,JE
            do i=IS,IE
               LONc = LONS(i,j)
               LATc = LATS(i,j)
               pressure = 0.5*(PE(i,j,k)+PE(i,j,k+1))
               call Rossby_Haurwitz(LONc,LATc, pressure, dummy_1, dummy_2, dummy_3, dummy_4, PS(i,j))
               U(i,j,k)  = dummy_1
               V(i,j,k)  = dummy_2
               PT(i,j,k) = dummy_3
               phis(i,j) = dummy_4
            enddo
         enddo
      enddo

      do k=KS,KE
         do j=JS,JE
            do i=IS,IE
               PKZ(i,j,k) = ( (PE(i,j,k+1)**kappa) - (PE(i,j,k)**kappa) ) /  &
                            ( kappa*(log(PE(i,j,k+1))-log(PE(i,j,k))) )

            enddo
         enddo
      enddo
      PT = PT/PKZ

    elseif (case_id == 5) then ! Mountain-Induced Rossby Wave

      do k=KS,KE
         do j=JS,JE
            do i=IS,IE
               LONc = LONS(i,j)
               LATc = LATS(i,j)
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

      do k=KS,KE
         do j=JS,JE
            do i=IS,IE
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
      do k=KS,KE
         height = (ztop - 0.5d0*dz) - (k)*dz  ! Layer middle height
         do j=JS,JE
            do i=IS,IE
               LONc = LONS(i,j)
               LATc = LATS(i,j)
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
         do j=JS,JE
            do i=IS,IE
               LONc = LONS(i,j)
               LATc = LATS(i,j)
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

      do k=KS,KE
         do j=JS,JE
            do i=IS,IE
               PKZ(i,j,k) = ( (PE(i,j,k+1)**kappa) - (PE(i,j,k)**kappa) ) /  &
                            ( kappa*(log(PE(i,j,k+1))-log(PE(i,j,k))) )

            enddo
         enddo
      enddo

      PT = PT/PKZ

    endif ! case_id

!--------------------
! Parse Tracers
!--------------------
   call ESMF_StateGet(IMPORT, 'TRADV' , TRADV_BUNDLE,   RC=STATUS)
   VERIFY_(STATUS)
   call ESMF_FieldBundleGet(TRADV_BUNDLE, fieldCount=NQ, RC=STATUS)
   VERIFY_(STATUS)

   allocate( TRACER(IS:IE, JS:JE, 1:KM), STAT=STATUS)
   VERIFY_(STATUS)
   allocate( state%vars%tracer(NQ), STAT=STATUS)
   VERIFY_(STATUS)

   if (FV3_STANDALONE /= 0) then
      !call ESMF_StateGet(IMPORT, 'TRADV' , TRADV_BUNDLE,   RC=STATUS)
      !VERIFY_(STATUS)

      call ESMF_GridCompGet(gc, grid=esmfGRID, rc=STATUS)
      VERIFY_(STATUS)

      !allocate( TRACER(IS:IE, JS:JE, 1:KM), STAT=STATUS)
      !VERIFY_(STATUS)

      TRACER(:,:,:)  = 0.0
      FIELDNAME = 'Q'
      call addTracer(STATE, TRADV_BUNDLE, TRACER, esmfGRID, FIELDNAME)

    if (case_tracers /= 1234) then

      do n=1,case_tracers
        TRACER(:,:,:)  = 0.0
        write(FIELDNAME, "('Q',i3.3)") n
        call addTracer(STATE, TRADV_BUNDLE, TRACER, esmfGRID, FIELDNAME)
      enddo

    else

!-----------------------------------------------------------------------
!     tracer q1
!-----------------------------------------------------------------------
      TRACER(:,:,:) = 0.0
      do k=KS,KE
         eta = 0.5*( (ak(k-1)+ak(k))/1.e5 + bk(k-1)+bk(k) )
         do j=JS,JE
            do i=IS,IE
               LONc = LONS(i,j)
               LATc = LATS(i,j)
               dummy_1 = tracer_q1_q2(LONc,LATc,eta,rot_ang,r0_6)
               TRACER(i,j,k) = dummy_1
            enddo
         enddo
      enddo
      FIELDNAME = 'Q1'
      call addTracer(STATE, TRADV_BUNDLE, TRACER, esmfGRID, FIELDNAME)

!-----------------------------------------------------------------------
!     tracer q2
!-----------------------------------------------------------------------
      do k=KS,KE
         eta = 0.5*( (ak(k-1)+ak(k))/1.e5 + bk(k-1)+bk(k) )
         do j=JS,JE
            do i=IS,IE
               LONc = LONS(i,j)
               LATc = LATS(i,j)
               dummy_1 = tracer_q1_q2(LONc,LATc,eta,rot_ang,r1_0)
               TRACER(i,j,k) = dummy_1
            enddo
         enddo
      enddo
      FIELDNAME = 'Q2'
      call addTracer(STATE, TRADV_BUNDLE, TRACER, esmfGRID, FIELDNAME)

!-----------------------------------------------------------------------
!     tracer q3
!-----------------------------------------------------------------------
      do k=KS,KE
         eta = 0.5*( (ak(k-1)+ak(k))/1.e5 + bk(k-1)+bk(k) )
         do j=JS,JE
            do i=IS,IE
               LONc = LONS(i,j)
               LATc = LATS(i,j)
               dummy_1 = tracer_q3(LONc,LATc,eta,rot_ang)
               TRACER(i,j,k) = dummy_1
            enddo
         enddo
      enddo
      FIELDNAME = 'Q3'
      call addTracer(STATE, TRADV_BUNDLE, TRACER, esmfGRID, FIELDNAME)

!-----------------------------------------------------------------------
!     tracer q4
!-----------------------------------------------------------------------
      TRACER(:,:,:)  = 1.0_r4
      FIELDNAME = 'Q4'
      call addTracer(STATE, TRADV_BUNDLE, TRACER, esmfGRID, FIELDNAME)
      VERIFY_(STATUS)

!-----------------------------------------------------------------------
!     tracer q5
!-----------------------------------------------------------------------
      if (allocated(Q5)) then
      TRACER(:,:,:)  = Q5(:,:,:)
      FIELDNAME = 'Q5'
      call addTracer(STATE, TRADV_BUNDLE, TRACER, esmfGRID, FIELDNAME)
      VERIFY_(STATUS)
      deallocate( Q5, STAT=STATUS)
      VERIFY_(STATUS)
      endif

!-----------------------------------------------------------------------
!     tracer q6
!-----------------------------------------------------------------------
      if (allocated(Q6)) then
      TRACER(:,:,:)  = Q6(:,:,:)
      FIELDNAME = 'Q6'
      call addTracer(STATE, TRADV_BUNDLE, TRACER, esmfGRID, FIELDNAME)
      VERIFY_(STATUS)
      deallocate( Q6, STAT=STATUS)
      VERIFY_(STATUS)
      endif

    endif

      deallocate( TRACER, STAT=STATUS)
      VERIFY_(STATUS)

    endif
    endif

    DEALLOCATE( PS )

    DYN_COLDSTART=.true.

    call SSI_CopyCoarseToFine(INTERNAL, U, 'U', STATE%f2c_SSI_arr_map,  _RC)
    call SSI_CopyCoarseToFine(INTERNAL, V, 'V', STATE%f2c_SSI_arr_map,  _RC)
    call SSI_CopyCoarseToFine(INTERNAL, PT, 'PT', STATE%f2c_SSI_arr_map,  _RC)
    call SSI_CopyCoarseToFine(INTERNAL, PE, 'PE', STATE%f2c_SSI_arr_map,  _RC)
    call SSI_CopyCoarseToFine(INTERNAL, PKZ, 'PKZ', STATE%f2c_SSI_arr_map,  _RC)
    call SSI_CopyCoarseToFine(IMPORT, phis, 'PHIS', STATE%f2c_SSI_arr_map,  _RC)
   call SSI_CopyCoarseToFine(INTERNAL, LONS, 'LONS', STATE%f2c_SSI_arr_map,  _RC)
   call SSI_CopyCoarseToFine(INTERNAL, LATS, 'LATS', STATE%f2c_SSI_arr_map,  _RC)

    RETURN_(ESMF_SUCCESS)
  end subroutine COLDSTART

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

      data a20_01 / 0.27381905404907E+05,  0.26590539035976E+05,  0.25752394878279E+05,  0.24865429808716E+05, &
                 0.23927536347865E+05,  0.22936541085572E+05,  0.21890203071294E+05,  0.20786212168493E+05, &
                 0.19622187372385E+05,  0.18395675090318E+05,  0.17104147384052E+05,  0.15745000173179E+05, &
                 0.14315551398919E+05,  0.12813039147516E+05,  0.11234619732416E+05,  0.95773657344247E+04, &
                 0.78382639990006E+04,  0.60142135898353E+04,  0.41020236978492E+04,  0.20984115047143E+04, &
                 0.00000000000000E+00 /

      data b20_01 / 0.00000000000000E+00,  0.28901070149364E-01,  0.59510487036309E-01,  0.91902866472543E-01, &
                 0.12615517459290E+00,  0.16234678535331E+00,  0.20055953931639E+00,  0.24087780374962E+00, &
                 0.28338853406205E+00,  0.32818133660555E+00,  0.37534853286773E+00,  0.42498522508382E+00, &
                 0.47718936329560E+00,  0.53206181388604E+00,  0.58970642961892E+00,  0.65023012121324E+00, &
                 0.71374293048299E+00,  0.78035810507338E+00,  0.85019217482527E+00,  0.92336502980036E+00, &
                 0.10000000000000E+01 /

      data a20_0178 / 0.32021324453921E+05,  0.31137565415634E+05,  0.30202026400316E+05,  0.29211673587770E+05, &
                      0.28163295404433E+05,  0.27053492108706E+05,  0.25878664766072E+05,  0.24635003578258E+05, &
                      0.23318475528610E+05,  0.21924811303582E+05,  0.20449491447964E+05,  0.18887731708932E+05, &
                      0.17234467521390E+05,  0.15484337584307E+05,  0.13631666474783E+05,  0.11670446243450E+05, &
                      0.95943169315531E+04,  0.73965459465018E+04,  0.50700062290314E+04,  0.26071531411601E+04, &
                      0.00000000000000E+00 /

      data b20_0178 / 0.00000000000000E+00,  0.27599078219223E-01,  0.56815203138214E-01,  0.87743118501982E-01, &
                      0.12048311914891E+00,  0.15514137625266E+00,  0.19183028162025E+00,  0.23066881216269E+00, &
                      0.27178291572025E+00,  0.31530591949337E+00,  0.36137896240390E+00,  0.41015145278854E+00, &
                      0.46178155290889E+00,  0.51643669184922E+00,  0.57429410846515E+00,  0.63554142614418E+00, &
                      0.70037726124166E+00,  0.76901186716541E+00,  0.84166781619770E+00,  0.91858072126555E+00, &
                      0.10000000000000E+01 /


      data a26 /  219.4067,   489.5209,   988.2418,  1805.2010,  2983.7240,  4462.3340,   &
                 6160.5870,  7851.2430,  7731.2710,  7590.1310,  7424.0860,   &
                 7228.7440,  6998.9330,  6728.5740,  6410.5090,  6036.3220,   &
                 5596.1110,  5078.2250,  4468.9600,  3752.1910,  2908.9490,   &
                  2084.739,   1334.443,    708.499,   252.1360,  0.0, 0.0     /

      data b26 / 0.0, 0.0, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,&
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

       data a137 &
         /1.000000, 2.000365, 3.102241, 4.666084, 6.827977, 9.746966, 13.605424, 18.608931, 24.985718, 32.985710,  &
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

       data b137 &
         /0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,  &
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
120   continue

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
126   continue

      CASE(40)
!--------------------------------------------------
! Pure sigma-coordinate with uniform spacing in "z"
!--------------------------------------------------
!
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
140   continue

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
!
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
160   continue

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
172   continue

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
137   continue

     CASE DEFAULT

        print*, 'Bad KM in FVdycoreCubed_GridComp:set_eta', km

     END SELECT

 end subroutine set_eta
#endif

subroutine addTracer_r8(state, bundle, var, grid, fieldname)
  type (DynState), pointer         :: STATE
  type (ESMF_FieldBundle)          :: BUNDLE
  real(r8), pointer                :: var(:,:,:)
  type (ESMF_Grid)                 :: GRID
  character(len=ESMF_MAXSTR)       :: FIELDNAME

  integer :: nq,rc,status

  character(len=ESMF_MAXSTR)       :: IAm='CoarseFV:addTracer_r8'

  character(len=ESMF_MAXSTR), allocatable :: fieldNames(:)
  integer :: dimen(3)

  call ESMF_FieldBundleGet(BUNDLE, fieldCount=NQ, RC=STATUS)
  VERIFY_(STATUS)

  allocate(fieldNames(NQ))
  call ESMF_FieldBundleGet(bundle, fieldNameList=fieldNames, RC=STATUS)
  VERIFY_(STATUS)
  call SSI_BundleCopyCoarseToFine(bundle, &
       var, fieldname, STATE%f2c_SSI_arr_map, rc=status)
  VERIFY_(STATUS)
  dimen = shape(var)
  do i = 1, NQ
     if (trim(fieldNames(i)) == trim(fieldname)) then
        state%vars%tracer(i)%is_r4 = .false.
        if(.not.associated(state%vars%tracer(i)%content)) then
           allocate(state%vars%tracer(i)%content(dimen(1),dimen(2),dimen(3)), stat=status)
           VERIFY_(status)
        endif
        call SSI_BundleCopyFineToCoarse(bundle, &
             state%vars%tracer(i)%content, fieldname, STATE%f2c_SSI_arr_map, rc=status)
        VERIFY_(STATUS)
        exit
    endif
  enddo

  STATE%GRID%NQ = NQ

  return
end subroutine addTracer_r8

subroutine addTracer_r4(state, bundle, var, grid, fieldname)
  type (DynState), pointer         :: STATE
  type (ESMF_FieldBundle)          :: BUNDLE
  real(r4), pointer                :: var(:,:,:)
  type (ESMF_Grid)                 :: GRID
  character(len=ESMF_MAXSTR)       :: FIELDNAME

  integer :: nq,rc,status

  character(len=ESMF_MAXSTR)       :: IAm='CoarseFV:addTracer_r4'
         
  character(len=ESMF_MAXSTR), allocatable :: fieldNames(:)
  integer :: dimen(3)

  call ESMF_FieldBundleGet(BUNDLE, fieldCount=NQ, RC=STATUS)
  VERIFY_(STATUS)

  allocate(fieldNames(NQ))
  call ESMF_FieldBundleGet(bundle, fieldNameList=fieldNames, RC=STATUS)
  VERIFY_(STATUS)
  call SSI_BundleCopyCoarseToFine(bundle, &
       var, fieldname, STATE%f2c_SSI_arr_map, rc=status)
  VERIFY_(STATUS)
  dimen = shape(var)
  do i = 1, NQ
     if (trim(fieldNames(i)) == trim(fieldname)) then
        state%vars%tracer(i)%is_r4 = .true.
        if(.not.associated(state%vars%tracer(i)%content_r4)) then
           allocate(state%vars%tracer(i)%content_r4(dimen(1),dimen(2),dimen(3)), stat=status)
           VERIFY_(status)
        endif
        call SSI_BundleCopyFineToCoarse(bundle, &
             state%vars%tracer(i)%content_r4, fieldname, STATE%f2c_SSI_arr_map, rc=status)
        VERIFY_(STATUS)
        exit
    endif
  enddo

  STATE%GRID%NQ = NQ

  return
end subroutine addTracer_r4

subroutine freeTracers(state)
  type (DynState) :: STATE

  if (associated(STATE%VARS%tracer)) then
     DEALLOCATE( STATE%VARS%tracer)   ! Comment out to output tracer to checkpoint file
        NULLIFY( STATE%VARS%tracer)
  end if

  return
end subroutine freeTracers

subroutine allocateTracers(state, import, rc)
  type (DynState)  :: state
  type(ESMF_State) :: import
  integer, optional, intent(out) :: rc

  type(ESMF_FieldBundle) :: bundle
  type(ESMF_Field      ) :: field
  type(ESMF_Array      ) :: array
  integer :: status
  integer :: n, nq
  type (ESMF_TypeKind_Flag)  :: kind
  character(len=ESMF_MAXSTR) :: IAm = "CoarseallocateTracers"
  character(len=ESMF_MAXSTR)          :: fieldname
  integer                          :: i1,in,j1,jn,im,jm,km


  i1 = state%grid%is
  in = state%grid%ie
  j1 = state%grid%js
  jn = state%grid%je
  km = state%grid%npz

  BUNDLE = bundleAdv

  call ESMF_FieldBundleGet ( BUNDLE, fieldCount=NQ, RC=STATUS )
  VERIFY_(STATUS)

  allocate(state%vars%tracer(nq), stat=status)
  VERIFY_(status)
  
  do n =  1, nq
     call ESMF_FieldBundleGet(bundle, fieldIndex=n, field=field, rc=status)
     VERIFY_(STATUS)
     call ESMF_FieldGet(FIELD, Array=Array, name=fieldname, RC=STATUS)
     VERIFY_(STATUS)
     call ESMF_ArrayGet(array,typekind=kind,rc=status)
     VERIFY_(STATUS)
     STATE%VARS%TRACER(N)%IS_R4  = (kind == ESMF_TYPEKIND_R4)   ! Is real*4?
     STATE%VARS%TRACER(N)%TNAME = fieldname

     if ( STATE%VARS%TRACER(N)%IS_R4 ) then
        allocate(STATE%VARS%TRACER(N)%content_r4(i1:in,j1:jn,km), stat=status)
        VERIFY_(STATUS)
     else
        allocate(STATE%VARS%TRACER(N)%content(i1:in,j1:jn,km),    stat=status)
        VERIFY_(STATUS)
     end if
  end do

  RETURN_(ESMF_SUCCESS)

end subroutine allocateTracers

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
    real(r8) :: arr_global(grid%npx,grid%ntiles*grid%npy,grid%npz)
    real(r8) :: rng(3,grid%npz)
    real(r8) :: GSUM

    real(kind=ESMF_KIND_R8)     :: locArr(grid%is:grid%ie,grid%js:grid%je)
    real(kind=ESMF_KIND_R8)     :: glbArr(grid%npx,grid%ntiles*grid%npy)

    istrt = grid%is
    iend  = grid%ie
    jstrt = grid%js
    jend  = grid%je
    kstrt = 1
    kend  = grid%npz
    im    = grid%npx
    jm    = grid%npy*grid%ntiles
    km    = grid%npz

  ! call write_parallel('GlobalSUm')
    do k=kstrt,kend
       locArr(:,:) = arr(:,:,k)
       call ArrayGather(locArr, glbArr, grid%grid)
       arr_global(:,:,k) = glbArr
    enddo

    IF (MAPL_AM_I_ROOT()) Then
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

  End Subroutine Write_Profile_R8

  Subroutine Write_Profile_R4(grid, arr, name, delp)
    type (DynGrid),   intent(IN) :: grid
    real(r4),         intent(IN) :: arr(grid%is:grid%ie,grid%js:grid%je,1:grid%npz)
    character(len=*), intent(IN) :: name
    real(r8), optional, intent(IN) :: delp(grid%is:grid%ie,grid%js:grid%je,1:grid%npz)

    integer  :: istrt,iend, jstrt,jend, kstrt,kend
    integer  :: im, jm, km, k
    real(r4) :: arr_global(grid%npx,grid%ntiles*grid%npy,grid%npz)
    real(r4) :: rng(3,grid%npz)
    real(r8) :: gsum_p
    real(r4) :: GSUM
    
    real(kind=ESMF_KIND_R8)     :: locArr(grid%is:grid%ie,grid%js:grid%je)
    real(kind=ESMF_KIND_R8)     :: glbArr(grid%npx,grid%ntiles*grid%npy)
      
    istrt = grid%is
    iend  = grid%ie
    jstrt = grid%js
    jend  = grid%je 
    kstrt = 1
    kend  = grid%npz
    im    = grid%npx
    jm    = grid%npy*grid%ntiles      
    km    = grid%npz
    
    do k=kstrt,kend
       locArr(:,:) = arr(:,:,k)       
       call ArrayGather(locArr, glbArr, grid%grid)
       arr_global(:,:,k) = glbArr
    enddo
    IF (MAPL_AM_I_ROOT()) Then
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
       arr_global(:,:,k) = glbArr
       locArr(:,:) = delp(:,:,k)
       call ArrayGather(locArr, glbArr, grid%grid)
       gsum_p = gsum_p + SUM(SUM(glbArr,DIM=1),DIM=1)
    enddo
    IF (MAPL_AM_I_ROOT()) Then
       GSUM     = SUM(SUM(SUM(arr_global,DIM=1),DIM=1),DIM=1)
       print*,'***********'
       Write(*,"('GlobalSum: ',e21.9)") GSUM/(grid%globalarea*gsum_p)
       print*,'***********'
       print*,' '
    End IF
    endif

  End Subroutine Write_Profile_R4

  function R8_TO_R4(dbl_var)
     real(REAL8), intent(IN) :: dbl_var(:,:)
     real(REAL4)  :: R8_TO_R4(LBOUND(dbl_var,1):UBOUND(dbl_var,1),&
                              LBOUND(dbl_var,2):UBOUND(dbl_var,2))
     integer :: i, j

     real(REAL8), parameter :: eps = 1.e-15_REAL8
     real(REAL8), parameter :: big = 1.e15_REAL8

     do j=LBOUND(dbl_var,2),UBOUND(dbl_var,2)
        do i=LBOUND(dbl_var,1),UBOUND(dbl_var,1)
           R8_TO_R4(i,j) = SIGN(MIN(big,MAX(eps,ABS(dbl_var(i,j)))),dbl_var(i,j))
        enddo
     enddo
  end function

  function R4_TO_R8(sngl_var)
     real(REAL4), intent(IN) :: sngl_var(:,:)
     real(REAL8)  :: R4_TO_R8(LBOUND(sngl_var,1):UBOUND(sngl_var,1),&
                              LBOUND(sngl_var,2):UBOUND(sngl_var,2))
     integer :: i, j

     real(REAL4), parameter :: eps = 1.e-15_REAL4
     real(REAL4), parameter :: big = 1.e15_REAL4

     do j=LBOUND(sngl_var,2),UBOUND(sngl_var,2)
        do i=LBOUND(sngl_var,1),UBOUND(sngl_var,1)
           R4_TO_R8(i,j) = SIGN(MIN(big,MAX(eps,ABS(sngl_var(i,j)))),sngl_var(i,j))
        enddo
     enddo
  end function

end Module CoarseFVdycoreCubed_GridComp
