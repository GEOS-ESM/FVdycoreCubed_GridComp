!  $id: DynCore_GridCompMod.F90,v 1.1.1.1 2007/05/29 12:26:20 atrayanov Exp $

#include "MAPL_Generic.h"

!#define SCALAR_WINDS
!#define INC_WINDS

!-----------------------------------------------------------------------
!              ESMA - Earth System Modeling Applications
!-----------------------------------------------------------------------

   Module FVdycoreCubed_GridComp

!BOP
!
! !MODULE: FVdycoreCubed_GridComp --- Dynamical Core Grid Component
!

! !USES:

   use ESMF                ! ESMF base class
   use MAPL                ! GEOS base class
   use m_set_eta,       only: set_eta

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
                           ADIABATIC, SW_DYNAMICS, AdvCore_Advection
   use m_topo_remap, only: dyn_topo_remap
   use CubeGridPrototype, only: register_grid_and_regridders
! Begin Coarse GC stuff
   use CoarseFVdycoreCubed_GridComp, only : coarse_setvm, &
                                            CoarseSetServices => SetServices, &
                                            DYN_wrap
! End Coarse GC stuff

! !PUBLIC MEMBER FUNCTIONS:

  implicit none
  private

  ! Include the MPI library definitons:
  include 'mpif.h'

  type(ESMF_FieldBundle), save :: bundleAdv
  integer :: NXQ = 0
  logical :: overwrite_Q = .true.

  public  SetServices      ! Register component methods

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

! Begin Coarse GC stuff
  type (ESMF_GridComp)             :: coarseGC
  type (ESMF_State)                :: coarseIM
  type (ESMF_State)                :: coarseEX
  type (ESMF_State)                :: coarseIN
  type (ESMF_VM)                   :: coarseVM
! End Coarse GC stuff

contains

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

   type(ESMF_GridComp), intent(inout) :: gc     ! gridded component
   integer, intent(out), optional     :: rc     ! return code


! !DESCRIPTION: Set services (register) for the FVCAM Dynamical Core
!               Grid Component.
!
!EOP
!----------------------------------------------------------------------

    integer                          :: FV3_STANDALONE
    integer                          :: status
    character(len=ESMF_MAXSTR)       :: IAm
    character(len=ESMF_MAXSTR)       :: COMP_NAME

    type (ESMF_Config)               :: CF
    type (ESMF_VM)                   :: VM

    type (MAPL_MetaComp),      pointer :: MAPL
    character (len=ESMF_MAXSTR)        :: LAYOUT_FILE
! Begin Coarse GC stuff
    integer, allocatable :: gcImg(:) ! holds fine GC image via "transfer" function
! End Coarse GC stuff

! Get the configuration from the component
!-----------------------------------------
    call ESMF_GridCompGet( GC, CONFIG = CF, RC=STATUS )
    call ESMF_GridCompGet( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // "SetServices"


    !call ESMF_VMGetCurrent(VM, rc=STATUS)
    call ESMF_GridCompGet( GC, VM=VM, RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_MemUtilsWrite(VM, trim(IAm)//': Begin', RC=STATUS )
    VERIFY_(STATUS)

!BOS

    call MAPL_AddImportSpec ( gc,                             &
         SHORT_NAME = 'DUDT',                                      &
         LONG_NAME  = 'eastward_wind_tendency',                    &
         UNITS      = 'm s-2',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                             &
         SHORT_NAME = 'DVDT',                                      &
         LONG_NAME  = 'northward_wind_tendency',                   &
         UNITS      = 'm s-2',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                             &
         SHORT_NAME = 'DWDT',                                      &
         LONG_NAME  = 'vertical_velocity_tendency',                &
         UNITS      = 'm s-2',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                             &
         SHORT_NAME = 'DTDT',                                      &
         LONG_NAME  = 'delta-p_weighted_temperature_tendency',     &
         UNITS      = 'Pa K s-1',                                  &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                             &
         SHORT_NAME = 'DQVANA',                                    &
         LONG_NAME  = 'specific_humidity_increment_from_analysis', &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                             &
         SHORT_NAME = 'DQLANA',                                    &
         LONG_NAME  = 'specific_humidity_liquid_increment_from_analysis', &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                             &
         SHORT_NAME = 'DQIANA',                                    &
         LONG_NAME  = 'specific_humidity_ice_increment_from_analysis', &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                             &
         SHORT_NAME = 'DQRANA',                                    &
         LONG_NAME  = 'specific_humidity_rain_increment_from_analysis', &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                             &
         SHORT_NAME = 'DQSANA',                                    &
         LONG_NAME  = 'specific_humidity_snow_increment_from_analysis', &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                             &
         SHORT_NAME = 'DQGANA',                                    &
         LONG_NAME  = 'specific_humidity_graupel_increment_from_analysis', &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                             &
         SHORT_NAME = 'DOXANA',                                    &
         LONG_NAME  = 'ozone_increment_from_analysis',             &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                             &
         SHORT_NAME = 'DPEDT',                                     &
         LONG_NAME  = 'edge_pressure_tendency',                    &
         UNITS      = 'Pa s-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                             &
         SHORT_NAME = 'PHIS',                                      &
         LONG_NAME  = 'surface_geopotential_height',               &
         UNITS      = 'm+2 s-2',                                 &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'VARFLT',                            &
        LONG_NAME          = 'variance_of_filtered_topography',   &
        UNITS              = 'm+2',                               &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART    = MAPL_RestartSkip,                            &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddImportSpec( gc,                              &
        SHORT_NAME = 'TRADV',                                        &
        LONG_NAME  = 'advected_quantities',                        &
        UNITS      = 'unknown',                                    &
        DATATYPE   = MAPL_BundleItem,               &
        RC=STATUS  )
    VERIFY_(STATUS)

! !EXPORT STATE:

    call MAPL_AddExportSpec ( gc,                                         &
         SHORT_NAME       = 'KE',                                         &
         LONG_NAME        = 'vertically_integrated_kinetic_energy',       &
         UNITS            = 'J m-2'  ,                                    &
         DIMS             = MAPL_DimsHorzOnly,                            &
         VLOCATION        = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'TAVE',                                                               &
         LONG_NAME  = 'vertically_averaged_dry_temperature',                                &
         UNITS      = 'K',                                                                  &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'UAVE',                                                               &
         LONG_NAME  = 'vertically_averaged_zonal_wind',                                     &
         UNITS      = 'm sec-1',                                                            &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         FIELD_TYPE = MAPL_VectorField,                                                     &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'KEPHY',                                                              &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_tendency_due_to_physics',       &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                &
       SHORT_NAME   = 'PEPHY',                                                   &
       LONG_NAME    = 'total_potential_energy_tendency_due_to_physics',          &
       UNITS        = 'W m-2',                                                   &
       DIMS         = MAPL_DimsHorzOnly,                                         &
       VLOCATION    = MAPL_VLocationNone,                              RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                &
       SHORT_NAME   = 'TEPHY',                                                   &
       LONG_NAME    = 'mountain_work_tendency_due_to_physics',                   &
       UNITS        = 'W m-2',                                                   &
       DIMS         = MAPL_DimsHorzOnly,                                         &
       VLOCATION    = MAPL_VLocationNone,                              RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                &
       SHORT_NAME         = 'KEANA',                                             &
       LONG_NAME          = 'total_kinetic_energy_tendency_due_to_analysis',     &
       UNITS              = 'W m-2',                                             &
       DIMS               = MAPL_DimsHorzOnly,                                   &
       VLOCATION          = MAPL_VLocationNone,                        RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                &
       SHORT_NAME         = 'PEANA',                                             &
       LONG_NAME          = 'total_potential_energy_tendency_due_to_analysis',   &
       UNITS              = 'W m-2',                                             &
       DIMS               = MAPL_DimsHorzOnly,                                   &
       VLOCATION          = MAPL_VLocationNone,                        RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                &
       SHORT_NAME         = 'TEANA',                                             &
       LONG_NAME          = 'mountain_work_tendency_due_to_analysis',            &
       UNITS              = 'W m-2',                                             &
       DIMS               = MAPL_DimsHorzOnly,                                   &
       VLOCATION          = MAPL_VLocationNone,                        RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                             &
         SHORT_NAME = 'KEHOT',                                                                &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_tendency_due_to_HOT',             &
         UNITS      = 'W m-2',                                                                &
         DIMS       = MAPL_DimsHorzOnly,                                                      &
         VLOCATION  = MAPL_VLocationNone,                                          RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                             &
         SHORT_NAME = 'KEDP',                                                                 &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_tendency_due_to_pressure_change', &
         UNITS      = 'W m-2',                                                                &
         DIMS       = MAPL_DimsHorzOnly,                                                      &
         VLOCATION  = MAPL_VLocationNone,                                          RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                             &
         SHORT_NAME = 'KEADV',                                                                &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_tendency_due_to_dynamics_advection',      &
         UNITS      = 'W m-2',                                                                &
         DIMS       = MAPL_DimsHorzOnly,                                                      &
         VLOCATION  = MAPL_VLocationNone,                                          RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                             &
         SHORT_NAME = 'KEPG',                                                                 &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_tendency_due_to_pressure_gradient',      &
         UNITS      = 'W m-2',                                                                &
         DIMS       = MAPL_DimsHorzOnly,                                                      &
         VLOCATION  = MAPL_VLocationNone,                                          RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'KEDYN',                                                              &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_tendency_due_to_dynamics',      &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'PEDYN',                                                              &
         LONG_NAME  = 'vertically_integrated_potential_energy_tendency_due_to_dynamics',    &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'TEDYN',                                                              &
         LONG_NAME  = 'mountain_work_tendency_due_to_dynamics',                             &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'KECDCOR',                                                            &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_tendency_due_to_cdcore',        &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'PECDCOR',                                                            &
         LONG_NAME  = 'vertically_integrated_potential_energy_tendency_due_to_cdcore',      &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'TECDCOR',                                                            &
         LONG_NAME  = 'mountain_work_tendency_due_to_cdcore',                               &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'QFIXER',                                                             &
         LONG_NAME  = 'vertically_integrated_potential_energy_tendency_due_to_CONSV',       &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'KEREMAP',                                                            &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_tendency_due_to_remap',         &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'PEREMAP',                                                            &
         LONG_NAME  = 'vertically_integrated_potential_energy_tendency_due_to_remap',       &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'TEREMAP',                                                            &
         LONG_NAME  = 'mountain_work_tendency_due_to_remap',                                &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'KEGEN',                                                              &
         LONG_NAME  = 'vertically_integrated_generation_of_kinetic_energy',                 &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'DKERESIN',                                                           &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_residual_from_inertial_terms',  &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'DKERESPG',                                                           &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_residual_from_PG_terms',        &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'DMDTANA',                                                            &
         LONG_NAME  = 'vertically_integrated_mass_tendency_due_to_analysis',                &
         UNITS      = 'kg m-2 s-1',                                                         &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'DOXDTANAINT',                                                        &
         LONG_NAME  = 'vertically_integrated_ozone_tendency_due_to_analysis',               &
         UNITS      = 'kg m-2 s-1',                                                         &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'DQVDTANAINT',                                                        &
         LONG_NAME  = 'vertically_integrated_water_vapor_tendency_due_to_analysis',         &
         UNITS      = 'kg m-2 s-1',                                                         &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'DQLDTANAINT',                                                        &
         LONG_NAME  = 'vertically_integrated_liquid_water_tendency_due_to_analysis',        &
         UNITS      = 'kg m-2 s-1',                                                         &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'DQIDTANAINT',                                                        &
         LONG_NAME  = 'vertically_integrated_ice_water_tendency_due_to_analysis',           &
         UNITS      = 'kg m-2 s-1',                                                         &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'DMDTDYN',                                                            &
         LONG_NAME  = 'vertically_integrated_mass_tendency_due_to_dynamics',                &
         UNITS      = 'kg m-2 s-1',                                                         &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'DOXDTDYNINT',                                                        &
         LONG_NAME  = 'vertically_integrated_ozone_tendency_due_to_dynamics',               &
         UNITS      = 'kg m-2 s-1',                                                         &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'DTHVDTDYNINT',                                                       &
         LONG_NAME  = 'vertically_integrated_THV_tendency_due_to_dynamics',                 &
         UNITS      = 'K kg m-2 s-1',                                                       &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'DTHVDTREMAP',                                                        &
         LONG_NAME  = 'vertically_integrated_THV_tendency_due_to_vertical_remapping',       &
         UNITS      = 'K kg m-2 s-1',                                                       &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'DTHVDTCONSV',                                                        &
         LONG_NAME  = 'vertically_integrated_THV_tendency_due_to_TE_conservation',          &
         UNITS      = 'K kg m-2 s-1',                                                       &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'DTHVDTPHYINT',                                                       &
         LONG_NAME  = 'vertically_integrated_THV_tendency_due_to_physics',                  &
         UNITS      = 'K kg m-2 s-1',                                                       &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'DTHVDTANAINT',                                                       &
         LONG_NAME  = 'vertically_integrated_THV_tendency_due_to_analysis',                 &
         UNITS      = 'K kg m-2 s-1',                                                       &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'DQVDTDYNINT',                                                        &
         LONG_NAME  = 'vertically_integrated_water_vapor_tendency_due_to_dynamics',         &
         UNITS      = 'kg m-2 s-1',                                                         &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'DQLDTDYNINT',                                                        &
         LONG_NAME  = 'vertically_integrated_liquid_water_tendency_due_to_dynamics',        &
         UNITS      = 'kg m-2 s-1',                                                         &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'DQIDTDYNINT',                                                        &
         LONG_NAME  = 'vertically_integrated_ice_water_tendency_due_to_dynamics',           &
         UNITS      = 'kg m-2 s-1',                                                         &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                    &
         SHORT_NAME = 'CONVKE',                                                      &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_convergence',            &
         UNITS      = 'W m-2',                                                       &
         DIMS       = MAPL_DimsHorzOnly,                                             &
         VLOCATION  = MAPL_VLocationCenter,                               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                    &
         SHORT_NAME = 'CONVTHV',                                                     &
         LONG_NAME  = 'vertically_integrated_thetav_convergence',                    &
         UNITS      = 'W m-2',                                                       &
         DIMS       = MAPL_DimsHorzOnly,                                             &
         VLOCATION  = MAPL_VLocationCenter,                               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                    &
         SHORT_NAME = 'CONVCPT',                                                     &
         LONG_NAME  = 'vertically_integrated_enthalpy_convergence',                  &
         UNITS      = 'W m-2',                                                       &
         DIMS       = MAPL_DimsHorzOnly,                                             &
         VLOCATION  = MAPL_VLocationCenter,                               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                    &
         SHORT_NAME = 'CONVPHI',                                                     &
         LONG_NAME  = 'vertically_integrated_geopotential_convergence',              &
         UNITS      = 'W m-2',                                                       &
         DIMS       = MAPL_DimsHorzOnly,                                             &
         VLOCATION  = MAPL_VLocationCenter,                               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'U',                                         &
         LONG_NAME  = 'eastward_wind',                             &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'V',                                         &
         LONG_NAME  = 'northward_wind',                            &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'T',                                         &
         LONG_NAME  = 'air_temperature',                           &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'PL',                                        &
         LONG_NAME  = 'mid_level_pressure',                        &
         UNITS      = 'Pa',                                        &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'ZLE',                                       &
         LONG_NAME  = 'edge_heights',                              &
         UNITS      = 'm',                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'ZLE0',                                      &
         LONG_NAME  = 'edge_heights_above_surface',                &
         UNITS      = 'm',                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'ZL0',                                       &
         LONG_NAME  = 'mid_layer_heights_above_surface',           &
         UNITS      = 'm',                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'ZL',                                        &
         LONG_NAME  = 'mid_layer_heights',                         &
         UNITS      = 'm',                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'S',                                         &
         LONG_NAME  = 'mid_layer_dry_static_energy',               &
         UNITS      = 'm',                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'PLE',                                       &
         LONG_NAME  = 'edge_pressure',                             &
         UNITS      = 'Pa',                                        &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'TH',                                        &
         LONG_NAME  = 'potential_temperature',                     &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'PLK',                                       &
         LONG_NAME  = 'mid-layer_p$^\kappa$',                         &
         UNITS      = 'Pa$^\kappa$',                                  &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'PKE',                                       &
         LONG_NAME  = 'edge_p$^\kappa$',                         &
         UNITS      = 'Pa$^\kappa$',                                  &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'W',                                         &
         LONG_NAME  = 'vertical_velocity',                         &
         UNITS      = 'm s-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'OMEGA',                                     &
         LONG_NAME  = 'vertical_pressure_velocity',                &
         UNITS      = 'Pa s-1',                                  &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'CX',                                        &
         LONG_NAME  = 'eastward_accumulated_courant_number',       &
         UNITS      = '',                                          &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'CY',                                        &
         LONG_NAME  = 'northward_accumulated_courant_number',      &
         UNITS      = '',                                          &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'CU',                                        &
         LONG_NAME  = 'eastward_accumulated_courant_number',       &
         UNITS      = '',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'CV',                                        &
         LONG_NAME  = 'northward_accumulated_courant_number',      &
         UNITS      = '',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'MX',                                        &
         LONG_NAME  = 'pressure_weighted_accumulated_eastward_mass_flux', &
         UNITS      = 'Pa m+2',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'MY',                                        &
         LONG_NAME  = 'pressure_weighted_accumulated_northward_mass_flux', &
         UNITS      = 'Pa m+2',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'MFX',                                       &
         LONG_NAME  = 'pressure_weighted_accumulated_eastward_mass_flux', &
         UNITS      = 'Pa m+2',                                    &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'MFY',                                       &
         LONG_NAME  = 'pressure_weighted_accumulated_northward_mass_flux', &
         UNITS      = 'Pa m+2',                                    &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'MFZ',                                       &
         LONG_NAME  = 'vertical_mass_flux',                        &
         UNITS      = 'kg m-2 s-1',                                &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'PV',                                        &
         LONG_NAME  = 'ertels_isentropic_potential_vorticity',     &
         UNITS      = 'm+2 kg-1 s-1',                            &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'EPV',                                       &
         LONG_NAME  = 'ertels_potential_vorticity',                &
         UNITS      = 'K m+2 kg-1 s-1',                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'Q',                                         &
         LONG_NAME  = 'specific_humidity',                         &
         UNITS      = '1',                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'QC',                                        &
         LONG_NAME  = 'specific_mass_of_condensate',                &
         UNITS      = '1',                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                       &
         SHORT_NAME = 'DUDTSUBZ',                                        &
         LONG_NAME  = 'tendency_of_eastward_wind_due_to_subgrid_dz',      &
         UNITS      = 'm/s/s',                                      &
         DIMS       = MAPL_DimsHorzVert,                                &
         FIELD_TYPE = MAPL_VectorField,                                 &
         VLOCATION  = MAPL_VLocationCenter,                  RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                       &
         SHORT_NAME = 'DVDTSUBZ',                                        &
         LONG_NAME  = 'tendency_of_northward_wind_due_to_subgrid_dz',     &
         UNITS      = 'm/s/s',                                      &
         DIMS       = MAPL_DimsHorzVert,                                &
         FIELD_TYPE = MAPL_VectorField,                                 &
         VLOCATION  = MAPL_VLocationCenter,                  RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                       &
         SHORT_NAME = 'DTDTSUBZ',                                        &
         LONG_NAME  = 'tendency_of_air_temperature_due_to_subgrid_dz',    &
         UNITS      = 'K s-1',                                        &
         DIMS       = MAPL_DimsHorzVert,                                &
         VLOCATION  = MAPL_VLocationCenter,                RC=STATUS    )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                       &
         SHORT_NAME = 'DWDTSUBZ',                                        &
         LONG_NAME  = 'tendency_of_vertical_velocity_due_to_subgrid_dz',     &
         UNITS      = 'm/s/s',                                      &
         DIMS       = MAPL_DimsHorzVert,                                &
         VLOCATION  = MAPL_VLocationCenter,                  RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'DTDT_RAY',                                  &
        LONG_NAME  = 'air_temperature_tendency_due_to_Rayleigh_friction',        &
        UNITS      = 'K s-1',                                  &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'DUDT_RAY',                                  &
        LONG_NAME  = 'tendency_of_eastward_wind_due_to_Rayleigh_friction',       &
        UNITS      = 'm s-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        FIELD_TYPE = MAPL_VectorField,                                 &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'DVDT_RAY',                                  &
        LONG_NAME  = 'tendency_of_northward_wind_due_to_Rayleigh_friction',      &
        UNITS      = 'm s-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        FIELD_TYPE = MAPL_VectorField,                                 &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'DWDT_RAY',                                  &
        LONG_NAME  = 'vertical_velocity_tendency_due_to_Rayleigh_friction',        &
        UNITS      = 'm/s/s',                                  &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                       &
         SHORT_NAME = 'DUDTANA',                                        &
         LONG_NAME  = 'tendency_of_eastward_wind_due_to_analysis',      &
         UNITS      = 'm/s/s',                                      &
         DIMS       = MAPL_DimsHorzVert,                                &
         FIELD_TYPE = MAPL_VectorField,                                 &
         VLOCATION  = MAPL_VLocationCenter,                  RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                       &
         SHORT_NAME = 'DVDTANA',                                        &
         LONG_NAME  = 'tendency_of_northward_wind_due_to_analysis',     &
         UNITS      = 'm/s/s',                                      &
         DIMS       = MAPL_DimsHorzVert,                                &
         FIELD_TYPE = MAPL_VectorField,                                 &
         VLOCATION  = MAPL_VLocationCenter,                  RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                       &
         SHORT_NAME = 'DTDTANA',                                        &
         LONG_NAME  = 'tendency_of_air_temperature_due_to_analysis',    &
         UNITS      = 'K s-1',                                        &
         DIMS       = MAPL_DimsHorzVert,                                &
         VLOCATION  = MAPL_VLocationCenter,                RC=STATUS    )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                       &
         SHORT_NAME = 'DDELPDTANA',                                     &
         LONG_NAME  = 'tendency_of_pressure_thickness_due_to_analysis', &
         UNITS      = 'K s-1',                                        &
         DIMS       = MAPL_DimsHorzVert,                                &
         VLOCATION  = MAPL_VLocationCenter,                RC=STATUS    )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DUDTDYN',                                   &
         LONG_NAME  = 'tendency_of_eastward_wind_due_to_dynamics', &
         UNITS      = 'm/s/s',                                 &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DVDTDYN',                                   &
         LONG_NAME  = 'tendency_of_northward_wind_due_to_dynamics',&
         UNITS      = 'm/s/s',                                 &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                    &
         SHORT_NAME = 'DTDTDYN',                                     &
         LONG_NAME  = 'tendency_of_air_temperature_due_to_dynamics', &
         UNITS      = 'K s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                             &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS    )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                      &
         SHORT_NAME = 'DQVDTDYN',                                      &
         LONG_NAME  = 'tendency_of_specific_humidity_due_to_dynamics', &
         UNITS      = 'kg/kg/s',                                     &
         DIMS       = MAPL_DimsHorzVert,                               &
         VLOCATION  = MAPL_VLocationCenter,                 RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                      &
         SHORT_NAME = 'DQIDTDYN',                                      &
         LONG_NAME  = 'tendency_of_ice_water_due_to_dynamics',         &
         UNITS      = 'kg/kg/s',                                     &
         DIMS       = MAPL_DimsHorzVert,                               &
         VLOCATION  = MAPL_VLocationCenter,                 RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                      &
         SHORT_NAME = 'DQLDTDYN',                                      &
         LONG_NAME  = 'tendency_of_liquid_water_due_to_dynamics',      &
         UNITS      = 'kg/kg/s',                                     &
         DIMS       = MAPL_DimsHorzVert,                               &
         VLOCATION  = MAPL_VLocationCenter,                 RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                      &
         SHORT_NAME = 'DOXDTDYN',                                      &
         LONG_NAME  = 'tendency_of_ozone_due_to_dynamics',             &
         UNITS      = 'mol mol-1 s-1',                                 &
         DIMS       = MAPL_DimsHorzVert,                               &
         VLOCATION  = MAPL_VLocationCenter,                 RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'PREF',                                      &
         LONG_NAME  = 'reference_air_pressure',                    &
         UNITS      = 'Pa',                                        &
         DIMS       = MAPL_DimsVertOnly,                           &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'AK',                                        &
         LONG_NAME  = 'hybrid_sigma_pressure_a',                   &
         UNITS      = '1',                                         &
         DIMS       = MAPL_DimsVertOnly,                           &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'BK',                                        &
         LONG_NAME  = 'hybrid_sigma_pressure_b',                   &
         UNITS      = '1',                                         &
         DIMS       = MAPL_DimsVertOnly,                           &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
       SHORT_NAME         = 'PHIS',                                &
       LONG_NAME          = 'surface_height',                      &
       UNITS              = 'm',                                   &
       DIMS               = MAPL_DimsHorzOnly,                     &
       VLOCATION          = MAPL_VLocationNone,          RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
       SHORT_NAME         = 'PS',                                  &
       LONG_NAME          = 'surface_pressure',                    &
       UNITS              = 'Pa',                                  &
       DIMS               = MAPL_DimsHorzOnly,                     &
       VLOCATION          = MAPL_VLocationNone,          RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
       SHORT_NAME         = 'TA',                                  &
       LONG_NAME          = 'surface_air_temperature',             &
       UNITS              = 'K',                                   &
       DIMS               = MAPL_DimsHorzOnly,                     &
       VLOCATION          = MAPL_VLocationNone,          RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
       SHORT_NAME         = 'QA',                                  &
       LONG_NAME          = 'surface_specific_humidity',           &
       UNITS              = '1',                                   &
       DIMS               = MAPL_DimsHorzOnly,                     &
       VLOCATION          = MAPL_VLocationNone,          RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
       SHORT_NAME         = 'US',                                  &
       LONG_NAME          = 'surface_eastward_wind',               &
       UNITS              = 'm s-1',                               &
       DIMS               = MAPL_DimsHorzOnly,                     &
       FIELD_TYPE         = MAPL_VectorField,                      &
       VLOCATION          = MAPL_VLocationNone,          RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
       SHORT_NAME         = 'VS',                                  &
       LONG_NAME          = 'surface_northward_wind',              &
       UNITS              = 'm s-1',                               &
       DIMS               = MAPL_DimsHorzOnly,                     &
       FIELD_TYPE         = MAPL_VectorField,                      &
       VLOCATION          = MAPL_VLocationNone,          RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
       SHORT_NAME         = 'SPEED',                               &
       LONG_NAME          = 'surface_wind_speed',                  &
       UNITS              = 'm s-1',                               &
       DIMS               = MAPL_DimsHorzOnly,                     &
       VLOCATION          = MAPL_VLocationNone,          RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
       SHORT_NAME         = 'WSPD_10M',                               &
       LONG_NAME          = 'wind_speed_at_10m',                  &
       UNITS              = 'm s-1',                               &
       DIMS               = MAPL_DimsHorzOnly,                     &
       VLOCATION          = MAPL_VLocationNone,          RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
       SHORT_NAME         = 'VVEL_UP_100_1000',                    &
       LONG_NAME          = 'max_vertical_velocity_up_between_100_1000_hPa', &
       UNITS              = 'm s-1',                               &
       DIMS               = MAPL_DimsHorzOnly,                     &
       VLOCATION          = MAPL_VLocationNone,          RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
       SHORT_NAME         = 'VVEL_DN_100_1000',                    &
       LONG_NAME          = 'max_vertical_velocity_down_between_100_1000_hPa', &
       UNITS              = 'm s-1',                               &
       DIMS               = MAPL_DimsHorzOnly,                     &
       VLOCATION          = MAPL_VLocationNone,          RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
       SHORT_NAME         = 'DZ',                                  &
       LONG_NAME          = 'surface_layer_height',                &
       UNITS              = 'm',                                   &
       DIMS               = MAPL_DimsHorzOnly,                     &
       VLOCATION          = MAPL_VLocationNone,          RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
       SHORT_NAME         = 'SLP',                                 &
       LONG_NAME          = 'sea_level_pressure',                  &
       UNITS              = 'Pa',                                  &
       DIMS               = MAPL_DimsHorzOnly,                     &
       VLOCATION          = MAPL_VLocationNone,          RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
       SHORT_NAME         = 'H1000',                               &
       LONG_NAME          = 'height_at_1000_mb',                   &
       UNITS              = 'm',                                   &
       DIMS               = MAPL_DimsHorzOnly,                     &
       VLOCATION          = MAPL_VLocationNone,          RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                        &
       SHORT_NAME         = 'TROPP_EPV',                                                 &
       LONG_NAME          = 'tropopause_pressure_based_on_EPV_estimate',                 &
       UNITS              = 'Pa',                                                        &
       DIMS               = MAPL_DimsHorzOnly,                                           &
       VLOCATION          = MAPL_VLocationNone,                                RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                        &
       SHORT_NAME         = 'TROPP_THERMAL',                                             &
       LONG_NAME          = 'tropopause_pressure_based_on_thermal_estimate',             &
       UNITS              = 'Pa',                                                        &
       DIMS               = MAPL_DimsHorzOnly,                                           &
       VLOCATION          = MAPL_VLocationNone,                                RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                        &
       SHORT_NAME         = 'TROPP_BLENDED',                                             &
       LONG_NAME          = 'tropopause_pressure_based_on_blended_estimate',             &
       UNITS              = 'Pa',                                                        &
       DIMS               = MAPL_DimsHorzOnly,                                           &
       VLOCATION          = MAPL_VLocationNone,                                RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                        &
       SHORT_NAME         = 'TROPK_BLENDED',                                             &
       LONG_NAME          = 'tropopause_index_based_on_blended_estimate',             &
       UNITS              = 'unitless',                                                  &
       DIMS               = MAPL_DimsHorzOnly,                                           &
       VLOCATION          = MAPL_VLocationNone,                                RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                        &
       SHORT_NAME         = 'TROPT',                                                     &
       LONG_NAME          = 'tropopause_temperature_using_blended_TROPP_estimate',       &
       UNITS              = 'K',                                                         &
       DIMS               = MAPL_DimsHorzOnly,                                           &
       VLOCATION          = MAPL_VLocationNone,                                RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                        &
       SHORT_NAME         = 'TROPQ',                                                     &
       LONG_NAME          = 'tropopause_specific_humidity_using_blended_TROPP_estimate', &
       UNITS              = 'kg/kg',                                                     &
       DIMS               = MAPL_DimsHorzOnly,                                           &
       VLOCATION          = MAPL_VLocationNone,                                RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'PLE0',                                      &
         LONG_NAME  = 'pressure_at_layer_edges_before_dynamics',   &
         UNITS      = 'Pa',                                        &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'PLE1',                                      &
         LONG_NAME  = 'pressure_at_layer_edges_after_dynamics',    &
         UNITS      = 'Pa',                                        &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DELP',                                      &
         LONG_NAME  = 'pressure_thickness',                        &
         UNITS      = 'Pa',                                        &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DELPTOP',                                      &
         LONG_NAME  = 'pressure_thickness_at_model_top',                        &
         UNITS      = 'Pa',                                        &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'U_AGRID',                                   &
         LONG_NAME  = 'eastward_wind_on_A-Grid',                   &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'V_AGRID',                                   &
         LONG_NAME  = 'northward_wind_on_A-Grid',                  &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)


    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'U_CGRID',                                   &
         LONG_NAME  = 'eastward_wind_on_C-Grid',                   &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'V_CGRID',                                   &
         LONG_NAME  = 'northward_wind_on_C-Grid',                  &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'U_DGRID',                                   &
         LONG_NAME  = 'eastward_wind_on_native_D-Grid',            &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'V_DGRID',                                   &
         LONG_NAME  = 'northward_wind_on_native_D-Grid',           &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'TV',                                        &
         LONG_NAME  = 'air_virtual_temperature',                   &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'THV',                                       &
         LONG_NAME  = 'scaled_virtual_potential_temperature',      &
         UNITS      = 'K/Pa$^\kappa$',                               &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                       &
         SHORT_NAME = 'DPLEDTDYN',                                       &
         LONG_NAME  = 'tendency_of_edge_pressure_due_to_dynamics', &
         UNITS      = 'Pa s-1',                                       &
         DIMS       = MAPL_DimsHorzVert,                                &
         VLOCATION  = MAPL_VLocationCenter,                RC=STATUS    )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                       &
         SHORT_NAME = 'DDELPDTDYN',                                     &
         LONG_NAME  = 'tendency_of_pressure_thickness_due_to_dynamics', &
         UNITS      = 'Pa s-1',                                       &
         DIMS       = MAPL_DimsHorzVert,                                &
         VLOCATION  = MAPL_VLocationCenter,                RC=STATUS    )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                       &
         SHORT_NAME = 'UKE',                                            &
         LONG_NAME  = 'eastward_flux_of_atmospheric_kinetic_energy',    &
         UNITS      = 'J m-1 s-1',                                      &
         DIMS       = MAPL_DimsHorzOnly,                                &
         FIELD_TYPE = MAPL_VectorField,                                 &
         VLOCATION  = MAPL_VLocationNone,                    RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                       &
         SHORT_NAME = 'VKE',                                            &
         LONG_NAME  = 'northward_flux_of_atmospheric_kinetic_energy',   &
         UNITS      = 'J m-1 s-1',                                      &
         DIMS       = MAPL_DimsHorzOnly,                                &
         FIELD_TYPE = MAPL_VectorField,                                 &
         VLOCATION  = MAPL_VLocationNone,                    RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'UCPT',                                      &
         LONG_NAME  = 'eastward_flux_of_atmospheric_enthalpy',     &
         UNITS      = 'J m-1 s-1',                                 &
         DIMS       = MAPL_DimsHorzOnly,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'VCPT',                                      &
         LONG_NAME  = 'northward_flux_of_atmospheric_enthalpy',    &
         UNITS      = 'J m-1 s-1',                                 &
         DIMS       = MAPL_DimsHorzOnly,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                       &
         SHORT_NAME = 'UPHI',                                           &
         LONG_NAME  = 'eastward_flux_of_atmospheric_potential_energy',  &
         UNITS      = 'J m-1 s-1',                                      &
         DIMS       = MAPL_DimsHorzOnly,                                &
         FIELD_TYPE = MAPL_VectorField,                                 &
         VLOCATION  = MAPL_VLocationNone,                    RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                       &
         SHORT_NAME = 'VPHI',                                           &
         LONG_NAME  = 'northward_flux_of_atmospheric_potential_energy', &
         UNITS      = 'J m-1 s-1',                                      &
         DIMS       = MAPL_DimsHorzOnly,                                &
         FIELD_TYPE = MAPL_VectorField,                                 &
         VLOCATION  = MAPL_VLocationNone,                    RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'UQV',                                       &
         LONG_NAME  = 'eastward_flux_of_atmospheric_water_vapor',  &
         UNITS      = 'kg m-1 s-1',                                &
         DIMS       = MAPL_DimsHorzOnly,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'VQV',                                       &
         LONG_NAME  = 'northward_flux_of_atmospheric_water_vapor', &
         UNITS      = 'kg m-1 s-1',                                &
         DIMS       = MAPL_DimsHorzOnly,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'UQL',                                       &
         LONG_NAME  = 'eastward_flux_of_atmospheric_liquid_water', &
         UNITS      = 'kg m-1 s-1',                                &
         DIMS       = MAPL_DimsHorzOnly,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'VQL',                                       &
         LONG_NAME  = 'northward_flux_of_atmospheric_liquid_water',&
         UNITS      = 'kg m-1 s-1',                                &
         DIMS       = MAPL_DimsHorzOnly,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'UQI',                                       &
         LONG_NAME  = 'eastward_flux_of_atmospheric_ice',          &
         UNITS      = 'kg m-1 s-1',                                &
         DIMS       = MAPL_DimsHorzOnly,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'VQI',                                       &
         LONG_NAME  = 'northward_flux_of_atmospheric_ice',         &
         UNITS      = 'kg m-1 s-1',                                &
         DIMS       = MAPL_DimsHorzOnly,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)


    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DKE',                                       &
         LONG_NAME  = 'tendency_of_atmosphere_kinetic_energy_content_due_to_dynamics',&
         UNITS      = 'W m-2',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DCPT',                                      &
         LONG_NAME  = 'tendency_of_atmosphere_dry_energy_content_due_to_dynamics',&
         UNITS      = 'W m-2',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'DPET',                                      &
         LONG_NAME  = 'tendency_of_atmosphere_topographic_potential_energy_due_to_dynamics',&
         UNITS      = 'W m-2',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'WRKT',                                      &
         LONG_NAME  = 'work_done_by_atmosphere_at_top',            &
         UNITS      = 'W m-2',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DQV',                                       &
         LONG_NAME  = 'tendency_of_atmosphere_water_vapor_content_due_to_dynamics',&
         UNITS      = 'kg m-2 s-1',                                &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DQL',                                       &
         LONG_NAME  = 'tendency_of_atmosphere_liquid_water_content_due_to_dynamics',&
         UNITS      = 'kg m-2 s-1',                                &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DQI',                                       &
         LONG_NAME  = 'tendency_of_atmosphere_ice_content_due_to_dynamics',&
         UNITS      = 'kg m-2 s-1',                                &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'CNV',                                       &
         LONG_NAME  = 'generation_of_atmosphere_kinetic_energy_content',&
         UNITS      = 'W m-2',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

#ifdef SKIP_TRACERS
     do ntracer=1,ntracers
        do nlev=1,nlevs
           write(myTracer, "('Q',i5.5,'_',i3.3)") ntracer-1, plevs(nlev)
           call MAPL_AddExportSpec ( gc,                             &
                SHORT_NAME = TRIM(myTracer),                              &
                LONG_NAME  = TRIM(myTracer),                             &
                UNITS      = '1',                                         &
                DIMS       = MAPL_DimsHorzOnly,                           &
                VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
           VERIFY_(STATUS)
        enddo
        write(myTracer, "('Q',i5.5)") ntracer-1
        call MAPL_AddExportSpec ( gc,                             &
             SHORT_NAME = TRIM(myTracer),                         &
             LONG_NAME  = TRIM(myTracer),                         &
             UNITS      = '1',                                    &
             DIMS       = MAPL_DimsHorzVert,                      &
             VLOCATION  = MAPL_VLocationCenter,               RC=STATUS  )
        VERIFY_(STATUS)
     enddo
#endif

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'UH25',                                      &
         LONG_NAME  = 'updraft_helicity_2_to_5_km',           &
         UNITS      = 'm+2 s-2',                                   &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'UH03',                                      &
         LONG_NAME  = 'updraft_helicity_0_to_3_km',           &
         UNITS      = 'm+2 s-2',                                   &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'SRH01',                                      &
         LONG_NAME  = 'storm_relative_helicity_0_to_1_km',           &
         UNITS      = 'm+2 s-2',                                   &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'SRH03',                                      &
         LONG_NAME  = 'storm_relative_helicity_0_to_3_km',           &
         UNITS      = 'm+2 s-2',                                   &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'SRH25',                                      &
         LONG_NAME  = 'storm_relative_helicity_2_to_5_km',           &
         UNITS      = 'm+2 s-2',                                   &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'VORT',                                      &
         LONG_NAME  = 'vorticity_at_mid_layer_heights',            &
         UNITS      = 's-1',                                       &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'VORT850',                                   &
         LONG_NAME  = 'vorticity_at_850_hPa',                      &
         UNITS      = 's-1',                                       &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'VORT700',                                   &
         LONG_NAME  = 'vorticity_at_700_hPa',                      &
         UNITS      = 's-1',                                       &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'VORT500',                                   &
         LONG_NAME  = 'vorticity_at_500_hPa',                      &
         UNITS      = 's-1',                                       &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'VORT200',                                   &
         LONG_NAME  = 'vorticity_at_200_hPa',                      &
         UNITS      = 's-1',                                       &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DIVG',                                      &
         LONG_NAME  = 'divergence_at_mid_layer_heights',           &
         UNITS      = 's-1',                                       &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'DIVG850',                                   &
         LONG_NAME  = 'divergence_at_850_hPa',                     &
         UNITS      = 's-1',                                       &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'DIVG700',                                   &
         LONG_NAME  = 'divergence_at_700_hPa',                     &
         UNITS      = 's-1',                                       &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'DIVG500',                                   &
         LONG_NAME  = 'divergence_at_500_hPa',                     &
         UNITS      = 's-1',                                       &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'DIVG200',                                   &
         LONG_NAME  = 'divergence_at_200_hPa',                     &
         UNITS      = 's-1',                                       &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'U850',                                      &
         LONG_NAME  = 'eastward_wind_at_850_hPa',                  &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'U700',                                      &
         LONG_NAME  = 'eastward_wind_at_700_hPa',                  &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'U500',                                      &
         LONG_NAME  = 'eastward_wind_at_500_hPa',                  &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'U250',                                      &
         LONG_NAME  = 'eastward_wind_at_250_hPa',                  &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'U200',                                      &
         LONG_NAME  = 'eastward_wind_at_200_hPa',                  &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'UTOP',                                      &
         LONG_NAME  = 'eastward_wind_at_model_top',                  &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'V850',                                      &
         LONG_NAME  = 'northward_wind_at_850_hPa',                 &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationNone,              &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'V700',                                      &
         LONG_NAME  = 'northward_wind_at_700_hPa',                 &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationNone,              &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'V500',                                      &
         LONG_NAME  = 'northward_wind_at_500_hPa',                 &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'V250',                                      &
         LONG_NAME  = 'northward_wind_at_250_hPa',                 &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'V200',                                      &
         LONG_NAME  = 'northward_wind_at_200_hPa',                 &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'VTOP',                                      &
         LONG_NAME  = 'northward_wind_at_model_top',                 &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'T850',                                      &
         LONG_NAME  = 'air_temperature_at_850_hPa',                &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'T700',                                      &
         LONG_NAME  = 'air_temperature_at_700_hPa',                &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'T500',                                      &
         LONG_NAME  = 'air_temperature_at_500_hPa',                &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'T300',                                      &
         LONG_NAME  = 'air_temperature_at_300_hPa',                &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,              &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'T250',                                      &
         LONG_NAME  = 'air_temperature_at_250_hPa',                &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,              &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'TTOP',                                      &
         LONG_NAME  = 'air_temperature_at_model_top',                &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,              &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'Q850',                                      &
         LONG_NAME  = 'specific_humidity_at_850_hPa',              &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'Q500',                                      &
         LONG_NAME  = 'specific_humidity_at_500_hPa',              &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'Q250',                                      &
         LONG_NAME  = 'specific_humidity_at_250_hPa',              &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'Z700',                                      &
         LONG_NAME  = 'geopotential_height_at_700_hPa',            &
         UNITS      = 'm',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,              &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'Z500',                                      &
         LONG_NAME  = 'geopotential_height_at_500_hPa',            &
         UNITS      = 'm',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,              &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'Z300',                                      &
         LONG_NAME  = 'geopotential_height_at_300_hPa',            &
         UNITS      = 'm',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,              &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'H850',                                      &
         LONG_NAME  = 'height_at_850_hPa',                         &
         UNITS      = 'm',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,              &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'H700',                                      &
         LONG_NAME  = 'height_at_700_hPa',                         &
         UNITS      = 'm',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,              &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'H500',                                      &
         LONG_NAME  = 'height_at_500_hPa',                         &
         UNITS      = 'm',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,              &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'H300',                                      &
         LONG_NAME  = 'height_at_300_hPa',                         &
         UNITS      = 'm',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,              &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'H250',                                      &
         LONG_NAME  = 'height_at_250_hPa',                         &
         UNITS      = 'm',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'OMEGA850',                                  &
         LONG_NAME  = 'omega_at_850_hPa',                          &
         UNITS      = 'Pa s-1',                                    &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'OMEGA700',                                  &
         LONG_NAME  = 'omega_at_700_hPa',                          &
         UNITS      = 'Pa s-1',                                    &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'OMEGA500',                                  &
         LONG_NAME  = 'omega_at_500_hPa',                          &
         UNITS      = 'Pa s-1',                                    &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'OMEGA200',                                  &
         LONG_NAME  = 'omega_at_200_hPa',                          &
         UNITS      = 'Pa s-1',                                    &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'OMEGA10',                                   &
         LONG_NAME  = 'omega_at_10_hPa',                           &
         UNITS      = 'Pa s-1',                                    &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'W850',                                      &
         LONG_NAME  = 'w_at_850_hPa',                              &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'W500',                                      &
         LONG_NAME  = 'w_at_500_hPa',                              &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'W200',                                      &
         LONG_NAME  = 'w_at_200_hPa',                              &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'W10',                                       &
         LONG_NAME  = 'w_at_10_hPa',                               &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)


    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'U50M',                                      &
         LONG_NAME  = 'eastward_wind_at_50_meters',                &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'V50M',                                      &
         LONG_NAME  = 'northward_wind_at_50_meters',               &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DXC',                                       &
         LONG_NAME  = 'cgrid_delta_x',                             &
         UNITS      = 'm'  ,                                       &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DYC',                                       &
         LONG_NAME  = 'cgrid_delta_y',                             &
         UNITS      = 'm'  ,                                       &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'AREA',                                      &
         LONG_NAME  = 'agrid_cell_area',                           &
         UNITS      = 'm+2'  ,                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                &
         SHORT_NAME = 'PT',                                        &
         LONG_NAME  = 'scaled_potential_temperature',              &
         UNITS      = 'K Pa$^{-\kappa}$',                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                &
         SHORT_NAME = 'PE',                                        &
         LONG_NAME  = 'air_pressure',                              &
         UNITS      = 'Pa',                                        &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'LONS',                                      &
         LONG_NAME  = 'Center_longitudes',                         &
         UNITS      = 'radians',                                   &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'LATS',                                      &
         LONG_NAME  = 'Center_latitudes',                          &
         UNITS      = 'radians',                                   &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
       SHORT_NAME         = 'DYNTIMER',                            &
       LONG_NAME          = 'timer_for_main_dynamics_run',         &
       UNITS              = 'seconds',                             &
       DIMS               = MAPL_DimsHorzOnly,                     &
       VLOCATION          = MAPL_VLocationNone,          RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
       SHORT_NAME         = 'PID',                                 &
       LONG_NAME          = 'process_id',                          &
       UNITS              = '',                                    &
       DIMS               = MAPL_DimsHorzOnly,                     &
       VLOCATION          = MAPL_VLocationNone,          RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                 &
         SHORT_NAME = 'QV_DYN_IN',                                 &
         LONG_NAME  = 'spec_humidity_at_begin_of_time_step',       &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                 &
         SHORT_NAME = 'T_DYN_IN',                                 &
         LONG_NAME  = 'temperature_at_begin_of_time_step',       &
         UNITS      = 'K',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                 &
         SHORT_NAME = 'U_DYN_IN',                                 &
         LONG_NAME  = 'u_wind_at_begin_of_time_step',       &
         UNITS      = 'm s-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                 &
         SHORT_NAME = 'V_DYN_IN',                                 &
         LONG_NAME  = 'v_wind_at_begin_of_time_step',       &
         UNITS      = 'm s-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                 &
         SHORT_NAME = 'PLE_DYN_IN',                                 &
         LONG_NAME  = 'edge_pressure_at_begin_of_time_step',       &
         UNITS      = 'Pa',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,             RC=STATUS  )
    VERIFY_(STATUS)

! !INTERNAL STATE:

!ALT: technically the first 2 records of "old" style FV restart have
!     6 ints: YYYY MM DD H M S
!     5 ints: I,J,K, KS (num true pressure levels), NQ (num tracers) headers

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'AK',                                        &
         LONG_NAME  = 'hybrid_sigma_pressure_a',                   &
         UNITS      = 'Pa',                                        &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsVertOnly,                           &
         RESTART    = MAPL_RestartRequired,                        &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'BK',                                        &
         LONG_NAME  = 'hybrid_sigma_pressure_b',                   &
         UNITS      = '1',                                         &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsVertOnly,                           &
         RESTART    = MAPL_RestartRequired,                        &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'U',                                         &
         LONG_NAME  = 'eastward_wind',                             &
         UNITS      = 'm s-1',                                     &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         RESTART    = MAPL_RestartRequired,                        &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'V',                                         &
         LONG_NAME  = 'northward_wind',                            &
         UNITS      = 'm s-1',                                     &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         RESTART    = MAPL_RestartRequired,                        &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'PT',                                        &
         LONG_NAME  = 'scaled_potential_temperature',              &
         UNITS      = 'K Pa$^{-\kappa}$',                          &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         RESTART    = MAPL_RestartRequired,                        &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'PE',                                        &
         LONG_NAME  = 'air_pressure',                              &
         UNITS      = 'Pa',                                        &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         RESTART    = MAPL_RestartRequired,                        &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'PKZ',                                       &
         LONG_NAME  = 'pressure_to_kappa',                         &
         UNITS      = 'Pa$^\kappa$',                               &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         RESTART    = MAPL_RestartRequired,                        &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'DZ',                                      &
         LONG_NAME  = 'height_thickness',                          &
         UNITS      = 'm',                                         &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'W',                                         &
         LONG_NAME  = 'vertical_velocity',                         &
         UNITS      = 'm s-1',                                     &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )

!AOO Add LONS and LATS to import to safe as field to be used
!at coarse side where MAPL state is not available
    call MAPL_AddInternalSpec( gc,                                 &
         SHORT_NAME = 'LONS',                                      &
         LONG_NAME  = 'Center_longitudes',                         &
         UNITS      = 'radians',                                   &
         DIMS       = MAPL_DimsHorzOnly,                           &
         RESTART    = MAPL_RestartSkip,                            &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddInternalSpec( gc,                                 &
         SHORT_NAME = 'LATS',                                      &
         LONG_NAME  = 'Center_latitudes',                          &
         UNITS      = 'radians',                                   &
         DIMS       = MAPL_DimsHorzOnly,                           &
         RESTART    = MAPL_RestartSkip,                            &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)
!EOS


! Set the Profiling timers
! ------------------------

    call MAPL_TimerAdd(GC,    name="INITIALIZE"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="RUN"           ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="RUN2"          ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="-DYN_INIT"     ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="--FMS_INIT"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="--FV_INIT"     ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="-DYN_ANA"      ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="-DYN_PROLOGUE" ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="-DYN_CORE"     ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="-DYN_EPILOGUE" ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="--FV_DYNAMICS" ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="--MASS_FIX"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="FINALIZE"      ,RC=STATUS)
    VERIFY_(STATUS)

! Register services for this component
! ------------------------------------

    call MAPL_GridCompSetEntryPoint ( gc, ESMF_METHOD_INITIALIZE,  Initialize, rc=status)
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( gc, ESMF_METHOD_RUN,   Run, rc=status)
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( gc, ESMF_METHOD_RUN,   RunAddIncs, rc=status)
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( gc, ESMF_METHOD_FINALIZE, Finalize, rc=status)
    VERIFY_(STATUS)
 !  call MAPL_GridCompSetEntryPoint ( gc, ESMF_SETREADRESTART, Coldstart, rc=status)
 !  VERIFY_(STATUS)

! Setup FMS/FV3
!--------------
    call MAPL_GetObjectFromGC (GC, MAPL,  RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetResource ( MAPL, LAYOUT_FILE, 'LAYOUT:', default='fvcore_layout.rc', rc=status )
    VERIFY_(STATUS)
    !call DynSetup(GC, state, rc=status)
    !VERIFY_(STATUS)

! Register prototype of cubed sphere grid and associated regridders
!------------------------------------------------------------------
    call register_grid_and_regridders()

! At this point check if FV is standalone and init the grid
!------------------------------------------------------
    call ESMF_ConfigGetAttribute ( CF, FV3_STANDALONE, Label="FV3_STANDALONE:", default=0, RC=STATUS)
    VERIFY_(STATUS)
    if (FV3_STANDALONE /=0) then
        call MAPL_GridCreate(GC, rc=status)
        VERIFY_(STATUS)
        call MAPL_AddExportSpec( gc,                              &
            SHORT_NAME = 'TRADVEX',                                    &
            LONG_NAME  = 'advected_quantities',                        &
            UNITS      = 'unknown',                                    &
            DATATYPE   = MAPL_BundleItem,               &
            RC=STATUS  )
        VERIFY_(STATUS)
    endif

    coarseGC = &
        ESMF_GridCompCreate(NAME="COARSE_DYN", config=CF, &
        RC=STATUS)
    VERIFY_(STATUS)

! Begin Coarse GC stuff
    gcImg = transfer(GC, gcImg)
    call ESMF_AttributeSet(coarseGC, name='GC_IMAGE', valueList=gcImg, rc=status)
    VERIFY_(STATUS)

    call ESMF_GridCompSetVM(coarseGC, userRoutine=coarse_setvm, rc=status)
    VERIFY_(STATUS)


    call ESMF_GridCompSetServices(coarseGC, userRoutine=CoarseSetServices, &
       rc=status)
    VERIFY_(STATUS)
! End Coarse GC stuff

! Generic SetServices
!--------------------

    call MAPL_GenericSetServices( GC, RC=STATUS )
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)

  end subroutine SetServices


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine Initialize ( gc, import, export, clock, rc )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: gc       ! composite gridded component
  type(ESMF_State),    intent(inout) :: import   ! import state
  type(ESMF_State),    intent(inout) :: export   ! export state
  type(ESMF_Clock),    intent(inout) :: clock    ! the clock

  integer, intent(out), OPTIONAL     :: rc       ! Error code:
                                                 ! = 0 all is well
                                                 ! otherwise, error
  type (ESMF_Config)                 :: cf

  type (MAPL_MetaComp),      pointer :: mapl

  character (len=ESMF_MAXSTR)        :: layout_file

  type (ESMF_Field)                  :: field
  real(r4), pointer                      :: pref(:), ak4(:), bk4(:)

  real(r8), pointer                  ::  ak(:)
  real(r8), pointer                  ::  bk(:)
  real(r8), pointer                  ::  ud(:,:,:)
  real(r8), pointer                  ::  vd(:,:,:)
  real(r8), pointer                  ::  pe(:,:,:)
  real(r8), pointer                  ::  pt(:,:,:)
  real(r8), pointer                  ::  pk(:,:,:)

  real(r4), pointer                  :: ple(:,:,:)
  real(r4), pointer                  ::   u(:,:,:)
  real(r4), pointer                  ::   v(:,:,:)
  real(r4), pointer                  ::   t(:,:,:)

! Begin Coarse GC stuff
  real(r4), pointer                  :: LATS(:,:), LONS(:,:)
  real(r4), pointer                  :: LATS_MAPL(:,:), LONS_MAPL(:,:)
! End Coarse GC stuff

  character(len=ESMF_MAXSTR)         :: ReplayMode
  real                               :: DNS_INTERVAL
  type (ESMF_TimeInterval)           :: Intv
  type (ESMF_Alarm)                  :: Alarm
  integer                            :: ColdRestart=0

  integer                            :: status
  character(len=ESMF_MAXSTR)         :: IAm
  character(len=ESMF_MAXSTR)         :: COMP_NAME

  type (ESMF_State)                  :: INTERNAL
  type (DynGrid),  pointer           :: DycoreGrid

  real(r4), pointer                      :: temp2d(:,:)

  integer                            :: ifirst
  integer                            :: ilast
  integer                            :: jfirst
  integer                            :: jlast
  integer                            :: km
  type(ESMF_FieldBundle)             :: tradv, tradvex
  integer                            :: i,numTracers,fv3_standalone
  type(ESMF_Grid)                    :: grid

! Begin
!------

    Iam = "Initialize"
    call ESMF_GridCompGet( GC, name=COMP_NAME, CONFIG=CF, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

! Call Generic Initialize
!------------------------

    call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, CLOCK,  RC=STATUS)
    VERIFY_(STATUS)

! Retrieve the pointer to the state
! ---------------------------------

    call MAPL_GetObjectFromGC (GC, MAPL,  RC=STATUS )
    VERIFY_(STATUS)

! Start the timers
!-----------------

    call MAPL_TimerOn(MAPL,"TOTAL")
    call MAPL_TimerOn(MAPL,"INITIALIZE")

! Check for ColdStart from the configuration
!--------------------------------------
    call MAPL_GetResource ( MAPL, ColdRestart, 'COLDSTART:', default=0, rc=status )
    VERIFY_(STATUS)
    if (ColdRestart /=0 ) then
      call Coldstart_thin( gc, import, export, clock, rc=STATUS )
      VERIFY_(STATUS)
    endif

    call MAPL_Get ( MAPL, INTERNAL_ESMF_STATE=INTERNAL, RC=STATUS )
    VERIFY_(STATUS)

! All fine PETs allocate EXPORT

    call MAPL_GetPointer(export, temp2d, 'DXC', ALLOC=.true., rc=status)
    VERIFY_(STATUS)

    call MAPL_GetPointer(export, temp2d, 'DYC', ALLOC=.true., rc=status)
    VERIFY_(STATUS)

    call MAPL_GetPointer(export, temp2d, 'AREA', ALLOC=.true., rc=status)
    VERIFY_(STATUS)


    call MAPL_GetPointer(EXPORT,PREF,'PREF',ALLOC=.true.,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,AK4 ,'AK'  ,ALLOC=.true.,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,BK4 ,'BK'  ,ALLOC=.true.,RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetPointer(INTERNAL, AK, 'AK', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, BK, 'BK', RC=STATUS)
    VERIFY_(STATUS)

    AK4 = AK
    BK4 = BK
    PREF = AK + BK * P00

    call MAPL_GetPointer(EXPORT,PLE,'PLE',ALLOC=.true.,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,U,  'U',  ALLOC=.true.,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,V,  'V',  ALLOC=.true.,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,T,  'T',  ALLOC=.true.,RC=STATUS)

! Initialize LATS and LONS into INTERNAL state to be retieved on coarse side
! needed for coldstart
    call MAPL_Get ( MAPL, lats = LATS_MAPL, lons = LONS_MAPL, RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, LATS, 'LATS', RC=STATUS)
    VERIFY_(STATUS)
    LATS = LATS_MAPL
    call MAPL_GetPointer(INTERNAL, LONS, 'LONS', RC=STATUS)
    VERIFY_(STATUS)
    LONS = LONS_MAPL

! Begin Coarse GC stuff
    call ESMF_GridCompGet( GC, grid=grid, RC=STATUS )
    VERIFY_(STATUS)
    call ESMF_GridCompSet( coarseGC, grid=grid, RC=STATUS )
    VERIFY_(STATUS)

    call ESMF_GridCompInitialize(coarseGC, importState=IMPORT, &
       exportState=EXPORT, clock=clock, _RC) ! run Initialize
! End Coarse GC stuff

! ======================================================================
!ALT: the next section addresses the problem when export variables have been
!     assigned values during Initialize. To prevent "connected" exports
!     being overwritten by DEFAULT in the Import spec in the other component
!     we label them as being "initailized by restart". A better solution
!     would be to move the computation to phase 2 of Initialize and
!     eliminate this section alltogether
! ======================================================================
    call ESMF_StateGet(EXPORT, 'PREF', FIELD, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_AttributeSet(field, NAME="MAPL_InitStatus", VALUE=MAPL_InitialRestart, RC=STATUS)
    VERIFY_(STATUS)

    call ESMF_StateGet(EXPORT, 'PLE', FIELD, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_AttributeSet(field, NAME="MAPL_InitStatus", VALUE=MAPL_InitialRestart, RC=STATUS)
    VERIFY_(STATUS)

    call ESMF_StateGet(EXPORT, 'U', FIELD, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_AttributeSet(field, NAME="MAPL_InitStatus", VALUE=MAPL_InitialRestart, RC=STATUS)
    VERIFY_(STATUS)

    call ESMF_StateGet(EXPORT, 'V', FIELD, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_AttributeSet(field, NAME="MAPL_InitStatus", VALUE=MAPL_InitialRestart, RC=STATUS)
    VERIFY_(STATUS)

    call ESMF_StateGet(EXPORT, 'T', FIELD, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_AttributeSet(field, NAME="MAPL_InitStatus", VALUE=MAPL_InitialRestart, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetResource ( MAPL, FV3_STANDALONE, Label="FV3_STANDALONE:", default=0, RC=STATUS)
    VERIFY_(STATUS)
    if (FV3_STANDALONE /=0) then
       call ESMF_StateGet(import,'TRADV',tradv,rc=status)
       VERIFY_(STATUS)
       call ESMF_StateGet(export,'TRADVEX',tradvex,rc=status)
       VERIFY_(STATUS)
       call ESMF_FieldBundleGet(tradv,fieldCount=numTracers,rc=status)
       VERIFY_(STATUS)
       do i=1,numTracers
         call ESMF_FieldBundleGet(tradv,fieldIndex=i,field=field,rc=status)
         VERIFY_(status)
         call MAPL_FieldBundleAdd(tradvex,field,rc=status)
         VERIFY_(status)
       enddo
    end if

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

    call MAPL_TimerOff(MAPL,"INITIALIZE")
    call MAPL_TimerOff(MAPL,"TOTAL")

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

  type(ESMF_GridComp), intent(inout) :: gc
  type (ESMF_State),   intent(inout) :: import
  type (ESMF_State),   intent(inout) :: export
  type (ESMF_Clock),   intent(inout) :: clock
  integer, intent(out), optional     :: rc

!EOP
  integer                    :: status
  type (ESMF_FieldBundle)                          :: bundle
  type (ESMF_Field)                                :: field
  type (ESMF_Config)                               :: cf
  type (ESMF_Grid)                                 :: ESMFGRID
  integer :: n

  type (MAPL_MetaComp), pointer :: mapl

  real(kind=4), pointer :: LATS(:,:)
  real(kind=4), pointer :: LONS(:,:)
  real(kind=4), pointer :: temp2d(:,:)

  logical, save                       :: firstime=.true.
  integer, save                       :: nq_saved = 0
  logical                             :: adjustTracers
  type(ESMF_Alarm)                    :: predictorAlarm
  type(ESMF_Grid)                     :: bgrid
  integer                             :: j,pos
  integer                             :: nqt
  logical                             :: tend
  logical                             :: exclude
  character(len=ESMF_MAXSTR)          :: tmpstring
  character(len=ESMF_MAXSTR)          :: fieldname
  character(len=ESMF_MAXSTR)          :: STRING
  character(len=ESMF_MAXSTR)          :: adjustTracerMode
  character(len=ESMF_MAXSTR) :: COMP_NAME
  character(len=ESMF_MAXSTR), allocatable :: xlist(:)
  character(len=ESMF_MAXSTR), allocatable :: biggerlist(:)
  integer, parameter                  :: XLIST_MAX = 60
  logical                             :: isPresent

  character(len=ESMF_MAXSTR) :: IAm

  Iam = "Run"

  call ESMF_GridCompGet( GC, name=COMP_NAME, CONFIG=CF, grid=ESMFGRID, RC=STATUS )
  VERIFY_(STATUS)
  Iam = trim(COMP_NAME) // trim(Iam)

  call ESMF_GridValidate(ESMFGRID,RC=STATUS)
  VERIFY_(STATUS)

  call MAPL_GetObjectFromGC (GC, MAPL,  RC=STATUS )
  VERIFY_(STATUS)

  call MAPL_Get( MAPL, LONS=LONS, LATS=LATS, RC=STATUS )
  VERIFY_(STATUS)

  call MAPL_TimerOn(MAPL,"TOTAL")
  call MAPL_TimerOn(MAPL,"RUN")

  call MAPL_GetPointer(EXPORT, temp2d, 'LONS', RC=STATUS)
  VERIFY_(STATUS)
  if( associated(temp2D) ) temp2d = LONS
  call MAPL_GetPointer(EXPORT, temp2d, 'LATS', RC=STATUS)
  VERIFY_(STATUS)
  if( associated(temp2D) ) temp2d = LATS

  call ESMF_GridCompRun(coarseGC, importState=IMPORT, &
      exportState=EXPORT, clock=clock, PHASE=1, rc=status)
  VERIFY_(STATUS)

  call MAPL_TimerOff(MAPL,"RUN")
  call MAPL_TimerOff(MAPL,"TOTAL")

  RETURN_(ESMF_SUCCESS)

end subroutine RUN

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

    type(ESMF_GridComp), intent(inout) :: gc
    type (ESMF_State),   intent(inout) :: import
    type (ESMF_State),   intent(inout) :: export
    type (ESMF_Clock),   intent(inout)    :: clock
    integer, intent(out), optional     :: rc

!EOP

! !Local Variables:

    type (MAPL_MetaComp), pointer :: genstate

    integer                                          :: status
    character(len=ESMF_MAXSTR) :: IAm

    Iam = "RunAddIncs"

    call MAPL_GetObjectFromGC (GC, GENSTATE,  RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_TimerOn(GENSTATE,"TOTAL")
    call MAPL_TimerOn(GENSTATE,"RUN2")

    call ESMF_GridCompRun(coarseGC, importState=IMPORT, &
        exportState=EXPORT, clock=clock, PHASE=2, rc=status)
    VERIFY_(STATUS)

    call MAPL_TimerOff(GENSTATE,"RUN2")
    call MAPL_TimerOff(GENSTATE,"TOTAL")

    RETURN_(ESMF_SUCCESS)

end subroutine RunAddIncs

!-----------------------------------------------------------------------
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

    type (ESMF_GridComp), intent(inout) :: gc
    type (ESMF_State),    intent(inout) :: import
    type (ESMF_State),    intent(inout) :: export
    type (ESMF_Clock),    intent(inout) :: clock
    integer, optional,    intent(  out) :: rc

!EOP

! Local variables
    type (DYN_wrap) :: wrap
    type (DynState), pointer  :: STATE

    character(len=ESMF_MAXSTR)        :: IAm
    character(len=ESMF_MAXSTR)        :: COMP_NAME
    integer                           :: status

    type (MAPL_MetaComp),     pointer :: MAPL
    type (ESMF_Config)                :: cf


! BEGIN

    Iam = "Finalize"
    call ESMF_GridCompGet( GC, name=COMP_NAME, config=cf, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

! Retrieve the pointer to the state
! ---------------------------------

    call MAPL_GetObjectFromGC (GC, MAPL,  RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_TimerOn(MAPL,"TOTAL")
    call MAPL_TimerOn(MAPL,"FINALIZE")

! Retrieve the pointer to the state
!----------------------------------

    call ESMF_GridCompFinalize(coarseGC, importState=IMPORT, &
        exportState=EXPORT, clock=clock, rc=status)
    VERIFY_(STATUS)

! Call Generic Finalize
!----------------------

    call MAPL_TimerOff(MAPL,"FINALIZE")
    call MAPL_TimerOff(MAPL,"TOTAL")

    call MAPL_GenericFinalize ( GC, IMPORT, EXPORT, CLOCK,  RC=STATUS)
    VERIFY_(STATUS)

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

!BOP

! !IROUTINE: Coldstart_Thin

! !DESCRIPTION:
!   Routine to coldstart from an isothermal state of rest.
!   The temperature can be specified in the config, otherwise
!   it is 300K. The surface pressure is assumed to be 1000 hPa.
!
! !INTERFACE:

subroutine Coldstart_Thin(gc, import, export, clock, rc)

! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: gc
    type(ESMF_State),    intent(inout) :: import
    type(ESMF_State),    intent(inout) :: export
    type (ESMF_Clock),   intent(inout) :: clock
    integer, intent(out), optional     :: rc

!EOP

    character(len=ESMF_MAXSTR)        :: IAm="FV:Coldstart"
    character(len=ESMF_MAXSTR)        :: COMP_NAME
    integer                           :: status

    type (MAPL_MetaComp),     pointer :: MAPL
    type (ESMF_State)                 :: INTERNAL

    type(ESMF_Config)                 :: CF

    integer                     :: case_id
    integer                     :: case_tracers

    integer :: FV3_STANDALONE
    integer :: n

! Tracer Stuff
    type (ESMF_Grid)                 :: esmfGRID
    type (ESMF_FieldBundle)          :: TRADV_BUNDLE
    character(len=ESMF_MAXSTR)       :: FIELDNAME

! Begin

    call ESMF_GridCompGet( GC, name=COMP_NAME, CONFIG=CF, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // trim(Iam)

! Retrieve the pointer to the state
! ---------------------------------

    call MAPL_GetObjectFromGC (GC, MAPL,  RC=STATUS )
    VERIFY_(STATUS)


    call MAPL_Get ( MAPL,                &
           INTERNAL_ESMF_STATE=INTERNAL, &
                               RC=STATUS )
    VERIFY_(STATUS)


! Check if running standalone model
    call ESMF_ConfigGetAttribute ( CF, FV3_STANDALONE, Label="FV3_STANDALONE:", default=0, RC=STATUS)
    VERIFY_(STATUS)

! 3D Baroclinic Test Cases

    call ESMF_ConfigGetAttribute( cf, case_id      , label='CASE_ID:'      , default=0 , rc = STATUS )
    VERIFY_(STATUS)
    call ESMF_ConfigGetAttribute( cf, case_tracers , label='CASE_TRACERS:' , default=1234, rc=STATUS)
    VERIFY_(STATUS)

!--------------------
! Parse Tracers
!--------------------
   if (FV3_STANDALONE /= 0) then
      call ESMF_StateGet(IMPORT, 'TRADV' , TRADV_BUNDLE,   RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_GridCompGet(gc, grid=esmfGRID, rc=STATUS)
      VERIFY_(STATUS)

      FIELDNAME = 'Q'
      call addTracer_thin(TRADV_BUNDLE, esmfGRID, FIELDNAME)

    if (case_tracers /= 1234) then

      do n=1,case_tracers
        write(FIELDNAME, "('Q',i3.3)") n
        call addTracer_thin(TRADV_BUNDLE, esmfGRID, FIELDNAME)
      enddo

    else

!-----------------------------------------------------------------------
!     tracer q1
!-----------------------------------------------------------------------
      FIELDNAME = 'Q1'
      call addTracer_thin(TRADV_BUNDLE, esmfGRID, FIELDNAME)

!-----------------------------------------------------------------------
!     tracer q2
!-----------------------------------------------------------------------
      FIELDNAME = 'Q2'
      call addTracer_thin(TRADV_BUNDLE, esmfGRID, FIELDNAME)

!-----------------------------------------------------------------------
!     tracer q3
!-----------------------------------------------------------------------
      FIELDNAME = 'Q3'
      call addTracer_thin(TRADV_BUNDLE, esmfGRID, FIELDNAME)

!-----------------------------------------------------------------------
!     tracer q4
!-----------------------------------------------------------------------
      FIELDNAME = 'Q4'
      call addTracer_thin(TRADV_BUNDLE, esmfGRID, FIELDNAME)

!-----------------------------------------------------------------------
!     tracer q5
!-----------------------------------------------------------------------
      if (case_id == 3) then
         FIELDNAME = 'Q5'
         call addTracer_thin(TRADV_BUNDLE, esmfGRID, FIELDNAME)

!-----------------------------------------------------------------------
!     tracer q6
!-----------------------------------------------------------------------
         FIELDNAME = 'Q6'
         call addTracer_thin(TRADV_BUNDLE, esmfGRID, FIELDNAME)
      endif

    endif
    endif

    RETURN_(ESMF_SUCCESS)
  end subroutine Coldstart_thin

subroutine addTracer_thin(bundle, grid, fieldname)
  type (ESMF_FieldBundle)          :: BUNDLE
  type (ESMF_Grid)                 :: GRID
  character(len=ESMF_MAXSTR)       :: FIELDNAME

  integer :: nq,rc,status
  type(DynTracers), pointer        :: t(:)

  character(len=ESMF_MAXSTR)       :: IAm='FV:addTracer_thin'

  type (ESMF_Field)                :: field

      call ESMF_FieldBundleGet(BUNDLE, fieldCount=NQ, RC=STATUS)
      VERIFY_(STATUS)

      NQ = NQ + 1

      field = MAPL_FieldCreateEmpty(name=trim(fieldname), grid=GRID, RC=STATUS)
      VERIFY_(STATUS)
      call ESMF_AttributeSet(field,name='VLOCATION',value=MAPL_VLocationCenter,rc=status)
      VERIFY_(STATUS)
      call ESMF_AttributeSet(field,name='DIMS',value=MAPL_DimsHorzVert,rc=status)
      VERIFY_(STATUS)
      call ESMF_AttributeSet(field,name='PRECISION',value=ESMF_KIND_R4,rc=status)
      VERIFY_(STATUS)
      call ESMF_AttributeSet(field,name='HALOWIDTH',value=0,rc=status)
      VERIFY_(STATUS)
      call ESMF_AttributeSet(field,name='DEFAULT_PROVIDED',value=.false.,rc=status)
      VERIFY_(STATUS)
      call MAPL_AllocateCoupling(field, rc=STATUS)
      VERIFY_(STATUS)
      call MAPL_FieldBundleAdd ( bundle, field, rc=STATUS )
      VERIFY_(STATUS)

      !STATE%GRID%NQ = NQ

  return
end subroutine addTracer_thin

end module FVdycoreCubed_GridComp
