#include "MAPL_Generic.h"

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
   use MAPL, only: MAPL_GridCompAddSpec, MAPL_GridCompSetEntryPoint
   use MAPL, only: VERTICAL_STAGGER_CENTER, VERTICAL_STAGGER_NONE, VERTICAL_STAGGER_EDGE
   use MAPL, only: MAPL_GridCompGetResource, MAPL_GridCompGet, MAPL_GridGet
   use MAPL, only: MAPL_GridCompSetGeometry
   use MAPL, only: MAPL_StateGetPointer
   use MAPL, only: MAPL_FieldBundleSameData, MAPL_FieldBundleAdd
   use MAPL, only: MAPL_STATEITEM_SERVICE
   use MAPL, only: MAPL_Verify, MAPL_Return, MAPL_Assert

   use m_set_eta, only: set_eta
   use mpp_mod, only: mpp_pe, mpp_root_pe
   use fv_arrays_mod, only: fv_atmos_type, FVPRC, REAL4, REAL8
   use fms_mod, only: fms_init, set_domain, nullify_domain
   use fv_control_mod, only: fv_init1, fv_init2, fv_end
   use fv_tracer2d_mod, only: offline_tracer_advection
   use fv_mp_mod, only: is, ie, js, je, is_master, tile
   use fv_grid_utils_mod, only: g_sum_r8
   use fv_diagnostics_mod, only: prt_maxmin, prt_minmax
   use FV_StateMod, only: FV_Atm, get_im_world_and_topology
   use FV_StateMod, only: AdvCoreTracers => T_TRACERS
   use FVdycoreCubed_GridComp, only: field_is_cloud_water_species
   use FVdycoreCubed_GridComp, only: get_short_name, is_name_in_list

   use pflogger, only: logger_t => logger

   implicit none
   private

   logical :: FV3_DynCoreIsRunning
   integer :: AdvCore_Advection
   logical, allocatable, save :: grids_on_my_pe(:)
   real(kind=FVPRC) :: dt

   ! Tracer I/O History stuff
   integer, parameter :: ntracers = 38

   !PUBLIC MEMBER FUNCTIONS:
   public SetServices

   !EOP

contains

   !BOP
   !IROUTINE: SetServices - Externally visible registration routine
   !INTERFACE:
   subroutine SetServices(gc, rc)
      !ARGUMENTS:
      type(ESMF_GridComp) :: gc
      integer, intent(out) :: rc

      !DESCRIPTION:
      ! User-supplied setservices routine.
      ! The register routine sets the subroutines to be called
      ! as the init, run, and finalize routines.  Note that those are
      ! private to the module.
      !EOP

      character(len=:), allocatable :: dycore
      character(len=ESMF_MAXSTR) :: my_tracer
      integer :: ndt, itracer
      integer :: status

#include "AdvCore_Import___.h"
#include "AdvCore_Export___.h"
      ! AdvCore provides advection SERVICE
      ! NOTE: SERVICE, irrespective of whether you are a provider or subscriber, adds the bundle
      ! to BOTH the export and import states
      call MAPL_GridCompAddSpec(gc, &
           state_intent=ESMF_STATEINTENT_EXPORT, &
           short_name="TRADV", &
           standard_name="advected_quantities", &
      ! TODO: pchakrab - we shouldn't need dims and vstagger for a bundle
           dims="xyz", &
           vstagger=VERTICAL_STAGGER_NONE, &
           units="unknown", &
           itemtype=MAPL_STATEITEM_SERVICE, _RC)

      ! 3D Tracers, for diagnostics
      do itracer = 1, ntracers
         write(my_tracer, "('TEST_TRACER',i5.5)") itracer - 1
         call MAPL_GridCompAddSpec(gc, &
              state_intent=ESMF_STATEINTENT_EXPORT, &
              short_name=trim(my_tracer), &
              standard_name=trim(my_tracer), &
              units='1', &
              dims="xyz", &
              vstagger=VERTICAL_STAGGER_CENTER, _RC)
      end do

      ! Register methods
      call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_INITIALIZE, Initialize, _RC)
      call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_RUN, Run, phase_name="run", _RC)
      call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_FINALIZE, Finalize, _RC)

      ! TODO: pchakrab - should we get RUN_DT from clock instead?
      call MAPL_GridCompGetResource(gc, "RUN_DT", ndt, default=0, _RC)
      dt = ndt

      ! Check if AdvCore is running without FV3. If not, initialize FMS/FV3
      call MAPL_GridCompGetResource(gc, "DYCORE", dycore, default="", _RC)
      select case (trim(dycore))
      case ("FV3")
         FV3_DynCoreIsRunning = .true.
         AdvCore_Advection = 0
      case ("FV3+ADV")
         FV3_DynCoreIsRunning = .true.
         AdvCore_Advection = 1
      case default
         FV3_DynCoreIsRunning = .false.
         AdvCore_Advection = 1
      end select
      if (.not. FV3_DynCoreIsRunning) call fv_setup(gc, _RC)
      call MAPL_GridCompGetResource(gc, "AdvCore_Advection", AdvCore_Advection, default=AdvCore_Advection, _RC)

      _RETURN(_SUCCESS)
   end subroutine SetServices

   !BOP
   !IROUTINE: Initialize - initialization routine
   !INTERFACE:
   subroutine Initialize(gc, import, export, clock, rc)
      type(ESMF_GridComp) :: gc
      type(ESMF_State) :: import
      type(ESMF_State) :: export
      type(ESMF_Clock) :: clock
      integer, intent(out) :: rc

      !DESCRIPTION:
      ! This initialization routine creates the import and export states,
      ! as well as the internal state, which is attached to the component.
      ! It also determines the distribution (and therefore the grid)
      ! and performs allocations of persistent data,
      !EOP

      !BOC
      real, pointer :: area(:, :)
      integer :: is, ie, js, je, status

      ! pchakrab - maybe use the FV_DynCoreIsRunning flag instead
      if (.not. FV3_DynCoreIsRunning) then
         call MAPL_GridCompSetGeometry(gc, _RC)
      end if
      ! pchakrab - below was the original logic for grid creation
      ! grid_created = .false.
      ! call ESMF_GridCompGet(gc, grid=grid, rc=status)
      ! if (status == ESMF_SUCCESS) then
      !    call ESMF_GridValidate(grid, rc=status)
      !    if (status == ESMF_SUCCESS) grid_created = .true.
      ! end if
      ! if (.not. grid_created) call MAPL_GridCompSetGeometry(gc, _RC)

      ! Compute Grid-Cell Area
      if (.not. FV3_DynCoreIsRunning) then
         is = FV_Atm(1)%bd%isc
         ie = FV_Atm(1)%bd%iec
         js = FV_Atm(1)%bd%jsc
         je = FV_Atm(1)%bd%jec
         call MAPL_StateGetPointer(export, area, "AREA", _RC)
         if (associated(area)) then
            area = FV_Atm(1)%gridstruct%area(is:ie, js:je)
         end if
      end if

      _RETURN(_SUCCESS)
      _UNUSED_DUMMY(import)
      _UNUSED_DUMMY(clock)
   end subroutine Initialize

   !BOP
   !IROUTINE: Run - run routine
   !INTERFACE:
   subroutine Run(gc, import, export, clock, rc)
      type(ESMF_GridComp) :: gc
      type(ESMF_State) :: import
      type(ESMF_State) :: export
      type(ESMF_Clock) :: clock
      integer, intent(out) :: rc

      !DESCRIPTION:
      ! The Run method advanced the advection one long time step, as
      ! specified in the configuration.  This may be broken down int a
      ! number of internal, small steps, also configurable.
      !EOP

      type(ESMF_Grid) :: esmfgrid
      type(ESMF_HConfig) :: hconfig
      type(ESMF_FieldBundle) :: tradv, tradv_import
      type(ESMF_Field) :: field
      type(ESMF_Array) :: array
      type(ESMF_TypeKind_Flag) :: tk
      type(ESMF_FieldBundle), save :: bundle_adv
      type(AdvCoreTracers), allocatable :: adv_tracers(:)
#include "AdvCore_DeclarePointer___.h"
      real(kind=FVPRC), allocatable, dimension(:, :, :) :: xPLE1
      real(kind=FVPRC), allocatable, dimension(:, :, :, :) :: tracers
      real(kind=REAL8), allocatable :: tmass0(:), tmass1(:)
      real, pointer :: temp3D(:, :, :)
      real(kind=REAL4), pointer :: tracer_r4(:, :, :)
      real(kind=REAL8), pointer :: tracer_r8(:, :, :)
      real(kind=FVPRC), allocatable :: DEBUG_ARRAY(:, :, :)
      real(kind=FVPRC), parameter :: fac1 = 1.0
      real :: tmp0, tmp1, tmp2
      logical :: same_tradv_data, rpt_mass, DEBUG_ADV, adjust_tracers
      logical, save :: first_run = .true.
      character(len=ESMF_MAXSTR) :: field_name, mytracer
      integer :: n, im, jm, lm, nq, QSPLIT, status
      class(logger_t), pointer :: logger

      _RETURN_UNLESS(AdvCore_Advection > 0)

      call MAPL_GridCompGet(gc, hconfig=hconfig, logger=logger, _RC)

#include "AdvCore_GetPointer___.h"
      xPLE1 = PLE1

      ! SERVICE adds the bundle containing tracers to be advected to both import and export states
      ! Here we copy the tracer data from the import bundle to the export
      call ESMF_StateGet(import, "TRADV", tradv_import, _RC)
      call ESMF_StateGet(export, "TRADV", tradv, _RC)
      ! Instead of copying, we ensure that bundle_imp and bundle point to the same data in the
      ! contained fields. This is important to check because we don't have a coupling mechanism yet
      same_tradv_data = MAPL_FieldBundleSameData(tradv_import, tradv, _RC)
      _ASSERT(same_tradv_data, "TRADV bundles in import and export do not point to the same data")

      adjust_tracers = should_adjust_tracers(gc, clock, _RC)
      if (adjust_tracers) then
         if (first_run) then
            first_run = .false.
            bundle_adv = get_adjusted_tracer_bundle(tradv, hconfig, _RC) ! save'd
         end if
         tradv = bundle_adv
      end if

      call ESMF_FieldBundleGet(tradv, fieldCount=nq, _RC)
      _RETURN_UNLESS(nq > 0)

      call logger%info("Advcore is Advecting the following %i tracers", nq)
      do n = 1, nq
         call ESMF_FieldBundleGet(tradv, fieldIndex=n, field=field, _RC)
         field_name = get_short_name(field, _RC)
         call logger%info(field_name)
      end do

      ! We allocate a list of tracers big enough to hold all items in the bundle
      call MAPL_GridCompGet(gc, grid=esmfgrid, num_levels=lm, _RC)
      call MAPL_GridGet(esmfgrid, im=im, jm=jm, _RC)
      allocate(tracers(im, jm, lm, nq), adv_tracers(nq), _STAT)

      ! Go through the bundle copying the fields into the tracer list
      ! This is essentially PULL_Q in DynCore
      do n = 1, nq
         call ESMF_FieldBundleGet(tradv, fieldIndex=n, field=field, _RC)
         call ESMF_FieldGet(field, array=array, _RC)
         call ESMF_ArrayGet(array, typekind=tk, _RC)
         adv_tracers(n)%is_r4 = (tk == ESMF_TYPEKIND_R4)
         adv_tracers(n)%tname = get_short_name(field, _RC)
         if (adv_tracers(n)%is_r4) then
            call ESMF_ArrayGet(array, farrayptr=tracer_r4, _RC)
            adv_tracers(n)%content_r4 => tracer_r4
            tracers(:, :, :, n) = tracer_r4
         else
            call ESMF_ArrayGet(array, farrayptr=tracer_r8, _RC)
            adv_tracers(n)%content => tracer_r8
            tracers(:, :, :, n) = tracer_r8
         end if
      end do

      ! Tracer Mass before advection
      call MAPL_GridCompGetResource(gc, "ADV_CORE_REPORT_TRACER_MASS", rpt_mass, default=.false., _RC)
      if (rpt_mass) then
         allocate(tmass0(nq), _STAT)
         call global_integral(tmass0, tracers, real(PLE0, kind=FVPRC), im, jm, lm, nq)
      end if

      ! Run FV3 advection
      call MAPL_GridCompGetResource(gc, "ADV_QSPLIT", QSPLIT, default=0, _RC)
      call offline_tracer_advection(tracers, &
           real(PLE0, kind=FVPRC), xPLE1, &
           real(MFX, kind=FVPRC), real(MFY, kind=FVPRC), &
           real(CX, kind=FVPRC), real(CY, kind=FVPRC), &
           FV_Atm(1)%gridstruct, FV_Atm(1)%flagstruct, FV_Atm(1)%bd, FV_Atm(1)%domain, &
           FV_Atm(1)%npx, FV_Atm(1)%npy, FV_Atm(1)%npz, nq, dt, QSPLIT)

      ! Tracer Mass after advection
      if (rpt_mass) then
         allocate(tmass1(nq), _STAT)
         call global_integral(tmass1, tracers, xPLE1, im, jm, lm, nq)
         ! Conserve Specific Mass of Constituents Keeping Mixing_Ratio Constant WRT_Dry_Air
         do n = 1, nq
            if (tmass1(n) > 0.0) then
               tmp0 = tmass1(n) - tmass0(n)
               tmp1 = tmp0 / tmass1(n)
               if (abs(tmp1) >= epsilon(1.0_REAL4)) then
                  tmp2 = tmp0 / tmass0(n)
                  call logger%info(trim(adv_tracers(n)%tname) // ": %f", tmp2)
                  ! tracers(:,:,:,n) = tracers(:,:,:,n) * tmass0(n)/tmass1(n)
               end if
            end if
         end do
      end if

      !--> Fill Export States
      !--> This section is used for diagnostics only.
      !--> It has no effect on CTM experiments.
      do n = 1, nq
         if (n <= min(ntracers, nq)) then
            write(mytracer, "('TEST_TRACER',i5.5)") n - 1
            call MAPL_StateGetPointer(export, temp3D, trim(mytracer), _RC)
            if (associated(temp3D)) temp3D = tracers(:, :, :, n)
         end if
      end do

      ! Clean negative tracers
      do n = 1, nq
         where (tracers(:, :, :, n) < tiny(0.0))
            tracers(:, :, :, n) = 0.0
         end where
      end do

      ! Save tracer data for subsequent runs
      do n = 1, nq
         if (adv_tracers(n)%is_r4) then
            adv_tracers(n)%content_r4 = tracers(:, :, :, n)
         else
            adv_tracers(n)%content = tracers(:, :, :, n)
         end if
      end do

      call MAPL_GridCompGetResource(gc, "DEBUG_ADV", DEBUG_ADV, default=.false., _RC)
      if (DEBUG_ADV) then
         prt_minmax = DEBUG_ADV
         if (mpp_pe() == 0) print*, ''
         if (mpp_pe() == 0) print*, '-------------- FV3 Tracer Debug After ADV --------------'
         allocate(DEBUG_ARRAY(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec, FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec, FV_Atm(1)%npz))
         do n = 1, nq
            DEBUG_ARRAY = tracers(:, :, :, n)
            call prt_maxmin( &
                 trim(adv_tracers(n)%tname), DEBUG_ARRAY, &
                 FV_Atm(1)%bd%isc, FV_Atm(1)%bd%iec, FV_Atm(1)%bd%jsc, FV_Atm(1)%bd%jec, &
                 0, FV_Atm(1)%npz, fac1)
         end do
         if (mpp_pe() == 0) print*, '-------------- FV3 Tracer Debug After ADV --------------'
         if (mpp_pe() == 0) print*, ''
         prt_minmax = .false.
      end if

      _RETURN(_SUCCESS)
   end subroutine Run

   !BOP
   !IROUTINE:  Finalize - user supplied finalize routine
   !INTERFACE:
   subroutine Finalize(gc, import, export, clock, rc)
      type(ESMF_GridComp) :: gc
      type(ESMF_State) :: import
      type(ESMF_State) :: export
      type(ESMF_Clock) :: clock
      integer, intent(out) :: rc

      !DESCRIPTION:
      !    Finalize merely destroys the FVadv object that was created in Initialize
      !    and releases the space for the persistent data .
      !EOP

      ! Clean up FV if AdvCore is running without FV3_DynCoreIsRunning
      if (.not. FV3_DynCoreIsRunning) then
         call fv_end(FV_Atm, grids_on_my_pe, .false.)
      end if

      _RETURN(_SUCCESS)
      _UNUSED_DUMMY(gc)
      _UNUSED_DUMMY(import)
      _UNUSED_DUMMY(export)
      _UNUSED_DUMMY(clock)
   end subroutine Finalize

   subroutine global_integral(QG, Q, PLE, im, jm, KM, nq)
      real(kind=REAL8), intent(out) :: QG(nq)
      real(kind=FVPRC), intent(in) :: Q(im, jm, KM, nq)
      real(kind=FVPRC), intent(in) :: PLE(im, jm, KM + 1)
      integer, intent(in) :: im, jm, KM, nq
      ! Locals
      integer :: k, n
      real(kind=REAL8), allocatable :: dp(:, :, :)
      real(kind=REAL8), allocatable :: qsum1(:, :)
      real(kind=REAL8) :: mass

      allocate(dp(im, jm, KM))
      allocate(qsum1(im, jm))

      ! Compute Pressure Thickness
      do k = 1, KM
         dp(:, :, k) = PLE(:, :, k + 1) - PLE(:, :, k)
      end do

      ! Compute Global Mass
      qsum1(:, :) = 0.d0
      do k = 1, KM
         qsum1(:, :) = qsum1(:, :) + dp(:, :, k)
      end do
      mass = g_sum_r8(FV_Atm(1)%domain, qsum1, is, ie, js, je, FV_Atm(1)%ng, FV_Atm(1)%gridstruct%area_64, 1)

      ! Loop over Tracers
      do n = 1, nq
         qsum1(:, :) = 0.d0
         do k = 1, KM
            qsum1(:, :) = qsum1(:, :) + Q(:, :, k, n) * dp(:, :, k)
         end do
         QG(n) = g_sum_r8(FV_Atm(1)%domain, qsum1, is, ie, js, je, FV_Atm(1)%ng, FV_Atm(1)%gridstruct%area_64, 1)
         if (mass > 0.0) QG(n) = QG(n) / mass
      end do

      deallocate(dp)
      deallocate(qsum1)
   end subroutine global_integral

   function should_adjust_tracers(gc, clock, rc) result(adjust_tracers)
      type(ESMF_GridComp), intent(inout) :: gc
      type(ESMF_Clock), intent(in) :: clock
      integer, intent(out) :: rc
      logical :: adjust_tracers

      class(logger_t), pointer :: logger
      character(len=:), allocatable :: adjust_tracer_mode
      type(ESMF_Alarm) :: predictor_alarm
      integer :: status

      adjust_tracers = .false.

      call MAPL_GridCompGetResource(gc, "EXCLUDE_ADVECTION_TRACERS", adjust_tracer_mode, default="ALWAYS", _RC)
      select case (trim(adjust_tracer_mode))
      case ("ALWAYS")
         adjust_tracers = .true.
      case ("PREDICTOR")
         call ESMF_ClockGetAlarm(clock, alarmName="PredictorAlarm", alarm=predictor_alarm, rc=status)
         if (status == ESMF_SUCCESS) then
            if (ESMF_AlarmIsRinging(predictor_alarm)) then
               adjust_tracers = .true.
            end if
         end if
      case ("NO")
         adjust_tracers = .false. ! proceed without warning
      case default
         call MAPL_GridCompGet(gc, logger=logger, _RC)
         call logger%info("Invalid value specified for EXCLUDE_ADVECTION_TRACERS, ignored")
         adjust_tracers = .false.
      end select

      _RETURN(_SUCCESS)
   end function should_adjust_tracers

   function get_adjusted_tracer_bundle(orig_bundle, hconfig, rc) result(adjusted_bundle)
      type(ESMF_FieldBundle), intent(in) :: orig_bundle
      type(ESMF_HConfig), intent(in) :: hconfig
      integer, intent(out) :: rc
      type(ESMF_FieldBundle) :: adjusted_bundle

      type(ESMF_Grid) :: grid
      type(ESMF_Field) :: field
      character(len=:), allocatable :: field_name
      character(len=ESMF_MAXSTR), allocatable :: xlist(:), biggerlist(:)
      integer :: i, nqt, num_excluded, status

      call ESMF_FieldBundleGet(orig_bundle, grid=grid, fieldCount=nqt, _RC)
      adjusted_bundle = ESMF_FieldBundleCreate(name="xTRADV", _RC) ! saved between runs
      call ESMF_FieldBundleSet(adjusted_bundle, grid=grid, _RC)
      xlist = ESMF_HConfigAsStringSeq(hconfig, keyString="EXCLUDE_ADVECTION_TRACERS_LIST", stringLen=ESMF_MAXSTR, _RC)
      num_excluded = 0
      if (allocated(xlist)) num_excluded = size(xlist)
      ! Exclude tracers that are to be advected by DynCore (cloud/water species)
      do i = 1, nqt
         call ESMF_FieldBundleGet(orig_bundle, fieldIndex=i, field=field, _RC)
         field_name = get_short_name(field, _RC)
         if (field_is_cloud_water_species(field_name)) then
            num_excluded = num_excluded + 1
            if (num_excluded > size(xlist)) then
               allocate(biggerlist(2*num_excluded), _STAT)
               biggerlist(1:num_excluded-1) = xlist
               call move_alloc(from=biggerlist, to=xlist)
            end if
            xlist(num_excluded) = trim(field_name)
         end if
      end do
      ! Add non-excluded tracers to adjusted bundle
      do i = 1, nqt
         call ESMF_FieldBundleGet(orig_bundle, fieldIndex=i, field=field, _RC)
         field_name = get_short_name(field, _RC)
         if (.not. is_name_in_list(field_name, xlist)) then
            call MAPL_FieldBundleAdd(adjusted_bundle, [field], _RC)
         end if
      end do

      _RETURN(_SUCCESS)
   end function get_adjusted_tracer_bundle

   subroutine fv_setup(gc, rc)
      type(ESMF_GridComp), intent(inout) :: gc
      integer, intent(out) :: rc

      type(ESMF_VM) :: vm
      type(ESMF_HConfig) :: hconfig
      integer :: p_split = 1
      integer :: comm, im_world, topology(2), num_levels, status

      call MAPL_GridCompGet(gc, hconfig=hconfig, num_levels=num_levels, _RC)

      call ESMF_VMGetCurrent(vm, _RC)
      call ESMF_VMGet(vm, mpiCommunicator=comm, _RC)
      call fms_init(comm)

      call fv_init1(FV_Atm, dt, grids_on_my_pe, p_split)
      call get_im_world_and_topology(hconfig, im_world, topology, _RC)
      associate(layout => FV_Atm(1)%layout, flags => FV_Atm(1)%flagstruct)
         ! Domain decomposition
         layout = topology
         if (flags%grid_type == 4) then
            layout(2) = layout(2) * 6
         end if
         ! Grid dimensions
         flags%npx = im_world
         flags%npy = im_world * 6
         flags%npz = num_levels
         ! TODO: pchakrab - check this! npy is always set to 6*npx!
         if (flags%npy == 6 * flags%npx) then
            flags%npy = flags%npx + 1
            flags%npx = flags%npx + 1
            flags%ntiles = 6
         else
            flags%npy = flags%npy + 1
            flags%npx = flags%npx + 1
            flags%ntiles = 1
         end if
      end associate
      call fv_init2(FV_Atm, dt, grids_on_my_pe, p_split)
   end subroutine fv_setup

end module AdvCore_GridCompMod

subroutine AdvCore_SetServices(gc, rc)
   use ESMF, only: ESMF_GridComp
   use AdvCore_GridCompMod, only: mySetServices => SetServices
   type(ESMF_GridComp) :: gc
   integer, intent(out) :: rc
   call mySetServices(gc, rc=rc)
end subroutine AdvCore_SetServices
