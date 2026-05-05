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
   use FV_StateMod, only: AdvCoreTracers => T_TRACERS
   use FV_StateMod, only: FV_Atm
   use FV_StateMod, only: get_im_world_and_topology

   use pflogger, only: logger_t => logger

   implicit none
   private

   logical :: FV3_DynCoreIsRunning = .false.
   integer :: AdvCore_Advection = 1
   logical, allocatable, save :: grids_on_my_pe(:)
   real(kind=FVPRC) :: dt

   ! Tracer I/O History stuff
   integer, parameter :: ntracers = 38
   logical, save :: first_run = .true.

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

      type(ESMF_VM) :: vm
      type(ESMF_HConfig) :: hconfig
      character(len=:), allocatable :: dycore, my_tracer
      integer :: comm, ndt, itracer, im_world, topology(2), num_levels
      integer :: p_split = 1
      integer :: status

#include "AdvCore_Import___.h"
      call MAPL_GridCompAddSpec(gc, &
           state_intent=ESMF_STATEINTENT_EXPORT, &
           short_name="TRADV", &
           standard_name="advected_quantities", &
      ! TODO: pchakrab - we shouldn't need dims and vstagger for a bundle
           dims="xyz", &
           vstagger=VERTICAL_STAGGER_NONE, &
           units="unknown", &
           itemtype=MAPL_STATEITEM_SERVICE, _RC)

#include "AdvCore_Export___.h"

      ! 3D Tracers
      ! TODO: pchakrab - do we need these ntracers
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

      ! Register methods with MAPL
      call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_INITIALIZE, Initialize, _RC)
      call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_RUN, Run, _RC)
      call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_FINALIZE, Finalize, _RC)

      ! Check if AdvCore is running without FV3_DynCoreIsRunning, if yes then setup the MAPL Grid
      call MAPL_GridCompGetResource(gc, "DYCORE", dycore, default="", _RC)
      ! TODO: pchakrab - use select case here
      if (trim(dycore) == "FV3") then
         FV3_DynCoreIsRunning = .true.
         AdvCore_Advection = 0
      end if
      if (trim(dycore) == "FV3+ADV") then
         FV3_DynCoreIsRunning = .true.
         AdvCore_Advection = 1
      end if
      call MAPL_GridCompGetResource(gc, "AdvCore_Advection", AdvCore_Advection, default=AdvCore_Advection, _RC)


      ! Start up FMS/MPP
      call ESMF_VMGet(vm, mpiCommunicator=comm, _RC)
      call fms_init(comm)

      call MAPL_GridCompGetResource(gc, "RUN_DT", ndt, default=0, _RC)
      dt = ndt

      if (.not. FV3_DynCoreIsRunning) then
         ! Make sure FV3 is setup
         call fv_init1(FV_Atm, dt, grids_on_my_pe, p_split)

         call MAPL_GridCompGet(gc, hconfig=hconfig, num_levels=num_levels, _RC)
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

         ! TODO: pchakrab - shouldn't we setup the geomtry here?

         call fv_init2(FV_Atm, dt, grids_on_my_pe, p_split)
      end if

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
      type(ESMF_Grid) :: grid
      real, pointer :: area(:, :)
      logical :: grid_created
      integer :: is, ie, js, je, status

      grid_created = .false.
      call ESMF_GridCompGet(gc, grid=grid, rc=status)
      if (status == ESMF_SUCCESS) then
         call ESMF_GridValidate(grid, rc=status)
         if (status == ESMF_SUCCESS) grid_created = .true.
      end if
      if (.not. grid_created) call MAPL_GridCompSetGeometry(gc, _RC)

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

      integer :: status
      type(ESMF_Grid) :: esmfgrid, bgrid
      type(ESMF_HConfig) :: hconfig
      type(ESMF_FieldBundle) :: tradv, tradv_import
      type(ESMF_Field) :: field
      type(ESMF_Array) :: array
      type(ESMF_TypeKind_Flag) :: kind
      type(ESMF_FieldBundle), save :: bundle_adv
      type(ESMF_Alarm) :: predictor_alarm

      ! Imports
      ! TODO: pchakrab - replace with AdvCore_DeclarePointer___.h once the typekind bug is fixed
      real(kind=REAL8), pointer, dimension(:, :, :) :: area
      real(kind=REAL8), pointer, dimension(:, :, :) :: CX
      real(kind=REAL8), pointer, dimension(:, :, :) :: CY
      real(kind=REAL8), pointer, dimension(:, :, :) :: MFX
      real(kind=REAL8), pointer, dimension(:, :, :) :: MFY
      real(kind=REAL8), pointer, dimension(:, :, :) :: PLE0
      real(kind=REAL8), pointer, dimension(:, :, :) :: PLE1
      ! TODO: pchakrab - this stays
      real(kind=FVPRC), allocatable, dimension(:, :, :) :: xPLE1
      real(kind=FVPRC), pointer, dimension(:, :, :, :) :: tracers
      real(kind=REAL8), allocatable :: tmass0(:), tmass1(:)
      type(AdvCoreTracers), pointer :: adv_tracers(:)

      integer :: im, jm, lm, nq, nqt, QSPLIT
      integer :: i, j, n
      integer, save :: nq_saved = 0
      logical :: exclude, same_tradv_data, rpt_mass, DEBUG_ADV, adjust_tracers

      real, pointer :: temp3D(:, :, :)
      real(kind=REAL4), pointer :: tracer_r4(:, :, :)
      real(kind=REAL8), pointer :: tracer_r8(:, :, :)

      character(len=ESMF_MAXSTR) :: field_name, mytracer
      character(len=:), allocatable :: adjust_tracer_mode
      character(len=ESMF_MAXSTR), allocatable :: xlist(:), biggerlist(:)

      real(kind=FVPRC), allocatable :: DEBUG_ARRAY(:, :, :)
      real(kind=FVPRC) :: fac1 = 1.0
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

      ! ALT: this section attempts to limit the amount of advected tracers
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
         call logger%info("Invalid value specified for EXCLUDE_ADVECTION_TRACERS, ignored")
         adjust_tracers = .false.
      end select
      if (adjust_tracers) then
         if (first_run) then
            xlist = ESMF_HConfigAsStringSeq(hconfig, keyString="EXCLUDE_ADVECTION_TRACERS_LIST", stringLen=ESMF_MAXSTR,&
                 & _RC)
            n = 0
            if (allocated(xlist)) n = size(xlist)
            ! Count the number of tracers
            call ESMF_FieldBundleGet(tradv, grid=bgrid, fieldCount=nqt, _RC)
            bundle_adv = ESMF_FieldBundleCreate(name="xTRADV", _RC)
            call ESMF_FieldBundleSet(bundle_adv, grid=bgrid, _RC)
            ! Loop over NQ in TRADV
            do i = 1, nqt
               ! Get field from TRADV and its name
               call ESMF_FieldBundleGet(tradv, fieldIndex=i, field=field, _RC)
               ! TODO: pchakrab - replace the call below with a call to get_short_name
               call ESMF_FieldGet(field, name=field_name, _RC)
               ! field_name = get_short_name(field, _RC)
               ! Exclude everything that is not cloud/water species
               if ((FV3_DynCoreIsRunning) .and. &
                    (field_name == 'Q') .or. &
                    (field_name == 'QLCN') .or. &
                    (field_name == 'QLLS') .or. &
                    (field_name == 'QICN') .or. &
                    (field_name == 'QILS') .or. &
                    (field_name == 'CLCN') .or. &
                    (field_name == 'CLLS') .or. &
                    (field_name == 'NCPL') .or. &
                    (field_name == 'NCPI') .or. &
                    (field_name == 'QRAIN') .or. &
                    (field_name == 'QSNOW') .or. &
                    (field_name == 'QGRAUPEL')) then
                  ! write(STRING,'(A,A)') "ADV is excluding ", TRIM(field_name)
                  ! call WRITE_PARALLEL( trim(STRING)   )
                  n = n + 1
                  if (n > size(xlist)) then
                     allocate(biggerlist(2 * n), _STAT)
                     biggerlist(1:n - 1) = xlist
                     call move_alloc(from=biggerlist, to=xlist)
                  end if
                  xlist(n) = trim(field_name)
               end if
               ! Loop over exclude_list
               exclude = .false.
               do j = 1, n
                  if (field_name == xlist(j)) then
                     exclude = .true.
                     exit
                  end if
               end do
               if (.not. exclude) then
                  call MAPL_FieldBundleAdd(bundle_adv, [field], _RC)
               end if
            end do

            if (allocated(xlist)) deallocate(xlist)
            if (allocated(biggerlist)) deallocate(biggerlist)

            first_run = .false.
         end if ! first_run
         tradv = bundle_adv
      end if ! adjust_tracers

      call MAPL_GridCompGet(gc, grid=esmfgrid, num_levels=lm, _RC)
      call MAPL_GridGet(esmfgrid, im=im, jm=jm, _RC)
      call ESMF_FieldBundleGet(tradv, fieldCount=nq, _RC)

      if (nq > 0) then
         ! We allocate a list of tracers big enough to hold all items in the bundle
         allocate(tracers(im, jm, lm, nq), _STAT)
         allocate(adv_tracers(nq), _STAT)

         if (nq /= nq_saved) then
            call logger%info("Advcore is Advecting the following %i tracers", nq)
         end if

         ! Go through the bundle copying the friendlies into the tracer list.
         do n = 1, nq
            call ESMF_FieldBundleGet(tradv, fieldIndex=n, field=field, _RC)
            ! TODO: pchakrab - replace the call below with a call to get_short_name
            call ESMF_FieldGet(field, array=array, name=field_name, _RC)
            ! field_name = get_short_name(field, _RC)
            call ESMF_ArrayGet(array, typekind=kind, _RC)
            adv_tracers(n)%is_r4 = (kind == ESMF_TYPEKIND_R4) ! Is real*4?
            adv_tracers(n)%tName = field_name
            if (nq /= nq_saved) call logger%info(trim('--' // field_name))
            if (adv_tracers(n)%is_r4) then
               call ESMF_ArrayGet(array, farrayptr=tracer_r4, _RC)
               adv_tracers(n)%content_r4 => tracer_r4
               tracers(:, :, :, n) = adv_tracers(n)%content_r4
            else
               call ESMF_ArrayGet(array, farrayptr=tracer_r8, _RC)
               adv_tracers(n)%content => tracer_r8
               tracers(:, :, :, n) = adv_tracers(n)%content
            end if
         end do

         if (nq /= nq_saved) then
            nq_saved = nq
         end if

         ! Get Tracer Mass before advection
         call MAPL_GridCompGetResource(gc, "ADV_CORE_REPORT_TRACER_MASS", rpt_mass, default=.false., _RC)
         if (rpt_mass) then
            allocate(tmass0(nq))
            call global_integral(tmass0, tracers, real(PLE0, kind=FVPRC), im, jm, lm, nq)
         end if

         ! Run FV3 advection
         call MAPL_GridCompGetResource(gc, "ADV_QSPLIT", QSPLIT, default=0, _RC)
         call offline_tracer_advection(tracers, &
              real(PLE0, kind=FVPRC), xPLE1, &
              real(MFX, kind=FVPRC), real(MFY, kind=FVPRC), &
              real(CX, kind=FVPRC), real(CY, kind=FVPRC), &
              FV_Atm(1)%gridstruct, FV_Atm(1)%flagstruct, &
              FV_Atm(1)%bd, FV_Atm(1)%domain, FV_Atm(1)%npx, FV_Atm(1)%npy, FV_Atm(1)%npz, &
              nq, dt, QSPLIT)

         ! Get Tracer Mass after advection
         if (rpt_mass) then
            allocate(tmass1(nq))
            call global_integral(tmass1, tracers, real(PLE1, kind=FVPRC), im, jm, lm, nq)
            ! Conserve Specific Mass of Constituents Keeping Mixing_Ratio Constant WRT_Dry_Air
            do n = 1, nq
               if (tmass1(n) > 0.0) then
                  if (ABS((tmass0(n) - tmass1(n)) / tmass1(n)) >= epsilon(1.0_REAL4)) then
                     if (is_master()) call logger%info(trim(adv_tracers(n)%tName) // ': ' // &
                          trim(adjustl(transfer((tmass1(n) - tmass0(n)) / tmass0(n), 'G21.14'))))
                     !!TRACERS(:,:,:,N) = TRACERS(:,:,:,N) * TMASS0(N)/TMASS1(N)
                  end if
                  125 format('Mass Conservation Adjustment in AdvCore:'2x, A, 2x, g21.14)
               end if
            end do
            deallocate(tmass0)
            deallocate(tmass1)
         end if

         ! Go through the bundle copying tracers back to the bundle.
         do n = 1, nq
            if (adv_tracers(n)%is_r4) then
               adv_tracers(n)%content_r4 = tracers(:, :, :, n)
            else
               adv_tracers(n)%content = tracers(:, :, :, n)
            end if

            !-----------------------------------------------
            !--> Fill Export States
            !--> This section is used for diagnostics only.
            !--> It has no effect on CTM experiments.
            !-----------------------------------------------
            if (n <= min(ntracers, nq)) then
               write(mytracer, "('TEST_TRACER',i5.5)") n - 1
               call MAPL_StateGetPointer(export, temp3D, trim(mytracer), _RC)
               if (associated(temp3D)) temp3D = tracers(:, :, :, n)
            end if
         end do

         ! Clean negative tracers and check
         call MAPL_GridCompGetResource(gc, "DEBUG_ADV", DEBUG_ADV, default=.false., _RC)
         if (DEBUG_ADV) then
            prt_minmax = DEBUG_ADV
            if (mpp_pe() == 0) print*, ''
            if (mpp_pe() == 0) print*, '-------------- FV3 Tracer Debug After ADV --------------'
            allocate(DEBUG_ARRAY(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec, FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec, FV_Atm(1)%npz))
         end if
         do n = 1, nq
            if (adv_tracers(n)%is_r4) then
               where (adv_tracers(n)%content_r4 < tiny(0.0))
                  adv_tracers(n)%content_r4 = 0.0
               end where
               if (DEBUG_ADV) DEBUG_ARRAY = adv_tracers(n)%content_r4
            else
               where (adv_tracers(n)%content < tiny(0.0))
                  adv_tracers(n)%content = 0.0
               end where
               if (DEBUG_ADV) DEBUG_ARRAY = adv_tracers(n)%content
            end if
            if (DEBUG_ADV) then
               call prt_maxmin(trim(adv_tracers(n)%tName), DEBUG_ARRAY, FV_Atm(1)%bd%isc, FV_Atm(1)%bd%iec, FV_Atm(1)&
                    %bd%jsc, FV_Atm(1)%bd%jec, 0, FV_Atm(1)%npz, fac1)
            end if
         end do
         if (DEBUG_ADV) then
            deallocate(DEBUG_ARRAY)
            if (mpp_pe() == 0) print*, '-------------- FV3 Tracer Debug After ADV --------------'
            if (mpp_pe() == 0) print*, ''
            prt_minmax = .false.
         end if

         ! Deallocate the list of tracers
         deallocate(tracers, adv_tracers, _STAT)

      end if ! NQ > 0

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

end module AdvCore_GridCompMod

subroutine AdvCore_SetServices(gc, rc)
   use ESMF, only: ESMF_GridComp
   use AdvCore_GridCompMod, only: mySetServices => SetServices
   type(ESMF_GridComp) :: gc
   integer, intent(out) :: rc
   call mySetServices(gc, rc=rc)
end subroutine AdvCore_SetServices