#define ASSERT_(a) if(.not. a) then; print *, __FILE__, __LINE__; call MPI_Abort(MPI_COMM_WORLD, 1, mpierr); endif

module fv_state_utilities_mod

  use fv_arrays_mod , only: fv_atmos_type

  implicit none

  private
  public :: set_resolution_dependent_geos_defaults

contains

  ! Copied from FV_StateMod::FV_Setup
  ! Create some resolution dependent defaults for FV3 in GEOS...
  ! These can be overrided in fv_core_nml in fvcore_layout.rc linked to input.nml
  subroutine set_resolution_dependent_geos_defaults(fix_mass, DT, FV_Atm)

    ! Arguments
    logical, intent(in) :: fix_mass
    real, intent(in) :: DT
    type(fv_atmos_type), intent(inout) :: FV_Atm(:)

    ! Start
    if (size(FV_Atm) /= 1) error stop

    ! Number of water species for FV3 determined later
    ! when reading the tracer bundle in fv_first_run
    FV_Atm(1)%flagstruct%nwat = 0
    ! Veritical resolution dependencies
    FV_Atm(1)%flagstruct%external_eta = .true.
    if (FV_Atm(1)%flagstruct%npz >= 70) then
       FV_Atm(1)%flagstruct%n_sponge = 9   ! ~0.2mb
       FV_Atm(1)%flagstruct%n_zfilter = 18 ! ~10mb
    endif
    if (FV_Atm(1)%flagstruct%npz >= 72) then
       FV_Atm(1)%flagstruct%n_sponge = 9   ! ~0.2mb
       FV_Atm(1)%flagstruct%n_zfilter = 25 ! ~10mb
    endif
    if (FV_Atm(1)%flagstruct%npz >= 90) then
       FV_Atm(1)%flagstruct%n_sponge = 9   ! ~0.2mb
       FV_Atm(1)%flagstruct%n_zfilter = 25 ! ~10mb
    endif
    if (FV_Atm(1)%flagstruct%npz >= 126) then
       FV_Atm(1)%flagstruct%n_sponge = 9   ! ~0.2mb
       FV_Atm(1)%flagstruct%n_zfilter = 23 ! ~10mb
    endif
    if (FV_Atm(1)%flagstruct%npz >= 132) then
       FV_Atm(1)%flagstruct%n_sponge = 9   ! ~0.2mb
       FV_Atm(1)%flagstruct%n_zfilter = 30 ! ~10mb
    endif
    if (FV_Atm(1)%flagstruct%npz >= 136) then
       FV_Atm(1)%flagstruct%n_sponge = 9   ! ~0.2mb
       FV_Atm(1)%flagstruct%n_zfilter = 30 ! ~10mb
    endif
    if (FV_Atm(1)%flagstruct%npz >= 180) then
       FV_Atm(1)%flagstruct%n_sponge = 18  ! ~0.2mb
       FV_Atm(1)%flagstruct%n_zfilter = 50 ! ~10mb
    endif
    FV_Atm(1)%flagstruct%tau = 0.
    FV_Atm(1)%flagstruct%rf_cutoff = 7.5e2
    FV_Atm(1)%flagstruct%d2_bg_k1 = 0.20
    FV_Atm(1)%flagstruct%d2_bg_k2 = 0.06
    FV_Atm(1)%flagstruct%remap_option = 0
    FV_Atm(1)%flagstruct%kord_tm =  9
    FV_Atm(1)%flagstruct%kord_mt =  9
    FV_Atm(1)%flagstruct%kord_wz =  9
    FV_Atm(1)%flagstruct%kord_tr =  9
    FV_Atm(1)%flagstruct%z_tracer = .true.
    ! Some default horizontal flags
    FV_Atm(1)%flagstruct%adjust_dry_mass = fix_mass
    FV_Atm(1)%flagstruct%consv_te = 1.
    FV_Atm(1)%flagstruct%consv_am = .false.
    FV_Atm(1)%flagstruct%fill = .true.
    FV_Atm(1)%flagstruct%dwind_2d = .false.
    FV_Atm(1)%flagstruct%delt_max = 0.002
    FV_Atm(1)%flagstruct%ke_bg = 0.0
    ! Some default damping options
    FV_Atm(1)%flagstruct%nord = 2
    FV_Atm(1)%flagstruct%dddmp = 0.2
    FV_Atm(1)%flagstruct%d4_bg = 0.12
    FV_Atm(1)%flagstruct%d2_bg = 0.0
    FV_Atm(1)%flagstruct%d_ext = 0.0
    ! Some default time-splitting options
    FV_Atm(1)%flagstruct%n_split = 0
    FV_Atm(1)%flagstruct%k_split = 1
    ! default NonHydrostatic settings (irrelavent to Hydrostatic)
    FV_Atm(1)%flagstruct%beta = 0.0
    FV_Atm(1)%flagstruct%a_imp = 1.0
    FV_Atm(1)%flagstruct%p_fac = 0.1
    ! Cubed-Sphere Global Resolution Specific adjustments
    if (FV_Atm(1)%flagstruct%ntiles == 6) then
       ! Cubed-sphere grid resolution and DT dependence
       !              based on ideal remapping DT
       if (FV_Atm(1)%flagstruct%npx >= 48) then
          FV_Atm(1)%flagstruct%k_split = CEILING(DT/1800.0  )
       endif
       if (FV_Atm(1)%flagstruct%npx >= 90) then
          FV_Atm(1)%flagstruct%k_split = CEILING(DT/ 900.0   )
       endif
       if (FV_Atm(1)%flagstruct%npx >= 180) then
          FV_Atm(1)%flagstruct%k_split = CEILING(DT/ 450.0   )
       endif
       if (FV_Atm(1)%flagstruct%npx >= 360) then
          FV_Atm(1)%flagstruct%k_split = CEILING(DT/ 225.0   )
       endif
       if (FV_Atm(1)%flagstruct%npx >= 720) then
          FV_Atm(1)%flagstruct%k_split = CEILING(DT/ 112.5   )
       endif
       if (FV_Atm(1)%flagstruct%npx >= 1440) then
          FV_Atm(1)%flagstruct%k_split = CEILING(DT/  37.5   )
       endif
       if (FV_Atm(1)%flagstruct%npx >= 2880) then
          FV_Atm(1)%flagstruct%k_split = CEILING(DT/  18.75  )
       endif
       if (FV_Atm(1)%flagstruct%npx >= 5760) then
          FV_Atm(1)%flagstruct%k_split = CEILING(DT/   9.375 )
       endif
       if (FV_Atm(1)%flagstruct%npz == 72) then
          FV_Atm(1)%flagstruct%fv_sg_adj = DT
          ! Monotonic Hydrostatic defaults
          FV_Atm(1)%flagstruct%hydrostatic = .true.
          FV_Atm(1)%flagstruct%make_nh = .false.
          FV_Atm(1)%flagstruct%vtdm4 = 0.0
          FV_Atm(1)%flagstruct%do_vort_damp = .false.
          FV_Atm(1)%flagstruct%d_con = 0.
          FV_Atm(1)%flagstruct%hord_mt =  10
          FV_Atm(1)%flagstruct%hord_vt =  10
          FV_Atm(1)%flagstruct%hord_tm =  10
          FV_Atm(1)%flagstruct%hord_dp =  10
          ! This is the best/fastest option for tracers
          FV_Atm(1)%flagstruct%hord_tr =  8
          ! NonMonotonic defaults for c360 (~25km) and finer
          if (FV_Atm(1)%flagstruct%npx >= 360) then
             ! This combination of horizontal advection schemes is critical
             ! for anomaly correlation NWP skill.
             ! Using all = 5 (like GFS) produces a substantial degredation in skill
             FV_Atm(1)%flagstruct%hord_mt =  5
             FV_Atm(1)%flagstruct%hord_vt =  6
             FV_Atm(1)%flagstruct%hord_tm =  6
             FV_Atm(1)%flagstruct%hord_dp = -6
             ! Must now include explicit vorticity damping
             FV_Atm(1)%flagstruct%d_con = 1.
             FV_Atm(1)%flagstruct%do_vort_damp = .true.
             FV_Atm(1)%flagstruct%vtdm4 = 0.02
          endif
       else
          FV_Atm(1)%flagstruct%fv_sg_adj = DT
          ! Monotonic Hydrostatic defaults
          FV_Atm(1)%flagstruct%hydrostatic = .false.
          FV_Atm(1)%flagstruct%make_nh = .false.
          FV_Atm(1)%flagstruct%vtdm4 = 0.0
          FV_Atm(1)%flagstruct%do_vort_damp = .false.
          FV_Atm(1)%flagstruct%d_con = 0.
          FV_Atm(1)%flagstruct%hord_mt =  10
          FV_Atm(1)%flagstruct%hord_vt =  10
          FV_Atm(1)%flagstruct%hord_tm =  10
          FV_Atm(1)%flagstruct%hord_dp =  10
          ! This is the best/fastest option for tracers
          FV_Atm(1)%flagstruct%hord_tr =  8
          ! NonMonotonic defaults for c360 (~25km) and finer
          if (FV_Atm(1)%flagstruct%npx >= 360) then
             FV_Atm(1)%flagstruct%hord_mt =  6
             FV_Atm(1)%flagstruct%hord_vt =  6
             FV_Atm(1)%flagstruct%hord_tm =  6
             FV_Atm(1)%flagstruct%hord_dp = -6
             ! Must now include explicit vorticity damping
             FV_Atm(1)%flagstruct%d_con = 1.
             FV_Atm(1)%flagstruct%do_vort_damp = .true.
             FV_Atm(1)%flagstruct%vtdm4 = 0.02
          endif
          ! continue to adjust vorticity damping with
          ! increasing resolution
          if (FV_Atm(1)%flagstruct%npx >= 1440) then
             FV_Atm(1)%flagstruct%vtdm4 = 0.04
          endif
          if (FV_Atm(1)%flagstruct%npx >= 2880) then
             FV_Atm(1)%flagstruct%vtdm4 = 0.06
          endif
          if (FV_Atm(1)%flagstruct%npx >= 5760) then
             FV_Atm(1)%flagstruct%vtdm4 = 0.08
          endif
       endif
    endif

  end subroutine set_resolution_dependent_geos_defaults

end module fv_state_utilities_mod
