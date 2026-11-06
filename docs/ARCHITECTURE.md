# FVdycoreCubed_GridComp: MAPL wrapper of the FV3 dynamical core

**Audience:** coding agents and developers working in the `FVdycoreCubed_GridComp` repository.
**Describes:** the MAPL2 line, branch `develop`.
**Verified against:** `FVdycoreCubed_GridComp` @ `4deed15` (2026-10-27).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `0f2e0fb`).
Other fixtures (for example GEOSadas, GEOSldas) may place this repository
elsewhere; paths inside this repository do not change.

Paths are relative to the `FVdycoreCubed_GridComp` repository root. Parent
composite (who calls this component, in which order):
`GEOSagcm_GridComp/ARCHITECTURE.md` (`../../../ARCHITECTURE.md` in the fixture), section on
Superdyn. Repository overview one level up:
`GEOSgcm_GridComp/docs/ARCHITECTURE.md` (`../../../../docs/ARCHITECTURE.md` in the fixture).

"Agcm" is `GEOS_AgcmGridComp` (the atmosphere composite); "Superdyn" is
`GEOS_SuperdynGridComp`. A `(file, routine)` citation means the statement was
checked there; *Inferred:* and *Unverified:* mark weaker statements.

## 1. Purpose and placement

This repository is the MAPL gridded-component layer around the GFDL FV3
(finite-volume cubed-sphere) dynamical core. It does not contain the core
itself: the core is the nested repository `@fvdycore` (section 9), which is
built into the same library (`CMakeLists.txt`, `SUBCOMPONENTS fvdycore`).

In the GEOSgcm fixture the repository is at
`src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSsuperdyn_GridComp/@FVdycoreCubed_GridComp`
(`components.yaml`, key `FVdycoreCubed_GridComp`, tag `v3.1.0`). Its parent is
`GEOS_SuperdynGridComp`, which imports `SetServices` from the Fortran module
`FVdycoreCubed_GridComp` (the `DYN` component, sections 3 to 6) and from
`AdvCore_GridCompMod` (the `ADV` component, section 7)
(`GEOS_SuperdynGridComp.F90`, `SetServices`).

## 2. Directory and file map

| Path | Content |
|---|---|
| `DynCore_GridCompMod.F90` | `DYN` component, about 9300 lines. Only `SetServices` is public |
| `FV_StateMod.F90` | `FV_Setup`, `FV_InitState`, `FV_Run`, `FV_Finalize`, `State_To_FV`, `FV_To_State`; module-level array `FV_Atm` |
| `AdvCore_GridCompMod.F90` | `ADV` component |
| `StandAlone_*.F90`, `c2c.F90`, `interp_restarts*.F90`, `rs_scale*.F90` | Standalone drivers and restart tools (section 10) |
| `Cube*.F90`, `LatLon2Cube.F90`, `AppGridCreate.F90`, `GetWeights*.F90`, `fv_regrid*.F90`, `CreateInterpWeights_GridCompMod.F90` | Grid creation and regridding helpers |
| `sw.f90`, `jw.f90`, `testcases_3_4_5_6_stand_alone.f90` | Idealized test cases used by the `Coldstart` path |
| `GEOS_FV3_Utilities.F90` | Utility module (`a2d2c`) |
| `.circleci/config.yml`, `.github/` | CI (jobs `ci/build`, `ci/run_fv3`, `ci/run_gcm`), `CODEOWNERS`, workflows |

Not listed in `CMakeLists.txt` (not built here): `c2c_rst.F90`,
`c2l_CFIO_offline.F90`, `CreateTopography.F90`, `CreateInterpWeights.F90` (a
program using the library module `CreateInterpWeights_GridCompMod`),
`dyn_interp_rst.F90`, `moist_interp_rst.F90`, `ESMF_GenerateCSGridDescription.F90`,
`.cppdefs`, `FVdycore_arch.mk`. *Inferred:* legacy or manually built tools.
GEOSsurface `generate_scrip_cube.F90` mentions `ESMF_GenerateCSGridDescription`.

## 3. The DYN component

### 3.1 How it is called

Superdyn (the parent) chooses the dycore through the resource `DYCORE:`
(default `FV3`): `FV3` adds `DYN` only, `FV3+ADV` adds `DYN` and `ADV`, and
`DATMO` adds a data-atmosphere `DYN` from another component
(`GEOS_SuperdynGridComp.F90`, `SetServices`). `DYN` registers two Run phases
(`DynCore_GridCompMod.F90`, `SetServices`):

| Phase | Routine | Called by | Does |
|---|---|---|---|
| Run 1 | `Run` | Superdyn `Run`, then (if `ADV` exists) `ADV` | Applies the analysis/bias tendencies, advances the FV3 state by one dynamics step, fills exports |
| Run 2 | `RunAddIncs` | Superdyn `RunAddIncs` | Adds the physics tendencies to the FV3 state, recomputes exports |

Order in the Agcm `Run` (`GEOS_AgcmGridComp.F90`):

1. `DO_UPDATE_ANA3D` fills the Superdyn imports `DUDT`, `DVDT`, `DPEDT`,
   `DTDT` by run type: zero for `FREERUN`, the stored bias tendency for
   `PREDICTOR` and `FORECAST`, the IAU increment for `CORRECTOR`. `DQVANA`,
   `DQLANA`, `DQIANA`, `DQRANA`, `DQSANA`, `DQGANA` are set in the per-tracer
   branch for `Q` (`DQVANA = Q`, `DQLANA = QLLS + QLCN`, ..., later differenced
   against the updated state), and `DOXANA` in the `DO3DT` branch.
2. Superdyn phase 1: `Run` copies its imports `DUDT`, `DVDT`, `DTDT`, `DWDT`
   into the `DYN` imports and runs `DYN` phase 1, which applies them (3.4).
3. Physics; then `DO_UPDATE_PHY` overwrites the Superdyn imports `DUDT`,
   `DVDT`, `DWDT`, `DPEDT`, `DTDT` with the Physics exports.
4. Superdyn phase 2: `RunAddIncs` copies the four imports again and runs
   `DYN` phase 2, which applies the physics tendencies (3.5).

Superdyn terminates the `DYN` imports `DUDT`, `DVDT`, `DTDT`, `DWDT` so that no
framework connection fills them (`GEOS_SuperdynGridComp.F90`, `SetServices`).
Superdyn has import specs only for those four names. The other `DYN` imports
(`DPEDT`, the `DQ*ANA` and `DOXANA` increments, `PHIS`, `VARFLT`, `TRADV`) have no
Superdyn spec; the Agcm terminates the `SDYN` imports `DUDT`, `DVDT`, `DWDT`,
`DTDT`, `DPEDT`, `DQ*ANA`, `DOXANA`, `PHIS`, `VARFLT` as manually filled
(`GEOS_AgcmGridComp.F90`, `SetServices`, `MAPL_TerminateImport` with
`CHILD = SDYN`). *Inferred:* the `DYN` imports reach the Agcm level by MAPL
promotion through Superdyn. `PHIS` and `VARFLT` are filled in the Agcm
`Initialize` (`GEOS_TopoGet`).

### 3.2 SetServices

(`DynCore_GridCompMod.F90`, `SetServices`)

- Imports (15 specs): `DUDT`, `DVDT`, `DWDT`, `DTDT`, `DPEDT` (tendencies
  set by the Agcm: bias, increment or physics), `DQVANA`, `DQLANA`, `DQIANA`, `DQRANA`, `DQSANA`, `DQGANA`,
  `DOXANA` (analysis increments), `PHIS`, `VARFLT` (topography); and `TRADV`,
  a bundle (`MAPL_BundleItem`, section 6). Only `VARFLT` sets
  `RESTART=MAPL_RestartSkip`; the other imports set no `RESTART` and take the
  MAPL2 default (`MAPL_RestartOptional`).
- Exports (about 236 specs): `CX`, `CY`, `MFX`, `MFY`, `PLE0`, `PLE1` (inputs of
  `ADV`), the state (`U`, `V`, `T`, `Q`, `W`, `PLE`, `ZLE`, `PS`, ...), grid and
  vertical coordinates (`PREF`, `AK`, `BK`, `AREA`, `DXC`, `DYC`), the
  `*_DYN_IN` set promoted to Physics, and tendency and energetics diagnostics.
- Internal state (section 5).
- Allocates the private state `DynState` (internal state `'DYNstate'`); calls
  `DynSetup` (alias for `FV_Setup`; reads `LAYOUT:`) and
  `register_grid_and_regridders` (`CubeGridPrototype`: bilinear regridders for
  the MAPL `regridder_manager`).
- Reads `FV3_STANDALONE:` (nonzero: own grid through `MAPL_GridCreate` and a
  `TRADVEX` export bundle), `DEBUG_DYN:`, `DEBUG_ADV:`, `DEBUG_TQ_ERRORS:`.
- Ends with `MAPL_GenericSetServices`.

### 3.3 Initialize

(`DynCore_GridCompMod.F90`, `Initialize`)

1. `MAPL_GenericInitialize` (reads the internal restart, section 5).
2. Reads `DO_ADD_INCS:` and `COLDSTART:`. A nonzero `COLDSTART` calls `Coldstart`,
   which builds an idealized state from the `CASE_*` keys, `T0:`, and optional
   `AK:`/`BK:` tables (or `set_eta`).
3. `DynInit` (alias for `FV_InitState`, `FV_StateMod.F90`): connects the
   `STATE%VARS` pointers to the MAPL internal-state arrays (`CREATE_VARS`),
   creates the alarm "Time2Run" at `RUN_DT:` and a six-hour dry-mass alarm,
   sets `FV_Atm%phis`, `varflt`, `ak`, `bk`, calls `State_To_FV`, and under
   `RUN_GTFV3` calls `geos_gtfv3_interface_f_init`.
4. Fills the exports `PREF`, `AK`, `BK`, `PLE`, `U`, `V`, `T`, `DXC`, `DYC`,
   `AREA`; marks `PREF`, `PLE`, `U`, `V`, `T` with `MAPL_InitialRestart`.
5. For `FV3_STANDALONE`, copies the `TRADV` fields to `TRADVEX`; for
   `REPLAY_MODE: Intermittent`, creates the alarm `INTERMITTENT`.

### 3.4 Run (phase 1)

(`DynCore_GridCompMod.F90`, `Run`)

1. Returns at once if the "Time2Run" alarm is not ringing (the dycore step may
   be longer than the parent's step, `RUN_DT:`).
2. Gets the `TRADV` bundle from the import state and selects the tracers the
   dycore advects (section 6); `PULL_Q` points `STATE%VARS%tracer(:)` at the
   selected bundle fields.
3. The "WMP REPLAY/ANA" section (skipped when `FV3_STANDALONE /= 0` or the
   `ReplayShutOff` clock alarm is ringing): intermittent replay (`dump_n_splash_`
   for `REPLAY_TYPE FULL`, `incremental_` otherwise, then `state_remap_`; reads
   the `REPLAY_*` keys, `NX:`, `NY:`). Then, if `DO_ADD_INCS` (default true)
   and not `ADIABATIC`, the tendencies the Agcm placed in the imports
   (`DUDT`, `DVDT`, `DPEDT`, `DTDT`: bias or IAU increment) are applied with
   `ADD_INCS` (`ANA_IS_WEIGHTED:` selects mass-weighted temperature
   tendencies). `DQVANA` etc. are used to rebuild the pre-increment moisture
   state (`QVOLD = QV - DQVANA`) for mass bookkeeping and diagnostics;
   `ADD_INCS` does not change the tracers.
4. Fills the `*_DYN_IN` exports (state at the start of the step).
5. Calls `DynRun` (alias for `FV_Run`, section 4) with `PLE0` the edge
   pressure before the step.
6. Fills the exports (`PLE0`, `PLE1`, `CX`, `CY`, `MFX`, `MFY`, `U`, `V`, `T`,
   `Q`, `PLE`, tendency diagnostics, ...). The mass fluxes come from
   `fillMassFluxes`, or from `computeMassFluxes` (time-centered) when
   `AdvCore_Advection == 2`. Then tropopause, EPV, vorticity and
   pressure-level diagnostics, and `freeTracers`.

### 3.5 RunAddIncs (phase 2)

(`DynCore_GridCompMod.F90`, `RunAddIncs`)

- Skipped if `SW_DYNAMICS`. Only `Run` checks the Time2Run alarm;
  `RunAddIncs` has no such guard.
- `DWDT` is copied into the `DYN` imports, but `ADD_INCS` uses it only in
  commented-out code.
- `PULL_Q`, then (if `DO_ADD_INCS` and not `ADIABATIC`) `ADD_INCS` adds the
  physics tendencies to the FV3 state, the same routine as in phase 1: the
  winds are converted from the A grid to the native D grid with
  `Agrid_To_Native` (called with `wind_increment_limiter = 800.d0/86400.d0`,
  see section 15); `PE` is advanced by `DT*DPEDT`; `PT` is re-weighted by old
  and new `DP` and `cp/cvm`; `getPKZ` recomputes `PKZ`. The tracers are not
  changed here.
- Recomputes the exports that depend on the new state (`U_DGRID`, `ZLE`, `PLE`,
  `PKE`, ..., `KEPHY`, `PEPHY`, `TEPHY`), then `freeTracers`.

### 3.6 Finalize

`DynFinalize` (alias for `FV_Finalize`: `fv_end`, and the gtFV3 finalize under
`RUN_GTFV3`), then `MAPL_GenericFinalize`.

## 4. FV_StateMod: the bridge to the FV3 core

(`FV_StateMod.F90`)

`FV_Atm(:)` is a module-level `save` array of the FV3 type `fv_atmos_type`
(`@fvdycore`, `fv_arrays_mod`). `AdvCore_GridCompMod` uses the same array
(`use FV_StateMod, only: FV_Atm`), so `DYN` and `ADV` share one FV3 state.

**`FV_Setup`** (called from `SetServices`): `fms_init`, `fv_init1`; reads the
grid size (`AGCM_IM:`, `AGCM_JM:`, `AGCM_LM:`; `npz==1` selects `SW_DYNAMICS`),
`NX:`, `NY:`, `RUN_DT:`; reads `DYCORE:` and sets `AdvCore_Advection` (0 for
`FV3`, 1 for `FV3+ADV`), which the resource `AdvCore_Advection:` may override
(2 selects time-centered mass fluxes); reads `ADIABATIC:`, `FV_OFF:`,
`fix_mass:` (default 1), `FV3_CONFIG:` (default `STOCK`; a string containing
`MONOTONIC` or `MERRA-2` selects hydrostatic, monotonic schemes) and
resolution-dependent defaults; then `fv_init2`.

**`FV_Run`** (called every dynamics step, first time with extra set-up):

1. First call, and only when not `ADIABATIC`: counts the water tracers
   (`Q`, `QLCN`, `QLLS`, `QICN`, `QILS` must be present: asserts
   `nwat_tracers == 5`; then `CLLS`, `CLCN`, `QRAIN`, `QSNOW`, `QGRAUPEL` are
   counted too) and, if `flagstruct%nwat` is 0, sets it to 1 (3 if
   non-hydrostatic; 6 when `nwat_tracers >= 10`) and `STATE%VARS%nwat`.
2. Sets `FV_Atm%phis` and `varflt` from the imports; packs the tracers into
   `FV_Atm%q` (section 6).
3. `State_To_FV` copies the MAPL arrays into the FV3 state; dry mass fixed to
   `PSDRY:`, mass check/fix under `MASSALARM`; first non-hydrostatic call:
   `Make_NH`, `na_init`, `adiabatic_init`.
4. `fv_dynamics` (or `geos_gtfv3_interface_f` when `run_gtfv3 /= 0`, section 11).
5. Fills the Rayleigh-friction tendencies (`DUDT_RAY`, ...); `fv_subgrid_z` when
   `fv_sg_adj > 0`; unpacks the tracers back; `FV_To_State`.

If `ADIABATIC:` is nonzero the tracers are passed through in order.

## 5. State and restarts

Internal state of `DYN` (`DynCore_GridCompMod.F90`, `SetServices`):

| Name | Notes |
|---|---|
| `AK`, `BK` | Vertical-only, edge; hybrid coordinate coefficients; restart required |
| `U`, `V`, `PT`, `PE`, `PKZ` | D-grid winds, scaled potential temperature (`scaled_potential_temperature`, K Pa^-kappa), edge pressure, `pk` exponent; double precision; restart required |
| `DZ`, `W` | Non-hydrostatic; double precision; default restart behavior (optional) |

`FV_InitState` wires the FV3 state arrays to these MAPL arrays
(`FV_StateMod.F90`, `FV_InitState`, `CREATE_VARS`), so after `MAPL_GenericInitialize`
reads a restart, the FV3 core sees those values.

Keys in the GEOSgcm fixture (`AGCM.rc.tmpl`): `DYN_INTERNAL_RESTART_FILE:
fvcore_internal_rst`, `DYN_INTERNAL_CHECKPOINT_FILE: fvcore_internal_checkpoint`,
`DYN_INTERNAL_HEADER: 1`. The dycore restart contains no tracers: tracers live in
the internal states of the components that own them (for example moist and
chemistry), and the dycore only borrows pointers to them through `TRADV`.
`freeTracers` only deallocates and nullifies the FV tracer pointers (a comment
there says "Comment out to output tracer to checkpoint file").

The `Coldstart` path (`COLDSTART:`, section 3.3) builds an idealized state
instead of reading a restart; `scripts/fv3.j` uses it (section 12).

## 6. Tracers: the `TRADV` bundle

`TRADV` is the dycore import bundle of all advected tracers. In the Agcm it is
connected from Physics to Superdyn by `MAPL_AddConnectivity` with
`SRC_NAME`/`DST_NAME` `'TRADV'` (`GEOS_AgcmGridComp.F90`, `SetServices`);
Superdyn declares no `TRADV` import spec of its own
(`GEOS_SuperdynGridComp.F90`, `SetServices`). *Inferred:* `DYN` (and `ADV`)
receive it through the MAPL promotion of unsatisfied child imports to the
parent.

Selection in `DYN` (`DynCore_GridCompMod.F90`, `Run`):

- `EXCLUDE_ADVECTION_TRACERS:` (default `ALWAYS`; `PREDICTOR` uses the clock
  alarm `PredictorAlarm`) and `EXCLUDE_ADVECTION_TRACERS_LIST:` name tracers to
  remove from the advected set.
- When `AdvCore_Advection >= 1` (`ADV` advects the other tracers), every
  tracer except `Q`, `QLCN`, `QLLS`, `QICN`, `QILS`, `CLCN`, `CLLS`, `NCPL`,
  `NCPI`, `QRAIN`, `QSNOW`, `QGRAUPEL` is added to the exclude list. The
  remaining fields are copied into a saved module bundle `bundleAdv` ('xTRADV').
- `PULL_Q` points `STATE%VARS%tracer(:)` at the bundle fields (no copy until
  `FV_Run` packs them).

Packing in `FV_Run` (`FV_StateMod.F90`):

- The cloud species `QLCN`/`QLLS` are summed into the liquid slot, `QICN`/`QILS`
  into the ice slot, `CLCN`/`CLLS` into the cloud-fraction slot. The convective
  part is carried as extra tracers (`idx_qlcn`, `idx_qicn`, `idx_clcn`);
  `QRAIN`, `QSNOW`, `QGRAUPEL` map to their own slots; all other tracers take the
  next free index; the final count must equal `ncnst` (asserted).
- After `fv_dynamics` the large-scale part is recovered as total minus
  convective, bounded below by zero (`MAX(.,0)`).

Consequence: with `DYCORE: FV3` every tracer in `TRADV` (minus the exclude
list) is advected inside `fv_dynamics`. With `FV3+ADV` only the species listed
above stay in `fv_dynamics`; all others are handled by `ADV` (section 7).

## 7. The ADV component

(`AdvCore_GridCompMod.F90`)

Offline tracer advection driven by the mass fluxes of `DYN`; present only with
`DYCORE: FV3+ADV` (or in the standalone drivers).

- Imports: `MFX`, `MFY`, `CX`, `CY`, `PLE0`, `PLE1` (double precision;
  connected from the `DYN` exports by Superdyn `SetServices`, `FV3+ADV` branch)
  and the bundle `TRADV` (the Agcm connection). Exports: `AREA` and 38
  diagnostic `TEST_TRACER*` fields.
- `SetServices` registers `Initialize`, `Run`, `Finalize`; reads `DYCORE:`,
  `AdvCore_Advection:`, `RUN_DT:`, `ADV_QSPLIT:`; calls `fms_init`. In standalone
  mode (`DYCORE` not `FV3`/`FV3+ADV`) it also calls `fv_init1`/`fv_init2`.
- `Run` (only if `AdvCore_Advection > 0`): copies the imports and builds the
  tracer set `xTRADV`. When the dycore runs (`FV3_DynCoreIsRunning`), it
  excludes the species that `DYN` keeps (`Q`, `QLCN`, `QLLS`, `QICN`, `QILS`,
  `CLCN`, `CLLS`, `NCPL`, `NCPI`, `QRAIN`, `QSNOW`, `QGRAUPEL`) and the exclude
  list. It then calls `offline_tracer_advection(TRACERS, PLE0, PLE1, MFX, MFY,
  CX, CY, ..., dt, QSPLIT)` and sets values below `tiny` to zero.

The tracer selection is written twice, once in `DYN` and once in `ADV`, with
opposite polarity (`DYN` excludes everything not in the water/cloud list; `ADV`
excludes the water/cloud list). *Inferred:* the two lists must be kept
consistent by hand (section 14); no code enforces the agreement.

## 8. Configuration keys (structure and scheme selection)

MAPL resource keys (`AGCM.rc`): `DYCORE:`, `AdvCore_Advection:` (components,
sections 3.1 and 7); `AGCM_IM:`, `AGCM_JM:`, `AGCM_LM:`, `NX:`, `NY:`, `RUN_DT:`
(grid, step); `FV3_CONFIG:`, `ADIABATIC:`, `FV_OFF:`, `ADJUST_SPLITS:`,
`RUN_GTFV3:` (scheme); `COLDSTART:` and `CASE_*`; `DO_ADD_INCS:`,
`ANA_IS_WEIGHTED:`, `REPLAY_MODE:`; `EXCLUDE_ADVECTION_TRACERS[_LIST]:`.

## 9. `@fvdycore`

`@fvdycore` is the nested repository holding the GFDL atmospheric cubed-sphere
core with GEOS modifications (`components.yaml`, key `fvdycore`, tag `geos/v3.1.0`, develop `geos/develop`).
Its `CMakeLists.txt` builds a library `fvdycore` from `model/`, `tools/` and
`geos_utils/`, added here as a sub-component. It is not part of this repository
(`.gitignore` ignores `/@fvdycore`) and is not documented here; see
`@fvdycore/README.md` (`../@fvdycore/README.md` in the fixture) (the upstream GFDL README).

## 10. Standalone drivers and tools

Executables built by `CMakeLists.txt`: `StandAlone_FV3_Dycore.x` (`MAPL_Cap`
with `DYN` as the root; used by `scripts/fv3.j`); `StandAlone_AdvCore.x`
(`ADV` as root); `StandAlone_DynAdvCore.x` (a small composite with children
`DYN` and `ADV`, connecting `MFX`, `MFY`, `CX`, `CY`, `PLE0`, `PLE1`; it runs
`DYN` phase 1 only, with no `PHASE` argument, then `ADV`); `rs_scale.x` (inserts
a dry mass into restarts); `c2c.x` (`c2c.F90`, program `gmao_regrid`: cube to
cube regridding); `interp_restarts.x` and `interp_restarts_bin.x` (interpolate
FV and moist restarts to a cubed sphere, with optional new vertical levels);
`fv3_driver.x` (`geos-gtfv3/driver/`, only if `BUILD_GEOS_GTFV3_INTERFACE`).

## 11. `geos-gtfv3/`

An optional bridge to a Python implementation of the FV3 core (the
NDSL/Pace "GeosDycoreWrapper", imported in `geos_gtfv3.py`). It is built only
when `BUILD_GEOS_GTFV3_INTERFACE` is on (default off; `CMakeLists.txt`), which
defines `RUN_GTFV3`, runs `geos-gtfv3/geos_gtfv3_interface.py` with `cffi` under
`mpirun -np 1` to produce `libgeos_gtfv3_interface_py.so`, and compiles
`geos_gtfv3_interface.f90` and `geos_gtfv3_interface.c`. Files:
`geos_gtfv3_interface.f90` (Fortran bindings `geos_gtfv3_interface_f_init`, `_f`,
`_f_finalize`), the `.c` C shim, `geos_gtfv3.py` (`geos_gtfv3_init`,
`geos_gtfv3`, `geos_gtfv3_finalize`, class `GEOSGTFV3`), further Python files
(`f_py_conversion.py`, `geos_gtfv3_debug.py`, `cuda_profiler.py`), `driver/` (a
standalone driver, with `driver/setenv/pyenv.sh`) and `README.md`
(`GTFV3_BACKEND`, default `gt:gpu`; `GTFV3_NAMELIST`, default `input.nml`).

At run time the resource `RUN_GTFV3:` selects the Python core in `FV_Run`
instead of `fv_dynamics`. The Fortran path passes only the first seven
tracers to Python (`q(:,:,:,1:7)`, `nq_tot=7`).

## 12. `scripts/`

`scripts/CMakeLists.txt` installs `fv3_setup` (configured with the precision,
hydrostatic and MPI-stack settings) into `bin/`, `fv3.j`, and a configured
`.FV3_VERSION` into `etc/`. `fv3_setup` (`tcsh`, interactive) prompts for
experiment id, resolution and levels and substitutes values into `fv3.j`, the
batch job. `fv3.j` writes `CAP.rc` (`ROOT_NAME: DYN`), `AGCM.rc`
(`FV3_STANDALONE: 1`, `DYCORE: FV3`, `COLDSTART: 1`, `CASE_ID: 1`,
`ADIABATIC: 1`), `HISTORY.rc` and `input.nml`, then runs
`StandAlone_FV3_Dycore.x`.

## 13. Build

`CMakeLists.txt` builds one library `FVdycoreCubed_GridComp` (`esma_add_library`,
`SUBCOMPONENTS fvdycore`), the executables of section 10, and
`add_subdirectory(scripts)`. Dependencies: `MAPL`, `GFTL_SHARED::gftl-shared`,
`GMAO_hermes` (source of `m_set_eta`, `m_topo_remap`), `GEOS_Shared`, `esmf`,
OpenMP. `FV_PRECISION` (`R4`, `R4R8`, `R8`; set in the fixture-level
`CMakeLists.txt`, default `R4R8`): `R4` and `R4R8` add `SINGLE_FV OVERLOAD_R4`;
`R4` links `FMS::fms_r4`; `R4R8` takes the `FMS::fms_r4` include directories and
links `$<LINK_ONLY:FMS::fms_r8>`; `R8` adds `FREAL8` flags and links
`FMS::fms_r8`. Other options: `BUILD_GEOS_GTFV3_INTERFACE`, `CRAY_POINTER`.

## 14. Known traps

- *Inferred:* the two tracer lists of section 7 must agree; no code checks it.
- `FV_Run` asserts exactly five water tracers when not `ADIABATIC` and requires
  `Q`, `QLCN`, `QLLS`, `QICN`, `QILS` in `TRADV`.
- Only `Run` checks the Time2Run alarm: when it is not ringing `Run` returns
  early and does not refresh the exports, while `RunAddIncs` runs on every
  call. *Inferred:* `RUN_DT:` is meant to be a multiple of the parent step.

## 15. Observations (noted, not fixed)

- `Run` reads `CONSV:` and `FILL:` into `LCONSV` and `LFILL`, which are not
  used afterwards (`DynCore_GridCompMod.F90`, `Run`). `FV_Setup` sets the
  FV3 flags itself (`flagstruct%fill`, `consv_te`).
- In `Run` the timers `"TOTAL"` and `"RUN"` are switched on before the early
  `return` for the Time2Run alarm and not switched off on that path.
- The module header comment says the component "supports only a 1D
  latitude-based decomposition" and mentions PILGRIM.
- The body of `FV_DA_Incs` is commented out. In `a2d3d` (`FV_StateMod.F90`)
  the code that uses the optional argument `wind_increment_limiter` is
  commented out, so the value `800.d0/86400.d0` passed by `ADD_INCS` has no
  effect; no other use of the argument was found in the repository.
- The export `PID` is set to 0 ("WMP need to get from MAPL gid").
- The gtFV3 path passes only seven tracers (section 11) while `FV_Atm%q` holds
  all packed tracers.
- The `Coldstart` entry point in `SetServices` is commented out; source files
  not referenced by `CMakeLists.txt` (section 2) have no build path here.
