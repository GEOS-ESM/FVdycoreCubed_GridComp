# FVdycoreCubed_GridComp: MAPL wrapper of the FV3 dynamical core

**Audience:** coding agents and developers working in the `FVdycoreCubed_GridComp` repository.
**Describes:** the MAPL3 line, branch `release/MAPL-v3`.
**Verified against:** `FVdycoreCubed_GridComp` @ `dadf5b5` (2026-11-06).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `f50622d`).
Other fixtures (for example GEOSadas, GEOSldas) may place this repository
elsewhere; paths inside this repository do not change.

> **Build status on this line:** built and ported to MAPL3.
> `GEOSsuperdyn_GridComp/CMakeLists.txt` lists `FVdycoreCubed_GridComp` as the
> only entry of `alldirs` (`GEOSdatmodyn_GridComp` is commented out: "not yet
> ported to MAPL3") and builds it as a `SUBCOMPONENT` of the library
> `GEOS_SuperdynGridComp`. Superdyn adds `DYN` (config file `dyn.yaml`) and
> runs it with the phases `run` and `run_add_incs`. `ADV` is compiled into this
> repository's library but Superdyn does not add it (section 7).

Paths are relative to the `FVdycoreCubed_GridComp` repository root. Parent
composite: `GEOSagcm_GridComp/ARCHITECTURE.md` (`../../../ARCHITECTURE.md` in
the fixture), section on Superdyn. Repository overview one level up:
`GEOSgcm_GridComp/docs/ARCHITECTURE.md` (`../../../../docs/ARCHITECTURE.md` in the fixture).

"Agcm" is `GEOS_AgcmGridComp`; "Superdyn" is `GEOS_SuperdynGridComp`. A
`(file, routine)` citation means the statement was checked there; *Inferred:*
and *Unverified:* mark weaker statements.

## 0. Differences from the MAPL2 line

- Specs are in `DynCore_StateSpecs.rc` and `AdvCore_StateSpecs.rc` (`mapl_acg`), not `MAPL_Add*Spec` calls; `DUDT`, `DVDT` are the vector import `D_UV_DT`.
- Entry points are phases `run` and `run_add_incs`; Superdyn copies no tendencies into `DYN`.
- Grid size and layout come from the config `geometry:` section, not `AGCM_IM`, `NX` and similar keys.
- `DYN` and `ADV` share one tracer predicate (`fv_shared.F90`).
- Intermittent replay, `TRADVEX`, the standalone executables and many regridding helper sources are gone; regression cases `adv-dyn`, `dyn-sa` exist.

## 1. Purpose and placement

This repository is the MAPL gridded-component layer around the GFDL FV3
cubed-sphere dynamical core. The core itself is the nested repository
`@fvdycore` (section 9), built into the same library (`CMakeLists.txt`,
`SUBCOMPONENTS fvdycore`).

In the GEOSgcm fixture the repository is at
`src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSsuperdyn_GridComp/@FVdycoreCubed_GridComp`
(`components.yaml`, key `FVdycoreCubed_GridComp`). Superdyn imports `SetServices`
from the Fortran module `FVdycoreCubed_GridComp` (the `DYN` component) and adds
`DYN` with `MAPL_GridCompAddChild(gc, "DYN", FV3_SetServices, "dyn.yaml")` when
`DYCORE:` (default `FV3`) is `FV3`; any other value fails
(`GEOS_SuperdynGridComp.F90`, `SetServices`, the `select case (trim(dycore))`). The `use` of
`AdvCore_GridCompMod` and the `FV3+ADV` and `DATMO` cases are commented out.

Each component has an external-procedure wrapper for loading from a shared
object: `FVdycoreCubed_SetServices` and `AdvCore_SetServices` (end of the
respective files); the regression cases name them in YAML
(`fvdycorecubed_setservices_`, `advcore_setservices_`, in
`libFVdycoreCubed_GridComp`).

## 2. Directory and file map

| Path | Content |
|---|---|
| `DynCore_GridCompMod.F90` | `DYN` component, about 5000 lines. Only `SetServices` is public |
| `DynCore_StateSpecs.rc`, `AdvCore_StateSpecs.rc` | Import, export and internal specs of `DYN` and `ADV` (`mapl_acg`) |
| `FV_StateMod.F90` | `FV_Setup`, `FV_InitState`, `FV_Run`, `FV_Finalize`, `State_To_FV`, `FV_To_State`; module-level array `FV_Atm` |
| `AdvCore_GridCompMod.F90` | `ADV` component |
| `fv_shared.F90` | Tracer-name helpers used by `DYN` and `ADV` |
| `rs_scale*.F90`, `interp_restarts.F90` | Restart tools built here (section 10) |
| `AppGridCreate.F90`, `CubeHalo.F90`, `GetWeights*.F90`, `fv_regrid*.F90`, `fv_regridding_utils.F90` | Grid creation and regridding helpers |
| `sw.f90`, `jw.f90`, `testcases_3_4_5_6_stand_alone.f90` | Idealized cases used by `Coldstart` |
| `GEOS_FV3_Utilities.F90` | Utility module (`a2d2c`) |
| `regression/` | CTest cases (section 12) |
| `.circleci/config.yml`, `.github/` | CI, `CODEOWNERS`, workflows |

Not listed in `CMakeLists.txt` (not built): `c2c_rst.F90`, `c2l_CFIO_offline.F90`,
`CreateTopography.F90`, `dyn_interp_rst.F90`, `moist_interp_rst.F90`,
`ESMF_GenerateCSGridDescription.F90`, the four `StandAlone_*.F90`, `.cppdefs`,
`FVdycore_arch.mk`. *Inferred:* legacy or manually built tools.

## 3. The DYN component

### 3.1 How it is called

Superdyn runs `DYN` as a child (`GEOS_SuperdynGridComp.F90`):

| Phase | Routine | Called by | Does |
|---|---|---|---|
| `run` | `run` | Superdyn `run`: `MAPL_GridCompRunChild(gc, "DYN", phase_name="run")` | Applies the analysis/bias tendencies, advances the FV3 state one dynamics step, fills exports |
| `run_add_incs` | `run_add_incs` | Superdyn `run_add_incs` | Adds the physics tendencies to the FV3 state, recomputes exports |

Superdyn `run` runs `ADV` only if its module variable `ADV` is not -1; it is
initialized to -1 and never assigned (`GEOS_SuperdynGridComp.F90`, module
variable `ADV` and `run`). Superdyn registers no imports and copies no tendencies into `DYN`:
`DYN` declares its own imports (3.2) and relies on connections to fill them.
`GEOS_AgcmGridComp.F90` (in `GEOSgcm_GridComp`, `GEOSagcm_GridComp/`) still calls
the Superdyn child with `ESMF_GridCompRun(..., PHASE=1)` in `Run`, the MAPL2
interface, and fills `DUDT`, `DVDT`, `DPEDT` with the `DO_UPDATE_ANA3D` calls
in `Run`. *Unverified:* what fills the `DYN` imports `D_UV_DT`, `DTDT`,
`DPEDT`, `DQ*ANA`, `DOXANA`, `PHIS`, `VARFLT` in the full GCM. In `adv-dyn`,
connections from data components fill `PHIS`, `VARFLT` and the six `DQ*ANA`
(`regression/adv-dyn/root.yaml`).

### 3.2 SetServices

(`DynCore_GridCompMod.F90`, `SetServices`)

- Specs come from `DynCore_StateSpecs.rc` through the generated includes
  `DynCore_Import___.h`, `DynCore_Export___.h`, `DynCore_Internal___.h`
  (`mapl_acg` in `CMakeLists.txt`).
- Imports (13 specs): `D_UV_DT` (vector wind tendency), `DWDT`, `DTDT`, `DPEDT`,
  `DQVANA`, `DQLANA`, `DQIANA`, `DQRANA`, `DQSANA`, `DQGANA`, `DOXANA`, `PHIS`,
  `VARFLT`. Only `VARFLT` sets `RESTART` (`SKIP`).
- `TRADV`: added in code as an export spec with `itemtype=MAPL_STATEITEM_SERVICE`;
  a service item appears in both the import and the export state (source
  comment). Section 6.
- Exports (206 rows in `DynCore_StateSpecs.rc`): `CX`, `CY`, `MFX`, `MFY` (R8),
  `PLE0`, `PLE1` (inputs of `ADV`), the state (`UV` vector, `T`, `Q`, `W`,
  `ZLE`, ...), grid and vertical coordinates, the `*_DYN_IN` set, tendency and
  energetics diagnostics.
- Internal state (section 5); private state `DynState` named `DYN_STATE`.
- Registers `Initialize`, `run`, `run_add_incs`, `Finalize`; the `Coldstart`
  entry point is commented out. `MAPL_GridCompSetGeometry` builds the grid from
  the config (section 8). Reads `FV3_STANDALONE:`, `DEBUG_DYN:`, `DEBUG_ADV:`,
  `DEBUG_TQ_ERRORS:`.

### 3.3 Initialize

(`DynCore_GridCompMod.F90`, `Initialize`)

1. `DynSetup` (alias for `FV_Setup`, section 4), timed `DYN_SETUP`.
2. Reads `DO_ADD_INCS:` and `COLDSTART:`. A nonzero `COLDSTART` calls `Coldstart`
   (idealized state from `CASE_*`, `T0:`, `IM:`, `JM:`).
3. `DynInit` (alias for `FV_InitState`): connects `STATE%VARS` to the MAPL
   internal-state arrays, creates the alarm "Time2Run" at `RUN_DT:` and a
   six-hour dry-mass alarm, sets `FV_Atm%phis`, `varflt`, `ak`, `bk`, calls
   `State_To_FV`, and under `RUN_GTFV3` `geos_gtfv3_interface_f_init`.
4. Fills the exports `PREF`, `UV` (if present), `T`, `DXC`, `DYC`, `AREA`.
5. Creates the sticky alarm "ReplayShutOff" (`REPLAY_SHUTOFF:`, default -3600;
   the absolute value in seconds is the ring interval).

There is no explicit `MAPL_GenericInitialize` call. *Inferred:* the framework
handles the internal-state restart.

### 3.4 run (phase `run`)

(`DynCore_GridCompMod.F90`, `run`)

1. Returns at once if the "Time2Run" alarm is not ringing (the dycore step may
   be longer than the parent's, `RUN_DT:`).
2. Gets `TRADV` from the import and export state and asserts the two bundles
   hold the same data (`MAPL_FieldBundleSameData`; a source comment says there
   is no coupling mechanism yet). If `EXCLUDE_ADVECTION_TRACERS` asks for it
   (`should_adjust_tracers`), the export bundle is replaced by an adjusted
   bundle `xTRADV` (section 6); `PULL_Q` points `STATE%VARS%tracer(:)` at the
   bundle fields.
3. The "WMP REPLAY/ANA" section runs when `FV3_STANDALONE` is false and the
   `ReplayShutOff` alarm is not ringing. Intermittent-replay code is not in this
   file (a comment: "we are not doing this anymore, right?"). If `DO_ADD_INCS`
   (default true) and not `ADIABATIC`, the imports `D_UV_DT`, `DPEDT`, `DTDT`
   are applied with `add_incs` (`ANA_IS_WEIGHTED:`, default `NO`, must be `YES`
   or `NO`, selects mass-weighted temperature tendencies). `DQVANA` etc. rebuild
   the pre-increment moisture state (`QVOLD = QV - DQVANA`) for diagnostics;
   `add_incs` does not change the tracers.
4. Fills the `*_DYN_IN` exports.
5. Calls `DynRun` (alias for `FV_Run`, section 4) with `PLE0` the edge pressure
   before the step.
6. Fills the exports (`PLE0`, `PLE1`, `CX`, `CY`, `MFX`, `MFY`, `UV`, `T`, `Q`,
   `PLE`, diagnostics), the mass fluxes from `fillMassFluxes` or, for
   `AdvCore_Advection == 2`, `computeMassFluxes`; then `free_tracers`.

### 3.5 run_add_incs (phase `run_add_incs`)

(`DynCore_GridCompMod.F90`, `run_add_incs`)

- Skipped if `SW_DYNAMICS`. Only `run` checks the Time2Run alarm.
- `PULL_Q`, then (if `DO_ADD_INCS` and not `ADIABATIC`) `add_incs` adds the
  physics tendencies (timer `PHYS_ADD_INCS`): `D_UV_DT` is converted from the A
  grid to the native D grid with `Agrid_To_Native` (`wind_increment_limiter =
  800.d0/86400.d0`, section 16); `PE` is advanced by `DT*DPEDT`; `PT` is
  re-weighted by old and new `DP` and `cp/cvm`; `getPKZ` recomputes `PKZ`. The
  tracers are not changed. The use of `DWDT` in `add_incs` is commented out.
- Recomputes the dependent exports (`U_DGRID`, `ZLE`, `PLE`, `PKE`, `KEPHY`,
  `PEPHY`, `TEPHY`, ...), then frees the tracer pointers.

### 3.6 Finalize

`DynFinalize` (alias for `FV_Finalize`: `fv_end`, and the gtFV3 finalize under
`RUN_GTFV3`).

## 4. FV_StateMod: the bridge to the FV3 core

(`FV_StateMod.F90`)

`FV_Atm(:)` is a module-level `save` array of the FV3 type `fv_atmos_type`
(`@fvdycore`, `fv_arrays_mod`). `AdvCore_GridCompMod` uses the same array
(`use FV_StateMod, only: FV_Atm`), so `DYN` and `ADV` share one FV3 state.

**`FV_Setup`** (called from `Initialize`): `fms_init`, `fv_init1`; takes the
grid size and layout from the component's config (`MAPL_GridCompGet` with
`hconfig`, `num_levels`; `setup_fv_dimensions_and_topology` reads
`geometry: esmf_geom:` `im_world`, `nx_face`, `ny_face`; `npz==1` selects
`SW_DYNAMICS`); reads `AGCM.STRETCH_FACTOR`, `FIXED_LATS`, `RUN_DT:`,
`DYCORE:` and sets `AdvCore_Advection` (0 for
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

Internal state of `DYN` (`DynCore_StateSpecs.rc`, `category: INTERNAL`):

| Name | Notes |
|---|---|
| `AK`, `BK` | Vertical-only, edge; `RESTART` required; `ADD2EXPORT` true |
| `U`, `V`, `PT`, `PE`, `PKZ` | D-grid winds, scaled potential temperature (K Pa^-kappa), edge pressure, `pk` exponent; R8; `RESTART` required. `PE` is also exported as `PLE` |
| `DZ`, `W` | Non-hydrostatic; R8; no `RESTART` entry |

`FV_InitState` wires the FV3 state arrays to these MAPL arrays
(`FV_StateMod.F90`), so once a restart is read the FV3 core sees the values.
The regression cases enable the internal checkpoint in the component YAML
(`mapl: misc: checkpoint: internal: true`). The dycore state holds no tracers:
the components that own them keep them, and `DYN` borrows pointers through
`TRADV`. `free_tracers` only deallocates and nullifies the FV tracer pointers
(a comment there says "Comment out to output tracer to checkpoint file").
`COLDSTART:` builds an idealized state instead of reading a restart; both
regression cases use it (`COLDSTART: true`, `CASE_ID: 1`).

## 6. Tracers: the `TRADV` bundle

`TRADV` is the dycore's bundle of all advected tracers, a service item in both
the import and export state of `DYN` (3.2). It is filled by connecting `TRADV`
of the owning component to `TRADV` of `DYN` (`regression/*/root.yaml`:
`DataMoist`, and for `adv-dyn` also `DataTracers` and `ADV`). Superdyn
re-exports `TRADV` from `DYN` (`MAPL_GridCompReexport`).

Selection in `DYN` (`run`, `should_adjust_tracers`, `get_adjusted_tracer_bundle`):

- `EXCLUDE_ADVECTION_TRACERS:` (default `ALWAYS`; `PREDICTOR` uses the clock
  alarm `PredictorAlarm`; other values are ignored) decides whether the bundle
  is adjusted. `EXCLUDE_ADVECTION_TRACERS_LIST` is a YAML sequence in the
  component config (`ESMF_HConfigAsStringSeq`; `[]` in the regression cases).
- When `AdvCore_Advection >= 1`, every tracer for which
  `field_is_cloud_water_species` (`fv_shared.F90`) is false is excluded. It is
  true for `Q`, `QLCN`, `QLLS`, `QICN`, `QILS`, `CLCN`, `CLLS`, `NCPL`, `NCPI`,
  `QRAIN`, `QSNOW`, `QGRAUPEL`. Names come from the field `standard_name`
  (`get_short_name`).
- The remaining fields go into a saved module bundle `bundle_adv` (`xTRADV`),
  and `PULL_Q` points `STATE%VARS%tracer(:)` at them.

Packing in `FV_Run` (`FV_StateMod.F90`): `QLCN`/`QLLS` are summed into the
liquid slot, `QICN`/`QILS` into ice, `CLCN`/`CLLS` into cloud fraction; the
convective parts are carried as extra tracers (`idx_qlcn`, `idx_qicn`,
`idx_clcn`); `QRAIN`, `QSNOW`, `QGRAUPEL` have their own slots; other tracers
take the next free index; the final count must equal `ncnst` (asserted). After
`fv_dynamics` the large-scale part is total minus convective, bounded below by
zero.

With `DYCORE: FV3` every tracer in `TRADV` (minus the exclude list) is advected
inside `fv_dynamics`; with `FV3+ADV` only the species above stay there and `ADV`
handles the rest (section 7).

## 7. The ADV component

(`AdvCore_GridCompMod.F90`, `AdvCore_StateSpecs.rc`)

Offline tracer advection driven by the mass fluxes of `DYN`. In
`GEOS_SuperdynGridComp.F90` the `ADV` `MAPL_GridCompAddChild` in `SetServices` is commented
out, but the `if (adjustl(DYCORE) == "FV3+ADV")` block that connects `MFX`,
`MFY`, `CX`, `CY`, `PLE0`, `PLE1` to `ADV` (in `SetServices`) is live and unreachable,
since any `DYCORE` other than `FV3` hits `_FAIL`. `ADV` is added by the
regression case `adv-dyn` (`root.yaml`, child `ADV` with `adv.yaml`).

- Imports (R8): `MFX`, `MFY`, `CX`, `CY`, `PLE0`, `PLE1`; `TRADV` is a service
  item as in `DYN`. Exports: `AREA` and 38 `TEST_TRACER*` diagnostics
  (`ntracers = 38`).
- `SetServices` registers `Initialize`, phase `run`, `Finalize`; reads `RUN_DT:`,
  `DYCORE:` and `AdvCore_Advection:` (`FV3` sets 0; `FV3+ADV` and any other
  value 1). If `DYCORE` is neither, the dycore is "not running" and
  `Initialize` calls `fv_setup` and `MAPL_GridCompSetGeometry` itself.
- `run` (only if `AdvCore_Advection > 0`) copies the imports, builds `xTRADV`
  with its own `get_adjusted_tracer_bundle`, reads `ADV_QSPLIT:` (default 0),
  `ADV_CORE_REPORT_TRACER_MASS:`, calls `offline_tracer_advection(...)` and sets
  values below `tiny` to zero.

`DYN` and `ADV` each have a function `get_adjusted_tracer_bundle`; both use
`field_is_cloud_water_species`, with opposite polarity (`DYN` excludes tracers
for which it is false, `ADV` those for which it is true). The species list is
defined once, in `fv_shared.F90`.

## 8. Configuration

Resource keys are read with `MAPL_GridCompGetResource` from the component's
configuration; in the regression cases that is its YAML file
(`regression/adv-dyn/dyn.yaml`, `regression/dyn-sa/dyn-sa.yaml`).

- Top-level keys: `DYCORE:`, `AdvCore_Advection:`; `RUN_DT:`; `FV3_CONFIG:`,
  `ADIABATIC:`, `FV_OFF:`, `ADJUST_SPLITS:`, `RUN_GTFV3:`, `fix_mass:`;
  `COLDSTART:`, `CASE_ID:`, `CASE_TRACERS:`; `FV3_STANDALONE:`; `DO_ADD_INCS:`,
  `ANA_IS_WEIGHTED:`, `REPLAY_SHUTOFF:`; `EXCLUDE_ADVECTION_TRACERS[_LIST]:`.
- `geometry:` with `esmf_geom:` (`class: CubedSphere`, `im_world`, `nx_face`,
  `ny_face`) and `vertical_grid:` (`num_levels`, `field_edge: PLE`,
  `field_center: PE`). `FV_Setup` reads these
  (`setup_fv_dimensions_and_topology`); `DYN` does not call `MAPL_GridCreate`.
- `mapl:` with `setServices` (`sharedObj: libFVdycoreCubed_GridComp`,
  `userRoutine: fvdycorecubed_setservices_`) and `misc` (checkpoint options).

Superdyn names the config file `dyn.yaml`; *Unverified:* the file used in the
full GCM was not found in the fixture.

## 9. `@fvdycore`

`@fvdycore` is the nested repository holding the GFDL atmospheric cubed-sphere
core with GEOS modifications (`components.yaml`, key `fvdycore`, branch `geos/release/MAPL-v3`, develop `geos/develop`).
Its `CMakeLists.txt` builds a library `fvdycore` from `model/`, `tools/` and
`geos_utils/`, added here as a sub-component. It is not part of this repository
(`.gitignore` ignores `/@fvdycore`) and is not documented here; see
`@fvdycore/README.md` (`../@fvdycore/README.md` in the fixture) (the upstream GFDL README).

## 10. Executables and tools

Built by `CMakeLists.txt`: `rs_scale.x` (inserts a dry mass into restarts),
`interp_restarts.x` (interpolates FV and moist restarts to a cubed sphere, with
optional new vertical levels), and `fv3_driver.x` (`geos-gtfv3/driver/`, only
with `BUILD_GEOS_GTFV3_INTERFACE`). The four `StandAlone_*.F90` sources remain
but have no `ecbuild_add_executable`; the standalone use of `DYN` is the
regression case `dyn-sa`.

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

## 12. Regression tests and CI

`regression/CMakeLists.txt` reads `regression/cases.txt` (`adv-dyn`, `dyn-sa`)
and calls `esma_add_regression_tests(NAME fv3 ...)`. Each case is a CTest with
label `REGRESSION`. `regression/run_case.cmake` runs `GEOS.x` on 6 ranks and
compares checkpoints against a baseline only when the compiler is `IntelLLVM`
(`regression/README.md`: `ctest -L REGRESSION`).

| Case | Content |
|---|---|
| `dyn-sa` | Children `DYNsa` (`dyn-sa.yaml`: `FV3_STANDALONE: true`, `COLDSTART: true`, `DYCORE: FV3`, `ADIABATIC: true`, 12 cells per face edge, 91 levels) and `DataMoist`; `TRADV` connected |
| `adv-dyn` | Children `DYN` (`DYCORE: FV3+ADV`), `ADV`, and data components for moist (`DataMoist`), extra tracers (`DataTracers`), `PHIS`/`VARFLT` (`DataPhisVarflt`) and `DQ*ANA` (`DataAna`); `MFX`, `MFY`, `CX`, `CY`, `PLE0`, `PLE1` connected from `DYN` to `ADV` |

The README says each case has a `cmpchk.py`; none is in the repository.
CI (`.circleci/config.yml`): `build-GEOSfvdycore-as-R4-on-<compiler>` (`ifx`,
`gfortran`; fixture branch `release/MAPL-v3`; `ctest -L REGRESSION`) and
`build-GEOSfvdycore-as-R8-on-ifx`.

## 13. `scripts/`

`scripts/CMakeLists.txt` installs `fv3_setup` (tcsh, interactive), `fv3.j` and
`.FV3_VERSION`. `fv3.j` writes MAPL2-format `CAP.rc`, `AGCM.rc`, `HISTORY.rc`
and runs `StandAlone_FV3_Dycore.x`, which is not built here (section 10).

## 14. Build

`CMakeLists.txt` builds one library `FVdycoreCubed_GridComp` (`esma_add_library`,
`SUBCOMPONENTS fvdycore`; `TYPE SHARED` unless `BUILD_GEOS_GTFV3_INTERFACE`), the
executables of section 10, `scripts` and `regression`. `mapl_acg` runs on
`DynCore_StateSpecs.rc` (`IMPORT_SPECS EXPORT_SPECS INTERNAL_SPECS`) and
`AdvCore_StateSpecs.rc` (`IMPORT_SPECS EXPORT_SPECS GET_POINTERS
DECLARE_POINTERS`). Dependencies: `MAPL`, `GFTL_SHARED::gftl-shared-v2`,
`GMAO_hermes` (`m_set_eta`, `m_topo_remap`), `GEOS_Shared`, `esmf`, OpenMP.
`FV_PRECISION` (`R4`, `R4R8`, `R8`; set in the fixture, default `R4R8`)
selects the FMS library (`fms_r4`, `fms_r8`) and the flags `SINGLE_FV
OVERLOAD_R4` or `FREAL8`. Other options:
`BUILD_GEOS_GTFV3_INTERFACE`, `CRAY_POINTER`.

## 15. Known traps

- `FV_Run` asserts exactly five water tracers when not `ADIABATIC` and requires
  `Q`, `QLCN`, `QLLS`, `QICN`, `QILS` in `TRADV`.
- Only `run` checks the Time2Run alarm: when it is not ringing `run` returns
  early and does not refresh the exports, while `run_add_incs` runs on every
  call. *Inferred:* `RUN_DT:` is meant to be a multiple of the parent step.
- `run` asserts the import and export `TRADV` bundles share data.

## 16. Observations (noted, not fixed)

- `run` reads `CONSV:` and `FILL:` into `LCONSV` and `LFILL`, which are not used
  afterwards. `FV_Setup` sets the FV3 flags itself.
- The module header comment (MAPL2-era) says the component "supports only a 1D
  latitude-based decomposition" and mentions PILGRIM; the comment block above
  `SetServices` still lists the imports `DUDT`, `DVDT`.
- The body of `FV_DA_Incs` is commented out. In `a2d3d` (`FV_StateMod.F90`) the
  code using the optional argument `wind_increment_limiter` is commented out,
  so the value `800.d0/86400.d0` passed by `add_incs` has no effect; no other use
  was found.
- The export `PID` is set to 0 ("WMP need to get from MAPL gid").
- The gtFV3 path passes only seven tracers (section 11) while `FV_Atm%q` holds
  all packed tracers.
- `setup_fv_dimensions_and_topology` sets `npy = 6*npx`, so its `else` branch is
  never taken (a `TODO` in the source asks about it).
- `FV3_STANDALONE` is read in `SetServices` and not used there; `run` and
  `Coldstart` read it again.
