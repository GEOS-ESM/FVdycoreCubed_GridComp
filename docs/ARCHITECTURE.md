# FVdycoreCubed_GridComp architecture

A map of the repository: what lives where, how the pieces fit, and what must not break. Paths are relative to the repository root. It names modules and routines rather than linking to lines; use them to search the code. It describes only things that change rarely; status, rules for changes, known traps and to-dos are in [`AGENTS.md`](../AGENTS.md).

## Bird's-eye view

This repository wraps the GFDL FV3 cubed-sphere dynamical core as MAPL components. In GEOS its main component is `DYN`, the child of Superdyn (`GEOS_SuperdynGridComp`, in `GEOSgcm_GridComp`); the Agcm runs Superdyn twice per step, once to advance the resolved flow (with analysis increments) and once to add the physics tendencies. The core itself is the nested repository `@fvdycore`; this repository holds the MAPL layer (`DYN`), the bridge to the core's data structures (`FV_StateMod`), an offline tracer-advection component (`ADV`), restart tools and an optional Python core.

The defining design choice: the dynamics state lives in MAPL's internal state, which `FV_StateMod` copies into the FV3 core before each step and back after it, while the tracers stay with the components that own them and are lent to the core through the `TRADV` service.

## Codemap

- `DynCore_GridCompMod.F90`, `DynCore_StateSpecs.rc`: `DYN`. `SetServices`, `Initialize`, phases `run` and `run_add_incs`, `Finalize`; `add_incs` applies tendencies; `Coldstart` builds an idealized state; tracer handling (`should_adjust_tracers`, `get_adjusted_tracer_bundle`, `PULL_Q`).
- `FV_StateMod.F90`: the bridge to FV3. `FV_Setup`, `FV_InitState`, `FV_Run`, `FV_Finalize`, `State_To_FV`, `FV_To_State`; the module array `FV_Atm` (the FV3 state).
- `AdvCore_GridCompMod.F90`, `AdvCore_StateSpecs.rc`: `ADV`, offline tracer advection from `DYN`'s mass fluxes (`offline_tracer_advection` in `@fvdycore`).
- `fv_shared.F90`: helpers both components use: `field_is_cloud_water_species`, `get_short_name`, `is_name_in_list`.
- `@fvdycore/`: the FV3 core (GFDL_atmos_cubed_sphere, branch `geos/release/MAPL-v3`), a separate nested repository built into this library.
- `geos-gtfv3/`: optional bridge to a Python FV3 core (built with `BUILD_GEOS_GTFV3_INTERFACE`), and its standalone driver `fv3_driver.x`.
- `rs_scale.F90`, `rs_scaleMod.F90`, `interp_restarts.F90` (with the grid and regridding helpers `AppGridCreate.F90`, `CubeHalo.F90`, `GetWeights*.F90`, `fv_regrid*.F90`): the restart tools `rs_scale.x` and `interp_restarts.x`.
- `sw.f90`, `jw.f90`, `testcases_3_4_5_6_stand_alone.f90`: idealized initial states for `Coldstart`.
- `regression/`: standalone regression cases. `scripts/`: the MAPL2-era `fv3_setup` and `fv3.j`.

## How it is put together

- **Placement.** Superdyn adds `DYN` with configuration `dyn.yaml` when `DYCORE` is `FV3` (any other value stops). It runs `DYN` phase `run` from its own `run` and phase `run_add_incs` from its own `run_add_incs`, and re-exports `DYN`'s state and `TRADV`. Each component also has an external `SetServices` wrapper (`fvdycorecubed_setservices_`, `advcore_setservices_`) for loading by name from `libFVdycoreCubed_GridComp`.
- **Grid.** `DYN` takes its cubed-sphere grid and levels from its configuration (`geometry:` with `esmf_geom: class: CubedSphere`, `im_world`, `nx_face`, `ny_face`, and `vertical_grid: num_levels`): `MAPL_GridCompSetGeometry` builds the ESMF grid and `FV_Setup` sizes FV3 from the same keys. One level means shallow-water dynamics.
- **State.** The internal state (`AK`, `BK`, D-grid `U`, `V`, `PT`, `PE`, `PKZ`, and non-hydrostatic `W`, `DZ`) holds FV3's prognostic state: `FV_InitState` points the bridge's `STATE%VARS` at it, `State_To_FV` copies it into `FV_Atm` before each step and `FV_To_State` copies it back. The MAPL internal restart is therefore the dycore restart.
- **`run`.** Returns unless its own alarm `Time2Run` (`RUN_DT`) rings. Otherwise it lends the tracers (below), applies the analysis or bias increments `D_UV_DT`, `DTDT`, `DPEDT` through `add_incs` (not with `FV3_STANDALONE`, after the `ReplayShutOff` alarm rings, with `ADIABATIC`, or without `DO_ADD_INCS`), calls `FV_Run`, and fills the exports, including the mass fluxes and Courant numbers `MFX`, `MFY`, `CX`, `CY` and the edge pressures `PLE0`, `PLE1` that `ADV` uses.
- **`run_add_incs`.** On every call (except for shallow-water dynamics): adds the physics tendencies through `add_incs` (winds converted from the A grid to the D grid, `PE`, `PT` and `PKZ` updated; tracers untouched), then recomputes the dependent exports.
- **`FV_Run`.** Packs the tracers into `FV_Atm%q`, copies the state in (`State_To_FV`), runs `fv_dynamics` (or the Python core), unpacks the tracers and copies the state back (`FV_To_State`). Water species are packed by name: liquid, ice and cloud fraction are each the sum of their convective and large-scale parts, with the convective part carried as an extra tracer and the large-scale part recovered after the step.
- **Tracers.** `TRADV` is a service `DYN` provides: each owner subscribes the tracers it wants advected, and `DYN` finds them in the bundle, by name through `standard_name` (`get_short_name`). `DYN` keeps no tracers of its own and frees its pointers after each phase.
- **`ADV`.** With `DYCORE: FV3+ADV`, `DYN` keeps only the cloud and water species (`field_is_cloud_water_species`), and `ADV`, which also provides `TRADV`, advects the rest with `DYN`'s mass fluxes. When `ADV`'s configuration has no `DYCORE`, it sets up FV3 and its grid itself, and its mass-flux imports come from elsewhere (a data component in `adv-sa`); with `DYCORE: FV3` it does nothing. Superdyn does not add `ADV` on this line.

## Contracts

### Imports of `DYN`

- From the Agcm (MAPL2 code, not compiled on this line): analysis or bias increments `D_UV_DT` (a vector), `DTDT`, `DPEDT`, and the `DQ*ANA` and `DOXANA` moisture and ozone increments (diagnostics only); the surface geopotential `PHIS` and the variance of the filtered topography `VARFLT`. `DWDT` is declared but not applied.
- Through `TRADV`: every advected tracer, from the components that own it (Moist, Chemistry, ...).

### Exports of `DYN`

- The dynamics state that Physics imports through Superdyn: `UV`, `T`, `PLE`, `ZLE`, `PREF`, `DELP`, `OMEGA`, `PS`, `AK`, `BK`, ...; pre-dynamics (`*_DYN_IN`) and tendency diagnostics; energetics.
- For `ADV`: `MFX`, `MFY`, `CX`, `CY`, `PLE0`, `PLE1` (`R8`).

### Services

- `DYN` provides `TRADV` (advection); `ADV` provides `TRADV` too. Neither subscribes to any service.

## Invariants

- Tracers are never in the dycore restart; their owners keep them.
- `DYN` and `ADV` share one FV3 state (`FV_Atm` is module data in `FV_StateMod`), and one definition of which tracers stay in the core (`field_is_cloud_water_species`, used with opposite polarity).
- Physics tendencies reach the dynamics state only in `run_add_incs`; analysis increments only in `run`.
- `run` works only when `Time2Run` rings; `run_add_incs` works on every call.

## Cross-cutting concerns

- **Generated specs:** `DYN` and `ADV` declare their specs in `*_StateSpecs.rc` (ACG). `TRADV` and `ADV`'s `TEST_TRACER*` diagnostics are added in `SetServices`.
- **Precision:** `FV_PRECISION` (`R4`, `R4R8`, `R8`), set by the fixture build, chooses the FMS library and the FV3 real kind; the internal state is `R8`.
- **Configuration:** resource keys are read from the component's configuration (`dyn.yaml` under Superdyn); the main ones are `DYCORE`, `AdvCore_Advection`, `RUN_DT`, `FV3_CONFIG`, `ADIABATIC`, `COLDSTART` and `CASE_*`, `FV3_STANDALONE`, `DO_ADD_INCS`, `ANA_IS_WEIGHTED`, `REPLAY_SHUTOFF`, `EXCLUDE_ADVECTION_TRACERS[_LIST]`, `RUN_GTFV3`.
- **Testing:** each regression case runs `DYN`, `ADV` or both for one hour (two or three steps) under a configurable Root, with data components that own the tracers (and, for `adv-sa`, provide the mass fluxes), and compares the checkpoints exactly against a stored baseline.

## Where to go next

- [`AGENTS.md`](../AGENTS.md): rules, status on this line, to-dos, known traps, build and test.
- Superdyn and the Agcm: `GEOSagcm_GridComp/AGENTS.md` and `GEOSagcm_GridComp/docs/ARCHITECTURE.md` in `GEOSgcm_GridComp` (`../../AGENTS.md` and `../../docs/ARCHITECTURE.md` in the fixture).
- `regression/README.md`: the regression cases.
- `@fvdycore/README.md`: the FV3 core (upstream GFDL).
