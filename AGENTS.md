# AGENTS.md

**Audience:** coding agents and developers working in the `FVdycoreCubed_GridComp` repository.
**Describes:** the MAPL3 line, branch `release/MAPL-v3`.
**Verified against:** `FVdycoreCubed_GridComp` @ `dadf5b5` (2026-11-06).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `f50622d`).
Other fixtures (for example GEOSadas, GEOSldas) may place this repository
elsewhere; paths inside this repository do not change.

## What this is

The MAPL gridded-component wrapper of the GFDL FV3 cubed-sphere dynamical core,
at `src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSsuperdyn_GridComp/@FVdycoreCubed_GridComp`
in the GEOSgcm fixture (`components.yaml`, key `FVdycoreCubed_GridComp`). Parent:
`GEOS_SuperdynGridComp`. The core is the nested repository `@fvdycore`; edit it
only there.

> **Build status on this line:** built and ported to MAPL3. Superdyn's
> `CMakeLists.txt` builds it as the only `SUBCOMPONENT`; Superdyn adds `DYN`
> (config `dyn.yaml`) and runs the phases `run` and `run_add_incs`. `ADV` is
> compiled but not added by Superdyn; the regression case `adv-dyn` runs it.

## Documents

| Document | Content |
|---|---|
| [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) | Read first: `DYN`, `FV_StateMod`, state, tracers, `ADV`, configuration, build, regression, traps |
| [`docs/decisions.md`](docs/decisions.md) | Previous architectural decisions |
| `../../ARCHITECTURE.md` | Parent: Agcm and Superdyn |
| `regression/README.md` | CTest cases `adv-dyn`, `dyn-sa` |

## Where things are

| Want to change | Go to |
|---|---|
| Imports, exports, internal state | `DynCore_StateSpecs.rc` (`mapl_acg`); `AdvCore_StateSpecs.rc` |
| Run phases, tracer selection, tendencies | `DynCore_GridCompMod.F90` (`run`, `run_add_incs`, `add_incs`) |
| FV3 call, tracer packing, grid size | `FV_StateMod.F90` (`FV_Run`, `FV_Setup`) |
| Offline advection; shared tracer predicate | `AdvCore_GridCompMod.F90`; `fv_shared.F90` |
| Python gtFV3 core; regression; standalone scripts | `geos-gtfv3/`; `regression/`; `scripts/` |

## Rules for changes

- Superdyn runs `DYN` phases `run` and `run_add_incs` and copies nothing; `DYN`
  reads its own imports (`D_UV_DT`, `DTDT`, `DPEDT`, ...). `run` applies the
  analysis/bias tendencies, `run_add_incs` the physics tendencies, both through
  `add_incs`.
- Specs are edited in the `.rc` files, not by `MAPL_Add*Spec` calls.
- Keep the `ADV` imports (`MFX`, `MFY`, `CX`, `CY`, `PLE0`, `PLE1`) consistent
  with `AdvCore_StateSpecs.rc` and the connections that fill them: the Superdyn
  `FV3+ADV` connection block (live but unreachable) and `regression/adv-dyn/root.yaml`.
- `DYN` and `ADV` use one predicate, `field_is_cloud_water_species`
  (`fv_shared.F90`), with opposite polarity.
- Tracers are not in the dycore restart; `DYN` borrows the `TRADV` fields.
- `DYN` and `ADV` share the module array `FV_Atm`.
- Files not in `CMakeLists.txt` are not compiled (architecture document,
  section 2). `FV_PRECISION` is set by the fixture build.
- Error handling in code you write or touch uses the MAPL3 macros (`@MAPL/include/MAPL_ErrLog.h`): `_RC` on calls, `_STAT` on `allocate`, `_ASSERT(cond, msg)` with a message, and `_RETURN(_SUCCESS)` at every successful exit. Convert `RC=STATUS`/`VERIFY_`/`ASSERT_`/`RETURN_` in any routine you change. `_RC` expands to two statements (`rc=status);_VERIFY(status`), so it cannot end a one-line `if`; use a block `if`. Do not add decorative `!----` separator lines.

## Build and test

Built through Superdyn. Tests: CTest cases `adv-dyn` and `dyn-sa` (label
`REGRESSION`) and the CI builds in `.circleci/config.yml`. Executables:
`rs_scale.x`, `interp_restarts.x`; the `StandAlone_*.F90` drivers are not built.
