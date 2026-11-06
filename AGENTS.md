# AGENTS.md

**Audience:** coding agents and developers working in the `FVdycoreCubed_GridComp` repository.
**Describes:** the MAPL2 line, branch `develop`.
**Verified against:** `FVdycoreCubed_GridComp` @ `4deed15` (2026-10-27).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `0f2e0fb`).
Other fixtures (for example GEOSadas, GEOSldas) may place this repository
elsewhere; paths inside this repository do not change.

## What this is

The MAPL gridded-component wrapper of the GFDL FV3 cubed-sphere dynamical core.
In the GEOSgcm fixture it is at
`src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSsuperdyn_GridComp/@FVdycoreCubed_GridComp`
(`components.yaml`, key `FVdycoreCubed_GridComp`). Its parent is
`GEOS_SuperdynGridComp` in `GEOSgcm_GridComp`. The core itself lives in the
nested repository `@fvdycore` (key `fvdycore`); edit it only in its own repository.

## Documents

Read `docs/ARCHITECTURE.md` first.

| Document | Content |
|---|---|
| [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) | The `DYN` component (SetServices, Initialize, Run phases 1 and 2), `FV_StateMod`, state and restarts, tracers (`TRADV`), `ADV`, configuration keys, `@fvdycore`, `geos-gtfv3`, `scripts`, build, traps |
| [`docs/decisions.md`](docs/decisions.md) | Previous architectural decisions (pattern, likely reason) |
| `../../ARCHITECTURE.md` | Parent: Agcm and Superdyn run order (in `GEOSgcm_GridComp`) |
| `README.md`, `CONTRIBUTING.md` | Project basics, contribution rules |

## Where things are

| Want to change | Go to |
|---|---|
| Component interface (imports, exports, internal state, run phases) | `DynCore_GridCompMod.F90` (module `FVdycoreCubed_GridComp`) |
| Call into the FV3 core, tracer packing, dry mass | `FV_StateMod.F90` (`FV_Run`, `FV_Setup`) |
| Offline tracer advection (`FV3+ADV`) | `AdvCore_GridCompMod.F90` |
| Python gtFV3 core | `geos-gtfv3/` |
| Standalone experiment setup, build options | `scripts/`, `CMakeLists.txt` |

## Rules for changes

- `DYN` has two Run phases; Superdyn copies `DUDT`, `DVDT`, `DTDT`, `DWDT` into the
  `DYN` imports before each. Phase 1 applies the analysis/bias tendencies the Agcm
  put there (`DO_UPDATE_ANA3D`); phase 2 applies the physics tendencies
  (`DO_UPDATE_PHY`). Both go through `ADD_INCS`.
- Keep exports used by `ADV` (`MFX`, `MFY`, `CX`, `CY`, `PLE0`, `PLE1`) consistent
  with `AdvCore_GridCompMod.F90` and Superdyn.
- *Inferred:* the tracer selection lists in `DYN` and `ADV` are separate and must
  agree; no code checks it (`docs/ARCHITECTURE.md` section 7).
- Tracers are not in the dycore restart; `DYN` borrows the `TRADV` bundle fields.
- `DYN` and `ADV` share the module array `FV_Atm` of `FV_StateMod`.
- Several `.F90` files are not in `CMakeLists.txt` (architecture document,
  section 2); they are not compiled here. `FV_PRECISION` is set by the fixture build.

## Build and test

Built as part of the fixture (GEOSgcm). Tests: this repository has no unit-test
directory; its checks are the CI jobs in `.circleci/config.yml` (`ci/build`,
`ci/run_fv3`, `ci/run_gcm`). Executables are listed in `docs/ARCHITECTURE.md`
section 10; standalone experiments: `scripts/fv3_setup`.
