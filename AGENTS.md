# FVdycoreCubed_GridComp: guide for agents

**Describes:** the MAPL3 line, branch `release/MAPL-v3`; verified against `FVdycoreCubed_GridComp` @ `4d7c4a3` (2026-11-10).
**Placement:** `src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSsuperdyn_GridComp/@FVdycoreCubed_GridComp` in the GEOSgcm fixture (`components.yaml`); other fixtures (GEOSfvdycore) may place it elsewhere. Paths below are relative to this repository.

> **Built and ported on this line.** Superdyn builds this repository as its only subcomponent and adds `DYN` (`dyn.yaml`). `ADV` is compiled but Superdyn does not add it; the regression cases `adv-dyn` and `adv-sa` run it. The FV3 core `@fvdycore` is a separate nested repository; edit it only there. The Agcm that runs Superdyn is not compiled on this line, so `DYN` runs only in the regression cases.

## Rules for agents

- Read [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) (what `DYN` and `ADV` are and how they fit: grid, state, run phases, tracers, contracts, invariants) before changing phases, tracer handling, exports or configuration. Update it when the structure, a contract or an invariant changes; it describes only what changes rarely, so keep status and to-dos out of it.
- Note bugs and oddities; do not fix them unless asked. Record a bug or an open question as a GitHub issue, not in a document; a document may link to the issue.
- Do not edit the generated `*___.h` files; edit `DynCore_StateSpecs.rc` or `AdvCore_StateSpecs.rc`.
- When porting code to MAPL3, follow "Porting to MAPL3" in the `GEOSgcm_GridComp` `AGENTS.md` (`../../../AGENTS.md` in the fixture).
- Keep this file current: update it when status, a rule or a trap changes, and update its "verified against" commit.
- Cite files and routines, never line numbers. Label unchecked statements *Inferred:* or *Unverified:*.
- Error handling in code you write or touch uses the MAPL3 macros (`@MAPL/include/MAPL_ErrLog.h`): `_RC` on calls, `_STAT` on `allocate`, `_ASSERT(cond, msg)` with a message, and `_RETURN(_SUCCESS)` at every successful exit. Convert `RC=STATUS`/`VERIFY_`/`ASSERT_`/`RETURN_` in any routine you change. `_RC` expands to two statements (`rc=status);_VERIFY(status`), so it cannot end a one-line `if`; use a block `if`. Do not add decorative `!----` separator lines.

## Status and to-dos

- Superdyn accepts only `DYCORE: FV3`; its `FV3+ADV` and `DATMO` cases are commented out, but its block connecting `DYN`'s mass fluxes to `ADV` is still there, unreachable.
- Every regression case is `ADIABATIC: true` and cold-started (`COLDSTART: true`, `CASE_ID: 1`): the water-tracer packing in `FV_Run`, `add_incs` and restart reading are not exercised.
- `run` asserts that the import and export `TRADV` bundles share data (`MAPL_FieldBundleSameData`), because MAPL3 has no coupling for a service bundle yet.
- Not ported or not built: `Coldstart` as an entry point (it is called from `Initialize` instead), the `StandAlone_*.F90` drivers and several other tools (`CMakeLists.txt` lists what is built); `scripts/fv3.j` still writes MAPL2 `.rc` files for the unbuilt `StandAlone_FV3_Dycore.x`.

## Known traps

- **Five water species are required.** Unless `ADIABATIC`, `FV_Run` asserts that the tracers it receives include each of `Q`, `QLCN`, `QLLS`, `QICN` and `QILS` exactly once, and the packing then asserts that every tracer found a slot.
- **Tracers are found by name, but bundle fields have MAPL-generated names.** Use `get_short_name` (`standard_name`), never the field name.
- **`run` may return early.** If `RUN_DT` is longer than the parent step, `run` does nothing on most calls and leaves its exports unchanged, while `run_add_incs` adds tendencies on every call.
- **`ADV`'s behaviour depends on its own `DYCORE`.** With `FV3+ADV` it advects using the FV3 set-up of `DYN` (`adv-dyn`); with `DYCORE` unset it sets FV3 and its grid up itself (`adv-sa`); with `FV3` it does nothing (`AdvCore_Advection` is 0). Set it to match whether `DYN` is present and, when it is, to the same value as `DYN`'s.
- **The tracer split has two sides.** A change to `field_is_cloud_water_species` changes both what `DYN` keeps and what `ADV` advects.
- **`EXCLUDE_ADVECTION_TRACERS_LIST` must be a YAML sequence** (`[]` for none), read with `ESMF_HConfigAsStringSeq`; the adjusted bundle is built once, on the first `run`, and reused.
- **An unknown `EXCLUDE_ADVECTION_TRACERS` value** (other than `ALWAYS` or `PREDICTOR`) is meant to be ignored, but the code logs through an unassociated `logger` in `should_adjust_tracers`.
- **The Python core gets only seven tracers** (`q(:,:,:,1:7)` in `FV_Run`).
- **`FV_PRECISION` changes the build** (`R4`, `R4R8`, `R8`): FMS library, FV3 real kind and flags differ. CI runs the regression tests only on the R4 builds.

## Build and test

- Built as part of a fixture (GEOSgcm, or GEOSfvdycore) through Superdyn; this repository has no standalone build. Executables: `rs_scale.x`, `interp_restarts.x`, and `fv3_driver.x` with `BUILD_GEOS_GTFV3_INTERFACE`.
- Regression tests (label `REGRESSION`): `regression/`, cases `adv-dyn`, `adv-sa`, `dyn-sa` (`cases.txt`); `regression/README.md` describes them. Each case directory has a `logging.yaml` (keep it in new cases). A case runs even when its baseline is missing; only the comparison is then skipped, and the comparison is done only for IntelLLVM (ifx) builds (`run_case.cmake`). `adv-sa` reads `adv-dyn`'s baseline exports as its input restart.
- Pull requests for this line target `release/MAPL-v3`, not `develop` (the MAPL2 line), and need a 0-diff label (`.github/workflows/enforce-labels.yml`). CI (`.circleci/config.yml`) builds the GEOSfvdycore fixture at `release/MAPL-v3` as R4 with ifx and gfortran (running the regression tests) and as R8 with ifx. Code owners: the FV3 gatekeepers, and the CMake team for `CMakeLists.txt` (`.github/CODEOWNERS`).
