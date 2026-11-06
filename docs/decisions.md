# FVdycoreCubed_GridComp: previous architectural decisions

**Audience:** coding agents and developers working in the `FVdycoreCubed_GridComp` repository.
**Describes:** the MAPL3 line, branch `release/MAPL-v3`.
**Verified against:** `FVdycoreCubed_GridComp` @ `dadf5b5` (2026-11-06).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `f50622d`).
Other fixtures (for example GEOSadas, GEOSldas) may place this repository
elsewhere; paths inside this repository do not change.

Only what is easy to see in the code now. Entries grow over time. Context:
[`ARCHITECTURE.md`](ARCHITECTURE.md).

### Dycore state aliases the MAPL internal state
- Pattern: `FV_InitState` points the FV3 state arrays at the MAPL internal-state arrays
  (`AK`, `BK`, `U`, `V`, `PT`, `PE`, `PKZ`, `DZ`, `W`) so a MAPL restart fills the dycore
  directly (`FV_StateMod.F90`, `FV_InitState`; `DynCore_GridCompMod.F90`, `SetServices`).  [Verified]
- Likely reason: reuse MAPL restart and checkpoint handling instead of the FV3 own I/O.  [Inferred]
- Question for maintainers: what is the role of the FV3 native restart I/O (`fv_io`) in the GEOS configuration?

### Tracers are borrowed, not owned
- Pattern: the dycore restart holds no tracers; `DYN` imports the `TRADV` bundle and points
  its tracer pointers at the bundle fields, then `FV_Run` packs and unpacks them around
  `fv_dynamics` (`DynCore_GridCompMod.F90`, `run`, `PULL_Q`; `FV_StateMod.F90`, `FV_Run`).  [Verified]
- Likely reason: tracers belong to the components that produce them and carry their own restarts.  [Inferred]
- Question for maintainers: is the `free_tracers` comment ("Comment out to output tracer to
  checkpoint file") a supported option?

### Two Run phases: dynamics, then physics increments
- Pattern: phase `run` applies the analysis/bias tendencies found in the imports
  and advances the state; phase `run_add_incs` adds the physics tendencies to the dycore state
  through the same `add_incs`, instead of the physics changing the state itself
  (`DynCore_GridCompMod.F90`, `SetServices`, `run`, `run_add_incs`).  [Verified]
- Likely reason: physics runs between the two phases on the dynamics state, and the dycore
  holds the prognostic variables on its native D-grid winds and edge pressure.  [Inferred]
- Question for maintainers: on this line Superdyn runs the two phases but copies nothing
  (`GEOS_SuperdynGridComp.F90`, `run`, `run_add_incs`); which connection fills `D_UV_DT`,
  `DTDT`, `DPEDT` in the full GCM?

### Offline tracer advection as a sibling component
- Pattern: with `DYCORE: FV3+ADV`, tracers outside a fixed water/cloud list are advected by
  the separate component `ADV` using the mass fluxes and Courant numbers exported by `DYN`
  (`MFX`, `MFY`, `CX`, `CY`, `PLE0`, `PLE1`) (`AdvCore_GridCompMod.F90`, `run`;
  `regression/adv-dyn/root.yaml`; the Superdyn `FV3+ADV` branch is commented out on this line).  [Verified]
- Likely reason: advect many chemistry/aerosol tracers with a different split or time step than
  the core, without the cost of carrying them through `fv_dynamics`.  [Inferred]
- Question for maintainers: what is the role of the time-centered mass flux option
  (`AdvCore_Advection: 2`, `computeMassFluxes`)?

### Water tracers packed by name
- Pattern: `FV_Run` packs `QLCN`/`QLLS`, `QICN`/`QILS`, `CLCN`/`CLLS` into the FV3 liquid, ice and
  cloud-fraction slots, carries the convective parts as extra tracers, and recovers the
  large-scale part afterwards (`FV_StateMod.F90`, `FV_Run`).  [Verified]
- Likely reason: the FV3 core knows a fixed set of water species while GEOS moist physics
  keeps convective and large-scale parts separate.  [Inferred]
- Question for maintainers: is the requirement of exactly five water tracers (`nwat_tracers == 5`)
  expected to stay?

### Component wrapper in this repository, core in a separate one
- Pattern: the MAPL layer lives here and the GFDL-derived core in the nested repository `@fvdycore`;
  the core is built into this repository's library (`CMakeLists.txt`, `SUBCOMPONENTS fvdycore`).  [Verified]
- Likely reason: track the upstream GFDL code separately from the GEOS interface.  [Inferred]
- Question for maintainers: which changes belong in `@fvdycore` and which in `FV_StateMod.F90`?

### Optional Python core
- Pattern: a build option (`BUILD_GEOS_GTFV3_INTERFACE`) and a runtime resource (`RUN_GTFV3:`)
  replace `fv_dynamics` by a call into a Python implementation through a C shim
  (`geos-gtfv3/`; `FV_StateMod.F90`, `FV_Run`).  [Verified]
- Likely reason: run the dycore on GPUs while the rest of the model stays in Fortran.  [Inferred]
- Question for maintainers: what is the role of the first-seven-tracers limit of the Python path?

### Tracer bundle as a MAPL3 service item
- Pattern: `DYN` and `ADV` declare `TRADV` as `MAPL_STATEITEM_SERVICE`, which puts the bundle in
  both the import and the export state, and `run` asserts the two share data
  (`DynCore_GridCompMod.F90`, `SetServices`, `run`).  [Verified]
- Likely reason: the source comment says the component is the provider and there is no coupling
  mechanism yet to copy the bundle.  [Verified: comment in `run`]
- Question for maintainers: is the same-data assertion meant to stay once a coupling mechanism exists?

### Tracer selection predicate shared by DYN and ADV
- Pattern: both components decide which tracers stay in the dycore with one function,
  `field_is_cloud_water_species` in `fv_shared.F90`, with opposite polarity
  (`DynCore_GridCompMod.F90`, `AdvCore_GridCompMod.F90`, `get_adjusted_tracer_bundle`).  [Verified]
- Likely reason: keep the two selections consistent.  [Inferred]
- Question for maintainers: should `FV_StateMod.F90` (the `nwat_tracers` count) use the same predicate?
