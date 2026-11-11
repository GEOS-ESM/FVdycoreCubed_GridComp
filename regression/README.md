# FVdycoreCubed Regression Tests

This directory contains CTest-based regression cases for `@FVdycoreCubed_GridComp`.

## Test Harness Overview

- Test case names are listed in `cases.txt`.
- `CMakeLists.txt` registers each case as a CTest with label `REGRESSION` (`esma_add_regression_tests` in `@cmake`, `NAME fv3`, `DATA_PATH FVdycoreCubed_GridComp`).
- Regression data is read from `$LOCAL_REGRESSION_DATA_DIR/FVdycoreCubed_GridComp/<case>`, or from `<build>/regression-data/FVdycoreCubed_GridComp/<case>` when `LOCAL_REGRESSION_DATA_DIR` is not set at configure time. The test `fv3/sync_data`, which every case requires, syncs that directory from S3 when ctest runs.
- `run_case.cmake` creates a temporary run directory, copies the case files and any input restarts (`checkpoints/<time>/` in the data directory) into it, and runs `GEOS.x` on 6 MPI ranks (`num_procs`). If the run fails, it prints `PET0.ESMF_LogFile` (when present) and the test fails.
- After the run, each `*.nc` file in the data directory's `checkpoints/last/` (the baseline) is compared exactly with the run's checkpoint of the same name, with `nccmp` (NaNs equal; `lons`, `lats`, `corner_lons` and `corner_lats` excluded). The comparison is done only for IntelLLVM (ifx) builds. If the baseline directory is missing, the case still runs and only the comparison is skipped, with a warning.
- Each case runs one hour of model time on a C12 cubed sphere with 91 levels. `adv-dyn` and `dyn-sa` cold-start `DYN` (`COLDSTART: true`, `CASE_ID: 1`, the Jablonowski-Williamson steady state) and run it adiabatically (`ADIABATIC: true`).

## Logging

Each case directory has a `logging.yaml`, loaded by `pflogger: logging.yaml` in `mapl.yaml`. It logs to the console only and sets the `MAPL` logger to ERROR. Without it, MAPL's default logging writes its warnings (once per variable spec per rank, since the cases run without a field dictionary) to `warnings_and_errors.log` under an MPI lock, which serializes the ranks and slows each test. Keep it in new cases.

## Available Cases

### `adv-dyn`

Coupled dynamics and offline advection: `DYN` (`DYCORE: FV3+ADV`) and `ADV`, with `DYN`'s mass fluxes, Courant numbers and edge pressures (`MFX`, `MFY`, `CX`, `CY`, `PLE0`, `PLE1`) connected to `ADV`. `DYN` advects `Q` (a water species) and `ADV` the other tracers. Two steps of 1800 s.

Key files:
- `mapl.yaml`, `cap_driver.yaml`, `cap_gridcomp.yaml`, `cap_restart.yaml`, `logging.yaml`, `history.yaml`
- `root.yaml`
- `dyn.yaml`, `adv.yaml`
- `data-moist.yaml`: `DataMoist` component; owns `Q`, advected by `DYN`
- `data-tracers.yaml`: `DataTracers` component; owns the synthetic tracers `Q001`, `Q002`, advected by `ADV`
- `data-phis-varflt.yaml`: `DataPhisVarflt` component; provides `PHIS` and `VARFLT` to `DYN`
- `data-ana.yaml`: `DataAna` component; provides the analysis increments `DQVANA`, `DQLANA`, `DQIANA`, `DQRANA`, `DQSANA`, `DQGANA` to `DYN`
- `input.nml`

### `adv-sa`

Standalone advection (`ADV` only). `DYCORE` is unset, so `ADV` runs `fv_setup` itself. The mass fluxes, Courant numbers and edge pressures (`MFX`, `MFY`, `CX`, `CY`, `PLE0`, `PLE1`) are read from the `DataForAdv` internal restart, which is a symlink to the `adv-dyn` baseline's `checkpoints/last/DYN_export.nc` (MAPL reads only the fields it declares). Two steps of 1800 s.

Key files:
- `mapl.yaml`, `cap_driver.yaml`, `cap_gridcomp.yaml`, `cap_restart.yaml`, `logging.yaml`, `history.yaml`
- `root.yaml`
- `adv.yaml`
- `data-for-adv.yaml`: `DataForAdv` component; provides the advection inputs from its internal restart
- `data-tracers.yaml`: `DataTracers` component; owns the synthetic tracers `Q001`, `Q002`, advected by `ADV`
- `input.nml`

Regression data: `checkpoints/1891-03-02T00:00:00/DataForAdv_internal.nc` (input restart, symlink; copied into the run directory) and `checkpoints/last/` (baseline).

### `dyn-sa`

Standalone dynamics: `DYNsa` (`DYCORE: FV3`, `FV3_STANDALONE: true`, so no analysis increments are applied) with `DataMoist` owning the advected tracers (`Q` and the other water species). Three steps of 1200 s.

Key files:
- `mapl.yaml`, `cap_driver.yaml`, `cap_gridcomp.yaml`, `cap_restart.yaml`, `logging.yaml`, `history.yaml`
- `root.yaml`
- `dyn-sa.yaml`
- `data-moist.yaml`
- `input.nml`

## Running the Tests

From your build directory:

```bash
# Ensure test targets are built
make -j8 build-tests

# Show discovered regression tests
ctest -N -L REGRESSION

# Run all FVdycoreCubed regression tests
ctest -L REGRESSION -j8 --output-on-failure

# Run one test case
ctest -R '^adv-dyn$' -j8 --output-on-failure
ctest -R '^adv-sa$' -j8 --output-on-failure
ctest -R '^dyn-sa$' -j8 --output-on-failure
```

## Troubleshooting

- If CTest reports `Error running <case>`, inspect the `PET0.ESMF_LogFile` printed in the command output.
- If GEOS shared libraries are not found at run time, run the tests through CTest so that the `LD_LIBRARY_PATH`/`DYLD_LIBRARY_PATH` test environment set by `esma_add_regression_tests` is applied.
- If no regression tests are found, check that this component was built with tests enabled and rerun `make -j8 build-tests`.
- If a case passes with a "Baseline directory ... not found" warning, nothing was compared: its data directory is missing. Check `LOCAL_REGRESSION_DATA_DIR` and the output of `fv3/sync_data`.
