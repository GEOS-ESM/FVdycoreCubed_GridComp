# FVdycoreCubed Regression Tests

This directory contains CTest-based regression cases for
`@FVdycoreCubed_GridComp`.

## Test Harness Overview

- Test case names are listed in `cases.txt`.
- `CMakeLists.txt` registers each case as a CTest and assigns label `REGRESSION`.
- `run_case.cmake` creates a temporary run directory, copies the selected case
  files, and runs `GEOS.x` with MPI.
- The harness currently runs each case with 6 MPI ranks (`num_procs = 6`).
- If the run fails, `run_case.cmake` emits `PET0.ESMF_LogFile` (when present)
  and exits with a CTest failure.

## Available Cases

### `adv-dyn`

Coupled advection + dynamics test.

 Key files:
- `mapl.yaml`, `cap_driver.yaml`, `cap_gridcomp.yaml`, `cap_restart.yaml`
- `root.yaml`
- `dyn.yaml`, `adv.yaml`
- `data-moist.yaml` — `DataMoist` component; owns `Q`, that is advected by DYN
- `data-tracers.yaml` — `DataTracers` component; owns synthetic tracers `Q001`, `Q002`, that are advected by ADV
- `data-ana.yaml`
- `input.nml`

### `dyn-sa`

Held-Suarez style standalone dynamics test.

 Key files:
- `mapl.yaml`, `cap_driver.yaml`, `cap_gridcomp.yaml`, `cap_restart.yaml`
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
ctest -R '^dyn-sa$' -j8 --output-on-failure
```

## Troubleshooting

- If CTest reports `Error running <case>`, inspect the emitted
  `PET0.ESMF_LogFile` in command output.
- If GEOS shared libraries are not found at runtime, confirm tests are run via
  CTest so the `LD_LIBRARY_PATH`/`DYLD_LIBRARY_PATH` test environment from
  `CMakeLists.txt` is applied.
- If no regression tests are found, verify this component was built with tests
  enabled and rerun `make -j8 build-tests`.

## Optional Baseline Comparison Scripts

Each case includes a `cmpchk.py` script that compares produced checkpoints
against external baseline datasets.

Notes:
- These scripts are not part of the CTest pass/fail criteria.
- Baseline paths are environment-specific and currently hard-coded in each
  `cmpchk.py`.
- Run manually from a completed case run directory after ensuring Python has
  `numpy` and `netCDF4` available.

Example:

```bash
python cmpchk.py
```
