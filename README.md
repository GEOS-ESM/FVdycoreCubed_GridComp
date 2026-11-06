# FV dynamical core and gridded component

## Contributing

Please check out our [contributing guidelines](CONTRIBUTING.md).

## Adaptive-split restarts

With `ADJUST_SPLITS: .true.`, the dynamics internal restart includes `K_SPLIT`,
`N_SPLIT`, `Q_SPLIT`, and `BASE_ACOUSTIC_DT`. These are one-element, double-precision
ungridded fields; the split counts are integer-valued. The baseline acoustic
timestep is computed after `FV_Setup` resolves the initial splits and remains
fixed while the active splits adjust.

The fields are optional and default to `-1` when absent. When bootstrapping an
older restart, enable MAPL's optional-field bootstrapping with
`MAPL_ENABLE_BOOTSTRAP: YES`. Missing values are initialized from `FV_Setup`,
with a log message; saved values, including zero split counts, are retained.
Split adjustment occurs once per dynamics timestep, not during restart
initialization.

With `ADJUST_SPLITS: .false.`, the additional fields are not registered or read,
and configured splitting and the existing restart layout are unchanged.

The focused split-restart tests use pFUnit and the CMake-configured Fortran
compiler, like MAPL's tests. In a GEOSgcm build with pFUnit available, run:

```sh
cmake --build <build-dir> --target build-tests
ctest --test-dir <build-dir> -R FVdycoreCubed.adaptive_split_restart --output-on-failure
```

The tests exercise the production state declaration and split routines with
stubbed MAPL storage and single-process collectives in both default-real
precisions. Full MAPL restart I/O and prognostic start/stop regression still
require the GEOSgcm fixture.

CI runs the `ESSENTIAL` tests in the existing `build_gcm` build tree via the
shared workflow's `run-tests: true` input (default: `false` for other callers).
The `tests` target builds test executables incrementally and runs only tests
labeled `ESSENTIAL`; it does not configure a second model build.

## License

All files are currently licensed under the Apache-2.0 license, see [`LICENSE`](LICENSE).

Previously, the code was licensed under the [NASA Open Source Agreement, Version 1.3](LICENSE-NOSA).
