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

The focused split-restart tests require Python 3 and `gfortran` and can be run
with `python3 -m unittest discover -s tests -p 'test_restart_splits.py' -v`.
They exercise the production split routines with stubbed MAPL storage and
single-process collectives in both default-real precisions; full prognostic
start/stop regression still requires the GEOSgcm fixture.

## License

All files are currently licensed under the Apache-2.0 license, see [`LICENSE`](LICENSE).

Previously, the code was licensed under the [NASA Open Source Agreement, Version 1.3](LICENSE-NOSA).
