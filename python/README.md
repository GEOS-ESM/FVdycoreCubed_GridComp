# pyFV3 integration

## Source structure

- `@pyFV3`: contains the model code ported to the NDSL middleware. This is shared code with NOAA and lives in GEOS as a `mepo` component.
- `interface`: all things to move data in and out between Fortran and Python
  - `cffi_lib`: all files required for the three way .so build
  - `pyFV3_interface.py`: entry point for Python interface
  - `pyFV3_wrapper.py`: bespoke wrapper around pyFV3 for GEOS (e.g. driver)

The only changes that are not included here are the required build changes in the parent CMakeLists.txt

## Fortran <-> Python interface

Using `cffi` we have a 3-way Fortran <-> C <-> Python interface with minimum overload.

Memory layouts are not compatible so a copy (and eventual upload to device) occurs on the Python side. Future accelerations will look into keeping a single memory layout for CPU based backend.

The configuration of the dycore is passed via the `FV_Flags` structure/user type. This can be extended but _cannot included anything but POD_.

## Control pyFV3 behavior

Input control to pyFV3 via environement variables

- `GEOS_PYFV3_BACKEND`: `NDSL` backend for `pyFV3`, defaults to `gt:gpu`
- `GEOS_PYFV3_SINGLE_RANK_OVERRIDE`: boolean to control a single rank override for debug

## Contact

Any questions contact the SMT team: <florian.g.deconinck@nasa.gov>
