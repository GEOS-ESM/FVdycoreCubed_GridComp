# interp_restarts.F90 cleanup notes

Working notes for cleaning up `interp_restarts.F90` (the program behind
`interp_restarts.x`, built by `CMakeLists.txt`). It is run by
`GEOS_Util/pre/remap_restart/remap_upper.py`.

- Base: `@FVdycoreCubed_GridComp` @ `69a5503` (tag `v3.1.1`, detached HEAD).
- Line numbers below refer to the **original** file at that commit unless
  stated otherwise.
- Related code read: `fv_regrid_c2c.F90` (`get_geos_ic`,
  `get_geos_cubed_ic`), `fv_regridding_utils.F90` (`fv_var`, `fv_rst`),
  `@fvdycore/model/fv_control.F90` (`fv_init1`, `fv_init2`, `run_setup`).

## Status

| Step | Description | Status |
|---|---|---|
| 1 | Remove unused code and dead paths (no behaviour change) | Done, not yet built or tested |
| 2 | Decide on `-scalers` | Open |
| 3 | Small, contained bug fixes | Open |
| 4 | Readability refactors | Open |
| 5 | Moist tracer index clash | Open, needs investigation |
| 6 | Same cleanup for `interp_restarts_bin.F90` | Open |
| 7 | Record remaining bugs under Observations in `docs/ARCHITECTURE.md` | Open |

## Step 1: done (uncommitted)

With whitespace ignored, the diff is 7 lines added and 65 removed; the rest is
re-indentation.

- Removed the `external :: getarg, iargc` / `integer iargc` block (34–37); the
  program uses `command_argument_count` / `get_command_argument`.
- Removed `NOTE` from the `use mpp_mod` list (11).
- Removed `tracer_names` (StringIntegerMap) and `use gFTL_StringIntegerMap`
  (23, 63); the map was filled in but never read.
- Removed `moist_variables` (61, 278); filled in but never read.
- First moist loop (270–289): removed the `iq0` / `lcnt_var` work that only fed
  `tracer_names`. The loop now only decrements `nmoist` for 2D variables.
  `iq0` and `lcnt_var` are still declared because the moist write loop uses them.
- Removed `subroutine rs_count` (741–756); never called.
- Removed `has_attr` from `reset_stretch_params` (771–776); it was assigned but
  never tested. The three `remove_attribute` calls are unchanged.
- Removed the empty `if (.not. Atm(1)%flagstruct%hydrostatic)` block
  (490–497) that held only commented-out `MAPL_NCIOAddVar` calls.
- Removed the always-true `file_exist(...)` checks around the fvcore (471) and
  moist (580) writes and dedented their contents. The program already stops at
  232 and 293 if either file is missing; those checks are unchanged.
- Renamed the leftover memory label `'interp_restarts: rs_count'` to
  `'interp_restarts: moist count'` (296). This only changes a log line.
- `-input_files` / `-output_files` parsing no longer reuses `nq` as a counter;
  the loops are now `do j=1,n` with `get_command_argument(i+j, ...)`, which
  gives the same result. `nq` is now only the tracer count.

Kept on purpose:

- `ks`, `ptop`, `pint`: outputs required by `set_eta`.
- The DZ/W write branch (`if (.not. Atm(1)%flagstruct%hydrostatic)`). Line 197
  sets `hydrostatic = .false.`, but `fv_init2` → `run_setup` reads
  `fv_core_nml` from `input.nml` and can override it, so the branch is not
  dead.

To finish step 1:

1. Build `interp_restarts.x`.
2. Run one `remap_upper` case with the old and new executables and compare the
   output restart files; they should match.
3. Create a branch (`git switch -c <branch>`), commit, and open a pull request
   in the FVdycoreCubed_GridComp repository.

## Step 2: decide what `-scalers` should do

`-scalers T|F` sets `scale_rst` (97, 151–159), which is never read.
`get_geos_ic` always sets `Atm(1)%flagstruct%adjust_dry_mass = .true.`
(`fv_regrid_c2c.F90:132`). `remap_upper.py` passes `-scalers F` when the input
`drymass` option is off, and the program ignores it without warning.

Options:

- **Make it work:** pass `scale_rst` into `get_geos_ic` and use it for
  `adjust_dry_mass`. Changes `fv_regrid_c2c.F90` and changes output for runs
  with `drymass` off.
- **Remove it:** delete the option from the program and from `remap_upper.py`,
  and document that dry mass is always adjusted.

This is a scientific decision; check with the `remap_restart` owners first.

## Step 3: small, contained bug fixes (proposed)

- `VERIFY_` (line 1) calls `MPI_Abort(A)` with one argument; the MPI interface
  is `MPI_Abort(comm, errorcode, ierror)`. It compiles only because `mpif.h`
  has no explicit interfaces. Use `MPI_Abort(MPI_COMM_WORLD, A, ierr)`.
- Memory leak: in the extra-file write loop, `allocate(r4_local(...,nlev))`
  (702) runs once per variable, but `deallocate(r4_local)` (723) runs once per
  file. `r4_local` is a pointer, so this leaks instead of failing.
- With no `-input_files` (which `remap_upper.py` does when there are no extra
  restarts), `size(extra_files)` (327) is taken on an unallocated array, and
  `size(rst_files)` (658, and in `get_geos_ic`) on a null pointer. Guard these
  checks.
- `vname` is `character(62)` (60); `fv_var%name` is `character(len=128)`.
- Bare `stop` when dimension 3 is neither `km` nor `km+1` (375, 386); replace
  with `mpp_error(FATAL, ...)` and a message.

## Step 4: readability refactors (proposed)

- Duplicated stretched-grid flag setup (199–208, 243–252) and
  `CubedSphereGridFactory` construction (409–416, 433–440): two small helpers.
- Duplicated "level count from dimension 3" logic for rank 3 and rank 4
  (368–376, 379–387): one helper.
- Moist write index chain (617–644): the `iwat < nwat` test is always true,
  because `get_geos_cubed_ic` sets `nwat` to the count of these exact names
  (`fv_regrid_c2c.F90:484–505`). The chain duplicates the reader chain at
  `fv_regrid_c2c.F90:508–529`; both should call one shared function (see
  step 5).
- fvcore write section: values are copied into a 4-byte buffer and converted
  back to `FVPRC` only so `prt_mxm` can print them. The two `PE` loops
  (531–537) can be merged, and the two `if (AmWriter)` AK/BK lines combined.
- Extra-file `MAPL_IOChangeRes` chain (671–687): use `.not. x` instead of
  `x .eqv. .false.`, or build the dimension list in one routine.
- `r4_local`, `r8_local`, `pt_local`, `r4_local2D` are pointers but are only
  ever allocated; they can be `allocatable`.

Test step 4 the same way as step 1 (compare output restarts).

## Step 5: moist tracer index clash (needs investigation)

Reader (`fv_regrid_c2c.F90:508–529`) and writer (617–644) give water species
fixed slots 1–8 (Q, QLLS, QLCN, QILS, QICN, QRAIN, QSNOW, QGRAUPEL), while other
tracers start at `nwat+1`. If a file has only some species (for example Q and
QRAIN but no QLLS), slots can collide or go past `ncnst`. Reproduce with a real
moist restart, then fix with one shared `moist_tracer_index(name)` function in
`fv_regridding_utils.F90`, used by both the reader and the writer. This touches
the regridding code, so make it a separate change.

## Other issues found (not yet scheduled)

- `MAPL_IOChangeRes` chain for extra files (671–687) has no fallback; other
  combinations (for example edge levels plus `unknown_dim1`) leave `OutCfg`
  undefined.
- Argument parsing stops the file list at any argument containing `-`
  (`index(astr,'-')`, 120, 134), so a file name with a hyphen ends the list
  early.
- `-im` / `-lm` have no defaults or checks; if missing, `npx` / `npz` are used
  uninitialized.
- `MPI_Info_create` (459) is never paired with `MPI_Info_free`.
- Most `rc=status` results are not checked; the extra-file `MAPL_VarWrite`
  calls pass no `rc` at all.
- For moist, every rank runs `MAPL_IOChangeRes`; for fvcore only writers do.

## Step 6: `interp_restarts_bin.F90`

Shares much of this code (`rs_count`, `rst_files` handling, the write loops)
and is chosen by `remap_upper.py` for binary restarts. Review it for the same
cleanup once steps 1–4 are settled here.
