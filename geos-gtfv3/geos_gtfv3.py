import f90nml
from datetime import datetime

from f_py_conversion import FortranPythonConversion
from pace.fv3core.initialization.geos_wrapper import GeosDycoreWrapper
from geos_gtfv3_debug import write_shape_or_value

#-GLOBAL-variables-start-
GTFV3_DYCORE = None
#-GLOBAL-variables-end-

def geos_gtfv3(
        comm,
        npx, npy, npz, ntiles,
        is_, ie, js, je, isd, ied, jsd, jed,
        bdt, nq_tot, ng, ptop, ks, layout_1, layout_2, adiabatic,
        u, v, w, delz,
        pt, delp, q,
        ps, pe, pk, peln, pkz,
        phis, q_con, omga,
        ua, va, uc, vc,
        ak, bk,
        mfx, mfy, cx, cy, diss_est):

    BACKEND = 'gt:gpu'

    rank = comm.Get_rank()
    RANK_PRINT = 0

    if rank == RANK_PRINT:
        print('P:', datetime.now().isoformat(timespec='milliseconds'),
              '--in top level geos-gtfv3, backend:', BACKEND, flush=True)

    # For Fortran<->NumPy conversion
    f_py = FortranPythonConversion(
        npx, npy, npz,
        is_, ie, js, je,
        isd, ied, jsd, jed,
        nq_tot)

    # Convert Fortran arrays to NumPy
    state_in = f_py.fortran_to_numpy(
        # input
        u, v, w, delz,
        pt, delp, q,
        ps, pe, pk, peln, pkz,
        phis, q_con, omga,
        ua, va, uc, vc,
        mfx, mfy, cx, cy, diss_est)
    # write_sum_of_vars(comm, state_in)
    # write_shape_or_value(comm, state_in)

    if rank == RANK_PRINT:
        print('P:', datetime.now().isoformat(timespec='milliseconds'),
              '--fortran->numpy, transpose, sp->dp, swap axes', flush=True)

    # Setup dycore object
    # This block gets executed only at the first time step
    global GTFV3_DYCORE
    if GTFV3_DYCORE is None:
        namelist = f90nml.read('input.nml')
        GTFV3_DYCORE = GeosDycoreWrapper(namelist, bdt, comm, BACKEND)
    assert GTFV3_DYCORE is not None, 'GTFV3_DYCORE is None'

    if rank == RANK_PRINT:
        print('P:', datetime.now().isoformat(timespec='milliseconds'),
              '--initialized dycore', flush=True)

    # Run gtFV3
    state_out = GTFV3_DYCORE(
        state_in['u'], state_in['v'], state_in['w'], state_in['delz'],
        state_in['pt'], state_in['delp'], state_in['q'],
        state_in['ps'], state_in['pe'], state_in['pk'], state_in['peln'], state_in['pkz'],
        state_in['phis'], state_in['q_con'], state_in['omga'],
        state_in['ua'], state_in['va'], state_in['uc'], state_in['vc'],
        state_in['mfxd'], state_in['mfyd'], state_in['cxd'], state_in['cyd'], state_in['diss_estd'])

    if rank == RANK_PRINT:
        print('P:', datetime.now().isoformat(timespec='milliseconds'),
              '--ran dycore', flush=True)

    # Convert NumPy arrays back to Fortran
    f_py.numpy_to_fortran(
        # input
        state_out,
        # output
        u, v, w, delz,
        pt, delp, q,
        ps, pe, pk, peln, pkz,
        phis, q_con, omga,
        ua, va, uc, vc,
        mfx, mfy, cx, cy, diss_est)

    if rank == RANK_PRINT:
       print('P:', datetime.now().isoformat(timespec='milliseconds'),
             '--numpy->fortran, transpose, dp->sp, swap axes', flush=True)
