import f90nml
from datetime import datetime

from f_py_conversion import fortran_state_to_numpy
from f_py_conversion import numpy_state_to_fortran

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

    BACKEND = 'dace:gpu'

    rank = comm.Get_rank()

    if rank == 0:
        print('P:', datetime.now().isoformat(timespec='milliseconds'),
              '--in top level geos-gtfv3, backend:', BACKEND, flush=True)

    # Convert Fortran arrays to NumPy
    state_in = fortran_state_to_numpy(
        npx, npy, npz,
        is_, ie, js, je, isd, ied, jsd, jed,
        bdt, nq_tot, ptop, ks,
        ak, bk,
        # input/output arrays
        u, v, w, delz,
        pt, delp, q,
        ps, pe, pk, peln, pkz,
        phis, q_con, omga,
        ua, va, uc, vc,
        mfx, mfy, cx, cy, diss_est)
    # write_sum_of_vars(comm, state_in)
    # write_shape_or_value(comm, state_in)

    if rank == 0:
        print('P:', datetime.now().isoformat(timespec='milliseconds'),
              '--created input data', flush=True)

    # Setup dycore object
    # This block gets executed only at the first time step
    global GTFV3_DYCORE
    if GTFV3_DYCORE is None:
        namelist = f90nml.read('input.nml')
        GTFV3_DYCORE = GeosDycoreWrapper(namelist, bdt, comm, BACKEND)
    assert GTFV3_DYCORE is not None, 'GTFV3_DYCORE is None'

    if rank == 0:
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

    # Convert NumPy arrays back to Fortran
    numpy_state_to_fortran(
        state_out,
        u, v, w, delz,
        pt, delp, q,
        ps, pe, pk, peln, pkz,
        phis, q_con, omga,
        ua, va, uc, vc,
        mfx, mfy, cx, cy, diss_est)

    if rank == 0:
       print('P:', datetime.now().isoformat(timespec='milliseconds'),
             '--converted NumPy arrays to Fortran', flush=True)
