from datetime import datetime
from f_py_conversion import fortran_state_to_numpy
from f_py_conversion import numpy_output_data_to_fortran
from geos_gtfv3_initialize import setup_namelist
from geos_gtfv3_initialize import setup_dycore
from geos_gtfv3_run import geos_gtfv3_run
from geos_gtfv3_debug import write_sum_of_vars

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

    assert nq_tot == 7, f'Expected 7, received: {nq_tot}'

    BACKEND = 'gtcuda'

    rank = comm.Get_rank()

    if rank == 0:
        print('P:', datetime.now().isoformat(timespec='milliseconds'),
              '--in top level geos-gtfv3, backend:', BACKEND, flush=True)

    # Initialize namelist
    # This function gets exercised only the first timestep
    dycore_config = setup_namelist(npx, npy, npz, ntiles, [layout_1, layout_2])
    assert dycore_config is not None

    if rank == 0:
        print('P:', datetime.now().isoformat(timespec='milliseconds'),
              '--setup namelist completed', flush=True)
        print(dycore_config)

    # Convert Fortran arrays to NumPy
    # This function gets exercised at every timestep
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

    if rank == 0:
        print('P:', datetime.now().isoformat(timespec='milliseconds'),
              '--created input data', flush=True)

    # Initialize dycore
    # This function gets exercised only the first timestep
    driver_object, dycore = setup_dycore(
        BACKEND, comm,
        npx, npy, npz,
        is_, ie, js, je, isd, ied, jsd, jed,
        state_in)
    assert dycore is not None, 'dycore is None'

    if rank == 0:
        print('P:', datetime.now().isoformat(timespec='milliseconds'),
              '--initialized dycore', flush=True)

    # Run gtFV3
    state_out = geos_gtfv3_run(
        comm,
        dycore_config,
        driver_object,
        dycore,
        state_in)
    # write_sum_of_vars(comm, state_out)

    # Convert NumPy arrays back to Fortran
    numpy_output_data_to_fortran(
        state_out, nq_tot,
        u, v, w, delz,
        pt, delp, q,
        ps, pe, pk, peln, pkz,
        phis, q_con, omga,
        ua, va, uc, vc,
        mfx, mfy, cx, cy, diss_est)

    if rank == 0:
       print('P:', datetime.now().isoformat(timespec='milliseconds'),
             '--converted NumPy arrays to Fortran', flush=True)
