import numpy as np

import fv3core
import fv3core.testing
import fv3gfs.util as util

from datetime import datetime

def geos_gtfv3_run(comm, spec, driver_object, dycore, fv3_input_data):
    '''
    Inputs:
        comm - MPI communicator
        spec - fv3core._config that has already been created
        driver_object
        fv3_input_data (u, v, w, etc.)
    Output:
        fv3_output_data
    '''

    rank = comm.Get_rank()

    if (rank == 0):
        print('P:', datetime.now().isoformat(timespec='milliseconds'),
              '--in top level function', flush=True)

    # Create state
    state = driver_object.state_from_inputs(fv3_input_data)
    if rank == 0:
        print('P:', datetime.now().isoformat(timespec='milliseconds'),
              '--created state', flush=True)

    # # Instantiate DryConvectiveAdjustment
    # if spec.namelist.fv_sg_adj > 0:
    #     fv_subgrid_z = fv3core.DryConvectiveAdjustment(
    #         spec.grid.grid_indexing,
    #         spec.namelist.nwat,
    #         spec.namelist.fv_sg_adj,
    #         spec.namelist.n_sponge,
    #         spec.namelist.hydrostatic,)
    #     if (spec.grid.rank == 0):
    #         print('P:', datetime.now().isoformat(timespec='milliseconds'),
    #               '--instantiated DryConvectiveAdjustment', flush=True)

    # Run a single timestep of DynamicalCore
    dycore.step_dynamics(
        state,
        fv3_input_data['consv_te'],
        fv3_input_data['do_adiabatic_init'],
        fv3_input_data['bdt'],
        fv3_input_data['ptop'],
        fv3_input_data['n_split'],
        fv3_input_data['ks'])
    if (spec.grid.rank == 0):
        print('P:', datetime.now().isoformat(timespec='milliseconds'),
              '--ran DynamicalCore::step_dynamics', flush=True)

    # Retrieve output data from state
    gtfv3_output_data = driver_object.outputs_from_state(state)
    if rank == 0:
        print('P:', datetime.now().isoformat(timespec='milliseconds'),
              '--retrieved output data from state', flush=True)

    return gtfv3_output_data
