import numpy as np

from datetime import datetime

def geos_gtfv3_run(comm, dycore_config, driver_object, dycore, state_in):
    '''
    Inputs:
        comm - MPI communicator
        dycore_config - namelist for dynamical core
        driver_object
        dycore - initialized DynamicalCore object
        state_in (u, v, w, etc.)
    Output:
        state_out
    '''

    rank = comm.Get_rank()

    # if rank == 0:
    #     print('P:', datetime.now().isoformat(timespec='milliseconds'),
    #           '--running step dynamics', flush=True)

    # Create state
    state = driver_object.state_from_inputs(state_in)
    # if rank == 0:
    #     print('P:', datetime.now().isoformat(timespec='milliseconds'),
    #           '--created state', flush=True)

    # # Instantiate DryConvectiveAdjustment
    # if spec.namelist.fv_sg_adj > 0:
    #     fv_subgrid_z = fv3core.DryConvectiveAdjustment(
    #         spec.grid.grid_indexing,
    #         spec.namelist.nwat,
    #         spec.namelist.fv_sg_adj,
    #         spec.namelist.n_sponge,
    #         spec.namelist.hydrostatic,)
    #     if rank == 0:
    #         print('P:', datetime.now().isoformat(timespec='milliseconds'),
    #               '--instantiated DryConvectiveAdjustment', flush=True)

    # Run a single timestep of DynamicalCore
    dycore.step_dynamics(
        state,
        conserve_total_energy = dycore_config.consv_te,
        do_adiabatic_init = True,
        timestep = state_in['bdt'],
        n_split = dycore_config.n_split)
    # if rank == 0:
    #     print('P:', datetime.now().isoformat(timespec='milliseconds'),
    #           '--ran DynamicalCore::step_dynamics', flush=True)

    # Retrieve output data from state
    state_out = driver_object.outputs_from_state(state)
    # if rank == 0:
    #     print('P:', datetime.now().isoformat(timespec='milliseconds'),
    #           '--retrieved output data from state', flush=True)

    return state_out
