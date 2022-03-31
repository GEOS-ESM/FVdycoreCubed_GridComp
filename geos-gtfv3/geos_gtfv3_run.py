import time
import numpy as np

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

    # Create state
    state = driver_object.state_from_inputs(state_in)

    # # Instantiate DryConvectiveAdjustment
    # if spec.namelist.fv_sg_adj > 0:
    #     fv_subgrid_z = fv3core.DryConvectiveAdjustment(
    #         spec.grid.grid_indexing,
    #         spec.namelist.nwat,
    #         spec.namelist.fv_sg_adj,
    #         spec.namelist.n_sponge,
    #         spec.namelist.hydrostatic,)

    # Run a single timestep of DynamicalCore
    start = time.time()
    dycore.step_dynamics(
        state,
        conserve_total_energy = dycore_config.consv_te,
        do_adiabatic_init = True,
        timestep = state_in['bdt'],
        n_split = dycore_config.n_split)
    end = time.time()
    print('[{:02d}] Time taken by dycore.step_dynamics: {:.5f}s'.format(rank, end-start), flush=True)

    # Retrieve output data from state
    state_out = driver_object.outputs_from_state(state)

    return state_out
