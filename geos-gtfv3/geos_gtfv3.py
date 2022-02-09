import numpy as np

import gt4py
import fv3core
import fv3core._config as spec
import fv3core.testing
import fv3gfs.util as util

from datetime import datetime

# import h5py

def geos_gtfv3(
        comm,
        npx, npy, npz,
        is_, ie, js, je,
        isd, ied, jsd, jed,
        bdt, nq_tot, ng, ptop, ks, layout_1, layout_2,
        adiabatic_int,
        hydrostatic_int, z_tracer_int, make_nh_int, fv_debug_int,
        reproduce_sum_int, do_sat_adj_int, do_vort_damp_int, rf_fast_int, fill_int,
        ncnst, n_split, k_split, fv_sg_adj, n_sponge, n_zfilter, nwat,
        hord_tr, hord_tm, hord_dp, hord_mt, hord_vt,
        nord, kord_tm, kord_tr, kord_wz, kord_mt,
        d_ext, beta, vtdm4, ke_bg, d_con, d2_bg, d2_bg_k1, d2_bg_k2,
        p_fac, a_imp, dddmp, d4_bg, rf_cutoff, tau, consv_te,
        fv3_input_data, grid_data,
        nested_int, stretched_grid_int, da_min, da_min_c):

    rank = comm.Get_rank()
    backend = 'gtx86'
    if (rank == 0):
        print('P:', datetime.now().isoformat(timespec='milliseconds'),
              '--in top level function with backend', backend, flush=True)

    nranks = comm.Get_size()
    for i in range(nranks):
        if i == rank:
            print('P: rank:', rank, flush=True)
            print('P: q:',
                  np.sum(fv3_input_data['qvapor']),
                  np.sum(fv3_input_data['qliquid']),
                  np.sum(fv3_input_data['qice']),
                  np.sum(fv3_input_data['qrain']),
                  np.sum(fv3_input_data['qsnow']),
                  np.sum(fv3_input_data['qgraupel']),
                  np.sum(fv3_input_data['qcld']),
                  flush=True)
        comm.Barrier()

    fv3core.set_backend(backend)
    fv3core.set_rebuild(False)
    fv3core.set_validate_args(True)
    # fv3core.utils.global_config.set_do_halo_exchange(True)

    spec.set_namelist_from_dict({
        'npx': npx, 'npy': npy, 'npz': npz,
        'layout': [layout_1, layout_2],
        'adiabatic': False if adiabatic_int == 0 else True,
        'hydrostatic': False if hydrostatic_int == 0 else True,
        'z_tracer': False if z_tracer_int == 0 else True,
        'make_nh': False if make_nh_int == 0 else True,
        'fv_debug': False if fv_debug_int == 0 else True,
        # reproduce_sum
        'do_sat_adj': False if do_sat_adj_int == 0 else True,
        'do_vort_damp': False if do_vort_damp_int == 0 else True,
        'rf_fast': False if rf_fast_int == 0 else True,
        'fill': False if fill_int == 0 else True,
        # 'ncnst': ncnst,
        'n_split': n_split,
        'k_split': k_split,
        'fv_sg_adj': fv_sg_adj,
        'n_sponge': n_sponge,
        'n_zfilter': n_zfilter,
        'nwat': nwat,
        'hord_tr': hord_tr,
        'hord_tm': hord_tm,
        'hord_dp': hord_dp,
        'hord_mt': hord_mt,
        'hord_vt': hord_vt,
        'nord': nord,
        'kord_tm': kord_tm,
        'kord_tr': kord_tr,
        'kord_wz': kord_wz,
        'kord_mt': kord_mt,
        'd_ext': d_ext,
        'beta': beta,
        'vtdm4': vtdm4,
        'ke_bg': ke_bg,
        'd_con': d_con,
        'd2_bg': d2_bg,
        'd2_bg_k1': d2_bg_k1,
        'd2_bg_k2': d2_bg_k2,
        'p_fac': p_fac,
        'a_imp': a_imp,
        'dddmp': dddmp,
        'd4_bg': d4_bg,
        'rf_cutoff': rf_cutoff,
        'tau': tau,
        'delt_max': 0.002}) # namelist

    layout = spec.namelist.layout
    partitioner = util.CubedSpherePartitioner(util.TilePartitioner(layout))
    communicator = util.CubedSphereCommunicator(comm, partitioner)

    # Create dicts of input and grid data
    fv3_input_data.update({
        'ak': grid_data['ak'], 'bk': grid_data['bk'],
        'do_adiabatic_init': True, # Has to be True
        'consv_te': consv_te, 'bdt': bdt, 'ptop': ptop,
        'n_split': n_split, # Need n_split here as well
        'ks': ks, 'comm': communicator}) # fv3_input_data
    grid_data.update({
        'npx': npx, 'npy': npy, 'npz': npz,
        'is_': is_, 'ie': ie, 'js': js, 'je': je,
        'isd': isd, 'ied': ied, 'jsd': jsd, 'jed': jed,
        'nested': False if nested_int == 0 else True,
        'stretched_grid': False if stretched_grid_int == 0 else True,
        'da_min': da_min, 'da_min_c': 0.0,}) # grid_data
    if rank == 0:
        print('P:', datetime.now().isoformat(timespec='milliseconds'),
              '--created input and grid data', flush=True)

    # Create grid
    grid = fv3core.testing.TranslateGrid(grid_data, comm.Get_rank()).python_grid()
    if rank == 0:
        print('P:', datetime.now().isoformat(timespec='milliseconds'),
              '--translated data', flush=True)
    spec.set_grid(grid)
    if rank == 0:
        print('P:', datetime.now().isoformat(timespec='milliseconds'),
              '--created grid', flush=True)

    # Create state
    driver_object = fv3core.testing.TranslateFVDynamics([grid])
    state = driver_object.state_from_inputs(fv3_input_data)
    if rank == 0:
        print('P:', datetime.now().isoformat(timespec='milliseconds'),
              '--created state', flush=True)

    # Instantiate DynamicalCore
    dycore = fv3core.DynamicalCore(
        communicator,
        grid_data=spec.grid.grid_data,
        grid_indexing=spec.grid.grid_indexing,
        damping_coefficients=spec.grid.damping_coefficients,
        config=spec.namelist.dynamical_core,
        ak=state["atmosphere_hybrid_a_coordinate"],
        bk=state["atmosphere_hybrid_b_coordinate"],
        phis=state["surface_geopotential"],)
    if rank == 0:
        print('P:', datetime.now().isoformat(timespec='milliseconds'),
              '--instantiated DynamicalCore', flush=True)

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
        # hf = h5py.File('output.h5', 'w')
        # for var in gtfv3_output_data:
        #     hf.create_dataset(var, data=gtfv3_output_data[var])
        # hf.close()

    nranks = comm.Get_size()
    for i in range(nranks):
        if i == rank:
            print('P: rank:', rank, flush=True)
            print('P: u:',
                  np.sum(gtfv3_output_data['u']),
                  np.sum(gtfv3_output_data['v']),
                  np.sum(gtfv3_output_data['w']),
                  np.sum(gtfv3_output_data['delz']),
                  flush=True)
            print('P: pt:',
                  np.sum(gtfv3_output_data['pt']),
                  np.sum(gtfv3_output_data['delp']),
                  np.sum(gtfv3_output_data['qvapor']),
                  np.sum(gtfv3_output_data['qliquid']),
                  np.sum(gtfv3_output_data['qice']),
                  np.sum(gtfv3_output_data['qrain']),
                  np.sum(gtfv3_output_data['qsnow']),
                  np.sum(gtfv3_output_data['qgraupel']),
                  np.sum(gtfv3_output_data['qcld']),
                  flush=True)
            print('P: ps:',
                  np.sum(gtfv3_output_data['ps']),
                  np.sum(gtfv3_output_data['pe']),
                  np.sum(gtfv3_output_data['pk']),
                  np.sum(gtfv3_output_data['peln']),
                  np.sum(gtfv3_output_data['pkz']),
                  flush=True)
            print('P: phis:',
                  np.sum(gtfv3_output_data['phis']),
                  np.sum(gtfv3_output_data['q_con']),
                  np.sum(gtfv3_output_data['omga']),
                  flush=True)
            print('P: ua:',
                  np.sum(gtfv3_output_data['ua']),
                  np.sum(gtfv3_output_data['va']),
                  np.sum(gtfv3_output_data['uc']),
                  np.sum(gtfv3_output_data['vc']),
                  flush=True)
            print('P: mfx:',
                  np.sum(gtfv3_output_data['mfxd']),
                  np.sum(gtfv3_output_data['mfyd']),
                  np.sum(gtfv3_output_data['cxd']),
                  np.sum(gtfv3_output_data['cyd']),
                  flush=True)
            print('P: diss_est:', np.sum(gtfv3_output_data['diss_estd']), flush=True)
        comm.Barrier()

    return gtfv3_output_data
