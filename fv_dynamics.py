from datetime import datetime

import cffi
import numpy as np
from mpi4py import MPI
# import h5py

import gt4py
from fort2py import fort_to_numpy

import fv3core
import fv3core._config as spec
import fv3core.testing
import fv3gfs.util as util

def fv_dynamics_top_level_function(
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
        u_ptr, v_ptr, w_ptr, delz_ptr,
        pt_ptr, delp_ptr, q_ptr,
        ps_ptr, pe_ptr, pk_ptr, peln_ptr, pkz_ptr,
        phis_ptr, q_con_ptr, omga_ptr,
        ua_ptr, va_ptr, uc_ptr, vc_ptr,
        ak_ptr, bk_ptr,
        mfxd_ptr, mfyd_ptr, cxd_ptr, cyd_ptr, diss_estd_ptr,
        dx_ptr, dy_ptr, dxa_ptr, dya_ptr, dxc_ptr, dyc_ptr,
        rdx_ptr, rdy_ptr, rdxa_ptr, rdya_ptr, rdxc_ptr, rdyc_ptr,
        cosa_ptr, cosa_s_ptr, sina_u_ptr, sina_v_ptr,
        cosa_u_ptr, cosa_v_ptr, rsin2_ptr, rsina_ptr, rsin_u_ptr, rsin_v_ptr,
        sin_sg_ptr, cos_sg_ptr,
        area_ptr, rarea_ptr, rarea_c_ptr, f0_ptr, fC_ptr,
        del6_u_ptr, del6_v_ptr, divg_u_ptr, divg_v_ptr,
        agrid_ptr, bgrid_ptr, a11_ptr, a12_ptr, a21_ptr, a22_ptr,
        edge_e_ptr, edge_w_ptr, edge_n_ptr, edge_s_ptr,
        nested_int, stretched_grid_int, da_min, da_min_c):

    rank = comm.Get_rank()
    backend = 'gtcuda'
    if (rank == 0):
        print('P:', datetime.now().isoformat(timespec='milliseconds'),
              '--in top level function with backend', backend, flush=True)

    origin = (ng, ng, 0) # fv3core.utils.gt4py_utils.origin
    origin0 = (0, 0, 0)

    # Convert Fortran arrays to NumPy
    ffi = cffi.FFI()
    u = fort_to_numpy(ffi, u_ptr, (ied-isd+1, jed+1-jsd+1, npz))
    v = fort_to_numpy(ffi, v_ptr, (ied+1-isd+1, jed-jsd+1, npz))
    w = fort_to_numpy(ffi, w_ptr, (ied-isd+1, jed-jsd+1, npz))
    delz = fort_to_numpy(ffi, delz_ptr, (ied-isd+1, jed-jsd+1, npz))
    pt = fort_to_numpy(ffi, pt_ptr, (ied-isd+1, jed-jsd+1, npz))
    delp = fort_to_numpy(ffi, delp_ptr, (ied-isd+1, jed-jsd+1, npz))
    q = fort_to_numpy(ffi, q_ptr, (ied-isd+1, jed-jsd+1, npz, nq_tot))
    ps = fort_to_numpy(ffi, ps_ptr, (ied-isd+1, jed-jsd+1))
    pe = fort_to_numpy(ffi, pe_ptr, (ie+1-(is_-1)+1, npz+1, je+1-(js-1)+1))
    pk = fort_to_numpy(ffi, pk_ptr, (ie-is_+1, je-js+1, npz+1))
    peln = fort_to_numpy(ffi, peln_ptr, (ie-is_+1, npz+1, je-js+1))
    pkz = fort_to_numpy(ffi, pkz_ptr, (ie-is_+1, je-js+1, npz))
    phis = fort_to_numpy(ffi, phis_ptr, (ied-isd+1, jed-jsd+1))
    q_con = fort_to_numpy(ffi, q_con_ptr, (ied-isd+1, jed-jsd+1, npz))
    omga = fort_to_numpy(ffi, omga_ptr, (ied-isd+1, jed-jsd+1, npz))
    ua = fort_to_numpy(ffi, ua_ptr, (ied-isd+1, jed-jsd+1, npz))
    va = fort_to_numpy(ffi, va_ptr, (ied-isd+1, jed-jsd+1, npz))
    uc = fort_to_numpy(ffi, uc_ptr, (ied+1-isd+1, jed-jsd+1, npz))
    vc = fort_to_numpy(ffi, va_ptr, (ied-isd+1, jed+1-jsd+1, npz))
    ak = fort_to_numpy(ffi, ak_ptr, (npz+1,))
    bk = fort_to_numpy(ffi, bk_ptr, (npz+1,))
    mfxd = fort_to_numpy(ffi, mfxd_ptr, (ie+1-is_+1, je-js+1, npz))
    mfyd = fort_to_numpy(ffi, mfyd_ptr, (ie-is_+1, je+1-js+1, npz))
    cxd = fort_to_numpy(ffi, cxd_ptr, (ie+1-is_+1, jed-jsd+1, npz))
    cyd = fort_to_numpy(ffi, cyd_ptr, (ied-isd+1, je+1-js+1, npz))
    diss_estd = fort_to_numpy(ffi, diss_estd_ptr, (ied-isd+1, jed-jsd+1, npz))

    dx = fort_to_numpy(ffi, dx_ptr, (ied-isd+1, jed+1-jsd+1))
    dy = fort_to_numpy(ffi, dy_ptr, (ied+1-isd+1, jed-jsd+1))
    dxa = fort_to_numpy(ffi, dxa_ptr, (ied-isd+1, jed-jsd+1))
    dya = fort_to_numpy(ffi, dya_ptr, (ied-isd+1, jed-jsd+1))
    dxc = fort_to_numpy(ffi, dxc_ptr, (ied+1-isd+1, jed-jsd+1))
    dyc = fort_to_numpy(ffi, dyc_ptr, (ied-isd+1, jed+1-jsd+1))
    rdx = fort_to_numpy(ffi, rdx_ptr, (ied-isd+1, jed+1-jsd+1))
    rdy = fort_to_numpy(ffi, rdy_ptr, (ied+1-isd+1, jed-jsd+1))
    rdxa = fort_to_numpy(ffi, rdxa_ptr, (ied-isd+1, jed-jsd+1))
    rdya = fort_to_numpy(ffi, rdya_ptr, (ied-isd+1, jed-jsd+1))
    rdxc = fort_to_numpy(ffi, rdxc_ptr, (ied+1-isd+1, jed-jsd+1))
    rdyc = fort_to_numpy(ffi, rdyc_ptr, (ied-isd+1, jed+1-jsd+1))
    cosa = fort_to_numpy(ffi, cosa_ptr, (ied+1-isd+1, jed+1-jsd+1))
    cosa_s = fort_to_numpy(ffi, cosa_s_ptr, (ied-isd+1, jed-jsd+1))
    sina_u = fort_to_numpy(ffi, sina_u_ptr, (ied+1-isd+1, jed-jsd+1))
    sina_v = fort_to_numpy(ffi, sina_v_ptr, (ied-isd+1, jed+1-jsd+1))
    cosa_u = fort_to_numpy(ffi, cosa_u_ptr, (ied+1-isd+1, jed-jsd+1))
    cosa_v = fort_to_numpy(ffi, cosa_v_ptr, (ied-isd+1, jed+1-jsd+1))
    rsin2 = fort_to_numpy(ffi, rsin2_ptr, (ied-isd+1, jed-jsd+1))
    rsina = fort_to_numpy(ffi, rsina_ptr, (ie+1-is_+1, je+1-js+1))
    rsin_u = fort_to_numpy(ffi, rsin_u_ptr, (ied+1-isd+1, jed-jsd+1))
    rsin_v = fort_to_numpy(ffi, rsin_v_ptr, (ied-isd+1, jed+1-jsd+1))
    sin_sg = fort_to_numpy(ffi, sin_sg_ptr, (ied-isd+1, jed-jsd+1, 9))
    cos_sg = fort_to_numpy(ffi, cos_sg_ptr, (ied-isd+1, jed-jsd+1, 9))
    area = fort_to_numpy(ffi, area_ptr, (ied-isd+1, jed-jsd+1))
    rarea = fort_to_numpy(ffi, rarea_ptr, (ied-isd+1, jed-jsd+1))
    rarea_c =  fort_to_numpy(ffi, rarea_c_ptr, (ied+1-isd+1,jed+1-jsd+1))
    f0 = fort_to_numpy(ffi, f0_ptr, (ied-isd+1, jed-jsd+1))
    fC = fort_to_numpy(ffi, fC_ptr, (ied+1-isd+1, jed+1-jsd+1))
    del6_u = fort_to_numpy(ffi, del6_u_ptr, (ied-isd+1, jed+1-jsd+1))
    del6_v = fort_to_numpy(ffi, del6_v_ptr, (ied+1-isd+1, jed-jsd+1))
    divg_u = fort_to_numpy(ffi, divg_u_ptr, (ied-isd+1, jed+1-jsd+1))
    divg_v = fort_to_numpy(ffi, divg_v_ptr, (ied+1-isd+1, jed-jsd+1))
    agrid = fort_to_numpy(ffi, agrid_ptr, (ied-isd+1, jed-jsd+1, 2))
    bgrid = fort_to_numpy(ffi, bgrid_ptr, (ied+1-isd+1, jed+1-jsd+1, 2))
    a11 = fort_to_numpy(ffi, a11_ptr, (ie+1-(is_-1)+1, je+1-(js-1)+1))
    a12 = fort_to_numpy(ffi, a12_ptr, (ie+1-(is_-1)+1, je+1-(js-1)+1))
    a21 = fort_to_numpy(ffi, a21_ptr, (ie+1-(is_-1)+1, je+1-(js-1)+1))
    a22 = fort_to_numpy(ffi, a22_ptr, (ie+1-(is_-1)+1, je+1-(js-1)+1))
    edge_e = fort_to_numpy(ffi, edge_e_ptr, (npy,))
    edge_w = fort_to_numpy(ffi, edge_w_ptr, (npy,))
    edge_n = fort_to_numpy(ffi, edge_n_ptr, (npx,))
    edge_s = fort_to_numpy(ffi, edge_s_ptr, (npx,))

    # if rank == 0:
    #     with np.printoptions(formatter={'float': lambda x: "{0:06.2f}".format(x)}, linewidth=10000):
    #         print('P: u/v/w/delz:', np.sum(u), np.sum(v), np.sum(w), np.sum(delz))
    nranks = comm.Get_size()
    for i in range(nranks):
        if i == rank:
            print('P: rank:', rank, flush=True)
            print('P: u:', np.sum(u), np.sum(v), np.sum(w), np.sum(delz), flush=True)
            print('P: pt:', np.sum(pt), np.sum(delp), np.sum(q), flush=True)
            print('P: ps:', np.sum(ps), np.sum(pe), np.sum(pk), np.sum(peln), np.sum(pkz), flush=True)
            print('P: phis:', np.sum(phis), np.sum(q_con), np.sum(omga), flush=True)
            print('P: ua:', np.sum(ua), np.sum(va), np.sum(uc), np.sum(vc), flush=True)
            print('P: mfx:', np.sum(mfxd), np.sum(mfyd), np.sum(cxd), np.sum(cyd), flush=True)
            print('P: diss_est:', np.sum(diss_estd), flush=True)
        comm.Barrier()

    if (rank == 0):
        print('P:', datetime.now().isoformat(timespec='milliseconds'),
              '--converted data (Fortran to NumPy)', flush=True)
        # hf = h5py.File('incoming-data.h5', 'w')
        # hf.create_dataset('u', data=u)
        # hf.create_dataset('v', data=v)
        # hf.create_dataset('w', data=w)
        # hf.create_dataset('delz', data=delz)
        # hf.create_dataset('pt', data=pt)
        # hf.create_dataset('delp', data=delp)
        # hf.create_dataset('q', data=q)
        # hf.create_dataset('ps', data=ps)
        # hf.create_dataset('pe', data=pe)
        # hf.create_dataset('pk', data=pk)
        # hf.create_dataset('peln', data=peln)
        # hf.create_dataset('pkz', data=pkz)
        # hf.create_dataset('phis', data=phis)
        # hf.create_dataset('q_con', data=q_con)
        # hf.create_dataset('omga', data=omga)
        # hf.create_dataset('ak', data=ak)
        # hf.create_dataset('bk', data=bk)
        # hf.close()

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

    # Grid dependent helper structure
    layout = spec.namelist.layout
    partitioner = util.CubedSpherePartitioner(util.TilePartitioner(layout))
    communicator = util.CubedSphereCommunicator(comm, partitioner)

    # Create a dict of input data
    input_data = {
        'u': u, 'v': v, 'w': w, 'delz': delz,
        'pt': pt, 'delp': delp, #'q': q,
        'ps': ps, 'pe': pe, 'pk': pk, 'pkz': pkz, 'peln': peln,
        'phis': phis, 'q_con': q_con, 'omga': omga,
        'ua': ua, 'va': va, 'uc': uc, 'vc': vc,
        'ak': ak, 'bk': bk,
        'mfxd': mfxd, 'mfyd': mfyd, 'cxd': cxd, 'cyd': cyd,
        'diss_estd': diss_estd,
        'qvapor': q[:,:,:,0],
        'qliquid': q[:,:,:,1],
        'qice': q[:,:,:,2],
        'qrain': q[:,:,:,3],
        'qsnow': q[:,:,:,4],
        'qgraupel': q[:,:,:,5],
        'qcld': q[:,:,:,6],
        # GEOS does not provide qo3mr and qsgs_tke - setting those to zero
        'qo3mr': np.zeros(q[:,:,:,0].shape),
        'qsgs_tke': np.zeros(q[:,:,:,0].shape),
        'do_adiabatic_init': True, # Has to be True
        'consv_te': consv_te,
        'bdt': bdt,
        'ptop': ptop,
        'n_split': n_split, # Need n_split here as well
        'ks': ks,
        'comm': communicator} # input_data
    # print('cxd/cyd shapes (input):', input_data['cxd'].shape, input_data['cyd'].shape)
    # print('u/v/w   shapes (input):',
    #       input_data['u'].shape, input_data['v'].shape, input_data['w'].shape)

    grid_data = {
        'npx': npx, 'npy': npy, 'npz': npz,
        'is_': is_, 'ie': ie, 'js': js, 'je': je,
        'isd': isd, 'ied': ied, 'jsd': jsd, 'jed': jed,
        'dx': dx, 'dy': dy, 'dxa': dxa, 'dya': dya, 'dxc': dxc, 'dyc': dyc,
        'rdx': rdx, 'rdy': rdy, 'rdxa': rdxa, 'rdya': rdya, 'rdxc': rdxc, 'rdyc': rdyc,
        'cosa' : cosa, 'cosa_s': cosa_s,
        'sina_u': sina_u, 'sina_v': sina_v, 'cosa_u': cosa_u, 'cosa_v': cosa_v,
        'rsin2': rsin2, 'rsina': rsina, 'rsin_u': rsin_u, 'rsin_v': rsin_v,
        'sin_sg1': sin_sg[:,:,0], 'sin_sg2': sin_sg[:,:,1],
        'sin_sg3': sin_sg[:,:,2], 'sin_sg4': sin_sg[:,:,3],
        'cos_sg1': cos_sg[:,:,0], 'cos_sg2': cos_sg[:,:,1],
        'cos_sg3': cos_sg[:,:,2], 'cos_sg4': cos_sg[:,:,3],
        'area': area, 'area_64': area, 'rarea': rarea, 'rarea_c': rarea_c,
        'f0': f0, 'fC': fC,
        'del6_u': del6_u, 'del6_v': del6_v, 'divg_u': divg_u, 'divg_v': divg_v,
        'agrid1': agrid[:,:,0], 'agrid2': agrid[:,:,1],
        'bgrid1': bgrid[:,:,0], 'bgrid2': bgrid[:,:,1],
        'a11': a11, 'a12': a12, 'a21': a21, 'a22': a22,
        'edge_e': edge_e, 'edge_w': edge_w, 'edge_n': edge_n, 'edge_s': edge_s,
        'nested': False if nested_int == 0 else True,
        'stretched_grid': False if stretched_grid_int == 0 else True,
        'da_min': da_min, 'da_min_c': 0.0,} # grid_data

    if (rank == 0):
        print('P:', datetime.now().isoformat(timespec='milliseconds'),
              '--wrapped data', flush=True)

    grid = fv3core.testing.TranslateGrid(grid_data, comm.Get_rank()).python_grid()
    if (rank == 0):
        print('P:', datetime.now().isoformat(timespec='milliseconds'),
              '--translated data', flush=True)
    spec.set_grid(grid)
    grid_vars_to_write = [
        'dx', 'dy', 'dxa', 'dya', 'dxc', 'dyc',
        'rdx', 'rdy', 'rdxa', 'rdya', 'rdxc', 'rdyc',
        'cosa', 'cosa_s', 'sina_u', 'sina_v', 'cosa_u', 'cosa_v',
        'rsin2', 'rsina', 'rsin_u', 'rsin_v',
        'sin_sg1', 'sin_sg2', 'sin_sg3', 'sin_sg4',
        'cos_sg1', 'cos_sg2', 'cos_sg3', 'cos_sg4',
        'area', 'rarea', 'rarea_c', 'f0', 'fC',
        'del6_u', 'del6_v', 'divg_u', 'divg_v',
        'agrid1', 'agrid2', 'bgrid1', 'bgrid2',
        'a11', 'a12', 'a21', 'a22',
        'edge_e', 'edge_w', 'edge_n', 'edge_s']

    if (spec.grid.rank == 0):
        print('P:', datetime.now().isoformat(timespec='milliseconds'),
              '--created grid', flush=True)
        # hf = h5py.File('grid-data.h5', 'w')
        # for var in grid_vars_to_write:
        #     # print('var:', var, grid.__dict__[var].shape, np.sum(grid.__dict__[var]))
        #     hf.create_dataset(var, data=grid.__dict__[var])
        # hf.close()
    driver_object = fv3core.testing.TranslateFVDynamics([grid])
    state = driver_object.state_from_inputs(input_data)
    # print('state.pe/state.q_con:', state['pe'].shape, state['q_con'].shape)
    if (spec.grid.rank == 0):
        print('P:', datetime.now().isoformat(timespec='milliseconds'),
              '--created state', flush=True)
        print(spec.namelist.dynamical_core, flush=True)
    vars_to_write = ['u', 'v', 'w', 'delz',
                     'pt', 'delp',
                     'ps', 'pe', 'pk', 'peln', 'pkz',
                     'phis', 'q_con', 'omga',
                     'ua', 'va', 'uc', 'vc']
    # if (rank == 0):
    #     hf = h5py.File('state-initial.h5', 'w')
    #     for var in vars_to_write:
    #         hf.create_dataset(var, data=state[var])
    #     hf.close()

    # print('pe after  translating:', state['pe'].shape, np.sum(state['pe']))
    # print(state.keys())
    # print('cxd/cyd shapes (state):', state['cxd'].shape, state['cyd'].shape)
    # print('u/v/w   shapes (state):', state['u'].shape, state['v'].shape, state['w'].shape)
    dycore = fv3core.DynamicalCore(
        communicator,
        grid_data=spec.grid.grid_data,
        grid_indexing=spec.grid.grid_indexing,
        damping_coefficients=spec.grid.damping_coefficients,
        config=spec.namelist.dynamical_core,
        ak=state["atmosphere_hybrid_a_coordinate"],
        bk=state["atmosphere_hybrid_b_coordinate"],
        phis=state["surface_geopotential"],)
    if (spec.grid.rank == 0):
        print('P:', datetime.now().isoformat(timespec='milliseconds'),
              '--instantiated DynamicalCore', flush=True)

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

    # First timestep
    dycore.step_dynamics(
        state,
        input_data['consv_te'],
        input_data['do_adiabatic_init'],
        input_data['bdt'],
        input_data['ptop'],
        input_data['n_split'],
        input_data['ks'])
    if (spec.grid.rank == 0):
        print('P:', datetime.now().isoformat(timespec='milliseconds'),
              '--ran DynamicalCore::step_dynamics', flush=True)

    # Retrieve output data from state
    output_data = driver_object.outputs_from_state(state)

    if (spec.grid.rank == 0):
        print('P:', datetime.now().isoformat(timespec='milliseconds'),
              '--retrieved output data from state', flush=True)

    # q needs special handling - no need to include qo3mr and qsgs_tks
    q[:,:,:,0] = output_data['qvapor']
    q[:,:,:,1] = output_data['qliquid']
    q[:,:,:,2] = output_data['qice']
    q[:,:,:,3] = output_data['qrain']
    q[:,:,:,4] = output_data['qsnow']
    q[:,:,:,5] = output_data['qgraupel']
    q[:,:,:,6] = output_data['qcld']

    nranks = comm.Get_size()
    for i in range(nranks):
        if i == rank:
            print('P: rank:', rank, flush=True)
            print('P: u:',
                  np.sum(output_data['u']),
                  np.sum(output_data['v']),
                  np.sum(output_data['w']),
                  np.sum(output_data['delz']),
                  flush=True)
            print('P: pt:', np.sum(output_data['pt']), np.sum(output_data['delp']), np.sum(q), q.shape, flush=True)
            print('P: ps:',
                  np.sum(output_data['ps']),
                  np.sum(output_data['pe']),
                  np.sum(output_data['pk']),
                  np.sum(output_data['peln']),
                  np.sum(output_data['pkz']),
                  flush=True)
            print('P: phis:',
                  np.sum(output_data['phis']),
                  np.sum(output_data['q_con']),
                  np.sum(output_data['omga']),
                  flush=True)
            print('P: ua:',
                  np.sum(output_data['ua']),
                  np.sum(output_data['va']),
                  np.sum(output_data['uc']),
                  np.sum(output_data['vc']),
                  flush=True)
            print('P: mfx:',
                  np.sum(output_data['mfxd']),
                  np.sum(output_data['mfyd']),
                  np.sum(output_data['cxd']),
                  np.sum(output_data['cyd']),
                  flush=True)
            print('P: diss_est:', np.sum(output_data['diss_estd']), flush=True)
        comm.Barrier()

    # if (rank == 0):
    #     hf = h5py.File('output.h5', 'w')
    #     for var in output_data:
    #         hf.create_dataset(var, data=output_data[var])
    #     hf.close()

    # u_out_f = output_data['u'].astype(np.float32).flatten(order='F')
    # ffi.memmove(u_ptr, u_out_f, 4*u_out_f.size)
    # v_out_f = output_data['v'].astype(np.float32).flatten(order='F')
    # ffi.memmove(v_ptr, v_out_f, 4*v_out_f.size)
    # w_out_f = output_data['w'].astype(np.float32).flatten(order='F')
    # ffi.memmove(w_ptr, w_out_f, 4*w_out_f.size)
    # delz_out_f = output_data['delz'].astype(np.float32).flatten(order='F')
    # ffi.memmove(delz_ptr, delz_out_f, 4*delz_out_f.size)

    # pt_out_f = output_data['pt'].astype(np.float32).flatten(order='F')
    # ffi.memmove(pt_ptr, pt_out_f, 4*pt_out_f.size)
    # delp_out_f = output_data['delp'].astype(np.float32).flatten(order='F')
    # ffi.memmove(delp_ptr, delp_out_f, 4*delp_out_f.size)
    # # q needs special handling - no need to include qo3mr and qsgs_tks
    # q[:,:,:,0] = output_data['qvapor']
    # q[:,:,:,1] = output_data['qliquid']
    # q[:,:,:,2] = output_data['qice']
    # q[:,:,:,3] = output_data['qrain']
    # q[:,:,:,4] = output_data['qsnow']
    # q[:,:,:,5] = output_data['qgraupel']
    # q[:,:,:,6] = output_data['qcld']
    # q_out_f = q.astype(np.float32).flatten(order='F')
    # ffi.memmove(q_ptr, q_out_f, 4*q_out_f.size)

    # ps_out_f = output_data['ps'].astype(np.float32).flatten(order='F')
    # ffi.memmove(ps_ptr, ps_out_f, 4*ps_out_f.size)
    # pe_out_f = output_data['pe'].astype(np.float32).flatten(order='F')
    # ffi.memmove(pe_ptr, pe_out_f, 4*pe_out_f.size)
    # pk_out_f = output_data['pk'].astype(np.float32).flatten(order='F')
    # ffi.memmove(pk_ptr, pk_out_f, 4*pk_out_f.size)
    # peln_out_f = output_data['peln'].astype(np.float32).flatten(order='F')
    # ffi.memmove(peln_ptr, peln_out_f, 4*peln_out_f.size)
    # pkz_out_f = output_data['pkz'].astype(np.float32).flatten(order='F')
    # ffi.memmove(pkz_ptr, pkz_out_f, 4*pkz_out_f.size)

    # phis_out_f = output_data['phis'].astype(np.float32).flatten(order='F')
    # ffi.memmove(phis_ptr, phis_out_f, 4*phis_out_f.size)
    # q_con_out_f = output_data['q_con'].astype(np.float32).flatten(order='F')
    # ffi.memmove(q_con_ptr, q_con_out_f, 4*q_con_out_f.size)
    # omga_out_f = output_data['omga'].astype(np.float32).flatten(order='F')
    # ffi.memmove(omga_ptr, omga_out_f, 4*omga_out_f.size)

    # ua_out_f = output_data['ua'].astype(np.float32).flatten(order='F')
    # ffi.memmove(ua_ptr, ua_out_f, 4*ua_out_f.size)
    # va_out_f = output_data['va'].astype(np.float32).flatten(order='F')
    # ffi.memmove(va_ptr, va_out_f, 4*va_out_f.size)
    # uc_out_f = output_data['uc'].astype(np.float32).flatten(order='F')
    # ffi.memmove(uc_ptr, uc_out_f, 4*uc_out_f.size)
    # vc_out_f = output_data['vc'].astype(np.float32).flatten(order='F')
    # ffi.memmove(vc_ptr, vc_out_f, 4*vc_out_f.size)

    # mfxd_out_f = output_data['mfxd'].astype(np.float32).flatten(order='F')
    # ffi.memmove(mfxd_ptr, mfxd_out_f, 4*mfxd_out_f.size)
    # mfyd_out_f = output_data['mfyd'].astype(np.float32).flatten(order='F')
    # ffi.memmove(mfyd_ptr, mfyd_out_f, 4*mfyd_out_f.size)
    # cxd_out_f = output_data['cxd'].astype(np.float32).flatten(order='F')
    # ffi.memmove(cxd_ptr, cxd_out_f, 4*cxd_out_f.size)
    # cyd_out_f = output_data['cyd'].astype(np.float32).flatten(order='F')
    # ffi.memmove(cyd_ptr, cyd_out_f, 4*cyd_out_f.size)
    # diss_estd_out_f = output_data['diss_estd'].astype(np.float32).flatten(order='F')
    # ffi.memmove(diss_estd_ptr, diss_estd_out_f, 4*diss_estd_out_f.size)

    if (spec.grid.rank == 0):
        print('P:', datetime.now().isoformat(timespec='milliseconds'),
              '--converted data (NumPy to Fortran)', flush=True)
