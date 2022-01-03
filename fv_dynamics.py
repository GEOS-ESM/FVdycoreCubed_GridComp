from datetime import datetime

import cffi
import numpy as np
from mpi4py import MPI
import h5py

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
        mfx_ptr, mfy_ptr, cx_ptr, cy_ptr, diss_est_ptr,
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
    backend = 'gtx86'
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
    # print('delz:', delz[6,12,13], delz[7,12,13], delz[6,13,13])
    # with np.printoptions(formatter={'float': lambda x: "{0:06.2f}".format(x)}, linewidth=10000):
    # print('u, v, w, delz:', np.sum(u), np.sum(v), np.sum(w), np.sum(delz))
    pt = fort_to_numpy(ffi, pt_ptr, (ied-isd+1, jed-jsd+1, npz))
    delp = fort_to_numpy(ffi, delp_ptr, (ied-isd+1, jed-jsd+1, npz))
    # print('delp:', delp[6,12,13], delp[7,12,13], delp[6,13,13])
    q = fort_to_numpy(ffi, q_ptr, (ied-isd+1, jed-jsd+1, npz, ncnst))
    # print('pt, delp, q:', np.sum(pt), np.sum(delp), np.sum(q))
    ps = fort_to_numpy(ffi, ps_ptr, (ied-isd+1, jed-jsd+1))
    pe = fort_to_numpy(ffi, pe_ptr, (ie+1-(is_-1)+1, npz+1, je+1-(js-1)+1))
    # pe = fort_to_numpy(ffi, pe_ptr, (ie+1-(is_-1)+1, je+1-(js-1)+1, npz+1))
    # if (rank == 0): print('shape(pe):', pe.shape)
    pk = fort_to_numpy(ffi, pk_ptr, (ie-is_+1, je-js+1, npz+1))
    peln = fort_to_numpy(ffi, peln_ptr, (ie-is_+1, npz+1, je-js+1))
    # peln = fort_to_numpy(ffi, peln_ptr, (ie-is_+1, je-js+1, npz+1))
    pkz = fort_to_numpy(ffi, pkz_ptr, (ie-is_+1, je-js+1, npz))
    # print('ps/pe/pk/peln/pkz:', np.sum(ps), np.sum(pe), np.sum(pk), np.sum(peln), np.sum(pkz))
    # print(ps.shape, pe.shape, pk.shape, peln.shape, pkz.shape)
    phis = fort_to_numpy(ffi, phis_ptr, (ied-isd+1, jed-jsd+1))
    # with np.printoptions(formatter={'float': lambda x: "{0:013.3f}".format(x)}, linewidth=10000):
    #      print(phis)
    q_con = fort_to_numpy(ffi, q_con_ptr, (ied-isd+1, jed-jsd+1, npz))
    # print('q_con:', q_con.shape, np.sum(q_con))
    omga = fort_to_numpy(ffi, omga_ptr, (ied-isd+1, jed-jsd+1, npz))
    # print('phis, q_con, omga:', np.sum(phis), np.sum(q_con), np.sum(omga))
    ua = fort_to_numpy(ffi, ua_ptr, (ied-isd+1, jed-jsd+1, npz))
    va = fort_to_numpy(ffi, va_ptr, (ied-isd+1, jed-jsd+1, npz))
    uc = fort_to_numpy(ffi, uc_ptr, (ied+1-isd+1, jed-jsd+1, npz))
    vc = fort_to_numpy(ffi, va_ptr, (ied-isd+1, jed+1-jsd+1, npz))
    # print('ua, va, uc, vc:', np.sum(ua), np.sum(va), np.sum(uc), np.sum(vc))
    ak = fort_to_numpy(ffi, ak_ptr, (npz+1,))
    bk = fort_to_numpy(ffi, bk_ptr, (npz+1,))
    # print('ak, bk:', np.sum(ak), np.sum(bk))
    mfx = fort_to_numpy(ffi, mfx_ptr, (ie+1-is_+1, je-js+1, npz))
    mfy = fort_to_numpy(ffi, mfy_ptr, (ie-is_+1, je+1-js+1, npz))
    cx = fort_to_numpy(ffi, cx_ptr, (ie+1-is_+1, jed-jsd+1, npz))
    cy = fort_to_numpy(ffi, cy_ptr, (ied-isd+1, je+1-js+1, npz))
    # print('mfx, mfy, cx, cy:', np.sum(mfx), np.sum(mfy), np.sum(cx), np.sum(cy))
    # ze0
    diss_est = fort_to_numpy(ffi, diss_est_ptr, (ied-isd+1, jed-jsd+1, npz))
    # print('diss_est:', np.sum(diss_est))
    dx = fort_to_numpy(ffi, dx_ptr, (ied-isd+1, jed+1-jsd+1))
    dy = fort_to_numpy(ffi, dy_ptr, (ied+1-isd+1, jed-jsd+1))
    dxa = fort_to_numpy(ffi, dxa_ptr, (ied-isd+1, jed-jsd+1))
    dya = fort_to_numpy(ffi, dya_ptr, (ied-isd+1, jed-jsd+1))
    dxc = fort_to_numpy(ffi, dxc_ptr, (ied+1-isd+1, jed-jsd+1))
    dyc = fort_to_numpy(ffi, dyc_ptr, (ied-isd+1, jed+1-jsd+1))
    # print('dx, dy, dxa, dya, dxc, dyc:',
    #       np.sum(dx), np.sum(dy), np.sum(dxa), np.sum(dya), np.sum(dxc), np.sum(dyc))
    rdx = fort_to_numpy(ffi, rdx_ptr, (ied-isd+1, jed+1-jsd+1))
    rdy = fort_to_numpy(ffi, rdy_ptr, (ied+1-isd+1, jed-jsd+1))
    rdxa = fort_to_numpy(ffi, rdxa_ptr, (ied-isd+1, jed-jsd+1))
    rdya = fort_to_numpy(ffi, rdya_ptr, (ied-isd+1, jed-jsd+1))
    rdxc = fort_to_numpy(ffi, rdxc_ptr, (ied+1-isd+1, jed-jsd+1))
    rdyc = fort_to_numpy(ffi, rdyc_ptr, (ied-isd+1, jed+1-jsd+1))
    # print('rdx, rdy, rdxa, rdya, rdxc, rdyc:',
    #       np.sum(rdx), np.sum(rdy), np.sum(rdxa), np.sum(rdya), np.sum(rdxc), np.sum(rdyc))
    cosa = fort_to_numpy(ffi, cosa_ptr, (ied+1-isd+1, jed+1-jsd+1))
    cosa_s = fort_to_numpy(ffi, cosa_s_ptr, (ied-isd+1, jed-jsd+1))
    # print('cosa, cosa_s:', np.sum(cosa), np.sum(cosa_s))
    sina_u = fort_to_numpy(ffi, sina_u_ptr, (ied+1-isd+1, jed-jsd+1))
    sina_v = fort_to_numpy(ffi, sina_v_ptr, (ied-isd+1, jed+1-jsd+1))
    cosa_u = fort_to_numpy(ffi, cosa_u_ptr, (ied+1-isd+1, jed-jsd+1))
    cosa_v = fort_to_numpy(ffi, cosa_v_ptr, (ied-isd+1, jed+1-jsd+1))
    # print('sina_u, sina_v, cosa_u, cosa_v:',
    #       np.sum(sina_u), np.sum(sina_v), np.sum(cosa_u), np.sum(cosa_v))
    rsin2 = fort_to_numpy(ffi, rsin2_ptr, (ied-isd+1, jed-jsd+1))
    rsina = fort_to_numpy(ffi, rsina_ptr, (ied+1-isd+1, jed+1-jsd+1))
    rsin_u = fort_to_numpy(ffi, rsin_u_ptr, (ied+1-isd+1, jed-jsd+1))
    rsin_v = fort_to_numpy(ffi, rsin_v_ptr, (ied-isd+1, jed+1-jsd+1))
    sin_sg = fort_to_numpy(ffi, sin_sg_ptr, (ied-isd+1, jed-jsd+1, 9))
    cos_sg = fort_to_numpy(ffi, cos_sg_ptr, (ied-isd+1, jed-jsd+1, 9))
    # print('rsin2, rsina, rsin_u, rsin_v, sin_sg, cos_sg:',
    #       np.sum(rsin2), np.sum(rsina), np.sum(rsin_u),
    #       np.sum(rsin_v), np.sum(sin_sg), np.sum(cos_sg))
    area = fort_to_numpy(ffi, area_ptr, (ied-isd+1, jed-jsd+1))
    rarea = fort_to_numpy(ffi, rarea_ptr, (ied-isd+1, jed-jsd+1))
    rarea_c =  fort_to_numpy(ffi, rarea_c_ptr, (ied+1-isd+1,jed+1-jsd+1))
    f0 = fort_to_numpy(ffi, f0_ptr, (ied-isd+1, jed-jsd+1))
    fC = fort_to_numpy(ffi, fC_ptr, (ied+1-isd+1, jed+1-jsd+1))
    # print('area, rarea, rarea_c, f0, fC:',
    #       np.sum(area), np.sum(rarea), np.sum(rarea_c), np.sum(f0), np.sum(fC))
    del6_u = fort_to_numpy(ffi, del6_u_ptr, (ied-isd+1, jed+1-jsd+1))
    del6_v = fort_to_numpy(ffi, del6_v_ptr, (ied+1-isd+1, jed-jsd+1))
    divg_u = fort_to_numpy(ffi, divg_u_ptr, (ied-isd+1, jed+1-jsd+1))
    divg_v = fort_to_numpy(ffi, divg_v_ptr, (ied+1-isd+1, jed-jsd+1))
    # print('del6_u, del6_v, divg_u, divg_v:',
    #       np.sum(del6_u), np.sum(del6_v), np.sum(divg_u), np.sum(divg_v))
    agrid = fort_to_numpy(ffi, agrid_ptr, (ied-isd+1, jed-jsd+1, 2))
    bgrid = fort_to_numpy(ffi, bgrid_ptr, (ied+1-isd+1, jed+1-jsd+1, 2))
    # print('agrid, bgrid:', np.sum(agrid), np.sum(bgrid))
    a11 = fort_to_numpy(ffi, a11_ptr, (ie+1-(is_-1)+1, je+1-(js-1)+1))
    a12 = fort_to_numpy(ffi, a12_ptr, (ie+1-(is_-1)+1, je+1-(js-1)+1))
    a21 = fort_to_numpy(ffi, a21_ptr, (ie+1-(is_-1)+1, je+1-(js-1)+1))
    a22 = fort_to_numpy(ffi, a22_ptr, (ie+1-(is_-1)+1, je+1-(js-1)+1))
    # print('a11, a12, a21, a22:', np.sum(a11), np.sum(a12), np.sum(a21), np.sum(a22))
    edge_e = fort_to_numpy(ffi, edge_e_ptr, (npy,))
    edge_w = fort_to_numpy(ffi, edge_w_ptr, (npy,))
    edge_n = fort_to_numpy(ffi, edge_n_ptr, (npx,))
    edge_s = fort_to_numpy(ffi, edge_s_ptr, (npx,))
    # print('edge_e/w/n/s:', np.sum(edge_e), np.sum(edge_w), np.sum(edge_n), np.sum(edge_s))
    # print('nested, stretched, da_min, da_min_c:', nested_int, stretched_grid_int,da_min,da_min_c)

    if (rank == 0):
        print('P:', datetime.now().isoformat(timespec='milliseconds'),
              '--converted data (Fortran to NumPy)', flush=True)
        hf = h5py.File('incoming-data.h5', 'w')
        hf.create_dataset('u', data=u)
        hf.create_dataset('v', data=v)
        hf.create_dataset('w', data=w)
        hf.create_dataset('delz', data=delz)
        hf.create_dataset('pt', data=pt)
        hf.create_dataset('delp', data=delp)
        hf.create_dataset('q', data=q)
        hf.create_dataset('ps', data=ps)
        hf.create_dataset('pe', data=pe)
        hf.create_dataset('pk', data=pk)
        hf.create_dataset('peln', data=peln)
        hf.create_dataset('pkz', data=pkz)
        hf.create_dataset('phis', data=phis)
        hf.create_dataset('q_con', data=q_con)
        hf.create_dataset('omga', data=omga)
        hf.create_dataset('ak', data=ak)
        hf.create_dataset('bk', data=bk)
        hf.close()

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
        'kord_tm': -9, #kord_tm,
        'kord_tr': kord_tr,
        'kord_wz': kord_wz,
        'kord_mt': kord_mt,
        'd_ext': d_ext,
        'beta': beta,
        # To get around the error:
        # File "/discover/nobackup/pchakrab/code/ext/vulcan/fv3core/fv3core/stencils/updatedzd.py", line 207, in __init__
        # raise NotImplementedError("damp <= 1e-5 in column_cols is untested")
        # NotImplementedError: damp <= 1e-5 in column_cols is untested
        'vtdm4': vtdm4,
        'ke_bg': ke_bg,
        'd_con': d_con,
        'd2_bg': d2_bg,
        'd2_bg_k1': d2_bg_k1,
        'd2_bg_k2': d2_bg_k2,
        'p_fac': p_fac,
        'a_imp': a_imp,
        'dddmp': 0.2,
        'd4_bg': 0.12,
        'rf_cutoff': 750.0,
        'tau': 0.0,
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
        'mfxd': mfx, 'mfyd': mfy, 'cxd': cx, 'cyd': cy,
        'diss_estd': diss_est,
        'qvapor': q[:,:,:,0],
        'qliquid': q[:,:,:,1],
        'qice': q[:,:,:,2],
        'qrain': q[:,:,:,3],
        'qsnow': q[:,:,:,4],
        'qgraupel': q[:,:,:,5],
        'qcld': q[:,:,:,6],
        'qo3mr': q[:,:,:,7],
        'qsgs_tke': q[:,:,:,8],
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
        # 'a11': np.zeros((ie+1-(is_-1)+1, je+1-(js-1)+1), np.float64, order='F'),
        # 'a12': np.zeros((ie+1-(is_-1)+1, je+1-(js-1)+1), np.float64, order='F'),
        # 'a21': np.zeros((ie+1-(is_-1)+1, je+1-(js-1)+1), np.float64, order='F'),
        # 'a22': np.zeros((ie+1-(is_-1)+1, je+1-(js-1)+1), np.float64, order='F'),
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
        hf = h5py.File('grid-data.h5', 'w')
        for var in grid_vars_to_write:
            # print('var:', var, grid.__dict__[var].shape, np.sum(grid.__dict__[var]))
            hf.create_dataset(var, data=grid.__dict__[var])
        hf.close()
    driver_object = fv3core.testing.TranslateFVDynamics([grid])
    state = driver_object.state_from_inputs(input_data)
    # print('state.pe/state.q_con:', state['pe'].shape, state['q_con'].shape)
    if (spec.grid.rank == 0):
        print('P:', datetime.now().isoformat(timespec='milliseconds'),
              '--created state', flush=True)
        # print(spec.namelist.dynamical_core, flush=True)
    vars_to_write = ['u', 'v', 'w', 'delz',
                     'pt', 'delp',
                     'ps', 'pe', 'pk', 'peln', 'pkz',
                     'phis', 'q_con', 'omga',
                     'ua', 'va', 'uc', 'vc']
    if (rank == 0):
        hf = h5py.File('state-initial.h5', 'w')
        for var in vars_to_write:
            hf.create_dataset(var, data=state[var])
        hf.close()

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

    if spec.namelist.fv_sg_adj > 0:
        fv_subgrid_z = fv3core.DryConvectiveAdjustment(
            spec.grid.grid_indexing,
            spec.namelist.nwat,
            spec.namelist.fv_sg_adj,
            spec.namelist.n_sponge,
            spec.namelist.hydrostatic,)
        if (spec.grid.rank == 0):
            print('P:', datetime.now().isoformat(timespec='milliseconds'),
                  '--instantiated DryConvectiveAdjustment', flush=True)

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

    if (rank == 0):
        hf = h5py.File('state-final.h5', 'w')
        for var in vars_to_write:
            hf.create_dataset(var, data=state[var])
        hf.close()
