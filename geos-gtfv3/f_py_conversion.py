import cffi
import numpy as np
from math import prod

ffi = cffi.FFI()
TYPEMAP = {
    'float': np.float32,
    'double': np.float64,
    'int': np.int32,}

def fort_to_numpy(fptr, dim):
    '''
    Input: Fortran data pointed to by fptr and of shape dim = (i, j, k)
    Output: C-ordered double precision NumPy data of shape (i, j, k)
    '''
    ftype = ffi.getctype(ffi.typeof(fptr).item)
    assert ftype in TYPEMAP
    return np.frombuffer(
        ffi.buffer(fptr, prod(dim)*ffi.sizeof(ftype)),
        TYPEMAP[ftype],
    ).reshape(tuple(reversed(dim))).transpose().astype(np.float64)



def fortran_input_data_to_numpy(
        npx, npy, npz,
        is_, ie, js, je, isd, ied, jsd, jed,
        ncnst, consv_te, bdt, ptop, n_split, ks,
        ak_ptr, bk_ptr,
        # input/output arrays
        u_ptr, v_ptr, w_ptr, delz_ptr,
        pt_ptr, delp_ptr, q_ptr,
        ps_ptr, pe_ptr, pk_ptr, peln_ptr, pkz_ptr,
        phis_ptr, q_con_ptr, omga_ptr,
        ua_ptr, va_ptr, uc_ptr, vc_ptr,
        mfxd_ptr, mfyd_ptr, cxd_ptr, cyd_ptr, diss_estd_ptr):
    '''
    Convert Fortran arrays pointed to by *_ptr to NumPy arrays
    Input: Pointers to Fortran arrays *_ptr
    Output: dict where dict[key] is a NumPy array plus consv_te etc.
    '''

    # 'q' requires special handling
    q = fort_to_numpy(q_ptr, (ied-isd+1, jed-jsd+1, npz, ncnst))

    inout_data = {
        'do_adiabatic_init': True, # REQUIRED (TODO???)
        'consv_te': consv_te, 'bdt': bdt, 'ptop': ptop,
        'n_split': n_split, 'ks': ks,

        'ak': fort_to_numpy(ak_ptr, (npz+1,)),
        'bk': fort_to_numpy(bk_ptr, (npz+1,)),

        'u': fort_to_numpy(u_ptr, (ied-isd+1, jed+1-jsd+1, npz)),
        'v': fort_to_numpy(v_ptr, (ied+1-isd+1, jed-jsd+1, npz)),
        'w': fort_to_numpy(w_ptr, (ied-isd+1, jed-jsd+1, npz)),
        'delz': fort_to_numpy(delz_ptr, (ied-isd+1, jed-jsd+1, npz)),

        'pt': fort_to_numpy(pt_ptr, (ied-isd+1, jed-jsd+1, npz)),
        'delp': fort_to_numpy(delp_ptr, (ied-isd+1, jed-jsd+1, npz)),
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

        'ps': fort_to_numpy(ps_ptr, (ied-isd+1, jed-jsd+1)),
        'pe': fort_to_numpy(pe_ptr, (ie+1-(is_-1)+1, npz+1, je+1-(js-1)+1)),
        'pk': fort_to_numpy(pk_ptr, (ie-is_+1, je-js+1, npz+1)),
        'pkz': fort_to_numpy(pkz_ptr, (ie-is_+1, je-js+1, npz)),
        'peln': fort_to_numpy(peln_ptr, (ie-is_+1, npz+1, je-js+1)),

        'phis': fort_to_numpy(phis_ptr, (ied-isd+1, jed-jsd+1)),
        'q_con': fort_to_numpy(q_con_ptr, (ied-isd+1, jed-jsd+1, npz)),
        'omga': fort_to_numpy(omga_ptr, (ied-isd+1, jed-jsd+1, npz)),

        'ua': fort_to_numpy(ua_ptr, (ied-isd+1, jed-jsd+1, npz)),
        'va': fort_to_numpy(va_ptr, (ied-isd+1, jed-jsd+1, npz)),
        'uc': fort_to_numpy(uc_ptr, (ied+1-isd+1, jed-jsd+1, npz)),
        'vc': fort_to_numpy(va_ptr, (ied-isd+1, jed+1-jsd+1, npz)),

        'mfxd': fort_to_numpy(mfxd_ptr, (ie+1-is_+1, je-js+1, npz)),
        'mfyd': fort_to_numpy(mfyd_ptr, (ie-is_+1, je+1-js+1, npz)),
        'cxd': fort_to_numpy(cxd_ptr, (ie+1-is_+1, jed-jsd+1, npz)),
        'cyd': fort_to_numpy(cyd_ptr, (ied-isd+1, je+1-js+1, npz)),
        'diss_estd': fort_to_numpy(diss_estd_ptr, (ied-isd+1, jed-jsd+1, npz)),} # input_data

    return inout_data

def fortran_grid_data_to_numpy(
        npx, npy, npz,
        is_, ie, js, je, isd, ied, jsd, jed,
        nested_int, stretched_grid_int, da_min, da_min_c,
        # input arrays - grid data
        dx_ptr, dy_ptr, dxa_ptr, dya_ptr, dxc_ptr, dyc_ptr,
        rdx_ptr, rdy_ptr, rdxa_ptr, rdya_ptr, rdxc_ptr, rdyc_ptr,
        cosa_ptr, cosa_s_ptr, sina_u_ptr, sina_v_ptr,
        cosa_u_ptr, cosa_v_ptr, rsin2_ptr, rsina_ptr, rsin_u_ptr, rsin_v_ptr,
        sin_sg_ptr, cos_sg_ptr,
        area_ptr, rarea_ptr, rarea_c_ptr, f0_ptr, fC_ptr,
        del6_u_ptr, del6_v_ptr, divg_u_ptr, divg_v_ptr,
        agrid_ptr, bgrid_ptr, a11_ptr, a12_ptr, a21_ptr, a22_ptr,
        edge_e_ptr, edge_w_ptr, edge_n_ptr, edge_s_ptr):
    '''
    Convert Fortran arrays pointed to by *_ptr to NumPy arrays
    Input: Pointers to Fortran arrays *_ptr
    Output: dict where dict[key] is a NumPy array
    '''

    # sin/cos_sg require special handling
    sin_sg = fort_to_numpy(sin_sg_ptr, (ied-isd+1, jed-jsd+1, 9))
    cos_sg = fort_to_numpy(cos_sg_ptr, (ied-isd+1, jed-jsd+1, 9))

    # a/bgrid require special handling
    agrid = fort_to_numpy(agrid_ptr, (ied-isd+1, jed-jsd+1, 2))
    bgrid = fort_to_numpy(bgrid_ptr, (ied+1-isd+1, jed+1-jsd+1, 2))

    grid_data = {
        'npx': npx, 'npy': npy, 'npz': npz,
        'is_': is_, 'ie': ie, 'js': js, 'je': je,
        'isd': isd, 'ied': ied, 'jsd': jsd, 'jed': jed,

        'nested': False if nested_int == 0 else True,
        'stretched_grid': False if stretched_grid_int == 0 else True,
        'da_min': da_min, 'da_min_c': 0.0,

        'dx': fort_to_numpy(dx_ptr, (ied-isd+1, jed+1-jsd+1)),
        'dy': fort_to_numpy(dy_ptr, (ied+1-isd+1, jed-jsd+1)),
        'dxa': fort_to_numpy(dxa_ptr, (ied-isd+1, jed-jsd+1)),
        'dya': fort_to_numpy(dya_ptr, (ied-isd+1, jed-jsd+1)),
        'dxc': fort_to_numpy(dxc_ptr, (ied+1-isd+1, jed-jsd+1)),
        'dyc': fort_to_numpy(dyc_ptr, (ied-isd+1, jed+1-jsd+1)),

        'rdx': fort_to_numpy(rdx_ptr, (ied-isd+1, jed+1-jsd+1)),
        'rdy': fort_to_numpy(rdy_ptr, (ied+1-isd+1, jed-jsd+1)),
        'rdxa': fort_to_numpy(rdxa_ptr, (ied-isd+1, jed-jsd+1)),
        'rdya': fort_to_numpy(rdya_ptr, (ied-isd+1, jed-jsd+1)),
        'rdxc': fort_to_numpy(rdxc_ptr, (ied+1-isd+1, jed-jsd+1)),
        'rdyc': fort_to_numpy(rdyc_ptr, (ied-isd+1, jed+1-jsd+1)),

        'cosa': fort_to_numpy(cosa_ptr, (ied+1-isd+1, jed+1-jsd+1)),
        'cosa_s': fort_to_numpy(cosa_s_ptr, (ied-isd+1, jed-jsd+1)),

        'sina_u': fort_to_numpy(sina_u_ptr, (ied+1-isd+1, jed-jsd+1)),
        'sina_v': fort_to_numpy(sina_v_ptr, (ied-isd+1, jed+1-jsd+1)),
        'cosa_u': fort_to_numpy(cosa_u_ptr, (ied+1-isd+1, jed-jsd+1)),
        'cosa_v': fort_to_numpy(cosa_v_ptr, (ied-isd+1, jed+1-jsd+1)),

        'rsin2': fort_to_numpy(rsin2_ptr, (ied-isd+1, jed-jsd+1)),
        'rsina': fort_to_numpy(rsina_ptr, (ie+1-is_+1, je+1-js+1)),
        'rsin_u': fort_to_numpy(rsin_u_ptr, (ied+1-isd+1, jed-jsd+1)),
        'rsin_v': fort_to_numpy(rsin_v_ptr, (ied-isd+1, jed+1-jsd+1)),

        'sin_sg1': sin_sg[:,:,0], 'sin_sg2': sin_sg[:,:,1],
        'sin_sg3': sin_sg[:,:,2], 'sin_sg4': sin_sg[:,:,3],
        'cos_sg1': cos_sg[:,:,0], 'cos_sg2': cos_sg[:,:,1],
        'cos_sg3': cos_sg[:,:,2], 'cos_sg4': cos_sg[:,:,3],

        'area': fort_to_numpy(area_ptr, (ied-isd+1, jed-jsd+1)),
        'area_64': fort_to_numpy(area_ptr, (ied-isd+1, jed-jsd+1)), # TODO: area_64 = area??
        'rarea': fort_to_numpy(rarea_ptr, (ied-isd+1, jed-jsd+1)),
        'rarea_c':  fort_to_numpy(rarea_c_ptr, (ied+1-isd+1,jed+1-jsd+1)),
        'f0': fort_to_numpy(f0_ptr, (ied-isd+1, jed-jsd+1)),
        'fC': fort_to_numpy(fC_ptr, (ied+1-isd+1, jed+1-jsd+1)),

        'del6_u': fort_to_numpy(del6_u_ptr, (ied-isd+1, jed+1-jsd+1)),
        'del6_v': fort_to_numpy(del6_v_ptr, (ied+1-isd+1, jed-jsd+1)),
        'divg_u': fort_to_numpy(divg_u_ptr, (ied-isd+1, jed+1-jsd+1)),
        'divg_v': fort_to_numpy(divg_v_ptr, (ied+1-isd+1, jed-jsd+1)),

        'agrid1': agrid[:,:,0], 'agrid2': agrid[:,:,1],
        'bgrid1': bgrid[:,:,0], 'bgrid2': bgrid[:,:,1],

        'a11': fort_to_numpy(a11_ptr, (ie+1-(is_-1)+1, je+1-(js-1)+1)),
        'a12': fort_to_numpy(a12_ptr, (ie+1-(is_-1)+1, je+1-(js-1)+1)),
        'a21': fort_to_numpy(a21_ptr, (ie+1-(is_-1)+1, je+1-(js-1)+1)),
        'a22': fort_to_numpy(a22_ptr, (ie+1-(is_-1)+1, je+1-(js-1)+1)),

        'edge_e': fort_to_numpy(edge_e_ptr, (npy,)),
        'edge_w': fort_to_numpy(edge_w_ptr, (npy,)),
        'edge_n': fort_to_numpy(edge_n_ptr, (npx,)),
        'edge_s': fort_to_numpy(edge_s_ptr, (npx,)),} # grid_data

    return grid_data

def numpy_output_data_to_fortran(
        output_data,
        u_ptr, v_ptr, w_ptr, delz_ptr,
        pt_ptr, delp_ptr, q_ptr,
        ps_ptr, pe_ptr, pk_ptr, peln_ptr, pkz_ptr,
        phis_ptr, q_con_ptr, omga_ptr,
        ua_ptr, va_ptr, uc_ptr, vc_ptr,
        mfxd_ptr, mfyd_ptr, cxd_ptr, cyd_ptr, diss_estd_ptr):

    # NumPy -> Fortran
    u_out_f = output_data['u'].astype(np.float32).flatten(order='F')
    ffi.memmove(u_ptr, u_out_f, 4*u_out_f.size)
    v_out_f = output_data['v'].astype(np.float32).flatten(order='F')
    ffi.memmove(v_ptr, v_out_f, 4*v_out_f.size)
    w_out_f = output_data['w'].astype(np.float32).flatten(order='F')
    ffi.memmove(w_ptr, w_out_f, 4*w_out_f.size)
    delz_out_f = output_data['delz'].astype(np.float32).flatten(order='F')
    ffi.memmove(delz_ptr, delz_out_f, 4*delz_out_f.size)

    pt_out_f = output_data['pt'].astype(np.float32).flatten(order='F')
    ffi.memmove(pt_ptr, pt_out_f, 4*pt_out_f.size)
    delp_out_f = output_data['delp'].astype(np.float32).flatten(order='F')
    ffi.memmove(delp_ptr, delp_out_f, 4*delp_out_f.size)
    # q needs special handling
    q = np.empty(list(output_data['qvapor'].shape)+[9])
    q[:,:,:,0] = output_data['qvapor']
    q[:,:,:,1] = output_data['qliquid']
    q[:,:,:,2] = output_data['qice']
    q[:,:,:,3] = output_data['qrain']
    q[:,:,:,4] = output_data['qsnow']
    q[:,:,:,5] = output_data['qgraupel']
    q[:,:,:,6] = output_data['qcld']
    # q[:,:,:,7] = output_data['qo3mr']
    # q[:,:,:,8] = output_data['qsgs_tke']
    q_out_f = q.astype(np.float32).flatten(order='F')
    ffi.memmove(q_ptr, q_out_f, 4*q_out_f.size)

    ps_out_f = output_data['ps'].astype(np.float32).flatten(order='F')
    ffi.memmove(ps_ptr, ps_out_f, 4*ps_out_f.size)
    pe_out_f = output_data['pe'].astype(np.float32).flatten(order='F')
    ffi.memmove(pe_ptr, pe_out_f, 4*pe_out_f.size)
    pk_out_f = output_data['pk'].astype(np.float32).flatten(order='F')
    ffi.memmove(pk_ptr, pk_out_f, 4*pk_out_f.size)
    peln_out_f = output_data['peln'].astype(np.float32).flatten(order='F')
    ffi.memmove(peln_ptr, peln_out_f, 4*peln_out_f.size)
    pkz_out_f = output_data['pkz'].astype(np.float32).flatten(order='F')
    ffi.memmove(pkz_ptr, pkz_out_f, 4*pkz_out_f.size)

    phis_out_f = output_data['phis'].astype(np.float32).flatten(order='F')
    ffi.memmove(phis_ptr, phis_out_f, 4*phis_out_f.size)
    q_con_out_f = output_data['q_con'].astype(np.float32).flatten(order='F')
    ffi.memmove(q_con_ptr, q_con_out_f, 4*q_con_out_f.size)
    omga_out_f = output_data['omga'].astype(np.float32).flatten(order='F')
    ffi.memmove(omga_ptr, omga_out_f, 4*omga_out_f.size)

    ua_out_f = output_data['ua'].astype(np.float32).flatten(order='F')
    ffi.memmove(ua_ptr, ua_out_f, 4*ua_out_f.size)
    va_out_f = output_data['va'].astype(np.float32).flatten(order='F')
    ffi.memmove(va_ptr, va_out_f, 4*va_out_f.size)
    uc_out_f = output_data['uc'].astype(np.float32).flatten(order='F')
    ffi.memmove(uc_ptr, uc_out_f, 4*uc_out_f.size)
    vc_out_f = output_data['vc'].astype(np.float32).flatten(order='F')
    ffi.memmove(vc_ptr, vc_out_f, 4*vc_out_f.size)

    mfxd_out_f = output_data['mfxd'].astype(np.float32).flatten(order='F')
    ffi.memmove(mfxd_ptr, mfxd_out_f, 4*mfxd_out_f.size)
    mfyd_out_f = output_data['mfyd'].astype(np.float32).flatten(order='F')
    ffi.memmove(mfyd_ptr, mfyd_out_f, 4*mfyd_out_f.size)
    cxd_out_f = output_data['cxd'].astype(np.float32).flatten(order='F')
    ffi.memmove(cxd_ptr, cxd_out_f, 4*cxd_out_f.size)
    cyd_out_f = output_data['cyd'].astype(np.float32).flatten(order='F')
    ffi.memmove(cyd_ptr, cyd_out_f, 4*cyd_out_f.size)
    diss_estd_out_f = output_data['diss_estd'].astype(np.float32).flatten(order='F')
    ffi.memmove(diss_estd_ptr, diss_estd_out_f, 4*diss_estd_out_f.size)
