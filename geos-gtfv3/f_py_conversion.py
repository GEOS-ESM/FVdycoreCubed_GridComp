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



# Create a class FortranPythonConversion with state consisting of
# is/ie/js/je, isd/ied/jsd/jed, nq_tot etc



def fortran_state_to_numpy(
        npx, npy, npz,
        is_, ie, js, je, isd, ied, jsd, jed,
        bdt, nq_tot, ptop, ks,
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
    q = fort_to_numpy(q_ptr, (ied-isd+1, jed-jsd+1, npz, nq_tot))

    numpy_state = {
        'do_adiabatic_init': True, # REQUIRED (TODO???)
        'bdt': bdt,
        'ptop': ptop,
        'ks': ks,

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

    return numpy_state



def numpy_output_data_to_fortran(
        output_data, nq_tot,
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
    q = np.empty(list(output_data['qvapor'].shape)+[nq_tot])
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
