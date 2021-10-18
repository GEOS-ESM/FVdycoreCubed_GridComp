import numpy as np
from mpi4py import MPI
from fort2py import fort_to_gt4py

def fv_dynamics_top_level_function(
        comm,
        npx, npy, npz, nq_tot, ng,
        is1, ie, js, je, # var is is renamed as is1
        isd, ied, jsd, jed, bdt,
        consv_te, fill_int, reproduce_sum_int,
        kappa, cp_air, zvir, ptop,
        ks, ncnst, n_split, q_split,
        u_ptr, v_ptr, w_ptr, delz_ptr, # _ptr => pointer to Fortran data
        hydrostatic_int,
        pt_ptr, delp_ptr, q_ptr,
        ps_ptr, pe_ptr, pk_ptr, peln_ptr, pkz_ptr,
        phis_ptr, q_con_ptr, omga_ptr,
        ua_ptr, va_ptr, uc_ptr, vc_ptr,
        ak_ptr, bk_ptr,
        mfx_ptr, mfy_ptr, cx_ptr, cy_ptr,
        # ze0_ptr,
        hybrid_z_int
):

    # Manage logicals
    fill = True if fill_int == 1 else False
    reproduce_sum = True if reproduce_sum_int == 1 else False
    hydrostatic = True if hydrostatic_int == 1 else False
    hybrid_z = True if hybrid_z_int == 1 else False

    print(
        'top level: I am process {} of {} on {}'.format(
            comm.Get_rank(),
            comm.Get_size(),
            MPI.Get_processor_name()),
        flush=True)
    comm.Barrier()

    print('domain dimensions:', npx, npy, npz)
    print('grid bounds:', isd, ied, jsd, jed)
    print('nq_tot, ng:', nq_tot, ng)
    print('consv_te, fill, reproduce_sum:', consv_te, fill, reproduce_sum)
    print('kappa, cp_air, zvir, ptop:', kappa, cp_air, zvir, ptop)
    print('ks, ncnst, n_split, q_split:', ks, ncnst, n_split, q_split)

    backend = 'numpy'

    # Convert Fortran arrays to GT4Py storage
    u = fort_to_gt4py(u_ptr, (ied-isd+1, jed+1-jsd+1, npz), (0, 0, 0), backend)
    v = fort_to_gt4py(v_ptr, (ied+1-isd+1, jed-jsd+1, npz), (0, 0, 0), backend)
    w = fort_to_gt4py(w_ptr, (1, 1, 1), (0, 0, 0), backend)
    delz = fort_to_gt4py(delz_ptr, (1, 1, 1), (0, 0, 0), backend)
    # with np.printoptions(formatter={'float': lambda x: "{0:06.2f}".format(x)}, linewidth=10000):
    print('u, v, w, delz:', np.sum(u), np.sum(v), np.sum(w), np.sum(delz))
    print('hydrostatic:', hydrostatic)
    pt = fort_to_gt4py(pt_ptr, (ied-isd+1, jed-jsd+1, npz), (0, 0, 0), backend)
    delp = fort_to_gt4py(delp_ptr, (ied-isd+1, jed-jsd+1, npz), (0, 0, 0), backend)
    pt = fort_to_gt4py(pt_ptr, (ied-isd+1, jed-jsd+1, npz), (0, 0, 0), backend)
    q = fort_to_gt4py(q_ptr, (ied-isd+1, jed-jsd+1, npz, ncnst), (0, 0, 0, 0), backend)
    print('pt, delp, q:', np.sum(pt), np.sum(delp), np.sum(q))
    ps = fort_to_gt4py(ps_ptr, (ied-isd+1, jed-jsd+1), (0, 0), backend)
    pe = fort_to_gt4py(pe_ptr, (ie+1-(is1-1)+1, npz+1, je+1-(js-1)+1), (0, 0, 0), backend)
    pk = fort_to_gt4py(pk_ptr, (ie-is1+1, je-js+1, npz+1), (0, 0, 0), backend)
    peln = fort_to_gt4py(peln_ptr, (ie-is1+1, npz+1, je-js+1), (0, 0, 0), backend)
    pkz = fort_to_gt4py(pkz_ptr, (ie-is1+1, je-js+1, npz), (0, 0, 0), backend)
    print('ps, pe, pk, peln, pkz:', np.sum(ps), np.sum(pe), np.sum(pk), np.sum(peln), np.sum(pkz))
    # print(ps.shape, pe.shape, pk.shape, peln.shape, pkz.shape)
    phis = fort_to_gt4py(phis_ptr, (ied-isd+1, jed-jsd+1), (0, 0), backend)
    q_con = fort_to_gt4py(q_con_ptr, (1, 1, 1), (0, 0, 0), backend)
    omga = fort_to_gt4py(omga_ptr, (ied-isd+1, jed-jsd+1, npz), (0, 0, 0), backend)
    print('phis, q_con, omga:', np.sum(phis), np.sum(q_con), np.sum(omga))
    ua = fort_to_gt4py(ua_ptr, (ied-isd+1, jed-jsd+1, npz), (0, 0, 0), backend)
    va = fort_to_gt4py(va_ptr, (ied-isd+1, jed-jsd+1, npz), (0, 0, 0), backend)
    uc = fort_to_gt4py(uc_ptr, (ied+1-isd+1, jed-jsd+1, npz), (0, 0, 0), backend)
    vc = fort_to_gt4py(va_ptr, (ied-isd+1, jed+1-jsd+1, npz), (0, 0, 0), backend)
    print('ua, va, uc, vc:', np.sum(ua), np.sum(va), np.sum(uc), np.sum(vc))
    ak = fort_to_gt4py(ak_ptr, (npz+1,), (0,), backend)
    bk = fort_to_gt4py(bk_ptr, (npz+1,), (0,), backend)
    mfx = fort_to_gt4py(mfx_ptr, (ie+1-is1+1, je-js+1, npz), (0, 0, 0), backend)
    mfy = fort_to_gt4py(mfy_ptr, (ie-is1+1, je+1-js+1, npz), (0, 0, 0), backend)
    cx = fort_to_gt4py(cx_ptr, (ie+1-is1+1, jed-jsd+1, npz), (0, 0, 0), backend)
    cy = fort_to_gt4py(cy_ptr, (ied-isd+1, je+1-js+1, npz), (0, 0, 0), backend)
    # ze0
    print('ak, bk:', np.sum(ak), np.sum(bk))
    print('mfx, mfy, cx, cy:', np.sum(mfx), np.sum(mfy), np.sum(cx), np.sum(cy))
    print('hybrid_z:', hybrid_z)
