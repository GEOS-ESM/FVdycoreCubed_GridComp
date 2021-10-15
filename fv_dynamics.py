import numpy as np
from mpi4py import MPI
from fort2py import fort_to_gt4py

def fv_dynamics_top_level_function(
        comm,
        npx, npy, npz, nq_tot, ng,
        isd, ied, jsd, jed, bdt,
        consv_te, fill, reproduce_sum,
        kappa, cp_air, zvir, ptop,
        ks, ncnst, n_split, q_split,
        u_p, v_p, w_p # _p => pointer to Fortran data
):
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
    u = fort_to_gt4py(u_p, (ied-isd+1, jed+1-jsd+1, npz), (0, 0, 0), backend)
    print(u.shape)
    with np.printoptions(formatter={'float': lambda x: "{0:05.2f}".format(x)}, linewidth=10000):
        print(u[0,:,:], flush=True)
