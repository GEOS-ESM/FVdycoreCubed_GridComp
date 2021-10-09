def fv_dynamics_top_level_function(
        comm_py, npx, npy, npz, nq_tot, ng,
        bdt, consv_te, fill, reproduce_sum, kappa,
        cp_air, zvir, ptop, ks, ncnst, n_split, q_split,
        u, v, w):
    print('fv_dynamics_top_level_function')
