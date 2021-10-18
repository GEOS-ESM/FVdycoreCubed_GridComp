#include <stdio.h>
#include "mpi.h"
#include "fvdynwrap.h"

void fv_dynamics_interface(
    // input
    MPI_Fint comm_f,
    int npx, int npy, int npz, int nq_tot, int ng,
    int is, int ie, int js, int je,
    int isd, int ied, int jsd, int jed,
    float bdt, float consv_te, int fill, int reproduce_sum,
    float kappa, float cp_air, float zvir, float ptop,
    int ks, int ncnst, int n_split, int q_split,
    // input/output
    float* u, float* v, float* w, float* delz,
    // input
    int hydrostatic,
    // input/output
    float* pt, float* delp, float* q,
    float* ps, float* pe, float* pk, float* peln, float* pkz,
    float* phis, float* q_con, float* omga,
    float* ua, float* va, float* uc, float* vc,
    /* // input */
    const float* ak, const float* bk,
    // input/output
    float* mfx, float* mfy, float* cx, float* cy,
    /* // TODO: Is this used anywhere? */
    /* float* ze0, */
    // input
    int hybrid_z
    ) {

    MPI_Comm comm_c = MPI_Comm_f2c(comm_f);
    printf("MPI communicator (F): %d\n", comm_f);
    printf("MPI communicator (C): %d\n", comm_c);
    fv_dynamics_py_wrapper(
        comm_c,
        npx, npy, npz, nq_tot, ng,
        is, ie, js, je,
        isd, ied, jsd, jed,
        bdt, consv_te, fill, reproduce_sum,
        kappa, cp_air, zvir, ptop,
        ks, ncnst, n_split, q_split,
        u, v, w, delz,
        hydrostatic,
        pt, delp, q,
        ps, pe, pk, peln, pkz,
        phis, q_con, omga,
        ua, va, uc, vc,
        ak, bk,
        mfx, mfy, cx, cy,
        // ze0,
        hybrid_z
        );
}
