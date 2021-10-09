#include <stdio.h>
#include "mpi.h"
#include "fvdynwrap.h"

void fv_dynamics_interface(
    // input
    MPI_Fint comm_f, int npx, int npy, int npz, int nq_tot, int ng,
    float bdt, float consv_te, int fill, int reproduce_sum, float kappa,
    float cp_air, float zvir, float ptop, int ks, int ncnst, int n_split, int q_split,
    // input/output
    float* u, float* v, float* w, float* delz,
    // input
    int hydrostatic,
    // input/output
    float* pt, float* delp, float* q,
    float* ps, float* pe, float* pk, float* peln, float* pkz, // can be re-computed from delp and ptop
    float* phis, float* q_con, float* omga,
    float* ua, float* va, float* uc, float* vc,
    // input
    const float* ak, const float* bk,
    // input/output
    float* mfx, float* mfy, float* cx, float* cy,
    // TODO: Is this used anywhere?
    float* ze0,
    // input
    int hybrid_z) {

    printf("MPI communicator: %d\n", comm_f);
    MPI_Comm comm_c = MPI_Comm_f2c(comm_f);
    fv_dynamics_py_wrapper(
        comm_c, npx, npy, npz, nq_tot, ng,
        bdt, consv_te, fill, reproduce_sum, kappa,
        cp_air, zvir, ptop, ks, ncnst, n_split, q_split,
        u, v, w);
}
