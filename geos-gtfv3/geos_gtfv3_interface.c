#include <stdio.h>
#include <time.h>
#include "mpi.h"
#include "geos_gtfv3_interface_py.h"

void geos_gtfv3_interface_c_init(
    MPI_Fint comm_f,
    int npx, int npy, int npz, int ntiles,
    int is, int ie, int js, int je,
    int isd, int ied, int jsd, int jed,
    float bdt, int nq_tot)
{
    MPI_Comm comm_c = MPI_Comm_f2c(comm_f);
    geos_gtfv3_interface_py_init(
        comm_c,
        npx, npy, npz, ntiles,
        is, ie, js, je, isd, ied, jsd, jed, bdt, nq_tot);
}

void geos_gtfv3_interface_c(
    // input
    MPI_Fint comm_f,
    int npx, int npy, int npz, int ntiles,
    int is, int ie, int js, int je,
    int isd, int ied, int jsd, int jed,
    float bdt, int nq_tot, int ng, float ptop, int ks, int layout_1, int layout_2,
    int adiabatic,

    // input/output
    float *u, float *v, float *w, float *delz,
    float *pt, float *delp, float *q,
    float *ps, float *pe, float *pk, float *peln, float *pkz,
    float *phis, float *q_con, float *omga, float *ua, float *va, float *uc, float *vc,

    // input
    const float *ak, const float *bk,

    // input/output
    float *mfx, float *mfy, float *cx, float *cy, float *diss_est)
{

    MPI_Comm comm_c = MPI_Comm_f2c(comm_f);
    /* printf("MPI communicator (F): %d\n", comm_f); */
    /* printf("MPI communicator (C): %d\n", comm_c); */

    // time_t now = time(&now);
    // char buf[20];
    // strftime(buf, sizeof(buf), "%FT%T", localtime(&now));
    // int rank;
    // MPI_Comm_rank(comm_c, &rank);
    // if (rank == 0) printf("C: %s.xxx --calling python interface\n", buf);

    geos_gtfv3_interface_py(
        comm_c,
        npx, npy, npz, ntiles,
        is, ie, js, je,
        isd, ied, jsd, jed,
        bdt, nq_tot, ng, ptop, ks, layout_1, layout_2,
        adiabatic,

        u, v, w, delz,
        pt, delp, q,
        ps, pe, pk, peln, pkz,
        phis, q_con, omga,
        ua, va, uc, vc,

        ak, bk,

        mfx, mfy, cx, cy, diss_est);
}

void geos_gtfv3_interface_c_finalize()
{
    geos_gtfv3_interface_py_finalize();
}
