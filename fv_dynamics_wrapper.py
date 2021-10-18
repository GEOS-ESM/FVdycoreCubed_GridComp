import cffi
from mpi4py import MPI

TMPFILEBASE = 'fvdynwrap'

ffi = cffi.FFI()

# MPI_Comm can be int or void*
if MPI._sizeof(MPI.Comm) == ffi.sizeof('int'):
    _mpi_comm_t = 'int'
else:
    _mpi_comm_t = 'void*'

source = '''
from {} import ffi
from mpi4py import MPI
from fv_dynamics import fv_dynamics_top_level_function
@ffi.def_extern()
def fv_dynamics_py_wrapper(
        comm_c,
        npx, npy, npz, nq_tot, ng,
        is1, ie, js, je,
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
        # ze0,
        hybrid_z):
    # comm_c -> comm_py
    comm_py = MPI.Intracomm() # new comm, internal MPI_Comm handle is MPI_COMM_NULL
    comm_ptr = MPI._addressof(comm_py)  # internal MPI_Comm handle
    comm_ptr = ffi.cast('{}*', comm_ptr)  # make it a CFFI pointer
    comm_ptr[0] = comm_c  # assign comm_c to comm_py's MPI_Comm handle
    print('fv_dynamics_py_wrapper')
    # Call top-level function for fv_dynamics
    fv_dynamics_top_level_function(
        comm_py, 
        npx, npy, npz, nq_tot, ng,
        is1, ie, js, je,
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
        # ze0,
        hybrid_z)
'''.format(TMPFILEBASE, _mpi_comm_t)

header = '''
extern void fv_dynamics_py_wrapper(
    {} comm_c,
    int npx, int npy, int npz, int nq_tot, int ng,
    int is1, int ie, int js, int je,
    int isd, int ied, int jsd, int jed,
    float bdt, float consv_te, int fill, int reproduce_sum,
    float kappa, float cp_air, float zvir, float ptop,
    int ks, int ncnst, int n_split, int q_split,
    float* u, float* v, float* w, float* delz,
    int hydrostatic,
    float* pt, float* delp, float* q,
    float* ps, float* pe, float* pk, float* peln, float* pkz,
    float* phis, float* q_con, float* omga,
    float* ua, float* va, float* uc, float* vc,
    const float* ak, const float* bk,
    float* mfx, float* mfy, float* cx, float* cy,
    // float* ze0,
    int hybrid_z);'''.format(_mpi_comm_t)
with open(TMPFILEBASE+'.h', 'w') as f:
    f.write(header)
ffi.embedding_api(header)

source_header = r'''#include "{}.h"'''.format(TMPFILEBASE)
ffi.set_source(TMPFILEBASE, source_header)

ffi.embedding_init_code(source)
ffi.compile(target='lib'+TMPFILEBASE+'.so', verbose=True)
