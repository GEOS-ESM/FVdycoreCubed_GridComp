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
        comm_c, npx, npy, npz, nq_tot, ng,
        bdt, consv_te, fill, reproduce_sum, kappa,
        cp_air, zvir, ptop, ks, ncnst, n_split, q_split,
        u, v, w):
    # comm_c -> comm_py
    comm_py = MPI.Intracomm() # new comm, internal MPI_Comm handle is MPI_COMM_NULL
    comm_ptr = MPI._addressof(comm_py)  # internal MPI_Comm handle
    comm_ptr = ffi.cast('{}*', comm_ptr)  # make it a CFFI pointer
    comm_ptr[0] = c_comm  # assign c_comm to comm_py's MPI_Comm handle
    # Call top-level function for fv_dynamics
    fv_dynamics_top_level_function(
        comm_py, npx, npy, npz, nq_tot, ng,
        bdt, consv_te, fill, reproduce_sum, kappa,
        cp_air, zvir, ptop, ks, ncnst, n_split, q_split,
        u, v, w)
'''.format(TMPFILEBASE, _mpi_comm_t)

header = '''
extern void fv_dynamics_py_wrapper(
    {} comm_c, int npx, int npy, int npz, int nq_tot, int ng,
    float bdt, float consv_te, int fill, int reproduce_sum, float kappa,
    float cp_air, float zvir, float ptop, int ks, int ncnst, int n_split, int q_split,
    float* u, float* v, float* w);'''.format(_mpi_comm_t)
with open(TMPFILEBASE+'.h', 'w') as f:
    f.write(header)
ffi.embedding_api(header)

source_header = r'''#include "{}.h"'''.format(TMPFILEBASE)
ffi.set_source(TMPFILEBASE, source_header)

ffi.embedding_init_code(source)
ffi.compile(target='lib'+TMPFILEBASE+'.so', verbose=True)
