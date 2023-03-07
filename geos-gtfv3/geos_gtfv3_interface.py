import cffi
from mpi4py import MPI

TMPFILEBASE = "geos_gtfv3_interface_py"

ffi = cffi.FFI()

# MPI_Comm can be int or void*
if MPI._sizeof(MPI.Comm) == ffi.sizeof("int"):
    _mpi_comm_t = "int"
else:
    _mpi_comm_t = "void*"

source = """
from {} import ffi
from datetime import datetime
from mpi4py import MPI
from geos_gtfv3 import geos_gtfv3, geos_gtfv3_finalize

@ffi.def_extern()
def geos_gtfv3_interface_py(
    comm_c,
    npx, npy, npz, ntiles,
    is_, ie, js, je, isd, ied, jsd, jed,
    bdt, nq_tot, ng, ptop, ks, layout_1, layout_2, adiabatic,
    u, v, w, delz,
    pt, delp, q,
    ps, pe, pk, peln, pkz,
    phis, q_con, omga,
    ua, va, uc, vc,
    ak, bk,
    mfx, mfy, cx, cy, diss_est):

    # comm_c -> comm_py
    comm_py = MPI.Intracomm() # new comm, internal MPI_Comm handle is MPI_COMM_NULL
    comm_ptr = MPI._addressof(comm_py)  # internal MPI_Comm handle
    comm_ptr = ffi.cast('{}*', comm_ptr)  # make it a CFFI pointer
    comm_ptr[0] = comm_c  # assign comm_c to comm_py's MPI_Comm handle

    # if comm_py.Get_rank() == 0:
    #     print('P:', datetime.now().isoformat(timespec='milliseconds'),
    #           '--in cffi interface', flush=True)

    geos_gtfv3(
        comm_py,
        npx, npy, npz, ntiles,
        is_, ie, js, je, isd, ied, jsd, jed,
        bdt, nq_tot, ng, ptop, ks, layout_1, layout_2, adiabatic,
        u, v, w, delz,
        pt, delp, q,
        ps, pe, pk, peln, pkz,
        phis, q_con, omga,
        ua, va, uc, vc,
        ak, bk,
        mfx, mfy, cx, cy, diss_est)

@ffi.def_extern()
def geos_gtfv3_interface_finalize_py():
    geos_gtfv3_finalize()

""".format(
    TMPFILEBASE, _mpi_comm_t
)

header = """
extern void geos_gtfv3_interface_py(
    {} comm_c,
    int npx, int npy, int npz, int ntiles,
    int is_, int ie, int js, int je, int isd, int ied, int jsd, int jed,
    float bdt, int nq_tot, int ng, float ptop, int ks, int layout_1, int layout_2,
    int adiabatic,
    // input/output
    float* u, float* v, float* w, float* delz,
    float* pt, float* delp, float* q,
    float* ps, float* pe, float* pk, float* peln, float* pkz,
    float* phis, float* q_con, float* omga, float* ua, float* va, float* uc, float* vc,
    // input
    const float* ak, const float* bk,
    // input/output
    float* mfx, float* mfy, float* cx, float* cy, float* diss_est);

extern void geos_gtfv3_interface_finalize_py();
""".format(
    _mpi_comm_t
)

with open(TMPFILEBASE + ".h", "w") as f:
    f.write(header)
ffi.embedding_api(header)

source_header = r'''#include "{}.h"'''.format(TMPFILEBASE)
ffi.set_source(TMPFILEBASE, source_header)

ffi.embedding_init_code(source)
ffi.compile(target="lib" + TMPFILEBASE + ".so", verbose=True)
