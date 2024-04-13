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
from geos_gtfv3 import geos_gtfv3_init, geos_gtfv3, geos_gtfv3_finalize
import traceback

@ffi.def_extern()
def geos_gtfv3_interface_py_init(
    fv_flags,
    comm_c,
    npx, npy, npz, ntiles,
    is_, ie, js, je, isd, ied, jsd, jed,
    bdt, nq_tot, ak, bk,
    ) -> int:

    # comm_c -> comm_py
    comm_py = MPI.Intracomm() # new comm, internal MPI_Comm handle is MPI_COMM_NULL
    comm_ptr = MPI._addressof(comm_py)  # internal MPI_Comm handle
    comm_ptr = ffi.cast('{}*', comm_ptr)  # make it a CFFI pointer
    comm_ptr[0] = comm_c  # assign comm_c to comm_py's MPI_Comm handle
    
    try:
        geos_gtfv3_init(
            fv_flags,
            comm_py,
            npx, npy, npz, ntiles,
            is_, ie, js, je, isd, ied, jsd, jed,
            bdt, nq_tot, ak, bk
            )
    except Exception as err:
        print("Error in Python:")
        print(traceback.format_exc())
        return -1
    return 0

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
def geos_gtfv3_interface_py_finalize():
    geos_gtfv3_finalize()

""".format(
    TMPFILEBASE, _mpi_comm_t, _mpi_comm_t
)

with open("fv_flags.h") as f:
    data = "".join([line for line in f if not line.startswith("#")])
    data = data.replace("CFFI_DLLEXPORT", "")
    ffi.embedding_api(data)

ffi.set_source(TMPFILEBASE, '#include "fv_flags.h"')

ffi.embedding_init_code(source)
ffi.compile(target="lib" + TMPFILEBASE + ".so", verbose=True)
