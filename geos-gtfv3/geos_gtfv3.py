import f90nml
from datetime import datetime
from f_py_conversion import FortranPythonConversion
from pace.fv3core.initialization.geos_wrapper import GeosDycoreWrapper
from cuda_profiler import CUDAProfiler
from mpi4py import MPI


class GEOSGTFV3:
    def __init__(
        self,
        namelist_path: str,
        bdt: float,
        comm: MPI.Intercomm,
        npx,
        npy,
        npz,
        is_,
        ie,
        js,
        je,
        isd,
        ied,
        jsd,
        jed,
        nq_tot,
        backend: str = "dace:gpu",
    ) -> None:
        self.namelist = f90nml.read(namelist_path)
        self.rank = comm.Get_rank()
        self.backend = backend
        self.dycore = GeosDycoreWrapper(self.namelist, bdt, comm, self.backend)
        # For Fortran<->NumPy conversion
        self.f_py = FortranPythonConversion(
            npx,
            npy,
            npz,
            is_,
            ie,
            js,
            je,
            isd,
            ied,
            jsd,
            jed,
            nq_tot,
        )

    def __call__(
        self,
        ng,
        ptop,
        ks,
        layout_1,
        layout_2,
        adiabatic,
        u,
        v,
        w,
        delz,
        pt,
        delp,
        q,
        ps,
        pe,
        pk,
        peln,
        pkz,
        phis,
        q_con,
        omga,
        ua,
        va,
        uc,
        vc,
        ak,
        bk,
        mfx,
        mfy,
        cx,
        cy,
        diss_est,
    ):
        RANK_PRINT = 0

        if self.rank == RANK_PRINT:
            print(
                "P:",
                datetime.now().isoformat(timespec="milliseconds"),
                "--in top level geos-gtfv3, backend:",
                self.backend,
                flush=True,
            )

        CUDAProfiler.start_cuda_profiler()
        with CUDAProfiler("fortran->py"):

            # Convert Fortran arrays to NumPy
            state_in = self.f_py.fortran_to_numpy(
                # input
                u,
                v,
                w,
                delz,
                pt,
                delp,
                q,
                ps,
                pe,
                pk,
                peln,
                pkz,
                phis,
                q_con,
                omga,
                ua,
                va,
                uc,
                vc,
                mfx,
                mfy,
                cx,
                cy,
                diss_est,
            )
            # write_sum_of_vars(comm, state_in)
            # write_shape_or_value(comm, state_in)

        if self.rank == RANK_PRINT:
            print(
                "P:",
                datetime.now().isoformat(timespec="milliseconds"),
                "--fortran->numpy, transpose, sp->dp, swap axes",
                flush=True,
            )

        if self.rank == RANK_PRINT:
            print(
                "P:",
                datetime.now().isoformat(timespec="milliseconds"),
                "--initialized dycore",
                flush=True,
            )

        # Run gtFV3
        with CUDAProfiler("Py dycore runtime"):
            state_out = self.dycore(
                state_in["u"],
                state_in["v"],
                state_in["w"],
                state_in["delz"],
                state_in["pt"],
                state_in["delp"],
                state_in["q"],
                state_in["ps"],
                state_in["pe"],
                state_in["pk"],
                state_in["peln"],
                state_in["pkz"],
                state_in["phis"],
                state_in["q_con"],
                state_in["omga"],
                state_in["ua"],
                state_in["va"],
                state_in["uc"],
                state_in["vc"],
                state_in["mfxd"],
                state_in["mfyd"],
                state_in["cxd"],
                state_in["cyd"],
                state_in["diss_estd"],
            )

        if self.rank == RANK_PRINT:
            print(
                "P:",
                datetime.now().isoformat(timespec="milliseconds"),
                "--ran dycore",
                flush=True,
            )

        # Convert NumPy arrays back to Fortran
        with CUDAProfiler("Py -> Fortran"):
            self.f_py.numpy_to_fortran(
                # input
                state_out,
                # output
                u,
                v,
                w,
                delz,
                pt,
                delp,
                q,
                ps,
                pe,
                pk,
                peln,
                pkz,
                phis,
                q_con,
                omga,
                ua,
                va,
                uc,
                vc,
                mfx,
                mfy,
                cx,
                cy,
                diss_est,
            )

        if self.rank == RANK_PRINT:
            print(
                "P:",
                datetime.now().isoformat(timespec="milliseconds"),
                "--numpy->fortran, transpose, dp->sp, swap axes",
                flush=True,
            )


# Below is the entry point to the interface
# ToDo: we should build the object outside of the sim loop from fortran
# potentially by writing a geos_gtfv3_setup_interface and caching the ptr Fortran side
GEOS_DYCORE = None


def geos_gtfv3(
    comm,
    npx,
    npy,
    npz,
    ntiles,
    is_,
    ie,
    js,
    je,
    isd,
    ied,
    jsd,
    jed,
    bdt,
    nq_tot,
    ng,
    ptop,
    ks,
    layout_1,
    layout_2,
    adiabatic,
    u,
    v,
    w,
    delz,
    pt,
    delp,
    q,
    ps,
    pe,
    pk,
    peln,
    pkz,
    phis,
    q_con,
    omga,
    ua,
    va,
    uc,
    vc,
    ak,
    bk,
    mfx,
    mfy,
    cx,
    cy,
    diss_est,
):
    BACKEND = "dace:gpu"
    NAMELIST_PATH = "input.nml"

    global GEOS_DYCORE
    if GEOS_DYCORE is None:
        GEOS_DYCORE = GEOSGTFV3(
            namelist_path=NAMELIST_PATH,
            bdt=bdt,
            comm=comm,
            npx=npx,
            npy=npy,
            npz=npz,
            is_=is_,
            ie=ie,
            js=js,
            je=je,
            isd=isd,
            ied=ied,
            jsd=jsd,
            jed=jed,
            nq_tot=nq_tot,
            backend=BACKEND,
        )
    GEOS_DYCORE(
        ng,
        ptop,
        ks,
        layout_1,
        layout_2,
        adiabatic,
        u,
        v,
        w,
        delz,
        pt,
        delp,
        q,
        ps,
        pe,
        pk,
        peln,
        pkz,
        phis,
        q_con,
        omga,
        ua,
        va,
        uc,
        vc,
        ak,
        bk,
        mfx,
        mfy,
        cx,
        cy,
        diss_est,
    )
