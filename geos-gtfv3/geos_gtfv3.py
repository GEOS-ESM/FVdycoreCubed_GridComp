import f90nml
from f_py_conversion import FortranPythonConversion
from pace.fv3core.initialization.geos_wrapper import GeosDycoreWrapper, MemorySpace
from cuda_profiler import CUDAProfiler, TimedCUDAProfiler
from mpi4py import MPI
from pace.util._optional_imports import cupy as cp
import numpy as np
from pace.dsl.gt4py_utils import is_gpu_backend
from typing import TYPE_CHECKING
import os

if TYPE_CHECKING:
    import cffi


class GEOSGTFV3:
    def __init__(
        self,
        namelist_path: str,
        bdt: float,
        comm: MPI.Intercomm,
        npx: int,
        npy: int,
        npz: int,
        is_: int,
        ie: int,
        js: int,
        je: int,
        isd: int,
        ied: int,
        jsd: int,
        jed: int,
        nq_tot: int,
        backend: str = "dace:gpu",
    ) -> None:
        self.namelist = f90nml.read(namelist_path)
        self.rank = comm.Get_rank()
        self.backend = backend
        # For Fortran<->NumPy conversion
        if is_gpu_backend(self.backend):
            numpy_module = cp
            fortran_mem_space = MemorySpace.DEVICE
        else:
            numpy_module = np
            fortran_mem_space = MemorySpace.HOST
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
            numpy_module,
        )
        # Setup Pace's dynamical core
        self.dycore = GeosDycoreWrapper(
            self.namelist, bdt, comm, self.backend, fortran_mem_space
        )

        self._timings = {}

    def finalize(self):
        import json

        with open("gtfv3_timings.json", "w") as f:
            json.dump(self._timings, f, indent=4)

    def __call__(
        self,
        ng,
        ptop,
        ks,
        layout_1,
        layout_2,
        adiabatic,
        u: "cffi.FFI.CData",
        v: "cffi.FFI.CData",
        w: "cffi.FFI.CData",
        delz: "cffi.FFI.CData",
        pt: "cffi.FFI.CData",
        delp: "cffi.FFI.CData",
        q: "cffi.FFI.CData",
        ps: "cffi.FFI.CData",
        pe: "cffi.FFI.CData",
        pk: "cffi.FFI.CData",
        peln: "cffi.FFI.CData",
        pkz: "cffi.FFI.CData",
        phis: "cffi.FFI.CData",
        q_con: "cffi.FFI.CData",
        omga: "cffi.FFI.CData",
        ua: "cffi.FFI.CData",
        va: "cffi.FFI.CData",
        uc: "cffi.FFI.CData",
        vc: "cffi.FFI.CData",
        ak: "cffi.FFI.CData",
        bk: "cffi.FFI.CData",
        mfx: "cffi.FFI.CData",
        mfy: "cffi.FFI.CData",
        cx: "cffi.FFI.CData",
        cy: "cffi.FFI.CData",
        diss_est: "cffi.FFI.CData",
    ):
        CUDAProfiler.start_cuda_profiler()
        with TimedCUDAProfiler("Fortran -> Python", self._timings):
            # Convert Fortran arrays to NumPy
            state_in = self.f_py.fortran_to_python(
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

        # Run gtFV3
        with TimedCUDAProfiler("Numerics", self._timings):
            state_out, self._timings = self.dycore(
                self._timings,
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

        # Convert NumPy arrays back to Fortran
        with TimedCUDAProfiler("Python -> Fortran", self._timings):
            self.f_py.python_to_fortran(
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


# Below is the entry point to the interface
# ToDo: we should build the object outside of the sim loop from fortran
# potentially by writing a geos_gtfv3_setup_interface and caching the ptr Fortran side
GEOS_DYCORE = None


def geos_gtfv3(
    comm: MPI.Intercomm,
    npx: int,
    npy: int,
    npz: int,
    ntiles: int,
    is_: int,
    ie: int,
    js: int,
    je: int,
    isd: int,
    ied: int,
    jsd: int,
    jed: int,
    bdt,
    nq_tot,
    ng,
    ptop,
    ks,
    layout_1,
    layout_2,
    adiabatic,
    u: "cffi.FFI.CData",
    v: "cffi.FFI.CData",
    w: "cffi.FFI.CData",
    delz: "cffi.FFI.CData",
    pt: "cffi.FFI.CData",
    delp: "cffi.FFI.CData",
    q: "cffi.FFI.CData",
    ps: "cffi.FFI.CData",
    pe: "cffi.FFI.CData",
    pk: "cffi.FFI.CData",
    peln: "cffi.FFI.CData",
    pkz: "cffi.FFI.CData",
    phis: "cffi.FFI.CData",
    q_con: "cffi.FFI.CData",
    omga: "cffi.FFI.CData",
    ua: "cffi.FFI.CData",
    va: "cffi.FFI.CData",
    uc: "cffi.FFI.CData",
    vc: "cffi.FFI.CData",
    ak: "cffi.FFI.CData",
    bk: "cffi.FFI.CData",
    mfx: "cffi.FFI.CData",
    mfy: "cffi.FFI.CData",
    cx: "cffi.FFI.CData",
    cy: "cffi.FFI.CData",
    diss_est: "cffi.FFI.CData",
):
    global GEOS_DYCORE
    if not GEOS_DYCORE:
        raise RuntimeError(
            "[GEOS WRAPPER] Bad init, did you call init?"
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


def geos_gtfv3_finalize():
    if GEOS_DYCORE is not None:
        GEOS_DYCORE.finalize()


def geos_gtfv3_init(
    comm: MPI.Intercomm,
    npx: int,
    npy: int,
    npz: int,
    ntiles: int,
    is_: int,
    ie: int,
    js: int,
    je: int,
    isd: int,
    ied: int,
    jsd: int,
    jed: int,
    bdt,
    nq_tot: int,
):
    # Read in the backend
    BACKEND = os.environ.get("GTFV3_BACKEND", "gt:gpu")

    # Read in the namelist
    NAMELIST_PATH = os.environ.get("GTFV3_NAMELIST", "input.nml")

    global GEOS_DYCORE
    if GEOS_DYCORE is not None:
        raise RuntimeError("[GEOS WRAPPER] Double init")
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
