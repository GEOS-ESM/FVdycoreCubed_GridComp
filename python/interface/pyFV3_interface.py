import os
from f_py_conversion import FortranPythonConversion
from cuda_profiler import CUDAProfiler, TimedCUDAProfiler
from mpi4py import MPI
from ndsl.optional_imports import cupy as cp
import numpy as np
from ndsl.dsl.gt4py_utils import is_gpu_backend
from typing import TYPE_CHECKING
from pyFV3_wrapper import GeosDycoreWrapper, MemorySpace
from fv_flags import FVFlags


if TYPE_CHECKING:
    import cffi


class PYFV3_WRAPPER:
    def __init__(
        self,
        fv_flags: FVFlags,
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
        tracer_count: int,
        ak_cdata: "cffi.FFI.CData",
        bk_cdata: "cffi.FFI.CData",
        backend: str = "dace:gpu",
    ) -> None:
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
            tracer_count,
            numpy_module,
        )

        # Input pressure levels
        ak = self.f_py._fortran_to_numpy(ak_cdata, [npz + 1])
        bk = self.f_py._fortran_to_numpy(bk_cdata, [npz + 1])

        # Setup pyFV3's dynamical core
        self.dycore = GeosDycoreWrapper(
            fv_flags=fv_flags,
            bdt=bdt,
            comm=comm,
            ak=ak,
            bk=bk,
            backend=self.backend,
            tracer_count=tracer_count,
            fortran_mem_space=fortran_mem_space,
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
WRAPPER = None


def pyfv3_run(
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
    global WRAPPER
    if not WRAPPER:
        raise RuntimeError("[GEOS WRAPPER] Bad init, did you call init?")
    WRAPPER(
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


def pyfv3_finalize():
    if WRAPPER is not None:
        WRAPPER.finalize()


def pyfv3_init(
    fv_flags: FVFlags,
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
    ak: "cffi.FFI.CData",
    bk: "cffi.FFI.CData",
):
    # Read in the backend
    BACKEND = os.environ.get("GEOS_PYFV3_BACKEND", "gt:gpu")

    global WRAPPER
    if WRAPPER is not None:
        raise RuntimeError("[GEOS WRAPPER] Double init")
    WRAPPER = PYFV3_WRAPPER(
        fv_flags=fv_flags,
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
        tracer_count=nq_tot,
        backend=BACKEND,
        ak_cdata=ak,
        bk_cdata=bk,
    )
