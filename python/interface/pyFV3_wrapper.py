import enum
import logging
import os
from datetime import timedelta
from typing import Dict, List, Tuple

import numpy as np
from gt4py.cartesian.config import build_settings as gt_build_settings, GT4PY_COMPILE_OPT_LEVEL
from mpi4py import MPI

import pyfv3
from ndsl import (
    CompilationConfig,
    CubedSphereCommunicator,
    CubedSpherePartitioner,
    DaceConfig,
    GridIndexing,
    PerformanceCollector,
    QuantityFactory,
    StencilConfig,
    StencilFactory,
    SubtileGridSizer,
    TilePartitioner,
    orchestrate,
    Backend,
)
import ndsl.constants
from ndsl.comm.comm_abc import Comm
from ndsl.dsl.dace.build import set_distributed_caches
from ndsl.dsl.typing import get_precision
from ndsl.grid import DampingCoefficients, GridData, MetricTerms
from ndsl.logging import ndsl_log, ndsl_log_on_rank_0
from ndsl.optional_imports import cupy as cp
from ndsl.utils import safe_assign_array
from fv_flags import FVFlags, FVFlags_to_DycoreConfig
from ndsl.comm.mpi import MPIComm

TRACERS_IN_FORTRAN = {
    "vapor": 0,
    "liquid": 1,
    "ice": 2,
    "rain": 3,
    "snow": 4,
    "graupel": 5,
    "cloud": 6,
}


class StencilBackendCompilerOverride:
    """Override the NDSL global stencil JIT to allow for 9-rank build
    on any setup.

    This is a workaround that requires to know _exactly_ when build is happening.
    Using this as a context manager, we leverage the DaCe build system to override
    the name and build the 9 codepaths required- while every other rank wait.

    This should be removed when we refactor the GT JIT to distribute building
    much more efficiently
    """

    def __init__(self, comm: MPI.Intracomm, config: DaceConfig):
        self.comm = comm
        self.config = config

        # Orchestration or mono-node is not concerned
        self.no_op = self.config.is_dace_orchestrated() or self.comm.Get_size() == 1

        # We abuse the DaCe build system
        if not self.no_op:
            set_distributed_caches(config, force_build=True)

        # We remove warnings from the stencils compiling when in critical and/or
        # error
        if ndsl_log.level > logging.WARNING:
            gt_build_settings["extra_compile_args"]["cxx"].append("-w")
            gt_build_settings["extra_compile_args"]["cuda"].append("-w")

    def __enter__(self):
        if self.no_op:
            return
        if self.config.do_compile:
            ndsl_log.info(f"Stencil backend compiles on {self.comm.Get_rank()}")
        else:
            ndsl_log.info(f"Stencil backend waits on {self.comm.Get_rank()}")
            self.comm.Barrier()

    def __exit__(self, type, value, traceback):
        if self.no_op:
            return
        if not self.config.do_compile:
            ndsl_log.info(f"Stencil backend read cache on {self.comm.Get_rank()}")
        else:
            ndsl_log.info(f"Stencil backend compiled on {self.comm.Get_rank()}")
            self.comm.Barrier()


@enum.unique
class MemorySpace(enum.Enum):
    CPU = 0
    GPU = 1


class InterfaceTransferType(enum.Enum):
    CPU_COPY = enum.auto()  # Copies because of layout mismatch
    CPU_ZERO_COPY = enum.auto()  # No copy - reuse the memory given (case of same layout)
    GPU_TRANSFER = enum.auto()  # Upload from central RAM to device, then download
    GPU_MAPPING = enum.auto()  # Paged memory mapped onto GPU memory space (HMM, ATS, ...)


class GeosDycoreWrapper:
    """
    Provides an interface for the Geos model to access the Pace dycore.
    Takes numpy arrays as inputs, returns a dictionary of numpy arrays as outputs
    """

    def __init__(
        self,
        fv_flags: FVFlags,
        ak: np.ndarray,
        bk: np.ndarray,
        phis: np.ndarray,
        bdt: float,
        comm: Comm,
        backend: Backend,
        tracer_count: int,
        fortran_mem_space: MemorySpace = MemorySpace.CPU,
    ):
        ndsl_log_on_rank_0.info("Setup pyfv3 (including grid)")
        # pyfv3 support only 0 and 7 tracers run
        if fv_flags.nwat not in [6, 0]:
            raise NotImplementedError(
                f"pyFV3 requires 6 or 0 water tracers to be advected, {fv_flags.nwat} given. Abort."
            )
        comm = MPIComm()

        # Make a custom performance collector for the GEOS wrapper
        self.perf_collector = PerformanceCollector("GEOS wrapper", comm)

        self.dycore_config = pyfv3._config.DynamicalCoreConfig()
        FVFlags_to_DycoreConfig(fv_flags, self.dycore_config)

        self.dycore_config.dt_atmos = int(bdt)
        assert self.dycore_config.dt_atmos != 0

        self.layout = self.dycore_config.layout
        partitioner = CubedSpherePartitioner(TilePartitioner(self.layout))
        self.communicator = CubedSphereCommunicator(
            comm,
            partitioner,
            timer=self.perf_collector.timestep_timer,
        )
        sizer = SubtileGridSizer.from_tile_params(
            nx_tile=self.dycore_config.npx - 1,  # NX/NY from config are cell-centers
            ny_tile=self.dycore_config.npy - 1,  # NX/NY from config are cell-centers
            nz=self.dycore_config.npz,
            n_halo=ndsl.constants.N_HALO_DEFAULT,
            layout=self.dycore_config.layout,
            tile_partitioner=partitioner.tile,
            tile_rank=self.communicator.tile.rank,
            backend=backend,
            pad_non_interface_dimensions=True,
        )
        quantity_factory = QuantityFactory(sizer=sizer, backend=backend)

        # set up the metric terms and grid data
        metric_terms = MetricTerms(
            quantity_factory=quantity_factory,
            communicator=self.communicator,
            ak=ak,
            bk=bk,
        )
        grid_data = GridData.new_from_metric_terms(metric_terms)

        stencil_config = StencilConfig(
            compilation_config=CompilationConfig(
                backend=backend,
                rebuild=False,
                validate_args=False,
            ),
        )

        # Build a DaCeConfig for orchestration.
        # This and all orchestration code are transparent when outside
        # configuration deactivate orchestration
        stencil_config.dace_config = DaceConfig(
            communicator=self.communicator,
            backend=stencil_config.backend,
            tile_nx=self.dycore_config.npx,
            tile_nz=self.dycore_config.npz,
        )

        self._grid_indexing = GridIndexing.from_sizer_and_communicator(sizer=sizer, comm=self.communicator)
        stencil_factory = StencilFactory(config=stencil_config, grid_indexing=self._grid_indexing)

        self.tracer_count = tracer_count
        if fv_flags.nwat == 6:
            tracer_names = TRACERS_IN_FORTRAN
        elif fv_flags.nwat == 0:
            tracer_names = {"vapor": 0}
        pyfv3.tracers.setup_fvtracers(quantity_factory, tracer_count, tracer_names)

        self.dycore_state = pyfv3.DycoreState.init_zeros(quantity_factory=quantity_factory)
        self.dycore_state.bdt = self.dycore_config.dt_atmos
        self.dycore_state.phis[:-1, :-1] = phis[:]

        damping_coefficients = DampingCoefficients.new_from_metric_terms(metric_terms)

        with StencilBackendCompilerOverride(MPI.COMM_WORLD, stencil_config.dace_config):
            self.dynamical_core = pyfv3.DynamicalCore(
                comm=self.communicator,
                grid_data=grid_data,
                stencil_factory=stencil_factory,
                quantity_factory=quantity_factory,
                damping_coefficients=damping_coefficients,
                config=self.dycore_config,
                timestep=timedelta(seconds=self.dycore_state.bdt),
                phis=self.dycore_state.phis,
                state=self.dycore_state,
            )

        self._fortran_mem_space = fortran_mem_space
        self._ndsl_mem_space = MemorySpace.GPU if backend.is_gpu_backend() else MemorySpace.CPU

        self.output_dict = self._allocate_output_dir()

        # Feedback information
        device_ordinal_info = (
            f"  Device PCI bus id: {cp.cuda.Device(0).pci_bus_id}\n" if backend.is_gpu_backend() else "N/A"
        )
        MPS_pipe_directory = os.getenv("CUDA_MPS_PIPE_DIRECTORY", None)
        MPS_is_on = (
            MPS_pipe_directory is not None
            and backend.is_gpu_backend()
            and os.path.exists(f"{MPS_pipe_directory}/log")
        )
        build_status = (
            f"Build : {os.environ['FV3_DACEMODE']}\n"
            if backend.is_orchestrated() and "FV3_DACEMODE" in os.environ
            else ""
        )
        ndsl_log_on_rank_0.info(
            "pyFV3 <> GEOS wrapper initialized (Rank 0):\n"
            f"           Bridge : {fortran_mem_space} <> {self._ndsl_mem_space}\n"
            f"          Backend : {backend}\n"
            f"            {build_status}"
            f"        Precision : {get_precision()} bit\n"
            f"     Optimization : -O{GT4PY_COMPILE_OPT_LEVEL}\n"
            f"     Local domain : {sizer.nx}x{sizer.ny}x{sizer.nz}"
            f"(halo: {sizer.n_halo})\n"
            f"           Layout : {partitioner.layout}\n"
            f"   Strides for 3D : {self.dycore_state.pt._data.strides}\n"
            f"       Device ord : {device_ordinal_info}\n"
            f"       Nvidia MPS : {MPS_is_on}\n"
        )

    def _run(self):
        with self.perf_collector.timestep_timer.clock("step_dynamics"):
            self.dynamical_core.step_dynamics(
                state=self.dycore_state,
                timer=self.perf_collector.timestep_timer,
            )

    def __call__(
        self,
        timings: dict[str, list[float | int]],
        u: np.ndarray,
        v: np.ndarray,
        w: np.ndarray,
        delz: np.ndarray,
        pt: np.ndarray,
        delp: np.ndarray,
        q: np.ndarray,
        ps: np.ndarray,
        pe: np.ndarray,
        pk: np.ndarray,
        peln: np.ndarray,
        pkz: np.ndarray,
        phis: np.ndarray,
        q_con: np.ndarray,
        omga: np.ndarray,
        ua: np.ndarray,
        va: np.ndarray,
        uc: np.ndarray,
        vc: np.ndarray,
        mfxd: np.ndarray,
        mfyd: np.ndarray,
        cxd: np.ndarray,
        cyd: np.ndarray,
        diss_estd: np.ndarray,
    ) -> Tuple[Dict[str, np.ndarray], Dict[str, List[float]]]:
        with self.perf_collector.timestep_timer.clock("numpy-to-dycore"):
            self.dycore_state = self._put_fortran_data_in_dycore(
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
                mfxd,
                mfyd,
                cxd,
                cyd,
                diss_estd,
            )
        # Enter orchestrated code - if applicable
        self._run()

        with self.perf_collector.timestep_timer.clock("dycore-to-numpy"):
            self.output_dict = self._prep_outputs_for_geos()

        # Collect performance of the timestep and write a json file for rank 0
        self.perf_collector.collect_performance()
        for k, v in self.perf_collector.times_per_step[0].items():
            if k not in timings.keys():
                timings[k] = [v]
            else:
                timings[k].append(v)
        self.perf_collector.clear()

        return self.output_dict, timings

    def _put_fortran_data_in_dycore(
        self,
        u: np.ndarray,
        v: np.ndarray,
        w: np.ndarray,
        delz: np.ndarray,
        pt: np.ndarray,
        delp: np.ndarray,
        q: np.ndarray,
        ps: np.ndarray,
        pe: np.ndarray,
        pk: np.ndarray,
        peln: np.ndarray,
        pkz: np.ndarray,
        phis: np.ndarray,
        q_con: np.ndarray,
        omga: np.ndarray,
        ua: np.ndarray,
        va: np.ndarray,
        uc: np.ndarray,
        vc: np.ndarray,
        mfxd: np.ndarray,
        mfyd: np.ndarray,
        cxd: np.ndarray,
        cyd: np.ndarray,
        diss_estd: np.ndarray,
    ) -> pyfv3.DycoreState:
        isc = self._grid_indexing.isc
        jsc = self._grid_indexing.jsc
        iec = self._grid_indexing.iec + 1
        jec = self._grid_indexing.jec + 1

        state = self.dycore_state

        # Assign compute domain:
        safe_assign_array(state.u.field[:], u[isc:iec, jsc : jec + 1, :])
        safe_assign_array(state.v.field[:], v[isc : iec + 1, jsc:jec, :])
        safe_assign_array(state.w.field[:], w[isc:iec, jsc:jec, :])
        safe_assign_array(state.ua.field[:], ua[isc:iec, jsc:jec, :])
        safe_assign_array(state.va.field[:], va[isc:iec, jsc:jec, :])
        safe_assign_array(state.uc.field[:], uc[isc : iec + 1, jsc:jec, :])
        safe_assign_array(state.vc.field[:], vc[isc:iec, jsc : jec + 1, :])

        safe_assign_array(state.delz.field[:], delz[isc:iec, jsc:jec, :])
        safe_assign_array(state.pt.field[:], pt[isc:iec, jsc:jec, :])
        safe_assign_array(state.delp.field[:], delp[isc:iec, jsc:jec, :])

        safe_assign_array(state.mfxd.field[:], mfxd)
        safe_assign_array(state.mfyd.field[:], mfyd)
        safe_assign_array(state.cxd.field[:], cxd[:, jsc:jec, :])
        safe_assign_array(state.cyd.field[:], cyd[isc:iec, :, :])

        safe_assign_array(state.ps.field[:], ps[isc:iec, jsc:jec])
        safe_assign_array(state.pe[isc - 1 : iec + 1, jsc - 1 : jec + 1, :], pe)
        safe_assign_array(state.pk.field[:], pk)
        safe_assign_array(state.peln.field[:], peln)
        safe_assign_array(state.pkz.field[:], pkz)
        safe_assign_array(state.phis.field[:], phis[isc:iec, jsc:jec])
        safe_assign_array(state.q_con.field[:], q_con[isc:iec, jsc:jec, :])
        safe_assign_array(state.omga.field[:], omga[isc:iec, jsc:jec, :])
        safe_assign_array(state.diss_estd.field[:], diss_estd[isc:iec, jsc:jec, :])

        # tracer quantities should be a 4d array in order:
        # vapor, liquid, ice, rain, snow, graupel, cloud
        # This codes works because Fortran as moved all those tracers at the top of the list
        # it will fail if they are not contiguous (second loop)
        state.tracers[isc:iec, jsc:jec, :-1, :] = q[isc:iec, jsc:jec, :, :]

        return state

    def _prep_outputs_for_geos(self) -> dict[str, np.ndarray]:
        output_dict = self.output_dict
        isc = self._grid_indexing.isc
        jsc = self._grid_indexing.jsc
        iec = self._grid_indexing.iec + 1
        jec = self._grid_indexing.jec + 1

        if self._fortran_mem_space != self._ndsl_mem_space:
            safe_assign_array(output_dict["u"], self.dycore_state.u[:-1, :, :-1])
            safe_assign_array(output_dict["v"], self.dycore_state.v[:, :-1, :-1])
            safe_assign_array(output_dict["w"], self.dycore_state.w[:-1, :-1, :-1])
            safe_assign_array(output_dict["ua"], self.dycore_state.ua[:-1, :-1, :-1])
            safe_assign_array(output_dict["va"], self.dycore_state.va[:-1, :-1, :-1])
            safe_assign_array(output_dict["uc"], self.dycore_state.uc[:, :-1, :-1])
            safe_assign_array(output_dict["vc"], self.dycore_state.vc[:-1, :, :-1])

            safe_assign_array(output_dict["delz"], self.dycore_state.delz[:-1, :-1, :-1])
            safe_assign_array(output_dict["pt"], self.dycore_state.pt[:-1, :-1, :-1])
            safe_assign_array(output_dict["delp"], self.dycore_state.delp[:-1, :-1, :-1])

            safe_assign_array(
                output_dict["mfxd"],
                self.dycore_state.mfxd[isc : iec + 1, jsc:jec, :-1],
            )
            safe_assign_array(
                output_dict["mfyd"],
                self.dycore_state.mfyd[isc:iec, jsc : jec + 1, :-1],
            )
            safe_assign_array(output_dict["cxd"], self.dycore_state.cxd[isc : iec + 1, :-1, :-1])
            safe_assign_array(output_dict["cyd"], self.dycore_state.cyd[:-1, jsc : jec + 1, :-1])

            safe_assign_array(output_dict["ps"], self.dycore_state.ps[:-1, :-1])
            safe_assign_array(
                output_dict["pe"],
                self.dycore_state.pe[isc - 1 : iec + 1, jsc - 1 : jec + 1, :],
            )
            safe_assign_array(output_dict["pk"], self.dycore_state.pk[isc:iec, jsc:jec, :])
            safe_assign_array(output_dict["peln"], self.dycore_state.peln[isc:iec, jsc:jec, :])
            safe_assign_array(output_dict["pkz"], self.dycore_state.pkz[isc:iec, jsc:jec, :-1])
            safe_assign_array(output_dict["phis"], self.dycore_state.phis[:-1, :-1])
            safe_assign_array(output_dict["q_con"], self.dycore_state.q_con[:-1, :-1, :-1])
            safe_assign_array(output_dict["omga"], self.dycore_state.omga[:-1, :-1, :-1])
            safe_assign_array(
                output_dict["diss_estd"],
                self.dycore_state.diss_estd[:-1, :-1, :-1],
            )

            # Tracers
            safe_assign_array(
                output_dict["tracers"],
                self.dycore_state.tracers.as_4D_array(),
            )
        else:
            output_dict["u"] = self.dycore_state.u[:-1, :, :-1]
            output_dict["v"] = self.dycore_state.v[:, :-1, :-1]
            output_dict["w"] = self.dycore_state.w[:-1, :-1, :-1]
            output_dict["ua"] = self.dycore_state.ua[:-1, :-1, :-1]
            output_dict["va"] = self.dycore_state.va[:-1, :-1, :-1]
            output_dict["uc"] = self.dycore_state.uc[:, :-1, :-1]
            output_dict["vc"] = self.dycore_state.vc[:-1, :, :-1]
            output_dict["delz"] = self.dycore_state.delz[:-1, :-1, :-1]
            output_dict["pt"] = self.dycore_state.pt[:-1, :-1, :-1]
            output_dict["delp"] = self.dycore_state.delp[:-1, :-1, :-1]
            output_dict["mfxd"] = self.dycore_state.mfxd[isc : iec + 1, jsc:jec, :-1]
            output_dict["mfyd"] = self.dycore_state.mfyd[isc:iec, jsc : jec + 1, :-1]
            output_dict["cxd"] = self.dycore_state.cxd[isc : iec + 1, :-1, :-1]
            output_dict["cyd"] = self.dycore_state.cyd[:-1, jsc : jec + 1, :-1]
            output_dict["ps"] = self.dycore_state.ps[:-1, :-1]
            output_dict["pe"] = self.dycore_state.pe[isc - 1 : iec + 1, jsc - 1 : jec + 1, :]
            output_dict["pk"] = self.dycore_state.pk[isc:iec, jsc:jec, :]
            output_dict["peln"] = self.dycore_state.peln[isc:iec, jsc:jec, :]
            output_dict["pkz"] = self.dycore_state.pkz[isc:iec, jsc:jec, :-1]
            output_dict["phis"] = self.dycore_state.phis[:-1, :-1]
            output_dict["q_con"] = self.dycore_state.q_con[:-1, :-1, :-1]
            output_dict["omga"] = self.dycore_state.omga[:-1, :-1, :-1]
            output_dict["diss_estd"] = self.dycore_state.diss_estd[:-1, :-1, :-1]

            # Tracers
            output_dict["tracers"] = self.dycore_state.tracers[:-1, :-1, :-1, :]
        return output_dict

    def _allocate_output_dir(self) -> dict[str, np.ndarray]:
        if self._fortran_mem_space == self._ndsl_mem_space:
            return {}

        nhalo = self._grid_indexing.n_halo
        shape_centered = self._grid_indexing.domain_full(add=(0, 0, 0))
        shape_x_interface = self._grid_indexing.domain_full(add=(1, 0, 0))
        shape_y_interface = self._grid_indexing.domain_full(add=(0, 1, 0))
        shape_2d = shape_centered[:-1]

        output_dict = {}
        output_dict["u"] = np.empty((shape_y_interface))
        output_dict["v"] = np.empty((shape_x_interface))
        output_dict["w"] = np.empty((shape_centered))
        output_dict["ua"] = np.empty((shape_centered))
        output_dict["va"] = np.empty((shape_centered))
        output_dict["uc"] = np.empty((shape_x_interface))
        output_dict["vc"] = np.empty((shape_y_interface))

        output_dict["delz"] = np.empty((shape_centered))
        output_dict["pt"] = np.empty((shape_centered))
        output_dict["delp"] = np.empty((shape_centered))

        output_dict["mfxd"] = np.empty((self._grid_indexing.domain_full(add=(1 - 2 * nhalo, -2 * nhalo, 0))))
        output_dict["mfyd"] = np.empty((self._grid_indexing.domain_full(add=(-2 * nhalo, 1 - 2 * nhalo, 0))))
        output_dict["cxd"] = np.empty((self._grid_indexing.domain_full(add=(1 - 2 * nhalo, 0, 0))))
        output_dict["cyd"] = np.empty((self._grid_indexing.domain_full(add=(0, 1 - 2 * nhalo, 0))))

        output_dict["ps"] = np.empty((shape_2d))
        output_dict["pe"] = np.empty((self._grid_indexing.domain_full(add=(2 - 2 * nhalo, 2 - 2 * nhalo, 1))))
        output_dict["pk"] = np.empty((self._grid_indexing.domain_full(add=(-2 * nhalo, -2 * nhalo, 1))))
        output_dict["peln"] = np.empty((self._grid_indexing.domain_full(add=(-2 * nhalo, -2 * nhalo, 1))))
        output_dict["pkz"] = np.empty((self._grid_indexing.domain_full(add=(-2 * nhalo, -2 * nhalo, 0))))
        output_dict["phis"] = np.empty((shape_2d))
        output_dict["q_con"] = np.empty((shape_centered))
        output_dict["omga"] = np.empty((shape_centered))
        output_dict["diss_estd"] = np.empty((shape_centered))

        output_dict["tracers"] = np.empty(
            (
                shape_centered[0],
                shape_centered[1],
                shape_centered[2],
                self.tracer_count,
            )
        )
        return output_dict
