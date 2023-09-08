import cffi
import numpy as np
from math import prod
from cuda_profiler import CUDAProfiler
from pace.dsl.typing import Float
from typing import Tuple, Optional, List, Dict, Union
from types import ModuleType
from pace.util._optional_imports import cupy as cp

DeviceArray = cp.ndarray if cp else None
PythonArray = Union[np.ndarray, (cp.ndarray if cp else None)]


class DummyStream:
    def __init__(self) -> None:
        pass

    def synchronize(self):
        pass


class FortranPythonConversion:
    """
    Convert Fortran arrays to NumPy and vice-versa
    """

    def __init__(
        self,
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
        num_tracers: int,
        numpy_module: ModuleType,
    ):
        # Python numpy-like module is given by the caller leaving
        # optional control of upload/download in the case
        # of GPU/CPU system
        self._target_np = numpy_module

        # Device parameters
        #  Pace targets gpu: we want the Pace layout to be on device
        self._python_targets_gpu = self._target_np == cp
        if self._python_targets_gpu:
            self._stream_A = cp.cuda.Stream(non_blocking=True)
            self._stream_B = cp.cuda.Stream(non_blocking=True)
        else:
            self._stream_A = DummyStream()
            self._stream_B = DummyStream()
        self._current_stream = self._stream_A

        # Layout & indexing
        self._npx, self._npy, self._npz = npx, npy, npz
        self._is, self._ie, self._js, self._je = is_, ie, js, je
        self._isd, self._ied, self._jsd, self._jed = isd, ied, jsd, jed
        assert num_tracers == 7, f"Expected 7 tracers, received: {num_tracers}"
        self._num_tracers = num_tracers

        # cffi init
        self._ffi = cffi.FFI()
        self._TYPEMAP = {
            "float": np.float32,
            "double": np.float64,
            "int": np.int32,
        }

    def _device_sync(self):
        """Synchronize the working CUDA streams"""
        self._stream_A.synchronize()
        self._stream_B.synchronize()

    def _fortran_to_numpy(
        self,
        fptr: "cffi.FFI.CData",
        dim: List[int],
    ):
        """
        Input: Fortran data pointed to by fptr and of shape dim = (i, j, k)
        Output: C-ordered double precision NumPy data of shape (i, j, k)
        """
        ftype = self._ffi.getctype(self._ffi.typeof(fptr).item)
        assert ftype in self._TYPEMAP
        return np.frombuffer(
            self._ffi.buffer(fptr, prod(dim) * self._ffi.sizeof(ftype)),
            self._TYPEMAP[ftype],
        )

    def _upload_and_transform(
        self,
        host_array: np.ndarray,
        dim: List[int],
        swap_axes: Optional[Tuple[int, int]] = None,
    ) -> DeviceArray:
        """Upload to device & transform to Pace compatible layout"""
        with self._current_stream:
            device_array = cp.asarray(host_array)
            final_array = self._transform_from_fortran_layout(
                device_array,
                dim,
                swap_axes,
            )
            self._current_stream = (
                self._stream_A
                if self._current_stream == self._stream_B
                else self._stream_B
            )
            return final_array

    def _transform_from_fortran_layout(
        self,
        array: PythonArray,
        dim: List[int],
        swap_axes: Optional[Tuple[int, int]] = None,
    ) -> PythonArray:
        """Transform from Fortran layout into a Pace compatible layout"""
        trf_array = array.reshape(tuple(reversed(dim))).transpose().astype(Float)
        if swap_axes:
            trf_array = self._target_np.swapaxes(
                trf_array,
                swap_axes[0],
                swap_axes[1],
            )
        return trf_array

    def _fortran_to_python_trf(
        self,
        fptr: np.ndarray,
        dim: List[int],
        swap_axes: Optional[Tuple[int, int]] = None,
    ) -> PythonArray:
        """Move fortran memory into python space"""
        np_array = self._fortran_to_numpy(fptr, dim)
        if self._python_targets_gpu:
            return self._upload_and_transform(np_array, dim, swap_axes)
        else:
            return self._transform_from_fortran_layout(
                np_array,
                dim,
                swap_axes,
            )

    def fortran_to_python(
        self,
        # input
        u_ptr: "cffi.FFI.CData",
        v_ptr: "cffi.FFI.CData",
        w_ptr: "cffi.FFI.CData",
        delz_ptr: "cffi.FFI.CData",
        pt_ptr: "cffi.FFI.CData",
        delp_ptr: "cffi.FFI.CData",
        q_ptr: "cffi.FFI.CData",
        ps_ptr: "cffi.FFI.CData",
        pe_ptr: "cffi.FFI.CData",
        pk_ptr: "cffi.FFI.CData",
        peln_ptr: "cffi.FFI.CData",
        pkz_ptr: "cffi.FFI.CData",
        phis_ptr: "cffi.FFI.CData",
        q_con_ptr: "cffi.FFI.CData",
        omga_ptr: "cffi.FFI.CData",
        ua_ptr: "cffi.FFI.CData",
        va_ptr: "cffi.FFI.CData",
        uc_ptr: "cffi.FFI.CData",
        vc_ptr: "cffi.FFI.CData",
        mfxd_ptr: "cffi.FFI.CData",
        mfyd_ptr: "cffi.FFI.CData",
        cxd_ptr: "cffi.FFI.CData",
        cyd_ptr: "cffi.FFI.CData",
        diss_estd_ptr: "cffi.FFI.CData",
    ):
        """
        Convert Fortran arrays pointed to by *_ptr to NumPy arrays
        Input: Pointers to Fortran arrays *_ptr
        Output: dict where dict[key] is a NumPy array
        """
        # Shorthands
        is_, ie, js, je = self._is, self._ie, self._js, self._je
        isd, ied, jsd, jed = self._isd, self._ied, self._jsd, self._jed
        npz, num_tracers = self._npz, self._num_tracers

        # q/pe/peln require special handling
        # pe/peln need to be have their axes swapped - (i, k, j) -> (i, j, k)
        q = self._fortran_to_python_trf(
            q_ptr, (ied - isd + 1, jed - jsd + 1, npz, num_tracers)
        )
        pe = self._fortran_to_python_trf(
            pe_ptr,
            (ie + 1 - (is_ - 1) + 1, npz + 1, je + 1 - (js - 1) + 1),
            swap_axes=(1, 2),
        )
        peln = self._fortran_to_python_trf(
            peln_ptr,
            (ie - is_ + 1, npz + 1, je - js + 1),
            swap_axes=(1, 2),
        )

        python_state = {
            "u": self._fortran_to_python_trf(
                u_ptr, (ied - isd + 1, jed + 1 - jsd + 1, npz)
            ),
            "v": self._fortran_to_python_trf(
                v_ptr, (ied + 1 - isd + 1, jed - jsd + 1, npz)
            ),
            "w": self._fortran_to_python_trf(
                w_ptr, (ied - isd + 1, jed - jsd + 1, npz)
            ),
            "delz": self._fortran_to_python_trf(
                delz_ptr, (ied - isd + 1, jed - jsd + 1, npz)
            ),
            "pt": self._fortran_to_python_trf(
                pt_ptr, (ied - isd + 1, jed - jsd + 1, npz)
            ),
            "delp": self._fortran_to_python_trf(
                delp_ptr, (ied - isd + 1, jed - jsd + 1, npz)
            ),
            "q": q,
            "ps": self._fortran_to_python_trf(ps_ptr, (ied - isd + 1, jed - jsd + 1)),
            "pe": pe,
            "pk": self._fortran_to_python_trf(
                pk_ptr, (ie - is_ + 1, je - js + 1, npz + 1)
            ),
            "pkz": self._fortran_to_python_trf(
                pkz_ptr, (ie - is_ + 1, je - js + 1, npz)
            ),
            "peln": peln,
            "phis": self._fortran_to_python_trf(
                phis_ptr, (ied - isd + 1, jed - jsd + 1)
            ),
            "q_con": self._fortran_to_python_trf(
                q_con_ptr, (ied - isd + 1, jed - jsd + 1, npz)
            ),
            "omga": self._fortran_to_python_trf(
                omga_ptr, (ied - isd + 1, jed - jsd + 1, npz)
            ),
            "ua": self._fortran_to_python_trf(
                ua_ptr, (ied - isd + 1, jed - jsd + 1, npz)
            ),
            "va": self._fortran_to_python_trf(
                va_ptr, (ied - isd + 1, jed - jsd + 1, npz)
            ),
            "uc": self._fortran_to_python_trf(
                uc_ptr, (ied + 1 - isd + 1, jed - jsd + 1, npz)
            ),
            "vc": self._fortran_to_python_trf(
                va_ptr, (ied - isd + 1, jed + 1 - jsd + 1, npz)
            ),
            "mfxd": self._fortran_to_python_trf(
                mfxd_ptr, (ie + 1 - is_ + 1, je - js + 1, npz)
            ),
            "mfyd": self._fortran_to_python_trf(
                mfyd_ptr, (ie - is_ + 1, je + 1 - js + 1, npz)
            ),
            "cxd": self._fortran_to_python_trf(
                cxd_ptr, (ie + 1 - is_ + 1, jed - jsd + 1, npz)
            ),
            "cyd": self._fortran_to_python_trf(
                cyd_ptr, (ied - isd + 1, je + 1 - js + 1, npz)
            ),
            "diss_estd": self._fortran_to_python_trf(
                diss_estd_ptr, (ied - isd + 1, jed - jsd + 1, npz)
            ),
        }

        self._device_sync()

        return python_state  # output

    def _transform_and_download(
        self,
        device_array: DeviceArray,
        dtype: type,
        swap_axes: Optional[Tuple[int, int]] = None,
    ) -> np.ndarray:
        with self._current_stream:
            if swap_axes:
                device_array = cp.swapaxes(
                    device_array,
                    swap_axes[0],
                    swap_axes[1],
                )
            host_array = cp.asnumpy(
                device_array.astype(dtype).flatten(order="F"),
            )
            self._current_stream = (
                self._stream_A
                if self._current_stream == self._stream_B
                else self._stream_B
            )
            return host_array

    def _transform_from_python_layout(
        self,
        array: PythonArray,
        dtype: type,
        swap_axes: Optional[Tuple[int, int]] = None,
    ) -> np.ndarray:
        """Copy back a numpy array in python layout to Fortran"""

        if self._python_targets_gpu:
            numpy_array = self._transform_and_download(array, dtype, swap_axes)
        else:
            if swap_axes:
                numpy_array = np.swapaxes(
                    array,
                    swap_axes[0],
                    swap_axes[1],
                )
            else:
                numpy_array = array
            numpy_array = numpy_array.astype(dtype).flatten(order="F")

        return numpy_array

    def _python_to_fortran_trf(
        self,
        array: PythonArray,
        fptr: "cffi.FFI.CData",
        ptr_offset: int = 0,
        swap_axes: Optional[Tuple[int, int]] = None,
    ) -> np.ndarray:
        """
        Input: Fortran data pointed to by fptr and of shape dim = (i, j, k)
        Output: C-ordered double precision NumPy data of shape (i, j, k)
        """
        ftype = self._ffi.getctype(self._ffi.typeof(fptr).item)
        assert ftype in self._TYPEMAP
        dtype = self._TYPEMAP[ftype]
        numpy_array = self._transform_from_python_layout(
            array,
            dtype,
            swap_axes,
        )
        self._ffi.memmove(fptr + ptr_offset, numpy_array, 4 * numpy_array.size)

    def python_to_fortran(
        self,
        # input
        python_state: Dict[str, PythonArray],
        # output
        u_ptr: "cffi.FFI.CData",
        v_ptr: "cffi.FFI.CData",
        w_ptr: "cffi.FFI.CData",
        delz_ptr: "cffi.FFI.CData",
        pt_ptr: "cffi.FFI.CData",
        delp_ptr: "cffi.FFI.CData",
        q_ptr: "cffi.FFI.CData",
        ps_ptr: "cffi.FFI.CData",
        pe_ptr: "cffi.FFI.CData",
        pk_ptr: "cffi.FFI.CData",
        peln_ptr: "cffi.FFI.CData",
        pkz_ptr: "cffi.FFI.CData",
        phis_ptr: "cffi.FFI.CData",
        q_con_ptr: "cffi.FFI.CData",
        omga_ptr: "cffi.FFI.CData",
        ua_ptr: "cffi.FFI.CData",
        va_ptr: "cffi.FFI.CData",
        uc_ptr: "cffi.FFI.CData",
        vc_ptr: "cffi.FFI.CData",
        mfxd_ptr: "cffi.FFI.CData",
        mfyd_ptr: "cffi.FFI.CData",
        cxd_ptr: "cffi.FFI.CData",
        cyd_ptr: "cffi.FFI.CData",
        diss_estd_ptr: "cffi.FFI.CData",
    ) -> None:
        """
        dp->sp, transpose, swap axes, numpy -> fortran
        """

        with CUDAProfiler("u/v/w/delz"):
            # u/v/w/delz
            self._python_to_fortran_trf(python_state["u"], u_ptr)
            self._python_to_fortran_trf(python_state["v"], v_ptr)
            self._python_to_fortran_trf(python_state["w"], w_ptr)
            self._python_to_fortran_trf(python_state["delz"], delz_ptr)

        with CUDAProfiler("pt/delp/q"):
            # pt/delp/q
            self._python_to_fortran_trf(python_state["pt"], pt_ptr)
            self._python_to_fortran_trf(python_state["delp"], delp_ptr)

        with CUDAProfiler("q"):
            # Dev Note: you should be able to unroll the below code in ptr + offset
            # since we are using Fortran layout (column-first)
            # q needs special handling
            # self.q = np.empty(list(python_state["qvapor"].shape) + [self._num_tracers])
            # self.q[:, :, :, 0] = python_state["qvapor"]
            # self.q[:, :, :, 1] = python_state["qliquid"]
            # self.q[:, :, :, 2] = python_state["qice"]
            # self.q[:, :, :, 3] = python_state["qrain"]
            # self.q[:, :, :, 4] = python_state["qsnow"]
            # self.q[:, :, :, 5] = python_state["qgraupel"]
            # self.q[:, :, :, 6] = python_state["qcld"]
            # self._python_to_fortran_trf(self.q, q_ptr)

            self._python_to_fortran_trf(python_state["qvapor"], q_ptr)
            offset = python_state["qvapor"].size
            self._python_to_fortran_trf(
                python_state["qliquid"],
                q_ptr,
                ptr_offset=offset,
            )
            offset += python_state["qliquid"].size
            self._python_to_fortran_trf(
                python_state["qice"],
                q_ptr,
                ptr_offset=offset,
            )
            offset += python_state["qice"].size
            self._python_to_fortran_trf(
                python_state["qrain"],
                q_ptr,
                ptr_offset=offset,
            )
            offset += python_state["qrain"].size
            self._python_to_fortran_trf(
                python_state["qsnow"],
                q_ptr,
                ptr_offset=offset,
            )
            offset += python_state["qsnow"].size
            self._python_to_fortran_trf(
                python_state["qgraupel"],
                q_ptr,
                ptr_offset=offset,
            )
            offset += python_state["qgraupel"].size
            self._python_to_fortran_trf(
                python_state["qcld"],
                q_ptr,
                ptr_offset=offset,
            )

        with CUDAProfiler("ps/pe/pk/peln/pkz"):
            # ps/pe/pk/peln/pkz
            self._python_to_fortran_trf(python_state["ps"], ps_ptr)
            self._python_to_fortran_trf(python_state["pk"], pk_ptr)
            self._python_to_fortran_trf(python_state["pkz"], pkz_ptr)
            # pe: (i, j, k) -> (i, k, j)
            self._python_to_fortran_trf(
                python_state["pe"],
                pe_ptr,
                swap_axes=(1, 2),
            )
            # peln: (i, j, k) -> (i, k, j)
            self._python_to_fortran_trf(
                python_state["peln"],
                peln_ptr,
                swap_axes=(1, 2),
            )

        with CUDAProfiler("phis/q_con/omga"):
            # phis/q_con/omga
            self._python_to_fortran_trf(python_state["phis"], phis_ptr)
            self._python_to_fortran_trf(python_state["q_con"], q_con_ptr)
            self._python_to_fortran_trf(python_state["omga"], omga_ptr)

        with CUDAProfiler("ua/va/uc/vc"):
            # ua/va/uc/vc
            self._python_to_fortran_trf(python_state["ua"], ua_ptr)
            self._python_to_fortran_trf(python_state["va"], va_ptr)
            self._python_to_fortran_trf(python_state["uc"], uc_ptr)
            self._python_to_fortran_trf(python_state["vc"], vc_ptr)

        with CUDAProfiler("mfx/mfy/cx/cy/diss_est"):
            # mfx/mfy/cx/cy/diss_est
            self._python_to_fortran_trf(python_state["mfxd"], mfxd_ptr)
            self._python_to_fortran_trf(python_state["mfyd"], mfyd_ptr)
            self._python_to_fortran_trf(python_state["cxd"], cxd_ptr)
            self._python_to_fortran_trf(python_state["cyd"], cyd_ptr)
            self._python_to_fortran_trf(
                python_state["diss_estd"],
                diss_estd_ptr,
            )

        self._device_sync()
