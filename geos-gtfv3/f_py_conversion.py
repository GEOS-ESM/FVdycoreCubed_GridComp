import cffi
import numpy as np
from math import prod
from cuda_profiler import CUDAProfiler
from pace.dsl.typing import Float


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
    ):

        self._npx, self._npy, self._npz = npx, npy, npz
        self._is, self._ie, self._js, self._je = is_, ie, js, je
        self._isd, self._ied, self._jsd, self._jed = isd, ied, jsd, jed
        self._num_tracers = num_tracers

        self._ffi = cffi.FFI()
        self._TYPEMAP = {
            "float": np.float32,
            "double": np.float64,
            "int": np.int32,
        }

    def _fortran_to_numpy_trf(self, fptr, dim):
        """
        Input: Fortran data pointed to by fptr and of shape dim = (i, j, k)
        Output: C-ordered double precision NumPy data of shape (i, j, k)
        """
        ftype = self._ffi.getctype(self._ffi.typeof(fptr).item)
        assert ftype in self._TYPEMAP
        if Float == self._TYPEMAP[ftype]:
            numpy_arr = (
                np.frombuffer(
                    self._ffi.buffer(fptr, prod(dim) * self._ffi.sizeof(ftype)),
                    self._TYPEMAP[ftype],
                )
                .reshape(tuple(reversed(dim)))
                .transpose()
            )
        else:
            numpy_arr = (
                np.frombuffer(
                    self._ffi.buffer(fptr, prod(dim) * self._ffi.sizeof(ftype)),
                    self._TYPEMAP[ftype],
                )
                .reshape(tuple(reversed(dim)))
                .transpose()
                .astype(Float)
            )
        return numpy_arr

    def fortran_to_numpy(
        self,
        # input
        u_ptr,
        v_ptr,
        w_ptr,
        delz_ptr,
        pt_ptr,
        delp_ptr,
        q_ptr,
        ps_ptr,
        pe_ptr,
        pk_ptr,
        peln_ptr,
        pkz_ptr,
        phis_ptr,
        q_con_ptr,
        omga_ptr,
        ua_ptr,
        va_ptr,
        uc_ptr,
        vc_ptr,
        mfxd_ptr,
        mfyd_ptr,
        cxd_ptr,
        cyd_ptr,
        diss_estd_ptr,
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
        q = self._fortran_to_numpy_trf(
            q_ptr, (ied - isd + 1, jed - jsd + 1, npz, num_tracers)
        )
        pe = np.swapaxes(
            self._fortran_to_numpy_trf(
                pe_ptr, (ie + 1 - (is_ - 1) + 1, npz + 1, je + 1 - (js - 1) + 1)
            ),
            1,
            2,
        )
        peln = np.swapaxes(
            self._fortran_to_numpy_trf(peln_ptr, (ie - is_ + 1, npz + 1, je - js + 1)),
            1,
            2,
        )

        numpy_state = {
            "u": self._fortran_to_numpy_trf(
                u_ptr, (ied - isd + 1, jed + 1 - jsd + 1, npz)
            ),
            "v": self._fortran_to_numpy_trf(
                v_ptr, (ied + 1 - isd + 1, jed - jsd + 1, npz)
            ),
            "w": self._fortran_to_numpy_trf(w_ptr, (ied - isd + 1, jed - jsd + 1, npz)),
            "delz": self._fortran_to_numpy_trf(
                delz_ptr, (ied - isd + 1, jed - jsd + 1, npz)
            ),
            "pt": self._fortran_to_numpy_trf(
                pt_ptr, (ied - isd + 1, jed - jsd + 1, npz)
            ),
            "delp": self._fortran_to_numpy_trf(
                delp_ptr, (ied - isd + 1, jed - jsd + 1, npz)
            ),
            "q": q,
            "ps": self._fortran_to_numpy_trf(ps_ptr, (ied - isd + 1, jed - jsd + 1)),
            "pe": pe,
            "pk": self._fortran_to_numpy_trf(
                pk_ptr, (ie - is_ + 1, je - js + 1, npz + 1)
            ),
            "pkz": self._fortran_to_numpy_trf(
                pkz_ptr, (ie - is_ + 1, je - js + 1, npz)
            ),
            "peln": peln,
            "phis": self._fortran_to_numpy_trf(
                phis_ptr, (ied - isd + 1, jed - jsd + 1)
            ),
            "q_con": self._fortran_to_numpy_trf(
                q_con_ptr, (ied - isd + 1, jed - jsd + 1, npz)
            ),
            "omga": self._fortran_to_numpy_trf(
                omga_ptr, (ied - isd + 1, jed - jsd + 1, npz)
            ),
            "ua": self._fortran_to_numpy_trf(
                ua_ptr, (ied - isd + 1, jed - jsd + 1, npz)
            ),
            "va": self._fortran_to_numpy_trf(
                va_ptr, (ied - isd + 1, jed - jsd + 1, npz)
            ),
            "uc": self._fortran_to_numpy_trf(
                uc_ptr, (ied + 1 - isd + 1, jed - jsd + 1, npz)
            ),
            "vc": self._fortran_to_numpy_trf(
                va_ptr, (ied - isd + 1, jed + 1 - jsd + 1, npz)
            ),
            "mfxd": self._fortran_to_numpy_trf(
                mfxd_ptr, (ie + 1 - is_ + 1, je - js + 1, npz)
            ),
            "mfyd": self._fortran_to_numpy_trf(
                mfyd_ptr, (ie - is_ + 1, je + 1 - js + 1, npz)
            ),
            "cxd": self._fortran_to_numpy_trf(
                cxd_ptr, (ie + 1 - is_ + 1, jed - jsd + 1, npz)
            ),
            "cyd": self._fortran_to_numpy_trf(
                cyd_ptr, (ied - isd + 1, je + 1 - js + 1, npz)
            ),
            "diss_estd": self._fortran_to_numpy_trf(
                diss_estd_ptr, (ied - isd + 1, jed - jsd + 1, npz)
            ),
        }

        return numpy_state  # output

    def _numpy_to_fortran_trf(self, numpy_array, fptr, ptr_offset=0):
        """
        Input: Fortran data pointed to by fptr and of shape dim = (i, j, k)
        Output: C-ordered double precision NumPy data of shape (i, j, k)
        """
        ftype = self._ffi.getctype(self._ffi.typeof(fptr).item)
        assert ftype in self._TYPEMAP
        if Float == self._TYPEMAP[ftype]:
            out_f = numpy_array.flatten(order="F")
            self._ffi.memmove(fptr + ptr_offset, out_f, 4 * out_f.size)
        else:
            out_f = numpy_array.astype(self._TYPEMAP[ftype]).flatten(order="F")
            self._ffi.memmove(fptr + ptr_offset, out_f, 4 * out_f.size)

    def numpy_to_fortran(
        self,
        # input
        numpy_state,
        # output
        u_ptr,
        v_ptr,
        w_ptr,
        delz_ptr,
        pt_ptr,
        delp_ptr,
        q_ptr,
        ps_ptr,
        pe_ptr,
        pk_ptr,
        peln_ptr,
        pkz_ptr,
        phis_ptr,
        q_con_ptr,
        omga_ptr,
        ua_ptr,
        va_ptr,
        uc_ptr,
        vc_ptr,
        mfxd_ptr,
        mfyd_ptr,
        cxd_ptr,
        cyd_ptr,
        diss_estd_ptr,
    ):
        """
        dp->sp, transpose, swap axes, numpy -> fortran
        """

        with CUDAProfiler("u/v/w/delz"):
            # u/v/w/delz
            self._numpy_to_fortran_trf(numpy_state["u"], u_ptr)
            self._numpy_to_fortran_trf(numpy_state["v"], v_ptr)
            self._numpy_to_fortran_trf(numpy_state["w"], w_ptr)
            self._numpy_to_fortran_trf(numpy_state["delz"], delz_ptr)

        with CUDAProfiler("pt/delp/q"):
            # pt/delp/q
            self._numpy_to_fortran_trf(numpy_state["pt"], pt_ptr)
            self._numpy_to_fortran_trf(numpy_state["delp"], delp_ptr)

        with CUDAProfiler("q"):
            # Dev Note: you should be able to unroll the below code in ptr + offset
            # since we are using Fortran layout (column-first)
            # q needs special handling
            # self.q = np.empty(list(numpy_state["qvapor"].shape) + [self._num_tracers])
            # self.q[:, :, :, 0] = numpy_state["qvapor"]
            # self.q[:, :, :, 1] = numpy_state["qliquid"]
            # self.q[:, :, :, 2] = numpy_state["qice"]
            # self.q[:, :, :, 3] = numpy_state["qrain"]
            # self.q[:, :, :, 4] = numpy_state["qsnow"]
            # self.q[:, :, :, 5] = numpy_state["qgraupel"]
            # self.q[:, :, :, 6] = numpy_state["qcld"]
            # self._numpy_to_fortran_trf(self.q, q_ptr)

            self._numpy_to_fortran_trf(numpy_state["qvapor"], q_ptr)
            offset = numpy_state["qvapor"].size
            self._numpy_to_fortran_trf(
                numpy_state["qliquid"],
                q_ptr,
                ptr_offset=offset,
            )
            offset += numpy_state["qliquid"].size
            self._numpy_to_fortran_trf(
                numpy_state["qice"],
                q_ptr,
                ptr_offset=offset,
            )
            offset += numpy_state["qice"].size
            self._numpy_to_fortran_trf(
                numpy_state["qrain"],
                q_ptr,
                ptr_offset=offset,
            )
            offset += numpy_state["qrain"].size
            self._numpy_to_fortran_trf(
                numpy_state["qsnow"],
                q_ptr,
                ptr_offset=offset,
            )
            offset += numpy_state["qsnow"].size
            self._numpy_to_fortran_trf(
                numpy_state["qgraupel"],
                q_ptr,
                ptr_offset=offset,
            )
            offset += numpy_state["qgraupel"].size
            self._numpy_to_fortran_trf(
                numpy_state["qcld"],
                q_ptr,
                ptr_offset=offset,
            )

        with CUDAProfiler("ps/pe/pk/peln/pkz"):
            # ps/pe/pk/peln/pkz
            self._numpy_to_fortran_trf(numpy_state["ps"], ps_ptr)
            self._numpy_to_fortran_trf(numpy_state["pk"], pk_ptr)
            self._numpy_to_fortran_trf(numpy_state["pkz"], pkz_ptr)
            #   pe
            pe_out_tmp = np.swapaxes(numpy_state["pe"], 1, 2)  # (i, j, k) -> (i, k, j)
            pe_out_f = pe_out_tmp.astype(np.float32).flatten(order="F")
            self._ffi.memmove(pe_ptr, pe_out_f, 4 * pe_out_f.size)
            #   peln
            peln_out_tmp = np.swapaxes(
                numpy_state["peln"], 1, 2
            )  # (i, j, k) -> (i, k, j)
            peln_out_f = peln_out_tmp.astype(np.float32).flatten(order="F")
            self._ffi.memmove(peln_ptr, peln_out_f, 4 * peln_out_f.size)

        with CUDAProfiler("phis/q_con/omga"):
            # phis/q_con/omga
            self._numpy_to_fortran_trf(numpy_state["phis"], phis_ptr)
            self._numpy_to_fortran_trf(numpy_state["q_con"], q_con_ptr)
            self._numpy_to_fortran_trf(numpy_state["omga"], omga_ptr)

        with CUDAProfiler("ua/va/uc/vc"):
            # ua/va/uc/vc
            self._numpy_to_fortran_trf(numpy_state["ua"], ua_ptr)
            self._numpy_to_fortran_trf(numpy_state["va"], va_ptr)
            self._numpy_to_fortran_trf(numpy_state["uc"], uc_ptr)
            self._numpy_to_fortran_trf(numpy_state["vc"], vc_ptr)

        with CUDAProfiler("mfx/mfy/cx/cy/diss_est"):
            # mfx/mfy/cx/cy/diss_est
            self._numpy_to_fortran_trf(numpy_state["mfxd"], mfxd_ptr)
            self._numpy_to_fortran_trf(numpy_state["mfyd"], mfyd_ptr)
            self._numpy_to_fortran_trf(numpy_state["cxd"], cxd_ptr)
            self._numpy_to_fortran_trf(numpy_state["cyd"], cyd_ptr)
            self._numpy_to_fortran_trf(numpy_state["diss_estd"], diss_estd_ptr)
