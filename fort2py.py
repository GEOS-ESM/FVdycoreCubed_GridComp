import numpy as np
import gt4py
from math import prod
import cffi

TYPEMAP = {
    'float': np.float32,
    'double': np.float64,
    'int': np.int32,}

def fort_to_numpy(ffi, ptr, dim):
    ftype = ffi.getctype(ffi.typeof(ptr).item)
    assert ftype in TYPEMAP
    return np.frombuffer(
        ffi.buffer(ptr, prod(dim)*ffi.sizeof(ftype)),
        TYPEMAP[ftype],
    ).reshape(tuple(reversed(dim))).transpose().astype(np.float64)

def fort_to_gt4py(ptr, dim, origin, backend):
    ffi = cffi.FFI()
    nparr = fort_to_numpy(ffi, ptr, dim).transpose()
    return gt4py.storage.from_array(nparr, backend, origin)
    # return gt4py.storage.wrap_cpu_array(nparr, backend, origin)
