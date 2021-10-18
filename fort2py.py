import numpy as np
import gt4py
from math import prod
import cffi

TYPEMAP = {
    'float': np.dtype('f4'),
    'int': np.dtype('i4'),}

def fort_to_numpy(ffi, ptr, dim):
    ftype = ffi.getctype(ffi.typeof(ptr).item)
    assert ftype in TYPEMAP
    return np.frombuffer(
        ffi.buffer(ptr, prod(dim)*ffi.sizeof(ftype)),
        TYPEMAP[ftype],
    ).reshape(tuple(reversed(dim)))

def fort_to_gt4py(ptr, dim, origin, backend):
    ffi = cffi.FFI()
    nparr = fort_to_numpy(ffi, ptr, dim)
    # return gt4py.storage.from_array(nparr, backend, origin)
    return gt4py.storage.wrap_cpu_array(nparr, backend, origin)
