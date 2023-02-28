import cffi
import numpy as np
from math import prod

class FortranPythonConversion:
    '''
    Convert Fortran arrays to NumPy and vice-versa
    '''

    def __init__(
        self,
        npx, npy, npz,
        is_, ie, js, je,
        isd, ied, jsd, jed,
        num_tracers):

        self._npx, self._npy, self._npz = npx, npy, npz
        self._is, self._ie, self._js, self._je = is_, ie, js, je
        self._isd, self._ied, self._jsd, self._jed = isd, ied, jsd, jed
        assert num_tracers == 7, f'Expected 7 tracers, received: {num_tracers}'
        self._num_tracers = num_tracers

        self._ffi = cffi.FFI()
        self._TYPEMAP = {
            'float': np.float32,
            'double': np.float64,
            'int': np.int32,}

    def _fort_to_numpy(self, fptr, dim):
        '''
        Input: Fortran data pointed to by fptr and of shape dim = (i, j, k)
        Output: C-ordered double precision NumPy data of shape (i, j, k)
        '''
        ftype = self._ffi.getctype(self._ffi.typeof(fptr).item)
        assert ftype in self._TYPEMAP
        return np.frombuffer(
            self._ffi.buffer(fptr, prod(dim)*self._ffi.sizeof(ftype)),
            self._TYPEMAP[ftype],
        ).reshape(tuple(reversed(dim))).transpose().astype(np.float64)

    def fortran_to_numpy(
        self,
        # input
        u_ptr, v_ptr, w_ptr, delz_ptr,
        pt_ptr, delp_ptr, q_ptr,
        ps_ptr, pe_ptr, pk_ptr, peln_ptr, pkz_ptr,
        phis_ptr, q_con_ptr, omga_ptr,
        ua_ptr, va_ptr, uc_ptr, vc_ptr,
        mfxd_ptr, mfyd_ptr, cxd_ptr, cyd_ptr, diss_estd_ptr):
        '''
        Convert Fortran arrays pointed to by *_ptr to NumPy arrays
        Input: Pointers to Fortran arrays *_ptr
        Output: dict where dict[key] is a NumPy array
        '''

        # Shorthands
        is_, ie, js, je = self._is, self._ie, self._js, self._je
        isd, ied, jsd, jed = self._isd, self._ied, self._jsd, self._jed
        npz, num_tracers = self._npz, self._num_tracers

        # q/pe/peln require special handling
        # pe/peln need to be have their axes swapped - (i, k, j) -> (i, j, k)
        q = self._fort_to_numpy(q_ptr, (ied-isd+1, jed-jsd+1, npz, num_tracers))
        pe = np.swapaxes(self._fort_to_numpy(pe_ptr, (ie+1-(is_-1)+1, npz+1, je+1-(js-1)+1)), 1, 2)
        peln = np.swapaxes(self._fort_to_numpy(peln_ptr, (ie-is_+1, npz+1, je-js+1)), 1, 2)
        
        numpy_state = {
            'u': self._fort_to_numpy(u_ptr, (ied-isd+1, jed+1-jsd+1, npz)),
            'v': self._fort_to_numpy(v_ptr, (ied+1-isd+1, jed-jsd+1, npz)),
            'w': self._fort_to_numpy(w_ptr, (ied-isd+1, jed-jsd+1, npz)),
            'delz': self._fort_to_numpy(delz_ptr, (ied-isd+1, jed-jsd+1, npz)),

            'pt': self._fort_to_numpy(pt_ptr, (ied-isd+1, jed-jsd+1, npz)),
            'delp': self._fort_to_numpy(delp_ptr, (ied-isd+1, jed-jsd+1, npz)),
            'q': q,

            'ps': self._fort_to_numpy(ps_ptr, (ied-isd+1, jed-jsd+1)),
            'pe': pe,
            'pk': self._fort_to_numpy(pk_ptr, (ie-is_+1, je-js+1, npz+1)),
            'pkz': self._fort_to_numpy(pkz_ptr, (ie-is_+1, je-js+1, npz)),
            'peln': peln,

            'phis': self._fort_to_numpy(phis_ptr, (ied-isd+1, jed-jsd+1)),
            'q_con': self._fort_to_numpy(q_con_ptr, (ied-isd+1, jed-jsd+1, npz)),
            'omga': self._fort_to_numpy(omga_ptr, (ied-isd+1, jed-jsd+1, npz)),

            'ua': self._fort_to_numpy(ua_ptr, (ied-isd+1, jed-jsd+1, npz)),
            'va': self._fort_to_numpy(va_ptr, (ied-isd+1, jed-jsd+1, npz)),
            'uc': self._fort_to_numpy(uc_ptr, (ied+1-isd+1, jed-jsd+1, npz)),
            'vc': self._fort_to_numpy(va_ptr, (ied-isd+1, jed+1-jsd+1, npz)),

            'mfxd': self._fort_to_numpy(mfxd_ptr, (ie+1-is_+1, je-js+1, npz)),
            'mfyd': self._fort_to_numpy(mfyd_ptr, (ie-is_+1, je+1-js+1, npz)),
            'cxd': self._fort_to_numpy(cxd_ptr, (ie+1-is_+1, jed-jsd+1, npz)),
            'cyd': self._fort_to_numpy(cyd_ptr, (ied-isd+1, je+1-js+1, npz)),
            'diss_estd': self._fort_to_numpy(diss_estd_ptr, (ied-isd+1, jed-jsd+1, npz)),}

        return numpy_state # output

    def numpy_to_fortran(
        self,
        # input
        numpy_state,
        # output
        u_ptr, v_ptr, w_ptr, delz_ptr,
        pt_ptr, delp_ptr, q_ptr,
        ps_ptr, pe_ptr, pk_ptr, peln_ptr, pkz_ptr,
        phis_ptr, q_con_ptr, omga_ptr,
        ua_ptr, va_ptr, uc_ptr, vc_ptr,
        mfxd_ptr, mfyd_ptr, cxd_ptr, cyd_ptr, diss_estd_ptr):
        '''
        dp->sp, transpose, swap axes, numpy -> fortran
        '''

        # u/v/w/delz
        u_out_f = numpy_state['u'].astype(np.float32).flatten(order='F')
        self._ffi.memmove(u_ptr, u_out_f, 4*u_out_f.size)
        v_out_f = numpy_state['v'].astype(np.float32).flatten(order='F')
        self._ffi.memmove(v_ptr, v_out_f, 4*v_out_f.size)
        w_out_f = numpy_state['w'].astype(np.float32).flatten(order='F')
        self._ffi.memmove(w_ptr, w_out_f, 4*w_out_f.size)
        delz_out_f = numpy_state['delz'].astype(np.float32).flatten(order='F')
        self._ffi.memmove(delz_ptr, delz_out_f, 4*delz_out_f.size)

        # pt/delp/q
        pt_out_f = numpy_state['pt'].astype(np.float32).flatten(order='F')
        self._ffi.memmove(pt_ptr, pt_out_f, 4*pt_out_f.size)
        delp_out_f = numpy_state['delp'].astype(np.float32).flatten(order='F')
        self._ffi.memmove(delp_ptr, delp_out_f, 4*delp_out_f.size)
        # q needs special handling
        q = np.empty(list(numpy_state['qvapor'].shape)+[self._num_tracers])
        q[:,:,:,0] = numpy_state['qvapor']
        q[:,:,:,1] = numpy_state['qliquid']
        q[:,:,:,2] = numpy_state['qice']
        q[:,:,:,3] = numpy_state['qrain']
        q[:,:,:,4] = numpy_state['qsnow']
        q[:,:,:,5] = numpy_state['qgraupel']
        q[:,:,:,6] = numpy_state['qcld']
        q_out_f = q.astype(np.float32).flatten(order='F')
        self._ffi.memmove(q_ptr, q_out_f, 4*q_out_f.size)

        # ps/pe/pk/peln/pkz
        ps_out_f = numpy_state['ps'].astype(np.float32).flatten(order='F')
        self._ffi.memmove(ps_ptr, ps_out_f, 4*ps_out_f.size)
        pe_out_tmp = np.swapaxes(numpy_state['pe'], 1, 2) # (i, j, k) -> (i, k, j)
        pe_out_f = pe_out_tmp.astype(np.float32).flatten(order='F')
        self._ffi.memmove(pe_ptr, pe_out_f, 4*pe_out_f.size)
        pk_out_f = numpy_state['pk'].astype(np.float32).flatten(order='F')
        self._ffi.memmove(pk_ptr, pk_out_f, 4*pk_out_f.size)
        peln_out_tmp = np.swapaxes(numpy_state['peln'], 1, 2) # (i, j, k) -> (i, k, j)
        peln_out_f = peln_out_tmp.astype(np.float32).flatten(order='F')
        self._ffi.memmove(peln_ptr, peln_out_f, 4*peln_out_f.size)
        pkz_out_f = numpy_state['pkz'].astype(np.float32).flatten(order='F')
        self._ffi.memmove(pkz_ptr, pkz_out_f, 4*pkz_out_f.size)

        # phis/q_con/omga
        phis_out_f = numpy_state['phis'].astype(np.float32).flatten(order='F')
        self._ffi.memmove(phis_ptr, phis_out_f, 4*phis_out_f.size)
        q_con_out_f = numpy_state['q_con'].astype(np.float32).flatten(order='F')
        self._ffi.memmove(q_con_ptr, q_con_out_f, 4*q_con_out_f.size)
        omga_out_f = numpy_state['omga'].astype(np.float32).flatten(order='F')
        self._ffi.memmove(omga_ptr, omga_out_f, 4*omga_out_f.size)

        # ua/va/uc/vc
        ua_out_f = numpy_state['ua'].astype(np.float32).flatten(order='F')
        self._ffi.memmove(ua_ptr, ua_out_f, 4*ua_out_f.size)
        va_out_f = numpy_state['va'].astype(np.float32).flatten(order='F')
        self._ffi.memmove(va_ptr, va_out_f, 4*va_out_f.size)
        uc_out_f = numpy_state['uc'].astype(np.float32).flatten(order='F')
        self._ffi.memmove(uc_ptr, uc_out_f, 4*uc_out_f.size)
        vc_out_f = numpy_state['vc'].astype(np.float32).flatten(order='F')
        self._ffi.memmove(vc_ptr, vc_out_f, 4*vc_out_f.size)

        # mfx/mfy/cx/cy/diss_est
        mfxd_out_f = numpy_state['mfxd'].astype(np.float32).flatten(order='F')
        self._ffi.memmove(mfxd_ptr, mfxd_out_f, 4*mfxd_out_f.size)
        mfyd_out_f = numpy_state['mfyd'].astype(np.float32).flatten(order='F')
        self._ffi.memmove(mfyd_ptr, mfyd_out_f, 4*mfyd_out_f.size)
        cxd_out_f = numpy_state['cxd'].astype(np.float32).flatten(order='F')
        self._ffi.memmove(cxd_ptr, cxd_out_f, 4*cxd_out_f.size)
        cyd_out_f = numpy_state['cyd'].astype(np.float32).flatten(order='F')
        self._ffi.memmove(cyd_ptr, cyd_out_f, 4*cyd_out_f.size)
        diss_estd_out_f = numpy_state['diss_estd'].astype(np.float32).flatten(order='F')
        self._ffi.memmove(diss_estd_ptr, diss_estd_out_f, 4*diss_estd_out_f.size)
