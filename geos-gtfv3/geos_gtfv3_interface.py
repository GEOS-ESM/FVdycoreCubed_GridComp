import cffi
from mpi4py import MPI

TMPFILEBASE = 'geos_gtfv3_interface_py'

ffi = cffi.FFI()

# MPI_Comm can be int or void*
if MPI._sizeof(MPI.Comm) == ffi.sizeof('int'):
    _mpi_comm_t = 'int'
else:
    _mpi_comm_t = 'void*'

source = '''
from {} import ffi
from datetime import datetime
from mpi4py import MPI
from geos_gtfv3 import geos_gtfv3

@ffi.def_extern()
def geos_gtfv3_interface_py(
    comm_c,
    npx, npy, npz,
    is_, ie, js, je,
    isd, ied, jsd, jed,
    bdt, nq_tot, ng, ptop, ks, layout_1, layout_2,
    adiabatic,
    hydrostatic, z_tracer, make_nh, fv_debug,
    reproduce_sum, do_sat_adj, do_vort_damp, rf_fast, fill,
    ntiles, ncnst, n_split, k_split, fv_sg_adj, n_sponge, n_zfilter, nwat,
    hord_tr, hord_tm, hord_dp, hord_mt, hord_vt,
    nord, kord_tm, kord_tr, kord_wz, kord_mt,
    d_ext, beta, vtdm4, ke_bg, d_con, d2_bg, d2_bg_k1, d2_bg_k2,
    p_fac, a_imp, dddmp, d4_bg, rf_cutoff, tau, consv_te,
    u, v, w, delz,
    pt, delp, q,
    ps, pe, pk, peln, pkz,
    phis, q_con, omga,
    ua, va, uc, vc,
    ak, bk,
    mfx, mfy, cx, cy, diss_est,
    dx, dy, dxa, dya, dxc, dyc, rdx, rdy, rdxa, rdya, rdxc, rdyc,
    cosa, cosa_s, sina_u, sina_v, cosa_u, cosa_v, rsin2, rsina, rsin_u, rsin_v,
    sin_sg, cos_sg,
    area, area_64, rarea, rarea_c, f0, fC, del6_u, del6_v, divg_u, divg_v,
    agrid, bgrid, a11, a12, a21, a22,
    edge_e, edge_w, edge_n, edge_s, nested, stretched_grid, da_min, da_min_c):

    # comm_c -> comm_py
    comm_py = MPI.Intracomm() # new comm, internal MPI_Comm handle is MPI_COMM_NULL
    comm_ptr = MPI._addressof(comm_py)  # internal MPI_Comm handle
    comm_ptr = ffi.cast('{}*', comm_ptr)  # make it a CFFI pointer
    comm_ptr[0] = comm_c  # assign comm_c to comm_py's MPI_Comm handle

    if comm_py.Get_rank() == 0:
        print('P:', datetime.now().isoformat(timespec='milliseconds'),
              '--in cffi interface', flush=True)

    geos_gtfv3(
        comm_py,
        npx, npy, npz,
        is_, ie, js, je,
        isd, ied, jsd, jed,
        bdt, nq_tot, ng, ptop, ks, layout_1, layout_2,
        adiabatic,
        hydrostatic, z_tracer, make_nh, fv_debug,
        reproduce_sum, do_sat_adj, do_vort_damp, rf_fast, fill,
        ntiles, ncnst, n_split, k_split, fv_sg_adj, n_sponge, n_zfilter, nwat,
        hord_tr, hord_tm, hord_dp, hord_mt, hord_vt,
        nord, kord_tm, kord_tr, kord_wz, kord_mt,
        d_ext, beta, vtdm4, ke_bg, d_con, d2_bg, d2_bg_k1, d2_bg_k2,
        p_fac, a_imp, dddmp, d4_bg, rf_cutoff, tau, consv_te,
        u, v, w, delz,
        pt, delp, q,
        ps, pe, pk, peln, pkz,
        phis, q_con, omga,
        ua, va, uc, vc,
        ak, bk,
        mfx, mfy, cx, cy, diss_est,
        dx, dy, dxa, dya, dxc, dyc, rdx, rdy, rdxa, rdya, rdxc, rdyc,
        cosa, cosa_s, sina_u, sina_v, cosa_u, cosa_v, rsin2, rsina, rsin_u, rsin_v,
        sin_sg, cos_sg,
        area, area_64, rarea, rarea_c, f0, fC, del6_u, del6_v, divg_u, divg_v,
        agrid, bgrid, a11, a12, a21, a22,
        edge_e, edge_w, edge_n, edge_s, nested, stretched_grid, da_min, da_min_c)
'''.format(TMPFILEBASE, _mpi_comm_t)

header = '''
extern void geos_gtfv3_interface_py(
    {} comm_c,
    int npx, int npy, int npz,
    int is_, int ie, int js, int je,
    int isd, int ied, int jsd, int jed,
    float bdt, int nq_tot, int ng, float ptop, int ks, int layout_1, int layout_2,
    int adiabatic,
    int hydrostatic, int z_tracer, int make_nh, int fv_debug,
    int reproduce_sum, int do_sat_adj, int do_vort_damp, int rf_fast, int fill,
    int ntiles, int ncnst, int n_split, int k_split,
    int fv_sg_adj, int n_sponge, int n_zfilter, int nwat,
    int hord_tr, int hord_tm, int hord_dp, int hord_mt, int hord_vt,
    int nord, int kord_tm, int kord_tr, int kord_wz, int kord_mt,
    float d_ext, float beta, float vtdm4, float ke_bg, float d_con,
    float d2_bg, float d2_bg_k1, float d2_bg_k2,
    float p_fac, float a_imp, float dddmp, float d4_bg, float rf_cutoff, float tau, float consv_te,
    // input/output
    float* u, float* v, float* w, float* delz,
    float* pt, float* delp, float* q,
    float* ps, float* pe, float* pk, float* peln, float* pkz,
    float* phis, float* q_con, float* omga, float* ua, float* va, float* uc, float* vc,
    // input
    const float* ak, const float* bk,
    // input/output
    float* mfx, float* mfy, float* cx, float* cy, float* diss_est,
    // input
    const float* dx, const float* dy,
    const float* dxa, const float* dya,
    const float* dxc, const float* dyc,
    const float* rdx, const float* rdy,
    const float* rdxa, const float* rdya,
    const float* rdxc, const float* rdyc,
    const float* cosa, const float* cosa_s, const float* sina_u, const float* sina_v,
    const float* cosa_u, const float* cosa_v,
    const float* rsin2, const float* rsina, const float* rsin_u, const float* rsin_v,
    const float* sin_sg, const float* cos_sg,
    const float* area, const double* area_64, const float* rarea, const float* rarea_c,
    const float* f0, const float* fC,
    const float* del6_u, const float* del6_v, const float* divg_u, const float* divg_v,
    const float* agrid, const float* bgrid,
    const float* a11, const float* a12, const float* a21, const float* a22,
    const double* edge_e, const double* edge_w, const double* edge_n, const double* edge_s,
    int nested, int stretched_grid, double da_min, double da_min_c);
'''.format(_mpi_comm_t)

with open(TMPFILEBASE+'.h', 'w') as f:
    f.write(header)
ffi.embedding_api(header)

source_header = r'''#include "{}.h"'''.format(TMPFILEBASE)
ffi.set_source(TMPFILEBASE, source_header)

ffi.embedding_init_code(source)
ffi.compile(target='lib'+TMPFILEBASE+'.so', verbose=True)
