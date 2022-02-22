from datetime import datetime
from f_py_conversion import fortran_input_data_to_numpy
from f_py_conversion import numpy_output_data_to_fortran
import geos_gtfv3_initialize
from geos_gtfv3_run import geos_gtfv3_run
from geos_gtfv3_debug import write_sum_of_untranslated_arrays

def geos_gtfv3(
        comm,
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
        area, rarea, rarea_c, f0, fC, del6_u, del6_v, divg_u, divg_v,
        agrid, bgrid, a11, a12, a21, a22,
        edge_e, edge_w, edge_n, edge_s, nested, stretched_grid, da_min, da_min_c):

    BACKEND = 'gtx86'

    rank = comm.Get_rank()

    if rank == 0:
        print('P:', datetime.now().isoformat(timespec='milliseconds'),
              '--in top level geos-gtfv3, backend:', BACKEND, flush=True)

    # Initialize namelist - set global var spec
    # This function gets exercised only the first timestep
    geos_gtfv3_initialize.initialize_namelist(
        npx, npy, npz, layout_1, layout_2,
        adiabatic, hydrostatic, z_tracer,
        do_sat_adj, do_vort_damp, rf_fast, fill,
        ntiles, n_split, k_split, fv_sg_adj, n_sponge, n_zfilter, nwat,
        hord_tr, hord_tm, hord_dp, hord_mt, hord_vt,
        nord, kord_tm, kord_tr, kord_wz, kord_mt,
        d_ext, beta, vtdm4, ke_bg, d_con, d2_bg, d2_bg_k1, d2_bg_k2,
        p_fac, a_imp, dddmp, d4_bg, rf_cutoff, tau)
    assert geos_gtfv3_initialize.spec is not None

    if rank == 0:
        print('P:', datetime.now().isoformat(timespec='milliseconds'),
              '--created namelist', flush=True)

    # Convert Fortran arrays to NumPy
    # This function gets exercised at every timestep
    fv3_input_data = fortran_input_data_to_numpy(
        npx, npy, npz,
        is_, ie, js, je, isd, ied, jsd, jed,
        ncnst, consv_te, bdt, ptop, n_split, ks,
        ak, bk,
        # input/output arrays
        u, v, w, delz,
        pt, delp, q,
        ps, pe, pk, peln, pkz,
        phis, q_con, omga,
        ua, va, uc, vc,
        mfx, mfy, cx, cy, diss_est)
    # write_sum_of_untranslated_arrays(comm, fv3_input_data)

    if rank == 0:
        print('P:', datetime.now().isoformat(timespec='milliseconds'),
              '--created input data', flush=True)

    # Initialize dycore - set global vars driver_object and dycore
    # This function gets exercised only the first timestep
    geos_gtfv3_initialize.initialize_dycore(
        BACKEND, comm,
        npx, npy, npz,
        is_, ie, js, je, isd, ied, jsd, jed,
        # grid data
        dx, dy, dxa, dya, dxc, dyc,
        rdx, rdy, rdxa, rdya, rdxc, rdyc,
        cosa, cosa_s, sina_u, sina_v,
        cosa_u, cosa_v, rsin2, rsina, rsin_u, rsin_v,
        sin_sg, cos_sg,
        area, rarea, rarea_c, f0, fC,
        del6_u, del6_v, divg_u, divg_v,
        agrid, bgrid, a11, a12, a21, a22,
        edge_e, edge_w, edge_n, edge_s,
        nested, stretched_grid, da_min, da_min_c,
        # input data
        fv3_input_data)
    assert geos_gtfv3_initialize.driver_object is not None, 'driver_object is None'
    assert geos_gtfv3_initialize.dycore is not None, 'dycore is None'

    if rank == 0:
        print('P:', datetime.now().isoformat(timespec='milliseconds'),
              '--initialized dycore', flush=True)

    # Run gtFV3
    gtfv3_output_data = geos_gtfv3_run(
        comm,
        geos_gtfv3_initialize.spec,
        geos_gtfv3_initialize.driver_object,
        geos_gtfv3_initialize.dycore,
        fv3_input_data)
    # write_sum_of_untranslated_arrays(comm, gtfv3_output_data)

    # Convert NumPy arrays back to Fortran
    numpy_output_data_to_fortran(
        gtfv3_output_data,
        u, v, w, delz,
        pt, delp, q,
        ps, pe, pk, peln, pkz,
        phis, q_con, omga,
        ua, va, uc, vc,
        mfx, mfy, cx, cy, diss_est)

    if rank == 0:
       print('P:', datetime.now().isoformat(timespec='milliseconds'),
             '--converted NumPy arrays to Fortran', flush=True)
