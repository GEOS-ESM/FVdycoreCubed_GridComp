from datetime import datetime

from f_py_conversion import fortran_grid_data_to_numpy

import gt4py
import fv3core
import fv3core._config as spec
import fv3core.testing
import fv3gfs.util as util


__initialized_dycore = False
__initialized_namelist = False
driver_object = None
dycore = None


def initialize_namelist(
        npx, npy, npz, layout_1, layout_2,
        adiabatic_int, hydrostatic_int, z_tracer_int,
        do_sat_adj_int, do_vort_damp_int, rf_fast_int, fill_int,
        n_split, k_split, fv_sg_adj, n_sponge, n_zfilter, nwat,
        hord_tr, hord_tm, hord_dp, hord_mt, hord_vt,
        nord, kord_tm, kord_tr, kord_wz, kord_mt,
        d_ext, beta, vtdm4, ke_bg, d_con, d2_bg, d2_bg_k1, d2_bg_k2,
        p_fac, a_imp, dddmp, d4_bg, rf_cutoff, tau):
    '''
    '''
    global __initialized_namelist

    if not __initialized_namelist:
        spec.set_namelist_from_dict({
            'npx': npx, 'npy': npy, 'npz': npz,
            'layout': [layout_1, layout_2],

            'adiabatic': False if adiabatic_int == 0 else True,
            'hydrostatic': False if hydrostatic_int == 0 else True,
            'z_tracer': False if z_tracer_int == 0 else True,
            'do_sat_adj': False if do_sat_adj_int == 0 else True,
            'do_vort_damp': False if do_vort_damp_int == 0 else True,
            'rf_fast': False if rf_fast_int == 0 else True,
            'fill': False if fill_int == 0 else True,

            'n_split': n_split, 'k_split': k_split, 'fv_sg_adj': fv_sg_adj,
            'n_sponge': n_sponge, 'n_zfilter': n_zfilter, 'nwat': nwat,
            'hord_tr': hord_tr, 'hord_tm': hord_tm, 'hord_dp': hord_dp,
            'hord_mt': hord_mt, 'hord_vt': hord_vt,
            'nord': nord,
            'kord_tm': kord_tm, 'kord_tr': kord_tr,
            'kord_wz': kord_wz, 'kord_mt': kord_mt,

            'd_ext': d_ext, 'beta': beta, 'vtdm4': vtdm4,
            'ke_bg': ke_bg, 'd_con': d_con,
            'd2_bg': d2_bg, 'd2_bg_k1': d2_bg_k1, 'd2_bg_k2': d2_bg_k2,
            'p_fac': p_fac, 'a_imp': a_imp, 'dddmp': dddmp, 'd4_bg': d4_bg,
            'rf_cutoff': rf_cutoff, 'tau': tau,
            'delt_max': 0.002, 'ntiles': 6}) # namelist
        __initialized_namelist = True


def initialize_dycore(
        backend, comm,
        npx, npy, npz,
        is_, ie, js, je, isd, ied, jsd, jed,
        # grid data
        dx_ptr, dy_ptr, dxa_ptr, dya_ptr, dxc_ptr, dyc_ptr,
        rdx_ptr, rdy_ptr, rdxa_ptr, rdya_ptr, rdxc_ptr, rdyc_ptr,
        cosa_ptr, cosa_s_ptr, sina_u_ptr, sina_v_ptr,
        cosa_u_ptr, cosa_v_ptr, rsin2_ptr, rsina_ptr, rsin_u_ptr, rsin_v_ptr,
        sin_sg_ptr, cos_sg_ptr,
        area_ptr, rarea_ptr, rarea_c_ptr, f0_ptr, fC_ptr,
        del6_u_ptr, del6_v_ptr, divg_u_ptr, divg_v_ptr,
        agrid_ptr, bgrid_ptr, a11_ptr, a12_ptr, a21_ptr, a22_ptr,
        edge_e_ptr, edge_w_ptr, edge_n_ptr, edge_s_ptr,
        nested_int, stretched_grid_int, da_min, da_min_c,
        # input data
        fv3_input_data):
    '''
    '''
    global __initialized_dycore, driver_object, dycore

    if not __initialized_namelist:
        raise ValueError('namelist has not been initialized')

    if not __initialized_dycore:

        fv3core.set_backend(backend)
        fv3core.set_rebuild(False)
        fv3core.set_validate_args(True)
        # fv3core.utils.global_config.set_do_halo_exchange(True)

        # Convert Fortran data to NumPy
        grid_data = fortran_grid_data_to_numpy(
            npx, npy, npz,
            is_, ie, js, je, isd, ied, jsd, jed,
            nested_int, stretched_grid_int, da_min, da_min_c,
            # input arrays - grid data
            dx_ptr, dy_ptr, dxa_ptr, dya_ptr, dxc_ptr, dyc_ptr,
            rdx_ptr, rdy_ptr, rdxa_ptr, rdya_ptr, rdxc_ptr, rdyc_ptr,
            cosa_ptr, cosa_s_ptr, sina_u_ptr, sina_v_ptr,
            cosa_u_ptr, cosa_v_ptr, rsin2_ptr, rsina_ptr, rsin_u_ptr, rsin_v_ptr,
            sin_sg_ptr, cos_sg_ptr,
            area_ptr, rarea_ptr, rarea_c_ptr, f0_ptr, fC_ptr,
            del6_u_ptr, del6_v_ptr, divg_u_ptr, divg_v_ptr,
            agrid_ptr, bgrid_ptr, a11_ptr, a12_ptr, a21_ptr, a22_ptr,
            edge_e_ptr, edge_w_ptr, edge_n_ptr, edge_s_ptr)

        # Create grid
        grid = fv3core.testing.TranslateGrid(grid_data, comm.Get_rank()).python_grid()
        spec.set_grid(grid)

        # Create driver_object
        driver_object = fv3core.testing.TranslateFVDynamics([grid])

        # Create state
        state = driver_object.state_from_inputs(fv3_input_data)

        # Communicator
        layout = spec.namelist.layout
        partitioner = util.CubedSpherePartitioner(util.TilePartitioner(layout))
        communicator = util.CubedSphereCommunicator(comm, partitioner)

        # Instantiate DynamicalCore
        dycore = fv3core.DynamicalCore(
            communicator,
            grid_data=spec.grid.grid_data,
            grid_indexing=spec.grid.grid_indexing,
            damping_coefficients=spec.grid.damping_coefficients,
            config=spec.namelist.dynamical_core,
            ak=state["atmosphere_hybrid_a_coordinate"],
            bk=state["atmosphere_hybrid_b_coordinate"],
            phis=state["surface_geopotential"],)

        # Initialization complete
        __initialized_dycore = True
