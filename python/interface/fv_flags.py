import dataclasses
from typing import Union
from pyfv3._config import (
    DynamicalCoreConfig,
    RiemannConfig,
    SatAdjustConfig,
    RemappingConfig,
    DGridShallowWaterLagrangianDynamicsConfig,
    AcousticDynamicsConfig,
)


@dataclasses.dataclass
class FVFlags:
    # Fortran flagstruct
    grid_type: int
    hord_mt: int
    kord_mt: int
    kord_wz: int
    hord_vt: int
    hord_tm: int
    hord_dp: int
    kord_tm: int
    hord_tr: int
    kord_tr: int
    scale_z: float
    w_max: float
    z_min: float
    lim_fac: float
    nord: int
    nord_tr: int
    dddmp: float
    d2_bg: float
    d4_bg: float
    vtdm4: float
    trdm2: float
    d2_bg_k1: float
    d2_bg_k2: float
    d2_divg_max_k1: float
    d2_divg_max_k2: float
    damp_k_k1: float
    damp_k_k2: float
    n_zs_filter: int
    nord_zs_filter: int
    full_zs_filter: int # [bool] but under ifx/icc those are int8_t
    RF_fast: int # [bool] but under ifx/icc those are int8_t
    Beljaars_TOFD: int # [bool] but under ifx/icc those are int8_t
    consv_am: int # [bool] but under ifx/icc those are int8_t
    do_sat_adj: int # [bool] but under ifx/icc those are int8_t
    do_f3d: int # [bool] but under ifx/icc those are int8_t
    no_dycore: int # [bool] but under ifx/icc those are int8_t
    convert_ke: int # [bool] but under ifx/icc those are int8_t
    do_vort_damp: int # [bool] but under ifx/icc those are int8_t
    use_old_omega: int # [bool] but under ifx/icc those are int8_t
    beta: float
    n_zfilter: int
    n_sponge: int
    d_ext: float
    nwat: int
    warm_start: int # [bool] but under ifx/icc those are int8_t
    inline_q: int # [bool] but under ifx/icc those are int8_t
    adiabatic: int # [bool] but under ifx/icc those are int8_t
    shift_fac: float
    do_schmidt: int # [bool] but under ifx/icc those are int8_t
    stretch_fac: float
    target_lat: float
    target_lon: float
    reset_eta: int # [bool] but under ifx/icc those are int8_t
    p_fac: float
    a_imp: float
    dz_min: float
    n_split: int
    m_split: int
    k_split: int
    use_logp: int # [bool] but under ifx/icc those are int8_t
    q_split: int
    print_freq: int
    write_3d_diags: int # [bool] but under ifx/icc those are int8_t
    npx: int
    npy: int
    npz: int
    npz_rst: int
    ncnst: int
    pnats: int
    dnats: int
    ntiles: int
    ndims: int
    nf_omega: int
    fv_sg_adj: int
    na_init: int
    nudge_dz: int # [bool] but under ifx/icc those are int8_t
    p_ref: float
    dry_mass: float
    nt_prog: int
    nt_phys: int
    tau_h2o: float
    delt_max: float
    d_con: float
    ke_bg: float
    consv_te: float
    tau: float
    rf_cutoff: float
    filter_phys: int # [bool] but under ifx/icc those are int8_t
    dwind_2d: int # [bool] but under ifx/icc those are int8_t
    breed_vortex_inline: int # [bool] but under ifx/icc those are int8_t
    range_warn: int # [bool] but under ifx/icc those are int8_t
    fill: int # [bool] but under ifx/icc those are int8_t
    fill_dp: int # [bool] but under ifx/icc those are int8_t
    fill_wz: int # [bool] but under ifx/icc those are int8_t
    check_negative: int # [bool] but under ifx/icc those are int8_t
    non_ortho: int # [bool] but under ifx/icc those are int8_t
    moist_phys: int # [bool] but under ifx/icc those are int8_t
    do_Held_Suarez: int # [bool] but under ifx/icc those are int8_t
    do_reed_physics: int # [bool] but under ifx/icc those are int8_t
    reed_cond_only: int # [bool] but under ifx/icc those are int8_t
    reproduce_sum: int # [bool] but under ifx/icc those are int8_t
    adjust_dry_mass: int # [bool] but under ifx/icc those are int8_t
    fv_debug: int # [bool] but under ifx/icc those are int8_t
    srf_init: int # [bool] but under ifx/icc those are int8_t
    mountain: int # [bool] but under ifx/icc those are int8_t
    old_divg_damp: int # [bool] but under ifx/icc those are int8_t
    remap_option: int
    gmao_remap: int
    z_tracer: int # [bool] but under ifx/icc those are int8_t
    fv_land: int # [bool] but under ifx/icc those are int8_t
    nudge: int # [bool] but under ifx/icc those are int8_t
    nudge_ic: int # [bool] but under ifx/icc those are int8_t
    ncep_ic: int # [bool] but under ifx/icc those are int8_t
    nggps_ic: int # [bool] but under ifx/icc those are int8_t
    ecmwf_ic: int # [bool] but under ifx/icc those are int8_t
    gfs_phil: int # [bool] but under ifx/icc those are int8_t
    agrid_vel_rst: int # [bool] but under ifx/icc those are int8_t
    use_new_ncep: int # [bool] but under ifx/icc those are int8_t
    use_ncep_phy: int # [bool] but under ifx/icc those are int8_t
    fv_diag_ic: int # [bool] but under ifx/icc those are int8_t
    external_ic: int # [bool] but under ifx/icc those are int8_t
    external_eta: int # [bool] but under ifx/icc those are int8_t
    read_increment: int # [bool] but under ifx/icc those are int8_t
    do_skeb: int # [bool] but under ifx/icc those are int8_t
    skeb_npass: int
    hydrostatic: int # [bool] but under ifx/icc those are int8_t
    phys_hydrostatic: int # [bool] but under ifx/icc those are int8_t
    use_hydro_pressure: int # [bool] but under ifx/icc those are int8_t
    do_uni_zfull: int # [bool] but under ifx/icc those are int8_t
    hybrid_z: int # [bool] but under ifx/icc those are int8_t
    Make_NH: int # [bool] but under ifx/icc those are int8_t
    make_hybrid_z: int # [bool] but under ifx/icc those are int8_t
    nudge_qv: int # [bool] but under ifx/icc those are int8_t
    add_noise: float
    a2b_ord: int
    c2l_ord: int
    dx_const: float
    dy_const: float
    deglat: float
    deglon_start: float
    adj_mass_vmr: int # [bool] but under ifx/icc those are int8_t
    compute_coords_locally: int # [bool] but under ifx/icc those are int8_t
    # Grid
    layout_x: int
    layout_y: int
    # Magic number needs to be last item
    mn_123456789: int


def _generic_config_bridge(
    py_config: Union[
        DynamicalCoreConfig,
        RiemannConfig,
        SatAdjustConfig,
        RemappingConfig,
        DGridShallowWaterLagrangianDynamicsConfig,
        AcousticDynamicsConfig,
    ],
    fv_config: FVFlags,
):
    keys = list(filter(lambda k: not k.startswith("__"), dir(type(py_config))))
    for k in keys:
        if hasattr(fv_config, k):
            v_fortran = getattr(fv_config, k)
            v_py = getattr(py_config, k)
            if isinstance(v_py, bool):
                setattr(py_config, k, v_fortran != 0)
            else:
                setattr(py_config, k, v_fortran)


def FVFlags_to_DycoreConfig(
    fv_config: FVFlags,
    py_config: DynamicalCoreConfig,
):
    if fv_config.mn_123456789 != 123456789:
        raise RuntimeError(
            "Magic number failed, pyFV3 interface is broken on the python side"
        )

    _generic_config_bridge(py_config, fv_config)
    py_config.layout = (
        getattr(fv_config, "layout_x"),
        getattr(fv_config, "layout_y"),
    )
