import dataclasses
from typing import Union
from pyFV3._config import (
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
    full_zs_filter: bool
    RF_fast: bool
    Beljaars_TOFD: bool
    consv_am: bool
    do_sat_adj: bool
    do_f3d: bool
    no_dycore: bool
    convert_ke: bool
    do_vort_damp: bool
    use_old_omega: bool
    beta: float
    n_zfilter: int
    n_sponge: int
    d_ext: float
    nwat: int
    warm_start: bool
    inline_q: bool
    adiabatic: bool
    shift_fac: float
    do_schmidt: bool
    stretch_fac: float
    target_lat: float
    target_lon: float
    reset_eta: bool
    p_fac: float
    a_imp: float
    dz_min: float
    n_split: int
    m_split: int
    k_split: int
    use_logp: bool
    q_split: int
    print_freq: int
    write_3d_diags: bool
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
    nudge_dz: bool
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
    filter_phys: bool
    dwind_2d: bool
    breed_vortex_inline: bool
    range_warn: bool
    fill: bool
    fill_dp: bool
    fill_wz: bool
    check_negative: bool
    non_ortho: bool
    moist_phys: bool
    do_Held_Suarez: bool
    do_reed_physics: bool
    reed_cond_only: bool
    reproduce_sum: bool
    adjust_dry_mass: bool
    fv_debug: bool
    srf_init: bool
    mountain: bool
    old_divg_damp: bool
    remap_option: int
    gmao_remap: int
    z_tracer: bool
    fv_land: bool
    nudge: bool
    nudge_ic: bool
    ncep_ic: bool
    nggps_ic: bool
    ecmwf_ic: bool
    gfs_phil: bool
    agrid_vel_rst: bool
    use_new_ncep: bool
    use_ncep_phy: bool
    fv_diag_ic: bool
    external_ic: bool
    external_eta: bool
    read_increment: bool
    do_skeb: bool
    skeb_npass: int
    hydrostatic: bool
    phys_hydrostatic: bool
    use_hydro_pressure: bool
    do_uni_zfull: bool
    hybrid_z: bool
    Make_NH: bool
    make_hybrid_z: bool
    nudge_qv: bool
    add_noise: float
    a2b_ord: int
    c2l_ord: int
    dx_const: float
    dy_const: float
    deglat: float
    deglon_start: float
    adj_mass_vmr: bool
    compute_coords_locally: bool
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
            setattr(py_config, k, getattr(fv_config, k))


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
