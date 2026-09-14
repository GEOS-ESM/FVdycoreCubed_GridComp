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
    full_zs_filter: bool | int
    RF_fast: bool | int
    Beljaars_TOFD: bool | int
    consv_am: bool | int
    do_sat_adj: bool | int
    do_f3d: bool | int
    no_dycore: bool | int
    convert_ke: bool | int
    do_vort_damp: bool | int
    use_old_omega: bool | int
    beta: float
    n_zfilter: int
    n_sponge: int
    d_ext: float
    nwat: int
    warm_start: bool | int
    inline_q: bool | int
    adiabatic: bool | int
    shift_fac: float
    do_schmidt: bool | int
    stretch_fac: float
    target_lat: float
    target_lon: float
    reset_eta: bool | int
    p_fac: float
    a_imp: float
    dz_min: float
    n_split: int
    m_split: int
    k_split: int
    use_logp: bool | int
    q_split: int
    print_freq: int
    write_3d_diags: bool | int
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
    nudge_dz: bool | int
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
    filter_phys: bool | int
    dwind_2d: bool | int
    breed_vortex_inline: bool | int
    range_warn: bool | int
    fill: bool | int
    fill_dp: bool | int
    fill_wz: bool | int
    check_negative: bool | int
    non_ortho: bool | int
    moist_phys: bool | int
    do_Held_Suarez: bool | int
    do_reed_physics: bool | int
    reed_cond_only: bool | int
    reproduce_sum: bool | int
    adjust_dry_mass: bool | int
    fv_debug: bool | int
    srf_init: bool | int
    mountain: bool | int
    old_divg_damp: bool | int
    remap_option: int
    gmao_remap: int
    z_tracer: bool | int
    fv_land: bool | int
    nudge: bool | int
    nudge_ic: bool | int
    ncep_ic: bool | int
    nggps_ic: bool | int
    ecmwf_ic: bool | int
    gfs_phil: bool | int
    agrid_vel_rst: bool | int
    use_new_ncep: bool | int
    use_ncep_phy: bool | int
    fv_diag_ic: bool | int
    external_ic: bool | int
    external_eta: bool | int
    read_increment: bool | int
    do_skeb: bool | int
    skeb_npass: int
    hydrostatic: bool | int
    phys_hydrostatic: bool | int
    use_hydro_pressure: bool | int
    do_uni_zfull: bool | int
    hybrid_z: bool | int
    Make_NH: bool | int
    make_hybrid_z: bool | int
    nudge_qv: bool | int
    add_noise: float
    a2b_ord: int
    c2l_ord: int
    dx_const: float
    dy_const: float
    deglat: float
    deglon_start: float
    adj_mass_vmr: bool | int
    compute_coords_locally: bool | int
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
