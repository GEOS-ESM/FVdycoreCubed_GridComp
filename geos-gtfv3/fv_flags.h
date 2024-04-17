#pragma once

/***
 * Dynamical core configuration from GEOS (namelist and extra flags)
 ***/

#include <stdbool.h>
#include <stdlib.h>

// Fortran FlagStruct
typedef struct
{
    int grid_type;
    int hord_mt;
    int kord_mt;
    int kord_wz;
    int hord_vt;
    int hord_tm;
    int hord_dp;
    int kord_tm;
    int hord_tr;
    int kord_tr;
    float scale_z;
    float w_max;
    float z_min;
    float lim_fac;
    int nord;
    int nord_tr;
    float dddmp;
    float d2_bg;
    float d4_bg;
    float vtdm4;
    float trdm2;
    float d2_bg_k1;
    float d2_bg_k2;
    float d2_divg_max_k1;
    float d2_divg_max_k2;
    float damp_k_k1;
    float damp_k_k2;
    int n_zs_filter;
    int nord_zs_filter;
    bool full_zs_filter;
    bool RF_fast;
    bool Beljaars_TOFD;
    bool consv_am;
    bool do_sat_adj;
    bool do_f3d;
    bool no_dycore;
    bool convert_ke;
    bool do_vort_damp;
    bool use_old_omega;
    float beta;
    int n_zfilter;
    int n_sponge;
    float d_ext;
    int nwat;
    bool warm_start;
    bool inline_q;
    bool adiabatic;
    float shift_fac;
    bool do_schmidt;
    float stretch_fac; // Fortran original type is real(kind=R_GRID) we maximize comp by using double.
    float target_lat;  // Fortran original type is real(kind=R_GRID) we maximize comp by using double.
    float target_lon;  // Fortran original type is real(kind=R_GRID) we maximize comp by using double.
    bool reset_eta;
    float p_fac;
    float a_imp;
    float dz_min;
    int n_split;
    int m_split;
    int k_split;
    bool use_logp;
    int q_split;
    int print_freq;
    bool write_3d_diags;
    int npx;
    int npy;
    int npz;
    int npz_rst;
    int ncnst;
    int pnats;
    int dnats;
    int ntiles;
    int ndims;
    int nf_omega;
    int fv_sg_adj;
    int na_init;
    bool nudge_dz;
    float p_ref;
    float dry_mass;
    int nt_prog;
    int nt_phys;
    float tau_h2o;
    float delt_max;
    float d_con;
    float ke_bg;
    float consv_te;
    float tau;
    float rf_cutoff;
    bool filter_phys;
    bool dwind_2d;
    bool breed_vortex_inline;
    bool range_warn;
    bool fill;
    bool fill_dp;
    bool fill_wz;
    bool check_negative;
    bool non_ortho;
    bool moist_phys;
    bool do_Held_Suarez;
    bool do_reed_physics;
    bool reed_cond_only;
    bool reproduce_sum;
    bool adjust_dry_mass;
    bool fv_debug;
    bool srf_init;
    bool mountain;
    bool old_divg_damp;
    int remap_option;
    int gmao_remap;
    bool z_tracer;
    bool fv_land;
    bool nudge;
    bool nudge_ic;
    bool ncep_ic;
    bool nggps_ic;
    bool ecmwf_ic;
    bool gfs_phil;
    bool agrid_vel_rst;
    bool use_new_ncep;
    bool use_ncep_phy;
    bool fv_diag_ic;
    bool external_ic;
    bool external_eta;
    bool read_increment;
    bool do_skeb;
    int skeb_npass;
    bool hydrostatic;
    bool phys_hydrostatic;
    bool use_hydro_pressure;
    bool do_uni_zfull;
    bool hybrid_z;
    bool Make_NH;
    bool make_hybrid_z;
    bool nudge_qv;
    float add_noise;
    int a2b_ord;
    int c2l_ord;
    float dx_const;
    float dy_const;
    float deglat;
    double deglon_start;
    bool adj_mass_vmr;
    bool compute_coords_locally;
    // Grid information
    int layout_x;
    int layout_y;
    // Magic number needs to be last item
    int mn_123456789;
} fv_flags_t;

typedef union
{
    int comm_int;
    void *comm_ptr;
} MPI_Comm_t;

extern int geos_gtfv3_interface_py_init(
    fv_flags_t *fv_flags,
    void *comm_c,
    int npx, int npy, int npz, int ntiles,
    int is_, int ie, int js, int je, int isd, int ied, int jsd, int jed,
    float bdt, int nq_tot,
    const float *ak, const float *bk);

extern int geos_gtfv3_interface_py(
    void *comm_c,
    int npx, int npy, int npz, int ntiles,
    int is_, int ie, int js, int je, int isd, int ied, int jsd, int jed,
    float bdt, int nq_tot, int ng, float ptop, int ks, int layout_1, int layout_2,
    int adiabatic,
    // input/output
    float *u, float *v, float *w, float *delz,
    float *pt, float *delp, float *q,
    float *ps, float *pe, float *pk, float *peln, float *pkz,
    float *phis, float *q_con, float *omga, float *ua, float *va, float *uc, float *vc,
    // input/output
    float *mfx, float *mfy, float *cx, float *cy, float *diss_est);

extern int geos_gtfv3_interface_py_finalize();
