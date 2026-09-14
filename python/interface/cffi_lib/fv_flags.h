#pragma once

/***
 * Dynamical core configuration from GEOS (namelist and extra flags)
 ***/

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
    unsigned char full_zs_filter;
    unsigned char rf_fast;
    unsigned char Beljaars_TOFD;
    unsigned char consv_am;
    unsigned char do_sat_adj;
    unsigned char do_f3d;
    unsigned char no_dycore;
    unsigned char convert_ke;
    unsigned char do_vort_damp;
    unsigned char use_old_omega;
    float beta;
    int n_zfilter;
    int n_sponge;
    float d_ext;
    int nwat;
    unsigned char warm_start;
    unsigned char inline_q;
    unsigned char adiabatic;
    float shift_fac;
    unsigned char do_schmidt;
    float stretch_fac; // Fortran original type is real(kind=R_GRID) we maximize comp by using double.
    float target_lat;  // Fortran original type is real(kind=R_GRID) we maximize comp by using double.
    float target_lon;  // Fortran original type is real(kind=R_GRID) we maximize comp by using double.
    unsigned char reset_eta;
    float p_fac;
    float a_imp;
    float dz_min;
    int n_split;
    int m_split;
    int k_split;
    unsigned char use_logp;
    int q_split;
    int print_freq;
    unsigned char write_3d_diags;
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
    unsigned char nudge_dz;
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
    unsigned char filter_phys;
    unsigned char dwind_2d;
    unsigned char breed_vortex_inline;
    unsigned char range_warn;
    unsigned char fill;
    unsigned char fill_dp;
    unsigned char fill_wz;
    unsigned char check_negative;
    unsigned char non_ortho;
    unsigned char moist_phys;
    unsigned char do_Held_Suarez;
    unsigned char do_reed_physics;
    unsigned char reed_cond_only;
    unsigned char reproduce_sum;
    unsigned char adjust_dry_mass;
    unsigned char fv_debug;
    unsigned char srf_init;
    unsigned char mountain;
    unsigned char old_divg_damp;
    int remap_option;
    int gmao_remap;
    unsigned char z_tracer;
    unsigned char fv_land;
    unsigned char nudge;
    unsigned char nudge_ic;
    unsigned char ncep_ic;
    unsigned char nggps_ic;
    unsigned char ecmwf_ic;
    unsigned char gfs_phil;
    unsigned char agrid_vel_rst;
    unsigned char use_new_ncep;
    unsigned char use_ncep_phy;
    unsigned char fv_diag_ic;
    unsigned char external_ic;
    unsigned char external_eta;
    unsigned char read_increment;
    unsigned char do_skeb;
    int skeb_npass;
    unsigned char hydrostatic;
    unsigned char phys_hydrostatic;
    unsigned char use_hydro_pressure;
    unsigned char do_uni_zfull;
    unsigned char hybrid_z;
    unsigned char Make_NH;
    unsigned char make_hybrid_z;
    unsigned char nudge_qv;
    float add_noise;
    int a2b_ord;
    int c2l_ord;
    float dx_const;
    float dy_const;
    float deglat;
    double deglon_start;
    unsigned char adj_mass_vmr;
    unsigned char compute_coords_locally;
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

extern int pyfv3_interface_py_init(
    fv_flags_t *fv_flags,
    void *comm_c,
    int npx, int npy, int npz, int ntiles,
    int is_, int ie, int js, int je, int isd, int ied, int jsd, int jed,
    float bdt, int nq_tot,
    const float *ak, const float *bk, const float *phis);

extern int pyfv3_interface_py_run(
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

extern int pyfv3_interface_py_finalize();
