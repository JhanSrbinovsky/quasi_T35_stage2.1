MODULE allocate_jules_arrays_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ALLOCATE_JULES_ARRAYS_MOD'

CONTAINS
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine ALLOCATE_JULES_ARRAYS
!
! Description: Routine that allocates memory to the JULES arrays
! This assume that the values in the jules_surface_types module have been set
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.
!
!   Code Owner: See Unified Model Code Owner's HTML page
!   This file belongs in section: Land

#if defined(UM_JULES)
SUBROUTINE allocate_jules_arrays(                                             &
  land_pts,nsurft,sm_levels,nice,nice_use                                     &
  )
#else
SUBROUTINE allocate_jules_arrays()

!Replacements for the argument list
  USE ancil_info,               ONLY:                                         &
    land_pts, nsurft
  USE jules_soil_mod,           ONLY:                                         &
    sm_levels
  USE jules_rivers_mod,         ONLY:                                         &
     tot_surf_runoff_gb, tot_sub_runoff_gb, acc_lake_evap_gb
  USE jules_sea_seaice_mod,     ONLY:                                         &
    nice, nice_use
#endif

!Common Non-science modules
  USE parkind1,                 ONLY:                                         &
    jprb, jpim
  USE yomhook,                  ONLY:                                         &
    lhook, dr_hook
  USE jules_print_mgr,          ONLY:                                         &
    jules_message, jules_print, PrNorm
  USE ereport_mod,              ONLY:                                         &
    ereport

!USE statements that apply to both UM and JULES
  USE ancil_info,               ONLY:                                         &
    ssi_index, sea_index, sice_index, sice_pts_ncat, fssi_ij,                 &
    sea_frac, sice_frac, sice_frac_ncat, sice_index_ncat,                     &
    l_lice_point, l_soil_point

  USE bvoc_vars,                ONLY:                                         &
    isoprene_gb, isoprene_pft, terpene_gb , terpene_pft,                      &
    methanol_gb, methanol_pft, acetone_gb, acetone_pft

  USE c_elevate,                ONLY:                                         &
    surf_hgt_surft, lw_down_elevcorr_surft

  USE crop_vars_mod,            ONLY:                                         &
    sthu_irr_gb, frac_irr_all, frac_irr_gb, frac_irr_old_gb, frac_irr_surft,  &
    plant_n_gb, smc_irr_gb, wt_ext_irr_surft, gs_irr_surft, dvimax_gb,        &
    gc_irr_surft, resfs_irr_surft, ext_irr_gb, wt_ext_irr_gb, fsmc_irr_gb,    &
    irrdaysdiag_gb, icntmax_gb, tl_1_day_av_gb, tl_1_day_av_use_gb,           &
    rn_1_day_av_gb, rn_1_day_av_use_gb, ndpy, nyav, prec_1_day_av_gb,         &
    prec_1_day_av_use_gb, irrig_water_gb, irrtiles

  USE c_z0h_z0m,                ONLY:                                         &
    z0h_z0m, z0h_z0m_classic

  USE fire_vars,                ONLY:                                         &
    burnt_area, burnt_area_ft,                                                &
    emitted_carbon, emitted_carbon_ft, emitted_carbon_DPM, emitted_carbon_RPM,&
    fire_em_CO2, fire_em_CO2_ft, fire_em_CO2_DPM, fire_em_CO2_RPM,            &
    fire_em_CO, fire_em_CO_ft, fire_em_CO_DPM, fire_em_CO_RPM,                &
    fire_em_CH4, fire_em_CH4_ft, fire_em_CH4_DPM, fire_em_CH4_RPM,            &
    fire_em_NOx, fire_em_NOx_ft, fire_em_NOx_DPM, fire_em_NOx_RPM,            &
    fire_em_SO2, fire_em_SO2_ft, fire_em_SO2_DPM, fire_em_SO2_RPM,            &
    fire_em_OC,  fire_em_OC_ft, fire_em_OC_DPM,  fire_em_OC_RPM,              &
    fire_em_BC, fire_em_BC_ft, fire_em_BC_DPM, fire_em_BC_RPM,                &
    pop_den, flash_rate

  USE fluxes,                   ONLY:                                         &
    anthrop_heat_surft, surf_ht_store_surft,                                  &
    sw_sicat, sw_rts_sicat, swup_rts_sicat, swdn_rts_sicat, alb_sicat,        &
    et_stom_gb, et_stom_pft

  !These should be spirited away to other modules, allowing it to be retired
  USE jules_mod,                ONLY:                                         &
    snowdep_surft, albobs_scaling_surft

  USE jules_radiation_mod,      ONLY:                                         &
    l_albedo_obs, l_spec_albedo

  USE jules_snow_mod,           ONLY:                                         &
    nsmax, cansnowtile

  USE jules_internal,           ONLY:                                         &
    unload_backgrnd_pft

  USE jules_soil_mod,           ONLY:                                         &
    l_bedrock, ns_deep

  USE jules_surface_mod,        ONLY:                                         &
    diff_frac

  USE jules_surface_types_mod,  ONLY:                                         &
    npft, nnvg, ntype

  USE jules_vegetation_mod,     ONLY:                                         &
    l_triffid, l_phenol, irr_crop, l_nitrogen, l_irrig_dmd

  USE nvegparm,                 ONLY:                                         &
    albsnc_nvg, albsnf_nvgu, albsnf_nvg, albsnf_nvgl, catch_nvg, emis_nvg,    &
    gs_nvg, infil_nvg, z0_nvg, ch_nvg, vf_nvg

  USE ozone_vars,               ONLY:                                         &
    o3_gb, flux_o3_pft, fo3_pft

  USE pftparm,                  ONLY:                                         &
    albsnc_max, albsnc_min, albsnf_maxu, albsnf_max, albsnf_maxl, alpha,      &
    alniru, alnir, alnirl, alparu, alpar, alparl, a_wl, a_ws, b_wl, catch0,   &
    c3, dcatch_dlai, dgl_dm, dgl_dt, dqcrit, dz0v_dh, emis_pft,eta_sl, fd,    &
    fsmc_of, f0, glmin, g_leaf_0, infil_f, kext, kpar, lai_alb_lim, neff,     &
    nl0, nr_nl, ns_nl, omegau, omega, omegal, omniru, omnir, omnirl, orient,  &
    r_grow, rootd_ft, sigl, tleaf_of, tlow, tupp, lma, nmass, vsl, vint, kn,  &
    knl, q10_leaf, fl_o3_ct, dfp_dcuo, ci_st, gpp_st, ief, tef, mef, aef,     &
    fef_co2, fef_co, fef_ch4, fef_nox, fef_so2, fef_oc, fef_bc, ccleaf_min,   &
    ccleaf_max, ccwood_min, ccwood_max, avg_ba,                               &
    nsw, nr, hw_sw, fsmc_p0, fsmc_mod, can_struct_a

  USE prognostics,              ONLY:                                         &
    lai_pft, canht_pft, smcl_gb, t_soil_gb, tsurf_elev_surft,                 &
    rgrain_surft, snow_surft, soot_ij,                                        &
    tstar_surft, canopy_surft, canopy_gb, cs_pool_gb, ti_sicat, z0msea_ij,    &
    gs_gb, gc_surft, smc_gb, di_ncat_sicat, k_sice_sicat, snow_grnd_surft,    &
    snow_mass_ij, snow_mass_sea_sicat, nsnow_surft, rho_snow_grnd_surft,      &
    snowdepth_surft, ds_surft, rgrainl_surft, sice_surft, sliq_surft,         &
    tsnow_surft, wood_prod_fast_gb, wood_prod_med_gb, wood_prod_slow_gb,      &
    frac_agr_prev_gb, frac_past_prev_gb,                                      &
    !This is required in the UM dump file but not for standalone.
    !It may be redundant information and therefore could be removed
    rho_snow_surft, tsoil_deep_gb, ns_pool_gb, n_inorg_gb

  USE p_s_parms,                ONLY:                                         &
    bexp_gb, sathh_gb, hcap_gb, hcon_gb, satcon_gb, smvccl_gb, smvcwt_gb,     &
    smvcst_gb, clay_gb

  USE switches_urban,           ONLY:                                         &
    l_moruses, l_urban2t

  USE theta_field_sizes,        ONLY:                                         &
    t_i_length, t_j_length

  USE trif,                     ONLY:                                         &
    crop, g_area, g_grow, g_root, g_wood, lai_max, lai_min, alloc_fast,       &
    alloc_med, alloc_slow, dpm_rpm_ratio, retran_r, retran_l

  USE trif_vars_mod,            ONLY:                                         &
    wp_fast_in_gb, wp_med_in_gb, wp_slow_in_gb, wp_fast_out_gb, wp_med_out_gb,&
    wp_slow_out_gb, lit_c_orig_pft, lit_c_ag_pft, n_leaf_pft, n_root_pft,     &
    n_stem_pft, resp_r_pft, resp_l_pft, lai_bal_pft, pc_s_pft, resp_s_diag_gb,&
    resp_s_pot_diag_gb, minl_n_gb, minl_n_pot_gb, immob_n_gb, immob_n_pot_gb, &
    fn_gb, leafC_pft, rootC_pft, stemC_pft, woodC_pft, droot_pft, dleaf_pft,  &
    dwood_pft, n_uptake_pft, n_demand_gb, n_uptake_growth_pft,                &
    n_demand_growth_pft, n_demand_spread_pft, n_uptake_spread_pft,            &
    n_demand_lit_pft, n_uptake_gb, n_demand_pft, exudates_pft, exudates_gb,   &
    dcveg_pft, dnveg_pft, dcveg_gb, dnveg_gb, n_veg_pft, n_veg_gb, n_loss_gb, &
    dpm_ratio_gb, root_litC_pft, leaf_litC_pft, wood_litC_pft, root_litN_pft, &
    leaf_litN_pft, wood_litN_pft, litterC_pft, litterN_pft, lit_N_pft,        &
    lit_N_t_gb, deposition_N_gb, n_fix_gb, n_fix_pft, n_leach_gb, n_gas_gb,   &
    fapar_diag_pft, fao_et0, cnsrv_carbon_veg2_gb, cnsrv_carbon_triffid_gb,   &
    cnsrv_veg_triffid_gb, cnsrv_soil_triffid_gb, cnsrv_prod_triffid_gb,       &
    root_abandon_pft, harvest_pft, harvest_gb, root_abandon_n_pft,            &
    harvest_n_pft, harvest_n_gb, n_fertiliser_pft, n_fertiliser_gb,           &
    frac_past_gb, lit_n_orig_pft, lit_n_ag_pft

  USE urban_param,              ONLY:                                         &
     hgt_gb, hwr_gb, wrr_gb, disp_gb, ztm_gb, albwl_gb, albrd_gb, emisw_gb,   &
     emisr_gb

  USE sf_diags_mod, ONLY: sf_diag

!Model-dependent USE statements
#if defined(UM_JULES)
  !Used for getting t_i_length and t_jlength and for FLAKE allocations
  USE atm_fields_bounds_mod,    ONLY:                                         &
    tdims

  !For FLAKE
  USE jules_surface_mod,        ONLY:                                         &
    l_flake_model

  !For FLAKE
  USE lake_mod,                 ONLY:                                         &
    surf_ht_flux_lake_ij, surf_ht_flux_lk_gb, sw_down_gb, coriolis_param_gb,  &
    u_s_lake_gb, lake_depth_gb, lake_fetch_gb, lake_albedo_gb, lake_t_snow_gb,&
    lake_t_ice_gb, lake_t_mean_gb, lake_t_mxl_gb, lake_shape_factor_gb,       &
    lake_h_snow_gb, lake_h_ice_gb, lake_h_mxl_gb, lake_t_sfc_gb, ts1_lake_gb, &
    nusselt_gb, g_dt_gb

  USE pftparm,                  ONLY:                                         &
    !dust_veg_scj is only used in the boundary layer scheme?
    dust_veg_scj

#else
  USE aero,                     ONLY:                                         &
    co2_3d_ij, rho_cd_modv1_ij, rho_aresist_ij, aresist_ij, resist_b_ij,      &
    rho_aresist_surft, aresist_surft, resist_b_surft, r_b_dust_ij,            &
    cd_std_dust_ij, u_s_std_surft

  USE ancil_info,               ONLY:                                         &
    land_pts, surft_index, soil_index, lice_index, ice_fract_ij,              &
    ice_fract_ncat_sicat,                                                     &
    ti_cat_sicat, pond_frac_cat_sicat, pond_depth_cat_sicat, sstfrz_ij,       &
    z1_uv_ij, z1_tq_ij, nsurft, frac_surft, dim_cs1, land_pts_trif, npft_trif,&
    surft_pts

  USE c_elevate,                ONLY:                                         &
    z_land_ij

  USE coastal,                  ONLY:                                         &
    tstar_land_ij, tstar_sea_ij, tstar_sice_ij, tstar_sice_sicat,             &
    tstar_ssi_ij, taux_land_ij, taux_ssi_ij, tauy_land_ij, tauy_ssi_ij,       &
    vshr_land_ij, vshr_ssi_ij, surf_ht_flux_land_ij, surf_ht_flux_sice_sicat

  USE cropparm,                 ONLY:                                         &
    beta1, beta2, beta3, gamma, delta, remob, cfrac_s, cfrac_r, cfrac_l,      &
    allo1, allo2, t_bse, t_opt, t_max, tt_emr, crit_pp, pp_sens, rt_dir,      &
    alpha1, alpha2, alpha3, mu, nu, yield_frac, initial_carbon, sen_dvi, t_mort

  USE crop_vars_mod,            ONLY:                                         &
    phot, dphotdt, dvi_cpft, croprootc_cpft => rootc_cpft, harvc_cpft,        &
    reservec_cpft, yield_diag_cpft, nonyield_diag_cpft, leafc_diag_cpft,      &
    stemc_diag_cpft, croplai_cpft, cropcanht_cpft, sow_date_cpft,             &
    tt_veg_cpft, tt_rep_cpft, harvest_trigger_cpft, harvest_counter_cpft

  USE dust_parameters_mod,      ONLY:                                         &
    ndiv

  USE fluxes,                   ONLY:                                         &
    alb_surft, e_sea_ij, ecan_ij, ecan_surft, ei_ij, ei_surft, esoil_ij,      &
    esoil_surft, ext_gb, fqw_1_ij, fqw_surft, fqw_sicat, fsmc_pft, ftl_1_ij,  &
    ftl_sicat, ftl_surft, h_sea_ij, hf_snow_melt_gb, land_albedo_ij,          &
    le_surft, melt_surft, sea_ice_htf_sicat, snomlt_sub_htf_gb, snow_melt_gb, &
    snowmelt_ij, sub_surf_roff_gb, surf_ht_flux_ij, surf_htf_surft,           &
    surf_roff_gb, radnet_surft, taux_1_ij, tauy_1_ij, tot_tfall_gb, tstar_ij, &
    sw_surft, emis_surft, alb_sicat, rflow_gb, rrun_gb,                       &
    snow_soil_htf, snow_smb_surft

  USE forcing,                  ONLY:                                         &
    qw_1_ij, tl_1_ij, u_0_ij, v_0_ij, u_1_ij, v_1_ij, pstar_ij, ls_rain_ij,   &
    con_rain_ij, ls_snow_ij, con_snow_ij, sw_down_ij, lw_down_ij, diff_rad_ij,&
    diurnal_temperature_range_ij

  USE jules_sea_seaice_mod,     ONLY:                                         &
    nice, nice_use

  USE jules_soil_mod,           ONLY:                                         &
    sm_levels

  USE jules_surface_types_mod,  ONLY:                                         &
    ncpft

  USE jules_vegetation_mod,     ONLY:                                         &
    l_crop

  USE max_dimensions,           ONLY:                                         &
    npft_max

  USE orog,                     ONLY:                                         &
    sil_orog_land_gb, ho2r2_orog_gb, h_blend_orog_ij, z0m_eff_ij

  USE p_s_parms,                ONLY:                                         &
    albsoil_gb, albobs_sw_gb, albobs_vis_gb, albobs_nir_gb, catch_surft,      &
    catch_snow_surft, cosz_gb, infil_surft, z0_surft, z0h_bare_surft,         &
    z0m_soil_gb, sthu_gb, sthf_gb, sthu_min_gb

  USE top_pdm,                  ONLY:                                         &
    a_fsat_gb, a_fwet_gb, c_fsat_gb, c_fwet_gb, drain_gb, dun_roff_gb,        &
    fexp_gb, fsat_gb, fch4_wetl_gb, fch4_wetl_cs_gb, fch4_wetl_npp_gb,        &
    fch4_wetl_resps_gb, fwetl_gb, gamtot_gb, qbase_gb, qbase_zw_gb, sthzw_gb, &
    ti_mean_gb, ti_sig_gb, zw_gb, inlandout_atm_gb

  USE trifctl,                  ONLY:                                         &
    lai_phen_pft, c_veg_pft, cv_gb, g_leaf_day_pft, g_leaf_dr_out_pft,        &
    lit_c_pft, lit_c_mn_gb, npp_dr_out_pft, resp_w_dr_out_pft,                &
    resp_s_dr_out_gb, frac_agr_gb, g_leaf_acc_pft, npp_acc_pft,               &
    resp_w_acc_pft, resp_s_acc_gb, g_leaf_phen_acc_pft, gpp_gb, npp_gb,       &
    resp_p_gb, g_leaf_pft, g_leaf_phen_pft, gpp_pft, npp_pft, resp_p_pft,     &
    resp_s_gb, resp_w_pft

  USE u_v_grid,                 ONLY:                                         &
    u_0_p_ij, v_0_p_ij, u_1_p_ij, v_1_p_ij, dtrdz_charney_grid_1_ij

#endif

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Allocates the model arrays using sizes determined during initialisation
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!
! Modules are ordered alphabetically. Please respect it if you add something
! new.
!
!-----------------------------------------------------------------------------

#if defined(UM_JULES)
! Input variables for dimensioning
INTEGER, INTENT(IN) ::                                                        &
  land_pts,                                                                   &
                       ! Number of land points
  nsurft,                                                                     &
                       ! Number of surface tiles
  sm_levels,                                                                  &
                       ! Number of soil layers
  nice,                                                                       &
                       ! Number of sea ice categories
  nice_use
                       ! Number of sea ice cats used in radiation and
                       !  explicit part of surface exchange
#endif

#if defined(UM_JULES)
! In the UM, we must define dim_cs1 as a local variable and calculate
! it
! This is because it is only calculated in atm_step_phys_init, which
! is no good for allocating stuff
INTEGER :: dim_cs1
#endif

!-----------------------------------------------------------------------
! Local variables for error trapping
!-----------------------------------------------------------------------
INTEGER ::                                                                    &
  error        = 0,                                                           &
                       ! Variable for trapping the error from each
                       ! individual call to allocate
  error_sum    = 0,                                                           &

                       ! Variable to track the sum of all errors
                       ! resulting from calls to allocate. Hence we
                       ! know that everything was successful if and
                       ! only if this is zero at the end

  temp_size,                                                                  &
  temp_tiles,                                                                 &
  temp_layers
                       ! For storing the size of array to allocate for variables
                       ! that are sometimes set to size 1.
                       ! Removes some duplicate allocate statements

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ALLOCATE_JULES_ARRAYS'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

#if defined(UM_JULES)
!DO WE NEED THIS? Only existed in the UM version.
!-----------------------------------------------------------------------
! Compute length of theta field in i and j directions.
! (The module keeps these values for future use, this is
!  earliest place in the code that they are needed.)
!-----------------------------------------------------------------------
t_i_length = tdims%i_end - tdims%i_start + 1
t_j_length = tdims%j_end - tdims%j_start + 1

! Define dim_cs1
IF( l_triffid ) THEN
  dim_cs1 = 4
ELSE
  dim_cs1 = 1
END IF
#endif

!------------------------------------------------------------------------------
!Begin the array allocation. We will follow the order that the USE statements
!have been done.
!------------------------------------------------------------------------------

!USE statements that apply to both UM and JULES

!  ====ancil_info module common====
  ALLOCATE( ssi_index(t_i_length*t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( sea_index(t_i_length*t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( sice_index(t_i_length*t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( sice_pts_ncat(nice), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( fssi_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( sea_frac(t_i_length*t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( sice_frac(t_i_length*t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( sice_frac_ncat(t_i_length*t_j_length,nice), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( sice_index_ncat(t_i_length*t_j_length,nice), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( l_lice_point(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( l_soil_point(land_pts), STAT=error )
  error_sum = error_sum + error

!  ====bvoc_vars module common====
! BVOC diagnostics
  ALLOCATE( isoprene_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( isoprene_pft(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( terpene_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( terpene_pft(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( methanol_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( methanol_pft(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( acetone_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( acetone_pft(land_pts,npft), STAT=error )
  error_sum = error_sum + error

!  ====c_elevate module common====
! Height above mean grid-box
  ALLOCATE( surf_hgt_surft(land_pts,nsurft), STAT=error )
  error_sum = error_sum + error

  ALLOCATE( lw_down_elevcorr_surft(land_pts,nsurft), STAT=error )
  error_sum = error_sum + error

!  ====crop_vars_mod module common====
! Irrigation variables
  ALLOCATE( sthu_irr_gb(land_pts,sm_levels), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( frac_irr_all(land_pts,1), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( frac_irr_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( frac_irr_old_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( frac_irr_surft(land_pts,nsurft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( plant_n_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( smc_irr_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( wt_ext_irr_surft(land_pts,sm_levels,nsurft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( gs_irr_surft(land_pts,nsurft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( dvimax_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( gc_irr_surft(land_pts,nsurft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( resfs_irr_surft(land_pts,nsurft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( ext_irr_gb(land_pts,sm_levels), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( wt_ext_irr_gb(land_pts,sm_levels), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( fsmc_irr_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( irrDaysDiag_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( irrig_water_gb(land_pts), STAT=error )
  error_sum = error_sum + error

  IF ( error_sum == 0 ) THEN
    sthu_irr_gb(:,:)           = 0.0
    frac_irr_all(:,:)          = 0.0
    frac_irr_gb(:)             = 0.0
    frac_irr_old_gb(:)         = 0.0
    frac_irr_surft(:,:)        = 0.0
    plant_n_gb(:)              = 0
    smc_irr_gb(:)              = 0.0
    wt_ext_irr_surft(:,:,:)    = 0.0
    gs_irr_surft(:,:)          = 0.0
    dvimax_gb(:)               = 0.0
    gc_irr_surft(:,:)          = 0.0
    resfs_irr_surft(:,:)       = 0.0
    ext_irr_gb(:,:)            = 0.0
    wt_ext_irr_gb(:,:)         = 0.0
    fsmc_irr_gb(:)             = 0.0
    irrDaysDiag_gb(:)          = 0.0
    irrig_water_gb(:)          = 0.0
  END IF

!-----------------------------------------------------------------------
! Allocate space for inferno diagnostic variables
!-----------------------------------------------------------------------
  ALLOCATE( burnt_area(land_pts), STAT = error )
  error_sum = error_sum + error
  ALLOCATE( burnt_area_ft(land_pts,npft), STAT = error )
  error_sum = error_sum + error

  ALLOCATE( emitted_carbon(land_pts), STAT = error )
  error_sum = error_sum + error
  ALLOCATE( emitted_carbon_ft(land_pts,npft), STAT = error )
  error_sum = error_sum + error
  ALLOCATE( emitted_carbon_DPM(land_pts), STAT = error )
  error_sum = error_sum + error
  ALLOCATE( emitted_carbon_RPM(land_pts), STAT = error )
  error_sum = error_sum + error

  ALLOCATE( fire_em_CO2(land_pts), STAT = error )
  error_sum = error_sum + error
  ALLOCATE( fire_em_CO2_ft(land_pts,npft), STAT = error )
  error_sum = error_sum + error
  ALLOCATE( fire_em_CO2_DPM(land_pts), STAT = error )
  error_sum = error_sum + error
  ALLOCATE( fire_em_CO2_RPM(land_pts), STAT = error )
  error_sum = error_sum + error

  ALLOCATE( fire_em_CO(land_pts), STAT = error )
  error_sum = error_sum + error
  ALLOCATE( fire_em_CO_ft(land_pts,npft), STAT = error )
  error_sum = error_sum + error
  ALLOCATE( fire_em_CO_DPM(land_pts), STAT = error )
  error_sum = error_sum + error
  ALLOCATE( fire_em_CO_RPM(land_pts), STAT = error )
  error_sum = error_sum + error

  ALLOCATE( fire_em_CH4(land_pts), STAT = error )
  error_sum = error_sum + error
  ALLOCATE( fire_em_CH4_ft(land_pts,npft), STAT = error )
  error_sum = error_sum + error
  ALLOCATE( fire_em_CH4_DPM(land_pts), STAT = error )
  error_sum = error_sum + error
  ALLOCATE( fire_em_CH4_RPM(land_pts), STAT = error )
  error_sum = error_sum + error

  ALLOCATE( fire_em_NOx(land_pts), STAT = error )
  error_sum = error_sum + error
  ALLOCATE( fire_em_NOx_ft(land_pts,npft), STAT = error )
  error_sum = error_sum + error
  ALLOCATE( fire_em_NOx_DPM(land_pts), STAT = error )
  error_sum = error_sum + error
  ALLOCATE( fire_em_NOx_RPM(land_pts), STAT = error )
  error_sum = error_sum + error

  ALLOCATE( fire_em_SO2(land_pts), STAT = error )
  error_sum = error_sum + error
  ALLOCATE( fire_em_SO2_ft(land_pts,npft), STAT = error )
  error_sum = error_sum + error
  ALLOCATE( fire_em_SO2_DPM(land_pts), STAT = error )
  error_sum = error_sum + error
  ALLOCATE( fire_em_SO2_RPM(land_pts), STAT = error )
  error_sum = error_sum + error

  ALLOCATE( fire_em_OC(land_pts), STAT = error )
  error_sum = error_sum + error
  ALLOCATE( fire_em_OC_ft(land_pts,npft), STAT = error )
  error_sum = error_sum + error
  ALLOCATE( fire_em_OC_DPM(land_pts), STAT = error )
  error_sum = error_sum + error
  ALLOCATE( fire_em_OC_RPM(land_pts), STAT = error )
  error_sum = error_sum + error

  ALLOCATE( fire_em_BC(land_pts), STAT = error )
  error_sum = error_sum + error
  ALLOCATE( fire_em_BC_ft(land_pts,npft), STAT = error )
  error_sum = error_sum + error
  ALLOCATE( fire_em_BC_DPM(land_pts), STAT = error )
  error_sum = error_sum + error
  ALLOCATE( fire_em_BC_RPM(land_pts), STAT = error )
  error_sum = error_sum + error

  ALLOCATE( pop_den(land_pts), STAT = error )
  error_sum = error_sum + error

  ALLOCATE( flash_rate(land_pts), STAT = error )
  error_sum = error_sum + error


! vars for Doell & Siebert crop calendar
  IF ( irr_crop == 1 ) THEN
    ALLOCATE( icntmax_gb(land_pts), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( tl_1_day_av_gb(land_pts), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( tl_1_day_av_use_gb(land_pts,ndpy,nyav), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( prec_1_day_av_gb(land_pts), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( prec_1_day_av_use_gb(land_pts,ndpy,nyav), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( rn_1_day_av_gb(land_pts), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( rn_1_day_av_use_gb(land_pts,ndpy,nyav), STAT=error )
    error_sum = error_sum + error

    IF ( error_sum == 0 ) THEN
      icntmax_gb(:)                = 0
      tl_1_day_av_gb(:)            = 0.0
      tl_1_day_av_use_gb(:,:,:)    = 0.0
      prec_1_day_av_gb(:)          = 0.0
      prec_1_day_av_use_gb(:,:,:)  = 0.0
      rn_1_day_av_gb(:)            = 0.0
      rn_1_day_av_use_gb(:,:,:)    = 0.0
    END IF
  END IF

!  ====c_z0h_z0m module, common
! Surface type variables.
  ALLOCATE( z0h_z0m(ntype), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( z0h_z0m_classic(ntype), STAT=error )
  error_sum = error_sum + error

!  ====fluxes module, common
  ALLOCATE( surf_ht_store_surft(land_pts,nsurft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( anthrop_heat_surft(land_pts,nsurft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( sw_rts_sicat(t_i_length*t_j_length, nice_use), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( swup_rts_sicat(t_i_length*t_j_length, nice_use), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( swdn_rts_sicat(t_i_length*t_j_length, nice_use), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( sw_sicat(t_i_length*t_j_length, nice_use), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( alb_sicat(t_i_length*t_j_length, nice_use, 4), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( et_stom_pft(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( et_stom_gb(land_pts), STAT=error )
  error_sum = error_sum + error

  IF ( error_sum == 0 ) THEN
    anthrop_heat_surft(:,:) = 0.0
    alb_sicat(:,:,:)        = 0.0
  END IF

!  ====jules_mod module common====
  ALLOCATE( snowdep_surft(land_pts,nsurft), STAT=error )
  error_sum = error_sum + error
  IF ( l_albedo_obs .AND. l_spec_albedo) THEN
    ALLOCATE( albobs_scaling_surft(land_pts,ntype,2), STAT=error )
    error_sum = error_sum + error
  ELSE IF ( l_albedo_obs .AND. (.NOT.l_spec_albedo) ) THEN
    ALLOCATE( albobs_scaling_surft(land_pts,ntype,1), STAT=error )
    error_sum = error_sum + error
  END IF

!========jules_radiation_mod module========
!Nothing to allocate

!========jules_snow_mod module========
  IF (ANY(cansnowtile(1:npft)) .EQV. .TRUE.) THEN
    ALLOCATE( unload_backgrnd_pft(land_pts,npft), STAT=error )
    error_sum = error_sum + error
  ELSE
    ALLOCATE( unload_backgrnd_pft(1,1), STAT=error )
    error_sum = error_sum + error
  END IF

!  ====jules_surface_mod module common====
  ALLOCATE( diff_frac(t_i_length*t_j_length), STAT=error )
  error_sum = error_sum + error

!========jules_surface_types_mod module========
!Nothing to allocate

!  ====jules_vegetation_mod module UM====
!Nothing to allocate

!  ====nvegparm module common====
! Non-veg surface type variables.
  ALLOCATE( albsnc_nvg(nnvg), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( albsnf_nvgu(nnvg), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( albsnf_nvg(nnvg), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( albsnf_nvgl(nnvg), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( catch_nvg(nnvg), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( emis_nvg(nnvg), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( gs_nvg(nnvg), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( infil_nvg(nnvg), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( z0_nvg(nnvg), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( ch_nvg(nnvg), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( vf_nvg(nnvg), STAT=error )
  error_sum = error_sum + error

!  ====ozone_vars module common====
! Ozone variables
  ALLOCATE( o3_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( flux_o3_pft(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( fo3_pft(land_pts,npft), STAT=error )
  error_sum = error_sum + error

  IF ( error_sum == 0 ) THEN
    o3_gb(:) = 0.0
  END IF

!  ====pftparm module common====
! Veg surface type variables.
  ALLOCATE( albsnc_max(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( albsnc_min(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( albsnf_maxu(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( albsnf_max(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( albsnf_maxl(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( alpha(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( alniru(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( alnir(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( alnirl(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( alparu(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( alpar(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( alparl(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( a_wl(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( a_ws(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( b_wl(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( catch0(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( c3(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( dcatch_dlai(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( dgl_dm(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( dgl_dt(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( dqcrit(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( dz0v_dh(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( emis_pft(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( eta_sl(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( fd(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( fsmc_of(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( f0(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( glmin(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( g_leaf_0(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( infil_f(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( kext(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( kpar(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( lai_alb_lim(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( neff(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( nl0(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( nr_nl(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( ns_nl(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( nsw(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( nr(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( hw_sw(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( can_struct_a(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( omegau(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( omega(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( omegal(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( omniru(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( omnir(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( omnirl(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( orient(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( r_grow(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( rootd_ft(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( fsmc_p0(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( fsmc_mod(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( sigl(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( tleaf_of(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( tlow(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( tupp(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( lma(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( nmass(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( vsl(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( vint(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( kn(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( knl(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( q10_leaf(npft), STAT=error )
  error_sum = error_sum + error

  IF( error_sum == 0 ) THEN
    fsmc_p0(:)    = 0.0
    fsmc_mod(:)   = 0
  END IF

! Ozone damage parameters
  ALLOCATE( fl_o3_ct(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( dfp_dcuo(npft), STAT=error )
  error_sum = error_sum + error
! BVOC emission parameters
  ALLOCATE( ci_st(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( gpp_st(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( ief(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( tef(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( mef(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( aef(npft), STAT=error )
  error_sum = error_sum + error
! INFERNO combustion parameters
  ALLOCATE( ccleaf_min(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( ccleaf_max(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( ccwood_min(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( ccwood_max(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( avg_ba(npft), STAT=error )
  error_sum = error_sum + error
! INFERNO emission parameters
  ALLOCATE( fef_co2(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( fef_co(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( fef_ch4(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( fef_nox(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( fef_so2(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( fef_oc(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( fef_bc(npft), STAT=error )
  error_sum = error_sum + error

!  ====prognostics module common====
  ALLOCATE( nsnow_surft(land_pts,nsurft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( rho_snow_grnd_surft(land_pts,nsurft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( snowdepth_surft(land_pts,nsurft), STAT=error )
  error_sum = error_sum + error

! Nitrogen scheme variables
  ALLOCATE( ns_pool_gb(land_pts,dim_cs1), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( n_inorg_gb(land_pts), STAT=error )
  error_sum = error_sum + error

  IF( error_sum == 0 ) THEN
    ns_pool_gb(:,:)   = 0.001
    n_inorg_gb(:)     = 0.0
  END IF

  !For multilayer snow variables, only allocate to full size if the scheme is
  !being used, ie nsmax > 0
  IF (nsmax > 0)THEN
    temp_size   = land_pts
    temp_tiles  = nsurft
    temp_layers = nsmax
  ELSE
    temp_size   = 1
    temp_tiles  = 1
    temp_layers = 1
  END IF

  ALLOCATE( ds_surft(temp_size,temp_tiles,temp_layers),      STAT=error )
  error_sum = error_sum + error
  ALLOCATE( rgrainl_surft(temp_size,temp_tiles,temp_layers), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( sice_surft(temp_size,temp_tiles,temp_layers),    STAT=error )
  error_sum = error_sum + error
  ALLOCATE( sliq_surft(temp_size,temp_tiles,temp_layers),    STAT=error )
  error_sum = error_sum + error
  ALLOCATE( tsnow_surft(temp_size,temp_tiles,temp_layers),   STAT=error )
  error_sum = error_sum + error

  !See comment in prognostics module
  ALLOCATE( rho_snow_surft(land_pts,nsurft,nsmax), STAT=error )
  error_sum = error_sum + error

! Allocate WP Pools
  IF ( l_triffid .OR. l_phenol ) THEN
    ALLOCATE( wood_prod_fast_gb(land_pts), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( wood_prod_med_gb(land_pts), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( wood_prod_slow_gb(land_pts), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( frac_agr_prev_gb(land_pts), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( frac_past_prev_gb(land_pts), STAT=error )
    error_sum = error_sum + error

    IF ( error_sum == 0 ) THEN
      wood_prod_fast_gb(:) = 0.0
      wood_prod_med_gb(:)  = 0.0
      wood_prod_slow_gb(:) = 0.0
      frac_agr_prev_gb(:)  = 0.0
      frac_past_prev_gb(:) = 0.0
    END IF
  END IF

! Only allocate the bedrock tsoil_deep_gb if bedrock is being used
  IF( l_bedrock ) THEN
    ALLOCATE( tsoil_deep_gb(land_pts,ns_deep), STAT=error )
    error_sum = error_sum + error
  END IF

!  ====p_s_parm module UM====
  ALLOCATE( bexp_gb(land_pts,sm_levels), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( sathh_gb(land_pts,sm_levels), STAT=error )
  error_sum = error_sum + error
  ALLOCATE(  hcap_gb(land_pts,sm_levels), STAT=error )
  error_sum = error_sum + error
  ALLOCATE(  hcon_gb(land_pts,0:sm_levels), STAT=error )
  error_sum = error_sum + error
  ALLOCATE(satcon_gb(land_pts,0:sm_levels), STAT=error )
  error_sum = error_sum + error
  ALLOCATE(smvccl_gb(land_pts,sm_levels), STAT=error )
  error_sum = error_sum + error
  ALLOCATE(smvcwt_gb(land_pts,sm_levels), STAT=error )
  error_sum = error_sum + error
  ALLOCATE(smvcst_gb(land_pts,sm_levels), STAT=error )
  error_sum = error_sum + error

!  ====switches_urban module common====
!Nothing to allocate

!========theta_field_sizes module========
!Nothing to allocate


!  ====trif module common====
! TRIFFID variables - only needed if TRIFFID and/or phenology is selected.
  IF ( l_triffid .OR. l_phenol ) THEN
    ALLOCATE( crop(npft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( g_area(npft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( g_grow(npft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( g_root(npft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( g_wood(npft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( lai_max(npft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( lai_min(npft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( alloc_fast(npft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( alloc_med(npft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( alloc_slow(npft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( dpm_rpm_ratio(npft), STAT=error )
    error_sum = error_sum + error
  END IF

!  ====trif_vars_mod module common====
  ALLOCATE( resp_l_pft(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( resp_r_pft(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( n_leaf_pft(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( n_root_pft(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( n_stem_pft(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( lai_bal_pft(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( pc_s_pft(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( fapar_diag_pft(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( fao_et0(land_pts), STAT=error )
  error_sum = error_sum + error

  IF ( error_sum == 0 ) THEN
    fao_et0(:) = 0.0
  END IF

  IF ( l_triffid .OR. l_phenol ) THEN
    ALLOCATE( wp_fast_in_gb(land_pts), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( wp_med_in_gb(land_pts), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( wp_slow_in_gb(land_pts), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( wp_fast_out_gb(land_pts), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( wp_med_out_gb(land_pts), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( wp_slow_out_gb(land_pts), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( lit_c_orig_pft(land_pts,npft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( lit_c_ag_pft(land_pts,npft), STAT=error )
    error_sum = error_sum + error

    ALLOCATE( cnsrv_carbon_veg2_gb(land_pts), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( cnsrv_carbon_triffid_gb(land_pts), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( cnsrv_veg_triffid_gb(land_pts), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( cnsrv_soil_triffid_gb(land_pts), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( cnsrv_prod_triffid_gb(land_pts), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( root_abandon_pft(land_pts,npft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( harvest_pft(land_pts,npft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( harvest_gb(land_pts), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( root_abandon_n_pft(land_pts,npft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( harvest_n_pft(land_pts,npft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( harvest_n_gb(land_pts), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( n_fertiliser_pft(land_pts,npft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( n_fertiliser_gb(land_pts), STAT=error )
    error_sum = error_sum + error

! Variables added for #7 (nitrogen scheme)
    ALLOCATE( resp_s_diag_gb(land_pts,dim_cs1+1), stat=error )
    error_sum = error_sum + error
    ALLOCATE( resp_s_pot_diag_gb(land_pts,dim_cs1+1), stat=error )
    error_sum = error_sum + error
    ALLOCATE ( n_veg_pft(land_pts,npft), stat=error)
    error_sum = error_sum + error
    ALLOCATE ( n_veg_gb(land_pts), stat=error)
    error_sum = error_sum + error
    ALLOCATE ( dnveg_pft(land_pts,npft), stat=error)
    error_sum = error_sum + error
    ALLOCATE ( dcveg_pft(land_pts,npft), stat=error)
    error_sum = error_sum + error
    ALLOCATE ( dnveg_gb(land_pts), stat=error)
    error_sum = error_sum + error
    ALLOCATE ( dcveg_gb(land_pts), stat=error)
    error_sum = error_sum + error
    ALLOCATE ( n_demand_gb(land_pts), stat=error)
    error_sum = error_sum + error
    ALLOCATE ( n_uptake_gb(land_pts), stat=error)
    error_sum = error_sum + error
    ALLOCATE ( deposition_N_gb(land_pts), stat=error)
    error_sum = error_sum + error
    ALLOCATE ( n_uptake_pft(land_pts,npft), stat=error)
    error_sum = error_sum + error
    ALLOCATE ( n_demand_pft(land_pts,npft), stat=error)
    error_sum = error_sum + error
    ALLOCATE ( dleaf_pft(land_pts,npft), stat=error)
    error_sum = error_sum + error
    ALLOCATE ( droot_pft(land_pts,npft), stat=error)
    error_sum = error_sum + error
    ALLOCATE ( dwood_pft(land_pts,npft), stat=error)
    error_sum = error_sum + error
    ALLOCATE ( n_demand_lit_pft(land_pts,npft), stat=error)
    error_sum = error_sum + error
    ALLOCATE ( n_demand_spread_pft(land_pts,npft), stat=error)
    error_sum = error_sum + error
    ALLOCATE ( n_demand_growth_pft(land_pts,npft), stat=error)
    error_sum = error_sum + error
    ALLOCATE ( n_uptake_growth_pft(land_pts,npft), stat=error)
    error_sum = error_sum + error
    ALLOCATE ( n_uptake_spread_pft(land_pts,npft), stat=error)
    error_sum = error_sum + error
    ALLOCATE ( n_fix_gb(land_pts), stat=error)
    error_sum = error_sum + error
    ALLOCATE ( n_fix_pft(land_pts,npft), stat=error)
    error_sum = error_sum + error
    ALLOCATE ( n_leach_gb(land_pts), stat=error)
    error_sum = error_sum + error
    ALLOCATE ( n_gas_gb(land_pts), stat=error)
    error_sum = error_sum + error  
    ALLOCATE ( exudates_pft(land_pts,npft), stat=error)
    error_sum = error_sum + error
    ALLOCATE ( exudates_gb(land_pts), stat=error)
    error_sum = error_sum + error
    ALLOCATE ( root_litC_pft(land_pts,npft), stat=error)
    error_sum = error_sum + error
    ALLOCATE ( leaf_litC_pft(land_pts,npft), stat=error)
    error_sum = error_sum + error
    ALLOCATE ( wood_litC_pft(land_pts,npft), stat=error)
    error_sum = error_sum + error
    ALLOCATE ( root_litN_pft(land_pts,npft), stat=error)
    error_sum = error_sum + error
    ALLOCATE ( leaf_litN_pft(land_pts,npft), stat=error)
    error_sum = error_sum + error
    ALLOCATE ( wood_litN_pft(land_pts,npft), stat=error)
    error_sum = error_sum + error
    ALLOCATE ( lit_N_pft(land_pts,npft), stat=error)
    error_sum = error_sum + error
    ALLOCATE ( lit_N_t_gb(land_pts), stat=error)
    error_sum = error_sum + error
    ALLOCATE( minl_n_pot_gb(land_pts,dim_cs1+1), stat=error )
    error_sum = error_sum + error
    ALLOCATE( immob_n_gb(land_pts,dim_cs1+1), stat=error )
    error_sum = error_sum + error
    ALLOCATE( immob_n_pot_gb(land_pts,dim_cs1+1), stat=error )
    error_sum = error_sum + error
    ALLOCATE( fn_gb(land_pts), stat=error )
    error_sum = error_sum + error
    ALLOCATE( minl_n_gb(land_pts,dim_cs1+1), stat=error )
    error_sum = error_sum + error
    ALLOCATE ( dpm_ratio_gb(land_pts), stat=error)
    error_sum = error_sum + error
    ALLOCATE ( n_loss_gb(land_pts), stat=error)
    error_sum = error_sum + error
    ALLOCATE ( leafC_pft(land_pts,npft), stat=error)
    error_sum = error_sum + error
    ALLOCATE ( rootC_pft(land_pts,npft), stat=error)
    error_sum = error_sum + error
    ALLOCATE ( woodC_pft(land_pts,npft), stat=error)
    error_sum = error_sum + error
    ALLOCATE ( litterC_pft(land_pts,npft), stat=error)
    error_sum = error_sum + error
    ALLOCATE ( stemC_pft(land_pts,npft), stat=error)
    error_sum = error_sum + error
    ALLOCATE ( litterN_pft(land_pts,npft), stat=error)
    error_sum = error_sum + error
    ALLOCATE ( retran_r(npft), stat=error)
    error_sum = error_sum + error
    ALLOCATE ( retran_l(npft), stat=error)
    error_sum = error_sum + error
    ALLOCATE ( lit_n_orig_pft(land_pts,npft), stat=error)
    error_sum = error_sum + error
    ALLOCATE ( lit_n_ag_pft(land_pts,npft), stat=error)
    error_sum = error_sum + error
    IF ( error_sum == 0 ) THEN
      n_loss_gb(:) = 0.0
      n_leach_gb(:) = 0.0
      n_gas_gb(:) = 0.0
      n_fix_pft(:,:) = 0.0
      n_fix_gb(:) = 0.0
    END IF
  END IF

!  ====urban_param common====
  ALLOCATE( wrr_gb(land_pts), STAT=error )
  error_sum = error_sum + error

!  ====urban_param module UM====
  IF ( l_urban2t .OR. l_moruses ) THEN
    CALL jules_print(                                                         &
        'allocate_jules_arrays',                                              &
        'Allocating URBAN-2T / MORUSES arrays',                               &
        level=PrNorm)
    temp_size = land_pts
  ELSE
    temp_size = 1
  END IF

  ALLOCATE( hgt_gb(temp_size), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( hwr_gb(temp_size), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( disp_gb(temp_size), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( ztm_gb(temp_size), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( albwl_gb(temp_size), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( albrd_gb(temp_size), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( emisw_gb(temp_size), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( emisr_gb(temp_size), STAT=error )
  error_sum = error_sum + error

  IF  ( error_sum == 0 ) THEN
    hgt_gb(:)   = 0.0
    hwr_gb(:)   = 0.0
    albwl_gb(:) = 0.0
    albrd_gb(:) = 0.0
    emisw_gb(:) = 0.0
    emisr_gb(:) = 0.0
    ztm_gb(:)   = 0.0
    disp_gb(:)  = 0.0
  END IF

!Model-dependent USE statements

#if defined(UM_JULES)
!  ====atm_fields_bounds_mod module UM
!Nothing to allocate

!  ====jules_surface_mod module UM====
!Nothing to allocate

!  ====lake_mod module UM====
  IF ( l_flake_model ) THEN
    ALLOCATE( surf_ht_flux_lake_ij(tdims%i_start:tdims%i_end                  &
                               ,tdims%j_start:tdims%j_end), STAT=error )
    error_sum = error_sum + error
    temp_size = land_pts
  ELSE
    ALLOCATE( surf_ht_flux_lake_ij(1,1), STAT=error )
    error_sum = error_sum + error
    temp_size = 1
  END IF

  ALLOCATE( surf_ht_flux_lk_gb(temp_size), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( sw_down_gb(temp_size), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( coriolis_param_gb(temp_size), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( u_s_lake_gb(temp_size), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( lake_depth_gb(temp_size), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( lake_fetch_gb(temp_size), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( lake_albedo_gb(temp_size), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( lake_t_snow_gb(temp_size), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( lake_t_ice_gb(temp_size), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( lake_t_mean_gb(temp_size), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( lake_t_mxl_gb(temp_size), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( lake_shape_factor_gb(temp_size), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( lake_h_snow_gb(temp_size), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( lake_h_ice_gb(temp_size), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( lake_h_mxl_gb(temp_size), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( lake_t_sfc_gb(temp_size), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( ts1_lake_gb(temp_size), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( nusselt_gb(temp_size), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( g_dt_gb(temp_size), STAT=error )
  error_sum = error_sum + error

!  ====pftparm module UM====
  ALLOCATE( dust_veg_scj(npft), STAT=error )
  error_sum = error_sum + error

#else
!  ====aero module JULES====
  ALLOCATE( co2_3d_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( rho_cd_modv1_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( rho_aresist_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( aresist_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( resist_b_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( rho_aresist_surft(land_pts,nsurft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( aresist_surft(land_pts,nsurft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( resist_b_surft(land_pts,nsurft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( r_b_dust_ij(t_i_length,t_j_length,ndiv), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( cd_std_dust_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( u_s_std_surft(land_pts,nsurft), STAT=error )
  error_sum = error_sum + error

!  ====ancil_info module JULES====
  ALLOCATE( surft_index(land_pts,ntype), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( soil_index(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( lice_index(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( ice_fract_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( ice_fract_ncat_sicat(t_i_length,t_j_length,nice), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( ti_cat_sicat(t_i_length,t_j_length,nice), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( pond_frac_cat_sicat(t_i_length,t_j_length,nice), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( pond_depth_cat_sicat(t_i_length,t_j_length,nice), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( sstfrz_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( z1_uv_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( z1_tq_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( frac_surft(land_pts,ntype), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( surft_pts(ntype), STAT=error )
  error_sum = error_sum + error

!  ====c_elevate module JULES====
! Land height
  ALLOCATE( z_land_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error

!  ====coastal module JULES====
! Coastal tiling variables
  ALLOCATE( tstar_land_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( tstar_sea_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( tstar_sice_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( tstar_sice_sicat(t_i_length,t_j_length,nice_use), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( tstar_ssi_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( taux_land_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( taux_ssi_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( tauy_land_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( tauy_ssi_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( vshr_land_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( vshr_ssi_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( surf_ht_flux_land_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( surf_ht_flux_sice_sicat(t_i_length,t_j_length,nice), STAT=error )
  error_sum = error_sum + error
!  ===Allocate coupled model diagnostics in sf_diag
  IF (sf_diag%l_lw_up_sice_weighted_cat) THEN
    ALLOCATE( sf_diag%lw_up_sice_weighted_cat(t_i_length,t_j_length,nice_use), &
             STAT=error )
    error_sum = error_sum + error
    sf_diag%lw_up_sice_weighted_cat = 0.0
  ELSE
    ALLOCATE( sf_diag%lw_up_sice_weighted_cat(1,1,1), STAT=error )
    error_sum = error_sum + error
  ENDIF
  IF (sf_diag%l_lw_up_sice_weighted) THEN
    ALLOCATE( sf_diag%lw_up_sice_weighted(t_i_length,t_j_length), &
             STAT=error )
    error_sum = error_sum + error
    sf_diag%lw_up_sice_weighted = 0.0
  ELSE
    ALLOCATE( sf_diag%lw_up_sice_weighted(1,1), STAT=error )
    error_sum = error_sum + error
  ENDIF
  IF (sf_diag%l_ftl_ice_sm) THEN
    ALLOCATE( sf_diag%ftl_ice_sm(t_i_length,t_j_length), &
             STAT=error )
    error_sum = error_sum + error
    sf_diag%ftl_ice_sm = 0.0
  ELSE
    ALLOCATE( sf_diag%ftl_ice_sm(1,1), STAT=error )
    error_sum = error_sum + error
  ENDIF
  IF (sf_diag%l_tstar_sice_weighted_cat) THEN
    ALLOCATE( sf_diag%tstar_sice_weighted_cat(t_i_length,t_j_length,nice_use), &
             STAT=error )
    error_sum = error_sum + error
    sf_diag%tstar_sice_weighted_cat = 0.0
  ELSE
    ALLOCATE( sf_diag%tstar_sice_weighted_cat(1,1,1), STAT=error )
    error_sum = error_sum + error
  ENDIF
  IF (sf_diag%l_tstar_sice_weighted) THEN
    ALLOCATE( sf_diag%tstar_sice_weighted(t_i_length,t_j_length), &
             STAT=error )
    error_sum = error_sum + error
    sf_diag%tstar_sice_weighted = 0.0
  ELSE
    ALLOCATE( sf_diag%tstar_sice_weighted(1,1), STAT=error )
    error_sum = error_sum + error
  ENDIF
  IF (sf_diag%l_ice_present_cat) THEN
    ALLOCATE( sf_diag%ice_present_cat(t_i_length,t_j_length,nice_use), &
             STAT=error )
    error_sum = error_sum + error
    sf_diag%ice_present_cat = 0.0
  ELSE
    ALLOCATE( sf_diag%ice_present_cat(1,1,1), STAT=error )
    error_sum = error_sum + error
  ENDIF
  IF (sf_diag%l_ice_present) THEN
    ALLOCATE( sf_diag%ice_present(t_i_length,t_j_length), &
             STAT=error )
    error_sum = error_sum + error
    sf_diag%ice_present = 0.0
  ELSE
    ALLOCATE( sf_diag%ice_present(1,1), STAT=error )
    error_sum = error_sum + error
  ENDIF
  IF (sf_diag%l_cd_ssi) THEN
    ALLOCATE( sf_diag%cd_ssi(t_i_length,t_j_length), STAT=error )
    error_sum = error_sum + error
    sf_diag%cd_ssi = 0.0
  ELSE
    ALLOCATE( sf_diag%cd_ssi(1,1), STAT=error )
    error_sum = error_sum + error
  END IF
!
!  ====cropparm module JULES====
!  ====crop_vars_mod module JULES====
!  Done together as they share the l_crop IF
  IF ( l_crop ) THEN
  ! Crop variables
    ALLOCATE( t_bse(ncpft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( t_opt(ncpft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( t_max(ncpft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( tt_emr(ncpft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( crit_pp(ncpft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( pp_sens(ncpft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( rt_dir(ncpft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( alpha1(ncpft), STAT=error )
    error_sum = error_sum + error

    ALLOCATE( alpha2(ncpft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( alpha3(ncpft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( beta1(ncpft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( beta2(ncpft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( beta3(ncpft), STAT=error )
    error_sum = error_sum + error

    ALLOCATE( gamma(ncpft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( delta(ncpft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( remob(ncpft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( cfrac_s(ncpft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( cfrac_r(ncpft), STAT=error )
    error_sum = error_sum + error

    ALLOCATE( cfrac_l(ncpft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( allo1(ncpft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( allo2(ncpft), STAT=error )
    error_sum = error_sum + error
    
    ALLOCATE( mu(ncpft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( nu(ncpft), STAT=error )
    error_sum = error_sum + error

    ALLOCATE( yield_frac(ncpft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( initial_carbon(ncpft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( sen_dvi(ncpft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( t_mort(ncpft), STAT=error )
    error_sum = error_sum + error 

    ALLOCATE( dvi_cpft(land_pts,ncpft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( croprootc_cpft(land_pts,ncpft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( harvc_cpft(land_pts,ncpft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( reservec_cpft(land_pts,ncpft), STAT=error )
    error_sum = error_sum + error

    ALLOCATE( yield_diag_cpft(land_pts,ncpft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( nonyield_diag_cpft(land_pts,ncpft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( harvest_trigger_cpft(land_pts,ncpft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( harvest_counter_cpft(land_pts,ncpft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( leafc_diag_cpft(land_pts,ncpft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( stemc_diag_cpft(land_pts,ncpft), STAT=error )
    error_sum = error_sum + error

    IF ( error_sum == 0 ) THEN
      yield_diag_cpft(:,:)    = 0.0
      nonyield_diag_cpft(:,:) = 0.0
      harvest_trigger_cpft(:,:)    = 0
      harvest_counter_cpft(:,:)    = 0
      stemc_diag_cpft(:,:)    = 0.0
      leafc_diag_cpft(:,:)    = 0.0
    END IF

    ALLOCATE( croplai_cpft(land_pts,ncpft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( cropcanht_cpft(land_pts,ncpft), STAT=error )
    error_sum = error_sum + error

    ALLOCATE( sow_date_cpft(land_pts,ncpft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( tt_veg_cpft(land_pts,ncpft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( tt_rep_cpft(land_pts,ncpft), STAT=error )
    error_sum = error_sum + error

  ELSE
    ALLOCATE( dvi_cpft(land_pts,ncpft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( croprootc_cpft(land_pts,ncpft), STAT=error )
    error_sum = error_sum + error
  END IF

  ALLOCATE( phot(t_i_length*t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( dphotdt(t_i_length*t_j_length), STAT=error )
  error_sum = error_sum + error

  IF ( l_irrig_dmd ) THEN
    ALLOCATE(irrtiles(npft_max), STAT=error )
    error_sum = error_sum + error
    IF ( error_sum == 0 ) THEN
      irrtiles(:)    = 0.0
    END IF
  END IF

!  ====diag_swchs module JULES====
!Nothing to allocate

!  ====dust_parameters_mod module JULES====
!Nothing to allocate

!  ====fluxes module JULES====
! Runoff components.
  ALLOCATE( sub_surf_roff_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( surf_roff_gb(land_pts), STAT=error )
  error_sum = error_sum + error
! (Forcing) fluxes
  ALLOCATE( alb_surft(land_pts,nsurft,4) )
  error_sum = error_sum + error
  ALLOCATE( tstar_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( e_sea_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( fqw_1_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( fsmc_pft(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( ftl_1_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( ftl_surft(land_pts,nsurft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( le_surft(land_pts,nsurft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( h_sea_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( taux_1_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( tauy_1_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( fqw_surft(land_pts,nsurft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( fqw_sicat(t_i_length,t_j_length,nice_use), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( ftl_sicat(t_i_length,t_j_length,nice_use), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( ecan_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( esoil_surft(land_pts,nsurft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( sea_ice_htf_sicat(t_i_length,t_j_length,nice), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( surf_ht_flux_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( surf_htf_surft(land_pts,nsurft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( snow_soil_htf(land_pts,nsurft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( snow_smb_surft(land_pts,nsurft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( sf_diag%sice_mlt_htf(t_i_length,t_j_length,nice), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( sf_diag%snomlt_surf_htf(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( land_albedo_ij(t_i_length,t_j_length,4) )
  error_sum = error_sum + error
  ALLOCATE( sf_diag%latent_heat(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( ei_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( ei_surft(land_pts,nsurft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( ecan_surft(land_pts,nsurft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( esoil_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( ext_gb(land_pts,sm_levels), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( snowmelt_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( hf_snow_melt_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( radnet_surft(land_pts,nsurft) )
  error_sum = error_sum + error
  ALLOCATE( sw_surft(land_pts,nsurft) )
  error_sum = error_sum + error
  ALLOCATE( emis_surft(land_pts,nsurft) )
  error_sum = error_sum + error
  ALLOCATE( snow_melt_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( snomlt_sub_htf_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( tot_tfall_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( melt_surft(land_pts,nsurft), STAT=error )
  error_sum = error_sum + error

  IF ( error_sum == 0 ) THEN
    alb_surft(:,:,:)     = 0.0
    ftl_surft(:,:)       = 0.0
    le_surft(:,:)        = 0.0
    fqw_surft(:,:)       = 0.0
    esoil_surft(:,:)     = 273.15
    ecan_surft(:,:)      = 0.0
    ei_surft(:,:)        = 0.0
    melt_surft(:,:)      = 0.0
    radnet_surft(:,:)    = 0.0
    sw_surft(:,:)        = 0.0
    emis_surft(:,:)      = 0.0
    surf_htf_surft(:,:)  = 0.0
    snow_soil_htf(:,:)   = 0.0
    snow_smb_surft(:,:)   = 0.0
    alb_sicat(:,:,:)     = 0.0
    surf_roff_gb(:)      = 0.0
    sub_surf_roff_gb(:)  = 0.0
  END IF

! Equivalent neutral winds. Only allocated at full size if required;
! otherwise allocated a minimal size to ensure a valid address is
! assigned in subroutine calls.
  IF (sf_diag%suv10m_n) THEN
    ALLOCATE( sf_diag%mu10m_n(t_i_length,t_j_length), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( sf_diag%mv10m_n(t_i_length,t_j_length), STAT=error )
    error_sum = error_sum + error
  ELSE
    ALLOCATE( sf_diag%mu10m_n(1,1), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( sf_diag%mv10m_n(1,1), STAT=error )
    error_sum = error_sum + error
  END IF

  IF (sf_diag%sfme) THEN
    ALLOCATE (sf_diag%fme(t_i_length,t_j_length), STAT=error)
    sf_diag%fme = 0.0
    error_sum = error_sum + error
  ELSE
    ALLOCATE (sf_diag%fme(1,1), STAT=error)
    error_sum = error_sum + error
  END IF

!  ====forcing module JULES====
! Forcing variables
  ALLOCATE( qw_1_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( tl_1_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( u_0_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( v_0_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( u_1_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( v_1_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( pstar_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( ls_rain_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( con_rain_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( ls_snow_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( con_snow_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( sw_down_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( lw_down_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( diff_rad_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( diurnal_temperature_range_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
!
!  ====jules_sea_seaice_mod module JULES====
!Nothing to allocate

!  ====jules_soil_mod module JULES====
!Nothing to allocate

!  ====jules_surface_types_mod module JULES====
!Nothing to allocate

!  ====jules_vegetation_mod module JULES====
!Nothing to allocate

!  ====orog module JULES====
! Orographic roughness variables
  ALLOCATE( sil_orog_land_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( ho2r2_orog_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( h_blend_orog_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( z0m_eff_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error

!  ====prognostics module JULES====
  ALLOCATE( lai_pft(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( canht_pft(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( smcl_gb(land_pts,sm_levels), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( t_soil_gb(land_pts,sm_levels), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( tsurf_elev_surft(land_pts,nsurft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( rgrain_surft(land_pts,nsurft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( snow_surft(land_pts,nsurft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( soot_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( tstar_surft(land_pts,nsurft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( canopy_surft(land_pts,nsurft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( canopy_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( cs_pool_gb(land_pts,dim_cs1), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( ti_sicat(t_i_length,t_j_length,nice), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( z0msea_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( gs_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( gc_surft(land_pts,nsurft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( smc_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( di_ncat_sicat(t_i_length,t_j_length,nice), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( k_sice_sicat(t_i_length,t_j_length,nice), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( snow_grnd_surft(land_pts,nsurft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( snow_mass_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( snow_mass_sea_sicat(t_i_length,t_j_length,nice_use), STAT=error )
  error_sum = error_sum + error

  IF ( error_sum == 0 ) THEN
    snow_grnd_surft(:,:) = 0.0
  END IF

!  ====p_s_parms module JULES====
! Plant and soil parameters
  ALLOCATE( albsoil_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( clay_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( albobs_sw_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( albobs_vis_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( albobs_nir_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( catch_surft(land_pts,nsurft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( catch_snow_surft(land_pts,nsurft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( cosz_gb(t_i_length*t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( infil_surft(land_pts,nsurft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( z0_surft(land_pts,nsurft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( z0h_bare_surft(land_pts,nsurft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( z0m_soil_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( sthu_gb(land_pts,sm_levels), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( sthf_gb(land_pts,sm_levels), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( sthu_min_gb(land_pts,sm_levels), STAT=error )
  error_sum = error_sum + error

!  ====screen module JULES====
! Screen variables
  ALLOCATE( sf_diag%q1p5m(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( sf_diag%q1p5m_surft(land_pts,nsurft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( sf_diag%t1p5m(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( sf_diag%t1p5m_surft(land_pts,nsurft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( sf_diag%u10m(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( sf_diag%v10m(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  IF (sf_diag%suv10m_n) THEN
    ALLOCATE( sf_diag%u10m_n(t_i_length,t_j_length), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( sf_diag%v10m_n(t_i_length,t_j_length), STAT=error )
    error_sum = error_sum + error
  ELSE
    ALLOCATE( sf_diag%u10m_n(1,1), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( sf_diag%v10m_n(1,1), STAT=error )
    error_sum = error_sum + error
  END IF
  IF (sf_diag%l_t10m .OR. sf_diag%l_q10m) THEN
    ALLOCATE( sf_diag%chr10m(t_i_length,t_j_length), STAT=error )
    error_sum = error_sum + error
    sf_diag%chr10m = 0.0
  ELSE
    ALLOCATE( sf_diag%chr10m(1,1), STAT=error )
    error_sum = error_sum + error
  END IF
  IF (sf_diag%l_t10m) THEN
    ALLOCATE( sf_diag%t10m(t_i_length,t_j_length), STAT=error )
    error_sum = error_sum + error
    sf_diag%t10m = 0.0
  ELSE
    ALLOCATE( sf_diag%t10m(1,1), STAT=error )
    error_sum = error_sum + error
  END IF
  IF (sf_diag%l_q10m) THEN
    ALLOCATE( sf_diag%q10m(t_i_length,t_j_length), STAT=error )
    error_sum = error_sum + error
    sf_diag%q10m = 0.0
  ELSE
    ALLOCATE( sf_diag%q10m(1,1), STAT=error )
    error_sum = error_sum + error
  END IF

!  ====top_pdm module JULES====
!       TOPMODEL and PDM variables
  ALLOCATE( a_fsat_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( a_fwet_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( c_fsat_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( c_fwet_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( drain_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( dun_roff_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( fch4_wetl_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( fch4_wetl_cs_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( fch4_wetl_npp_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( fch4_wetl_resps_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( fexp_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( fsat_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( fwetl_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( gamtot_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( qbase_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( qbase_zw_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( sthzw_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( ti_mean_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( ti_sig_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( zw_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( inlandout_atm_gb(land_pts), STAT=error )
  error_sum = error_sum + error

!  ====trifctl module JULES====
! Triffid variables
  ALLOCATE( g_leaf_acc_pft(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( npp_acc_pft(land_pts_TRIF,npft_trif), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( resp_w_acc_pft(land_pts_TRIF,npft_trif), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( resp_s_acc_gb(land_pts_TRIF,DIM_CS1), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( g_leaf_phen_acc_pft(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( gpp_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( npp_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( resp_p_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( g_leaf_pft(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( g_leaf_phen_pft(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( gpp_pft(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( npp_pft(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( resp_p_pft(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( resp_s_gb(land_pts,dim_cs1), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( resp_w_pft(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( lai_phen_pft(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( c_veg_pft(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( cv_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( g_leaf_day_pft(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( g_leaf_dr_out_pft(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( lit_c_pft(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( lit_c_mn_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( npp_dr_out_pft(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( resp_w_dr_out_pft(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( resp_s_dr_out_gb(land_pts,5), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( frac_agr_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( frac_past_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  IF ( error_sum == 0 ) THEN
    frac_past_gb(:) = 0.0
    frac_agr_gb(:) = 0.0
  END IF

!  ====u_v_grid module JULES====
! Grid-change variables
  ALLOCATE( u_0_p_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( v_0_p_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( u_1_p_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( v_1_p_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( dtrdz_charney_grid_1_ij(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error

!  ====jules_rivers module JULES====
  ALLOCATE( tot_surf_runoff_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( tot_sub_runoff_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( acc_lake_evap_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( rflow_gb(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( rrun_gb(land_pts), STAT=error )
  error_sum = error_sum + error

#endif


!-----------------------------------------------------------------------
! Write out an error if there was one
!-----------------------------------------------------------------------
! Check for error.
  IF ( error_sum /= 0 )                                                       &
    CALL ereport("allocate_jules_arrays", 101,                                &
                 "Error allocating JULES model arrays")


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE allocate_jules_arrays
END MODULE allocate_jules_arrays_mod
