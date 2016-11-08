MODULE surf_couple_extra_mod
! *****************************COPYRIGHT****************************************
! (c) Crown copyright, Met Office 2014. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms and
! conditions set out therein.
!
! [Met Office Ref SC0237]
! *****************************COPYRIGHT****************************************
IMPLICIT NONE

PRIVATE

PUBLIC :: surf_couple_extra

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SURF_COUPLE_EXTRA_MOD'

CONTAINS

!===============================================================================
! Public subroutine
!===============================================================================
SUBROUTINE surf_couple_extra(                                                 &
#if !defined(UM_JULES)
   ! Arguments used by JULES-standalone
   !Driving data and associated INTENT(IN)
   ls_rain, con_rain, ls_snow, con_snow, tl_1, lw_down, qw_1, u_1, v_1, pstar,&
   !Fluxes INTENT(IN)
   ei_surft, surf_htf_surft, ecan_surft, ext, sw_surft,                       &
   !Misc INTENT(IN)
   a_step, tile_frac, hcons,                                                  &
   !Fluxes INTENT(INOUT)
   melt_surft,                                                                &
   !Fluxes INTENT(OUT)
   snomlt_surf_htf, snowmelt, snomlt_sub_htf, sub_surf_roff, surf_roff,       &
   tot_tfall, snow_melt, rrun, rflow, snow_soil_htf, snow_smb_surft            &
#else
   !Arguments for the UM-----------------------------------------
   !IN
   land_pts, row_length, rows, river_row_length, river_rows, land_index,      &
   ls_rain, con_rain, ls_snow, con_snow, surf_ht_flux_land,                   &
   cca_2d, smlt, nsurft, surft_pts, surft_index, tile_frac,                   &
   ei_surft, surf_htf_surft,                                                  &
   tstar_surft, hcons,                                                        &
   lice_pts, lice_index, soil_pts, soil_index, ext,                           &
   stf_sub_surf_roff,                                                         &
   ecan_surft, fexp_gb, gamtot_gb, ti_mean_gb, ti_sig_gb, cs_ch4,             &
   a_fsat_gb, c_fsat_gb, a_fwet_gb, c_fwet_gb,                                &
   ntype, fqw_surft,                                                          &
   halo_i, halo_j, model_levels,                                              &
   delta_lambda, delta_phi, xx_cos_theta_latitude, i_river_vn,                &
   aocpl_row_length, aocpl_p_rows, xpa, xua, xva, ypa, yua, yva,              &
   g_p_field, g_r_field, n_proc, global_row_length, global_rows,              &
   global_river_row_length, global_river_rows, flandg, river_vel, river_mcoef,&
   trivdir, trivseq, r_area, slope, flowobs1, r_inext, r_jnext, r_land,       &
   substore, surfstore, flowin, bflowin, smvcst, smvcwt,                      &
   a_step, n_rows, offx, offy, n_procx, n_procy, g_rows, g_row_length,        &
   at_extremity, frac_disturb, satcon, soil_clay_ij, resp_s_gb, npp_gb,       &
   z0m_soil_gb,                                                               &
   !INOUT
   a_steps_since_riv, melt_surft, t_soil_gb, tsurf_elev_surft,                &
   rgrain_surft, snow_grnd_surft, snow_surft,                                 &
   soil_layer_moisture, sthf_gb, sthu_gb, canopy_surft, fsat_gb, fwetl_gb,    &
   zw_gb, sthzw_gb,                                                           &
   snow_depth, snowmelt, ls_rainfrac_land,                                    &
   tot_surf_runoff, tot_sub_runoff, acc_lake_evap, twatstor,                  &
   asteps_since_triffid, g_leaf_acc_pft, g_leaf_phen_acc_pft, npp_acc_pft,    &
   resp_s_acc_gb,                                                             &
   resp_w_acc_pft, cs_pool_gb, frac_surft, lai_pft, canht_pft,                &
   catch_snow_surft, catch_surft, infil_surft,                                &
   !OUT
   inlandout_atm_gb, surf_ht_flux_ld, snow_melt,                              &
   snomlt_sub_htf, dhf_surf_minus_soil,                                       &
   canopy_gb, snomlt_surf_htf, smc_gb, sub_surf_roff,                         &
   surf_roff,                                                                 &
   tot_tfall, z0_surft, z0h_bare_surft, snow_soil_htf, snow_smb_surft,        &
   land_sea_mask                                                              &
#endif
   )

!Module imports

!Common modules
!Import interfaces to subroutines called
USE hydrol_mod,               ONLY:                                           &
  hydrol
USE snow_mod,                 ONLY:                                           &
  snow
USE river_control_mod,        ONLY:                                           &
  river_control

!Variables in modules that don't change names between UM and JULES
USE atm_fields_bounds_mod,    ONLY:                                           &
  tdims, tdims_s, pdims, pdims_s
USE jules_soil_mod,           ONLY:                                           &
  sm_levels, l_soil_sat_down, confrac
USE theta_field_sizes,        ONLY:                                           &
  t_i_length, t_j_length
USE jules_vegetation_mod,     ONLY:                                           &
  l_crop, l_triffid, l_trif_eq, l_phenol, phenol_period, triffid_period,      &
  can_model, l_irrig_dmd, irr_crop, l_irrig_limit, l_inferno, ignition_method,&
  l_fao_ref_evapotranspiration
USE trif_vars_mod, ONLY:                                                      &
  cnsrv_carbon_veg2_gb, cnsrv_veg_triffid_gb, cnsrv_soil_triffid_gb,          &
  cnsrv_prod_triffid_gb, cnsrv_carbon_triffid_gb, fao_et0
USE jules_hydrology_mod,      ONLY:                                           &
  l_hydrology, l_pdm, l_top, l_var_rainfrac
USE jules_surface_types_mod,  ONLY:                                           &
  npft, ncpft, nnpft
USE trif_vars_mod,            ONLY: frac_past_gb

! Variables required only in UM-mode
#if defined(UM_JULES)
  USE veg_control_mod,          ONLY:                                         &
     veg_control
  USE diagnostics_riv_mod,      ONLY:                                         &
     diagnostics_riv
  USE prognostics,              ONLY:                                         &
     nsnow_surft, snowdepth_surft, rgrainl_surft, rho_snow_surft,             &
     sice_surft, sliq_surft, tsnow_surft, rho_snow_grnd_surft, ds_surft

  USE atm_fields_real_mod,      ONLY: disturb_veg_prev,                       &
     wood_prod_fast_d1, wood_prod_med_d1, wood_prod_slow_d1

  USE timestep_mod,             ONLY:                                         &
     timestep
  USE jules_radiation_mod,      ONLY:                                         &
     l_snow_albedo
  USE um_parcore,               ONLY:                                         &
     mype

!CABLE_LSM:
USE hydrol_mod_cable,        ONLY : hydrol_cable
USE lsm_switches_mod,        ONLY : lsm_id
USE timestep_mod,            ONLY : timestep_number 

  USE p_s_parms,                ONLY:                                         &
     ! To avoid clash with non _levs versions
     clapp_levs  => bexp_gb,                                                  &
     sathh_levs  => sathh_gb,                                                 &
     hcap_levs   => hcap_gb,                                                  &
     hcon_levs   => hcon_gb,                                                  &
     satcon_levs => satcon_gb,                                                &
     smvccl_levs => smvccl_gb,                                                &
     smvcwt_levs => smvcwt_gb,                                                &
     smvcst_levs => smvcst_gb
  USE atm_step_local,           ONLY:                                         &
     land_pts_trif, npft_trif, STASHwork19, STASHwork8, STASHwork26

!Variables in modules that change name between JULES and UM
  USE atm_step_local,           ONLY:                                         &
     dim_cs1, dim_cs2
  USE river_inputs_mod,         ONLY:                                         &
     l_inland, l_rivers

!Error reporting
  USE umPrintMgr
  USE ereport_mod,              ONLY:                                         &
     ereport

!For detecting the SCM, rather than using an ifdef
  USE model_domain_mod,         ONLY:                                         &
     model_type, mt_single_column

  USE stash_array_mod,          ONLY: sf

! JULES-standalone only
#else

!Modules that change name between JULES and UM
  USE jules_rivers_mod,         ONLY:                                         &
     l_rivers
  USE ancil_info,               ONLY:                                         &
    lice_pts, lice_index, soil_pts, soil_index,                               &
    dim_cs1, frac_surft, land_pts, nsurft, land_index, surft_pts, surft_index,&
    row_length, rows
  USE switches,                 ONLY:                                         &
    l_inland

!JULES-standalone only
  USE diag_swchs,               ONLY:                                         &
    stf_sub_surf_roff,                                                        &
    smlt => stf_hf_snow_melt,                                                 &
    srflow, srrun
  USE trifctl,                  ONLY:                                         &
    asteps_since_triffid, g_leaf_acc_pft, npp_acc_pft, g_leaf_phen_acc_pft,   &
    resp_s_acc_gb, resp_w_acc_pft, g_leaf_dr_out_pft, npp_dr_out_pft,         &
    resp_w_dr_out_pft, resp_s_dr_out_gb, c_veg_pft, cv_gb, lit_c_pft,         &
    lit_c_mn_gb, g_leaf_day_pft, g_leaf_phen_pft, lai_phen_pft, frac_agr_gb,  &
    resp_s_gb, npp_gb
  USE p_s_parms,                ONLY:                                         &
    clapp_levs => bexp_gb,                                                    &
    hcon_levs => hcon_gb,                                                     &
    sathh_levs => sathh_gb,                                                   &
    smvcwt_levs => smvcwt_gb,                                                 &
    catch_snow_surft, sthu_gb, sthf_gb, catch_surft, infil_surft,             &
    z0_surft,clay_gb,                                                                 &
    hcap_levs => hcap_gb,                                                     &
    smvcst_levs => smvcst_gb,                                                 &
    satcon_levs => satcon_gb,                                                 &
    z0h_bare_surft, z0m_soil_gb
  USE model_time_mod,           ONLY:                                         &
    timestep_len
  USE top_pdm,                  ONLY:                                         &
    inlandout_atm_gb, fexp_gb, gamtot_gb, ti_mean_gb, ti_sig_gb, dun_roff_gb, &
    drain_gb, fsat_gb, fwetl_gb, qbase_gb, qbase_zw_gb, zw_gb, sthzw_gb,      &
    a_fsat_gb, c_fsat_gb, a_fwet_gb, c_fwet_gb, fch4_wetl_gb, fch4_wetl_cs_gb,&
    fch4_wetl_npp_gb, fch4_wetl_resps_gb
  USE prognostics,              ONLY:                                         &
    canopy_surft, canopy_gb, smc_gb, tstar_surft, rgrain_surft, rgrainl_surft,&
    rho_snow_grnd_surft, sice_surft, sliq_surft, snow_grnd_surft, snow_surft, &
    tsnow_surft, ds_surft, snow_mass_ij,                                      &
    soil_layer_moisture => smcl_gb,                                           &
    t_soil_gb, snowdepth_surft, nsnow_surft, cs_pool_gb, canht_pft, lai_pft,  &
    rho_snow_surft, tsurf_elev_surft
  USE datetime_mod,             ONLY:                                         &
    secs_in_day

!For science that isn't in the UM at this point (or at all)
  !Crops
  USE zenith_mod,               ONLY:                                         &
    photoperiod
  USE crop_vars_mod,            ONLY:                                         &
    phot, dphotdt, dvi_cpft, rootc_cpft, harvc_cpft, reservec_cpft,           &
    croplai_cpft, cropcanht_cpft, dvimax_gb, frac_irr_surft,                  &
    plant_n_gb, nday_crop, sthu_irr_gb
  USE sf_diags_mod,             ONLY:                                         &
    sf_diag
  USE trifctl,                  ONLY:                                         &
    npp_pft
  USE p_s_parms,                ONLY:                                         &
    smvccl_gb
  USE fire_mod,                 ONLY:                                         &
    fire_prog, fire_diag, l_fire
  USE fire_timestep_mod,        ONLY:                                         &
    fire_timestep
  USE metstats_mod,             ONLY:                                         &
    metstats_prog, metstats_input, l_metstats
  USE metstats_timestep_mod,    ONLY:                                         &
    metstats_timestep
  USE model_time_mod,           ONLY:                                         &
    current_time
  USE jules_rivers_trip_mod,    ONLY:                                         &
    adjust_routestore
! Interactive fires (INFERNO)
  USE inferno_io_mod,           ONLY:                                         &
    inferno_io

  USE fao_evapotranspiration,   ONLY:                                         &
    fao_ref_evapotranspiration  
  USE forcing,                  ONLY:                                         &
    pstar_ij, sw_down_ij, lw_down_ij
  USE fluxes,                   ONLY:                                         &
    surf_ht_flux_ij, tstar_ij
  USE model_interface_mod,      ONLY:                                         &
    tiles_to_gbm

#endif

!Dr Hook
USE parkind1,                 ONLY: jprb, jpim
USE yomhook,                  ONLY: lhook, dr_hook



  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Private subroutine for accessing the JULES science routines that are called
!   after the implicit code
!
! Current Code Owner: Richard Gilham
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

!Subroutine Arguments

!-----------------------------------------------------------------------------
! Integer variables required for UM_JULES definitions
#if defined(UM_JULES)
INTEGER, INTENT(IN) ::                                                        &
  land_pts,                                                                   &
  row_length,                                                                 &
  rows,                                                                       &
  river_row_length,                                                           &
  river_rows,                                                                 &
  land_index(land_pts)
REAL :: clay_gb(land_pts)
#endif

!-----------------------------------------------------------------------------
! Common variables to both UM_JULES and STANDALONE
!Driving data and associated INTENT(IN)
REAL, INTENT(INOUT) ::                                                        &
  ls_rain(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
  con_rain(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
REAL, INTENT(IN) ::                                                           &
  ls_snow(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
  con_snow(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)

!Fluxes INTENT(IN)
REAL, INTENT(IN) ::                                                           &
  ei_surft(land_pts,nsurft),                                                  &
       !Sublimation of snow (kg/m2/s)
  surf_htf_surft(land_pts,nsurft),                                            &
       !Surface heat flux (W/m2)
  ecan_surft(land_pts,nsurft),                                                &
       !Canopy evaporation from land tiles (kg/m2/s).
  ext(land_pts,sm_levels)
       !Extraction of water from each soil layer (kg/m2/s).

!Misc INTENT(IN)
INTEGER, INTENT(IN) ::                                                        &
  a_step
REAL, INTENT(IN)  ::                                                          &
  tile_frac(land_pts,nsurft)

!Fluxes INTENT(INOUT)
REAL, INTENT(INOUT) ::                                                        &
  melt_surft(land_pts,nsurft),                                                &
        !Surface or canopy snowmelt rate (kg/m2/s)
        !On output, this is the total melt rate for the tile
        !(i.e. sum of  melt on canopy and ground).
  snomlt_sub_htf(land_pts)
        !Sub-canopy snowmelt heat flux (W/m2)

REAL, INTENT(INOUT) ::                                                        &
   hcons(land_pts)

#if defined(UM_JULES)
! Snowmelt initialised earlier for UM runs
REAL, INTENT(INOUT) :: snowmelt(row_length,rows)
REAL, INTENT(OUT)   :: dhf_surf_minus_soil(land_pts)
        !heat flux difference across the FLake snowpack (W/m2)
#else
! Snowmelt purely output in stand-alone JULES
REAL, INTENT(  OUT) :: snowmelt(row_length,rows)
REAL :: dhf_surf_minus_soil(land_pts)
        !heat flux difference across the FLake snowpack (W/m2)
#endif


!Fluxes INTENT(OUT)
REAL, INTENT(OUT)::                                                           &
  sub_surf_roff(land_pts),                                                    &
        !Sub-surface runoff (kg/m2/s).
  surf_roff(land_pts),                                                        &
        !Surface runoff (kg/m2/s).
  tot_tfall(land_pts),                                                        &
        !Total throughfall (kg/m2/s).
  snomlt_surf_htf(row_length,rows),                                           &
        !Gridbox snowmelt heat flux (W/m2)
  snow_soil_htf(land_pts,nsurft),                                             &
        !Tiled snow->soil heat flux (W/m2) 
  snow_smb_surft(land_pts,nsurft),                                             &
        !Tiled snow mass balance (kg/m2/s) 
  snow_melt(land_pts)

!Local constants
INTEGER ::                                                                    &
   i,j,l,n,                                                                   &
        !Various counters
   p_field
        !Number of model points

!Local variables
INTEGER ::                                                                    &
   phenol_call,                                                               &
        !indicates whether phenology is to be called
  triffid_call,                                                               &
        !indicates whether TRIFFID is to be called
  crop_call,                                                                  &
        !indicates whether crop model is to be called
  crop_period,                                                                &
        !crops have a daily calling period
  nstep_trif
        !Number of atmospheric timesteps between calls to TRIFFID
        !vegetation model

PARAMETER( crop_period = 1 )
        ! Crop code hard wired to run daily : crop_period = 1

REAL ::                                                                       &
  a_boxareas(row_length,rows),                                                &
  ls_rain_land(land_pts),                                                     &
  con_rain_land(land_pts),                                                    &
  ls_snow_land(land_pts),                                                     &
  con_snow_land(land_pts),                                                    &
  con_rainfrac_land(land_pts),                                                &
  inlandout_atmos(row_length,rows),                                           &
  lying_snow(land_pts)

!-----------------------------------------------------------------------------
!JULES standalone-only arguments
#if !defined(UM_JULES)

REAL, INTENT(IN) ::                                                           &
   tl_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &
   lw_down(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),              &
   qw_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &
   u_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                  &
   v_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                  &
   pstar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                &
   sw_surft(land_pts,nsurft)
         !Surface net SW radiation on land tiles (W/m2)

!Fluxes INTENT(OUT)
REAL, INTENT(OUT)::                                                           &
  rrun(land_pts),                                                             &
         !Surface runoff after river routing (kg/m2/s)
  rflow(land_pts)
         !River runoff (kg/m2/s)

! Local variables - JULES standalone only
REAL ::                                                                       &
   timestep,                                                                  &
         ! Model timestep (s)
  surf_ht_flux_ld(land_pts),                                                  &
         ! Surface heat flux on land (W/m2)
  ls_rainfrac_land(land_pts),                                                 &
  cs_ch4(land_pts),                                                           &
         ! soil carbon used in wetland CH4 emissions model if TRIFFID
         ! is switched off
  trad(land_pts)    
         ! gridbox effective radiative temperature (assuming emissivity=1)   
!-----------------------------------------------------------------------------
!UM-only arguments
#else

REAL, INTENT(IN) ::                                                           &
  surf_ht_flux_land(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),     &
  cca_2d(row_length,rows),                                                    &
  tstar_surft(land_pts,nsurft)

LOGICAL, INTENT(IN) ::                                                        &
  smlt

INTEGER, INTENT(IN) ::                                                        &
  nsurft,                                                                     &
  surft_pts(nsurft),                                                          &
  surft_index(land_pts,nsurft)

!Arguments added for hydrol
INTEGER, INTENT(IN) ::                                                        &
  lice_pts,                                                                   &
  lice_index(land_pts),                                                       &
  soil_pts,                                                                   &
  soil_index(land_pts)

REAL, INTENT(IN) ::                                                           &
  fexp_gb(land_pts),                                                          &
  gamtot_gb(land_pts),                                                        &
  ti_mean_gb(land_pts),                                                       &
  ti_sig_gb(land_pts),                                                        &
  resp_s_gb(land_pts,dim_cs1),                                                &
  npp_gb(land_pts),                                                           &
  a_fsat_gb(land_pts),                                                        &
  c_fsat_gb(land_pts),                                                        &
  a_fwet_gb(land_pts),                                                        &
  c_fwet_gb(land_pts)

LOGICAL, INTENT(IN) ::                                                        &
  stf_sub_surf_roff

!Arguments added for rivers
INTEGER, INTENT(IN) ::                                                        &
  ntype,                                                                      &
  i_river_vn,                                                                 &
  aocpl_row_length,                                                           &
  aocpl_p_rows,                                                               &
  g_p_field,                                                                  &
  g_r_field,                                                                  &
  n_proc,                                                                     &
  global_row_length,                                                          &
  global_rows,                                                                &
  global_river_row_length,                                                    &
  global_river_rows,                                                          &
  halo_i,                                                                     &
  halo_j,                                                                     &
  model_levels

REAL, INTENT(IN) ::                                                           &
  fqw_surft(land_pts,nsurft),                                                 &
  delta_lambda,                                                               &
  delta_phi,                                                                  &
  xx_cos_theta_latitude(tdims_s%i_start:tdims_s%i_end,                        &
                        tdims_s%j_start:tdims_s%j_end),                       &
  xpa(aocpl_row_length+1),                                                    &
  xua(0:aocpl_row_length),                                                    &
  xva(aocpl_row_length+1),                                                    &
  ypa(aocpl_p_rows),                                                          &
  yua(aocpl_p_rows),                                                          &
  yva(0:aocpl_p_rows),                                                        &
  flandg(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),        &
  river_vel,                                                                  &
  river_mcoef,                                                                &
  trivdir(river_row_length, river_rows),                                      &
  trivseq(river_row_length, river_rows),                                      &
  r_area(row_length, rows),                                                   &
  slope(row_length, rows),                                                    &
  flowobs1(row_length, rows),                                                 &
  r_inext(row_length, rows),                                                  &
  r_jnext(row_length, rows),                                                  &
  r_land(row_length, rows),                                                   &
  substore(row_length, rows),                                                 &
  surfstore(row_length, rows),                                                &
  flowin(row_length, rows),                                                   &
  bflowin(row_length, rows),                                                  &
  smvcst(land_pts),                                                           &
  smvcwt(land_pts)

INTEGER, INTENT(IN) ::                                                        &
  n_rows,                                                                     &
  offx,                                                                       &
  offy,                                                                       &
  n_procx,                                                                    &
  n_procy,                                                                    &
  g_rows (0:n_proc-1),                                                        &
  g_row_length (0:n_proc-1)

REAL, INTENT(IN) ::                                                           &
  frac_disturb(land_pts),                                                     &
  satcon(land_pts),                                                           &
  soil_clay_ij(row_length,rows),                                              &
  z0m_soil_gb(land_pts)

LOGICAL, INTENT(IN) ::                                                        &
  at_extremity(4)

INTEGER, INTENT(INOUT)  ::                                                    &
  a_steps_since_riv

REAL, INTENT(INOUT) ::                                                        &
  snow_depth(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),            &
  ls_rainfrac_land(land_pts),                                                 &
  catch_snow_surft(land_pts,nsurft),                                          &
  t_soil_gb(land_pts,sm_levels),                                              &
  tsurf_elev_surft(land_pts,nsurft),                                          &
  rgrain_surft(land_pts,nsurft),                                              &
  snow_grnd_surft(land_pts,nsurft),                                           &
  snow_surft(land_pts,nsurft),                                                &
  soil_layer_moisture(land_pts,sm_levels),                                    &
  sthf_gb(land_pts,sm_levels)

!Arguments added for hydrol
REAL, INTENT(INOUT) ::                                                        &
  infil_surft(land_pts,nsurft),                                               &
  catch_surft(land_pts,nsurft),                                               &
  sthu_gb(land_pts,sm_levels),                                                &
  canopy_surft(land_pts,nsurft),                                              &
  fsat_gb(land_pts),                                                          &
  fwetl_gb(land_pts),                                                         &
  zw_gb(land_pts),                                                            &
  sthzw_gb(land_pts)

!Arguments added for rivers
REAL, INTENT(INOUT) ::                                                        &
  tot_surf_runoff(land_pts),                                                  &
  tot_sub_runoff(land_pts),                                                   &
  acc_lake_evap(row_length,rows),                                             &
  twatstor(river_row_length, river_rows)

!Arguments added for veg
INTEGER, INTENT(INOUT) ::                                                     &
  asteps_since_triffid

REAL, INTENT(INOUT) ::                                                        &
  g_leaf_acc_pft(land_pts,npft),                                              &
  g_leaf_phen_acc_pft(land_pts,npft),                                         &
  npp_acc_pft(land_pts_trif,npft_trif),                                       &
  resp_s_acc_gb(land_pts_trif,dim_cs1),                                       &
  resp_w_acc_pft(land_pts_trif,npft_trif),                                    &
  cs_pool_gb(land_pts,dim_cs1),                                               &
  frac_surft(land_pts,ntype),                                                 &
  lai_pft(land_pts,npft),                                                     &
  canht_pft(land_pts,npft)

REAL, INTENT(OUT) ::                                                          &
  inlandout_atm_gb(land_pts),                                                 &
  surf_ht_flux_ld(land_pts)

REAL, INTENT(OUT) ::                                                          &
  canopy_gb(land_pts),                                                        &
  smc_gb(land_pts)

!Arguments added for veg
REAL, INTENT(OUT) ::                                                          &
  z0_surft(land_pts,nsurft),                                                  &
  z0h_bare_surft(land_pts,nsurft)

LOGICAL, INTENT(IN) ::                                                        &
  land_sea_mask(row_length, rows)

!Local variables
REAL ::                                                                       &
  riverout(row_length, rows),                                                 &
  box_outflow(river_row_length, river_rows),                                  &
  box_inflow(river_row_length, river_rows),                                   &
  inlandout_riv(river_row_length,river_rows),                                 &
  infil(land_pts),                                                            &
  dun_roff_gb(land_pts),                                                      &
  qbase_gb(land_pts),                                                         &
  qbase_zw_gb(land_pts),                                                      &
  drain_gb(land_pts),                                                         &
  cs_ch4(land_pts),                                                           &
  fch4_wetl_gb(land_pts),                                                     &
  fch4_wetl_cs_gb(land_pts),                                                  &
  fch4_wetl_npp_gb(land_pts),                                                 &
  fch4_wetl_resps_gb(land_pts)

INTEGER ::                                                                    &
  nstep_trip,                                                                 &
  gather_pe_trip

LOGICAL ::                                                                    &
  first_routing,                                                              &
  invert_atmos,                                                               &
  trip_call

!Local variables to pass from science routines to STASH routines
REAL ::                                                                       &
  c_veg_pft(land_pts,npft),                                                   &
  cv_gb(land_pts),                                                            &
  lit_c_pft(land_pts,npft),                                                   &
  lit_c_mn_gb(land_pts),                                                      &
  g_leaf_day_pft(land_pts,npft),                                              &
  g_leaf_phen_pft(land_pts,npft),                                             &
  lai_phen_pft(land_pts,npft),                                                &
  g_leaf_dr_out_pft(land_pts,npft),                                           &
  npp_dr_out_pft(land_pts,npft),                                              &
  resp_w_dr_out_pft(land_pts,npft),                                           &
  resp_s_dr_out_gb(land_pts,dim_cs1+1)

! Need INVERT_OCEAN (hence N->S configuration of ATMOS grid)
! for river routing
LOGICAL, PARAMETER :: invert_ocean=.FALSE.

#endif

!-----------------------------------------------------------------------------

CHARACTER(LEN=256)            :: message
INTEGER                       :: errorstatus
CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'SURF_COUPLE_EXTRA'

!Dr Hook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------------
!End of header
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

#if !defined(UM_JULES)
!Set up a some commonly used values
p_field  = t_i_length * t_j_length
timestep = REAL(timestep_len)

!initialise
trad(:) = 0.0
#endif

  SELECT CASE( lsm_id )
    CASE( 'jules' )
! ----------------------------------------------------------------------
! Section HYD.1 Compress fields to land points, then call hydrology.
! ----------------------------------------------------------------------

IF (l_hydrology .AND. land_pts /= 0 ) THEN
! Inland basin outflow is added to soil moisture at each timestep.
! This flux changes its value only when the river routing scheme
! has been called in the previous timestep

#if defined(UM_JULES)
  IF(l_rivers)THEN
    !Pass inland flow to soil moisture every timestep
    IF  (.NOT. l_inland) THEN
      inlandout_atmos = 0.0
      inlandout_riv   = 0.0
    END IF
  END IF  ! l_rivers
#endif

  !Compress fields to land points
  DO l = 1, land_pts
    j = (land_index(l)-1)/row_length + 1
    i = land_index(l) - (j-1)*row_length
    ls_rain_land(l)    = ls_rain(i,j)
    con_rain_land(l)   = con_rain(i,j)
    con_snow_land(l)   = con_snow(i,j)
    ls_snow_land(l)    = ls_snow(i,j)

#if defined(UM_JULES)
    surf_ht_flux_ld(l) = surf_ht_flux_land(i,j)
    lying_snow(l)      = snow_depth(i,j)
    snow_melt(l)       = snowmelt(i,j)
#endif

    !pass jules the modelled rain fractions
    IF (l_var_rainfrac) THEN

#if defined(UM_JULES)
      con_rainfrac_land(l) = MIN(cca_2d(i,j),0.5)  !As in diagnostics_conv
#endif
      !provide some safety checking for convective rain with no CCA
      IF (con_rainfrac_land(l) == 0.0 .AND. con_rain_land(l) > 0.0) THEN
        con_rainfrac_land(l) = confrac
      !and for very small CCA amounts
      ELSE IF (con_rainfrac_land(l) < 0.01 .AND. con_rain_land(l) > 0.0) THEN
        con_rainfrac_land(l) = 0.01
      END IF

      !provide some safety checking for ls rain with no rainfrac
      IF (ls_rainfrac_land(l) == 0.0 .AND. ls_rain_land(l) > 0.0) THEN
        ls_rainfrac_land(l) = 0.5
      !and for very small rainfrac amounts
      ELSE IF (ls_rainfrac_land(l) < 0.01 .AND. ls_rain_land(l) > 0.0) THEN
        ls_rainfrac_land(l) = 0.01
      END IF

    ELSE
      !use original default values
      con_rainfrac_land(l) = confrac
      ls_rainfrac_land(l) = 1.0
    END IF
  END DO !land_pts

!-------------------------------------------------------------------------------
!   Snow processes.
!-------------------------------------------------------------------------------
  CALL snow ( land_pts,timestep,smlt,nsurft,surft_pts,                        &
              surft_index,catch_snow_surft,con_snow_land,con_rain_land,       &
              tile_frac,ls_snow_land,ls_rain_land,                            &
              ei_surft,hcap_levs(:,1),hcons,melt_surft,                       &
              soil_layer_moisture(:,1),sthf_gb(:,1),surf_htf_surft,           &
              t_soil_gb(:,1),tsurf_elev_surft,                                &
              tstar_surft,smvcst_levs(:,1),rgrain_surft,rgrainl_surft,        &
              rho_snow_grnd_surft,                                            &
              sice_surft,sliq_surft,snow_grnd_surft,snow_surft,               &
              snowdepth_surft, tsnow_surft,nsnow_surft,ds_surft,              &
              snomlt_surf_htf,lying_snow,rho_snow_surft,snomlt_sub_htf,       &
              snow_melt,snow_soil_htf,surf_ht_flux_ld,snow_smb_surft,         &
              dhf_surf_minus_soil )

!-------------------------------------------------------------------------------
!   Land hydrology.
!-------------------------------------------------------------------------------

! Calculate soil carbon for use in the wetland CH4 scheme only
! (only used if TRIFFID is switched off):
#if !defined(UM_JULES)
  DO j=1,soil_pts
    i=soil_index(j)
      cs_ch4(i)=cs_pool_gb(i,1)
  ENDDO
#endif

  CALL hydrol (                                                               &
               lice_pts,lice_index,soil_pts,soil_index,nsnow_surft,           &
               land_pts,sm_levels,clapp_levs,catch_surft,con_rain_land,       &
               ecan_surft,ext,hcap_levs,hcon_levs,ls_rain_land,               &
               con_rainfrac_land, ls_rainfrac_land,                           &
               satcon_levs,sathh_levs,snowdepth_surft,snow_soil_htf,          &
               surf_ht_flux_ld,timestep,                                      &
               smvcst_levs,smvcwt_levs,canopy_surft,                          &
               stf_sub_surf_roff,soil_layer_moisture,sthf_gb,sthu_gb,         &
               t_soil_gb,tsurf_elev_surft,canopy_gb,smc_gb,snow_melt,         &
               sub_surf_roff,surf_roff,tot_tfall,                             &
               inlandout_atm_gb,l_inland,nsurft,surft_pts,surft_index,        &
               infil_surft,melt_surft,tile_frac,                              &
               l_top,l_pdm,fexp_gb,gamtot_gb,ti_mean_gb,ti_sig_gb,            &
               cs_ch4,cs_pool_gb,                                             &
               dun_roff_gb,drain_gb,fsat_gb,fwetl_gb,qbase_gb,qbase_zw_gb,    &
               zw_gb,sthzw_gb,a_fsat_gb,c_fsat_gb,a_fwet_gb,c_fwet_gb,        &
               resp_s_gb,npp_gb,fch4_wetl_gb,                                 &
               fch4_wetl_cs_gb,fch4_wetl_npp_gb,fch4_wetl_resps_gb,           &
               dim_cs1,l_soil_sat_down,l_triffid )

!-------------------------------------------------------------------------------
!   Reset snowmelt over land points.
!-------------------------------------------------------------------------------
  !Copy land points output back to full fields array.
  DO l = 1, land_pts
    j=(land_index(l)-1)/row_length + 1
    i=land_index(l) - (j-1)*row_length
#if defined(UM_JULES)
    snow_depth(i,j) = lying_snow(l)
#else
    snow_mass_ij(i,j) = lying_snow(l)
#endif
    snowmelt(i,j) = snow_melt(l)
  END DO

END IF ! ( l_hydrology .AND. land_pts /= 0 )

!-------------------------------------------------------------------
! RIVER ROUTING
!-------------------------------------------------------------------
IF ( l_rivers ) THEN
#if defined(UM_JULES)
  CALL river_control(                                                         &
    !LOGICAL, INTENT(IN)
    invert_ocean,                                                             &
    !INTEGER, INTENT(IN)
    n_proc, land_pts, row_length, rows, river_row_length, river_rows,         &
    land_index, ntype, i_river_vn, aocpl_row_length, aocpl_p_rows, g_p_field, &
    g_r_field, mype, global_row_length, global_rows, global_river_row_length, &
    global_river_rows, halo_i, halo_j, model_levels, nsurft,                  &
    !REAL, INTENT(IN)
    fqw_surft, delta_lambda, delta_phi, xx_cos_theta_latitude,                &
    xpa, xua, xva, ypa, yua, yva, flandg, river_vel, river_mcoef, trivdir,    &
    trivseq, r_area, slope, flowobs1, r_inext, r_jnext, r_land, substore,     &
    surfstore, flowin, bflowin, smvcst_levs, smvcwt_levs,                     &
    surf_roff, sub_surf_roff, frac_surft,                                     &
    !INTEGER, INTENT(INOUT)
    a_steps_since_riv,                                                        &
    !REAL, INTENT(INOUT)
    tot_surf_runoff, tot_sub_runoff, acc_lake_evap, twatstor,                 &
    soil_layer_moisture, sthu_gb,                                             &
    !LOGICAL, INTENT(OUT)
    trip_call,                                                                &
    !REAL, INTENT(OUT)
    inlandout_atm_gb, inlandout_atmos, inlandout_riv, riverout, box_outflow,  &
    box_inflow                                                                &
    )
#else
  CALL river_control( land_pts,sub_surf_roff                                  &
                             ,surf_roff,srflow,srrun,rflow,rrun)
#endif
END IF ! l_rivers (ATMOS)

#if !defined(UM_JULES)
!-------------------------------------------------------------------------------
!   Calculate irrigation demand and planting dates if required
!-------------------------------------------------------------------------------
IF( l_irrig_dmd ) THEN

! calculate maximum dvi per grid cell
  IF ( irr_crop == 2 ) THEN
    DO l=1,land_pts
      dvimax_gb(l) = MAXVAL(dvi_cpft(l,:))
    END DO
  END IF

  CALL irrig_dmd(land_pts, sm_levels, frac_irr_surft, a_step, plant_n_gb,     &
                 sthf_gb, smvccl_gb, smvcst_levs, smvcwt_levs, sthzw_gb,      &
                 sthu_irr_gb, sthu_gb,                                        &
                 soil_layer_moisture, irr_crop, dvimax_gb)

  IF( irr_crop == 1 ) THEN
    CALL calc_crop_date(land_index, land_pts, t_i_length, t_j_length, nsurft, &
                        frac_surft, sw_surft, tstar_surft, lw_down, tl_1,     &
                        con_rain, ls_rain, con_snow, ls_snow,                 &
                        plant_n_gb, nday_crop)
  END IF

  IF ( l_irrig_limit ) THEN
    CALL adjust_routestore()
  END IF
ELSE
! if .not. l_irrig_dmd, set sthu_irr_gb to 0.0 in case it is still reported
  sthu_irr_gb(:,:) = 0.0
END IF ! l_irrig_dmd


!-------------------------------------------------------------------------------
! Run crop code if required
!-------------------------------------------------------------------------------
IF ( l_crop ) THEN
  crop_call = MOD ( FLOAT(a_step),                                            &
                    REAL(crop_period) * REAL(secs_in_day) / timestep )

  DO n=1,ncpft
    DO l=1,land_pts
      npp_acc_pft(l,nnpft+n) = npp_acc_pft(l,nnpft+n)                         &
                            + (npp_pft(l,nnpft+n) * timestep)
    END DO
  END DO

  CALL photoperiod(p_field, phot, dphotdt)

  CALL crop(p_field, land_pts, land_index, a_step,                            &
            crop_call, sm_levels, frac_surft, phot, dphotdt,                  &
            sf_diag%t1p5m_surft, t_soil_gb, sthu_gb, satcon_levs, smvccl_gb,  &
            smvcst_levs, npp_acc_pft,                                         &
            canht_pft, lai_pft, dvi_cpft, rootc_cpft, harvc_cpft,             &
            reservec_cpft, croplai_cpft, cropcanht_cpft,                      &
            catch_surft, infil_surft, z0_surft)
END IF  ! l_crop

!------------------------------------------------------------------------------
!   Update metstats for this timestep
!-------------------------------------------------------------------------------
IF ( l_metstats ) THEN
  !Compress variables to land points as metstats has no knowledge of
  !i and j. Also use a TYPE to keep the argument list short
  DO l = 1, land_pts
    j = ( land_index(l)-1 ) / t_i_length + 1
    i = land_index(l) - (j-1) * t_i_length
    metstats_input(l)%temp     = tl_1(i,j)
    metstats_input(l)%spec_hum = qw_1(i,j)
    metstats_input(l)%wind_u   = u_1(i,j)
    metstats_input(l)%wind_v   = v_1(i,j)
    metstats_input(l)%ls_rain  = ls_rain(i,j)
    metstats_input(l)%con_rain = con_rain(i,j)
    metstats_input(l)%ls_snow  = ls_snow(i,j)
    metstats_input(l)%con_snow = con_snow(i,j)
    metstats_input(l)%press    = pstar(i,j)
  END DO

  CALL metstats_timestep(metstats_input,metstats_prog,                        &
       !Things that really ought to come in via USE but won't work with the UM
                         current_time%time, timestep,land_pts)
END IF

!-------------------------------------------------------------------------------
!   Call to fire module
!-------------------------------------------------------------------------------
IF ( l_fire ) THEN
  CALL fire_timestep(metstats_prog, smc_gb, fire_prog, fire_diag,             &
       !Things that really ought to come in via USE but won't work with the UM
                     current_time%time, current_time%month, timestep, land_pts)
END IF
#endif

#if !defined(UM_JULES)
!-------------------------------------------------------------------------------
!   Call to INFERNO (interactive fire module)
!-------------------------------------------------------------------------------
IF ( l_inferno ) THEN
   CALL inferno_io( sf_diag%t1p5m_surft, sf_diag%q1p5m_surft, pstar, sthu_gb, &
                    sm_levels,                                                &
                    frac_surft, dim_cs1, cs_ch4, canht_pft, ls_rain, con_rain,&
                    land_pts, ignition_method)
END IF
#endif

! ----------------------------------------------------------------------
! Section 19 -- VEGETATION DYNAMICS
! ----------------------------------------------------------------------

! initialize carbon conservation diagnostics
! otherwise they can be non-zero on non-triffid timesteps
IF ( l_triffid ) THEN
  DO l = 1, land_pts
    cnsrv_carbon_veg2_gb(l) = 0.0
    cnsrv_carbon_triffid_gb(l) = 0.0
    cnsrv_veg_triffid_gb(l) = 0.0
    cnsrv_soil_triffid_gb(l) = 0.0
    cnsrv_prod_triffid_gb(l) = 0.0
  END DO
ENDIF

#if defined(UM_JULES)
!change 2d to 1d soil clay content for soil respiration
IF ( l_triffid ) THEN
  DO l = 1, land_pts
    j = (land_index(l)-1)/row_length + 1
    i = land_index(l) - (j-1)*row_length
    clay_gb(l) = soil_clay_ij(i,j)
  ENDDO
ENDIF
#endif

!-------------------------------------------------------------------------
!   If leaf phenology is activated, check whether the atmosphere model
!   has run an integer number of phenology calling periods.
!-------------------------------------------------------------------------
#if defined(UM_JULES)

IF (l_phenol .OR. l_triffid) THEN
  CALL veg_control(                                                           &
    land_pts, land_index, nsurft, can_model,                                  &
    a_step, asteps_since_triffid,                                             &
    land_pts_trif, npft_trif,                                                 &
    phenol_period, triffid_period, row_length, rows,                          &
    l_phenol, l_triffid, l_trif_eq,                                           &
    timestep, frac_disturb, frac_past_gb, satcon_levs,                        &
    g_leaf_acc_pft, g_leaf_phen_acc_pft, npp_acc_pft,                         &
    resp_s_acc_gb, resp_w_acc_pft,                                            &
    cs_pool_gb, frac_surft, lai_pft, clay_gb, z0m_soil_gb, canht_pft,    &
    catch_snow_surft, catch_surft, infil_surft, z0_surft, z0h_bare_surft,     &
    c_veg_pft, cv_gb, lit_c_pft, lit_c_mn_gb, g_leaf_day_pft, g_leaf_phen_pft,&
    lai_phen_pft, g_leaf_dr_out_pft, npp_dr_out_pft, resp_w_dr_out_pft,       &
    resp_s_dr_out_gb                                                          &
     )
END IF
#else
phenol_call=1
triffid_call=1
IF ( l_phenol ) phenol_call = MOD ( FLOAT(a_step),                            &
            REAL(phenol_period) * REAL(secs_in_day) / timestep )

IF ( l_triffid ) THEN
  nstep_trif = INT( REAL(secs_in_day) * REAL(triffid_period) / timestep )
  IF ( asteps_since_triffid == nstep_trif ) triffid_call = 0
ENDIF

IF ( triffid_call == 0 ) THEN
!-------------------------------------------------------------------------------
!     Run includes dynamic vegetation
!-------------------------------------------------------------------------------
  CALL veg2( land_pts, land_index, nsurft, can_model                          &
            ,a_step, asteps_since_triffid                                     &
            ,phenol_period, triffid_period, l_phenol, l_triffid, l_trif_eq    &
            ,timestep, frac_agr_gb, frac_past_gb, satcon_levs                 &
            ,g_leaf_acc_pft, g_leaf_phen_acc_pft, npp_acc_pft                 &
            ,resp_s_acc_gb, resp_w_acc_pft                                    &
            ,cs_pool_gb, frac_surft, lai_pft, clay_gb, z0m_soil_gb            &
            ,canht_pft, catch_snow_surft, catch_surft, infil_surft, z0_surft  &
            ,z0h_bare_surft, c_veg_pft, cv_gb, lit_c_pft, lit_c_mn_gb         &
            ,g_leaf_day_pft, g_leaf_phen_pft, lai_phen_pft, g_leaf_dr_out_pft &
            ,npp_dr_out_pft, resp_w_dr_out_pft, resp_s_dr_out_gb )

ELSE

  IF ( phenol_call == 0 )                                                     &
!-------------------------------------------------------------------------------
!       Run includes phenology,  but not dynamic vegetation
!       therefore call veg1 rather than veg2
!-------------------------------------------------------------------------------
    CALL veg1( land_pts, nsurft, can_model, a_step, phenol_period, l_phenol   &
              ,timestep, satcon_levs, z0m_soil_gb, g_leaf_acc_pft, frac_surft &
              ,lai_pft,  canht_pft                                            &
              ,catch_snow_surft, catch_surft, infil_surft, z0_surft           &
              ,z0h_bare_surft                                                 &
              ,g_leaf_day_pft, g_leaf_phen_pft, g_leaf_phen_acc_pft           &
              ,lai_phen_pft )

ENDIF  !  triffid_call

IF (l_fao_ref_evapotranspiration) THEN
   trad = ( tiles_to_gbm(tstar_surft**4) )**0.25
   CALL fao_ref_evapotranspiration(soil_pts, soil_index,                      &
     land_pts, land_index, sf_diag%t1p5m,                                     &
     sw_down_ij, lw_down_ij, surf_ht_flux_ij, sf_diag%u10m,                   &
     sf_diag%v10m, sf_diag%q1p5m, pstar_ij, trad, fao_et0)
END IF

#endif

!------------------------------------------------------------------------------
!Call STASH diagnostic routines - UM-only
!------------------------------------------------------------------------------
#if defined(UM_JULES)
IF (model_type /= mt_single_column) THEN
  IF (l_hydrology .AND. sf(0,8) ) THEN
!DEPENDS ON: diagnostics_hyd
    CALL diagnostics_hyd(                                                     &
      row_length, rows, model_levels,                                         &
      n_rows, global_row_length, global_rows,                                 &
      halo_i, halo_j, offx, offy, mype,                                       &
      n_proc, n_procx, n_procy,                                               &
      g_rows, g_row_length,                                                   &
      at_extremity,                                                           &
      land_pts, sm_levels,                                                    &
      !Put inland basin outflow in call to diagriv
      land_index,inlandout_atm_gb,                                            &
      smc_gb, surf_roff, sub_surf_roff,                                       &
      lying_snow, snow_melt,                                                  &
      canopy_gb,t_soil_gb,                                                    &
      tsurf_elev_surft,snow_soil_htf,snow_smb_surft,                          &
      soil_layer_moisture,                                                    &
      nsurft, snomlt_surf_htf, sthu_gb, sthf_gb,                              &
      tot_tfall, snow_surft, melt_surft,                                      &
      rgrain_surft, land_sea_mask,                                            &
      dun_roff_gb, drain_gb, qbase_gb, qbase_zw_gb,                           &
      fch4_wetl_gb,fch4_wetl_cs_gb,fch4_wetl_npp_gb,                          &
      fch4_wetl_resps_gb,                                                     &
      fexp_gb,gamtot_gb,ti_mean_gb,ti_sig_gb,                                 &
      fsat_gb,fwetl_gb,zw_gb,sthzw_gb,                                        &
      timestep,                                                               &
      STASHwork8                                                              &
      )
  END IF

  IF ( l_rivers .AND. trip_call .AND. sf(0,26) ) THEN

    CALL diagnostics_riv(                                                     &
      row_length, rows,                                                       &
      river_row_length, river_rows,                                           &
      at_extremity,                                                           &
      at_extremity,                                                           &
      riverout,                                                               &
      box_outflow, box_inflow,                                                &
      !Put inland basin outflow in call to diagriv
      twatstor,inlandout_riv,                                                 &
      STASHwork26                                                             &
      )
  END IF

  IF (sf(0,19)) THEN
! DEPENDS ON: diagnostics_veg
    CALL diagnostics_veg(                                                     &
      row_length, rows, n_rows,                                               &
      global_row_length, global_rows,                                         &
      dim_cs1, dim_cs2,                                                       &
      halo_i, halo_j, offx, offy, mype,                                       &
      n_proc, n_procx, n_procy,                                               &
      g_rows, g_row_length,                                                   &
      at_extremity,                                                           &
      land_pts,                                                               &
      land_index,                                                             &
      ntype,npft,                                                             &
      c_veg_pft,cv_gb,g_leaf_phen_pft,                                        &
      lit_c_pft,lit_c_mn_gb,g_leaf_day_pft,                                   &
      lai_phen_pft,g_leaf_dr_out_pft,npp_dr_out_pft,                          &
      resp_w_dr_out_pft,resp_s_dr_out_gb,frac_disturb,disturb_veg_prev,       &
      wood_prod_fast_d1, wood_prod_med_d1, wood_prod_slow_d1,                 &
      frac_surft,lai_pft,canht_pft,cs_pool_gb,                                &
      STASHwork19                                                             &
      )
  END IF
END IF
#endif
     CASE( 'cable' )
! ----------------------------------------------------------------------
! Section HYD.1 Compress fields to land points, then call hydrology.
! ----------------------------------------------------------------------

IF (l_hydrology .AND. land_pts /= 0 ) THEN
! Inland basin outflow is added to soil moisture at each timestep.
! This flux changes its value only when the river routing scheme
! has been called in the previous timestep

#if defined(UM_JULES)
  IF(l_rivers)THEN
    !Pass inland flow to soil moisture every timestep
    IF  (.NOT. l_inland) THEN
      inlandout_atmos = 0.0
      inlandout_riv   = 0.0
    END IF
  END IF  ! l_rivers
#endif

  !Compress fields to land points
  DO l = 1, land_pts
    j = (land_index(l)-1)/row_length + 1
    i = land_index(l) - (j-1)*row_length
    ls_rain_land(l)    = ls_rain(i,j)
    con_rain_land(l)   = con_rain(i,j)
    con_snow_land(l)   = con_snow(i,j)
    ls_snow_land(l)    = ls_snow(i,j)

#if defined(UM_JULES)
    surf_ht_flux_ld(l) = surf_ht_flux_land(i,j)
    lying_snow(l)      = snow_depth(i,j)
    snow_melt(l)       = snowmelt(i,j)
#endif

    !pass jules the modelled rain fractions
    IF (l_var_rainfrac) THEN

#if defined(UM_JULES)
      con_rainfrac_land(l) = MIN(cca_2d(i,j),0.5)  !As in diagnostics_conv
#endif
      !provide some safety checking for convective rain with no CCA
      IF (con_rainfrac_land(l) == 0.0 .AND. con_rain_land(l) > 0.0) THEN
        con_rainfrac_land(l) = confrac
      !and for very small CCA amounts
      ELSE IF (con_rainfrac_land(l) < 0.01 .AND. con_rain_land(l) > 0.0) THEN
        con_rainfrac_land(l) = 0.01
      END IF

      !provide some safety checking for ls rain with no rainfrac
      IF (ls_rainfrac_land(l) == 0.0 .AND. ls_rain_land(l) > 0.0) THEN
        ls_rainfrac_land(l) = 0.5
      !and for very small rainfrac amounts
      ELSE IF (ls_rainfrac_land(l) < 0.01 .AND. ls_rain_land(l) > 0.0) THEN
        ls_rainfrac_land(l) = 0.01
      END IF

    ELSE
      !use original default values
      con_rainfrac_land(l) = confrac
      ls_rainfrac_land(l) = 1.0
    END IF
  END DO !land_pts

!-------------------------------------------------------------------------------
!   Snow processes.
!-------------------------------------------------------------------------------
  CALL snow ( land_pts,timestep,smlt,nsurft,surft_pts,                        &
              surft_index,catch_snow_surft,con_snow_land,con_rain_land,       &
              tile_frac,ls_snow_land,ls_rain_land,                            &
              ei_surft,hcap_levs(:,1),hcons,melt_surft,                       &
              soil_layer_moisture(:,1),sthf_gb(:,1),surf_htf_surft,           &
              t_soil_gb(:,1),tsurf_elev_surft,                                &
              tstar_surft,smvcst_levs(:,1),rgrain_surft,rgrainl_surft,        &
              rho_snow_grnd_surft,                                            &
              sice_surft,sliq_surft,snow_grnd_surft,snow_surft,               &
              snowdepth_surft, tsnow_surft,nsnow_surft,ds_surft,              &
              snomlt_surf_htf,lying_snow,rho_snow_surft,snomlt_sub_htf,       &
              snow_melt,snow_soil_htf,surf_ht_flux_ld,snow_smb_surft,         &
              dhf_surf_minus_soil )

!-------------------------------------------------------------------------------
!   Land hydrology.
!-------------------------------------------------------------------------------

! Calculate soil carbon for use in the wetland CH4 scheme only
! (only used if TRIFFID is switched off):
#if !defined(UM_JULES)
  DO j=1,soil_pts
    i=soil_index(j)
      cs_ch4(i)=cs_pool_gb(i,1)
  ENDDO
#endif

  CALL hydrol_cable (                                                         &
               lice_pts,lice_index,soil_pts,soil_index,nsnow_surft,           &
               land_pts,sm_levels,clapp_levs,catch_surft,con_rain_land,       &
               ecan_surft,ext,hcap_levs,hcon_levs,ls_rain_land,               &
               con_rainfrac_land, ls_rainfrac_land,                           &
               satcon_levs,sathh_levs,snowdepth_surft,snow_soil_htf,          &
               surf_ht_flux_ld,timestep,                                      &
               smvcst_levs,smvcwt_levs,canopy_surft,                          &
               stf_sub_surf_roff,soil_layer_moisture,sthf_gb,sthu_gb,         &
               t_soil_gb,tsurf_elev_surft,canopy_gb,smc_gb,snow_melt,         &
               sub_surf_roff,surf_roff,tot_tfall,                             &
               inlandout_atm_gb,l_inland,nsurft,surft_pts,surft_index,        &
               infil_surft,melt_surft,tile_frac,                              &
               l_top,l_pdm,fexp_gb,gamtot_gb,ti_mean_gb,ti_sig_gb,            &
               cs_ch4,cs_pool_gb,                                             &
               dun_roff_gb,drain_gb,fsat_gb,fwetl_gb,qbase_gb,qbase_zw_gb,    &
               zw_gb,sthzw_gb,a_fsat_gb,c_fsat_gb,a_fwet_gb,c_fwet_gb,        &
               resp_s_gb,npp_gb,fch4_wetl_gb,                                 &
               fch4_wetl_cs_gb,fch4_wetl_npp_gb,fch4_wetl_resps_gb,           &
               dim_cs1,l_soil_sat_down,l_triffid,                             &
               mype, timestep_number  & !CABLE_LSM:                           & 
              )  

!-------------------------------------------------------------------------------
!   Reset snowmelt over land points.
!-------------------------------------------------------------------------------
  !Copy land points output back to full fields array.
  DO l = 1, land_pts
    j=(land_index(l)-1)/row_length + 1
    i=land_index(l) - (j-1)*row_length
#if defined(UM_JULES)
    snow_depth(i,j) = lying_snow(l)
#else
    snow_mass_ij(i,j) = lying_snow(l)
#endif
    snowmelt(i,j) = snow_melt(l)
  END DO

END IF ! ( l_hydrology .AND. land_pts /= 0 )

!-------------------------------------------------------------------
! RIVER ROUTING
!-------------------------------------------------------------------
IF ( l_rivers ) THEN
#if defined(UM_JULES)
  CALL river_control(                                                         &
    !LOGICAL, INTENT(IN)
    invert_ocean,                                                             &
    !INTEGER, INTENT(IN)
    n_proc, land_pts, row_length, rows, river_row_length, river_rows,         &
    land_index, ntype, i_river_vn, aocpl_row_length, aocpl_p_rows, g_p_field, &
    g_r_field, mype, global_row_length, global_rows, global_river_row_length, &
    global_river_rows, halo_i, halo_j, model_levels, nsurft,                  &
    !REAL, INTENT(IN)
    fqw_surft, delta_lambda, delta_phi, xx_cos_theta_latitude,                &
    xpa, xua, xva, ypa, yua, yva, flandg, river_vel, river_mcoef, trivdir,    &
    trivseq, r_area, slope, flowobs1, r_inext, r_jnext, r_land, substore,     &
    surfstore, flowin, bflowin, smvcst_levs, smvcwt_levs,                     &
    surf_roff, sub_surf_roff, frac_surft,                                     &
    !INTEGER, INTENT(INOUT)
    a_steps_since_riv,                                                        &
    !REAL, INTENT(INOUT)
    tot_surf_runoff, tot_sub_runoff, acc_lake_evap, twatstor,                 &
    soil_layer_moisture, sthu_gb,                                             &
    !LOGICAL, INTENT(OUT)
    trip_call,                                                                &
    !REAL, INTENT(OUT)
    inlandout_atm_gb, inlandout_atmos, inlandout_riv, riverout, box_outflow,  &
    box_inflow                                                                &
    )
#else
  CALL river_control( land_pts,sub_surf_roff                                  &
                             ,surf_roff,srflow,srrun,rflow,rrun)
#endif
END IF ! l_rivers (ATMOS)

#if !defined(UM_JULES)
!-------------------------------------------------------------------------------
!   Calculate irrigation demand and planting dates if required
!-------------------------------------------------------------------------------
IF( l_irrig_dmd ) THEN

! calculate maximum dvi per grid cell
  IF ( irr_crop == 2 ) THEN
    DO l=1,land_pts
      dvimax_gb(l) = MAXVAL(dvi_cpft(l,:))
    END DO
  END IF

  CALL irrig_dmd(land_pts, sm_levels, frac_irr_surft, a_step, plant_n_gb,     &
                 sthf_gb, smvccl_gb, smvcst_levs, smvcwt_levs, sthzw_gb,      &
                 sthu_irr_gb, sthu_gb,                                        &
                 soil_layer_moisture, irr_crop, dvimax_gb)

  IF( irr_crop == 1 ) THEN
    CALL calc_crop_date(land_index, land_pts, t_i_length, t_j_length, nsurft, &
                        frac_surft, sw_surft, tstar_surft, lw_down, tl_1,     &
                        con_rain, ls_rain, con_snow, ls_snow,                 &
                        plant_n_gb, nday_crop)
  END IF

  IF ( l_irrig_limit ) THEN
    CALL adjust_routestore()
  END IF
ELSE
! if .not. l_irrig_dmd, set sthu_irr_gb to 0.0 in case it is still reported
  sthu_irr_gb(:,:) = 0.0
END IF ! l_irrig_dmd


!-------------------------------------------------------------------------------
! Run crop code if required
!-------------------------------------------------------------------------------
IF ( l_crop ) THEN
  crop_call = MOD ( FLOAT(a_step),                                            &
                    REAL(crop_period) * REAL(secs_in_day) / timestep )

  DO n=1,ncpft
    DO l=1,land_pts
      npp_acc_pft(l,nnpft+n) = npp_acc_pft(l,nnpft+n)                         &
                            + (npp_pft(l,nnpft+n) * timestep)
    END DO
  END DO

  CALL photoperiod(p_field, phot, dphotdt)

  CALL crop(p_field, land_pts, land_index, a_step,                            &
            crop_call, sm_levels, frac_surft, phot, dphotdt,                  &
            sf_diag%t1p5m_surft, t_soil_gb, sthu_gb, satcon_levs, smvccl_gb,  &
            smvcst_levs, npp_acc_pft,                                         &
            canht_pft, lai_pft, dvi_cpft, rootc_cpft, harvc_cpft,             &
            reservec_cpft, croplai_cpft, cropcanht_cpft,                      &
            catch_surft, infil_surft, z0_surft)
END IF  ! l_crop

!------------------------------------------------------------------------------
!   Update metstats for this timestep
!-------------------------------------------------------------------------------
IF ( l_metstats ) THEN
  !Compress variables to land points as metstats has no knowledge of
  !i and j. Also use a TYPE to keep the argument list short
  DO l = 1, land_pts
    j = ( land_index(l)-1 ) / t_i_length + 1
    i = land_index(l) - (j-1) * t_i_length
    metstats_input(l)%temp     = tl_1(i,j)
    metstats_input(l)%spec_hum = qw_1(i,j)
    metstats_input(l)%wind_u   = u_1(i,j)
    metstats_input(l)%wind_v   = v_1(i,j)
    metstats_input(l)%ls_rain  = ls_rain(i,j)
    metstats_input(l)%con_rain = con_rain(i,j)
    metstats_input(l)%ls_snow  = ls_snow(i,j)
    metstats_input(l)%con_snow = con_snow(i,j)
    metstats_input(l)%press    = pstar(i,j)
  END DO

  CALL metstats_timestep(metstats_input,metstats_prog,                        &
       !Things that really ought to come in via USE but won't work with the UM
                         current_time%time, timestep,land_pts)
END IF

!-------------------------------------------------------------------------------
!   Call to fire module
!-------------------------------------------------------------------------------
IF ( l_fire ) THEN
  CALL fire_timestep(metstats_prog, smc_gb, fire_prog, fire_diag,             &
       !Things that really ought to come in via USE but won't work with the UM
                     current_time%time, current_time%month, timestep, land_pts)
END IF
#endif

#if !defined(UM_JULES)
!-------------------------------------------------------------------------------
!   Call to INFERNO (interactive fire module)
!-------------------------------------------------------------------------------
IF ( l_inferno ) THEN
   CALL inferno_io( sf_diag%t1p5m_surft, sf_diag%q1p5m_surft, pstar, sthu_gb, &
                    sm_levels,                                                &
                    frac_surft, dim_cs1, cs_ch4, canht_pft, ls_rain, con_rain,&
                    land_pts, ignition_method)
END IF
#endif

! ----------------------------------------------------------------------
! Section 19 -- VEGETATION DYNAMICS
! ----------------------------------------------------------------------

! initialize carbon conservation diagnostics
! otherwise they can be non-zero on non-triffid timesteps
IF ( l_triffid ) THEN
  DO l = 1, land_pts
    cnsrv_carbon_veg2_gb(l) = 0.0
    cnsrv_carbon_triffid_gb(l) = 0.0
    cnsrv_veg_triffid_gb(l) = 0.0
    cnsrv_soil_triffid_gb(l) = 0.0
    cnsrv_prod_triffid_gb(l) = 0.0
  END DO
ENDIF

#if defined(UM_JULES)
!change 2d to 1d soil clay content for soil respiration
IF ( l_triffid ) THEN
  DO l = 1, land_pts
    j = (land_index(l)-1)/row_length + 1
    i = land_index(l) - (j-1)*row_length
    clay_gb(l) = soil_clay_ij(i,j)
  ENDDO
ENDIF
#endif

!-------------------------------------------------------------------------
!   If leaf phenology is activated, check whether the atmosphere model
!   has run an integer number of phenology calling periods.
!-------------------------------------------------------------------------
#if defined(UM_JULES)

IF (l_phenol .OR. l_triffid) THEN
  CALL veg_control(                                                           &
    land_pts, land_index, nsurft, can_model,                                  &
    a_step, asteps_since_triffid,                                             &
    land_pts_trif, npft_trif,                                                 &
    phenol_period, triffid_period, row_length, rows,                          &
    l_phenol, l_triffid, l_trif_eq,                                           &
    timestep, frac_disturb, frac_past_gb, satcon_levs,                        &
    g_leaf_acc_pft, g_leaf_phen_acc_pft, npp_acc_pft,                         &
    resp_s_acc_gb, resp_w_acc_pft,                                            &
    cs_pool_gb, frac_surft, lai_pft, clay_gb, z0m_soil_gb, canht_pft,    &
    catch_snow_surft, catch_surft, infil_surft, z0_surft, z0h_bare_surft,     &
    c_veg_pft, cv_gb, lit_c_pft, lit_c_mn_gb, g_leaf_day_pft, g_leaf_phen_pft,&
    lai_phen_pft, g_leaf_dr_out_pft, npp_dr_out_pft, resp_w_dr_out_pft,       &
    resp_s_dr_out_gb                                                          &
     )
END IF
#else
phenol_call=1
triffid_call=1
IF ( l_phenol ) phenol_call = MOD ( FLOAT(a_step),                            &
            REAL(phenol_period) * REAL(secs_in_day) / timestep )

IF ( l_triffid ) THEN
  nstep_trif = INT( REAL(secs_in_day) * REAL(triffid_period) / timestep )
  IF ( asteps_since_triffid == nstep_trif ) triffid_call = 0
ENDIF

IF ( triffid_call == 0 ) THEN
!-------------------------------------------------------------------------------
!     Run includes dynamic vegetation
!-------------------------------------------------------------------------------
  CALL veg2( land_pts, land_index, nsurft, can_model                          &
            ,a_step, asteps_since_triffid                                     &
            ,phenol_period, triffid_period, l_phenol, l_triffid, l_trif_eq    &
            ,timestep, frac_agr_gb, frac_past_gb, satcon_levs                 &
            ,g_leaf_acc_pft, g_leaf_phen_acc_pft, npp_acc_pft                 &
            ,resp_s_acc_gb, resp_w_acc_pft                                    &
            ,cs_pool_gb, frac_surft, lai_pft, clay_gb, z0m_soil_gb            &
            ,canht_pft, catch_snow_surft, catch_surft, infil_surft, z0_surft  &
            ,z0h_bare_surft, c_veg_pft, cv_gb, lit_c_pft, lit_c_mn_gb         &
            ,g_leaf_day_pft, g_leaf_phen_pft, lai_phen_pft, g_leaf_dr_out_pft &
            ,npp_dr_out_pft, resp_w_dr_out_pft, resp_s_dr_out_gb )

ELSE

  IF ( phenol_call == 0 )                                                     &
!-------------------------------------------------------------------------------
!       Run includes phenology,  but not dynamic vegetation
!       therefore call veg1 rather than veg2
!-------------------------------------------------------------------------------
    CALL veg1( land_pts, nsurft, can_model, a_step, phenol_period, l_phenol   &
              ,timestep, satcon_levs, z0m_soil_gb, g_leaf_acc_pft, frac_surft &
              ,lai_pft,  canht_pft                                            &
              ,catch_snow_surft, catch_surft, infil_surft, z0_surft           &
              ,z0h_bare_surft                                                 &
              ,g_leaf_day_pft, g_leaf_phen_pft, g_leaf_phen_acc_pft           &
              ,lai_phen_pft )

ENDIF  !  triffid_call

IF (l_fao_ref_evapotranspiration) THEN
   trad = ( tiles_to_gbm(tstar_surft**4) )**0.25
   CALL fao_ref_evapotranspiration(soil_pts, soil_index,                      &
     land_pts, land_index, sf_diag%t1p5m,                                     &
     sw_down_ij, lw_down_ij, surf_ht_flux_ij, sf_diag%u10m,                   &
     sf_diag%v10m, sf_diag%q1p5m, pstar_ij, trad, fao_et0)
END IF

#endif

!------------------------------------------------------------------------------
!Call STASH diagnostic routines - UM-only
!------------------------------------------------------------------------------
#if defined(UM_JULES)
IF (model_type /= mt_single_column) THEN
  IF (l_hydrology .AND. sf(0,8) ) THEN
!DEPENDS ON: diagnostics_hyd
    CALL diagnostics_hyd(                                                     &
      row_length, rows, model_levels,                                         &
      n_rows, global_row_length, global_rows,                                 &
      halo_i, halo_j, offx, offy, mype,                                       &
      n_proc, n_procx, n_procy,                                               &
      g_rows, g_row_length,                                                   &
      at_extremity,                                                           &
      land_pts, sm_levels,                                                    &
      !Put inland basin outflow in call to diagriv
      land_index,inlandout_atm_gb,                                            &
      smc_gb, surf_roff, sub_surf_roff,                                       &
      lying_snow, snow_melt,                                                  &
      canopy_gb,t_soil_gb,                                                    &
      tsurf_elev_surft,snow_soil_htf,snow_smb_surft,                          &
      soil_layer_moisture,                                                    &
      nsurft, snomlt_surf_htf, sthu_gb, sthf_gb,                              &
      tot_tfall, snow_surft, melt_surft,                                      &
      rgrain_surft, land_sea_mask,                                            &
      dun_roff_gb, drain_gb, qbase_gb, qbase_zw_gb,                           &
      fch4_wetl_gb,fch4_wetl_cs_gb,fch4_wetl_npp_gb,                          &
      fch4_wetl_resps_gb,                                                     &
      fexp_gb,gamtot_gb,ti_mean_gb,ti_sig_gb,                                 &
      fsat_gb,fwetl_gb,zw_gb,sthzw_gb,                                        &
      timestep,                                                               &
      STASHwork8                                                              &
      )
  END IF

  IF ( l_rivers .AND. trip_call .AND. sf(0,26) ) THEN

    CALL diagnostics_riv(                                                     &
      row_length, rows,                                                       &
      river_row_length, river_rows,                                           &
      at_extremity,                                                           &
      at_extremity,                                                           &
      riverout,                                                               &
      box_outflow, box_inflow,                                                &
      !Put inland basin outflow in call to diagriv
      twatstor,inlandout_riv,                                                 &
      STASHwork26                                                             &
      )
  END IF

  IF (sf(0,19)) THEN
! DEPENDS ON: diagnostics_veg
    CALL diagnostics_veg(                                                     &
      row_length, rows, n_rows,                                               &
      global_row_length, global_rows,                                         &
      dim_cs1, dim_cs2,                                                       &
      halo_i, halo_j, offx, offy, mype,                                       &
      n_proc, n_procx, n_procy,                                               &
      g_rows, g_row_length,                                                   &
      at_extremity,                                                           &
      land_pts,                                                               &
      land_index,                                                             &
      ntype,npft,                                                             &
      c_veg_pft,cv_gb,g_leaf_phen_pft,                                        &
      lit_c_pft,lit_c_mn_gb,g_leaf_day_pft,                                   &
      lai_phen_pft,g_leaf_dr_out_pft,npp_dr_out_pft,                          &
      resp_w_dr_out_pft,resp_s_dr_out_gb,frac_disturb,disturb_veg_prev,       &
      wood_prod_fast_d1, wood_prod_med_d1, wood_prod_slow_d1,                 &
      frac_surft,lai_pft,canht_pft,cs_pool_gb,                                &
      STASHwork19                                                             &
      )
  END IF
END IF
#endif
        
    CASE DEFAULT
      CALL ereport('surf_couple_explicit', 101, 'Unrecognised surface scheme')

  END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE surf_couple_extra
END MODULE surf_couple_extra_mod
