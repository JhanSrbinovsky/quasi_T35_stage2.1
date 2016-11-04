! *****************************COPYRIGHT****************************************
! (c) Crown copyright, Met Office 2014. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms and
! conditions set out therein.
!
! [Met Office Ref SC0237]
! *****************************COPYRIGHT****************************************

MODULE surf_couple_explicit_mod

IMPLICIT NONE

PRIVATE
PUBLIC :: surf_couple_explicit

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SURF_COUPLE_EXPLICIT_MOD'

CONTAINS

!===============================================================================
! Public subroutine
!===============================================================================
SUBROUTINE surf_couple_explicit(                                              &
    !Arguments used by JULES-standalone
    !Misc INTENT(IN) DONE
    bq_1, bt_1, zh, photosynth_act_rad,                                       &
    curr_year, curr_day_number, curr_hour, curr_minute, curr_second,          &
    !Forcing INTENT(IN) DONE
    qw, tl, pstar, lw_down,                                                   &
    !Fluxes INTENT(IN) DONE
    sw_surft, tstar,                                                          &
    !Diagnostics, INTENT(INOUT)
    sf_diag,                                                                  &
    !Fluxes INTENT(OUT) DONE
    fqw_1,ftl_1,ftl_surft, fqw_surft, fqw_ice, ftl_ice, fsmc, emis_surft,     &
    !Misc INTENT(OUT)
    radnet_sice, rhokm_1, rhokm_land, rhokm_ssi,                              &
    !Out of explicit and into implicit only INTENT(OUT)
    cdr10m, cdr10m_n, cd10m_n,                                                &
    alpha1, alpha1_sice, ashtf_prime, ashtf_prime_surft, epot_surft,          &
    fraca, resfs, resft, rhokh, rhokh_surft, rhokh_sice, rhokh_sea,           &
    dtstar_surft, dtstar,                                                     &
    z0hssi, z0h_surft, z0mssi, z0m_surft, chr1p5m, chr1p5m_sice, canhc_surft, &
    wt_ext_surft, flake,                                                      &
    !Out of explicit and into extra only INTENT(OUT)
    hcons,                                                                    &
    !Out of explicit and into implicit and extra INTENT(OUT)
    tile_frac                                                                 &
#if defined(UM_JULES)
    ,                                                                         &
    !Additional arguments for the UM-----------------------------------------
    !JULES prognostics module
    canopy_surft, snow_surft, k_sice_sicat, cs_pool_gb, canht_pft, lai_pft,   &
    t_soil_gb, tsurf_elev_surft, ti_sicat, ti_cat_sicat,                      &
    tstar_surft, z0msea_ij, smc_gb, gc_surft, gs_gb,                          &
    !JULES ancil_info module
    land_pts, z1_uv_ij, z1_tq_ij, land_index, nsurft, ice_fract_ncat_sicat,   &
    frac_surft, surft_index, surft_pts,                                       &
    !JULES coastal module
    fland, flandg, tstar_sea_ij, tstar_sice_sicat, vshr_land_ij, vshr_ssi_ij, &
    !JULES aero module
    co2_3d_ij, rho_aresist_ij, aresist_ij, resist_b_ij, rho_aresist_surft,    &
    aresist_surft, resist_b_surft, r_b_dust_ij, cd_std_dust_ij, u_s_std_surft,&
    !JULES trifctl module
    asteps_since_triffid, g_leaf_acc_pft, npp_acc_pft, resp_w_acc_pft,        &
    resp_s_acc_gb, gpp_gb, npp_gb, resp_p_gb, g_leaf_pft, gpp_pft, npp_pft,   &
    resp_p_pft, resp_s_gb, resp_w_pft,                                        &
    !JULES p_s_parms module
    catch_surft, catch_snow_surft, hcon_gb, smvccl_gb, smvcst_gb, smvcwt_gb,  &
    sthf_gb, sthu_gb, z0_surft,                                               &
    z0h_bare_surft, z0m_soil_gb, albsoil_gb, cosz_gb, soil_clay_ij,           &
    !JULES orog module
    ho2r2_orog_gb, sil_orog_land_gb, h_blend_orog_ij, z0m_eff_ij,             &
    !JULES u_v_grid module
    u_1_p_ij, v_1_p_ij, u_0_p_ij, v_0_p_ij,                                   &
    !JULES switches module
    l_spec_z0,                                                                &
    !JULES c_elevate module
    z_land_ij,                                                                &
    !Not in a JULES module
    numcycles, cycleno, z1_uv_top, z1_tq_top, ddmfx,                          &
    l_aero_classic, z0m_scm, z0h_scm, recip_l_mo_sea, rib, rib_surft,         &
    flandfac, fseafac, wt_ext,fb_surf, u_s, t1_sd, q1_sd, rhostar,            &
    z0m_gb, z0h_eff, vshr, resp_s_tot, emis_soil                              &
#endif
    )

!Module imports

!Common modules
USE ereport_mod, ONLY : ereport

USE jules_sea_seaice_mod,     ONLY:                                           &
  nice, nice_use
USE jules_soil_mod,           ONLY:                                           &
  sm_levels
USE jules_vegetation_mod,     ONLY:                                           &
  l_triffid, l_phenol, l_q10
USE jules_surface_types_mod,  ONLY:                                           &
  npft, ntype
USE sf_diags_mod, ONLY: strnewsfdiag
USE ozone_vars,               ONLY:                                           &
  o3_gb
USE atm_fields_bounds_mod,    ONLY:                                           &
  pdims_s, pdims, tdims
USE dust_param,               ONLY:                                           &
  ndiv
USE trif_vars_mod,            ONLY:                                           &
  resp_l_pft, resp_r_pft, n_leaf_pft, n_stem_pft, n_root_pft, lai_bal_pft

!Modules that change name between JULES and UM
#if defined(UM_JULES)
  USE atm_step_local,         ONLY:                                           &
    land_pts_trif, npft_trif,dim_cs1, dim_cs2, co2_dim_len,co2_dim_row
  USE carbon_options_mod,     ONLY:                                           &
    l_co2_interactive
  USE dust_parameters_mod,    ONLY:                                           &
    l_dust, l_dust_diag
!  USE s_main_force,           ONLY: & !A SCM module currently doesn't work
!    l_spec_z0
  USE rad_input_mod,          ONLY:                                           &
    co2_mmr
  USE gen_phys_inputs_mod,    ONLY:                                           &
    lq_mix_bl => l_mr_physics2
#else
  USE ancil_info,             ONLY:                                           &
    land_pts_trif, npft_trif,dim_cs1, dim_cs2, co2_dim_len,co2_dim_row,       &
    land_pts, nsurft
  USE prognostics,            ONLY:                                           &
    canopy_surft, snow_surft, k_sice_sicat, cs_pool_gb, canht_pft, lai_pft,   &
    t_soil_gb, tsurf_elev_surft, ti_sicat, tstar_surft,                       &
    z0msea_ij, gs_gb, smc_gb, gc_surft
  USE ancil_info,             ONLY:                                           &
    land_pts, npft_trif,dim_cs1, dim_cs2, co2_dim_len,co2_dim_row,            &
    z1_uv_ij, z1_tq_ij, land_index, nsurft, ice_fract_ncat_sicat,             &
    ti_cat_sicat, frac_surft, surft_index, surft_pts
  USE coastal,                ONLY:                                           &
    fland, flandg, tstar_sea_ij, tstar_sice_sicat, vshr_land_ij, vshr_ssi_ij
  USE aero,                   ONLY:                                           &
    co2_mmr, co2_3d_ij, rho_aresist_ij, aresist_ij, resist_b_ij,              &
    rho_aresist_surft, aresist_surft, resist_b_surft, r_b_dust_ij,            &
    cd_std_dust_ij, u_s_std_surft
  USE trifctl,                ONLY:                                           &
    asteps_since_triffid, g_leaf_acc_pft,npp_acc_pft,resp_w_acc_pft,          &
    resp_s_acc_gb, gpp_gb, npp_gb, resp_p_gb, g_leaf_pft, gpp_pft, npp_pft,   &
    resp_p_pft, resp_s_gb, resp_w_pft
  USE p_s_parms,              ONLY:                                           &
    catch_surft,catch_snow_surft,hcon_gb,smvccl_gb,smvcst_gb,smvcwt_gb,       &
    sthf_gb,sthu_gb,z0_surft,z0h_bare_surft, z0m_soil_gb, albsoil_gb, cosz_gb,&
    clay_gb
  USE orog,                   ONLY:                                           &
    ho2r2_orog_gb, sil_orog_land_gb, h_blend_orog_ij, z0m_eff_ij
  USE u_v_grid,               ONLY:                                           &
    u_0_p_ij, u_1_p_ij, v_0_p_ij, v_1_p_ij
  USE switches,               ONLY:                                           &
    l_co2_interactive, lq_mix_bl, l_dust, l_spec_z0
  USE c_elevate,              ONLY:                                           &
    z_land_ij
#endif

!CABLE_LSM:
USE lsm_switches_mod,        ONLY: lsm_id
USE um_parcore,              ONLY : mype
USE timestep_mod,            ONLY : timestep_number 
USE p_s_parms,                ONLY :                                          &
  bexp_gb, sathh_gb, satcon_gb
USE atm_fields_real_mod, ONLY : soil_alb!, lw_down
USE trignometric_mod, ONLY :                                                  &
  sin_theta_latitude, true_latitude, true_longitude
USE cable_gather_um_data_decs, ONLY :                                         & 
  ls_rain_cable, ls_snow_cable, surf_down_sw_cable
!CABLE_LSM: End

!Dr Hook
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Coupling routine between the UM or JULES system code and land surface
!   explicit science routines. Calls the appropriate LSM-specific code.
!
! Current Code Owner: Richard Gilham
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

! Subroutine arguments

!UM-only arguments
#if defined(UM_JULES)
!JULES prognostics module
REAL, INTENT(IN) ::                                                           &
  canopy_surft(land_pts,nsurft),                                              &
  snow_surft(land_pts,nsurft),                                                &
  k_sice_sicat(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use), &
  cs_pool_gb(land_pts,dim_cs1),                                               &
  canht_pft(land_pts,npft),                                                   &
  lai_pft(land_pts,npft),                                                     &
  t_soil_gb(land_pts,sm_levels),                                              &
  tsurf_elev_surft(land_pts,nsurft),                                          &
  ti_sicat(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),              &
  ti_cat_sicat(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice),     &
  tstar_surft(land_pts,nsurft)
REAL, INTENT(INOUT) ::                                                        &
  z0msea_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),             &
  gs_gb(land_pts)
REAL, INTENT(OUT) ::                                                          &
  smc_gb(land_pts),                                                           &
  gc_surft(land_pts,nsurft)
!JULES ancil_info module
!JULES ancil_info module
INTEGER, INTENT(IN) ::                                                        &
  land_pts,                                                                   &
  land_index(land_pts),                                                       &
  nsurft
REAL, INTENT(IN) ::                                                           &
  z1_uv_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),              &
  z1_tq_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
INTEGER, INTENT(OUT)  ::                                                      &
  surft_index(land_pts,nsurft),                                               &
  surft_pts(nsurft)
REAL, INTENT(IN) ::                                                           &
  ice_fract_ncat_sicat(tdims%i_start:tdims%i_end,                             &
                       tdims%j_start:tdims%j_end,nice_use),                   &
  frac_surft(land_pts,ntype)
!JULES coastal module
REAL, INTENT(IN) ::                                                           &
  fland(land_pts),                                                            &
  flandg(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),        &
  tstar_sea_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),          &
  tstar_sice_sicat(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use)
REAL, INTENT(OUT) ::                                                          &
  vshr_land_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),          &
  vshr_ssi_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
!JULES aero module
REAL, INTENT(IN) ::                                                           &
  co2_3d_ij(co2_dim_len,co2_dim_row)
REAL, INTENT(OUT) ::                                                          &
  rho_aresist_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),        &
  aresist_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),            &
  resist_b_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),           &
  rho_aresist_surft(land_pts,nsurft),                                         &
  aresist_surft(land_pts,nsurft),                                             &
  resist_b_surft(land_pts,nsurft),                                            &
  r_b_dust_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,ndiv),      &
  cd_std_dust_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),        &
  u_s_std_surft(land_pts,nsurft)
!JULES trifctl module
INTEGER, INTENT(IN)::                                                         &
  asteps_since_triffid
REAL, INTENT(INOUT) ::                                                        &
  g_leaf_acc_pft(land_pts,npft),                                              &
  npp_acc_pft(land_pts_trif,npft_trif),                                       &
  resp_w_acc_pft(land_pts_trif,npft_trif),                                    &
  resp_s_acc_gb(land_pts_trif,dim_cs1)
REAL, INTENT(OUT) ::                                                          &
  gpp_gb(land_pts),                                                           &
  npp_gb(land_pts),                                                           &
  resp_p_gb(land_pts),                                                        &
  g_leaf_pft(land_pts,npft),                                                  &
  gpp_pft(land_pts,npft),                                                     &
  npp_pft(land_pts,npft),                                                     &
  resp_p_pft(land_pts,npft),                                                  &
  resp_s_gb(land_pts,dim_cs1),                                                &
  resp_w_pft(land_pts,npft)
!JULES p_s_parms module
REAL, INTENT(IN) ::                                                           &
  catch_surft(land_pts,nsurft),                                               &
  catch_snow_surft(land_pts,nsurft),                                          &
  hcon_gb(land_pts),                                                          &
  smvccl_gb(land_pts,sm_levels),                                              &
  smvcst_gb(land_pts,sm_levels),                                              &
  smvcwt_gb(land_pts,sm_levels),                                              &
  sthf_gb(land_pts,sm_levels),                                                &
  sthu_gb(land_pts,sm_levels),                                                &
  z0_surft(land_pts,nsurft),                                                  &
  z0h_bare_surft(land_pts,nsurft),                                            &
  z0m_soil_gb(land_pts),                                                      &
  albsoil_gb(land_pts),                                                       &
  cosz_gb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
  soil_clay_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
!JULES orog module
REAL, INTENT(IN) ::                                                           &
  ho2r2_orog_gb(land_pts),                                                    &
  sil_orog_land_gb(land_pts),                                                 &
  h_blend_orog_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
REAL, INTENT(OUT) ::                                                          &
  z0m_eff_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
!JULES u_v_grid module
REAL, INTENT(IN) ::                                                           &
  u_1_p_ij(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),      &
  v_1_p_ij(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),      &
  u_0_p_ij(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),      &
  v_0_p_ij(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end)
!JULES switches module
LOGICAL, INTENT(IN) ::                                                        &
  l_spec_z0
!JULES c_elevate module
REAL, INTENT(IN) ::                                                           &
  z_land_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
!Not in a JULES module
INTEGER, INTENT(IN) :: numcycles, cycleno
REAL, INTENT(IN) ::                                                           &
  !These variables are INTENT(IN) to sf_expl, but not used with the
  !current configuration of standalone JULES (initialised to 0 below)
  z1_uv_top(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),             &
    ! Height of top of lowest uv-layer
  z1_tq_top(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),             &
    ! Height of top of lowest Tq-layer
  ddmfx(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
    ! Convective downdraught mass-flux at cloud base
LOGICAL, INTENT(IN) :: l_aero_classic
REAL, INTENT(IN) ::                                                           &
  !These variables are required for prescribed roughness lengths in
  !SCM mode in UM - not used standalone
  z0m_scm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
    !Fixed Sea-surface roughness length for momentum (m).(SCM)
  z0h_scm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
    !Fixed Sea-surface roughness length for heat (m). (SCM)
REAL, INTENT(OUT) ::                                                          &
  recip_l_mo_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),        &
    !Reciprocal of the surface Obukhov  length at sea points. (m-1).
  rib(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                   &
    !Mean bulk Richardson number for lowest layer.
  rib_surft(LAND_PTS,nsurft),                                                 &
    !RIB for land tiles.
  flandfac(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),      &
  fseafac(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),       &
  wt_ext(land_pts,sm_levels),                                                 &
    !cumulative fraction of transp'n
  fb_surf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
    !Surface flux buoyancy over density (m^2/s^3)
  u_s(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                   &
    !Surface friction velocity (m/s)
  t1_sd(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &
    !Standard deviation of turbulent fluctuations of layer 1 temp; used in
    !initiating convection.
  q1_sd(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &
    !Standard deviation of turbulent flux of layer 1 humidity; used in
    !initiating convection.
  rhostar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
    !Surface air density
  z0m_gb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                &
    !Gridbox mean roughness length for momentum (m).
  z0h_eff(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
    !Effective grid-box roughness length for heat, moisture (m)
  vshr(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                  &
    !Magnitude of surface-to-lowest atm level wind shear (m per s).
  resp_s_tot(dim_cs2),                                                        &
    !Total soil respiration (kg C/m2/s).
  emis_soil(land_pts)
    !Emissivity of underlying soil
REAL ::                                                                       &
  clay_gb(land_pts)
    !Clay content of underlying soil
#else
  INTEGER, PARAMETER :: numcycles = 1 !Number of cycles (iterations) for
                                      !iterative SISL.
  INTEGER, PARAMETER :: cycleno = 1   !Iteration no
REAL::                                                                        &
  !These variables are INTENT(IN) to sf_expl, but not used with the
  !current configuration of standalone JULES (initialised to 0 below)
  z1_uv_top(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),             &
    ! Height of top of lowest uv-layer
  z1_tq_top(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),             &
    ! Height of top of lowest Tq-layer
  ddmfx(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
    ! Convective downdraught mass-flux at cloud base
LOGICAL :: l_dust_diag
  !In standalone, this switch essentially does the same job as l_dust
LOGICAL, PARAMETER :: l_aero_classic = .FALSE. ! switch for CLASSIC aerosol
                                               !scheme - NEVER USED!!!
REAL::                                                                        &
  !These variables are required for prescribed roughness lengths in
  !SCM mode in UM - not used standalone
  z0m_scm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
    !Fixed Sea-surface roughness length for momentum (m).(SCM)
  z0h_scm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
    !Fixed Sea-surface roughness length for heat (m). (SCM)
  recip_l_mo_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),        &
    !Reciprocal of the surface Obukhov  length at sea points. (m-1).
  rib(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                   &
    !Mean bulk Richardson number for lowest layer.
  rib_surft(LAND_PTS,nsurft),                                                 &
    !RIB for land tiles.
  flandfac(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),      &
  fseafac(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),       &
  wt_ext(land_pts,sm_levels),                                                 &
    !cumulative fraction of transp'n
  fb_surf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
    !Surface flux buoyancy over density (m^2/s^3)
  u_s(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                   &
    !Surface friction velocity (m/s)
  t1_sd(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &
    !Standard deviation of turbulent fluctuations of layer 1 temp; used in
    !initiating convection.
  q1_sd(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &
    !Standard deviation of turbulent flux of layer 1 humidity; used in
    !initiating convection.
  rhostar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
    !Surface air density
  z0m_gb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                &
    !Gridbox mean roughness length for momentum (m).
  z0h_eff(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
    !Effective grid-box roughness length for heat, moisture (m)
  vshr(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                  &
    !Magnitude of surface-to-lowest atm level wind shear (m per s).
  resp_s_tot(dim_cs2),                                                        &
    !Total soil respiration (kg C/m2/s).
  emis_soil(land_pts)
#endif

!Misc INTENT(IN)
REAL, INTENT(IN) ::                                                           &
  bq_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                  &
    !A buoyancy parameter (beta q tilde).
  bt_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                  &
    !A buoyancy parameter (beta T tilde).
  zh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                    &
    !Height above surface of top of boundary layer (metres).
  photosynth_act_rad(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
    !Net downward shortwave radiation in band 1 (w/m2).
INTEGER, INTENT(IN) ::                                                        &
  curr_year, curr_day_number, curr_hour, curr_minute, curr_second

!Forcing INTENT(IN)
REAL, INTENT(IN) ::                                                           &
  qw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                    &
    !Total water content
  tl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                    &
    !Ice/liquid water temperature
  pstar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                 &
    !Surface pressure (Pascals).
  lw_down(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
    !Surface downward LW radiation (W/m2).

!Fluxes INTENT(IN)
REAL, INTENT(IN) ::                                                           &
  sw_surft(land_pts,nsurft),                                                  &
    !Surface net SW radiation on land tiles (W/m2).
  tstar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
    !GBM surface temperature (K).

!diagnostic array
TYPE (strnewsfdiag), INTENT(INOUT) :: sf_diag

!Fluxes INTENT(OUT)
REAL, INTENT(OUT) ::                                                          &
  fqw_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &
    !Moisture flux between layers (kg per square metre per sec).
  ftl_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &
    !FTL(,K) contains net turbulent sensible heat flux into layer K from below;
    !so FTL(,1) is the surface sensible heat, H.(W/m2)
  ftl_surft(land_pts,nsurft),                                                 &
    !Surface FTL for land tiles
  fqw_surft(land_pts,nsurft),                                                 &
    !Surface FQW for land tiles
  fqw_ice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),      &
    !Surface FQW for sea-ice
  ftl_ice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),      &
    !Surface FTL for sea-ice
  fsmc(land_pts,npft),                                                        &
    !Moisture availability factor.
  emis_surft(land_pts,nsurft)                                                
    !Emissivity for land tiles

!Misc INTENT(OUT)
REAL, INTENT(OUT) ::                                                          &
  radnet_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),  &
    !Surface net radiation on sea-ice (W/m2)
  rhokm_1(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),       &
    !Exchange coefficients for momentum on P-grid
  rhokm_land(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),    &
  rhokm_ssi(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end)

!Out of explicit and into implicit only INTENT(OUT)
REAL, INTENT(OUT) ::                                                          &
  cdr10m(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),        &
!   Interpolation coefficient for diagnosis of 10 m winds
  cdr10m_n(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),      &
!   Interpolation coefficient for diagnosis of neutral 10 m winds
  cd10m_n(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),       &
!   Neutral drag coefficient for calculation of pseudostress
  alpha1(land_pts,nsurft),                                                    &
    !Mean gradient of saturated specific humidity with respect to temperature
    !between the bottom model layer and tile surfaces
  alpha1_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),  &
    !ALPHA1 for sea-ice.
  ashtf_prime(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),  &
    !Coefficient to calculate surface heat flux into sea-ice.
  ashtf_prime_surft(land_pts,nsurft),                                         &
    !Coefficient to calculate surface heat flux into land tiles.
  epot_surft(land_pts,nsurft),                                                &
    !Local EPOT for land tiles.
  fraca(land_pts,nsurft),                                                     &
    !Fraction of surface moisture flux with only aerodynamic resistance for
    !snow-free land tiles.
  resfs(land_pts,nsurft),                                                     &
    !Combined soil, stomatal and aerodynamic resistance factor for fraction
    !(1-FRACA) of snow-free land tiles.
  resft(land_pts,nsurft),                                                     &
    !Total resistance factor. FRACA+(1-FRACA)*RESFS for snow-free land, 1 for
    !snow.
  rhokh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &
    !Grid-box surface exchange coefficients
  rhokh_surft(land_pts,nsurft),                                               &
    !Surface exchange coefficients for land tiles
  rhokh_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),   &
    !Surface exchange coefficients for sea-ice
  rhokh_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),             &
    !Surface exchange coefficients for sea
  dtstar_surft(land_pts,nsurft),                                              &
    !Change in TSTAR over timestep for land tiles
  dtstar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),       &
    !Change is TSTAR over timestep for sea-ice
  z0hssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                &
    !Roughness length for heat and moisture over sea (m).
  z0h_surft(land_pts,nsurft),                                                 &
    !Tile roughness lengths for heat and moisture (m).
  z0mssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                &
    !Roughness length for momentum over sea (m).
  z0m_surft(land_pts,nsurft),                                                 &
    !Tile roughness lengths for momentum.
  chr1p5m(land_pts,nsurft),                                                   &
    !Ratio of coefffs for calculation of 1.5m temp for land tiles.
  chr1p5m_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),          &
    !CHR1P5M for sea and sea-ice (leads ignored).
  canhc_surft(land_pts,nsurft),                                               &
    !Areal heat capacity of canopy for land tiles (J/K/m2).
  wt_ext_surft(land_pts,sm_levels,nsurft),                                    &
    !Fraction of evapotranspiration which is extracted from each soil layer
    !by each tile.
  flake(land_pts,nsurft)
    !Lake fraction.

!Out of explicit and into extra only INTENT(OUT)
REAL, INTENT(OUT) ::                                                          &
  hcons(land_pts)
    !Soil thermal conductivity including water and ice

!Out of explicit and into implicit and extra INTENT(OUT)
REAL, INTENT(OUT) ::                                                          &
  tile_frac(land_pts,nsurft)
    !Tile fractions including snow cover in the ice tile.

INTEGER ::                                                                    &
   i,j,l         !Various counters

!-----------------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------------
! Temp until the model switching is implemented
! Having a parameter until then should hopefully help the compiler eliminate
! dead code

  !Dr Hook variables
  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='SURF_COUPLE_EXPLICIT'

!-----------------------------------------------------------------------------
!End of header
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!change 2d UM clay to 1d jules clay content for soil respiration
#if defined(UM_JULES)
IF ( l_triffid ) THEN
  DO l = 1, land_pts
    j = (land_index(l)-1)/pdims%i_end + 1
    i = land_index(l) - (j-1)*pdims%i_end
    clay_gb(l) = soil_clay_ij(i,j)
  ENDDO
ENDIF
#endif

#if !defined(UM_JULES)
  z0m_scm(:,:)       = 0.0
  z0h_scm(:,:)       = 0.0
  z1_uv_top          = 0.0
  z1_tq_top          = 0.0
  ddmfx(:,:)         = 0.0
  l_dust_diag = l_dust
!-----------------------------------------------------------------------
! Set message passing variables for call to sf_expl
!-----------------------------------------------------------------------
! Copy wind fields into halo versions (note that halos are 0)
!u_0_px(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end) = u_0_p_ij(:,:)
!v_0_px(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end) = v_0_p_ij(:,:)
!u_1_p_ijx(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end) = u_1_p_ij(:,:)
!v_1_px(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end) = v_1_p_ij(:,:)
#endif
  SELECT CASE( lsm_id )
    CASE( 'jules' )
! DEPENDS ON: sf_expl_l
      CALL sf_expl_l (                                                        &
        !IN date-related values
        curr_year, curr_day_number, curr_hour, curr_minute, curr_second,      &
        !IN values defining field dimensions and subset to be processed :
        land_pts, nice, nice_use,                                             &
        !IN  parameters for iterative SISL scheme
        numcycles, cycleno,                                                   &
        !IN parameters required from boundary-layer scheme :
        bq_1,bt_1,z1_uv_ij,z1_uv_top,z1_tq_ij,z1_tq_top,qw,tl,                &
        !IN soil/vegetation/land surface data :
        land_index,nsurft,sm_levels,canopy_surft,catch_surft,catch_snow_surft,&
        hcon_gb,ho2r2_orog_gb, fland,flandg,                                  &
        snow_surft,sil_orog_land_gb,smvccl_gb,smvcst_gb,smvcwt_gb,sthf_gb,    &
        sthu_gb,z0_surft, z0h_bare_surft, z0m_soil_gb,                        &
        !IN sea/sea-ice data :
        ice_fract_ncat_sicat,k_sice_sicat,                                    &
        !IN everything not covered so far :
        pstar,lw_down,sw_surft,zh,ddmfx,                                      &
        co2_mmr,co2_3d_ij,l_co2_interactive,l_phenol,l_triffid,               &
        asteps_since_triffid,cs_pool_gb,frac_surft,canht_pft,                 &
        photosynth_act_rad, lai_pft,                                          &
        lq_mix_bl,t_soil_gb,tsurf_elev_surft,                                 &
        ti_sicat,ti_cat_sicat,tstar,tstar_sea_ij,                             &
        tstar_sice_sicat,                                                     &
        tstar_surft,z_land_ij,albsoil_gb,cosz_gb,                             &
        l_aero_classic,l_dust,l_dust_diag,clay_gb,o3_gb,                      &
        !IN idealised and SCM things
        l_spec_z0, z0m_scm, z0h_scm,                                          &
        !IN variables for message passing
        u_1_p_ij, v_1_p_ij, u_0_p_ij, v_0_p_ij,                               &
        !INOUT diagnostics
        sf_diag,                                                              &
        !INOUT data :
        l_q10,z0msea_ij,gs_gb,g_leaf_acc_pft,npp_acc_pft,resp_w_acc_pft,      &
        resp_s_acc_gb,                                                        &
        !OUT Diagnostic not requiring STASH flags :
        recip_l_mo_sea,fqw_1,ftl_1,ftl_surft,                                 &
        radnet_sice,rhokm_1,rib,rib_surft,                                    &
        !OUT variables for message passing
        flandfac, fseafac, rhokm_land, rhokm_ssi,                             &
        cdr10m, cdr10m_n, cd10m_n,                                            &
        !OUT diagnostics required for soil moisture nudging scheme :
        wt_ext,                                                               &
        !OUT data required for tracer mixing :
        rho_aresist_ij,aresist_ij,resist_b_ij,                                &
        rho_aresist_surft,aresist_surft,resist_b_surft,                       &
        !OUT data required for mineral dust scheme
        r_b_dust_ij,cd_std_dust_ij,u_s_std_surft,                             &
        !OUT data required elsewhere in UM system :
        fb_surf,u_s,t1_sd,q1_sd,                                              &
        !OUT data required elsewhere in boundary layer or surface code
        alpha1,alpha1_sice,ashtf_prime,ashtf_prime_surft,fqw_surft,           &
        epot_surft,fqw_ice,ftl_ice,fraca,rhostar,resfs,resft,                 &
        rhokh,rhokh_surft,rhokh_sice,rhokh_sea,dtstar_surft,dtstar,           &
        h_blend_orog_ij,z0hssi,z0h_surft,z0h_eff,z0m_gb,z0mssi,z0m_surft,     &
        z0m_eff_ij,chr1p5m,chr1p5m_sice,smc_gb,hcons,vshr,vshr_land_ij,       &
        vshr_ssi_ij, gpp_gb,npp_gb,resp_p_gb,g_leaf_pft,gpp_pft,npp_pft,      &
        resp_p_pft,resp_s_gb,resp_s_tot,resp_l_pft,resp_r_pft,resp_w_pft,     &
        n_leaf_pft,n_root_pft,n_stem_pft,lai_bal_pft,                         &
        gc_surft,canhc_surft,wt_ext_surft,flake,                              &
        surft_index,surft_pts,tile_frac,fsmc,emis_surft,emis_soil             &
        )

    CASE( 'cable' )

! DEPENDS ON: sf_expl_l_cable
      CALL sf_expl_l_cable (                                                  &
        !IN date-related values
        curr_year, curr_day_number, curr_hour, curr_minute, curr_second,      &
        !IN values defining field dimensions and subset to be processed :
        land_pts, nice, nice_use,                                             &
        !IN  parameters for iterative SISL scheme
        numcycles, cycleno,                                                   &
        !IN parameters required from boundary-layer scheme :
        bq_1,bt_1,z1_uv_ij,z1_uv_top,z1_tq_ij,z1_tq_top,qw,tl,                &
        !IN soil/vegetation/land surface data :
        land_index,nsurft,sm_levels,canopy_surft,catch_surft,catch_snow_surft,&
        hcon_gb,ho2r2_orog_gb, fland,flandg,                                  &
        snow_surft,sil_orog_land_gb,smvccl_gb,smvcst_gb,smvcwt_gb,sthf_gb,    &
        sthu_gb,z0_surft, z0h_bare_surft, z0m_soil_gb,                        &
        !IN sea/sea-ice data :
        ice_fract_ncat_sicat,k_sice_sicat,                                    &
        !IN everything not covered so far :
        pstar,lw_down,sw_surft,zh,ddmfx,                                      &
        co2_mmr,co2_3d_ij,l_co2_interactive,l_phenol,l_triffid,               &
        asteps_since_triffid,cs_pool_gb,frac_surft,canht_pft,                 &
        photosynth_act_rad, lai_pft,                                          &
        lq_mix_bl,t_soil_gb,tsurf_elev_surft,                                 &
        ti_sicat,ti_cat_sicat,tstar,tstar_sea_ij,                             &
        tstar_sice_sicat,                                                     &
        tstar_surft,z_land_ij,albsoil_gb,cosz_gb,                             &
        l_aero_classic,l_dust,l_dust_diag,clay_gb,o3_gb,                      &
        !IN idealised and SCM things
        l_spec_z0, z0m_scm, z0h_scm,                                          &
        !IN variables for message passing
        u_1_p_ij, v_1_p_ij, u_0_p_ij, v_0_p_ij,                               &
        !INOUT diagnostics
        sf_diag,                                                              &
        !INOUT data :
        l_q10,z0msea_ij,gs_gb,g_leaf_acc_pft,npp_acc_pft,resp_w_acc_pft,      &
        resp_s_acc_gb,                                                        &
        !OUT Diagnostic not requiring STASH flags :
        recip_l_mo_sea,fqw_1,ftl_1,ftl_surft,                                 &
        radnet_sice,rhokm_1,rib,rib_surft,                                    &
        !OUT variables for message passing
        flandfac, fseafac, rhokm_land, rhokm_ssi,                             &
        cdr10m, cdr10m_n, cd10m_n,                                            &
        !OUT diagnostics required for soil moisture nudging scheme :
        wt_ext,                                                               &
        !OUT data required for tracer mixing :
        rho_aresist_ij,aresist_ij,resist_b_ij,                                &
        rho_aresist_surft,aresist_surft,resist_b_surft,                       &
        !OUT data required for mineral dust scheme
        r_b_dust_ij,cd_std_dust_ij,u_s_std_surft,                             &
        !OUT data required elsewhere in UM system :
        fb_surf,u_s,t1_sd,q1_sd,                                              &
        !OUT data required elsewhere in boundary layer or surface code
        alpha1,alpha1_sice,ashtf_prime,ashtf_prime_surft,fqw_surft,           &
        epot_surft,fqw_ice,ftl_ice,fraca,rhostar,resfs,resft,                 &
        rhokh,rhokh_surft,rhokh_sice,rhokh_sea,dtstar_surft,dtstar,           &
        h_blend_orog_ij,z0hssi,z0h_surft,z0h_eff,z0m_gb,z0mssi,z0m_surft,     &
        z0m_eff_ij,chr1p5m,chr1p5m_sice,smc_gb,hcons,vshr,vshr_land_ij,       &
        vshr_ssi_ij, gpp_gb,npp_gb,resp_p_gb,g_leaf_pft,gpp_pft,npp_pft,      &
        resp_p_pft,resp_s_gb,resp_s_tot,resp_l_pft,resp_r_pft,resp_w_pft,     &
        n_leaf_pft,n_root_pft,n_stem_pft,lai_bal_pft,                         &
        gc_surft,canhc_surft,wt_ext_surft,flake,                              &
        surft_index,surft_pts,tile_frac,fsmc,emis_surft,emis_soil,            &
!CABLE_LSM:
        mype, timestep_number, cycleno, numcycles,                            &
        true_latitude, true_longitude,                                        &
        bexp_gb, hcon_gb,satcon_gb, sathh_gb,                                 &
        soil_alb,                                                             & 
        surf_down_sw_cable, ls_rain_cable, ls_snow_cable,                     &
        cosz_gb, sin_theta_latitude                                           &
!CABLE_LSM: End
        )


    CASE DEFAULT
      CALL ereport('surf_couple_explicit', 101, 'Unrecognised surface scheme')

  END SELECT

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END SUBROUTINE surf_couple_explicit

END MODULE surf_couple_explicit_mod

