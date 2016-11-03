! *****************************COPYRIGHT****************************************
! (c) Crown copyright, Met Office 2014. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms and
! conditions set out therein.
!
! [Met Office Ref SC0237]
! *****************************COPYRIGHT****************************************

MODULE surf_couple_radiation_mod

IMPLICIT NONE

PRIVATE
PUBLIC :: surf_couple_radiation

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SURF_COUPLE_RADIATION_MOD'

CONTAINS

!===============================================================================
! Public subroutine
!===============================================================================
SUBROUTINE surf_couple_radiation(                                             &
  !Fluxes INTENT(IN)
  tstar,                                                                      &
  !Misc INTENT(IN)
  ws10m, chloro,                                                              &
  n_band, max_n_swbands, wavelength_short, wavelength_long,                   &
  !Misc INTENT(OUT)
  sea_ice_albedo,                                                             &
  !Fluxes INTENT(OUT)
  alb_surft, land_albedo                                                      &
#if defined(UM_JULES)
  ,                                                                           &
  !UM-only args: INTENT(IN)
  pond_frac_cat_sicat, pond_depth_cat_sicat,                                  &
  !(ancil_info mod)
  nsurft, land_pts, land_index, surft_pts, surft_index,                       &
  row_length, rows, ice_fract_ij, ice_fract_ncat_sicat, frac_surft,           &
  !(p_s_parms mod)
  cosz_gb, albobs_sw_gb, albobs_vis_gb, albobs_nir_gb,                        &
  z0_surft, albsoil_gb,                                                       &
  !(coastal mod)
  flandg, tstar_sice_sicat,                                                   &
  !(prognostics mod)
  snow_mass_ij, snow_mass_sea_sicat, di_ncat_sicat, lai_pft, canht_pft,       &
  rgrain_surft, snow_surft, soot_ij, tstar_surft, ho2r2_orog_gb,              &
  !UM-only args: INTENT(OUT)
  albobs_sc, open_sea_albedo                                                  &
#endif
  )

!Module imports.

!Common modules
USE ereport_mod, ONLY : ereport

USE jules_sea_seaice_mod,     ONLY:                                           &
  nice, nice_use,                                                             &
  alpham, alphac, alphab, dtice, dt_bare, dalb_bare_wet,                      &
  pen_rad_frac, sw_beta,                                                      &
  albicev_cice, albicei_cice, albsnowv_cice, albsnowi_cice,                   &
  albpondv_cice, albpondi_cice,                                               &
  ahmax, dalb_mlt_cice, dalb_mlts_v_cice, dalb_mlts_i_cice,                   &
  dt_bare_cice, dt_snow_cice, pen_rad_frac_cice, sw_beta_cice,                &
  snowpatch
USE jules_surface_mod,        ONLY:                                           &
  l_aggregate
USE jules_surface_types_mod,  ONLY:                                           &
  ntype, npft
USE jules_radiation_mod,      ONLY:                                           &
  l_snow_albedo

!Potential troublemaker
USE theta_field_sizes,        ONLY:                                           &
  t_i_length, t_j_length


!JULES-only modules
#if !defined(UM_JULES)
  USE ancil_info,               ONLY:                                         &
    ice_fract_ij, ice_fract_ncat_sicat, pond_frac_cat_sicat,                  &
    pond_depth_cat_sicat, land_pts, land_index, nsurft, surft_pts,            &
    surft_index, frac_surft, row_length, rows
  USE p_s_parms,                ONLY:                                         &
    cosz_gb, albsoil_gb, albobs_sw_gb, albobs_vis_gb, albobs_nir_gb, z0_surft
  USE coastal,                  ONLY:                                         &
    flandg, tstar_sice_sicat
  USE orog,                     ONLY:                                         &
    ho2r2_orog_gb
  USE prognostics,              ONLY:                                         &
    snow_mass_ij, snow_mass_sea_sicat, di_ncat_sicat, lai_pft, canht_pft,     &
    rgrain_surft, snow_surft, soot_ij, tstar_surft
#endif

!CABLE_LSM:
USE lsm_switches_mod,        ONLY: lsm_id
USE um_parcore,              ONLY : mype
USE timestep_mod,            ONLY : timestep_number 

!Dr Hook
USE parkind1,                 ONLY:                                           &
  jprb, jpim

USE yomhook,                  ONLY:                                           &
  lhook, dr_hook

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Coupling routine between the UM or JULES system code and land surface
!   radiation science routines. Calls the appropriate LSM-specific code.
!
!   Some variables exist in modules only in JULES, others only in the UM
!   Options in order of preference
!   -UM and JULES share the same module names.
!   -UM and JULES have different module names and USE statements go on an ifdef
!   -The UM flavour of the variable does not live in a module. Pass in using a
!    an ifdef'ed argument list
!
!   If there are lots of ifs, then we could cosider splitting lsm_couple into
!   jules_couple and lsm_couple
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
  !UM-only args: INTENT(IN)
  !(ancil_info mod)
  INTEGER, INTENT(IN)::                                                       &
    nsurft, land_pts, land_index(land_pts), surft_pts(nsurft),                &
    surft_index(land_pts,nsurft), row_length, rows

  REAL, INTENT(IN) ::                                                         &
    ice_fract_ij(row_length, rows),                                           &
    ice_fract_ncat_sicat(row_length, rows,nice_use),                          &
    frac_surft(land_pts, ntype)

  !(p_s_parms mod)
  REAL, INTENT(IN) ::                                                         &
    cosz_gb(row_length, rows),                                                &
    albobs_sw_gb(land_pts),                                                   &
    albobs_vis_gb(land_pts),                                                  &
    albobs_nir_gb(land_pts),                                                  &
    z0_surft(land_pts,nsurft),                                                &
    albsoil_gb(land_pts)

  !(orog mod)
  REAL, INTENT(IN) ::                                                         &
    ho2r2_orog_gb(land_pts)

  !(coastal mod)
  REAL, INTENT(IN) ::                                                         &
  flandg(row_length, rows),                                                   &
  tstar_sice_sicat(row_length, rows)

  !(prognostics mod)
  REAL, INTENT(IN) ::                                                         &
    snow_mass_ij(row_length, rows),                                           &
    snow_mass_sea_sicat(row_length, rows, nice_use),                          &
    di_ncat_sicat(row_length, rows, nice_use),                                &
    canht_pft(land_pts, npft),                                                &
    lai_pft(land_pts, npft),                                                  &
    rgrain_surft(land_pts, nsurft),                                           &
    snow_surft(land_pts, nsurft),                                             &
    soot_ij(row_length, rows),                                                &
    tstar_surft(land_pts, nsurft)

  !UM-only args: INTENT(OUT)
  REAL, INTENT(OUT) ::                                                        &
    albobs_sc(t_i_length,t_j_length,ntype,2),                                 &
      !albedo scaling factors to obs
    open_sea_albedo(row_length,rows,2,max_n_swbands)
      !Surface albedo for Open Sea (direct and diffuse components, for each
      !band, with zeros for safety where no value applies)

   REAL, INTENT(IN) ::                                                        &
     pond_frac_cat_sicat(row_length, rows, nice),                             &
!     Meltpond fraction on sea ice categories
     pond_depth_cat_sicat(row_length, rows, nice)
!     Meltpond depth on sea ice categories (m)

#else
  REAL ::                                                                     &
      albobs_sc(t_i_length,t_j_length,ntype,2),                               &
        !albedo scaling factors to obs
      open_sea_albedo(row_length,rows,2,max_n_swbands)
        !Surface albedo for Open Sea (direct and diffuse components, for each
        !band, with zeros for safety where no value applies)
#endif

!Fluxes INTENT(IN)
REAL, INTENT(IN) ::                                                           &
  tstar(row_length,rows)            !Surface temperature

!Misc INTENT(IN)
REAL, INTENT(IN) ::                                                           &
  ws10m(row_length,rows),                                                     &
                                    !10m wind speed
  chloro(row_length,rows)           !nr surface chlorophyll content

INTEGER, INTENT(IN) ::                                                        &
  n_band,                                                                     &
  max_n_swbands

REAL, INTENT(IN)    ::                                                        &
  wavelength_short(n_band),                                                   &
  wavelength_long(n_band)

!Misc INTENT(OUT)
REAL, INTENT(OUT) ::                                                          &
  sea_ice_albedo(row_length,rows,4)   !Surface Albedo for sea ice
                                      ! (*,1) - direct beam visible
                                      ! (*,2) - diffuse visible
                                      ! (*,3) - direct beam near-ir
                                      ! (*,4) - diffuse near-ir

!Fluxes INTENT(OUT)
REAL, INTENT(OUT) ::                                                          &
  alb_surft(land_pts,nsurft,4),                                               &
                                      !Albedos for surface tiles.
                                      ! (*,*,1) - Direct beam visible
                                      ! (*,*,2) - Diffuse visible
                                      ! (*,*,3) - Direct beam near-IR
                                      ! (*,*,4) - Diffuse near-IR
  land_albedo(t_i_length,t_j_length,4)!GBM albedos.

!-----------------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------------

!Dr Hook variables
  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='SURF_COUPLE_RADIATION'

!-----------------------------------------------------------------------------
!End of header

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  SELECT CASE( lsm_id )
    CASE( 'jules' )
! DEPENDS ON: ftsa
      CALL ftsa (                                                             &
        !INTENT(IN)
        !input fields
        flandg, ice_fract_ij, ice_fract_ncat_sicat, tstar, tstar_sice_sicat,  &
        cosz_gb, ws10m, chloro,                                               &
        snow_mass_ij, snow_mass_sea_sicat, di_ncat_sicat,                     &
        pond_frac_cat_sicat, pond_depth_cat_sicat,                            &
        !max and min sea ice albedo specifications
        alpham, alphac, alphab, dtice,                                        &
        dt_bare, dalb_bare_wet, pen_rad_frac, sw_beta,                        &
        ! parameters for CICE multi-band albedo scheme:
        albicev_cice, albicei_cice, albsnowv_cice, albsnowi_cice,             &
        albpondv_cice, albpondi_cice,                                         &
        ahmax, dalb_mlt_cice, dalb_mlts_v_cice, dalb_mlts_i_cice,             &
        dt_bare_cice, dt_snow_cice,                                           &
        pen_rad_frac_cice, sw_beta_cice, snowpatch,                           &
        !size and control variables
        row_length*rows, max_n_swbands,                                       &
        n_band, nice, nice_use,                                               &
        !spectral boundaries
        wavelength_short,                                                     &
        wavelength_long,                                                      &
        !INTENT(OUT)
        !output arguments
        sea_ice_albedo,                                                       &
        open_sea_albedo)

! DEPENDS ON: tile_albedo
      CALL tile_albedo (                                                      &
        t_i_length*t_j_length, land_pts, land_index, nsurft, surft_pts,       &
        surft_index,                                                          &
        l_aggregate,                                                          &
        l_snow_albedo, albsoil_gb, albobs_sw_gb, albobs_vis_gb, albobs_nir_gb,&
        cosz_gb,frac_surft,                                                   &
        lai_pft,canht_pft,rgrain_surft, snow_surft, soot_ij, tstar_surft,     &
        z0_surft,ho2r2_orog_gb,                                               &
        alb_surft, land_albedo, albobs_sc)


    CASE( 'cable' )
! DEPENDS ON: ftsa
      CALL ftsa (                                                             &
        !INTENT(IN)
        !input fields
        flandg, ice_fract_ij, ice_fract_ncat_sicat, tstar, tstar_sice_sicat,  &
        cosz_gb, ws10m, chloro,                                               &
        snow_mass_ij, snow_mass_sea_sicat, di_ncat_sicat,                     &
        pond_frac_cat_sicat, pond_depth_cat_sicat,                            &
        !max and min sea ice albedo specifications
        alpham, alphac, alphab, dtice,                                        &
        dt_bare, dalb_bare_wet, pen_rad_frac, sw_beta,                        &
        ! parameters for CICE multi-band albedo scheme:
        albicev_cice, albicei_cice, albsnowv_cice, albsnowi_cice,             &
        albpondv_cice, albpondi_cice,                                         &
        ahmax, dalb_mlt_cice, dalb_mlts_v_cice, dalb_mlts_i_cice,             &
        dt_bare_cice, dt_snow_cice,                                           &
        pen_rad_frac_cice, sw_beta_cice, snowpatch,                           &
        !size and control variables
        row_length*rows, max_n_swbands,                                       &
        n_band, nice, nice_use,                                               &
        !spectral boundaries
        wavelength_short,                                                     &
        wavelength_long,                                                      &
        !INTENT(OUT)
        !output arguments
        sea_ice_albedo,                                                       &
        open_sea_albedo)

! DEPENDS ON: tile_albedo_cable
      CALL tile_albedo_cable (                                                &
        t_i_length*t_j_length, land_pts, land_index, nsurft, surft_pts,       &
        surft_index,                                                          &
        l_aggregate,                                                          &
        l_snow_albedo, albsoil_gb, albobs_sw_gb, albobs_vis_gb, albobs_nir_gb,&
        cosz_gb,frac_surft,                                                   &
        lai_pft,canht_pft,rgrain_surft, snow_surft, soot_ij, tstar_surft,     &
        z0_surft,ho2r2_orog_gb,                                               &
        alb_surft, land_albedo, albobs_sc,                                    &
        !CABLE_LSM:
        mype, timestep_number                                                 &  
        )


    CASE DEFAULT
      CALL ereport('surf_couple_radiation', 101, 'Unrecognised surface scheme')

  END SELECT

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END SUBROUTINE surf_couple_radiation
END MODULE surf_couple_radiation_mod

