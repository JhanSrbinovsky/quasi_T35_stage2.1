#if !defined(UM_JULES)

SUBROUTINE init_vars_tmp()

  USE Ancil_info, ONLY: ice_fract_ncat_sicat, ti_cat_sicat, lice_pts, lice_index
  use jules_surface_mod
  USE Top_pdm, ONLY : inlandout_atm_gb
  USE Forcing, ONLY : u_0_ij,v_0_ij
  USE Prognostics
  USE Aero
  USE Orog
  USE Trifctl, ONLY: cv_gb, g_leaf_acc_pft, g_leaf_phen_acc_pft, npp_acc_pft, &
                     resp_s_acc_gb, resp_w_acc_pft, c_veg_pft
  USE Coastal
  USE C_Densty
  USE jules_sea_seaice_mod
  USE C_kappai
  USE C_Rough
  USE p_s_parms, ONLY : satcon_gb

  USE jules_sea_seaice_mod, ONLY : l_ssice_albedo

  USE update_mod, ONLY : l_imogen

  USE logging_mod, ONLY : log_info, log_debug, log_warn, log_error, log_fatal

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Initialises various variables that may change their initialisation in
!   future versions
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
  INTEGER :: i,l  ! Loop counters


!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
! Initialise accumulated fluxes for TRIFFID and phenology.
! This is not necessary if these are read from a restart file - but at
! present they're not.
!-----------------------------------------------------------------------
  g_leaf_acc_pft(:,:)       = 0.0
  g_leaf_phen_acc_pft(:,:)  = 0.0
  npp_acc_pft(:,:)          = 0.0
  resp_s_acc_gb(:,:)        = 0.0
  resp_w_acc_pft(:,:)       = 0.0

!-----------------------------------------------------------------------
! Set saturated hydraulic conductivity to zero at land ice points.
!-----------------------------------------------------------------------
  IF ( lice_pts > 0 ) THEN
    CALL log_info("init_vars_tmp",                                            &
                  "Setting satcon to zero at land ice points")
    DO i=1,lice_pts
      l = lice_index(i)
      satcon_gb(l,:) = 0.0
    ENDDO
  ENDIF

!-----------------------------------------------------------------------
! Set surface velocity to be zero
!-----------------------------------------------------------------------
  u_0_ij(:,:) = 0.0
  v_0_ij(:,:) = 0.0

!-----------------------------------------------------------------------
! Set CO2 variables
!-----------------------------------------------------------------------
  co2_3d_ij(:,:) = 0.0

!-----------------------------------------------------------------------
! Set coastal tiling variables
!-----------------------------------------------------------------------
  tstar_sea_ij(:,:)       = 280.0
  tstar_sice_sicat(:,:,:) = 270.0

!-----------------------------------------------------------------------
! Set orographic roughness variables
!-----------------------------------------------------------------------
  h_blend_orog_ij(:,:) = 0.0
  sil_orog_land_gb(:)  = 0.0
  ho2r2_orog_gb(:)     = 0.0

!-----------------------------------------------------------------------
! Set up prognostics which are not currently in dump
!-----------------------------------------------------------------------
  ti_cat_sicat(:,:,:)         = 270.0
  z0msea_ij(:,:)              = z0hsea
  snow_mass_ij(:,:)           = 0.0
  ice_fract_ncat_sicat(:,:,:) = 0.0
  di_ncat_sicat(:,:,:)        = 0.0
  snow_mass_sea_sicat(:,:,:)  = 0.0
  k_sice_sicat(:,:,:)         = 2.0 * kappai / de

!-----------------------------------------------------------------------
! Set up sea-ice parameter variables
!-----------------------------------------------------------------------
  IF(L_SSICE_ALBEDO) THEN
    alpham=0.65
    alphac=0.80
    alphab=0.57
    dtice=2.00
  ELSE
    alpham=0.50
    alphac=0.80
    alphab=-1.00
    dtice=10.00
  ENDIF

  soot_ij(:,:) = 0.0

!------------------------------------------------------------------
! Temporary initialisation of variables added during UM integration
!------------------------------------------------------------------

  inlandout_atm_gb(:) = 0.0

!-----------------------------------------------------------------------------
! Unless IMOGEN is on, cv will not have been initialised
! Unless TRIFFID is on, it will also not be used...
!-----------------------------------------------------------------------------
  IF ( .NOT. l_imogen ) cv_gb(:) = 0.0

! Initialise c_veg to 0 as not doing so can cause issues with IMOGEN outputs
  c_veg_pft(:,:) = 0.0

  RETURN

END SUBROUTINE init_vars_tmp
#endif
