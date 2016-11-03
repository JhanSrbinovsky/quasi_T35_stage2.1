! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE p_s_parms

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Module containing plant and soil variables (plus a few others)
!-----------------------------------------------------------------------------

  REAL, ALLOCATABLE ::                                                        &
    !Used in both the UM and JULES
    bexp_gb(:,:),                                                             &
      !  Exponent for soil moisture characteristic functions
      !    Clapp-Hornberger model: b is the Clapp-Hornberger exponent
      !    van Genuchten model: b=1/(n-1)  (metres)
    sathh_gb(:,:),                                                            &
      !  Parameter for soil moisture characteristic functions
      !  Clapp-Hornberger model: sathh is the saturated soil water pressure (m)
      !  van Genuchten model: sathh=1/alpha
    hcap_gb(:,:),                                                             &
      !  Soil heat capacity (J/K/m3)
    hcon_gb(:,:),                                                             &
      !  Soil thermal conductivity (W/m/K)
    satcon_gb(:,:),                                                           &
      !  Saturated hydraulic conductivity (kg/m2/s)
    smvccl_gb(:,:),                                                           &
      !  Critical volumetric SMC (cubic m per cubic m of soil)
    smvcst_gb(:,:),                                                           &
      !  Volumetric saturation point (m3/m3 of soil)
    smvcwt_gb(:,:),                                                           &
      !  Volumetric wilting point (cubic m per cubic m of soil)
    clay_gb(:)
      !  Fraction of clay

    !JULES only
#if !defined(UM_JULES)
  REAL, ALLOCATABLE ::                                                        &
    albsoil_gb(:),                                                            &
      !  Soil albedo
    albobs_sw_gb(:),                                                          &
      !  Obs SW albedo
    albobs_vis_gb(:),                                                         &
      !  Obs VIS albedo
    albobs_nir_gb(:),                                                         &
      !  Obs NIR albedo
    catch_surft(:,:),                                                         &
      !  Surface/canopy water capacity of snow-free land tiles (kg/m2)
    catch_snow_surft(:,:),                                                    &
      !  Snow interception capacity (kg/m2)
    cosz_gb(:),                                                               &
      !  Cosine of the zenith angle
    infil_surft(:,:),                                                         &
      !  Maximum possible surface infiltration for tiles (kg/m2/s)
    z0_surft(:,:),                                                            &
      !  Surface roughness on tiles (m).
    z0h_bare_surft(:,:),                                                      &
      !  Surface thermal roughness on tiles before allowance for snow
      !  cover (m).
    z0m_soil_gb(:),                                                           &
      !  Bare soil roughness, for momentum (m).
    sthu_gb(:,:),                                                             &
      ! Unfrozen soil moisture content of the layers as a fraction of saturation
    sthf_gb(:,:),                                                             &
      !  Frozen soil moisture content of the layers as a fraction of saturation.
    sthu_min_gb(:,:)
      ! Minimum unfrozen water content for each layer. Used to normalise
      ! thaw depth calculation based on unfrozen water content fraction.
#endif

END MODULE p_s_parms

