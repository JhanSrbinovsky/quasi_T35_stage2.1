! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE crop_vars_mod

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Module holding various variables for the crop model
!
! Current Code Owner: Tom Osborne
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

  REAL, ALLOCATABLE ::                                                        &
!-----------------------------------------------------------------------------
! Variables
!-----------------------------------------------------------------------------
    sow_date_cpft(:,:),                                                       &
        ! Sowing date of each crop functional type
    tt_veg_cpft(:,:),                                                         &
        ! Thermal requirement of stage 1 for crop pfts
    tt_rep_cpft(:,:),                                                         &
        ! Thermal requirement of stage 2 for crop pfts
    phot(:),                                                                  &
      !  Photoperiod (hours) for crop model
    dphotdt(:),                                                               &
      !  Rate of change of phot for crops
!-----------------------------------------------------------------------------
! Prognostics
!-----------------------------------------------------------------------------
    dvi_cpft(:,:),                                                            &
                !  Development index for crop tiles
    rootc_cpft(:,:),                                                          &
                !  Root carbon pool for crop tiles
    harvc_cpft(:,:),                                                          &
                !  Carbon in 'harvest parts' pool for crop tiles
    reservec_cpft(:,:),                                                       &
                !  carbon in stem reserves pool for crop tiles
    croplai_cpft(:,:),                                                        &
                !  Leaf area index for crop tiles
    cropcanht_cpft(:,:),                                                      &
                !  Canopy height for crop tiles
    sthu_irr_gb(:,:),                                                         &
                !  Unfrozen soil wetness over irrigation
!-----------------------------------------------------------------------------
! Diagnostics
!-----------------------------------------------------------------------------
    yield_diag_cpft(:,:),                                                     &
        ! Harvested carbon
    stemc_diag_cpft(:,:),                                                     &
        ! Stem carbon
    leafc_diag_cpft(:,:),                                                     &
        ! Leaf carbon
    irrig_water_gb(:),                                                        &
        ! Addition of irrigation water to soil (kg/m2/s).
    nonyield_diag_cpft(:,:)
        ! Carbon leaving the crop model which is not yield

  INTEGER, ALLOCATABLE ::                                                     &
    harvest_trigger_cpft(:,:),                                                &
        ! Trigger condition for the harvest in this timestep
    harvest_counter_cpft(:,:)
        ! 1 if timestep contains a harvest, 0 otherwise

!-----------------------------------------------------------------------------
! Irrigation variables
!-----------------------------------------------------------------------------
  REAL, ALLOCATABLE ::                                                        &
    frac_irr_all(:,:),                                                        &
    frac_irr_gb(:),                                                           &
        ! Irrigation fraction
    frac_irr_old_gb(:)
        ! Previous irrigated fraction in grid box

  LOGICAL :: frac_irr_all_tiles
        ! Switch to assign irrigation fraction to ALL tiles, or to pre-defined
        ! tiles ONLY
        ! set TRUE to reproduce original results

  INTEGER :: nirrtile
        ! Nr of tiles that can have irrigated fraction
  INTEGER, ALLOCATABLE :: irrtiles(:)
        ! Tiles that can have irrigated fraction
        ! Only used when frac_irr_all_tiles = .FALSE.

! Constants
  INTEGER, PARAMETER ::                                                       &
    NDPY = 365,                                                               &
        ! No. of days per year
    NYAV = 3,                                                                 &
        ! No. of years averaged for crop plant estimates
    nday_crop = 150
        ! Cropping date

  INTEGER :: iyear_old, startyr, startmon, startday, starttime
        ! Variables to store required datetime information

  INTEGER ,ALLOCATABLE ::                                                     &
    plant_n_gb(:),                                                            &
        ! best plant date for non-rice
    icntmax_gb(:)
        ! counter for start date for non-rice

! Variables for calculating average temperature, precip and radiation
  REAL, ALLOCATABLE ::                                                        &
    tl_1_day_av_gb(:),                                                        &
    tl_1_day_av_use_gb(:,:,:),                                                &
    prec_1_day_av_gb(:),                                                      &
    prec_1_day_av_use_gb(:,:,:),                                              &
    rn_1_day_av_gb(:),                                                        &
    rn_1_day_av_use_gb(:,:,:)

  REAL, ALLOCATABLE ::                                                        &
    frac_irr_surft(:,:),                                                      &
        ! Irrigation fraction in tile
    smc_irr_gb(:),                                                            &
        ! Available moisture in the soil profile in irrig frac (mm)
    wt_ext_irr_surft(:,:,:),                                                  &
        ! Fraction of transpiration over irrigated fraction which is
        ! extracted from each soil layer by each tile
    gs_irr_surft(:,:),                                                        &
        ! Conductance for irrigated surface types
    dvimax_gb(:),                                                             &
    gc_irr_surft(:,:),                                                        &
        ! Interactive canopy resistance
    resfs_irr_surft(:,:),                                                     &
        ! Combined soil, stomatal and aerodynamic resistance factor for
        ! fraction 1-FRACA
    ext_irr_gb(:,:),                                                          &
        ! Extraction of water from each soil layer over irrigation
    wt_ext_irr_gb(:,:),                                                       &
        ! Fraction of transpiration extracted from each soil layer
        ! over irrigated fraction.
    fsmc_irr_gb(:),                                                           &
        ! Soil moisture availability factor over irrigated fraction
    irrDaysDiag_gb(:)
        ! Number of days on which irrigation is applied

END MODULE crop_vars_mod
