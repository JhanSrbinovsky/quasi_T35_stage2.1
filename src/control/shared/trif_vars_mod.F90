! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE trif_vars_mod

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Module holding various shared variables for TRIFFID
!
! Current Code Owner: Andy Wiltshire
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Diagnostics
!-----------------------------------------------------------------------------

  REAL, ALLOCATABLE :: wp_fast_in_gb(:)
                        ! INOUT Fast-turnover wood product C pool input.
  REAL, ALLOCATABLE :: wp_med_in_gb(:)
                        ! INOUT Fast-turnover wood product C pool input.
  REAL, ALLOCATABLE :: wp_slow_in_gb(:)
                        ! INOUT Fast-turnover wood product C pool input.
  REAL, ALLOCATABLE :: wp_fast_out_gb(:)
                        ! INOUT Fast-turnover wood product C pool output.
  REAL, ALLOCATABLE :: wp_med_out_gb(:)
                        ! INOUT Fast-turnover wood product C pool output.
  REAL, ALLOCATABLE :: wp_slow_out_gb(:)
                        ! INOUT Fast-turnover wood product C pool output.
  REAL, ALLOCATABLE :: lit_c_orig_pft(:,:)
  REAL, ALLOCATABLE :: lit_c_ag_pft(:,:)
  REAL, ALLOCATABLE :: lit_n_orig_pft(:,:)
  REAL, ALLOCATABLE :: lit_n_ag_pft(:,:)

  REAL,   ALLOCATABLE :: root_abandon_pft(:,:)
                                  ! Root carbon moved to soil carbon during
                                  ! landuse change.
                                  !    (kg C/(m2 PFT)/yr)
  REAL,   ALLOCATABLE :: harvest_pft(:,:)
                                  ! Carbon harvested from crops.
                                  !    (kg C/(m2 PFT)/yr)
  REAL,   ALLOCATABLE :: harvest_gb(:)
                                  ! Carbon harvested from crops - gridbox mean.
                                  !    (kg C/(m2 land)/yr)

  REAL,   ALLOCATABLE :: root_abandon_n_pft(:,:)
                                  ! Root nitrogen moved to soil nitrogen during
                                  ! landuse change.
                                  !    (kg N/(m2 PFT)/yr)
  REAL,   ALLOCATABLE :: harvest_n_pft(:,:)
                                  ! Nitrogen harvested from crops.
                                  !    (kg N/(m2 PFT)/yr)
  REAL,   ALLOCATABLE :: harvest_n_gb(:)
                                  ! Nitrogen harvested from crops - gridbox mean
                                  !    (kg N/(m2 land)/yr)
  REAL,   ALLOCATABLE :: n_fertiliser_pft(:,:)
                                  ! Nitrogen available to crop PFTs in addition
                                  ! to soil nitrogen (kg N/(m2 PFT)/yr)
  REAL,   ALLOCATABLE :: n_fertiliser_gb(:)
                                  ! N avaialable to crop PFTs in addition to
                                  ! soil N - gridbox mean (kg N/(m2 land)/yr)


  REAL, ALLOCATABLE :: n_leaf_pft(:,:)
                        ! Leaf N content scaled by LAI in sf_stom
                        ! (kgN/m2)
  REAL, ALLOCATABLE :: n_root_pft(:,:)
                        ! Root N content scaled by LAI_BAL in sf_stom
                        ! (kgN/m2)
  REAL, ALLOCATABLE :: n_stem_pft(:,:)
                        ! Stem N content scaled by LAI_BAL in sf_stom
                        ! (kgN/m2)
  REAL, ALLOCATABLE :: resp_r_pft(:,:)
                        ! Root maintenance respiration (kg C/m2/s)
  REAL, ALLOCATABLE :: resp_l_pft(:,:)
                        ! Leaf maintenance respiration (kg C/m2/s)
  REAL, ALLOCATABLE :: lai_bal_pft(:,:)
                        ! balanced lai from sf_stom
  REAL,   ALLOCATABLE :: frac_past_gb(:)
!                                    ! Fraction of pasture
  REAL, ALLOCATABLE :: pc_s_pft(:,:)
                        ! Carbon available for spreading a PFT.
                        !    (kg C/m2/yr)
  REAL, ALLOCATABLE :: fapar_diag_pft(:,:)
                        ! Fraction of Absorbed Photosynthetically Active
                        ! Radiation diagnostic
  REAL, ALLOCATABLE :: fao_et0(:)
                        ! FAO Penman-Monteith evapotranspiration for 
                        ! reference crop in kg m-2 s-1
  REAL,   ALLOCATABLE :: cnsrv_carbon_veg2_gb(:)
                             ! IN Diagnostic of error in land carbon
                             ! conservation in the veg2 routine
  REAL,   ALLOCATABLE :: cnsrv_carbon_triffid_gb(:)
                             ! IN Diagnostic of error in land carbon
                             ! conservation in the triffid routine
  REAL,   ALLOCATABLE :: cnsrv_veg_triffid_gb(:)
                             ! IN Diagnostic of error in vegetation carbon
                             !    conservation in the triffid routine
  REAL,   ALLOCATABLE :: cnsrv_soil_triffid_gb(:)
                             ! IN Diagnostic of error in soil carbon
                             !    conservation in the triffid routine
  REAL,   ALLOCATABLE :: cnsrv_prod_triffid_gb(:)
                             ! IN Diagnostic of error in wood product carbon
                             !    conservation in the triffid routine

!-------------------------------------------------------------------------------
! Variables added for ticket #7,#127 (nitrogen scheme)
!-------------------------------------------------------------------------------
  REAL, ALLOCATABLE :: deposition_n_gb(:)
                                  ! Nitrogen deposition (kgN/m2/s)

  REAL, ALLOCATABLE :: leafC_pft(:,:)
                        ! Leaf Carbon on PFTs (kgC/m2)
  REAL, ALLOCATABLE :: rootC_pft(:,:)
                        ! Root Carbon on PFTs (kgC/m2)
  REAL, ALLOCATABLE :: stemC_pft(:,:)
                        ! Stem Carbon on PFTs (kgC/m2)
  REAL, ALLOCATABLE :: woodC_pft(:,:)
                        ! Wood Carbon on PFTs (kgC/m2)
  REAL, ALLOCATABLE :: droot_pft(:,:)
                        ! Increment in Leaf Carbon on PFTs (kgC/m2/trif period)
  REAL, ALLOCATABLE :: dleaf_pft(:,:)
                        ! Increment in Leaf Carbon on PFTs (kgC/m2/trif period)
  REAL, ALLOCATABLE :: dwood_pft(:,:)
                        ! Increment in Wood Carbon on PFTs (kgC/m2/trif period)

  REAL, ALLOCATABLE :: root_litC_pft(:,:)
                        ! Root Litter C Turnover on PFTS (kgC/m2/360 day)
  REAL, ALLOCATABLE :: leaf_litC_pft(:,:)
                        ! Leaf Litter C Turnover on PFTS (kgC/m2/360 day)
  REAL, ALLOCATABLE :: wood_litC_pft(:,:)
                        ! Wood Litter C Turnover on PFTS (kgC/m2/360 day)

  REAL, ALLOCATABLE :: root_litN_pft(:,:)
                        ! Root Litter N Turnover on PFTS (kgN/m2/360 day)
  REAL, ALLOCATABLE :: leaf_litN_pft(:,:)
                        ! Leaf Litter N Turnover on PFTS (kgN/m2/360 day)
  REAL, ALLOCATABLE :: wood_litN_pft(:,:)
                        ! Wood Litter N Turnover on PFTS (kgN/m2/360 day)


  REAL, ALLOCATABLE :: litterC_pft(:,:)
                        ! Local litter production (kgC/m2/360 day)
  REAL, ALLOCATABLE :: litterN_pft(:,:)
                        ! Local litter production (kgN/m2/360 day)

  REAL, ALLOCATABLE :: lit_n_pft(:,:)
                        ! Total litter production (kgN/m2/360 day)

  REAL, ALLOCATABLE :: n_uptake_growth_pft(:,:)
                        ! Vegetation N uptake for growth on PFTS
                        !(kgN/m2/360 days)
  REAL, ALLOCATABLE :: n_demand_growth_pft(:,:)
                        ! Vegetation N demand for growth on PFTS
                        !(kgN/m2/360 days)
  REAL, ALLOCATABLE :: n_demand_lit_pft(:,:)
                        ! Vegetation N demand for balanced literature
                        !production on PFTS (kgN/m2/360 days)
  REAL, ALLOCATABLE :: n_demand_spread_pft(:,:)
                        ! Vegetation N demand for spreading on
                        !PFTS(kgN/m2/360 days)
  REAL, ALLOCATABLE :: n_uptake_spread_pft(:,:)
                        ! Vegetation N uptake for spreading in PFTS
                        !(kgN/m2/360 days)
  REAL, ALLOCATABLE :: n_uptake_pft(:,:)
                        ! Vegetation N uptake on PFTS (kgN/m2/360 days)
  REAL, ALLOCATABLE :: n_demand_pft(:,:)
                        ! Vegetation N demand on PFTS(kgN/m2/360 days)
  REAL, ALLOCATABLE :: n_uptake_gb(:)
                        ! Vegetation N uptake (kgN/m2/360 days)
  REAL, ALLOCATABLE :: n_demand_gb(:)
                        ! Vegetation N demand (kgN/m2/360 days)
  REAL, ALLOCATABLE :: n_veg_pft(:,:)
                        ! Veg N on PFTS (kgN/m2)
  REAL, ALLOCATABLE :: n_veg_gb(:)
                        ! Veg N (kgN/m2)
  REAL, ALLOCATABLE :: n_loss_gb(:)
                        ! Other N Gaseous Loss (kgN/m2/360 days)
  REAL, ALLOCATABLE :: n_fix_pft(:,:)
                        ! Fixed N on PFTS (kgN/m2/360 days)
  REAL, ALLOCATABLE :: n_fix_gb(:)
                        ! Fixed N (kgN/m2/360 days)
  REAL, ALLOCATABLE :: n_leach_gb(:)
                        ! Leached N (kgN/m2/s)
  REAL, ALLOCATABLE :: n_gas_gb(:)
                        ! Mineralised N Gas Emmissions (kgN/m2/360 days)
  REAL, ALLOCATABLE :: dpm_ratio_gb(:)
                        ! Ratio of Decomposatble Plant Material to Resistant (:)
  REAL, ALLOCATABLE :: dnveg_pft(:,:)
                        ! Increment in Veg N on PFTs (kgN/m2/trif period)
  REAL, ALLOCATABLE :: dcveg_pft(:,:)
                        ! Increment in Veg C on PFTs (kgC/m2/trif period)
  REAL, ALLOCATABLE :: dnveg_gb(:)
                        ! Increment in Veg N (kgN/m2/trif period)
  REAL, ALLOCATABLE :: dcveg_gb(:)
                        ! Increment in Veg C (kgC/m2/trif period)

  REAL, ALLOCATABLE :: immob_n_gb(:,:)
                        ! Immobilised N on soil pools (kgN/m2/360 days)
  REAL, ALLOCATABLE :: immob_n_pot_gb(:,:)
                        ! Unlimited Immobilised N on soil pools
                        !(kgN/m2/360 days)

  REAL, ALLOCATABLE :: minl_n_gb(:,:)
                        ! Mineralised N on soil pools (kgN/m2/360 days)
  REAL, ALLOCATABLE :: minl_n_pot_gb(:,:)
                        ! Unlimited Mineralised N on soil pools
                        !(kgN/m2/360 days)

  REAL, ALLOCATABLE :: resp_s_diag_gb(:,:)
                        ! Diagnosed soil respiration after TRIFFID on soil pools
                        !(kgC/m2/360 days)
  REAL, ALLOCATABLE :: resp_s_pot_diag_gb(:,:)
                        ! Diagnosed unlimited soil respiration after TRIFFID
                        !(kgC/m2/360 days)

  REAL, ALLOCATABLE :: fn_gb(:)
                        ! Nitrogen decomposition factor (:)

  REAL, ALLOCATABLE :: lit_n_t_gb(:)
                        ! Total N litter flux (kgC/m2/360 days)
  REAL, ALLOCATABLE :: nl_bal_gb(:)
                        ! Canopy N concentration at balanced LAI (kgN/KgC)
  REAL, ALLOCATABLE :: exudates_pft(:,:)
                        ! Unallocated C due to lack on N added to soil
                        !respiration on PFTS ((kgC/m2/360 days)
  REAL, ALLOCATABLE :: exudates_gb(:)
                        ! Unallocated C due to lack on N added to soil
                        !respiration ((kgC/m2/360 days)

END MODULE trif_vars_mod
