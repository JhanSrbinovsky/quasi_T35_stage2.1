#if !defined(UM_JULES)
! Module containing the variables for Topmodel and PDM.

! DBC Arguably sthzw and zw should be stored in PROGNOSTICS since they
! are indeed prognostics. fsat needs to persist between timesteps (but
! can be initialised (recalculated) from soil moisture).


MODULE top_pdm

  IMPLICIT NONE

  REAL, DIMENSION(:), ALLOCATABLE :: fexp_gb
    ! Decay factor in Sat. Conductivity in water table layer
  REAL, DIMENSION(:), ALLOCATABLE :: gamtot_gb
    ! Integrated complete Gamma function
    ! DBC gamtot doesn't need to be in a module in this version, but left there
    !for now for compatability.
  REAL, DIMENSION(:), ALLOCATABLE :: ti_mean_gb
    ! Mean topographic index
  REAL, DIMENSION(:), ALLOCATABLE :: ti_sig_gb
    ! Standard dev. of topographic index
  REAL, DIMENSION(:), ALLOCATABLE :: fsat_gb
    ! Surface saturation fraction
  REAL, DIMENSION(:), ALLOCATABLE :: fwetl_gb
    ! Wetland fraction
  REAL, DIMENSION(:), ALLOCATABLE :: zw_gb
    ! Water table depth (m)
  REAL, DIMENSION(:), ALLOCATABLE :: drain_gb
    ! Drainage out of bottom (nshyd) soil layer (kg/m2/s)
  REAL, DIMENSION(:), ALLOCATABLE :: dun_roff_gb
    ! Dunne part of sfc runoff (kg/m2/s)
  REAL, DIMENSION(:), ALLOCATABLE :: qbase_gb
    ! Base flow (kg/m2/s)
  REAL, DIMENSION(:), ALLOCATABLE :: qbase_zw_gb
    ! Base flow from ZW layer (kg/m2/s)
  REAL, DIMENSION(:), ALLOCATABLE :: fch4_wetl_gb
    ! Scaled wetland methane flux with
    ! substrate used in atmospheric chemistry.
    ! (Set by switch l_wetland_ch4_npp:
    ! =FALSE (default) FCH4_WETL=FCH4_WETL_CS
    ! =TRUE FCH4_WETL=FCH4_WETL_NPP)
    ! (Note different units: 10^-9 kg C/m2/s).
  REAL, DIMENSION(:), ALLOCATABLE :: fch4_wetl_cs_gb
    ! Scaled wetland methane flux using
    ! soil carbon as substrate (kg C/m2/s).
  REAL, DIMENSION(:), ALLOCATABLE :: fch4_wetl_npp_gb
    ! Scaled wetland methane flux using
    ! NPP as substrate (kg C/m2/s).
  REAL, DIMENSION(:), ALLOCATABLE :: fch4_wetl_resps_gb
    ! Scaled wetland methane flux using
    ! soil respiration as substrate (kg C/m2/s).

  REAL, ALLOCATABLE :: inlandout_atm_gb(:)
    ! TRIP inland basin outflow (for land points only)(kg/m2/s)
  REAL, ALLOCATABLE :: sthzw_gb(:)
    ! soil moist fraction in deep (water table) layer.
  REAL, ALLOCATABLE :: a_fsat_gb(:)
    ! Fitting parameter for Fsat in LSH model
  REAL, ALLOCATABLE :: c_fsat_gb(:)
    ! Fitting parameter for Fsat in LSH model
  REAL, ALLOCATABLE :: a_fwet_gb(:)
    ! Fitting parameter for Fwet in LSH model
  REAL, ALLOCATABLE :: c_fwet_gb(:)
    ! Fitting parameter for Fwet in LSH model

END MODULE top_pdm
#endif
