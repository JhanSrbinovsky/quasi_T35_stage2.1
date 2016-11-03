! *****************************COPYRIGHT****************************************
! (c) Crown copyright, Met Office 2016. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms and
! conditions set out therein.
!
! [Met Office Ref SC0237]
! *****************************COPYRIGHT****************************************
!
! Code Description:
!   Language: FORTRAN 90
!
! Current Code Owner: Stephane Mangeon
!

MODULE inferno_io_mod

USE parkind1,                       ONLY : jpim

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = "INFERNO_IO_MOD"

CONTAINS

SUBROUTINE inferno_io( t1p5m_tile, q1p5m_tile, pstar, sthu,                   &
  sm_levels, frac, dim_cs1, cs, canht, ls_rain, con_rain,                     &
  land_pts, ignition_method)

USE inferno_mod,  ONLY :                                                      &
  calc_flam, calc_ignitions, calc_burnt_area,                                 &
  calc_emitted_carbon, calc_emitted_carbon_soil, calc_emission

USE yomhook,      ONLY : lhook, dr_hook

USE fire_vars,                ONLY:                                           &
  burnt_area, burnt_area_ft,                                                  &
  emitted_carbon, emitted_carbon_ft,                                          &
  emitted_carbon_DPM, emitted_carbon_RPM,                                     &
  fire_em_CO2, fire_em_CO2_ft,                                                &
  fire_em_CO2_DPM, fire_em_CO2_RPM,                                           &
  fire_em_CO, fire_em_CO_ft,                                                  &
  fire_em_CO_DPM, fire_em_CO_RPM,                                             &
  fire_em_CH4, fire_em_CH4_ft,                                                &
  fire_em_CH4_DPM, fire_em_CH4_RPM,                                           &
  fire_em_NOx, fire_em_NOx_ft,                                                &
  fire_em_NOx_DPM, fire_em_NOx_RPM,                                           &
  fire_em_SO2, fire_em_SO2_ft,                                                &
  fire_em_SO2_DPM, fire_em_SO2_RPM,                                           &
  fire_em_OC,  fire_em_OC_ft,                                                 &
  fire_em_OC_DPM,  fire_em_OC_RPM,                                            &
  fire_em_BC, fire_em_BC_ft,                                                  &
  fire_em_BC_DPM, fire_em_BC_RPM,                                             &
  pop_den, flash_rate

USE pftparm,                  ONLY:                                           &
  fef_co2, fef_co, fef_ch4, fef_nox, fef_so2, fef_oc, fef_bc,                 &
    ! PFT Emission factors from namelists
  ccleaf_min, ccleaf_max, ccwood_min, ccwood_max,                             &
    ! PFT Combustion completeness from namelists
  avg_ba
    ! Average Burned Area per PFT

USE jules_surface_types_mod,        ONLY : npft
USE ancil_info,                     ONLY : nsurft
USE parkind1,                       ONLY : jprb, jpim

! To compute leaf, root and wood carbon
USE pftparm,          ONLY : sigl, a_ws, eta_sl, a_wl, b_wl

IMPLICIT NONE
!
! Description:
!   Called every model timestep, this subroutine updates INFERNO's
!   driving variables and calls the scientific routines.
!
! Current Code Owner: Stephane Mangeon, Imperial College London
!
! Code Description:
!   Language: Fortran 90.
!

INTEGER, INTENT(IN) ::                                                        &
  land_pts, sm_levels, dim_cs1, ignition_method

REAL,    INTENT(IN) ::                                                        &
  t1p5m_tile(land_pts,nsurft),                                                &
  q1p5m_tile(land_pts,nsurft),                                                &
  pstar(land_pts),                                                            &
  sthu(land_pts,sm_levels),                                                   &
  frac(land_pts,nsurft),                                                      &
  cs(land_pts,dim_cs1),                                                       &
  canht(land_pts,npft),                                                       &
  ls_rain(land_pts),                                                          &
  con_rain(land_pts)

  ! Local temporary variables used in the interactive fire code
REAL                        ::                                                &
  inferno_temp(land_pts),                                                     &
    ! The temperature (K)
  inferno_rhum(land_pts),                                                     &
    ! The Relative Humidity (%)
  inferno_sm(land_pts),                                                       &
    ! The Soil Moisture (Fraction of saturation)
  inferno_rain(land_pts),                                                     &
    ! The total rainfall (kg/m2/s)
  inferno_fuel(land_pts),                                                     &
    ! The fuel density (fine litter and leaves - kg/m3)
  qsat(land_pts),                                                             &
    ! Saturation humidity
  ignitions(land_pts),                                                        &
    ! The number of ignitions (#/m2/s)
  lai_bal_inf(land_pts,npft),                                                 &
    ! The balanced lai used to compute carbon pools
  leaf_inf(land_pts,npft),                                                    &
    ! The leaf carbon
  wood_inf(land_pts,npft),                                                    &
    ! The wood carbon
  dpm_fuel(land_pts),                                                         &
    ! The amount of DPM that is available to burn (kgC.m-2)
  rpm_fuel(land_pts),                                                         &
    ! The amount of RPM that is available to burn (kgC.m-2)
  flam(land_pts)
    ! The Flammability

REAL ,   PARAMETER      ::                                                    &
  fef_co2_dpm = 1637.0, fef_co_dpm  = 89.0,                                   &
  fef_ch4_dpm = 3.92, fef_nox_dpm = 2.51,                                     &
  fef_so2_dpm = 0.40,                                                         &
  fef_oc_dpm  = 8.2,   fef_bc_dpm  = 0.56,                                    &
    ! HARDCODED Emission factors for DPM in g kg-1
  fef_co2_rpm = 1489.0, fef_co_rpm  = 127.0,                                  &
  fef_ch4_rpm = 5.96,  fef_nox_rpm = 0.90,                                    &
  fef_so2_rpm = 0.40,                                                         &
  fef_oc_rpm  = 8.2,   fef_bc_rpm  = 0.56,                                    &
      ! HARDCODED Emission factors for RPM in g kg-1
  pmtofuel    = 0.7,                                                          &
    ! Plant Material that is available as fuel (on the surface)
  fuel_low    = 0.02,  fuel_high   = 0.2
    !Fuel availability high/low threshold

INTEGER :: i, l ! counters for loops

REAL(KIND=jprb)               :: zhook_handle
CHARACTER (LEN=*),  PARAMETER :: RoutineName = "INFERNO_IO"

EXTERNAL qsat_wat_jls ! Saturation vapour pressure routine

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-------------------------------------------------
! Initialisation
!-------------------------------------------------
! Driving variables
inferno_temp(:) = 0.0
inferno_rhum(:) = 0.0
inferno_sm(:)   = 0.0
inferno_rain(:) = 0.0
inferno_fuel(:) = 0.0
! Work variables
qsat(:)         = 0.0
flam(:)         = 0.0
lai_bal_inf(:,:)= 0.0
leaf_inf(:,:)   = 0.0
wood_inf(:,:)   = 0.0
ignitions(:)    = 0.0

! INFERNO diagnostic variables
burnt_area(:)         = 0.0
burnt_area_ft(:,:)    = 0.0
emitted_carbon(:)     = 0.0
emitted_carbon_ft(:,:)= 0.0
emitted_carbon_DPM(:) = 0.0
emitted_carbon_RPM(:) = 0.0
fire_em_CO2(:)        = 0.0
fire_em_CO2_ft(:,:)   = 0.0
fire_em_CO2_DPM(:)    = 0.0
fire_em_CO2_RPM(:)    = 0.0
fire_em_CO(:)         = 0.0
fire_em_CO_ft(:,:)    = 0.0
fire_em_CO_DPM(:)     = 0.0
fire_em_CO_RPM(:)     = 0.0
fire_em_CH4(:)        = 0.0
fire_em_CH4_ft(:,:)   = 0.0
fire_em_CH4_DPM(:)    = 0.0
fire_em_CH4_RPM(:)    = 0.0
fire_em_NOx(:)        = 0.0
fire_em_NOx_ft(:,:)   = 0.0
fire_em_NOx_DPM(:)    = 0.0
fire_em_NOx_RPM(:)    = 0.0
fire_em_SO2(:)        = 0.0
fire_em_SO2_ft(:,:)   = 0.0
fire_em_SO2_DPM(:)    = 0.0
fire_em_SO2_RPM(:)    = 0.0
fire_em_OC(:)         = 0.0
fire_em_OC_ft(:,:)    = 0.0
fire_em_OC_DPM(:)     = 0.0
fire_em_OC_RPM(:)     = 0.0
fire_em_BC(:)         = 0.0
fire_em_BC_ft(:,:)    = 0.0
fire_em_BC_DPM(:)     = 0.0
fire_em_BC_RPM(:)     = 0.0

!---------------------------------------------------------------
! Get the available DPM and RPM using a scaling parameter
!---------------------------------------------------------------
dpm_fuel(:) = pmtofuel * cs(:,1)
rpm_fuel(:) = pmtofuel * cs(:,2)

!---------------------------------------------------------------
! Get the inferno meteorological variables for the whole gridbox
!---------------------------------------------------------------

! Soil Humidity (inferno_sm)
inferno_sm(:)  = (sthu(:,1))

! Rainfall (inferno_rain)
inferno_rain(:)  = ls_rain(:) + con_rain(:)

!---------------------------------------------------------------
! Diagnose the balanced-growth leaf area index and the leaf,
! wood, root and total vegetation carbon.
!---------------------------------------------------------------
DO i = 1, npft
  DO l = 1, land_pts
    lai_bal_inf(l,i) = MAX((a_ws(i)*eta_sl(i)*canht(l,i)                      &
                        /a_wl(i))**(1.0/(b_wl(i)-1.0)),0.0)

    leaf_inf(l,i)    = MAX(sigl(i)*lai_bal_inf(l,i),0.0)

    wood_inf(l,i)    = MAX(a_wl(i)*(lai_bal_inf(l,i)**b_wl(i)),0.0)
  END DO
END DO

!---------------------------------------------------------------
! Fire calculations - per PFT
!---------------------------------------------------------------
DO i = 1, npft
  ! Calculate the fuel density
  ! We use normalised Leaf Carbon + the available DPM
  inferno_fuel(:) = (leaf_inf(:,i)+dpm_fuel-fuel_low)                         &
                    /(fuel_high-fuel_low)

  WHERE(inferno_fuel < 0.0) inferno_fuel = 0.0

  WHERE(inferno_fuel > 1.0) inferno_fuel = 1.0

  inferno_temp(:) = t1p5m_tile(:,i)

  DO l = 1, land_pts
    !------------------------------------------------------------
    ! Conditional statements to make sure we are dealing with
    ! reasonable weather. Note initialisation to 0 already done.
    ! If the driving variables are singularities, we assume
    ! no burnt area.
    !------------------------------------------------------------

    ! Temperatures constrained akin to qsat (from the WMO)
    IF ((inferno_temp(l)>338.15).OR.(inferno_temp(l)<183.15)) CYCLE

    ! The maximum rain rate ever observed is 38mm in one minute,
    ! here we assume 0.5mm/s stops fires altogether
    IF ((inferno_rain(l)>0.5)   .OR.(inferno_rain(l)<0.0   )) CYCLE

    ! Fuel Density is an index constrained to 0-1
    IF ((inferno_fuel(l)>1.0)   .OR.(inferno_fuel(l)<0.0   )) CYCLE

    ! Soil moisture is a fraction of saturation
    IF ((inferno_sm(l)  >1.0)   .OR.(inferno_sm(l)  <0.0   )) CYCLE

    ! Get the tile relative humidity using saturation routine
    CALL qsat_wat(qsat(l), inferno_temp(l), pstar(l), 1)

    inferno_rhum(l)  = q1p5m_tile(l,i) / qsat(l) * 100.0

    ! Relative Humidity should be constrained to 0-100
    IF ((inferno_rhum(l)>100.0) .OR.(inferno_rhum(l)<0.0   )) CYCLE

    ! If all these checks are passes, start fire calculations
    CALL calc_ignitions(                                                      &
      !Point intent(IN)
      pop_den(l), flash_rate(l), ignition_method,                             &
      !Point intent(OUT)
      ignitions(l))

    CALL calc_flam(                                                           &
      !Point Intent(IN)
      inferno_temp(l), inferno_rhum(l), inferno_fuel(l),                      &
      inferno_sm(l), inferno_rain(l),                                         &
      !Point Intent(INOUT)
      flam(l))

    CALL calc_burnt_area(                                                     &
      !Point INTENT(IN)
      flam(l), ignitions(l), avg_ba(i),                                       &
      !Point INTENT(OUT)
      burnt_area_ft(l,i))
  END DO

  CALL calc_emitted_carbon(                                                   &
    !Array INTENT(IN)
    land_pts, burnt_area_ft(:,i), inferno_sm(:),                              &
    leaf_inf(:,i), wood_inf(:,i),                                             &
    ccleaf_min(i), ccleaf_max(i),                                             &
    ccwood_min(i), ccwood_max(i),                                             &
    !Array INTENT(OUT)
    emitted_carbon_ft(:,i))

  CALL calc_emission(                                                         &
    !Array INTENT(IN)
    land_pts, emitted_carbon_ft(:,i),                                         &
    fef_co2(i), fef_co(i), fef_ch4(i), fef_nox(i), fef_so2(i),                &
    fef_oc(i), fef_bc(i),                                                     &
    !Array INTENT(OUT)
    fire_em_CO2_ft(:,i), fire_em_CO_ft(:,i),                                  &
    fire_em_CH4_ft(:,i), fire_em_NOx_ft(:,i),                                 &
    fire_em_SO2_ft(:,i),                                                      &
    fire_em_OC_ft(:,i), fire_em_BC_ft(:,i) )

  ! We add pft-specific variables to the gridbox totals
  burnt_area(:)     = burnt_area(:)                                           &
                      + frac(:,i) * burnt_area_ft(:,i)
  emitted_carbon(:) = emitted_carbon(:)                                       &
                      + frac(:,i) * emitted_carbon_ft(:,i)
  fire_em_CO2(:)    = fire_em_CO2(:)                                          &
                      + frac(:,i) * fire_em_CO2_ft(:,i)
  fire_em_CO(:)     = fire_em_CO(:)                                           &
                      + frac(:,i) * fire_em_CO_ft(:,i)
  fire_em_CH4(:)    = fire_em_CH4(:)                                          &
                      + frac(:,i) * fire_em_CH4_ft(:,i)
  fire_em_NOx(:)    = fire_em_NOx(:)                                          &
                      + frac(:,i) * fire_em_NOx_ft(:,i)
  fire_em_SO2(:)    = fire_em_SO2(:)                                          &
                      + frac(:,i) * fire_em_SO2_ft(:,i)
  fire_em_OC(:)     = fire_em_OC(:)                                           &
                      + frac(:,i) * fire_em_OC_ft(:,i)
  fire_em_BC(:)     = fire_em_BC(:)                                           &
                      + frac(:,i) * fire_em_BC_ft(:,i)
END DO

!---------------------------------------------------------------
! In addition we diagnose the soil carbon (DPM and RPM only).
! However, this is not added to the gridbox totals as it was
! observed to lead to unrealistic emission with the currently
! recommended configuration.
!---------------------------------------------------------------

CALL calc_emitted_carbon_soil(                                                &
  ! Array INTENT(IN)
  land_pts, burnt_area, dpm_fuel, rpm_fuel,                                   &
  inferno_sm, dim_cs1,                                                        &
  ! Array INTENT(OUT)
  emitted_carbon_DPM, emitted_carbon_RPM )

! Decomposable Plant Material
CALL calc_emission(                                                           &
  !Array INTENT(IN)
  land_pts, emitted_carbon_DPM(:),                                            &
  fef_co2_dpm, fef_co_dpm, fef_ch4_dpm,                                       &
  fef_nox_dpm, fef_so2_dpm, fef_oc_dpm, fef_bc_dpm,                           &
  !Array INTENT(OUT)
  fire_em_CO2_DPM(:), fire_em_CO_DPM(:), fire_em_CH4_DPM(:),                  &
  fire_em_NOx_DPM(:), fire_em_SO2_DPM(:), fire_em_OC_DPM(:),                  &
  fire_em_BC_DPM(:))

! Resistant Plant Material
CALL calc_emission(                                                           &
  !Array INTENT(IN)
  land_pts, emitted_carbon_RPM(:),                                            &
  fef_co2_rpm, fef_co_rpm, fef_ch4_rpm,                                       &
  fef_nox_rpm, fef_so2_rpm, fef_oc_rpm, fef_bc_rpm,                           &
  !Array INTENT(OUT)
  fire_em_CO2_RPM(:), fire_em_CO_RPM(:), fire_em_CH4_RPM(:),                  &
  fire_em_NOx_RPM(:), fire_em_SO2_RPM(:), fire_em_OC_RPM(:),                  &
  fire_em_BC_RPM(:))

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE inferno_io
END MODULE inferno_io_mod
