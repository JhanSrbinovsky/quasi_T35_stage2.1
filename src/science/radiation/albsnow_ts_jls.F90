! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Routine to calculate spectral snow albedos for JULES, using a
! two-stream model.

! *********************************************************************
SUBROUTINE albsnow_ts (p_field,land_field,land_index,                         &
                       snow_pts, snow_index,                                  &
                       cosz,albudir,albudif,                                  &
                       rgrain,snowmass,soot,                                  &
                       alb_snow)

USE jules_snow_mod, ONLY : cnr_g, cnr_om, dce

USE water_constants_mod, ONLY: rho_ice

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

! Subroutine arguments
!   Scalar arguments with intent(in):
INTEGER, INTENT(IN) :: p_field
                             ! Total number of grid points.
INTEGER, INTENT(IN) :: land_field
                             ! Number of land points.
INTEGER, INTENT(IN) :: snow_pts
                             ! Number of snow-covered points

!   Array arguments with intent(in):
INTEGER, INTENT(IN) :: land_index(land_field)
                             ! Index of land points.
INTEGER, INTENT(IN) :: snow_index(land_field)
                             ! Index of land points with snow cover

REAL, INTENT(IN) :: albudir(land_field,2)
                             ! Direct albedos of the substrate
REAL, INTENT(IN) :: albudif(land_field,2)
                             ! Diffuse albedos of the substrate
REAL, INTENT(IN) :: cosz(p_field)
                             ! Zenith cosine.
REAL, INTENT(IN) :: rgrain(land_field)
                             ! Snow grain size (microns).
REAL, INTENT(IN) :: snowmass(land_field)
                             ! Mass loading of snow
REAL, INTENT(IN) :: soot(p_field)
                              ! Snow soot content (kg/kg).

!   Array arguments with intent(out):
REAL, INTENT(OUT) :: alb_snow(land_field,4)
!                              Snow albedo.
!                              (:,1) - Direct beam visible
!                              (:,2) - Diffuse visible
!                              (:,3) - Direct beam near-IR
!                              (:,4) - Diffuse near-IR
!
!
! Local scalars:
REAL :: deff
                             ! Effective dimension of snow crystals
REAL :: tau
                             ! Optical depth of snow pack
REAL :: omega
                             ! Albedo fo single scattering of snowpack
REAL :: asym
                             ! Asymmetry of snow crystals
REAL :: frwd
                             ! Forward scattering in the snowpack
REAL :: beta
                             ! Diffuse backward scattering coefficient
REAL :: beta_s
                             ! Solar backward scattering coefficient
REAL, PARAMETER :: diff = 2.0
                             ! Diffusivity factor
REAL :: alpha_1
                             ! Two-stream coefficient
REAL :: alpha_2
                             ! Two-stream coefficient
REAL :: alpha_3
                             ! Two-stream coefficient
REAL :: alpha_4
                             ! Two-stream coefficient
REAL :: lmb
                             ! Two-stream coefficient
REAL :: gamma2
                             ! Two-stream coefficient
REAL :: pp
                             ! Two-stream coefficient
REAL :: sigma_1
                             ! Two-stream coefficient
REAL :: sigma_2
                             ! Two-stream coefficient
REAL :: rdd
                             ! Diffuse reflection coefficient
REAL :: tdd
                             ! Diffuse transmission coefficient
REAL :: tss
                             ! Solar transmission coefficient
REAL :: rsd
                             ! Direct-diffuse reflection coefficient
REAL :: tsd
                             ! Direct-diffuse transmission coefficient


INTEGER                                                                       &
 band,i,j,l                  ! Loop counters.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ALBSNOW_TS'


IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
!
! 2-Band scheme, VIS and NIR.
!
DO band=1, 2
  DO j=1,snow_pts
    l = snow_index(j)
    i = land_index(l)
!
!
!   Provisionally set albedos for snow layer to albedos for
!   ground to cover the case of no snow.
    alb_snow(l,2*band-1) = albudir(l,band)
    alb_snow(l,2*band)   = albudif(l,band)
!
!
!   Recalculate where snow is present.
    IF ( (snowmass(l) > 0.0) .AND. (cosz(i) > EPSILON(1.0)) ) THEN
!
!     Radiative effective dimension, recalling that rgrain
!     is in microns.
      deff = 2.0 * rgrain(l) * 1.0e-6
!
!     Single scattering properties
      tau = 3.0 * snowmass(l) / (rho_ice * deff)
      asym = cnr_g(1,band) + cnr_g(2,band) * LOG(deff)
      omega = 1.0 / (1.0 + cnr_om(1,1,band) + LOG(deff) *                     &
                    (cnr_om(2,1,band) + LOG(deff) * cnr_om(3,1,band)) +       &
                    deff * cnr_om(2,2,band) * soot(i) / dce)
!
!     Rescale
      frwd  = asym * asym
      asym  = (asym - frwd) / (1.0 - frwd)
      tau   = (1.0 - omega * frwd) * tau
      omega = (1.0 - frwd) * omega / (1.0 - omega * frwd)
!
!     Two stream coefficients (PIFM85)
      beta   = (4.0 - frwd - 3.0 * asym) / 8.0
      beta_s = 0.5 - 0.75 * asym * cosz(i)
      alpha_1 = diff * (1.0 - omega * (1.0 - beta))
      alpha_2 = diff * omega * beta
      alpha_3 = omega * beta_s
      alpha_4 = omega * (1.0 - beta_s)
!
!     Reflection and Transmission coefficients
      lmb    = SQRT(alpha_1*alpha_1 - alpha_2*alpha_2)
      gamma2 = alpha_2 / (alpha_1 + lmb)
      pp     = EXP( -lmb * tau )
      rdd = gamma2 * (1.0 - pp*pp) / (1.0 - (gamma2*pp)**2)
      tdd = pp * (1.0 - gamma2*gamma2) / (1.0 - (gamma2*pp)**2)
      tss = EXP( -tau / cosz(i) )
!     Calculate explicit solar transmission and reflection
!     unless very close to the solar beam.
      IF ( ABS(lmb * cosz(i) - 1.0) > 100.0 * EPSILON(1.0) ) THEN
        sigma_1 = ( cosz(i) *                                                 &
          (alpha_1 * alpha_3 + alpha_2 * alpha_4) - alpha_3) /                &
          ( (lmb * cosz(i))**2 - 1.0)
        sigma_2 = ( cosz(i) *                                                 &
          (alpha_1 * alpha_4 + alpha_2 * alpha_3) + alpha_4) /                &
          ( (lmb * cosz(i))**2 - 1.0)
        rsd = sigma_1 * (1.0 - tdd) - sigma_2 * rdd
        tsd = sigma_2 * (tss - tdd) - sigma_1 * rdd * tss
      ELSE
        rsd = rdd
        tsd = tdd
      END IF
!
!     Overall reflection coefficients
      alb_snow(l,2*band-1) = rsd + tdd *                                      &
        (albudif(l,band) * tsd + albudir(l,band) * tss) /                     &
        (1.0 - albudif(l,band) * rdd)
      alb_snow(l,2*band)   = rdd +                                            &
        albudif(l,band) * tdd * tdd /                                         &
        (1.0 - albudif(l,band) * rdd)
      alb_snow(l,2*band-1) = MAX(0.0, MIN(alb_snow(l,2*band-1), 1.0))
      alb_snow(l,2*band) = MAX(0.0, MIN(alb_snow(l,2*band), 1.0))
!
    END IF
  END DO
END DO


IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE albsnow_ts
