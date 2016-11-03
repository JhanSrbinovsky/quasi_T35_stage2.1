! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine DPM_RPM ------------------------------------------------
!
! Purpose : Calculates the DPM_RPM ratio of litter input for use
!            in the RothC soil carbon sub-model
!
! -------------------------------------------------------------------

SUBROUTINE dpm_rpm(land_pts,trif_pts,trif_index,                              &
                     lit_c,dpm_ratio)

USE jules_surface_types_mod

USE trif, ONLY : dpm_rpm_ratio

USE ereport_mod, ONLY : ereport

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

INTEGER                                                                       &
 land_pts                                                                     &
                          ! IN Total number of land points.
,trif_pts                                                                     &
                            ! IN Number of points on which
!                                 !    TRIFFID may operate
,trif_index(land_pts)                                                         &
                          ! IN Indices of land points on
!                                 !    which TRIFFID may operate
,n,l,t                        ! WORK Loop counters

REAL                                                                          &
 lit_c(land_pts,npft)                                                         &
                          ! IN  PFT carbon litter
!                                 !    (kg C/m2/360days).
,dpm_ratio(land_pts)      ! OUT Gridbox mean DPM ratio of litter
                          !     input DPM ratio of litter input

REAL                                                                          &
 total_litter                                                                 &
                          ! WORK  gridbox total lit_c
,dpm(land_pts)            ! WORK  litter input to DPM

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='DPM_RPM'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

! calculate DPM litter input and hence ratio
DO t=1,trif_pts
  l=trif_index(t)
  dpm(l) = 0.0
  total_litter = 0.0

  DO n = 1,npft
    dpm(l) = dpm(l) + lit_c(l,n)*dpm_rpm_ratio(n)/                            &
                                 (1+dpm_rpm_ratio(n))
    total_litter = total_litter + lit_c(l,n)
  END DO
  dpm_ratio(l) = dpm(l) / total_litter

  dpm_ratio(L) = MAX(0.0, dpm_ratio(l))
  dpm_ratio(L) = MIN(1.0, dpm_ratio(l))

END DO

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE dpm_rpm
