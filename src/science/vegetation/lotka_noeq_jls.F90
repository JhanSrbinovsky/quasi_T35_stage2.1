! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine LOTKA --------------------------------------------------
!
! Purpose : Updates fractional coverage of each functional type.
!           Based on the Lotka-Volterra equations of interspecies
!           competition.
!
! -------------------------------------------------------------------
SUBROUTINE lotka_noeq (land_pts,trif_pts,trif_index                           &
,                 c_veg,frac_agric,frac_agr_prev                              &
,                 GAMMA,lai,ht,pc_s                                           &
,                 frac,dfrac,frac_na,dfrac_na)

USE jules_surface_types_mod, ONLY :                                           &
!      imported scalars with intent(in)
  npft,nnpft,ntype,soil,ice,urban,lake

USE descent
USE pftparm
USE trif
USE jules_vegetation_mod, ONLY: frac_min,frac_seed,pow

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
,k,l,m,n,t,i,j                                                                &
                            ! WORK Loop counters.
,dom(land_pts,nnpft)        ! WORK Dominance hierachy.

REAL                                                                          &
 c_veg(land_pts,npft)                                                         &
                            ! IN Carbon content of vegetation
                            !    (kg C/m2).
,frac_agric(land_pts)                                                         &
,frac_agr_prev(land_pts)                                                      &
                            ! IN Fraction of agriculture.
,GAMMA                                                                        &
                            ! IN Inverse timestep (/360days).
,lai(land_pts,npft)                                                           &
                            ! IN Leaf area index.
,ht(land_pts,npft)                                                            &
                            ! IN canopy height (m)
,pc_s(land_pts,npft)                                                          &
                            ! IN Net carbon flux available for
                            !    spreading (kg C/m2/360days).
,frac(land_pts,ntype)                                                         &
                            ! INOUT Fractional cover of each
!                                 !       Functional Type.
,frac_na(land_pts,ntype)                                                      &
                                  ! INOUT Fractional cover of each
!                                 !       Functional Type.
,dfrac(land_pts,npft)                                                         &
                            ! OUT Increment to the areal fraction
!                                 !     during the timestep (/timestep).
,dfrac_na(land_pts,npft)                                                      &
                                  ! OUT Increment to the areal fraction
!                                 !     during the timestep (/timestep).
!                                 !     prior to the land use change
,com(land_pts,nnpft,nnpft)                                                    &
                            ! WORK Coefficients representing
!                                 !      the influence of one type
!                                 !      (second argument) on another
!                                 !      (first argument).
,nosoil(land_pts)                                                             &
                            ! WORK Fractional area not available
!                                 !      to vegetation.
,space(land_pts,nnpft)                                                        &
                            ! WORK Space available for invasion.
,spacemax(land_pts)                                                           &
                            ! WORK Total space in the grid cell

,fracn
                            ! WORK Fraction used in the spreading
!                                 !      calculation.

REAL                          :: temp_ht(npft)
                            ! WORK Store height for assigning dominance.
REAL                          :: tallest
                            ! Work Height of the tallest PFT from temp_ht


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='LOTKA_NOEQ'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

!----------------------------------------------------------------------
! Take working copies of dfrac and frac (removed dfrac_na, it gets
! overwritten below anyway
!----------------------------------------------------------------------
DO n = 1,ntype
  DO l = 1,land_pts
    frac_na(l,n) = frac(l,n)
  END DO
END DO

!----------------------------------------------------------------------
! Define competition coefficients and the dominance hierachy
!----------------------------------------------------------------------
DO n = 1,nnpft
  DO m = 1,nnpft
    DO t = 1,trif_pts
      l = trif_index(t)
      com(l,n,m) = 1.0
    END DO
  END DO
END DO

!NEW: COM(I,J) BASED ON HEIGHT OF PFT I & J
DO t = 1,trif_pts
  l=trif_index(t)

  !DOM IS THE DOMINANCE HIERARCHY BASED ON HEIGHTS
  !WARNING: If two or more PFTs
  !have exactly the same height, the first PFT in the list
  !will get the higher dominance.
  temp_ht(1:nnpft) = ht(l,1:nnpft)
  DO i = 1,nnpft
    tallest = MAXVAL(temp_ht)
    DO j = 1,nnpft
      IF (ABS(temp_ht(j)-tallest) < EPSILON(1.0)) THEN
        dom(l,i) = j
        temp_ht(j) = -999.
        EXIT
      END IF
    END DO
  END DO

  DO j = 1,nnpft
     n = dom(l,j)
     com(l,n,n) = 1.0
     IF (j < nnpft) THEN
       DO k = j+1,nnpft
          m = dom(l,k)
          com(l,m,n) = 1.0
          com(l,n,m) = 0.0
       END DO
     END IF
  END DO
END DO

!----------------------------------------------------------------------
! Calculate the space available for the expansion of each FT
!----------------------------------------------------------------------
DO t = 1,trif_pts
  l = trif_index(t)
  nosoil(l) = frac(l,urban) + frac(l,lake) + frac(l,ice)
END DO

!----------------------------------------------------------------------
! Exclude non-crop types from agricultural regions
!----------------------------------------------------------------------
DO k = 1,nnpft
  DO t = 1,trif_pts
    l = trif_index(t)
    n = dom(l,k)
    space(l,n) = 1.0 - nosoil(l) - frac_agric(l)*(1-crop(n))                  &
                     - frac_min*(nnpft-k)
  END DO
END DO

DO n = 1,nnpft
  DO m = 1,nnpft
    DO t=1,trif_pts
      l=trif_index(t)
      space(l,n) = space(l,n) - (com(l,n,m)*frac(l,m))
    END DO
  END DO
END DO

!----------------------------------------------------------------------
! We are removing the implicit timestepping approach to equilibrium.
! Divide the update equation by FRAC to eliminate the (unstable)
! bare soil solution.
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Update the areal fractions (Compete code now here:)
!----------------------------------------------------------------------
spacemax(:) = 0.0
DO t = 1,trif_pts
  l = trif_index(t)
  DO n = 1,nnpft
    dfrac(l,n) = 0.0
  END DO
  spacemax(l) = 1.0 - nosoil(l) - frac_min*(nnpft-1)
END DO

DO t = 1,trif_pts
  l = trif_index(t)
  DO j = 1,nnpft
    n = dom(l,j)
    fracn = frac(l,n)
    fracn = MAX(fracn,frac_seed)
    dfrac(l,n) = (pc_s(l,n)*space(l,n)/c_veg(l,n) - g_area(n)) / (GAMMA/fracn)
    frac(l,n) = frac(l,n) + dfrac(l,n)

    IF (frac(l,n) > (spacemax(l)-(1.0-crop(n))*frac_agric(l))) THEN
      dfrac(l,n) = dfrac(l,n) +                                               &
                   (spacemax(l)-(1.0-crop(n))*frac_agric(l))-frac(l,n)
      frac(l,n) = spacemax(l)-(1.0-crop(n))*frac_agric(l)
    END IF

    IF (frac(l,n) < frac_min) THEN
      dfrac(l,n) = dfrac(l,n) + (frac_min-frac(l,n))
      frac(l,n) = frac_min
    END IF

    spacemax(l) = spacemax(l)-frac(l,n)+frac_min

  END DO
END DO

!----------------------------------------------------------------------
! Diagnose the new bare soil fraction
!----------------------------------------------------------------------
DO t = 1,trif_pts
  l = trif_index(t)
  frac(l,soil) = spacemax(l)
  IF (frac(l,soil) < frac_min) THEN !AVOID NEGATIVE SOIL FRACTIONS
    frac(l,soil) = frac_min
  END IF
END DO


!----------------------------------------------------------------------
! For the dynamic land use we require a double call to COMPETE
! in order to be able to clearly diagnose what change is caused by the
! land use change.
!----------------------------------------------------------------------

!----------------------------------------------------------------------
! Exclude non-crop types from agricultural regions
!----------------------------------------------------------------------
DO k = 1,nnpft
  DO t = 1,trif_pts
    l = trif_index(t)
    n = dom(l,k)
    space(l,n) = 1.0 - nosoil(l) - frac_agr_prev(l)*(1-crop(n))               &
                     - frac_min*(nnpft-k)
  END DO
END DO

DO n = 1,nnpft
  DO m = 1,nnpft
    DO t = 1,trif_pts
       l = trif_index(t)
       space(l,n) = space(l,n) - (com(l,n,m)*frac_na(l,m))
     END DO
   END DO
END DO

!SECOND CALL TO COMPETE TO ACCOUNT FOR CHANGES FROM LU
!----------------------------------------------------------------------
! We are removing the implicit timestepping approach to equilibrium.
! Divide the update equation by FRAC to eliminate the (unstable)
! bare soil solution.
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Update the areal fractions (Compete code now here:)
!----------------------------------------------------------------------
spacemax(:) = 0.0
DO t = 1,trif_pts
  l = trif_index(t)
  DO n = 1,nnpft
    dfrac_na(l,n) = 0.0
  END DO
  spacemax(l) = 1 - nosoil(l) - frac_min*(nnpft-1)
END DO

DO t = 1,trif_pts
  l = trif_index(t)
  DO j = 1,nnpft
    n = dom(l,j)
    fracn = frac_na(l,n)
    fracn = MAX(fracn,frac_seed)

    dfrac_na(l,n) = (pc_s(l,n)*space(l,n)/c_veg(l,n) - g_area(n))             &
                    / (GAMMA/fracn)
    frac_na(l,n) = frac_na(l,n) + dfrac_na(l,n)

    IF (frac_na(l,n) < frac_min) THEN
        dfrac_na(l,n) = dfrac_na(l,n) + (frac_min-frac_na(l,n))
        frac_na(l,n) = frac_min
    ELSE IF (frac_na(l,n) > (spacemax(l)-(1.0-crop(n))*frac_agr_prev(l))) THEN
     IF ((spacemax(l)-(1.0-crop(n))*frac_agr_prev(l)) < frac_min ) THEN
        dfrac_na(l,n) = dfrac_na(l,n) + (frac_min - frac_na(l,n))
        frac_na(l,n) = frac_min
     ELSE
        dfrac_na(l,n) = dfrac_na(l,n) + (spacemax(l)-                         &
                    (1.0-crop(n))*frac_agr_prev(l))-frac_na(l,n)
        frac_na(l,n) = spacemax(l)-(1.0-crop(n))*frac_agr_prev(l)
     END IF
    END IF

    spacemax(l) = spacemax(l) - frac_na(l,n) + frac_min

  END DO
END DO

!----------------------------------------------------------------------
! Diagnose the new bare soil fraction
!----------------------------------------------------------------------
DO t = 1,trif_pts
  l = trif_index(t)
  frac_na(l,soil) = spacemax(l)
  IF (frac_na(l,soil) < frac_min) THEN !AVOID NEGATIVE SOIL FRACTIONS
    frac_na(l,soil) = frac_min
  END IF
END DO


IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE lotka_noeq
