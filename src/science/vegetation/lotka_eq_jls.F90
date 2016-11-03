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
SUBROUTINE lotka_eq (land_pts,trif_pts,trif_index                             &
,                 c_veg,frac_prev,frac_agric,ht,pc_s                          &
,                 frac_eq)

USE jules_surface_types_mod, ONLY :                                           &
!      imported scalars with intent(in)
  npft,nnpft,ntype,soil,urban,lake,ice

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
,dom(land_pts,nnpft)         ! WORK Dominance hierachy.

REAL                                                                          &
 c_veg(land_pts,npft)                                                         &
                            ! IN Carbon content of vegetation
                            !    (kg C/m2).
,frac_prev(land_pts,npft)                                                     &
                            ! IN Previous time step fractional cover
,frac_agric(land_pts)                                                         &
                            ! IN Fraction of agriculture.
,ht(land_pts,npft)                                                            &
!                             ! IN canopy height (m)
,pc_s(land_pts,npft)                                                          &
                            ! IN Net carbon flux available for
                            !    spreading (kg C/m2/360days).
,frac_eq(land_pts,ntype)                                                      &
                            ! INOUT Fractional cover of each
!                                 !       Functional Type.
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

,spacemax(land_pts)
                            ! WORK Total space in the grid cell

REAL                          :: temp_ht(nnpft)
                            ! WORK Store height for assigning dominance.
REAL                          :: tallest
                            ! Work Height of the tallest PFT from temp_ht

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='LOTKA_EQ'

!----------------------------------------------------------------------
! Define competition coefficients and the dominance hierachy
!----------------------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

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
  l = trif_index(t)

  !DOM IS THE DOMINANCE HIERARCHY BASED ON HEIGHTS
  !WARNING:If two or more PFTs
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
     com(l,n,n) = 0.0
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
  nosoil(l) = frac_prev(l,lake) + frac_prev(l,urban) + frac_prev(l,ice)
END DO

!----------------------------------------------------------------------
! Exclude non-crop types from agricultural regions
!----------------------------------------------------------------------
DO j = 1,nnpft
  DO t = 1,trif_pts
    l = trif_index(t)
    n = dom(l,j)
    space(l,n) = 1.0 - nosoil(l) - frac_agric(l)*(1-crop(n))                  &
                     - frac_min*(nnpft-j)
  END DO
END DO

DO n = 1,nnpft
  DO m = 1,nnpft
    DO t = 1,trif_pts
      l = trif_index(t)
      space(l,n) = space(l,n) - com(l,n,m)*frac_eq(l,m)
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
  spacemax(l) = 1. - nosoil(l) - frac_min*(nnpft-1)
END DO

DO t = 1,trif_pts
  l = trif_index(t)
  DO j = 1,nnpft

    n = dom(l,j)
    IF (pc_s(l,n) < 0) THEN
      frac_eq(l,n) = frac_min
    ELSE
      frac_eq(l,n) = space(l,n) - (g_area(n)*c_veg(l,n)/pc_s(l,n))
    END IF

    IF (frac_eq(l,n) < frac_min) THEN
        frac_eq(l,n) = frac_min
    ELSE IF (frac_eq(l,n) > spacemax(l)) THEN
        frac_eq(l,n) = spacemax(l)
    END IF

    spacemax(l) = spacemax(l) - frac_eq(l,n) + frac_min
  END DO
END DO

!----------------------------------------------------------------------
! Diagnose the new bare soil fraction
!----------------------------------------------------------------------
DO t = 1,trif_pts
  l = trif_index(t)
  frac_eq(l,soil) = spacemax(l)
  IF (frac_eq(l,soil) < frac_min) THEN !AVOID NEGATIVE SOIL FRACTIONS
    frac_eq(l,soil) = frac_min
  END IF
END DO

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE lotka_eq
