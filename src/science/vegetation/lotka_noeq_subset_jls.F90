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
SUBROUTINE lotka_noeq_subset (land_pts,trif_pts,trif_index                    &
,                 c_veg,frac_agric,frac_agr_prev                              &
,                 GAMMA,ht,pc_s                                               &
,                 frac,dfrac,frac_na,dfrac_na,nsub,crop_lotka)

USE jules_surface_types_mod, ONLY :                                           &
!      imported scalars with intent(in)
  npft,nnpft,ntype,soil,ice,urban,lake

USE trif, ONLY: crop, g_area
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
,k,l,m,n,o,t,i,j                                                              &
                            ! WORK Loop counters.
,dom(land_pts,nsub)        ! WORK Dominance hierachy.

REAL                                                                          &
 c_veg(land_pts,npft)                                                         &
                            ! IN Carbon content of vegetation
                            !    (kg C/m2).
,frac_agric(land_pts)                                                         &
,frac_agr_prev(land_pts)                                                      &
                            ! IN Fraction of agriculture.
,GAMMA                                                                        &
                            ! IN Inverse timestep (/360days).
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
,com(land_pts,nsub,nsub)                                                      &
                            ! WORK Coefficients representing
!                                 !      the influence of one type
!                                 !      (second argument) on another
!                                 !      (first argument).
,nosoil(land_pts)                                                             &
                            ! WORK Fractional area not available
!                                 !      to vegetation.
,space(land_pts,nsub)                                                         &
                            ! WORK Space available for invasion.
,spacemax(land_pts)                                                           &
                            ! WORK Total space in the grid cell
,fracn
                            ! WORK Fractions used in the spreading
!                                 !      calculation.

REAL                          :: temp_ht(nsub)
                            ! WORK Store height for assigning dominance.
REAL                          :: tallest
                            ! Work Height of the tallest PFT from temp_ht

INTEGER                                                                       &
 nsub                                                                         &
                            ! IN number of PFTs to use
,crop_lotka
                            ! IN switch to use crop PFTs (=1) or
                            !    non-crop PFTs (not = 1)
REAL                                                                          &
pft_idx(nsub)
                            ! WORK index of the PFTs to use
                            !      (either crop or non-crop PFTS)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='LOTKA_NOEQ_SUBSET'
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

!----------------------------------------------------------------------
! Define the subset of PFTs to use
!----------------------------------------------------------------------
m=1
DO n=1,npft
  IF ( ABS(crop(n)-crop_lotka) < EPSILON(1.0) ) THEN
    pft_idx(m)=n
    m=m+1
  END IF
END DO

!----------------------------------------------------------------------
! Take working copies of dfrac and frac (removed dfrac_na, it gets
! overwritten below anyway)
!----------------------------------------------------------------------
DO n = 1,nsub
  o=pft_idx(n)
  DO l = 1,land_pts
    frac_na(l,o) = frac(l,o)
  END DO
END DO

!----------------------------------------------------------------------
! Define competition coefficients and the dominance hierachy
!----------------------------------------------------------------------
DO n = 1,nsub
  DO m = 1,nsub
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

  DO n = 1,nsub
    o=pft_idx(n)
    temp_ht(n) = ht(l,o)
  END DO

  DO i = 1,nsub
    tallest = MAXVAL(temp_ht)
    DO j = 1,nsub
      IF (ABS(temp_ht(j)-tallest) < EPSILON(1.0)) THEN
        dom(l,i) = j
        temp_ht(j) = -999.
        EXIT
      END IF
    END DO
  END DO

  DO j = 1,nsub
     n = dom(l,j)
     com(l,n,n) = 1.0
     IF (j < nsub) THEN
       DO k = j+1,nsub
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
DO k = 1,nsub
  DO t = 1,trif_pts
    l = trif_index(t)
    n = dom(l,k)
    IF ( crop_lotka >= 1 ) THEN
      space(l,n) = frac_agric(l)-frac_min*(nsub-k)
    ELSE
      space(l,n) = 1.0 - nosoil(l) - frac_agric(l) - frac_min*(nsub-k)
    END IF
  END DO
END DO

DO n = 1,nsub
  DO m = 1,nsub
    o=pft_idx(m)
    DO t=1,trif_pts
      l=trif_index(t)
      space(l,n) = space(l,n) - (com(l,n,m)*frac(l,o))
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
  DO n = 1,nsub
    o=pft_idx(n)
    dfrac(l,o) = 0.0
  END DO
  IF ( crop_lotka >= 1 ) THEN
    spacemax(l) = frac_agric(l) - frac_min*(nsub-1)
  ELSE
    spacemax(l) = 1 - nosoil(l) - frac_agric(l) - frac_min*(nsub-1)
  END IF
END DO

DO t = 1,trif_pts
  l = trif_index(t)
  DO j = 1,nsub
    n = dom(l,j)
    o=pft_idx(n)
    fracn = frac(l,o)
    fracn = MAX(fracn,frac_seed)
    dfrac(l,o) = (pc_s(l,o)*space(l,n)/c_veg(l,o) - g_area(o)) / (GAMMA/fracn)
    frac(l,o) = frac(l,o) + dfrac(l,o)
    IF (frac(l,o) <  frac_min) THEN
        dfrac(l,o) = dfrac(l,o) + (frac_min-frac(l,o))
        frac(l,o) = frac_min
    ELSE IF (frac(l,o) >  (spacemax(l))) THEN
     IF ((spacemax(l)) <  frac_min ) THEN
        frac(l,o) = frac_min
        dfrac(l,o) = frac(l,o) - fracn
     ELSE
        dfrac(l,o) = dfrac(l,o) + spacemax(l) - frac(l,o)
        frac(l,o) = spacemax(l)
     END IF
    END IF
    spacemax(l) = spacemax(l)-frac(l,o)+frac_min

  END DO
END DO

!----------------------------------------------------------------------
! Diagnose the new bare soil fraction
!----------------------------------------------------------------------
DO t = 1,trif_pts
  l = trif_index(t)
  frac(l,soil) = 1.0 - nosoil(l)
  DO n = 1,nnpft
    frac(l,soil) = frac(l,soil) - frac(l,n)
  END DO
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
DO k = 1,nsub
  DO t = 1,trif_pts
    l = trif_index(t)
    n = dom(l,k)
    IF ( crop_lotka >= 1 ) THEN
      space(l,n) = frac_agr_prev(l)-frac_min*(nsub-k)
    ELSE
      space(l,n) = 1.0-nosoil(l)-frac_agr_prev(l)-frac_min*(nsub-k)
    END IF
  END DO
END DO

DO n = 1,nsub
  DO m = 1,nsub
    o=pft_idx(m)
    DO t = 1,trif_pts
       l = trif_index(t)
       space(l,n) = space(l,n) - (com(l,n,m)*frac_na(l,o))
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
  DO n = 1,nsub
    o=pft_idx(n)
    dfrac_na(l,o) = 0.0
  END DO
  IF ( crop_lotka >= 1 ) THEN
    spacemax(l) = frac_agr_prev(l) - frac_min*(nsub-1)
  ELSE
    spacemax(l) = 1 - nosoil(l) -frac_agr_prev(l)- frac_min*(nsub-1)
  END IF
END DO

DO t = 1,trif_pts
  l = trif_index(t)
  DO j = 1,nsub
    n = dom(l,j)
    o=pft_idx(n)
    fracn = frac_na(l,o)
    fracn = MAX(fracn,frac_seed)

    dfrac_na(l,o) = (pc_s(l,o)*space(l,n)/c_veg(l,o) - g_area(o))             &
                     / (GAMMA/fracn)
    frac_na(l,o) = frac_na(l,o) + dfrac_na(l,o)

    IF (frac_na(l,o) < frac_min) THEN
        dfrac_na(l,o) = dfrac_na(l,o) + (frac_min-frac_na(l,o))
        frac_na(l,o) = frac_min
    ELSE IF (frac_na(l,o) > (spacemax(l))) THEN
     IF ( spacemax(l) <  frac_min ) THEN
        frac_na(l,o) = frac_min
        dfrac_na(l,o) = frac_na(l,o) - fracn
     ELSE
        dfrac_na(l,o) = dfrac_na(l,o) + spacemax(l)-frac_na(l,o)
        frac_na(l,o) = spacemax(l)
     END IF
    END IF

    spacemax(l) = spacemax(l) - frac_na(l,o) + frac_min

  END DO
END DO

!----------------------------------------------------------------------
! Diagnose the new bare soil fraction
!----------------------------------------------------------------------
DO t = 1,trif_pts
  l = trif_index(t)
  frac_na(l,soil) = 1.0 - nosoil(l)
  DO n = 1,nnpft
    frac_na(l,soil) = frac_na(l,soil) - frac_na(l,n)
  END DO
  IF (frac_na(l,soil) < frac_min) THEN !AVOID NEGATIVE SOIL FRACTIONS
    frac_na(l,soil) = frac_min
  END IF
END DO

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE lotka_noeq_subset
