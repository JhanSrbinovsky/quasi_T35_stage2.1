! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!! Subroutine WOODPROD -----------------------------------------------
!!!
!!! Purpose : Update wood products pool, driven by forest clearance
!!!           from land use change
!!!
!!! Allocation rules between the pools are:
!!!
!!!             QUICK            MED      SLOW
!!!             1 yr + burn     10 yr     100 yr
!!!     NL      60%             30%       10%
!!!     BL      60%             40%       0%
!!!     SH      80%             20%       0%
!!!
!!!  These values are from McGuire et al 2001 "Carbon balance of the
!!!  terrestrial biosphere in the 20th century:.." Global Biogeochemical
!!!  Cycles, Vol 15(1): 183-206, table 3, recalculated to make use of
!!!  TRIFFID simulated ratios of above ground : below ground carbon storage.
!!!
!!! -------------------------------------------------------------------
SUBROUTINE woodprod (trif_pts,land_pts,trif_index                             &
     ,                    GAMMA                                               &
     ,                    lit_c_ag                                            &
     ,                    wood_prod_fast,wood_prod_med                        &
     ,                    wood_prod_slow                                      &
     ,                    wp_fast_in,wp_med_in,wp_slow_in                     &
     ,                    wp_fast_out,wp_med_out,wp_slow_out                  &
)

USE jules_surface_types_mod

USE trif, ONLY : alloc_fast, alloc_med, alloc_slow

IMPLICIT NONE

INTEGER                                                                       &
trif_pts                                                                      &
     ,land_pts                                                                &
     ,trif_index(land_pts)                                                    &
                            ! IN Indices of land points on
!                                 !    which TRIFFID may operate
     ,j,k,l,t                    ! Loop counters

REAL                                                                          &
GAMMA                                                                         &
                            ! IN Inverse timestep (/360days).
     ,wood_prod_fast(land_pts)                                                &
                            ! INOUT Fast-turnover wood product C pool.
     ,wood_prod_med(land_pts)                                                 &
                            ! INOUT Medium-turnover wood product C pool.
     ,wood_prod_slow(land_pts)                                                &
                            ! INOUT Slow-turnover wood product C pool.
     ,wp_fast_in(land_pts)                                                    &
                            ! OUT Fast-turnover wood product C pool input.
     ,wp_med_in(land_pts)                                                     &
                            ! OUT Fast-turnover wood product C pool input.
     ,wp_slow_in(land_pts)                                                    &
                            ! OUT Fast-turnover wood product C pool input.
     ,wp_fast_out(land_pts)                                                   &
                            ! OUT Fast-turnover wood product C pool output.
     ,wp_med_out(land_pts)                                                    &
                            ! OUT Fast-turnover wood product C pool output.
     ,wp_slow_out(land_pts)                                                   &
                            ! OUT Fast-turnover wood product C pool output.
     ,lit_c_ag(land_pts,nnpft)                                                &
!                                   IN litter carbon
!                                      (kg C/m2/360days).
     ,denom                                                                   &
                            ! WORK Denominator of update
!                                 !      equation.
     ,numer                                                                   &
!                            ! WORK Numerator of the update
!                                 !      equation.
     ,wpresp(3)
!                      !Wood products psuedo respiration rates


! LOSS RATE CONSTANTS OF THE POOLS
     DATA wpresp /1.0,1.e-1,1.e-2/


  DO t=1,trif_pts
    l=trif_index(t)

! Calculate fluxes

! Wood products respiration (loss)

     wp_fast_out(l)=wood_prod_fast(l)*wpresp(1)
     wp_med_out(l) =wood_prod_med(l) *wpresp(2)
     wp_slow_out(l)=wood_prod_slow(l)*wpresp(3)

! Gain

    DO k=1,npft
        wp_fast_in(l)=wp_fast_in(l)                                           &
                       + (lit_c_ag(l,k)*alloc_fast(k))
        wp_med_in(l)=wp_med_in(l)                                             &
                       + (lit_c_ag(l,k)*alloc_med(k))
        wp_slow_in(l)=wp_slow_in(l)                                           &
                       + (lit_c_ag(l,k)*alloc_slow(k))
    END DO


! Update Wood product pools

    wood_prod_fast(l)=wood_prod_fast(l)-(wp_fast_out(l))/GAMMA
    wood_prod_med(l)=wood_prod_med(l)-(wp_med_out(l))/GAMMA
    wood_prod_slow(l)=wood_prod_slow(l)-(wp_slow_out(l))/GAMMA

    wood_prod_fast(l)=MAX(wood_prod_fast(l),0.0)
    wood_prod_med(l)=MAX(wood_prod_med(l),0.0)
    wood_prod_slow(l)=MAX(wood_prod_slow(l),0.0)

    wood_prod_fast(l)=wood_prod_fast(l)                                       &
                       + (wp_fast_in(l))/GAMMA
    wood_prod_med(l)=wood_prod_med(l)                                         &
                       + (wp_med_in(l))/GAMMA
    wood_prod_slow(l)=wood_prod_slow(l)                                       &
                       + (wp_slow_in(l))/GAMMA

  END DO

RETURN
END SUBROUTINE woodprod
