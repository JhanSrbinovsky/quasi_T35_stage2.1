! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! *********************************************************************
! SUBROUTINE sparm
!
! Routine to calculate the gridbox mean land surface parameters from
! the areal fractions of the surface types and the structural
! properties of the plant functional types.
!
! *********************************************************************
SUBROUTINE sparm (land_pts, nsurft, can_model, l_aggregate,                   &
                  surft_pts, surft_index,                                     &
                  frac, ht, lai, satcon, z0m_soil,                            &
                  catch_snow, catch_surft, infil_surft,                       &
                  z0_surft, z0h_surft_bare)

USE jules_surface_types_mod, ONLY :                                           &
  lake, npft, ntype, urban_canyon, urban_roof, soil

USE blend_h
USE nvegparm
USE pftparm
use jules_surface_mod, ONLY : i_aggregate_opt, l_vary_z0m_soil
USE jules_snow_mod, ONLY : cansnowtile,snowloadlai
USE c_z0h_z0m, ONLY : z0h_z0m, z0h_z0m_classic
USE urban_param, ONLY : ztm_gb
USE switches_urban, ONLY : l_urban2t, l_moruses
USE dust_param, ONLY : z0_soil
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

LOGICAL                                                                       &
 l_aggregate       ! IN Logical to set aggregate surface scheme

INTEGER                                                                       &
 land_pts                                                                     &
                       ! IN Number of land points to be processed.
,nsurft                                                                       &
                       ! IN Number of surface tiles.
,can_model                                                                    &
                              ! IN Swith for thermal vegetation
,surft_pts(ntype)                                                             &
                              ! IN Number of land points which
!                                   !    include the nth surface type.
,surft_index(land_pts,ntype)   ! IN Indices of land points which
!                                   !    include the nth surface type.

REAL, INTENT(IN) ::                                                           &
 frac(land_pts,ntype),                                                        &
                 ! IN Fractional cover of each surface type.
 ht(land_pts,npft),                                                           &
                 ! IN Vegetation height (m).
 lai(land_pts,npft),                                                          &
                 ! IN Leaf area index.
 satcon(land_pts),                                                            &
                 ! IN Saturated hydraulic conductivity (kg/m2/s).
 z0m_soil(land_pts)
                 ! IN z0m of bare soil (m)

REAL                                                                          &
 catch_surft(land_pts,nsurft)                                                 &
                              ! OUT Canopy capacity for each tile
!                                   !     (kg/m2).
,infil_surft(land_pts,nsurft)                                                 &
                              ! OUT Maximum surface infiltration
!                                   !     rate for each tile (kg/m2/s).
,z0_surft(land_pts,nsurft)                                                    &
                              ! OUT Roughness length for each
!                                   !     tile (m).
,z0h_surft_bare(land_pts,nsurft)                                              &
                              ! OUT Snow-free thermal roughness length
!                                   !     for each tile (m).
,catch_snow(land_pts,nsurft)  ! OUT Snow capacity for tile (kg/m2)

REAL                                                                          &
 catch(land_pts)                                                              &
                            ! WORK GBM canopy capacity (kg/m2).
,catch_t(land_pts,ntype)                                                      &
                            ! WORK Capacities for types.
,fz0(land_pts)                                                                &
                            ! WORK Aggregation function of Z0.
,fz0h(land_pts)                                                               &
                            ! WORK Aggregation function of Z0H.
,infil(land_pts)                                                              &
                            ! WORK GBM infiltration rate(kg/m2/s).
,infil_t(land_pts,ntype)                                                      &
                            ! WORK Infiltration rates for types.
,z0(land_pts)                                                                 &
                            ! WORK GBM roughness length (m).
,z0h(land_pts)                                                                &
                            ! WORK GBM thermal roughness length (m).
,z0_t(land_pts,ntype)       ! WORK Roughness lengths for types.

INTEGER                                                                       &
 j,l,n                      ! WORK Loop counters

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SPARM'


!----------------------------------------------------------------------
! Set parameters for vegetated surface types
!----------------------------------------------------------------------
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
DO n=1,npft
! DEPENDS ON: pft_sparm
  CALL pft_sparm (land_pts,n,surft_index(:,n),surft_pts(n)                    &
,                 ht(:,n),lai(:,n),satcon                                     &
,                 catch_t(:,n),infil_t(:,n),z0_t(:,n))
END DO

! cansnowtile is only used when l_aggregate is .FALSE. - needs to be
! consistent with logic where cansnowtile is set.
IF (can_model  ==  4 .AND. .NOT. l_aggregate) THEN
  DO n=1,npft
    IF ( cansnowtile(n) ) THEN
      DO j=1,surft_pts(n)
        l = surft_index(j,n)
        catch_snow(l,n) = snowloadlai*lai(l,n)
      END DO
    END IF
  END DO
END IF

!----------------------------------------------------------------------
! Set parameters for non-vegetated surface types
!----------------------------------------------------------------------
DO n=npft+1,ntype
  DO j=1,surft_pts(n)
    l = surft_index(j,n)
    catch_t(l,n) = catch_nvg(n-npft)
    infil_t(l,n) = infil_nvg(n-npft)*satcon(l)
  END DO
END DO

IF (l_vary_z0m_soil) THEN
  DO n=npft+1,ntype
    IF (n == soil) THEN
      ! for soil, get the z0m from the input data:
      DO l=1,land_pts
        z0_t(l,n) = z0m_soil(l)
      END DO
    ELSE
      ! for non soil types, get the z0m from the input namelist:
      DO j=1,surft_pts(n)
        l = surft_index(j,n)
        z0_t(l,n) = z0_nvg(n-npft)
      END DO
    END IF
  END DO
ELSE
  ! get all the non veg types:
  DO n=npft+1,ntype
    DO j=1,surft_pts(n)
      l = surft_index(j,n)
      z0_t(l,n) = z0_nvg(n-npft)
    END DO
  END DO
END IF

! MORUSES Set canyon & roof roughness length
IF ( l_moruses ) THEN
  n = urban_canyon
  DO j=1,surft_pts(n)
    l = surft_index(j,n)
      z0_t(l,n)          = ztm_gb(l)
      z0_t(l,urban_roof) = ztm_gb(l)
    END DO
  END IF

IF ( l_aggregate ) THEN
!----------------------------------------------------------------------
! Form means and copy to tile arrays if required for aggregate tiles
!----------------------------------------------------------------------
  DO l=1,land_pts
    catch(l) = 0.0
    fz0(l) = 0.0
    fz0h(l) = 0.0
    infil(l) = 0.0
    z0(l) = 0.0
  END DO

  DO n=1,ntype
    DO j=1,surft_pts(n)
      l = surft_index(j,n)
      fz0(l) = fz0(l) + frac(l,n) / (LOG(lb / z0_t(l,n)))**2
!     Explicit aggregation of z0h if required.
      IF (i_aggregate_opt == 1)                                               &
        fz0h(l) = fz0h(l) + frac(l,n) /                                       &
          ( LOG(lb / z0_t(l,n)) * LOG(lb / (z0h_z0m(n) * z0_t(l,n))) )
    END DO
  END DO

  DO l=1,land_pts
    z0(l) = lb * EXP(-SQRT(1. / fz0(l)))
!   Explicit aggregation of z0h if required.
    IF (i_aggregate_opt == 1)                                                 &
      z0h(l) = lb * EXP(-1.0 / (fz0h(l) * LOG(lb / z0(l))) )
  END DO

  DO n=1,ntype
    DO j=1,surft_pts(n)
      l = surft_index(j,n)
      catch(l) = catch(l) + frac(l,n) * catch_t(l,n)
      infil(l) = infil(l) + frac(l,n) * infil_t(l,n)
    END DO
  END DO

  DO l=1,land_pts
!         Canopy capacity is average over non-lake surface types
    catch_surft(l,1) = 0.
    IF ( lake > 0 ) THEN
      IF ( frac(l,lake) < 1. )                                                &
        catch_surft(l,1) = catch(l) / (1. - frac(l,lake))
    END IF
    infil_surft(l,1) = infil(l)
    z0_surft(l,1) = z0(l)
    IF (i_aggregate_opt == 1) z0h_surft_bare(l,1) = z0h(l)
  END DO

ELSE
!----------------------------------------------------------------------
! Copy surface-type arrays to tiles if separate tiles used
!----------------------------------------------------------------------
  DO n=1,ntype
    DO j=1,surft_pts(n)
      l = surft_index(j,n)
      catch_surft(l,n) = catch_t(l,n)
      infil_surft(l,n) = infil_t(l,n)
      z0_surft(l,n) = z0_t(l,n)
    END DO
  END DO

END IF
!----------------------------------------------------------------------
! Set bare soil roughness for use in 1 tile dust scheme
!----------------------------------------------------------------------
! The following line has been added to readlsta so to fix CRUNs
! and is now a duplication.
! When l_vary_z0m_soil is true, the variable is also passed down to the
! dust emission, so this isn't used.

z0_soil=z0_nvg(soil-npft)

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE sparm
