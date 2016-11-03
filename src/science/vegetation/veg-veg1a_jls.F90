! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Version 1A of vegetation section: models leaf phenology

! Subroutine Interface:
SUBROUTINE veg1(                                                              &
               land_pts,nsurft,can_model                                      &
,              a_step,phenol_period,l_phenol                                  &
,              atimestep, satcon, z0m_soil                                    &
,              g_leaf_ac,frac,lai,ht                                          &
,              catch_s,catch_t,infil_t,z0_t, z0h_t                            &
,              g_leaf_day,g_leaf_phen,g_leaf_phen_ac,lai_phen                 &
               )

USE jules_surface_types_mod, ONLY : npft, ntype, nnpft

use jules_surface_mod, ONLY : l_aggregate

use jules_vegetation_mod, ONLY : l_gleaf_fix

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE jules_print_mgr, ONLY :                                                   &
    jules_message,                                                            &
    jules_print
IMPLICIT NONE

! Description:
!   Updates Leaf Area Index for Plant Functional Types (PFTs) and uses
!   this to derive new vegetation parameters for PFTs along with gridbox
!   mean values where appropriate.

! Method:
!   Calls PHENOL which models phenology and updates Leaf Area Index
!   (LAI), then passes new LAI into SPARM along with canopy height
!   and fractional cover of Plant Functional Types.  SPARM uses this to
!   derive the vegetation parameters for each PFT, and also derives
!   gridbox means where this is required.

INTEGER                                                                       &
 land_pts                                                                     &
                       ! IN Number of land points to be processed.
,nsurft                                                                       &
                       ! IN Number of land-surface tiles.
,can_model                                                                    &
                              ! IN Swith for thermal vegetation
,a_step                                                                       &
                       ! IN Atmospheric timestep number.
,phenol_period         ! IN Phenology period (days).

INTEGER                                                                       &
 j,l,n                      ! WORK loop counters.

LOGICAL                                                                       &
 l_phenol                     ! IN .T. for interactive leaf
!                                   !    phenology.

REAL                                                                          &
 atimestep                                                                    &
                              ! IN Atmospheric timestep (s).
,satcon(land_pts)                                                             &
                              ! IN Saturated hydraulic
!                                   !    conductivity (kg/m2/s).
,z0m_soil(land_pts)                                                           &
                              ! IN z0m of bare soil
,g_leaf_ac(land_pts,npft)                                                     &
                              ! INOUT Accumulated leaf turnover
!                                   !       rate.
,frac(land_pts,ntype)                                                         &
                              ! INOUT Fractions of surface types.
,lai(land_pts,npft)                                                           &
                              ! INOUT LAI of plant functional
!                                   !       types.
,ht(land_pts,npft)                                                            &
                              ! INOUT Height of plant functional
!                                   !       types (m).
,catch_s(land_pts,nsurft)                                                     &
                              ! OUT Snow capacity for tiles
!                                   !     (kg/m2).
,catch_t(land_pts,nsurft)                                                     &
                              ! OUT Canopy capacity for tiles
!                                   !     (kg/m2).
,infil_t(land_pts,nsurft)                                                     &
                              ! OUT Maximum surface infiltration
!                                   !     rate for tiles (kg/m2/s).
,g_leaf_day(land_pts,npft)                                                    &
                              ! OUT Mean leaf turnover rate for
!                                   !      input to PHENOL (/360days).
,g_leaf_phen(land_pts,npft)                                                   &
                              ! OUT Mean leaf turnover rate over
!                                   !      phenology period (/360days).
,g_leaf_phen_ac(land_pts,npft)                                                &
                              ! OUT Accumulated mean leaf turnover rate
                                    !over phenology period (/360days)
,lai_phen(land_pts,npft)                                                      &
                              ! OUT LAI of PFTs after phenology.
!                                   !     Required as separate variable
!                                   !     for top-level argument list
!                                   !     matching with VEG_IC2A.
,z0_t(land_pts,nsurft)       &! OUT Roughness length for tiles (m)
,z0h_t(land_pts,nsurft)       ! OUT Thermal roughness length for
                              !     tiles (m)

INTEGER                                                                       &
 nstep_phen                                                                   &
                              ! WORK Number of atmospheric
!                                   !      timesteps between calls to
!                                   !      PHENOL.
,surft_pts(ntype)                                                             &
                              ! WORK Number of land points which
!                                   !      include the nth surface type.
,surft_index(land_pts,ntype)   ! WORK Indices of land points which
!                                   !      include the nth surface type.

REAL                                                                          &
 dtime_phen
                              ! WORK The phenology timestep (yr).

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='VEG1'

!-----------------------------------------------------------------------
! Initialisations
!-----------------------------------------------------------------------
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

DO n=1,nsurft
  DO l=1,land_pts
    catch_s(l,n)=0.0
    catch_t(l,n)=0.0
    infil_t(l,n)=0.0
    z0_t(l,n)=0.0
  END DO
END DO

DO n=1,npft
  DO l=1,land_pts
    g_leaf_phen(l,n)=0.0
    g_leaf_day(l,n)=0.0
  END DO
END DO

!-----------------------------------------------------------------------
! Calculate the number of atmospheric timesteps between calls to PHENOL
! and TRIFFID.
!-----------------------------------------------------------------------
nstep_phen=INT(86400.0*phenol_period/atimestep)


!-----------------------------------------------------------------------
! Create the surft_index array of land points with each surface type
!-----------------------------------------------------------------------
! DEPENDS ON: tilepts
CALL tilepts(land_pts,frac,surft_pts,surft_index)

IF (l_phenol .AND. MOD(a_step,nstep_phen) == 0) THEN

!-----------------------------------------------------------------------
! Calculate the phenology timestep in years.
!-----------------------------------------------------------------------
  dtime_phen=FLOAT(phenol_period)/360.0


  DO n=1,nnpft

!-----------------------------------------------------------------------
! Calculate the mean turnover rate and update the leaf phenological
! state.
!-----------------------------------------------------------------------
    DO j=1,surft_pts(n)
      l=surft_index(j,n)
      g_leaf_day(l,n)=g_leaf_ac(l,n)/dtime_phen
    END DO

    WRITE(jules_message,*) 'Calling phenology'
    CALL jules_print('veg-veg1a_jls',jules_message)

! DEPENDS ON: phenol
    CALL phenol (land_pts,surft_pts(n),surft_index(:,n),n,                    &
                 g_leaf_day(:,n),ht(:,n),dtime_phen,                          &
                 g_leaf_phen(:,n),lai(:,n))

    WRITE(jules_message,*) 'Phenology completed normally'
    CALL jules_print('veg-veg1a_jls',jules_message)

    DO l=1,land_pts
      lai_phen(l,n)=lai(l,n)
    END DO

    IF( l_gleaf_fix ) THEN
      !gleaf fix copied from veg-veg2a
      !----------------------------------------------------------------
      ! Increment the leaf turnover rate for driving TRIFFID and reset
      ! the accumulation over atmospheric model timesteps to zero.
      !----------------------------------------------------------------
      DO j=1,surft_pts(n)
        l=surft_index(j,n)
        g_leaf_phen_ac(l,n)=g_leaf_phen_ac(l,n)                               &
                           + (g_leaf_phen(l,n)*dtime_phen)
      END DO
    END IF

!-----------------------------------------------------------------------
! Reset the accumulation over atmospheric model timesteps to zero.
!-----------------------------------------------------------------------
    DO l=1,land_pts
      g_leaf_ac(l,n)=0.0
    END DO
  END DO
END IF

!-----------------------------------------------------------------------
! Calculate gridbox mean vegetation parameters from fractions of
! surface functional types
!-----------------------------------------------------------------------
! DEPENDS ON: sparm
CALL sparm (land_pts, nsurft, can_model, l_aggregate,                         &
            surft_pts, surft_index,                                           &
            frac, ht, lai, satcon, z0m_soil, catch_s, catch_t,                &
            infil_t, z0_t, z0h_t)

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE veg1
