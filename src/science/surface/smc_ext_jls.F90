! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE smc_ext_mod

  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SMC_EXT_MOD'

CONTAINS

SUBROUTINE smc_ext (npnts,nshyd,surft_pts,surft_index,ft                      &
,                   f_root,sthu,v_crit,v_sat,v_wilt                           &
,                   wt_ext,fsmc)
!---------------------------------------------------------------------

! Description:
!     Calculates the soil moisture availability factor and
!     the fraction of the transpiration which is extracted from each
!     soil layer.

! Documentation : UM Documentation Paper 25
!---------------------------------------------------------------------
USE pftparm, ONLY : fsmc_mod

  IMPLICIT NONE

! Subroutine arguments
  INTEGER, INTENT(IN) ::                                                        &
   npnts                                                                        &
                        ! Number of gridpoints.
  ,nshyd                                                                        &
                        ! Number of soil moisture layers.
  ,surft_pts                                                                    &
                        ! Number of points containing the
  !                     !    given surface type.
  ,surft_index(npnts)                                                           &    
                        ! Indices on the land grid of the
  !                     !    points containing the given
  !                     !    surface type.                          
  ,ft                   ! Plant functional type.

  REAL, INTENT(IN) ::                                                           &
   f_root(npnts,nshyd)                                                          &
                        ! Fraction of roots in each soil
  !                     !    layer.
  ,sthu(npnts,nshyd)                                                            &
                        ! Unfrozen soil moisture content of
  !                     !    each layer as a fraction of
  !                     !    saturation.
  ,v_crit(npnts,nshyd)                                                          &
                        ! Volumetric soil moisture
  !                     !    concentration at field capacity
  !                     !    (m3 H2O/m3 soil).
  ,v_sat(npnts,nshyd)                                                           &
                        ! Volumetric soil moisture
  !                     !    concentration at saturation
  !                     !    (m3 H2O/m3 soil).
  ,v_wilt(npnts,nshyd)  ! Volumetric soil moisture
  !                     !    concentration below which
  !                     !    stomata close (m3 H2O/m3 soil).

  REAL, INTENT(OUT) ::                                                          &
   wt_ext(npnts,nshyd)  ! Cummulative fraction of transpiration
  !                     !    extracted from each soil layer
  !                     !    (kg/m2/s).

  REAL, INTENT(OUT) ::                                                          &
   fsmc(npnts)          ! Soil moisture availability
  !                     !    factor.

  ! work
  INTEGER ::                                                                    &
   i,j,n                ! Loop counters

  REAL ::                                                                       &
   fsmc_l(npnts,nshyd)                                                          &
  !                     ! Soil moisture availability
  !                     !    factor for each soil layer.
  ,v_layer(npnts,nshyd) !
  !                     ! Volumetric soil moisture
  !                     !    concentration of layer
  !                     !    (m3 H2O/m3 soil).

  REAL ::                                                                       & 
   v_root_zone(npnts)                                                           &              
  !                     ! Volumetric soil moisture
  !                     !    concentration of root zone
  !                     !    (m3 H2O/m3 soil).                                    
  ,v_wilt_root_zone(npnts)                                                      &                       
  !                     ! Volumetric soil moisture
  !                     !    concentration of root zone below which
  !                     !    stomata close (m3 H2O/m3 soil).                   
  ,v_crit_root_zone(npnts)
  !                     ! Volumetric soil moisture
  !                     !    concentration at field capacity
  !                     !    (m3 H2O/m3 soil).  

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='SMC_EXT'

  !----------------------------------------------------------------------
  ! Initialisations
  !----------------------------------------------------------------------
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  fsmc(:) = 0.0
  v_layer(:,:) = sthu(:,:) * v_sat(:,:)
  wt_ext(:,:) = 0.0
  
  IF ( fsmc_mod(ft) == 1 ) THEN
    v_root_zone  = calc_weighted_mean(npnts, nshyd, surft_pts, surft_index,             &
                              sthu * v_sat, f_root)
    
    v_wilt_root_zone = calc_weighted_mean(npnts, nshyd, surft_pts, surft_index,         &
                              v_wilt, f_root)
    v_crit_root_zone = calc_weighted_mean(npnts, nshyd, surft_pts, surft_index,         &
                              v_crit, f_root)

    DO j=1,surft_pts
      i=surft_index(j)
 
      fsmc(i) = fsmc_layer(ft, v_root_zone(i), v_wilt_root_zone(i), v_crit_root_zone(i), 1.0)
    END DO

    DO n=1,nshyd
      DO j=1,surft_pts
        i=surft_index(j)
        wt_ext(i,n) = f_root(i,n) * v_layer(i,n)
      END DO
    END DO

    wt_ext = calc_norm_weights(npnts, nshyd, surft_pts, surft_index,             &
                               wt_ext)

  ELSE
    !----------------------------------------------------------------------
    ! Calculate the soil moisture availability factor for each layer and
    ! weight with the root fraction to calculate the total availability
    ! factor.
    !----------------------------------------------------------------------
    DO n=1,nshyd
!CDIR NODEP
      DO j=1,surft_pts
        i=surft_index(j)
        fsmc_l(i,n) = fsmc_layer(ft, sthu(i,n), v_wilt(i,n), v_crit(i,n), v_sat(i,n))
      END DO
    END DO

    fsmc = calc_weighted_mean(npnts, nshyd, surft_pts, surft_index,            &
                              fsmc_l, f_root)

    !----------------------------------------------------------------------
    ! Calculate the fraction of the transpiration which is extracted from
    ! each soil layer.
    !----------------------------------------------------------------------
    DO n=1,nshyd
      DO j=1,surft_pts
        i=surft_index(j)
        IF (fsmc(i) > 0.0)                                                     &
          wt_ext(i,n) = f_root(i,n) * fsmc_l(i,n) / fsmc(i)
      END DO
    END DO

  END IF

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END SUBROUTINE smc_ext

FUNCTION calc_norm_weights(npnts, nshyd, surft_pts, surft_index,             &
                           weights) RESULT (norm_weights)
!-----------------------------------------------------------------------------
! Description:
!   Normalises an array of weights over the soil levels dimension.
!   The weights should be >= 0.
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

  IMPLICIT NONE

! function arguments
  INTEGER, INTENT(IN) :: npnts              ! Number of gridpoints.
  INTEGER, INTENT(IN) :: nshyd              ! Number of soil moisture layers.
  INTEGER, INTENT(IN) :: surft_pts          ! Number of points containing the
!                                           !    given surface type.
  INTEGER, INTENT(IN) :: surft_index(npnts) ! Indices on the land grid of the
!                                           !    points containing the given
!                                           !    surface type.

  REAL, INTENT(IN) :: weights(npnts,nshyd)  ! Unnormalised weights. 
!                                           ! Should all be >= 0. 

! work
  REAL :: norm_factor                       ! Factor used in the normalisation
  INTEGER :: i,j,n                          ! Loop counters

! returns
  REAL :: norm_weights(npnts,nshyd)         ! Normalised weights 

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='CALC_NORM_WEIGHTS'
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  norm_weights(:, :) = 0.0 
  norm_factor = 0.0 

  DO n=1,nshyd
    DO j=1,surft_pts
      i=surft_index(j)

      norm_factor = SUM(weights(i,:))

      IF (norm_factor > 0.0) THEN
        norm_weights(i,n) = weights(i,n) / norm_factor
      END IF
    END DO
  END DO

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END FUNCTION calc_norm_weights

FUNCTION calc_weighted_mean(npnts, nshyd, surft_pts, surft_index,             &
                            var, norm_weights) RESULT (weighted_mean)
!-----------------------------------------------------------------------------
! Description:
!   Calculates the weighted mean of a variable over soil levels
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

  IMPLICIT NONE

! function arguments
  INTEGER, INTENT(IN) :: npnts              ! Number of gridpoints.
  INTEGER, INTENT(IN) :: nshyd              ! Number of soil moisture layers.
  INTEGER, INTENT(IN) :: surft_pts          ! Number of points containing the
!                                           !    given surface type.
  INTEGER, INTENT(IN) :: surft_index(npnts) ! Indices on the land grid of the
!                                           !    points containing the given
!                                           !    surface type.

  REAL, INTENT(IN) :: var(npnts,nshyd)      ! Variable to take the weighted mean of
  REAL, INTENT(IN) :: norm_weights(npnts,nshyd)  
!                                           ! Normalised weights to use in the mean

! work
  INTEGER :: i,j,n                          ! Loop counters

! returns
  REAL :: weighted_mean(npnts)     

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='CALC_WEIGHTED_MEAN'
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  weighted_mean(:)=0.0 

  DO n=1,nshyd
!CDIR NODEP
    DO j=1,surft_pts
      i=surft_index(j)

      weighted_mean(i) = weighted_mean(i) + norm_weights(i,n) * var(i,n)

    END DO
  END DO

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END FUNCTION calc_weighted_mean

FUNCTION fsmc_layer(pft, v, v_wilt, v_crit, extra_factor) RESULT (fsmc_l)
!-----------------------------------------------------------------------------
! Description:
!   Calculates the soil water availability factor of a soil layer
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
USE pftparm, ONLY : fsmc_p0

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: pft  ! Plant functional type.
  REAL, INTENT(IN) :: v, extra_factor  ! v*extra_factor is the volumetric soil moisture
!                             !    concentration in layer
!                             !    (m3 H2O/m3 soil).  
  REAL, INTENT(IN) :: v_wilt  ! Volumetric soil moisture
!                             !    concentration below which
!                             !    stomata close (m3 H2O/m3 soil).
  REAL, INTENT(IN) :: v_crit  ! Volumetric soil moisture
!                             !    concentration at field capacity
!                             !    (m3 H2O/m3 soil).

! returns
  REAL :: fsmc_l  ! soil water availability factor for this soil layer

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='FSMC_LAYER'
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
  
  IF ( abs(v_crit - v_wilt) > 0.0 ) THEN
    fsmc_l = (v * extra_factor - v_wilt) / (v_crit - v_wilt) / ( 1.0 - fsmc_p0(pft))
  ELSE
    fsmc_l = 0.0
  END IF
  
  fsmc_l = MAX(fsmc_l,0.0)
  fsmc_l = MIN(fsmc_l,1.0)
 
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END FUNCTION fsmc_layer

END MODULE smc_ext_mod
