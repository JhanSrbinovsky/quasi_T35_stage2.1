! *****************************COPYRIGHT*******************************

! (c) [University of Edinburgh] [2009]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the JULES collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC237]

! *****************************COPYRIGHT*******************************
!  SUBROUTINE LAYERSNOW-----------------------------------------------

! Description:
!     Divide snow pack into layers if it exceeds a minimum depth.

! Subroutine Interface:
MODULE layersnow_mod
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='LAYERSNOW_MOD'
CONTAINS
SUBROUTINE layersnow ( land_pts,surft_pts,surft_index,snowdepth,              &
                       nsnow,ds )

USE jules_snow_mod, ONLY :                                                    &
 nsmax                                                                        &
                    !  Maximum possible number of snow layers
,dzsnow             !  Prescribed snow layer depths (m)

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Scalar arguments with intent(in)
INTEGER, INTENT(IN) ::                                                        &
 land_pts                                                                     &
                      !  Number of land points
,surft_pts             !  Number of tile points

! Array arguments with intent(in)
INTEGER, INTENT(IN) ::                                                        &
 surft_index(land_pts)    !  Index of tile points

REAL, INTENT(IN) ::                                                           &
 snowdepth(land_pts)     !  Snow depth (m)

! Array arguments with intent(out)
INTEGER, INTENT(OUT) ::                                                       &
 nsnow(land_pts)         !  Number of snow layers

REAL, INTENT(OUT) ::                                                          &
 ds(land_pts,nsmax)     !  Snow layer thicknesses (m)

! Local scalars
INTEGER ::                                                                    &
 k                                                                            &
                     !  Tile point index
,l                                                                            &
                     !  Land point index
,n                   !  Snow layer index

REAL ::                                                                       &
 remains             !  Remaining depth of snow for other layers

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='LAYERSNOW'

!-------------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialise number of layers. This value is not changed for the zero-layer
! model.
nsnow(:) = 0

IF ( nsmax > 0 ) THEN

  ! Initialise layer depths (this value will persist at locations where tile
  ! frac=0).
  ds(:,:) = 0.0

  DO k=1,surft_pts
    l = surft_index(k)

!   Only divide into layers if depth is >= a threshold.
    IF ( snowdepth(l) >= dzsnow(1) ) THEN
      remains = snowdepth(l)

      DO n=1,nsmax
        ds(l,n) = dzsnow(n)
        remains = remains - dzsnow(n)
        IF ( remains <= dzsnow(n) .OR. n == nsmax ) THEN
          ds(l,n) = ds(l,n) + remains
          EXIT
        END IF
      END DO

!     Set number of layers.
      nsnow(l) = n
    END IF    !  >dzSnow(1)

  END DO

END IF  !  nsmax
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE layersnow
END MODULE layersnow_mod
