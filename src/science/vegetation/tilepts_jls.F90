! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! SUBROUTINE tilepts
!
! Purpose:
! Counts the number of points containing each surface type and creates
! a surft_index array specifying the location of these points on the land
! grid.

! Subroutine Interface:

SUBROUTINE tilepts(land_pts,frac,surft_pts,surft_index)

USE jules_surface_mod, ONLY : all_tiles

USE ancil_info, ONLY : l_lice_point

USE jules_surface_types_mod, ONLY : ntype, ice, elev_ice

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

INTEGER, INTENT(IN) :: land_pts  !num land points to process

REAL, INTENT(IN) :: frac(land_pts,ntype)  !fractions of surface types

INTEGER, INTENT(OUT) ::  surft_pts(ntype)                                     &
                                          ! Number of land points which
                                          ! include the nth surface type
,                        surft_index(land_pts,ntype)
                                          ! Indices of land points which
                                          ! include the nth surface type

LOGICAL :: use_tile  ! Indicates if we will model the tile for the current
                     ! land point
LOGICAL :: an_ice_tile  ! Identify if we are on one of several ice tiles

INTEGER :: n,l,c   !local counters: type, land pts, count

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='TILEPTS'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------
! Create the surft_index array of land points with each surface type
!-----------------------------------------------------------------------
surft_pts (:) = 0
surft_index(:,:) = 0

DO n=1,ntype
  c=0
!CDIR NODEP
  DO l=1,land_pts
    IF( all_tiles == 0 ) THEN
! If all_tiles is off, we model tiles with frac > 0 only
      use_tile = ( frac(l,n) > 0.0 )
    ELSE
! If all_tiles is on, we model:
!   * Only the ice tile on land ice points (l_lice_point .TRUE.)
!   * All tiles except the ice tile on non-land ice points
!  (this is more complicated where land ice can have multiple tiles)

! First, identify whether this n is a land ice tile
      an_ice_tile=.FALSE.
      IF ( (n == ice) .OR.                                 &
           (maxval(elev_ice) > 0 .AND.                     &
            n >= minval(elev_ice,(elev_ice.gt.0)) .AND.    &
            n<= maxval(elev_ice)) ) an_ice_tile=.TRUE.

! Default to TRUE
      use_tile = .TRUE.

! Correct state for whether this tile is on an land ice gridbox
      IF (.NOT. l_lice_point(l)) THEN
        IF (an_ice_tile) use_tile = .FALSE.
      ELSE
        IF (.NOT.an_ice_tile) use_tile = .FALSE.
      END IF

    END IF

    IF( use_tile ) THEN
      c = c + 1
      surft_index(c,n) = l
    END IF
  END DO
  surft_pts(n) = c
END DO

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE tilepts
