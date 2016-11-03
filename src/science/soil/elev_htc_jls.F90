! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE ELEV_HTC------------------------------------------------

! Description:
!     Updates tiled bedrock temperatures for land ice. No external 
!     subroutines are called. Modified from ice_htc @ julesvn4.1 
!     using only the explicit scheme for a single layer

! Subroutine Interface:
MODULE elev_htc_mod

  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ELEV_HTC_MOD'

  CONTAINS

  SUBROUTINE elev_htc (                                             &
   npnts,lice_pts,lice_index,nsurft,                                &
   dz,snow_soil_htf,timestep,                                       &
   tsurf_elev_surft                                                       &
  )
  
  USE jules_snow_mod, ONLY :                                        &
  !      imported scalars with intent(in)
   snow_hcap
  
  USE parkind1, ONLY: jprb, jpim
  USE yomhook, ONLY: lhook, dr_hook
  IMPLICIT NONE
  
  ! Subroutine arguments
  !   Scalar arguments with intent(IN) :
  INTEGER                                                           &
   lice_pts                                                         &
                        ! IN Number of land ice points.
  ,npnts                                                            &
                        ! IN Number of land points
  ,nsurft               ! IN Number of tiles
  
  REAL                                                              &
   timestep             ! IN Model timestep (s).
  
  
  !   Array arguments with intent(IN) :
  INTEGER                                                           &
   lice_index(npnts)    ! IN Array of ice points.
  
  REAL                                                              &
   dz                                                               &
                        ! IN Thicknesses of the bedrock layer (m).
  ,snow_soil_htf(npnts,nsurft) ! IN Net downward surface heat flux (W/m2).
  
  !   Array arguments with intent(INOUT) :
  REAL                                                              &
   tsurf_elev_surft(npnts,nsurft)   ! INOUT Sub-surface temperatures (K).
  
  ! Local scalars:
  INTEGER                                                           &
   i,j,n                ! WORK Loop counters.
  
  
  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='ELEV_HTC'
  
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  ! Original explicit scheme
  !
  !--------------------------------------------------------------------
  ! Update the sub-surface temperatures
  !--------------------------------------------------------------------
  !CDIR NOVECTOR
  !   CDIR$ IVDEP here would force vectorization but has been seen to change 
  !         results in similar code (eg ice_htc)
    DO j=1,lice_pts
      i=lice_index(j)
      DO n=1,nsurft
  
        tsurf_elev_surft(i,n)=tsurf_elev_surft(i,n)       &
         +1.0/(snow_hcap*dz)*                             &
         snow_soil_htf(i,n)*timestep
  
      END DO
    END DO
  
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  
  RETURN
  
  END SUBROUTINE elev_htc
END MODULE elev_htc_mod
