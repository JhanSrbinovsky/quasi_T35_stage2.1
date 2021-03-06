#if defined(UM_JULES)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: River Routing

MODULE INIT_A2T_4A_mod
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='INIT_A2T_4A_MOD'

CONTAINS
SUBROUTINE INIT_A2T_4A(A_REALHD, XPA, XUA, XVA, YPA, YUA, YVA)

!  Routine: INIT_A2T_4A ---------------------------------------------
!
!  Purpose: Initialises the Lat and Long values of ATMOS for
!           regridding to TRIP river routing grid (based on INITA2O)
!           NB this will need to be extended to pick up the values
!           for the river routing grid from the ancil header of one
!           file later so will leave in redundant code for guidance.
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.3 programming standards.
!  -------------------------------------------------------------------

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE atm_fields_bounds_mod
USE UM_ParParams
USE nlsizes_namelist_mod, ONLY: A_LEN_REALHD, aocpl_row_length, aocpl_p_rows

IMPLICIT NONE
! Description: Gridline coordinates for interpolation and area-averaging
! between atmosphere and river-routing grids (Part of TYPAOCPL.h)

REAL, INTENT(IN)  :: A_REALHD(A_LEN_REALHD)   !real header
REAL, INTENT(OUT) :: XPA(AOCPL_ROW_LENGTH+1)  !Atmosphere TP longitude coordina
REAL, INTENT(OUT) :: XUA(0:AOCPL_ROW_LENGTH)  !Atmosphere U longitude coordinat
REAL, INTENT(OUT) :: XVA(AOCPL_ROW_LENGTH+1)  !Atmosphere V longitude coordinat
REAL, INTENT(OUT) :: YPA(AOCPL_P_ROWS)        !Atmosphere TP latitude coordinat
REAL, INTENT(OUT) :: YUA(AOCPL_P_ROWS)        !Atmosphere U latitude coordinate
REAL, INTENT(OUT) :: YVA(0:AOCPL_P_ROWS)      !Atmosphere V latitude coordinate

INTEGER :: I,J

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='INIT_A2T_4A'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! This would be the "correct" ENDGame grid, assuming that array
! index 1 refers to
! the first entry (and 0 to the extensions - in ENDGame u&v index 0
! referes to the
! first entry)
!
!! Pick up the ATMOS x coords from the dump header
!      XUA(0)=A_REALHD(4)-A_REALHD(1)
!      DO I=1,AOCPL_ROW_LENGTH
!        XPA(I)=A_REALHD(4)+(I-0.5)*A_REALHD(1)
!        XVA(I)=A_REALHD(4)+(I-0.5)*A_REALHD(1)
!        XUA(I)=A_REALHD(4)+(I-1)*A_REALHD(1)
!      ENDDO
!      XPA(AOCPL_ROW_LENGTH+1)=A_REALHD(4)+(AOCPL_ROW_LENGTH+0.5)*A_REALHD(1)
!      XVA(AOCPL_ROW_LENGTH+1)=A_REALHD(4)+(AOCPL_ROW_LENGTH+0.5)*A_REALHD(1)
!
!! Pick up the ATMOS y coords from the dump header
!      DO J=1,AOCPL_P_ROWS
!! JCT  S->N now  YTA(J)=A_REALHD(3)-(J-1)*A_REALHD(2)
!        YPA(J)=A_REALHD(3)+(J-0.5)*A_REALHD(2)
!        YUA(J)=A_REALHD(3)+(J-0.5)*A_REALHD(2)
!      ENDDO
!      DO J=0,AOCPL_P_ROWS
!! JCT  S->N now  YUA(J)=A_REALHD(3)-(J-0.5)*A_REALHD(2)
!        YVA(J)=A_REALHD(3)+(J)*A_REALHD(2)
!      ENDDO

! river routing will assume that a p-point gridbox is defined such that
! the first entry is at XPA(1)=0, YPA(1)=-90 and that the gridbox is defined
!
!    XUA(0),YVA(0)  XUA(1),YVA(0)
!          +--------------+
!          |              |
!          |      x  XPA(1)=0, YPA(1)=-90
!          |              |
!          +--------------+
!    XUA(0),YVA(1)  XUA(1),YVA(1)
!
! for ENDGame that is different. In ENDGame the first p-point is defined at
! XPA(1)=dx/2, YPA(1)=-90+dx/2, XUA(0)=0

XUA(0)=A_REALHD(4)
DO I=1,AOCPL_ROW_LENGTH
  XPA(I)=A_REALHD(4)+(I-0.5)*A_REALHD(1)
  XVA(I)=A_REALHD(4)+(I-0.5)*A_REALHD(1)
  XUA(I)=A_REALHD(4)+I*A_REALHD(1)
END DO
XPA(AOCPL_ROW_LENGTH+1)=A_REALHD(4)+(AOCPL_ROW_LENGTH+0.5)*A_REALHD(1)
XVA(AOCPL_ROW_LENGTH+1)=A_REALHD(4)+(AOCPL_ROW_LENGTH+0.5)*A_REALHD(1)

DO J=1,AOCPL_P_ROWS
    YPA(J)=A_REALHD(3)+(J-0.5)*A_REALHD(2)
    YUA(J)=A_REALHD(3)+(J-0.5)*A_REALHD(2)
END DO
DO J=0,AOCPL_P_ROWS
  YVA(J)=A_REALHD(3)+(J)*A_REALHD(2)
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
!----------------------------------------------------------------------
END SUBROUTINE INIT_A2T_4A
END MODULE INIT_A2T_4A_mod
#endif
