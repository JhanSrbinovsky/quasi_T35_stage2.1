! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module with UM setting of
! params used in methane flux from wetlands.

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.

MODULE c_ch4

IMPLICIT NONE

! Required for large-scale hydrology L_TOP code.
! Used in calculation of methane flux from wetlands.
!
REAL,PARAMETER :: T0_CH4 = 273.15             ! T0 value
REAL,PARAMETER :: Q10_CH4_CS = 3.7            ! Q10 value for soil
!                                             !   carbon substrate
REAL,PARAMETER :: Q10_CH4_NPP = 1.5           ! Q10 value for NPP
!                                             !   substrate
REAL,PARAMETER :: Q10_CH4_RESPS = 1.5         ! Q10 value for soil
!                                             !   repsiration substrate

REAL,PARAMETER :: CONST_CH4_NPP  = 9.99E-3    ! Scale factor for NPP substrate
REAL,PARAMETER :: CONST_CH4_RESPS  = 4.36E-3  ! Scale factor for
!                                             !  soil respiration substrate

REAL:: CONST_CH4_CS            ! Scale factor for soil carbon substrate
!                              !    (set in ch4_wetl_jls.F90)

END MODULE c_ch4
