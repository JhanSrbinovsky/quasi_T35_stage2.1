#if !defined(UM_JULES)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE precision_mod

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   This module defines the appropriate kinds for 32 and 64 bit integers and
!   reals
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

  INTEGER, PARAMETER ::                                                       &
    real32 = selected_real_kind(6, 37),                                       &
        ! Kind for 32-bit reals
    real64 = selected_real_kind(15, 307),                                     &
        ! Kind for 64-bit reals
    int32 = selected_int_kind(9),                                             &
        ! Kind for 32-bit integers
    int64 = selected_int_kind(15)
        ! Kind for 64-bit integers

END MODULE precision_mod
#endif
