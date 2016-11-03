#if !defined(UM_JULES) && defined(MPI_DUMMY)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

MODULE mpi

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   This module contains just the variable definitions from a "real" mpi
!   module that JULES needs to compile against
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

  INTEGER, PARAMETER :: MPI_COMM_WORLD = 1
  INTEGER, PARAMETER :: MPI_REAL = 1
  INTEGER, PARAMETER :: MPI_INTEGER = 1
  INTEGER, PARAMETER :: MPI_LOGICAL = 1
  INTEGER, PARAMETER :: MPI_INFO_NULL = 1
  INTEGER, PARAMETER :: MPI_LAND = 1

! Setting this to 1 means use the default INTEGER size
  INTEGER, PARAMETER :: MPI_ADDRESS_KIND = 1

END MODULE mpi
#endif
