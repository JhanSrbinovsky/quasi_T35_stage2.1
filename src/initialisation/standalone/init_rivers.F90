#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
MODULE init_rivers_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE init_rivers(nml_dir)
!-----------------------------------------------------------------------------
! Description:
!   Reads in the river routing namelist items and checks them for consistency
!   Initialises river routing parameters
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
  USE mpi, ONLY : MPI_COMM_WORLD

  USE io_constants, ONLY : NAMELIST_unit

  USE string_utils_mod, ONLY : to_string

  USE jules_rivers_mod, ONLY : jules_rivers, check_jules_rivers,              &
                               l_rivers, rivers_type, rivers_timestep,        &
                               rivers_umtrip, rivers_umrfm

  USE logging_mod, ONLY : log_info, log_debug, log_warn, log_error, log_fatal

  IMPLICIT NONE


! Arguments
  CHARACTER(len=*), INTENT(IN) :: nml_dir  ! The directory containing the
                                           ! namelists

! Work variables
  INTEGER :: error  ! Error indicator
  INTEGER :: ntasks ! Parallel mode indicator

!-----------------------------------------------------------------------------
! Read river routing namelist
!----------------------------------------------------------------------------
  CALL log_info("init_rivers", "Reading JULES_RIVERS namelist...")

! Open the river routing parameters namelist file
  OPEN(NAMELIST_UNIT, FILE=(TRIM(nml_dir) // '/' // 'jules_rivers.nml'),      &
                 STATUS='old', POSITION='rewind', ACTION='read', IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_rivers",                                             &
                   "Error opening namelist file jules_rivers.nml " //         &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

  READ(NAMELIST_UNIT, nml=jules_rivers, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_rivers",                                             &
                   "Error reading namelist JULES_RIVERS " //                  &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

! Close the namelist file
  CLOSE(NAMELIST_UNIT, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_rivers",                                             &
                   "Error closing namelist file jules_rivers.nml " //         &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

  CALL check_jules_rivers()

  IF (l_rivers) THEN
    CALL log_info("init_rivers",                                              &
                  TRIM(rivers_type) // " river routing is selected")
    CALL log_info("init_rivers",                                              &
                  "River routing timestep = " //                              &
                   TRIM(to_string(rivers_timestep)))

! Check if river routing selected in parallel mode - not currently available
! option in JULES standalone
    SELECT CASE (rivers_type)
    CASE(rivers_umrfm, rivers_umtrip)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ntasks,error)
      IF (ntasks > 1) THEN
        CALL log_fatal("init_rivers",                                         &
                       "River routing not currently available when " //       &
                       "running in parallel mode with this rivers_type.")
      END IF
    END SELECT
  ELSE
    CALL log_info("init_rivers", "No river routing selected")
  END IF

 END SUBROUTINE init_rivers

END MODULE init_rivers_mod
#endif
