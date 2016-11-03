#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


MODULE logging_mod

  USE mpi

  USE ISO_FORTRAN_ENV, ONLY : OUTPUT_UNIT, ERROR_UNIT

  IMPLICIT NONE

! Log levels - these can be combined using bitwise operators to indicate
! any combination of log levels
  INTEGER, PARAMETER ::                                                       &
    LOG_LEVEL_INFO  = 1,                                                      &
    LOG_LEVEL_DEBUG = 2,                                                      &
    LOG_LEVEL_WARN  = 4,                                                      &
    LOG_LEVEL_ERROR = 8,                                                      &
    LOG_LEVEL_FATAL = 16

! Determines what log levels are printed to log_unit - this is a bitwise
! combination of values from above.
! The default (31) is to print everything
  INTEGER, PARAMETER :: LOG_PRINT_LEVEL = 31

! Determines what log levels cause a program to stop - this is a bitwise
! combination of values from above.
! Fatal errors will always cause the program to stop, by definition.
! The default (0) is to stop only for fatal errors
! Setting this to 15 (i.e. stop for everything, even info) or 14 (i.e. stop
! for everything except info) are useful options for debugging
  INTEGER, PARAMETER :: LOG_STOP_LEVEL = 0

! The maximum line length for log entries. Any message larger than this is
! truncated
  INTEGER, PARAMETER :: LOG_MAX_LINE_LEN = 500

! The calculated number of tasks id and task prefix
! These are cached the first time they are calculated to avoid calling the
! MPI routines every time a log message is written
  INTEGER :: ntasks = 0
  CHARACTER(len=20) :: task_prefix = ""


! Visibilities
  PRIVATE
  PUBLIC :: log_info, log_debug, log_warn, log_error, log_fatal

CONTAINS

  SUBROUTINE write_to_log(log_level, proc_name, message)

#if defined(COMPILER_INTEL)
    USE IFCORE, ONLY : TRACEBACKQQ
#endif

    USE string_utils_mod, ONLY : to_string

!-----------------------------------------------------------------------------
! Description:
!   Logs the given message at the given level and performs any necessary
!   action
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
    INTEGER, INTENT(IN) :: log_level
      ! The level at which to log the given message
    CHARACTER(len=*), INTENT(IN) :: proc_name
      ! The name of the originating routine/function
    CHARACTER(len=*), INTENT(IN) :: message
      ! The message to log

! Work variables
    CHARACTER(len=LOG_MAX_LINE_LEN) :: full_message
    INTEGER :: task_id, error

!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Construct the full log message
!-----------------------------------------------------------------------------
! Combine the name of the originating procedure with the message
    full_message = TRIM(proc_name) // ': ' // message

! Prepend a fragment to the message depending on the log level
    SELECT CASE ( log_level )
      CASE ( LOG_LEVEL_INFO )
        full_message = "[INFO] " // full_message

      CASE ( LOG_LEVEL_DEBUG )
        full_message = "[DEBUG] " // full_message

      CASE ( LOG_LEVEL_WARN )
        full_message = "[WARNING] " // full_message

      CASE ( LOG_LEVEL_ERROR )
        full_message = "[ERROR] " // full_message

      CASE ( LOG_LEVEL_FATAL )
        full_message = "[FATAL ERROR] " // full_message

      CASE DEFAULT
! This should never happen since the only public access to write_to_log is
! through the log_* routines defined below
        CALL log_fatal("write_to_log", "Unknown log level")
    END SELECT

! Calculate the number of tasks and task prefix if we have not already
    IF( ntasks <= 0 ) THEN
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, error)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, task_id, error)
      IF( ntasks > 1 ) THEN
        task_prefix = "{MPI Task " // TRIM(to_string(task_id)) // "}"
      END IF
    END IF

! Prepend the task prefix
    full_message = ADJUSTL(TRIM(task_prefix) // " " // full_message)

!-----------------------------------------------------------------------------
! Use a bitwise and to check if we want to print log messages for the
! given log level
!-----------------------------------------------------------------------------
    IF ( IAND(LOG_PRINT_LEVEL, log_level) > 0 ) THEN
      IF( log_level <= LOG_LEVEL_DEBUG ) THEN
! Print info and debug to stdout
        WRITE(OUTPUT_UNIT, "(A)") TRIM(full_message)
      ELSE
! Print errors to stdout
        WRITE(ERROR_UNIT, "(A)") TRIM(full_message)
      END IF
    END IF

!-----------------------------------------------------------------------------
! Check if we need to stop the program
!-----------------------------------------------------------------------------
    IF ( IAND(LOG_STOP_LEVEL, log_level) > 0 .OR.                             &
         log_level == LOG_LEVEL_FATAL ) THEN
! If stopping abnormally, attempt to print a stack trace
#if defined(COMPILER_GFORTRAN)
! Currently, the default Met Office version of gfortran is too old to support
! BACKTRACE, as is the default gfortran version in some common OS package
! systems
! We'll add this back in when gfortran 4.8 is fully available
!      CALL BACKTRACE()
#elif defined(COMPILER_INTEL)
      CALL TRACEBACKQQ(USER_EXIT_CODE = -1)
#endif

! Abort MPI with a non-zero exit code
      CALL MPI_ABORT(MPI_COMM_WORLD, 1, error)
    END IF

  END SUBROUTINE write_to_log


  SUBROUTINE log_info(proc_name, message)

!-----------------------------------------------------------------------------
! Description:
!   Utility function to write a message to the log at level info
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
    CHARACTER(len=*), INTENT(IN) :: proc_name
    CHARACTER(len=*), INTENT(IN) :: message

    CALL write_to_log(LOG_LEVEL_INFO, proc_name, message)

  END SUBROUTINE log_info


  SUBROUTINE log_debug(proc_name, message)

!-----------------------------------------------------------------------------
! Description:
!   Utility function to write a message to the log at level debug
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
    CHARACTER(len=*), INTENT(IN) :: proc_name
    CHARACTER(len=*), INTENT(IN) :: message

    CALL write_to_log(LOG_LEVEL_DEBUG, proc_name, message)

  END SUBROUTINE log_debug


  SUBROUTINE log_warn(proc_name, message)

!-----------------------------------------------------------------------------
! Description:
!   Utility function to write a message to the log at level warn
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
    CHARACTER(len=*), INTENT(IN) :: proc_name
    CHARACTER(len=*), INTENT(IN) :: message

    CALL write_to_log(LOG_LEVEL_WARN, proc_name, message)

  END SUBROUTINE log_warn


  SUBROUTINE log_error(proc_name, message)

!-----------------------------------------------------------------------------
! Description:
!   Utility function to write a message to the log at level error
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
    CHARACTER(len=*), INTENT(IN) :: proc_name
    CHARACTER(len=*), INTENT(IN) :: message

    CALL write_to_log(LOG_LEVEL_ERROR, proc_name, message)

  END SUBROUTINE log_error


  SUBROUTINE log_fatal(proc_name, message)

!-----------------------------------------------------------------------------
! Description:
!   Utility function to write a message to the log at level fatal
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
    CHARACTER(len=*), INTENT(IN) :: proc_name
    CHARACTER(len=*), INTENT(IN) :: message

    CALL write_to_log(LOG_LEVEL_FATAL, proc_name, message)

  END SUBROUTINE log_fatal

END MODULE logging_mod
#endif
