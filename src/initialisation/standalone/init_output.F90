#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE init_output(nml_dir)

  USE io_constants, ONLY : MAX_SDF_NAME_LEN, MAX_VAR_FILE, NAMELIST_UNIT

  USE datetime_mod, ONLY : DATETIME_STR_LEN, datetime, datetime_from_string

  USE string_utils_mod, ONLY : to_string

  USE model_time_mod, ONLY : timestep_len, main_run_start, main_run_end

  USE model_interface_mod, ONLY : IDENTIFIER_LEN

  USE output_mod, ONLY : output_dir, run_id, register_output_profile

  USE theta_field_sizes, ONLY : t_i_length, t_j_length

  USE jules_snow_mod, ONLY : nsmax

  USE jules_vegetation_mod, ONLY : l_fapar_diag, l_fao_ref_evapotranspiration

  USE logging_mod, ONLY : log_info, log_debug, log_warn, log_error, log_fatal

  USE errormessagelength_mod, ONLY: errormessagelength

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads in information about what output has been requested and sets up
!   the output profiles
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Arguments
  CHARACTER(len=*), INTENT(IN) :: nml_dir  ! The directory containing the
                                           ! namelists

! Work variables
  LOGICAL :: dir_exists  ! Used to check existence of output directory

  TYPE(datetime) :: output_start_dt, output_end_dt
                         ! Datetime objects created from strings given in
                         ! output profiles

  INTEGER :: nvars_in  ! The number of variables that were input for the
                       ! profile currently being processed, before filtering
  LOGICAL :: remove_var  ! Indicates if the current variable should be removed
                         ! from the list as it is not allowed

  INTEGER :: i,j  ! Index variables

  INTEGER :: error  ! Error indicator
  CHARACTER(LEN=errormessagelength) :: iomessage

!-----------------------------------------------------------------------------
! Definition of the jules_output namelist
!-----------------------------------------------------------------------------
  INTEGER :: nprofiles  ! The number of output profiles

  NAMELIST /jules_output/ output_dir, run_id, nprofiles

!-----------------------------------------------------------------------------
! Definition of the JULES_OUTPUT_PROFILE namelist
!-----------------------------------------------------------------------------
  CHARACTER(len=MAX_SDF_NAME_LEN) :: profile_name  ! The name of the profile
  LOGICAL :: output_initial ! T - this profile should output initial data
                            !     for each section it is outputting
                            ! F - this profile should not output initial data
  LOGICAL :: output_spinup  ! T - generate output during spinup
                            ! F - don't generate output during spinup
  LOGICAL :: output_main_run  ! T - generate output for the part of
                              !     the main run specified by
                              !     output_start and output_end
                              ! F - don't generate output for any
                              !     of the main run
  CHARACTER(len=DATETIME_STR_LEN) :: output_start, output_end
                                        ! The start and end times for output
                                        ! during the main run as strings
                                        ! If not given, output_start is assumed
                                        ! to be the start of the main run and
                                        ! output_end is assumed to be the end
                                        ! of the main run
  INTEGER :: output_period
                                        ! The period of output - assumed to be
                                        ! every timestep
                                        ! This can be a special period
  INTEGER :: sample_period
                                        ! The period for sampling data - assumed
                                        ! to be every timestep
                                        ! Must be a multiple of the model
                                        ! timestep
  INTEGER :: file_period                ! The period for new files - unless this
                                        ! is PERIOD_YEAR or PERIOD_MONTH, it is
                                        ! ignored

  INTEGER :: nvars  ! The number of variables in the output profile
  CHARACTER(len=IDENTIFIER_LEN) :: var(MAX_VAR_FILE)
                        ! The model identifiers for variables in the output
                        ! profile
  CHARACTER(len=IDENTIFIER_LEN) :: var_name(MAX_VAR_FILE)
                        ! The name to use for variables in output files
  CHARACTER(len=1) :: output_type(MAX_VAR_FILE)
                        ! The type of output to use - this should be one of
                        ! 'S', 'A', 'M', 'N' or 'X'

  NAMELIST /jules_output_profile/ profile_name, output_initial, output_spinup,&
                                  output_main_run, output_start, output_end,  &
                                  output_period, sample_period, file_period,  &
                                  nvars, var, var_name, output_type


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Initialise
!-----------------------------------------------------------------------------
  nprofiles = 0

!-----------------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------------
  OPEN(NAMELIST_UNIT, FILE=(TRIM(nml_dir) // '/' // 'output.nml'),            &
                 STATUS='old', POSITION='rewind', ACTION='read', IOSTAT=error,&
                 IOMSG=iomessage)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_output",                                             &
                   "Error opening namelist file output.nml " //               &
                   "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //       &
                   TRIM(iomessage) // ")")

!-----------------------------------------------------------------------------
! Read the JULES_OUTPUT namelist
!-----------------------------------------------------------------------------
  CALL log_info("init_output", "Reading JULES_OUTPUT namelist...")
  READ(NAMELIST_UNIT, nml=jules_output, IOSTAT=error, IOMSG=iomessage)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_output",                                             &
                   "Error reading namelist JULES_OUTPUT " //                  &
                   "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //       &
                   TRIM(iomessage) // ")")

!-----------------------------------------------------------------------------
! Process the information from the JULES_OUTPUT namelist
!-----------------------------------------------------------------------------
! First check that a run_id was given
  IF ( LEN_TRIM(run_id) == 0 )                                                &
    CALL log_fatal("init_output", "No run_id given")

! Check that the output directory exists
! The Intel compiler requires a different form of the statement for directories,
! which we swap in with an ifdef
#if defined(COMPILER_INTEL)
  INQUIRE(DIRECTORY=output_dir, EXIST=dir_exists)
#else
  INQUIRE(FILE=output_dir, EXIST=dir_exists)
#endif
  IF ( .NOT. dir_exists )                                                     &
    CALL log_fatal("init_output",                                             &
                   "Output directory does not exist - " // TRIM(output_dir))

! Warn if no output is requested
! We do this after the directory check as we will always produce dumps
  IF ( nprofiles < 1 ) THEN
    CALL log_warn("init_output",                                              &
                  "No output profiles given - output will not be " //         &
                  "generated for this run")
    RETURN
  END IF

!-----------------------------------------------------------------------------
! Read and process information about each output profile in turn
!-----------------------------------------------------------------------------
  DO i = 1,nprofiles
! Set namelist values to their defaults before reading the next profile
    profile_name    = ""
    output_initial  = .FALSE.
    output_spinup   = .FALSE.
    output_main_run = .FALSE.
    output_start    = ""
    output_end      = ""
    output_period   = timestep_len
    sample_period   = timestep_len
    file_period     = 0
    nvars           = 0
    var(:)          = ''
    var_name(:)     = ''
    output_type(:)  = 'S'

! Read the namelist
    CALL log_info("init_output", "Reading JULES_OUTPUT_PROFILE namelist...")
    READ(NAMELIST_UNIT, nml=jules_output_profile, IOSTAT=error,               &
         IOMSG=iomessage)
    IF ( error /= 0 )                                                         &
      CALL log_fatal("init_output",                                           &
                     "Error reading namelist JULES_OUTPUT_PROFILE " //        &
                     "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //     &
                     TRIM(iomessage) // ")")

! Filter any variables not allowed from the list, warning about their removal
    nvars_in = nvars
    nvars = 0
    DO j = 1,nvars_in
! Assume that we are going to keep the variable until we decide otherwise
      remove_var = .FALSE.

! At the moment, the only variables not allowed are the multilayer snow
! variables, and those derived from them, when nsmax < 1
      IF ( nsmax < 1 ) THEN
        SELECT CASE ( var(j) )
          CASE ( 'snow_ice_gb', 'snow_ice_tile',                              &
                 'snow_liq_gb', 'snow_liq_tile',                              &
                 'rgrainl', 'snow_ds', 'snow_ice', 'snow_liq', 'tsnow' )
            remove_var = .TRUE.

          CASE DEFAULT
            remove_var = .FALSE.

        END SELECT
      END IF

      IF ( remove_var ) THEN
! If the variable needs to be removed, issue a warning that we are doing it
        CALL log_error("init_output",                                         &
                       "Variable " // TRIM(var(j)) // " is not available " // &
                       "for output with the current configuration - " //      &
                       "removing from output")
      ELSE
! Otherwise add the variable to the collapsed list - we can do this since
! nvars <= j so we are not overwriting unprocessed data
        nvars = nvars + 1
        var(nvars)         = var(j)
        var_name(nvars)    = var_name(j)
        output_type(nvars) = output_type(j)
      END IF
    END DO

! If this profile is not providing any variables, we can skip it
    IF ( nvars < 1 ) THEN
      CALL log_error("init_output",                                           &
                     "Profile " // TRIM(profile_name) // " is not " //        &
                     "outputting any variables - ignoring")
      CYCLE
    END IF

! Set up the datetime objects for output start and end
! If not given, they default to start and end of main run respectively
    IF ( LEN_TRIM(output_start) == 0 ) THEN
      output_start_dt = main_run_start
    ELSE
      output_start_dt = datetime_from_string(output_start)
    END IF
    IF ( LEN_TRIM(output_end) == 0 ) THEN
      output_end_dt = main_run_end
    ELSE
      output_end_dt = datetime_from_string(output_end)
    END IF

    IF ( ANY(var(1:nvars) == 'fapar') ) THEN
      l_fapar_diag = .TRUE.
    END IF

    IF ( ANY(var(1:nvars) == 'fao_et0') ) THEN
      l_fao_ref_evapotranspiration = .TRUE.
    END IF

! Register the output profile - this performs more error checking
    CALL register_output_profile(profile_name, output_initial,                &
                                 output_spinup, output_main_run,              &
                                 output_start_dt, output_end_dt,              &
                                 output_period, sample_period, file_period,   &
                                 var(1:nvars), var_name(1:nvars),             &
                                 output_type(1:nvars))

  END DO

  CLOSE(NAMELIST_UNIT, IOSTAT=error, IOMSG=iomessage)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_output",                                             &
                   "Error closing namelist file output.nml " //               &
                   "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //       &
                   TRIM(iomessage) // ")")

  RETURN

END SUBROUTINE init_output
#endif
