#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE init_top()

  USE io_constants, ONLY : MAX_FILE_NAME_LEN, MAX_SDF_NAME_LEN, NAMELIST_UNIT

  USE string_utils_mod, ONLY : to_string

  USE templating_mod, ONLY : tpl_has_var_name, tpl_substitute_var

  USE model_interface_mod, ONLY : IDENTIFIER_LEN, populate_var, get_var_id

  USE jules_hydrology_mod, ONLY : l_top

  USE dump_mod, ONLY: ancil_dump_read
  
  USE errormessagelength_mod, ONLY: errormessagelength

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads in the TOPMODEL properties and checks them for consistency
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
  INTEGER, PARAMETER :: MAX_TOPMODEL_VARS = 3
         ! The maximum possible number of TOPMODEL variables that can be given

  INTEGER :: nvars_required      ! The number of variables that are
                                 ! required in this configuration
  CHARACTER(len=IDENTIFIER_LEN) :: required_vars(MAX_TOPMODEL_VARS)
                                 ! The variable identifiers of the required
                                 ! variables

  INTEGER :: nvars_file       ! The number of variables that will be set
                              ! from the given file (template?)

  INTEGER :: error  ! Error indicator

  INTEGER :: i  ! Loop counter

!-----------------------------------------------------------------------------
! Definition of the jules_top namelist - this combines local variables with
! some from c_topog
!-----------------------------------------------------------------------------
  LOGICAL :: read_from_dump
  CHARACTER(len=MAX_FILE_NAME_LEN) :: file
                        ! The name of the file (or variable name template) to
                        ! use for variables that need to be filled from file

  INTEGER :: nvars      ! The number of variables in this section
  CHARACTER(len=IDENTIFIER_LEN) :: var(MAX_TOPMODEL_VARS)
                        ! The variable identifiers of the variables
  LOGICAL :: use_file(MAX_TOPMODEL_VARS)
                        !   T - the variable uses the file
                        !   F - the variable is set using a constant value
  CHARACTER(len=MAX_SDF_NAME_LEN) :: var_name(MAX_TOPMODEL_VARS)
                        ! The name of each variable in the file
  CHARACTER(len=MAX_SDF_NAME_LEN) :: tpl_name(MAX_TOPMODEL_VARS)
                        ! The name to substitute in a template for each
                        ! variable
  CHARACTER(LEN=errormessagelength) :: iomessage
                        ! I/O error message string
  REAL :: const_val(MAX_TOPMODEL_VARS)
                        ! The constant value to use for each variable if
                        ! use_file = F for that variable
  NAMELIST /jules_top/ read_from_dump, file, nvars, var, use_file, var_name, &
                       tpl_name, const_val


!-----------------------------------------------------------------------------

! If TOPMODEL is not on, we have nothing to do
  IF ( .NOT. l_top ) RETURN


!-----------------------------------------------------------------------------
! Initialise
!-----------------------------------------------------------------------------
  read_from_dump = .FALSE.
  nvars_required = 0
  nvars_file     = 0
  nvars          = 0
  use_file(:)    = .TRUE.  ! Default is for every variable to be read from file
  file           = ''      ! Empty file name
  tpl_name(:)    = ''      ! Empty template string

!-----------------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------------
  CALL log_info("init_top", "Reading JULES_TOP namelist...")

  READ(NAMELIST_UNIT, nml=jules_top, IOSTAT=error, IOMSG=iomessage)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_top",                                                &
                   "Error reading namelist JULES_TOP " //                     &
                   "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //       &
                   TRIM(iomessage) // ")")

  ancil_dump_read%top = read_from_dump

  IF ( .NOT. ancil_dump_read%top) THEN !we read from the ancil file
    !---------------------------------------------------------------------------
    ! Set up TOPMODEL properties using namelist values
    !---------------------------------------------------------------------------
    ! Set up the required variables
    ! All the TOPMODEL variables are always required
    nvars_required = MAX_TOPMODEL_VARS
    ! ti_skew may be available in later releases
    required_vars(:) = (/ 'fexp   ', 'ti_mean', 'ti_sig ' /)


    !---------------------------------------------------------------------------
    ! Check that all the required variables are there
    !---------------------------------------------------------------------------
    DO i = 1,nvars_required
      IF ( .NOT. ANY(var(1:nvars) == required_vars(i)) )                      &
        CALL log_fatal("init_top",                                            &
                       "No value given for required variable '" //            &
                       TRIM(required_vars(i)) // "'")
    END DO

    !---------------------------------------------------------------------------
    ! Check which variables we will be using and partition them into variables
    ! set to constant values and variables set from file
    !---------------------------------------------------------------------------
    DO i = 1,nvars
    !---------------------------------------------------------------------------
    ! If the variable is one of the required vars, then we will be using it
    !---------------------------------------------------------------------------
      IF ( ANY(required_vars(1:nvars_required) == var(i)) ) THEN
        IF ( use_file(i) ) THEN
          CALL log_info("init_top",                                           &
                        "'" // TRIM(var(i)) // "' will be read from file")

          ! If the variable will be filled from file, register it here
          nvars_file = nvars_file + 1
          ! Since nvars_file <= i (so we will not overwrite unprocessed values)
          ! and we do not need the values from these arrays for any non-file
          !variables from now on, we can just compress them down onto variables
          !that are in the file
          var(nvars_file) = var(i)
          var_name(nvars_file) = var_name(i)
          tpl_name(nvars_file) = tpl_name(i)
        ELSE
          ! If the variable is being set as a constant, just populate it here
          CALL log_info("init_top",                                           &
                        "'" // TRIM(var(i)) // "' will be set to a " //       &
                        "constant = " // to_string(const_val(i)))

          CALL populate_var(get_var_id(var(i)), CONST_VAL=const_val(i))
        END IF
      ELSE
        ! If the variable is not a required variable, warn about not using it
        CALL log_warn("init_top",                                             &
                      "Provided variable '" // TRIM(var(i)) //                &
                      "' is not required, so will be ignored")
      END IF
    END DO

    !---------------------------------------------------------------------------
    ! Set variables from file
    !---------------------------------------------------------------------------
      IF ( nvars_file > 0 ) THEN
    !   Check that a file name was provided.
      IF ( LEN_TRIM(file) == 0 )                                              &
        CALL log_fatal("init_top", "No file name provided")

      IF ( tpl_has_var_name(file) ) THEN
    ! We are using a file name template, so loop through the variables setting
    ! one from each file
        DO i = 1,nvars_file
    ! If using a variable name template, check that a template string was
    !provided for the current variable
          IF ( LEN_TRIM(tpl_name(i)) == 0 )                                   &
            CALL log_fatal("init_top",                                        &
                           "No variable name template substitution " //       &
                           "(tpl_name) provided for " // TRIM(var(i)))

          CALL fill_variables_from_file(                                      &
            tpl_substitute_var(file, tpl_name(i)),                            &
            (/ var(i) /), (/ var_name(i) /)                                   &
          )
        END DO
      ELSE
    ! We are not using a file name template, so set all variables from the same
    ! file
        CALL fill_variables_from_file(                                        &
          file, var(1:nvars_file), var_name(1:nvars_file)                     &
        )
      END IF
    END IF

  ELSE !We read from the dump file
    CALL log_info("init_top",                                                 &
                  "topmodel ancils will be read from the dump file.  " //     &
                  "Namelist values ignored")

  END IF !.NOT. ancil_dump_read%frac

  RETURN

END SUBROUTINE init_top
#endif
