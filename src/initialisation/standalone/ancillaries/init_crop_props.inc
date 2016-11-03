#if !defined(UM_JULES)

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Reads in JULES_CROP_PROPS namelist if necessary

SUBROUTINE init_crop_props()

  USE io_constants, ONLY : MAX_FILE_NAME_LEN, MAX_SDF_NAME_LEN, NAMELIST_UNIT
  
  USE templating_mod, ONLY : tpl_has_var_name, tpl_substitute_var

  USE model_interface_mod, ONLY : IDENTIFIER_LEN, populate_var, get_var_id
  
  USE string_utils_mod, ONLY : to_string

  USE jules_vegetation_mod, ONLY : l_crop, l_prescsow

  USE dump_mod, ONLY: ancil_dump_read
  
  USE errormessagelength_mod, ONLY: errormessagelength
  
  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!    Reads in TTveg, TTrep
!    Reads in the sowing date for each of the crop pfts if l_prescsow=T
!
! Current Code Owner:  Tom Osborne
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
  INTEGER, PARAMETER :: MAX_CROP_VARS = 3
         ! The maximum possible number of crop variables that can be given
         
  INTEGER :: nvars_required      ! The number of variables that are
                                 ! required in this configuration
         
  CHARACTER(len=IDENTIFIER_LEN), ALLOCATABLE :: required_vars(:)
                                 ! The variable identifiers of the required
                                 ! variables

  INTEGER :: nvars_file       ! The number of variables that will be set
                              ! from the given file (template?)                                
                                 
  INTEGER :: error  ! Error indicator

  INTEGER :: i  ! Loop counter
  
!-----------------------------------------------------------------------------
! Definition of the jules_crop_props namelist
!-----------------------------------------------------------------------------
  LOGICAL :: read_from_dump
  CHARACTER(len=MAX_FILE_NAME_LEN) :: file
                        ! The name of the file (or variable name template) to
                        ! use for variables that need to be filled from file
  INTEGER :: nvars      ! The number of variables in this section
  CHARACTER(len=IDENTIFIER_LEN) :: var(MAX_CROP_VARS)
                        ! The variable identifiers of the variables
  CHARACTER(LEN=errormessagelength) :: iomessage
                        ! Error message string for I/O errors
  LOGICAL :: use_file(MAX_CROP_VARS)
                        !   T - the variable uses the file
                        !   F - the variable is set using a constant value
  CHARACTER(len=MAX_SDF_NAME_LEN) :: var_name(MAX_CROP_VARS)
                        ! The name of each variable in the file
  CHARACTER(len=MAX_SDF_NAME_LEN) :: tpl_name(MAX_CROP_VARS)
                        ! The name to substitute in a template for each
                        ! variable
  REAL :: const_val(MAX_CROP_VARS)
                        ! The constant value to use for each variable if
                        ! use_file = F for that variable
  NAMELIST /jules_crop_props/                                                 &
                       read_from_dump, file, nvars, var, use_file, var_name,  &
                       tpl_name, const_val

!-----------------------------------------------------------------------------

! Nothing to do if crop model is not on
  IF ( .NOT. l_crop ) RETURN
  
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
  CALL log_info("init_crop_props", "Reading JULES_CROP_PROPS namelist...")

  READ(NAMELIST_UNIT, nml=jules_crop_props, IOSTAT=error, IOMSG=iomessage)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_crop_props",                                         &
                  "Error reading namelist JULES_CROP_PROPS " //               &
                  "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //        &
                  TRIM(iomessage) // ")")

  ancil_dump_read%crop_props = read_from_dump

  IF ( .NOT. ancil_dump_read%crop_props) THEN !we read from the ancil file

    !---------------------------------------------------------------------------
    ! Set up crop properties using namelist values
    !---------------------------------------------------------------------------
    ! Set up the required variables
    IF (l_prescsow) THEN
      nvars_required = MAX_CROP_VARS
      ALLOCATE(required_vars(nvars_required))
      required_vars(:) = (/ 'cropsowdate', 'cropttveg  ', 'cropttrep  ' /)
    ELSE
      nvars_required = MAX_CROP_VARS - 1
      ALLOCATE(required_vars(nvars_required))
      required_vars(:) = (/ 'cropttveg', 'cropttrep' /)
    END IF

    !---------------------------------------------------------------------------
    ! Check that all the required variables are there
    !---------------------------------------------------------------------------
    DO i = 1,nvars_required
      IF ( .NOT. ANY(var(1:nvars) == required_vars(i)) )                      &
        CALL log_fatal("init_crop_props",                                     &
                       "No value given for required variable '" //            &
                       TRIM(required_vars(i)) // "'")
    END DO

    !---------------------------------------------------------------------------
    ! Check which variables we will be using and partition them into variables
    ! set to constant values and variables set from file
    !---------------------------------------------------------------------------
    DO i = 1,nvars
      !-------------------------------------------------------------------------
      ! If the variable is one of the required vars, then we will be using it
      !-------------------------------------------------------------------------
      IF ( ANY(required_vars(1:nvars_required) == var(i)) ) THEN
        IF ( use_file(i) ) THEN
          CALL log_info("init_crop_props",                                    &
                        "'" // TRIM(var(i)) // "' will be read from file")

          ! If the variable will be filled from file, register it here
          nvars_file = nvars_file + 1
          ! Since nvars_file <= i (so we will not overwrite unprocessed values)
          ! and we do not need the values from these arrays for any non-file
          ! variables from now on, we can just compress them down onto variables
          !that are in the file
          var(nvars_file) = var(i)
          var_name(nvars_file) = var_name(i)
          tpl_name(nvars_file) = tpl_name(i)
        ELSE
          ! If the variable is being set as a constant, just populate it here
          CALL log_info("init_crop_props",                                    &
                        "'" // TRIM(var(i)) // "' will be set to a " //       &
                        "constant = " // to_string(const_val(i)))

          CALL populate_var(get_var_id(var(i)), CONST_VAL=const_val(i))
        END IF
      ELSE
        ! If the variable is not a required variable, warn about not using it
        CALL log_warn("init_crop_props",                                      &
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
        CALL log_fatal("init_crop_props", "No file name provided")

      IF ( tpl_has_var_name(file) ) THEN
        ! We are using a file name template, so loop through the variables
        !setting one from each file
        DO i = 1,nvars_file
          ! If using a variable name template, check that a template string was
          !provided for the current variable
          IF ( LEN_TRIM(tpl_name(i)) == 0 )                                   &
            CALL log_fatal("init_crop_props",                                 &
                           "No variable name template substitution " //       &
                           "(tpl_name) provided for " // TRIM(var(i)))
  
          CALL fill_variables_from_file(                                      &
            tpl_substitute_var(file, tpl_name(i)),                            &
            (/ var(i) /), (/ var_name(i) /)                                   &
          )
        END DO
      ELSE
        ! We are not using a file name template, so set all variables from the
        ! same file
        CALL fill_variables_from_file(                                        &
          file, var(1:nvars_file), var_name(1:nvars_file)                     &
        )
      END IF
    END IF

    DEALLOCATE(required_vars)

  ELSE !We read from the dump file
    CALL log_info("init_crop_props",                                          &
                  "crop properties will be read from the dump file.  " //     &
                  "Namelist values ignored")

  END IF !.NOT. ancil_dump_read%crop_props
  
END SUBROUTINE init_crop_props
#endif