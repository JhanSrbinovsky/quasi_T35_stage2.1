#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE init_soil_props()

  USE io_constants, ONLY : MAX_SDF_NAME_LEN, MAX_FILE_NAME_LEN, NAMELIST_UNIT

  USE string_utils_mod, ONLY : to_string

  USE templating_mod, ONLY : tpl_has_var_name, tpl_substitute_var

  USE model_interface_mod, ONLY : IDENTIFIER_LEN, populate_var, get_var_id

  USE dump_mod, ONLY: ancil_dump_read
  
  USE errormessagelength_mod, ONLY: errormessagelength

  USE jules_vegetation_mod, ONLY : l_triffid

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads in the soil properties and checks them for consistency
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
  INTEGER, PARAMETER :: MAX_SOIL_VARS = 10 ! The maximum possible number of
                                           ! soil variables that can be given

  INTEGER :: nvars_required      ! The number of soil variables that are
                                 ! required in this configuration
  INTEGER :: nvars_optional      ! The number of optional soil variables 
                                 ! in this configuration
  CHARACTER(len=IDENTIFIER_LEN) :: required_vars(MAX_SOIL_VARS)
                                 ! The variable identifiers of the required
                                 ! variables
  CHARACTER(len=IDENTIFIER_LEN) :: optional_vars(MAX_SOIL_VARS)
                                 ! The variable identifiers of any optional
                                 ! variables

  INTEGER :: nvars_file       ! The number of variables that will be set
                              ! from the given file (template?)
! Variables passed to fill_variables_from_file
  CHARACTER(len=IDENTIFIER_LEN) :: file_var(MAX_SOIL_VARS)
                        ! The variable identifiers of the variables to set
                        ! from file
  CHARACTER(len=MAX_SDF_NAME_LEN) :: file_var_name(MAX_SOIL_VARS)
                        ! The name of each variable in the file
  CHARACTER(len=MAX_SDF_NAME_LEN) :: file_tpl_name(MAX_SOIL_VARS)
                        ! The name to substitute in a template for each
                        ! variable

  INTEGER :: i    ! Index variables

  INTEGER :: error  ! Error indicator

!-----------------------------------------------------------------------------
! Definition of the jules_soil_props namelist
!-----------------------------------------------------------------------------
  LOGICAL :: read_from_dump
  LOGICAL :: const_z            ! T - the same properties are used for each
                                !     soil layer
                                ! F - properties for each layer are read from
                                !     file
  CHARACTER(len=MAX_FILE_NAME_LEN) :: file
                        ! The name of the file (or variable name template) to
                        ! use for variables that need to be filled from file

  INTEGER :: nvars      ! The number of variables in this section
  CHARACTER(len=IDENTIFIER_LEN) :: var(MAX_SOIL_VARS)
                        ! The variable identifiers of the variables
  LOGICAL :: use_file(MAX_SOIL_VARS)
                        !   T - the variable uses the file
                        !   F - the variable is set using a constant value
  CHARACTER(len=MAX_SDF_NAME_LEN) :: var_name(MAX_SOIL_VARS)
                        ! The name of each variable in the file
  CHARACTER(len=MAX_SDF_NAME_LEN) :: tpl_name(MAX_SOIL_VARS)
                        ! The name to substitute in a template for each
                        ! variable
  CHARACTER(LEN=errormessagelength) :: iomessage
                        ! I/O error message string
  REAL :: const_val(MAX_SOIL_VARS)
                        ! The constant value to use for each variable if
                        ! use_file = F for that variable
  NAMELIST /jules_soil_props/ read_from_dump, const_z, file, nvars, var,      &
                              use_file, var_name, tpl_name, const_val

!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Initialise
!-----------------------------------------------------------------------------
  read_from_dump = .FALSE.
  nvars_required = 0
  nvars_optional = 0  
  nvars_file     = 0
  const_z        = .FALSE. ! Default is to read a value for each soil level
  nvars          = 0
  use_file(:)    = .TRUE.  ! Default is for every variable to be read from file
  file           = ''      ! Empty file name.

!-----------------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------------
  CALL log_info("init_soil", "Reading JULES_SOIL_PROPS namelist...")

! Read the soil properties namelist
  READ(NAMELIST_UNIT, nml=jules_soil_props, IOSTAT=error, IOMSG=iomessage)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_soil",                                               &
                   "Error reading namelist JULES_SOIL_PROPS " //              &
                   "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //       &
                   TRIM(iomessage) // ")")

  ancil_dump_read%soil_props = read_from_dump

  IF ( .NOT. ancil_dump_read%soil_props) THEN !we read from the ancil file

    !---------------------------------------------------------------------------
    ! Set up soil properties using namelist values
    !---------------------------------------------------------------------------
    ! Set up the required variables

    IF ( l_triffid ) THEN  
     !In older versions the soil clay content was set to zero instead of being read in
     !If it is not available in ancullaries give stern warning and set it to zero
      IF ( .NOT. ANY(var(1:nvars) == 'clay   ') )  THEN
        CALL log_warn("init_soil",                                            &
                       "No value given for soil clay content " //             &
                       "soil clay content is required with TRIFFID ON" //     &
                       "It is set to 0.0 as for previous versions" //         &
                       "This is WRONG - please try and find values for clay") 
        nvars_required = MAX_SOIL_VARS - 1
        required_vars(1:nvars_required) = (/ 'b      ',                       & 
                            'sathh  ', 'satcon ',                             &
                            'sm_sat ', 'sm_crit', 'sm_wilt',                  &
                            'hcap   ', 'hcon   ', 'albsoil' /)
        nvars_optional = 1
        optional_vars(1:nvars_optional) = (/ 'clay   ' /)
        CALL populate_var(get_var_id(optional_vars(1)), CONST_VAL=0.0)
      ELSE
        ! All the soil variables are required
        nvars_required = MAX_SOIL_VARS
        required_vars(1:nvars_required) = (/ 'b      ',                       &
                            'sathh  ', 'satcon ',                             &
                            'sm_sat ', 'sm_crit', 'sm_wilt',                  &
                            'hcap   ', 'hcon   ', 'albsoil',                  &
                            'clay   ' /)
       nvars_optional = 0
     END IF
   ELSE
      ! All the soil variables are required except for clay content
      nvars_required = MAX_SOIL_VARS - 1
      required_vars(1:nvars_required) = (/ 'b      ',                         & 
                            'sathh  ', 'satcon ',                             &
                            'sm_sat ', 'sm_crit', 'sm_wilt',                  &
                            'hcap   ', 'hcon   ', 'albsoil' /)
      nvars_optional = 1
      optional_vars(1:nvars_optional) = (/ 'clay   ' /)
   END IF

    !---------------------------------------------------------------------------
    ! Check that all the required variables are there
    !---------------------------------------------------------------------------
    DO i = 1,nvars_required
      IF ( .NOT. ANY(var(1:nvars) == required_vars(i)) )                      &
        CALL log_fatal("init_soil",                                           &
                       "No value given for required variable '" //            &
                       TRIM(required_vars(i)) // "'")
    END DO

    !-----------------------------------------------------------------------------
    ! Constant Z (i.e. spacially varying but constant through vertical levels) is
    ! implemented by having a separate input variable in model_interface_mod
    ! called <var>_const_z that has no vertical levels
    !
    ! Hence, once the previous check is done, we add _const_z to both required
    ! and provided variable identifiers if asked for
    !-----------------------------------------------------------------------------
    IF ( const_z ) THEN
      DO i = 1,nvars
        !albsoil, clay and albobs_*
        !are soil variables that do not have multiple levels
        IF ( var(i) /= 'albsoil' .AND. var(i) /= 'albobs_sw' .AND.            &
             var(i) /= 'albobs_vis' .AND. var(i) /= 'albobs_nir'  .AND.       &
             var(i) /= 'clay' )                                               &
          var(i) = TRIM(var(i)) // "_const_z"
      END DO

      DO i = 1,nvars_required
        !albsoil, clay and albobs_*
        !are soil variables that do not have multiple levels
        IF ( var(i) /= 'albsoil' .AND. var(i) /= 'albobs_sw' .AND.            &
             var(i) /= 'albobs_vis' .AND. var(i) /= 'albobs_nir'  .AND.       &
             var(i) /= 'clay' )                                               &
          required_vars(i) = TRIM(required_vars(i)) // "_const_z"
      END DO
      IF( nvars_optional > 0 ) THEN
        DO i = 1,nvars_optional
          !clay does not have multiple levels
          IF ( var(i) /= 'clay' )                                               &
            optional_vars(i) = TRIM(optional_vars(i)) // "_const_z"
        END DO
      END IF
    END IF

    !---------------------------------------------------------------------------
    ! Check which variables we will be using and partition them into variables
    ! set to constant values and variables set from file
    !---------------------------------------------------------------------------
    DO i = 1,nvars
    !---------------------------------------------------------------------------
    ! If the variable is one of the required vars, then we will be using it
    !---------------------------------------------------------------------------
      IF ( ANY(required_vars(1:nvars_required) == TRIM(var(i)))               &
         .OR. (nvars_optional > 0 .AND.                                       &
         ANY(optional_vars(1:nvars_optional) == TRIM(var(i))) )) THEN
        IF ( use_file(i) ) THEN
          CALL log_info("init_soil",                                          &
                        "'" // TRIM(var(i)) // "' will be read from file")

          !If the variable will be filled from file, register it here
          nvars_file = nvars_file + 1
          file_var(nvars_file) = var(i)
          file_var_name(nvars_file) = var_name(i)
          file_tpl_name(nvars_file) = tpl_name(i)
        ELSE
          !If the variable is being set as a constant, just populate it here
          CALL log_info("init_soil",                                          &
                        "'" // TRIM(var(i)) // "' will be set to a " //       &
                        "constant = " // to_string(const_val(i)))

          CALL populate_var(get_var_id(var(i)), CONST_VAL=const_val(i))
        END IF
      ELSE
  ! If the variable is not a required variable, warn about not using it
        CALL log_warn("init_soil",                                            &
                      "Provided variable '" // TRIM(var(i)) //                &
                      "' is not required, so will be ignored")
      END IF
    END DO

    !---------------------------------------------------------------------------
    ! Set variables from file
    !---------------------------------------------------------------------------
    IF ( nvars_file > 0 ) THEN
      !Check that a file name was provided
      IF ( LEN_TRIM(file) == 0 )                                              &
        CALL log_fatal("init_soil", "No file name provided")

      IF ( tpl_has_var_name(file) ) THEN
        !We are using a file name template, so loop through the variables
        !setting one from each file
        DO i = 1,nvars_file
          CALL fill_variables_from_file(                                      &
            tpl_substitute_var(file, file_tpl_name(i)),                       &
            (/ file_var(i) /), (/ file_var_name(i) /)                         &
          )
        END DO
      ELSE
        !We are not using a file name template, so set all variables from the
        !same file
        CALL fill_variables_from_file(                                        &
          file,file_var(1:nvars_file), file_var_name(1:nvars_file)            &
        )
      END IF
    END IF

  ELSE !We read from the dump file
    CALL log_info("init_soil_props",                                          &
                  "soil properties will be read from the dump file.  " //     &
                  "Namelist values ignored")

  END IF !.NOT. ancil_dump_read%soil_props

  RETURN

END SUBROUTINE init_soil_props
#endif
