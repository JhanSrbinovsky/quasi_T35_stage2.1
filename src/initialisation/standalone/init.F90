#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE init(nml_dir)

  USE time_varying_input_mod, ONLY : seek_all_to_current_datetime

  USE model_time_mod, ONLY : is_spinup

  USE init_grid_mod, ONLY : init_grid

  USE init_ancillaries_mod, ONLY : init_ancillaries

  USE init_params_mod, ONLY : init_params

  USE initial_conditions_mod, ONLY : init_ic

  USE spinup_mod, ONLY : spinup_init

  USE dump_mod, ONLY : write_dump

  USE logging_mod, ONLY : log_info, log_debug, log_warn, log_error, log_fatal

  USE init_vegetation_mod, ONLY: init_vegetation

  USE init_hydrology_mod, ONLY : init_hydrology

  USE init_radiation_mod, ONLY : init_radiation

  USE init_snow_mod, ONLY : init_snow

  USE init_surface_types_mod, ONLY : init_surface_types

  USE init_soil_mod, ONLY : init_soil

  USE init_surface_mod, ONLY : init_surface

  USE init_rivers_mod, ONLY : init_rivers

  USE metstats_mod, ONLY : l_metstats
  USE metstats_init_mod, ONLY : metstats_init

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   The main initialisation routine - initialises the model by calling
!   specialised routines
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


!-----------------------------------------------------------------------------


! Initialise the surface types
  CALL init_surface_types(nml_dir)

! Initialise the surface options
  CALL init_surface(nml_dir)

! Initialise radiation options
  CALL init_radiation(nml_dir)

! Initialise hydrology options
  CALL init_hydrology(nml_dir)

! Initialise soil options
  CALL init_soil(nml_dir)

! Initialise vegetation options
  CALL init_vegetation(nml_dir)

! Initialise snow options
  CALL init_snow(nml_dir)

! Initialise river routing parameters, ancils and grid
  CALL init_rivers(nml_dir)

! Intialise the times for the run
  CALL init_time(nml_dir)

! Initialise the input, model and output grids
  CALL init_grid(nml_dir)

! Initialise the model ancils
  CALL init_ancillaries(nml_dir)

! Initialise model parameters
  CALL init_params(nml_dir)

! Initialise urban schemes
  CALL init_urban(nml_dir)

! Initialise fire module
  CALL init_fire(nml_dir)

! Initialise the metstats module (used by fire)
  IF (l_metstats) CALL metstats_init()

! Initialise meteorological forcing
  CALL init_drive(nml_dir)

! Initialise IMOGEN
  CALL init_imogen(nml_dir)

! Initialise other prescribed data
  CALL init_prescribed_data(nml_dir)

! Initialise the model prognostics
  CALL init_ic(nml_dir)

! Initialise output
  CALL init_output(nml_dir)


!-----------------------------------------------------------------------------
! Other initialisation that does not depend on user input
!-----------------------------------------------------------------------------
! Temporary (i.e. for this version) initialisation of variables.
  CALL init_vars_tmp()

! Set index arrays and initialise other variables.
  CALL init_parms()

! Seek the input files to the start of the run
  CALL seek_all_to_current_datetime()

! Save initial state if spinning up. Arrays are allocated here.
  IF ( is_spinup ) CALL spinup_init()

! Write an initial dump
  CALL write_dump()


  CALL log_info("init", "Initialisation is complete")

  RETURN

END SUBROUTINE init
#endif
