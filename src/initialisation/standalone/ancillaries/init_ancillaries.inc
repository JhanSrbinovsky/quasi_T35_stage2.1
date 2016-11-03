#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE init_ancillaries(nml_dir)

  USE io_constants, ONLY : NAMELIST_UNIT

  USE string_utils_mod, ONLY : to_string

  USE ancil_info, ONLY : frac_surft

  USE aero, ONLY : co2_mmr_mod => co2_mmr

  USE dump_mod, ONLY: ancil_dump_read
  
  USE errormessagelength_mod, ONLY: errormessagelength

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads in information about the model ancillaries
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

  CHARACTER(LEN=errormessagelength) :: iomessage

! Work variables
  INTEGER :: error  ! Error indicator

  LOGICAL :: read_from_dump
  REAL    :: co2_mmr
  NAMELIST /jules_co2/ read_from_dump, co2_mmr


!-----------------------------------------------------------------------------


! Open the ancillaries namelist file
  OPEN(NAMELIST_UNIT, FILE=(TRIM(nml_dir) // '/' // 'ancillaries.nml'),       &
                 STATUS='old', POSITION='rewind', ACTION='read', IOSTAT=error,&
                 IOMSG=iomessage)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_ancillaries",                                        &
                   "Error opening namelist file ancillaries.nml " //          &
                   "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //       &
                   TRIM(iomessage) // ")")

! Defer to specialist routines to process each namelist
  CALL init_frac()
  CALL init_soil_props()
  CALL init_top()
  CALL init_agric()
  CALL init_crop_props()
  CALL init_irrig()
  CALL init_rivers_props()

! Read the JULES_CO2 namelist
! This is so simple it doesn't really need its own subroutine
  read_from_dump = .FALSE.

  !Copy across the default value set up in the module.
  co2_mmr = co2_mmr_mod

  CALL log_info("init_ancillaries", "Reading JULES_CO2 namelist...")
  READ(NAMELIST_UNIT, nml=jules_co2, IOSTAT=error, IOMSG=iomessage)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_ancillaries",                                        &
                   "Error reading namelist JULES_CO2 " //                     &
                   "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //       &
                   TRIM(iomessage) // ")")

  read_from_dump = ancil_dump_read%co2

  IF ( .NOT. ancil_dump_read%co2) THEN
    co2_mmr_mod = co2_mmr
  ELSE !We read from the dump file
    CALL log_info("init_frac",                                                &
                  "co2_mmr will be read from the dump file.  " //             &
                  "Namelist value ignored")
  END IF !.NOT. ancil_dump_read%frac

! Close the namelist file
  CLOSE(NAMELIST_UNIT, IOSTAT=error, IOMSG=iomessage)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_ancillaries",                                        &
                   "Error closing namelist file ancillaries.nml " //          &
                   "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //       &
                   TRIM(iomessage) // ")")
  
  RETURN

END SUBROUTINE init_ancillaries
#endif
