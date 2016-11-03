#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
MODULE init_soil_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE init_soil(nml_dir)
!-----------------------------------------------------------------------------
! Description:
!   Reads in the soil namelist items and checks them for consistency
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
  USE io_constants, ONLY : NAMELIST_UNIT

  USE string_utils_mod, ONLY : to_string

  USE jules_soil_mod

  USE jules_surface_mod, ONLY: l_elev_land_ice

  USE logging_mod, ONLY : log_info, log_error, log_fatal

  USE errormessagelength_mod, ONLY: errormessagelength

  IMPLICIT NONE

! Arguments
  CHARACTER(len=*), INTENT(IN) :: nml_dir  ! The directory containing the
                                           ! namelists
! Work variables
  INTEGER :: error  ! Error indicator
  CHARACTER(LEN=errormessagelength) :: iomessage

!-----------------------------------------------------------------------------
! First, we read the soil namelist
!-----------------------------------------------------------------------------
  CALL log_info("init_soil", "Reading JULES_SOIL namelist...")

  OPEN(NAMELIST_UNIT, FILE=(TRIM(nml_dir) // '/' // 'jules_soil.nml'),        &
       STATUS='old', POSITION='rewind', ACTION='read', IOSTAT=error,          &
       IOMSG=iomessage)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_soil",                                               &
                   "Error opening namelist file jules_soil.nml " //           &
                   "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //       &
                   TRIM(iomessage) // ")")

  READ(NAMELIST_UNIT, nml=jules_soil, IOSTAT=error, IOMSG=iomessage)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_soil",                                               &
                   "Error reading namelist JULES_SOIL " //                    &
                   "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //       &
                   TRIM(iomessage) // ")")

  CLOSE(NAMELIST_UNIT, IOSTAT=error, IOMSG=iomessage)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_soil",                                               &
                   "Error closing namelist file jules_soil.nml " //           &
                   "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //       &
                   TRIM(iomessage) // ")")

  CALL check_jules_soil(sm_levels)

!-----------------------------------------------------------------------------
! Print some human friendly summary information about the selected options
!-----------------------------------------------------------------------------
  CALL log_info("init_soil",                                                  &
                "Using " // TRIM(to_string(sm_levels)) // " soil levels")

  CALL log_info("init_soil",                                                  &
                "Soil levels: " // TRIM(to_string(dzsoil)))

  IF ( l_elev_land_ice )                                                      &
  CALL log_info("init_soil",                                                  &
              "Tiled ice subsurface depth: " // TRIM(to_string(dzsoil_elev)))

  IF ( l_vg_soil )                                                            &
    CALL log_info("init_soil", "van Genuchten model will be used")

  IF ( l_soil_sat_down ) THEN
    CALL log_info("init_soil",                                                &
                  "l_soil_sat_down = T - excess water is pushed down")
  ELSE
    CALL log_info("init_soil",                                                &
                  "l_soil_sat_down = F - excess water is pushed up")
  ENDIF

  IF ( soilHc_method == 1 ) THEN
    CALL log_info("init_soil",                                                &
                  "soilHc_method = 1 - following Cox et al (1999)")
  ELSE IF ( soilHc_method == 2 ) THEN
    CALL log_info("init_soil",                                                &
                  "soilHc_method = 2 - following simplified Johansen (1975)")
  ELSE
    CALL log_info("init_soil",                                                &
                  "soilHc_method = 3 - Chadburn et al. (2015)")
  ENDIF

  IF ( l_bedrock )                                                            &
    CALL log_info("init_soil", "Bedrock will be included at base of soil")

  RETURN
END SUBROUTINE init_soil

END MODULE init_soil_mod
#endif
