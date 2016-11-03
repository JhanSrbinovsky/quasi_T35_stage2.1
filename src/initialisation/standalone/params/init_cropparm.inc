#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

SUBROUTINE init_cropparm(nml_dir)

  USE missing_data_mod, ONLY : rmdi

  USE io_constants, ONLY : NAMELIST_UNIT

  USE string_utils_mod, ONLY : to_string

  USE jules_surface_types_mod, ONLY : ncpft, nnpft
  
  USE cropparm
  
  USE cropparm_io
  
  USE jules_vegetation_mod, ONLY : l_crop
  
  USE errormessagelength_mod, ONLY: errormessagelength

  IMPLICIT NONE
  
!-----------------------------------------------------------------------------
! Description:
!   Initialises the crop PFT parameters
!
! Current Code Owner: Tom Osborne
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!
!-----------------------------------------------------------------------------
! Arguments
  CHARACTER(len=*), INTENT(IN) :: nml_dir  ! The directory containing the
                                           ! namelists

! Work variables
  INTEGER :: error  ! Error indicator
  CHARACTER(LEN=errormessagelength) :: iomessage

  INTEGER :: i ! WORK Loop counter

!----------------------------------------------------------------------------

! Nothing to do if crop model is not selected
  IF ( .NOT. l_crop ) RETURN

!-----------------------------------------------------------------------------
! Initialise namelist values before reading them
!-----------------------------------------------------------------------------

  T_BSE_io(:)     = RMDI
  T_OPT_io(:)     = RMDI
  T_MAX_io(:)     = RMDI
  TT_EMR_io(:)    = RMDI
  
  CRIT_PP_io(:)   = RMDI 
  PP_SENS_io(:)   = RMDI
  RT_DIR_io(:)    = RMDI
  ALPHA1_io(:)    = RMDI
  
  ALPHA2_io(:)    = RMDI
  ALPHA3_io(:)    = RMDI
  BETA1_io(:)     = RMDI
  BETA2_io(:)     = RMDI
  BETA3_io(:)     = RMDI
  
  GAMMA_io(:)     = RMDI
  DELTA_io(:)     = RMDI
  REMOB_io(:)     = RMDI
  CFRAC_S_io(:)   = RMDI
  CFRAC_R_io(:)   = RMDI
  
  CFRAC_L_io(:)   = RMDI
  ALLO1_io(:)     = RMDI
  ALLO2_io(:)     = RMDI
  
  mu_io(:)             = RMDI
  nu_io(:)             = RMDI
  yield_frac_io(:)     = RMDI
  initial_carbon_io(:) = RMDI
  sen_dvi_io(:)        = RMDI
  t_mort_io(:)         = RMDI
  
!-----------------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------------
  CALL log_info("init_cropparm", "Reading JULES_CROPPARM namelist...")

! Open the crop pft parameters namelist file
  OPEN(NAMELIST_UNIT, FILE=(TRIM(nml_dir) // '/' // 'crop_params.nml'),       &
                 STATUS='old', POSITION='rewind', ACTION='read', IOSTAT=error,&
                 IOMSG=iomessage)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_cropparm",                                           &
                   "Error opening namelist file crop_params.nml " //          &
                   "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //       &
                   TRIM(iomessage) // ")")

  READ(NAMELIST_UNIT, nml=jules_cropparm, IOSTAT=error, IOMSG=iomessage)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_cropparm",                                           &
                   "Error reading namelist JULES_CROPPARM " //                &
                   "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //       &
                   TRIM(iomessage) // ")")

! Close the namelist file
  CLOSE(NAMELIST_UNIT, IOSTAT=error, IOMSG=iomessage)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_cropparm",                                           &
                   "Error closing namelist file crop_params.nml " //          &
                   "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //       &
                   TRIM(iomessage) // ")")


!-----------------------------------------------------------------------------
! Process the namelist values
!-----------------------------------------------------------------------------
! Copy values from fixed length arrays used in namelist into allocated arrays  

  T_BSE(:)   = T_BSE_io(1:ncpft)
  T_OPT(:)   = T_OPT_io(1:ncpft)
  T_MAX(:)   = T_MAX_io(1:ncpft)
  TT_EMR(:)  = TT_EMR_io(1:ncpft)
  CRIT_PP(:) = CRIT_PP_io(1:ncpft)
  PP_SENS(:) = PP_SENS_io(1:ncpft)
  RT_DIR(:)  = RT_DIR_io(1:ncpft)
  ALPHA1(:)  = ALPHA1_io(1:ncpft)
  ALPHA2(:)  = ALPHA2_io(1:ncpft)
  ALPHA3(:)  = ALPHA3_io(1:ncpft)
  BETA1(:)   = BETA1_io(1:ncpft)
  BETA2(:)   = BETA2_io(1:ncpft)
  BETA3(:)   = BETA3_io(1:ncpft)
  GAMMA(:)   = GAMMA_io(1:ncpft)
  DELTA(:)   = DELTA_io(1:ncpft)
  REMOB(:)   = REMOB_io(1:ncpft)
  CFRAC_S(:) = CFRAC_S_io(1:ncpft)
  CFRAC_R(:) = CFRAC_R_io(1:ncpft)
  CFRAC_L(:) = CFRAC_L_io(1:ncpft)
  ALLO1(:)   = ALLO1_io(1:ncpft)
  ALLO2(:)   = ALLO2_io(1:ncpft)
  mu(:)             = mu_io(1:ncpft)
  nu(:)             = nu_io(1:ncpft)
  yield_frac(:)     = yield_frac_io(1:ncpft)  
  initial_carbon(:) = initial_carbon_io(1:ncpft)
  sen_dvi(:)        = sen_dvi_io(1:ncpft)
  t_mort(:)         = t_mort_io(1:ncpft)

!-----------------------------------------------------------------------------
! Check that all required variables were present in the namelist.
! The namelist variables were initialised to RMDI.
!-----------------------------------------------------------------------------  
  error = 0  

  IF ( ANY( ABS( T_BSE(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_cropparm", "No value for T_BSE")
  END IF
  IF ( ANY( ABS( T_OPT(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_cropparm", "No value for T_OPT")
  END IF
  IF ( ANY( ABS( T_MAX(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_cropparm", "No value for T_MAX")
  END IF
  IF ( ANY( ABS( TT_EMR(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_cropparm", "No value for TT_EMR")
  END IF  
  
  
  IF ( ANY( ABS( CRIT_PP(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_cropparm", "No value for CRIT_PP")
  END IF
  IF ( ANY( ABS( PP_SENS(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_cropparm", "No value for PP_SENS")
  END IF
  IF ( ANY( ABS( RT_DIR(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_cropparm", "No value for RT_DIR")
  END IF  
  IF ( ANY( ABS( ALPHA1(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_cropparm", "No value for ALPHA1")
  END IF  
  

  IF ( ANY( ABS( ALPHA2(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_cropparm", "No value for ALPHA2")
  END IF
  IF ( ANY( ABS( ALPHA3(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_cropparm", "No value for ALPHA3")
  END IF
  IF ( ANY( ABS( BETA1(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_cropparm", "No value for BETA1")
  END IF
  IF ( ANY( ABS( BETA2(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_cropparm", "No value for BETA2")
  END IF  
  IF ( ANY( ABS( BETA3(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_cropparm", "No value for BETA3")
  END IF
  

  IF ( ANY( ABS( GAMMA(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_cropparm", "No value for GAMMA")
  END IF
  IF ( ANY( ABS( DELTA(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_cropparm", "No value for DELTA")
  END IF
  IF ( ANY( ABS( REMOB(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_cropparm", "No value for REMOB")
  END IF
  IF ( ANY( ABS( CFRAC_S(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_cropparm", "No value for CFRAC_S")
  END IF  
  IF ( ANY( ABS( CFRAC_R(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_cropparm", "No value for CFRAC_R")
  END IF  


  IF ( ANY( ABS( CFRAC_L(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_cropparm", "No value for CFRAC_L")
  END IF
  IF ( ANY( ABS( ALLO1(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_cropparm", "No value for ALLO1")
  END IF
  IF ( ANY( ABS( ALLO2(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_cropparm", "No value for ALLO2")
  END IF
  

  IF ( ANY( ABS( mu(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_cropparm", "No value for mu")
  END IF
  IF ( ANY( ABS( nu(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_cropparm", "No value for nu")
  END IF

  IF ( ANY( ABS( yield_frac(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_cropparm", "No value for yield_frac")
  END IF
  IF ( ANY( ABS( initial_carbon(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_cropparm", "No value for initial_carbon")
  END IF
  IF ( ANY( ABS( sen_dvi(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_cropparm", "No value for sen_dvi")
  END IF
  IF ( ANY( ABS( t_mort(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_cropparm", "No value for t_mort")
  END IF


  IF ( error /= 0 )                                                           &                                       
    CALL log_fatal("init_cropparm",                                           &
                   "Variable(s) missing from namelist - see earlier " //      &
                   "error message(s)")
                   
  RETURN
  
END SUBROUTINE init_cropparm
#endif

