#if !defined(UM_JULES)
! This code is designed to make a quick cross-check to ensure that
! the flags in IMOGEN are set properly for the currently allowed configu
! The currently allowed configurations should fit with the available doc
! C. Huntingford (4th March 2004)
  SUBROUTINE IMOGEN_CHECK(                                                    &
    C_EMISSIONS,INCLUDE_CO2,INCLUDE_NON_CO2,LAND_FEED,                        &
    OCEAN_FEED,ANLG,ANOM                                                      &
  )

    USE logging_mod, ONLY : log_fatal

    USE jules_print_mgr, ONLY :                                               &
      jules_message,                                                          &
      jules_print
    IMPLICIT NONE

    LOGICAL ::                                                                &
      C_EMISSIONS,                                                            &
             !IN If true, means CO2 concentration is calcula
      INCLUDE_CO2,                                                            &
             !IN Are adjustments to CO2 values allowed?
      INCLUDE_NON_CO2,                                                        &
             !IN Are adjustments to non-CO2 values allowed?
      LAND_FEED,                                                              &
             !IN Are land feedbacks allowed on atmospheric C
      OCEAN_FEED,                                                             &
             !IN Are ocean feedbacks allowed on atmospheric
      CHECK_FLAG,                                                             &
             !IN Confirms that configuration is OK
      ANLG,                                                                   &
             !IN True if the analogue model is used
      ANOM   !IN True if the anomalies are used



    CHECK_FLAG=.FALSE.

! Now make checks on structures which will hopefully all be available at
! To make this code easy to read, simple "if"
! statements comments are used for each instance and case printed.

! Spin-up
    IF(.NOT.ANOM) THEN
      WRITE(jules_message,*) 'Run is a spin-up'
      CALL jules_print('imogen_check',jules_message)
      WRITE(jules_message,*) 'Are flags set OK for any future transient run?'
      CALL jules_print('imogen_check',jules_message)
      CHECK_FLAG=.TRUE.
    ENDIF

! User prescribed anomalies
    IF(ANOM .AND. (.NOT.ANLG)) THEN
      WRITE(jules_message,*) 'Run is reading in user anomalies'
      CALL jules_print('imogen_check',jules_message)
      CHECK_FLAG=.TRUE.
    ENDIF

! Currently coded analogue model possibilities
    IF(ANOM .AND. ANLG) THEN
      IF((.NOT.INCLUDE_CO2) .AND. INCLUDE_NON_CO2) THEN
        WRITE(jules_message,*) 'Run is for prescribed non-CO2 gases only'
        CALL jules_print('imogen_check',jules_message)
        CHECK_FLAG=.TRUE.
      ENDIF

      IF(INCLUDE_CO2) THEN
        IF(.NOT.C_EMISSIONS) THEN
          IF(INCLUDE_NON_CO2) THEN
            WRITE(jules_message,*) 'Run is for prescribed CO2 and non-CO2 gases'
            CALL jules_print('imogen_check',jules_message)
            CHECK_FLAG=.TRUE.
          ELSE
            WRITE(jules_message,*) 'Run is for prescribed CO2 only'
            CALL jules_print('imogen_check',jules_message)
            CHECK_FLAG=.TRUE.
          ENDIF
        ENDIF

        IF(C_EMISSIONS) THEN
          IF(.NOT.LAND_FEED .AND. .NOT.OCEAN_FEED) THEN
            IF(.NOT.INCLUDE_NON_CO2) THEN
              WRITE(jules_message,*) 'Run is for CO2 emissions only'
              CALL jules_print('imogen_check',jules_message)
              WRITE(jules_message,*) 'There are NO land/ocean feedbacks'
              CALL jules_print('imogen_check',jules_message)
              CHECK_FLAG=.TRUE.
            ELSE
              WRITE(jules_message,*) 'Run is for CO2 emissions along with ' //&
                      'prescribed values of non-CO2 gases'
              CALL jules_print('imogen_check',jules_message)
              WRITE(jules_message,*) 'There are NO land/ocean feedbacks'
              CALL jules_print('imogen_check',jules_message)
              CHECK_FLAG=.TRUE.
            ENDIF
          ENDIF

          IF(LAND_FEED.AND.(.NOT.OCEAN_FEED)) THEN
            IF(.NOT.INCLUDE_NON_CO2) THEN
              WRITE(jules_message,*) 'Run is for CO2 emissions only'
              CALL jules_print('imogen_check',jules_message)
              WRITE(jules_message,*) 'There are land but no ocean feedbacks'
              CALL jules_print('imogen_check',jules_message)
              CHECK_FLAG=.TRUE.
            ELSE
              WRITE(jules_message,*) 'Run is for CO2 emissions along with ' //&
                      'prescribed values of non-CO2 gases'
              CALL jules_print('imogen_check',jules_message)
              WRITE(jules_message,*) 'There are land but no ocean feedbacks'
              CALL jules_print('imogen_check',jules_message)
              CHECK_FLAG=.TRUE.
            ENDIF
          ENDIF

          IF(.NOT.LAND_FEED .AND. OCEAN_FEED) THEN
            IF(.NOT.INCLUDE_NON_CO2) THEN
              WRITE(jules_message,*) 'Run is for CO2 emissions only'
              CALL jules_print('imogen_check',jules_message)
              WRITE(jules_message,*) 'There are ocean but no land feedbacks'
              CALL jules_print('imogen_check',jules_message)
              CHECK_FLAG=.TRUE.
            ELSE
              WRITE(jules_message,*) 'Run is for CO2 emissions along with ' //&
                      'prescribed values of non-CO2 gases'
              CALL jules_print('imogen_check',jules_message)
              WRITE(jules_message,*) 'There are ocean but no land feedbacks'
              CALL jules_print('imogen_check',jules_message)
              CHECK_FLAG=.TRUE.
            ENDIF
          ENDIF

          IF(LAND_FEED .AND. OCEAN_FEED) THEN
            IF(.NOT.INCLUDE_NON_CO2) THEN
              WRITE(jules_message,*) 'Run is for CO2 emissions only'
              CALL jules_print('imogen_check',jules_message)
              WRITE(jules_message,*) 'There are ocean and land feedbacks'
              CALL jules_print('imogen_check',jules_message)
              CHECK_FLAG=.TRUE.
            ELSE
              WRITE(jules_message,*) 'Run is for CO2 emissions along with ' //&
                      'prescribed values of non-CO2 gases'
              CALL jules_print('imogen_check',jules_message)
              WRITE(jules_message,*) 'There are ocean and land feedbacks'
              CALL jules_print('imogen_check',jules_message)
              CHECK_FLAG=.TRUE.
            ENDIF
          ENDIF
! End of check over different combinations with prescribed carbon
! emissions.
        ENDIF
! End of check if CO2 is included.
      ENDIF
! End of check if analogue model used.
    ENDIF

    IF(.NOT.CHECK_FLAG)                                                       &
      CALL log_fatal("IMOGEN_CHECK",                                          &
                     'Combination not yet allowed')

    RETURN

  END SUBROUTINE IMOGEN_CHECK
#endif
