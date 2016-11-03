#if !defined(UM_JULES)
SUBROUTINE imogen_update_carb

  USE model_time_mod, ONLY : current_time, main_run_start,                    &
         main_run_end

  USE ancil_info, ONLY : dim_cs1

  USE parallel_mod, ONLY : is_master_task, gather_land_field,                 &
        scatter_land_field

  USE model_grid_mod, ONLY : global_land_pts

  USE ancil_info, ONLY : land_pts

  USE trifctl, ONLY : cv_gb

  USE prognostics, ONLY : cs_pool_gb

  USE imogen_run, ONLY :                                                      &
    INCLUDE_CO2,C_EMISSIONS,ANOM,ANLG,CO2_INIT_PPMV,LAND_FEED,                &
    NYR_EMISS,OCEAN_FEED

  USE imogen_clim, ONLY :                                                     &
    DCTOT,LAT

  USE imogen_constants, ONLY :                                                &
    CONV,OCEAN_AREA,NFARRAY

  USE imogen_progs, ONLY :                                                    &
    CO2_PPMV,CO2_START_PPMV,DTEMP_O,CO2_CHANGE_PPMV,FA_OCEAN

  USE imogen_anlg_vals, ONLY :                                                &
    T_OCEAN_INIT

  USE imogen_io_vars, ONLY :                                                  &
    YR_EMISS,C_EMISS

  USE logging_mod, ONLY : log_fatal

  USE jules_print_mgr, ONLY :                                                 &
    jules_message,                                                            &
    jules_print
  IMPLICIT NONE

  INTEGER ::                                                                  &
    EMISS_TALLY ! Checks that datafile of emissions includes
                ! year of interest

  INTEGER :: i,n ! loop counters

  REAL ::                                                                     &
    C_EMISS_LOCAL  ! Local value of C_EMISS

  REAL ::                                                                     &
    D_LAND_ATMOS,                                                             &
                ! Change in atmospheric CO2 concentration due to
                ! feedbacks (ppm/year).
    D_OCEAN_ATMOS,                                                            &
                ! Change in atmospheric CO2 concentration due to
                ! feedbacks (ppm/year)
    D_LAND_ATMOS_FIELD(LAND_PTS)
                !d_land_atmos scattered on to each processor
  REAL, ALLOCATABLE ::                                                        &
    global_data_dctot(:),                                                     &
                !  dcot on full grid
    global_data_lat(:),                                                       &
                !  lat on full grid
    global_data_tmp(:)  ! d_land_atmos on full grid
                        ! same value for each cell
!-----------------------------------------------------------------


!-----------------------------------------------------------------
! Calculate land atmosphere carbon exchange.
!-----------------------------------------------------------------
  IF(LAND_FEED) THEN
    DO I=1,LAND_PTS
      DCTOT(I) = cv_gb(I) - DCTOT(I)
      DO N=1,dim_cs1
        DCTOT(I) = DCTOT(I) + cs_pool_gb(I,N)
      ENDDO
    ENDDO

!need to gather dctot and lat to global grid
    IF ( is_master_task() ) ALLOCATE(global_data_dctot(global_land_pts))
    IF ( is_master_task() ) ALLOCATE(global_data_lat(global_land_pts))
    IF ( is_master_task() ) ALLOCATE(global_data_tmp(global_land_pts))
    CALL gather_land_field(DCTOT, global_data_dctot)
    CALL gather_land_field(LAT, global_data_lat)

    IF ( is_master_task() ) THEN
      CALL DIFFCARB_LAND(global_land_pts,D_LAND_ATMOS,                        &
                global_data_lat,global_data_dctot,CONV)
      global_data_tmp(:) = D_LAND_ATMOS
    ENDIF

!need to scatter D_LAND_ATMOS to each processor
    CALL scatter_land_field(global_data_tmp,D_LAND_ATMOS_FIELD)
    D_LAND_ATMOS = D_LAND_ATMOS_FIELD(1)

    IF ( ALLOCATED(global_data_dctot) ) DEALLOCATE(global_data_dctot)
    IF ( ALLOCATED(global_data_lat) ) DEALLOCATE(global_data_lat)
    IF ( ALLOCATED(global_data_tmp) ) DEALLOCATE(global_data_tmp)

  ENDIF

!-----------------------------------------------------------------
! Now do the carbon cycling update.
!-----------------------------------------------------------------
  IF(INCLUDE_CO2 .AND. C_EMISSIONS .AND. ANOM .AND. ANLG) THEN

! Include anthropogenic carbon emissions.
    EMISS_TALLY=0
    DO N = 1,NYR_EMISS
      IF(YR_EMISS(N) == current_time%year) THEN
        C_EMISS_LOCAL = C_EMISS(N)
        CO2_PPMV = CO2_PPMV + CONV * C_EMISS_LOCAL
        EMISS_TALLY = EMISS_TALLY + 1
! We have found the right year so we can exit the loop
        EXIT
      ENDIF
    ENDDO

    IF(EMISS_TALLY /= 1)                                                      &
      CALL log_fatal("imogen_update_carb",                                    &
                     'IMOGEN: Emission dataset does not match run')

    IF(LAND_FEED) THEN
! Update with land feedbacks if required
      CO2_PPMV = CO2_PPMV + D_LAND_ATMOS
    ENDIF

    IF(OCEAN_FEED) THEN
! Update with ocean feedbacks if required
      CALL OCEAN_CO2(                                                         &
        current_time%year-main_run_start%year+1,                              &
        1,CO2_PPMV,CO2_INIT_PPMV,DTEMP_O(1),                                  &
        FA_OCEAN,OCEAN_AREA,CO2_CHANGE_PPMV,                                  &
        main_run_end%year-main_run_start%year,                                &
        T_OCEAN_INIT,NFARRAY,D_OCEAN_ATMOS                                    &
      )
      CO2_PPMV = CO2_PPMV + D_OCEAN_ATMOS
    ENDIF
  ENDIF

  IF(INCLUDE_CO2) CO2_CHANGE_PPMV = CO2_PPMV - CO2_START_PPMV


  RETURN

END SUBROUTINE imogen_update_carb
#endif
