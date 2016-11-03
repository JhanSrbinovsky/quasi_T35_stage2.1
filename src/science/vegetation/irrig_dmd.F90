! *****************************COPYRIGHT****************************************
! (c) Crown copyright, Met Office 2011. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use
! and distribution under the JULES collaboration agreement, subject
! to the terms and conditions set out therein.
!
! [Met Office Ref SC0237]
! *****************************COPYRIGHT****************************************

! Calculate irrigation demand for each irrigated grid box.

SUBROUTINE irrig_dmd ( land_pts,sm_levels,frac_irr                            &
                       ,a_step,plant_n                                        &
                       ,sthf,smvccl,smvcst,smvcwt,sthzw                       &
                       ,sthu_irr,sthu,smcl                                    &
                       ,irr_crop,dvimax )

! Description:
!   Calculates irrigation demand over unfrozen soils as the amount of
!     water needed to alleviate soil water deficit
!
!
! Method:
!
! Notes hadrd 2011-07-20:
!   This code was originally included within control.f90
!     but is now put in a separate file
!   Instead of called every timestep from control, this could only be called
!   once per day
!
! Current Code Owner: code originally developed by Nic Gedney
!
! Code Description:
!   Language: Fortran 90.
!   This code is partially modified to JULES coding standards v1.
!   Rutger Dankers, July 2011
!
!-------------------------------------------------------------------------------


  USE c_densty, ONLY : rho_water

  USE jules_soil_mod, ONLY : dzsoil

  USE jules_vegetation_mod, ONLY : l_irrig_dmd, l_irrig_limit

  USE time_info_mod, ONLY : l_360, l_leap, current_model_time

  USE datetime_utils_mod, ONLY : day_of_year, days_in_year

  !USE model_time_mod, ONLY : timestep_len, current_time

  USE crop_vars_mod, ONLY : nday_crop, irrDaysDiag_gb, irrig_water_gb

  USE timestep_mod, ONLY : timestep

  USE ereport_mod, ONLY : ereport

  USE jules_hydrology_mod, ONLY : zw_max

  USE jules_rivers_mod, ONLY : rivers_sto_per_m2_on_landpts,                  &
      rivers_adj_on_landpts

!-------------------------------------------------------------------------------

  IMPLICIT NONE

!   Scalar arguments with intent(IN) :
  INTEGER, INTENT(IN) ::                                                      &
    land_pts                                                                  &
                      ! Number of land ice points
    ,sm_levels                                                                &
                      ! No. of soil moisture levels
    ,a_step                                                                   &
                       ! Atmospheric timestep number
                       ! necessary when called from control.f90 each timestep
    ,irr_crop
                      ! Switch for irrigation cropping model

!   Array arguments with intent(IN) :

  INTEGER, INTENT(IN) ::                                                      &
    plant_n(land_pts)
                      ! best plant date for non-rice

  REAL, INTENT(IN) ::                                                         &
    frac_irr(land_pts)                                                        &
                      ! irrigation fraction for this year

    ,sthf(land_pts,sm_levels)                                                 &
                      ! Frozen soil moisture content of
                      !    each layer as a fraction of
                      !    saturation
    ,smvccl(land_pts,sm_levels)                                               &
                      ! Critical volumetric SMC
                      !    (cubic m per cubic m of soil)
    ,smvcst(land_pts,sm_levels)                                               &
                      ! Volumetric saturation point
                      !    (m3/m3 of soil)
    ,smvcwt(land_pts,sm_levels)                                               &
                      ! Volumetric wilting point
                      !    (m3/m3 of soil)
    ,dvimax(land_pts)
                      ! Maximum DVI for crop types
! hadrd - may be replaced with dynamic irrigation fraction?


!   Array arguments with intent(INOUT) :
  REAL, INTENT(INOUT) ::                                                      &
    sthu_irr(land_pts,sm_levels)                                              &
                      ! Unfrozen soil moisture content of
                      !    each layer as a fraction of
                      !    saturation in irrigated fraction
    ,sthu(land_pts,sm_levels)                                                 &
                      ! Unfrozen soil moisture content of
                      !    each layer as a fraction of
                      !    saturation.
    ,smcl(land_pts,sm_levels)                                                 &
                      ! Soil moisture content of each
                      !    layer (kg/m2)
    ,sthzw(land_pts)
                      ! Soil moist fraction in deep layer.

!   LOCAL scalars :
  LOGICAL :: l_julian_crop
  LOGICAL , PARAMETER :: irr_zw_all = .TRUE.
                      ! Withdraw all available water from the
                      ! groundwater store (T), or withdraw
                      ! until the wilting point only (F)

  INTEGER ::                                                                  &
    l,n                                                                       &
                      ! Loop counters
    ,tspd                                                                     &
                      ! timesteps per day
    ,days_in_yr                                                               &
                      ! Total number of days in current year
    ,day_of_yr                                                                &
                      ! Current day of year
    ,day                                                                      &
                      ! Current day
    ,month                                                                    &
                      ! Current month
    ,year
                      ! Current year
  INTEGER, PARAMETER :: secs_in_day = 86400
!   LOCAL arrays :
  LOGICAL :: l_irrigated(land_pts)
                      ! flag for cells being irrigated

  REAL ::                                                                     &
    irrigwater_levels_day(land_pts,sm_levels)                                 &
                     ! Addition of irrigation water to
                     ! each soil layer (kg/m2/day)
    ,sthu_tmp                                                                 &
                     ! Dummy for STHU_IRR
    ,smclsatzw(land_pts)                                                      &
                     ! Moisture content in deep layer
                     !      at saturation (kg/m2)
    ,smclwiltzw(land_pts)                                                     &
                     ! Moisture content in deep layer
                     !      at wilting point (kg/m2)
    ,smclzw(land_pts)                                                         &
                      ! Actual moisture content
                      ! in deep layer (kg/m2)
    ,smclzw_rest(land_pts)                                                    &
                      ! Deep moisture content after
                      ! extraction of irrigwater_levels_day (kg/m2)

    ,zdepth(0:sm_levels)                                                      &
                      ! Lower soil layer boundary depth (m).
    ,irrig_sum                                                                &
                      ! Total irrigwater_levels_day over soil layers
                      ! (kg/m2)
    ,irrig_lim                                                                &
                      ! Total irrigwater_levels_day constrained by
                      ! water in deep soil (kg/m2)
    ,irrig_max                                                                &
                      ! Maximum available water for
                      ! irrigation from deep soil (kg/m2)
    ,irrig_dif                                                                &
                      ! Amount by which irrigwater_levels_day is
                      ! decreased (kg/m2)
    ,irrig_riv(land_pts)
                      ! total irrigation water extracted
                      ! from river routing storage
                      ! (kg/m2)

!-------------------------------------------------------------------------------
! End of header
!-------------------------------------------------------------------------------
  CALL current_model_time(year, month, day)

  days_in_yr = days_in_year(year, l_360, l_leap)
  day_of_yr  = day_of_year(year, month, day, l_360, l_leap)

  tspd = secs_in_day / timestep ! #timesteps per day
  l_irrigated(:) = .FALSE.

! Initialise water for irrigation to be zero:
  irrigwater_levels_day(:,:) = 0.0

! Initialise variables for constraining irrigation supply
  IF ( l_irrig_limit ) THEN

  ! Calculate depth of soil column
    zdepth(:) = 0.0
    DO n=1,sm_levels
      zdepth(n) = zdepth(n-1) + dzsoil(n)
    END DO

  ! Calculate water availability in deep GW store
    smclsatzw(:) = rho_water*smvcst(:,sm_levels) * (zw_max-zdepth(sm_levels))
    smclwiltzw(:) = rho_water*smvcwt(:,sm_levels) * (zw_max-zdepth(sm_levels))
    smclzw(:) = sthzw(:) * smclsatzw(:)
    ! initialise other variables
    smclzw_rest(:) = smclzw(:)
    rivers_adj_on_landpts(:) = 1.0
    irrig_riv(:) = 0.0

  ENDIF ! l_irrig_limit

!-------------------------------------------------------------------------------
! Calculate irrigation demand
!-------------------------------------------------------------------------------

! Calculate demand for irrigated water at start of every day:
! hadrd - can be left out if subroutine is called just once per day from
! control.f90
  IF ( MOD(a_step, tspd) == 0 .AND. ( a_step > 1 ) ) THEN
! Loop over land points
    DO l=1, land_pts

! Determine if irrigation should be applied
      l_julian_crop = .FALSE.

      SELECT CASE ( irr_crop )
        CASE ( 1 )
          ! original crop model based on Doell & Siebert, 2002
          ! khalla- If today is before harvest date (plant_n + 150) and after
          ! plant_n
          IF ( (day_of_yr - plant_n(l) < nday_crop) .AND.                     &
               (day_of_yr - plant_n(l) >= 0) )                                &
            l_julian_crop = .TRUE.

          ! khalla- If plant_n is less than 150 days before end of year and
          !         today is early in year, before harvest date
          IF ( (plant_n(l) > days_in_yr - nday_crop) .AND.                    &
               (day_of_yr <= plant_n(l) + nday_crop - days_in_yr) )           &
            l_julian_crop = .TRUE.

          IF ( plant_n(l) <= 0 ) l_julian_crop = .FALSE.
                                                    ! crop code not yet called

        CASE ( 2 )
          ! use JULES-crop development index
          ! for the time being, the maximum dvi across all crop tiles is used to
          ! trigger irrigation. However, it may be better to vary the irrigated
          ! fraction according to which tiles have a suitable dvi
            IF ( dvimax(l) > -1.0 .AND. dvimax(l) < 2.0 )                     &
              l_julian_crop = .TRUE.

        CASE DEFAULT
          CALL ereport('irrig_dmd', 101, 'Invalid value for irr_crop')

      END SELECT


! Irrigate if crop is planted and irrigation fraction > 0
      IF ( l_julian_crop .AND. frac_irr(l) > 0.0 ) THEN

! Irrigate top two soil layers
        DO n=1,2       !SM_LEVELS
! Only irrigate if there is no frozen soil:
          IF ( sthf(l,n) <= 0.0 ) THEN
            IF ( sthu_irr(l,n) < (smvccl(l,n) / smvcst(l,n)) ) THEN
              sthu_tmp = sthu_irr(l,n)
              sthu_irr(l,n) = smvccl(l,n) / smvcst(l,n)

! Ensure that irrigated soil moisture is less than saturation:
              sthu_irr(l,n) = MIN( sthu_irr(l,n), 1.0 - sthf(l,n) )

! Ensure that gridbox mean soil moisture is less than saturation:
              sthu_irr(l,n) = MIN( sthu_irr(l,n), (( 1.0 - sthf(l,n)          &
                            - sthu(l,n)) / frac_irr(l) + sthu_tmp) )
              irrigwater_levels_day(l,n) = (sthu_irr(l,n) - sthu_tmp)         &
                              * rho_water * dzsoil(n)                         &
                              * smvcst(l,n) * frac_irr(l)
              sthu(l,n) = sthu(l,n) + frac_irr(l) * (sthu_irr(l,n) - sthu_tmp)

              IF ( sthu(l,n) + sthf(l,n) > 1.0 ) THEN
                CALL ereport('irrig_dmd', 101, 'GDM super saturation')
              END IF

! Update gridbox moisture content
              smcl(l,n) = (sthu(l,n) + sthf(l,n))                             &
                        * rho_water * dzsoil(n) * smvcst(l,n)

            END IF   ! ( STHU_IRR(L,N) <  (SMVCCL(L,N) / SMVCST(L,N)) )
          END IF   ! ( STHF(L,N) <= 0.0 )

        END DO ! SM_LEVELS

!-------------------------------------------------------------------------------
! Constrain irrigation supply if required
!-------------------------------------------------------------------------------

        IF ( l_irrig_limit ) THEN

            irrig_sum = 0.0
            irrig_max = 0.0
            irrig_lim = 0.0
            irrig_dif = 0.0

            ! total irrigation water added
            ! note that irrigwater_levels_day is already multiplied by frac_irr
            ! i.e. units are kg/m2 for entire grid box
             DO n=1,sm_levels
               irrig_sum = irrig_sum + irrigwater_levels_day(l,n)
             END DO

             IF ( irrig_sum > 0.0 ) THEN
              ! amount of water added should not be more than available
              ! from deep groundwater store
              IF ( .NOT. irr_zw_all ) THEN
                  ! option (1) withdraw water until the wilting point
                irrig_max = MAX( smclzw(l) - smclwiltzw(l), 0.0)
              ELSE
                  ! option (2) withdraw all available water
                irrig_max = MAX( smclzw(l), 0.0)
              END IF
              irrig_lim = MIN(irrig_max, irrig_sum)
              smclzw_rest(l) = smclzw(l) - irrig_lim !min(irrig_max,irrig_sum)

              ! re-calculate soil moisture fraction in deep layer
              sthzw(l) = MAX( smclzw_rest(l) / smclsatzw(l), 0.0 )

              ! if irrigation is constrained, try extracting from river
              ! routestore
              IF ( rivers_sto_per_m2_on_landpts(l) > 0.0) THEN
                irrig_riv(l) = MAX(irrig_sum - irrig_lim, 0.0)
                irrig_riv(l) = MIN(rivers_sto_per_m2_on_landpts(l),irrig_riv(l))
                rivers_adj_on_landpts(l) =                                    &
                  (rivers_sto_per_m2_on_landpts(l) - irrig_riv(l)) /          &
                  rivers_sto_per_m2_on_landpts(l)
                irrig_lim = irrig_lim + irrig_riv(l)
              END IF ! rivers_sto_per_m2_on_landpts gt 0


              ! re-calculate soil moisture
              ! if irrigation is constrained, some or all of the extra water
              ! that was added should be subtracted again
              DO n=1,sm_levels
                ! amount of water to be subtracted
                irrig_dif = ( 1.0 - irrig_lim/irrig_sum ) *                   &
                              irrigwater_levels_day(l,n)
                ! update soil moisture variables
                irrigwater_levels_day(l,n) =                                  &
                  MAX( irrigwater_levels_day(l,n) - irrig_dif, 0.0)
                sthu_irr(l,n) = sthu_irr(l,n) - irrig_dif /                   &
                                (rho_water *dzsoil(n) *                       &
                                smvcst(l,n) * frac_irr(l))
                smcl(l,n) = smcl(l,n) - irrig_dif
                sthu(l,n) = smcl(l,n) / (rho_water *dzsoil(n) * smvcst(l,n))  &
                            - sthf(l,n)
              END DO ! n
            END IF ! irrigsum gt 0
        END IF ! l_irrig_limit

      END IF   ! ( L_JULIAN_CROP .AND. FRAC_IRR(L) >  0.0 )

      irrig_water_gb(l) = SUM(irrigwater_levels_day(l,:)) / FLOAT(secs_in_day)
      IF ( irrig_water_gb(l) > EPSILON(0.0) ) l_irrigated(l) = .TRUE.

    END DO ! LAND_PTS

  END IF ! (MOD(A_STEP,TSPD)

  ! number of days on which irrigation is applied
  WHERE ( l_irrigated ) irrDaysDiag_gb = irrDaysDiag_gb+1

  RETURN

END SUBROUTINE irrig_dmd
