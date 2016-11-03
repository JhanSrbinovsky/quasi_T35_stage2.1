! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

 SUBROUTINE crop(p_field, land_pts, land_index, a_step, crop_call, sm_levels, &
                 frac, phot, dphotdt, t_surft, t_soil, sthu, satcon, smvccl,  &
                 smvcst, npp_ft_acc, canht, lai, dvi, rootc, harvc, reservec, &
                 croplai, cropcanht, catch_t, infil_t, z0_t)

  USE cropparm, ONLY : t_mort, yield_frac
  USE crop_vars_mod, ONLY : sow_date_cpft, tt_veg_cpft, tt_rep_cpft,          &
                            yield_diag_cpft, stemc_diag_cpft,                 &
                            leafc_diag_cpft, nonyield_diag_cpft,              &
                            harvest_trigger_cpft, harvest_counter_cpft

  USE jules_surface_types_mod, ONLY : ncpft, ntype, nnpft, npft

  USE crop_utils_mod, ONLY : reset_crop, cropharvc_min,                       &
                             leafc_from_prognostics, stemc_from_prognostics,  &
                             lai_from_leafc, canht_from_stemc

  USE datetime_utils_mod, ONLY : day_of_year, days_in_year

  USE time_info_mod, ONLY : l_360, l_leap, current_model_time

  USE jules_vegetation_mod, ONLY : l_prescsow

  USE jules_print_mgr, ONLY : jules_print

  USE ancil_info, ONLY : nsurft

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Calls crop routines and harvests crop when DVI exceeds 2
!
! Method:
!   See JULES-crop technical documentation (Tom Osborne, Josh Hooker).
!
! Current Code Owner: Tom Osborne
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

  INTEGER, INTENT(IN) ::                                                      &
    p_field,                                                                  &
      ! Number of model points
    land_pts,                                                                 &
      ! Number of land points to be processed.
    land_index(land_pts),                                                     &
      ! Index of land points
    a_step,                                                                   &
      ! Atmospheric timestep number
    crop_call,                                                                &
      ! Indicates when to call sow and partition (daily timestep)
    sm_levels
      ! Number of soil levels

  REAL, INTENT(IN) ::                                                         &
    frac(land_pts,ntype),                                                     &
      ! Fractional coverage of tiles
    phot(p_field),                                                            &
      ! Photoperiod (hours)
    dphotdt(p_field),                                                         &
      ! Change in photoperiod in hours between this day and the previous day
    t_surft(land_pts,npft),                                                   &
      ! Temperature of tile
    t_soil(land_pts,sm_levels),                                               &
      ! Soil temperatures
    sthu(land_pts,sm_levels),                                                 &
      ! Unfrozen soil moisture content of each layer as a fraction of saturation
    satcon(land_pts),                                                         &
      ! Saturated hydraulic conductivity (kg/m2/s)
    smvccl(land_pts,sm_levels),                                               &
      ! Critical volumetric SMC (cubic m per cubic m of soil)
    smvcst(land_pts,sm_levels)
      ! Volumetric saturation point (m3/m3 of soil)

  REAL, INTENT(INOUT) ::                                                      &
    npp_ft_acc(land_pts,npft),                                                &
      ! Accumulated npp
    canht(land_pts,npft),                                                     &
      ! Canopy height of tiles
    lai(land_pts,npft),                                                       &
      ! Lai of tiles
    dvi(land_pts,ncpft),                                                      &
      ! Development index of crop tiles
    rootc(land_pts,ncpft),                                                    &
      ! Root carbon of crop tiles
    harvc(land_pts,ncpft),                                                    &
      ! Harvest carbon of crop tiles
    reservec(land_pts,ncpft)
      ! Carbon in stem reserve pool

  REAL, INTENT(OUT) ::                                                        &
    croplai(land_pts,ncpft),                                                  &
      ! Lai of crop tiles
    cropcanht(land_pts,ncpft),                                                &
      ! Canopy height of crop tiles
    catch_t(land_pts,nsurft),                                                 &
      ! Canopy capacity (kg/m2).
    infil_t(land_pts,nsurft),                                                 &
      ! Maximum surface infiltration rate (kg/m2/s).
    z0_t(land_pts,nsurft)
      ! Roughness length (m)

! Local variables
  INTEGER :: n                            ! crop tile number
  INTEGER :: i,j,l                        ! counters
  INTEGER :: surft_pts(ntype)              ! number of land points which
                                          ! include the nth surface type
  INTEGER :: surft_index(land_pts,ntype)   ! indices of land points which
                                          ! include the nth surface type

  INTEGER ::  day        ! Current day
  INTEGER ::  month      ! Current month
  INTEGER ::  year       ! Current year
  INTEGER ::  day_of_yr  ! Current day of year
  INTEGER ::  days_in_yr ! Total number of days in current year
  INTEGER ::  day_before_sowdate ! Day of year for day just before sowing date

  REAL :: nonyield_diag_dummy             ! dummy variable used when making
       ! sure that crop is consistent at initialisation if it's not sown

!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Create the surft_index array of land points with each surface type
!-----------------------------------------------------------------------------
! NB: tilepts loops over all land surface types
  CALL tilepts(land_pts,frac,surft_pts,surft_index)

!-----------------------------------------------------------------------------
! Initialise two of the crop diagnostics at start of run
!-----------------------------------------------------------------------------
  IF ( a_step == 1 ) THEN
    DO  n=1,ncpft
      DO j=1,surft_pts(nnpft+n)
        l=surft_index(j,nnpft+n)
        IF ( dvi(l,n) < 0.0 ) THEN
          CALL reset_crop(n, nonyield_diag_dummy, dvi(l,n),                   &
                          lai(l,n+nnpft), canht(l,n+nnpft),                   &
                          rootc(l,n), harvc(l,n), reservec(l,n),              &
                          stemc_diag_cpft(l,n), leafc_diag_cpft(l,n))
        ELSE
          !------------------------------------------------
          ! Derive stem and leaf pools from prognostics
          !------------------------------------------------
          stemc_diag_cpft(l,n) = stemc_from_prognostics(n, canht(l,n+nnpft))
          leafc_diag_cpft(l,n) = leafc_from_prognostics(n, dvi(l,n),          &
                                                          lai(l,n+nnpft))
        END IF
      END DO
    END DO
  END IF

  IF ( l_prescsow ) THEN
    ! Get the current model year, month and day
    CALL current_model_time(year, month, day)

    day_of_yr  = day_of_year(year, month, day, l_360, l_leap)
    days_in_yr = days_in_year(year, l_360, l_leap)
  END IF 

!-----------------------------------------------------------------------------
! Loop over crop tiles and sow/grow/harvest the crop
!-----------------------------------------------------------------------------
  DO  n=1,ncpft

!-----------------------------------------------------------------------------
! Loop over crop tile points
!-----------------------------------------------------------------------------
    DO j=1,surft_pts(nnpft+n)

      l = surft_index(j,nnpft+n)
      i = land_index(l)

      !reset yield and non-yield carbon
      yield_diag_cpft(l,n) = 0.0
      harvest_trigger_cpft(l,n) = 0
      harvest_counter_cpft(l,n) = 0
      nonyield_diag_cpft(l,n) = 0.0

      IF ( nint(dvi(l,n)) == -2 .AND. crop_call == 0 ) THEN
          CALL sow(n, sm_levels, t_soil(l,:), sthu(l,:), smvcst(l,:),         &
                   smvccl(l,:), dphotdt(i), sow_date_cpft(l,n), dvi(l,n))
      END IF

      IF ( dvi(l,n) >= -1.0 .AND. dvi(l,n) < 0.0 ) THEN
          CALL emerge(n, t_surft(l,n+nnpft), dvi(l,n),                        &
                      rootc(l,n), harvc(l,n),                                 &
                      stemc_diag_cpft(l,n), leafc_diag_cpft(l,n))
      END IF

      IF ( dvi(l,n) >= 0.0 ) THEN

        CALL develop(n, t_surft(l,n+nnpft), phot(i), tt_veg_cpft(l,n),        &
        tt_rep_cpft(l,n), dvi(l,n))

        IF (crop_call == 0) THEN
          CALL partition(n, npp_ft_acc(l,n+nnpft),                            &
                         dvi(l,n), rootc(l,n), harvc(l,n),                    &
                         reservec(l,n), nonyield_diag_cpft(l,n),              &
                         stemc_diag_cpft(l,n), leafc_diag_cpft(l,n),          &
                         harvest_counter_cpft(l,n), harvest_trigger_cpft(l,n))
          npp_ft_acc(l,n+nnpft) = 0.0
        END IF

        !------------------------------------------------
        ! Derive new prognostics from carbon pools and DVI
        !------------------------------------------------
        lai(l,n+nnpft) = lai_from_leafc(n, dvi(l,n), leafc_diag_cpft(l,n))
        canht(l,n+nnpft) = canht_from_stemc(n, stemc_diag_cpft(l,n))

        IF ( dvi(l,n) >= 2.0 ) THEN
          harvest_counter_cpft(l,n) = 1
          harvest_trigger_cpft(l,n) = 1
        ELSE IF( lai(l,n+nnpft) >= 15.0 ) THEN
          harvest_counter_cpft(l,n) = 1
          harvest_trigger_cpft(l,n) = 2
        ELSE IF ( dvi(l,n) > 1.0 .AND. t_soil(l,2) < t_mort(n) )  THEN
          harvest_counter_cpft(l,n) = 1
          harvest_trigger_cpft(l,n) = 3
        ! harvest_trigger_cpft = 4 is set in partition.F90
        ELSE IF ( l_prescsow ) THEN
            day_before_sowdate = NINT(sow_date_cpft(l,n)) - 1
            IF ( day_before_sowdate == 0 ) day_before_sowdate = days_in_yr
          IF (day_of_yr == day_before_sowdate) THEN
            harvest_counter_cpft(l,n) = 1
            harvest_trigger_cpft(l,n) = 5
          END IF
        END IF        

        IF ( harvest_counter_cpft(l,n) == 1 ) THEN
          yield_diag_cpft(l,n)  = yield_frac(n) * harvc(l,n) - cropharvc_min
          harvc(l,n) = harvc(l,n) - yield_diag_cpft(l,n)
          dvi(l,n) = -2.0
        END IF

        IF ( nint(dvi(l,n)) == -2 ) THEN
          CALL reset_crop(n, nonyield_diag_cpft(l,n), dvi(l,n),               &
                            lai(l,n+nnpft), canht(l,n+nnpft),                 &
                            rootc(l,n), harvc(l,n), reservec(l,n),            &
                            stemc_diag_cpft(l,n), leafc_diag_cpft(l,n))
        END IF

      END IF  !DVI gt 0.0

      croplai(l,n)   = lai(l,n+nnpft)
      cropcanht(l,n) = canht(l,n+nnpft)

      IF (crop_call == 0) THEN
        nonyield_diag_cpft(l,n) = nonyield_diag_cpft(l,n) +                   &
                                    npp_ft_acc(l,n+nnpft)
        npp_ft_acc(l,n+nnpft) = 0.0
      END IF

    END DO  ! surft_pts

!--------------------------------------------
! Update canopy capacity and roughness length
!--------------------------------------------
    CALL pft_sparm(land_pts, n+nnpft, surft_index(:,n+nnpft),                 &
                   surft_pts(n+nnpft), canht(:,n+nnpft), lai(:,n+nnpft),      &
                   satcon(:), catch_t(:,n+nnpft), infil_t(:,n+nnpft),         &
                   z0_t(:,n+nnpft))

  END DO

END SUBROUTINE crop
