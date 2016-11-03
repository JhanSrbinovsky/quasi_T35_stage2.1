#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE update_derived_variables()

  USE datetime_mod, ONLY : SECS_IN_DAY

  USE model_time_mod, ONLY : current_time, timestep_len

  USE jules_surface_mod, ONLY : l_point_data, l_aggregate

  USE jules_vegetation_mod, ONLY : can_model, l_irrig_dmd

  USE jules_soil_mod, ONLY : sm_levels

  USE ancil_info, ONLY : land_pts, land_index, nsurft, surft_pts, surft_index,&
                         frac_surft

  USE theta_field_sizes, ONLY : t_i_length

  USE forcing, ONLY:                                                          &
    pstar_ij, u_1_ij, v_1_ij, u_0_ij, v_0_ij, con_rain_ij, con_snow_ij,       &
    ls_rain_ij, ls_snow_ij, sw_down_ij, lw_down_ij, qw_1_ij, tl_1_ij

  USE u_v_grid, ONLY: u_0_p_ij, v_0_p_ij, u_1_p_ij, v_1_p_ij


  USE imogen_drive_vars

  USE prognostics, ONLY : canht_pft, lai_pft

  USE p_s_parms, ONLY : satcon_gb, catch_snow_surft, catch_surft, infil_surft,&
                        z0_surft, z0h_bare_surft, z0m_soil_gb, sthu_gb

  !USE crop_vars_mod, ONLY : readFileIrr, frac_irr_gb, frac_irr_surft, sthu_irr, dvi

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Updates variables that are derived from those given in time-varying files
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
  INTEGER :: insd  ! Timestep in day - used to index IMOGEN arrays

  INTEGER :: i,j,l  ! Index variables


!------------------------------------------------------------------------------


  IF( l_imogen ) THEN
!-------------------------------------------------------------------------------
! If IMOGEN is enabled, copy the correct timestep of the current climatology
! into the driving variables
!-------------------------------------------------------------------------------
! Get the timestep in the day that we are on
    insd = (current_time%time / timestep_len) + 1

    DO l = 1,land_pts
      j = (land_index(l) - 1) / t_i_length + 1
      i = land_index(l) - (j - 1) * t_i_length

      pstar_ij(i,j)    = PSTAR_OUT(l,current_time%month,current_time%day,insd)
      u_1_ij(i,j)      = WIND_OUT(l,current_time%month,current_time%day,insd)
      v_1_ij(i,j)      = 0.0
      u_0_ij(i,j)      = 0.0
      v_0_ij(i,j)      = 0.0
      con_rain_ij(i,j) =                                          &
        CONV_RAIN_OUT(l,current_time%month,current_time%day,insd) &
        / REAL(SECS_IN_DAY)
      con_snow_ij(i,j) =                                          &
        CONV_SNOW_OUT(l,current_time%month,current_time%day,insd) &
        / REAL(SECS_IN_DAY)
      ls_rain_ij(i,j)  =                                          &
        LS_RAIN_OUT(l,current_time%month,current_time%day,insd)   &
        / REAL(SECS_IN_DAY)
      ls_snow_ij(i,j)  =                                          &
        LS_SNOW_OUT(l,current_time%month,current_time%day,insd)   &
        / REAL(SECS_IN_DAY)
      sw_down_ij(i,j)  = SW_OUT(l,current_time%month,current_time%day,insd)
      lw_down_ij(i,j)  = LW_OUT(l,current_time%month,current_time%day,insd)
      qw_1_ij(i,j)     = QHUM_OUT(l,current_time%month,current_time%day,insd)
      tl_1_ij(i,j)     = T_OUT(l,current_time%month,current_time%day,insd)
    ENDDO

  ELSE
!-----------------------------------------------------------------------------
! Otherwise, update the main driving variables based on what was given in the
! input file(s) - the update of those variables has already taken place
!-----------------------------------------------------------------------------
! Initialise wind based on io_wind_speed
! This is just a case of initialising variables that are not set via file
! u_1 is always set from file, v_1 is only set from file if using both
! wind components
    IF ( io_wind_speed ) v_1_ij(:,:) = 0.0
    u_0_ij(:,:) = 0.0
    v_0_ij(:,:) = 0.0

    IF (l_daily_disagg) THEN
!-----------------------------------------------------------------------------
! If we are using the disaggregator, we need to disaggregate the precip at
! the start of every day, then impose a diurnal cycle
!-----------------------------------------------------------------------------
      IF ( current_time%time == 0 ) THEN
          CALL update_precip_variables()
          CALL fill_disaggregated_precip_arrays()
      END IF

      CALL impose_diurnal_cycle()
    ELSE
!-----------------------------------------------------------------------------
! Otherwise, just update the precip variables based on the io_precip_type
!-----------------------------------------------------------------------------
      CALL update_precip_variables()
    END IF

!-----------------------------------------------------------------------------
! Radiation variables are updated in CONTROL, since the new albedos need to be
! known before the update takes place
!-----------------------------------------------------------------------------
  END IF

!-------------------------------------------------------------------------------
! Copy information to U, V and T grids (assume that att grids are the same)
!-------------------------------------------------------------------------------
  u_0_p_ij(:,:) = u_0_ij(:,:)
  v_0_p_ij(:,:) = v_0_ij(:,:)
  u_1_p_ij(:,:) = u_1_ij(:,:)
  v_1_p_ij(:,:) = v_1_ij(:,:)


!-----------------------------------------------------------------------------
! Update any variables dependent on variables being prescribed that are not
! driving variables
!-----------------------------------------------------------------------------
  IF ( have_prescribed_veg ) THEN
    WHERE ( lai_pft(:,:) <= 0.01 )
      lai_pft = 0.01
    END WHERE

    CALL sparm(land_pts, nsurft, can_model, l_aggregate, surft_pts,           &
               surft_index, frac_surft, canht_pft, lai_pft, satcon_gb, z0m_soil_gb,       &
               catch_snow_surft, catch_surft, infil_surft, z0_surft, z0h_bare_surft)
  END IF

!-----------------------------------------------------------------------------
! Update irrigation fractions
!-----------------------------------------------------------------------------
  IF ( l_irrig_dmd ) THEN
    CALL assign_irrig_fraction()
    CALL update_irrig_variables()
  END IF

  RETURN

END SUBROUTINE update_derived_variables

#endif
