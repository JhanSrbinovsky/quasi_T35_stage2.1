#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE init_parms()

  USE theta_field_sizes, ONLY : t_i_length, t_j_length

  use jules_surface_mod, ONLY : l_aggregate

  USE jules_sea_seaice_mod, ONLY: nice

  USE jules_vegetation_mod, ONLY : can_model

  USE ancil_info, ONLY : ssi_index, fssi_ij, ice_fract_ij,                    &
                         ice_fract_ncat_sicat,                                &
                         sea_frac, sice_frac, sice_frac_ncat, sea_index,      &
                         sice_index, sice_frac, sea_frac, sice_pts_ncat,      &
                         sice_index_ncat, land_pts, nsurft, sea_pts,          &
                         sice_pts, ssi_pts, surft_pts, surft_index,           &
                         frac_surft, land_index

  USE coastal, ONLY: tstar_sice_ij, tstar_land_ij, tstar_ssi_ij, tstar_sea_ij,&
                     flandg, tstar_sice_sicat

  USE prognostics, ONLY : tstar_surft, canht_pft, lai_pft

  USE p_s_parms, ONLY : catch_surft, catch_snow_surft, infil_surft, satcon_gb,&
                        z0_surft, z0h_bare_surft, z0m_soil_gb

  USE fluxes, ONLY : tstar_ij

  USE u_v_grid, ONLY : dtrdz_charney_grid_1_ij

  USE logging_mod, ONLY : log_info, log_debug, log_warn, log_error, log_fatal

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Initialises various variables that may change their initialisation in
!   future versions
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
  INTEGER :: error, error_sum  ! Error indicator

  INTEGER :: i,j,l,n  ! Loop counters


!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Calculate surface parameters.
!-----------------------------------------------------------------------------
  CALL sparm(land_pts, nsurft, can_model, l_aggregate, surft_pts, surft_index,&
             frac_surft, canht_pft, lai_pft, satcon_gb, z0m_soil_gb,          &
             catch_snow_surft, catch_surft, infil_surft, z0_surft,            &
             z0h_bare_surft)

!-----------------------------------------------------------------------
! Set up index for sea and sea-ice
!-----------------------------------------------------------------------
  SSI_PTS = 0
  SSI_INDEX(:) = 0
  DO j = 1,t_j_length
    DO i = 1,t_i_length
      SSI_PTS = SSI_PTS + 1
      IF ( FLANDG(i,j) < 1.0 ) THEN
        SSI_INDEX(SSI_PTS) = (j - 1) * t_i_length + i
      ENDIF
      fssi_ij(i,j)=1.0 - FLANDG(i,j)
    ENDDO
  ENDDO

!-----------------------------------------------------------------------
! Set sea ice fraction.
!-----------------------------------------------------------------------
  DO j = 1,t_j_length
    DO i = 1,t_i_length
      ice_fract_ij(I,J) = 0.0
      tstar_sice_ij(I,J) = 0.0
      DO N=1,NICE
        ice_fract_ij(I,J) = ice_fract_ij(I,J) + ice_fract_ncat_sicat(I,J,N)
      ENDDO
      IF (ice_fract_ij(I,J) > 0.0) THEN
        DO N=1,NICE  !assuming nice=nice_use here
          tstar_sice_ij(I,J) = tstar_sice_ij(I,J)                             &
                    + ice_fract_ncat_sicat(I,J,N) * tstar_sice_sicat(I,J,N) / &
                                                          ice_fract_ij(I,J)
        END DO
      END IF
    ENDDO
  ENDDO

!-----------------------------------------------------------------------
! Initialise sea and sea-ice indices
!-----------------------------------------------------------------------
  SEA_PTS  = 0
  SICE_PTS = 0
  SEA_INDEX(:)  = 0
  SICE_INDEX(:) = 0
  SICE_FRAC(:) = 0.0
  SEA_FRAC(:)  = 0.0
  DO L = 1,SSI_PTS
    J = (SSI_INDEX(L) - 1) / t_i_length + 1
    I = SSI_INDEX(L) - (J - 1) * t_i_length
    IF ( SSI_INDEX(L) > 0 ) THEN
      IF ( ice_fract_ij(I,J) > 0.0 ) THEN
        SICE_PTS = SICE_PTS + 1
        SICE_INDEX(SICE_PTS) = L
        SICE_FRAC(L) = ice_fract_ij(I,J)
      ENDIF
      IF( ice_fract_ij(I,J) < 1.0 ) THEN
        SEA_PTS = SEA_PTS + 1
        SEA_INDEX(SEA_PTS) = L
        SEA_FRAC(L) = 1.0 - SICE_FRAC(L)
      ENDIF
    ENDIF
  ENDDO

  SICE_PTS_NCAT(:) = 0
  SICE_INDEX_NCAT(:,:) = 0
  SICE_FRAC_NCAT(:,:) = 0.0
  DO N = 1,NICE
    DO L = 1,SSI_PTS
      J = (SSI_INDEX(L) - 1) / t_i_length + 1
      I = SSI_INDEX(L) - (J - 1) * t_i_length
      IF ( SSI_INDEX(L) > 0 ) THEN
        IF ( ice_fract_ncat_sicat(I,J,N) > 0.0 ) THEN
          SICE_PTS_NCAT(N) = SICE_PTS_NCAT(N) + 1
          SICE_INDEX_NCAT(SICE_PTS_NCAT(N),N) = L
          SICE_FRAC_NCAT(L,N) = ice_fract_ncat_sicat(I,J,N)
        ENDIF
      ENDIF
    ENDDO
  ENDDO

!-----------------------------------------------------------------------
! Set up gridbox "prognostics".
!-----------------------------------------------------------------------

  tstar_ij(:,:)      = 0.0
  tstar_land_ij(:,:) = 0.0
  tstar_ssi_ij(:,:)  = 0.0

  DO L = 1,LAND_PTS
    J = (LAND_INDEX(L) - 1) / t_i_length + 1
    I = LAND_INDEX(L) - (J - 1) * t_i_length
    IF ( L_AGGREGATE ) THEN
      tstar_land_ij(I,J) = TSTAR_surft(L,1)
    ELSE
      DO N = 1,nsurft
        tstar_land_ij(I,J) = tstar_land_ij(I,J) + frac_surft(L,N) *           &
                                                               TSTAR_surft(L,N)
      ENDDO
    ENDIF
  ENDDO

  tstar_ssi_ij(:,:) = (1.0 - ice_fract_ij(:,:)) * tstar_sea_ij(:,:)           &
                 + ice_fract_ij(:,:) * tstar_sice_ij(:,:)
  tstar_ij(:,:) = flandg(:,:) * tstar_land_ij(:,:)                            &
             + (1.0 - flandg(:,:)) * tstar_ssi_ij(:,:)

!-----------------------------------------------------------------------------
! Set up information on U, V and T grids (assume that att grids are the same)
!-----------------------------------------------------------------------------
  dtrdz_charney_grid_1_ij(:,:) = 0.0

  RETURN

END SUBROUTINE init_parms
#endif
