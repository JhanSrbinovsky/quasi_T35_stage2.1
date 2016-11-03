#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE topmodel_init()

  USE ancil_info, ONLY : land_pts, soil_pts, soil_index

  USE p_s_parms, ONLY : bexp_gb, satcon_gb, sathh_gb, sthf_gb, sthu_gb,       &
                        smvcst_gb

  USE jules_soil_mod, ONLY : dzsoil, sm_levels

  USE top_pdm, ONLY : fexp_gb, gamtot_gb, ti_mean_gb, ti_sig_gb, zw_gb,       &
                      fsat_gb, fwetl_gb, qbase_gb, sthzw_gb

  USE jules_hydrology_mod, ONLY : l_wetland_unfrozen

  USE c_densty, ONLY : rho_water  !  density of pure water (kg/m3)

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Finish initialising TOPMODEL by calculating surface saturated fraction
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
  LOGICAL, PARAMETER :: l_gamtot = .FALSE.  !  Switch for calculation of gamtot.
                                            !  .FALSE. so that calc_fsat
                                            !  calculates fsat_gb.

  INTEGER :: i, j, n                    ! Index variables

  REAL :: qbase_l(land_pts,sm_levels+1) ! Base flow from each layer (kg/m2/s).
  REAL :: top_crit(land_pts)            ! Critical topographic index required
                                        !   to calculate the surface saturation
                                        !   fraction.
  REAL :: zdepth(0:sm_levels)           ! Lower soil layer boundary depth (m).
  REAL :: wutot(land_pts)               ! UNFROZEN to TOTAL fraction at ZW.
  REAL :: dumwutot(land_pts)            ! Dummy UNFROZEN to TOTAL fraction at ZW
                                        !(always set to 1).
  REAL :: ksz(land_pts,0:sm_levels)     ! Saturated hydraulic conductivity for
                                        !   each layer (kg/m2/s).
  REAL :: dumsthf(land_pts,sm_levels)   ! Dummy Frozen soil moisture content of
                                        !   each layer as a fraction of
                                        !   saturation (always set to 0).
  REAL :: smcl(land_pts,sm_levels)      ! Soil moisture content
                                        !   of each layer (kg/m2).
  REAL :: smclsat(land_pts,sm_levels)   ! Saturated soil moisture content
                                        !   of each layer (kg/m2).
  REAL :: zw_inund(land_pts)            ! Adjusted Water table depth (m).

!-------------------------------------------------------------------------------
  zdepth(0) = 0.0
  DO n = 1,sm_levels
    zdepth(n) = zdepth(n - 1) + dzsoil(n)
  END DO

! Set values that are retained at non-soil points.
  fsat_gb(:) = 0.0
  fwetl_gb(:) = 0.0
  dumwutot(:) = 1.0

  IF ( soil_pts /= 0 ) THEN
    DO j = 1,soil_pts
      i = soil_index(j)
      DO n = 0,sm_levels
        ksz(i,n) = satcon_gb(i,n)
      END DO
      DO n = 1,sm_levels
        smclsat(i,n) = rho_water*dzsoil(n)*smvcst_gb(i,n)
        smcl(i,n) = (sthu_gb(i,n)+sthf_gb(i,n))*smclsat(i,n)
      ENDDO
    END DO

    IF(L_wetland_unfrozen)THEN
      CALL calc_zw_inund(land_pts,sm_levels,soil_pts,soil_index,zdepth        &
                        ,bexp_gb(:,1),sathh_gb(:,1),smclsat,smcl,sthu_gb,     &
                        sthzw_gb,zw_gb,zw_inund,wutot)
      dumsthf(:,:) = 0.0
    ELSE
      zw_inund(:) = zw_gb(:)
      dumsthf(:,:) = sthf_gb(:,:)
    ENDIF

!   We need top_crit - get this from calc_baseflow.
    CALL calc_baseflow_jules(soil_pts, soil_index, land_pts, sm_levels,       &
                             zdepth, ksz, bexp_gb, fexp_gb, ti_mean_gb,       &
                             zw_inund, dumsthf,                               &
                             top_crit, qbase_gb, qbase_l)

!   Call calc_fsat with 1st argument (l_gamtot)=.FALSE. so as to get fsat_gb.
    CALL calc_fsat(l_gamtot, soil_pts, soil_index, land_pts, ti_mean_gb,      &
                   ti_sig_gb,                                                 &
                   dumwutot, top_crit, gamtot_gb, fsat_gb, fwetl_gb)

    IF(L_wetland_unfrozen)THEN
      DO j = 1,soil_pts
        i = soil_index(j)
        fsat_gb(i)=wutot(i)*fsat_gb(i)
        fwetl_gb(i)=wutot(i)*fwetl_gb(i)
      END DO
    ENDIF

  END IF  !  soil_pts

  RETURN

END SUBROUTINE topmodel_init
#endif
