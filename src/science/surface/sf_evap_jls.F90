! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  SUBROUTINE SF_EVAP------------------------------------------------

!  Purpose: Calculate surface evaporation and sublimation amounts
!           (without applying them to the surface stores).

!  Suitable for single column usage.

!  Documentation: UMDP 24

!--------------------------------------------------------------------
SUBROUTINE sf_evap (                                                          &
 land_pts,nsurft,                                                             &
 land_index,surft_index,surft_pts,nshyd,fland,                                &
 ashtf_prime_surft,canopy,dtrdz_1,flake,fraca,snow_surft,resfs,               &
 resft,rhokh_1,tile_frac,smc,wt_ext_surft,timestep,GAMMA,                     &
 fqw_1,fqw_surft,ftl_1,ftl_surft,tstar_surft,                                 &
 ecan,ecan_surft,elake_surft,esoil,esoil_surft,ei_surft,ext)

USE atm_fields_bounds_mod
USE theta_field_sizes, ONLY : t_i_length

USE c_r_cp
USE c_lheat
USE c_0_dg_c

USE jules_vegetation_mod, ONLY : l_irrig_dmd

USE crop_vars_mod, ONLY: frac_irr_gb, frac_irr_surft, wt_ext_irr_surft        &
                         ,ext_irr_gb, smc_irr_gb, resfs_irr_surft

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

INTEGER                                                                       &
 land_pts                                                                     &
                       ! IN Number of land points to be processed.
,nsurft                                                                       &
                       ! IN Number of tiles per land point.
,land_index(land_pts)                                                         &
                       ! IN Index of land points.
,surft_index(land_pts,nsurft)                                                 &
!                            ! IN Index of tile points.
,surft_pts(nsurft)                                                            &
                       ! IN Number of tile points.
,nshyd                 ! IN Number of soil moisture levels.


REAL                                                                          &
 fland(land_pts)                                                              &
                       ! IN Fraction of gridbox which is land.
,ashtf_prime_surft(land_pts,nsurft)                                           &
!                            ! IN Adjusted SEB coefficient
,canopy(land_pts,nsurft)                                                      &
!                            ! IN Surface/canopy water on land
!                            !    tiles (kg/m2).
,dtrdz_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                 &
!                            ! IN -g.dt/dp for surface layer
,flake(land_pts,nsurft)                                                       &
                       ! IN Lake fraction.
,fraca(land_pts,nsurft)                                                       &
                       ! IN Fraction of surface moisture flux
!                            !    with only aerodynamic resistance
!                            !    for land tiles.
,snow_surft(land_pts,nsurft)                                                  &
!                            ! IN Lying snow amount on tiles (kg/m2).
,resfs(land_pts,nsurft)                                                       &
                       ! IN Combined soil, stomatal and aerodynam.
!                            !    resistance factor for fraction 1-FRACA
!                            !    of land tiles.
,resft(land_pts,nsurft)                                                       &
                       ! IN Total resistance factor
!                            !    FRACA+(1-FRACA)*RESFS.
,rhokh_1(land_pts,nsurft)                                                     &
!                            ! IN Surface exchange coefficients.
,tile_frac(land_pts,nsurft)                                                   &
!                            ! IN Tile fractions.
,smc(land_pts)                                                                &
                       ! IN Available soil moisture (kg/m2).
,wt_ext_surft(land_pts,nshyd,nsurft)                                          &
!                            ! IN Fraction of transpiration
!                            !    extracted from each soil layer
!                            !    by each tile.
,timestep                                                                     &
                       ! IN Timestep in seconds.
,GAMMA                 ! IN implicit weight in level 1

REAL                                                                          &
 fqw_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                       ! INOUT Surface moisture flux (kg/m2/s).
,fqw_surft(land_pts,nsurft)                                                   &
!                            ! INOUT Local FQW_1 for tiles.
,ftl_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                       ! INOUT Surface sensible heat flux (W/m2).
,ftl_surft(land_pts,nsurft)                                                   &
!                            ! INOUT Local FTL_1 for tiles.
,tstar_surft(land_pts,nsurft)
!                            ! INOUT Tile surface temperatures (K).

REAL                                                                          &
 ecan(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                    &
                       ! OUT Gridbox mean evaporation from canopy/
!                            !     surface store (kg per sq m per s).
!                            !     Zero over sea and sea-ice.
,ecan_surft(land_pts,nsurft)                                                  &
!                            ! OUT ECAN for land tiles.
,elake_surft(land_pts,nsurft)                                                 &
!                            ! OUT Lake evaporation.
,esoil(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                       ! OUT Gridbox mean evapotranspiration from
!                            !     soil moisture (kg per sq m per s).
!                            !     Zero over sea and sea-ice.
,esoil_surft(land_pts,nsurft)                                                 &
!                            ! OUT ESOIL for land tiles.
,ei_surft(land_pts,nsurft)                                                    &
!                            ! OUT Sublimation from snow or land-ice
!                            !     (kg per sq m per s).
,ext(land_pts,nshyd)   ! OUT Extraction of water from each
!                            !     soil layer (kg/m2/s).

                        !    of land tiles.

REAL                                                                          &
 esoil_irr(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)               &
!                            ! WORK Gridbox mean evapotranspiration from
!                            !     soil moisture (kg per sq m per s).
,esoil_nir(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)               &
!                            ! WORK Gridbox mean evapotranspiration from
!                            !     soil moisture (kg per sq m per s).
,esoil_irr_surft(land_pts,nsurft)                                             &
!                            ! WORK Evapotranspiration from soil
!                            !     moisture through irrigated
!                            !     fraction of land tiles (kg/m2/s).
,esoil_nir_surft(land_pts,nsurft)                                             &
!                            ! WORK Evapotranspiration from soil
!                            !     moisture through non-irrigated
!                            !     fraction of land tiles (kg/m2/s).
,smc_nir(land_pts)
!                            ! WORK Available soil moisture (kg/m2).
!                            !     fraction (kg/m2/s).


REAL                                                                          &
 dfqw(land_pts)                                                               &
                       ! Increment in GBM moisture flux.
,dftl(land_pts)                                                               &
                       ! Increment in GBM sensible heat flux.
,e_surft_old(land_pts,nsurft)                                                 &
!                            ! Surface moisture flux before adjustment.
,le_surft_old(land_pts,nsurft)
!                            ! Surf latent heat flux before adjustment.

REAL                                                                          &
 diff_lat_htf                                                                 &
                       ! Increment in local latent heat flux.
,diff_sens_htf                                                                &
                       ! Increment in local sensible heat flux.
,dtstar                                                                       &
                       ! Increment in local surface temperature.
,edt                                                                          &
                       ! Moisture flux x timestep
,rhokh1_prime          ! Modified forward time-weighted
                       ! transfer coefficient.

INTEGER                                                                       &
 i,j                                                                          &
             ! Loop counter (horizontal field index).
,k                                                                            &
             ! Loop counter (land, snow or land-ice field index).
,m                                                                            &
             ! Loop counter (soil level index).
,l                                                                            &
             ! Loop counter (land point field index).
,n           ! Loop counter (tile index).

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SF_EVAP'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

DO n=1,nsurft
  DO k=1,surft_pts(n)
    l = surft_index(k,n)
    e_surft_old(l,n) = fqw_surft(l,n)
    IF (snow_surft(l,n)  >   0.) THEN
      le_surft_old(l,n) = (lc + lf)*fqw_surft(l,n)
    ELSE
      le_surft_old(l,n) = lc*fqw_surft(l,n)
    END IF
  END DO
END DO

DO n=1,nsurft
  DO l=1,land_pts
    ecan_surft(l,n) = 0.
    esoil_surft(l,n) = 0.
    IF (l_irrig_dmd) THEN
      esoil_nir_surft(l,n) = 0.0
      esoil_irr_surft(l,n) = 0.0
    END IF
    elake_surft(l,n) = 0.
    ei_surft(l,n) = 0.
  END DO
END DO

!-----------------------------------------------------------------------
! Sublimation from snow-covered land tiles
!-----------------------------------------------------------------------
DO n=1,nsurft
  DO k=1,surft_pts(n)
    l = surft_index(k,n)
    IF (snow_surft(l,n)  >   0.) THEN
      ei_surft(l,n) =  fqw_surft(l,n)
      edt = ei_surft(l,n)*timestep
      IF ( edt  >   snow_surft(l,n) )                                         &
        ei_surft(l,n) = snow_surft(l,n) / timestep
        fqw_surft(l,n) = fqw_surft(l,n) -  ei_surft(l,n)
    END IF
  END DO
END DO

!-----------------------------------------------------------------------
! Surface evaporation from and condensation onto snow-free land
!-----------------------------------------------------------------------
DO j=tdims%j_start,tdims%j_end
  DO i=tdims%i_start,tdims%i_end
    ecan(i,j) = 0.
    esoil(i,j) = 0.
    IF (l_irrig_dmd) THEN
      esoil_nir(i,j) = 0.0
      esoil_irr(i,j) = 0.0
    END IF
  END DO
END DO

DO n=1,nsurft
  DO k=1,surft_pts(n)
    l = surft_index(k,n)
    j=(land_index(l)-1)/t_i_length + 1
    i = land_index(l) - (j-1)*t_i_length
    IF ( fqw_surft(l,n)  >   0.0 ) THEN
      ecan_surft(l,n) = (1. - flake(l,n)) *                                   &
                       fraca(l,n) * fqw_surft(l,n) / resft(l,n)
      esoil_surft(l,n) = (1. - flake(l,n)) *                                  &
                        (1. - fraca(l,n))*resfs(l,n)*fqw_surft(l,n)           &
                                                      / resft(l,n)
      elake_surft(l,n) = flake(l,n)*fqw_surft(l,n) / resft(l,n)
      edt = ecan_surft(l,n)*timestep
      IF ( edt  >   canopy(l,n) ) THEN
        esoil_surft(l,n) =  (1. - flake(l,n)) *                               &
                           (1. - fraca(l,n)*canopy(l,n)/edt) *                &
                               resfs(l,n)*fqw_surft(l,n)/resft(l,n)
        ecan_surft(l,n) = canopy(l,n) / timestep
      END IF
    ELSE IF (snow_surft(l,n) <= 0.) THEN
      IF (tstar_surft(l,n) >= tm) THEN
        ecan_surft(l,n) = (1. - flake(l,n))*fqw_surft(l,n)
        elake_surft(l,n) = flake(l,n)*fqw_surft(l,n)
      ELSE
        ei_surft(l,n) =  fqw_surft(l,n)
      END IF
    END IF
    ecan(i,j) = ecan(i,j) + tile_frac(l,n)*ecan_surft(l,n)
    esoil(i,j) = esoil(i,j) + tile_frac(l,n)*esoil_surft(l,n)

    IF ( l_irrig_dmd ) THEN

! split esoil_surft into irrigated and non-irrigated fraction
! for irrigated fraction, esoil_surft scaled according to resfs_irr_surft/resfs
      IF ( resfs(l,n) > 0.0 ) THEN
        esoil_irr_surft(l,n) = esoil_surft(l,n) *                             &
                              resfs_irr_surft(l,n) / resfs(l,n)
      ELSE
        esoil_irr_surft(l,n) = esoil_surft(l,n)
      END IF

! for non-irrigated fraction, assume total esoil_surft is
! linear combination of irrigated and non-irrigated fraction
      IF ( frac_irr_surft(l,n) <= EPSILON(1.0) ) THEN
! but when irrigated fraction is (close to) 0
        esoil_nir_surft(l,n) = esoil_surft(l,n)
      ELSEIF ( frac_irr_surft(l,n) < 1.0 ) THEN
        esoil_nir_surft(l,n) = ( esoil_surft(l,n) -                           &
                                frac_irr_surft(l,n) * esoil_irr_surft(l,n) ) /&
                              (1.0 - frac_irr_surft(l,n) )
      ELSE
        esoil_nir_surft(l,n) = 0.0
      END IF

! hadrd - because ESOIL_IRR_surft is scaled according to resfs_irr_surft/resfs,
! the contribution of the irrigated fraction *in this tile* may (in some cases)
! exceed the total tile ESOIL_surft, and ESOIL_NIR_surft may become negative.
! Including the line below prevents this, and gives an ESOIL_NIR_surft that
! (at the 1-2 grid boxes where this has been checked) is closer to
! ESOIL_surft when irrigation is switched off. However, it may be worth
! pointing out that ESOIL_NIR_surft can still turn out to be less than
! ESOIL_surft without irrigation.
      esoil_nir_surft(l,n) = MAX( esoil_nir_surft(l,n), 0.0 )

! hadrd - (TILE_FRAC*FRAC_IRR_surft) gives value for total grid box, not the
! irrigated fraction only to convert to irrigated fraction of total gridbox,
! multiply by (TILE_FRAC*FRAC_IRR_surft)/FRAC_IRR
      IF ( frac_irr_gb(l) > EPSILON(1.0) ) THEN
        esoil_irr(i,j) = esoil_irr(i,j) +                                     &
                       tile_frac(l,n) * frac_irr_surft(l,n) *                 &
                       esoil_irr_surft(l,n) / frac_irr_gb(l)
      ELSE
        esoil_irr(i,j) = 0.0
      END IF

      IF ( frac_irr_gb(l) > EPSILON(1.0) ) THEN

! hadrd - same for non-irrigated fraction
        IF ( 1.0 - frac_irr_gb(l) > EPSILON(1.0) ) THEN
          esoil_nir(i,j) = esoil_nir(i,j) +                                   &
                        tile_frac(l,n) * (1.0 - frac_irr_surft(l,n)) *        &
                        esoil_nir_surft(l,n)                                  &
                        / (1.0 - frac_irr_gb(l))
        ELSE
          esoil_nir(i,j) = 0.0
        END IF
      ELSE
        esoil_nir(i,j) = esoil(i,j)
      END IF

! total smc is linear combination of smc in irrig and non-irrig fraction
! hadrd - this could/should be done outside loop over tiles
      IF ( smc_irr_gb(l) < 0.0 ) smc_irr_gb(l) = 0.0
        smc_nir(l) = smc_irr_gb(l)
      IF ( frac_irr_gb(l) < 1.0 ) THEN
        smc_nir(l) = ( smc(l) - frac_irr_gb(l) * smc_irr_gb(l) ) /            &
                    (1.0 - frac_irr_gb(l))
! hadrd - added case where total grid box frac_irr_gb is 1.0
      ELSE
        smc_nir(l) = 0.0
      END IF

      IF ( frac_irr_gb(l) <= EPSILON(1.0) ) THEN
        smc_nir(l) = smc(l)
      END IF

      IF ( smc_nir(l) < 0.0 ) smc_nir(l) = 0.0 ! this may not be necessary?

    END IF ! l_irrig_dmd

  END DO
END DO

!-----------------------------------------------------------------------
! Soil evapotranspiration
!-----------------------------------------------------------------------
DO l=1,land_pts
  j=(land_index(l)-1)/t_i_length + 1
  i = land_index(l) - (j-1)*t_i_length


  IF ( l_irrig_dmd ) THEN

! non-irrigated fraction
    edt = esoil_nir(i,j) * timestep
    IF ( edt > smc_nir(l) ) THEN
      DO n=1,nsurft
        esoil_nir_surft(l,n) = smc_nir(l) * esoil_nir_surft(l,n) / edt
      END DO
      esoil_nir(i,j) = smc_nir(l) / timestep
    END IF

! irrigated fraction
    edt = esoil_irr(i,j) * timestep
    IF ( edt > smc_irr_gb(l) ) THEN
      DO n=1,nsurft
        esoil_irr_surft(l,n) = smc_irr_gb(l) * esoil_irr_surft(l,n) / edt
      END DO
      esoil_irr(i,j) = smc_irr_gb(l) / timestep
    END IF

! combine irrigated and non-irrigated fractions
    DO n=1,nsurft
      esoil_surft(l,n) = esoil_nir_surft(l,n)
      IF ( frac_irr_surft(l,n) > EPSILON(1.0) ) THEN
        esoil_surft(l,n) = frac_irr_surft(l,n) * esoil_irr_surft(l,n) +       &
                          (1.0 - frac_irr_surft(l,n) ) * esoil_nir_surft(l,n)
      END IF

      esoil(i,j) = esoil_nir(i,j)
      IF ( frac_irr_gb(l) > EPSILON(1.0) ) THEN
        esoil(i,j) = frac_irr_gb(l) * esoil_irr(i,j) +                        &
                     (1.0 - frac_irr_gb(l) ) * esoil_nir(i,j)
      END IF

    END DO

   ELSE ! l_irrig_dmd is F

      ! revert to original code
     edt = esoil(i,j)*timestep
     IF ( edt  >   smc(l) ) THEN
       DO n=1,nsurft
         esoil_surft(l,n) = smc(l)*esoil_surft(l,n) / edt
       END DO
       esoil(i,j) = smc(l) / timestep
     END IF

  END IF ! l_irrig_dmd


END DO

!-----------------------------------------------------------------------
! Extraction of water from each layer
!-----------------------------------------------------------------------

DO m=1,nshyd
  DO l=1,land_pts
    ext(l,m) = 0.0
    ext_irr_gb(l,m) = 0.0
  END DO
END DO

DO m=1,nshyd
  DO n=1,nsurft
    DO k=1,surft_pts(n)
      l = surft_index(k,n)
      ext(l,m) = ext(l,m) + tile_frac(l,n)*wt_ext_surft(l,m,n)                &
                                          *esoil_surft(l,n)

      IF (l_irrig_dmd) THEN
! irrigated fraction
        IF ( frac_irr_gb(l) > EPSILON(1.0) ) THEN
          ext_irr_gb(l,m) = ext_irr_gb(l,m) + tile_frac(l,n) *                &
                            wt_ext_irr_surft(l,m,n)                           &
                            * esoil_irr_surft(l,n) *                          &
                            frac_irr_surft(l,n)/frac_irr_gb(l)
! hadrd - with original code, FRAC_IRR_surft(L,N)/FRAC_IRR(L) is always 1.0
! since frac_irr_gb is the same in all tiles
        END IF
      END IF
    END DO
  END DO
END DO

!-----------------------------------------------------------------------
! Calculate increments to surface heat fluxes, moisture fluxes and
! temperatures
!-----------------------------------------------------------------------
DO l=1,land_pts
  dftl(l) = 0.
  dfqw(l) = 0.
END DO

DO n=1,nsurft
  DO k=1,surft_pts(n)
    l = surft_index(k,n)
    j=(land_index(l)-1)/t_i_length + 1
    i = land_index(l) - (j-1)*t_i_length
    rhokh1_prime = 1. / ( 1. / rhokh_1(l,n)                                   &
                       + GAMMA*dtrdz_1(i,j) )
    diff_lat_htf = (lc + lf)*ei_surft(l,n) + lc*ecan_surft(l,n)               &
                    + lc*esoil_surft(l,n) + lc*elake_surft(l,n)               &
                    - le_surft_old(l,n)
    dtstar = - diff_lat_htf /                                                 &
                ( cp*rhokh1_prime + ashtf_prime_surft(l,n) )
    diff_sens_htf = cp * rhokh1_prime * dtstar
    ftl_surft(l,n) = ftl_surft(l,n) + diff_sens_htf
    tstar_surft(l,n) = tstar_surft(l,n) + dtstar
    dftl(l) = dftl(l) + tile_frac(l,n)*diff_sens_htf
    dfqw(l) = dfqw(l) + tile_frac(l,n)*( ecan_surft(l,n) +                    &
                  esoil_surft(l,n) + ei_surft(l,n) + elake_surft(l,n)         &
                  - e_surft_old(l,n) )
  END DO
END DO

!-----------------------------------------------------------------------
! Update level 1 temperature and humidity and GBM heat and moisture
! fluxes due to limited moisture availability
!-----------------------------------------------------------------------
DO l=1,land_pts
  j=(land_index(l)-1)/t_i_length + 1
  i = land_index(l) - (j-1)*t_i_length
  ftl_1(i,j) = ftl_1(i,j) + fland(l)*dftl(l)
  fqw_1(i,j) = fqw_1(i,j) + fland(l)*dfqw(l)
END DO

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE sf_evap
