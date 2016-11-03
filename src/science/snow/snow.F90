! *****************************COPYRIGHT*******************************

! (c) [University of Edinburgh] [2009]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the JULES collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC237]

! *****************************COPYRIGHT*******************************
!  SUBROUTINE SNOW-----------------------------------------------------

! Description:
!     Calling routine for snow module

! Subroutine Interface:
MODULE snow_mod
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SNOW_MOD'
CONTAINS
SUBROUTINE snow ( land_pts,timestep,stf_hf_snow_melt,nsurft,                  &
                  surft_pts,surft_index,catch_snow,con_snow,                  &
                  con_rain,tile_frac,ls_snow,ls_rain,ei_surft,                &
                  hcaps,hcons,melt_surft,                                     &
                  smcl1,sthf1,surf_htf_surft,                                 &
                  tsoil1,tsurf_elev_surft,tstar_surft,                        &
                  v_sat1,rgrain,rgrainl,rho_snow_grnd,sice,                   &
                  sliq,snow_grnd,snow_surft,snowdepth,                        &
                  tsnow,nsnow,ds,hf_snow_melt,lying_snow,                     &
                  rho_snow,snomlt_sub_htf,snow_melt,                          &
                  snow_soil_htf, surf_ht_flux_ld, snow_smb_surft,             &
                  dhf_surf_minus_soil )

USE canopysnow_mod,  ONLY: canopysnow
USE compactsnow_mod, ONLY: compactsnow
USE layersnow_mod,   ONLY: layersnow
USE relayersnow_mod, ONLY: relayersnow
USE snowgrain_mod,   ONLY: snowgrain
USE snowpack_mod,    ONLY: snowpack
USE snowtherm_mod,   ONLY: snowtherm

USE c_lheat, ONLY :                                                           &
!  imported scalar parameters
 lf                 !  latent heat of fusion of water at 0degc

USE jules_snow_mod, ONLY :                                                    &
 nsmax                                                                        &
                    !  Maximum possible number of snow layers
,l_snow_infilt                                                                &
                    !  Include infiltration of rain into the snowpack
,r0                                                                           &
                    !  Grain size for fresh snow (microns)
,cansnowtile        !  switch for canopy snow model

USE jules_surface_types_mod, ONLY : lake

USE jules_radiation_mod, ONLY : l_snow_albedo, l_embedded_snow

USE jules_internal, ONLY : unload_backgrnd_pft
                    !  Background unloading rate (s-1)


USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Scalar arguments with intent(in)
INTEGER, INTENT(IN) ::                                                        &
 land_pts              ! Total number of land points

REAL, INTENT(IN) ::                                                           &
 timestep              ! Timestep length (s)

LOGICAL, INTENT(IN) ::                                                        &
 stf_hf_snow_melt      ! Stash flag for snowmelt heat flux

! Array arguments with intent(in)
INTEGER, INTENT(IN) ::                                                        &
 nsurft                                                                       &
                       !  Number of land tiles
,surft_pts(nsurft)                                                            &
                       ! Number of tile points
,surft_index(land_pts,nsurft)
                       ! Index of tile points

REAL, INTENT(IN) ::                                                           &
 catch_snow(land_pts,nsurft)                                                  &
                       ! canopy snow capacity (kg/m2)
,con_snow(land_pts)                                                           &
                       ! Convective snowfall rate (kg/m2/s)
,tile_frac(land_pts,nsurft)                                                   &
                       ! Tile fractions
,ls_snow(land_pts)                                                            &
                       ! Large-scale snowfall rate (kg/m2/s)
,ei_surft(land_pts,nsurft)                                                    &
                       ! Sublimation of snow (kg/m2/s)
,hcaps(land_pts)                                                              &
                       ! Soil heat capacity of top layer(J/K/m3).
,hcons(land_pts)                                                              &
                       ! Thermal conductivity of top soil layer,
!                            ! including water and ice (W/m/K)
,smcl1(land_pts)                                                              &
                       ! Moisture content of surface soil
!                            ! layer (kg/m2)
,sthf1(land_pts)                                                              &
                       ! Frozen soil moisture content of
!                            ! surface layer as a fraction of saturation.
,surf_htf_surft(land_pts,nsurft)                                              &
                       ! Surface heat flux (W/m2)
,tstar_surft(land_pts,nsurft)                                                 &
                       ! Tile surface temperature (K)
,v_sat1(land_pts)      ! Surface soil layer volumetric
!                            ! moisture concentration at saturation

! Array arguments with intent(inout)
REAL, INTENT(INOUT) ::                                                        &
 con_rain(land_pts)                                                           &
                       ! Convective rainfall rate (kg/m2/s)
,ls_rain(land_pts)
                       ! Large-scale rainfall rate (kg/m2/s)
REAL, INTENT(INOUT) ::                                                        &
 melt_surft(land_pts,nsurft)                                                  &
                       ! Surface or canopy snowmelt rate (kg/m2/s)
                       ! On output, this is the total melt rate
                       ! for the tile (i.e. sum of  melt on canopy
                       ! and ground).
,tsoil1(land_pts)                                                             &
                       ! Soil surface layer temperature (K)
,tsurf_elev_surft(land_pts,nsurft)                                            &
                       ! Temperature of elevated subsurface tiles (K)
,rgrain(land_pts,nsurft)                                                      &
                       ! Snow surface grain size (microns)
,rgrainl(land_pts,nsurft,nsmax)                                               &
                       ! Snow layer grain size (microns)
,rho_snow_grnd(land_pts,nsurft)                                               &
                       ! Snowpack bulk density (kg/m3)
,sice(land_pts,nsurft,nsmax)                                                  &
                       ! Ice content of snow layers (kg/m2)
,sliq(land_pts,nsurft,nsmax)                                                  &
                       ! Liquid content of snow layers (kg/m2)
,snow_grnd(land_pts,nsurft)                                                   &
                       ! Snow beneath canopy (kg/m2)
,snow_surft(land_pts,nsurft)                                                  &
                       ! Snow mass (Kg/m2)
,snowdepth(land_pts,nsurft)                                                   &
                       ! Snow depth (m)
,tsnow(land_pts,nsurft,nsmax)
                       ! Snow layer temperatures (K)

! Array arguments with intent(out)
INTEGER, INTENT(OUT) ::                                                       &
 nsnow(land_pts,nsurft)   ! Number of snow layers

REAL, INTENT(OUT) ::                                                          &
 ds(land_pts,nsurft,nsmax)                                                    &
                          ! Snow layer thicknesses (m)
,hf_snow_melt(land_pts)                                                       &
                          ! Gridbox snowmelt heat flux (W/m2)
,lying_snow(land_pts)                                                         &
                          ! Gridbox snow mass (kg m-2)
,rho_snow(land_pts,nsurft,nsmax)                                              &
                          ! Snow layer densities (kg/m3)
,snomlt_sub_htf(land_pts)                                                     &
                          ! Sub-canopy snowmelt heat flux (W/m2)
,snow_melt(land_pts)                                                          &
                          ! Gridbox snowmelt (Kg/m2/s)
,surf_ht_flux_ld(land_pts)                                                    &   
                          ! Surface heat flux on land (W/m2).
,snow_soil_htf(land_pts,nsurft)                                               &
                          ! Heat flux into the uppermost
                          ! subsurface layer (W/m2)
                          ! i.e. snow to ground, or into
                          ! snow/soil composite layer
,snow_smb_surft(land_pts,nsurft)                                              &
                          ! Rate of change of mass on tiles
,dhf_surf_minus_soil(land_pts)
                          ! heat flux difference
                          ! across the FLake snowpack (W/m2)

! Local scalars
INTEGER ::                                                                    &
 i                                                                            &
                      ! land point index and loop counter
,j                                                                            &
                      ! tile pts loop counter
,k                                                                            &
                      ! tile number
,n                                                                            &
                      ! tile loop counter
,n_can                ! Actual or dummy pointer to array
                      ! defined only on PFTs

! Local arrays
REAL ::                                                                       &
 csnow(land_pts,nsmax)                                                        &
                      ! Areal heat capacity of layers (J/K/m2)
,ksnow(land_pts,nsmax)                                                        &
                      ! Thermal conductivity of layers (W/m/K)
,rho0(land_pts)                                                               &
                      ! Density of fresh snow (kg/m3)
!                         ! Where NSNOW=0, rho0 is the density
!                         ! of the snowpack.
,snowfall(land_pts)                                                           &
                      ! Snowfall reaching the ground in timestep
                      ! (kg/m2) - includes any canopy unloading
,infiltration(land_pts)                                                       &
                      ! Infiltration of rainfall into snow pack
                      ! on the current tile (kg/m2)
,infil_rate_con_gbm(land_pts)                                                 &
                      ! Grid-box mean rate of infiltration of
                      ! convective rainfall into
                      ! snowpack (kg/m2/s)
,infil_rate_ls_gbm(land_pts)                                                  &
                      ! Grid-box mean rate of infiltration of
                      ! large-scale rainfall into
                      ! snowpack (kg/m2/s)
,snowmass(land_pts)                                                           &
                      ! Snow mass on the ground (Kg/m2)
,rgrain0(land_pts)                                                            &
                      ! Fresh snow grain size (microns)
,sice0(land_pts)                                                              &
                      ! Ice content of fresh snow (kg/m2)
!                           ! Where NSNOW=0, SICE0 is the mass of
!                           ! the snowpack.
,snow_can(land_pts,nsurft)                                                    &
                      ! Canopy snow load (Kg/m2)
,tsnow0(land_pts)                                                             & 
                      ! Temperature of fresh snow (K)
,snow_surft_old(land_pts,nsurft)
                      ! Copy of snow_surft

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SNOW'

!-----------------------------------------------------------------------
! Initialise gridbox variables.
!-----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

lying_snow(:) = 0.0
snomlt_sub_htf(:) = 0.0
snow_melt(:) = 0.0
infil_rate_con_gbm(:) = 0.0
infil_rate_ls_gbm(:)  = 0.0

snow_surft_old(:,:)=snow_surft(:,:)

! initialise FLake tile flux divergence
dhf_surf_minus_soil(:)  = 0.0

DO n=1,nsurft
!-----------------------------------------------------------------------
! Set snow mass variables
!-----------------------------------------------------------------------

  IF ( cansnowtile(n) ) THEN
!-----------------------------------------------------------------------
!         With the canopy snow model, snow_surft is held on the canopy,
!         while the mass of the snowpack (on the ground) is snow_grnd.
!-----------------------------------------------------------------------
    DO k=1,surft_pts(n)
      i = surft_index(k,n)
      snow_can(i,n) = snow_surft(i,n)
      snowmass(i) = snow_grnd(i,n)

!-----------------------------------------------------------------------
! Subtract sublimation and melt of canopy snow.
!-----------------------------------------------------------------------
      snow_can(i,n) = snow_can(i,n) -                                         &
                     ( ei_surft(i,n) + melt_surft(i,n) ) * timestep
    END DO

  ELSE

!-----------------------------------------------------------------------
! Without the snow canopy model, all the snow is in the snowpack
!-----------------------------------------------------------------------
    DO k=1,surft_pts(n)
      i = surft_index(k,n)
      snowmass(i) = snow_surft(i,n)
    END DO

  END IF

!-----------------------------------------------------------------------
! Canopy interception, throughfall and unloading of snow
!-----------------------------------------------------------------------
! The canopy model applies only on PFTs, so set a dummy pointer
! in the array on unused tiles.
  IF (cansnowtile(n)) THEN
    n_can = n
  ELSE
    n_can = 1
  END IF

  CALL canopysnow ( land_pts,surft_pts(n),timestep,cansnowtile(n),            &
                    surft_index(:,n),catch_snow(:,n),con_snow,                &
                    ls_snow,unload_backgrnd_pft(:,n_can),                     &
                    melt_surft(:,n), snow_can(:,n),snowfall )

!-----------------------------------------------------------------------
! Divide snow pack into layers
!-----------------------------------------------------------------------
CALL layersnow ( land_pts,surft_pts(n),surft_index(:,n),                    &
                   snowdepth(:,n),nsnow(:,n),ds(:,n,:) )

!-----------------------------------------------------------------------
! Thermal properties of snow layers
!-----------------------------------------------------------------------
  IF ( nsmax > 0 )                                                            &
    CALL snowtherm ( land_pts,surft_pts(n),nsnow(:,n),                        &
                     surft_index(:,n),ds(:,n,:),sice(:,n,:),                  &
                     sliq(:,n,:),csnow,ksnow )

!-----------------------------------------------------------------------
! Snow thermodynamics and hydrology
!-----------------------------------------------------------------------
  IF (l_snow_infilt) THEN
!   Where there is snow on the ground direct rainfall into infiltration.
    WHERE (nsnow(:,n) > 0)
      infiltration = ( ls_rain + con_rain ) * timestep
      infil_rate_con_gbm = infil_rate_con_gbm +                               &
        tile_frac(:,n) * con_rain
      infil_rate_ls_gbm = infil_rate_ls_gbm +                                 &
        tile_frac(:,n) * ls_rain
    ELSEWHERE
      infiltration = 0.0
    END WHERE
  ELSE
    infiltration(:) = 0.0
  END IF

  CALL snowpack ( land_pts,surft_pts(n),timestep,cansnowtile(n),              &
                  nsnow(:,n),surft_index(:,n),csnow,ei_surft(:,n),            &
                  hcaps,hcons,infiltration,                                   &
                  ksnow,rho_snow_grnd(:,n),smcl1,                             &
                  snowfall,sthf1,surf_htf_surft(:,n),                         &
                  tile_frac(:,n),v_sat1,ds(:,n,:),melt_surft(:,n),            &
                  sice(:,n,:),sliq(:,n,:),snomlt_sub_htf,                     &
                  snowdepth(:,n),snowmass,tsnow(:,n,:),                       &
                  tsoil1,tsurf_elev_surft(:,n),snow_soil_htf(:,n),            &
                  rho_snow(:,n,:), rho0,sice0,tsnow0 )

!-----------------------------------------------------------------------
! Growth of snow grains
!-----------------------------------------------------------------------
  IF ( l_snow_albedo .OR. l_embedded_snow ) THEN
    CALL snowgrain ( land_pts,surft_pts(n),timestep,nsnow(:,n),               &
                     surft_index(:,n),sice(:,n,:),snowfall,                   &
                     snowmass,tsnow(:,n,:),tstar_surft(:,n),                  &
                     rgrain(:,n),rgrainl(:,n,:),rgrain0 )
  ELSE
!   Deafult initialization required for bit-comparison in the UM.
    rgrain0(:) = r0
  END IF

  IF ( nsmax > 0 ) THEN
!-----------------------------------------------------------------------
! Mechanical compaction of snow
!-----------------------------------------------------------------------
    CALL compactsnow ( land_pts,surft_pts(n),timestep,nsnow(:,n),             &
                       surft_index(:,n),sice(:,n,:),sliq(:,n,:),              &
                       tsnow(:,n,:),rho_snow(:,n,:),ds(:,n,:) )

!-----------------------------------------------------------------------
! Redivide snowpack after changes in depth, conserving mass and energy
!-----------------------------------------------------------------------
    CALL relayersnow ( land_pts,surft_pts(n),surft_index(:,n),                &
                       rgrain0,rho0,sice0,snowfall,                           &
                       snowmass,tsnow0,nsnow(:,n),ds(:,n,:),                  &
                       rgrain(:,n),rgrainl(:,n,:),sice(:,n,:),                &
                       rho_snow_grnd(:,n),sliq(:,n,:),                        &
                       tsnow(:,n,:),rho_snow(:,n,:),                          &
                       snowdepth(:,n) )
  END IF  !  NSMAX>0

!-----------------------------------------------------------------------
! Copy into final snow mass variables
!-----------------------------------------------------------------------
  IF ( cansnowtile(n) ) THEN
    DO k=1,surft_pts(n)
      i = surft_index(k,n)
      snow_grnd(i,n) = snowmass(i)
      snow_surft(i,n) = snow_can(i,n)
    END DO
  ELSE
    DO k=1,surft_pts(n)
      i = surft_index(k,n)
      snow_surft(i,n) = snowmass(i)
    END DO
  END IF

!-----------------------------------------------------------------------
! Increment gridbox lying snow and snow melt.
!-----------------------------------------------------------------------
  DO k=1,surft_pts(n)
    i = surft_index(k,n)
    lying_snow(i) = lying_snow(i) +                                           &
                       tile_frac(i,n) * snow_surft(i,n)
!     Add snow beneath canopy.
    IF ( cansnowtile(n) )                                                     &
      lying_snow(i) = lying_snow(i) +                                         &
                         tile_frac(i,n) * snow_grnd(i,n)

!     Snow melt.
    snow_melt(i) = snow_melt(i) + tile_frac(i,n) * melt_surft(i,n)
  END DO

END DO  !  tiles

!-----------------------------------------------------------------------
! Calculate the total snowmelt heat flux.
!-----------------------------------------------------------------------
IF ( stf_hf_snow_melt ) THEN
  DO i=1,land_pts
    hf_snow_melt(i) = lf * snow_melt(i)
  END DO
END IF

!-----------------------------------------------------------------------
! Calculate gridbox surface heat flux over land.
!-----------------------------------------------------------------------
surf_ht_flux_ld(:) = 0.0
DO n=1,nsurft
  DO j=1,surft_pts(n)
    i = surft_index(j,n)
    surf_ht_flux_ld(i) = surf_ht_flux_ld(i) +                                 &
                          tile_frac(i,n) * snow_soil_htf(i,n)
    IF (n == lake) THEN
      dhf_surf_minus_soil(i) = surf_htf_surft(i,n) - snow_soil_htf(i,n)
    END IF
  END DO
END DO
!----------------------------------------------------------------------
! Remove the rainfall redirected into precipitation.
!----------------------------------------------------------------------
IF (l_snow_infilt) THEN
  con_rain(:) = con_rain(:) - infil_rate_con_gbm(:)
  ls_rain(:)  = ls_rain(:)  - infil_rate_ls_gbm(:)
END IF

snow_smb_surft(:,:)=(snow_surft(:,:)-snow_surft_old(:,:))/timestep

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE snow
END MODULE snow_mod
