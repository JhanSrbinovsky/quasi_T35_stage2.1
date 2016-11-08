! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Routine to calculate albedos of land-surface tiles and gridbox-mean
! albedo for JULES.


SUBROUTINE tile_albedo_cable (                                                &
 p_field,land_field,land_index,nsurft,surft_pts,                              &
 surft_index,l_aggregate,l_snow_albedo,albsoil,                               &
 albobs_sw, albobs_vis, albobs_nir,                                           &
 cosz,frac,lai_in,canht,rgrain,snow_surft,soot,tstar_surft,                   &
 z0_surft,ho2r2_orog,alb_surft,land_albedo,albobs_sc,                         &
 mype, timestep_number, cycleno, numcycles, row_length, rows,                 &
 fland, surf_down_sw, soil_alb, land_alb )

USE jules_print_mgr, ONLY :                                                   &
    jules_message,                                                            &
    jules_print

USE jules_surface_types_mod, ONLY : npft, ntype, lake, soil, ice,             &
                                    urban_canyon, elev_ice

USE jules_snow_mod,  ONLY : kland, maskd, tcland,                             &
                            rho_snow_const, rho_snow_fresh,                   &
                            cansnowtile,l_snowdep_surf,                       &
                            can_clump, lai_alb_lim_sn,                        &
                            n_lai_exposed,                                    &
                            amax, aicemax, rho_firn_albedo

USE c_0_dg_c
USE nvegparm
USE pftparm
USE water_constants_mod, ONLY : rho_ice

USE jules_snow_mod,  ONLY : kland, maskd, tcland, rho_snow_const,             &
                            cansnowtile,l_snowdep_surf

USE jules_surface_mod,   ONLY : l_point_data,                                 &
                                l_flake_model,                                &
                                l_elev_land_ice

USE jules_radiation_mod, ONLY : l_spec_albedo, l_spec_alb_bs,                 &
                                l_albedo_obs,                                 &
                                l_embedded_snow, l_mask_snow_orog

USE jules_vegetation_mod, ONLY : can_model

USE switches_urban, ONLY : l_urban2t, l_moruses_albedo
USE urban_param, ONLY : albwl_gb, albrd_gb, hwr_gb
USE lake_mod,    ONLY : albedo_whiteice_ref                                   &
                       ,albedo_blueice_ref                                    &
                       ,c_albice_MR                                           &
                       ,lake_h_ice_gb

USE ancil_info, ONLY: l_soil_point, l_lice_point

USE ereport_mod, ONLY : ereport
USE jules_mod, ONLY : snowdep_surft, albobs_scaling_surft

USE prognostics, ONLY : snowdepth_surft, rho_snow_grnd_surft,                 &
                        nsnow_surft, sice_surft, sliq_surft, ds_surft

!CABLE_LSM:
USE cable_rad_main_mod, ONLY : cable_rad_main

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE missing_data_mod, ONLY : rmdi
IMPLICIT NONE

! Subroutine arguments
!   Scalar arguments with intent(in):
INTEGER, INTENT(IN) ::                                                        &
 p_field                                                                      &
                             ! Total number of grid points.
,land_field                                                                   &
                             ! No. of land points.
,nsurft                      ! Number of surface tiles.

LOGICAL, INTENT(IN) ::                                                        &
 l_aggregate                                                                  &
                             ! IN Logical to set aggregate
!                                  !    surface scheme
,l_snow_albedo               ! .TRUE. for prognostic snow albedo.

!   Array arguments with intent(in):
INTEGER, INTENT(IN) ::                                                        &
 land_index(land_field)                                                       &
                             ! Index of land points.
,surft_pts(ntype)                                                             &
                             ! Number of land points which
!                                  ! include the nth surface type.
,surft_index(land_field,ntype)
                             ! Indices of land points which
!                                  ! include the nth surface type.

REAL, INTENT(IN) ::                                                           &
 rgrain(land_field,nsurft)
                             ! Snow grain size on tiles
!                                  ! (microns).
REAL, INTENT(IN) ::                                                           &
 albsoil(land_field)                                                          &
                             ! Soil albedo.
,albobs_sw(land_field)                                                        &
                             ! Observed snow-free sw albedo.
,albobs_vis(land_field)                                                       &
                             ! Observed snow-free vis albedo.
,albobs_nir(land_field)                                                       &
                             ! Observed snow-free nir albedo.
,cosz(p_field)                                                                &
                             ! Cosine of the zenith angle.
,frac(land_field,ntype)                                                       &
                             ! Fractional cover of each
!                                  ! surface type.
,lai_in(land_field,npft)                                                      &
                             ! Leaf area index.
,canht(land_field,npft)                                                       &
                             ! Canopy height
,snow_surft(land_field,nsurft)                                                &
                             ! Canopy snow on tiles (kg/m2)
,soot(p_field)                                                                &
                             ! Snow soot content (kg/kg).
,tstar_surft(land_field,nsurft)                                               &
                             ! Tile surface temperatures (K).
,z0_surft(land_field,nsurft)                                                  &
                             ! Surface roughness on tiles (m).
,ho2r2_orog(land_field)
                             ! Standard deviation of surface
                             ! orography

!CABLE_LSM:
INTEGER :: mype, timestep_number
integer :: cycleno, numcycles                ! ENDGAME cycle number
integer :: row_length, rows ! grid
real :: surf_down_sw(row_length,rows,4) ! 4-band ShortWave forcing
real :: fland(land_field)                 ! land fraction 
real :: soil_alb(land_field)              ! Snow-free, bare soil albedo: 
real :: land_alb(row_length,rows)         ! Mean land_albedo

LOGICAL, SAVE :: first_call = .true.
!CABLE_LSM: End

!   Array arguments with intent(out):
REAL, INTENT(OUT) ::                                                          &
 alb_surft(land_field,nsurft,4)                                               &
                              !Albedos for surface tiles.
!                                  !   (*,*,1) - Direct beam visible
!                                  !   (*,*,2) - Diffuse visible
!                                  !   (*,*,3) - Direct beam near-IR
!                                  !   (*,*,4) - Diffuse near-IR
,land_albedo(p_field,4)                                                       &
                              ! GBM albedos.
,albobs_sc(p_field,nsurft,2)  ! albedo scaling to obs in VIS and NIR
                              ! for diagnostics output by the UM

! Local arrays:
REAL ::                                                                       &
 albsnc(land_field,ntype)                                                     &
                             ! Snow-covered albedo of surf types.
,albsnf(land_field,ntype)                                                     &
                             ! Snow-free albedo of surf types.
,albsfsc(land_field,2)                                                        &
                             ! local scaling factor to match obs
,albobs_surft(land_field,ntype)                                               &
                             ! Albedo of the tiles (full veg) after
                             ! being scaled to obs
,alb_type(land_field,ntype,4)                                                 &
                             ! Albedos of surface types.
,snow_mass_alb(land_field)                                                    &
                             ! Mass of snow in calculation of
                             ! albedo
,alb_snow(land_field,ntype,4)                                                 &
                             ! Snow albedos.
,alb_snow_surft(land_field,4)                                                 &
                             ! Tiled snow albedo: used with
                             ! embedded snow.
,fsnow(land_field)                                                            &
                             ! Weighting factor for albedo.
,lai(land_field,npft)                                                         &
                             ! Adjusted leaf area index.
,snowd(land_field)                                                            &
                             ! Work variable (snow depth).
,tstar(land_field)                                                            &
                             ! Copy of TSTAR_SURFT.
,z0(land_field)                                                               &
                              ! Copy of Z0_SURFT.
,albsfm_sw(land_field)                                                        &
                             ! Model grid-box mean sf sw albedo
,albsfm_vis(land_field)                                                       &
                             ! Model grid-box mean sf vis albedo
,albsfm_nir(land_field)
                             ! Model grid-box mean sf nir albedo

INTEGER ::                  snow_pts(ntype)
                             ! Number of points with snow cover
INTEGER ::                  snow_index(land_field,ntype)
                             ! List of points with snow cover

INTEGER, PARAMETER ::       ilayers_dummy=1

INTEGER :: n_surft_pft
                            ! Index in a tiled array indicating
                            ! which tile contains the PFT being
                            ! considered at that time.
INTEGER :: n_surft_nvg
                            ! Index in a tiled array indicating
                            ! which tile contains the unvegetated
                            ! surface being considered at that time.

INTEGER ::                  albpft_call
!                            ! Option for call to albpft, matching
!                            ! the internal one


REAL ::                                                                       &
 fapar_dir_dummy(land_field,npft,ilayers_dummy)                               &
!                                 ! Profile of absorbed PAR -
!                                 ! Direct beam - DUMMY
,fapar_dif_dummy(land_field,npft,ilayers_dummy)                               &
!                                 ! Profile of absorbed PAR -
!                                 ! Diffuse beam -DUMMY
,fapar_dir2dif_dummy(land_field,npft,ilayers_dummy)                           &
!                                 ! DUMMY
,fapar_dif2dif_dummy(land_field,npft,ilayers_dummy)                           &
!                                 ! DUMMY
,fapar_dir2dir_dummy(land_field,npft,ilayers_dummy)                           &
!                                 ! DUMMY
,fsun_dummy(land_field,npft,ilayers_dummy)
!                                 ! DUMMY
REAL :: albudir(land_field,2,ntype)
!                                 ! Direct albedo of underlying
!                                 ! surface
REAL :: albudif(land_field,2,ntype)
!                                 ! Diffuse albedo of underlying
!                                 ! surface

! Local scalars:
REAL ::                                                                       &
 dsa                                                                          &
                             ! Deep-snow albedo.
,flit                                                                         &
                             ! Weighting factor for albedo.
,mask_orog                                                                    &
                             ! Orographic masking factor
,snowdepth_eff
                             ! Effective snow depth

REAL ::                                                                       &
 rho_snow_surf                                                                &
,snow_alb_vis_as                                                              &
,snow_alb_nir_as                                                              &
,ssum

INTEGER ::                                                                    &
 band,i,j,k,l,n              ! Loop counters

LOGICAL ::                                                                    &
 pointflag                   ! Switch for treatment of snow

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='TILE_ALBEDO'

INTEGER ::              errcode            ! Error code
CHARACTER (LEN = 80) :: errmsg             ! Error message text

! declare parameters used in ccontrol
!#include "chsunits.h"
!#include "csubmmax.h"
! We need this file for CAN_MODEL
!#include "ccontrol.h"


IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
DO n=1,nsurft
  DO band=1,4
    DO l=1,land_field
      alb_surft(l,n,band) = 0.
    END DO
  END DO
END DO
DO n=1,ntype
  DO band=1,4
    DO l=1,land_field
      alb_type(l,n,band) = 0.
      alb_snow(l,n,band) = 0.
    END DO
  END DO
END DO

! Equivalent snowdepth_surft for surface calculations.
snowdep_surft(:,:) = snowdepth_surft(:,:)
DO n=1,nsurft
  IF ( (can_model==4).AND.cansnowtile(n).AND.l_snowdep_surf ) THEN
    DO l=1,land_field
      snowdep_surft(l,n)=snow_surft(l,n)/rho_snow_const
    END DO
  END IF
END DO

! Impose minimum LAI for bare vegetation
DO n=1,npft
! Although this loop is over PFTs, with aggregation we
! have only one tile and snow loading, so point to the
! tile associated with the PFT.
  IF (.NOT.l_aggregate) THEN
    n_surft_pft=n
  ELSE
    n_surft_pft=1
  END IF
  DO j=1,surft_pts(n)
    l = surft_index(j,n)
    IF (snowdep_surft(l,n_surft_pft) > EPSILON(0.0)) THEN
      lai(l,n) = MAX( lai_in(l,n), lai_alb_lim_sn(n) )
    ELSE
      lai(l,n) = MAX( lai_in(l,n), lai_alb_lim(n) )
    END IF
  END DO
END DO


! scaling factor to get the model albedo to agree with the obs
! (initialise it as it goes into the module for use elsewhere)
IF (l_albedo_obs) THEN
  albobs_scaling_surft = 1.0
END IF
albobs_sc = rmdi

IF (l_spec_albedo) THEN
!----------------------------------------------------------------------
! Spectral albedo scheme, can have prognostic or diagnostic snow albedo
!----------------------------------------------------------------------

! Set albedos of vegetated surface types
! The logical argument (getProfile in albpft) is FALSE to indicate
! that profiles through the canopy should not be calculated.
! The 0 indicates that no scaling to obs is required
! Set the underlying albedo to the values for bare soil.
  DO n=1,npft
    DO l=1,land_field
      albudif(l,1,n) = albsoil(l)
      albudif(l,2,n) = albsoil(l)
      albudir(l,1,n) = albsoil(l)
      albudir(l,2,n) = albsoil(l)
    END DO
  END DO
! DEPENDS ON: albpft
  CALL albpft       (p_field,land_field,                                      &
                     land_index,surft_index,surft_pts,                        &
                     ilayers_dummy,.FALSE.,0,                                 &
                     albudir,albudif,cosz,lai,alb_type,                       &
                     fapar_dir_dummy,fapar_dif_dummy,                         &
                     fapar_dir2dif_dummy,fapar_dif2dif_dummy,                 &
                     fapar_dir2dir_dummy,fsun_dummy )

! Set albedos of non-vegetated surface types
  DO band=1,4
    DO n=npft+1,ntype
      DO j=1,surft_pts(n)
        l = surft_index(j,n)
        alb_type(l,n,band) = albsnf_nvg(n-npft)
        IF ( albsnf_nvg(n-npft) <  0. )                                       &
                                               ! Soil tile
          alb_type(l,n,band) = albsoil(l)
      END DO
    END DO
  END DO

! MORUSES: overwrite alb_type with values from canyonalb. Loop around band
! could be taken out as no snow involved at the moment.
  IF ( l_moruses_albedo ) THEN
    n = urban_canyon
    DO j=1,surft_pts(n)
      l = surft_index(j,n)
      i = land_index(l)
      DO band = 1,4
        ! DEPENDS ON: canyonalb
        CALL canyonalb(cosz(i),hwr_gb(l),albwl_gb(l),albrd_gb(l),             &
           alb_type(l,n,band))
      END DO
    END DO
  END IF

  IF ( l_albedo_obs ) THEN
! Average the model snow-free diffuse albedos (2 for VIS and 4 for NIR)
! over the grid box, by the tile fraction:
    albsfm_vis(:) = 0.0
    albsfm_nir(:) = 0.0
    DO n=1,ntype
      DO j=1,surft_pts(n)
        l = surft_index(j,n)
        albsfm_vis(l) = albsfm_vis(l) + ( alb_type(l,n,2) * frac(l,n) )
        albsfm_nir(l) = albsfm_nir(l) + ( alb_type(l,n,4) * frac(l,n) )
      END DO
    END DO

! Work out the scaling factor that needs to be applied to make the model
! snow-free albedo agree with the obs:
    DO l=1,land_field
      albsfsc(l,1) = albobs_vis(l) / albsfm_vis(l)
      albsfsc(l,2) = albobs_nir(l) / albsfm_nir(l)
      ! If point is land ice then do not do anything:
      IF ( l_lice_point(l) ) albsfsc(l,:) = 1.0
    END DO

! Recalculate the above albedos, but with scaled input albedos, within limits
! and store the scaling:
!
! starting with the non-veg (only need to do calculations twice, as the
! diffuse and direct albedos are the same):
    DO band=1,2
      DO n=npft+1,ntype
        DO j=1,surft_pts(n)
          l = surft_index(j,n)
          IF ( n == soil ) THEN
            alb_type(l,n,2*band-1) =                                          &
               MIN( MAX(albsoil(l) * albsfsc(l,band),                         &
               albsnf_nvgl(n-npft)), albsnf_nvgu(n-npft))
            albobs_scaling_surft(l,n,band) = alb_type(l,n,2*band-1)/ albsoil(l)
          ELSE IF ( n == urban_canyon .AND. l_moruses_albedo ) THEN
            ! MORUSES previously overwrites albsnf_nvg/alb_type with call to
            ! canyonalb, therefore scaling will currently overwrite canyonalb
          ELSE
            alb_type(l,n,2*band-1) =                                          &
               MIN(MAX(albsnf_nvg(n-npft)*albsfsc(l,band),                    &
               albsnf_nvgl(n-npft)), albsnf_nvgu(n-npft))
            albobs_scaling_surft(l,n,band) =                                  &
               alb_type(l,n,2*band-1)/albsnf_nvg(n-npft)
          END IF
        END DO
      END DO
    END DO
! now fill in the diffuse albedo from the direct value:
    DO band=1,2
      DO n=npft+1,ntype
        DO j=1,surft_pts(n)
          l = surft_index(j,n)
          alb_type(l,n,2*band) = alb_type(l,n,2*band-1)
        END DO
      END DO
    END DO

!   now do the veg tiles, by calling albpft again and telling it to
!   scale the vegetation scattering and reflectivity parameters:
!
!   Put the intended scaling into albobs_scaling_surft for each PFT
!   and then albpft will correct it for the non-linearity
!   (which is what the 1 indicates):
    DO band=1,2
      DO n=1,npft
        DO l=1,land_field
          albobs_scaling_surft(l,n,band) = albsfsc(l,band)
        END DO
      END DO
    END DO
!
!   Scale the soil albedo using the provisional scaling
    DO n=1,npft
      DO l=1,land_field
        albudif(l,1,n) = albsoil(l) * albobs_scaling_surft(l,soil,1)
        albudif(l,2,n) = albsoil(l) * albobs_scaling_surft(l,soil,2)
        albudir(l,1,n) = albsoil(l) * albobs_scaling_surft(l,soil,1)
        albudir(l,2,n) = albsoil(l) * albobs_scaling_surft(l,soil,2)
      END DO
    END DO
!
!   DEPENDS ON: albpft
    CALL albpft       (p_field,land_field,                                    &
                       land_index,surft_index,surft_pts,                      &
                       ilayers_dummy,.FALSE.,1,                               &
                       albudir,albudif,cosz,lai,alb_type,                     &
                       fapar_dir_dummy,fapar_dif_dummy,                       &
                       fapar_dir2dif_dummy,fapar_dif2dif_dummy,               &
                       fapar_dir2dir_dummy,fsun_dummy )
!
!
  END IF ! End of test on l_albedo_obs
!
! ----------------------------------------------------------------
! If snow is embedded in the canopy calculate the albedo for the
! snow on the ground, recalculate for the exposed LAI and for
! any snow retained on the canopy. This is not compatible with
! the aggregation of tiled properties.
! ----------------------------------------------------------------
  IF (l_embedded_snow) THEN
!
    DO n=1, npft
!
!     Point to the tile associated with the PFT.
      IF (.NOT.l_aggregate) THEN
        n_surft_pft=n
      ELSE
        n_surft_pft=1
      END IF
!
!     Gather snow-covered points.
      snow_pts(n)=0
      DO j=1,surft_pts(n)
        l = surft_index(j,n)
!       Allow a numerical tolerance for checking snow depth.
        IF (snowdepth_surft(l,n_surft_pft) >                                  &
            128 * EPSILON(snowdepth_surft)) THEN
          snow_pts(n) = snow_pts(n) + 1
          snow_index(snow_pts(n),n) = l
        END IF
      END DO
!
!     Calculate the albedos for the snow on the ground.
      IF ( (can_model==4).AND.cansnowtile(n) ) THEN
        DO j=1,snow_pts(n)
          l = snow_index(j,n)
!         Provisional check on the density. In a future upgrade
!         reconfiguration will ensure that the density of snow on
!         the ground is properly initialized.
          snow_mass_alb(l) = MIN(rho_ice,                                     &
            MAX(rho_snow_grnd_surft(l,n_surft_pft),rho_snow_fresh) ) *        &
            snowdepth_surft(l,n_surft_pft)
        END DO
      ELSE
        DO j=1,snow_pts(n)
          l = snow_index(j,n)
          snow_mass_alb(l) = snow_surft(l,n_surft_pft)
        END DO
      END IF
!     The albedo of the underlying surface is now the scaled value for
!     bare soil.
! DEPENDS ON: albsnow_ts
      CALL albsnow_ts(p_field,land_field,land_index,                          &
                   snow_pts(n), snow_index(:,n),                              &
                   cosz,albudir(:,:,n),albudif(:,:,n),                        &
                   rgrain(:,n_surft_pft),snow_mass_alb,soot,                  &
                   alb_snow_surft)
!     Prepare to calculate the albedo for the exposed canopy. Set
!     the appropriate exposed LAI and set the snow albedo on
!     the underlying surface.
      DO j=1,snow_pts(n)
        l = snow_index(j,n)
        snowdepth_eff = MAX(0.0, snowdepth_surft(l,n_surft_pft))
!
        IF (l_mask_snow_orog) THEN
!         Include orographic masking, based on Roesch et al. (2006),
!         J. Clim., vol. 19, p.3828. Note that they record the snow
!         loading using SWE in m. Here we take a typical snow
!         density as 200 kgm-3.
          snowdepth_eff = snowdepth_eff *                                     &
                ( snowdepth_eff /                                             &
                ( snowdepth_eff + EPSILON(1.0) +                              &
                                  0.00075*ho2r2_orog(l)) )
        END IF
!
        IF (snowdepth_eff < canht(l,n)) THEN
          lai(l,n) = lai(l, n) *                                              &
            (1.0 - snowdepth_eff / canht(l,n) ) ** n_lai_exposed(n)
        ELSE
          lai(l,n) = 0.0
        END IF
!
        albudir(l,1,n) = alb_snow_surft(l,1)
        albudif(l,1,n) = alb_snow_surft(l,2)
        albudir(l,2,n) = alb_snow_surft(l,3)
        albudif(l,2,n) = alb_snow_surft(l,4)
      END DO
    END DO
!
! ----------------------------------------------------------------
!   Calculate the albedo of the canopy using the exposed LAI.
!   embeded snow is assumed).
! ----------------------------------------------------------------
    IF (l_albedo_obs) THEN
!     We have already calculated the scaling without snow. Iterating
!     again would scale the transmission and reflection coefficients
!     by a snow albedo and would be wrong.
      albpft_call = 2
    ELSE
!     Just calculate the albedo directly with the exposed LAI.
      albpft_call = 0
    END IF
!   This call is made only at snow points, with the underlying albedo
!   set for snow.
!   DEPENDS ON: albpft
    CALL albpft       (p_field,land_field,                                    &
                       land_index,snow_index,snow_pts,                        &
                       ilayers_dummy,.FALSE.,albpft_call,                     &
                       albudir,albudif,cosz,lai,alb_type,                     &
                       fapar_dir_dummy,fapar_dif_dummy,                       &
                       fapar_dir2dif_dummy,fapar_dif2dif_dummy,               &
                       fapar_dir2dir_dummy,fsun_dummy )
!
! ----------------------------------------------------------------
! If embedded snow is assumed with a canopy model, allow for
! snow on the canopy. Then, assign final albedos.
! ----------------------------------------------------------------
!
    DO n=1, npft
!
!     Point to the tile associated with the PFT.
      IF (.NOT.l_aggregate) THEN
        n_surft_pft=n
      ELSE
        n_surft_pft=1
      END IF
!
!     We need only reset alb_type if there is a canopy model --
!     at this point alb_type holds the albedo up to the top of
!     the canopy, but excluding any snow thereon. Snow is clumped
!     on to a fraction 1/can_clump of the canopy and the albedo
!     is weighted accordingly.
      IF ( (can_model==4).AND.cansnowtile(n) ) THEN
        DO j=1,snow_pts(n)
          l = snow_index(j,n)
          albudir(l,1,n) = alb_type(l,n,1)
          albudif(l,1,n) = alb_type(l,n,2)
          albudir(l,2,n) = alb_type(l,n,3)
          albudif(l,2,n) = alb_type(l,n,4)
          snow_mass_alb(l) = snow_surft(l,n_surft_pft) * can_clump(n)
        END DO
! DEPENDS ON: albsnow_ts
        CALL albsnow_ts(p_field,land_field,land_index,                        &
                     snow_pts(n), snow_index(:,n),                            &
                     cosz,albudir(:,:,n),albudif(:,:,n),                      &
                     rgrain(:,n_surft_pft),snow_mass_alb,soot,                &
                     alb_snow_surft)
!
!       Pass the final albedos back into alb_type.
        DO band=1,4
          DO j=1,snow_pts(n)
            l = snow_index(j,n)
            alb_type(l,n,band) = ( alb_snow_surft(l,band) +                   &
              (can_clump(n) - 1.0) * alb_type(l,n,band) ) /                   &
              can_clump(n)
          END DO
        END DO
!
      END IF
!
    END DO
!
!   Finish off with the non-vegetated tiles.
    DO n=npft+1, ntype
!
!     Point to the tile associated with the surface type.
      IF (.NOT.l_aggregate) THEN
        n_surft_nvg=n
      ELSE
        n_surft_nvg=1
      END IF
!
!     Gather snow-covered points.
      snow_pts(n)=0
      DO j=1,surft_pts(n)
        l = surft_index(j,n)
!       Allow a numerical tolerance for checking snow depth.
        IF (snowdepth_surft(l,n_surft_nvg) >                                  &
            128 * EPSILON(snowdepth_surft)) THEN
          snow_pts(n) = snow_pts(n) + 1
          snow_index(snow_pts(n),n) = l
        END IF
      END DO
!
!     The (adjusted) snow-free albedo will have been set on
!     unvegetated tiles
      DO j=1,snow_pts(n)
        l = snow_index(j,n)
        albudir(l,1,n) = alb_type(l,n,1)
        albudif(l,1,n) = alb_type(l,n,2)
        albudir(l,2,n) = alb_type(l,n,3)
        albudif(l,2,n) = alb_type(l,n,4)
      END DO
!
! DEPENDS ON: albsnow_ts
      CALL albsnow_ts(p_field,land_field,land_index,                          &
                   snow_pts(n), snow_index(:,n),                              &
                   cosz,albudir(:,:,n),albudif(:,:,n),                        &
                   rgrain(:,n_surft_nvg),snow_surft(:,n),soot,                &
                   alb_snow_surft)

!-----------------------------------------------------------------------------
! For land ice surfaces where deep, dense snow may be emulating firn/bare ice
! scattering physics as in albedo_ts less valid. As in MAR, scale albedo above
! threshold with surface density (Gruell and Konzellmann '94) using ~ top 10cm 
!-----------------------------------------------------------------------------
      IF (l_elev_land_ice) THEN
        DO j=1,snow_pts(n)
          l = snow_index(j,n)
          IF (l_lice_point(l) .AND. nsnow_surft(l,n) > 0) THEN

            ssum=0.
            DO k=1,nsnow_surft(l,n)
              ssum=ssum+ds_surft(l,n,k)
              IF (ssum > 0.1) EXIT
            END DO
            k=min(k,nsnow_surft(l,n))

            IF (sum(ds_surft(l,n,1:k)) > 1e-3) THEN
              rho_snow_surf=( sum(sice_surft(l,n,1:k))+sum(sliq_surft(l,n,1:k)) )   &
                            / sum(ds_surft(l,n,1:k))
            ELSE
              rho_snow_surf=rho_snow_const
            END IF

            IF (rho_snow_surf > rho_firn_albedo) THEN
              snow_alb_vis_as=aicemax(1) + ( rho_snow_surf - rho_ice) * &
                             (  (amax(1) - aicemax(1))  /               &
                                (rho_snow_const - rho_ice) )

              snow_alb_nir_as=aicemax(2) + ( rho_snow_surf - rho_ice) * &
                           (  (amax(2) - aicemax(2))  /                 &
                              (rho_snow_const - rho_ice) )

              alb_snow_surft(l,1) = MIN(alb_snow_surft(l,1),snow_alb_vis_as)
              alb_snow_surft(l,2) = MIN(alb_snow_surft(l,2),snow_alb_vis_as)
              alb_snow_surft(l,3) = MIN(alb_snow_surft(l,3),snow_alb_nir_as)
              alb_snow_surft(l,4) = MIN(alb_snow_surft(l,4),snow_alb_nir_as)
            ELSE
              alb_snow_surft(l,1) = MIN(alb_snow_surft(l,1),amax(1))
              alb_snow_surft(l,2) = MIN(alb_snow_surft(l,2),amax(1))
              alb_snow_surft(l,3) = MIN(alb_snow_surft(l,3),amax(2))
              alb_snow_surft(l,4) = MIN(alb_snow_surft(l,4),amax(2))
            END IF
          END IF
        END DO
      END IF

      IF (l_mask_snow_orog) THEN
!
!       Include orographic masking, based on Roesch et al. (2006),
!       J. Clim., vol. 19, p.3828, where we take a typical snow
!       density as 200 kgm-3. Note that they measure SWE in m.
!       Note also that we are relying on the previous initialization
!       of alb_type to the background value, which will not be
!       adjusted if snow is present.
        DO j=1, snow_pts(n)
          l = snow_index(j,n)
          snowdepth_eff = MAX(0.0, snowdepth_surft(l,n_surft_nvg))
          mask_orog = snowdepth_eff /                                         &
                ( snowdepth_eff + EPSILON(1.0) +                              &
                                0.00075*ho2r2_orog(l) )
          DO band=1,4
            alb_type(l,n,band) = alb_type(l,n,band) + mask_orog *             &
              (alb_snow_surft(l,band) - alb_type(l,n,band))
          END DO
        END DO
!
      ELSE
!
        DO band=1,4
          DO j=1, snow_pts(n)
            l = snow_index(j,n)
            alb_type(l,n,band) = alb_snow_surft(l,band)
          END DO
        END DO
!
      END IF
!
    END DO
!
  END IF

IF (l_albedo_obs) THEN
! ! Check on 0 > albedo > 1
    DO band=1,4
      DO n=1,ntype
      !DO l=1,land_field
        DO j=1,surft_pts(n)
          l = surft_index(j,n)
          IF ( (alb_type(l,n,band) < 0.) .OR. (alb_type(l,n,band) > 1.) ) THEN
            WRITE(jules_message,'(a,i5,F24.12)')                              &
                'ERROR: Albedo < 0.0 or > 1.0: ',l,alb_type(l,n,band)
            CALL jules_print('tile_albeno_jls',jules_message)
            WRITE(jules_message,'(a,I3,a,I3)') 'Tile number', n, ', Band', band
            CALL jules_print('tile_albeno_jls',jules_message)
            WRITE(jules_message,'(a,F24.12)') 'Scaling', albsfsc(l,band)
            CALL jules_print('tile_albeno_jls',jules_message)
            IF (band <= 2) THEN
              WRITE(jules_message,'(a,F24.12,F24.12,F24.12)')                 &
                  'obs, model, soil:',albobs_vis(l), albsfm_vis(l), albsoil(l)
              CALL jules_print('tile_albeno_jls',jules_message)
            ELSE
              WRITE(jules_message,'(a,F24.12,F24.12,F24.12)')                 &
                  'obs, model, soil:',albobs_nir(l), albsfm_nir(l), albsoil(l)
              CALL jules_print('tile_albeno_jls',jules_message)
            END IF
            errcode = 1
            errmsg  = 'Unphysical albedos being created'
            CALL ereport ('tile_albedo',errcode,errmsg)
          END IF
        END DO
      END DO
    END DO

  END IF ! ends l_albedo_obs for the spectral albedo scheme

  IF ( l_spec_alb_bs ) THEN
  ! If the spectral albedos are being used with a single ('blue' sky) value for
  ! the direct ('black' sky) and diffuse ('white' sky) beams, then change the
  ! direct albedos to the diffuse values:
    DO band=1,2
      ! at the moment only the veg tiles have a different direct/diffuse value:
      DO n=1,npft
        DO j=1,surft_pts(n)
          l = surft_index(j,n)
          alb_type(l,n,2*band-1) = alb_type(l,n,2*band)
        END DO
      END DO
    END DO
  END IF

! Re-set albedos of frozen lakes if FLake is used
! using the algorithm from the flake interface program
! (reference Mironov & Ritter 2004).
  IF ( l_flake_model.AND.(.NOT.l_aggregate) ) THEN
    n = lake
    DO j=1,surft_pts(n)
      l = surft_index(j,n)
      IF ((lake_h_ice_gb(l) > 0.0).AND.(tstar_surft(l,n) <= tm)) THEN
        alb_type(l,n,:) =   albedo_whiteice_ref                               &
                          + EXP( -c_albice_MR*(tm-tstar_surft(l,n)) / tm )    &
                          * (albedo_blueice_ref - albedo_whiteice_ref)
      END IF
    END DO
  END IF

! ---------------------------------------------------------------------
! Assign snow albedos using the default scheme, essentially assuming
! that the snow is on the canopy.
! ---------------------------------------------------------------------
  IF (.NOT.l_embedded_snow) THEN
    IF (l_snow_albedo) THEN
!----------------------------------------------------------------------
! Spectral albedo scheme with prognostic snow albedo
!----------------------------------------------------------------------
! Calculate snow albedos
! DEPENDS ON: albsnow
      CALL albsnow(p_field,land_field,land_index,                             &
                   nsurft,surft_index,surft_pts,l_aggregate,                  &
                   cosz,rgrain,snowdep_surft,soot,alb_snow)

!-----------------------------------------------------------------------------
! For land ice sheet surfaces where deep, dense snow may be emulating firn/bare ice
! scattering physics as in albedo_ts less valid. As in MAR, scale albedo above
! threshold with surface density (Gruell and Konzellmann '94) using ~ top 10cm 
!-----------------------------------------------------------------------------
      IF (l_elev_land_ice .AND. maxval(elev_ice) > 0) THEN
        DO n=minval(elev_ice,(elev_ice > 0)),maxval(elev_ice)
        DO j=1,surft_pts(n)
          l = surft_index(j,n)
          IF (l_lice_point(l) .AND. nsnow_surft(l,n) > 0) THEN

            ssum=0.
            DO k=1,nsnow_surft(l,n)
              ssum=ssum+ds_surft(l,n,k)
              if (ssum > 0.1) EXIT
            END DO
            k=MIN(k,nsnow_surft(l,n))

            IF (sum(ds_surft(l,n,1:k)) > 1e-3) THEN
              rho_snow_surf=( sum(sice_surft(l,n,1:k))+sum(sliq_surft(l,n,1:k)) )   &
                            / sum(ds_surft(l,n,1:k))
            ELSE
              rho_snow_surf=rho_snow_const
            END IF

            IF (rho_snow_surf > rho_firn_albedo) THEN
              snow_alb_vis_as=aicemax(1) + ( rho_snow_surf - rho_ice) * &
                             (  (amax(1) - aicemax(1))  /               &
                                (rho_snow_const - rho_ice) )

              snow_alb_nir_as=aicemax(2) + ( rho_snow_surf - rho_ice) * &
                           (  (amax(2) - aicemax(2))  /                 &
                              (rho_snow_const - rho_ice) )

              alb_snow(l,n,1) = MIN(alb_snow(l,n,1),snow_alb_vis_as)
              alb_snow(l,n,2) = MIN(alb_snow(l,n,2),snow_alb_vis_as)
              alb_snow(l,n,3) = MIN(alb_snow(l,n,3),snow_alb_nir_as)
              alb_snow(l,n,4) = MIN(alb_snow(l,n,4),snow_alb_nir_as)
            ELSE
              alb_snow(l,n,1) = MIN(alb_snow(l,n,1),amax(1))
              alb_snow(l,n,2) = MIN(alb_snow(l,n,2),amax(1))
              alb_snow(l,n,3) = MIN(alb_snow(l,n,3),amax(2))
              alb_snow(l,n,4) = MIN(alb_snow(l,n,4),amax(2))
            END IF
          END IF
        END DO
        END DO
      END IF

! Adjust surface type albedos for snow cover
      DO l=1,land_field
        snowd(l) = snowdep_surft(l,1)
        z0(l) = z0_surft(l,1)
      END DO
      DO n=1,ntype
        IF (.NOT. l_aggregate) THEN
          DO j=1,surft_pts(n)
            l = surft_index(j,n)
            snowd(l) = snowdep_surft(l,n)
            z0(l) = z0_surft(l,n)
          END DO
        END IF
! Calculate snow albedo weighting factor.
        fsnow(:) = 0.0
        IF ( l_point_data ) THEN
          DO j=1,surft_pts(n)
            l = surft_index(j,n)
            IF ( snowd(l) > 0.) fsnow(l) = 1.0 - EXP( -50.0*snowd(l) )
          END DO
        ELSE
          DO j=1,surft_pts(n)
            l = surft_index(j,n)
            IF ( snowd(l) > 0.) fsnow(l) = snowd(l) /                         &
                 ( snowd(l) + 10.*z0(l) )
          END DO
        END IF
! Calculate weighted tile albedo.
        DO j=1,surft_pts(n)
          l = surft_index(j,n)
          DO band=1,4
            alb_type(l,n,band) = fsnow(l)*alb_snow(l,n,band)                  &
                                + (1. - fsnow(l))*alb_type(l,n,band)
          END DO

        END DO
      END DO   !  ntype

    ELSE
!----------------------------------------------------------------------
! spectral albedo scheme with diagnosed snow albedo
!----------------------------------------------------------------------
! Adjust surface type albedos for snow cover
      DO l=1,land_field
        tstar(l) = tstar_surft(l,1)
        snowd(l) = snowdep_surft(l,1)
      END DO
! Set albedos of snow covered vegetated surface types:
      DO n=1,npft
        DO j=1,surft_pts(n)
          l = surft_index(j,n)
          flit = 1.0 - EXP(-kext(n)*lai(l,n))
          albsnc(l,n) = albsnc_min(n)*(1 - flit) + albsnc_max(n)*flit
        END DO
      END DO
! and the snow covered non-veg types:
      DO n=npft+1,ntype
        DO j=1,surft_pts(n)
          l = surft_index(j,n)
          albsnc(l,n) = albsnc_nvg(n-npft)
        END DO
      END DO
! now apply the snow covered albedos to the snow free ones in alb_type
      DO band=1,4
        DO n=1,ntype
          IF (.NOT. l_aggregate) THEN
            DO j=1,surft_pts(n)
              l = surft_index(j,n)
              tstar(l) = tstar_surft(l,n)
              snowd(l) = snowdep_surft(l,n)
            END DO
          END IF
          DO j=1,surft_pts(n)
            l = surft_index(j,n)
            IF ( tstar(l)  <   tcland ) THEN
              dsa = albsnc(l,n)
            ELSE IF ( tstar(l)  <   tm ) THEN
              dsa = albsnc(l,n) + kland*(alb_type(l,n,band) - albsnc(l,n)) *  &
                                     (tstar(l) - tcland)
            ELSE
              dsa = albsnc(l,n) + kland*(alb_type(l,n,band) - albsnc(l,n)) *  &
                                     (tm - tcland)
            END IF
            alb_type(l,n,band) = alb_type(l,n,band) +                         &
                                   (dsa-alb_type(l,n,band)) *                 &
                                   ( 1. - EXP(-maskd*snowd(l)) )
          END DO
        END DO
      END DO

    END IF ! ends test on snow scheme for spectral albedo
!
  END IF ! ends test on embedded snow.

ELSE
!----------------------------------------------------------------------
! Non-spectral albedo scheme with diagnosed snow albedo
!----------------------------------------------------------------------
  IF (l_snow_albedo) THEN
    errcode = 1
    errmsg  = 'l_snow_albedo is dependent on l_spec_albedo'
    CALL ereport ('tile_albedo',errcode,errmsg)
  END IF

! Set albedos of vegetated surface types
  DO n=1,npft
    DO j=1,surft_pts(n)
      l = surft_index(j,n)
      flit = 1.0 - EXP(-kext(n)*lai(l,n))
      albsnc(l,n) = albsnc_min(n)*(1 - flit) + albsnc_max(n)*flit
      albsnf(l,n) = albsoil(l)*(1 - flit) + albsnf_max(n)*flit
    END DO
  END DO

! Set albedos of non-vegetated surface types
  DO n=npft+1,ntype
    DO j=1,surft_pts(n)
      l = surft_index(j,n)
      albsnc(l,n) = albsnc_nvg(n-npft)
      albsnf(l,n) = albsnf_nvg(n-npft)
      IF ( albsnf_nvg(n-npft) <  0. ) albsnf(l,n) = albsoil(l)
    END DO
  END DO

! MORUSES: Set canyon albedos
  IF ( l_moruses_albedo ) THEN
    n = urban_canyon
    DO j=1,surft_pts(n)
      l = surft_index(j,n)
      i = land_index(l)
      ! DEPENDS ON: canyonalb
      CALL canyonalb( cosz(i), hwr_gb(l), albwl_gb(l), albrd_gb(l),           &
         albsnf(l,n) )
    END DO
  END IF

  IF ( l_albedo_obs ) THEN
! Average the model snow-free albedo over the grid box, by the tile fraction:
    albsfm_sw(:) = 0.0
    DO n=1,ntype
      DO j=1,surft_pts(n)
        l = surft_index(j,n)
        albsfm_sw(l) = albsfm_sw(l) + ( albsnf(l,n) * frac(l,n) )
      END DO
    END DO

! Work out the scaling factor that needs to be applied to make the model
! albedo agree with the obs:
    DO l=1,land_field
      albsfsc(l,1) = albobs_sw(l) / albsfm_sw(l)
      ! If point is land ice then do not do anything:
      IF ( l_lice_point(l) ) albsfsc(l,1) = 1.0
    END DO

! Recalculate the above albedos, but with scaled input albedos, within limits:
! scale the bare soil albedo first:
    DO l=1,land_field
      albobs_surft(l,soil)=albsoil(l)*albsfsc(l,1)
    END DO
! apply the limits:
    DO l=1,land_field
      albobs_surft(l,soil) =MIN(                                              &
                                MAX(albobs_surft(l,soil),                     &
                                    albsnf_nvgl(soil-npft)),                  &
                                albsnf_nvgu(soil-npft))
    END DO
! and make a note of what the scaling ended up being, after the limits:
    DO l=1,land_field
      albobs_scaling_surft(l,soil,1) = albobs_surft(l,soil) / albsoil(l)
    END DO

! now do the veg tiles
    DO n=1,npft
       !DO l=1,land_field
       DO j=1,surft_pts(n)
         l = surft_index(j,n)
! scale the fully veg albedo, and apply limits:
         albobs_surft(l,n) = MIN(                                             &
                                 MAX(albsnf_max(n)*albsfsc(l,1),              &
                                     albsnf_maxl(n)),                         &
                                 albsnf_maxu(n))
! work out the albedo of the tile, with bare soil under the actual partial veg:
         flit = 1.0 - EXP(-kext(n)*lai(l,n))
         albsnf(l,n) = (1 - flit)*albobs_surft(l,soil) + flit*albobs_surft(l,n)
! and store the scaling:
         albobs_scaling_surft(l,n,1) =  albobs_surft(l,n) / albsnf_max(n)
      END DO
    END DO

! non-vegetated surface types:
    DO n=npft+1,ntype
      IF ( n == soil ) THEN
        DO j=1,surft_pts(n)
          l = surft_index(j,n)
          albsnf(l,n) = albobs_surft(l,soil)
        END DO
      ELSE IF ( n == urban_canyon .AND. l_moruses_albedo ) THEN
        ! MORUSES previously overwrites albsnf_nvg/albsnf with call to
        ! canyonalb, therefore scaling will currently overwrite canyonalb
      ELSE
        DO j=1,surft_pts(n)
          l = surft_index(j,n)
! apply scaling and limits:
          albsnf(l,n) = MIN( MAX(albsnf_nvg(n-npft)*albsfsc(l,1),             &
                         albsnf_nvgl(n-npft)), albsnf_nvgu(n-npft))
          albobs_surft(l,n)= albsnf(l,n)
! and store the scaling
          albobs_scaling_surft(l,n,1) =  albobs_surft(l,n) / albsnf_nvg(n-npft)
        END DO
      END IF
    END DO

! ! Check on 0 > albedo > 1
   DO n=1,ntype
    !DO l=1,land_field
      DO j=1,surft_pts(n)
        l = surft_index(j,n)
        IF ( (albsnf(l,n) < 0.) .OR. (albsnf(l,n) > 1.) ) THEN
          WRITE(jules_message,'(a,F24.12)')                                   &
             'ERROR: Albedo < 0.0 or > 1.0: ', albsnf(l,n)
          CALL jules_print('tile_albeno_jls',jules_message)
          WRITE(jules_message,'(a,I3)') 'Tile number', n
          CALL jules_print('tile_albeno_jls',jules_message)
          WRITE(jules_message,'(a,F24.12)') 'Scaling', albsfsc(l,1)
          CALL jules_print('tile_albeno_jls',jules_message)
          WRITE(jules_message,'(a,F24.12,F24.12,F24.12)')                     &
             'obs, model, soil:', albobs_sw(l), albsfm_sw(l), albsoil(l)
          CALL jules_print('tile_albeno_jls',jules_message)
          errcode = 1
          errmsg  = 'Unphysical albedos being created'
          CALL ereport ('tile_albedo',errcode,errmsg)
        END IF
      END DO
    END DO

  END IF ! ends l_albedo_obs for the non-spectral albedo scheme

! Re-set albedos of frozen lakes if FLake is used
! using the algorithm from the flake interface program
! (reference Mironov & Ritter 2004).
  IF ( l_flake_model.AND.(.NOT.l_aggregate) ) THEN
    n = lake
    DO j=1,surft_pts(n)
      l = surft_index(j,n)
      IF ((lake_h_ice_gb(l) > 0.0).AND.(tstar_surft(l,n) <= tm)) THEN
        albsnf(l,n) =   albedo_whiteice_ref                                   &
                      + EXP( -c_albice_MR*(tm-tstar_surft(l,n)) / tm )        &
                      * (albedo_blueice_ref - albedo_whiteice_ref)
      END IF
    END DO
  END IF

! Adjust surface type albedos for snow cover
  DO l=1,land_field
    tstar(l) = tstar_surft(l,1)
    snowd(l) = snowdep_surft(l,1)
  END DO
  DO n=1,ntype
    IF (.NOT. l_aggregate) THEN
      DO j=1,surft_pts(n)
        l = surft_index(j,n)
        tstar(l) = tstar_surft(l,n)
        snowd(l) = snowdep_surft(l,n)
      END DO
    END IF
    DO j=1,surft_pts(n)
      l = surft_index(j,n)
      IF ( tstar(l)  <   tcland ) THEN
        dsa = albsnc(l,n)
      ELSE IF ( tstar(l)  <   tm ) THEN
        dsa = albsnc(l,n) + kland*(albsnf(l,n) - albsnc(l,n))                 &
                                 *(tstar(l) - tcland)
      ELSE
        dsa = albsnc(l,n) + kland*(albsnf(l,n) - albsnc(l,n))                 &
                                 *(tm - tcland)
      END IF
      alb_type(l,n,1) = albsnf(l,n) + (dsa - albsnf(l,n)) *                   &
                          ( 1. - EXP(-maskd*snowd(l)) )
    END DO
  END DO

! Copy albedo to all bands
  DO band=2,4
    DO n=1,ntype
      DO j=1,surft_pts(n)
        l = surft_index(j,n)
        alb_type(l,n,band) = alb_type(l,n,1)
      END DO
    END DO
  END DO

END IF       ! Spectral or non-spectral albedo schemes

!----------------------------------------------------------------------
! Calculate GBM surface albedo
!----------------------------------------------------------------------

DO band=1,4
  DO i=1,p_field
    land_albedo(i,band) = 0.
  END DO
  DO n=1,ntype
    DO j=1,surft_pts(n)
      l = surft_index(j,n)
      i = land_index(l)
      land_albedo(i,band) = land_albedo(i,band) +                             &
                            frac(l,n)*alb_type(l,n,band)
    END DO
  END DO
END DO

!----------------------------------------------------------------------
! Copy albedos as required for aggregate or distinct tiles
!----------------------------------------------------------------------

IF (l_aggregate) THEN
  DO band=1,4
    DO l=1,land_field
      i = land_index(l)
      alb_surft(l,1,band) = land_albedo(i,band)
    END DO
  END DO
ELSE
  DO band=1,4
    DO n=1,ntype
      DO j=1,surft_pts(n)
        l = surft_index(j,n)
        alb_surft(l,n,band) = alb_type(l,n,band)
      END DO
    END DO
  END DO
END IF

!----------------------------------------------------------------------
! Copy the albedo scaling to obs for output as diagnostic.
!----------------------------------------------------------------------
IF (l_albedo_obs) THEN
  IF (l_aggregate) THEN
! Output the unmodified scaling factor albsfsc
    IF (l_spec_albedo) THEN
! There are 2 scalings, one for VIS and NIR
       DO band=1,2
         DO l=1,land_field
           i = land_index(l)
           albobs_sc(i,1,band) = albsfsc(l,band)
         END DO
       END DO
    ELSE
! Both diagnostics will see the same scaling
       DO l=1,land_field
         i = land_index(l)
         albobs_sc(i,1,1) = albsfsc(l,1)
         albobs_sc(i,1,2) = albsfsc(l,1)
       END DO
    END IF
  ELSE
! Output the modified scaling factor on tiles, albobs_scaling_surft
    IF (l_spec_albedo) THEN
! There are 2 scalings, one for VIS and NIR
       DO band=1,2
         DO n=1,nsurft
           DO l=1,land_field
             i = land_index(l)
             albobs_sc(i,n,band) = albobs_scaling_surft(l,n,band)
           END DO
         END DO
       END DO
    ELSE
! Both diagnostics will see the same scaling
       DO n=1,nsurft
         DO l=1,land_field
           i = land_index(l)
           albobs_sc(i,n,1) = albobs_scaling_surft(l,n,1)
           albobs_sc(i,n,2) = albobs_scaling_surft(l,n,1)
         END DO
       END DO
    END IF
  END IF
ELSE
! no obs scaling used, so the scaling is 1.0
  albobs_sc(:,:,:)=1.0
END IF

!CABLE_LSM:
if(.NOT.first_call) then
        CALL cable_rad_main( mype, timestep_number,                            & 
          cycleno, numcycles, row_length, rows, land_field, ntype,             &
          frac, fland,                                                    &
          surf_down_sw,                                                        &
          cosz,                                                              &
          snow_surft,                                                           &
          soil_alb,         &
          land_albedo,      &
          alb_surft,        &
          land_alb )
endif
first_call=.FALSE. 


IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE tile_albedo_cable
