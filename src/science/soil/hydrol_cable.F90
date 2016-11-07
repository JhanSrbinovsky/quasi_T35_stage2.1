! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE HYDROL-------------------------------------------------

! Description:
!     Surface hydrology module which also updates the
!     sub-surface temperatures. Includes soil water phase
!     changes and the effect of soil water and ice on the
!     thermal and hydraulic characteristics of the soil.
!     This version is for use with MOSES II land surface
!     scheme.

! Documentation : UM Documentation Paper 25

! Subroutine Interface:
MODULE hydrol_mod
  CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='HYDROL_MOD'

CONTAINS
SUBROUTINE hydrol (                                                           &
                   lice_pts,lice_index,soil_pts,soil_index,                   &
                   nsnow,                                                     &
                   npnts,nshyd,b,can_cpy,con_rain,                            &
                   e_canopy,ext,hcap,hcon,ls_rain,                            &
                   con_rainfrac, ls_rainfrac,                                 &
                   satcon,sathh,snowdepth,                                    &
                   snow_soil_htf,surf_ht_flux,timestep,                       &
                   v_sat,v_wilt,                                              &
                   can_wcnt,                                                  &
                   stf_sub_surf_roff,smcl,sthf,sthu,tsoil,tsurf_elev_surft,   &
                   can_wcnt_gb,smc,                                           &
                   snow_melt,                                                 &
                   sub_surf_roff,surf_roff,tot_tfall,                         &
! add new inland basin variable
                   inlandout_atm,l_inland,                                    &
! Additional variables for MOSES II
                   nsurft,surft_pts,surft_index,                              &
                   infil_surft,                                               &
                   melt_surft,tile_frac,                                      &
! Additional variables required for large-scale hydrology:
                   l_top,l_pdm,fexp,gamtot,ti_mean,ti_sig,                    &
                   cs_ch4,cs,                                                 &
                   dun_roff,drain,fsat,fwetl,qbase,qbase_zw,                  &
                   zw,sthzw,a_fsat,c_fsat,a_fwet,c_fwet,                      &
                   resp_s,npp,fch4_wetl,                                      &
                   fch4_wetl_cs,fch4_wetl_npp,fch4_wetl_resps,                &
                   dim_cs1,l_soil_sat_down,l_triffid)

USE jules_hydrology_mod, ONLY :                                               &
 l_wetland_unfrozen                                                           &
,ti_max                                                                       &
,zw_max

USE elev_htc_mod, ONLY : elev_htc

USE jules_soil_mod, ONLY :                                                    &
 dzsoil,                                                                      &
              !  Thicknesses of the soil layers (m)
 dzsoil_elev

USE crop_vars_mod, ONLY : sthu_irr_gb, frac_irr_gb, ext_irr_gb

USE jules_vegetation_mod, ONLY : l_irrig_dmd,l_nitrogen

USE prognostics, ONLY : n_inorg_gb

USE trif_vars_mod, ONLY : n_leach_gb

USE jules_surface_mod, ONLY : l_elev_land_ice

USE c_densty, ONLY :                                                          &
   rho_water  !  density of pure water (kg/m3)

USE jules_surface_mod, ONLY: sorp

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine arguments
!   Scalar arguments with intent(IN) :
INTEGER, INTENT(IN) ::                                                        &
 lice_pts                                                                     &
                     ! IN Number of land ice points.
,npnts                                                                        &
                     ! IN Number of gridpoints.
,nshyd                                                                        &
                     ! IN Number of soil moisture levels.
,soil_pts                                                                     &
                     ! IN Number of soil points.
,nsurft                                                                       &
                     ! IN Number of tiles
,dim_cs1             ! IN Number of soil carbon pools

REAL, INTENT(IN) ::                                                           &
 timestep            ! IN Model timestep (s).

LOGICAL , INTENT(IN) ::                                                       &
 stf_sub_surf_roff                                                            &
                     ! IN Stash flag for sub-surface runoff.
,l_top                                                                        &
             ! IN Flag for TOPMODEL-based hydrology.
,l_pdm                                                                        &
             ! IN Flag for PDM hydrology.
,l_soil_sat_down                                                              &
             ! IN Switch controlling direction of movement of
             !    soil moisture in excess of saturation
,l_triffid   ! IN Switch to use TRIFFID.


!   Array arguments with intent(IN) :
INTEGER, INTENT(IN) ::                                                        &
 lice_index(npnts)                                                            &
                     ! IN Array of land ice points.
,soil_index(npnts)                                                            &
                     ! IN Array of soil points.
,nsnow(npnts,nsurft) ! IN Number of snow layers

REAL, INTENT(IN) ::                                                           &
 b(npnts,nshyd)                                                               &
                     ! IN Clapp-Hornberger exponent.
,can_cpy(npnts,nsurft)                                                        &
                      !IN Canopy/surface capacity of
!                          !    land tiles (kg/m2).
,con_rain(npnts)                                                              &
                     ! IN Convective rain (kg/m2/s).
,e_canopy(npnts,nsurft)                                                       &
!                          ! IN Canopy evaporation from
!                          !    land tiles (kg/m2/s).
,ext(npnts,nshyd)                                                             &
                     ! IN Extraction of water from each soil
!                          !    layer (kg/m2/s).
,hcap(npnts,nshyd)                                                            &
                     ! IN Soil heat capacity (J/K/m3).
,hcon(npnts,0:nshyd)                                                          &
                     ! IN Soil thermal conductivity (W/m/K).
,ls_rain(npnts)                                                               &
                     ! IN Large-scale rain (kg/m2/s).
,con_rainfrac(npnts)                                                          &
                     ! IN Convective rain fraction
,ls_rainfrac(npnts)                                                           &
                     ! IN large scale rain fraction

,satcon(npnts,0:nshyd)                                                        &
                     ! IN Saturated hydraulic conductivity
!                          !    (kg/m2/s).
,sathh(npnts,nshyd)                                                           &
                     ! IN Saturated soil water pressure (m).
,snow_melt(npnts)                                                             &
                     ! IN Snowmelt (kg/m2/s).
,snowdepth(npnts,nsurft)                                                      &
                     ! Snow depth (on ground) (m)
,snow_soil_htf(npnts,nsurft)                                                  &
                     ! IN Tiled snowpack-> soil heat flux
,surf_ht_flux(npnts)                                                          &
                     ! IN Net downward surface heat flux (W/m2)
,v_sat(npnts,nshyd)                                                           &
                     ! IN Volumetric soil moisture
!                          !    concentration at saturation
!                          !    (m3 H2O/m3 soil).
,v_wilt(npnts,nshyd)                                                          &
                     ! IN Volumetric soil moisture
!                          !    concentration below which
!                          !    stomata close (m3 H2O/m3 soil).
,fexp(npnts)                                                                  &
                     ! IN Decay factor in Sat. Conductivity
!                          !    in water table layer.
,gamtot(npnts)                                                                &
                     ! IN Integrated complete Gamma function.
,ti_mean(npnts)                                                               &
                     ! IN Mean topographic index.
,ti_sig(npnts)                                                                &
                     ! IN Standard dev. of topographic index.
,cs_ch4(npnts)                                                                &
                     ! IN Soil carbon used in CH4 wetlands
!                          !   if TRIFFID is switched off
!                          !    (kg C/m2).
,cs(npnts,dim_cs1)                                                            &
                     ! IN Soil carbon (kg C/m2).
!                          !   For RothC (dim_cs1=4), the pools
!                          !   are DPM, RPM, biomass and humus.
,resp_s(npnts,dim_cs1)                                                        &
                     ! IN Soil respiration in pools (kg C/m2/s).
,npp(npnts)                                                                   &
                     ! IN Gridbox mean net primary
!                          !    productivity (kg C/m2/s).
,a_fsat(npnts)                                                                &
                     ! IN Fitting parameter for Fsat in LSH model
,c_fsat(npnts)                                                                &
                     ! IN Fitting parameter for Fsat in LSH model
,a_fwet(npnts)                                                                &
                     ! IN Fitting parameter for Fwet in LSH model
,c_fwet(npnts)
                     ! IN Fitting parameter for Fwet in LSH model


!   Array arguments with intent(INOUT) :

REAL, INTENT(INOUT) ::                                                        &
 can_wcnt(npnts,nsurft)                                                       &
!                          ! INOUT Canopy water content for
!                          !       land tiles (kg/m2).
,smcl(npnts,nshyd)                                                            &
                     ! INOUT Soil moisture content of each
!                          !       layer (kg/m2).
,sthf(npnts,nshyd)                                                            &
                     ! INOUT Frozen soil moisture content of
!                          !       each layer as a fraction of
!                          !       saturation.
,sthu(npnts,nshyd)                                                            &
                     ! INOUT Unfrozen soil moisture content of
!                          !       each layer as a fraction of
!                          !       saturation.
,tsoil(npnts,nshyd)                                                           &
                     ! INOUT Sub-surface temperatures (K).
,tsurf_elev_surft(npnts,nsurft)                                               &
                     ! INOUT Tiled sub-surface temperatures (K).
,fsat(npnts)                                                                  &
                      ! INOUT Surface saturation fraction.
,fwetl(npnts)                                                                 &
                      ! INOUT Wetland fraction.
,zw(npnts)                                                                    &
                      ! INOUT Water table depth (m).
,sthzw(npnts)         ! INOUT soil moist fract. in deep-zw layer.

!   Array arguments with intent(OUT) :
REAL, INTENT(OUT) ::                                                          &
 can_wcnt_gb(npnts)                                                           &
                      ! OUT Gridbox canopy water content (kg/m2).
,smc(npnts)                                                                   &
                      ! OUT Available soil moisture in a layer at
                      !     the surface (kg/m2)
,sub_surf_roff(npnts)                                                         &
                      ! OUT Sub-surface runoff (kg/m2/s).
,surf_roff(npnts)                                                             &
                      ! OUT Surface runoff (kg/m2/s).
,tot_tfall(npnts)                                                             &
                      ! OUT Total throughfall (kg/m2/s).
,dun_roff(npnts)                                                              &
                      ! OUT Dunne part of sfc runoff (kg/m2/s).
,qbase(npnts)                                                                 &
                      ! OUT Base flow (kg/m2/s).
,qbase_zw(npnts)                                                              &
                      ! OUT Base flow from ZW layer (kg/m2/s).
,drain(npnts)                                                                 &
                      ! OUT Drainage out of nshyd'th level (kg/m2/s).
,fch4_wetl(npnts)                                                             &
!                     ! OUT Scaled wetland methane flux.
                      !     (default substrate) for use in
                      !     atmos chemistry model (10^-9 kg C/m2/s)
,fch4_wetl_cs(npnts)                                                          &
                      ! OUT Scaled methane flux (soil carbon
!                          ! substrate) (kg C/m2/s)
,fch4_wetl_npp(npnts)                                                         &
                      ! OUT Scaled methane flux (npp
!                          ! substrate) (kg C/m2/s)
,fch4_wetl_resps(npnts)
                      ! OUT Scaled methane flux (soil respiration
!                          ! substrate) (kg C/m2/s)


! Additional variables for MOSES II
INTEGER, INTENT(IN) ::                                                        &
 surft_pts(nsurft)                                                            &
                     ! IN Number of tile points.
,surft_index(npnts,nsurft)
!                          ! IN Index of tile points.

REAL, INTENT(IN) ::                                                           &
 infil_surft(npnts,nsurft)                                                    &
!                          ! IN Maximum surface infiltration
,melt_surft(npnts,nsurft)                                                     &
!                          ! IN Snowmelt on tiles (kg/m2/s).
,tile_frac(npnts,nsurft)                                                      &
                     ! IN Tile fractions.

! Declare variable for inland basin outflow
,inlandout_atm(npnts)            ! IN TRIP INLAND BASIN
!                       OUTFLOW FOR LAND POINTS ONLY,kg/m2/s=mm
LOGICAL, INTENT(IN) ::                                                        &
 l_inland                   ! IN True if re-routing inland
                            !   basin flow to soil moisture

! Local scalars:
INTEGER                                                                       &
 i,j                                                                          &
                      ! WORK Loop counters.
,n                    ! WORK Tile loop counter.

! Local arrays:

REAL                                                                          &
 dsmc_dt(npnts)                                                               &
                      ! WORK Rate of change of soil moisture
!                           !      due to water falling onto the
!                           !      surface after surface runoff
!                           !      (kg/m2/s).
,w_flux(npnts,0:nshyd)                                                        &
                      ! WORK Fluxes of water between layers
!                           !      (kg/m2/s).
,ksz(npnts,0:nshyd)                                                           &
                      ! WORK Saturated hydraulic
!                           !      conductivity in layer (kg/m2/s).
,qbase_unfr(npnts)                                                            &
                      ! WORK Base flow in unfrozon soil (kg/m2/s).
!                           !      in inundation calculation (m).
,qbase_l(npnts,nshyd+1)                                                       &
!                           ! WORK Base flow from each level (kg/m2/s).
,qbase_l_unfr(npnts,nshyd+1)                                                  &
!                     ! WORK As qbase_l but for unfrozon soil (kg/m2/s).
,top_crit(npnts)                                                              &
                      ! WORK Critical TI when ZW <=0.0
,dumtop_crit(npnts)                                                           &
                      ! WORK Dummy for top_crit
,dumsthf(npnts,nshyd)                                                         &
                      ! WORK Dummy Frozen soil moisture content of
!                           ! each layer as a fraction of
!                           ! saturation (always set to 0).
,zdepth(0:nshyd)                                                              &
                      ! WORK Lower soil layer boundary depth (m).
,tsoil_d(npnts)                                                               &
                      ! WORK Soil temperature in the top metre
,zw_inund(npnts)                                                              &
                      ! WORK Water table depth used
,wutot(npnts)                                                                 &
!                     ! WORK Ratio of unfrozen to total soil
!                           ! moisture at ZW.
,dumwutot(npnts)
!                     ! WORK Dummy of wutot - always set to 1.0

! WORK variables required for irrigation code
REAL                                                                          &
w_flux_irr(npnts,0:nshyd)                                                     &
                      ! WORK The fluxes of water between layers
!                          !     in irrigated fraction (kg/m2/s).
,w_flux_nir(npnts,0:nshyd)                                                    &
                      ! WORK The fluxes of water between layers
!                          !     in non-irrigated fraction (kg/m2/s).
,smcl_irr(npnts,nshyd)                                                        &
                      ! WORK Total soil moisture contents of each
!                          !       layer in irrigated fraction (kg/m2).
,smcl_nir(npnts,nshyd)                                                        &
                      ! WORK Total soil moisture contents of each
!                          !       layer in non-irrigated fraction (kg/m2).
,sthu_nir(npnts,nshyd)                                                        &
                      ! WORK Unfrozen soil moisture content of
!                          !    each layer as a fraction of
!                          !    saturation in irrigated fraction.
,ext_nir(npnts,nshyd)                                                         &
                      ! WORK Extraction of water from each soil
!                          !    layer in non-irrigated fraction (kg/m2/s).
,smclsat(npnts,nshyd)                                                         &
                      ! WORK The saturation moisture content of
!                           !     each layer (kg/m2).
,smclzw(npnts)                                                                &
                      ! WORK moisture content in deep layer.
,smclsatzw(npnts)
                      ! WORK moisture content in deep layer

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='HYDROL'

! End of header--------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialise w_flux variables that are used in irrigation code
  w_flux(:,:) = 0.0
  w_flux_irr(:,:) = 0.0 ! to prevent random values reported over areas
                      ! that are not included as soil points (i.e. ice points)

!----------------------------------------------------------------------
! Set up variables required for LSH scheme:
!----------------------------------------------------------------------
zdepth(:)=0.0

DO n=1,nshyd
   zdepth(n)=zdepth(n-1)+dzsoil(n)
END DO
!-----------------------------------------------------------------------
! Calculate throughfall and surface runoff, and update the canopy water
! content
!-----------------------------------------------------------------------
! DEPENDS ON: surf_hyd
CALL surf_hyd (npnts,nsurft,surft_pts,surft_index,                            &
               can_cpy,e_canopy,tile_frac,infil_surft,con_rain,               &
               ls_rain, con_rainfrac, ls_rainfrac, melt_surft,                &
               snow_melt,timestep,                                            &
               can_wcnt,can_wcnt_gb,dsmc_dt,                                  &
               l_top,l_pdm,nshyd,soil_pts,soil_index,                         &
               surf_roff,tot_tfall,                                           &
               dun_roff,fsat,v_sat,sthu,sthf)

!-----------------------------------------------------------------------
! Specify the reduction of hydraulic conductivity with depth:
! Initialiase base flow to zero:
!-----------------------------------------------------------------------

DO n=0,nshyd
!CDIR NODEP
  DO j=1,soil_pts
    i=soil_index(j)
    ksz(i,n)=satcon(i,n)
  END DO
END DO

DO n=1,nshyd
!CDIR NODEP
  DO j=1,soil_pts
    smclsat(soil_index(j),n)=rho_water*dzsoil(n)*v_sat(soil_index(j),n)
    qbase_l(soil_index(j),n)=0.0
    qbase_l_unfr(soil_index(j),n)=0.0
    dumsthf(soil_index(j),n)=0.0
  END DO
END DO

DO i=1,npnts
  qbase(i)=0.0
  qbase_zw(i)=0.0
  wutot(i)=0.0
  drain(i)=0.0
  qbase_unfr(i)=0.0
  dumwutot(i)=0.0
  zw_inund(i)=0.0
END DO

IF(l_top)THEN
  IF (soil_pts /= 0) THEN
! DEPENDS ON: calc_baseflow_jules
    CALL calc_baseflow_jules(                                                 &
        soil_pts,soil_index,npnts,nshyd                                       &
       ,zdepth,ksz                                                            &
       ,b,fexp,ti_mean,zw,sthf                                                &
       ,dumtop_crit,qbase,qbase_l                                             &
        )

    IF(L_wetland_unfrozen)THEN
! DEPENDS ON: calc_zw_inund
      CALL calc_zw_inund(npnts,nshyd,soil_pts,soil_index,zdepth               &
       ,b,sathh,smclsat,smcl,sthu,sthzw                                       &
       ,zw,zw_inund,wutot                                                     &
        )

! Now call again to get the unfrozen equivalents to calculate fsat and fwet:
      CALL calc_baseflow_jules(                                               &
        soil_pts,soil_index,npnts,nshyd                                       &
       ,zdepth,ksz                                                            &
       ,b,fexp,ti_mean,zw_inund,dumsthf                                       &
       ,top_crit,qbase_unfr,qbase_l_unfr                                      &
        )
    ELSE
      top_crit(:)=dumtop_crit(:)
    END IF

  END IF
END IF


IF(l_inland)THEN

 DO i=1,npnts

! Add inland basin outflow to change in soil moisture store

   dsmc_dt(i)=dsmc_dt(i)+inlandout_atm(i)

 END DO
END IF


!-----------------------------------------------------------------------
! Update the layer soil moisture contents and calculate the
! gravitational drainage.
!-----------------------------------------------------------------------
IF (soil_pts /= 0) THEN

IF (l_irrig_dmd) THEN

!-----------------------------------------------------------------------
! If l_irrig_dmd = TRUE, call soil_hyd separately for irrigated and
! non-irrigated fraction
! afterwards, call soil_hyd ONLY to update water table/drainge with
! gridbox total w_flux, smcl
!-----------------------------------------------------------------------

! split into irrigated / non-irrigated fraction
  DO n=1,nshyd
    DO j=1,soil_pts
      i=soil_index(j)
      smclsat(i,n) = rho_water*dzsoil(n)*v_sat(i,n)
! hadrd - gridbox sthu is assumed to be combination of sthu of
! non-irrigated fraction
! and sthu_irr_gb, i.e.
! sthu = frac_irr_gb*sthu_irr_gb + (1-frac_irr_gb)*sthu_nir
      IF ( frac_irr_gb(i) < 1 ) THEN
        sthu_nir(i,n) = ( sthu(i,n) - frac_irr_gb(i)*sthu_irr_gb(i,n) )       &
                        / (1.0 - frac_irr_gb(i))
        ext_nir(i,n) = ( ext(i,n) - frac_irr_gb(i)*ext_irr_gb(i,n) )          &
                        / (1.0 - frac_irr_gb(i))
      ELSE
        sthu_nir(i,n) = sthu(i,n)
        ext_nir(i,n) = ext(i,n)
      END IF
      smcl_irr(i,n) = smcl(i,n) + ( sthu_irr_gb(i,n)-sthu(i,n) ) * smclsat(i,n)
      smcl_nir(i,n) = smcl(i,n) + ( sthu_nir(i,n)-sthu(i,n) )    * smclsat(i,n)
    END DO
  END DO

! first call soil_hyd for non-irrigated fraction
! DEPENDS ON: soil_hyd
  CALL soil_hyd (npnts,nshyd,soil_pts,soil_index,b,dzsoil,                    &
                ext_nir,dsmc_dt,satcon,ksz,sathh,timestep,v_sat,              &
                sub_surf_roff,smcl_nir,sthu_nir,surf_roff,w_flux_nir,         &
                stf_sub_surf_roff,                                            &
                zw,sthzw,zdepth,qbase,qbase_l,                                &
                dun_roff,drain,l_top,l_soil_sat_down,                         &
                smclzw,smclsatzw,smclsat)

! next call soil_hyd for irrigated fraction
! DEPENDS ON: soil_hyd
  CALL soil_hyd (npnts,nshyd,soil_pts,soil_index,b,dzsoil,                    &
                ext_irr_gb,dsmc_dt,satcon,ksz,sathh,timestep,v_sat,           &
                sub_surf_roff,smcl_irr,sthu_irr_gb,surf_roff,w_flux_irr,      &
                stf_sub_surf_roff,                                            &
                zw,sthzw,zdepth,qbase,qbase_l,                                &
                dun_roff,drain,l_top,l_soil_sat_down,                         &
                smclzw,smclsatzw,smclsat)

! re-calculate total grid box soil moisture
! hadrd - perhaps this could be done more efficiently with WHERE
  DO n=0,nshyd
    DO j=1,soil_pts
      i=soil_index(j)

      ! Ensure sensible values if irrigation fraction is very small:
      IF ( frac_irr_gb(i) <= EPSILON(1.0) ) THEN
        w_flux_irr(i,n) = 0.0
        w_flux(i,n) = w_flux_nir(i,n)
        IF ( n >= 1 ) THEN
          sthu_irr_gb(i,n) = sthu_nir(i,n)
          smcl(i,n) = smcl_nir(i,n)
          sthu(i,n) = sthu_nir(i,n)
        END IF
      ELSE
        w_flux(i,n) = frac_irr_gb(i) * w_flux_irr(i,n)                        &
                    + ( 1.0 - frac_irr_gb(i) ) * w_flux_nir(i,n)
        IF ( n >= 1 ) THEN
          smcl(i,n) = frac_irr_gb(i) * smcl_irr(i,n)                          &
                      + ( 1.0 - frac_irr_gb(i) ) * smcl_nir(i,n)
          sthu(i,n) = frac_irr_gb(i) * sthu_irr_gb(i,n)                       &
                      + ( 1.0 - frac_irr_gb(i) ) * sthu_nir(i,n)
        END IF
      END IF
    END DO
  END DO


! DEPENDS ON: soil_hyd_update
  CALL soil_hyd_update (npnts,nshyd,soil_pts,soil_index,                      &
                    dzsoil,v_sat,zdepth,smclzw,sthzw,smclsat,smclsatzw)

! DEPENDS ON: soil_hyd_wt
  CALL soil_hyd_wt (npnts,nshyd,soil_pts,soil_index,                          &
              b,dsmc_dt,sathh,timestep,v_sat,                                 &
              sub_surf_roff,smcl,surf_roff,w_flux,                            &
              stf_sub_surf_roff,                                              &
              zw,sthzw,qbase,qbase_l,                                         &
              drain,l_top,smclzw,smclsatzw,smclsat)

ELSE ! if .NOT. l_irrig_dmd

! DEPENDS ON: soil_hyd
  CALL soil_hyd (npnts,nshyd,soil_pts,soil_index,b,dzsoil,                    &
                ext,dsmc_dt,satcon,ksz,sathh,timestep,v_sat,                  &
                sub_surf_roff,smcl,sthu,surf_roff,w_flux,                     &
                stf_sub_surf_roff,                                            &
                zw,sthzw,zdepth,qbase,qbase_l,                                &
                dun_roff,drain,l_top,l_soil_sat_down,                         &
                smclzw,smclsatzw,smclsat)

! DEPENDS ON: soil_hyd_wt
  CALL soil_hyd_wt (npnts,nshyd,soil_pts,soil_index,                          &
              b,dsmc_dt,sathh,timestep,v_sat,                                 &
              sub_surf_roff,smcl,surf_roff,w_flux,                            &
              stf_sub_surf_roff,                                              &
              zw,sthzw,qbase,qbase_l,                                         &
              drain,l_top,smclzw,smclsatzw,smclsat)


END IF !l_irrig_dmd

!-----------------------------------------------------------------------
! Calculate surface saturation and wetland fractions:
!-----------------------------------------------------------------------
  IF(l_top)THEN

    DO i=1,npnts
      fsat(i)=0.0
      fwetl(i)=0.0
! Zero soil porosity over land ice:
      IF(v_sat(i,nshyd) <= 0.0) zw(i)=zw_max
    END DO

    DO j=1,soil_pts
      i=soil_index(j)
      qbase_zw(i)=qbase_l(i,nshyd+1)

!Now use fit for fsat and fwet:
      IF(L_wetland_unfrozen)THEN
        fsat(i)=wutot(i)*a_fsat(i)*EXP(-c_fsat(i)*top_crit(i))
        fwetl(i)=wutot(i)*a_fwet(i)*EXP(-c_fwet(i)*top_crit(i))
      ELSE
        fsat(i)=a_fsat(i)*EXP(-c_fsat(i)*top_crit(i))
        fwetl(i)=a_fwet(i)*EXP(-c_fwet(i)*top_crit(i))
      END IF

      IF(top_crit(i) >= ti_max)THEN
        fsat(i)=0.0
        fwetl(i)=0.0
      END IF
    END DO
  END IF  !  l_top

ELSE ! soil pts

!---------------------------------------------------------------------
! If required by STASH flag and there are no soil points,
! set sub-surface runoff to zero.
!---------------------------------------------------------------------

  IF(stf_sub_surf_roff) THEN
    DO i=1,npnts
      sub_surf_roff(i)=0.0
    END DO
  END IF

END IF  !  soil_pts

!-----------------------------------------------------------------------
! Update the soil temperatures and the frozen moisture fractions
!-----------------------------------------------------------------------
IF (soil_pts /= 0) THEN
! DEPENDS ON: soil_htc
  CALL soil_htc (npnts,nshyd,nsurft,soil_pts,soil_index,                      &
                 surft_pts,surft_index,nsnow,                                 &
                 b,dzsoil,tile_frac,hcap,hcon,                                &
                 sathh,surf_ht_flux,timestep,v_sat,                           &
                 w_flux,smcl,snowdepth,sthu,sthf,tsoil)
END IF

!-----------------------------------------------------------------------
! Update the sub-surface temperatures for land ice
!-----------------------------------------------------------------------
IF (lice_pts /= 0) THEN

  IF (.NOT. l_elev_land_ice) THEN
! DEPENDS ON: ice_htc
  CALL ice_htc (npnts,nshyd,lice_pts,lice_index,dzsoil,                       &
                surf_ht_flux,timestep,                                        &
                tsoil)
  ELSE
    CALL elev_htc (npnts,lice_pts,lice_index,nsurft,                          &
                   dzsoil_elev,snow_soil_htf,timestep,                        &
                   tsurf_elev_surft)
  END IF

END IF

!-----------------------------------------------------------------------
! Diagnose the available soil moisture in a layer at the surface.
!-----------------------------------------------------------------------
! DEPENDS ON: soilmc
CALL soilmc ( npnts,nshyd,soil_pts,soil_index,                                &
              dzsoil,sthu,v_sat,v_wilt,smc )

!-----------------------------------------------------------------------
! Calculate mean soil temperature and scaled CH4 flux:
!-----------------------------------------------------------------------

DO i=1,npnts
  fch4_wetl(i)=0.0
  fch4_wetl_cs(i)=0.0
  fch4_wetl_npp(i)=0.0
  fch4_wetl_resps(i)=0.0
END DO
IF(l_top)THEN
  IF (soil_pts /= 0) THEN
! DEPENDS ON: soilt
    CALL soilt(npnts,nshyd,soil_pts,soil_index                                &
             ,dzsoil,tsoil,tsoil_d)
! DEPENDS ON: ch4_wetl
    CALL ch4_wetl(npnts,soil_pts,dim_cs1,soil_index,l_triffid                 &
      ,tsoil_d,cs_ch4,cs,resp_s,npp,fwetl,fch4_wetl                           &
      ,fch4_wetl_cs,fch4_wetl_npp,fch4_wetl_resps                             &
       )
  END IF
END IF

!-----------------------------------------------------------------------
! Calculate Nitrogen Leaching
!-----------------------------------------------------------------------

IF (l_nitrogen) THEN
  n_leach_gb(:)=0.
  DO i=1,npnts
    n_leach_gb(i)=0.
    IF (SUM(smcl(i,1:3)) > EPSILON(0.)) THEN
      n_leach_gb(i) = sub_surf_roff(i) * (n_inorg_gb(i) / sorp) / SUM(smcl(i,1:3))
      ! Layers 1-3 (top 1m)
    END IF
    n_leach_gb(i) = MIN(n_leach_gb(i),(n_inorg_gb(i) / sorp) / timestep)
             ! H20 flux                ! N concentration Kg [N]/kg [H20]
    N_inorg_gb(i) = N_inorg_gb(i) - n_leach_gb(i) * timestep
  ENDDO
ENDIF

!-----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE hydrol
END MODULE hydrol_mod
