! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Description:
! Subroutine to calculate gridbox mean values of surface conductance
! and carbon fluxes. Also returns net primary productivity, leaf
! turnover and wood respiration of each plant functional type for
! driving TRIFFID.

! *********************************************************************
SUBROUTINE physiol (land_pts,land_index                                       &
,                   nshyd,nsurft,surft_pts,surft_index,l_aggregate            &
,                   dim_cs1, dim_cs2                                          &
,                   co2,co2_3d,co2_dim_len                                    &
,                   co2_dim_row,l_co2_interactive                             &
,                   l_triffid,l_q10                                           &
,                   can_model,cs,frac,fland,ht,ipar,lai,pstar,q1              &
,                   sthu,sthf,tsoil,tstar_surft                               &
,                   v_crit,v_sat,v_wilt,wind,z0_surft,z1,o3                   &
,                   canhc_surft,vfrac_surft,emis_surft,emis_soil              &
,                   flake,g_leaf,gs,gs_surft                                  &
,                   gpp,gpp_pft,npp,npp_pft,resp_p,resp_p_pft                 &
,                   resp_s,resp_l_pft,resp_r_pft,resp_w_pft,n_leaf            &
,                   n_root,n_stem,lai_bal,smct,wt_ext_surft,fsmc,wt_ext       &
,                   albsoil,cos_zenith_angle                                  &
,                   can_rad_mod,ilayers,flux_o3_pft,fo3_pft                   &
,                   et_stom_gb,et_stom_pft)

USE sf_stom_mod, ONLY: sf_stom

USE smc_ext_mod, ONLY: smc_ext

USE atm_fields_bounds_mod
USE theta_field_sizes, ONLY : t_i_length, t_j_length

USE c_densty

USE ancil_info, ONLY: l_soil_point

USE jules_surface_types_mod, ONLY :                                           &
   lake, npft, nnpft, ntype, soil, urban_canyon, urban_roof

USE nvegparm

USE jules_soil_mod, ONLY : dzsoil
USE jules_surface_mod, ONLY : diff_frac
USE jules_vegetation_mod, ONLY : l_crop, l_irrig_dmd, l_soil_resp_lev2

USE pftparm

USE urban_param, ONLY : emisr_gb, emisw_gb, hwr_gb, emiss,                    &
   diffus_road, diffus_wall, diffus_roof, cap_road, cap_wall,                 &
   cap_roof, dz_roof_p, omega_day

USE jules_radiation_mod, ONLY : l_spec_albedo, l_albedo_obs,                  &
                                l_spec_alb_bs

USE jules_mod, ONLY : albobs_scaling_surft

USE switches_urban, ONLY :                                                    &
   l_urban2t, l_moruses_emissivity, l_moruses_storage,                        &
   l_moruses_storage_thin

USE bvoc_vars,                ONLY:                                           &
  isoprene_gb, isoprene_pft, terpene_gb , terpene_pft,                        &
  methanol_gb, methanol_pft, acetone_gb, acetone_pft

USE trif_vars_mod, ONLY: fapar_diag_pft

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE cropparm, ONLY : rt_dir, cfrac_r
USE crop_vars_mod, ONLY : rootc_cpft, gs_irr_surft, smc_irr_gb,               &
   sthu_irr_gb, frac_irr_gb, frac_irr_surft, wt_ext_irr_surft, gc_irr_surft

USE jules_print_mgr, ONLY :                                                   &
    jules_message,                                                            &
    jules_print,                                                              &
    PrNorm

USE ereport_mod, ONLY: ereport

IMPLICIT NONE

INTEGER                                                                       &
 land_pts                                                                     &
                            ! IN Number of land points to be
!                                 !    processed.
,land_index(land_pts)                                                         &
                            ! IN Index of land points on the
!                                 !    P-grid.
,co2_dim_len                                                                  &
                            ! IN Length of a CO2 field row.
,co2_dim_row                                                                  &
                            ! IN Number of CO2 field rows.
,nshyd                                                                        &
                            ! IN Number of soil moisture
!                                 !    levels.
,nsurft                                                                       &
                            ! IN Number of surface tiles.
,surft_pts(ntype)                                                             &
                            ! IN Number of land points which
!                                 !    include the nth surface type.
,surft_index(land_pts,ntype)                                                  &
                            ! IN Indices of land points which
!                                 !    include the nth surface type.
,can_model                                                                    &
                            ! IN Swith for thermal vegetation
!                                 !    canopy
,dim_cs1, dim_cs2           ! soil carbon dimensions

LOGICAL                                                                       &
        l_co2_interactive                                                     &
                            ! switch for 3D CO2 field
,       l_aggregate                                                           &
                            ! IN Logical to set aggregate
                            !    surface scheme
,       l_triffid                                                             &
                            ! TRUE if using TRIFFID
,       l_q10               ! TRUE if using Q10 for soil resp

LOGICAL ::                                                                    &
   firstcall = .TRUE.

INTEGER                                                                       &
 can_rad_mod                                                                  &
!                                !Switch for canopy radiation model
,ilayers                                                                      &
!                                !No of layers in canopy radiation model
,albpft_call                     ! Flag for albpft, scaling to obs

REAL                                                                          &
 alb_type_dummy(land_pts,ntype,4)                                             &
!                                 ! WORK Dummy argument for albedo
!                                 ! subroutine
,fapar_dir(land_pts,npft,ilayers)                                             &
!                                 ! WORK Profile of absorbed PAR -
!                                 ! Direct beam
,fapar_dif(land_pts,npft,ilayers)                                             &
!                                 ! WORK Profile of absorbed PAR -
!                                 ! Diffuse beam
,faparv(land_pts,ilayers)                                                     &
!                                 ! WORK Profile of absorbed PAR.
,fapar_dif2dif(land_pts,npft,ilayers)                                         &
!                                 ! WORK Profile of absorbed PAR
!                                 ! - Diffuse (fraction/LAI)
,fapar_dir2dif(land_pts,npft,ilayers)                                         &
!                                 ! WORK Profile of scattered PAR
!                                 ! -direct to diffuse (fraction/LAI)
,fapar_dir2dir(land_pts,npft,ilayers)                                         &
!                                 ! WORK Profile of absorbed only
!                                 ! direct PAR (fraction/LAI)
,fapar_shd(land_pts,ilayers)                                                  &
!                                 ! WORK fraction of total absorbed PAR
!                                 ! by shaded leaves
,fapar_sun(land_pts,ilayers)                                                  &
!                                 ! WORK fraction of total absorbed PAR
!                                 ! by sunlit leaves
,fsun(land_pts,npft,ilayers)
!                                 ! WORK fraction of leaves that are
!                                 ! sunlit



REAL                                                                          &
 co2                                                                          &
                            ! IN Atmospheric CO2 concentration
,co2_3d(co2_dim_len,co2_dim_row)                                              &
!                                 ! IN 3D atmos CO2 concentration
!                                 !    (kg CO2/kg air).
,cs(land_pts,dim_cs1)                                                         &
                           ! IN Soil carbon (kg C/m2).
,veg_frac(dim_cs2)                                                            &
                           ! WORK vegetated fraction of gridbox
,frac(land_pts,ntype)                                                         &
                            ! IN Surface type fractions.
,fland(land_pts)                                                              &
                            ! IN Land fraction on land tiles
,ht(land_pts,npft)                                                            &
                            ! IN Canopy height (m).
,ipar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                    &
                            ! IN Incident PAR (W/m2).
,lai(land_pts,npft)                                                           &
                            ! IN Leaf area index.
,pstar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                            ! IN Surface pressure (Pa).
,q1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                      &
                            ! IN Specific humidity at level 1
!                                 !    (kg H2O/kg air).
,sthu(land_pts,nshyd)                                                         &
                            ! IN Soil moisture content in each
!                                 !    layer as a fraction of saturation
,sthf(land_pts,nshyd)                                                         &
                            ! IN Frozen Soil moisture content in each
!                                 !    layer as a fraction of saturation
,tsoil(land_pts,nshyd)                                                        &
                            ! IN Soil temperature (K).
,tstar_surft(land_pts,nsurft)                                                 &
                            ! IN Tile surface temperatures (K).
,v_crit(land_pts,nshyd)                                                       &
                            ! IN Volumetric soil moisture
!                                 !    concentration at field capacity
!                                 !    (m3 H2O/m3 soil).
,v_sat(land_pts,nshyd)                                                        &
                            ! IN Volumetric soil moisture
!                                 !    concentration at saturation
!                                 !    (m3 H2O/m3 soil).
,v_wilt(land_pts,nshyd)                                                       &
                            ! IN Volumetric soil moisture
!                                 !    concentration below which
!                                 !    stomata close (m3 H2O/m3 soil).

,wind(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                    &
                            ! IN Windspeed (m/s).
,z0_surft(land_pts,nsurft)                                                    &
                            ! IN Tile roughness lengths (m).
,z1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                      &
                            ! IN Windspeed reference height(m).
,o3(land_pts)                                                                 &
                            ! IN Surface ozone concentration (ppb).

,gs(land_pts)                                                                 &
                            ! INOUT Gridbox mean surface
!                                 !       conductance (m/s).
,wt_ext(land_pts,nshyd)                                                       &
                            ! OUT Gridbox-mean WT_EXT.
!                                       NB This is only non-zero if
!                                          l_aggregate=TRUE.
,albsoil(land_pts)                                                            &
!                                  ! Soil albedo.
,cos_zenith_angle(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
!                                  ! Cosine of the zenith angle

REAL                                                                          &
 canhc_surft(land_pts,nsurft)                                                 &
                            ! OUT Areal heat capacity of canopy
!                                 !     for land tiles (J/K/m2).
,flake(land_pts,nsurft)                                                       &
                            ! OUT Lake fraction.
,g_leaf(land_pts,npft)                                                        &
                            ! OUT Leaf turnover rate (/360days).
,gs_surft(land_pts,nsurft)                                                    &
                            ! OUT Surface conductance for
!                                 !     land tiles (m/s).
,gpp(land_pts)                                                                &
                            ! OUT Gridbox mean gross primary
!                                 !     productivity (kg C/m2/s).
,gpp_pft(land_pts,npft)                                                       &
                            ! OUT Gross primary productivity
!                                 !     (kg C/m2/s).
,npp(land_pts)                                                                &
                            ! OUT Gridbox mean net primary
!                                 !     productivity (kg C/m2/s).
,npp_pft(land_pts,npft)                                                       &
                            ! OUT Net primary productivity
!                                 !     (kg C/m2/s).
,resp_p(land_pts)                                                             &
                            ! OUT Gridbox mean plant respiration
!                                 !     (kg C/m2/s).
,resp_p_pft(land_pts,npft)                                                    &
                            ! OUT Plant respiration (kg C/m2/s).
,resp_s(land_pts,dim_cs1)                                                     &
                           ! OUT Soil respiration (kg C/m2/s).
,resp_l_pft(land_pts,npft)                                                    &
                            ! OUT Leaf maintenance respiration
!                                 !     (kg C/m2/s).
,resp_r_pft(land_pts,npft)                                                    &
                            ! OUT Root maintenance respiration
!                                 !     (kg C/m2/s).
,resp_w_pft(land_pts,npft)                                                    &
                            ! OUT Wood maintenance respiration
!                                 !     (kg C/m2/s).
,n_leaf(land_pts,npft)                                                        &
                            ! OUT Leaf N content scaled by LAI
!                                 !     (kg N/m2).
,n_root(land_pts,npft)                                                        &
                            ! OUT Root N content scaled by LAI_bal
!                                 !     (kg N/m2).
,n_stem(land_pts,npft)                                                        &
                            ! OUT Stem N content scaled by LAI_bal
!                                 !     (kg N/m2).
,lai_bal(land_pts,npft)                                                       &
                            ! OUT LAI_bal
,smct(land_pts)                                                               &
                            ! OUT Available moisture in the
!                                 !     soil profile (mm).
,vfrac_surft(land_pts,nsurft)                                                 &
                            ! OUT Fractional canopy coverage for
!                                 !     land tiles.
,emis_surft(land_pts,nsurft)                                                  &
                            ! OUT Emissivity for land tiles
,emis_soil(land_pts)                                                          &
                            ! OUT Emissivity of underlying soil
,wt_ext_surft(land_pts,nshyd,nsurft)                                          &
!                                 ! OUT Fraction of evapotranspiration
!                                 !     which is extracted from each
!                                 !     soil layer by each tile.
,fsmc(land_pts,npft)                                                          &
                            ! OUT Moisture availability factor.
,flux_o3_pft(land_pts,npft)                                                   &
                            ! OUT Flux of O3 to stomata (nmol O3/m2/s).
,fo3_pft(land_pts,npft)                                                       &
                            ! OUT Ozone exposure factor.
,et_stom_gb(land_pts)                                                         &
                            ! OUT Transpiration through stom (W/m2)
,et_stom_pft(land_pts,npft)

REAL                                                                          &
fsmc_irr(land_pts,npft)                                                       &
!                            ! WORK Soil moisture availability
!                                 !     factor over irrigated fraction.
,sthu_nir(land_pts,nshyd)                                                     &
,sthu_surft(land_pts,nshyd)                                                   &
,wt_ext_irr(land_pts,nshyd)                                                   &
!                            ! WORK Fraction of transpiration over
!                                 !      irrigated fraction extracted
!                                 !      from each soil layer.
,wt_ext_irr_type(land_pts,nshyd,ntype)                                        &
!                            ! WORK Fraction of transpiration over
!                                 !     irrigated area extracted from
!                                 !     each soil layer (kg/m2/s).
,gs_irr_type(land_pts,ntype)                                                  &
!                            ! WORK Conductance for irrigated surface types

,gsoil_irr(land_pts)         ! WORK Bare soil conductance over
!                                 !      irrigated fraction.

INTEGER                                                                       &
open_index(land_pts)                                                          &
!                            ! WORK Index of land points
!                                 !      with open stomata.
,open_pts                    ! WORK Number of land points
!                                 !      with open stomata.


!  External routines called :-
EXTERNAL root_frac,raero,soil_evap,                                   &
 leaf_lit,cancap,microbe


REAL                                                                          &
 canhc(land_pts)                                                              &
                            ! WORK Canopy heat capacity (J/K/m2).
,ch_type(land_pts,ntype)                                                      &
                            ! WORK CANHC for surface types.
,root_param(land_pts,npft)                                                    &
                            ! WORK Parameter for F_ROOT
,f_root(land_pts,nshyd)                                                       &
                            ! WORK Fraction of roots in each soil
!                                 !      layer.
,gsoil(land_pts)                                                              &
                            ! WORK Bare soil conductance.
,gs_type(land_pts,ntype)                                                      &
                            ! WORK Conductance for surface types.
,pstar_land(land_pts)                                                         &
                            ! WORK Surface pressure (Pa).
,ipar_land(land_pts)                                                          &
                            ! WORK Incident PAR (W/m2).
,q1_land(land_pts)                                                            &
                            ! WORK ecific humidity at level 1
,ra(land_pts)                                                                 &
                            ! WORK Aerodynamic resistance (s/m).
,rib(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                     &
                            ! WORK Bulk Richardson Number.
,tstar(land_pts)                                                              &
                            ! WORK Surface temperature (K).
,vfrac(land_pts)                                                              &
                            ! WORK Fractional canopy coverage.
,vf_type(land_pts,ntype)                                                      &
                            ! WORK VFRAC for surface types.
,wt_ext_type(land_pts,nshyd,ntype)                                            &
!                                 ! WORK WT_EXT for surface types.
,fsoil(land_pts,npft)                                                         &
                            ! WORK Fraction of ground below canopy
!                                 !      contributing to evaporation.
,fsoil_tot(land_pts)                                                          &
                            ! WORK Total fraction of soil
!                                 !      contributing to evaporation
!                                 !      ( = bare soil fraction
!                                 !    + fraction of ground below canopy
,z0(land_pts)                                                                 &
                            ! WORK Roughness length (m).
,emisc(land_pts)                                                              &
!                                 ! WORK Emissivity of canyon
,dz_wall                                                                      &
!                                 ! WORK Damping depth: wall
,dz_road                                                                      &
!                                 ! WORK Damping depth: road
,dz_roof                                                                      &
!                                 ! WORK Damping depth: roof
,dz_roof_c                  ! WORK Damping depth: roof

REAL :: albudir(land_pts,2,npft)
!                                 ! Direct albedo of underlying
!                                 ! surface
REAL :: albudif(land_pts,2,npft)
!                                 ! Diffuse albedo of underlying
!                                 ! surface
!
INTEGER                                                                       &
 i,j,k,l,m,n                ! Loop indices

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PHYSIOL'

CHARACTER (LEN=80) :: errmsg

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------
! Initialisations
!-----------------------------------------------------------------------

f_root(:,:)           = 0.0
gpp(:)                = 0.0
npp(:)                = 0.0
resp_p(:)             = 0.0
smct(:)               = 0.0
ra(:)                 = 0.0
canhc(:)              = 0.0
vfrac(:)              = 0.0
veg_frac(:)           = 0.0
wt_ext(:,:)           = 0.0
rib(:,:)              = 0.0
g_leaf(:,:)           = 0.0
gpp_pft(:,:)           = 0.0
npp_pft(:,:)           = 0.0
resp_p_pft(:,:)        = 0.0
resp_w_pft(:,:)        = 0.0
resp_l_pft(:,:)        = 0.0
resp_r_pft(:,:)        = 0.0
fsmc(:,:)             = 0.0
resp_s(:,:)           = 0.0
faparv(:,:)           = 0.0
wt_ext_type(:,:,:)    = 0.0
alb_type_dummy(:,:,:) = 0.0
fapar_dir(:,:,:)      = 0.0
fapar_dif(:,:,:)      = 0.0

isoprene_gb(:)         = 0.0
isoprene_pft(:,:)      = 0.0
terpene_gb(:)          = 0.0
terpene_pft(:,:)       = 0.0
methanol_gb(:)         = 0.0
methanol_pft(:,:)      = 0.0
acetone_gb(:)          = 0.0
acetone_pft(:,:)       = 0.0
et_stom_gb(:)          = 0.0
et_stom_pft(:,:)       = 0.0

root_param(:,:)       = 0.0

! irrigation variables
  smc_irr_gb(:) = 0.0
  wt_ext_irr(:,:)=0.0
  wt_ext_irr_type(:,:,:)=0.0
  fsmc_irr(:,:)=0.0

!Compress some variables to land points only for st_stom
!$OMP PARALLEL IF(land_pts > 1) DEFAULT(NONE) PRIVATE(i, j, l, n)             &
!$OMP SHARED(land_pts, land_index, t_i_length, pstar_land, pstar, ipar_land,  &
!$OMP        ipar, q1_land, q1, ntype, gs_type, gs, l_irrig_dmd, gs_irr_type, &
!$OMP        fsoil_tot, frac, gsoil, v_crit, soil, npft, sthu, v_sat,         &
!$OMP        gsoil_irr, gs_nvg, sthu_irr_gb)
!$OMP DO SCHEDULE(STATIC)
DO l=1,land_pts
  j             = (land_index(l)-1)/t_i_length + 1
  i             = land_index(l) - (j-1)*t_i_length
  pstar_land(l) = pstar(i,j)
  ipar_land(l)  = ipar(i,j)
  q1_land(l)    = q1(i,j)
END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO n=1,ntype
  DO l=1,land_pts
    gs_type(l,n)=gs(l)
    IF (l_irrig_dmd) THEN
      gs_irr_type(l,n)=gs(l)
    END IF
  END DO
END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO l=1,land_pts
  fsoil_tot(l) = frac(l,soil)
  gsoil(l) = 0.
  IF (v_crit(l,1) > 0.)                                                       &
    gsoil(l) = gs_nvg(soil - npft) *                                          &
               (sthu(l,1) * v_sat(l,1) / v_crit(l,1))**2
  ! repeated for irrigation
    gsoil_irr(l) = 0.0
  IF (l_irrig_dmd) THEN
    IF (v_crit(l,1) > 0.0)                                                    &
      gsoil_irr(l) = gs_nvg(soil - npft) *                                    &
                   (sthu_irr_gb(l,1) * v_sat(l,1) / v_crit(l,1))**2
  END IF
END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
gs(:)=0.0

!-----------------------------------------------------------------------
! Calculate light absorption by the plant canopy
!-----------------------------------------------------------------------
IF ( can_rad_mod /= 1 ) THEN
! The logical argument (getProfile in albpft) is TRUE to indicate
! that profiles through the canopy should be calculated.
  IF (l_albedo_obs .AND. l_spec_albedo .AND.                                  &
        ( .NOT. l_spec_alb_bs ) ) THEN
! scale the input PFT scattering and reflectivity parameters to match
! albedo obs, using pre-corrected scalings:
    albpft_call=2
  ELSE
! do not scale to match albedo obs:
    albpft_call=0
  END IF
!
IF ( albpft_call > 0 ) THEN
! Set all underlying albedos to that of bare soil.
! The necessary logical conditions for the assignment will have been
! checked when the albedo was calculated.
  DO n=1,npft
    DO l=1,land_pts
      albudif(l,1,n) = albsoil(l) * albobs_scaling_surft(l,soil,1)
      albudif(l,2,n) = albsoil(l) * albobs_scaling_surft(l,soil,2)
      albudir(l,1,n) = albsoil(l) * albobs_scaling_surft(l,soil,1)
      albudir(l,2,n) = albsoil(l) * albobs_scaling_surft(l,soil,2)
    END DO
  END DO
ELSE
  DO n=1,npft
    DO l=1,land_pts
      albudif(l,1,n) = albsoil(l)
      albudif(l,2,n) = albsoil(l)
      albudir(l,1,n) = albsoil(l)
      albudir(l,2,n) = albsoil(l)
    END DO
  END DO
END IF
!
! DEPENDS ON: albpft
  CALL albpft     (t_i_length * t_j_length,land_pts,land_index,               &
                   surft_index,surft_pts,ilayers,.TRUE.,albpft_call,          &
                   albudir,albudif,cos_zenith_angle,                          &
                   lai,alb_type_dummy,                                        &
                   fapar_dir,fapar_dif,fapar_dir2dif,                         &
                   fapar_dif2dif,fapar_dir2dir,fsun)
END IF



!-----------------------------------------------------------------------
! Loop over Plant Functional Types to calculate the available moisture
! and the values of canopy conductance, the carbon fluxes and the leaf
! turnover rate
!-----------------------------------------------------------------------
DO n=1,npft

  IF ( l_aggregate ) THEN
!$OMP PARALLEL DO IF(land_pts > 1) DEFAULT(NONE) PRIVATE(l) SHARED(land_pts,   &
!$OMP&             tstar, tstar_surft, z0, z0_surft) SCHEDULE(STATIC)
    DO l=1,land_pts
      tstar(l) = tstar_surft(l,1)
      z0(l) = z0_surft(l,1)
    END DO
!$OMP END PARALLEL DO
  ELSE
!$OMP PARALLEL DO IF(land_pts > 1) DEFAULT(NONE) PRIVATE(l) SHARED(land_pts,   &
!$OMP&            n, tstar, tstar_surft, z0, z0_surft) SCHEDULE(STATIC)
    DO l=1,land_pts
      tstar(l) = tstar_surft(l,n)
      z0(l) = z0_surft(l,n)
    END DO
!$OMP END PARALLEL DO
  END IF

  root_param(:,n)=rootd_ft(n)

  IF ( l_crop .AND. n > nnpft ) THEN
!$OMP PARALLEL DO IF(surft_pts(n) > 1) DEFAULT(NONE) PRIVATE(i, j)             &
!$OMP&            SHARED(cfrac_r, n, nnpft, rootc_cpft, root_param, rootd_ft,  &
!$OMP&                   rt_dir, surft_index, surft_pts) SCHEDULE(STATIC)
    DO j=1,surft_pts(n)
      i = surft_index(j,n)

      root_param(i,n) = rootd_ft(n) *                                         &
         ( rootc_cpft(i,n-nnpft) / cfrac_r(n-nnpft) ) ** rt_dir(n-nnpft)
    END DO
!$OMP END PARALLEL DO
  END IF

  DO j=1,surft_pts(n)
    i = surft_index(j,n)
! DEPENDS ON: root_frac
    CALL root_frac(n,nshyd,dzsoil,root_param(i,n),f_root(i,:))
  END DO

  !IF (l_irrig_dmd) THEN
! Soil moisture *in this tile* is a combination of soil moisture
! in irrigated and non-irrigated fraction of grid box, dependent
! on irrigated fraction *in this tile*
! However, in the irrigated fraction *in this tile*, soil moisture
! is the same as in the overall gridbox irrigated fraction
!$OMP PARALLEL DO IF(land_pts > 1) DEFAULT(NONE) PRIVATE(l, k)                 &
!$OMP&            SHARED(frac_irr_gb, frac_irr_surft, land_pts, l_irrig_dmd, n,&
!$OMP&                   nshyd, sthu, sthu_irr_gb, sthu_nir, sthu_surft)       &
!$OMP&            SCHEDULE(STATIC)
    DO l=1,land_pts
      DO k=1,nshyd
        sthu_surft(l,k) = sthu(l,k)

        IF ( l_irrig_dmd ) THEN
          sthu_nir(l,k) = sthu(l,k)
          IF ( frac_irr_gb(l) < 1.0 ) THEN
            sthu_nir(l,k) = ( sthu(l,k) - frac_irr_gb(l) * sthu_irr_gb(l,k) ) &
                          / (1.0 - frac_irr_gb(l))
          ELSE
            sthu_nir(l,k) = sthu_irr_gb(l,k)
          END IF
          sthu_surft(l,k) = frac_irr_surft(l,n) * sthu_irr_gb(l,k) +          &
                          (1.0 - frac_irr_surft(l,n)) * sthu_nir(l,k)
        END IF
      END DO
    END DO
!$OMP END PARALLEL DO
  !END IF

  CALL smc_ext (land_pts,nshyd,surft_pts(n),surft_index(:,n), n               &
,               f_root,sthu_surft,v_crit,v_sat,v_wilt                         &
,               wt_ext_type(:,:,n),fsmc(:,n))

  IF (l_irrig_dmd) THEN
    CALL smc_ext (land_pts,nshyd,surft_pts(n),surft_index(:,n), n             &
,               f_root,sthu_irr_gb,v_crit,v_sat,v_wilt                        &
,               wt_ext_irr_type(:,:,n),fsmc_irr(:,n))
  END IF

! DEPENDS ON: raero
  CALL raero (land_pts,land_index,surft_pts(n),surft_index(:,n)               &
,             rib,wind,z0,z0,z1,ra)
!-----------------------------------------------------------------------
! Calculate light absorption by the plant canopy
!-----------------------------------------------------------------------

  IF (can_rad_mod /= 1) THEN

!$OMP PARALLEL DO IF(ilayers > 1) DEFAULT(NONE) PRIVATE(i, k, l, m)           &
!$OMP SHARED(ilayers, surft_pts, surft_index, land_index, faparv, diff_frac,  &
!$OMP        fapar_dir, fapar_dif, n)                                         &
!$OMP SCHEDULE(STATIC)
    DO m=1,ilayers
      DO k=1,surft_pts(n)
        l = surft_index(k,n)
        i = land_index(l)

        faparv(l,m) = (1.-diff_frac(i)) * fapar_dir(l,n,m)                    &
                    + diff_frac(i) * fapar_dif(l,n,m)
      END DO !  points
    END DO  !  layer
!$OMP END PARALLEL DO

    IF ( can_rad_mod == 5 .OR. can_rad_mod == 6 ) THEN
!$OMP PARALLEL DO IF(ilayers > 1) DEFAULT(NONE) PRIVATE(i, k, l, m)           &
!$OMP SHARED(ilayers, surft_pts, surft_index, land_index, fapar_shd,          &
!$OMP        diff_frac, fapar_dif2dif, fapar_dir2dir, fsun, fapar_sun,        &
!$OMP        fapar_dir2dif, n)                                                &
!$OMP SCHEDULE(STATIC)
      DO m=1,ilayers
        DO k=1,surft_pts(n)
          l = surft_index(k,n)
          i = land_index(l)

          fapar_shd(l,m) = diff_frac(i) * fapar_dif2dif(l,n,m)                &
               + ( 1.0 - diff_frac(i) ) * fapar_dir2dif(l,n,m)

          IF ( fsun(l,n,m) > EPSILON(fsun) ) THEN
            fapar_sun(l,m) = fapar_shd(l,m)                                   &
                + ( 1. - diff_frac(i) ) * fapar_dir2dir(l,n,m )               &
                / fsun(l,n,m)
          ELSE
            fapar_sun(l,m) = 0.0
          END IF

        END DO !  points
      END DO  !  layer
!$OMP END PARALLEL DO
    END IF  !  can_rad_mod=5/6

  END IF  !  can_rad_mod

  CALL sf_stom (land_pts,land_index                                           &
,               surft_pts(n),surft_index(:,n),n                               &
,               co2,co2_3d,co2_dim_len                                        &
,               co2_dim_row,l_co2_interactive                                 &
,               fsmc(:,n),ht(:,n),ipar_land,lai(:,n)                          &
,               ht(:,n),pstar_land                                            &
,               q1_land,ra,tstar,o3                                           &
,               can_rad_mod,ilayers,faparv                                    &
,               gpp_pft(:,n),npp_pft(:,n),resp_p_pft(:,n)                     &
,               resp_l_pft(:,n),resp_r_pft(:,n),resp_w_pft(:,n)               &
,               n_leaf(:,n),n_root(:,n),n_stem(:,n)                           &
,               lai_bal(:,n),et_stom_pft(:,n)                                 &
,               gs_type(:,n)                                                  &
,               fapar_sun,fapar_shd,fsun(:,n,:)                               &
,               flux_o3_pft(:,n),fo3_pft(:,n),fapar_diag_pft(:,n)             &
,               isoprene_pft(:,n),terpene_pft(:,n)                            &
,               methanol_pft(:,n),acetone_pft(:,n)                            &
,               open_index,open_pts                                           &
                )
  IF (l_irrig_dmd) THEN
! adjust conductance for irrigated fraction
!$OMP PARALLEL DO IF(surft_pts(n) > 1) DEFAULT(NONE) PRIVATE(k, l)             &
!$OMP&            SHARED(gs_irr_type, gs_type, n, surft_index, surft_pts)      &
!$OMP&            SCHEDULE(STATIC)
    DO k=1,surft_pts(n)
      l = surft_index(k,n)
      gs_irr_type(l,n) = gs_type(l,n)
    END DO
!$OMP END PARALLEL DO

    DO j=1,open_pts
      l = surft_index(open_index(j),n)

      IF ( frac_irr_surft(l,n) > 0.0 ) THEN
        IF ( fsmc(l,n) > 0.0 ) THEN
          gs_irr_type(l,n) = gs_type(l,n) * fsmc_irr(l,n) / fsmc(l,n)

        ELSE
! hadrd - this should only happen if sthu is 0.0, see smc_ext
          WRITE (errmsg,*) 'tile:', n, 'point:',l, 'fsmc:', fsmc(l,n),        &
              'sthu', sthu(l,n)
          CALL ereport("physiol", 101,                                        &
              "***warning fsmc <= 0 in physiol.F90***" //                     &
               errmsg )
        END IF
      END IF
    END DO
  END IF

! DEPENDS ON: soil_evap
  CALL soil_evap (land_pts,nshyd,surft_pts(n),surft_index(:,n)                &
,                 gsoil,lai(:,n),gs_type(:,n),wt_ext_type(:,:,n)              &
,                 fsoil(:,n)                                                  &
,                 gsoil_irr,gs_irr_type(:,n)                                  &
,                 wt_ext_irr_type(:,:,n)                                      &
                  )

! DEPENDS ON: leaf_lit
  CALL leaf_lit (land_pts,surft_pts(n),surft_index(:,n)                       &
,                n,fsmc(:,n),tstar,g_leaf(:,n))

! DEPENDS ON: cancap
  CALL cancap (land_pts,surft_pts(n),surft_index(:,n),can_model,n             &
,              ht(:,n),lai(:,n),ch_type(:,n),vf_type(:,n))

!$OMP PARALLEL DO IF(land_pts > 1) DEFAULT(NONE) PRIVATE(l) SHARED(frac,       &
!$OMP&            fsoil, fsoil_tot, land_pts, n) SCHEDULE(STATIC)
  DO l=1,land_pts
    fsoil_tot(l) = fsoil_tot(l) + frac(l,n)*fsoil(l,n)
  END DO
!$OMP END PARALLEL DO

END DO

!----------------------------------------------------------------------
! Non-vegetated surface types
!----------------------------------------------------------------------
DO n=npft+1,ntype
!$OMP PARALLEL DO IF(surft_pts(n) > 1) DEFAULT(NONE) PRIVATE(l, m)             &
!$OMP&            SHARED(gs_irr_type, gs_nvg, gs_type, l_irrig_dmd, n, npft,   &
!$OMP&                   surft_index, surft_pts) SCHEDULE(STATIC)
  DO m=1,surft_pts(n)
    l=surft_index(m,n)
    gs_type(l,n) = gs_nvg(n-npft)
    IF (l_irrig_dmd) THEN
      gs_irr_type(l,n) = gs_nvg(n-npft) ! irrigation
    END IF
  END DO
!$OMP END PARALLEL DO
END DO

! Copy soil conductance and add bare soil fraction to extraction from
! surface layer
n = soil
!$OMP PARALLEL DO IF (surft_pts(n) > 1) DEFAULT(NONE) PRIVATE(l, m)            &
!$OMP&            SHARED(gs_irr_type, gs_type, gsoil, gsoil_irr, l_irrig_dmd,  &
!$OMP&                   n, surft_index, surft_pts, wt_ext_type,               &
!$OMP&                   wt_ext_irr_type) SCHEDULE(STATIC)
DO m=1,surft_pts(n)
  l=surft_index(m,n)
  gs_type(l,n) = gsoil(l)
  wt_ext_type(l,1,n) = 1.0
  IF (l_irrig_dmd) THEN
    gs_irr_type(l,n) = gsoil_irr(l) ! irrigation
    wt_ext_irr_type(l,1,n) = 1.0 ! irrigation
  END IF
END DO
!$OMP END PARALLEL DO

!----------------------------------------------------------------------
! Canopy heat capacity and coverage for non-vegetated surfaces
!----------------------------------------------------------------------
DO n=npft+1,ntype
!$OMP PARALLEL DO IF(surft_pts(n) > 1) DEFAULT(NONE) PRIVATE(l, m)             &
!$OMP&            SHARED(ch_nvg, ch_type, n, npft, surft_index, surft_pts,     &
!$OMP&                   vf_nvg, vf_type) SCHEDULE(STATIC)
  DO m=1,surft_pts(n)
    l=surft_index(m,n)
    ch_type(l,n) = ch_nvg(n-npft)
    vf_type(l,n) = vf_nvg(n-npft)
  END DO
!$OMP END PARALLEL DO
END DO

! MORUSES - Update ch_type (ch_nvg) with calculated values for urban fabric.
! Roof is done at the same time as canyon; cannot have one without the other.
! Differs from original coding for efficiency reasons; taken out of previous
! loop.
IF ( l_moruses_storage ) THEN

  dz_wall   = ( 2.0 * diffus_wall / omega_day )**( 1.0 / 2.0 )
  dz_road   = ( 2.0 * diffus_road / omega_day )**( 1.0 / 2.0 )
  dz_roof_c = ( 2.0 * diffus_roof / omega_day )**( 1.0 / 2.0 )
  ! dz_roof is very thin to represent insulation
  IF ( l_moruses_storage_thin ) THEN
    dz_roof = MIN( dz_roof_c, dz_roof_p )
  ELSE
    dz_roof=dz_roof_c
  END IF

  n = urban_canyon
!$OMP PARALLEL DO IF(surft_pts(n) > 1) DEFAULT(NONE) PRIVATE(l, m)             &
!$OMP&            SHARED(ch_type, dz_road, dz_wall, dz_roof, hwr_gb, l_aggregate, &
!$OMP&                   n, surft_index, surft_pts, urban_roof, vf_type)       &
!$OMP&            SCHEDULE(STATIC)
  DO m = 1, surft_pts(n)
    l = surft_index(m,n)

    ! Canyon
    ch_type(l,n) = cap_road * dz_road +                                       &
       ( 2.0 * hwr_gb(l) * cap_wall * dz_wall )

    ! Where there's a roof there's a canyon
    ch_type(l,urban_roof) = cap_roof * dz_roof

  END DO
!$OMP END PARALLEL DO

END IF

!----------------------------------------------------------------------
! Surface emissivity
!----------------------------------------------------------------------

! Initialise emis_surft
emis_surft(:,:) = 0.0

! URBAN-2T & MORUSES: Set canyon emissivity

IF ( l_moruses_emissivity ) THEN
  n = urban_canyon
  IF ( firstcall) THEN
    CALL jules_print('physiol_jls',                                           &
        'MORUSES canyon emissivity calculated',level=PrNorm)
  END IF
  DO m = 1, surft_pts(n)
    l = surft_index(m,n)
    ! DEPENDS ON: urbanemis
    CALL urbanemis(hwr_gb(l), emisr_gb(l), emiss, emisw_gb(l),                &
       emisc(l))
  END DO
END IF

! Calculate EMIS_surft

IF ( l_aggregate ) THEN

  DO n=1,npft
!$OMP PARALLEL DO IF(surft_pts(n) > 1) DEFAULT(NONE) PRIVATE(l, m)             &
!$OMP&          SHARED(emis_pft, emis_surft, frac, n, surft_index, surft_pts)  &
!$OMP&          SCHEDULE(STATIC)
    DO m=1,surft_pts(n)
      l = surft_index(m,n)
      emis_surft(l,1) = emis_surft(l,1) + frac(l,n) * emis_pft(n)
    END DO
!$OMP END PARALLEL DO
  END DO

! MORUSES: Use emisc instead of EMIS_NVG if canyon tile is updated by MORUSES

! Implementation inelegant. If emis_nvg becomes dependent on land_pts then
! could overwrite emis_nvg with SUBROUTINE urbanemis and do away with the
! if statements. Alternatively could create a 2d work array e.g.
! emis_nvg_wk(land_pts,ntype) or re-write to have just
! emis_surft(land_pts,ntype) then copy aggregated value to
! emis_surft(land_pts,1). Similar to albedo in tile_albedo.

  IF ( l_moruses_emissivity ) THEN

! MORUSES

    DO n=npft+1,ntype
      IF ( n == urban_canyon ) THEN
        DO m=1,surft_pts(n)
          l = surft_index(m,n)
          IF ( firstcall ) THEN
            WRITE(jules_message,*) 'EMIS_surft: Emissivity of canyon used',   &
                n,m,l
            CALL jules_print('physiol_jls',jules_message,level=PrNorm)
          END IF
          emis_surft(l,1) = emis_surft(l,1)                                   &
             + frac(l,n) * emisc(l)
        END DO
      ELSE
!$OMP PARALLEL DO IF(surft_pts(n) > 1) DEFAULT(NONE) PRIVATE(l, m)             &
!$OMP&             SHARED(emis_nvg, emis_surft, frac, n, npft, surft_index,    &
!$OMP&                    surft_pts) SCHEDULE(STATIC)
        DO m=1,surft_pts(n)
          l = surft_index(m,n)
          emis_surft(l,1) = emis_surft(l,1)                                   &
             + frac(l,n) * emis_nvg(n - npft)
        END DO
!$OMP END PARALLEL DO
      END IF
    END DO

    firstcall = .FALSE.

  ELSE ! .NOT. l_moruses_emissivity
    DO n=npft+1,ntype
!$OMP PARALLEL DO IF(surft_pts(n) > 1) DEFAULT(NONE) PRIVATE(l, m)             &
!$OMP&             SHARED(emis_nvg, emis_surft, frac, n, npft, surft_index,    &
!$OMP&                    surft_pts) SCHEDULE(STATIC)
      DO m=1,surft_pts(n)
        l = surft_index(m,n)
        emis_surft(l,1) = emis_surft(l,1)                                     &
           + frac(l,n) * emis_nvg(n - npft)
      END DO
!$OMP END PARALLEL DO
    END DO
  END IF

ELSE ! .NOT. L_AGGREGATE

  DO n=1,npft

!$OMP PARALLEL DO IF(surft_pts(n) > 1) DEFAULT(NONE) PRIVATE(l, m)             &
!$OMP&             SHARED(emis_pft, emis_surft, n, surft_index, surft_pts)     &
!$OMP&             SCHEDULE(STATIC)
    DO m=1,surft_pts(n)
      l = surft_index(m,n)
      emis_surft(l,n) = emis_pft(n)
    END DO
!$OMP END PARALLEL DO
  END DO

  DO n=npft+1,ntype
    DO m=1,surft_pts(n)
      l = surft_index(m,n)
      emis_surft(l,n) = emis_nvg(n - npft)
    END DO
  END DO

! MORUSES: Overwrite EMIS_surft for urban canyon
  IF ( l_moruses_emissivity ) THEN
    n = urban_canyon
!$OMP PARALLEL DO IF(surft_pts(n) > 1) DEFAULT(NONE) PRIVATE(l, m)             &
!$OMP             SHARED(emisc, emis_surft, n, surft_index, &
!$OMP                    surft_pts) SCHEDULE(STATIC)
    DO m = 1, surft_pts(n)
      l = surft_index(m,n)
      emis_surft(l,n) = emisc(l)
    END DO
!$OMP END PARALLEL DO
    firstcall = .FALSE.
  END IF

END IF ! L_AGGREGATE

DO l=1,land_pts
  emis_soil(l) = emis_nvg(soil - npft)
END DO

!----------------------------------------------------------------------
! Calculate the rate of soil respiration
!----------------------------------------------------------------------
! set VEG_FRAC according to whether it is full or dummy field
IF (l_triffid) THEN
!$OMP PARALLEL DO IF(land_pts > 1) DEFAULT(NONE) PRIVATE(l) SHARED(frac,       &
!$OMP             land_pts, npft, veg_frac) SCHEDULE(STATIC)
  DO l=1,land_pts
    veg_frac(l) = SUM(frac(l,1:npft))
  END DO
!$OMP END PARALLEL DO
END IF
IF (l_soil_resp_lev2) THEN
! DEPENDS ON: microbe
  CALL microbe (land_pts,dim_cs1,dim_cs2,l_triffid,l_q10,cs,                  &
               sthu(:,2)+sthf(:,2),v_sat(:,2),v_wilt(:,2),tsoil(:,2),resp_s,       &
               veg_frac)
ELSE
! DEPENDS ON: microbe
  CALL microbe (land_pts,dim_cs1,dim_cs2,l_triffid,l_q10,cs,                  &
               sthu(:,1),v_sat(:,1),v_wilt(:,1),tsoil(:,1),resp_s,                 &
               veg_frac)
ENDIF

!----------------------------------------------------------------------
! Form gridbox mean values
!----------------------------------------------------------------------

DO n=1,ntype
!$OMP PARALLEL DO IF(surft_pts(n) > 1) DEFAULT(NONE) PRIVATE(l, m)             &
!$OMP             SHARED(frac, gs, gs_type, n, surft_index, surft_pts)         &
!$OMP             SCHEDULE(STATIC)
  DO m=1,surft_pts(n)
    l=surft_index(m,n)
    gs(l) = gs(l) + frac(l,n)*gs_type(l,n)
  END DO
!$OMP END PARALLEL DO
END DO

IF ( l_aggregate ) THEN
  DO n=1,ntype
!$OMP PARALLEL DO IF(surft_pts(n) > 1) DEFAULT(NONE) PRIVATE(l, m)             &
!$OMP             SHARED(canhc, ch_type, frac, l_irrig_dmd, n, nshyd,          &
!$OMP                    surft_index, surft_pts, vfrac, vf_type, wt_ext,       &
!$OMP                    wt_ext_irr, wt_ext_type, wt_ext_irr_type)             &
!$OMP             SCHEDULE(STATIC)
    DO m=1,surft_pts(n)
      l=surft_index(m,n)
      canhc(l) = canhc(l) + frac(l,n)*ch_type(l,n)
      vfrac(l) = vfrac(l) + frac(l,n)*vf_type(l,n)
      DO k=1,nshyd
        wt_ext(l,k) = wt_ext(l,k) + frac(l,n)*wt_ext_type(l,k,n)
        IF (l_irrig_dmd) THEN
          wt_ext_irr(l,k) = wt_ext_irr(l,k)                                   &
                          + frac(l,n) * wt_ext_irr_type(l,k,n)
        END IF
      END DO
    END DO
!$OMP END PARALLEL DO
  END DO
!$OMP PARALLEL DO IF(land_pts > 1) DEFAULT(NONE) PRIVATE(k, l)                 &
!$OMP             SHARED(canhc, canhc_surft, flake, frac, gs, gs_surft, lake,  &
!$OMP                    land_pts, nshyd, vfrac, vfrac_surft, wt_ext,          &
!$OMP                    wt_ext_surft) SCHEDULE(STATIC)
  DO l=1,land_pts
    IF ( lake > 0 ) THEN
      flake(l,1) = frac(l,lake)
    ELSE
      flake(l,1) = 0.
    END IF
    gs_surft(l,1) = 0.
    IF (flake(l,1) <  1.)                                                     &
      gs_surft(l,1) = gs(l) / (1. - flake(l,1))
    canhc_surft(l,1) = canhc(l)
    vfrac_surft(l,1) = vfrac(l)
    DO k=1,nshyd
      wt_ext_surft(l,k,1) = wt_ext(l,k)
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE
  gs_surft(:,:)=0.0
  canhc_surft(:,:)=0.0
  vfrac_surft(:,:)=0.0
  IF (l_irrig_dmd) THEN
    gs_irr_surft(:,:)=0.0 ! irrigation
  END IF
  DO n=1,ntype
!$OMP PARALLEL DO IF(surft_pts(n) > 1) DEFAULT(NONE) PRIVATE(k, l, m)         &
!$OMP SHARED(surft_pts, surft_index, flake, gs_surft, gs_type, l_irrig_dmd,   &
!$OMP        gs_irr_surft, gs_irr_type, canhc_surft, ch_type, vfrac_surft,    &
!$OMP        vf_type, nshyd, wt_ext, frac, wt_ext_type, wt_ext_surft,         &
!$OMP        wt_ext_irr, wt_ext_irr_type, wt_ext_irr_surft, n)                &
!$OMP             SCHEDULE(STATIC)
    DO m=1,surft_pts(n)
      l=surft_index(m,n)
      flake(l,n) = 0.
      gs_surft(l,n) = gs_type(l,n)
      IF (l_irrig_dmd) THEN
        gs_irr_surft(l,n) = gs_irr_type(l,n) ! irrigation
      END IF
      canhc_surft(l,n) = ch_type(l,n)
      vfrac_surft(l,n) = vf_type(l,n)

      DO k=1,nshyd
        wt_ext(l,k) = wt_ext(l,k) + frac(l,n)*wt_ext_type(l,k,n)
        wt_ext_surft(l,k,n) = wt_ext_type(l,k,n)
        IF (l_irrig_dmd) THEN
          wt_ext_irr(l,k) = wt_ext_irr(l,k) + frac(l,n) * wt_ext_irr_type(l,k,n)
                                                          ! irrigation
          wt_ext_irr_surft(l,k,n) = wt_ext_irr_type(l,k,n)
        END IF
      END DO
    END DO
!$OMP END PARALLEL DO
  END DO
  IF ( lake > 0 ) THEN
    n = lake    ! Lake tile
!$OMP PARALLEL DO IF(surft_pts(n) > 1) DEFAULT(NONE) PRIVATE(l, m)             &
!$OMP             SHARED(flake, n, surft_index, surft_pts) SCHEDULE(STATIC)
    DO m=1,surft_pts(n)
      l=surft_index(m,n)
      flake(l,n) = 1.
    END DO
!$OMP END PARALLEL DO
  END IF
END IF

DO n=1,npft
!$OMP PARALLEL DO IF(surft_pts(n) > 1) DEFAULT(NONE) PRIVATE(l,m)              &
!$OMP SHARED(surft_pts, surft_index, n, gpp, frac, gpp_pft, npp, npp_pft,      &
!$OMP        resp_p, resp_p_pft, isoprene_gb, isoprene_pft, terpene_gb,        &
!$OMP        terpene_pft, methanol_gb, methanol_pft, acetone_gb, acetone_pft,  &
!$OMP        et_stom_gb, et_stom_pft)                                          &
!$OMP             SCHEDULE(STATIC)
  DO m=1,surft_pts(n)
    l=surft_index(m,n)

    gpp(l)      = gpp(l) + frac(l,n) * gpp_pft(l,n)
    npp(l)      = npp(l) + frac(l,n) * npp_pft(l,n)
    resp_p(l)   = resp_p(l) + frac(l,n) * resp_p_pft(l,n)

    isoprene_gb(l) = isoprene_gb(l) + frac(l,n) * isoprene_pft(l,n)
    terpene_gb(l)  = terpene_gb(l)  + frac(l,n) * terpene_pft(l,n)
    methanol_gb(l) = methanol_gb(l) + frac(l,n) * methanol_pft(l,n)
    acetone_gb(l)  = acetone_gb(l)  + frac(l,n) * acetone_pft(l,n)
    et_stom_gb(l)  = et_stom_gb(l)  + frac(l,n) + et_stom_pft(l,n)

  END DO
!$OMP END PARALLEL DO
END DO

!----------------------------------------------------------------------
! Diagnose the available moisture in the soil profile
!----------------------------------------------------------------------
!!AJM START
! Available water for plant transpiration
DO n=1,nshyd
!$OMP PARALLEL DO IF(land_pts > 1) DEFAULT(NONE) PRIVATE(l) SHARED(dzsoil,     &
!$OMP             land_pts, n, smct, sthu, v_sat, v_wilt, wt_ext)              &
!$OMP             SCHEDULE(STATIC)
  DO l=1,land_pts
    smct(l) = smct(l) + MAX( 0. ,                                             &
         wt_ext(l,n)*rho_water*dzsoil(n)*                                     &
               (sthu(l,n)*v_sat(l,n)-v_wilt(l,n)))
  END DO
!$OMP END PARALLEL DO
END DO

! Add available water for evaporation from bare soil
!$OMP PARALLEL DO IF(land_pts > 1) DEFAULT(NONE) PRIVATE(l) SHARED(dzsoil,     &
!$OMP             fsoil_tot, land_pts, smct, sthu, v_sat)                      &
!$OMP             SCHEDULE(STATIC)
DO l=1,land_pts
   smct(l) = (1.0-fsoil_tot(l))*smct(l) +                                     &
                   fsoil_tot(l)*rho_water*dzsoil(1)*                          &
                       MAX(0.0,sthu(l,1))*v_sat(l,1)
END DO
!$OMP END PARALLEL DO

IF (l_irrig_dmd) THEN
! Available water for plant transpiration in irrigated fraction.
  DO n=1,nshyd
!$OMP PARALLEL DO IF(land_pts > 1) DEFAULT(NONE) PRIVATE(l) SHARED(dzsoil,     &
!$OMP             land_pts, n, smc_irr_gb, sthu_irr_gb, v_sat, v_wilt,         &
!$OMP             wt_ext_irr) SCHEDULE(STATIC)
    DO l=1,land_pts
      smc_irr_gb(l) = smc_irr_gb(l) + MAX( 0.0 ,                              &
           wt_ext_irr(l,n) * rho_water * dzsoil(n) *                          &
           ( sthu_irr_gb(l,n) * v_sat(l,n) - v_wilt(l,n) ) )
    END DO
!$OMP END PARALLEL DO
  END DO

! Add available water for evaporation from bare soil in irrig frac.
!$OMP PARALLEL DO IF(land_pts > 1) DEFAULT(NONE) PRIVATE(l) SHARED(dzsoil,     &
!$OMP             fsoil_tot, land_pts, smc_irr_gb, sthu_irr_gb, v_sat)         &
!$OMP             SCHEDULE(STATIC)
  DO l=1,land_pts
     smc_irr_gb(l) = (1.0-fsoil_tot(l))*smc_irr_gb(l) +                       &
                   fsoil_tot(l)*rho_water*dzsoil(1)*                          &
                       MAX(0.0,sthu_irr_gb(l,1))*v_sat(l,1)
  END DO
!$OMP END PARALLEL DO

  gc_irr_surft(:,:)=gs_irr_surft(:,:)

END IF


IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE physiol
