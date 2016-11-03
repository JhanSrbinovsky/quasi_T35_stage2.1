! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine TRIFFID ------------------------------------------------
!
!                     Top-down
!                     Representation of
!                     Interactive
!                     Foliage and
!                     Flora
!                     Including
!                     Dynamics
!
! Purpose : Simulates changes in vegetation structure, areal
!           coverage and the carbon contents of vegetation and soil.
!           can be used to advance these variables dynamically
!           (GAMMA=1/TIMESTEP) or to iterate towards  equilibrium
!           (GAMMA --> 0.0, FORW=1.0).
!
! -------------------------------------------------------------------
SUBROUTINE triffid (land_pts,trif_pts,trif_index,forw,GAMMA                   &
,                   frac_agric,frac_past,g_leaf,npp,resp_s                    &
,                   resp_w,cs,frac,ht,lai,resp_frac                           &
,                   c_veg,cv,lit_c,lit_c_t                                    &
,                   wood_prod_total_gb,wp_total_out_gb)

USE jules_surface_types_mod

USE jules_surface_mod, ONLY: cmass, n_inorg_turnover

USE pftparm
USE trif

USE prognostics, ONLY: wood_prod_fast_gb, wood_prod_med_gb,                   &
                       wood_prod_slow_gb, frac_agr_prev_gb,                   &
                       frac_past_prev_gb, n_inorg_gb

USE trif_vars_mod, ONLY: cnsrv_veg_triffid_gb, cnsrv_soil_triffid_gb,         &
                         cnsrv_prod_triffid_gb, cnsrv_carbon_triffid_gb,      &
                         wp_fast_in_gb, wp_med_in_gb, wp_slow_in_gb,          &
                         wp_fast_out_gb, wp_med_out_gb, wp_slow_out_gb,       &
                         lit_c_orig_pft, lit_c_ag_pft, lit_n_orig_pft,        &
                         lit_n_ag_pft, pc_s_pft, leafc_pft, rootc_pft,        &
                         stemc_pft, woodc_pft, droot_pft, dleaf_pft,          &
                         dwood_pft, root_litc_pft, leaf_litc_pft,             &
                         wood_litc_pft, root_litn_pft, leaf_litn_pft,         &
                         wood_litn_pft, litterc_pft, lit_n_t_gb, lit_n_pft,   &
                         n_fix_pft, n_fix_gb, n_gas_gb, n_uptake_pft,         &
                         n_uptake_gb, n_demand_pft, n_demand_gb, exudates_pft,&
                         exudates_gb, littern_pft, n_uptake_growth_pft,       &
                         n_demand_growth_pft, n_demand_spread_pft,            &
                         n_uptake_spread_pft, n_demand_lit_pft, n_veg_gb,     &
                         n_veg_pft, dcveg_pft, dcveg_gb, dnveg_pft,           &
                         lai_bal_pft, dnveg_gb, deposition_n_gb, n_loss_gb,   &
                         harvest_pft,  root_abandon_pft,  harvest_gb,         &
                         harvest_n_pft, harvest_n_gb, n_fertiliser_pft,       &
                         n_fertiliser_gb, root_abandon_n_pft


USE jules_vegetation_mod, ONLY: l_veg_compete, l_ht_compete,                  &
                                l_trait_phys,l_trif_eq,l_landuse,             &
                                can_rad_mod, l_nitrogen,                      &
                                l_trif_crop

USE CN_utils_mod, ONLY: calc_n_comps_triffid

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

INTEGER                                                                       &
 land_pts                                                                     &
                            ! IN Total number of land points.
,trif_pts                                                                     &
                            ! IN Number of points on which
!                                 !    TRIFFID may operate
,trif_index(land_pts)                                                         &
                            ! IN Indices of land points on
!                                 !    which TRIFFID may operate
,l,n,t,k                    ! WORK Loop counters

INTEGER                                                                       &
 nsub                                                                         &
                            ! WORK number of PFTs to use in call
                            !      to lotka_noeq_subset
                            !      (if l_trif_crop is true)
,crop_lotka
                            ! WORK indicates to l_trif_crop if
                            !      crop or natural PFTs are
                            !      competing. 1 = crops as in the
                            !      "crop" variable

REAL                                                                          &
 forw                                                                         &
                            ! IN Forward timestep weighting.
,GAMMA                                                                        &
                            ! IN Inverse timestep (/360days).
,frac_agric(land_pts)                                                         &
                            ! IN Fraction of agriculture.
,frac_past(land_pts)                                                          &
                            ! IN Fraction of pasture.
,g_leaf(land_pts,npft)                                                        &
                            ! IN Turnover rate for leaf and
!                                 !    fine root biomass (/360days).
,npp(land_pts,npft)                                                           &
                            ! INOUT Net primary productivity
!                                 !       (kg C/m2/360days).
,resp_s(land_pts,5)                                                           &
                            ! INOUT Soil respiration
!                                 !       (kg C/m2/360days).
,resp_w(land_pts,npft)                                                        &
                            ! INOUT Wood maintenance respiration
!                                 !       (kg C/m2/360days).
,cs(land_pts,4)                                                               &
                            ! INOUT Soil carbon (kg C/m2).
,frac(land_pts,ntype)                                                         &
                            ! INOUT Fractional cover of each
!                                 !       Functional Type.
,frac_na(land_pts,ntype)                                                      &
,ht(land_pts,npft)                                                            &
                            ! INOUT Vegetation height (m).
,lai(land_pts,npft)                                                           &
                            ! INOUT Leaf area index.
,c_veg(land_pts,npft)                                                         &
                            ! OUT Total carbon content of
!                                 !     the vegetation (kg C/m2).
,cv(land_pts)                                                                 &
                            ! OUT Gridbox mean vegetation
!                                 !     carbon (kg C/m2).
,lit_c(land_pts,npft)                                                         &
                            ! OUT Carbon Litter (kg C/m2/360days).

,lit_c_t(land_pts)          ! OUT Gridbox mean carbon litter
!                                 !     (kg C/m2/360days).

REAL                                                                          &
 dfrac(land_pts,npft)                                                         &
                            ! WORK Change in areal fraction
!                                 !      during the timestep
!                                 !      (/timestep).
,dfrac_na(land_pts,npft)                                                      &
,frac_flux                                                                    &
                            ! WORK PFT fraction to be used
!                                 !      in the calculation of
!                                 !      the gridbox mean fluxes.
,frac_flux_nat                                                                &
                                  ! WORK PFT fraction to be used
!                                 !      in the calculation of
!                                 !      the gridbox mean fluxes.
,resp_frac(land_pts)                                                          &
                            ! WORK  respired fraction of RESP_S
,leaf(land_pts,npft)                                                          &
                            ! WORK Leaf biomass for balanced LAI
                                  !      (kg C/m2).
!pc_s_pft is a diagnostic now (in trif_vars_mod)
!,pc_s_pft(land_pts,npft)                                              &
!                            ! WORK Net carbon flux available
!                                 !      for spreading
!                                 !      (kg C/m2/yr).
,ns_gb(land_pts)                                                              &
                            ! WORK Total Soil N
                                  !    (kg N/m2)
,phen(land_pts,npft)                                                          &
                            ! WORK Phenological state.
,root(land_pts,npft)                                                          &
                            ! WORK Root biomass (kg C/m2).
,wood(land_pts,npft)                                                          &
                            ! WORK Woody biomass (kg C/m2).
,cratio                                                                       &
                            ! WORK Ratio above/below veg carbon
,nratio                                                                       &
                            ! WORK Ratio above/below veg nitrogen
,frac_prev(land_pts,ntype)                                                    &
                            ! WORK Copy of frac for equilibrium
!                                 !       TRIFFID
,n_unmet(land_pts,npft)                                                       &
                            ! WORK N this unmet through exudate limit
,n_unmet_gb(land_pts)                                                         &
                            ! WORK N this unmet through exudate limit
,deposition(land_pts)                                                         &
                            ! WORK N deposition (kgN m-2 360days-1)
,n_avail(land_pts,npft)                                                       &
                            ! WORK inorganic N on tiles (kgN m-2)
,n_leaf                                                                       &
                            ! WORK Leaf N pool (kg2 m-2)
,n_root                                                                       &
                            ! WORK Root N pool (kg2 m-2)
,n_stem                                                                       &
                            ! WORK Stem N pool (kg2 m-2)
,neg_n(land_pts)
                            ! OUT Negative N required to prevent ns<0.
                                  !    (kg N)

REAL                                                                          &
subset_space(land_pts)                                                        &
                            ! WORK Used if l_trif_crop=.T.
                            ! For the natural vegetation type this
                            ! is the area where vegetation is prevented
                            ! from growing becasue the area has been
                            ! reserved for crops/pasture/forestry. For
                            ! other vegetation types this is the area
                            ! reserved for them to grow in, e.g. crops.
,subset_space_prev(land_pts)
                            ! WORK the value of subset_space at the
                            ! previous TRIFFID timestep.

REAL                                                                          &
cnsrv_veg_flux(land_pts)                                                      &
                            ! WORK net carbon flux into vegetation
,cnsrv_soil_flux(land_pts)                                                    &
                            ! WORK net carbon flux into soil
,cnsrv_prod_flux(land_pts)                                                    &
                            ! WORK net carbon flux into wood products
,cnsrv_carbon_flux(land_pts)                                                  &
                            ! WORK net carbon flux into land
,implicit_resp_correction(land_pts)
                            ! WORK respiration carried to next triffid
                            !      timestep to account for applying
                            !      minimum soil carbon constraint
REAL::                                                                        &
wood_prod_total_gb(land_pts)                                                  &
                            ! OUT total of the wood product pools
,wp_total_out_gb(land_pts)
                            ! OUT total flux out of wood products
INTEGER, PARAMETER :: veg_types = 3
                              ! Thr number of vegetation types.
                              ! PFTs only compete with other PFTs
                              ! of the same vegetation type.
                              ! 1=natural vegetation
                              ! 2= nat veg, crops
                              ! 3= nat veg, crops, pasture
                              ! 4=nat, crop, pasture, forestry

REAL, PARAMETER :: harvest_rate = 0.3
                            ! Fraction of litter diverted as harvest
                            ! This fraction ends up in the product
                            !  pools instead of in the soil

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='TRIFFID'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

!----------------------------------------------------------------------
! Check carbon conservation (1/3)
! Calculate carbon stores and net carbon fluxes at start of
! routine. The wood product decay flux has not yet been calculated and
! will be added later. The litter and wood product decay flux have not
! yet been calculated and will be added later.
!----------------------------------------------------------------------
  DO t=1,trif_pts
    l=trif_index(t)

    DO n=1,nnpft
       lai_bal_pft(l,n) = (a_ws(n)*eta_sl(n)*ht(l,n)                          &
                       /a_wl(n))**(1.0/(b_wl(n)-1))
       IF(l_trait_phys) THEN
         leaf(l,n) = cmass*lma(n)*lai_bal_pft(l,n)
       ELSE
         leaf(l,n) = sigl(n)*lai_bal_pft(l,n)
       END IF
       root(l,n) = leaf(l,n)
       wood(l,n) = a_wl(n)*(lai_bal_pft(l,n)**b_wl(n))
       c_veg(l,n) = leaf(l,n) + root(l,n) + wood(l,n)
    END DO
    cnsrv_veg_triffid_gb(l)=SUM(frac(l,1:npft)*c_veg(l,1:npft))

    cnsrv_soil_triffid_gb(l)=SUM(cs(l,:))

    cnsrv_prod_triffid_gb(l)=wood_prod_fast_gb(l)+wood_prod_med_gb(l)+        &
                          wood_prod_slow_gb(l)

    cnsrv_carbon_triffid_gb(l)=cnsrv_veg_triffid_gb(l)+                       &
                 cnsrv_soil_triffid_gb(l)+cnsrv_prod_triffid_gb(l)

    cnsrv_veg_flux(l)=SUM(npp(l,1:npft)*frac(l,1:npft))/gamma

    cnsrv_soil_flux(l)=-(1.0-resp_frac(l))*SUM(resp_s(l,1:4))/gamma

    cnsrv_prod_flux(l)=0.0

    cnsrv_carbon_flux(l)=SUM(npp(l,1:npft)*frac(l,1:npft))/gamma-             &
                     (1.0-resp_frac(l))*SUM(resp_s(l,1:4))/gamma
  END DO

!----------------------------------------------------------------------
! Loop through Functional Types
!----------------------------------------------------------------------
! Make sure wood product fluxes are initialised to zero

n_loss_gb(:)=0.

DO l=1,land_pts

   wp_fast_out_gb(l)=0.0
   wp_med_out_gb(l)=0.0
   wp_slow_out_gb(l)=0.0

   wp_fast_in_gb(l)=0.0
   wp_med_in_gb(l)=0.0
   wp_slow_in_gb(l)=0.0

   IF (.NOT.(l_landuse)) frac_agr_prev_gb(l)=frac_agric(l)
   IF (.NOT.(l_landuse)) frac_past_prev_gb(l)=frac_past(l)
END DO

n_veg_pft(:,:)=0.
dnveg_pft(:,:)=0.
dnveg_gb(:)=0.
dcveg_gb(:)=0.
n_demand_gb(:)=0.
n_uptake_gb(:)=0.
n_fix_gb(:)=0.
litterC_pft(:,:)=0.
litterN_pft(:,:)=0.

DO t=1,trif_pts
  l=trif_index(t)
  deposition(l) = MAX(deposition_N_gb(l) * 60. * 60. * 24. * 360.,0.)
    ! convert from kgN/m2/s to kgN/m2/360days
  DO n=1,npft
    n_fix_pft(l,n) = MAX(npp(l,n) * 1.6/1.e3,0.)
                    ! Fixation rate 1.6 gN(kg C)-1
    N_fix_gb(l) = N_fix_gb(l) + frac(l,n) * n_fix_pft(l,n)
  ENDDO
!----------------------------------------------------------------------
! Make deposition available for vegetation uptake
!----------------------------------------------------------------------
  n_inorg_gb(l) = n_inorg_gb(l) + (deposition(l) + N_fix_gb(l)) / gamma
ENDDO

!----------------------------------------------------------------------
! Loop through TRIFFID points
!----------------------------------------------------------------------
DO n=1,npft
  DO t=1,trif_pts
    l=trif_index(t)

    dnveg_pft(l,n)=0.

!----------------------------------------------------------------------
! Diagnose the balanced-growth leaf area index and the associated leaf,
! wood, root and total vegetation carbon
!----------------------------------------------------------------------
    lai_bal_pft(l,n) = (a_ws(n)*eta_sl(n)*ht(l,n)                             &
              /a_wl(n))**(1.0/(b_wl(n)-1))
    IF(l_trait_phys) THEN
      leaf(l,n) = cmass*lma(n)*lai_bal_pft(l,n)
    ELSE
      leaf(l,n) = sigl(n)*lai_bal_pft(l,n)
    END IF
    root(l,n) = leaf(l,n)
    wood(l,n) = a_wl(n) * (lai_bal_pft(l,n)**b_wl(n))
    c_veg(l,n) = leaf(l,n) + root(l,n) + wood(l,n)
!----------------------------------------------------------------------
! Diagnose the phenological state
!----------------------------------------------------------------------
    phen(l,n) = lai(l,n)/lai_bal_pft(l,n)
    phen(l,n) = MIN(1.0, phen(l,n))
    phen(l,n) = MAX(0.0, phen(l,n))

!----------------------------------------------------------------------
! Diagnose Nitrogen Pools
!----------------------------------------------------------------------
    CALL calc_n_comps_triffid(l, n, phen(l,n), lai_bal_pft(l,n), wood(l,n),    &
                              root(l,n), n_leaf, n_root, n_stem)

    n_veg_pft(l,n) = n_leaf + n_root + n_stem

    dnveg_pft(l,n) = n_veg_pft(l,n)
    dnveg_gb(l) = dnveg_gb(l) - frac(l,n) * n_veg_pft(l,n)
    dcveg_gb(l) = dcveg_gb(l) - frac(l,n) * c_veg(l,n)
  END DO

!----------------------------------------------------------------------
! Assume Inorganic N equally spread amongst PFTs and bare ground
!----------------------------------------------------------------------
  n_uptake_pft(:,n) = 0.0
  n_demand_pft(:,n) = 0.0
  n_avail(:,n) = n_inorg_gb(:)
!----------------------------------------------------------------------
! Update vegetation carbon contents
!----------------------------------------------------------------------
! DEPENDS ON: vegcarb
  CALL vegcarb (land_pts,trif_pts,trif_index,n,forw                           &
,               GAMMA,g_leaf(:,n),npp(:,n),resp_w(:,n),phen(:,n)              &
,               leaf(:,n),root(:,n),wood(:,n),n_avail(:,n)                    &
,               dleaf_pft(:,n),droot_pft(:,n),dwood_pft(:,n)                  &
,               dcveg_pft(:,n),pc_s_pft(:,n)                                  &
,               n_demand_lit_pft(:,n),n_demand_spread_pft(:,n)                &
,               n_demand_growth_pft(:,n),n_uptake_growth_pft(:,n)             &
,               n_uptake_spread_pft(:,n),exudates_pft(:,n)                    &
,               n_fertiliser_pft(:,n))

!----------------------------------------------------------------------
! Diagnose N demand and N uptake
! Note N_uptake spreading updated in competition
! N uptake could be greater than n_inorg_gb at this stage due to need
! maintain min frac
!----------------------------------------------------------------------

  n_uptake_pft(:,n) = n_uptake_growth_pft(:,n) + n_demand_lit_pft(:,n) +      &
                  n_uptake_spread_pft(:,n)
  n_demand_pft(:,n) = n_demand_growth_pft(:,n) + n_demand_lit_pft(:,n) +      &
                  n_demand_spread_pft(:,n)

END DO

!-----------------------------------------------------------------------
! Diagnose the new value of Canopy Height, Leaf Area Index and Total
! Vegetation Carbon
!-----------------------------------------------------------------------
DO n=1,nnpft

  DO t=1,trif_pts
    l=trif_index(t)

    ht(l,n) = wood(l,n) / (a_ws(n) * eta_sl(n))                               &
            * (a_wl(n)/wood(l,n))**(1.0/b_wl(n))
!    lma replaces sigl
    IF (l_trait_phys) THEN
      lai_bal_pft(l,n) = leaf(l,n) / (lma(n)*cmass)
    ELSE
      lai_bal_pft(l,n) = leaf(l,n) / sigl(n)
    END IF
    lai(l,n) = phen(l,n) * lai_bal_pft(l,n)
    c_veg(l,n) = leaf(l,n) + root(l,n) + wood(l,n)
!----------------------------------------------------------------------
! Diagnose Nitrogen Pools
!----------------------------------------------------------------------
    CALL calc_n_comps_triffid(l, n, phen(l,n), lai_bal_pft(l,n), wood(l,n),   &
                              root(l,n), n_leaf, n_root, n_stem)

    n_veg_pft(l,n) = n_leaf + n_root + n_stem
    dnveg_pft(l,n) = n_veg_pft(l,n) - dnveg_pft(l,n)
!----------------------------------------------------------------------
! Reduce NPP by exudates. Important for Litterfall calculation
! Note this is the diagnostic from TRIFFID
!----------------------------------------------------------------------
    npp(l,n) = npp(l,n) - exudates_pft(l,n)
  END DO
END DO


!----------------------------------------------------------------------
! If l_veg_compete=T, then there is the choice between the old co-competition
! (lotka_jls.F90)
! or the new purely height-based competition (l_ht_compete=T, recommended).
! If the user is not using the original 5 PFTs (BT, NT, C3, C4, SH),
! l_ht_compete must be True.
! If l_ht_compete=T, there are separate subroutines for dynamic (lotka_noeq)
! and equilibrium mode (lotka_eq).
!----------------------------------------------------------------------

IF ( l_veg_compete ) THEN
  IF ( l_trif_crop ) THEN
! Loop over vegetation types: natural, crop, pasture (,forestry?)
    DO k=1,veg_types
      crop_lotka=k-1.0
      IF ( crop_lotka == 0.0 ) THEN
        DO l=1,land_pts
          subset_space(l)=frac_agric(l)+frac_past(l)
          subset_space_prev(l)=frac_agr_prev_gb(l)+frac_past_prev_gb(l)
        END DO
      ELSE IF ( crop_lotka == 1.0 ) THEN
        DO l=1,land_pts
          subset_space(l)=frac_agric(l)
          subset_space_prev(l)=frac_agr_prev_gb(l)
        END DO
      ELSE IF ( crop_lotka == 2.0 ) THEN
        DO l=1,land_pts
          subset_space(l)=frac_past(l)
          subset_space_prev(l)=frac_past_prev_gb(l)
        END DO
      END IF
      nsub=0
      DO n=1,npft
        IF ( crop(n) == crop_lotka ) THEN
          nsub=nsub+1
        END IF
      END DO
      IF ( nsub > 0 ) THEN
! DEPENDS ON: lotka_noeq_subset
        CALL lotka_noeq_subset( land_pts,trif_pts,trif_index                  &
,                 c_veg,subset_space,subset_space_prev                        &
,                 GAMMA                                                       &
,                 ht,pc_s_pft,frac,dfrac,frac_na,dfrac_na                     &
,                 nsub,crop_lotka)
      END IF
    END DO
  ELSE
    IF ( l_ht_compete ) THEN
      IF ( l_trif_eq ) THEN
! DEPENDS ON: lotka_eq
        frac_prev(:,:) = frac(:,:)   !For diagnosing bare soil
        CALL lotka_eq (land_pts,trif_pts,trif_index                           &
,                   c_veg,frac_prev,frac_agric,ht,pc_s_pft                    &
,                   frac)

        frac_na(:,:) =frac(:,:)
        dfrac_na(:,:)=0.0
        dfrac(:,:)   =0.0
      ELSE
! DEPENDS ON: lotka_noeq
        CALL lotka_noeq( land_pts,trif_pts,trif_index                         &
,                      c_veg,frac_agric,frac_agr_prev_gb                      &
,                      GAMMA,lai_bal_pft                                      &
,                      ht,pc_s_pft,frac,dfrac,frac_na,dfrac_na)
      END IF
    ELSE
! DEPENDS ON: lotka
      CALL lotka (land_pts,trif_pts,trif_index                                &
,                 c_veg,forw,frac_agric,frac_agr_prev_gb                      &
,                 GAMMA,lai_bal_pft,pc_s_pft                                  &
,                 frac,dfrac,frac_na,dfrac_na)
    END IF
  END IF
ELSE
  dfrac(:,:) = 0.0
  dfrac_na(:,:) = 0.0
  frac_na(:,:) = frac(:,:)
END IF

!----------------------------------------------------------------------
! Calculate the fertiliser requirement at the gridbox level
!----------------------------------------------------------------------
DO t=1,trif_pts
  l=trif_index(t)
  n_fertiliser_gb(l)=0.
    DO n=1,npft
      frac_flux = frac(l,n) - (1.0-forw) * dfrac(l,n)
      n_fertiliser_gb(l) = n_fertiliser_gb(l) +                               &
                           frac_flux * n_fertiliser_pft(l,n)
    END DO
END DO

!----------------------------------------------------------------------
! Limit the N uptake to the size of pool. Additional demand is to
! maintain frac_min and LAI_min will be taken as negative litter N
!----------------------------------------------------------------------

DO t=1,trif_pts
  l=trif_index(t)

  n_unmet_gb(l)=0.0
  n_demand_gb(l)=0.0
  n_uptake_gb(l)=0.0
  exudates_gb(l)=0.0

  n_inorg_gb(l)=n_inorg_gb(l) + (n_fertiliser_gb(l) / gamma)

    DO n=1,npft
      n_unmet(l,n)=0.

      frac_flux = frac(l,n) - (1.0-forw) * dfrac(l,n)

!----------------------------------------------------------------------
! If required N uptake to keep frac_min is in excess of n_inorg_gb on tiles
! force the model to meet this through negative N uptake at the PFT level
!----------------------------------------------------------------------
      IF (l_nitrogen .AND. n_uptake_pft(l,n)/gamma > n_inorg_gb(l)) THEN
        n_unmet(l,n) = n_uptake_pft(l,n) - (n_inorg_gb(l)*gamma)
        n_uptake_pft(l,n) = n_inorg_gb(l)*gamma
      END IF

!Update diagnostics on fracs prior to competition
    n_demand_gb(l) = n_demand_gb(l) + frac_flux * n_demand_pft(l,n)
    n_uptake_gb(l) = n_uptake_gb(l) + frac_flux * n_uptake_pft(l,n)
    exudates_gb(l) = exudates_gb(l) + frac_flux * exudates_pft(l,n)
    n_unmet_gb(l)  = n_unmet_gb(l)  + frac_flux * n_unmet(l,n)
  END DO

!----------------------------------------------------------------------
! Update Inorganic N with plant uptake
!----------------------------------------------------------------------
   n_inorg_gb(l)=n_inorg_gb(l)-n_uptake_gb(l)/GAMMA
END DO

!----------------------------------------------------------------------
! Diagnose the litterfall from the carbon balance of each vegetation
! type
!----------------------------------------------------------------------

DO t=1,trif_pts
  l=trif_index(t)
  lit_c_t(l) = 0.0
  lit_N_t_gb(l) = 0.0

  DO n=1,nnpft
    lit_c_ag_pft(l,n) = 0.0
    lit_c_orig_pft(l,n) = 0.0

    frac_flux=frac(l,n)-(1.0-forw)*dfrac(l,n)
    frac_flux_nat=frac_na(l,n)-(1.0-forw)*dfrac_na(l,n)
    lit_c_orig_pft(l,n) = npp(l,n)                                            &
        -GAMMA/frac_flux*(c_veg(l,n)*frac(l,n)                                &
        -(c_veg(l,n)-dcveg_pft(l,n))*(frac(l,n)-dfrac(l,n)))

    lit_n_orig_pft(l,n) = n_uptake_pft(l,n)                                   &
        -GAMMA/frac_flux*(n_veg_pft(l,n)*frac(l,n)                            &
        -(n_veg_pft(l,n)-dnveg_pft(l,n))*(frac(l,n)-dfrac(l,n)))

! ----------------------------------------------------------------------
! Diagnose the frac fluxes from 'natural' and frac_agric changes
! IF its a crop type all the change in frac should be considered 'natural'
    IF ( (crop(n) == 1) .AND. (.NOT. l_trif_crop) ) THEN
! ............................................................
! CASE 1 - CROPS
      lit_c(l,n) = lit_c_orig_pft(l,n) * frac_flux
      lit_c_ag_pft(l,n) = 0.0
      lit_n_pft(l,n) = lit_n_orig_pft(l,n) * frac_flux
      lit_n_ag_pft(l,n) = 0.0

    ELSE
! DFRAC>0.0 - expansion of frac in presence of Agric
      IF (frac(l,n) <  frac_na(l,n)) THEN
! ............................................................
! CASE 2 - CONTRACTION DUE TO LAND USE
        lit_c(l,n) = npp(l,n)                                                 &
             -(GAMMA/frac_flux_nat)*(c_veg(l,n)*frac_na(l,n)                  &
             -(c_veg(l,n)-dcveg_pft(l,n))*(frac_na(l,n)-dfrac_na(l,n)))


        lit_c_ag_pft(l,n) = npp(l,n)                                          &
             -(GAMMA/(frac_flux))*(c_veg(l,n)*frac(l,n)                       &
             -(c_veg(l,n)-dcveg_pft(l,n))*(frac(l,n)-dfrac(l,n)))

        lit_c(l,n)=lit_c(l,n)*frac_flux_nat
        lit_c_ag_pft(l,n)=lit_c_ag_pft(l,n)*frac_flux

        lit_c_ag_pft(l,n)=lit_c_ag_pft(l,n)-lit_c(l,n)

        lit_n_pft(l,n) = n_uptake_pft(l,n)                                    &
             -(GAMMA/frac_flux_nat)*(n_veg_pft(l,n)*frac_na(l,n)              &
             -(n_veg_pft(l,n)-dnveg_pft(l,n))*(frac_na(l,n)-dfrac_na(l,n)))
        lit_n_pft(l,n)=lit_n_pft(l,n)*frac_flux_nat
        lit_n_ag_pft(l,n)=lit_n_orig_pft(l,n)*frac_flux
        lit_n_ag_pft(l,n)=lit_n_ag_pft(l,n)-lit_n_pft(l,n)

        IF (lit_c_ag_pft(l,n) <  0.0) THEN
          lit_c_ag_pft(l,n)=0.0
          lit_c(l,n)=lit_c_orig_pft(l,n) * frac_flux

          lit_n_ag_pft(l,n)=0.0
          lit_n_pft(l,n)=lit_n_orig_pft(l,n) * frac_flux
        END IF
      ELSE
! ............................................................
! CASE 3 - NON-CROP EXPANSION

        lit_c(l,n) = lit_c_orig_pft(l,n) * frac_flux
        lit_c_ag_pft(l,n) = 0.0

        lit_n_pft(l,n) = lit_n_orig_pft(l,n) * frac_flux
        lit_n_ag_pft(l,n) = 0.0
      END IF
    END IF

! FATE OF FRAC AGFOREST CLEARANCE CARBON

    IF (l_landuse) THEN

      cratio=root(l,n)/(leaf(l,n)+wood(l,n)+root(l,n))
!          CR=1.0  ! i.e. it all goes into the soil
!          CR=0.0  ! i.e. it all goes to wood products

      nratio = n_root/n_veg_pft(l,n)

      IF ( (l_trif_crop) ) THEN
        root_abandon_pft(l,n)=lit_c_ag_pft(l,n)*cratio/frac_flux
        lit_c_t(l) = lit_c_t(l)                                               &
          +lit_c(l,n)                                                         &
          +root_abandon_pft(l,n) * frac_flux

        root_abandon_n_pft(l,n)=lit_n_ag_pft(l,n)*nratio/frac_flux
        lit_n_t_gb(l) = lit_n_t_gb(l)                                         &
          +lit_n_pft(l,n)                                                     &
          +root_abandon_n_pft(l,n) * frac_flux

! calculate the remaining carbon that goes to the 1-year turnover pool for
! wood products (not crop harvest).
        lit_c_ag_pft(l,n)=lit_c_ag_pft(l,n)-                                  &
                         root_abandon_pft(l,n) * frac_flux

        lit_n_ag_pft(l,n)=lit_n_ag_pft(l,n)-                                  &
                         root_abandon_n_pft(l,n) * frac_flux
      ELSE
        lit_c_t(l) = lit_c_t(l)                                               &
          +lit_c(l,n)                                                         &
          +lit_c_ag_pft(l,n)*cratio*(1.0-crop(n))

        lit_n_t_gb(l) = lit_n_t_gb(l)                                         &
          +lit_n_pft(l,n)                                                     &
          +lit_n_ag_pft(l,n)*nratio*(1.0-crop(n))

        root_abandon_pft(l,n)=                                              &
                  lit_c_ag_pft(l,n)*cratio*(1.0-crop(n))/frac_flux

        root_abandon_n_pft(l,n)=                                              &
                  lit_n_ag_pft(l,n)*nratio*(1.0-crop(n))/frac_flux

! calculate the remaining carbon that goes to the 1-year turnover pool for
! wood products (not crop harvest).
        lit_c_ag_pft(l,n)=lit_c_ag_pft(l,n)*(1.0-cratio)*(1.0-crop(n))

        lit_n_ag_pft(l,n)=lit_n_ag_pft(l,n)*(1.0-nratio)*(1.0-crop(n))
      END IF
    ELSE
      cratio=1.0
      nratio=1.0

      IF ( (l_trif_crop) ) THEN
        root_abandon_pft(l,n)=lit_c_ag_pft(l,n)*cratio/frac_flux
        lit_c_t(l) = lit_c_t(l)                                               &
          +lit_c(l,n)                                                         &
          +root_abandon_pft(l,n) * frac_flux

        root_abandon_n_pft(l,n)=lit_n_ag_pft(l,n)*nratio/frac_flux
        lit_n_t_gb(l) = lit_n_t_gb(l)                                         &
          +lit_n_pft(l,n)                                                     &
          +root_abandon_n_pft(l,n) * frac_flux
      ELSE
        lit_c_t(l) = lit_c_t(l)                                               &
          +lit_c(l,n)                                                         &
          +lit_c_ag_pft(l,n)*cratio*(1.0-crop(n))
        root_abandon_pft(l,n)=                                                &
                  lit_c_ag_pft(l,n)*cratio*(1.0-crop(n))/frac_flux

        lit_n_t_gb(l) = lit_n_t_gb(l)                                         &
          +lit_n_pft(l,n)                                                     &
          +lit_n_ag_pft(l,n)*nratio*(1.0-crop(n))
        root_abandon_n_pft(l,n)=                                              &
                  lit_n_ag_pft(l,n)*nratio*(1.0-crop(n))/frac_flux
      END IF

      lit_c_ag_pft(l,n)=0.0
      lit_n_ag_pft(l,n)=0.0
    END IF
  END DO
END DO


!----------------------------------------------------------------------
! Harvest crop carbon
!----------------------------------------------------------------------
IF ( l_trif_crop ) THEN
  DO t=1,trif_pts
    l=trif_index(t)
    harvest_gb(l)=0.0
    harvest_n_gb(l)=0.0
    DO n=1,npft
      IF ( (crop(n) == 1.0) .AND. (lit_c(l,n) > 0.0) ) THEN
        frac_flux = frac(l,n) - (1.0-forw)*dfrac(l,n)
        harvest_pft(l,n) = harvest_rate*lit_c(l,n)
        lit_c_ag_pft(l,n) = lit_c_ag_pft(l,n) + harvest_pft(l,n)
        lit_c(l,n) = lit_c(l,n) - harvest_pft(l,n)
        lit_c_t(l) = lit_c_t(l) - harvest_pft(l,n)
        harvest_gb(l)=harvest_gb(l)+harvest_pft(l,n)
!convert from [kg/(m2 of land)/(360 days)] to [kg/(m2 of PFT)/(360 days)]
        harvest_pft(l,n) = harvest_pft(l,n)/frac_flux
      ELSE
        harvest_pft(l,n) = 0.0
      END IF

!Harvest Nitrogen as well
      IF ( (crop(n) == 1.0) .AND. (lit_n_pft(l,n) > 0.0) ) THEN
        frac_flux = frac(l,n) - (1.0-forw)*dfrac(l,n)
        harvest_n_pft(l,n) = harvest_rate*lit_n_pft(l,n)
!         lit_n_ag_pft(l,n) = lit_n_ag_pft(l,n) + harvest_n_pft(l,n)
        lit_n_pft(l,n) = lit_n_pft(l,n) - harvest_n_pft(l,n)
        lit_n_t_gb(l) = lit_n_t_gb(l) - harvest_n_pft(l,n)
        harvest_n_gb(l)=harvest_n_gb(l)+harvest_n_pft(l,n)
!convert from [kg/(m2 of land)/(360 days)] to [kg/(m2 of PFT)/(360 days)]
        harvest_n_pft(l,n) = harvest_n_pft(l,n)/frac_flux
      ELSE
        harvest_n_pft(l,n) = 0.0
      END IF
    END DO
  END DO
END IF


!----------------------------------------------------------------------
! Check carbon conservation (2/3)
! Add the litter fluxes to the net carbon fluxes
!----------------------------------------------------------------------
  DO t=1,trif_pts
    l=trif_index(t)
    cnsrv_veg_flux(l)=cnsrv_veg_flux(l)-                                      &
                      (SUM(lit_c_ag_pft(l,:))+lit_c_t(l))/GAMMA

    cnsrv_soil_flux(l)=cnsrv_soil_flux(l)+lit_c_t(l)/GAMMA

    cnsrv_prod_flux(l)=cnsrv_prod_flux(l)+                                    &
                          SUM(lit_c_ag_pft(l,:))/GAMMA
  END DO


!----------------------------------------------------------------------
! Call SOIL_C to update the soil carbon and nitrogen content
!----------------------------------------------------------------------
! DEPENDS ON: soilcarb

CALL soilcarb (land_pts,trif_pts,trif_index                                   &
,              lit_c,frac,resp_frac                                           &
,              forw,GAMMA,lit_c_t,lit_n_t_gb,resp_s,cs,ns_gb,neg_n            &
,              implicit_resp_correction)

!----------------------------------------------------------------------
! Turnover Inorganic N pool
!----------------------------------------------------------------------

DO t=1,trif_pts
  l=trif_index(t)

  n_loss_gb(l) = max(0.,N_inorg_turnover * n_inorg_gb(l))
  N_inorg_gb(l) = N_inorg_gb(l) - n_loss_gb(l) / GAMMA

! This condition only occurs if soil pools are at min N content due to either
! pool exhaustion or unmet seed N demand

  IF (N_inorg_gb(l) < 0.) THEN
    n_loss_gb(l) = n_loss_gb(l) + (N_inorg_gb(l) * GAMMA)
    N_inorg_gb(l) = 0.
  END IF

END DO

!----------------------------------------------------------------------
! Diagnose the gridbox mean vegetation carbon and nitrogen
!----------------------------------------------------------------------
DO t=1,trif_pts
  l=trif_index(t)

  leafC_pft(L,:)   = leaf(L,:)
  rootc_pft(L,:)   = root(L,:)
  woodc_pft(L,:)   = wood(L,:)
  litterC_pft(l,:) = root_litC_pft(l,:) +                                     &
                     leaf_litc_pft(l,:) + wood_litc_pft(l,:)
  litterN_pft(l,:) = root_litn_pft(l,:) +                                     &
                     leaf_litn_pft(l,:) + wood_litn_pft(l,:)

  cv(l) = 0.0
  n_veg_gb(l) = 0.0
  DO n=1,nnpft
    cv(l) = cv(l) + frac(l,n) * c_veg(l,n)
    n_veg_gb(l) = n_veg_gb(l) + frac(l,n) * n_veg_pft(l,n)
  END DO
  dnveg_gb(l) = dnveg_gb(l) + n_veg_gb(l)
  dcveg_gb(l) = dcveg_gb(l) + cv(l)

END DO

DO l=1,land_pts
  frac_agr_prev_gb(l)=frac_agric(l)
  frac_past_prev_gb(l)=frac_past(l)
END DO

!----------------------------------------------------------------------
! Update the Wood Product Pools
!----------------------------------------------------------------------

! DEPENDS ON: woodprod
IF (forw == 0.0) THEN
     CALL woodprod ( trif_pts,land_pts,trif_index                             &
     ,               GAMMA                                                    &
     ,               lit_c_ag_pft                                             &
     ,               wood_prod_fast_gb,wood_prod_med_gb                       &
     ,               wood_prod_slow_gb                                        &
     ,               wp_fast_in_gb,wp_med_in_gb,wp_slow_in_gb                 &
     ,               wp_fast_out_gb,wp_med_out_gb,wp_slow_out_gb )

ELSE

     wood_prod_fast_gb(:)=0.0
     wood_prod_med_gb(:)=0.0
     wood_prod_slow_gb(:)=0.0
     wp_fast_in_gb(:)=0.0
     wp_med_in_gb(:)=0.0
     wp_slow_in_gb(:)=0.0
     wp_fast_out_gb(:)=0.0
     wp_med_out_gb(:)=0.0
     wp_slow_out_gb(:)=0.0

END IF

!----------------------------------------------------------------------
! Change the units of PFT litter fluxes for output
! Change units to kg/(m2 of PFT)/(360 days)
!            from kg/(m2 of land)/(360 days)
!----------------------------------------------------------------------
DO t=1,trif_pts
  l=trif_index(t)
  DO n=1,nnpft
    frac_flux=frac(l,n)-(1.0-forw)*dfrac(l,n)
    IF (frac(l,n) <  frac_na(l,n)) THEN
      frac_flux_nat=frac_na(l,n)-(1.0-forw)*dfrac_na(l,n)
      lit_c(l,n)=lit_c(l,n)/frac_flux_nat
      lit_c_ag_pft(l,n)=lit_c_ag_pft(l,n)/frac_flux
    ELSE
      lit_c(l,n)=lit_c(l,n)/frac_flux
    ENDIF
  END DO
END DO

!----------------------------------------------------------------------
! Check carbon conservation (3/3)
! Add the wood product decay flux to the net carbon fluxes.
! Calculate the error in carbon conservation as:
! error=(final store)-(initial store)-(net flux)
!      ="extra carbon in store not accounted for by fluxes"
!----------------------------------------------------------------------
  DO t=1,trif_pts
    l=trif_index(t)

    wp_total_out_gb(l)=                                                       &
             (wp_fast_out_gb(l)+wp_med_out_gb(l)+wp_slow_out_gb(l))/GAMMA

    wood_prod_total_gb(l) =                                                   &
             wood_prod_fast_gb(l)+wood_prod_med_gb(l)+wood_prod_slow_gb(l)

    cnsrv_prod_flux(l)=cnsrv_prod_flux(l)-wp_total_out_gb(l)

    cnsrv_carbon_flux(l)=cnsrv_carbon_flux(l)-wp_total_out_gb(l)

    cnsrv_veg_triffid_gb(l)=cv(l)-cnsrv_veg_triffid_gb(l)-                    &
                         cnsrv_veg_flux(l)+                                   &
                         exudates_gb(l)/GAMMA

    cnsrv_soil_triffid_gb(l)=SUM(cs(l,:))-cnsrv_soil_triffid_gb(l)-           &
                          cnsrv_soil_flux(l)-                                 &
                          implicit_resp_correction(l)

    cnsrv_prod_triffid_gb(l)=wood_prod_total_gb(l)-                           &
                          cnsrv_prod_triffid_gb(l)-                           &
                          cnsrv_prod_flux(l)

    cnsrv_carbon_triffid_gb(l)=(cv(l)+SUM(cs(l,:))+                           &
                             wood_prod_total_gb(l) )-                         &
                            cnsrv_carbon_triffid_gb(l)-                       &
                            cnsrv_carbon_flux(l)-                             &
                            implicit_resp_correction(l)+                      &
                            exudates_gb(l)/GAMMA
  END DO

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE triffid
