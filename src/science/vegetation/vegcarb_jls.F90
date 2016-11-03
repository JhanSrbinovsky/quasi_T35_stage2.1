! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine VEGCARB ------------------------------------------------
!
! Purpose : Updates carbon contents of the vegetation.
!
! -------------------------------------------------------------------
 SUBROUTINE vegcarb (land_pts,trif_pts,trif_index                             &
,                    n,forw,GAMMA,g_leaf,npp,resp_w,phen                      &
,                    leaf,root,wood,n_avail,dleaf,droot,dwood                 &
,                    dcveg,pc_s,n_demand_lit,n_demand_spread                  &
,                    n_demand_growth,n_uptake_growth                          &
,                    n_uptake_spread,exudate,n_fertiliser)

USE pftparm
USE trif
USE jules_surface_mod, ONLY: cmass
USE jules_vegetation_mod, ONLY: l_trait_phys,l_Nitrogen,                      &
                                can_rad_mod, l_trif_crop
USE trif_vars_mod, ONLY:  root_litc_pft,leaf_litc_pft,wood_litc_pft,          &
                          root_litn_pft,leaf_litn_pft,wood_litn_pft

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
,n                                                                            &
                            ! IN Plant functional type.
,l,t                        ! WORK Loop counters

REAL                                                                          &
 forw                                                                         &
                            ! IN Forward timestep weighting.
,GAMMA                                                                        &
                            ! IN Inverse timestep (/360days).
,g_leaf(land_pts)                                                             &
                            ! IN Turnover rate for leaf and
!                                 !    fine root biomass (/360days).
,phen(land_pts)                                                               &
                            ! IN Phenological state.
,npp(land_pts)                                                                &
                            ! INOUT Net primary productivity
!                                 !       (kg C/m2/360days).
,resp_w(land_pts)                                                             &
                            ! INOUT Wood maintenance respiration
!                                 !       (kg C/m2/360days).
,leaf(land_pts)                                                               &
                            ! INOUT Leaf biomass (kg C/m2).
,root(land_pts)                                                               &
                            ! INOUT Root biomass (kg C/m2).
,wood(land_pts)                                                               &
                            ! INOUT Woody biomass (kg C/m2).
,n_avail(land_pts)                                                            &
                            ! IN Available Nitrogen (kg n/m2).
,dleaf(land_pts)                                                              &
                            ! OUT Leaf biomass (kg C/m2).
,droot(land_pts)                                                              &
                            ! OUT Root biomass (kg C/m2).
,dwood(land_pts)                                                              &
                            ! OUT Woody biomass (kg C/m2).
,dcveg(land_pts)                                                              &
                            ! OUT Change in vegetation carbon
!                                 !     during the timestep
!                                 !     (kg C/m2/timestep).
,n_uptake_growth(land_pts)                                                    &
                            ! OUT(kg N/m2/360day).
,n_uptake_spread(land_pts)                                                    &
                            ! OUT(kg N/m2/360day).
,n_avail_growth(land_pts)                                                     &
                            ! OUT(kg N/m2/360day).
,n_demand_growth(land_pts)                                                    &
                            ! OUT(kg N/m2/360day).
,n_fertiliser(land_pts)                                                       &
                            ! OUT(kg N/m2/360day).
,n_demand_lit(land_pts)                                                       &
                            ! WORK(kg N/m2/360day).
,exudate(land_pts)                                                            &
                            ! OUT(kg C/m2/360day).
,pc_s(land_pts)                                                               &
                            ! OUT Net carbon flux available
!                                 !     for spreading (kg C/m2/360days).
,dfpar_dlai                                                                   &
                            ! WORK Rate of change of FPAR
!                                 !      with leaf area index.
,dlai                                                                         &
                            ! WORK Increment to the leaf area
!                                 !      index.
,dlamg_dlai,dlit_dlai                                                         &
                            ! WORK Required for calculation
!                                 !      of the equilibrium increments.
,dnpp_dlai(land_pts)                                                          &
                            ! WORK Rate of change of NPP
!                                 !      with leaf area index
!                                 !      (kg C/m2/360days/LAI).
,dpc_dlai(land_pts)                                                           &
                            ! WORK Rate of change of PC
!                                 !      with leaf area index
!                                 !      (kg C/m2/360days/LAI).
,dpcg_dlai(land_pts)                                                          &
                            ! WORK Rate of change of PC_G
!                                 !      with leaf area index
!                                 !      (kg C/m2/360days/LAI).
,drespw_dlai                                                                  &
                            ! WORK Rate of change of RESP_W
!                                 !      with leaf area index
,fpar                                                                         &
                            ! WORK PAR interception factor.
,lai(land_pts)                                                                &
                            ! WORK Balanced leaf area index.
,lambda_g                                                                     &
                            ! WORK Fraction of NPP available
!                                 !      for spreading.
,lit_c_l(land_pts)                                                            &
                            ! WORK Local rate of Carbon Litter
!                                 !      production (kg C/m2/360days).
,pc(land_pts)                                                                 &
                            ! WORK Net carbon flux available
!                                 !      to vegetation (kg C/m2/360days)
,pc_g(land_pts)                                                               &
                            ! WORK Net carbon flux available
!                                 !      for growth (kg C/m2/360days).
,n_demand_spread(land_pts)                                                    &
                            ! WORK(kg N/m2/360day).
,n_demand_phenology(land_pts)                                                 &
                            ! WORK(kg N/m2/360day).
,n_root,n_leaf,n_stem,n_leaf0,n_leaf1                                         &
                            ! WORK N_pools (kg N)
,plant_cn(land_pts)                                                           &
                            ! WORK Plant C:N Ratio
,dphen(land_pts)                                                              &
                            ! WORK Change in Phen due to phenology
                                  ! since last TRIFFID call
,g_leaf_N                                                                     &
                            ! WORK leaf N turnover during growth
,root_cn                                                                      &
                            ! WORK Root C:N Ratio
,stem_cn                                                                      &
                            ! WORK Stem C:N Ratio
,leaf_cn                                                                      &
                            ! WORK Leaf C:N Ratio
,nl_bal
                            ! WORK Balanced LAI Leaf N
                                  ! Concentration kgN/kgC

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='VEGCARB'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

DO t=1,trif_pts
  l=trif_index(t)

! lma replaces sigl
  IF ( l_trait_phys ) THEN
    lai(l) = leaf(l)/(cmass*lma(n))
  ELSE
    lai(l) = leaf(l)/sigl(n)
  END IF

   IF (l_nitrogen) THEN
     ! Diagnose dphen for leaf turnover
     dphen(l)=g_leaf(l)/gamma

     g_leaf_n=0.
     IF (g_leaf(l) <  0.) THEN
       ! Background turnover only
       g_leaf(l)=g_leaf_0(n)*phen(l)
       g_leaf_n=g_leaf(l)
     ENDIF

   ENDIF

!----------------------------------------------------------------------
! Calculate the local production rate for carbon litter
!----------------------------------------------------------------------
  lit_c_l(l) = g_leaf(l)*leaf(l)+g_root(n)*root(l)                            &
               + g_wood(n)*wood(l)

! leaf N at full leaf
  CALL calc_n_comps_triffid(l, n, 1., lai(l), wood(l), root(l),               &
                            n_leaf1, n_root, n_stem)
! leaf N without phenological change
  CALL calc_n_comps_triffid(l, n, phen(l) + dphen(l), lai(l), wood(l), root(l),&
                            n_leaf0, n_root, n_stem)
! leaf N after phenology
  CALL calc_n_comps_triffid(l, n, phen(l), lai(l), wood(l), root(l),          &
                            n_leaf, n_root, n_stem)

  root_cn = root(l)/n_root
  stem_cn = wood(l)/n_stem
  leaf_cn = leaf(l)/n_leaf1
  plant_cn=(leaf(l)+root(l)+wood(l))/(n_leaf+n_root+n_stem)

  root_litc_pft(l,n)=g_root(n)*root(l)
  leaf_litc_pft(l,n)=g_leaf(l)*leaf(l)
  wood_litc_pft(l,n)=g_wood(n)*wood(l)

  root_litn_pft(l,n)=(1.-retran_r(n))*root_litc_pft(l,n)/root_cn
  leaf_litn_pft(l,n)=(1.-retran_l(n))*leaf_litc_pft(l,n)/leaf_cn
  wood_litn_pft(l,n)=wood_litc_pft(l,n)/stem_cn


!----------------------------------------------------------------------
! Diagnose the net local carbon flux into the vegetation
!----------------------------------------------------------------------
  pc(l) = npp(l) - lit_c_l(l)
!----------------------------------------------------------------------
! Variables required for the implicit and equilibrium calculations
!----------------------------------------------------------------------
  dlit_dlai = (g_leaf(l)*leaf(l)+g_root(n)*root(l))/lai(l)                    &
            + b_wl(n)*g_wood(n)*wood(l)/lai(l)

  fpar = (1 - EXP(-kpar(n)*lai(l)))/kpar(n)
  dfpar_dlai = EXP(-kpar(n)*lai(l))

  dnpp_dlai(l) = npp(l)*dfpar_dlai/fpar                                       &
               + (1-r_grow(n))*resp_w(l)                                      &
               *(dfpar_dlai/fpar-b_wl(n)/lai(l))

  lambda_g = 1 - (lai(l) - lai_min(n))                                        &
                /(lai_max(n) - lai_min(n))

  dlamg_dlai = -1.0/(lai_max(n) - lai_min(n))

  pc_g(l) = lambda_g * npp(l) - lit_c_l(l)
  dpcg_dlai(l) = lambda_g*dnpp_dlai(l)                                        &
               + dlamg_dlai*npp(l)                                            &
               - dlit_dlai
  dpc_dlai(l) = dnpp_dlai(l) - dlit_dlai

!----------------------------------------------------------------------
! Nitrogen required for litter loss
! Reduces N available for growth, possibly leading to negative growth
! in the absence of enough N
!----------------------------------------------------------------------

  n_demand_lit(l)=(root_litn_pft(l,n) + leaf_litn_pft(l,n) + wood_litn_pft(l,n))
  n_demand_phenology(l)=(n_leaf - n_leaf0) * GAMMA
  n_demand_lit(l) = n_demand_lit(l) + n_demand_phenology(l)

  n_avail(l)=n_avail(l)-n_demand_lit(l)/GAMMA
!----------------------------------------------------------------------
! Split available N between growing and spreading
!----------------------------------------------------------------------
  n_avail_growth(l)=lambda_g * n_avail(l)

END DO

!----------------------------------------------------------------------
! Update vegetation carbon contents
! Split available N between growing and spreading
! Excess C is considered exudate
!----------------------------------------------------------------------
DO t=1,trif_pts
  l=trif_index(t)
  dcveg(l) = leaf(l)+root(l)+wood(l)
END DO

! DEPENDS ON: plant_growth_n
CALL plant_growth_n (land_pts,trif_pts,trif_index                             &
,                  n,dpcg_dlai,forw,GAMMA,pc_g,lit_c_l,n_avail_growth         &
,                  phen,leaf,root,wood,dleaf,droot,dwood                      &
,                  n_demand_growth,n_uptake_growth,exudate,n_fertiliser)

DO t=1,trif_pts
  l=trif_index(t)
  dcveg(l) = leaf(l)+root(l)+wood(l)-dcveg(l)
  n_avail(l)=n_avail(l)-n_uptake_growth(l)/GAMMA
END DO

!----------------------------------------------------------------------
! Diagnose the carbon available for spreading and apply implicit
! corrections to the driving fluxes.
!----------------------------------------------------------------------
DO t=1,trif_pts
  l=trif_index(t)

!   lma replaces sigl
  IF ( l_trait_phys ) THEN
    dlai = leaf(l)/(cmass*lma(n)) - lai(l)
  ELSE
    dlai = leaf(l)/sigl(n) - lai(l)
  END IF
  pc_s(l) = pc(l) + forw*dpc_dlai(l)*dlai - dcveg(l)*GAMMA - exudate(l)

  CALL calc_n_comps_triffid(l, n, phen (l), lai(l), wood(l), root(l),         &
                            n_leaf, n_root, n_stem)

  root_cn = root(l)/n_root
  stem_cn = wood(l)/n_stem
  leaf_cn = leaf(l)/n_leaf
  plant_cn=(leaf(l)+root(l)+wood(l))/(n_leaf+n_root+n_stem)

!----------------------------------------------------------------------
! N demand for spreading a function of whole plant C:N ratio
!----------------------------------------------------------------------
  n_demand_spread(l)=pc_s(l)/plant_cn(l)
  n_uptake_spread(l)=n_demand_spread(l)

!----------------------------------------------------------------------
! Limit spreading if N limiting
!----------------------------------------------------------------------
  IF ( (l_trif_crop .AND. crop(n) == 1) .AND.                                 &
       (n_demand_spread(l)/GAMMA > n_avail(l)) ) THEN
    n_fertiliser(l)=n_fertiliser(l) +                                         &
                    ((n_demand_spread(l)/GAMMA - n_avail(l)) * GAMMA)
    n_avail(l)=n_demand_spread(l)/GAMMA
  END IF

  IF (l_nitrogen .and. n_demand_spread(l)/GAMMA > n_avail(l)) THEN
    n_uptake_spread(l)=n_avail(l)*gamma
    exudate(l)=exudate(l)+pc_s(l)
    pc_s(l)=n_uptake_spread(l)*plant_cn(l)
    exudate(l)=exudate(l)-pc_s(l)
  END IF

  drespw_dlai = resp_w(l)*b_wl(n)/lai(l)

  npp(l) = npp(l) + forw*dnpp_dlai(l)*dlai
  resp_w(l) = resp_w(l) + forw*drespw_dlai*dlai

END DO

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE vegcarb
