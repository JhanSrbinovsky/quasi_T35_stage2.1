! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine plant_growth_n -------------------------------------------------
!
! Purpose : Increments leaf, root and wood carbon.
! Andy Wiltshire March 2014
! -------------------------------------------------------------------
SUBROUTINE plant_growth_n (land_pts,trif_pts,trif_index                       &
,                  n,dpcg_dlai,forw,GAMMA,pc_g,lit_c_l,n_avail                &
,                  phen,leaf,root,wood,dleaf,droot,dwood                      &
,                  n_demand_growth,n_uptake_growth,exudate                    &
,                  n_fertiliser)

USE descent
USE pftparm
USE trif
USE jules_vegetation_mod, ONLY: l_trait_phys,l_nitrogen,can_rad_mod,          &
                                l_trif_crop
USE jules_surface_mod, ONLY: cmass
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
                      ! IN Vegetation type.
,l,t                        ! WORK Loop counters

REAL                                                                          &
 dpcg_dlai(land_pts)                                                          &
                      ! IN Rate of change of PC_G with
!                                 !    leaf area index
!                                 !    (kg C/m2/360days/LAI).
,forw                                                                         &
                      ! IN Forward timestep weighting.
,GAMMA                                                                        &
                      ! IN Inverse timestep (/360days).
,pc_g(land_pts)                                                               &
                      ! IN Net carbon flux available
!                                 !    for growth (kg C/m2/360days).
,lit_c_l(land_pts)                                                            &
                      ! IN Local rate of Carbon Litter
!                                 !      production (kg C/m2/360days).
,phen(land_pts)                                                               &
                            ! IN Phenological state.
,n_avail(land_pts)                                                            &
                      ! INOUT Plant Available Nitrogen (kg N/m2)
,leaf(land_pts)                                                               &
                      ! INOUT Leaf biomass (kg C/m2).
,root(land_pts)                                                               &
                      ! INOUT Root biomass (kg C/m2).
,wood(land_pts)                                                               &
                      ! INOUT Woody biomass (kg C/m2).
,exudate(land_pts)                                                            &
                      ! INOUT Exudate kgC/m2/360day
,n_uptake_growth(land_pts)                                                    &
                      ! OUT(kg N/m2/360day).
,n_demand_growth(land_pts)                                                    &
                      ! OUT(kg N/m2/360day).
,n_fertiliser(land_pts)
                      ! OUT(kg N/m2/360day).

REAL                                                                          &
 denom                                                                        &
                      ! WORK Denominator of update
!                                 !      equation.
,dleaf,droot,dwood                                                            &
                      ! WORK Increments to leaf, root
!                                 !      and woody biomass (kg C/m2).
,dl_dw                                                                        &
                      ! WORK Rate of change of leaf
!                                 !      carbon with wood carbon.
,dlai_dw                                                                      &
                      ! WORK Rate of change of leaf area
!                                 !      index with wood carbon
!                                 !      (LAI m2/kg C).
,dr_dw                                                                        &
                      ! WORK Rate of change of root
!                                 !      carbon with wood carbon.
,numer                                                                        &
                      ! WORK Numerator of the update
!                                 !      equation.
,wood_max                                                                     &
                      ! WORK Maximum wood carbon (kg C/m2).
,wood_min                                                                     &
                      ! WORK Minimum wood carbon (kg C/m2).
,n_veg(land_pts),n_veg_old(land_pts),                                         &
                      ! WORK Vegetation N content (kg n/m2)
lai,                                                                          &
                      ! WORK Balanced LAI
dleaf_pot,                                                                    &
                      ! WORK Unlimited inc leaf C (kgC m-2)
dwood_pot,                                                                    &
                      ! WORK Unlimited inc wood leaf C (kgC m-2)
droot_pot,                                                                    &
                      ! WORK Unlimited inc root leaf C (kgC m-2)
nl_bal,                                                                       &
                      ! WORK Balanced LAI leaf N Concentration
                          ! kgN/kgC
n_leaf,                                                                       &
                      ! WORK leaf N (kgN m-2)
n_root,                                                                       &
                      ! WORK root N (kgN m-2)
n_stem ,                                                                      &
                      ! WORK stem N (kgN m-2)
nit_pot,                                                                      &
                      ! WORK Unlimited N demand (kgN m-2)
x1,                                                                           &
                      ! WORK minimum bracket of the root
x2,                                                                           &
                      ! WORK maximum bracket of the root
rtbis,                                                                        &
                      ! WORK root of the bisection
dx,                                                                           &
                      ! WORK bisection
fmid,                                                                         &
                      ! WORK bisection
xmid,                                                                         &
                      ! WORK bisection
nl_bal_pot,                                                                   &
                      ! WORK Balanced LAI leaf N Concentration
                          ! kgN/kgC
n_plant_pot,                                                                  &
                      ! WORK unlimited plant N (kgN m-2)
n_root_pot,                                                                   &
                      ! WORK unlimited root N (kgN m-2)
n_stem_pot,                                                                   &
                      ! WORK unlimited stem N (kgN m-2)
n_leaf_pot,                                                                   &
                      ! WORK unlimited leaf N (kgN m-2)
lai_pot
                      ! WORK unlimited Balanced LAI

INTEGER                                                                       &
it  !number of tries towards solution

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PLANT_GROWTH_N'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

exudate(:)=0.0
n_fertiliser(:)=0.0

DO t=1,trif_pts
  l=trif_index(t)

!----------------------------------------------------------------------
! Calculate the unlimited increment to the wood carbon
!----------------------------------------------------------------------
  dl_dw = leaf(l)/(b_wl(n)*wood(l))
  dr_dw = dl_dw
  IF(l_trait_phys) THEN
    dlai_dw = dl_dw/(cmass * lma(n))
           !gC/gleaf * gleaf/m2
  ELSE
    dlai_dw = dl_dw/sigl(n)
  END IF

  numer = pc_g(l)
  denom = (1+dl_dw+dr_dw)*GAMMA-forw*dlai_dw*dpcg_dlai(l)
  denom = MAX(denom,denom_min)

  dwood_pot = numer/denom
!----------------------------------------------------------------------
! Ensure that the local leaf area index does not drop below its
! minimum value or exceed its maximum value.
!----------------------------------------------------------------------
  wood_min = a_wl(n)*lai_min(n)**b_wl(n)
  wood_max = a_wl(n)*lai_max(n)**b_wl(n)

  dwood_pot = MAX((wood_min-wood(l)),dwood_pot)
  dwood_pot = MIN((wood_max-wood(l)),dwood_pot)

!----------------------------------------------------------------------
! Diagnose the unlimited increments to leaf and root carbon
!----------------------------------------------------------------------
  IF( l_trait_phys ) THEN
    dleaf_pot = (cmass*lma(n)) *                                              &
                ((wood(l)+dwood_pot)/a_wl(n))**(1.0/b_wl(n)) - leaf(l)
    lai=leaf(l)/(cmass*lma(n))   ! Balanced LAI
  ELSE
    dleaf_pot = sigl(n)*((wood(l)+dwood_pot)/a_wl(n))**(1.0/b_wl(n))          &
               -leaf(l)
    lai=leaf(l)/sigl(n)          ! Balanced LAI
  END IF
  droot_pot = dleaf_pot

  CALL calc_n_comps_triffid(l, n, phen(l), lai, wood(l), root(l),             &
                            n_leaf, n_root, n_stem)

  n_veg(l)=n_leaf+n_root+n_stem

  IF( l_trait_phys ) THEN
    lai_pot = (leaf(l)+dleaf_pot) / (cmass*lma(n))
  ELSE
    lai_pot = (leaf(l)+dleaf_pot)/ sigl(n)
  END IF
  CALL calc_n_comps_triffid(l, n,phen(l), lai_pot, wood(l)+dwood_pot,         &
                            root(l)+droot_pot, n_leaf_pot,                    &
                            n_root_pot, n_stem_pot)

  n_plant_pot=n_leaf_pot+n_root_pot+n_stem_pot

  nit_pot= n_plant_pot-n_veg(l)
  n_demand_growth(l)=nit_pot * GAMMA
  n_uptake_growth(l)=n_demand_growth(l)

  IF ( (l_trif_crop .AND. crop(n) == 1) .AND. (nit_pot > n_avail(l)) ) THEN
    n_fertiliser(l) = nit_pot-n_avail(l)
    n_avail(l)      = nit_pot
    n_fertiliser(l) = n_fertiliser(l) * GAMMA
  END IF

!---------------------------------------------------------------------
! So, if demand is greater than availbility we need to limit growth.
! Veg can only lose as much C as local lit_c - all NPP used as exudates
! More C can be lost through contration
! N Limitation does not force the plant to use C as exudates
!---------------------------------------------------------------------
    IF (nit_pot >  n_avail(l) .AND. l_nitrogen) THEN

    x1=-1.*lit_c_l(l)
    x2=pc_g(l)  !maximum bracket of the root
    rtbis=x1    !root of the bisection
    dx=x2-x1

    it=0  !number of tries towards solution
    fmid=EPSILON(1.0)+1.

     DO WHILE (ABS(fmid) >  EPSILON(1.0) .AND. it <= 20)
!   Abort search if >20 iterations needed to find solution
       it=it+1

       dx=dx*0.5
       xmid=rtbis+dx

       numer = xmid
       denom = (1+dl_dw+dr_dw)*GAMMA-forw*dlai_dw*dpcg_dlai(l)
       denom = MAX(denom,denom_min)

       dwood = numer/denom
!---------------------------------------------------------------------
! Ensure that the local leaf area index does not drop below its
! minimum value or exceed its maximum value.
!---------------------------------------------------------------------
       wood_min = a_wl(n)*lai_min(n)**b_wl(n)
       wood_max = a_wl(n)*lai_max(n)**b_wl(n)

       dwood = MAX((wood_min-wood(l)),dwood)
       dwood = MIN((wood_max-wood(l)),dwood)

!---------------------------------------------------------------------
! Diagnose the increments to leaf and root carbon
!---------------------------------------------------------------------
       IF (l_trait_phys) THEN
         dleaf = cmass*lma(n)*                                                &
                 ((wood(l)+dwood)/a_wl(n))**(1.0/b_wl(n)) -leaf(l)
         lai_pot = (leaf(l)+dleaf)/(cmass*lma(n))
       ELSE
         dleaf = sigl(n)*((wood(l)+dwood)/a_wl(n))**(1.0/b_wl(n)) -leaf(l)
         lai_pot = (leaf(l)+dleaf)/sigl(n)
       END IF
       droot = dleaf

!calculate the Nitrogen required to satisfy the demand in new growth
       lai_pot = (leaf(l)+dleaf)/sigl(n)

       CALL calc_n_comps_triffid(l, n, phen(l), lai_pot, wood(l)+dwood,       &
                                 root(l)+droot, n_leaf_pot,                   &
                                 n_root_pot, n_stem_pot)

       n_plant_pot =n_leaf_pot+n_root_pot+n_stem_pot
       nit_pot = n_plant_pot-n_veg(l)
       fmid = nit_pot - n_avail(l)

       IF (fmid <  0.0) rtbis=xmid

    END DO !bisection

!---------------------------------------------------------------------
! Update carbon contents
!---------------------------------------------------------------------
      leaf(l) = leaf(l)+dleaf
      root(l) = root(l)+droot
      wood(l) = wood(l)+dwood

! put the excess C into exudates. This is added directly to Soil Respiration
      exudate(l)=(dleaf_pot+droot_pot+dwood_pot) -                            &
      (dleaf+droot+dwood)

    ELSE
      dleaf=dleaf_pot
      droot=droot_pot
      dwood=dwood_pot

!----------------------------------------------------------------------
! Update carbon contents
!----------------------------------------------------------------------
      leaf(l) = leaf(l)+dleaf_pot
      root(l) = root(l)+droot_pot
      wood(l) = wood(l)+dwood_pot
      exudate(l)=0.0

    END IF

  IF( l_trait_phys ) THEN
    lai=leaf(l)/(cmass * lma(n))
  ELSE
    lai=leaf(l)/sigl(n)
  END IF
  CALL calc_n_comps_triffid(l, n, phen(l), lai, wood(l),                      &
                            root(l), n_leaf,                                  &
                            n_root, n_stem)

  n_veg_old(l)=n_veg(l)
  n_veg(l)=n_leaf+n_root+n_stem

  n_avail(l)=n_avail(l) -(n_veg(l)-n_veg_old(l))
  n_uptake_growth(l)=(n_veg(l)-n_veg_old(l))*GAMMA

  exudate(l)=exudate(l)*GAMMA
END DO

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE plant_growth_n

