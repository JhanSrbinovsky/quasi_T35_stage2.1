! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine SOILCARB -----------------------------------------------
!
! Purpose : Updates carbon and nitrogen contents of the soil.
!
! -------------------------------------------------------------------
SUBROUTINE soilcarb (land_pts,trif_pts,trif_index                             &
,                    lit_c,frac,resp_frac                                     &
,                    forw,GAMMA,lit_c_t,lit_n_t,resp_s_pot,                   &
                     cs,ns_gb,neg_n,implicit_resp_correction)

USE jules_surface_types_mod
USE prognostics, ONLY: ns_pool_gb, n_inorg_gb
USE jules_vegetation_mod, ONLY: l_nitrogen
USE trif_vars_mod, ONLY: minl_n_gb, minl_n_pot_gb, immob_n_gb, immob_n_pot_gb,&
                         fn_gb, resp_s_diag_gb, resp_s_pot_diag_gb,           &
                         dpm_ratio_gb, n_gas_gb
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE jules_vegetation_mod, ONLY: l_nitrogen
USE jules_soil_mod, ONLY : cs_min


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
,l,t                        ! WORK Loop counters

REAL                                                                          &
 forw                                                                         &
                            ! IN Forward timestep weighting.
,GAMMA                                                                        &
                            ! IN Inverse timestep (/360days).
,lit_c_t(land_pts)                                                            &
                            ! IN Total carbon litter
!                                 !    (kg C/m2/360days).
,lit_n_t(land_pts)                                                            &
                            ! IN Total nitrogen litter
!                                 !    (kg N/m2/360days).
,lit_c(land_pts,npft)                                                         &
                            ! IN Carbon Litter (kg C/m2/360days).
,resp_s_pot(land_pts,5)                                                       &
                            ! INOUT Soil respiration
!                                 !    (kg C/m2/360days).
,neg_n(land_pts)                                                              &
                            ! OUT Negative N required to prevent ns<0.
                                  !    (kg N)
,resp_s(land_pts,5)                                                           &
                            ! OUT Soil respiration
!                                 !    (kg C/m2/360days).
,ns_gb(land_pts)                                                              &
                            ! OUT Total Soil N
                                  !    (kg N/m2)
,cs(land_pts,4)                                                               &
                            ! INOUT Soil carbon (kg C/m2).
!                                 !    The 4 soil C pools are
!                                 !    DPM, RPM, biomass and humus.
,dcs(land_pts,4)                                                              &
                            ! WORK Increment to the soil carbon
!                                 !      (kg C/m2).
,dpc_dcs(land_pts,4)                                                          &
                            ! WORK Rate of change of PC with
!                                 !      soil carbon (/360days).
,pc(land_pts,4)                                                               &
                            ! WORK Net carbon accumulation in
!                                 !      the soil (kg C/m2/360days).
,pn(land_pts,4)                                                               &
                            ! WORK Net nitrogen accumulation in
!                                 !      the soil (kg N/m2/360days).
,frac(land_pts,ntype)                                                         &
                            ! INOUT Fractional cover of each
,resp_frac(land_pts)                                                          &
                            ! respired fraction of RESP_S
,cn(land_pts,5)
                            ! WORK C:N ratios of pools

REAL implicit_resp_correction(land_pts)
                            ! OUT respiration carried to next triffid
                            !     timestep to account for applying
                            !     minimum soil carbon constraint

REAL, PARAMETER               :: soil_cn   = 10.0
REAL, PARAMETER               :: lit_cn    = 300.0
REAL, PARAMETER               :: nminl_gas = 0.01

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SOILCARB'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

!----------------------------------------------------------------------
! Note: The 4 soil carbon pools are  1 decomposable plant material,
! 2 resistant plant material, 3 biomass, 4 humus.
!----------------------------------------------------------------------

n_gas_gb(:)=0.
fn_gb(:)=1.

DO l=1,land_pts

  resp_s_pot(l,5) = resp_s_pot(l,1) + resp_s_pot(l,2) +                       &
                resp_s_pot(l,3) + resp_s_pot(l,4)

  cn(l,1) = cs(l,1) / ns_pool_gb(l,1)   ! C:N ratio of DPM - Prognostic
  cn(l,2) = cs(l,2) / ns_pool_gb(l,2)   ! C:N ratio of DPM - Prognostic
  cn(l,3) = soil_cn           ! C:N ratio of BIO - Fixed
  cn(l,4) = soil_cn           ! C:N ratio of HUM - Fixed

  cn(l,1) = MIN(MAX(cn(l,1),1.e-6),1000.)
  cn(l,2) = MIN(MAX(cn(l,2),1.e-6),1000.)
  cn(l,3) = MIN(MAX(cn(l,3),1.e-6),1000.)
  cn(l,4) = MIN(MAX(cn(l,4),1.e-6),1000.)

  cn(l,5) = SUM(cs(l,:)) / SUM(ns_pool_gb(l,:))

  minl_n_pot_gb(l,1) = resp_s_pot(l,1) / cn(l,1)
  minl_n_pot_gb(l,2) = resp_s_pot(l,2) / cn(l,2)
  minl_n_pot_gb(l,3) = resp_s_pot(l,3) / cn(l,3)
  minl_n_pot_gb(l,4) = resp_s_pot(l,4) / cn(l,4)

  minl_n_pot_gb(l,5) = minl_n_pot_gb(l,1) + minl_n_pot_gb(l,2) +              &
                    minl_n_pot_gb(l,3) + minl_n_pot_gb(l,4)

  immob_n_pot_gb(l,1) = (0.46 * resp_frac(l) * resp_s_pot(l,1) / cn(l,3))     &
                   + (0.54 * resp_frac(l) * resp_s_pot(l,1) / cn(l,4))
  immob_n_pot_gb(l,2) = (0.46 * resp_frac(l) * resp_s_pot(l,2) / cn(l,3))     &
                   + (0.54 * resp_frac(l) * resp_s_pot(l,2) / cn(l,4))
  immob_n_pot_gb(l,3) = (0.46 * resp_frac(l) * resp_s_pot(l,3) / cn(l,3))     &
                   + (0.54 * resp_frac(l) * resp_s_pot(l,3) / cn(l,4))
  immob_n_pot_gb(l,4) = (0.46 * resp_frac(l) * resp_s_pot(l,4) / cn(l,3))     &
                   + (0.54 * resp_frac(l) * resp_s_pot(l,4) / cn(l,4))

  immob_n_pot_gb(l,5) = immob_n_pot_gb(l,1) + immob_n_pot_gb(l,2) +           &
                     immob_n_pot_gb(l,3) + immob_n_pot_gb(l,4)

!----------------------------------------------------------------------
! If soil N demand greater than N_available then limit decay through fn_gb
!----------------------------------------------------------------------
  fn_gb(l)=1.0
  IF( l_nitrogen .AND.                                                        &
     ((immob_n_pot_gb(l,5) - minl_n_pot_gb(l,5)) / gamma > n_inorg_gb(l)) ) THEN

   fn_gb(l) = (((minl_n_pot_gb(l,3) + minl_n_pot_gb(l,4) -                    &
              immob_n_pot_gb(l,3) - immob_n_pot_gb(l,4)) / GAMMA) +           &
              n_inorg_gb(l)) /                                                &
           ((immob_n_pot_gb(l,1) - minl_n_pot_gb(l,1)) / gamma +              &
            (immob_n_pot_gb(l,2) - minl_n_pot_gb(l,2)) / gamma)

   fn_gb(l) = MIN(MAX(fn_gb(l), 0.0), 1.0)
  END IF

  minl_n_gb(l,1) = fn_gb(l) * minl_n_pot_gb(l,1)
  minl_n_gb(l,2) = fn_gb(l) * minl_n_pot_gb(l,2)
  minl_n_gb(l,3) = minl_n_pot_gb(l,3)
  minl_n_gb(l,4) = minl_n_pot_gb(l,4)

  immob_n_gb(l,1) = fn_gb(l) * immob_n_pot_gb(l,1)
  immob_n_gb(l,2) = fn_gb(l)*immob_n_pot_gb(l,2)
  immob_n_gb(l,3) = immob_n_pot_gb(l,3)
  immob_n_gb(l,4) = immob_n_pot_gb(l,4)

  resp_s(l,1) = fn_gb(l) * resp_s_pot(l,1)
  resp_s(l,2) = fn_gb(l) * resp_s_pot(l,2)
  resp_s(l,3) = resp_s_pot(l,3)
  resp_s(l,4) = resp_s_pot(l,4)

  resp_s(l,5)  = resp_s(l,1) + resp_s(l,2) +                                  &
                 resp_s(l,3) + resp_s(l,4)
  minl_n_gb(l,5)  = minl_n_gb(l,1) + minl_n_gb(l,2) +                         &
                 minl_n_gb(l,3) + minl_n_gb(l,4)
  immob_n_gb(l,5) = immob_n_gb(l,1) + immob_n_gb(l,2) +                       &
                 immob_n_gb(l,3) + immob_n_gb(l,4)

END DO

! calculate DPM:RPM ratio of input litter Carbon
! DEPENDS ON: dpm_rpm
CALL dpm_rpm(land_pts,trif_pts,trif_index,                                    &
             lit_c,dpm_ratio_gb)

DO t=1,trif_pts
  l = trif_index(t)
!----------------------------------------------------------------------
! Diagnose the net local nitrogen flux into the soil
!----------------------------------------------------------------------
  pn(l,1) = dpm_ratio_gb(l) * lit_n_t(l) - minl_n_gb(l,1)
  pn(l,2) = (1.-dpm_ratio_gb(l)) * lit_n_t(l) - minl_n_gb(l,2)
  pn(l,3) = 0.46*immob_n_gb(l,5) - minl_n_gb(l,3)
  pn(l,4) = 0.54*immob_n_gb(l,5) - minl_n_gb(l,4)

  ns_pool_gb(l,1) = ns_pool_gb(l,1) + pn(l,1) / gamma
  ns_pool_gb(l,2) = ns_pool_gb(l,2) + pn(l,2) / gamma
  ns_pool_gb(l,3) = ns_pool_gb(l,3) + pn(l,3) / gamma
  ns_pool_gb(l,4) = ns_pool_gb(l,4) + pn(l,4) / gamma

  neg_n(l)=0.
  IF (ns_pool_gb(l,1) <  cs_min / lit_cn) THEN
    neg_n(l)  = neg_n(l) + ns_pool_gb(l,1) - cs_min / lit_cn
  END IF
  IF (ns_pool_gb(l,2) <  cs_min / lit_cn) THEN
    neg_n(l)  = neg_n(l) + ns_pool_gb(l,2) - cs_min / lit_cn
  END IF
  IF (ns_pool_gb(l,3) <  cs_min / soil_cn) THEN
    neg_n(l) = neg_n(l) + ns_pool_gb(l,3) - cs_min / soil_cn
  END IF
  IF (ns_pool_gb(l,4) <  cs_min / soil_cn) THEN
    neg_n(l) = neg_n(l) + ns_pool_gb(l,4) - cs_min / soil_cn
  END IF

  ns_pool_gb(l,1) = MAX(ns_pool_gb(l,1),cs_min / lit_cn)
  ns_pool_gb(l,2) = MAX(ns_pool_gb(l,2),cs_min / lit_cn)
  ns_pool_gb(l,3) = MAX(ns_pool_gb(l,3),cs_min / soil_cn)
  ns_pool_gb(l,4) = MAX(ns_pool_gb(l,4),cs_min / soil_cn)

! increase immobilisation to account for mininum n content
  immob_n_gb(l,5) = immob_n_gb(l,5) - (neg_n(l) * gamma)

! calculate mineralised gas emissions
  n_gas_gb(l) = nminl_gas * MAX((minl_n_gb(l,5) - immob_n_gb(l,5)),0.)

  !----------------------------------------------------------------------
! Update inorganic N
!----------------------------------------------------------------------
  n_inorg_gb(l) = n_inorg_gb(l) +                                             &
                  (minl_n_gb(l,5) - immob_n_gb(l,5) - n_gas_gb(l))            &
                  / gamma
  ns_gb(l)      = ns_pool_gb(l,1) + ns_pool_gb(l,2) +                         &
                  ns_pool_gb(l,3) + ns_pool_gb(l,4)

ENDDO

DO t=1,trif_pts
  l=trif_index(t)

!----------------------------------------------------------------------
! Diagnose the net local carbon flux into the soil
!----------------------------------------------------------------------
  pc(l,1) = dpm_ratio_gb(l) * lit_c_t(l) - resp_s(l,1)
  pc(l,2) = (1. - dpm_ratio_gb(l)) * lit_c_t(l) - resp_s(l,2)
  pc(l,3) = 0.46 * resp_frac(l) * resp_s(l,5) - resp_s(l,3)
  pc(l,4) = 0.54 * resp_frac(l) * resp_s(l,5) - resp_s(l,4)

!----------------------------------------------------------------------
! Variables required for the implicit and equilibrium calculations
!----------------------------------------------------------------------
  dpc_dcs(l,1) = resp_s(l,1) / cs(l,1)
  dpc_dcs(l,2) = resp_s(l,2) / cs(l,2)
  dpc_dcs(l,3) = resp_s(l,3) / cs(l,3)
  dpc_dcs(l,4) = resp_s(l,4) / cs(l,4)

!----------------------------------------------------------------------
! Save current value of soil carbon
!----------------------------------------------------------------------
  dcs(l,1) = cs(l,1)
  dcs(l,2) = cs(l,2)
  dcs(l,3) = cs(l,3)
  dcs(l,4) = cs(l,4)

END DO

!----------------------------------------------------------------------
! Update soil carbon
!----------------------------------------------------------------------
! DEPENDS ON: decay
CALL decay (land_pts,trif_pts,trif_index                                      &
,           dpc_dcs,forw,GAMMA,pc,cs)

!----------------------------------------------------------------------
! Apply implicit correction to the soil respiration rate.
!----------------------------------------------------------------------
DO t=1,trif_pts
  l=trif_index(t)

  dcs(l,1) = cs(l,1) - dcs(l,1)
  dcs(l,2) = cs(l,2) - dcs(l,2)
  dcs(l,3) = cs(l,3) - dcs(l,3)
  dcs(l,4) = cs(l,4) - dcs(l,4)

  resp_s(l,1) = resp_s(l,1) + forw * dpc_dcs(l,1) * dcs(l,1)
  resp_s(l,2) = resp_s(l,2) + forw * dpc_dcs(l,2) * dcs(l,2)
  resp_s(l,3) = resp_s(l,3) + forw * dpc_dcs(l,3) * dcs(l,3)
  resp_s(l,4) = resp_s(l,4) + forw * dpc_dcs(l,4) * dcs(l,4)

  implicit_resp_correction(l)=                                                &
                  (gamma*SUM(dcs(l,1:4))-SUM(pc(l,1:4))) / gamma
END DO


! sum total respiration
DO l=1,land_pts
  resp_s(l,5) = resp_s(l,1) + resp_s(l,2) +                                   &
                resp_s(l,3) + resp_s(l,4)

  resp_s_diag_gb(l,:) = resp_s(l,:)
  resp_s_pot_diag_gb(l,:) = resp_s_pot(l,:)

  resp_s_pot(l,1) = resp_s(l,1)
  resp_s_pot(l,2) = resp_s(l,2)
  resp_s_pot(l,3) = resp_s(l,3)
  resp_s_pot(l,4) = resp_s(l,4)

  resp_s_pot(l,5) = resp_s_pot(l,1) + resp_s_pot(l,2) +                       &
                resp_s_pot(l,3) + resp_s_pot(l,4)

END DO

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE soilcarb
