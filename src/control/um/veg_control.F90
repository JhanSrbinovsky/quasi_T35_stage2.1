#if defined(UM_JULES)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Top-level control routine for vegetation section
!
MODULE veg_control_mod
  CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='VEG_CONTROL_MOD'

CONTAINS
SUBROUTINE veg_control(                                                       &
  land_pts, land_index, nsurft, can_model,                                    &
  a_step, asteps_since_triffid,                                               &
  land_pts_trif, npft_trif,                                                   &
  phenol_period, triffid_period, row_length, rows,                            &
  l_phenol, l_triffid, l_trif_eq,                                             &
  atimestep, frac_disturb, frac_past, satcon,                                 &
  g_leaf_acc, g_leaf_phen_acc, npp_acc,                                       &
  resp_s_acc, resp_w_acc,                                                     &
  cs, frac, lai, clay_gb, z0m_soil, ht,                                       &
  catch_s, catch_t, infil_t, z0_t, z0h_t,                                     &
  c_veg, cv, lit_c, lit_c_mn, g_leaf_day, g_leaf_phen,                        &
  lai_phen, g_leaf_dr_out, npp_dr_out, resp_w_dr_out,                         &
  resp_s_dr_out                                                               &
   )

!Module imports

!Common
  USE jules_surface_types_mod
  USE jules_vegetation_mod,     ONLY:                                         &
    i_veg_vn, i_veg_vn_1b, i_veg_vn_2b

!Modules that change name between JULES and UM
#if defined(UM_JULES)
  USE atm_step_local,           ONLY:                                         &
    dim_cs1
#else
  USE ancil_info,               ONLY:                                         &
    dim_cs1
#endif

  USE ereport_mod,              ONLY:                                         &
    ereport

  USE timestep_mod,             ONLY:                                         &
    timestep

  !Dr Hook
  USE yomhook,  ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  USE errormessagelength_mod, ONLY: errormessagelength

  IMPLICIT NONE

! Description: Controls calling of the vegetation code.
!              In due course the final UM-isms will be removed to allow this
!              to be called from JULES-standalone control too.

INTEGER, INTENT(IN) ::                                                        &
  land_pts,                                                                   &
  land_index(land_pts),                                                       &
  nsurft,                                                                     &
  can_model,                                                                  &
  a_step,                                                                     &
  phenol_period,                                                              &
  triffid_period,                                                             &
  row_length,                                                                 &
  rows,                                                                       &
  land_pts_trif,                                                              &
  npft_trif

REAL, INTENT(IN) ::                                                           &
  atimestep,                                                                  &
  frac_disturb(land_pts),                                                     &
  satcon(land_pts),                                                           &
  z0m_soil(land_pts)

LOGICAL, INTENT(IN) ::                                                        &
  l_phenol, l_triffid, l_trif_eq

INTEGER, INTENT(INOUT) ::                                                     &
  asteps_since_triffid

REAL, INTENT(INOUT) ::                                                        &
  g_leaf_acc(land_pts,npft),                                                  &
  g_leaf_phen_acc(land_pts,npft),                                             &
  npp_acc(land_pts_trif,npft_trif),                                           &
  resp_s_acc(land_pts_trif,dim_cs1),                                          &
  resp_w_acc(land_pts_trif,npft_trif),                                        &
  cs(land_pts,dim_cs1),                                                       &
  frac(land_pts,ntype),                                                       &
  lai(land_pts,npft),                                                         &
  frac_past(land_pts),                                                        &
  ht(land_pts,npft),                                                          &
  clay_gb(land_pts)

REAL, INTENT(OUT) ::                                                          &
  catch_s(land_pts,nsurft),                                                   &
  catch_t(land_pts,nsurft),                                                   &
  infil_t(land_pts,nsurft),                                                   &
  z0_t(land_pts,nsurft),                                                      &
  z0h_t(land_pts,nsurft),                                                     &
  c_veg(land_pts,npft),                                                       &
  cv(land_pts),                                                               &
  lit_c(land_pts,npft),                                                       &
  lit_c_mn(land_pts),                                                         &
  g_leaf_day(land_pts,npft),                                                  &
  g_leaf_phen(land_pts,npft),                                                 &
  lai_phen(land_pts,npft),                                                    &
  g_leaf_dr_out(land_pts,npft),                                               &
  npp_dr_out(land_pts,npft),                                                  &
  resp_w_dr_out(land_pts,npft),                                               &
  resp_s_dr_out(land_pts,dim_cs1+1)

  INTEGER :: phenol_call, triffid_call, nstep_trif

  CHARACTER(LEN=errormessagelength)       :: cmessage
  INTEGER                  :: errorstatus

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='VEG_CONTROL'

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  !-------------------------------------------------------------------------
  !   If leaf phenology is activated, check whether the atmosphere model
  !   has run an integer number of phenology calling periods.
  !-------------------------------------------------------------------------
  phenol_call=1
  IF (l_phenol) THEN
    phenol_call = MOD(FLOAT(a_step),(FLOAT(phenol_period)* (86400.0/timestep)))
  END IF

  triffid_call=1
  IF (l_triffid) THEN
    nstep_trif = INT(86400.0*triffid_period/timestep)
    IF (asteps_since_triffid == nstep_trif) THEN
      triffid_call = 0
    END IF
  END IF

  !Call to veg 1 or 2 as appriopriate
  IF ((phenol_call == 0).OR.(triffid_call == 0)) THEN
    SELECT CASE ( i_veg_vn )
      CASE ( i_veg_vn_2b )
! DEPENDS ON: veg2
        CALL veg2(                                                            &
          land_pts, land_index, nsurft, can_model,                            &
          a_step, asteps_since_triffid,                                       &
          phenol_period, triffid_period,                                      &
          l_phenol, l_triffid, l_trif_eq,                                     &
          atimestep, frac_disturb, frac_past, satcon,                         &
          g_leaf_acc, g_leaf_phen_acc, npp_acc,                               &
          resp_s_acc, resp_w_acc,                                             &
          cs, frac, lai, clay_gb, z0m_soil, ht,                               &
          catch_s, catch_t, infil_t, z0_t, z0h_t,                             &
          c_veg, cv, lit_c, lit_c_mn, g_leaf_day, g_leaf_phen,                &
          lai_phen, g_leaf_dr_out, npp_dr_out, resp_w_dr_out,                 &
          resp_s_dr_out                                                       &
          )

      CASE ( i_veg_vn_1b )
! DEPENDS ON: veg1
        CALL veg1(                                                            &
          land_pts,nsurft,can_model,                                          &
          a_step,phenol_period,l_phenol,                                      &
          atimestep,satcon, z0m_soil,                                         &
          g_leaf_acc,frac,lai,ht,                                             &
          catch_s,catch_t,infil_t,z0_t,z0h_t,                                 &
          g_leaf_day,g_leaf_phen,g_leaf_phen_acc,lai_phen                     &
          )

      CASE DEFAULT ! i_veg_vn
        errorstatus = 10
        WRITE (cmessage,'(A,A,I6)') 'Vegetation scheme version value',        &
               'i_veg_vn = ',i_veg_vn
        CALL Ereport ('VEG_CTL', errorstatus, cmessage)

    END SELECT ! i_veg_vn
  END IF

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END SUBROUTINE veg_control
END MODULE veg_control_mod
#endif
