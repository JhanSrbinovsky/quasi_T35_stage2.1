! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine JULES_INIT ----------------------------------
!
! Description: Initialisation of JULES.
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.
!
!   Code Owner: See Unified Model Code Owner's HTML page
!   This file belongs in section: Land

MODULE jules_init_mod
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='JULES_INIT_MOD'

CONTAINS

! For the scm the call takes a few extra arguments due to differences in
! module from the main UM
SUBROUTINE jules_init(land_index                                              &
#if defined(SCMA)
                      ,land_field, ntiles, sm_levels                          &
#endif
                      )

!Sort out scm-dependent USE statements
#if defined(SCMA)
USE s_main_force,             ONLY: rho_snow_grnd => i_rho_snow_grnd,         &
                                    snowdepth     => i_snowdepth,             &
                                    nsnow,                                    &
                                    ds            => i_ds,                    &
                                    sice          => i_sice,                  &
                                    sliq          => i_sliq
#else
USE nlsizes_namelist_mod,     ONLY: ntiles, land_field, sm_levels

USE atm_fields_mod,           ONLY: rho_snow_grnd,                            &
                                    snowdepth,                                &
                                    nsnow,                                    &
                                    ds,                                       &
                                    sice,                                     &
                                    sliq
#endif

!USE in relevant variables from UM modules
USE atm_fields_mod,           ONLY: clapp_horn, sat_soilw_suction,            &
                                    sat_soil_cond, therm_cap, therm_cond,     &
                                    vol_smc_crit, vol_smc_wilt, vol_smc_sat


USE dyn_coriolis_mod,         ONLY: f3_at_u
USE theta_field_sizes,        ONLY: t_i_length

!USE in relevant variables from JULES modules
USE p_s_parms,                ONLY: bexp_gb, sathh_gb, hcap_gb, hcon_gb,      &
                                    satcon_gb, smvccl_gb,                     &
                                    smvcwt_gb, smvcst_gb
USE prognostics,            ONLY: snowdepth_jules     => snowdepth_surft,     &
                                  rho_snow_grnd_jules => rho_snow_grnd_surft, &
                                  nsnow_jules         => nsnow_surft,         &
                                  ds_jules            => ds_surft,            &
                                  sice_jules          => sice_surft,          &
                                  sliq_jules          => sliq_surft
USE jules_surface_mod,        ONLY: l_aggregate, l_flake_model
USE jules_snow_mod,           ONLY: nsmax
USE lake_mod,                 ONLY: coriolis_param_gb, nusselt_gb, nusselt_0

USE parkind1,                 ONLY: jprb, jpim
USE yomhook,                  ONLY: lhook, dr_hook
USE jules_print_mgr,          ONLY: jules_message, jules_print

IMPLICIT NONE


! Arguments for scm
#if defined(SCMA)
INTEGER, INTENT(IN) :: land_field
INTEGER, INTENT(IN) :: ntiles
INTEGER, INTENT(IN) :: sm_levels
#endif
! Subroutine argument (have to order this way to ensure land_field is defined)
INTEGER, INTENT(IN) :: land_index    (MAX(1,land_field))


! WORK variables:
INTEGER :: i,j,l,n

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*),   PARAMETER :: RoutineName='JULES_INIT'

! START OF EXECUTABLE CODE

! Dimension the JULES fields.
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

#if ! defined(SCMA)
! Initialise JULES arrays.
! Soil properties on soil moisture levels.
DO l = 1, land_field
  DO n = 1, sm_levels
    bexp_gb(l,n)      = clapp_horn(l)
    sathh_gb(l,n)  = sat_soilw_suction(l)
    hcap_gb(l,n)   = therm_cap(l)
    smvccl_gb(l,n) = vol_smc_crit(l)
    smvcwt_gb(l,n) = vol_smc_wilt(l)
    smvcst_gb(l,n) = vol_smc_sat(l)
  END DO
  DO n = 0, sm_levels
    hcon_gb(l,n)   = therm_cond(l)
    satcon_gb(l,n) = sat_soil_cond(l)
  END DO
END DO
#endif

! snowdepth needed in AP1 for JULES radiation
!---------------------------------------------
  snowdepth_jules = snowdepth
  rho_snow_grnd_jules = rho_snow_grnd
  nsnow_jules         = nsnow
  IF (nsmax > 0)THEN
    ds_jules            = ds
    sice_jules          = sice
    sliq_jules          = sliq
  END IF

! FLake model
!--------------
IF (     l_flake_model                                                        &
    .AND.(.NOT.l_aggregate)) THEN

! initialise the Nusselt number
    nusselt_gb(:) = nusselt_0

    DO l=1,land_field

        j=(land_index(l)-1)/t_i_length + 1
        i = land_index(l) - (j-1)*t_i_length

! set the Coriolis parameter : ABSOLUTE VALUE
!
! To get the value at theta points,
! average the adjacent values at u points.
!
        coriolis_param_gb(l) = ABS( (f3_at_u(i,j)+f3_at_u(i-1,j))/2.0 )

    END DO

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE jules_init
END MODULE jules_init_mod

