#if defined(UM_JULES)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Calls routines to initialize veg parameters and accumulated C fluxes
!
! Subroutine Interface
MODULE init_veg_mod
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='INIT_VEG_MOD'

CONTAINS
SUBROUTINE init_veg(a_step, triffid_period_arg, nstep_since_triffid)

!Use in relevant subroutines
USE init_urban_mod,           ONLY: init_urban
USE init_min_mod,             ONLY: init_min
USE init_acc_mod,             ONLY: init_acc

!USE in JULES modules
USE jules_surface_mod,        ONLY: l_aggregate
USE jules_surface_types_mod,  ONLY: npft, ntype
USE jules_vegetation_mod,     ONLY: l_nrun_mid_trif, l_triffid, l_phenol,     &
                                    triffid_period, can_model
USE switches_urban,           ONLY: l_urban2t

!USE in UM modules
USE nlsizes_namelist_mod,     ONLY: ntiles, land_field
USE submodel_mod,             ONLY: atmos_im
USE conversions_mod,          ONLY: rsec_per_day
USE model_time_mod,           ONLY: secs_per_stepim
USE atm_fields_mod,           ONLY: frac_typ, soil_carb1,                     &
                                    hgt, hwr, wrr, disp, ztm, albwl, albrd,   &
                                    emisw, emisr,                             &
                                    canht_pft, lai_pft, sat_soil_cond,        &
                                    z0m_soil, catch_snow, catch_tile,         &
                                    infil_tile, z0_tile, z0h_tile,            &
                                    npp_pft_acc, g_phlf_pft_acc,              &
                                    rsp_w_pft_acc, rsp_s_acc1, g_lf_pft_acc
USE umPrintMgr,               ONLY: umPrint, umMessage
USE yomhook,                  ONLY: lhook, dr_hook
USE parkind1,                 ONLY: jprb, jpim

IMPLICIT NONE

!
! Description:
!   Initializes vegetation parameters from fractions of surface types
!   and initializes accumulated carbon fluxes to zero if a new TRIFFID
!   calling period is starting.
!
! Method:
!   Calls routine SPARM to initialize vegetation parameters.
!   Calls routine INIT_ACC to initialize accumulated carbon fluxes.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Vegetation
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!

!Arguments
INTEGER, INTENT(IN)     :: a_step  !Current timestep in atmosphere model
INTEGER, INTENT(INOUT)  :: triffid_period_arg
INTEGER, INTENT(INOUT)  :: nstep_since_triffid

! #include "typduma.h"

!Local variables
INTEGER :: tile_pts(ntype)
  !Number of land points which include the nth surface type
INTEGER :: tile_index(land_field,ntype)
  !Indices of land points which include the nth surface type
INTEGER :: nstep_trif
  !Number of atmospheric timesteps between calls to TRIFFID.
INTEGER :: i,l,n
  !Counters

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='INIT_VEG'

!-----------------------------------------------------------------------
! If TRIFFID on, call INIT_MIN to ensure PFT fractions are GE minimum
! fraction except where vegetation excluded by ice, water or urban
!-----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (l_triffid) THEN
  CALL init_min(land_field,frac_typ,soil_carb1)
END IF

IF ( l_urban2t ) THEN
  ! This hijacks TILE 9 to be an urban_roof tile and so changes D1(JFRAC_TYP)
  ! It MUST therefore be called before call to TILEPTS.
  CALL init_urban(land_field, frac_typ, hgt, hwr, wrr, disp, ztm, albwl,      &
                  albrd, emisw, emisr)
END IF

!-----------------------------------------------------------------------
! Call TILEPTS to initialise TILE_PTS and TILE_INDEX
!-----------------------------------------------------------------------
! DEPENDS ON: tilepts
CALL tilepts(land_field,frac_typ,tile_pts,tile_index)

!-----------------------------------------------------------------------
! Initialise tiled and gridbox mean vegetation parameters
!-----------------------------------------------------------------------
! DEPENDS ON: sparm
CALL sparm (land_field, ntiles, can_model, l_aggregate, tile_pts, tile_index, &
            frac_typ, canht_pft, lai_pft, sat_soil_cond, z0m_soil, catch_snow,&
            catch_tile, infil_tile, z0_tile, z0h_tile)

IF (l_triffid) THEN
  !-----------------------------------------------------------------------
  ! If this is an NRUN and re-start from mid-way through a TRIFFID calling
  ! period has not been requested: (i) initialise accumulation prognostics
  ! to zero, (ii) set TRIFFID_PERIOD in integer header, and
  ! (iii) initialise ASTEPS_SINCE_TRIFFID integer header to zero.
  ! If mid-period restart is requested then leave the accumulated fields
  ! unchanged, and if a new calling period is specified then reset
  ! calling period header the new value provided that the number of
  ! atmosphere timesteps since the last call to TRIFFID does not exceed
  ! the new calling period .
  !-----------------------------------------------------------------------
  IF (a_step == 0) THEN
    IF (l_nrun_mid_trif) THEN

      IF (triffid_period /= triffid_period_arg) THEN
        nstep_trif = INT(rsec_per_day*triffid_period/secs_per_stepim(atmos_im))
        IF (nstep_since_triffid >  nstep_trif) THEN
          WRITE(umMessage,*) '**ERROR IN TRIFFID** YOU HAVE SELECTED TO'
          CALL umPrint(umMessage,src='init_veg')
          WRITE(umMessage,*) 'START MID-WAY THROUGH A TRIFFID CALLING'
          CALL umPrint(umMessage,src='init_veg')
          WRITE(umMessage,*) 'PERIOD BUT YOUR INITIAL DUMP CONTAINS'
          CALL umPrint(umMessage,src='init_veg')
          WRITE(umMessage,*) 'PROGNOSTICS ACCUMULATED OVER A PERIOD'
          CALL umPrint(umMessage,src='init_veg')
          WRITE(umMessage,*) 'LONGER THAN THE NEW CALLING PERIOD'
          CALL umPrint(umMessage,src='init_veg')
        ELSE
          triffid_period_arg = triffid_period
        END IF
      END IF

    ELSE

      CALL init_acc(land_field, npp_pft_acc, g_phlf_pft_acc, rsp_w_pft_acc,   &
                    rsp_s_acc1)
      triffid_period_arg = triffid_period
      nstep_since_triffid = 0

    END IF
  END IF
END IF

IF (l_phenol) THEN
  ! Initialise accumulated leaf turnover rate to zero
  g_lf_pft_acc(:,:) = 0.0
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE init_veg
END MODULE init_veg_mod
#endif
