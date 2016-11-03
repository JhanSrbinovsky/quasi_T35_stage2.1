#if !defined(UM_JULES)
! *****************************COPYRIGHT********************************************
! (c) Crown copyright, Met Office 2011. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use
! and distribution under the JULES collaboration agreement, subject
! to the terms and conditions set out therein.
!
! [Met Office Ref SC0237]
! *****************************COPYRIGHT********************************************

SUBROUTINE update_irrig_variables ()

 USE logging_mod, ONLY : log_info

 USE model_time_mod, ONLY : current_time

 USE datetime_mod, ONLY : datetime_to_string

 USE ancil_info, ONLY : land_pts, frac_surft

 USE p_s_parms, ONLY : sthu_gb

 USE jules_surface_types_mod, ONLY : ncpft, nnpft

 USE jules_vegetation_mod, ONLY : irr_crop

 USE crop_vars_mod, ONLY :                                                    &
   frac_irr_surft, frac_irr_gb,                                               &
   frac_irr_old_gb, sthu_irr_gb,                                              &
   dvi_cpft

 USE jules_soil_mod, ONLY : sm_levels

 IMPLICIT NONE

  INTEGER ::                                                                  &
   l,n

  CHARACTER (len=50) ::                                                       &
   frac_irr_str

! Local arrays
  REAL ::                                                                     &
   sthu_nir(land_pts,sm_levels)
                        ! soil moisture content in non-irrigated fraction as
                        ! a fraction of saturation

 !-------------------------------------------------------------------------------
!   When using JULES-crop model, adjust irrigation fraction
!   to those tiles that have suitable development index
!-------------------------------------------------------------------------------
!   When irr_crop=2, JULES-crop is used to determine irrigation period.
!   However, within a grid box the development index can be different for each
!   crop tile while irrigation is applied to the (single) irrigated fraction of
!   that grid box.
!   Here frac_irr_surft is set to zero in those tiles where the dvi_cpft is below
!   the threshold and irrigation is only applied to those tiles with suitable
!   dvi_cpft, i.e. the actual irrigation fraction grows or shrinks according to
!   which crop tiles have the right dvi_cpft. Note that frac_irr_gb can also be
!   allocated to non-crop tiles, these will (currently) only be irrigated when
!   any crop has a dvi_cpft above the threshold - see subroutine irrig_dmd.

  IF ( irr_crop == 2 ) THEN
    DO l=1,land_pts
      DO n=1,ncpft
        IF ( dvi_cpft(l,n) <= -1 ) THEN
!         no suitable dvi_cpft - set frac_irr_gb for this crop tile to zero
!             and subtract from overall irrigation fraction
          frac_irr_gb(l) = frac_irr_gb(l)                                     &
                      - frac_surft(l,n + nnpft) * frac_irr_surft(l,n + nnpft)
          frac_irr_gb(l) = MAX( frac_irr_gb(l), 0.0 )
          frac_irr_surft(l,n + nnpft) = 0.0
        END IF
      END DO
    END DO
  END IF

!-------------------------------------------------------------------------------
!   Update soil moisture content in irrigated fraction (sthu_irr_gb)
!-------------------------------------------------------------------------------
!   If frac_irr expands, the 'added fraction' does not have the moisture content
!   of the irrigated fraction, but of the non-irrigated fraction (which is not
!   a prognostic). Total gridbox moisture content should remain the same
!   Conversely, when frac_irr shrinks, the moisture content in the non-irrig
!   fraction should become higher, but in the (remaining) irrigated fraction
!   (and total gridbox) it remains the same, so no need to cater for that here

  DO l=1,land_pts
    IF ( frac_irr_gb(l) > frac_irr_old_gb(l) ) THEN
      DO n=1,sm_levels
        sthu_nir(l,n) = (sthu_gb(l,n) - frac_irr_old_gb(l) * sthu_irr_gb(l,n))&
                      / (1.0 - frac_irr_old_gb(l))
        sthu_irr_gb(l,n) = ( frac_irr_old_gb(l) * sthu_irr_gb(l,n) +          &
                          (frac_irr_gb(l)-frac_irr_old_gb(l)) * sthu_nir(l,n))&
                      / frac_irr_gb(l)
      END DO
    END IF
  END DO

!------------------------------------------------------------------------------
! Write to stdout irrigation fraction for each gridbox.
!-------------------------------------------------------------------------------
  WRITE(frac_irr_str, "('Range of frac irrigation =',8f5.2,' to ',8f5.2)")    &
                                            MINVAL(frac_irr_gb),              &
                                            MAXVAL(frac_irr_gb)

  IF (current_time%day == 1) THEN
    CALL log_info("update_ancil_irrig",                                       &
                  '### NB The ranges below include any ice points. ###')
    CALL log_info("update_ancil_irrig",                                       &
                  datetime_to_string(current_time) // frac_irr_str)
  END IF


END SUBROUTINE update_irrig_variables
#endif
