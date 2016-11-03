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

SUBROUTINE assign_irrig_fraction ()

  USE logging_mod, ONLY : log_fatal

  USE string_utils_mod, ONLY : to_string

  USE ancil_info, ONLY :                                             &
  land_pts, nsurft, frac_surft

  USE jules_surface_types_mod, ONLY : ntype

  USE crop_vars_mod, ONLY :                                          &
!  imported logicals with intent(in)
   frac_irr_all_tiles,                                               &
!  imported arrays with intent(in)
   frac_irr_all,irrtiles,                                            &
   frac_irr_gb, frac_irr_surft,                                          &
   frac_irr_old_gb,                                                     &
   nirrtile

  !USE jules_vegetation_mod, ONLY : irr_crop

  IMPLICIT NONE

  INTEGER ::                                                         &
   m,l,k,n
                         ! Loop indices

  INTEGER ::   alltiles(nsurft)

! Local arrays
  REAL ::                                                            &
   frac_irr_rest(land_pts)                                           !&
                        ! Remaining irrigation fraction to assign to tiles
   !frac_irr_old_gb(land_pts)
                        ! previous irrigated fraction in grid box

!-------------------------------------------------------------------------------
!   Update irrigation ancillary
!-------------------------------------------------------------------------------

! Loop over land points and assign irrigation fraction
! note frac_irr_all has shape (land_points, year) or (land_points, 1)
  DO l=1,land_pts
    frac_irr_old_gb(l) = frac_irr_gb(l)
    frac_irr_gb(l)=frac_irr_all(l,1)

!     Set a minimum threshold for irrigation:
    IF( frac_irr_gb(l) < 1.0E-3 ) frac_irr_gb(l)= 0.0
  END DO

!-------------------------------------------------------------------------------
!   Assign irrigation fraction to each tile
!   hadrd - this was originally done in physiol.f90
!-------------------------------------------------------------------------------
  IF ( .NOT. frac_irr_all_tiles ) THEN
! assign irrigated fraction to pre-defined tiles ONLY
! in decreasing order of importance
    frac_irr_surft(:,:) = 0.0
    frac_irr_rest(:) = frac_irr_gb(:)

! create index of all tiles (veg+non-veg)
    DO m=1,nsurft
      alltiles(m) = m
    END DO

    DO l=1,land_pts
      DO k=1,nirrtile ! loop over irrigated pfts

! for each tile, check if the index corresponds with index in irrtiles
        DO n=1,ntype
          IF ( alltiles(n) == irrtiles(k) ) THEN
! assign (remaining) irrigated fraction to this tile
            frac_irr_surft(l,n) = min( frac_irr_rest(l), frac_surft(l,n) )
            frac_irr_rest(l) = frac_irr_rest(l) - frac_irr_surft(l,n)
! check for negative remaining frac_irr_gb
            IF (frac_irr_rest(l) < 0.0) THEN
              CALL log_fatal("update_ancil_irrig",                            &
                "Error in assigning irrigated fraction at point " //          &
                to_string(l) // " - " //                                      &
                "irrigated fraction: " // to_string(frac_irr_gb(l)) // ", " //   &
                "tile fraction: " // to_string(frac_surft(l,n)) // ", " //          &
                "irrigated fraction assigned:"// to_string(frac_irr_surft(l,n)))
            END IF
          END IF
        END DO ! ntype
      END DO ! irrtilenames
    END DO ! land_pts

    ! to prevent negative fractions (this can happen due to numerical inaccuracy)
    frac_irr_surft(:,:) = max(frac_irr_surft(:,:),0.0)

! in sf_evap, frac_irr_surft is used as a multiplier to frac
    WHERE ( frac_surft > 0.0 ) frac_irr_surft = frac_irr_surft / frac_surft

  ELSE ! if frac_irr_all_tiles = TRUE
! assign irrigation fraction to all tiles
! to reproduce the original results
    frac_irr_surft(:,:) = 0.0
    DO n=1,nsurft
      WHERE ( frac_irr_gb > epsilon(1.0) ) frac_irr_surft(:,n) = frac_irr_gb
    END DO

  END IF ! .not. frac_irr_all_tiles

END SUBROUTINE assign_irrig_fraction
#endif
