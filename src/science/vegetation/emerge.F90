! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

SUBROUTINE emerge(n, t_surft, dvi, rootc, harvc, stemc, leafc)

  USE cropparm, ONLY : t_bse, tt_emr, initial_carbon

  USE timestep_mod, ONLY : timestep

  USE crop_utils_mod, ONLY : carbon_fraction_from_dvi,                        &
                             croprootc_min, cropharvc_min,                    &
                             calc_cropstemc_min, calc_cropleafc_min

  USE jules_print_mgr, ONLY : jules_print

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Increases DVI based on thermal time, until crop emerges
!   DVI increases from -1 to 0
!
! Method:
!   Crop should already have been sown, but not yet emerged when this
!   subroutine is called.
!   Subroutine calculates thermal time (effective temperature)
!   from tile temperature and uses this to increase crop
!   development index (DVI). Once dvi >= 0.0, crop has emergedm and is started
!   off with a bit of carbon.
!
!   See JULES-crop technical documentation (Tom Osborne, Josh Hooker).
!
! Current Code Owner: Tom Osborne
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

  INTEGER, INTENT(IN) :: n      ! crop tile number

  REAL, INTENT(IN)    :: t_surft ! temperature (T1P5M) on tile

  REAL, INTENT(INOUT) :: dvi    ! crop development index

  REAL, INTENT(INOUT) :: rootc  ! root carbon pool
  REAL, INTENT(INOUT) :: harvc  ! harvest carbon pool
  REAL, INTENT(INOUT) :: leafc  ! leaf carbon pool
  REAL, INTENT(INOUT) :: stemc  ! stem carbon pool

! Local variables

  REAL :: teff    ! effective temperature

  REAL :: f_root  ! fraction to roots
  REAL :: f_stem  ! fraction to stem
  REAL :: f_leaf  ! fraction to leaf
  REAL :: f_harv  ! fraction to harvested parts


!-----------------------------------------------------------------------------


  teff = t_surft - t_bse(n)

  teff = MAX(teff, 0.0) ! Presume that development can't go back

  teff = teff * ( timestep / 86400.0 ) ! times fraction of full day

  dvi = dvi + ( teff / tt_emr(n) )

  IF (dvi >= 0.0) THEN ! Crop has emerged

    CALL carbon_fraction_from_dvi(n, dvi, f_root, f_stem, f_leaf, f_harv)

    !----------------------------
    !  Start crop with a small amount of carbon
    !----------------------------

    rootc = croprootc_min         + initial_carbon(n) * f_root
    leafc = calc_cropleafc_min(n) + initial_carbon(n) * f_leaf
    stemc = calc_cropstemc_min(n) + initial_carbon(n) * f_stem
    harvc = cropharvc_min         + initial_carbon(n) * f_harv

  END IF

END SUBROUTINE emerge
