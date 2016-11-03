#if !defined(UM_JULES)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  *** JULES version of science_fixes_mod ***

! Description:
!   Module holding logical that enables full implementation of
!   emissivity for sea and sea-ice with the UM.
!
! Code Owner: See Unified Model Code Owner's HTML page

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3

MODULE science_fixes_mod

IMPLICIT NONE

LOGICAL :: l_emis_ssi_full = .TRUE.

! Set this logical to .TRUE. to correct the updating of the surface
! temperature in the implicit solver.
LOGICAL :: l_dtcanfix = .TRUE.

! Fixes the calculation of surface exchange in coastal grid-boxes when
! coastal tiling is switched on. Should have no effect in standalone
! JULES, but default to TRUE anyway
LOGICAL :: l_fix_ctile_orog = .TRUE.

! Fixes how ustar is included in the exchange coefficient for dust deposition. 
! Has no effect in standalone JULES, but default to TRUE anyway
LOGICAL :: l_fix_ustar_dust = .TRUE.

END MODULE science_fixes_mod
#endif
