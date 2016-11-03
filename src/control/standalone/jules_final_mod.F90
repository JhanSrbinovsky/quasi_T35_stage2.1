#if !defined(UM_JULES)
!******************************COPYRIGHT**************************************
! (c) Centre for Ecology and Hydrology 2016. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms
! and conditions set out therein.
!
! [Met Office Ref SC0237] 
!******************************COPYRIGHT**************************************

MODULE jules_final_mod

IMPLICIT NONE

CONTAINS

!#############################################################################

SUBROUTINE jules_final

  USE jules_vegetation_mod, ONLY :                                            &
!  imported scalars that are not altered
    can_rad_mod

  USE logging_mod, ONLY :                                                     &
!  imported procedures
    log_warn


  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Writes final warning messages at the end of a JULES run.
!
! Current Code Owner: Douglas Clark (CEH)
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

! Local scalar parameters.
  CHARACTER(LEN=*), PARAMETER ::                                              &
    proc_name = 'jules_final'   !  Name of this procedure.

!-----------------------------------------------------------------------------
! Warn if a deprecated value of can_rad_mod is used.
! Using multiple calls to make this more prominent!
!-----------------------------------------------------------------------------
  SELECT CASE ( can_rad_mod )
    CASE ( 2, 3 )
      CALL log_warn(proc_name,'#############################################')
      CALL log_warn(proc_name,'#############################################')
      CALL log_warn(proc_name, 'can_rad_mod=2 and 3 are deprecated and ' //   &
                    'will be removed in future. Use 1 or 6 instead.')
      CALL log_warn(proc_name,'#############################################')
      CALL log_warn(proc_name,'#############################################')
  END SELECT

END SUBROUTINE jules_final

!#############################################################################

END MODULE jules_final_mod
#endif
