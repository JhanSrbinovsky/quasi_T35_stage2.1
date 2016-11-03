! *****************************COPYRIGHT****************************************
! (c) Crown copyright, Met Office 2012. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use
! and distribution under the JULES collaboration agreement, subject
! to the terms and conditions set out therein.
!
! [Met Office Ref SC0237]
! *****************************COPYRIGHT****************************************

SUBROUTINE calc_mon_no(year, day, L_CAL360, imonth)

!  < Module imports >
! None

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Returns the month number. Needed for the UM as month number doesn't exist
!
! Method:
!
! Current Code Owner: Richard Gilham (Met Office)
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

! Subroutine arguments
  INTEGER, INTENT(IN)  :: year,day
  LOGICAL, INTENT(IN)  :: L_CAL360

  INTEGER, INTENT(OUT) :: imonth

! Local constants
! None

! Local variables
  INTEGER :: end_jan, end_feb, end_mar, end_apr, end_may, end_jun,            &
               end_jul, end_aug, end_sep, end_oct, end_nov, end_dec

  LOGICAL :: is_leap_year

!End of header
!-----------------------------------------------------------------------------

    is_leap_year = ( MOD(year, 4) == 0 .AND. MOD(year, 100) /= 0 ) .OR.       &
                   ( MOD(year, 400) == 0 )

    IF (L_CAL360) THEN
      end_jan = 30
      end_feb = end_jan + 30
      end_mar = end_feb + 30
      end_apr = end_mar + 30
      end_may = end_apr + 30
      end_jun = end_may + 30
      end_jul = end_jun + 30
      end_aug = end_jul + 30
      end_sep = end_aug + 30
      end_oct = end_sep + 30
      end_nov = end_oct + 30
      end_dec = end_nov + 30
    ELSE
      end_jan = 31
      IF (is_leap_year) THEN
        end_feb = end_jan + 29
      ELSE
        end_feb = end_jan + 28
      END IF
      end_mar = end_feb + 31
      end_apr = end_mar + 30
      end_may = end_apr + 31
      end_jun = end_may + 30
      end_jul = end_jun + 31
      end_aug = end_jul + 31
      end_sep = end_aug + 30
      end_oct = end_sep + 31
      end_nov = end_oct + 30
      end_dec = end_nov + 31
    END IF

    IF      (day >= 1      .AND. day <= end_jan) THEN
      imonth = 1
    ELSE IF (day > end_jan .AND. day <= end_feb) THEN
      imonth = 2
    ELSE IF (day > end_feb .AND. day <= end_mar) THEN
      imonth = 3
    ELSE IF (day > end_mar .AND. day <= end_apr) THEN
      imonth = 4
    ELSE IF (day > end_apr .AND. day <= end_may) THEN
      imonth = 5
    ELSE IF (day > end_may .AND. day <= end_jun) THEN
      imonth = 6
    ELSE IF (day > end_jun .AND. day <= end_jul) THEN
      imonth = 7
    ELSE IF (day > end_jul .AND. day <= end_aug) THEN
      imonth = 8
    ELSE IF (day > end_aug .AND. day <= end_sep) THEN
      imonth = 9
    ELSE IF (day > end_sep .AND. day <= end_oct) THEN
      imonth = 10
    ELSE IF (day > end_oct .AND. day <= end_nov) THEN
      imonth = 11
    ELSE IF (day > end_nov .AND. day <= end_dec) THEN
      imonth = 12
    END IF

    RETURN
END SUBROUTINE calc_mon_no
