!  ------------------COMDECK C_PI---------------------------------------
!  
!   4.0 19/09/95  New value for PI. Old value incorrect
!                 from 12th decimal place. D. Robinson
!   5.1 7/03/00   Fixed/Free format P.Selwood
!  

      ! Pi
      Real, Parameter :: Pi                 = 3.14159265358979323846

      ! Conversion factor degrees to radians
      Real, Parameter :: Pi_Over_180        = Pi/180.0

      ! Conversion factor radians to degrees
      Real, Parameter :: Recip_Pi_Over_180  = 180.0/Pi

! ----------------------------------------------------------------------
