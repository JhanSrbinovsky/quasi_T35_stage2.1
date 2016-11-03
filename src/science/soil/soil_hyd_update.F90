! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE SOIL_HYD_UPDATE-----------------------------------------------


SUBROUTINE soil_hyd_update(npnts,nshyd,soil_pts,soil_index,                   &
                 dz,v_sat,zdepth,smclzw,sthzw,smclsat,smclsatzw)

USE jules_hydrology_mod, ONLY :                                               &
!      imported scalar parameters
   zw_max

USE c_densty, ONLY :                                                          &
!      imported scalar parameters
   rho_water  !  density of pure water (kg/m3)

IMPLICIT NONE

! Subroutine arguments
!   Scalar arguments with intent(IN) :
INTEGER, INTENT(IN) ::                                                        &
npnts                                                                         &
                      ! IN Number of gridpoints.
,nshyd                                                                        &
                      ! IN Number of soil moisture levels.
,soil_pts
                      ! IN Number of soil points.
REAL, INTENT(IN) ::                                                           &
dz(nshyd)
                      ! IN Thicknesses of the soil layers (m).

!   Array arguments with intent(IN) :
INTEGER, INTENT(IN) ::                                                        &
soil_index(npnts)
                      ! IN Array of soil points.

REAL, INTENT(IN) ::                                                           &
v_sat(npnts,nshyd)                                                            &
                      ! IN Volumetric soil moisture
!                           !    concentration at saturation
!                           !    (m3 H2O/m3 soil).
,zdepth(0:nshyd)
                      ! IN Soil layer depth at lower boundary (m).
REAL, INTENT(IN) ::                                                           &
sthzw(npnts)         ! IN soil moist fraction in deep layer.

!   Array arguments with intent(OUT) :
REAL, INTENT(OUT) ::                                                          &
smclsat(npnts,nshyd)
                      ! OUT The saturation moisture content of
                           !     each layer (kg/m2).
REAL, INTENT(OUT) ::                                                          &
smclsatzw(npnts)                                                              &
                      ! OUT moisture content in deep layer
                           !      at saturation.
,smclzw(npnts)
                      ! OUT moisture content in deep layer.

! Local scalars:
INTEGER                                                                       &
 i,j,n                ! WORK Loop counters.


!-----------------------------------------------------------------------------

    DO n=1,nshyd
      DO j=1,soil_pts
        i=soil_index(j)
        smclsat(i,n)=rho_water*dz(n)*v_sat(i,n)
      END DO
    END DO

    DO j=1,soil_pts
      i=soil_index(j)
      smclsatzw(i)=rho_water*v_sat(i,nshyd)                                   &
          *(zw_max-zdepth(nshyd))
      smclzw(i)=sthzw(i)*smclsatzw(i)
    END DO

RETURN
END SUBROUTINE soil_hyd_update
