! SUBROUTINE BEDROCK -------------------------------------------
! Description:
!     Updates deep soil temperatures in deep layers


! Subroutine Interface:
SUBROUTINE bedrock ( npnts,soil_pts,soil_index,timestep,                      &
    tsoil,hcsoil,dzsoil,hflux_in                                              &
)

USE prognostics, ONLY : tsoil_deep_gb

USE jules_soil_mod, ONLY : ns_deep, hcapdeep, hcondeep, dzdeep

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

INTEGER :: npnts                                                              &
             ! IN number of land points
          ,soil_pts                                                           &
             ! IN number of soil points
          ,soil_index(npnts)                                                  &
             ! IN index of soil points
          ,i,j,n  ! WORK loop counter

REAL :: tsoil(npnts)                                                          &
             ! IN soil temp at base of column (Celsius)
       ,hcsoil(npnts)                                                         &
             ! IN heat conductivity of base soil layer
       ,dzsoil                                                                &
             ! IN thickness of base soil layer
       ,timestep
             ! IN model timestep

REAL :: hflux_in(npnts)
             ! OUT heat flux from base of soil column into bedrock layers

REAL :: hctop(npnts)                                                          &
             ! WORK interpolated heat conductivity where bedrock joins soil
       ,dztop                                                                 &
             ! WORK interpolated layer thickness for heat transfer
!            !      between base of soil and top of bedrock.
       ,dtsd(npnts,ns_deep)                                                   &
             ! WORK temperature increment
       ,tsoil_k(npnts)
             ! WORK temperature of base soil layer in Kelvin

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='BEDROCK'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)


! Calculate dztop (doesn't change spatially)
dztop = 0.5*(dzdeep + dzsoil)

! Convert base soil temp to Kelvin
DO j=1,soil_pts
  i = soil_index(j)
  tsoil_k(i) = tsoil(i) + 273.15
END DO

DO j=1,soil_pts
  i = soil_index(j)

  !-------------------------------------------------------------------------
  ! Calculate the thermal conductivity at the top of the bedrock.
  !-------------------------------------------------------------------------
  hctop(i) = (dzdeep*hcsoil(i) + dzsoil*hcondeep)/(dzdeep+dzsoil)

  ! Calculate the heat flux from the soil above.
  hflux_in(i) = hctop(i)*(tsoil_k(i)-tsoil_deep_gb(i,1))/dztop

  !-------------------------------------------------------------------------
  ! Calculate the temperature increments to the bedrock layers.
  !-------------------------------------------------------------------------
  IF (ns_deep > 1) THEN
    !! bottom:
    dtsd(i,ns_deep) = hcondeep*timestep*(tsoil_deep_gb(i,ns_deep-1)-          &
        tsoil_deep_gb(i,ns_deep))/(hcapdeep*dzdeep**2)
    !! top:
    dtsd(i,1) = timestep*( hcondeep*(tsoil_deep_gb(i,2)-                      &
        tsoil_deep_gb(i,1))/dzdeep + hflux_in(i) )/(hcapdeep*dzdeep)
  ELSE
    dtsd(i,1) = timestep*hflux_in(i)/(hcapdeep*dzdeep)
  END IF

  IF (ns_deep > 2) THEN
    DO n=2,ns_deep-1
      dtsd(i,n) = timestep*hcondeep*(tsoil_deep_gb(i,n+1)+                    &
        tsoil_deep_gb(i,n-1)-2*tsoil_deep_gb(i,n))/(hcapdeep*dzdeep**2)
    END DO
  END IF
  !-------------------------------------------------------------------------
  ! Update the layer temperatures
  !-------------------------------------------------------------------------
  DO n=1,ns_deep
    tsoil_deep_gb(i,n) = MAX(tsoil_deep_gb(i,n) + dtsd(i,n),0.0)
    tsoil_deep_gb(i,n) = MIN(tsoil_deep_gb(i,n),10000.0)
  END DO

END DO


IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE bedrock
