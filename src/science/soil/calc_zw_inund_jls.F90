#if !defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!   Subroutine CALC_ZW_INUND-------------------------------------------

!   Purpose: To calculate the mean water table depth from the soil
!            moisture deficit as described in Koster et al., 2000.,
!            using the Newton-Raphson method.

! Documentation: UNIFIED MODEL DOCUMENTATION PAPER NO 25

SUBROUTINE calc_zw_inund(npnts,nshyd,soil_pts,soil_index,zdepth               &
  ,bexp,sathh,smclsat,smcl,sthu,sthzw                                         &
  ,zw,zw_inund,wutot)

USE c_densty
USE jules_hydrology_mod, ONLY : zw_max

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine arguments
!   Scalar arguments with intent(IN) :
INTEGER                                                                       &
 npnts                                                                        &
                      ! IN Number of gridpoints.
,nshyd                                                                        &
                      ! IN Number of soil moisture levels.
,soil_pts             ! IN Number of soil points.

!   Array arguments with intent(IN) :
INTEGER                                                                       &
 soil_index(npnts)    ! IN Array of soil points.

REAL                                                                          &
 bexp(npnts)                                                                  &
!                     ! IN Clapp-Hornberger exponent.
,sathh(npnts)                                                                 &
!                     ! IN Saturated soil water pressure (m).
,smcl(npnts,nshyd)                                                            &
!                     ! IN Total soil moisture contents
!                           !      of each layer (kg/m2).
,smclsat(npnts,nshyd)                                                         &
                      ! IN Soil moisture contents of each
!                           !      layer at saturation (kg/m2).
,sthu(npnts,nshyd)                                                            &
!                     ! IN Unfrozen soil moisture content of
!                           !    each layer as a fraction of
!                           !    saturation.
,sthzw(npnts)                                                                 &
!                     ! IN Fraction soil moisture content
!                           !       in deep layer.
,zdepth(0:nshyd)
!                     ! IN Soil layer depth at lower boundary (m).


!   Array arguments with intent(INOUT) :
REAL                                                                          &
 zw(npnts)            ! INOUT Water table depth (m).
REAL                                                                          &
 zw_inund(npnts)                                                              &
!                     ! OUT Adjusted water table depth (m).
,wutot(npnts)         ! OUT UNFROZEN to TOTAL fraction at ZW.

REAL                                                                          &
 unfroz_frac(npnts,nshyd+1)                                                   &
!                     ! WORK fractional space available for
!                     !      unfrozen soil moisture in soil layer.
,zwest_l(nshyd)                                                               &
!                     ! WORK water estimate with Newton-Raphson
!                     !      iteration.
,psi(nshyd)                                                                   &
!                     ! WORK soil pressure in soil layer.
,sthf(npnts,nshyd)                                                            &
!                     ! WORK Frozen soil moisture content of
!                     !    each layer as a fraction of
!                     !    saturation.
,sthuk                                                                        &
!                     ! WORK fractional unfrozen soil moisture
!                     !      interpolated to water table.
,sthfk                                                                        &
!                     ! WORK fractional frozen soil moisture
!                     !      interpolated to water table.
,sthzwu                                                                       &
!                     ! WORK Estimated fractional unfrozen soil
!                     !      moisture content in water table layer.
!                     !      (assuming sthf in this layer=sthf(nshyd)
,avgzd(nshyd+1)
!                     ! WORK Soil layer depth in centre (m).


REAL                                                                          &
 zw_old                                                                       &
!                     ! WORK zw from last timestep
,psisat                                                                       &
!                     ! WORK Saturated soil water pressure (m)
!                     !      (negative).
,fac1                                                                         &
!                     ! WORK Soil moisture related variable.
,fac2                                                                         &
!                     ! WORK Soil moisture related variable.
,facb                                                                         &
!                     ! WORK fac^bexp.
,unfroz_frac_avg                                                              &
!                     ! WORK Weighted average of unfroz_frac.
,zwest_mean                                                                   &
!                     ! WORK Weighted average of zwest_l.
,zwdef
!                     ! WORK Water table depth related parameter.

INTEGER                                                                       &
 i,j,n                                                                        &
            ! Loop counters
,it                                                                           &
            ! Loop counters
,niter                                                                        &
            ! Number of iterations
,maxinter
PARAMETER(niter=3)
!PARAMETER(niter=5)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CALC_ZW_INUND'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

!----------------------------------------------------------------------
! Calculate depth of centre point of each layer:
!----------------------------------------------------------------------
DO n=1,nshyd
  avgzd(n)=0.5*(zdepth(n-1)+zdepth(n))
END DO
avgzd(nshyd+1)=0.5*(zw_max-zdepth(nshyd))+zdepth(nshyd)

!----------------------------------------------------------------------
! Calculate the available space for unfrozen soil moisture in each soil layer:
!----------------------------------------------------------------------
DO j=1,soil_pts
i=soil_index(j)

  DO n=1,nshyd
    IF(smclsat(i,n) > epsilon(1.0))THEN
      sthf(i,n)=smcl(i,n)/smclsat(i,n)-sthu(i,n)
    ELSE
      sthf(i,n)=0.0
    END IF

    unfroz_frac(i,n)=1.0-sthf(i,n)
    unfroz_frac(i,n)=max(unfroz_frac(i,n),0.0)
    unfroz_frac(i,n)=min(unfroz_frac(i,n),1.0)
  END DO

  unfroz_frac(i,nshyd+1)=1.0-sthf(i,nshyd)
  unfroz_frac(i,nshyd+1)=max(unfroz_frac(i,nshyd+1),0.0)
  unfroz_frac(i,nshyd+1)=min(unfroz_frac(i,nshyd+1),1.0)
END DO

!----------------------------------------------------------------------
! Add code to allow for the fact that the water table isnt well estimated
! using the method in calc_zw_jls.F90 when one layer is significantly
! frozen (see documentation).
! Instead use the same fundamental equations above but for each
! layer containing or above the provisional water table
! and weight accordingly:
!----------------------------------------------------------------------
DO j=1,soil_pts
  i=soil_index(j)

  zw_inund(i)=zw(i)
  wutot(i)=1.0
  unfroz_frac_avg=0.0
  DO n=1,nshyd
    unfroz_frac_avg=unfroz_frac_avg+unfroz_frac(i,n)*(zdepth(n)-zdepth(n-1))
  END DO
  unfroz_frac_avg=unfroz_frac_avg+unfroz_frac(i,nshyd+1)*(zw_max-zdepth(nshyd))
  unfroz_frac_avg=unfroz_frac_avg/zw_max

  IF(1.0-unfroz_frac_avg > EPSILON(1.0).and.zw(i) < zw_max)THEN
    zw_old=zw(i)
    DO n=1,nshyd
      zwest_l(n)=0.0
    END DO
    zwest_mean=0.0
    zwdef=zw_old
    IF(zw_old > zdepth(nshyd))zwdef=zdepth(nshyd)

!-----------------------------------------------------------------------
! interpolate thetaf to zw and use this in psisat
!-----------------------------------------------------------------------
    IF(zw(i) < avgzd(1))sthfk=sthf(i,1)
    DO n=1,nshyd-1
      IF(zw(i) >= avgzd(n).and.zw(i) < avgzd(n+1))THEN
        sthfk=(sthf(i,n)*(avgzd(n+1)-zw(i))                                   &
        +sthf(i,n+1)*(zw(i)-avgzd(n))  )                                      &
        /(avgzd(n+1)-avgzd(n))
      END IF
    END DO
    IF(zw(i) >= avgzd(nshyd))sthfk=sthf(i,nshyd)

!-----------------------------------------------------------------------
! Now include frozen soil into calc of PSI_sat:
!-----------------------------------------------------------------------
    fac1=(1.0-sthfk)
    fac1=min(fac1,1.0)
    fac1=max(fac1,0.0)
    facb=fac1**bexp(i)
    facb=max(facb,epsilon(1.0))
    psisat=-sathh(i)/facb

    DO n=1,nshyd

      fac2=sthu(i,n)
      fac2=min(fac2,1.0)
      fac2=max(fac2,0.0)
      facb=fac2**bexp(i)
      facb=max(facb,epsilon(1.0))
      psi(n)=-sathh(i)/facb

      zwest_l(n)=psisat-psi(n)+0.5*(zdepth(n)+zdepth(n-1))
      zwest_l(n)=max(zwest_l(n),0.0)
      zwest_l(n)=min(zwest_l(n),zw_max)

      IF(zwdef <  zdepth(n).AND.zwdef >  zdepth(n-1))                         &
        zwest_mean=zwest_mean+zwest_l(n)*(zwdef-zdepth(n-1))

      IF(zwdef >= zdepth(n))                                                  &
       zwest_mean=zwest_mean+zwest_l(n)*(zdepth(n)-zdepth(n-1))

     END DO
     zwest_mean=zwest_mean/zwdef

     zwest_mean=max(zwest_mean,0.0)
     zwest_mean=min(zwest_mean,zw_max)
     IF(zwest_mean > zw_old)zw_inund(i)=unfroz_frac_avg*zw_old                &
                               +(1.0-unfroz_frac_avg)*zwest_mean

  END IF
END DO


!-----------------------------------------------------------------------
! wutot needs to be based on corrected zw due to partially
! frozen soil moisture fraction:
! interpolate thetas to zw_inund and use this in wutot
!-----------------------------------------------------------------------
DO j=1,soil_pts
  i=soil_index(j)

     IF(zw_inund(i) < avgzd(1))sthfk=sthf(i,1)
     IF(zw_inund(i) < avgzd(1))sthuk=sthu(i,1)

     DO n=1,nshyd-1

        IF(zw_inund(i) >= avgzd(n).and.zw_inund(i) < avgzd(n+1))THEN
          sthfk=(sthf(i,1)*(avgzd(n+1)-zw_inund(i))                           &
          +sthf(i,n+1)*(zw_inund(i)-avgzd(n))  )                              &
          /(avgzd(n+1)-avgzd(n))
          sthuk=(  sthu(i,n)*(avgzd(n+1)-zw_inund(i))                         &
          +sthu(i,n+1)*(zw_inund(i)-avgzd(n))  )                              &
          /(avgzd(n+1)-avgzd(n))
        END IF

     END DO

     sthzwu=sthzw(i)-sthf(i,nshyd)
     sthzwu=max(sthzwu,0.0)

     IF(zw_inund(i) >= avgzd(nshyd))sthfk=sthf(i,nshyd)
     IF(zw_inund(i) >= avgzd(nshyd).and.zw_inund(i) < avgzd(nshyd+1))THEN
       sthuk=(  sthu(i,nshyd)*(avgzd(nshyd+1)-zw_inund(i))                    &
       +sthzwu*(zw_inund(i)-avgzd(nshyd))  )                                  &
       /(avgzd(nshyd+1)-avgzd(nshyd))
     END IF
     IF(zw_inund(i) >= avgzd(nshyd+1))sthuk=sthzwu

     IF(sthuk+sthfk > epsilon(1.0))                                           &
        wutot(i)=sthuk/(sthuk+sthfk)

   END DO

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE calc_zw_inund

#endif
