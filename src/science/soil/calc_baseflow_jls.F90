! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE CALC_BASEFLOW------------------------------------------

! Description:
!     Calculates subsurface runoff (aka baseflow).

MODULE calc_baseflow_jls_mod
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CALC_BASEFLOW_JLS_MOD'

CONTAINS
! Subroutine Interface:
SUBROUTINE calc_baseflow(                                                     &
 soil_pts,soil_index,npnts,nshyd                                              &
,zdepth,ksz                                                                   &
,b,fexp,ti_mean,zw,sthf                                                       &
,top_crit,qbase,qbase_l                                                       &
 )

USE jules_hydrology_mod, ONLY : ti_max, zw_max
USE jules_soil_mod, ONLY : dzsoil

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE jules_print_mgr, ONLY :                                                   &
    jules_message,                                                            &
    jules_print
IMPLICIT NONE

INTEGER                                                                       &
 npnts                                                                        &
                     ! IN Number of gridpoints.
,nshyd                                                                        &
                     ! IN Number of soil moisture levels.
,soil_pts            ! IN Number of soil points.

!   Array arguments with intent(IN) :
INTEGER                                                                       &
 soil_index(npnts)   ! IN Array of soil points.

REAL                                                                          &
 b(npnts)                                                                     &
                     ! IN Clapp-Hornberger exponent.
,fexp(npnts)                                                                  &
                     ! IN Decay factor in Sat. Conductivity
!                          !    in water table layer.
,ti_mean(npnts)                                                               &
                     ! IN Mean topographic index.
,zw(npnts)                                                                    &
                      ! IN   Water table depth (m).
,ksz(npnts,0:nshyd)                                                           &
                      ! IN Saturated hydraulic condy
!                           !      for each layer (kg/m2/s).
,ksfz(npnts,nshyd+1)                                                          &
                      ! WORK Function of sat. hydraulic
!                           !      conductivity, frozen soil
!                           !      and mean topographic index.
,sthf(npnts,nshyd)
                     ! INOUT Frozen soil moisture content of
!                          !       each layer as a fraction of
!                          !       saturation.

REAL                                                                          &
 qbase(npnts)                                                                 &
                      ! OUT Total base flow (kg/m2/s).
,qbase_l(npnts,nshyd+1)                                                       &
!                           ! OUT Base flow from each layer (kg/m2/s).
,top_crit(npnts)
                      ! OUT Critical topographic index required
!                           !     to calc surface saturation frac.

REAL                                                                          &
 qbase_max(npnts)                                                             &
                      ! WORK Max possible base flow (kg/m2/s).
,qbase_max_l(npnts,nshyd+1)                                                   &
!                           ! WORK Max possible base flow
!                           !     from each level (kg/m2/s).
,zdepth(0:nshyd)                                                              &
                      ! WORK Lower soil layer boundary depth (m).
,qbase_min(npnts)     ! WORK Residual base flow at zw_max

! Local scalars:
INTEGER                                                                       &
 i,j                                                                          &
                      ! WORK Loop counters.
,n                                                                            &
                      ! WORK Tile loop counter.
,errorstatus          ! WORK Error status for ereport.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CALC_BASEFLOW'

!-----------------------------------------------------------------------
! Initialise TOP_CRIT to maximum value.
! Initialise baseflow components to zero.
!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO j=1,soil_pts
  i=soil_index(j)
  top_crit(i)=ti_max
  qbase_max(i)=0.0
  DO n=1,nshyd
    qbase_max_l(i,n)=0.0
  END DO
END DO

DO j=1,soil_pts
  i=soil_index(j)

!-----------------------------------------------------------------------
! Calculate layer dependent variable with is dependent on
! effective saturated conductivity:
!-----------------------------------------------------------------------
  DO n=1,nshyd
    IF(sthf(i,n) <  1.0)THEN
      ksfz(i,n)=0.5*(ksz(i,n-1)+ksz(i,n))                                     &
               *(1.0-sthf(i,n))**(2.0*b(i)+3.0)*EXP(-ti_mean(i))
    ELSE
      ksfz(i,n)=0.0
    END IF
  END DO
  IF(sthf(i,nshyd) <  1.0)THEN
    ksfz(i,nshyd+1)=ksz(i,nshyd)/fexp(i)                                      &
                   *(1.0-sthf(i,nshyd))**(2.0*b(i)+3.0)                       &
                   *EXP(-ti_mean(i))
  ELSE
    ksfz(i,nshyd+1)=0.0
  END IF


  IF(EXP(-fexp(i)*(zw_max-zdepth(nshyd))) >  0.05)THEN
    WRITE(jules_message,*)'CB_L1A: maximum water table depth is too low!!!'
    CALL jules_print('calc_baseflow_jls',jules_message)
    WRITE(jules_message,*)'at soil point',i,fexp(i)
    CALL jules_print('calc_baseflow_jls',jules_message)
    WRITE(jules_message,*)EXP(-fexp(i)*(zw_max-zdepth(nshyd)))
    CALL jules_print('calc_baseflow_jls',jules_message)
    WRITE(jules_message,*)'try zw_max>',-LOG(0.05)/fexp(i)+zdepth(nshyd)
    CALL jules_print('calc_baseflow_jls',jules_message)
    errorstatus=1000
!            CALL EREPORT('CALC_BASEFLOW', ERRORSTATUS,                  &
!     &        'ERROR ZW_MAX is TOO SMALL')
  END IF

  !-----------------------------------------------------------------------
  ! Calculate base flow between maximum allowed water table depth and
  ! "infinity":
  !-----------------------------------------------------------------------
  qbase_min(i)=ksfz(i,nshyd+1) * EXP(-fexp(i)*(zw_max-zdepth(nshyd)))

  !-----------------------------------------------------------------------
  ! Calculate maximum possible and actual base flow for each layer:
  !-----------------------------------------------------------------------

  DO n=1,nshyd

    qbase_max_l(i,n)=ksfz(i,n)*(zdepth(n)-zdepth(n-1))

    IF(zw(i) <= zdepth(n-1)) THEN
      qbase_l(i,n)=qbase_max_l(i,n)
    END IF

    IF(zw(i) <  zdepth(n).AND.zw(i) >  zdepth(n-1)) THEN
      qbase_l(i,n)=ksfz(i,n)*(zdepth(n)-zw(i))
    END IF

    IF(n == 1.AND.zw(i) <  zdepth(n)) THEN
      qbase_l(i,n)=ksfz(i,n)*(zdepth(n)-zw(i))
    END IF

  END DO

  qbase_max_l(i,nshyd+1)=ksfz(i,nshyd+1) - qbase_min(i)

  IF(zw(i) <= zdepth(nshyd))THEN
    qbase_l(i,nshyd+1)=qbase_max_l(i,nshyd+1)
  ELSE
    qbase_l(i,nshyd+1)=ksfz(i,nshyd+1)                                        &
                      *EXP(-fexp(i)*(zw(i)-zdepth(nshyd)))                    &
                      -qbase_min(i)
  END IF

!-----------------------------------------------------------------------
! Calculate total possible and actual base flow:
!-----------------------------------------------------------------------
  DO n=1,nshyd+1
    qbase_l(i,n)=MAX(0.0,qbase_l(i,n))
    qbase(i)=qbase(i)+qbase_l(i,n)
    qbase_max(i)=qbase_max(i)+qbase_max_l(i,n)
  END DO

!-----------------------------------------------------------------------
! Calculate critical topographic index.
!-----------------------------------------------------------------------
  IF (qbase(i) >  qbase_max(i)) THEN
    qbase(i)=qbase_max(i)
  END IF

!Check that QBASE_MAX(I)/QBASE(I) will not underflow.
  IF(qbase_max(i) >  EPSILON(qbase_max(i)) .AND.                              &
     qbase(i) >  qbase_max(i)*(EPSILON(qbase(i)))) THEN
    top_crit(i) = LOG(qbase_max(i)/qbase(i))
  END IF

!The checks from the version ported from the UM repo are slightly different.
!Copied here just in case the above checks ever fail

!   !Check that QBASE_MAX(I)/QBASE(I) will not underflow.
!   SELECT CASE (model_type)
!   CASE (mt_single_column)
!
!     IF ( qbase_max(i) > 1e-20 .AND.                            &
!          qbase(i)     > qbase_max(i) * 1e-20) THEN
!       top_crit(i) = LOG( qbase_max(i)/qbase(i) )
!     END IF
!
!   CASE DEFAULT
!     IF ( qbase_max(i) > EPSILON(qbase_max(i)) .AND.            &
!          qbase(i)     > qbase_max(i)*(EPSILON(qbase(i)))) THEN
!       top_crit(i) = LOG( qbase_max(i) / qbase(i) )
!     END IF
!
!   END SELECT

  END DO

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END SUBROUTINE calc_baseflow
END MODULE calc_baseflow_jls_mod
