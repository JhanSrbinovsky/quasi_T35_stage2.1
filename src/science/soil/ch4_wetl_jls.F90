! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE CH4_WETL-----------------------------------------------

! Description:
!     Calculates methane emissions from wetland area.

! Subroutine Interface:
SUBROUTINE ch4_wetl(npnts,soil_pts,dim_cs1,soil_index,l_triffid               &
  ,tsoil_d,cs_ch4,cs,resp_s,npp,f_wetl,fch4_wetl                              &
  ,fch4_wetl_cs,fch4_wetl_npp,fch4_wetl_resps                                 &
   )

USE c_0_dg_c, ONLY :                                                          &
!      imported scalar parameters
   tm

USE c_ch4, ONLY :                                                             &
!      imported scalar parameters
   t0_ch4,const_ch4_cs,q10_ch4_cs,const_ch4_npp,q10_ch4_npp                   &
  ,const_ch4_resps,q10_ch4_resps

use jules_surface_mod, ONLY :                                                 &
!      imported arrays with intent(in)
   kaps_roth

USE jules_hydrology_mod, ONLY: l_wetland_ch4_npp

#if defined(UM_JULES)
USE jules_vegetation_mod, ONLY:i_veg_vn,i_veg_vn_1b,i_veg_vn_2b
#endif

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

INTEGER, INTENT(IN) ::                                                        &
 npnts                                                                        &
                     ! IN Number of gridpoints.
,soil_pts                                                                     &
                     ! IN Number of soil points.
,dim_cs1                                                                      &
                     ! IN Number of soil carbon pools
,soil_index(npnts)   ! IN Array of soil points.

LOGICAL, INTENT(IN) ::                                                        &
 l_triffid           ! TRUE if using TRIFFID

REAL, INTENT(IN) ::                                                           &
 tsoil_d(npnts)                                                               &
                     ! IN Diagnosed soil temp to 1 metre (K).
,cs_ch4(npnts)                                                                &
                     ! IN Soil carbon used in CH4 wetlands
!                          !   if TRIFFID is switched off
!                          !   (kg C/m2).
,cs(npnts,dim_cs1)                                                            &
                     ! IN Soil carbon
!                          !   For RothC (dim_cs1=4), the pools
!                          !   are DPM, RPM, biomass and humus.
!                          !   (kg C/m2).
,resp_s(npnts,dim_cs1)                                                        &
                     ! IN Soil respiration in pools (kg C/m2/s).
,npp(npnts)                                                                   &
                     ! IN Gridbox mean net primary
!                          !    productivity (kg C/m2/s).
,f_wetl(npnts)       ! IN Wetland fraction

REAL, INTENT(INOUT) ::                                                        &
 fch4_wetl(npnts)                                                             &
                     ! OUT Scaled methane flux (default substrate
!                          !    used in atmos chem) (10^-9 kg C/m2/s)
,fch4_wetl_cs(npnts)                                                          &
!                     ! OUT Scaled methane flux (soil carbon
!                          !    substrate) (kg C/m2/s)
,fch4_wetl_npp(npnts)                                                         &
                     ! OUT Scaled methane flux (npp
!                          !    substrate) (kg C/m2/s)
,fch4_wetl_resps(npnts)
                     ! OUT Scaled methane flux (soil respiration
!                          !    substrate) (kg C/m2/s)

REAL ::                                                                       &
 q10t_ch4_cs                                                                  &
                     ! Q10 value at T
,q10t_ch4_npp                                                                 &
                     ! Q10 value at T
,q10t_ch4_resps                                                               &
                     ! Q10 value at T
,const_tdep_cs                                                                &
!                    ! T and Q10(0) dependent function
,const_tdep_npp                                                               &
!                    ! T and Q10(0) dependent function
,const_tdep_resps
!                    ! T and Q10(0) dependent function

INTEGER ::                                                                    &
 i,j,k

REAL ::                                                                       &
 sumkaps             ! Sum of kaps_roth values.

!     Local aray variables.
REAL ::                                                                       &
 cs_eff(npnts)                                                                &
                     ! Effective soil carbon (kg C/m2).
,resp_s_tot(npnts)
                     ! Soil respiration total (kg C/m2/s).

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CH4_WETL'

!-----------------------------------------------------------------------
!     Calculate an effective soil carbon for wetland methane emission.
!-----------------------------------------------------------------------
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

#if defined(UM_JULES)
!----------------------------------------------------------------------------
! set value of CONST_CH4_CS dependent upon whether running  1B or 2B veg scheme
! ideally this should set in future as a JULES input via namelist
! ---------------------------------------------------------------------------
IF (i_veg_vn == i_veg_vn_1b)  CONST_CH4_CS        = 5.41E-12
IF (i_veg_vn == i_veg_vn_2b)  CONST_CH4_CS        = 5.41E-10
#else
CONST_CH4_CS = 7.41E-12
#endif


IF (l_triffid) THEN
  sumkaps = 0.0
  DO k=1,dim_cs1
    sumkaps = sumkaps + kaps_roth(k)
  END DO
!       Weight each pool by specific respiration rate.
  DO j=1,soil_pts
    i=soil_index(j)
    cs_eff(i) = 0.0
    DO k=1,dim_cs1
      cs_eff(i) = cs_eff(i) + cs(i,k) * kaps_roth(k)
    END DO
    cs_eff(i) = cs_eff(i) / sumkaps
  END DO
ELSE
!       Use the single soil carbon pool.
  DO j=1,soil_pts
    i=soil_index(j)
    cs_eff(i) = cs_ch4(i)
  END DO
END IF  !  l_triffid

DO j=1,soil_pts
  i=soil_index(j)
  resp_s_tot(i)=0.0
  DO k=1,dim_cs1
    IF(resp_s(i,k) > 0.0)                                                     &
      resp_s_tot(i)=resp_s_tot(i)+resp_s(i,k)
  ENDDO
ENDDO
!-----------------------------------------------------------------------
!     Calculate scaled wetland methane emission.
!-----------------------------------------------------------------------
const_tdep_cs = t0_ch4 * LOG(q10_ch4_cs)
const_tdep_npp = t0_ch4 * LOG(q10_ch4_npp)
const_tdep_resps = t0_ch4 * LOG(q10_ch4_resps)
DO j=1,soil_pts
  i=soil_index(j)
  IF( tsoil_d(i) > tm .AND. f_wetl(i) > 0.0 )THEN
    q10t_ch4_cs = EXP(const_tdep_cs/tsoil_d(i))
    q10t_ch4_npp = EXP(const_tdep_npp/tsoil_d(i))
    q10t_ch4_resps = EXP(const_tdep_resps/tsoil_d(i))
    fch4_wetl_cs(i) = const_ch4_cs*cs_eff(i)*                                 &
                   f_wetl(i)*q10t_ch4_cs**(0.1*(tsoil_d(i)-t0_ch4))
    IF(npp(i) > 0.0)THEN
      fch4_wetl_npp(i) = const_ch4_npp*npp(i)*                                &
                   f_wetl(i)*q10t_ch4_npp**(0.1*(tsoil_d(i)-t0_ch4))
    END IF
    IF(resp_s_tot(i) > 0.0)THEN
      fch4_wetl_resps(i) = const_ch4_resps*resp_s_tot(i)*                     &
                   f_wetl(i)*q10t_ch4_resps**(0.1*(tsoil_d(i)-t0_ch4))
    END IF

  END IF
END DO

IF(l_wetland_ch4_npp)THEN
  fch4_wetl(:)=1.e9*fch4_wetl_npp(:)
ELSE
  fch4_wetl(:)=1.e9*fch4_wetl_cs(:)
END IF

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ch4_wetl
