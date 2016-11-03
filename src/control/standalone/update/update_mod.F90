#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


MODULE update_mod

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Module variables
!-----------------------------------------------------------------------------
  LOGICAL ::                                                                  &
    l_imogen = .FALSE.
      ! Switch for using IMOGEN to provide driving data
      ! This is read in JULES DRIVE

  INTEGER ::                                                                  &
    io_precip_type,                                                           &
      !  Flag indicating how precipitation is input
      !      1 = total precipitation is read in
      !      2 = values for total rainfall and total snowfall are read in
      !      3 = values for large-scale rainfall, convective rainfall and
      !          total snowfall are read in
      !      4 = values for convective rainfall, large-scale rainfall,
      !          convective snowfall and large-scale snowfall are read in
    io_rad_type,                                                              &
      !  Flag indicating how radiation is input
      !      1 downward fluxes provided
      !      2 net (all wavelength) longwave flux and downward shortwave flux
      !        are provided
      !      3 net downward fluxes are provided
    precip_disagg_method = 2
      !  Flag indicating how precip will be disaggregated
      !      1 no disaggregation
      !      2 imogen method
      !      3 as imogen method except no redistribution when precip is
      !        over max_precip_rate
      !      4 as above, but wet timesteps are distributed randomly through
      !        day, rather than in sequence

  LOGICAL ::                                                                  &
    io_wind_speed,                                                            &
      !   T means that the windspeed is input
      !   F means 2 components of wind are input
    use_diff_rad,                                                             &
      !   T means diffuse radiation is input
      !   F means diffuse radiation is set to a constant fraction
      !     (diffFracConst)
    l_daily_disagg = .FALSE.,                                                 &
      !   T means expect daily driving data and disaggregate it
      !   F means diaggregationonly happens if l_imogen=T
    l_disagg_const_rh = .FALSE.
      !   Switch controlling sub-daily disaggregation of humidity.
      !   T means keep relative humidity constant.
      !   F means keep specific humidity constant.

  REAL ::                                                                     &
    t_for_snow      = 274.0,                                                  &
      !   air temperature (K) at or below which precipitation is assumed to
      !   be snow
    t_for_con_rain  = 373.15,                                                 &
      !   air temperature (K) at or above which rainfall is taken to
      !   be convective (rather than large-scale) in nature
    diff_frac_const = 0.0,                                                    &
      !   a constant value for fraction of radiation that is diffuse
    dur_ls_rain     = 3600.0,                                                 &
      !   duration of large scale rain event (seconds)
    dur_conv_rain   = 21600.0,                                                &
      !   duration of convective rain event (seconds)
    dur_ls_snow     = 3600.0,                                                 &
      !   duration of large scale snow event (seconds)
    dur_conv_snow   = 3600.0
      !   duration of convective snow event (seconds)


  LOGICAL ::                                                                  &
    have_prescribed_veg = .FALSE.
      !   T - a vegetation variable has been prescribed that requires a call
      !       to sparm after update
      !   F - no such variables have been prescribed

CONTAINS


! Fortran INCLUDE statements would be preferred, but (at least) the pgf90
! compiler objects to their use when the included file contains pre-processor
! directives. At present, such directives are used to exclude files from
! the UM build, so are required. This may change in the future.
#include "update_derived_variables.inc"
#include "update_precip_variables.inc"
#include "fill_disaggregated_precip_arrays.inc"
#include "impose_diurnal_cycle.inc"
#include "calc_downward_rad.inc"
#include "assign_irrig_fraction.inc"
#include "update_irrig_variables.inc"

END MODULE update_mod
#endif
