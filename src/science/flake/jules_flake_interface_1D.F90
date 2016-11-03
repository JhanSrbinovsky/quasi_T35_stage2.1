#if defined(UM_JULES)
!-----------------------------------------------------------
! FLake is freely available under the terms of the MIT license.
!
! Copyright (c) 2016
!
! Permission is hereby granted, free of charge, to any person obtaining
! a copy of this software and associated documentation files (the
! "Software"), to deal in the Software without restriction, including
! without limitation the rights to use, copy, modify, merge, publish,
! distribute, sublicense, and/or sell copies of the Software, and to
! permit persons to whom the Software is furnished to do so, subject to
! the following conditions:
!
! The above copyright notice and this permission notice shall be
! included in all copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
! MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
! NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
! LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
! OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
! WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!-----------------------------------------------------------

SUBROUTINE flake_interface (                                                  &
!    .       .       .     !..................!
!                           # points in arrays
!    .       .       .     !..................!
                               npts                                           &
                              ,tile_pts                                       &
                              ,tile_index                                     &
!    .       .       .     !..................!
!                           heat, light and momentum fluxes
!    .       .       .     !..................!
                             , u_star_w_inout                                 &
                             , Q_w_inout                                      &
                             , dwsw_inout                                     &
!    .       .       .     !..................!
!                           "fixed" info
!    .       .       .     !..................!
                             , depth_w_inout                                  &
                             , fetch_inout                                    &
                             , par_Coriolis_inout                             &
                             , del_time_inout                                 &
!    .       .       .     !..................!
!                           lake properties
!    .       .       .     !..................!
                             , albedo_inout                                   &
                             , T_snow_inout                                   &
                             , T_ice_inout                                    &
                             , T_mnw_inout                                    &
                             , T_wML_inout                                    &
                             , C_T_inout                                      &
                             , h_snow_inout                                   &
                             , h_ice_inout                                    &
                             , h_ML_inout                                     &
                             , T_sfc_inout                                    &
                             , TS1_inout                                      &
                             , G_DT_inout                                     &
!    .       .       .     !..................!
!                           trap counters
!    .       .       .     !..................!
                             , trap_frozen                                    &
                             , trap_unfrozen )

!------------------------------------------------------------------------------
!
! Description:
!
!  The FLake interface is
!  a communication routine between "flake_driver"
!  and a prediction system that uses FLake.
!  It assigns the FLake variables at the previous time step
!  to their input values given by the driving model,
!  calls a number of routines to compute the heat and radiation fluxes,
!  calls "flake_driver",
!  and returns the updated FLake variables to the driving model.
!  The "flake_interface" does not contain any Flake physics.
!  It only serves as a convenient means to organize calls of "flake_driver"
!  and of external routines that compute heat and radiation fluxes.
!  The interface may (should) be changed so that to provide
!  the most convenient use of FLake.
!  Within a 3D atmospheric prediction system,
!  "flake_driver" may be called in a DO loop within "flake_interface"
!  for each grid-point where a lake is present.
!  In this way, the driving atmospheric model should call "flake_interface"
!  only once, passing the FLake variables to "flake_interface" as 2D fields.
!
!  Lines embraced with "!_tmp" contain temporary parts of the code.
!  These should be removed prior to using FLake in applications.
!  Lines embraced/marked with "!_dev" may be replaced
!  as improved parameterizations are developed and tested.
!  Lines embraced/marked with "!_dm" are DM's comments
!  that may be helpful to a user.

!
!
! Current Code Owner: DWD, Dmitrii Mironov
!  Phone:  +49-69-8062 2705
!  Fax:    +49-69-8062 3721
!  E-mail: dmitrii.mironov@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.00       2005/11/17 Dmitrii Mironov
!  Initial release
! 2.00       2007/8/24  Gabriel Rooney
!!!----------------------------------------------------------------------------
!!! Altering this routine quite a bit to simplify the calling,
!!! and remove a bit of functionality which either doesn't work yet
!!! e.g. snow, or for which there is insufficient information
!!! e.g. bottom sediment stuff (switched off in FLAKE_CONFIGURE.f).
!!! See FLake documentation, and also "useful hints" for more info.
!!!
!!!
!!! G.Rooney 24/8/2007
!!!
!!! Another note:
!!! The KINDs ireal and iinteger make it not possible to pass variables
!!! e.g. TIMESTEP directly into the corresponding variables here
!!! (at least I assume this was the cause of the problem).
!!! Defining ordinary REALs and then setting one equal to the other seems
!!! to work, for whatever reason.
!!! G.Rooney 3/9/2007
!!!----------------------------------------------------------------------------

!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Declarations:
!
! Modules used:

USE data_parameters , ONLY :                                                  &
    ireals,                  & ! KIND-type parameter for real variables
    iintegers                  ! KIND-type parameter for "normal" integer variables

USE flake_derivedtypes         ! Definitions of several derived TYPEs

USE flake_parameters , ONLY :                                                 &
  tpl_T_f                     , & ! Fresh water freezing point [K]
  tpl_rho_w_r                 , & ! Maximum density of fresh water [kg m^{-3}]
  tpl_T_r                     , & ! Temperature of maximum density of fresh water [K]
  h_Snow_min_flk              , & ! Minimum snow thickness [m]
  h_Ice_min_flk                   ! Minimum ice thickness [m]

USE flake_paramoptic_ref       ! Reference values of the optical characteristics
                               ! of the lake water, lake ice and snow

USE flake_albedo_ref           ! Reference values the albedo for the lake water, lake ice and snow

USE flake           , ONLY :                                                  &
  T_snow_p_flk, T_snow_n_flk  , & ! Temperature at the air-snow interface [K]
  T_ice_p_flk, T_ice_n_flk    , & ! Temperature at the snow-ice or air-ice interface [K]
  T_mnw_p_flk, T_mnw_n_flk    , & ! Mean temperature of the water column [K]
  T_wML_p_flk, T_wML_n_flk    , & ! Mixed-layer temperature [K]
  T_bot_p_flk, T_bot_n_flk    , & ! Temperature at the water-bottom sediment interface [K]
  T_B1_p_flk, T_B1_n_flk      , & ! Temperature at the bottom of the upper layer of the sediments [K]
  C_T_p_flk, C_T_n_flk        , & ! Shape factor (thermocline)
  h_snow_p_flk, h_snow_n_flk  , & ! Snow thickness [m]
  h_ice_p_flk, h_ice_n_flk    , & ! Ice thickness [m]
  h_ML_p_flk, h_ML_n_flk      , & ! Thickness of the mixed-layer [m]
  H_B1_p_flk, H_B1_n_flk      , & ! Thickness of the upper layer of bottom sediments [m]
                                  !
  Q_snow_flk                  , & ! Heat flux through the air-snow interface [W m^{-2}]
  Q_ice_flk                   , & ! Heat flux through the snow-ice or air-ice interface [W m^{-2}]
  Q_w_flk                     , & ! Heat flux through the ice-water or air-water interface [W m^{-2}]
!
  Q_bot_flk                   , & ! Heat flux through the water-bottom sediment interface [W m^{-2}]
  I_atm_flk                   , & ! Radiation flux at the lower boundary of the atmosphere [W m^{-2}],
                                  ! i.e. the incident radiation flux with no regard for the surface albedo
  I_snow_flk                  , & ! Radiation flux through the air-snow interface [W m^{-2}]
  I_ice_flk                   , & ! Radiation flux through the snow-ice or air-ice interface [W m^{-2}]
  I_w_flk                     , & ! Radiation flux through the ice-water or air-water interface [W m^{-2}]
  I_h_flk                     , & ! Radiation flux through the mixed-layer-thermocline interface [W m^{-2}]
  I_bot_flk                   , & ! Radiation flux through the water-bottom sediment interface [W m^{-2}]
  I_intm_0_h_flk              , & ! Mean radiation flux over the mixed layer [W m^{-1}]
  I_intm_h_D_flk              , & ! Mean radiation flux over the thermocline [W m^{-1}]
  Q_star_flk                  , & ! A generalized heat flux scale [W m^{-2}]
  u_star_w_flk                , & ! Friction velocity in the surface layer of lake water [m s^{-1}]
!
  w_star_sfc_flk              , & ! Convective velocity scale, using a generalized heat flux scale [m s^{-1}]
  dMsnowdt_flk                    ! The rate of snow accumulation [kg m^{-2} s^{-1}]

!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations

! Interface variables

INTEGER ::   npts     & ! number of points in passed-down arrays
           , tile_pts   ! number of points with non-zero lake fraction

INTEGER, DIMENSION (npts) ::  tile_index  ! index of points with non-zero lake fraction

REAL                   ::  del_time_inout     & ! The model time step [s]
                         , point_albedo         ! albedo at a point, to fill in lake_albedo array

REAL, DIMENSION (npts) ::  Q_w_inout          & ! Heat flux through the ice-water or air-water interface [W m^{-2}]
                         , u_star_w_inout     & ! Friction velocity in the surface layer of lake water [m s^{-1}]
                         , dwsw_inout         & ! Downwelling shortwave irradiance [W m^{-2}]
                         , depth_w_inout      & ! lake depth [m]
                         , fetch_inout        & ! Typical wind fetch [m]
                         , par_Coriolis_inout & ! Coriolis parameter [s^{-1}]
                         , albedo_inout       & ! lake albedo (snow, ice or water)
                         , T_snow_inout       & ! Temperature at the air-snow interface [K]
                         , T_ice_inout        & ! Temperature at the snow-ice or air-ice interface [K]
                         , T_mnw_inout        & ! Mean temperature of the water column [K]
                         , T_wML_inout        & ! Mixed-layer temperature [K]
                         , C_T_inout          & ! Shape factor (thermocline)
                         , h_snow_inout       & ! Snow thickness [m]
                         , h_ice_inout        & ! Ice thickness [m]
                         , h_ML_inout         & ! Thickness of the mixed-layer [m]
                         , T_sfc_inout        & ! Surface temperature [K]
                         , TS1_inout          & ! Lake temperature at 1st 'soil' level [K]
                         , G_DT_inout           ! Ground heat flux over delta T [W m-2 K-1]

!  Input (procedure arguments)

REAL (KIND = ireals) ::                                                       &
  dMsnowdt_in                       , & ! The rate of snow accumulation [kg m^{-2} s^{-1}]
  I_atm_in                          , & ! Solar radiation flux at the surface [W m^{-2}]
  Q_atm_lw_in                       , & ! Long-wave radiation flux from the atmosphere [W m^{-2}]
  height_u_in                       , & ! Height above the lake surface where the wind speed is measured [m]
  height_tq_in                      , & ! Height where temperature and humidity are measured [m]
  U_a_in                            , & ! Wind speed at z=height_u_in [m s^{-1}]
  T_a_in                            , & ! Air temperature at z=height_tq_in [K]
  q_a_in                            , & ! Air specific humidity at z=height_tq_in
  P_a_in                                ! Surface air pressure [N m^{-2} = kg m^{-1} s^{-2}]

REAL (KIND = ireals) ::                                                       &
  depth_w                           , & ! The lake depth [m]
  fetch                             , & ! Typical wind fetch [m]
  depth_bs                          , & ! Depth of the thermally active layer of the bottom sediments [m]
  T_bs                              , & ! Temperature at the outer edge of
                                        ! the thermally active layer of the bottom sediments [K]
  par_Coriolis                      , & ! The Coriolis parameter [s^{-1}]
  del_time                              ! The model time step [s]

REAL (KIND = ireals) ::                                                       &
  T_snow_in                        , & ! Temperature at the air-snow interface [K]
  T_ice_in                         , & ! Temperature at the snow-ice or air-ice interface [K]
  T_mnw_in                         , & ! Mean temperature of the water column [K]
  T_wML_in                         , & ! Mixed-layer temperature [K]
  T_bot_in                         , & ! Temperature at the water-bottom sediment interface [K]
  T_B1_in                          , & ! Temperature at the bottom of the upper layer of the sediments [K]
  C_T_in                           , & ! Shape factor (thermocline)
  h_snow_in                        , & ! Snow thickness [m]
  h_ice_in                         , & ! Ice thickness [m]
  h_ML_in                          , & ! Thickness of the mixed-layer [m]
  H_B1_in                          , & ! Thickness of the upper layer of bottom sediments [m]
  T_sfc_p                              ! Surface temperature at the previous time step [K]

!  Input/Output (procedure arguments)

REAL (KIND = ireals) ::                                                       &
  albedo_water                        , & ! Water surface albedo with respect to the solar radiation
  albedo_ice                          , & ! Ice surface albedo with respect to the solar radiation
  albedo_snow                             ! Snow surface albedo with respect to the solar radiation

TYPE (opticpar_medium) ::                                                     &
  opticpar_water                       , & ! Optical characteristics of water
  opticpar_ice                         , & ! Optical characteristics of ice
  opticpar_snow                            ! Optical characteristics of snow

!  Output (procedure arguments)

REAL (KIND = ireals) ::                                                       &
  T_snow_out                        , & ! Temperature at the air-snow interface [K]
  T_ice_out                         , & ! Temperature at the snow-ice or air-ice interface [K]
  T_mnw_out                         , & ! Mean temperature of the water column [K]
  T_wML_out                         , & ! Mixed-layer temperature [K]
  T_bot_out                         , & ! Temperature at the water-bottom sediment interface [K]
  T_B1_out                          , & ! Temperature at the bottom of the upper layer of the sediments [K]
  C_T_out                           , & ! Shape factor (thermocline)
  h_snow_out                        , & ! Snow thickness [m]
  h_ice_out                         , & ! Ice thickness [m]
  h_ML_out                          , & ! Thickness of the mixed-layer [m]
  H_B1_out                          , & ! Thickness of the upper layer of bottom sediments [m]
  T_sfc_n                               ! Updated surface temperature [K]


INTEGER ::   ipts, k      ! loop counters
INTEGER ::   trap_frozen, trap_unfrozen ! trap counters

INTEGER iopt ! for optic implied-do loop

TYPE (opticpar_medium), PARAMETER ::                                          &
  opticpar_windermere = opticpar_medium(1,                     & ! no of bands
    (/1._ireals, (0._ireals,iopt=2,nband_optic_max)/),         & ! fractions
    (/0.36_ireals, (1.E+10_ireals,iopt=2,nband_optic_max)/))     ! coeffs

REAL, PARAMETER :: less_small = 0.02 ! less small than c_small


!==============================================================================
!  Start calculations
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!  loop over elements of input arrays
!------------------------------------------------------------------------------

DO k=1,tile_pts
  ipts = tile_index(k)

!------------------------------------------------------------------------------
!  Set initial values from calling variables
!------------------------------------------------------------------------------
del_time     = del_time_inout
u_star_w_flk = u_star_w_inout(ipts)
Q_w_flk      = Q_w_inout(ipts)
I_atm_flk    = dwsw_inout(ipts)
depth_w      = depth_w_inout(ipts)
fetch        = fetch_inout(ipts)
par_Coriolis = par_Coriolis_inout(ipts)
T_snow_p_flk = T_snow_inout(ipts)                                             &
               + less_small ! *** c_small in FLake is too small
!                           ! *** for the melting if-test to work
!                           ! *** so nudging up the surface snow T
!                           ! *** is equivalent to using a bigger c_small
T_ice_p_flk  = T_ice_inout(ipts)
T_mnw_p_flk  = T_mnw_inout(ipts)
T_wML_p_flk  = T_wML_inout(ipts)
C_T_p_flk    = C_T_inout(ipts)
if (h_snow_inout(ipts) >  h_Snow_min_flk) then
  h_snow_p_flk = h_snow_inout(ipts)
else
  h_snow_p_flk = 0.0
endif
h_ice_p_flk  = h_ice_inout(ipts)
h_ML_p_flk   = h_ML_inout(ipts)
T_sfc_p      = T_sfc_inout(ipts)

!------------------------------------------------------------------------------
!  Set albedo to be the same as that used in JULES
!------------------------------------------------------------------------------

  albedo_water = albedo_inout(ipts)
  albedo_ice   = albedo_inout(ipts)
  albedo_snow  = albedo_inout(ipts)

!------------------------------------------------------------------------------
!  Set optical characteristics of the lake water, lake ice and snow
!------------------------------------------------------------------------------

opticpar_water = opticpar_water_ref
opticpar_ice   = opticpar_ice_opaque   ! Opaque ice
opticpar_snow  = opticpar_snow_opaque  ! Opaque snow

!------------------------------------------------------------------------------
!  Turn off the snow / zero snowfall within FLake
!------------------------------------------------------------------------------

dMsnowdt_flk =   0.0     ! zero value turns off snow module

!------------------------------------------------------------------------------
!  fix the switched-off bottom sediment variables
!------------------------------------------------------------------------------
T_bot_p_flk = tpl_T_r
T_B1_p_flk  = tpl_T_r
T_bs        = tpl_T_r
H_B1_p_flk  =   0.0
depth_bs    =   0.0

!------------------------------------------------------------------------------
!  Compute solar radiation fluxes (positive downward)
!------------------------------------------------------------------------------
! DEPENDS ON: flake_radflux
 CALL flake_radflux ( depth_w, albedo_water, albedo_ice, albedo_snow,         &
                      opticpar_water, opticpar_ice, opticpar_snow )

!------------------------------------------------------------------------------
!  determine fluxes and albedo
!------------------------------------------------------------------------------
IF(h_ice_p_flk >= h_Ice_min_flk) THEN            ! Ice exists
  IF(h_snow_p_flk >= h_Snow_min_flk) THEN        ! There is snow above the ice
    Q_snow_flk = Q_w_flk
    Q_ice_flk  = 0._ireals
    Q_w_flk    = 0._ireals
    point_albedo = albedo_snow
  ELSE                                           ! No snow above the ice
    Q_snow_flk = 0._ireals
    Q_ice_flk  = Q_w_flk
    Q_w_flk    = 0._ireals
    point_albedo = albedo_ice
  END IF
ELSE                                             ! No ice-snow cover
    Q_snow_flk = 0._ireals
    Q_ice_flk  = 0._ireals
    point_albedo = albedo_water
END IF

!------------------------------------------------------------------------------
!  Advance FLake variables
!------------------------------------------------------------------------------

! DEPENDS ON: flake_driver
CALL flake_driver ( depth_w, depth_bs, T_bs, par_Coriolis,                    &
                    opticpar_water%extincoef_optic(1),                        &
                    del_time, T_sfc_p, T_sfc_n )

!------------------------------------------------------------------------------
!  Set output values
!------------------------------------------------------------------------------

T_snow_out = T_snow_n_flk
T_ice_out  = T_ice_n_flk
T_mnw_out  = T_mnw_n_flk
T_wML_out  = T_wML_n_flk
T_bot_out  = T_bot_n_flk
T_B1_out   = T_B1_n_flk
C_T_out    = C_T_n_flk
h_snow_out = h_snow_n_flk
h_ice_out  = h_ice_n_flk
h_ML_out   = h_ML_n_flk
H_B1_out   = H_B1_n_flk

!------------------------------------------------------------------------------
!  reset calling variables from output values
!------------------------------------------------------------------------------

T_snow_inout(ipts) = T_snow_out
T_ice_inout(ipts)  = T_ice_out
T_mnw_inout(ipts)  = T_mnw_out
T_wML_inout(ipts)  = T_wML_out
C_T_inout(ipts)    = C_T_out
h_snow_inout(ipts) = h_snow_out
h_ice_inout(ipts)  = h_ice_out
h_ML_inout(ipts)   = h_ML_out
T_sfc_inout(ipts)  = T_sfc_n
albedo_inout(ipts) = point_albedo

IF ( h_ice_out >  0.0 ) THEN
  IF(ABS(tpl_T_f - TS1_inout( ipts )) < EPSILON(0.0)) THEN
    trap_frozen = trap_frozen + 1
    G_DT_inout(ipts) = 1.0
  ELSE
    G_DT_inout(ipts) = abs( (   (1-point_albedo)*dwsw_inout(ipts)             &
                              + Q_w_inout(ipts) )                             &
                            / ( tpl_T_f - TS1_inout( ipts ) ) )
  END IF
ELSE
  IF(ABS(T_sfc_n - TS1_inout( ipts )) < EPSILON(0.0)) THEN
    trap_unfrozen = trap_unfrozen + 1
    G_DT_inout(ipts) = 1.0
  ELSE
    G_DT_inout(ipts) = abs( (   (1-point_albedo)*dwsw_inout(ipts)             &
                              + Q_w_inout(ipts) )                             &
                            / ( T_sfc_n - TS1_inout( ipts ) ) )
  END IF
END IF
!------------------------------------------------------------------------------
!  end loop over npts
!------------------------------------------------------------------------------
  ENDDO


!------------------------------------------------------------------------------
!  End calculations
!==============================================================================

END SUBROUTINE flake_interface
#endif
