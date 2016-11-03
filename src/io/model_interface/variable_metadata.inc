#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!-----------------------------------------------------------------------------
! DATA statement for variable metadata
!-----------------------------------------------------------------------------


DATA metadata/                                                                &
!-----------------------------------------------------------------------------
! Metadata for latitude
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'latitude',                                                               &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox latitude",                                                       &
! Units
    "degrees"                                                                 &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for longitude
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'longitude',                                                              &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox longitude",                                                      &
! Units
    "degrees"                                                                 &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for land_fraction
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'land_fraction',                                                          &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for surf_hgt
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'surf_hgt',                                                               &
! Variable type
    VAR_TYPE_SURFT,                                                            &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for z_land_land
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'z_land_land',                                                            &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for frac
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'frac',                                                                   &
! Variable type
    VAR_TYPE_TYPE,                                                            &
! Long name
    "Fractional cover of each surface type",                                  &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for albsoil
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'albsoil',                                                                &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for soil clay content
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'clay',                                                                   &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Soil clay content",                                                      &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for albobs_sw
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'albobs_sw',                                                              &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for albobs_vis
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'albobs_vis',                                                             &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for albobs_nir
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'albobs_nir',                                                             &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for b
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'b',                                                                      &
! Variable type
    VAR_TYPE_SOIL,                                                            &
! Long name
    "Gridbox Brooks-Corey exponent for each soil layer",                      &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for b_const_z
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'b_const_z',                                                              &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for sathh
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'sathh',                                                                  &
! Variable type
    VAR_TYPE_SOIL,                                                            &
! Long name
    "Gridbox saturated soil water pressure for each soil layer",              &
! Units
    "m"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for sathh_const_z
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'sathh_const_z',                                                          &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for satcon
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'satcon',                                                                 &
! Variable type
    VAR_TYPE_SOIL,                                                            &
! Long name
    "Gridbox saturated hydraulic conductivity for each soil layer",           &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for satcon_const_z
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'satcon_const_z',                                                         &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for sm_sat
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'sm_sat',                                                                 &
! Variable type
    VAR_TYPE_SOIL,                                                            &
! Long name
    "Gridbox volumetric moisture content at saturation for each soil layer",  &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for sm_sat_const_z
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'sm_sat_const_z',                                                         &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for sm_crit
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'sm_crit',                                                                &
! Variable type
    VAR_TYPE_SOIL,                                                            &
! Long name
    "Gridbox volumetric moisture content at critical point for each soil layer", &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for sm_crit_const_z
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'sm_crit_const_z',                                                        &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for sm_wilt
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'sm_wilt',                                                                &
! Variable type
    VAR_TYPE_SOIL,                                                            &
! Long name
    "Gridbox  volumetric moisture content at wilting point for each soil layer", &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for sm_wilt_const_z
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'sm_wilt_const_z',                                                        &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for hcap
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'hcap',                                                                   &
! Variable type
    VAR_TYPE_SOIL,                                                            &
! Long name
    "Gridbox dry soil heat capacity for each soil layer",                     &
! Units
    "J K-1 m-3"                                                               &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for hcap_const_z
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'hcap_const_z',                                                           &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for hcon
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'hcon',                                                                   &
! Variable type
    VAR_TYPE_SOIL,                                                            &
! Long name
    "Gridbox dry soil thermal conductivity for each soil layer",              &
! Units
    "W m-1 K-1"                                                               &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for hcon_const_z
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'hcon_const_z',                                                           &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for fexp
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fexp',                                                                   &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for ti_mean
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'ti_mean',                                                                &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for ti_sig
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'ti_sig',                                                                 &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for wrr
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'wrr',                                                                    &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for hwr
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'hwr',                                                                    &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for hgt
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'hgt',                                                                    &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for ztm
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'ztm',                                                                    &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for disp
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'disp',                                                                   &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for albwl
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'albwl',                                                                  &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for albrd
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'albrd',                                                                  &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for emisw
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'emisw',                                                                  &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for emisr
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'emisr',                                                                  &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for canopy
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'canopy',                                                                 &
! Variable type
    VAR_TYPE_SURFT,                                                            &
! Long name
    "Tile surface/canopy water for snow-free land tiles",                     &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for cs
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'cs',                                                                     &
! Variable type
    VAR_TYPE_SCPOOL,                                                          &
! Long name is dynamic - set in get_var_attrs
    "",                                                                       &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for ns
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'ns',                                                                     &
! Variable type
    VAR_TYPE_SCPOOL,                                                          &
! Long name
    "Gridbox soil Nitrogen in each pool (DPM,RPM,bio,hum)",                   &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for immob_n
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'immob_n',                                                                &
! Variable type
    VAR_TYPE_SCPOOL,                                                          &
! Long name
    "Gridbox soil Nitrogen immobilisation in each pool (DPM,RPM,bio,hum)",    &
! Units
    "kg m-2 per 360days"                                                      &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for immob_n_pot
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'immob_n_pot',                                                            &
! Variable type
    VAR_TYPE_SCPOOL,                                                          &
! Long name
    "Gridbox soil potential Nitrogen immobilisation in each pool (DPM,RPM,bio,hum)", &
! Units
    "kg m-2 per 360days"                                                      &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for minl_n
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'minl_n',                                                                 &
! Variable type
    VAR_TYPE_SCPOOL,                                                          &
! Long name
    "Gridbox soil Nitrogen mineralisation in each pool (DPM,RPM,bio,hum)",    &
! Units
    "kg m-2 per 360days"                                                      &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for minl_n_pot
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'minl_n_pot',                                                             &
! Variable type
    VAR_TYPE_SCPOOL,                                                          &
! Long name
    "Gridbox soil potential Nitrogen mineralisation in each pool (DPM,RPM,bio,hum)", &
! Units
    "kg m-2 per 360days"                                                      &
  ),                                                                          &

!-----------------------------------------------------------------------------
! Metadata for resp_s_diag
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'resp_s_diag',                                                            &
! Variable type
    VAR_TYPE_SCPOOL,                                                          &
! Long name
    "Gridbox soil respiration carbon flux in each pool (DPM,RPM,bio,hum)",    &
! Units
    "kg m-2 per 360days"                                                      &
  ),                                                                          &

!-----------------------------------------------------------------------------
! Metadata for resp_s_pot_diag
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'resp_s_pot_diag',                                                        &
! Variable type
    VAR_TYPE_SCPOOL,                                                          &
! Long name
    "Gridbox potential soil respiration carbon flux in each pool (DPM,RPM,bio,hum)", &
! Units
    "kg m-2 per 360days"                                                      &
  ),                                                                          &

!-----------------------------------------------------------------------------
! Metadata for soil_cn
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'soil_CN',                                                                &
! Variable type
    VAR_TYPE_SCPOOL,                                                          &
! Long name
    "Gridbox Soil C:N in each pool (DPM,RPM,bio,hum)",                        &
! Units
    ":"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for soil_cn_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'soil_CN_gb',                                                             &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox Soil C:N",                                                       &
! Units
    ":"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for resp_s_pot_diag_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'resp_s_pot_diag_gb',                                                     &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox mean potential soil respiration",                                &
! Units
    "kg m-2 per 360days"                                                      &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for resp_s_diag_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'resp_s_diag_gb',                                                         &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox mean soil respiration",                                          &
! Units
    "kg m-2 per 360days"                                                      &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for min_pot_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'minl_n_pot_gb',                                                          &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox mean potential soil N mineralisation",                           &
! Units
    "kg m-2 per 360days"                                                      &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for min_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'minl_n_gb',                                                              &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox mean potential soil N mineralisation",                           &
! Units
    "kg m-2 per 360days"                                                      &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for immob_pot_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'immob_n_pot_gb',                                                         &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox mean potential soil N immobilisation",                           &
! Units
    "kg m-2 per 360days"                                                      &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for immob_n_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'immob_n_gb',                                                             &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox mean potential soil N immobilisation",                           &
! Units
    "kg m-2 per 360days"                                                      &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for ns_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'ns_gb',                                                                  &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox soil Nitrogen(total)",                                           &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for fN
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fN',                                                                     &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox Nitrogen availability rate modifier",                            &
! Units
    ":"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for gs
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'gs',                                                                     &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox surface conductance to evaporation",                             &
! Units
    "m s-1"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_tile
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_tile',                                                              &
! Variable type
    VAR_TYPE_SURFT,                                                            &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for sthuf
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'sthuf',                                                                  &
! Variable type
    VAR_TYPE_SOIL,                                                            &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for t_soil
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    't_soil',                                                                 &
! Variable type
    VAR_TYPE_SOIL,                                                            &
! Long name
    "Gridbox sub-surface temperature of each layer",                          &
! Units
    "K"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for tsoil_deep
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'tsoil_deep',                                                             &
! Variable type
    VAR_TYPE_BEDROCK,                                                         &
! Long name
    "Temperature in bedrock layers",                                          &
! Units
    "K"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for tstar_tile
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'tstar_tile',                                                             &
! Variable type
    VAR_TYPE_SURFT,                                                            &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for lai
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'lai',                                                                    &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT leaf area index",                                                    &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for canht
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'canht',                                                                  &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT canopy height",                                                      &
! Units
    "m"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for sthzw
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'sthzw',                                                                  &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Soil wetness in deep (water table/TOPMODEL) layer",                      &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for zw
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'zw',                                                                     &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox mean depth to water table",                                      &
! Units
    "m"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for rgrain
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'rgrain',                                                                 &
! Variable type
    VAR_TYPE_SURFT,                                                            &
! Long name
    "Tile snow surface grain size",                                           &
! Units
    "microns"                                                                 &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for cv
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'cv',                                                                     &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox mean vegetation carbon at end of model timestep.",               &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for rho_snow
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'rho_snow',                                                               &
! Variable type
    VAR_TYPE_SURFT,                                                            &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_depth
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_depth',                                                             &
! Variable type
    VAR_TYPE_SURFT,                                                            &
! Long name
    "Tile snow depth (on ground)",                                            &
! Units
    "m"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_grnd
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_grnd',                                                              &
! Variable type
    VAR_TYPE_SURFT,                                                            &
! Long name
    "Tile snow on ground below canopy (snow_grnd)",                           &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for nsnow
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'nsnow',                                                                  &
! Variable type
    VAR_TYPE_SURFT,                                                            &
! Long name
   "Number of snow layers",                                                   &
! Units
   "1"                                                                        &
 ),                                                                           &
!-----------------------------------------------------------------------------
! Metadata for snow_ds
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_ds',                                                                &
! Variable type
    VAR_TYPE_SNOW,                                                            &
! Long name
    "Depth of snow layers for each tile",                                     &
! Units
    "m"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_ice
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_ice',                                                               &
! Variable type
    VAR_TYPE_SNOW,                                                            &
! Long name
    "Ice mass in snow layers for each tile",                                  &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_liq
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_liq',                                                               &
! Variable type
    VAR_TYPE_SNOW,                                                            &
! Long name
    "Liquid mass in snow layers for each tile",                               &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for tsnow
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'tsnow',                                                                  &
! Variable type
    VAR_TYPE_SNOW,                                                            &
! Long name
    "Temperature of snow layers for each tile",                               &
! Units
    "K"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for rgrainl
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'rgrainl',                                                                &
! Variable type
    VAR_TYPE_SNOW,                                                            &
! Long name
    "Grain size in snow layers for each tile",                                &
! Units
    "microns"                                                                 &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for pstar
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'pstar',                                                                  &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox surface pressure",                                               &
! Units
    "Pa"                                                                      &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for q
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'q',                                                                      &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for t
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    't',                                                                      &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for rad_net
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'rad_net',                                                                &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Surface net radiation of land points",                                   &
! Units
    "W m-2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for lw_net
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'lw_net',                                                                 &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox net downward longwave radiation at surface",                     &
! Units
    "W m-2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for sw_net
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'sw_net',                                                                 &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox net downward shortwave radiation at surface",                    &
! Units
    "W m-2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for lw_down
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'lw_down',                                                                &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox surface downward LW radiation",                                  &
! Units
    "W m-2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for sw_down
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'sw_down',                                                                &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox surface downward SW radiation",                                  &
! Units
    "W m-2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for diff_rad
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'diff_rad',                                                               &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for precip
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'precip',                                                                 &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox precipitation rate",                                             &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for tot_rain
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'tot_rain',                                                               &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for tot_snow
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'tot_snow',                                                               &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for con_rain
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'con_rain',                                                               &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox convective rainfall",                                            &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for ls_rain
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'ls_rain',                                                                &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox large-scale rainfall",                                           &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for con_snow
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'con_snow',                                                               &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox convective snowfall",                                            &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for ls_snow
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'ls_snow',                                                                &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox large-scale snowfall",                                           &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for wind
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'wind',                                                                   &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox wind speed",                                                     &
! Units
    "m s-1"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for u
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'u',                                                                      &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for v
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'v',                                                                      &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for z1_tq_in
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'z1_tq_in',                                                               &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for diurnal temperature range
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'dt_range',                                                               &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Diurnal temperature range",                                              &
! Units
    "K"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for albedo_land
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'albedo_land',                                                            &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox albedo (as used for net shortwave calculation)",                 &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for canopy_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'canopy_gb',                                                              &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox canopy water content",                                           &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for cs_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'cs_gb',                                                                  &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox soil carbon (total)",                                            &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for depth_frozen
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'depth_frozen',                                                           &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox depth of frozen ground at surface (temperature definition)",     &
! Units
    "m"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for depth_frozen_sthf
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'depth_frozen_sthf',                                                      &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox depth of frozen ground at surface (soil moisture definition)",   &
! Units
    "m"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for depth_unfrozen
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'depth_unfrozen',                                                         &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox depth of unfrozen ground at surface (temeprature definition)",   &
! Units
    "m"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for depth_unfrozen_sthf
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'depth_unfrozen_sthf',                                                    &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox depth of unfrozen ground at surface (soil moisture definition)", &
! Units
    "m"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for drain
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'drain',                                                                  &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Drainage from bottom (nshyd) soil layer",                                &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for elake
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'elake',                                                                  &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox mean evaporation from lakes",                                    &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for emis_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'emis_gb',                                                                &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox emissivity",                                                     &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for fch4_wetl
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fch4_wetl',                                                              &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Scaled methane flux from wetland fraction (for use in atmos chem)",      &
! Units
    "10^-9 kg m-2 s-1"                                                        &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for fch4_wetl_cs
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fch4_wetl_cs',                                                           &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Scaled methane flux from wetland fraction (cs substrate)",               &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for fch4_wetl_npp
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fch4_wetl_npp',                                                          &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Scaled methane flux from wetland fraction (npp substrate)",              &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for fch4_wetl_resps
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fch4_wetl_resps',                                                        &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Scaled methane flux from wetland fraction (soil resp substrate)",        &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for fsat
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fsat',                                                                   &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Surface saturated fraction",                                             &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for fsmc_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fsmc_gb',                                                                &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox soil moisture availability factor (beta)",                       &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for fwetl
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fwetl',                                                                  &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Wetland fraction at end of model timestep.",                             &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for gpp_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'gpp_gb',                                                                 &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox gross primary productivity",                                     &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for hf_snow_melt
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'hf_snow_melt',                                                           &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox snowmelt heat flux",                                             &
! Units
    "W m-2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for land_index
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'land_index',                                                             &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Index (gridbox number) of land points",                                  &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for lice_index
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'lice_index',                                                             &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Index (gridbox number) of land ice points",                              &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for lit_c_mean
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'lit_c_mean',                                                             &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox mean carbon litter",                                             &
! Units
    "kg m-2 per 360days"                                                      &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for lw_up
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'lw_up',                                                                  &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox surface upward LW radiation of land points",                     &
! Units
    "W m-2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for npp_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'npp_gb',                                                                 &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox net primary productivity",                                       &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for qbase
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'qbase',                                                                  &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Baseflow (lateral subsurface runoff)",                                   &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for qbase_zw
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'qbase_zw',                                                               &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Baseflow from deep layer",                                               &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for resp_p_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'resp_p_gb',                                                              &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox plant respiration",                                              &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for resp_s_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'resp_s_gb',                                                              &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox soil respiration (total)",                                       &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for resp_s_dr_out
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'resp_s_dr_out',                                                          &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox mean soil respiration for driving TRIFFID",                      &
! Units
    "kg m-2 per 360days"                                                      &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for n_root
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'n_root',                                                                 &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT Root N content scaled by LAI_BAL in sf_stom",                        &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for n_stem
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'n_stem',                                                                 &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT Stem N content scaled by LAI_BAL in sf_stom",                        &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for n_leaf
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'n_leaf',                                                                 &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT Leaf N content scaled by LAI in sf_stom",                            &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for pc_s
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'pc_s',                                                                   &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT net carbon available for spreading",                                 &
! Units
    "kg m-2 per 360days"                                                      &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for runoff
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'runoff',                                                                 &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox runoff rate",                                                    &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for sat_excess_roff
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'sat_excess_roff',                                                        &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Saturation excess surface ('Dunne') runoff",                             &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for smc_avail_top
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'smc_avail_top',                                                          &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name is dynamic - set in get_var_attrs
    "",                                                                       &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for smc_avail_tot
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'smc_avail_tot',                                                          &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox available moisture in soil column",                              &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for smc_tot
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'smc_tot',                                                                &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox total soil moisture in column",                                  &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snomlt_sub_htf
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snomlt_sub_htf',                                                         &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox sub-canopy snowmelt heat flux",                                  &
! Units
    "W m-2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_can_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_can_gb',                                                            &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox snow on canopy",                                                 &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_depth_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_depth_gb',                                                          &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Snow depth (on ground)",                                                 &
! Units
    "m"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_frac
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_frac',                                                              &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox snow-covered fraction of land points",                           &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_grnd_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_grnd_gb',                                                           &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox snow below canopy (snow_grnd)",                                  &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_ice_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_ice_gb',                                                            &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox ice content of snow on ground",                                  &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_liq_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_liq_gb',                                                            &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox liquid content of snow on ground",                               &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_melt_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_melt_gb',                                                           &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox rate of snowmelt",                                               &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for soil_index
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'soil_index',                                                             &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Index (gridbox number) of soil points",                                  &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for sub_surf_roff
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'sub_surf_roff',                                                          &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox sub-surface runoff",                                             &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for surf_roff
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'surf_roff',                                                              &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox surface runoff",                                                 &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for swet_liq_tot
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'swet_liq_tot',                                                           &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox unfrozen soil moisture as fraction of saturation",               &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for swet_tot
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'swet_tot',                                                               &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox soil moisture as fraction of saturation",                        &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for tfall
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'tfall',                                                                  &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox throughfall",                                                    &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for trad
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'trad',                                                                   &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox effective radiative temperature (assuming emissivity=1)",        &
! Units
    "K"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for c_veg
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'c_veg',                                                                  &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT total carbon content of the vegetation at the end of model timestep.", &
! Units
    "kg m-2"                                                                  &
  ),                                                                       &
!-----------------------------------------------------------------------------
! Metadata for fapar
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fapar',                                                                  &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT Fraction of Absorbed Photosynthetically Active Radiation",           &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for fao_et0
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fao_et0',                                                                &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "FAO Penman-Monteith evapotranspiration for reference crop",              &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for flux_o3_stom
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'flux_o3_stom',                                                           &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "Flux of O3 to stomata",                                                  &
! Units
    "mol m-2 s-1"                                                             &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for fsmc
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fsmc',                                                                   &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT soil moisture availability factor (beta)",                           &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for g_leaf
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'g_leaf',                                                                 &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT leaf turnover rate",                                                 &
! Units
    "per 360days"                                                             &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for g_leaf_day
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'g_leaf_day',                                                             &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT mean leaf turnover rate for input to PHENOL",                        &
! Units
    "per 360days"                                                             &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for g_leaf_dr_out
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'g_leaf_dr_out',                                                          &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT mean leaf turnover rate for driving TRIFFID",                        &
! Units
    "per 360days"                                                             &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for g_leaf_phen
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'g_leaf_phen',                                                            &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT mean leaf turnover rate over phenology period",                      &
! Units
    "per 360days"                                                             &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for gpp
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'gpp',                                                                    &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT gross primary productivity",                                         &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for lai_phen
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'lai_phen',                                                               &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT leaf area index after phenology",                                    &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for lit_c
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'lit_c',                                                                  &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT carbon Litter",                                                      &
! Units
    "kg m-2 per 360days"                                                      &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for npp_dr_out
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'npp_dr_out',                                                             &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT mean NPP for driving TRIFFID",                                       &
! Units
    "kg m-2 per 360days"                                                      &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for npp
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'npp',                                                                    &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT net primary productivity",                                           &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for o3_exp_fac
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'o3_exp_fac',                                                             &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "Ozone exposure factor",                                                  &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for resp_p
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'resp_p',                                                                 &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT plant respiration",                                                  &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for resp_w_dr_out
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'resp_w_dr_out',                                                          &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT mean wood respiration carbon flux for driving TRIFFID",              &
! Units
    "kg m-2 per 360days"                                                      &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for resp_l
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'resp_l',                                                                 &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT leaf maintenance respiration carbon flux (in sf_stom)",              &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for resp_r
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'resp_r',                                                                 &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT root maintenance respiration carbon flux (in sf_stom)",              &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for resp_w
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'resp_w',                                                                 &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT wood maintenance respiration carbon flux (in sf_stom)",              &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for resp_s
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'resp_s',                                                                 &
! Variable type
    VAR_TYPE_SCPOOL,                                                          &
! Long name is dynamic - set in get_var_attrs
    "",                                                                       &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for cosz
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'cosz',                                                                   &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Cosine of the zenith angle",                                             &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for diff_frac
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'diff_frac',                                                              &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox fraction of radiation that is diffuse",                          &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for ecan_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'ecan_gb',                                                                &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox mean evaporation from canopy/surface store",                     &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for ei_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'ei_gb',                                                                  &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox sublimation from lying snow or sea-ice",                         &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for et_stom
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'et_stom',                                                                &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT APPROXIMATE stomatal transpiration",                                 &
! Units
    "W m-2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for et_stom_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'et_stom_gb',                                                             &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox APPROXIMATE stomatal transpiration",                             &
! Units
    "W m-2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for esoil_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'esoil_gb',                                                               &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox surface evapotranspiration from soil moisture store",            &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for fqw_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fqw_gb',                                                                 &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox moisture flux from surface",                                     &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for ftl_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'ftl_gb',                                                                 &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox surface sensible heat flux",                                     &
! Units
    "W m-2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for land_albedo_1
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'land_albedo_1',                                                          &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox albedo for waveband 1",                                          &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for land_albedo_2
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'land_albedo_2',                                                          &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox albedo for waveband 2",                                          &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for land_albedo_3
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'land_albedo_3',                                                          &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox albedo for waveband 3",                                          &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for land_albedo_4
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'land_albedo_4',                                                          &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox albedo for waveband 4",                                          &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for latent_heat
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'latent_heat',                                                            &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox surface latent heat flux",                                       &
! Units
    "W m-2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for q1p5m_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'q1p5m_gb',                                                               &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox specific humidity at 1.5m height",                               &
! Units
    "kg kg-1"                                                                 &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for qw1
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'qw1',                                                                    &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox specific humidity (total water content)",                        &
! Units
    "kg kg-1"                                                                 &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for rainfall
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'rainfall',                                                               &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox rainfall rate",                                                  &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snomlt_surf_htf
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snomlt_surf_htf',                                                        &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox heat flux used for surface melting of snow",                     &
! Units
    "W m-2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snowfall
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snowfall',                                                               &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox snowfall rate",                                                  &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_mass_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_mass_gb',                                                           &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox snowmass (on land this is lying_snow=snow_tile+snow_grnd)",      &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for surf_ht_flux_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'surf_ht_flux_gb',                                                        &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox net downward heat flux at surface over land and sea-ice fraction", &
! Units
    "W m-2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for t1p5m_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    't1p5m_gb',                                                               &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox temperature at 1.5m height",                                     &
! Units
    "K"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for taux1
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'taux1',                                                                  &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox westerly component of surface wind stress",                      &
! Units
    "N m-2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for tauy1
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'tauy1',                                                                  &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox southerly component of surface wind stress",                     &
! Units
    "N m-2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for tl1
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'tl1',                                                                    &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox ice/liquid water temperature",                                   &
! Units
    "K"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for tstar_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'tstar_gb',                                                               &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox surface temperature",                                            &
! Units
    "K"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for u1
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'u1',                                                                     &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox westerly wind component",                                        &
! Units
    "m s-1"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for u10m
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'u10m',                                                                   &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox westerly wind component at 10 m height",                         &
! Units
    "m s-1"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for v1
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'v1',                                                                     &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox southerly wind component",                                       &
! Units
    "m s-1"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for v10m
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'v10m',                                                                   &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox southerly wind component at 10 m height",                        &
! Units
    "m s-1"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for ext
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'ext',                                                                    &
! Variable type
    VAR_TYPE_SOIL,                                                            &
! Long name
    "Gridbox extraction of water from each soil layer",                       &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for smcl
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'smcl',                                                                   &
! Variable type
    VAR_TYPE_SOIL,                                                            &
! Long name
    "Gridbox moisture content of each soil layer",                            &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for soil_wet
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'soil_wet',                                                               &
! Variable type
    VAR_TYPE_SOIL,                                                            &
! Long name
    "Gridbox total moisture content of each soil layer, as fraction of saturation", &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for sthf
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'sthf',                                                                   &
! Variable type
    VAR_TYPE_SOIL,                                                            &
! Long name
    "Gridbox frozen moisture content of each soil layer as a fraction of saturation", &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for sthu
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'sthu',                                                                   &
! Variable type
    VAR_TYPE_SOIL,                                                            &
! Long name
    "Gridbox unfrozen moisture content of each soil layer as a fraction of saturation", &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for alb_tile_1
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'alb_tile_1',                                                             &
! Variable type
    VAR_TYPE_SURFT,                                                            &
! Long name
    "Tile land albedo, waveband 1",                                           &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for alb_tile_2
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'alb_tile_2',                                                             &
! Variable type
    VAR_TYPE_SURFT,                                                            &
! Long name
    "Tile land albedo, waveband 2",                                           &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for alb_tile_3
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'alb_tile_3',                                                             &
! Variable type
    VAR_TYPE_SURFT,                                                            &
! Long name
    "Tile land albedo, waveband 3",                                           &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for alb_tile_4
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'alb_tile_4',                                                             &
! Variable type
    VAR_TYPE_SURFT,                                                            &
! Long name
    "Tile land albedo, waveband 4",                                           &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for anthrop_heat
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'anthrop_heat',                                                           &
! Variable type
    VAR_TYPE_SURFT,                                                            &
! Long name
    "Anthropogenic heat flux for each tile",                                  &
! Units
    "W m-2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for catch
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'catch',                                                                  &
! Variable type
    VAR_TYPE_SURFT,                                                            &
! Long name
    "Tile surface/canopy water capacity of snow-free land tiles",             &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for ecan
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'ecan',                                                                   &
! Variable type
    VAR_TYPE_SURFT,                                                            &
! Long name
    "Tile evaporation from canopy/surface store for snow-free land tiles",    &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for ei
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'ei',                                                                     &
! Variable type
    VAR_TYPE_SURFT,                                                            &
! Long name
    "Tile sublimation from lying snow for land tiles",                        &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for emis
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'emis',                                                                   &
! Variable type
    VAR_TYPE_SURFT,                                                            &
! Long name
    "Tile emissivity",                                                        &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for esoil
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'esoil',                                                                  &
! Variable type
    VAR_TYPE_SURFT,                                                            &
! Long name
    "Tile surface evapotranspiration from soil moisture store for snow-free land tiles", &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for fqw
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fqw',                                                                    &
! Variable type
    VAR_TYPE_SURFT,                                                            &
! Long name
    "Tile surface moisture flux for land tiles",                              &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for ftl
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'ftl',                                                                    &
! Variable type
    VAR_TYPE_SURFT,                                                            &
! Long name
    "Tile surface sensible heat flux for land tiles",                         &
! Units
    "W m-2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for gc
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'gc',                                                                     &
! Variable type
    VAR_TYPE_SURFT,                                                            &
! Long name
    "Tile surface conductance to evaporation for land tiles",                 &
! Units
    "m s-1"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for le
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'le',                                                                     &
! Variable type
    VAR_TYPE_SURFT,                                                            &
! Long name
    "Tile surface latent heat flux for land tiles",                           &
! Units
    "W m-2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for q1p5m
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'q1p5m',                                                                  &
! Variable type
    VAR_TYPE_SURFT,                                                            &
! Long name
    "Tile specific humidity at 1.5m over land tiles",                         &
! Units
    "kg kg-1"                                                                 &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for rad_net_tile
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'rad_net_tile',                                                           &
! Variable type
    VAR_TYPE_SURFT,                                                            &
! Long name
    "Tile surface net radiation",                                             &
! Units
    "W m-2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_can_melt
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_can_melt',                                                          &
! Variable type
    VAR_TYPE_SURFT,                                                            &
! Long name
    "Tile melt of canopy snow",                                               &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_can
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_can',                                                               &
! Variable type
    VAR_TYPE_SURFT,                                                            &
! Long name
    "Tile canopy snow",                                                       &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_grnd_rho
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_grnd_rho',                                                          &
! Variable type
    VAR_TYPE_SURFT,                                                            &
! Long name
    "Tile density of snow on ground",                                         &
! Units
    "kg m-3"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_ground
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_ground',                                                            &
! Variable type
    VAR_TYPE_SURFT,                                                            &
! Long name
    "Tile snow on ground (snow_tile or snow_grnd)",                           &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_ice_tile
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_ice_tile',                                                          &
! Variable type
    VAR_TYPE_SURFT,                                                            &
! Long name
    "Tile total ice content in snow on ground",                               &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_liq_tile
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_liq_tile',                                                          &
! Variable type
    VAR_TYPE_SURFT,                                                            &
! Long name
    "Tile total liquid content in snow on ground",                            &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_mass
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_mass',                                                              &
! Variable type
    VAR_TYPE_SURFT,                                                            &
! Long name
    "Tile lying snow (snow_tile+snow_grnd)",                                  &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_melt
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_melt',                                                              &
! Variable type
    VAR_TYPE_SURFT,                                                            &
! Long name
    "Tile snow melt (melt_tile)",                                             &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for surf_ht_flux
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'surf_ht_flux',                                                           &
! Variable type
    VAR_TYPE_SURFT,                                                            &
! Long name
    "Downward heat flux for each tile",                                       &
! Units
    "W m-2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_soil_htf
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_soil_htf',                                                          &
! Variable type
    VAR_TYPE_SURFT,                                                            &
! Long name
    "Downward heat flux after snowpack to subsurface",                        &
! Units
    "W m-2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_smb_surft
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_smb_surft',                                                          &
! Variable type
    VAR_TYPE_SURFT,                                                            &
! Long name
    "Rate of change of snowpack mass on tiles",                               &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for surf_ht_store
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'surf_ht_store',                                                          &
! Variable type
    VAR_TYPE_SURFT,                                                            &
! Long name
    "C*(dT/dt) for each tile",                                                &
! Units
    "W m-2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for t1p5m
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    't1p5m',                                                                  &
! Variable type
    VAR_TYPE_SURFT,                                                            &
! Long name
    "Tile temperature at 1.5m over land tiles",                               &
! Units
    "K"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for tstar
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'tstar',                                                                  &
! Variable type
    VAR_TYPE_SURFT,                                                            &
! Long name
    "Tile surface temperature",                                               &
! Units
    "K"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for tsurf_elev_surft
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'tsurf_elev_surft',                                                       &
! Variable type
    VAR_TYPE_SURFT,                                                           &
! Long name
    "Temperature of elevated subsurface tiles",                               &
! Units
    "K"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for z0
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'z0',                                                                     &
! Variable type
    VAR_TYPE_SURFT,                                                            &
! Long name
    "Tile surface roughness",                                                 &
! Units
    "m"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for z0h
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'z0h',                                                                    &
! Variable type
    VAR_TYPE_SURFT,                                                            &
! Long name
    "Tile surface thermal roughness",                                         &
! Units
    "m"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for tile_index
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'tile_index',                                                             &
! Variable type
    VAR_TYPE_TYPE,                                                            &
! Long name
    "Index (gridbox number) of land points with each surface type",           &
! Units
    "m"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for co2_mmr
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'co2_mmr',                                                                &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Concentration of atmospheric CO2, expressed as a mass mixing ratio.",    &
! Units
    ""                                                                        &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for ozone
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'ozone',                                                                  &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for isoprene_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'isoprene_gb',                                                            &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox mean isoprene emission flux",                                    &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for isoprene
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'isoprene',                                                               &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT isoprene emission flux",                                             &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for terpene_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'terpene_gb',                                                             &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox mean (mono-)terpene emission flux",                              &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for terpene
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'terpene',                                                                &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT (mono-)terpene emission flux",                                       &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for methanol_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'methanol_gb',                                                            &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox mean methanol emission flux",                                    &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for methanol
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'methanol',                                                               &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT methanol emission flux",                                             &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for acetone_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'acetone_gb',                                                             &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox mean acetone emission flux",                                     &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for acetone
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'acetone',                                                                &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT acetone emission flux",                                              &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for croprootc
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'croprootc',                                                              &
! Variable type
    VAR_TYPE_CPFT,                                                            &
! Long name
    "Root carbon of crop at the end of model timestep",                       &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for cropharvc
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'cropharvc',                                                              &
! Variable type
    VAR_TYPE_CPFT,                                                            &
! Long name
    "Carbon in harvested parts of crop at the end of model timestep",         &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for cropreservec
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'cropreservec',                                                           &
! Variable type
    VAR_TYPE_CPFT,                                                            &
! Long name
    "Carbon in crop stem reserve pool at the end of model timestep",          &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for cropdvi
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'cropdvi',                                                                &
! Variable type
    VAR_TYPE_CPFT,                                                            &
! Long name
    "Developmental index of crops at the end of model timestep",              &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for cropyield
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'cropyield',                                                              &
! Variable type
    VAR_TYPE_CPFT,                                                            &
! Long name
    "Yield carbon at the end of model timestep",                              &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for harvest_trigger
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'harvest_trigger',                                                        &
! Variable type
    VAR_TYPE_CPFT,                                                            &
! Long name
    "Trigger condition for crop harvest",                                     &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for harvest_counter
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'harvest_counter',                                                        &
! Variable type
    VAR_TYPE_CPFT,                                                            &
! Long name
    "Counter for crop harvest",                                               &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for cropstemc
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'cropstemc',                                                              &
! Variable type
    VAR_TYPE_CPFT,                                                            &
! Long name
    "Carbon in stem parts of crop at the end of model timestep",              &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for cropleafc
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'cropleafc',                                                              &
! Variable type
    VAR_TYPE_CPFT,                                                            &
! Long name
    "Carbon in leaf parts of crop at the end of model timestep",              &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for croplai
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'croplai',                                                                &
! Variable type
    VAR_TYPE_CPFT,                                                            &
! Long name
    "Leaf area index of crop at the end of model timestep",                   &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for cropcanht
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'cropcanht',                                                              &
! Variable type
    VAR_TYPE_CPFT,                                                            &
! Long name
    "Canopy height of crop at the end of model timestep",                     &
! Units
    "m"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for cropsowdate
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'cropsowdate',                                                            &
! Variable type
    VAR_TYPE_CPFT,                                                            &
! Long name
    "Day of year when crop is sown",                                          &
! Units
    "day of year"                                                             &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for cropttveg
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'cropttveg',                                                              &
! Variable type
    VAR_TYPE_CPFT,                                                            &
! Long name
    "Thermal time between emergence and flowering",                           &
! Units
    "degree days"                                                             &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for cropttrep
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'cropttrep',                                                              &
! Variable type
    VAR_TYPE_CPFT,                                                            &
! Long name
    "Thermal time between flowering and maturity/harvest",                    &
! Units
    "degree days"                                                             &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Neutral winds
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! Metadata for u10m_n
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'u10m_n',                                                                 &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox equivalent neutral westerly wind component at 10 m height",      &
! Units
    "m s-1"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for v10m_n
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'v10m_n',                                                                 &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox equivalent neutral southerly wind component at 10 m height",     &
! Units
    "m s-1"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for mu10m_n
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'mu10m_n',                                                                &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox westerly pseudostress referenced to 10 m height",                &
! Units
    "m2 s-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for mv10m_n
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'mv10m_n',                                                                &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox southerly pseudostress referenced to 10 m height",               &
! Units
    "m2 s-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for sthu_irr
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'sthu_irr',                                                               &
! Variable type
    VAR_TYPE_SOIL,                                                            &
! Long name
    "Gridbox wetness of each soil layer over irrigation",                     &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for frac_irrig
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'frac_irrig',                                                             &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Irrigated fraction",                                                     &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for irrDaysDiag
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'irrDaysDiag',                                                            &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Number of days on which irrigation is applied",                          &
! Units
    "days"                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for irrig_water
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'irrig_water',                                                            &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Irrigation water demand",                                                &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for frac_agr
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'frac_agr',                                                               &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox agricultural/crop fraction",                                     &
! Units
    "-"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for frac_agr_prev
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'frac_agr_prev',                                                          &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox agricultural/crop fraction from previous land use update",       &
! Units
    "-"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for frac_past
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'frac_past',                                                              &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox pasture fraction",                                               &
! Units
    "-"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for frac_past_prev
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'frac_past_prev',                                                         &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox pasture fraction from previous land use update",                 &
! Units
    "-"                                                                       &
  ),                                                                          &    
!-----------------------------------------------------------------------------
! Metadata for Wood_prod_fast
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'wood_prod_fast',                                                         &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Wood Products Fast C Pool",                                              &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for Wood_prod_med
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'wood_prod_med',                                                          &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Wood Products Medium C Pool",                                            &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for Wood_prod_slow
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'wood_prod_slow',                                                         &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Wood Products Slow C Pool",                                              &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for WP_fast_in
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'WP_fast_in',                                                             &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Wood Products Fast C Pool Input",                                        &
! Units
    "kg m-2 per 360 days"                                                     &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for WP_med_in
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'WP_med_in',                                                              &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Wood Products Medium C Pool Input",                                      &
! Units
    "kg m-2 per 360 days"                                                     &
  ),                                                                          &

!-----------------------------------------------------------------------------
! Metadata for WP_slow_in
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'WP_slow_in',                                                             &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Wood Products Slow C Pool Input",                                        &
! Units
    "kg m-2 per 360 days"                                                     &
  ),                                                                          &

!-----------------------------------------------------------------------------
! Metadata for WP_fast_out
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'WP_fast_out',                                                            &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Wood Products Fast C Pool Output",                                       &
! Units
    "kg m-2 per 360 days"                                                     &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for WP_med_out
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'WP_med_out',                                                             &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Wood Products Medium C Pool Output",                                     &
! Units
    "kg m-2 per 360 days"                                                     &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for WP_slow_out
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'WP_slow_out',                                                            &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Wood Products Slow C Pool Output",                                       &
! Units
    "kg m-2 per 360 days"                                                     &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for cnsrv_carbon_veg2
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'cnsrv_carbon_veg2',                                                      &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "error in land carbon conservation in veg2 routine",                      &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for cnsrv_carbon_triffid
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'cnsrv_carbon_triffid',                                                   &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "error in land carbon conservation in triffid routine",                   &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for cnsrv_veg_triffid
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'cnsrv_veg_triffid',                                                      &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "error in vegetation carbon conservation in triffid routine",             &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for cnsrv_soil_triffid
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'cnsrv_soil_triffid',                                                     &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "error in soil carbon conservation in triffid routine",                   &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for cnsrv_prod_triffid
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'cnsrv_prod_triffid',                                                     &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "error in wood product carbon conservation in triffid routine",           &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for lit_c_orig
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'lit_c_orig',                                                             &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT carbon Litter including LU",                                         &
! Units
    "kg m-2 per 360days"                                                      &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for lit_n_orig
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'lit_n_orig',                                                             &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT nitrogen Litter including LU",                                       &
! Units
    "kg N  m-2 per 360days"                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for n_inorg
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'n_inorg',                                                                &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Inorganic N Pool",                                                       &
! Units
    "kg m-2"                                                                  &
  ),                                                                          & 
!-----------------------------------------------------------------------------
! Metadata for n_loss
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'n_loss',                                                                 &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Nitrogen loss term (fixed fraction of N_INORG)",                         &
! Units
    "kg m-2 360 days"                                                         &
  ),                                                                          & 
!-----------------------------------------------------------------------------
! Metadata for dpm_ratio
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'dpm_ratio',                                                              &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "DPM:RPM ratio of overall litter input (one value per GB)",               &
! Units
    "fraction (limited to between 0 and 1)"                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for lit_c_ag
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'lit_c_ag',                                                               &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT carbon Litter from LU/Agriculture",                                  &
! Units
    "kg m-2 per 360days"                                                      &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for root_abandon
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'root_abandon',                                                           &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "carbon flux from roots abandoned during landuse change to soil",         &
! Units
    "kg C  m-2 per 360days"                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for root_abandon_n
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'root_abandon_n',                                                           &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "nitrogen flux from roots abandoned during landuse change to soil",         &
! Units
    "kg N  m-2 per 360days"                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for harvest
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'harvest',                                                                &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "flux of carbon to product pools due to harvest",                         &
! Units
    "kg C  m-2 per 360days"                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for harvest_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'harvest_gb',                                                             &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "gridbox mean flux of carbon to product pools due to harvest",            &
! Units
    "kg C  m-2 per 360days"                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for harvest_n
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'harvest_n',                                                              &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "flux of nitrogen to atmosphere due to harvest",                          &
! Units
    "kg N  m-2 per 360days"                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for harvest_n_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'harvest_n_gb',                                                           &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "gridbox mean flux of nitrogen to product pools due to harvest",          &
! Units
    "kg N  m-2 per 360days"                                                   &
  ),                                                                          &  
!-----------------------------------------------------------------------------
! Metadata for n_fertiliser
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'n_fertiliser',                                                           &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "N addition from fertiliser",                                             &
! Units
    "kg N  m-2 per 360days"                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for n_fertiliser_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'n_fertiliser_gb',                                                        &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox N addition from fertiliser",                                     &
! Units
    "kg N  m-2 per 360days"                                                   &
  ),                                                                          &  
 !-----------------------------------------------------------------------------
! Metadata for lit_n_ag
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'lit_n_ag',                                                               &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT nitrogen loss due to LU/Agriculture",                                &
! Units
    "kg N  m-2 per 360days"                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for river routing latitude (2d)
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'latitude_2d',                                                            &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Routing gridbox latitude",                                               &
! Units
    "degrees"                                                                 &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for river routing longitude (2d)
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'longitude_2d',                                                            &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Routing gridbox latitude",                                               &
! Units
    "degrees"                                                                 &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for river direction
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'direction',                                                              &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "River routing direction flag",                                           &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for river sequence
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'sequence',                                                               &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "River routing sequence flag",                                            &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for river drainage area
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'area',                                                                   &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "River routing drainage area",                                            &
! Units
    "10^4 km2"                                                                &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for river routing latitude (1d)
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'rivers_ygrid',                                                           &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Routing gridbox latitude",                                               &
! Units
    "degrees"                                                                 &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for river routing longitude (1d)
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'rivers_xgrid',                                                           &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Routing gridbox longitude",                                              &
! Units
    "degrees"                                                                 &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for river flow
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'rflow',                                                              &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Routing gridbox flow",                                                   &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for river routing grid runoff
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'rrun',                                                                   &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Routing gridbox runoff",                                                 &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for Metstat temp_max_00h_r
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'temp_max_00h_r',                                                         &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Maximum temperature to 00hrs local- running total",                      &
! Units
    "Kelvin"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for Metstat temp_ave_00h_r
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'temp_ave_00h_r',                                                         &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Average temperature to 00hrs local- running total",                      &
! Units
    "Kelvin"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for Metstat prec_tot_00h_r
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'prec_tot_00h_r',                                                         &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Total precipitation to 00hrs local- running total",                      &
! Units
    "kg/ms2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for Metstat prec_tot_12h_r
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'prec_tot_12h_r',                                                         &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Total precipitation to 12hrs local- running total",                      &
! Units
    "kg/m2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for Metstat rhum_min_00h_r
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'rhum_min_00h_r',                                                         &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Minimum relative humidity to 00hrs local- running total",                &
! Units
    "%"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for Metstat dewp_ave_00h_r
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'dewp_ave_00h_r',                                                         &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Average dewpoint to 00 hrs local- running total",                        &
! Units
    "Kelvin"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for Metstat wind_ave_00h_r
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'wind_ave_00h_r',                                                         &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Average wind speed to 00hrs local- running total",                       &
! Units
    "m/s"                                                                     &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for Metstat temp_max_00h
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'temp_max_00h',                                                           &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Maximum temperature to 00hrs local",                                     &
! Units
    "Kelvin"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for Metstat temp_ave_00h
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'temp_ave_00h',                                                           &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Average temperature to 00hrs local",                                     &
! Units
    "Kelvin"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for Metstat temp_pnt_12h
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'temp_pnt_12h',                                                           &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Timestep temperature at last 12hrs local",                               &
! Units
    "Kelvin"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for Metstat prec_tot_00h
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'prec_tot_00h',                                                           &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Total precipitation to 00hrs local",                                     &
! Units
    "kg/ms2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for Metstat prec_tot_12h
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'prec_tot_12h',                                                           &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Total precipitation to 12hrs local",                                     &
! Units
    "kg/m2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for Metstat rhum_min_00h
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'rhum_min_00h',                                                           &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Minimum relative humidity to 00hrs local",                               &
! Units
    "%"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for Metstat rhum_pnt_12h
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'rhum_pnt_12h',                                                           &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Timestep relative humidity at last 12hrs local",                         &
! Units
    "%"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for Metstat dewp_ave_00h
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'dewp_ave_00h',                                                           &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Average dewpoint to 00 hrs local",                                       &
! Units
    "Kelvin"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for Metstat wind_ave_00h
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'wind_ave_00h',                                                           &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Average wind speed to 00hrs local",                                      &
! Units
    "m/s"                                                                     &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for Metstat wind_pnt_12h
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'wind_pnt_12h',                                                           &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Timestep wind speed at last 12hrs local",                                &
! Units
    "m/s"                                                                     &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for fire_mcarthur
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fire_mcarthur',                                                          &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "McArthur FFDI",                                                          &
! Units
    "None"                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for fire_mcarthur_r_dr
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fire_mcarthur_r_dr',                                                     &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "McArthur FFDI- Drought R term",                                          &
! Units
    "None"                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for fire_mcarthur_n_dr
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fire_mcarthur_n_dr',                                                     &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "McArthur FFDI- Drought N term",                                          &
! Units
     "None"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for fire_canadian
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fire_canadian',                                                          &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Canadian FWI",                                                           &
! Units
      "None"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for fire_canadian_ffmc
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fire_canadian_ffmc',                                                     &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Canadian FWI- FFMC",                                                     &
! Units
    "None"                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for fire_canadian_ffmc_mois
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fire_canadian_ffmc_mois',                                                &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Canadian FWI- FFMC Moisture",                                            &
! Units
    "None"                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for fire_canadian_dmc
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fire_canadian_dmc',                                                      &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Canadian FWI- DMC",                                                      &
! Units
    "None"                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for fire_canadian_dc
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fire_canadian_dc',                                                       &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Canadian FWI- DC",                                                       &
! Units
    "None"                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for fire_canadian_isi
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fire_canadian_isi',                                                      &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Canadian FWI- ISI",                                                      &
! Units
    "None"                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for fire_canadian_bui
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fire_canadian_bui',                                                      &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Canadian FWI- BUI",                                                      &
! Units
    "None"                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for fire_nesterov
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fire_nesterov',                                                          &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Nesterov Index",                                                         &
! Units
    "None"                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for leafC
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'leafC',                                                                  &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "Carbon in LEAF biomass on PFTs",                                         &
! Units
    "kg/m2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for rootC
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'rootC',                                                                  &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "Carbon in ROOT biomass on PFTs",                                         &
! Units
    "kg/m2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for stemC
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'stemC',                                                                  &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "Carbon in STEM biomass on PFTs",                                         &
! Units
    "kg/m2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for woodC
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'woodC',                                                                  &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "Carbon in WOODY biomass on PFTs",                                        &
! Units
    "kg/m2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for dleaf
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'dleaf',                                                                  &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "Increment to carbon in LEAF biomass on PFTs",                            &
! Units
    "kg/m2/TRIFFID PERIOD"                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for droot
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'droot',                                                                  &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "Increment to carbon in ROOT biomass on PFTs",                            &
! Units
    "kg/m2/TRIFFID PERIOD"                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for dwood
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'dwood',                                                                  &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "Increment to carbon in WOODY biomass on PFTs",                           &
! Units
    "kg/m2/TRIFFID PERIOD"                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for leaf_litC
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'leaf_litC',                                                              &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "Litter CARBON due to leaf turnover",                                     &
! Units
    "kg/m2/360 days"                                                          &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for root_litC
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'root_litC',                                                              &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "Litter CARBON due to root turnover",                                     &
! Units
    "kg/m2/360 days"                                                          &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for wood_litC
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'wood_litC',                                                              &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "Litter CARBON due to wood turnover",                                     &
! Units
    "kg/m2/360 days"                                                          &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for leaf_litN
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'leaf_litN',                                                              &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "Litter NITROGEN due to leaf turnover",                                   &
! Units
    "kg/m2/360 days"                                                          &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for exudates
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'exudates',                                                               &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "EXUDATES on PFTs: excess C not assimilable into plant due lack of N availability",  &
! Units
    "kg/m2/360 days"                                                          &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for exudates_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'exudates_gb',                                                            &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "EXUDATES (GBM): excess C not assimilable into plant due lack of N availability",  &
! Units
    "kg/m2/360 days"                                                          &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for lai_bal
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'lai_bal',                                                                &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "LAI in balanced growth state (seasonal maximum LAI)",                    &
! Units
    "m2/m2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for litterN
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'litterN',                                                                &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "Local Littter N production",                                             &
! Units
    "kg/,2/360 days"                                                          &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for n_uptake_growth
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'n_uptake_growth',                                                        &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "N TAKEN UP FOR GROWTH of existing plant biomass",                        &
! Units
    "kg/m2/360 days"                                                          &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for n_demand_growth
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'n_demand_growth',                                                        &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "Demand for N for GROWTH of existing plant biomass",                      &
! Units
    "kg/m2/360 days"                                                          &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for n_gas
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'n_gas',                                                                  &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Mineralised N Gas Emissions",                                            &
! Units
    "kg/m2/360 days"                                                          &
  )           ,                                                               &
!-----------------------------------------------------------------------------
! Metadata for n_leach
!-----------------------------------------------------------------------------
var_metadata(                                                                 &
! String identifier
    'n_leach',                                                                &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Leached N Losses",                                                       &
! Units
    "kg/m2/S"                                                                 &
  )           ,                                                               &
!-----------------------------------------------------------------------------
! Metadata for n_fix_gb
!-----------------------------------------------------------------------------
var_metadata(                                                                 &
! String identifier
    'n_fix_gb',                                                               &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Fixed N",                                                                &
! Units
    "kg/m2/360 days"                                                          &
  )           ,                                                               &
!-----------------------------------------------------------------------------
! Metadata for n_fix
!-----------------------------------------------------------------------------
var_metadata(                                                                 &
! String identifier
    'n_fix',                                                                  &
! Variable type
    var_type_pft,                                                             &
! Long name
    "Fixed N on PFTs",                                                        &
! Units
    "kg/m2/360 days"                                                          &
  )           ,                                                               &
!-----------------------------------------------------------------------------
! Metadata for n_demand_spread
!-----------------------------------------------------------------------------
var_metadata(                                                                 &
! String identifier
    'n_demand_spread',                                                        &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "Demand for N for SPREADING plant across gridbox",                        &
! Units
    "kg/m2/360 days"                                                          &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for n_demand_lit
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'n_demand_lit',                                                           &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "N demand of litter: N lost in leaf, wood, root biomass",                 &
! Units
    "kg/m2/360 days"                                                          &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for deposition_n
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'deposition_n',                                                           &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Nitrogen deposition",                                                    &
! Units
    "kg/m2/s"                                                                 &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for root_litN
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'root_litN',                                                              &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "N lost as LITTER due to ROOT turnover",                                  &
! Units
    "kg/m2/360 days"                                                          &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for n_veg
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'n_veg',                                                                  &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "Plant NITROGEN content (PFTs): N_LEAF+N_ROOT+N_WOOD from C equivalents", &
! Units
    "kg/m2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for n_veg_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'n_veg_gb',                                                               &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Plant NITROGEN content (GBM): N_LEAF+N_ROOT+N_WOOD from C equivalents",  &
! Units
    "kg/m2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for dcveg
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'dcveg',                                                                  &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "Increment in vegetation C content (PFTs)",                               &
! Units
    "kg/m2/TRIFFID TIMESTEP"                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for dcveg_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'dcveg_gb',                                                               &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Increment in vegetation C content (GBM)",                                &
! Units
    "kg/m2/TRIFFID TIMESTEP"                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for dnveg
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'dnveg',                                                                  &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "Increment in vegetation N content (PFTs)",                               &
! Units
    "kg/m2/TRIFFID TIMESTEP"                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for dnveg_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'dnveg_gb',                                                               &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Increment in vegetation N content (GBM)",                                &
! Units
    "kg/m2/TRIFFID TIMESTEP"                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for n_demand_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'n_demand_gb',                                                            &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Demand for N in total (GBM)",                                            &
! Units
    "kg/m2/360 DAYS"                                                          &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for n_demand
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'n_demand',                                                               &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "Demand for N in total (PFTs)",                                           &
! Units
    "kg/m2/360 DAYS"                                                          &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for wood_litN
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifierN
    'wood_litN',                                                              &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "N lost as LITTER due to WOOD turnover",                                  &
! Units
    "kg/m2/360 days"                                                          &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for litterC
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'litterC',                                                                &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "LITTER CARBON as the sum of biomass lost through LEAF, ROOT and WOOD turnover",&
! Units
    "kg/m2/360 days"                                                          &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for lit_n
!-----------------------------------------------------------------------------
  var_metadata(                                                              &
! String identifier
    'lit_n',                                                                 &
! Variable type
    VAR_TYPE_PFT,                                                            &
! Long name
    "NITROGEN LITTER (PFTs)",                                                &
! Units
    "kg/m2/360 days"                                                         &
  ),                                                                         &
!-----------------------------------------------------------------------------
! Metadata for lit_n_t
!-----------------------------------------------------------------------------
  var_metadata(                                                              &
! String identifier
    'lit_n_t',                                                               &
! Variable type
    VAR_TYPE_SURFACE,                                                        &
! Long name
    "Total NITROGEN LITTER",                                                 &
! Units
    "kg/m2/360 days"                                                         &
  ),                                                                         &
!-----------------------------------------------------------------------------
! Metadata for n_uptake
!-----------------------------------------------------------------------------
  var_metadata(                                                              &
! String identifier
    'n_uptake',                                                              &
! Variable type
    VAR_TYPE_PFT,                                                            &
! Long name
    "Total NITROGEN UPTAKE by plants (PFTs)",                                &
! Units
    "kg/m2/360 days"                                                         &
  ),                                                                         &
!-----------------------------------------------------------------------------
! Metadata for n_uptake_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                              &
! String identifier
    'n_uptake_gb',                                                           &
! Variable type
    VAR_TYPE_SURFACE,                                                        &
! Long name
    "Total NITROGEN UPTAKE by plants (GBM)",                                 &
! Units
    "kg/m2/360 days"                                                         &
  ),                                                                         &  
!-----------------------------------------------------------------------------
! Metadata for burnt_area_gb (INFERNO)
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'burnt_area_gb',                                                          &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox mean burnt area fraction",                                       &
! Units
    "s-1"                                                                     &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for burnt_area
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'burnt_area',                                                             &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT burnt area fraction",                                                &
! Units
    "s-1"                                                                     &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for emitted_carbon_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'emitted_carbon_gb',                                                        &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox mean emitted carbon",                                              &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for emitted_carbon
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'emitted_carbon',                                                           &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT emitted carbon",                                                       &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for emitted_carbon_DPM
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'emitted_carbon_gb',                                                        &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "DPM emitted carbon",                                                       &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for emitted_carbon_RPM
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'emitted_carbon_RPM',                                                       &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "RPM emitted carbon",                                                       &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for fire_em_CO2_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fire_em_CO2_gb',                                                         &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox mean CO2 emission from fires",                                   &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for fire_em_CO2
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fire_em_CO2',                                                            &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT CO2 emission from fires",                                            &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          & 
!-----------------------------------------------------------------------------
! Metadata for fire_em_CO2_DPM
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fire_em_CO2_DPM',                                                         &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "DPM CO2 emission from fires",                                            &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for fire_em_CO2_RPM
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fire_em_CO2_RPM',                                                        &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "RPM CO2 emission from fires",                                            &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for fire_em_CO_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fire_em_CO_gb',                                                          &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox mean CO emission from fires",                                    &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &       
!-----------------------------------------------------------------------------
! Metadata for fire_em_CO
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fire_em_CO',                                                             &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT CO emission from fires",                                             &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for fire_em_CO_DPM
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fire_em_CO_DPM',                                                         &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "DPM CO emission from fires",                                             &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &       
!-----------------------------------------------------------------------------
! Metadata for fire_em_CO_RPM
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fire_em_CO_RPM',                                                         &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "RPM CO emission from fires",                                             &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &       
!-----------------------------------------------------------------------------
! Metadata for fire_em_CH4_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fire_em_CH4_gb',                                                         &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox mean CH4 emission from fires",                                   &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &   
!-----------------------------------------------------------------------------
! Metadata for fire_em_CH4
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fire_em_CH4',                                                            &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "Gridbox mean CH4 emission from fires",                                   &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &   
!-----------------------------------------------------------------------------
! Metadata for fire_em_CH4_DPM
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fire_em_CH4_DPM',                                                        &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "DPM CH4 emission from fires",                                            &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &   
!-----------------------------------------------------------------------------
! Metadata for fire_em_CH4_RPM
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fire_em_CH4_RPM',                                                        &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "RPM CH4 emission from fires",                                            &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &   
!-----------------------------------------------------------------------------
! Metadata for fire_em_NOx_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fire_em_NOx_gb',                                                         &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox mean NOx emission from fires",                                   &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &   
!-----------------------------------------------------------------------------
! Metadata for fire_em_NOx
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fire_em_NOx',                                                            &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "Gridbox mean NOx emission from fires",                                   &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &   
!-----------------------------------------------------------------------------
! Metadata for fire_em_NOx_DPM
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fire_em_NOx_DPM',                                                        &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "DPM NOx emission from fires",                                            &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &   
!-----------------------------------------------------------------------------
! Metadata for fire_em_NOx_RPM
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fire_em_NOx_RPM',                                                        &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "RPM NOx emission from fires",                                            &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &   
!-----------------------------------------------------------------------------
! Metadata for fire_em_SO2_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fire_em_SO2_gb',                                                         &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox mean SO2 emission from fires",                                   &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &   
!-----------------------------------------------------------------------------
! Metadata for fire_em_SO2
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fire_em_SO2',                                                            &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "Gridbox mean SO2 emission from fires",                                   &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &   
!-----------------------------------------------------------------------------
! Metadata for fire_em_SO2_DPM
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fire_em_SO2_DPM',                                                        &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "DPM SO2 emission from fires",                                            &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &   
!-----------------------------------------------------------------------------
! Metadata for fire_em_SO2_RPM
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fire_em_SO2_RPM',                                                        &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "RPM SO2 emission from fires",                                            &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &   
!-----------------------------------------------------------------------------
! Metadata for fire_em_OC_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fire_em_OC_gb',                                                          &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox mean OC emission from fires",                                    &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &       
!-----------------------------------------------------------------------------
! Metadata for fire_em_OC
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fire_em_OC',                                                             &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT OC emission from fires",                                             &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for fire_em_OC_DPM
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fire_em_OC_DPM',                                                         &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "DPM OC emission from fires",                                             &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &       
!-----------------------------------------------------------------------------
! Metadata for fire_em_OC_RPM
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fire_em_OC_RPM',                                                         &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "RPM OC emission from fires",                                             &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &       
!-----------------------------------------------------------------------------
! Metadata for fire_em_BC_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fire_em_BC_gb',                                                          &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox mean BC emission from fires",                                    &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &       
!-----------------------------------------------------------------------------
! Metadata for fire_em_BC
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fire_em_BC',                                                             &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT BC emission from fires",                                             &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for fire_em_BC_DPM
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fire_em_BC_DPM',                                                         &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "DPM BC emission from fires",                                             &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for fire_em_BC_RPM
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fire_em_BC_RPM',                                                         &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "RPM BC emission from fires",                                             &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for pop_den
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'pop_den',                                                                &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Population Density",                                                     &
! Units
    "ppl/km2"                                                                 &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for flash_rate
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'flash_rate',                                                             &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Cloud to Ground lightning flash rate",                                   &
! Units
    "flashes/km2"                                                             &
  )                                                                           &
/
#endif
