module cable_gather_UM_data_decs
   implicit none

   REAL, DIMENSION(:,:), allocatable :: &
    ls_rain_cable,    &
    ls_snow_cable,    &
    conv_rain_cable,    &
    conv_snow_cable

   REAL, DIMENSION(:,:,:), allocatable :: &
    surf_down_sw_cable
      
Contains


subroutine cable_set_ls_precip( row_length, rows,                         &
                ls_rain, ls_snow )

  integer :: row_length,rows
  
  real, dimension(row_length,rows) :: &
    ls_rain,    &
    ls_snow

    if(.NOT. allocated(ls_rain_cable) ) &
      allocate( ls_rain_cable(row_length,rows),    &
                ls_snow_cable(row_length,rows) )
    
    ls_rain_cable = ls_rain
    ls_snow_cable = ls_snow

End subroutine cable_set_ls_precip

subroutine cable_set_conv_precip( row_length, rows,                         &
                conv_rain, conv_snow )

  integer :: row_length,rows
  
  real, dimension(row_length,rows) :: &
    conv_rain,    &
    conv_snow
    
    if(.NOT. allocated(conv_rain_cable) ) &
      allocate(conv_rain_cable(row_length,rows),    &
                conv_snow_cable(row_length,rows) )
    
    conv_rain_cable = conv_rain
    conv_snow_cable = conv_snow

End subroutine cable_set_conv_precip

subroutine cable_set_shortwave( row_length, rows, nbands,                     &
                surf_down_sw )

  integer :: row_length,rows, nbands
  
  real, dimension(row_length,rows,nbands) :: &
    surf_down_sw

    if(.NOT. allocated(conv_rain_cable) ) &
      allocate(surf_down_sw_cable(row_length,rows, nbands) )
    
    surf_down_sw_cable = surf_down_sw

End subroutine cable_set_shortwave

End module cable_gather_UM_data_decs



