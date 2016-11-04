!==============================================================================
! This source code is part of the 
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CABLE Academic User Licence Agreement 
! (the "Licence").
! You may not use this file except in compliance with the Licence.
! A copy of the Licence and registration form can be obtained from 
! http://www.cawcr.gov.au/projects/access/cable
! You need to register and read the Licence agreement before use.
! Please contact cable_help@nf.nci.org.au for any questions on 
! registration and the Licence.
!
! Unless required by applicable law or agreed to in writing, 
! software distributed under the Licence is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the Licence for the specific language governing permissions and 
! limitations under the Licence.
! ==============================================================================
!
! Purpose:
!
! Called from: JULES: surf_couple_
!
! Contact: Jhan.Srbinovsky@csiro.au
!
! History: Developed for CABLE-JULES coupling in UM 10.5
!
!
! ==============================================================================

module cable_explicit_main_mod
  
contains

SUBROUTINE cable_explicit_main(                                                &
            mype, timestep_number, cycleno, numcycles,                         &
            ! grid, model, dimensions. PFT frac per landpoint    
            row_length, rows, land_pts, ntiles,                                &
            npft, sm_levels,                                                   &
            latitude, longitude,                                               &
            land_index, tile_frac, tile_pts, tile_index,                       &
            ! Soil parameters **jhan: could be undef.@some point  issue here
            bexp, hcon, satcon, sathh,                                         &
            smvcst, smvcwt, smvccl, albsoil,                                   & 
            ! packs/ unpacked ssnow% snowd 
            snow_tile, lw_down, cosine_zenith_angle,                           &
            surf_down_sw, ls_rain, ls_snow,                                    &          
            tl_1, qw_1, vshr_land, pstar, z1_tq, z1_uv, canopy_tile,           &
            Fland, CO2_MMR, sthu, canht_ft, lai_ft ,                           &
            sin_theta_latitude, dzsoil, FTL_TILE, FQW_TILE,                    &
            TSTAR_TILE, U_S, U_S_STD_TILE, CD_TILE, CH_TILE,                   &
            RADNET_TILE, FRACA, rESFS, RESFT, Z0H_TILE, Z0M_TILE,              &
            RECIP_L_MO_TILE, EPOT_TILE ) 

  implicit none
 
  !--- IN ARGS FROM sf_exch_cable, passed from surf_couple_explicit() down ----
   INTEGER ::                                                      & 
     mype, timestep_number, cycleno, numcycles
   INTEGER ::                                                      & 
      row_length, rows, & ! UM grid resolution
      land_pts,         & ! # of land points being processed
      ntiles,           & ! # of tiles 
      npft,             & ! # of plant functional types
      sm_levels          ! # of soil layers 

   REAL,  DIMENSION(row_length,rows) ::                             &
      latitude,   &
      longitude
   
   ! index of land points being processed
   INTEGER, DIMENSION(land_pts) :: land_index 

   REAL, DIMENSION(land_pts, ntiles) ::                            &
      tile_frac
   
   ! # of land points on each tile
   INTEGER,  DIMENSION(ntiles) :: tile_pts 

   INTEGER,  DIMENSION(land_pts, ntiles) ::                         & 
      tile_index   ! index of tile points being processed

   !___UM soil/snow/radiation/met vars
   REAL,  DIMENSION(land_pts) :: & 
      bexp,    & ! => parameter b in Campbell equation 
      hcon,    & ! Soil thermal conductivity (W/m/K).
      satcon,  & ! hydraulic conductivity @ saturation [mm/s]
      sathh,   &
      smvcst,  &
      smvcwt,  &
      smvccl,  &
      albsoil

   REAL,  DIMENSION(land_pts, ntiles) ::                         &
      snow_tile
   
   REAL,  DIMENSION(row_length,rows) ::                             &
      lw_down,    &
      cosine_zenith_angle
 
   REAL, DIMENSION(row_length, rows, 4) ::                         &
      surf_down_sw 
   
   REAL,  DIMENSION(row_length,rows) ::                             &
      ls_rain,    &
      ls_snow

   REAL,  DIMENSION(row_length,rows) ::                             &
      tl_1,       &
      qw_1,       &  
      vshr_land,  &
      pstar,      &
      z1_tq,      &
      z1_uv

   REAL, DIMENSION(land_pts, ntiles) ::                             &
      canopy_tile
   
   REAL,  DIMENSION(land_pts) :: & 
      fland
   
   REAL :: co2_mmr

   REAL, DIMENSION(land_pts, sm_levels) ::                         &
      sthu 
   
   REAL, DIMENSION(land_pts, npft) ::                              &
      canht_ft, lai_ft 
   
   REAL :: sin_theta_latitude(row_length,rows) 
     
   REAL,  DIMENSION(sm_levels) :: dzsoil

   REAL, DIMENSION(land_pts,ntiles) :: &
      FTL_TILE,   &  ! Surface FTL for land tiles     
      FQW_TILE       ! Surface FQW for land tiles     
   
   !___return temp and roughness, friction velocities/drags/ etc

   REAL, DIMENSION(land_pts,ntiles) :: &
      TSTAR_TILE
   
   REAL, DIMENSION(row_length,rows)  :: &
      U_S               ! Surface friction velocity (m/s)

   REAL, DIMENSION(land_pts,ntiles) :: &
      U_S_STD_TILE,     & ! Surface friction velocity
      CD_TILE,    &     ! Drag coefficient
      CH_TILE           ! Transfer coefficient for heat & moisture

   !___return miscelaneous 
   REAL,  DIMENSION(land_pts,ntiles) :: &
      RADNET_TILE,    & ! Surface net radiation
      FRACA,          & ! Fraction of surface moisture
      RESFS,          & ! Combined soil, stomatal & aerodynamic resistance
                        ! factor for fraction (1-FRACA) of snow-free land tiles
      RESFT,          & ! Total resistance factor.
                        ! FRACA+(1-FRACA)*RESFS for snow-free l_tile_pts,        
                        ! 1 for snow.    
      Z0H_TILE,       &
      Z0M_TILE,       &
      RECIP_L_MO_TILE,& ! Reciprocal of the Monin-Obukhov length for tiles (m^-1)
      EPOT_TILE

!REAL,  DIMENSION(land_pts, ntiles) ::                            & 
!snow_flg3l   ! 3 layer snow flag
!
!REAL, DIMENSION(land_pts, ntiles) ::                            &
!   snow_rho1l, &
!   snow_age
!
!REAL,  DIMENSION(land_pts, ntiles,3) ::                       &
!   snow_cond
!
!REAL, DIMENSION(land_pts, ntiles,3) ::                          &
!   snow_rho3l,    &
!   snow_depth3l,  &
!   snow_mass3l,   &
!   snow_tmp3l
!
!REAL, DIMENSION(land_pts, ntiles, sm_levels) :: & 
!   sthu_tile, &
!   sthf_tile, &
!   smcl_tile, &
!   tsoil_tile
 
!real :: SOIL_ORDER_casa(:)
!real :: LAI_casa(:,:)
!real :: PHENPHASE_casa(:,:)
!real :: C_pool_casa(:,:,:)
!real :: N_pool_casa(:,:,:)
!real :: N_dep_casa(:)
!real :: N_fix_casa(:)
!real :: P_pool_casa(:,:,:)
!real :: P_dust_casa(:)
!real :: P_weath_casa(:)

  
! rml 2/7/13 Extra atmospheric co2 variables
!LOGICAL, INTENT(IN) :: L_CO2_INTERACTIVE
!INTEGER, INTENT(IN) ::                              &
!  CO2_DIM_LEN                                      &
!  ,CO2_DIM_ROW
!REAL, INTENT(IN) :: CO2_3D(CO2_DIM_LEN,CO2_DIM_ROW)  ! co2 mass mixing ratio

  !--- End IN ARGS ----------------------------------------------------------

   !___UM parameters: water density
   REAL, parameter :: rho_water = 1000.

   !___true IF vegetation (tile) fraction is greater than 0
   LOGICAL, DIMENSION(land_pts, ntiles) :: L_tile_pts
  
!   !r825 adds CASA vars here
!   REAL, DIMENSION(land_pts,ntiles,10) :: &
!      CPOOL_TILE,    & ! Carbon Pools
!      NPOOL_TILE       ! Nitrogen Pools
!
!   REAL, DIMENSION(land_pts,ntiles,12) :: &
!      PPOOL_TILE       ! Phosphorus Pools
!
!   REAL, DIMENSION(land_pts) :: &
!      SOIL_ORDER,    & ! Soil Order (1 to 12)
!      NIDEP,         & ! Nitrogen Deposition
!      NIFIX,         & ! Nitrogen Fixation
!      PWEA,          & ! Phosphorus from Weathering
!      PDUST            ! Phosphorus from Dust
!
!   REAL, DIMENSION(land_pts,ntiles) :: &
!      GLAI, &          ! Leaf Area Index for Prognostics LAI
!      PHENPHASE        ! Phenology Phase for Casa-CNP
!                                  
!   REAL, DIMENSION(land_pts,ntiles) :: &
!      NPP_FT_ACC,     &
!      RESP_W_FT_ACC
      


  !--- declare local vars ------------------------------------------------------ 

  character(len=*), parameter :: subr_name = "cable_explicit_main"
  logical, save :: first_call = .true.
  
  !--- End header -------------------------------------------------------------
  
  if(mype==0) then
    !For now this is hardwired in explicit only but can be ubiquitous
    if(first_call)  & 
      write (6, *) "CABLE_LSM: You have reached CABLE. Placeholders in the &
     &UM/JULES have been extended. surf_couple_*.F90 subroutines switch between &
     &JULES/CABLE based on the character value of lsm_id. By default this is &
     &set to *jules*, however the namelist lsm_switches can be configured in &
     &app/um/rose-app.conf and re-defined as *cable*. In which case you will end &
     &up here. The FCM configuratiion has been modified to build the entire &
     &CABLE model, however only these dummy routines currently exist, echoing this &
     &message on the first call and returing to run JULES."
    
    write (6, *) "CABLE_LSM:Subr: ", subr_name,  "@ timestep: ", timestep_number 
  
  endif
     
  !----------------------------------------------------------------------------
  !--- Organize report writing for CABLE.                         -------------
  !--- Progress log and IN args @ timestep X,Y,Z                  -------------
  !----------------------------------------------------------------------------
  
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  
  !----------------------------------------------------------------------------
  !--- CALL _driver to run specific and necessary components of CABLE with IN -
  !--- args PACKED to force CABLE
  !----------------------------------------------------------------------------
  
  
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  
  !----------------------------------------------------------------------------
  !--- CALL _driver to run specific and necessary components of CABLE with IN -
  !--- args PACKED to force CABLE
  !----------------------------------------------------------------------------
  
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  
  !----------------------------------------------------------------------------
  !--- Organize report writing for CABLE.                         -------------
  !--- OUT args @ timestep X,Y,Z                                  -------------
  !----------------------------------------------------------------------------

  !jhan: call checks as required by namelis      
  
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  
  first_call = .false.        

return

End subroutine cable_explicit_main
  
End module cable_explicit_main_mod











































