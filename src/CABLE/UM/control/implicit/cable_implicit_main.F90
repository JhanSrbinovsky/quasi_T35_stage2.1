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
! Called from: JULES: surf_couple_ pathway
!
! Contact: Jhan.Srbinovsky@csiro.au
!
! History: Developed for CABLE-JULES coupling in UM 10.5
!
!
! ==============================================================================

module cable_implicit_main_mod
  
contains

SUBROUTINE cable_implicit_main(                                                &        
              mype, timestep_number,cycleno, numcycles,                        & 
              row_length,rows, land_pts, ntiles, npft, sm_levels,              &
              ls_rain, ls_snow,                                    &
              conv_rain, conv_snow!, &
              !canopy_gb, smcl &
            )  
  
  implicit none
 
  !--- IN ARGS FROM sf_impl2_cable, passed from surf_couple_implicit() down ----
  integer :: mype, timestep_number, cycleno, numcycles
  integer :: row_length,rows, land_pts, ntiles, npft, sm_levels
!jjhan
!_rain & _ snow still have  to be defined. 1st decs in atm_step(ls_) &
!atmos_phys2(conv_) 
  REAL, DIMENSION(row_length,rows) :: &
    ls_rain,    &
    ls_snow,    &
    conv_rain,    &
    conv_snow

REAL ::                                                                       &
 canopy_gb(land_pts), &
 smcl(land_pts, sm_levels) 
  !--- End IN ARGS  -----------------------------------------------------------

  !--- declare local vars ------------------------------------------------------ 

  character(len=*), parameter :: subr_name = "cable_implicit_main"
  logical, save :: first_call = .true.
  
  !--- End header -------------------------------------------------------------
  
  if(mype==0) then
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

End subroutine cable_implicit_main
  
End module cable_implicit_main_mod











































