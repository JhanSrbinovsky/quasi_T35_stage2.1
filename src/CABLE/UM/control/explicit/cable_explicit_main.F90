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

module cable_explicit_main_mod
  
contains

SUBROUTINE cable_explicit_main( mype, timestep_number ) 
  
  implicit none
 
  !--- IN ARGS FROM sf_exch_cable, passed from surf_couple_explicit() down ----
  integer :: mype, timestep_number
  !--- End IN ARGS  -----------------------------------------------------------

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











































