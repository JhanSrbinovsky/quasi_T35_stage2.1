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

module cable_hyd_main_mod
  
contains

SUBROUTINE cable_hyd_main( ) 

  !USE UM modules 
  USE timestep_mod,   ONLY : timestep_number  ! number
  
  !--- USE CABLE modules 
  USE cable_decs_mod, ONLY : CABLE_file
  USE cable_logs_mod, ONLY :                                                   &
    first_CABLE_call,     & ! LOGICAL var: CABLE wide first call  
    KBL_progress,         & ! TYPEd details of progress log 
    cable_logs,           & ! SUBR to set CABLE wide stuff
    init_nml_prefs        & ! SUBR to init. nml prefs
  
  USE cable_progs_mod, ONLY : cable !TYPEd CABLE progs, set in UM IO

  USE cable_common_module, ONLY : knode_gl,        & ! processor number
                                  ktau_gl,         & ! number
                                  kwidth_gl          ! width in S 
  
  USE cable_write_logs_mod, ONLY : write_progress_log

  implicit none
 
  !--- IN ARGS FROM surf_couple_hydiation() ------------------------------------
  
  !--- End IN ARGS FROM surf_couple_hydiation() --------------------------------


  !--- declare local vars ------------------------------------------------------ 

  character(len=*), parameter :: subr_name = "cable_hyd_main"
  logical, save :: first_call = .true.
  
  !--- End header -------------------------------------------------------------
  
  ktau_gl   = timestep_number
  
  !----------------------------------------------------------------------------
  !--- Organize report writing for CABLE.                         -------------
  !--- Progress log and IN args @ timestep X,Y,Z                  -------------
  !----------------------------------------------------------------------------
  
  ! get unique "KBL_progress" wherever this is called from 
  if( first_CABLE_call ) call cable_logs()
  
  !jhan: check this bc not always 0
  ! IF we are interested in logging CALLs to different parts of CABLE 
  if( knode_gl==0 ) then
    if( KBL_progress%want ) &
      call write_progress_log( KBL_progress, ktau_gl, subr_name, cycleno )
  endif

  !jhan: call checks as required by namelis      

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

End subroutine cable_hyd_main
  
End module cable_hyd_main_mod











































