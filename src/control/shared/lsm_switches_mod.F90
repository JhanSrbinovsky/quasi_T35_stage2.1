! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE lsm_switches_mod

!-----------------------------------------------------------------------------
! Description:
!   CABLE_LSM:
!   Declare default and read namelist overide for LSM switch between 
!  JULES & CABLE
!
! Current Code Owner: Jhan Srbinovsky
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!
! NB: Checks have not been implemented at this time
!-----------------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Switches
!-----------------------------------------------------------------------------
  INTEGER, parameter :: n_chars_lsm_id = 5
  character(len=n_chars_lsm_id):: lsm_id='jules' ! default case

!-----------------------------------------------------------------------------
! Single namelist definition for UM and standalone
!-----------------------------------------------------------------------------
  NAMELIST /lsm_switches/ lsm_id


  CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='LSM_SWITCHES_MOD'

CONTAINS


#if defined(UM_JULES)
  SUBROUTINE read_nml_lsm_switches(unitnumber)

! Description:
!  Read the JULES_SOIL namelist

    USE setup_namelist, ONLY : setup_nml_type
    USE check_iostat_mod, ONLY: check_iostat
    USE UM_parcore,       ONLY: mype

    USE parkind1, ONLY: jprb, jpim
    USE yomhook, ONLY: lhook, dr_hook

    USE errormessagelength_mod, ONLY: errormessagelength

    IMPLICIT NONE

! Subroutine arguments
    INTEGER, INTENT(IN) :: unitnumber

    INTEGER :: my_comm
    INTEGER :: mpl_nml_type
    INTEGER :: ErrorStatus
    INTEGER :: icode
    CHARACTER(LEN=errormessagelength) :: iomessage
    REAL(KIND=jprb)               :: zhook_handle

    CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_LSM_SWITCHES'
    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

    ! set number of each type of variable in my_namelist type
    INTEGER, PARAMETER :: no_of_types = 1
    INTEGER, PARAMETER :: n_chars =n_chars_lsm_id 

    TYPE my_namelist
      SEQUENCE
      character(len=n_chars_lsm_id):: lsm_id   
    END TYPE my_namelist

    TYPE (my_namelist) :: my_nml

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    CALL gc_get_communicator(my_comm, icode)

    CALL setup_nml_type( no_of_types, mpl_nml_type, n_chars_in=n_chars)
    
    IF (mype == 0) THEN
      !rewind unitnumber

      READ (unitnumber, NML=lsm_switches, IOSTAT=errorstatus,              &
            IOMSG=iomessage)
      
      CALL check_iostat(errorstatus, "namelist LSM switches", iomessage)
      
      my_nml % lsm_id          = lsm_id    

    END IF

    CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

    IF (mype /= 0) THEN
      lsm_id          = my_nml % lsm_id   
    END IF

    CALL mpl_type_free(mpl_nml_type,icode)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END SUBROUTINE read_nml_lsm_switches
#endif

END MODULE lsm_switches_mod
