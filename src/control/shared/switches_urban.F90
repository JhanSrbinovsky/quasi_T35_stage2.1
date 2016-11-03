MODULE switches_urban

! Description:
!   Module containing switches for the parametrisations of urban scheme MORUSES
!
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! -----------------------------------------------------------------------------

  IMPLICIT NONE

  LOGICAL ::                                                                  &
! l_urban2T and l_moruses, when initially set, are mutually exclusive as they
! are used to set the other switches to their default values in init_urban.
!
#if !defined(UM_JULES)
     l_urban_empirical      = .TRUE.,  & ! Empirical relationships for urban
                                         ! geometry (WRR & HWR)
#endif
     l_moruses_macdonald    = .TRUE.,  & ! MacDonald formulation for
                                         ! displacement height and effective
                                         ! roughness length for momentum
     l_urban2T              = .FALSE., & ! Original URBAN-2T switch
     l_moruses              = .FALSE., & ! Indicates any MORUSES switch used

! Independent parametristaion switches
     l_moruses_albedo       = .TRUE.,  & ! SW canyon albedo
     l_moruses_emissivity   = .FALSE., & ! LW canyon emissivity
     l_moruses_rough        = .TRUE.,  & ! Heat transfer
     l_moruses_storage      = .TRUE.,  & ! Storage
     l_moruses_storage_thin = .TRUE.     ! Storage thin roof

!-----------------------------------------------------------------------
! Set up a namelist to allow switches to be set
! In the UM, l_urban2T is in the namelist, l_moruses is set by inspecting
!   other MORUSES switches (i.e. if l_moruses_albedo = T then l_moruses = T
!
! In standalone JULES, l_moruses is in the namelist, and l_urban2T is
!   set by looking at whether 2 urban tiles have been specified
!
! Either way, after initialisation, we have l_moruses = T => l_urban2T = T,
!   but not necessarily the other way round
!-----------------------------------------------------------------------
  NAMELIST /urban_switches/                                                   &
#if defined(UM_JULES)
     l_urban2T,                                                               &
#else
     l_moruses,l_urban_empirical,                                             &
#endif
     l_moruses_albedo,l_moruses_emissivity,l_moruses_rough,l_moruses_storage, &
     l_moruses_storage_thin,l_moruses_macdonald

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SWITCHES_URBAN'

CONTAINS
  SUBROUTINE print_nlist_urban_switches()
    USE jules_print_mgr, ONLY : jules_print
    IMPLICIT NONE
    CHARACTER(LEN=50000) :: lineBuffer

    CALL jules_print('switches_urban',                                        &
        'Contents of namelist urban_switches')

    WRITE(lineBuffer,*)' l_urban2T = ',l_urban2T
    CALL jules_print('switches_urban',lineBuffer)
    WRITE(lineBuffer,*)' l_moruses_albedo = ',l_moruses_albedo
    CALL jules_print('switches_urban',lineBuffer)
    WRITE(lineBuffer,*)' l_moruses_emissivity = ',l_moruses_emissivity
    CALL jules_print('switches_urban',lineBuffer)
    WRITE(lineBuffer,*)' l_moruses_rough = ',l_moruses_rough
    CALL jules_print('switches_urban',lineBuffer)
    WRITE(lineBuffer,*)' l_moruses_storage = ',l_moruses_storage
    CALL jules_print('switches_urban',lineBuffer)
    WRITE(lineBuffer,*)' l_moruses_storage_thin = ',l_moruses_storage_thin
    CALL jules_print('switches_urban',lineBuffer)
    WRITE(lineBuffer,*)' l_moruses_macdonald = ',l_moruses_macdonald
    CALL jules_print('switches_urban',lineBuffer)

    CALL jules_print('switches_urban',                                        &
        '- - - - - - end of namelist - - - - - -')

  END SUBROUTINE print_nlist_urban_switches

#if defined(UM_JULES)

  SUBROUTINE read_nml_urban_switches (unitnumber)

! Description:
!  Read the URBAN_SWITCHES namelist

    USE setup_namelist,   ONLY : setup_nml_type
    USE check_iostat_mod, ONLY : check_iostat
    USE UM_parcore,       ONLY : mype

    USE parkind1,         ONLY : jprb, jpim
    USE yomhook,          ONLY : lhook, dr_hook

    USE errormessagelength_mod, ONLY: errormessagelength

    IMPLICIT NONE

! Subroutine arguments
    INTEGER, INTENT(IN) :: unitnumber

    INTEGER :: my_comm
    INTEGER :: mpl_nml_type
    INTEGER :: ErrorStatus
    INTEGER :: icode
    CHARACTER(LEN=errormessagelength) :: iomessage
    REAL(KIND=jprb) :: zhook_handle

    CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_URBAN_SWITCHES'
    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

    ! set number of each type of variable in my_namelist type
    INTEGER, PARAMETER :: no_of_types = 1
    INTEGER, PARAMETER :: n_log = 7

    TYPE my_namelist
      SEQUENCE
      LOGICAL :: l_urban2T
      LOGICAL :: l_moruses_albedo
      LOGICAL :: l_moruses_emissivity
      LOGICAL :: l_moruses_rough
      LOGICAL :: l_moruses_storage
      LOGICAL :: l_moruses_storage_thin
      LOGICAL :: l_moruses_macdonald
    END TYPE my_namelist

    TYPE (my_namelist) :: my_nml

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    CALL gc_get_communicator(my_comm, icode)

    CALL setup_nml_type(no_of_types, mpl_nml_type, n_log_in=n_log)

    IF (mype == 0) THEN

      READ (UNIT=unitnumber, NML=urban_switches, IOSTAT=errorstatus,          &
            IOMSG=iomessage)
      CALL check_iostat(errorstatus, "namelist urban_switches", iomessage)

      my_nml % l_urban2T              = l_urban2T
      my_nml % l_moruses_albedo       = l_moruses_albedo
      my_nml % l_moruses_emissivity   = l_moruses_emissivity
      my_nml % l_moruses_rough        = l_moruses_rough
      my_nml % l_moruses_storage      = l_moruses_storage
      my_nml % l_moruses_storage_thin = l_moruses_storage_thin
      my_nml % l_moruses_macdonald    = l_moruses_macdonald
    END IF

    CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

    IF (mype /= 0) THEN

      l_urban2T              = my_nml % l_urban2T
      l_moruses_albedo       = my_nml % l_moruses_albedo
      l_moruses_emissivity   = my_nml % l_moruses_emissivity
      l_moruses_rough        = my_nml % l_moruses_rough
      l_moruses_storage      = my_nml % l_moruses_storage
      l_moruses_storage_thin = my_nml % l_moruses_storage_thin
      l_moruses_macdonald    = my_nml % l_moruses_macdonald

    END IF

    CALL mpl_type_free(mpl_nml_type,icode)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END SUBROUTINE read_nml_urban_switches
#endif

END MODULE switches_urban
