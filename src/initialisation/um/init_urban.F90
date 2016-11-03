#if defined(UM_JULES)
! *****************************COPYRIGHT*******************************
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE init_urban_mod
  CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='INIT_URBAN_MOD'

CONTAINS
SUBROUTINE init_urban (land_pts, frac, hgt_p, hwr_p, wrr_p, disp_p, ztm_p,    &
                       albwl_p, albrd_p, emisw_p, emisr_p)
! Description:
!   Routine to initialize URBAN parameters.
!   This is intended for testing only as it initialises the urban arrays using
!     set values
!
!!!  6.1   10/01/07   First written. Peter Clark and Aurore Porson
!
! Code Owner: See Unified Model Code Owner's HTML page
! This file belongs in section: Land

  USE urban_param,              ONLY:                                         &
     hgt_gb, hwr_gb, wrr_gb, disp_gb, ztm_gb, albwl_gb, albrd_gb, emisw_gb,   &
     emisr_gb,a, cdz, kappa2, z0m_mat

  USE switches_urban, ONLY :                                                  &
     l_urban2T, l_moruses, l_moruses_albedo, l_moruses_emissivity,            &
     l_moruses_rough, l_moruses_storage, l_moruses_storage_thin,              &
     l_moruses_macdonald

  USE jules_surface_types_mod, ONLY :                                         &
     urban, ice, urban_canyon, urban_roof, ntype

  USE parkind1, ONLY: jprb, jpim
  USE yomhook, ONLY: lhook, dr_hook
  USE jules_print_mgr, ONLY :                                                 &
    jules_message,                                                            &
    jules_print,                                                              &
    PrNorm
  IMPLICIT NONE

! Arguments:
  INTEGER, INTENT(IN) ::                                                      &
     land_pts                      ! Number of land points to be processed.
                                   ! include the nth surface type.

  REAL, INTENT(INOUT) ::                                                      &
     frac(land_pts,ntype),       & ! Fractional cover of each surface type.
     hgt_p(land_pts),            & ! From d1: Building height
     hwr_p(land_pts),            & ! From d1: Height to width
     wrr_p(land_pts),            & ! From d1: Width ratio
     albwl_p(land_pts),          & ! From d1: Wall albedo
     albrd_p(land_pts),          & ! From d1: Road albedo
     emisw_p(land_pts),          & ! From d1: Wall emmissivity
     emisr_p(land_pts),          & ! From d1: Road emmissivity
     disp_p(land_pts),           & ! From d1: Displacemnet height
     ztm_p(land_pts)               ! From d1: Roughness length

! Local declarations:
  REAL ::                                                                     &
     sc_hwr(land_pts),       & ! working variable
     d_h(land_pts),          & ! working variable
     lambdaf, lambdap          ! Frontal and planar area index

  INTEGER ::                                                                  &
     l                         ! WORK Loop counters

  REAL :: urban_fraction        ! Temporary store

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='INIT_URBAN'


!----------------------------------------------------------------------
! Set parameters for urban morphology
!----------------------------------------------------------------------

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialising arrays
  hgt_gb(:)   = 0.0
  hwr_gb(:)   = 0.0
  wrr_gb(:)   = 0.0
  albwl_gb(:) = 0.0
  albrd_gb(:) = 0.0
  emisw_gb(:) = 0.0
  emisr_gb(:) = 0.0
  ztm_gb(:)   = 0.0
  disp_gb(:)  = 0.0

! Check logic

! Initialise MORUSES switch
  l_moruses = .FALSE.

  IF ( l_moruses_albedo .OR. l_moruses_emissivity                             &
     .OR. l_moruses_rough .OR. l_moruses_storage ) THEN
    ! Turn on l_moruses if any of the independent switches used
    l_moruses = .TRUE. ! This does not mean that all moruses switches are true
                       ! This is used to set z0 to ztm in sparm.
    l_urban2T = .TRUE. ! MORUSES must be used with URBAN-2T
  END IF

  CALL jules_print('init_urban','init_urban',level=PrNorm)
  CALL jules_print('init_urban','',level=PrNorm)
  CALL jules_print('init_urban','Urban switches used',level=PrNorm)
  WRITE(jules_message,*) 'l_urban2T             ', l_urban2T
  CALL jules_print('init_urban',jules_message,level=PrNorm)
  WRITE(jules_message,*) 'l_moruses             ', l_moruses
  CALL jules_print('init_urban',jules_message,level=PrNorm)
  WRITE(jules_message,*) 'l_moruses_albedo      ', l_moruses_albedo
  CALL jules_print('init_urban',jules_message,level=PrNorm)
  WRITE(jules_message,*) 'l_moruses_emissivity  ', l_moruses_emissivity
  CALL jules_print('init_urban',jules_message,level=PrNorm)
  WRITE(jules_message,*) 'l_moruses_rough       ', l_moruses_rough
  CALL jules_print('init_urban',jules_message,level=PrNorm)
  WRITE(jules_message,*) 'l_moruses_storage     ', l_moruses_storage
  CALL jules_print('init_urban',jules_message,level=PrNorm)
  WRITE(jules_message,*) 'l_moruses_storage_thin', l_moruses_storage_thin
  CALL jules_print('init_urban',jules_message,level=PrNorm)
  WRITE(jules_message,*) 'l_moruses_macdonald   ', l_moruses_macdonald
  CALL jules_print('init_urban',jules_message,level=PrNorm)
  CALL jules_print('init_urban','',level=PrNorm)

  ! Issue warnings on logic
  IF ( l_moruses_storage_thin .AND. .NOT. l_moruses_storage ) THEN
    CALL jules_print('init_urban','',level=PrNorm)
    CALL jules_print('init_urban',                                            &
        'WARNING: MORUSES storage parametrisation not used.',level=PrNorm)
    CALL jules_print('init_urban',                                            &
        'l_moruses_storage      = .FALSE. when',level=PrNorm)
    CALL jules_print('init_urban',                                            &
        'l_moruses_storage_thin = .TRUE.',level=PrNorm)
  END IF
  IF ( l_moruses_macdonald .AND. .NOT. l_moruses ) THEN
    CALL jules_print('init_urban','',level=PrNorm)
    CALL jules_print('init_urban',                                            &
        'WARNING: MORUSES is not switched on.',level=PrNorm)
    CALL jules_print('init_urban',                                            &
        'l_moruses              = .FALSE. when',level=PrNorm)
    CALL jules_print('init_urban',                                            &
        'l_moruses_macdonald    = .TRUE.',level=PrNorm)
  END IF

  IF ( l_urban2T ) THEN

    CALL jules_print('init_urban',                                            &
        "Setting URBAN-2T parameters",level=PrNorm)

! Fill allocatable arrays so that these can be passed in a module.
    wrr_gb(:)= wrr_p(:)

    IF ( l_moruses ) THEN
      CALL jules_print('init_urban',"Setting MORUSES parameters",level=PrNorm)

      hwr_gb(:)   = hwr_p(:)
      hgt_gb(:)   = hgt_p(:)
      albwl_gb(:) = albwl_p(:)
      albrd_gb(:) = albrd_p(:)
      emisw_gb(:) = emisw_p(:)
      emisr_gb(:) = emisr_p(:)

      IF ( l_moruses_macdonald ) THEN
        !       Macdonald Formulation
        CALL jules_print('init_urban','Using MacDonald formulation',          &
            level=PrNorm)
        sc_hwr(:) = 0.5 * ( hwr_gb(:) / (2.0 * ATAN(1.0)) )
        d_h(:)    = 1.0 - wrr_gb(:) * ( a**(wrr_gb(:) - 1.0) )
        DO l = 1, land_pts
!       Urban present (set to 1.0 otherwise) > 0 also as a check
          IF ( wrr_gb(l) > 0.0 .AND. wrr_gb(l) < 1.0 ) THEN
            disp_gb(l)   = d_h(l) * hgt_gb(l)
            ztm_gb(l)    = (cdz * (1.0 - d_h(l)) *                            &
                            sc_hwr(l) * wrr_gb(l) / kappa2)**(-0.5)
            ztm_gb(l)    = (1.0 - d_h(l))*EXP(-ztm_gb(l))
            ztm_gb(l)    = ztm_gb(l) * hgt_gb(l)
            ztm_gb(l)    = MAX(ztm_gb(l),z0m_mat)
          END IF
        END DO
        ztm_p(:)  = ztm_gb(:)
        disp_p(:) = disp_gb(:)
      ELSE
        CALL jules_print('init_urban','',level=PrNorm)
        CALL jules_print('init_urban',                                        &
            'WARNING: Ancillary file used instead of MacDonald formulation',  &
            level=PrNorm)
        ztm_gb(:)  = ztm_p(:)
        disp_gb(:) = disp_p(:)
      END IF

    END IF ! l_moruses

  ELSE    ! .NOT. l_urban2T

    CALL jules_print('init_urban','URBAN-2T OR MORUSES not used',level=PrNorm)
    CALL jules_print('init_urban',                                            &
        'All associated parameters initialised to zero',level=PrNorm)

  END IF   ! l_urban2T

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END SUBROUTINE init_urban
END MODULE init_urban_mod
#endif
