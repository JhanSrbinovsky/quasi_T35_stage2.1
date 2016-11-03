! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Subroutine UPDATE_VEG
!
!  Purpose: The routine is entered when any of the ancillary
!           fields have to be updated. It then checks to see if
!           leaf area index and/or canopy height have been updated.
!           If this is the case, then the subroutine SPARM is called
!           to ensure that all other vegetation parameters are
!           consistent.
!
!           Code Owner: Please refer to the UM file CodeOwners.txt
!           This file belongs in section: jules

MODULE update_veg_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE ::                                       &
  ModuleName='UPDATE_VEG_MOD'

CONTAINS

SUBROUTINE update_veg ()

!Use in relevant subroutines

!USE in JULES modules
USE jules_surface_mod, ONLY: l_aggregate
USE jules_vegetation_mod, ONLY: can_model
USE jules_surface_types_mod

!USE in UM modules (check!)
USE atm_fields_bounds_mod
USE atm_fields_mod,           ONLY: frac_typ,                                 &
                                    canht_pft, lai_pft, sat_soil_cond,        &
                                    z0m_soil, catch_snow, catch_tile,         &
                                    infil_tile, z0_tile, z0h_tile

!Dispose of this lot after we ditch D1?
USE nlsizes_namelist_mod, ONLY:                                               &
    land_field, ntiles
USE UM_ParParams
USE um_stashcode_mod, ONLY:                                                   &
    stashcode_lai, stashcode_canopy_height, stashcode_z0m_soil
USE ancil_file_mod, ONLY: n_anc_upd
USE cancila_mod, ONLY: fieldcode, update

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

!     Local variables
INTEGER :: tile_pts(ntype)              ! No of land pts which include surface
                                        !  type ntype
INTEGER :: tile_index(land_field,ntype) ! Indices of land points which
                                        !  include the surface type ntype

INTEGER :: field
INTEGER :: field_lai, field_canht, field_z0msoil
LOGICAL :: update_lai, update_canht, update_z0msoil

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UPDATE_VEG'

! Update vegetation parameters if required

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

field_lai   = 0
update_lai = .FALSE.
field_canht = 0
update_canht = .FALSE.
field_z0msoil = 0
update_z0msoil = .FALSE.
! check to see which fields are present as ancillary files to do updates,
! and whether or not they are updated at this time.
DO field =  1,n_anc_upd
  IF (fieldcode(3,field) == stashcode_lai) THEN
    ! field is present:
    field_lai = field
    ! is it updated?:
    update_lai = update(field_lai)
  END IF
  IF (fieldcode(3,field) == stashcode_canopy_height) THEN
    field_canht = field
    update_canht = update(field_canht)
  END IF
  IF (fieldcode(3,field) == stashcode_z0m_soil) THEN
    field_z0msoil = field
    update_z0msoil = update(field_z0msoil)
  END IF
END DO

IF ((field_lai > 0 .AND. field_canht > 0) .OR. field_z0msoil > 0)  THEN

  IF (update_lai .OR. update_canht .OR. update_z0msoil ) THEN
    !-----------------------------------------------------------------------
    ! Call TILEPTS to initialise TILE_PTS and TILE_INDEX
    !-----------------------------------------------------------------------

    ! DEPENDS ON: tilepts
    CALL tilepts(land_field,frac_typ,tile_pts,tile_index)
    !-----------------------------------------------------------------------
    ! Initialise tiled and gridbox mean vegetation parameters
    !-----------------------------------------------------------------------

    ! DEPENDS ON: sparm
      CALL sparm (land_field, ntiles, can_model, l_aggregate, tile_pts,       &
                  tile_index, frac_typ, canht_pft, lai_pft, sat_soil_cond,    &
                  z0m_soil, catch_snow, catch_tile, infil_tile, z0_tile,      &
                  z0h_tile)
  END IF

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE update_veg
END MODULE update_veg_mod
