! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module setting Tile ID numbers. The UM uses these to identify which land
!   surface tiles are present

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.

MODULE land_tile_ids

  USE max_dimensions,   ONLY : &
      ntype_max,               &
      snow_layers_max,         &
      elev_ice_tile_max

  USE missing_data_mod, ONLY : imdi

  IMPLICIT NONE

!-----------------------------------------------------------------------

! Land surface types :
!       1 - Broadleaf Tree
!     101 - Broadleaf Tree deciduous
!     102 - Broadleaf Tree evergreen tropical
!     103 - Broadleaf Tree evergreen temperate
!       2 - Needleleaf Tree
!     201 - Needleleaf Tree deciduous
!     202 - Needleleaf Tree evergreen
!       3 - C3 Grass
!     301 - C3 Crop
!     302 - C3 Pasture
!       4 - C4 Grass
!     401 - C4 Crop
!     402 - C4 Pasture
!       5 - Shrub
!     501 - Shrub deciduous
!     502 - Shrub evergreen
!       6 - Urban
!     601 - Urban roof
!     602 - Urban canyon
!       7 - Water
!       8 - Soil
!       9 - Ice
! 901-925 - Elev_Ice(1-25)


!-----------------------------------------------------------------------

  INTEGER, PRIVATE :: i ! Loop counter

! Surface type definitions
  INTEGER, PARAMETER ::                                                       &
      id_brd_leaf          =  1                                               &
                  ! ID number of surface type 'Broadleaf Tree'
     ,id_brd_leaf_dec      =  101                                             &
                  ! ID number of surface type 'Broadleaf Tree deciduous'
     ,id_brd_leaf_eg_trop  =  102                                             &
                  ! ID number of surface type 'Broadleaf Tree eg tropical'
     ,id_brd_leaf_eg_temp  =  103                                             &
                  ! ID number of surface type 'Broadleaf Tree eg temperate'
     ,id_ndl_leaf          =  2                                               &
                  ! ID number of surface type 'Needleleaf Tree'
     ,id_ndl_leaf_dec      =  201                                             &
                  ! ID number of surface type 'Needleleaf Tree deciduous'
     ,id_ndl_leaf_eg       =  202                                             &
                  ! ID number of surface type 'Needleleaf Tree evergreen'
     ,id_C3grass           =  3                                               &
                  ! ID number of surface type 'C3 Grass'
     ,id_C3crop            =  301                                             &
                  ! ID number of surface type 'C3 Crop'
     ,id_C3pasture         =  302                                             &
                  ! ID number of surface type 'C3 Pasture'
     ,id_C4grass           =  4                                               &
                  ! ID number of surface type 'C4 Grass'
     ,id_C4crop            =  401                                             &
                  ! ID number of surface type 'C4 Crop'
     ,id_C4pasture         =  402                                             &
                  ! ID number of surface type 'C4 Pasture'
     ,id_shrub             =  5                                               &
                  ! ID number of surface type 'Shrub'
     ,id_shrub_dec         =  501                                             &
                  ! ID number of surface type 'Shrub deciduous'
     ,id_shrub_eg          =  502                                             &
                  ! ID number of surface type 'Shrub evergreen'
     ,id_urban             =  6                                               &
                  ! ID number of surface type 'Urban' (l_urban2T == .FALSE.)
     ,id_lake              =  7                                               &
                  ! ID number of surface type 'Lake'
     ,id_soil              =  8                                               &
                  ! ID number of surface type 'Soil'
     ,id_ice               =  9                                               &
                  ! ID number of surface type 'Ice'
     ,id_elev_ice(25)      =  (/ (i, i = 901, 925) /)                         &
                  ! ID numbers of surface type 'Elev_Ice'

! For use when l_urban2T == .TRUE. (Urban 2-tile scheme on)

     ,id_urban_canyon      = 601                                              &
                  ! ID number of surface type 'Urban roof'
     ,id_urban_roof        = 602                                              &
                  ! ID number of surface type 'Urban canyon'
     ,id_usr_type(90)      = (/ (i, i = 10, 99) /)
                  ! ID numbers for user defined surface type (10-99)

  INTEGER :: surface_type_ids(ntype_max) = 0
                  ! Array which maps pseudo levels to tile types

  INTEGER :: ml_snow_type_ids(ntype_max*snow_layers_max) = 0
                  ! Array which maps pseudo levels to tile types

  INTEGER :: tile_ids_in(ntype_max*snow_layers_max) = imdi
                  ! Tile IDs in the input header

CONTAINS

  SUBROUTINE set_tile_id_arrays ( )

    USE ereport_mod,             ONLY : ereport

    USE jules_snow_mod,          ONLY : nsmax

    USE jules_surface_types_mod, ONLY : &
        ntype,                          &
        brd_leaf,                       &
        brd_leaf_dec,                   &
        brd_leaf_eg_trop,               &
        brd_leaf_eg_temp,               &
        ndl_leaf,                       &
        ndl_leaf_dec,                   &
        ndl_leaf_eg,                    &
        c3_grass,                       &
        c3_crop,                        &
        c3_pasture,                     &
        c4_grass,                       &
        c4_crop,                        &
        c4_pasture,                     &
        shrub,                          &
        shrub_dec,                      &
        shrub_eg,                       &
        urban,                          &
        urban_canyon,                   &
        urban_roof,                     &
        lake,                           &
        soil,                           &
        ice,                            &
        usr_type,                       &
        elev_ice

    USE um_parvars,              ONLY : mype

    USE umPrintMgr,              ONLY : &
        umPrint,                        &
        umMessage,                      &
        printstatus,                    &
        prstatus_diag

    IMPLICIT NONE

    INTEGER      :: i, j ! Loop counter
    INTEGER      :: pseudo, errorstatus
    CHARACTER (LEN=20) :: umFormat
    CHARACTER (LEN=80) :: cmessage
    CHARACTER (LEN=18), PARAMETER :: routinename='set_tile_id_arrays'

! Set up surface type (pl_code = 9) ids array

    IF (brd_leaf         > 0)                                                 &
                    surface_type_ids(brd_leaf         ) = id_brd_leaf
    IF (brd_leaf_dec     > 0)                                                 &
                    surface_type_ids(brd_leaf_dec     ) = id_brd_leaf_dec
    IF (brd_leaf_eg_trop > 0)                                                 &
                    surface_type_ids(brd_leaf_eg_trop ) = id_brd_leaf_eg_trop
    IF (brd_leaf_eg_temp > 0)                                                 &
                    surface_type_ids(brd_leaf_eg_temp ) = id_brd_leaf_eg_temp
    IF (ndl_leaf         > 0)                                                 &
                    surface_type_ids(ndl_leaf         ) = id_ndl_leaf
    IF (ndl_leaf_dec     > 0)                                                 &
                    surface_type_ids(ndl_leaf_dec     ) = id_ndl_leaf_dec
    IF (ndl_leaf_eg      > 0)                                                 &
                    surface_type_ids(ndl_leaf_eg      ) = id_ndl_leaf_eg
    IF (c3_grass         > 0)                                                 &
                    surface_type_ids(c3_grass         ) = id_c3grass
    IF (c3_crop          > 0)                                                 &
                    surface_type_ids(c3_crop          ) = id_c3crop
    IF (c3_pasture       > 0)                                                 &
                    surface_type_ids(c3_pasture       ) = id_c3pasture
    IF (c4_grass         > 0)                                                 &
                    surface_type_ids(c4_grass         ) = id_c4grass
    IF (c4_crop          > 0)                                                 &
                    surface_type_ids(c4_crop          ) = id_c4crop
    IF (c4_pasture       > 0)                                                 &
                    surface_type_ids(c4_pasture       ) = id_c4pasture
    IF (shrub            > 0)                                                 &
                    surface_type_ids(shrub            ) = id_shrub
    IF (shrub_dec        > 0)                                                 &
                    surface_type_ids(shrub_dec        ) = id_shrub_dec
    IF (shrub_eg         > 0)                                                 &
                    surface_type_ids(shrub_eg         ) = id_shrub_eg
    IF (urban            > 0)                                                 &
                    surface_type_ids(urban            ) = id_urban
    IF (lake             > 0)                                                 &
                    surface_type_ids(lake             ) = id_lake
    IF (soil             > 0)                                                 &
                    surface_type_ids(soil             ) = id_soil
    IF (ice              > 0)                                                 &
                    surface_type_ids(ice              ) = id_ice
    IF (urban_canyon     > 0)                                                 &
                    surface_type_ids(urban_canyon     ) = id_urban_canyon
    IF (urban_roof       > 0)                                                 &
                    surface_type_ids(urban_roof       ) = id_urban_roof

    DO i = 1, ntype
      IF (usr_type(i)    > 0)                                                 &
                    surface_type_ids(usr_type(i)      ) = id_usr_type(i)
    END DO

    IF ( SIZE(id_elev_ice) /= elev_ice_tile_max ) THEN
      errorstatus = 6973
      WRITE(cmessage, '(A)') 'SIZE(id_elev_ice) must match elev_ice_tile_max'
      CALL ereport ( routinename, errorstatus, cmessage )
    END IF

    DO i = 1, elev_ice_tile_max
      IF (elev_ice(i)    > 0)                                                 &
                    surface_type_ids(elev_ice(i)      ) = id_elev_ice(i)
    END DO

! Set up multilayer snow type (pl_code = 11) ids array
    IF ( nsmax > 0 ) THEN
      DO i = 1, ntype
        DO j = 1, nsmax
          pseudo = i + (j-1) * ntype
          ml_snow_type_ids(pseudo) = surface_type_ids(i)*1000+j
        END DO
      END DO
    END IF

    IF (mype == 0 .AND. printstatus >= prstatus_diag ) THEN
      WRITE(umFormat,'(a3,i3,a8)') '(a,',ntype,'(1x,i6))'
      WRITE(umMessage,umFormat) 'Surface types present =',                    &
         surface_type_ids(1:ntype)
      CALL umPrint(umMessage,src='set_tile_id_arrays')
    END IF

! Check that tile IDs are unique
    DO i = 1, ntype
      IF ( COUNT( surface_type_ids(:) == surface_type_ids(i) ) /= 1 ) THEN
        errorstatus = 6173
        WRITE(cmessage, '(A,I6)')                                             &
           'Surface type ID not unique :', surface_type_ids(i)
        CALL ereport ( routinename, errorstatus, cmessage )
      END IF
    END DO

  END SUBROUTINE set_tile_id_arrays

END MODULE land_tile_ids
