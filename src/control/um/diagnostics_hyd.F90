#if defined(UM_JULES)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine diagnostics_hyd
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: land surface

      SUBROUTINE diagnostics_hyd(                                             &
                             row_length, rows, model_levels                   &
      ,                      n_rows, global_row_length, global_rows           &
      ,                      halo_i, halo_j, off_x, off_y, me                 &
      ,                      n_proc, n_procx, n_procy                         &
      ,                      g_rows, g_row_length                             &
      ,                      at_extremity                                     &
      ,                      land_points, dsm_levels                          &
!  Add inland basin outflow to arguments
      ,                      land_index,inlandout_atm                         &
      ,                      smc, surf_roff, sub_surf_roff                    &
      ,                      snow_depth_land, snow_melt, canopy, t_soil       &
      ,                      tsurf_elev_surft, snow_soil_htf, snow_smb_surft  &
      ,                      soil_layer_moisture                              &
      ,                      nsurft, snomlt_surf_htf, sthu, sthf              &
      ,                      tot_tfall, snow_surft, melt_surft                &
      ,                      rgrain, land_sea_mask                            &
      ,                      dun_roff, drain, qbase, qbase_zw                 &
      ,                      fch4_wetl,fch4_wetl_cs,fch4_wetl_npp             &
      ,                      fch4_wetl_resps                                  &
      ,                      fexp, gamtot, ti_mean, ti_sig                    &
      ,                      fsat, fwetl, zw, sthzw                           &
      ,                      timestep, stashwork                        )

! Description:
!   Calculates hydrology-related diagnostics (held in STASH section 8).
!
! Method:
!   Required level lists and logical switches are determined by the
!   calling routine from STASH requests and STASHflags.
!   Intercepted arrays and diagnostic arrays are input from the
!   hydrology routines called previously. Each diagnostic is simply
!   copied into the stashwork array to be passed on to STASH for
!   output processing.
!
!   Diagnostics currently available (in order calculated):
!   Item  Description
!    208  Soil moisture content
!     23  Snow depth
!    201  Snow melt
!    234  Surface run-off
!    235  Sub-surface run-off
!    225  Soil temperature
!    223  Soil layer moisture
!    204  Surface run-off (accumulated)
!    205  Sub-surface run-off (accumulated)
!    209  Canopy water content
!    202  Land snow melt heat flux
!    231  Land snow melt rate
!    233  Canopy throughfall rate
!    236  Snow amount on tiles
!    237  Snow melt rate on tiles
!    238  Snow grain size on tiles
!    229  Unfrozen soil moisture fraction
!    230  Frozen soil moisture fraction
!    239  Baseflow
!    240  Dunne Runoff
!    241  Baseflow from water table layer
!    242  Wetland methane flux
!    252  Drainage out of "nshyd"th model layer
!    245  Inland basin outflow on atmos grid
!    252  Drainage out of bottom "nshyd"th soil layer (currently layer 4).

      USE mask_compression, ONLY:  expand_from_mask
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE ereport_mod, ONLY : ereport

      USE submodel_mod, ONLY: internal_model_index, atmos_im
      USE stash_array_mod, ONLY:                                              &
          sf, si, stlist, stindex, len_stlist, stash_pseudo_levels,           &
          num_stash_pseudo, stash_levels, num_stash_levels

      USE missing_data_mod, ONLY: rmdi

      USE errormessagelength_mod, ONLY: errormessagelength

      IMPLICIT NONE

      LOGICAL                                                                 &
        at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

! Arguments with Intent IN. ie: Input variables.

      INTEGER                                                                 &
        row_length                                                            &
                         ! number of points on a row
      , rows                                                                  &
                         ! number of rows in a theta field
      , n_rows                                                                &
                         ! number of rows in a v field
      , model_levels                                                          &
                         ! number of model levels
      , number_format
                         ! switch controlling number format diagnostics
                         ! are written out in. See PP_WRITE for details.


      INTEGER                                                                 &
        global_row_length                                                     &
                            !IN. NUMBER OF points on a global row
      , global_rows                                                           &
                            !IN. NUMBER OF global rows
      , me                                                                    &
                            !IN. Processor number
      , halo_i                                                                &
                            !IN. size of large halo in x direction
      , halo_j                                                                &
                            !IN. size of large halo in y direction
      , off_x                                                                 &
                            !IN. size of small halo in x direction
      , off_y                                                                 &
                            !IN. size of small halo in y direction
      , n_proc                                                                &
      , n_procx                                                               &
      , n_procy                                                               &
      , g_rows(0:n_proc-1)                                                    &
      , g_row_length(0:n_proc-1)                                              &
      , land_points                                                           &
                    ! No.of land points being processed, can be 0.
      , dsm_levels                                                            &
      , nsurft      ! No. of land-surface tiles ( MOSES II )

      REAL                                                                    &
        lat_rot_NP                                                            &
      , long_rot_NP


      REAL                                                                    &
        timestep

! Primary Arrays used in all models
      INTEGER                                                                 &
        land_index(land_points)      ! set from land_sea_mask

      REAL                                                                    &
        snow_depth_land(land_points)                                          &
                                     !
      , snow_melt(land_points)                                                &
                                ! snowmelt (kg/m2/s)
      , canopy (land_points)                                                  &
                              ! canopy water content (kg/m2)
      , smc(land_points)                                                      &
                            ! available soil moisture in the
!                                 rootzone (kg/m2).
      , surf_roff(land_points)                                                &
                                ! surface runoff (kg/m2/s).
      , sub_surf_roff(land_points)                                            &
                                    ! sub-surface runoff
! Declare inland basin outflow variable
      , inlandout_atm(land_points)                                            &
                                    !inland basin outflow

      , t_soil(land_points,dsm_levels)                                        &
      , tsurf_elev_surft(land_points,nsurft)                                  &
      , snow_soil_htf(land_points,nsurft)                                     &
      , snow_smb_surft(land_points,nsurft)                                    &
      , soil_layer_moisture(land_points,dsm_levels)                           &
      , snomlt_surf_htf(row_length, rows)                                     &
      , sthu(land_points,dsm_levels)                                          &
                                       ! Unfrozen soil moisture
!                content of each layer as a fraction of saturation
      , sthf(land_points,dsm_levels)                                          &
                                       ! Frozen soil moisture content
!                of each layer as a fraction of saturation
      , tot_tfall(land_points)                                                &
                                       ! total throughfall (kg/m2/s)
      , snow_surft(land_points,nsurft)                                        &
                                       ! Lying snow on tiles (kg/m2)
      , melt_surft(land_points,nsurft)                                        &
                                       ! Snowmelt on tiles (kg/m2/s)
      , rgrain(land_points,nsurft)                                            &
                                       ! Snow grain size (microns)

! Additional variables required for large-scale hydrology:
      , qbase(land_points)                                                    &
                                    ! Baseflow (kg/m2/s).
      , dun_roff(land_points)                                                 &
                                    ! Dunne runoff (kg/m2/s).
      , qbase_zw(land_points)                                                 &
                                    ! Baseflow from Zw layer (kg/m2/s).
      , drain(land_points)                                                    &
                                    ! Drainage out of "nshyd" later
!                                   ! (kg/m2/s).
      , fch4_wetl(land_points)                                                &
                                    ! Wetland methane flux
!                                   ! (default substrate) for atmos
!                                   ! chemistry (10^-9 kg C/m2/s).
      ,fch4_wetl_cs(land_points)                                              &
                                    ! Wetland methane flux (soil carbon
!                                   ! substrate) (kg C/m2/s)
      ,fch4_wetl_npp(land_points)                                             &
                                    ! Wetland methane flux (npp
!                                   ! substrate) (kg C/m2/s)
      ,fch4_wetl_resps(land_points)                                           &
                                    ! Wetland methane flux (soil respiration
!                                   ! substrate) (kg C/m2/s)
      , TI_MEAN(LAND_POINTS)                                                  &
                                    ! Mean topographic index.
      , TI_SIG(LAND_POINTS)                                                   &
                                    ! Std. dev. of topographic index.
      , FEXP(LAND_POINTS)                                                     &
                                    ! Decay factor in Sat. Conductivity
!                                   !   in water table layer.
      , GAMTOT(LAND_POINTS)                                                   &
                                    ! Integrated complete Gamma
!                                   !   function.
      , FSAT(LAND_POINTS)                                                     &
                                    ! Surface saturation fraction.
      , FWETL(LAND_POINTS)                                                    &
                                    ! Wetland fraction.
      , ZW(LAND_POINTS)                                                       &
                                    ! Water table depth (m).
      , STHZW(LAND_POINTS)          ! Saturation fraction in water
!                                   !   table layer.


      LOGICAL                                                                 &
        land_sea_mask(row_length, rows)

      INTEGER                                                                 &
       PSLEVEL                                                                &
                     !  loop counter for pseudolevels
      ,PSLEVEL_OUT                                                            &
                     !  index for pseudolevels sent to STASH
      ,LEVEL                                                                  &
                     !  loop counter for levels
      ,LEVEL_OUT     !  index for levels sent to STASH

      LOGICAL                                                                 &
       PLLTILE(nsurft)                                                        &
                          ! pseudolevel list for surface types
      ,LIST(DSM_LEVELS)   ! level list for soil levels


! Local variables

      INTEGER                                                                 &
        i, j, k, l                                                            &
      ,    icode                ! Return code  =0 Normal exit  >1 Error

      CHARACTER(LEN=errormessagelength)  cmessage
      CHARACTER(LEN=*) RoutineName
      PARAMETER ( RoutineName='DIAGNOSTICS_HYD')

      INTEGER                                                                 &
        im_index        ! internal model index

      REAL                                                                    &
        interp_data(row_length,rows)                                          &
      , interp_data_3(row_length,rows,dsm_levels)

! Diagnostic variables
       REAL                                                                   &
        stashwork(*)    ! STASH workspace

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! ----------------------------------------------------------------------
! Section 1.  Initialisation.
! ----------------------------------------------------------------------


      IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
      icode = 0 ! Initialise error status
      im_index = internal_model_index(atmos_im)

! ----------------------------------------------------------------------
! Soil moisture content
! ----------------------------------------------------------------------
! Item 8 208  smc

      IF (sf(208,8)) THEN

         DO j = 1, rows
            DO i = 1, row_length
               interp_data(i,j) = rmdi
            END DO
         END DO

         DO l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = smc(l)
         END DO

! DEPENDS ON: copydiag
         CALL copydiag (stashwork(si(208,8,im_index)),interp_data,            &
              row_length,rows,0,0,0,0, at_extremity,                          &
              atmos_im,8,208,                                                 &
              icode,cmessage)

         IF (icode  >   0) THEN
            cmessage="Error in copydiag( item 208)"
            GOTO 9999
         END IF

      END IF


! ----------------------------------------------------------------------
! Snow depth
! ----------------------------------------------------------------------
! Item 8 023 snow_depth_land

      IF (sf(023,8)) THEN

         DO j = 1, rows
            DO i = 1, row_length
               interp_data(i,j) = rmdi
            END DO
         END DO

         DO l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = snow_depth_land(l)
         END DO

! DEPENDS ON: copydiag
         CALL copydiag(stashwork(si(023,8,im_index)),interp_data,             &
              row_length,rows,0,0,0,0, at_extremity,                          &
              atmos_im,8,023,                                                 &
              icode,cmessage)

         IF (icode  >   0) THEN
            cmessage="Error in copydiag( item 023)"
            GOTO 9999
         END IF

      END IF

! ----------------------------------------------------------------------
! Snow melt
! ----------------------------------------------------------------------
! Item 201

      IF (sf(201,8)) THEN

         DO j = 1, rows
            DO i = 1, row_length
               interp_data(i,j) = rmdi
            END DO
         END DO

         DO l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = snow_melt(l) * timestep
         END DO

! DEPENDS ON: copydiag
         CALL copydiag(stashwork(si(201,8,im_index)),interp_data,             &
              row_length,rows,0,0,0,0, at_extremity,                          &
              atmos_im,8,201,                                                 &
              icode,cmessage)

         IF (icode  >   0) THEN
            cmessage="Error in copydiag( item 201)"
            GOTO 9999
         END IF

      END IF


! ----------------------------------------------------------------------
! Surface run-off.
! ----------------------------------------------------------------------
! Item 234  surf_roff

      IF (sf(234,8)) THEN

         DO j = 1, rows
            DO i = 1, row_length
               interp_data(i,j) = rmdi
            END DO
         END DO

         DO l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = surf_roff(l)
         END DO

! DEPENDS ON: copydiag
         CALL copydiag(stashwork(si(234,8,im_index)),interp_data,             &
              row_length,rows,0,0,0,0, at_extremity,                          &
              atmos_im,8,234,                                                 &
              icode,cmessage)

         IF (icode  >   0) THEN
            cmessage="Error in copydiag( item 234)"
            GOTO 9999
         END IF

      END IF

! ----------------------------------------------------------------------
! Sub-surface run-off.
! ----------------------------------------------------------------------
! Item 235  sub_surf_roff

      IF (sf(235,8)) THEN

         DO j = 1, rows
            DO i = 1, row_length
               interp_data(i,j) = rmdi
            END DO
         END DO

         DO l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = sub_surf_roff(l)
         END DO

! DEPENDS ON: copydiag
         CALL copydiag(stashwork(si(235,8,im_index)),interp_data,             &
              row_length,rows,0,0,0,0, at_extremity,                          &
              atmos_im,8,235,                                                 &
              icode,cmessage)

         IF (icode  >   0) THEN
            cmessage="Error in copydiag( item 235)"
            GOTO 9999
         END IF

      END IF


! ----------------------------------------------------------------------
!  Soil temperature
! ----------------------------------------------------------------------
! Item 8 225 t_soil

      IF (sf(225,8)) THEN

         DO k = 1, dsm_levels
            DO j= 1, rows
               DO i = 1, row_length
                  interp_data_3(i,j,k) = rmdi
               END DO
            END DO

            DO l = 1, land_points
               j=(land_index(l)-1)/row_length + 1
               i=land_index(l) - (j-1)*row_length
               interp_data_3(i,j,k) = t_soil(l,k)
            END DO

         END DO

! DEPENDS ON: copydiag_3d
         CALL copydiag_3d(stashwork(si(225,8,im_index)),interp_data_3,        &
              row_length,rows,dsm_levels,0,0,0,0,                             &
              at_extremity,                                                   &
              stlist(1,stindex(1,225,8,im_index)),len_stlist,                 &
              stash_levels,num_stash_levels+1,                                &
              atmos_im,8,225,                                                 &
              icode,cmessage)



         IF (icode  >   0) THEN
            cmessage="Error in copydiag_3d( item 225)"
            GOTO 9999
         END IF

      END IF

! ----------------------------------------------------------------------
! Temperature of elevated sub-surface
! ----------------------------------------------------------------------
! Item 576 tsurf_elev_surft

      IF (sf(576,8)) THEN
! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(nsurft,LEN_STLIST,                         &
             STLIST(1,STINDEX(1,576,8,im_index)),                       &
             PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              &
             ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage= RoutineName //                                    &
             " : error in set_pseudo_list(item 576 = tsurf_elev_surft)"
            GOTO 9999
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,nsurft
          IF (PLLTILE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
            CALL expand_from_mask (                                     &
                STASHWORK(SI(576,8,im_index)+(PSLEVEL_OUT-1)            &
                *row_length*rows),tsurf_elev_surft(1,PSLEVEL_OUT),      &
                land_sea_mask,row_length*rows,land_points)
          END IF
        END DO
      END IF

! ----------------------------------------------------------------------
! Heat flux to elevated sub-surface
! ----------------------------------------------------------------------
! Item 577 snow_soil_htf

      IF (sf(577,8)) THEN
! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(nsurft,LEN_STLIST,                         &
             STLIST(1,STINDEX(1,577,8,im_index)),                       &
             PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              &
             ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage= RoutineName //                                    &
                 " : error in set_pseudo_list(item 577 = snow_soil_htf)"
            GOTO 9999
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,nsurft
          IF (PLLTILE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
            CALL expand_from_mask (                                     &
                STASHWORK(SI(577,8,im_index)+(PSLEVEL_OUT-1)            &
                *row_length*rows),snow_soil_htf(1,PSLEVEL_OUT),         &
                land_sea_mask,row_length*rows,land_points)
          END IF
        END DO
      END IF
!
! ----------------------------------------------------------------------
! Tiled snowpack mass balance 
! ----------------------------------------------------------------------
! Item 578 snow_smb_surft

      IF (sf(578,8)) THEN
! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(nsurft,LEN_STLIST,                         &
             STLIST(1,STINDEX(1,578,8,im_index)),                       &
             PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              &
             ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage= RoutineName //                                    &
                " : error in set_pseudo_list(item 578 = snow_smb_surft)"
            GOTO 9999
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,nsurft
          IF (PLLTILE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
            CALL expand_from_mask (                                     &
                STASHWORK(SI(578,8,im_index)+(PSLEVEL_OUT-1)            &
                *row_length*rows),snow_smb_surft(1,PSLEVEL_OUT),         &
                land_sea_mask,row_length*rows,land_points)
          END IF
        END DO
      END IF
!
! ----------------------------------------------------------------------
! Soil layer moisture
! ----------------------------------------------------------------------
! Item 8 223 soil_layer_moisture

      IF (sf(223,8)) THEN

         DO k = 1, dsm_levels
            DO j= 1, rows
               DO i = 1, row_length
                  interp_data_3(i,j,k) = rmdi
               END DO
            END DO
             DO l = 1, land_points
               j=(land_index(l)-1)/row_length + 1
               i=land_index(l) - (j-1)*row_length
               interp_data_3(i,j,k) =                                         &
                    soil_layer_moisture(l,k)
             END DO
         END DO

! DEPENDS ON: copydiag_3d
         CALL copydiag_3d(stashwork(si(223,8,im_index)),interp_data_3,        &
              row_length,rows,dsm_levels,0,0,0,0,                             &
              at_extremity,                                                   &
              stlist(1,stindex(1,223,8,im_index)),len_stlist,                 &
              stash_levels,num_stash_levels+1,                                &
              atmos_im,8,223,                                                 &
              icode,cmessage)

         IF (icode  >   0) THEN
            cmessage="Error in copydiag( item 223)"
            GOTO 9999
         END IF

      END IF

! ----------------------------------------------------------------------
! Surface run-off (accumulated)
! ----------------------------------------------------------------------
! Item 204  surf_roff

      IF (sf(204,8)) THEN

         DO j = 1, rows
            DO i = 1, row_length
               interp_data(i,j) = rmdi
            END DO
         END DO

         DO l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = surf_roff(l) * timestep
         END DO

! DEPENDS ON: copydiag
         CALL copydiag(stashwork(si(204,8,im_index)),interp_data,             &
              row_length,rows,0,0,0,0, at_extremity,                          &
              atmos_im,8,204,                                                 &
              icode,cmessage)

         IF (icode  >   0) THEN
            cmessage="Error in copydiag( item 204 )"
            GOTO 9999
         END IF

      END IF

! ----------------------------------------------------------------------
! Sub-surface run-off (accumulated)
! ----------------------------------------------------------------------
! Item 205  sub_surf_roff

      IF (sf(205,8)) THEN

         DO j = 1, rows
            DO i = 1, row_length
               interp_data(i,j) = rmdi
            END DO
         END DO

         DO l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = sub_surf_roff(l)                               &
           *  timestep
         END DO

! DEPENDS ON: copydiag
         CALL copydiag(stashwork(si(205,8,im_index)),interp_data,             &
              row_length,rows,0,0,0,0, at_extremity,                          &
              atmos_im,8,205,                                                 &
              icode,cmessage)

         IF (icode  >   0) THEN
            cmessage="Error in copydiag( item 205)"
            GOTO 9999
         ENDIF
      ENDIF

! ----------------------------------------------------------------------
! Canopy water content
! ----------------------------------------------------------------------
! Item 209

      IF (sf(209,8)) THEN

         DO j = 1, rows
            DO i = 1, row_length
               interp_data(i,j) = rmdi
            END DO
         END DO

         DO l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = canopy(l)
         END DO

! DEPENDS ON: copydiag
         CALL copydiag(stashwork(si(209,8,im_index)),interp_data,             &
              row_length,rows,0,0,0,0, at_extremity,                          &
              atmos_im,8,209,                                                 &
              icode,cmessage)

         IF (icode  >   0) THEN
            cmessage="Error in copydiag( item 209)"
            GOTO 9999
         ENDIF
      ENDIF


! ----------------------------------------------------------------------
! Land snow melt heat flux (W/m2)
! ----------------------------------------------------------------------
! Item 202 SNOMLT_SURF_HTF

      IF (sf(202,8)) THEN
        CALL expand_from_mask (                                               &
            STASHWORK(SI(202,8,im_index)),                                    &
            snomlt_surf_htf,                                                  &
            land_sea_mask,row_length*rows,land_points)
      END IF


! ----------------------------------------------------------------------
! Land snow melt heat rate (Kg/m2/s)
! ----------------------------------------------------------------------
! Item 231 snow_melt

      IF (sf(231,8)) THEN
        CALL expand_from_mask (                                               &
            STASHWORK(SI(231,8,im_index)),                                    &
            snow_melt,                                                        &
            land_sea_mask,row_length*rows,land_points)
      END IF


! ----------------------------------------------------------------------
! Canopy throughfall rate (Kg/m2/s)
! ----------------------------------------------------------------------
! Item 233 tot_tfall

      IF (sf(233,8)) THEN
        CALL expand_from_mask (                                               &
            STASHWORK(SI(233,8,im_index)),                                    &
            tot_tfall,                                                        &
            land_sea_mask,row_length*rows,land_points)
      END IF


! ----------------------------------------------------------------------
! Snow amount on tiles (Kg/m2)
! ----------------------------------------------------------------------
! Item 236 snow_tile

      IF (sf(236,8)) THEN
! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(nsurft,LEN_STLIST,                               &
             STLIST(1,STINDEX(1,236,8,im_index)),                             &
             PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,                    &
             ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                         &
          "daghyd  : error in set_pseudo_list(item 236 = snow_tile)"
            GOTO 9999
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,nsurft
          IF (PLLTILE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
            CALL expand_from_mask (                                           &
                STASHWORK(SI(236,8,im_index)+(PSLEVEL_OUT-1)                  &
                *row_length*rows),snow_surft(1,PSLEVEL_OUT),                  &
                land_sea_mask,row_length*rows,land_points)
          END IF
        END DO
      END IF



! ----------------------------------------------------------------------
! Snow melt rate on tiles (Kg/m2)
! ----------------------------------------------------------------------
! Item 237 melt_tile

      IF (sf(237,8)) THEN
! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(nsurft,LEN_STLIST,                               &
             STLIST(1,STINDEX(1,237,8,im_index)),                             &
             PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,                    &
             ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                         &
          "daghyd  : error in set_pseudo_list(item 237 = melt_tile)"
            GOTO 9999
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,nsurft
          IF (PLLTILE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
            CALL expand_from_mask (                                           &
                STASHWORK(SI(237,8,im_index)+(PSLEVEL_OUT-1)                  &
                *row_length*rows),melt_surft(1,PSLEVEL_OUT),                  &
                land_sea_mask,row_length*rows,land_points)
          END IF
        END DO
      END IF



! ----------------------------------------------------------------------
! Snow grain size on tiles
! ----------------------------------------------------------------------
! Item 238 rgrain

      IF (sf(238,8)) THEN
! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(nsurft,LEN_STLIST,                               &
             STLIST(1,STINDEX(1,238,8,im_index)),                             &
             PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,                    &
             ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                         &
          "daghyd  : error in set_pseudo_list(item 238 = rgrain)"
            GOTO 9999
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,nsurft
          IF (PLLTILE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
            CALL expand_from_mask (                                           &
                STASHWORK(SI(238,8,im_index)+(PSLEVEL_OUT-1)                  &
                *row_length*rows),rgrain(1,PSLEVEL_OUT),                      &
                land_sea_mask,row_length*rows,land_points)
          END IF
        END DO
      END IF



! ----------------------------------------------------------------------
! Unfrozen soil moisture fraction
! ----------------------------------------------------------------------
! Item 229 sthu

      IF (sf(229,8)) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(DSM_LEVELS,LEN_STLIST,                           &
             STLIST(1,STINDEX(1,229,8,im_index)),                             &
             LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                         &
          "daghyd  : error in set_levels_list(item 229 = sthu)"
            GOTO 9999
        END IF
        LEVEL_OUT=0
        DO LEVEL=1,DSM_LEVELS
          IF(LIST(LEVEL)) THEN
            LEVEL_OUT=LEVEL_OUT+1
            CALL expand_from_mask (                                           &
                STASHWORK(SI(229,8,im_index)+(LEVEL_OUT-1)                    &
                *row_length*rows),sthu(1,LEVEL_OUT),                          &
                land_sea_mask,row_length*rows,land_points)
          END IF
        END DO
      END IF


! ----------------------------------------------------------------------
! Frozen soil moisture fraction
! ----------------------------------------------------------------------
! Item 230 sthf

      IF (sf(230,8)) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(DSM_LEVELS,LEN_STLIST,                           &
             STLIST(1,STINDEX(1,230,8,im_index)),                             &
             LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                         &
          "daghyd  : error in set_levels_list(item 230 = sthf)"
            GOTO 9999
        END IF
        LEVEL_OUT=0
        DO LEVEL=1,DSM_LEVELS
          IF(LIST(LEVEL)) THEN
            LEVEL_OUT=LEVEL_OUT+1
            CALL expand_from_mask (                                           &
                STASHWORK(SI(230,8,im_index)+(LEVEL_OUT-1)                    &
                *row_length*rows),sthf(1,LEVEL_OUT),                          &
                land_sea_mask,row_length*rows,land_points)
          END IF
        END DO
      END IF


! ----------------------------------------------------------------------
! Baseflow
! ----------------------------------------------------------------------
! Item 239  baseflow

      IF (sf(239,8)) THEN

         DO j = 1, rows
            DO i = 1, row_length
               interp_data(i,j) = rmdi
            END DO
         END DO

         DO l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = qbase(l)
         END DO

! DEPENDS ON: copydiag
         CALL copydiag(stashwork(si(239,8,im_index)),interp_data,             &
              row_length,rows,0,0,0,0, at_extremity,                          &
              atmos_im,8,239,                                                 &
              icode,cmessage)

         IF (icode  >   0) THEN
            cmessage="Error in copydiag( item 239)"
            GOTO 9999
         END IF

      END IF

! ----------------------------------------------------------------------
! Dunne Runoff
! ----------------------------------------------------------------------
! Item 240 Dunne Runoff

      IF (sf(240,8)) THEN

         DO j = 1, rows
            DO i = 1, row_length
               interp_data(i,j) = rmdi
            END DO
         END DO

         DO l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = dun_roff(l)
         END DO

! DEPENDS ON: copydiag
         CALL copydiag(stashwork(si(240,8,im_index)),interp_data,             &
              row_length,rows,0,0,0,0, at_extremity,                          &
              atmos_im,8,240,                                                 &
              icode,cmessage)

         IF (icode  >   0) THEN
            cmessage="Error in copydiag( item 240)"
            GOTO 9999
         END IF

      END IF

! ----------------------------------------------------------------------
! Baseflow from water table layer
! ----------------------------------------------------------------------
! Item 241  baseflow from water table layer

      IF (sf(241,8)) THEN

         DO j = 1, rows
            DO i = 1, row_length
               interp_data(i,j) = rmdi
            END DO
         END DO

         DO l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = qbase_zw(l)
         END DO

! DEPENDS ON: copydiag
         CALL copydiag(stashwork(si(241,8,im_index)),interp_data,             &
              row_length,rows,0,0,0,0, at_extremity,                          &
              atmos_im,8,241,                                                 &
              icode,cmessage)

         IF (icode  >   0) THEN
            cmessage="Error in copydiag( item 241)"
            GOTO 9999
         END IF

      END IF

! ----------------------------------------------------------------------
! Wetland methane flux used as input into UKCA
! ----------------------------------------------------------------------
! Item 242  Wetland methane flux used as input into UKCA

      IF (sf(242,8)) THEN

         DO j = 1, rows
            DO i = 1, row_length
               interp_data(i,j) = rmdi
            END DO
         END DO

         DO l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = fch4_wetl(l)
         END DO

! DEPENDS ON: copydiag
         CALL copydiag(stashwork(si(242,8,im_index)),interp_data,             &
              row_length,rows,0,0,0,0, at_extremity,                          &
              atmos_im,8,242,                                                 &
              icode,cmessage)

         IF (icode  >   0) THEN
            cmessage="Error in copydiag( item 242)"
            GOTO 9999
         END IF

      END IF
      IF (sf(243,8)) THEN
         DO j = 1, rows
            DO i = 1, row_length
               interp_data(i,j) = rmdi
            END DO
         END DO
         DO l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = ti_mean(l)
         END DO
! DEPENDS ON: copydiag
         CALL copydiag(stashwork(si(243,8,im_index)),interp_data,             &
              row_length,rows,0,0,0,0, at_extremity,                          &
              atmos_im,8,243,                                                 &
              icode,cmessage)
         IF (icode  >   0) THEN
            cmessage="Error in copydiag( item 243)"
            GOTO 9999
         END IF
      END IF
      IF (sf(244,8)) THEN
         DO j = 1, rows
            DO i = 1, row_length
               interp_data(i,j) = rmdi
            END DO
         END DO
         DO l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = ti_sig(l)
         END DO
! DEPENDS ON: copydiag
         CALL copydiag(stashwork(si(244,8,im_index)),interp_data,             &
              row_length,rows,0,0,0,0, at_extremity,                          &
              atmos_im,8,244,                                                 &
              icode,cmessage)
         IF (icode  >   0) THEN
            cmessage="Error in copydiag( item 244)"
            GOTO 9999
         END IF
      END IF
      IF (sf(251,8)) THEN
         DO j = 1, rows
            DO i = 1, row_length
               interp_data(i,j) = rmdi
            END DO
         END DO
         DO l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = fexp(l)
         END DO
! DEPENDS ON: copydiag
         CALL copydiag(stashwork(si(251,8,im_index)),interp_data,             &
              row_length,rows,0,0,0,0, at_extremity,                          &
              atmos_im,8,251,                                                 &
              icode,cmessage)
         IF (icode  >   0) THEN
            cmessage="Error in copydiag( item 251)"
            GOTO 9999
         END IF
      END IF
      IF (sf(246,8)) THEN
         DO j = 1, rows
            DO i = 1, row_length
               interp_data(i,j) = rmdi
            END DO
         END DO
         DO l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = gamtot(l)
         END DO
! DEPENDS ON: copydiag
         CALL copydiag(stashwork(si(246,8,im_index)),interp_data,             &
              row_length,rows,0,0,0,0, at_extremity,                          &
              atmos_im,8,246,                                                 &
              icode,cmessage)
         IF (icode  >   0) THEN
            cmessage="Error in copydiag( item 246)"
            GOTO 9999
         END IF
      END IF
      IF (sf(247,8)) THEN
         DO j = 1, rows
            DO i = 1, row_length
               interp_data(i,j) = rmdi
            END DO
         END DO
         DO l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = fsat(l)
         END DO
! DEPENDS ON: copydiag
         CALL copydiag(stashwork(si(247,8,im_index)),interp_data,             &
              row_length,rows,0,0,0,0, at_extremity,                          &
              atmos_im,8,247,                                                 &
              icode,cmessage)
         IF (icode  >   0) THEN
            cmessage="Error in copydiag( item 247)"
            GOTO 9999
         END IF
      END IF
      IF (sf(248,8)) THEN
         DO j = 1, rows
            DO i = 1, row_length
               interp_data(i,j) = rmdi
            END DO
         END DO
         DO l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = fwetl(l)
         END DO
! DEPENDS ON: copydiag
         CALL copydiag(stashwork(si(248,8,im_index)),interp_data,             &
              row_length,rows,0,0,0,0, at_extremity,                          &
              atmos_im,8,248,                                                 &
              icode,cmessage)
         IF (icode  >   0) THEN
            cmessage="Error in copydiag( item 248)"
            GOTO 9999
         END IF
      END IF
      IF (sf(249,8)) THEN
         DO j = 1, rows
            DO i = 1, row_length
               interp_data(i,j) = rmdi
            END DO
         END DO
         DO l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = zw(l)
         END DO
! DEPENDS ON: copydiag
         CALL copydiag(stashwork(si(249,8,im_index)),interp_data,             &
              row_length,rows,0,0,0,0, at_extremity,                          &
              atmos_im,8,249,                                                 &
              icode,cmessage)
         IF (icode  >   0) THEN
            cmessage="Error in copydiag( item 249)"
            GOTO 9999
         END IF
      END IF
      IF (sf(250,8)) THEN
         DO j = 1, rows
            DO i = 1, row_length
               interp_data(i,j) = rmdi
            END DO
         END DO
         DO l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = sthzw(l)
         END DO
! DEPENDS ON: copydiag
         CALL copydiag(stashwork(si(250,8,im_index)),interp_data,             &
              row_length,rows,0,0,0,0, at_extremity,                          &
              atmos_im,8,250,                                                 &
              icode,cmessage)
         IF (icode  >   0) THEN
            cmessage="Error in copydiag( item 250)"
            GOTO 9999
         END IF
      END IF
! ----------------------------------------------------------------------
! Drainage out of  "nshyd"th model layer
! ----------------------------------------------------------------------
! Item 252  drainage

      IF (sf(252,8)) THEN

         DO j = 1, rows
            DO i = 1, row_length
               interp_data(i,j) = rmdi
            END DO
         END DO

         DO l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = drain(l)
         END DO

! DEPENDS ON: copydiag
         CALL copydiag(stashwork(si(252,8,im_index)),interp_data,             &
              row_length,rows,0,0,0,0, at_extremity,                          &
              atmos_im,8,252,                                                 &
              icode,cmessage)

         IF (icode  >   0) THEN
            cmessage="Error in copydiag( item 252)"
            GOTO 9999
         END IF

      END IF

! ----------------------------------------------------------------------
! Output inland basin outflow on atmosphere grid

! --------------------------------------------------------------------
! Inland basin outflow (atmos grid)
! --------------------------------------------------------------------
! Item 245  Inland basin outflow

      IF (sf(245,8)) THEN

         DO j = 1, rows
            DO i = 1, row_length
               interp_data(i,j) = rmdi
            END DO
         END DO

         DO l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = inlandout_atm(l)
         END DO

! DEPENDS ON: copydiag
         CALL copydiag(stashwork(si(245,8,im_index)),interp_data,             &
              row_length,rows,0,0,0,0, at_extremity,                          &
              atmos_im,8,245,                                                 &
              icode,cmessage)

         IF (icode  >   0) THEN
            cmessage="Error in copydiag( item 245)"
            GOTO 9999
         END IF

      END IF
! ----------------------------------------------------------------------
! Wetland methane flux using soil carbon-based substrate
! ----------------------------------------------------------------------
! Item 260  Wetland methane flux using soil carbon-based substrate

      If (sf(260,8)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do

         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = fch4_wetl_cs(l)
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(260,8,im_index)),interp_data,             &
              row_length,rows,0,0,0,0, at_extremity,                          &
              atmos_im,8,260,                                                 &
              icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 260)"
            goto 9999
         End if

      End if
! ----------------------------------------------------------------------
! Wetland methane flux using NPP-based substrate
! ----------------------------------------------------------------------
! Item 261  Wetland methane flux using NPP-based substrate

      If (sf(261,8)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do

         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = fch4_wetl_npp(l)
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(261,8,im_index)),interp_data,             &
              row_length,rows,0,0,0,0, at_extremity,                          &
              atmos_im,8,261,                                                 &
              icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 261)"
            goto 9999
         End if

      End if
! ----------------------------------------------------------------------
! Wetland methane flux using soil-respiration-based substrate
! ----------------------------------------------------------------------
! Item 262  Wetland methane flux using soil-respiration-based substrate

      If (sf(262,8)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do

         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = fch4_wetl_resps(l)
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(262,8,im_index)),interp_data,             &
              row_length,rows,0,0,0,0, at_extremity,                          &
              atmos_im,8,262,                                                 &
              icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 262)"
            goto 9999
         End if

      End if
! -------------------------------------------
 9999 CONTINUE
      IF (icode /= 0) THEN

        CALL Ereport(RoutineName,icode,Cmessage)
      ENDIF

      IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE diagnostics_hyd
#endif
