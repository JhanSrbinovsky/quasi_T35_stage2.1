#if defined(UM_JULES)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine diagnostics_veg ----------------------------------------
!
! Purpose : Calculates diagnostics for dynamic vegetation and
!           outputs them.
!
! -----------------------------------------------------------------
!
! Subroutine diagnostics_veg
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: vegetation

      SUBROUTINE diagnostics_veg(                                             &
                             row_length, rows, n_rows                         &
      ,                      global_row_length, global_rows                   &
      ,                      DIM_CS1, DIM_CS2                                 &
      ,                      halo_i, halo_j, off_x, off_y, me                 &
      ,                      n_proc, n_procx, n_procy                         &
      ,                      g_rows, g_row_length                             &
      ,                      at_extremity                                     &
      ,                      land_pts                                         &
      ,                      land_index                                       &
      ,                      ntype,npft                                       &
      ,                      c_veg,cv,g_leaf_phen                             &
      ,                      lit_c,lit_c_mn,g_leaf_day                        &
      ,                      lai_phen,g_leaf_dr_out,npp_dr_out                &
      ,                      resp_w_dr_out,resp_s_dr_out,frac_disturb         &
      ,                      disturb_veg_prev,wood_prod_fast_d1               &
      ,                      wood_prod_med_d1,wood_prod_slow_d1               &
      ,                      frac,lai,ht,cs                                   &
      ,                      STASHwork                                  )


! Purpose:
!          Calculates diagnostics and outputs them.
!
      USE submodel_mod, ONLY: internal_model_index, atmos_im
      USE stash_array_mod, ONLY:                                              &
          sf, si, stlist, stindex, len_stlist, stash_pseudo_levels,           &
          num_stash_pseudo, stash_levels, num_stash_levels
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE ereport_mod, ONLY : ereport
      USE missing_data_mod, ONLY: rmdi
      USE errormessagelength_mod, ONLY: errormessagelength
      USE trif_vars_mod, ONLY: WP_fast_in_gb,  WP_med_in_gb,  WP_slow_in_gb,  &
                               WP_fast_out_gb, WP_med_out_gb, WP_slow_out_gb

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
      , land_pts                                                              &
                 ! No.of land points being processed, can be 0.
      , ntype                                                                 &
                    ! Max. No. of land surface tiles
      , npft                                                                  &
                    ! No. of plant functional types
      , DIM_CS1, DIM_CS2     ! soil carbon dimensions

      REAL                                                                    &
        lat_rot_NP                                                            &
      , long_rot_NP


! Primary Arrays used in all models
      INTEGER                                                                 &
        land_index(land_pts)      ! set from land_sea_mask


      REAL                                                                    &
       C_VEG(land_pts,NPFT)                                                   &
                            ! Total carbon content of vegetation
!                              ! (kg C/m2).
      ,CV(land_pts)                                                           &
                            ! Gridbox mean veg carbon (kg C/m2).
      ,LIT_C(land_pts,NPFT)                                                   &
                            ! Carbon Litter (kg C/m2/360days).
      ,LIT_C_MN(land_pts)                                                     &
                            ! Gridbox mean carbon litter
!                              ! (kg C/m2/360days)
      ,G_LEAF_DAY(land_pts,NPFT)                                              &
                                   ! Mean leaf turnover rate for
!                                     ! input to PHENOL (/360days).
      ,G_LEAF_PHEN(land_pts,NPFT)                                             &
                                   ! Mean leaf turnover rate over
!                                     ! phenology period (/360days).
      ,G_LEAF_DR_OUT(land_pts,NPFT)                                           &
                                   ! Mean leaf turnover rate for
!                                     ! driving TRIFFID (/360days).
      ,LAI_PHEN(land_pts,NPFT)                                                &
                                   ! LAI of PFTs after phenology.
      ,NPP_DR_OUT(land_pts,NPFT)                                              &
                                   ! Mean NPP for driving TRIFFID
!                                     ! (kg C/m2/360days).
      ,RESP_W_DR_OUT(land_pts,NPFT)                                           &
                                   ! Mean wood respiration for
!                                     ! driving TRIFFID
!                                     ! (kg C/m2/360days).
      ,RESP_S_DR_OUT(land_pts,DIM_CS1+1)                                      &
                                         ! Mean soil respiration for
!                                     ! driving TRIFFID
!                                     ! (kg C/m2/360days).
      ,FRAC_DISTURB(land_pts)                                                 &
                                   ! Fraction of gridbox in which
!                                     !    vegetation is disturbed.
      ,DISTURB_VEG_PREV(land_pts)                                             &
                                   ! Previous disturbed fraction
      ,WOOD_PROD_FAST_D1(land_pts)                                            &
                                   ! C flux from wood product pool
      ,WOOD_PROD_MED_D1(land_pts)                                             &
                                   ! C flux from wood product pool
      ,WOOD_PROD_SLOW_D1(land_pts)                                            &
                                   ! C flux from wood product pool
      ,FRAC(land_pts,NTYPE)                                                   &
                                   ! Fractions of surface types.
      ,LAI(land_pts,NPFT)                                                     &
                                   ! LAI of plant functional
!                                     !       types.
      ,HT(land_pts,NPFT)                                                      &
                                   ! Height of plant functional
!                                     !       types (m).
      ,CS(land_pts,DIM_CS1)   ! Soil carbon content
!                                     !       (kg C/m2).

! Diagnostic variables
       REAL                                                                   &
        STASHwork(*)    ! STASH workspace


! Local variables

      LOGICAL                                                                 &
       PLLTYPE(NTYPE)                                                         &
                          ! pseudolevel list for surface types
      ,PLLPFT(NPFT)       ! pseudolevel list for PFTs

      INTEGER                                                                 &
       PSLEVEL                                                                &
                     !  loop counter for pseudolevels
      ,PSLEVEL_OUT   !  index for pseudolevels sent to STASH

      INTEGER                                                                 &
        i, j, k, l                                                            &
      ,    icode                ! Return code  =0 Normal exit  >1 Error

      CHARACTER(LEN=errormessagelength)  cmessage
      CHARACTER(LEN=*) RoutineName
      PARAMETER ( RoutineName='DIAGNOSTICS_VEG')

      INTEGER                                                                 &
        im_index        ! internal model index

      REAL                                                                    &
        interp_data(row_length,rows)

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

!  ITEM 1: VEGETATION CARBON ON PLANT FUNCTIONAL TYPES

      IF (sf(1,19)) THEN

! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NPFT,LEN_STLIST,                                 &
             STLIST(1,STINDEX(1,1,19,im_index)),                              &
             PLLPFT,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,                     &
             ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                         &
          "dagveg  : error in set_pseudo_list(item 1 = c_veg)"
            GOTO 9999
        END IF
        PSLEVEL_OUT=0

        DO PSLEVEL=1,NPFT
          IF (PLLPFT(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
            DO j = 1, rows
              DO i = 1, row_length
                interp_data(i,j) = rmdi
              END DO
            END DO

            DO l = 1, land_pts
              j=(land_index(l)-1)/row_length + 1
              i=land_index(l) - (j-1)*row_length
              interp_data(i,j) = c_veg(l,pslevel_out)
            END DO

! DEPENDS ON: copydiag
            CALL copydiag (STASHwork(si(1,19,im_index)+(pslevel_out-1)        &
                 *row_length*rows),interp_data,                               &
                 row_length,rows,0,0,0,0, at_extremity,                       &
                 atmos_im,19,1,                                               &
                 icode,cmessage)

            IF (icode  >   0) THEN
               cmessage="Error in copydiag( item 1901)"
               GOTO 9999
            END IF
          ENDIF
        ENDDO

      END IF     !   sf(1,19)


!  ITEM 2: GRIDBOX MEAN VEGETATION CARBON

      IF (sf(2,19)) THEN
        DO j = 1, rows
          DO i = 1, row_length
            interp_data(i,j) = rmdi
          END DO
        END DO

        DO l = 1, land_pts
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          interp_data(i,j) = cv(l)
        END DO

! DEPENDS ON: copydiag
        CALL copydiag (STASHwork(si(2,19,im_index)),interp_data,              &
             row_length,rows,0,0,0,0, at_extremity,                           &
             atmos_im,19,2,                                                   &
             icode,cmessage)

        IF (icode  >   0) THEN
           cmessage="Error in copydiag( item 1902)"
           GOTO 9999
        END IF

      END IF     !   sf(2,19)

!  ITEM 3: PHENOLOGICAL LEAF TURNOVER RATE PFTS

      IF (sf(3,19)) THEN

! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NPFT,LEN_STLIST,                                 &
             STLIST(1,STINDEX(1,3,19,im_index)),                              &
             PLLPFT,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,                     &
             ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                         &
          "dagveg  : error in set_pseudo_list(item 3 = g_leaf_phen)"
            GOTO 9999
        END IF
        PSLEVEL_OUT=0

        DO PSLEVEL=1,NPFT
          IF (PLLPFT(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
            DO j = 1, rows
              DO i = 1, row_length
                interp_data(i,j) = rmdi
              END DO
            END DO

            DO l = 1, land_pts
              j=(land_index(l)-1)/row_length + 1
              i=land_index(l) - (j-1)*row_length
              interp_data(i,j) = g_leaf_phen(l,pslevel_out)
            END DO

! DEPENDS ON: copydiag
            CALL copydiag (STASHwork(si(3,19,im_index)+(pslevel_out-1)        &
                 *row_length*rows),interp_data,                               &
                 row_length,rows,0,0,0,0, at_extremity,                       &
                 atmos_im,19,3,                                               &
                 icode,cmessage)

            IF (icode  >   0) THEN
               cmessage="Error in copydiag( item 1903)"
               GOTO 9999
            END IF
          ENDIF
        ENDDO

      END IF     !   sf(3,19)


!  ITEM 4: LITTER CARBON ON PLANT FUNCTIONAL TYPES

      IF (sf(4,19)) THEN

! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NPFT,LEN_STLIST,                                 &
             STLIST(1,STINDEX(1,4,19,im_index)),                              &
             PLLPFT,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,                     &
             ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                         &
          "dagveg  : error in set_pseudo_list(item 4 = lit_c)"
            GOTO 9999
        END IF
        PSLEVEL_OUT=0

        DO PSLEVEL=1,NPFT
          IF (PLLPFT(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
            DO j = 1, rows
              DO i = 1, row_length
                interp_data(i,j) = rmdi
              END DO
            END DO

            DO l = 1, land_pts
              j=(land_index(l)-1)/row_length + 1
              i=land_index(l) - (j-1)*row_length
              interp_data(i,j) = lit_c(l,pslevel_out)
            END DO

! DEPENDS ON: copydiag
            CALL copydiag (STASHwork(si(4,19,im_index)+(pslevel_out-1)        &
                 *row_length*rows),interp_data,                               &
                 row_length,rows,0,0,0,0, at_extremity,                       &
                 atmos_im,19,4,                                               &
                 icode,cmessage)

            IF (icode  >   0) THEN
               cmessage="Error in copydiag( item 1904)"
               GOTO 9999
            END IF
          ENDIF
        ENDDO

      END IF     !   sf(4,19)


!  ITEM 5: GRIDBOX MEAN LITTER CARBON

      IF (sf(5,19)) THEN
        DO j = 1, rows
          DO i = 1, row_length
            interp_data(i,j) = rmdi
          END DO
        END DO

        DO l = 1, land_pts
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          interp_data(i,j) = lit_c_mn(l)
        END DO

! DEPENDS ON: copydiag
        CALL copydiag (STASHwork(si(5,19,im_index)),interp_data,              &
             row_length,rows,0,0,0,0, at_extremity,                           &
             atmos_im,19,5,                                                   &
             icode,cmessage)

        IF (icode  >   0) THEN
           cmessage="Error in copydiag( item 1905)"
           GOTO 9999
        END IF

      END IF     !   sf(5,19)

!  ITEM 6: MEAN LEAF TURNOVER RATE ON PFTS FOR PHENOLOGY

      IF (sf(6,19)) THEN

! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NPFT,LEN_STLIST,                                 &
             STLIST(1,STINDEX(1,6,19,im_index)),                              &
             PLLPFT,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,                     &
             ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                         &
          "dagveg  : error in set_pseudo_list(item 6 = g_leaf_day)"
            GOTO 9999
        END IF
        PSLEVEL_OUT=0

        DO PSLEVEL=1,NPFT
          IF (PLLPFT(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
            DO j = 1, rows
              DO i = 1, row_length
                interp_data(i,j) = rmdi
              END DO
            END DO

            DO l = 1, land_pts
              j=(land_index(l)-1)/row_length + 1
              i=land_index(l) - (j-1)*row_length
              interp_data(i,j) = g_leaf_day(l,pslevel_out)
            END DO

! DEPENDS ON: copydiag
            CALL copydiag (STASHwork(si(6,19,im_index)+(pslevel_out-1)        &
                 *row_length*rows),interp_data,                               &
                 row_length,rows,0,0,0,0, at_extremity,                       &
                 atmos_im,19,6,                                               &
                 icode,cmessage)

            IF (icode  >   0) THEN
               cmessage="Error in copydiag( item 1906)"
               GOTO 9999
            END IF
          ENDIF
        ENDDO

      END IF     !   sf(6,19)


!  ITEM 7: LEAF AREA INDEX ON PLANT FUNCTIONAL TYPES AFTER PHENOLOGY

      IF (sf(7,19)) THEN

! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NPFT,LEN_STLIST,                                 &
             STLIST(1,STINDEX(1,7,19,im_index)),                              &
             PLLPFT,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,                     &
             ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                         &
          "dagveg  : error in set_pseudo_list(item 7 = lai_phen)"
            GOTO 9999
        END IF
        PSLEVEL_OUT=0

        DO PSLEVEL=1,NPFT
          IF (PLLPFT(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
            DO j = 1, rows
              DO i = 1, row_length
                interp_data(i,j) = rmdi
              END DO
            END DO

            DO l = 1, land_pts
              j=(land_index(l)-1)/row_length + 1
              i=land_index(l) - (j-1)*row_length
              interp_data(i,j) = lai_phen(l,pslevel_out)
            END DO

! DEPENDS ON: copydiag
            CALL copydiag (STASHwork(si(7,19,im_index)+(pslevel_out-1)        &
                 *row_length*rows),interp_data,                               &
                 row_length,rows,0,0,0,0, at_extremity,                       &
                 atmos_im,19,7,                                               &
                 icode,cmessage)

            IF (icode  >   0) THEN
               cmessage="Error in copydiag( item 1907)"
               GOTO 9999
            END IF
          ENDIF
        ENDDO

      END IF     !   sf(7,19)


!  ITEM 8: MEAN LEAF TURNOVER RATE ON PFTS FOR TRIFFID

      IF (sf(8,19)) THEN

! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NPFT,LEN_STLIST,                                 &
             STLIST(1,STINDEX(1,8,19,im_index)),                              &
             PLLPFT,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,                     &
             ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                         &
          "dagveg  : error in set_pseudo_list(item 8 = g_leaf_dr_out)"
            GOTO 9999
        END IF
        PSLEVEL_OUT=0

        DO PSLEVEL=1,NPFT
          IF (PLLPFT(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
            DO j = 1, rows
              DO i = 1, row_length
                interp_data(i,j) = rmdi
              END DO
            END DO

            DO l = 1, land_pts
              j=(land_index(l)-1)/row_length + 1
              i=land_index(l) - (j-1)*row_length
              interp_data(i,j) = g_leaf_dr_out(l,pslevel_out)
            END DO

! DEPENDS ON: copydiag
            CALL copydiag (STASHwork(si(8,19,im_index)+(pslevel_out-1)        &
                 *row_length*rows),interp_data,                               &
                 row_length,rows,0,0,0,0, at_extremity,                       &
                 atmos_im,19,8,                                               &
                 icode,cmessage)

            IF (icode  >   0) THEN
               cmessage="Error in copydiag( item 1908)"
               GOTO 9999
            END IF
          ENDIF
        ENDDO

      END IF     !   sf(8,19)


!  ITEM 9: MEAN NPP ON PFTS FOR TRIFFID

      IF (sf(9,19)) THEN

! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NPFT,LEN_STLIST,                                 &
             STLIST(1,STINDEX(1,9,19,im_index)),                              &
             PLLPFT,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,                     &
             ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                         &
          "dagveg  : error in set_pseudo_list(item 9 = npp_dr_out)"
            GOTO 9999
        END IF
        PSLEVEL_OUT=0

        DO PSLEVEL=1,NPFT
          IF (PLLPFT(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
            DO j = 1, rows
              DO i = 1, row_length
                interp_data(i,j) = rmdi
              END DO
            END DO

            DO l = 1, land_pts
              j=(land_index(l)-1)/row_length + 1
              i=land_index(l) - (j-1)*row_length
              interp_data(i,j) = npp_dr_out(l,pslevel_out)
            END DO

! DEPENDS ON: copydiag
            CALL copydiag (STASHwork(si(9,19,im_index)+(pslevel_out-1)        &
                 *row_length*rows),interp_data,                               &
                 row_length,rows,0,0,0,0, at_extremity,                       &
                 atmos_im,19,9,                                               &
                 icode,cmessage)

            IF (icode  >   0) THEN
               cmessage="Error in copydiag( item 1909)"
               GOTO 9999
            END IF
          ENDIF
        ENDDO

      END IF     !   sf(9,19)


!  ITEM 10: MEAN WOOD RESPIRATION ON PFTS FOR TRIFFID

      IF (sf(10,19)) THEN

! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NPFT,LEN_STLIST,                                 &
             STLIST(1,STINDEX(1,10,19,im_index)),                             &
             PLLPFT,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,                     &
             ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                         &
          "dagveg  : error in set_pseudo_list(item 10 = resp_w_dr_out)"
            GOTO 9999
        END IF
        PSLEVEL_OUT=0

        DO PSLEVEL=1,NPFT
          IF (PLLPFT(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
            DO j = 1, rows
              DO i = 1, row_length
                interp_data(i,j) = rmdi
              END DO
            END DO

            DO l = 1, land_pts
              j=(land_index(l)-1)/row_length + 1
              i=land_index(l) - (j-1)*row_length
              interp_data(i,j) = resp_w_dr_out(l,pslevel_out)
            END DO

! DEPENDS ON: copydiag
            CALL copydiag (STASHwork(si(10,19,im_index)+(pslevel_out-1)       &
                 *row_length*rows),interp_data,                               &
                 row_length,rows,0,0,0,0, at_extremity,                       &
                 atmos_im,19,10,                                              &
                 icode,cmessage)

            IF (icode  >   0) THEN
               cmessage="Error in copydiag( item 1910)"
               GOTO 9999
            END IF
          ENDIF
        ENDDO

      END IF     !   sf(10,19)


!  ITEM 11: MEAN SOIL RESPIRATION FOR TRIFFID

      IF (sf(11,19)) THEN
        DO j = 1, rows
          DO i = 1, row_length
            interp_data(i,j) = rmdi
          END DO
        END DO

        DO l = 1, land_pts
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          interp_data(i,j) = resp_s_dr_out(l,5)
        END DO

! DEPENDS ON: copydiag
        CALL copydiag (STASHwork(si(11,19,im_index)),interp_data,             &
             row_length,rows,0,0,0,0, at_extremity,                           &
             atmos_im,19,11,                                                  &
             icode,cmessage)

        IF (icode  >   0) THEN
           cmessage="Error in copydiag( item 1911)"
           GOTO 9999
        END IF

      END IF     !   sf(11,19)

!  ITEM 12: DISTURBED FRACTION OF VEGETATION

      IF (sf(12,19)) THEN
        DO j = 1, rows
          DO i = 1, row_length
            interp_data(i,j) = rmdi
          END DO
        END DO

        DO l = 1, land_pts
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          interp_data(i,j) = frac_disturb(l)
        END DO

! DEPENDS ON: copydiag
        CALL copydiag (STASHwork(si(12,19,im_index)),interp_data,             &
             row_length,rows,0,0,0,0, at_extremity,                           &
             atmos_im,19,12,                                                  &
             icode,cmessage)

        IF (icode  >   0) THEN
           cmessage="Error in copydiag( item 1912)"
           GOTO 9999
        END IF

      END IF     !   sf(12,19)

!  ITEM 13: SURFACE TYPE FRACTIONS AFTER TRIFFID

      IF (sf(13,19)) THEN

! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NTYPE,LEN_STLIST,                                &
             STLIST(1,STINDEX(1,13,19,im_index)),                             &
             PLLTYPE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,                    &
             ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                         &
          "dagveg  : error in set_pseudo_list(item 13 = frac)"
            GOTO 9999
        END IF
        PSLEVEL_OUT=0

        DO PSLEVEL=1,NTYPE
          IF (PLLTYPE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
            DO j = 1, rows
              DO i = 1, row_length
                interp_data(i,j) = rmdi
              END DO
            END DO

            DO l = 1, land_pts
              j=(land_index(l)-1)/row_length + 1
              i=land_index(l) - (j-1)*row_length
              interp_data(i,j) = frac(l,pslevel_out)
            END DO

! DEPENDS ON: copydiag
            CALL copydiag (STASHwork(si(13,19,im_index)+(pslevel_out-1)       &
                 *row_length*rows),interp_data,                               &
                 row_length,rows,0,0,0,0, at_extremity,                       &
                 atmos_im,19,13,                                              &
                 icode,cmessage)

            IF (icode  >   0) THEN
               cmessage="Error in copydiag( item 1913)"
               GOTO 9999
            END IF
          ENDIF
        ENDDO

      END IF     !   sf(13,19)


!  ITEM 14: LEAF AREA INDEX ON PLANT FUNCTIONAL TYPES AFTER TRIFFID

      IF (sf(14,19)) THEN

! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NPFT,LEN_STLIST,                                 &
             STLIST(1,STINDEX(1,14,19,im_index)),                             &
             PLLPFT,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,                     &
             ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                         &
          "dagveg  : error in set_pseudo_list(item 14 = lai)"
            GOTO 9999
        END IF
        PSLEVEL_OUT=0

        DO PSLEVEL=1,NPFT
          IF (PLLPFT(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
            DO j = 1, rows
              DO i = 1, row_length
                interp_data(i,j) = rmdi
              END DO
            END DO

            DO l = 1, land_pts
              j=(land_index(l)-1)/row_length + 1
              i=land_index(l) - (j-1)*row_length
              interp_data(i,j) = lai(l,pslevel_out)
            END DO

! DEPENDS ON: copydiag
            CALL copydiag (STASHwork(si(14,19,im_index)+(pslevel_out-1)       &
                 *row_length*rows),interp_data,                               &
                 row_length,rows,0,0,0,0, at_extremity,                       &
                 atmos_im,19,14,                                              &
                 icode,cmessage)

            IF (icode  >   0) THEN
               cmessage="Error in copydiag( item 1914)"
               GOTO 9999
            END IF
          ENDIF
        ENDDO

      END IF     !   sf(14,19)


!  ITEM 15: CANOPY HEIGHT ON PLANT FUNCTIONAL TYPES AFTER TRIFFID

      IF (sf(15,19)) THEN

! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NPFT,LEN_STLIST,                                 &
             STLIST(1,STINDEX(1,15,19,im_index)),                             &
             PLLPFT,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,                     &
             ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                         &
          "dagveg  : error in set_pseudo_list(item 15 = ht)"
            GOTO 9999
        END IF
        PSLEVEL_OUT=0

        DO PSLEVEL=1,NPFT
          IF (PLLPFT(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
            DO j = 1, rows
              DO i = 1, row_length
                interp_data(i,j) = rmdi
              END DO
            END DO

            DO l = 1, land_pts
              j=(land_index(l)-1)/row_length + 1
              i=land_index(l) - (j-1)*row_length
              interp_data(i,j) = ht(l,pslevel_out)
            END DO

! DEPENDS ON: copydiag
            CALL copydiag (STASHwork(si(15,19,im_index)+(pslevel_out-1)       &
                 *row_length*rows),interp_data,                               &
                 row_length,rows,0,0,0,0, at_extremity,                       &
                 atmos_im,19,15,                                              &
                 icode,cmessage)

            IF (icode  >   0) THEN
               cmessage="Error in copydiag( item 1915)"
               GOTO 9999
            END IF
          ENDIF
        ENDDO

      END IF     !   sf(15,19)

!  ITEM 16: SOIL CARBON CONTENT AFTER TRIFFID

      IF (sf(16,19)) THEN
        DO j = 1, rows
          DO i = 1, row_length
            interp_data(i,j) = rmdi
          END DO
        END DO

        DO l = 1, land_pts
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          interp_data(i,j) = cs(l,1) + cs(l,2) + cs(l,3) + cs(l,4)
        END DO

! DEPENDS ON: copydiag
        CALL copydiag (STASHwork(si(16,19,im_index)),interp_data,             &
             row_length,rows,0,0,0,0, at_extremity,                           &
             atmos_im,19,16,                                                  &
             icode,cmessage)

        IF (icode  >   0) THEN
           cmessage="Error in copydiag( item 1916)"
           GOTO 9999
        END IF

      END IF     !   sf(16,19)


!  ITEM 17-20: MEAN SOIL RESPIRATION FOR TRIFFID, INDIVID. POOLS
! 17: DPM
      IF (sf(17,19)) THEN
        DO j = 1, rows
          DO i = 1, row_length
            interp_data(i,j) = rmdi
          END DO
        END DO

        DO l = 1, land_pts
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          interp_data(i,j) = resp_s_dr_out(l,1)
        END DO

! DEPENDS ON: copydiag
        CALL copydiag (STASHwork(si(17,19,im_index)),interp_data,             &
             row_length,rows,0,0,0,0, at_extremity,                           &
             atmos_im,19,17,                                                  &
             icode,cmessage)

        IF (icode  >   0) THEN
           cmessage="Error in copydiag( item 1917)"
           GOTO 9999
        END IF
      END IF

! 18: RPM
      IF (sf(18,19)) THEN
        DO j = 1, rows
          DO i = 1, row_length
            interp_data(i,j) = rmdi
          END DO
        END DO

        DO l = 1, land_pts
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          interp_data(i,j) = resp_s_dr_out(l,2)
        END DO

! DEPENDS ON: copydiag
        CALL copydiag (STASHwork(si(18,19,im_index)),interp_data,             &
             row_length,rows,0,0,0,0, at_extremity,                           &
             atmos_im,19,18,                                                  &
             icode,cmessage)

        IF (icode  >   0) THEN
           cmessage="Error in copydiag( item 1918)"
           GOTO 9999
        END IF
      END IF

! 19: BIO
      IF (sf(19,19)) THEN
        DO j = 1, rows
          DO i = 1, row_length
            interp_data(i,j) = rmdi
          END DO
        END DO

        DO l = 1, land_pts
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          interp_data(i,j) = resp_s_dr_out(l,3)
        END DO

! DEPENDS ON: copydiag
        CALL copydiag (STASHwork(si(19,19,im_index)),interp_data,             &
             row_length,rows,0,0,0,0, at_extremity,                           &
             atmos_im,19,19,                                                  &
             icode,cmessage)

        IF (icode  >   0) THEN
           cmessage="Error in copydiag( item 1919)"
           GOTO 9999
        END IF
      END IF

! 20: HUM
      IF (sf(20,19)) THEN
        DO j = 1, rows
          DO i = 1, row_length
            interp_data(i,j) = rmdi
          END DO
        END DO

        DO l = 1, land_pts
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          interp_data(i,j) = resp_s_dr_out(l,4)
        END DO

! DEPENDS ON: copydiag
        CALL copydiag (STASHwork(si(20,19,im_index)),interp_data,             &
             row_length,rows,0,0,0,0, at_extremity,                           &
             atmos_im,19,20,                                                  &
             icode,cmessage)

        IF (icode  >   0) THEN
           cmessage="Error in copydiag( item 1920)"
           GOTO 9999
        END IF
      END IF

!  ITEM 21-24: SOIL CARBON CONTENT AFTER TRIFFID, INDIVID. POOLS
! 21: DPM
      IF (sf(21,19)) THEN
        DO j = 1, rows
          DO i = 1, row_length
            interp_data(i,j) = rmdi
          END DO
        END DO

        DO l = 1, land_pts
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          interp_data(i,j) = cs(l,1)
        END DO

! DEPENDS ON: copydiag
        CALL copydiag (STASHwork(si(21,19,im_index)),interp_data,             &
             row_length,rows,0,0,0,0, at_extremity,                           &
             atmos_im,19,21,                                                  &
             icode,cmessage)

        IF (icode  >   0) THEN
           cmessage="Error in copydiag( item 1921)"
           GOTO 9999
        END IF
      END IF

! 22: RPM
      IF (sf(22,19)) THEN
        DO j = 1, rows
          DO i = 1, row_length
            interp_data(i,j) = rmdi
          END DO
        END DO

        DO l = 1, land_pts
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          interp_data(i,j) = cs(l,2)
        END DO

! DEPENDS ON: copydiag
        CALL copydiag (STASHwork(si(22,19,im_index)),interp_data,             &
             row_length,rows,0,0,0,0, at_extremity,                           &
             atmos_im,19,22,                                                  &
             icode,cmessage)

        IF (icode  >   0) THEN
           cmessage="Error in copydiag( item 1922)"
           GOTO 9999
        END IF
      END IF

! 23: BIO
      IF (sf(23,19)) THEN
        DO j = 1, rows
          DO i = 1, row_length
            interp_data(i,j) = rmdi
          END DO
        END DO

        DO l = 1, land_pts
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          interp_data(i,j) = cs(l,3)
        END DO

! DEPENDS ON: copydiag
        CALL copydiag (STASHwork(si(23,19,im_index)),interp_data,             &
             row_length,rows,0,0,0,0, at_extremity,                           &
             atmos_im,19,23,                                                  &
             icode,cmessage)

        IF (icode  >   0) THEN
           cmessage="Error in copydiag( item 1923)"
           GOTO 9999
        END IF
      END IF

! 24: HUM
      IF (sf(24,19)) THEN
        DO j = 1, rows
          DO i = 1, row_length
            interp_data(i,j) = rmdi
          END DO
        END DO

        DO l = 1, land_pts
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          interp_data(i,j) = cs(l,4)
        END DO

! DEPENDS ON: copydiag
        CALL copydiag (STASHwork(si(24,19,im_index)),interp_data,             &
             row_length,rows,0,0,0,0, at_extremity,                           &
             atmos_im,19,24,                                                  &
             icode,cmessage)

        IF (icode  >   0) THEN
           cmessage="Error in copydiag( item 1924)"
           GOTO 9999
        END IF
      END IF


! 19031: Previous agricultural fraction
! Agricultural fraction from previous call to TRIFFID.

      IF (sf(31,19)) THEN
        DO j = 1, rows
          DO i = 1, row_length
            interp_data(i,j) = rmdi
          END DO
        END DO

        DO l = 1, land_pts
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          interp_data(i,j) = disturb_veg_prev(l)
        END DO

! DEPENDS ON: copydiag
        CALL copydiag (STASHwork(si(31,19,im_index)),interp_data,             &
             row_length,rows,0,0,0,0, at_extremity,                           &
             atmos_im,19,31,                                                  &
             icode,cmessage)

        IF (icode  >   0) THEN
           cmessage="Error in copydiag( item 19031)"
           GOTO 9999
        END IF

      END IF     !   sf(31,19)


! 19032: Wood product pool Carbon (FAST turnover rate pool), kgC/m2
! Fast-turnover wood product carbon pool from cleared vegetation

      IF (sf(32,19)) THEN
        DO j = 1, rows
          DO i = 1, row_length
            interp_data(i,j) = rmdi
          END DO
        END DO

        DO l = 1, land_pts
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          interp_data(i,j) = wood_prod_fast_d1(l)
        END DO

! DEPENDS ON: copydiag
        CALL copydiag (STASHwork(si(32,19,im_index)),interp_data,             &
             row_length,rows,0,0,0,0, at_extremity,                           &
             atmos_im,19,32,                                                  &
             icode,cmessage)

        IF (icode  >   0) THEN
           cmessage="Error in copydiag( item 19032)"
           GOTO 9999
        END IF

      END IF     !   sf(32,19)


! 19033: Wood product pool Carbon (MEDIUM turnover rate pool), kgC/m2
! Medium-turnover wood product carbon pool from cleared vegetation

      IF (sf(33,19)) THEN
        DO j = 1, rows
          DO i = 1, row_length
            interp_data(i,j) = rmdi
          END DO
        END DO

        DO l = 1, land_pts
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          interp_data(i,j) = wood_prod_med_d1(l)
        END DO

! DEPENDS ON: copydiag
        CALL copydiag (STASHwork(si(33,19,im_index)),interp_data,             &
             row_length,rows,0,0,0,0, at_extremity,                           &
             atmos_im,19,33,                                                  &
             icode,cmessage)

        IF (icode  >   0) THEN
           cmessage="Error in copydiag( item 19033)"
           GOTO 9999
        END IF

      END IF     !   sf(33,19)

! 19034: Wood product pool Carbon (SLOW turnover rate pool), kgC/m2
! Slow-turnover wood product carbon pool from cleared vegetation

      IF (sf(34,19)) THEN
        DO j = 1, rows
          DO i = 1, row_length
            interp_data(i,j) = rmdi
          END DO
        END DO

        DO l = 1, land_pts
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          interp_data(i,j) = wood_prod_slow_d1(l)
        END DO

! DEPENDS ON: copydiag
        CALL copydiag (STASHwork(si(34,19,im_index)),interp_data,             &
             row_length,rows,0,0,0,0, at_extremity,                           &
             atmos_im,19,34,                                                  &
             icode,cmessage)

        IF (icode  >   0) THEN
           cmessage="Error in copydiag( item 19034)"
           GOTO 9999
        END IF

      END IF     !   sf(34,19)


! 19036:  LIT C FLUX TO FAST POOL KGC/M2/YR
! Carbon flux into the fast-turnover wood product pool from cleared vegetation

      IF (SF(36,19)) THEN
        Do j = 1, rows
          Do i = 1, row_length
            interp_data(i,j) = rmdi
          End Do
        End Do

        Do l = 1, land_pts
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          interp_data(i,j) = wp_fast_in_gb(l)
        End Do

! DEPENDS ON: copydiag
        Call copydiag (STASHwork(si(36,19,im_index)),interp_data,             &
             row_length,rows,0,0,0,0, at_extremity,                           &
             atmos_im,19,36,                                                  &
!#include <argppx/argppx.h>
             icode,cmessage)

        If (icode  >   0) then
           cmessage="Error in copydiag( item 19036)"
           goto 9999
        End if

      END IF     !   sf(36,19)

! 19037:  LIT C FLUX TO MED POOL KGC/M2/YR
! Carbon flux into the medium-turnover wood product pool from cleared vegetation

      IF (SF(37,19)) THEN
        Do j = 1, rows
          Do i = 1, row_length
            interp_data(i,j) = rmdi
          End Do
        End Do

        Do l = 1, land_pts
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          interp_data(i,j) = wp_med_in_gb(l)
        End Do

! DEPENDS ON: copydiag
        Call copydiag (STASHwork(si(37,19,im_index)),interp_data,             &
             row_length,rows,0,0,0,0, at_extremity,                           &
             atmos_im,19,37,                                                  &
!#include <argppx/argppx.h>
             icode,cmessage)

        If (icode  >   0) then
           cmessage="Error in copydiag( item 19037)"
           goto 9999
        End if

      END IF     !   sf(37,19)


! 19038:  LIT C FLUX TO SLOW POOL KGC/M2/YR
! Carbon flux into the slow-turnover wood product pool from cleared vegetation

      IF (SF(38,19)) THEN
        Do j = 1, rows
          Do i = 1, row_length
            interp_data(i,j) = rmdi
          End Do
        End Do

        Do l = 1, land_pts
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          interp_data(i,j) = wp_slow_in_gb(l)
        End Do

! DEPENDS ON: copydiag
        Call copydiag (STASHwork(si(38,19,im_index)),interp_data,             &
             row_length,rows,0,0,0,0, at_extremity,                           &
             atmos_im,19,38,                                                  &
!#include <argppx/argppx.h>
             icode,cmessage)

        If (icode  >   0) then
           cmessage="Error in copydiag( item 19038)"
           goto 9999
        End if

      END IF     !   sf(38,19)


! 19039: FAST WP POOL DECOMP C FLUX KGC/M2/YR
! CO2 flux from decomposition of the fast-turnover wood product pool

      IF (SF(39,19)) THEN
        Do j = 1, rows
          Do i = 1, row_length
            interp_data(i,j) = rmdi
          End Do
        End Do

        Do l = 1, land_pts
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          interp_data(i,j) = wp_fast_out_gb(l)
        End Do

! DEPENDS ON: copydiag
        Call copydiag (STASHwork(si(39,19,im_index)),interp_data,             &
             row_length,rows,0,0,0,0, at_extremity,                           &
             atmos_im,19,39,                                                  &
!#include <argppx/argppx.h>
             icode,cmessage)

        If (icode  >   0) then
           cmessage="Error in copydiag( item 19039)"
           goto 9999
        End if

      END IF     !   sf(39,19)


! 19040: MED WP POOL DECOMP C FLUX KGC/M2/YR
! CO2 flux from decomposition of the medium-turnover wood product pool

      IF (SF(40,19)) THEN
        Do j = 1, rows
          Do i = 1, row_length
            interp_data(i,j) = rmdi
          End Do
        End Do

        Do l = 1, land_pts
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          interp_data(i,j) = wp_med_out_gb(l)
        End Do

! DEPENDS ON: copydiag
        Call copydiag (STASHwork(si(40,19,im_index)),interp_data,             &
             row_length,rows,0,0,0,0, at_extremity,                           &
             atmos_im,19,40,                                                  &
!#include <argppx/argppx.h>
             icode,cmessage)

        If (icode  >   0) then
           cmessage="Error in copydiag( item 19040)"
           goto 9999
        End if

      END IF     !   sf(40,19)


! 19041: SLOW WP POOL DECOMP C FLUX KGC/M2/YR
! CO2 flux from decomposition of the slow-turnover wood product pool

      IF (SF(41,19)) THEN
        Do j = 1, rows
          Do i = 1, row_length
            interp_data(i,j) = rmdi
          End Do
        End Do

        Do l = 1, land_pts
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          interp_data(i,j) = wp_slow_out_gb(l)
        End Do

! DEPENDS ON: copydiag
        Call copydiag (STASHwork(si(41,19,im_index)),interp_data,             &
             row_length,rows,0,0,0,0, at_extremity,                           &
             atmos_im,19,41,                                                  &
!#include <argppx/argppx.h>
             icode,cmessage)

        If (icode  >   0) then
           cmessage="Error in copydiag( item 19041)"
           goto 9999
        End if

      END IF     !   sf(41,19)


! 19042: TOTAL WP POOL DECOMP C FLUX KGC/M2/YR
! Total CO2 flux from decomposition of all three wood product pools combined

      IF (SF(42,19)) THEN
        Do j = 1, rows
          Do i = 1, row_length
            interp_data(i,j) = rmdi
          End Do
        End Do

        Do l = 1, land_pts
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          interp_data(i,j)=wp_fast_out_gb(l)+wp_med_out_gb(l)+wp_slow_out_gb(l)
        End Do

! DEPENDS ON: copydiag
        Call copydiag (STASHwork(si(42,19,im_index)),interp_data,             &
             row_length,rows,0,0,0,0, at_extremity,                           &
             atmos_im,19,42,                                                  &
!#include <argppx/argppx.h>
             icode,cmessage)

        If (icode  >   0) then
           cmessage="Error in copydiag( item 19042)"
           goto 9999
        End if

      END IF     !   sf(42,19)


 9999 CONTINUE
      IF (icode /= 0) THEN

        CALL Ereport(RoutineName,icode,Cmessage)
      ENDIF

      IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE diagnostics_veg
#endif
