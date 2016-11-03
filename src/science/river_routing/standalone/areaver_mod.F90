#if !defined(UM_JULES)
!###############################################################################
!###############################################################################

MODULE areaver_mod

!-----------------------------------------------------------------------------
! Description:
!   Contains regular lat/lon regridding functions for standalone running
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

CONTAINS

!###############################################################################
!###############################################################################
  SUBROUTINE PRE_AREAVER(GAPS_LAMBDA_SRCE,LAMBDA_SRCE                         &
      ,GAPS_PHI_SRCE,PHI_SRCE,CYCLIC_SRCE,LROW_SRCE,WANT,MASK_SRCE            &
      ,GAPS_LAMBDA_TARG,LAMBDA_TARG,GAPS_PHI_TARG,PHI_TARG                    &
      ,CYCLIC_TARG,SPHERICAL                                                  &
      ,MAXL,COUNT_TARG,BASE_TARG,INDEX_SRCE,WEIGHT,ICODE,CMESSAGE)

!     Subroutine PRE_AREAVER
!
!     Calculate weights for area-averaging data on the source grid to
!     data on the target grid.
!
!    J.Gregory  <- programmer of some or all of previous code or changes
!
!    Model            Modification history from model version 3.0:
!   version  Date
!     3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
!                     portability.  Author Tracey Smith.
!     3.3   3.9.93    Correct out-of-bounds errors for IXL and IYT
!     4.1   22.5.96   J.M.Gregory  Add argument SPHERICAL to produce
!                     correct weights for a spherical surface.
!     5.3   22.10.01  Change for south -> north ordering. C.F.Durman.
!     5.5  28.02.03  Add defined RIVERS. C. Bunton
!     6.0  12.09.03  Change DEF from A20 to A26. D. Robinson.
!
!     Standard, paper 4 version 4 (14.12.90)
!
!
!     The source grid and target grid are each specified by a lambda set
!     and a phi set of coordinates delimiting the boxes. These sets are
!     supplied in 1D arrays aa_bb for coordinate aa=LAMBDA,PHI on grid
!     bb=SRCE,TARG.  The number of gaps is specified by GAPS_aa_bb,
!     which is equal to the number of lines IF (CY! IC_bb) and aa is
!     LAMBDA, otherwise one less. (By "gap" we mean the interval
!     spanning a box in one coordinate only. The total number of boxes,
!     or grid points, is the product of the numbers of gaps in the two
!     coordinates.) Whether the axes are cyclic is not known until
!     run-time, so the dimensions of the arrays LAMBDA_bb are not known
!     at compile-time, and they are dimensioned assumed-size. There are
!     no restrictions on the base meridian of lambda, and it does not
!     have to be the same for source and target grids. The lambda
!     coordinates should be in increasing order (running from left to
!     right), the phi increasing (top to bottom). The coordinates must
!     be given in degrees for cyclic axes, because a range of 360 is
!     assumed, or IF (SPHERICAL), when trigonometric functions are
!     used. IF (SPHERICAL), the weights computed are for a spherical
!     surface, assuming that LAMBDA is longitude and PHI latitude.
!     Otherwise, LAMBDA and PHI are treated as Cartesian axes on a
!     plane.
!
!     The array MASK_SRCE is the land/sea mask for the source grid. The
!     logical value indicating points to be used should be supplied in
!     WANT. The first dimension of MASK_SRCE should be supplied in
!     LROW_SRCE. Specifying this separately allows for the possibility
!     of extra rows and columns in MASK_SRCE which are to be ignored.
!
!     The arrays COUNT_TARG and BASE_TARG should be dimensioned in the
!     calling program to the number of boxes in the target array.
!
!     The arrays INDEX_SRCE and WEIGHT are returned in the form which
!     the area-averaging routine DO_AREAVER expects. They are continuous
!     lists comprising consecutive groups of entries. There is a group
!     for each target point, for which the number of entries is spec-
!     ified by COUNT_TARG, and the groups appear in the normal order of
!     grid points. The size required for INDEX_SRCE and WEIGHT depends
!     on how many source boxes go into each target box, on average, and
!     is not known at compile-time. The maximum that could be needed is
!     (GAPS_LAMBDA_SRCE+GAPS_LAMBDA_TARG)*(GAPS_PHI_SRCE+GAPS_PHI_TARG)
!     and the size to which the arrays are actually dimensioned should
!     be supplied in MAXL. The size used is returned in MAXL. It is the
!     responsibility of the calling routine to provide enough space.

!-----------------------------------------------------------------------
    USE c_pi, ONLY :                                                          &
!      imported scalar parameters
       pi_over_180
!-----------------------------------------------------------------------
!
    IMPLICIT NONE

    INTEGER ::                                                                &
       GAPS_LAMBDA_SRCE       &!IN number of lambda gaps in source grid
       ,GAPS_PHI_SRCE          &!IN number of phi gaps in source grid
       ,GAPS_LAMBDA_TARG       &!IN number of lambda gaps in target grid
        ,GAPS_PHI_TARG          &!IN number of phi gaps in target grid
       ,LROW_SRCE              &!IN first dimension of MASK_SRCE
       ,MAXL                   &!INOUT maximum entries in output lists
        ,COUNT_TARG(GAPS_LAMBDA_TARG,GAPS_PHI_TARG)                           &
       !                              !OUT no. of source boxes in target box
        ,BASE_TARG(GAPS_LAMBDA_TARG,GAPS_PHI_TARG)                            &
       !                              !OUT first index in list for target box
        ,INDEX_SRCE(MAXL)       &!OUT list of source box indices
       ,ICODE                   !OUT return code

    LOGICAL ::                                                                &
         CYCLIC_SRCE            &!IN source grid is cyclic
       ,CYCLIC_TARG            &!IN target grid is cyclic
       ,WANT                   &!IN indicator of wanted points in mask
       ,MASK_SRCE(LROW_SRCE,*) &!IN land/sea mask for source grid
       ,SPHERICAL               !IN calculate weights for a sphere

    REAL ::                                                                   &
       LAMBDA_SRCE(*)         &!IN source lambda line coordinates
       ,PHI_SRCE(*)            &!IN source phi line coordinates
       ,LAMBDA_TARG(*)         &!IN target lambda line coordinates
       ,PHI_TARG(*)            &!IN target phi line coordinates
       ,WEIGHT(MAXL)            !OUT list of weights for source boxes

    CHARACTER*(80) :: CMESSAGE !OUT error message

    INTEGER ::                                                                &
       LINES_LAMBDA_SRCE      &! number of source lambda lines
       ,LINES_PHI_SRCE         &! number of source phi lines
       ,LINES_LAMBDA_TARG      &! number of target lambda lines
       ,LINES_PHI_TARG         &! number of target phi lines
       ,COUNT_LAMBDA(GAPS_LAMBDA_TARG)                                        &
       !                              ! number of source lambda gaps per target
       ,COUNT_PHI(GAPS_PHI_TARG)                                              &
       !                              ! number of source phi gaps per target
       ,INDEX_LAMBDA(GAPS_LAMBDA_SRCE+GAPS_LAMBDA_TARG)                       &
       !                              ! source lambda gap indices
       ,INDEX_PHI(GAPS_PHI_SRCE+GAPS_PHI_TARG)                                &
       !                              ! source phi gap indices
       ,IX1,IY1,IX2,IY2        &! working SRCE/TARG LAMBDA/PHI indices
       ,IX1N,IX1W              &! working indices
       ,IXL(GAPS_LAMBDA_TARG+1) &! source line on the left of target line
       ,IX2N                   &! target gap to the right of IX2
       ,IYT(GAPS_PHI_TARG+1)   &! source line above target line
       ,IXP,IYP,IP             &! pointers into lists
       ,IX,IY,I                 ! loop indices

    REAL  ::                                                                  &
       BASLAM                 &! minimum lambda for TEMP coordinates
       ,BTARG                  &! width of target gap
       ,BLO,BHI                &! limits of gap
       ,TEMP_SRCE(GAPS_LAMBDA_SRCE+1)                                         &
       !                              ! adjusted version of LAMBDA_SRCE
       ,TEMP_TARG(GAPS_LAMBDA_TARG+1)                                         &
       !                              ! adjusted version of LAMBDA_TARG
       ,FRAC_LAMBDA(GAPS_LAMBDA_SRCE+GAPS_LAMBDA_TARG)                        &
       !                              ! fractions of target lambda gaps
       ,FRAC_PHI(GAPS_PHI_SRCE+GAPS_PHI_TARG)                                 &
       !                              ! fractions of target phi gaps
       ,SUM                     ! sum of weights

!     1  Set up the lambda coordinates to make them easier to handle.
!
!     1.1  Produce in TEMP_SRCE a monotonically increasing set of angles
!     equivalent to LAMBDA_SRCE i.e. equal under modulo 360.

    IF (CYCLIC_SRCE) THEN
      LINES_LAMBDA_SRCE=GAPS_LAMBDA_SRCE
    ELSE
      LINES_LAMBDA_SRCE=GAPS_LAMBDA_SRCE+1
    ENDIF
    BASLAM=LAMBDA_SRCE(1)
    DO IX1=1,LINES_LAMBDA_SRCE
      IF (LAMBDA_SRCE(IX1) <  BASLAM) THEN
        TEMP_SRCE(IX1)=LAMBDA_SRCE(IX1)+360.
      ELSE
        TEMP_SRCE(IX1)=LAMBDA_SRCE(IX1)
      ENDIF
    ENDDO

!     1.2  Produce in TEMP_TARG a set of angles equivalent to
!     LAMBDA_TARG i.e. equal under modulo 360, but all in the range
!     BASLAM to BASLAM+360, where BASLAM=min(TEMP_LAMBDA).

    IF (CYCLIC_TARG) THEN
      LINES_LAMBDA_TARG=GAPS_LAMBDA_TARG
    ELSE
      LINES_LAMBDA_TARG=GAPS_LAMBDA_TARG+1
    ENDIF
    DO IX2=1,LINES_LAMBDA_TARG
      TEMP_TARG(IX2)=MOD(LAMBDA_TARG(IX2)-BASLAM,360.)
      IF (TEMP_TARG(IX2) <  0.) TEMP_TARG(IX2)=TEMP_TARG(IX2)+360.
      TEMP_TARG(IX2)=TEMP_TARG(IX2)+BASLAM
    ENDDO

!     2  For each target lambda line, find the index of the next source
!     lambda line to the left.

    DO IX2=1,LINES_LAMBDA_TARG
      DO IX1=1,LINES_LAMBDA_SRCE
        IF (TEMP_TARG(IX2) >= TEMP_SRCE(IX1)) IXL(IX2)=IX1
      ENDDO
    ENDDO

!     3  Find which source lambda gaps cover each target gap and the
!     fractions they contribute.
!
!     At this point IXL(target_line) gives the index of the next source
!     lambda line to the left of the target lambda line, wrapping round
!     if the source grid is cyclic. This is also the index of the source
!     gap in which the target line falls. Similarly, the index of the
!     target line is also that of the target gap of which it is the
!     left-hand limit. Therefore also IXL(target_gap+1, wrapping round
!     if reqd.), is the index of the source gap which contains the
!     right-hand limit of the target gap. For each target gap, we loop
!     over all source gaps and find the fraction covered by each. Record
!     the fraction and the source index in cumulative lists. If the
!     source grid is not cyclic, parts of the target gap lying outside
!     the source grid are neglected.

    IXP=0
    DO IX2=1,GAPS_LAMBDA_TARG
      IX=0
      IX2N=MOD(IX2,LINES_LAMBDA_TARG)+1
      BTARG=TEMP_TARG(IX2N)-TEMP_TARG(IX2)
      IF (BTARG <  0.) THEN
        BTARG=BTARG+360.
        IX1N=IXL(IX2N)+LINES_LAMBDA_SRCE
      ELSE
        IX1N=IXL(IX2N)
      ENDIF
      DO IX1W=IXL(IX2),IX1N
        IX1=MOD(IX1W-1,LINES_LAMBDA_SRCE)+1
        IF (CYCLIC_SRCE.OR.IX1 /= LINES_LAMBDA_SRCE) THEN
          IF (IX1W == IXL(IX2)) THEN
            BLO=TEMP_TARG(IX2)
          ELSE
            BLO=TEMP_SRCE(IX1)
          ENDIF
          IF (IX1W == IX1N) THEN
            BHI=TEMP_TARG(IX2N)
          ELSE
            BHI=TEMP_SRCE(MOD(IX1,LINES_LAMBDA_SRCE)+1)
          ENDIF
          IF (BHI <  BLO) BHI=BHI+360.
          IF (ABS(BHI-BLO) >  1E-7*ABS(BTARG)) THEN
            IX=IX+1
            INDEX_LAMBDA(IXP+IX)=IX1
            FRAC_LAMBDA(IXP+IX)=(BHI-BLO)/BTARG
          ENDIF
        ENDIF
      ENDDO
      COUNT_LAMBDA(IX2)=IX
      IXP=IXP+COUNT_LAMBDA(IX2)
    ENDDO

!     4  For each target phi line, find the index of the next source phi
!     line above. Comments as for section 2, without wrap-round. Note
!     that this assumes that the atmosphere gridbox ordering is from
!     south to north; a north to south atmosphere would require LE
!     instead of GE in the following IF test.

    LINES_PHI_SRCE=GAPS_PHI_SRCE+1
    LINES_PHI_TARG=GAPS_PHI_TARG+1
    DO IY2=1,LINES_PHI_TARG
      IYT(IY2)=0
      DO IY1=1,LINES_PHI_SRCE
        IF (PHI_TARG(IY2) >= PHI_SRCE(IY1)) IYT(IY2)=IY1
      ENDDO
    ENDDO

!     5  Find which source phi gaps cover each target gap and the
!     fractions they contribute. Comments as for section 3, without
!     wrap-round.

    IYP=0
    DO IY2=1,GAPS_PHI_TARG
      IY=0
      IF (SPHERICAL) THEN
!     Contain angle between +-90. There is no real area outside
!     these limits on a sphere.
        BTARG=SIN(MAX(MIN(PHI_TARG(IY2+1),90.),-90.)*PI_OVER_180)             &
           -SIN(MAX(MIN(PHI_TARG(IY2),90.),-90.)*PI_OVER_180)
      ELSE
        BTARG=PHI_TARG(IY2+1)-PHI_TARG(IY2)
      ENDIF
      DO IY1=IYT(IY2),IYT(IY2+1)
        IF (.NOT.(IY1 == 0.OR.IY1 == LINES_PHI_SRCE)) THEN
          IF (IY1 == IYT(IY2)) THEN
            BLO=PHI_TARG(IY2)
          ELSE
            BLO=PHI_SRCE(IY1)
          ENDIF
          IF (IY1 == IYT(IY2+1)) THEN
            BHI=PHI_TARG(IY2+1)
          ELSE
            BHI=PHI_SRCE(IY1+1)
          ENDIF
          IF (SPHERICAL) THEN
            BLO=MAX(MIN(BLO,90.),-90.)
            BHI=MAX(MIN(BHI,90.),-90.)
          ENDIF
          IF (ABS(BHI-BLO) >  1E-7*ABS(BTARG)) THEN
            IY=IY+1
            INDEX_PHI(IYP+IY)=IY1
!             Both numerator and denominator in the following are -ve.
            IF (SPHERICAL) THEN
              FRAC_PHI(IYP+IY) = (SIN(BHI*PI_OVER_180)-SIN(BLO*PI_OVER_180))  &
                                  /BTARG
            ELSE
              FRAC_PHI(IYP+IY)=(BHI-BLO)/BTARG
            ENDIF
          ENDIF
        ENDIF
      ENDDO
      COUNT_PHI(IY2)=IY
      IYP=IYP+COUNT_PHI(IY2)
    ENDDO

!     6  For each target box, loop over contributing source boxes and
!     calculate the weights for each one, ignoring land boxes.
!
!     After the first pass for each target box, go back and normalise
!     the weights to compensate for land source boxes and any outside
!     the source grid. Record the source box index and the weight in
!     cumulative lists.

    IP=0
    IYP=0
    DO IY2=1,GAPS_PHI_TARG
      IXP=0
      DO IX2=1,GAPS_LAMBDA_TARG
        I=0
        SUM=0
        DO IY=IYP+1,IYP+COUNT_PHI(IY2)
          DO IX=IXP+1,IXP+COUNT_LAMBDA(IX2)
            IF ( MASK_SRCE(INDEX_LAMBDA(IX),INDEX_PHI(IY)) .EQV. WANT) THEN
              I=I+1
              INDEX_SRCE(IP+I) = INDEX_LAMBDA(IX)                             &
                 + (INDEX_PHI(IY)-1)*GAPS_LAMBDA_SRCE
              WEIGHT(IP+I)=FRAC_LAMBDA(IX)*FRAC_PHI(IY)
              SUM=SUM+WEIGHT(IP+I)
            ENDIF
          ENDDO
        ENDDO
        COUNT_TARG(IX2,IY2)=I
        BASE_TARG(IX2,IY2)=IP
        DO I=1,COUNT_TARG(IX2,IY2)
          WEIGHT(IP+I)=WEIGHT(IP+I)/SUM
        ENDDO
        IP=IP+COUNT_TARG(IX2,IY2)
        IXP=IXP+COUNT_LAMBDA(IX2)
      ENDDO
      IYP=IYP+COUNT_PHI(IY2)
    ENDDO
    MAXL=IP

    ICODE=0
    CMESSAGE=' '
    RETURN

  END SUBROUTINE pre_areaver
!###############################################################################
!###############################################################################
  SUBROUTINE DO_AREAVER( GAPS_LAMBDA_SRCE,GAPS_PHI_SRCE,LROW_SRCE             &
      ,INVERT_SRCE,DATA_SRCE,GAPS_LAMBDA_TARG,GAPS_PHI_TARG,COUNT_TARG        &
      ,BASE_TARG,LROW_TARG,WANT,MASK_TARG,INDEX_SRCE,WEIGHT,ADJUST            &
      ,DATA_TARG,ADJUST_TARG,ICODE,CMESSAGE )

!     Subroutine DO_AREAVER -------------------------------------------
!
!   Purpose:
!
!     Perform area-averaging to transform data from the source grid to
!     the target grid, or adjust the values on the source grid to have
!     the area-averages supplied on the target grid. The latter mode
!     is intended for adjusting values obtained by interpolating from
!     "target" to "source" in order to conserve the area-averages.
!     This mode should be used ONLY if each source box belongs in
!     exactly one target box. ADJUST=0 selects normal area-averaging,
!     ADJUST=1 selects adjustment by addition (use this mode for fields
!     which may have either sign), ADJUST=2 selects adjustment by
!     multiplication (for fields which are positive-definite or
!     negative-definite).
!
!     For two-way conservative coupling, ADJUST=3 makes an adjustment
!     field for fields which may have either sign, ADJUST=4 makes an
!     adjustment field for fields which are positive-definite or
!     negative-definite, ADJUST=5 performs conservative adjustment
!     by addition (use this mode for fields which may have either sign)
!     and ADJUST=6 selects conservative adjustment by multiplication
!     (for fields which are positive-definite or negative-definite).
!
!     The shape of the source and target grids are specified by their
!     dimensions GAPS_aa_bb, which give the number of gaps in the
!     aa=LAMBDA,PHI coordinate in the bb=SRCE,TARG grid. (The product
!     of GAPS_LAMBDA_bb and GAPS_PHI_bb is the number of boxes in the
!     bb grid.)
!
!     The input and output data are supplied as 2D arrays DATA_SRCE and
!     DATA_TARG, whose first dimensions should also be supplied. Speci-
!     fying these sizes separately from the actual dimensions of the
!     grids allows for columns and rows in the arrays to be ignored.
!     A target land/sea mask should be supplied in MASK_TARG, with the
!     value indicating wanted points specified in WANT. Points which
!     are unwanted or which lie outside the source grid are not altered
!     in DATA_TARG. DATA_SRCE can optionally be supplied with its rows
!     in reverse order (i.e. with the first row corresponding to
!     minimum LAMBDA).
!
!     The arrays COUNT_TARG, BASE_TARG, INDEX_SRCE and WEIGHT should be
!     supplied as returned by PRE_AREAVER q.v.
!
!     Programming Standard, paper 4 version 4 (14.12.90)
!
!    Model            Modification history from model version 5.2:
!   Version  date
!    5.3  07.11.01  Extended for use with 2-way conservative coupling
!                   scheme by introducing new output argument
!                   ADJUST_TARG and 4 new operation selectors.
!    5.5  11.04.03  Corrected bug in calculation of adjustment field
!                   for case = 3.
!  5.5  28.02.03  Add defined A20_1A (river routing). C. Bunton
!  6.0  12.09.03  Change DEF from A20 to A26. D. Robinson
!
!   Logical components covered :
!
!   Project task :
!
!   External documentation: Unified Model documentation paper No:
!                           Version:
!
    IMPLICIT NONE

    INTEGER ::                                                                &
       GAPS_LAMBDA_SRCE       &!IN number lambda gaps in source grid
       ,GAPS_PHI_SRCE          &!IN number phi gaps in source grid
       ,LROW_SRCE              &!IN first dimension of source arrays
       ,GAPS_LAMBDA_TARG       &!IN number lambda gaps in target grid
       ,GAPS_PHI_TARG          &!IN number phi gaps in target grid
       ,LROW_TARG              &!IN first dimension of target arrays
       ,COUNT_TARG(GAPS_LAMBDA_TARG,GAPS_PHI_TARG)                            &
       !                              !IN no. of source boxes in target box
       ,BASE_TARG(GAPS_LAMBDA_TARG,GAPS_PHI_TARG)                             &
       !                              !IN first index in list for target box
       ,INDEX_SRCE(*)          &!IN list of source box indices
       ,ADJUST                 &!IN selects normal or adjust mode
       ,ICODE                   !OUT return code

    LOGICAL ::                                                                &
       INVERT_SRCE            &!IN DATA_SRCE rows in reverse order
       ,WANT                   &!IN indicator of wanted points in mask
       ,MASK_TARG(LROW_TARG,*)  !IN land/sea mask for target grid

!     NB alternative intents below apply for normal/adjust mode
      REAL ::                                                                 &
         DATA_SRCE(LROW_SRCE,*) &!IN/INOUT data on source grid
         ,WEIGHT(*)              &!IN list of weights for source boxes
         ,DATA_TARG(LROW_TARG,*) &!INOUT/IN data on target grid
         ,ADJUST_TARG(LROW_TARG,*) !OUT factors by which DATA_SRCE
      !must be adjusted to give DATA_TARG

      CHARACTER ::                                                            &
         CMESSAGE*(*)            !OUT error message

      INTEGER ::                                                              &
         IP                     &! pointer into lists
         ,I                      &! loop index
         ,IX1(GAPS_LAMBDA_SRCE*GAPS_PHI_SRCE)                                 &
!                              ! working SRCE LAMBDA indices
         ,IY1(GAPS_LAMBDA_SRCE*GAPS_PHI_SRCE)                                 &
!                              ! working SRCE PHI indices
         ,IX2,IY2                 ! working TARG LAMBDA/PHI indices

      REAL ::                                                                 &
         TEMP_TARG              &! workspace for area-average
         ,DELTA                  &! additive adjustment
         ,RATIO                   ! multiplicative adjustment

!-------------------------------------------------------------------------------

!     Loop over all target boxes and calculate values as required.
!
!     The weights and source box indices are recorded in continuous
!     lists. COUNT_TARG indicates how many consecutive entries in these
!     lists apply to each target box.
!
      DO IY2=1,GAPS_PHI_TARG
        DO IX2=1,GAPS_LAMBDA_TARG
          IF (MASK_TARG(IX2,IY2).EQV.WANT) THEN
            IF (COUNT_TARG(IX2,IY2) /= 0) THEN
              TEMP_TARG=0.
              DO I=1,COUNT_TARG(IX2,IY2)
                IP=BASE_TARG(IX2,IY2)+I
                IX1(I)=MOD(INDEX_SRCE(IP)-1,GAPS_LAMBDA_SRCE)+1
                IY1(I)=(INDEX_SRCE(IP)-1)/GAPS_LAMBDA_SRCE+1
                IF (INVERT_SRCE) IY1(I)=GAPS_PHI_SRCE-IY1(I)+1
                TEMP_TARG=TEMP_TARG+WEIGHT(IP)*DATA_SRCE(IX1(I),IY1(I))
              ENDDO
            ELSE
              IF(ADJUST == 5)THEN
                TEMP_TARG=0.0
              ELSEIF(ADJUST == 6)THEN
                TEMP_TARG=1.0
              ELSE
                TEMP_TARG=DATA_TARG(IX2,IY2)
              ENDIF
            ENDIF
            IF (ADJUST == 0) THEN
              DATA_TARG(IX2,IY2)=TEMP_TARG
            ELSEIF (ADJUST == 1) THEN
              DELTA=DATA_TARG(IX2,IY2)-TEMP_TARG
              DO I=1,COUNT_TARG(IX2,IY2)
                DATA_SRCE(IX1(I),IY1(I))=DATA_SRCE(IX1(I),IY1(I))+DELTA
              ENDDO
            ELSEIF (ADJUST == 2.AND.TEMP_TARG /= 0.) THEN
              RATIO=DATA_TARG(IX2,IY2)/TEMP_TARG
              DO I=1,COUNT_TARG(IX2,IY2)
                DATA_SRCE(IX1(I),IY1(I))=DATA_SRCE(IX1(I),IY1(I))*RATIO
              ENDDO
            ELSEIF (ADJUST == 3) THEN
              ADJUST_TARG(IX2,IY2)=TEMP_TARG-DATA_TARG(IX2,IY2)
            ELSEIF (ADJUST == 4) THEN
              IF (TEMP_TARG == 0) THEN
                ADJUST_TARG(IX2,IY2)=1.0
              ELSE
                ADJUST_TARG(IX2,IY2)=DATA_TARG(IX2,IY2)/TEMP_TARG
              ENDIF
            ELSEIF (ADJUST == 5) THEN
              DATA_TARG(IX2,IY2)=DATA_TARG(IX2,IY2)+TEMP_TARG
            ELSEIF (ADJUST == 6) THEN
              DATA_TARG(IX2,IY2)=DATA_TARG(IX2,IY2)*TEMP_TARG
            ENDIF
          ENDIF
        ENDDO
      ENDDO

      ICODE=0
      CMESSAGE=' '
      RETURN

    END SUBROUTINE do_areaver
!###############################################################################
!###############################################################################
  END MODULE areaver_mod
#endif
