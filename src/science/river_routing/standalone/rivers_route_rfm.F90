! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    MODULE jules_rivers_rfm_mod------------------------------------------------
!
! Description:
!     Driver and science routines for calculating river flow routing
!     using the RFM kinematic wave model
!     see Bell et al. 2007 Hydrol. Earth Sys. Sci. 11. 532-549
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

MODULE jules_rivers_rfm_mod

CONTAINS

!###############################################################################
! subroutine rivers_drive_rfm
! Driver routine for runoff routing by the RFM (kinematic wave) model.

SUBROUTINE rivers_drive_rfm( sub_runoffin, surf_runoffin,                     &
                             runoff_out, rivflow )

!-------------------------------------------------------------------------------
!
! Description:
!   Perform the routing of surface and sub-surface runoff using RFM
!
!   This routine regrids the total surface runoff to the river grid and passes
!   it to the RFM routines to be routed.
!
!-------------------------------------------------------------------------------
! Modules used:

  USE jules_rivers_mod, ONLY :                                                &
!  imported scalars with intent(in)
    np_rivers,nx_rivers,ny_rivers                                             &
    ,nx_grid, ny_grid                                                         &
    ,rivers_dlat,rivers_dlon                                                  &
    ,rivers_regrid, rivers_reglatlon                                          &
!  imported arrays with intent(in)
    ,rivers_lat2d, rivers_lon2d, rivers_boxareas

  USE rivers_utils, ONLY :                                                    &
!  imported procedures
     rivers_route_regrid, rivers_route_regrid_invert

  USE jules_print_mgr, ONLY :                                                 &
    jules_message,                                                            &
    jules_print

  IMPLICIT NONE

!-------------------------------------------------------------------------------

! Array arguments with intent(in)
  REAL, INTENT(IN) :: sub_runoffin(nx_grid, ny_grid)
                               ! Average rate of sub surface runoff since
                               ! last rivers call on model grid in kg m-2 s-1
  REAL, INTENT(IN) :: surf_runoffin(nx_grid, ny_grid)
                               ! Average rate of surface runoff since last
                               ! rivers call on model grid in kg m-2 s-1

! Array arguments with intent(out)
  REAL, INTENT(OUT) :: runoff_out(nx_grid, ny_grid)
                               ! Total runoff diagnostic on model grid
                               ! in kg m-2 s-1
  REAL, INTENT(OUT) :: rivflow(nx_grid, ny_grid)
                               ! River flow diagnostic on model grid
                               ! in kg m-2 s-1

! Local scalar variables

  INTEGER :: i, j              !  work
  INTEGER :: l                 !  loop counter (land point)
  INTEGER :: ip                !  loop counter
  INTEGER :: ix                !  work
  INTEGER :: iy                !  work
  INTEGER :: maxlin            !  an array size

! Local array variables.

  INTEGER :: x(np_rivers)      !  x index on the rivers grid of
                               !  active rivers points
  INTEGER :: xnext(np_rivers)  !  x location of the next downstream point
                               !     <0 indicates next point is off grid
  INTEGER :: y(np_rivers)      !  y index on the rivers grid of
                               !  active rivers points
  INTEGER :: ynext(np_rivers)  !  y location of the next downstream point
                               !  If xnext(i)<0, ynext(i) gives the direction
                               !  of the next point, in terms of index in
                               !  flow_dir_set.

  REAL :: outflow(nx_rivers,ny_rivers)
                               !  rate of channel flow leaving gridbox
                               !  (kg m-2 s-1)
  REAL :: baseflow(nx_rivers,ny_rivers)
                               !  rate of channel base flow leaving
                               !  gridbox (kg m-2 s-1)
  REAL :: surf_runoff_rivers(nx_rivers, ny_rivers)
                               !  average rate of surface runoff since last
                               !  rivers call (kg m-2 s-1) on rivers grid.
  REAL :: sub_runoff_rivers(nx_rivers, ny_rivers)
                               !  average rate of sub-surface runoff since last
                               !  rivers call (kg m-2 s-1) on rivers grid.
  REAL :: tot_runoff_rivers(nx_rivers, ny_rivers)
                               ! runoff diagnostic on rivers grid (kg m-2 s-1)

!-------------------------------------------------------------------------------
! Initialisation
!-------------------------------------------------------------------------------

  outflow(:,:) = 0.0
  baseflow(:,:) = 0.0
  rivflow(:,:) = 0.0
  runoff_out(:,:) = 0.0

!-------------------------------------------------------------------------------
! Regrid surface and subsurface runoff from land grid to rivers grid
!-------------------------------------------------------------------------------

  sub_runoff_rivers(:,:) = 0.0
  surf_runoff_rivers(:,:) = 0.0

  IF ( rivers_regrid ) THEN

!   Calculate an array size used in regridding.
    maxlin = ( nx_grid + nx_rivers )*( ny_grid + ny_rivers )
    CALL rivers_route_regrid( maxlin,sub_runoffin,sub_runoff_rivers )
    CALL rivers_route_regrid( maxlin,surf_runoffin,surf_runoff_rivers )

  ELSE

!   "Main" and rivers grids are identical (including points in same order).
!   No need to regrid.

    sub_runoff_rivers(:,:) = sub_runoffin(:,:)
    surf_runoff_rivers(:,:) = surf_runoffin(:,:)

  END IF

!-------------------------------------------------------------------------------
! Call the routing routine.
!-------------------------------------------------------------------------------

   CALL rivers_rfm( surf_runoff_rivers, sub_runoff_rivers, rivers_boxareas,   &
                   outflow, baseflow)
!   CALL rivers_rfm_vector( surf_runoff_rivers, sub_runoff_rivers, &
!                 outflow, baseflow)

!-------------------------------------------------------------------------------
!   Regrid from rivers to land grid
!-------------------------------------------------------------------------------

! Sum the surface and subsurface runoffs
  tot_runoff_rivers(:,:) = surf_runoff_rivers(:,:) + sub_runoff_rivers(:,:)

  IF ( rivers_regrid ) THEN

    CALL rivers_route_regrid_invert( maxlin,outflow,rivflow )

    ! Sanity check regridding routines - expect input = output
    CALL rivers_route_regrid_invert( maxlin, tot_runoff_rivers, runoff_out )

  ELSE

    rivflow(:,:) = outflow(:,:)
    runoff_out(:,:) = tot_runoff_rivers(:,:)

  END IF

END SUBROUTINE rivers_drive_rfm

!###############################################################################
! subroutine rivers_rfm
!
!-----------------------------------------------------------------------------
! Description:
!   Perform the routing of surface and sub-surface runoff for Regional
!   model.
!   Calculates river outflow (kg m-2 s-1) and baseflow (kg m-2 s-1) for the
!   RFM kinematic wave river routing model.
!
! Method:
!   See Bell et al. 2007 Hydrol. Earth Sys. Sci. 11. 532-549
!   This subroutine routes surface and subsurface runoff
!   for the RCM. The routing procedure is performed on the
!   RCM grid, so no regridding is needed. The routing
!   procedure is currently based on the Kinematic Wave
!   routing model used at CEH Wallingford.
!
! Author: V.A.Bell, CEH Wallingford, 21.08.03
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!
! Code History:
!    Modified from river2a.f by sjd 13/05/05
!    Modified from UM routine riv_rout-river2a.F90 by hl 13/04/14
!
!  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 5.5:
! VERSION  DATE
!   6.0   12/09/03  Change DEF from A20 to A26. D. Robinson
!   6.0   12.09.03  Routing code added. V.A.Bell
!   x.x   02/05/12  Additional implementation and developments S. Dadson
!   x.x   06/01/15  Formal implementation within JULES code base H. Lewis
!
! NOTE ON UNITS:
!   This routine, based on Simon Dadson's work, includes a modification to the
!   'standard' RFM routines which assume a regular x/y-grid based
!   implementation to account for potential variable grid box areas (e.g. from
!   lat/lon grid).
!   Stores are calculated in units of m x m2 rather than mm, and flows are
!   initially calculated in units of m3/s.
!   For consistency with other routines (for now), the output is converted
!   again to a flux density kg/m2/s

  SUBROUTINE rivers_rfm( rfm_surf_ingrid, rfm_sub_surf_ingrid,                &
                         rfm_rivers_area, outflow, baseflow )

  USE jules_rivers_mod, ONLY :                                                &
!  imported scalar parameters
       rivers_timestep,rivers_type, nx_rivers, ny_rivers, rivers_dx           &
       ,rivers_dlat,rivers_count, rivers_step, rivers_first                   &
!  imported scalars with intent (inout)
!       ,l_rfmFirst, l_rfm_overbank &
!  imported arrays with intent(in)
       ,rivers_lat2d, rivers_lon2d, rfm_flowobs1                              &
!  imported arrays with intent(inout)
       ,rfm_inext, rfm_jnext, rfm_iarea, rfm_land                             &
       ,rfm_substore, rfm_surfstore, rfm_flowin, rfm_bflowin                  &
       ,rfm_rivflow, rfm_baseflow

  USE timestep_mod, ONLY : timestep

  USE planet_constants_mod, ONLY : planet_radius  ! the Earth's radius (m)

  USE c_pi, ONLY : pi_over_180

  USE c_mm_m, ONLY : mm_2_m

  USE jules_riversparm, ONLY : cland, criver, cbland, cbriver,                &
                               runoff_factor, retl, retr, slfac, a_thresh

  USE jules_print_mgr, ONLY :                                                 &
     jules_message,                                                           &
     jules_print

!-------------------------------------------------------------------------------

  IMPLICIT NONE

! IN Arguments

! Array arguments with intent(in)
  REAL, INTENT(IN) :: rfm_surf_ingrid(nx_rivers, ny_rivers)
         !  average rate of surface runoff since last rivers call (kg m-2 s-1)
  REAL, INTENT(IN) :: rfm_sub_surf_ingrid(nx_rivers, ny_rivers)
         !  average rate of sub-surface runoff since last call (kg m-2 s-1)
  REAL, INTENT(IN) :: rfm_rivers_area(nx_rivers, ny_rivers)
         ! grid box area in m2

  REAL, INTENT(OUT) :: outflow(nx_rivers,ny_rivers)
         !  rate of channel flow leaving gridbox (kg m-2 s-1)
  REAL, INTENT(OUT) :: baseflow(nx_rivers,ny_rivers)
         !  rate of channel base flow leaving gridbox (kg m-2 s-1)

! internal variables
  INTEGER ::                                                                  &
       landtype                                                               &
         !local for land type
       ,in,jn                                                                 &
         !local co-ords of downstream point
       ,i,j
         !co-ordinate counters in do loops

  REAL ::                                                                     &
     landtheta, rivertheta                                                    &
         !surface wave speeds
     ,sublandtheta, subrivertheta                                             &
         !sub-surface wave speeds
     ,returnflow                                                              &
         !returnflow
     ,flowobs1_m3s                                                            &
         !intial river flow [m3/s]
     ,dt                                                                      &
         !rfm model t/s
     ,dx
         !grid spacing in metres

  REAL ::                                                                     &
     substore_n(nx_rivers,ny_rivers)                                          &
         ! subsurface store at next timestep(m3)
     ,surfstore_n(nx_rivers,ny_rivers)                                        &
         ! surface store at next timestep (m3)
     ,flowin_n(nx_rivers,ny_rivers)                                           &
         ! surface lateral inflow next time (m3)
     ,bflowin_n(nx_rivers,ny_rivers)                                          &
         ! sub-surface lat inflow next time (m3)
     ,rarea(nx_rivers,ny_rivers)                                              &
         ! accumulated area (iarea)
     ,surf_roff(nx_rivers,ny_rivers)                                          &
         !INTERNAL surf_runoff
     ,sub_surf_roff(nx_rivers,ny_rivers)
         !INTERNAL sub_surf_runoff

  INTEGER :: nland, nsea, nriv

!-------------------------------------------------------------------------------
! Set up rivers parameters
!-------------------------------------------------------------------------------

  dt = rivers_timestep * timestep                     ! rivers model timestep(s)

  IF (rivers_dx <= 0) THEN
    dx = planet_radius * (ABS(rivers_dlat)*pi_over_180)  ! horiz. gridsize (m)
  ELSE
    dx = rivers_dx
  END IF

  rivertheta = criver * dt / dx           ! surface wave speeds
  landtheta = cland * dt / dx
  sublandtheta = cbland * dt / dx
  subrivertheta = cbriver * dt / dx

  ! Check condition for numerical stability
  IF (landtheta > 1.0 .OR. sublandtheta > 1.0 .OR.                            &
      rivertheta > 1.0 .OR. subrivertheta > 1.0) THEN
    WRITE(jules_message,*) 'ERROR: ' //                                       &
       'Unstable theta reduce wavespeed in RFM routing. Setting to zero.'
    CALL jules_print('rivers_rfm',jules_message)
    rivertheta = 0.0
    landtheta = 0.0
    sublandtheta = 0.0
    subrivertheta = 0.0
  END IF

!-------------------------------------------------------------------------------
! Initialise variables at first timestep
!-------------------------------------------------------------------------------
!! HL: TO FIX - how much of this could/should be done at initialisation steps??

  nland = 0
  nsea = 0
  nriv = 0

  IF (rivers_first) THEN

    ! From the cumulative catchment areas dataset, determine which
    ! squares are land or river
    !------------------------------------------------------------------
    ! rfm_land(i,j) is set at the first entry, and the values must be
    ! remembered for subsequent timesteps
    !-------------------------------------------------------------------

    rivers_first = .FALSE.

    DO j = 1, ny_rivers
      DO i = 1, nx_rivers

        IF (rfm_iarea(i,j) < 0) THEN
          rfm_land(i,j) = 0     !Sea
          nsea = nsea+1
        ELSE IF (rfm_iarea(i,j) > a_thresh) THEN
          rfm_land(i,j) = 1     !river
          nriv = nriv + 1
        ELSE
          rfm_land(i,j) = 2     !land
          nland = nland+1
        END IF
        rarea(i,j) = REAL(rfm_iarea(i,j)) + 1.0
        ! include current point, so add 1

        ! set to sea if top or left edge drains outside model domain
        IF (rfm_inext(i,j) < 1 .OR. rfm_jnext(i,j) < 1) THEN
          rfm_land(i,j) = 0  !Sea
        END IF

        ! set to sea if bottom or right edge drains outside model domain
        IF (rfm_inext(i,j) > nx_rivers .OR. rfm_jnext(i,j) > ny_rivers) THEN
          rfm_land(i,j) = 0  !Sea
        END IF

      END DO !end of rivers loop, i
    END DO !end of rivers loop, j

    ! Initialise surface and sub-surface stores to zero
    rfm_surfstore(:,:) = 0.0
    rfm_substore(:,:) = 0.0

    ! Initialise surface and sub-surface inflows
    rfm_flowin(:,:) = 0.0
    rfm_bflowin(:,:) = 0.0
    rfm_rivflow(:,:) = 0.0

!! TO FIX - update to enable initialisation of storage based on flowobs1
!!          N.B. need to check correct implementation of units
    DO j=1,ny_rivers
      DO i=1,nx_rivers

! Initialise surface and sub-surface stores using flow observations if available
        IF ( rfm_flowobs1(i,j) > 0.0 ) THEN
          flowobs1_m3s = rfm_flowobs1(i,j) * rfm_rivers_area(i,j) / mm_2_m
          rfm_surfstore(i,j) = flowobs1_m3s * dt / rivertheta
          rfm_substore(i,j) = flowobs1_m3s * dt / subrivertheta
        END IF
      END DO
    END DO

  END IF   ! end riversFirst
!! HL: END TO FIX

!-------------------------------------------------------------------------------
! Processing for each timestep
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Convert runoff from kg/sec to mm (over the model timestep)
! assume dx in metres and dt in seconds
!-------------------------------------------------------------------------------

  DO j = 1,ny_rivers
    DO i = 1,nx_rivers

      IF ( rfm_surf_ingrid(i,j) >= 0.0 .AND.                                  &
              rfm_sub_surf_ingrid(i,j) >= 0.0 ) THEN

        ! multiply by grid box area (in m) and divide by 1000 to convert
        ! depth in mm to depth in m (sjd 02/05/12)

        surf_roff(i,j) = runoff_factor * rfm_surf_ingrid(i,j) *               &
                            dt * rfm_rivers_area(i,j) / mm_2_m
        sub_surf_roff(i,j) = runoff_factor * rfm_sub_surf_ingrid(i,j) *       &
                            dt * rfm_rivers_area(i,j) / mm_2_m

      ELSE

        ! ignore no data (-1.00e+20) and other negative values
        surf_roff(i,j) = 0.0
        sub_surf_roff(i,j) = 0.0

      END IF

    END DO
  END DO

!-------------------------------------------------------------------------------
! Initialise accumulated inflows and stores for the next timestep
!-------------------------------------------------------------------------------

  flowin_n(:,:) = 0.0
  bflowin_n(:,:) = 0.0
  surfstore_n(:,:) = 0.0
  substore_n(:,:) = 0.0

!-------------------------------------------------------------------------------
! Rivers runoff using simple kinematic wave model
!-------------------------------------------------------------------------------

  DO j = 1,ny_rivers
    DO i = 1,nx_rivers

      in = rfm_inext(i,j)
      jn = rfm_jnext(i,j)
      landtype = rfm_land(i,j)

      IF (landtype == 2) THEN  !land

        ! land surface
        surfstore_n(i,j) = (1.0 - landtheta) * rfm_surfstore(i,j) +           &
                               rfm_flowin(i,j) + surf_roff(i,j)

        ! land subsurface
        substore_n(i,j) = (1.0 - sublandtheta) * rfm_substore(i,j) +          &
                               rfm_bflowin(i,j) + sub_surf_roff(i,j)

        ! return flow
        IF (retl > 0) THEN ! changed to allow for -ve retr/l (sjd: 4/6/09)
          returnflow = AMAX1( ABS( substore_n(i,j) * retl ), 0.0 )
        ELSE
          returnflow = -1.0 * AMAX1( ABS( surfstore_n(i,j) * retl ), 0.0 )
        END IF

        substore_n(i,j) = substore_n(i,j) - returnflow
        surfstore_n(i,j) = surfstore_n(i,j) + returnflow

        flowin_n(in,jn) = flowin_n(in,jn) + landtheta * rfm_surfstore(i,j)
        bflowin_n(in,jn) = bflowin_n(in,jn) + sublandtheta * rfm_substore(i,j)

      ELSE IF (landtype == 1) THEN  !river

        ! river subsurface
        substore_n(i,j) = (1.0 - subrivertheta) * rfm_substore(i,j) +         &
                            rfm_bflowin(i,j) + sub_surf_roff(i,j)

        ! river surface
        surfstore_n(i,j) = (1.0 - rivertheta) * rfm_surfstore(i,j) +          &
                              rfm_flowin(i,j) + surf_roff(i,j)

        ! return flow
        IF (retr > 0) THEN ! changed to allow for -ve retr/l (sjd: 4/6/09)
          returnflow = AMAX1( ABS( substore_n(i,j) * retr ), 0.0 )
        ELSE
          returnflow = -1.0 * AMAX1( ABS( surfstore_n(i,j) * retr ), 0.0 )
        END IF
        substore_n(i,j) = substore_n(i,j) - returnflow
        surfstore_n(i,j) = surfstore_n(i,j) + returnflow

        flowin_n(in,jn) = flowin_n(in,jn) + rivertheta * rfm_surfstore(i,j)
        bflowin_n(in,jn) = bflowin_n(in,jn) + subrivertheta * rfm_substore(i,j)

      END IF ! land or river

    END DO !end of rivers loop, i
  END DO !end of rivers loop, j

!-------------------------------------------------------------------------------
! Housekeeping for next timestep
!-------------------------------------------------------------------------------

  DO j = 1, ny_rivers
    DO i = 1, nx_rivers

      ! keep inflows for next timestep
      rfm_flowin(i,j) = flowin_n(i,j)
      rfm_bflowin(i,j) = bflowin_n(i,j)

      ! calculate flow in rivers (m3/s)
      IF (rfm_land(i,j) == 1) THEN

        rfm_rivflow(i,j) = rivertheta / dt * rfm_surfstore(i,j)
        rfm_baseflow(i,j) = subrivertheta / dt * rfm_substore(i,j)

        in = rfm_inext(i,j)
        jn = rfm_jnext(i,j)

        ! calculate flow into the sea if the next point is sea
        !(assume it's equal to river flow in the adjacent land pt)
        IF (rfm_land(in,jn) == 0) THEN
          rfm_rivflow(in,jn) = rfm_rivflow(i,j)
        END IF

      END IF

      !calculate flow for land points (m3/s)
      IF (rfm_land(i,j) == 2) THEN

        rfm_rivflow(i,j) = landtheta / dt * rfm_surfstore(i,j)
        rfm_baseflow(i,j) = sublandtheta / dt * rfm_substore(i,j)

        in = rfm_inext(i,j)
        jn = rfm_jnext(i,j)

        ! calculate flow into the sea if the next point is sea
        !(assume it's equal to river flow in the adjacent land pt)
        IF (rfm_land(in,jn) == 0) THEN
          rfm_rivflow(in,jn) = rfm_rivflow(i,j)
        END IF

      END IF

      ! keep rivers stores for next timestep  (m3)
      rfm_surfstore(i,j) = surfstore_n(i,j)
      rfm_substore(i,j) = substore_n(i,j)

    END DO !end rivers loop, i
  END DO !end rivers loop, j

!-------------------------------------------------------------------------------
! Return flows in flux density units kg/m2/s
!-------------------------------------------------------------------------------
!! HL: N.B. Add option to select what output units required for river flow

  outflow(:,:) = rfm_rivflow(:,:) * mm_2_m / rfm_rivers_area(:,:)
  baseflow(:,:) = rfm_baseflow(:,:) * mm_2_m / rfm_rivers_area(:,:)

END SUBROUTINE rivers_rfm

!###############################################################################
!###############################################################################

!###############################################################################
! subroutine rivers_rfm_vector
!
!-----------------------------------------------------------------------------
! Description:
!   Perform the routing of surface and sub-surface runoff for Regional
!   model.
!   Calculates river outflow (kg m-2 s-1) and baseflow (kg m-2 s-1) for the
!   RFM kinematic wave river routing model.
!
! Method:
!   See Bell et al. 2007 Hydrol. Earth Sys. Sci. 11. 532-549
!   This subroutine routes surface and subsurface runoff
!   for the RCM. The routing procedure is performed on the
!   RCM grid, so no regridding is needed. The routing
!   procedure is currently based on the Kinematic Wave
!   routing model used at CEH Wallingford.
!
! Author: V.A.Bell, CEH Wallingford, 21.08.03
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!
! Code History:
!    Modified from river2a.f by sjd 13/05/05
!    Modified from UM routine riv_rout-river2a.F90 by hl 13/04/14
!    Updated to route runoff on riv_pts vector only by hl 24/04/14
!
!  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 5.5:
! VERSION  DATE
!   6.0   12/09/03  Change DEF from A20 to A26. D. Robinson
!   6.0   12.09.03  Routing code added. V.A.Bell
!   x.x   02/05/12  Additional implementation and developments S. Dadson
!   x.x   06/01/15  Formal implementation within JULES code base H. Lewis
!
! NOTE ON UNITS:
!   This routine, based on Simon Dadson's work, includes a modification to the
!   'standard' RFM routines which assume a regular x/y-grid based
!   implementation to account for potential variable grid box areas (e.g. from
!   lat/lon grid).
!   Stores are calculated in units of m x m2 rather than mm, and flows are
!   initially calculated in units of m3/s.
!   For consistency with other routines (for now), the output is converted
!   again to a flux density kg/m2/s

  SUBROUTINE rivers_rfm_vector( rfm_surf_ingrid, rfm_sub_surf_ingrid,         &
                                outflow, baseflow )

  USE jules_rivers_mod, ONLY :                                                &
!  imported scalar parameters
       rivers_timestep,rivers_type, nx_rivers, ny_rivers, np_rivers           &
       ,rivers_dlat,rivers_count, rivers_step, rivers_first, rivers_dx        &
!  imported scalars with intent (inout)
!       ,l_rfmFirst, l_rfm_overbank &
!  imported arrays with intent(in)
       ,rivers_lat2d, rivers_lon2d, rfm_flowobs1_rp                           &
!  imported arrays with intent(inout)
       ,rivers_next_rp, rivers_index_rp                                       &
       ,rfm_iarea_rp, rfm_land_rp                                             &
       ,rfm_substore_rp, rfm_surfstore_rp, rfm_flowin_rp, rfm_bflowin_rp      &
       ,rfm_rivflow_rp, rfm_baseflow_rp, rivers_boxareas_rp

  USE timestep_mod, ONLY : timestep

  USE planet_constants_mod, ONLY : planet_radius  ! the Earth's radius (m)

  USE c_pi, ONLY : pi_over_180

  USE c_mm_m, ONLY : mm_2_m

  USE jules_riversparm, ONLY : cland, criver, cbland, cbriver,                &
                               runoff_factor, retl, retr, slfac, a_thresh

  USE rivers_utils, ONLY :                                                    &
!  imported procedures
     rivers_get_xy_pos

  USE jules_print_mgr, ONLY :                                                 &
     jules_message,                                                           &
     jules_print

!-------------------------------------------------------------------------------

  IMPLICIT NONE

! IN Arguments

! Array arguments with intent(in)
  REAL, INTENT(IN) :: rfm_surf_ingrid(nx_rivers, ny_rivers)
         !  average rate of surface runoff since last rivers call (kg m-2 s-1)
  REAL, INTENT(IN) :: rfm_sub_surf_ingrid(nx_rivers, ny_rivers)
         !  average rate of sub-surface runoff since last call (kg m-2 s-1)
!  REAL, INTENT(IN) :: rfm_rivers_area(nx_rivers, ny_rivers)
!         ! grid box area in m2

  REAL, INTENT(OUT) :: outflow(nx_rivers,ny_rivers)
         !  rate of channel flow leaving gridbox (kg m-2 s-1)
  REAL, INTENT(OUT) :: baseflow(nx_rivers,ny_rivers)
         !  rate of channel base flow leaving gridbox (kg m-2 s-1)

! internal variables
  INTEGER ::                                                                  &
       landtype                                                               &
         !local for land type
       ,imin,imax                                                             &
         !local for max/min vector indices
       ,in,jn,rn                                                              &
         !local co-ords of downstream point
       ,i,j,ip,rp
         !co-ordinate counters in do loops

  REAL ::                                                                     &
     landtheta, rivertheta                                                    &
         !surface wave speeds
     ,sublandtheta, subrivertheta                                             &
         !sub-surface wave speeds
     ,returnflow                                                              &
         !returnflow
     ,flowobs1_m3s                                                            &
         !intial river flow [m3/s]
     ,dt                                                                      &
         !rfm model t/s
     ,dx
         !grid spacing in metres

  REAL ::                                                                     &
     substore_n(np_rivers)                                                    &
         ! subsurface store at next timestep(m3)
     ,surfstore_n(np_rivers)                                                  &
         ! surface store at next timestep (m3)
     ,flowin_n(nx_rivers*ny_rivers)                                           &
         ! surface lateral inflow next time (m3)
     ,bflowin_n(nx_rivers*ny_rivers)                                          &
         ! sub-surface lat inflow next time (m3)
     ,rarea(np_rivers)                                                        &
         ! accumulated area (iarea)
     ,surf_roff(np_rivers)                                                    &
         !INTERNAL surf_runoff
     ,sub_surf_roff(np_rivers)
         !INTERNAL sub_surf_runoff

!  REAL, ALLOCATABLE ::                                                        &
!     substore_n(:)                                                     &
!         ! subsurface store at next timestep(m3)
!     ,surfstore_n(:)                                                   &
!         ! surface store at next timestep (m3)
!     ,flowin_n(:)                                                      &
!         ! surface lateral inflow next time (m3)
!     ,bflowin_n(:)                                                     &
         ! sub-surface lat inflow next time (m3)
!     ,rarea(np_rivers)                                                        &
         ! accumulated area (iarea)
!     ,surf_roff(:)                                                     &
         !INTERNAL surf_runoff
!     ,sub_surf_roff(:)
         !INTERNAL sub_surf_runoff

  INTEGER :: nland, nsea, nriv

!-------------------------------------------------------------------------------
! Set up rivers parameters
!-------------------------------------------------------------------------------

  dt = rivers_timestep * timestep                     ! rivers model timestep(s)

  IF (rivers_dx <= 0) THEN
    dx = planet_radius * (ABS(rivers_dlat)*pi_over_180)  ! horiz. gridsize (m)
  ELSE
    dx = rivers_dx
  END IF

  rivertheta = criver * dt / dx           ! surface wave speeds
  landtheta = cland * dt / dx
  sublandtheta = cbland * dt / dx
  subrivertheta = cbriver * dt / dx

  ! Check condition for numerical stability
  IF (landtheta > 1.0 .OR. sublandtheta > 1.0 .OR.                            &
      rivertheta > 1.0 .OR. subrivertheta > 1.0) THEN
    WRITE(jules_message,*) 'ERROR: ' //                                       &
       'Unstable theta reduce wavespeed in RFM routing. Setting to zero.'
    CALL jules_print('rivers_rfm',jules_message)
    rivertheta = 0.0
    landtheta = 0.0
    sublandtheta = 0.0
    subrivertheta = 0.0
  END IF

!-------------------------------------------------------------------------------
! Initialise variables at first timestep
!-------------------------------------------------------------------------------
!! HL: TO FIX - how much of this could/should be done at initialisation steps??

  nland = 0
  nsea = 0
  nriv = 0

  IF (rivers_first) THEN

    ! From the cumulative catchment areas dataset, determine which
    ! squares are land or river
    !------------------------------------------------------------------
    ! rfm_land(i,j) is set at the first entry, and the values must be
    ! remembered for subsequent timesteps
    !-------------------------------------------------------------------

    rivers_first = .FALSE.
    imin = MINVAL(rivers_index_rp)
    imax = MAXVAL(rivers_index_rp)

    DO ip = 1, np_rivers

      rp = rivers_index_rp(ip)

      IF (rfm_iarea_rp(ip) < 0) THEN
        rfm_land_rp(rp) = 0     !Sea
        nsea = nsea+1
      ELSE IF (rfm_iarea_rp(ip) > a_thresh) THEN
        rfm_land_rp(rp) = 1     !river
        nriv = nriv + 1
      ELSE
        rfm_land_rp(rp) = 2     !land
        nland = nland+1
      END IF
!      rarea(ip) = REAL(rfm_iarea_rp(ip)) + 1.0
      ! include current point, so add 1

      ! set to sea if top or left edge drains outside model domain
      IF (rivers_next_rp(ip) < imin .OR. rivers_next_rp(ip) > imax) THEN
        rfm_land_rp(rp) = 0  !Sea
      END IF

! Initialise surface and sub-surface stores using flow observations if available
      IF ( rfm_flowobs1_rp(ip) > 0.0 ) THEN
        flowobs1_m3s = rfm_flowobs1_rp(ip) * rivers_boxareas_rp(ip) / mm_2_m
        rfm_surfstore_rp(ip) = flowobs1_m3s * dt / rivertheta
        rfm_substore_rp(ip) = flowobs1_m3s * dt / subrivertheta
      END IF
    END DO

!! HL - need to avoid this if initialising with prognostics
    ! Initialise surface and sub-surface stores to zero
    rfm_surfstore_rp(:) = 0.0  !  HL only if not prognostics
    rfm_substore_rp(:) = 0.0   !  HL only if not prognostics

    ! Initialise surface and sub-surface inflows
    rfm_flowin_rp(:) = 0.0       !  HL only if not prognostics
    rfm_bflowin_rp(:) = 0.0      !  HL only if not prognostics
    rfm_rivflow_rp(:) = 0.0
    rfm_baseflow_rp(:) = 0.0

  END IF   ! end riversFirst
!! HL: END TO FIX

!-------------------------------------------------------------------------------
! Processing for each timestep
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Convert runoff from kg/sec to mm (over the model timestep)
! assume dx in metres and dt in seconds
!-------------------------------------------------------------------------------

  DO ip = 1,np_rivers

    CALL rivers_get_xy_pos(rivers_index_rp(ip),nx_rivers,ny_rivers,i,j)

    IF ( rfm_surf_ingrid(i,j) >= 0.0 .AND.                                    &
              rfm_sub_surf_ingrid(i,j) >= 0.0 ) THEN

      ! multiply by grid box area (in m) and divide by 1000 to convert
      ! depth in mm to depth in m (sjd 02/05/12)

      surf_roff(ip) = runoff_factor * rfm_surf_ingrid(i,j) *                  &
                            dt * rivers_boxareas_rp(ip) / mm_2_m
      sub_surf_roff(ip) = runoff_factor * rfm_sub_surf_ingrid(i,j) *          &
                            dt * rivers_boxareas_rp(ip) / mm_2_m

    ELSE

      ! ignore no data (-1.00e+20) and other negative values
      surf_roff(ip) = 0.0
      sub_surf_roff(ip) = 0.0

    END IF

  END DO

!-------------------------------------------------------------------------------
! Initialise accumulated inflows and stores for the next timestep
!-------------------------------------------------------------------------------

  flowin_n(:) = 0.0
  bflowin_n(:) = 0.0
  surfstore_n(:) = 0.0
  substore_n(:) = 0.0

!-------------------------------------------------------------------------------
! Rivers runoff using simple kinematic wave model
!-------------------------------------------------------------------------------

  DO ip = 1,np_rivers

    rn = rivers_next_rp(ip)
    rp = rivers_index_rp(ip)

    landtype = rfm_land_rp(rp)

    IF (landtype == 2) THEN  !land

      ! land surface
      surfstore_n(ip) = (1.0 - landtheta) * rfm_surfstore_rp(ip) +            &
                               rfm_flowin_rp(ip) + surf_roff(ip)

      ! land subsurface
      substore_n(ip) = (1.0 - sublandtheta) * rfm_substore_rp(ip) +           &
                               rfm_bflowin_rp(ip) + sub_surf_roff(ip)

      ! return flow
      IF (retl > 0) THEN ! changed to allow for -ve retr/l (sjd: 4/6/09)
        returnflow = AMAX1( ABS( substore_n(ip) * retl ), 0.0 )
      ELSE
        returnflow = -1.0 * AMAX1( ABS( surfstore_n(ip) * retl ), 0.0 )
      END IF

      substore_n(ip) = substore_n(ip) - returnflow
      surfstore_n(ip) = surfstore_n(ip) + returnflow

      flowin_n(rn) = flowin_n(rn) + landtheta * rfm_surfstore_rp(ip)
      bflowin_n(rn) = bflowin_n(rn) + sublandtheta * rfm_substore_rp(ip)

    ELSE IF (landtype == 1) THEN  !river

      ! river subsurface
      substore_n(ip) = (1.0 - subrivertheta) * rfm_substore_rp(ip) +          &
                            rfm_bflowin_rp(ip) + sub_surf_roff(ip)

      ! river surface
      surfstore_n(ip) = (1.0 - rivertheta) * rfm_surfstore_rp(ip) +           &
                              rfm_flowin_rp(ip) + surf_roff(ip)

      ! return flow
      IF (retr > 0) THEN ! changed to allow for -ve retr/l (sjd: 4/6/09)
        returnflow = AMAX1( ABS( substore_n(ip) * retr ), 0.0 )
      ELSE
        returnflow = -1.0 * AMAX1( ABS( surfstore_n(ip) * retr ), 0.0 )
      END IF
      substore_n(ip) = substore_n(ip) - returnflow
      surfstore_n(ip) = surfstore_n(ip) + returnflow

      flowin_n(rn) = flowin_n(rn) + rivertheta * rfm_surfstore_rp(ip)
      bflowin_n(rn) = bflowin_n(rn) + subrivertheta * rfm_substore_rp(ip)

    END IF ! land or river

  END DO !end of rivers loop, ip

!-------------------------------------------------------------------------------
! Housekeeping for next timestep
!-------------------------------------------------------------------------------
  DO ip = 1, np_rivers

    rp = rivers_index_rp(ip)
    rn = rivers_next_rp(ip)

    ! keep inflows for next timestep
    rfm_flowin_rp(ip) = flowin_n(rp)
    rfm_bflowin_rp(ip) = bflowin_n(rp)

    ! calculate flow in rivers (m3/s)
    IF (rfm_land_rp(rp) == 1) THEN

      rfm_rivflow_rp(rp) = rivertheta / dt * rfm_surfstore_rp(ip)
      rfm_baseflow_rp(rp) = subrivertheta / dt * rfm_substore_rp(ip)

      ! accumulate flow into the sea if the next point is sea
      !(assume it's equal to river flow in the adjacent land pt)
      IF (rfm_land_rp(rn) == 0) THEN
        rfm_rivflow_rp(rn) = rfm_rivflow_rp(rn) + rfm_rivflow_rp(rp)
      END IF

    END IF

    !calculate flow for land points (m3/s)
    IF (rfm_land_rp(rp) == 2) THEN

      rfm_rivflow_rp(rp) = landtheta / dt * rfm_surfstore_rp(ip)
      rfm_baseflow_rp(rp) = sublandtheta / dt * rfm_substore_rp(ip)

      ! accumulate flow into the sea if the next point is sea
      !(assume it's equal to river flow in the adjacent land pt)
      IF (rfm_land_rp(rn) == 0) THEN
        rfm_rivflow_rp(rn) = rfm_rivflow_rp(rn) + rfm_rivflow_rp(rp)
      END IF

    END IF

    ! keep rivers stores for next timestep  (m3)
    rfm_surfstore_rp(ip) = surfstore_n(ip)
    rfm_substore_rp(ip) = substore_n(ip)

  END DO !end rivers loop, ip

!-------------------------------------------------------------------------------
! Return flows in flux density units kg/m2/s
!-------------------------------------------------------------------------------
!! HL: N.B. Add option to select what output units required for river flow
  DO ip=1,np_rivers
    CALL rivers_get_xy_pos(rivers_index_rp(ip),nx_rivers,ny_rivers,i,j)
    rp = rivers_index_rp(ip)
    outflow(i,j) = rfm_rivflow_rp(rp) * mm_2_m / rivers_boxareas_rp(ip)
    baseflow(i,j) = rfm_baseflow_rp(rp) * mm_2_m / rivers_boxareas_rp(ip)
  END DO

END SUBROUTINE rivers_rfm_vector

!###############################################################################
!###############################################################################



END MODULE jules_rivers_rfm_mod
