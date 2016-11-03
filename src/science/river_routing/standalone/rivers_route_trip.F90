! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    MODULE jules_rivers_trip_mod-----------------------------------------------

! Description:
!     Driver and science routines for calculating river flow routing
!     using the TRIP (linear) model
!     see Oki et al 1999 J.Met.Soc.Japan, 77, 235-255.
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

MODULE jules_rivers_trip_mod

CONTAINS

!###############################################################################
! subroutine rivers_drive_trip
! Driver routine for runoff routing by the TRIP (linear) model.

SUBROUTINE rivers_drive_trip( sub_runoffin, surf_runoffin, runoff_out, rivflow )

!-------------------------------------------------------------------------------
!
! Description:
!   Perform the routing of total surface runoff using TRIP
!
!   This routine regrids the total surface runoff to the river grid and passes
!   it to the TRIP routines to be routed.
!
!-------------------------------------------------------------------------------
! Modules used:

  USE jules_rivers_mod, ONLY :                                                &
!  imported scalars with intent(in)
     np_rivers,nx_rivers,ny_rivers,rivers_regrid                              &
     ,rivers_dlat, rivers_dlon                                                &
    ,nx_grid, ny_grid                                                         &
!  imported arrays with intent(in)
    ,rivers_index_rp, rivers_lat2d, rivers_lon2d                              &
    ,rivers_boxareas                                                          &
!  imported arrays with intent(inout)
    ,rivers_next_rp

  USE rivers_utils, ONLY :                                                    &
!  imported procedures
     rivers_get_xy_pos, get_rivers_len,                                       &
     rivers_route_regrid, rivers_route_regrid_invert

!-------------------------------------------------------------------------------
  USE jules_print_mgr, ONLY :                                                 &
    jules_message,                                                            &
    jules_print

  IMPLICIT NONE

! Array arguments with intent(in)

  REAL, INTENT(IN) :: sub_runoffin(nx_grid, ny_grid)
                               ! Sub surface runoff on model grid in kg m-2 s-1
  REAL, INTENT(IN) :: surf_runoffin(nx_grid, ny_grid)
                               ! Surface runoff on model grid in kg m-2 s-1

! Array arguments with intent(out)

  REAL, INTENT(OUT) :: runoff_out(nx_grid, ny_grid)
                               ! Runoff diagnostic on model grid
  REAL, INTENT(OUT) :: rivflow(nx_grid, ny_grid)
                               ! River flow diagnostic on model grid

! Local scalar variables.

  INTEGER :: ip                ! loop counters
  INTEGER :: maxLin            ! an array size

! Local array variables.

  INTEGER :: x(np_rivers)      !  x index on active rivers points
  INTEGER :: xnext(np_rivers)  !  x location of the next downstream point
!                                 <0 indicates next point is off grid
  INTEGER :: y(np_rivers)      !  y index on active rivers points
  INTEGER :: ynext(np_rivers)  !  y location of the next downstream point
                               !  If xnext(i)<0, yNext(i) gives the direction
                               !  of the next point, in terms of index in
                               !  flow_dir_set.

  REAL :: runoff_grid(nx_grid, ny_grid)
                               !  average rate of runoff since last
                               !  rivers call (kg m-2 s-1), on a grid
                               !  that is a superset of (can be equal
                               !  to) the model grid
  REAL :: outflow(nx_rivers,ny_rivers)
                               !  rate of channel flow leaving gridbox
                               !  (kg m-2 s-1)
  REAL :: outflow_flux(nx_rivers,ny_rivers)
                               !  rate of channel flow leaving gridbox (kg s-1)
  REAL :: runoff_rivers(nx_rivers,ny_rivers)
                               !  average rate of runoff, on river grid
                               !  (kg m-2 s-1)
  REAL :: rflux_rivers(nx_rivers,ny_rivers)
                               ! average rate of runoff flux (kg s-1)
  REAL :: riverslength(nx_rivers,ny_rivers)
                               !  distance to next downstream gridpoint (m)
                               !  This is only calculated at "active" river
                               !  points.
  LOGICAL :: rivers_mask(nx_rivers,ny_rivers)
                               !  TRUE at land points on river grid, else FALSE

!-------------------------------------------------------------------------------
! Sum the surface and subsurface runoffs
!-------------------------------------------------------------------------------

  runoff_grid(:,:) = (surf_runoffin(:,:) + sub_runoffin(:,:)) ! * frac(:,:)

!-------------------------------------------------------------------------------
! Regrid runoff from land grid to river routing grid, and convert from
! flux density (kg m-2 s-1) to flux (kg s-1).
!-------------------------------------------------------------------------------

  runoff_rivers(:,:) = 0.0
  rflux_rivers(:,:) = 0.0

  IF ( rivers_regrid ) THEN

!   Calculate an array size used in regridding
    maxlin = ( nx_grid + nx_rivers )*( ny_grid + ny_rivers )
    CALL rivers_route_regrid( maxlin, runoff_grid, runoff_rivers )

!   Convert from flux density (kg m-2 s-1) to flux (kg s-1).
    rflux_rivers(:,:) = runoff_rivers(:,:) * rivers_boxareas(:,:)

  ELSE

!   "Main" and routing grids are identical (including points in same order).
!   No need to regrid, simply convert units from flux density (kg m-2 s-1) to
!   flux (kg s-1).
    rflux_rivers(:,:) = runoff_grid(:,:) * rivers_boxareas(:,:)

  END IF

!-------------------------------------------------------------------------------
! Calculate location on rivers grid of each land point, and find the point
! immediately downstream. Also set a mask to show land points on the rivers
! grid.
!-------------------------------------------------------------------------------

! Initialise.
  xnext(:) = 0
  ynext(:) = 0
  rivers_mask(:,:) = .FALSE.

  DO ip=1,np_rivers

!   Get location of this point.
    CALL rivers_get_xy_pos(rivers_index_rp(ip),nx_rivers,ny_rivers,x(ip),y(ip))

!   Set mask to TRUE at this land point.
    rivers_mask(x(ip),y(ip)) = .TRUE.

!   Get location of next downstream point.
    IF ( rivers_next_rp(ip) > 0 ) THEN
!     Downstream point is on the grid.
      CALL rivers_get_xy_pos(rivers_next_rp(ip),nx_rivers,ny_rivers,          &
                              xnext(ip),ynext(ip) )
    ELSE
!!    ELSEIF ( riversNext(ip) == -9 ) THEN
!     There is no defined flow direction, and so no downstream point.
      xnext(ip) = -9
      ynext(ip) = -9
!!    ELSE
!     Downstream point is off edge of grid.
!     rivers_next_rp has been set to -1*index of direction in flow_dir_set.
!!      xnext(ip) = -1 * rivers_next_rp(ip)
!!      ynext(ip) = xnext(ip)
    END IF

  END DO

!-------------------------------------------------------------------------------
! Calculate distance between grid points.
!-------------------------------------------------------------------------------
  CALL get_rivers_len( x,xnext,y,ynext,rivers_lat2d,rivers_lon2d,riverslength )

!-------------------------------------------------------------------------------
! Call the river routing routine.
!-------------------------------------------------------------------------------
  CALL rivers_trip( riverslength, rflux_rivers, outflow_flux )

!-------------------------------------------------------------------------------
! Add runoff on TRIP sea points to the river flow variable.
! Regridding may have resulted in some land runoff appearing in TRIP sea
! gridboxes, so we need to account for this water. Even without regridding,
! a land point that is sea in TRIP will not have been added to flow.
!-------------------------------------------------------------------------------
!!  WHERE( .NOT. rivers_mask(:,:) ) outflow(:,:) = rflux_rivers(:,:)

!-------------------------------------------------------------------------------
!   Regrid from rivers to land grid
!-------------------------------------------------------------------------------

  IF ( rivers_regrid ) THEN

    outflow(:,:) = outflow_flux(:,:)/rivers_boxareas(:,:)
    CALL rivers_route_regrid_invert( maxlin,outflow,rivflow )

    ! Sanity check regridding routines - expect input = output
    runoff_rivers(:,:) = rflux_rivers(:,:)/rivers_boxareas(:,:)
    CALL rivers_route_regrid_invert( maxlin, runoff_rivers, runoff_out )

  ELSE

    rivflow(:,:) = outflow_flux(:,:)/rivers_boxareas(:,:)
    runoff_out(:,:) = rflux_rivers(:,:)/rivers_boxareas(:,:)

  END IF

END SUBROUTINE rivers_drive_trip
!###############################################################################

!###############################################################################
! subroutine rivers_trip
!
!-----------------------------------------------------------------------------
! Description:
!   Calculates river outflow (kg s-1) and updates channel storage for the
!   TRIP river routing model.
!
! Method:
!   See Oki et al. 1999, J.Met.Soc.Japan, 77, 235-255.
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

SUBROUTINE rivers_trip( riverslength,runoff,outflow )

  USE jules_rivers_mod, ONLY :                                                &
!  imported scalars with intent(in)
     np_rivers,nx_rivers,ny_rivers,                                           &
     rivers_timestep,nseqmax                                                  &
!  imported arrays with intent(in)
    ,rivers_index_rp,rivers_next_rp,rivers_seq_rp,rivers_sto_rp

  USE jules_riversparm, ONLY :                                                &
!  imported parameters
     rivers_meander, rivers_speed

  USE rivers_utils, ONLY :                                                    &
!  imported procedures
     rivers_get_xy_pos

  USE timestep_mod, ONLY:                                                     &
     timestep

  USE jules_print_mgr, ONLY :                                                 &
    jules_message,                                                            &
    jules_print

  IMPLICIT NONE

!-------------------------------------------------------------------------------

  REAL, INTENT(IN) ::                                                         &
!  arrays with intent(in)
     riverslength(nx_rivers,ny_rivers)                                        &
     !  distance between gridpoints (m)
     ,runoff(nx_rivers,ny_rivers)
     !  rate of runoff generation in each gridbox (kg s-1)

  REAL, INTENT(OUT) ::                                                        &
!  arrays with intent(out)
     outflow(nx_rivers,ny_rivers)
     !  rate of channel flow leaving gridbox (kg s-1)

  INTEGER ::                                                                  &
!  local scalars
     i,ip,ix,iy,jx,jy,iseq       !  work/loop counters

  REAL ::                                                                     &
!  local scalars
     coeff                                                                    &
     !  coefficient in the routing model (s-1)
     ,exp_coeffdt                                                             &
     !  working variable exp[c*dt]
     ,dt                                                                      &
     !  timestep of routing model (s)
     ,store_old
     !  channel storage (kg)

  REAL ::                                                                     &
!  local arrays
    inflow(nx_rivers,ny_rivers)
    !  rate of channel flow entering gridbox (kg s-1)

!-------------------------------------------------------------------------------

  dt = REAL(rivers_timestep) * timestep

! Initialise inflow with runoff generated over each gridbox.
  inflow(:,:) = runoff(:,:)

! Initialise outlow (for clarity at non-rivers points).
  outflow(:,:) = 0.0

!-------------------------------------------------------------------------------
! Loop over active rivers points.
!-------------------------------------------------------------------------------

  DO iseq = 1, nseqmax

    DO ip=1,np_rivers

!   Get index (location in rivers vector) of the point to consider.
      IF (NINT(rivers_seq_rp(ip)) == iseq) THEN

!   Get location of this point in rivers grid.
        CALL rivers_get_xy_pos( rivers_index_rp(ip),nx_rivers,ny_rivers,ix,iy )

!-------------------------------------------------------------------------------
!   Calculate the coefficient "c" of the model.
!   c=u/(d*r), where u is effective flow speed,
!   d is distance between gridpoints, and r is meander ratio.
!-------------------------------------------------------------------------------
        coeff = rivers_speed / ( riverslength(ix,iy) * rivers_meander )
        exp_coeffdt = EXP(-(coeff*dt))

!-------------------------------------------------------------------------------
!   Save value of channel storage at start of timestep.
!-------------------------------------------------------------------------------
        store_old = rivers_sto_rp(ip)

!-------------------------------------------------------------------------------
!   Calculate channel storage at end of timestep.
!   Eqn.4 of Oki et al, 1999, J.Met.Soc.Japan, 77, 235-255.
!-------------------------------------------------------------------------------
        rivers_sto_rp(ip) = store_old * exp_coeffdt                           &
                         + (1.0 - exp_coeffdt) * inflow(ix,iy) / coeff

!-------------------------------------------------------------------------------
!   Calculate outflow as inflow minus change in storage.
!-------------------------------------------------------------------------------
        outflow(ix,iy) = inflow(ix,iy) + (store_old - rivers_sto_rp(ip)) / dt

!-------------------------------------------------------------------------------
!   Add outflow to inflow of next downstream point.
!-------------------------------------------------------------------------------
        IF ( rivers_next_rp(ip) > 0 ) THEN
!     Get location in grid of next downstream point.
          CALL rivers_get_xy_pos(rivers_next_rp(ip),nx_rivers,ny_rivers,jx,jy)
          inflow(jx,jy) = inflow(jx,jy) + outflow(ix,iy)
        END IF

      END IF

    END DO !points
  END DO !river sequences

END SUBROUTINE rivers_trip

!###############################################################################

SUBROUTINE regrid_routestore

  USE jules_rivers_mod, ONLY :                                                &
!  imported scalars with intent(in)
    np_rivers, nx_rivers, ny_rivers, nx_grid, ny_grid, rivers_regrid          &
!  imported array with intent(in)
    ,rivers_index_rp, rivers_sto_rp                                           &
    ,rivers_boxareas                                                          &
!  imported array with intent(out)
    ,rivers_sto_per_m2_on_landpts

  USE parallel_mod, ONLY : is_master_task,                                    &
    gather_land_field, scatter_land_field

  USE rivers_utils, ONLY :                                                    &
!  imported procedures
     rivers_get_xy_pos, rivers_map_from_model_grid                            &
    , rivers_route_regrid_invert

  USE model_grid_mod, ONLY : latitude_of_land_pts, longitude_of_land_pts,     &
    global_land_pts

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Uses to the river storage on the river points array (a prognostic) to
!   fill the river storage on land points array, which is used to calculate
!   the water available for irrigation.
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

! Local arrays
  REAL, ALLOCATABLE :: global_lat_of_land_pts(:)
            ! latitude of global model land points
  REAL, ALLOCATABLE :: global_lon_of_land_pts(:)
            ! longitude of global model land points
  REAL, ALLOCATABLE :: global_rivers_sto_per_m2_on_landpts(:)
            ! Water storage (kg m-2) on global model land points
  REAL :: rivers_sto_per_m2_rgrid(nx_rivers, ny_rivers)
            ! water storage on routing grid in kg m-2
  REAL :: rivers_sto_rgrid(nx_rivers, ny_rivers)
            ! water storage on routing grid in kg
  REAL :: rivers_sto_per_m2_model_grid(nx_grid, ny_grid)
            ! water storage on model grid
  REAL :: dummy_rivflow(nx_grid, ny_grid)
            ! dummy variable used as argument in rivers_map_from_model_grid
  REAL :: dummy_rflow(global_land_pts)
            ! dummy variable used as argument in rivers_map_from_model_grid
  INTEGER :: maxlin ! an array size
  INTEGER :: ip, ix, iy ! array indices

  rivers_sto_per_m2_on_landpts(:) = 0.0

  IF ( is_master_task() ) THEN
    ALLOCATE(global_lat_of_land_pts(global_land_pts))
    ALLOCATE(global_lon_of_land_pts(global_land_pts))
    ALLOCATE(global_rivers_sto_per_m2_on_landpts(global_land_pts))
  ELSE
    ALLOCATE(global_lat_of_land_pts(1))
    ALLOCATE(global_lon_of_land_pts(1))
    ALLOCATE(global_rivers_sto_per_m2_on_landpts(1))
  END IF

  CALL gather_land_field(latitude_of_land_pts, global_lat_of_land_pts)
  CALL gather_land_field(longitude_of_land_pts, global_lon_of_land_pts)

  global_rivers_sto_per_m2_on_landpts(:) = 0.0

  IF ( is_master_task() ) THEN
    rivers_sto_rgrid(:,:) = -1.0
    rivers_sto_per_m2_model_grid(:,:) = 0.0

    DO ip=1,np_rivers
      CALL rivers_get_xy_pos(rivers_index_rp(ip), nx_rivers, ny_rivers,       &
                            ix, iy)
      rivers_sto_rgrid(ix,iy) = rivers_sto_rp(ip)
    END DO

    rivers_sto_per_m2_rgrid(:,:) = rivers_sto_rgrid(:,:) / rivers_boxareas(:,:)

    IF ( rivers_regrid ) THEN

  !   Calculate an array size used in regridding
      maxlin = ( nx_grid + nx_rivers )*( ny_grid + ny_rivers )

      CALL rivers_route_regrid_invert( maxlin, rivers_sto_per_m2_rgrid,       &
                                       rivers_sto_per_m2_model_grid )

    ELSE

!   "Main" and routing grids are identical (including points in same order).
      rivers_sto_per_m2_model_grid(:,:) = rivers_sto_per_m2_rgrid(:,:)

    END IF

    CALL rivers_map_from_model_grid( global_lat_of_land_pts,                  &
                    global_lon_of_land_pts,rivers_sto_per_m2_model_grid,      &
                    dummy_rivflow, global_rivers_sto_per_m2_on_landpts,       &
                    dummy_rflow )
  END IF

  CALL scatter_land_field(global_rivers_sto_per_m2_on_landpts,                &
                          rivers_sto_per_m2_on_landpts)

  DEALLOCATE(global_rivers_sto_per_m2_on_landpts)
  DEALLOCATE(global_lon_of_land_pts)
  DEALLOCATE(global_lat_of_land_pts)

END SUBROUTINE regrid_routestore

!###############################################################################

SUBROUTINE adjust_routestore

!-----------------------------------------------------------------------------
! Description:
!   Uses to the river adjustment factor on land points (which accounts for
!   water extracted from the river storage for irrigation)
!   to adjust the river storage on river points array (a prognostic).
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

  USE jules_rivers_mod, ONLY :                                                &
!  imported scalars with intent(in)
    np_rivers, nx_rivers, ny_rivers, nx_grid, ny_grid, rivers_regrid          &
!  imported array with intent(in)
    ,rivers_index_rp, rivers_sto_rp, rivers_adj_on_landpts

  USE rivers_utils, ONLY :                                                    &
!  imported procedures
     rivers_get_xy_pos, rivers_map_to_model_grid                              &
    , rivers_route_regrid

  USE parallel_mod, ONLY : is_master_task, gather_land_field

  USE model_grid_mod, ONLY : latitude_of_land_pts, longitude_of_land_pts,     &
    global_land_pts

  IMPLICIT NONE

  REAL, ALLOCATABLE :: global_lat_of_land_pts(:)
            ! latitude of global model land points
  REAL, ALLOCATABLE :: global_lon_of_land_pts(:)
            ! longitude of global model land points
  REAL, ALLOCATABLE :: global_rivers_adj_on_landpts(:)
            ! rivers adjustment factor on global model land points

  REAL :: rivers_adj_model_grid(nx_grid, ny_grid)
            ! rivers adjustment factor on model grid
  REAL :: rivers_adj_rgrid(nx_rivers, ny_rivers)
            ! rivers adjustment factor on routing grid
  REAL :: dummy_surf_runoffin(nx_grid, ny_grid)
            ! dummy variable used as argument in rivers_map_to_model_grid
  REAL :: dummy_tot_surf_runoff(global_land_pts)
            ! dummy variable used as argument in rivers_map_to_model_grid

  INTEGER :: maxlin ! an array size
  INTEGER :: ip, ix, iy ! array indices

  IF ( is_master_task() ) THEN
    ALLOCATE(global_lat_of_land_pts(global_land_pts))
    ALLOCATE(global_lon_of_land_pts(global_land_pts))
    ALLOCATE(global_rivers_adj_on_landpts(global_land_pts))
  ELSE
    ALLOCATE(global_lat_of_land_pts(1))
    ALLOCATE(global_lon_of_land_pts(1))
    ALLOCATE(global_rivers_adj_on_landpts(1))
  END IF

  CALL gather_land_field(latitude_of_land_pts, global_lat_of_land_pts)
  CALL gather_land_field(longitude_of_land_pts, global_lon_of_land_pts)
  CALL gather_land_field(rivers_adj_on_landpts, global_rivers_adj_on_landpts)

  IF ( is_master_task() ) THEN
    rivers_adj_model_grid(:,:) = 1.0
    rivers_adj_rgrid(:,:) = 1.0

    CALL rivers_map_to_model_grid( global_lat_of_land_pts,                    &
                                   global_lon_of_land_pts,                    &
                                   global_rivers_adj_on_landpts,              &
                                   dummy_tot_surf_runoff,                     &
                                   rivers_adj_model_grid, dummy_surf_runoffin )

    IF ( rivers_regrid ) THEN

!   Calculate an array size used in regridding
      maxlin = ( nx_grid + nx_rivers )*( ny_grid + ny_rivers )

      CALL rivers_route_regrid( maxlin, rivers_adj_model_grid, rivers_adj_rgrid)

    ELSE

!   "Main" and routing grids are identical (including points in same order).
      rivers_adj_rgrid(:,:) = rivers_adj_model_grid(:,:)

    END IF

!-------------------------------------------------------------------------------
! Apply the correction factor to the route storage, which is a prognostic
!-------------------------------------------------------------------------------

    DO ip=1,np_rivers
      CALL rivers_get_xy_pos(rivers_index_rp(ip), nx_rivers, ny_rivers,       &
                           ix, iy)
      rivers_sto_rp(ip) = rivers_adj_rgrid(ix,iy) * rivers_sto_rp(ip)
    END DO

  END IF

  DEALLOCATE(global_rivers_adj_on_landpts)
  DEALLOCATE(global_lon_of_land_pts)
  DEALLOCATE(global_lat_of_land_pts)

END SUBROUTINE adjust_routestore

!###############################################################################

END MODULE jules_rivers_trip_mod

