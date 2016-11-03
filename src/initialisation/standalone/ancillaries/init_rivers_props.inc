#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

SUBROUTINE init_rivers_props()

  USE mpi, ONLY : MPI_COMM_WORLD

  USE io_constants, ONLY : MAX_SDF_NAME_LEN, MAX_FILE_NAME_LEN, NAMELIST_UNIT
 
  USE grid_utils_mod, ONLY : grid_create, grid_info

  USE input_mod, ONLY :                                                       &
      od_grid => grid,                                                        &
      dummy_grid => grid,                                                     &
      fill_variables_from_file, use_subgrid

  USE missing_data_mod, ONLY: imdi, rmdi

  USE ancil_info, ONLY : land_pts

  USE string_utils_mod, ONLY : to_string

  USE templating_mod, ONLY : tpl_has_var_name, tpl_substitute_var
 
  USE model_interface_mod, ONLY : IDENTIFIER_LEN, populate_var, get_var_id

  USE model_grid_mod, ONLY : latitude, longitude, global_land_pts, model_grid

  USE c_pi, ONLY : pi_over_180
  
  USE planet_constants_mod, ONLY : planet_radius  ! the Earth's radius (m)
      
  USE jules_rivers_mod, ONLY :                                                &
       l_rivers,                                                              &
       rivers_rfm, rivers_trip, rivers_umrfm, rivers_umtrip,                  &
       flow_dir_delta, rivers_type, np_rivers, nx_rivers, ny_rivers,          &
       rivers_dlat, rivers_dlon, rivers_lat1, rivers_lon1, nx_grid, ny_grid,  &
       nseqmax, reg_lon1, reg_lat1, reg_dlon, reg_dlat, rivers_dx,            &
       rivers_reglatlon, rivers_regrid, l_rivers,                             &
       rivers_index_rp, rivers_next_rp, rivers_seq_rp, rivers_dir_rp,         &
       rivers_sto_rp, rivers_dra_rp, rivers_lat_rp, rivers_lon_rp,            &
       rivers_seq, rivers_dir, rivers_sto, rivers_dra, rivers_boxareas,       &
       rivers_lat2d, rivers_lon2d, rivers_xgrid, rivers_ygrid,                &
       rfm_flowobs1, rfm_surfstore, rfm_substore, rfm_flowin, rfm_bflowin,    &
       rfm_rivflow, rfm_baseflow, rfm_inext, rfm_jnext, rfm_iarea, rfm_land,  &
       rfm_flowobs1_rp, rfm_surfstore_rp, rfm_substore_rp, rfm_flowin_rp,     &
       rfm_bflowin_rp, rfm_rivflow_rp, rfm_baseflow_rp,                       &
       rfm_iarea_rp, rfm_land_rp, rivers_boxareas_rp,                         &
       rivers_sto_per_m2_on_landpts, rivers_adj_on_landpts
 
  USE rivers_utils, ONLY: rivers_get_xy_pos, rivers_earth_area
  
  USE jules_vegetation_mod, ONLY: l_irrig_limit

  USE jules_print_mgr, ONLY :                                                 &
    jules_message,                                                            &
    jules_print

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!    Reads in river routing parameters
!
! Current Code Owner:  Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

! Work variables

  INTEGER :: error  ! Error indicator
  
  INTEGER :: ntasks ! Parallel mode indicator

  INTEGER, PARAMETER :: MAX_RIVERS_VARS = 6 ! The maximum possible number of
                                            ! routing variables 

  INTEGER :: nvars_required      ! The number of routing variables that are
                                 ! required in this configuration
  INTEGER :: nvars_optional      ! The number of optional routing variables 
                                 ! in this configuration
  CHARACTER(len=IDENTIFIER_LEN) :: required_vars(MAX_RIVERS_VARS)
                                 ! The variable identifiers of the required
                                 ! variables
  CHARACTER(len=IDENTIFIER_LEN) :: optional_vars(MAX_RIVERS_VARS)
                                 ! The variable identifiers of any optional
                                 ! variables

  INTEGER :: nvars_file       ! The number of variables that will be set
                              ! from the given file (template?)

! Variables passed to fill_variables_from_file

  CHARACTER(len=IDENTIFIER_LEN) :: file_var(MAX_RIVERS_VARS)
                        ! The variable identifiers of the variables to set
                        ! from file
  CHARACTER(len=MAX_SDF_NAME_LEN) :: file_var_name(MAX_RIVERS_VARS)
                        ! The name of each variable in the file
  CHARACTER(len=MAX_SDF_NAME_LEN) :: file_tpl_name(MAX_RIVERS_VARS)
                        ! The name to substitute in a template for each
                        ! variable

  CHARACTER(len=MAX_FILE_NAME_LEN) :: file
                        ! The name of the file (or variable name template) to
                        ! use for variables that need to be filled from file

  INTEGER :: nvars      ! The number of variables in this section
  CHARACTER(len=IDENTIFIER_LEN) :: var(MAX_RIVERS_VARS)
                        ! The variable identifiers of the variables
  LOGICAL :: use_file(MAX_RIVERS_VARS)
                        !   T - the variable uses the file
                        !   F - the variable is set using a constant value
  CHARACTER(len=MAX_SDF_NAME_LEN) :: var_name(MAX_RIVERS_VARS)
                        ! The name of each variable in the file
  CHARACTER(len=MAX_SDF_NAME_LEN) :: tpl_name(MAX_RIVERS_VARS)
                        ! The name to substitute in a template for each
                        ! variable
  REAL :: const_val(MAX_RIVERS_VARS)
                        ! The constant value to use for each variable if
                        ! use_file = F for that variable

  INTEGER :: i,ip,ix,iy,l  ! Index variables
  INTEGER :: inext, jnext
  INTEGER :: i1, i2, ilat, ilon, step
  REAL :: reg_lat2, reg_lon2, min_lon

  LOGICAL :: use_sub_local
  LOGICAL :: set_dra

! The names and sizes of the x and y dimensions

  CHARACTER(len=MAX_SDF_NAME_LEN) :: x_dim_name, y_dim_name
  INTEGER :: nx, ny
  
  TYPE(grid_info) :: local_grid, rivers_grid

  NAMELIST /jules_rivers_props/                                               &
     x_dim_name, y_dim_name, nx, ny, file,                                    &
     nvars, var, use_file, var_name, tpl_name, const_val,                     &
     rivers_regrid, rivers_reglatlon,                                         &
     nx_grid, ny_grid, reg_lat1, reg_lon1, reg_dlon, reg_dlat, rivers_dx

!-----------------------------------------------------------------------------
! Notes on river routing grid assumptions (see also jules_rivers_mod):
!
!        The JULES river routing implementation assumes that river routing 
!        variables are defined on a river grid with cell (1,1) in the lower
!        left corner (i.e. W-E, S-N). Ancillary information may be read in
!        on alternative grids, but this routine will translate these points
!        to the anticipated ordering.
!
!        The JULES river routing algorithms also require information on the
!        dimensions of the full model grid, which may not be typically defined 
!        for some standalone applications where only land_pts are specified.
!        Information on the full model grid is therefore also read in and 
!        checked here.
!
!-----------------------------------------------------------------------------

! If rivers are not on, there is nothing to do
  IF( .NOT. l_rivers ) RETURN

!-----------------------------------------------------------------------------
! Initialise
!-----------------------------------------------------------------------------
  
  nvars_required = 0
  nvars_optional = 0
  nvars_file     = 0
  nvars          = 0
  use_file(:)    = .TRUE. ! Default is for every variable to be read from file
  file           = ''     ! Empty file name.
  
  IF ( global_land_pts <= 1 ) THEN
     l_rivers = .FALSE.
     CALL log_warn("init_rivers_props",                                        &
                    "River routing not appropriate for single point runs. " // &
                    "River routing disabled.")

  ELSE

!-----------------------------------------------------------------------------
! Read river routing namelist
!-----------------------------------------------------------------------------
  CALL log_info("init_rivers_props", "Reading JULES_RIVERS_PROPS namelist...")

  READ(NAMELIST_UNIT, nml=jules_rivers_props, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_rivers_props",                                       &
                   "Error reading namelist JULES_RIVERS_PROPS " //            &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

!-----------------------------------------------------------------------------
! Define regular routing grid lat/lon
!-----------------------------------------------------------------------------

  use_sub_local = use_subgrid
  use_subgrid = .FALSE.

! Temporarily copy saved grid (set for full model grid) to a local variable 
! before overwriting to define the river routing grid
  
  local_grid = grid_create(dummy_grid%is_1d,dummy_grid%dim_name,              &
                    dummy_grid%nx,dummy_grid%x_name,dummy_grid%nx,            & 
                    dummy_grid%y_name,dummy_grid%ny)

! Read latitude dimension from file

  od_grid = grid_create(grid_is_1d=.TRUE.,                                    &
                        dim_name=TRIM(y_dim_name),npoints=ny)
  ALLOCATE(rivers_ygrid(ny), STAT=error)
  CALL fill_variables_from_file(file,                                         &
                                (/ 'rivers_ygrid' /),(/ TRIM(y_dim_name) /))

! Read longitude dimension from file

  od_grid = grid_create(grid_is_1d=.TRUE.,                                    &
                        dim_name=TRIM(x_dim_name),npoints=nx)
  ALLOCATE(rivers_xgrid(nx), STAT=error)
  CALL fill_variables_from_file(file,                                         &
                                (/ 'rivers_xgrid' /),(/ TRIM(x_dim_name) /))
 
!-----------------------------------------------------------------------------
! Build the routing grid object from the namelist values
!-----------------------------------------------------------------------------

  dummy_grid = grid_create(.FALSE., "", 0, x_dim_name, nx, y_dim_name, ny)
  rivers_grid = grid_create(.FALSE., "", 0, x_dim_name, nx, y_dim_name, ny)

  nx_rivers = rivers_grid%nx
  ny_rivers = rivers_grid%ny

  ! Set all arrays S to N by default

  IF ( rivers_ygrid(1) > MINVAL(rivers_ygrid) ) THEN
    i1 = ny_rivers
    i2 = 1
    step = -1
  ELSE
    i1 = 1
    i2 = ny_rivers
    step = 1
  ENDIF

!-----------------------------------------------------------------------------
! Set up routing properties using namelist values
!-----------------------------------------------------------------------------

! Set up the required and optional routing variables

  SELECT CASE ( rivers_type )
    
    CASE ( rivers_trip, rivers_umtrip )
    
      nvars_required = 2
      required_vars(1:nvars_required) = (/ 'direction', 'sequence ' /)
      nvars_optional = 0

    CASE ( rivers_rfm, rivers_umrfm )

      IF ( rivers_reglatlon ) THEN
         nvars_required = 1
         required_vars(1:nvars_required) = (/ 'direction' /)
         nvars_optional = 4
         optional_vars(1:nvars_optional) =                                    &
           (/ 'latitude_2d ', 'longitude_2d',                                 &
              'area        ', 'sequence    ' /)
      ELSE
         nvars_required = 3
         required_vars(1:nvars_required) = (/ 'direction   ',                 &
                                     'latitude_2d ', 'longitude_2d' /)        
         nvars_optional = 2
         optional_vars(1:nvars_optional) =                                    &
           (/ 'area        ', 'sequence    ' /)
      END IF

    CASE default
      CALL log_fatal("init_rivers_props",                                     &
                     "Do not recognise rivers_type= '" //                     &
                     TRIM(rivers_type) // "'")
    END SELECT
    
!-----------------------------------------------------------------------------
! Check that all the required variables are there
!-----------------------------------------------------------------------------
  DO i = 1,nvars_required
    IF ( .NOT. ANY(var(1:nvars) == required_vars(i)) )                        &
      CALL log_fatal("init_rivers_props",                                     &
                     "No value given for required variable '" //              &
                     TRIM(required_vars(i)) // "'")
  END DO

!-----------------------------------------------------------------------------
! Allocate routing-specific arrays
!-----------------------------------------------------------------------------
  
  IF ( nx_rivers > 0 .AND. ny_rivers > 0) THEN

    ALLOCATE( rivers_dir(nx_rivers, ny_rivers), STAT=error )
    ALLOCATE( rivers_seq(nx_rivers, ny_rivers), STAT=error )
    ALLOCATE( rivers_sto(nx_rivers, ny_rivers), STAT=error )
    ALLOCATE( rivers_dra(nx_rivers, ny_rivers), STAT=error )

    ALLOCATE( rivers_lat2d(nx_rivers, ny_rivers), STAT=error )
    ALLOCATE( rivers_lon2d(nx_rivers, ny_rivers), STAT=error )

    ALLOCATE( rivers_boxareas(nx_rivers, ny_rivers), STAT=error )

! Initialise to impossible values.

    rivers_dir(:,:)   = rmdi    
    rivers_seq(:,:)   = rmdi 
    rivers_sto(:,:)   = rmdi
    rivers_dra(:,:)   = rmdi
    rivers_lat2d(:,:) = rmdi
    rivers_lon2d(:,:) = rmdi
    rivers_boxareas(:,:) = rmdi
 
! Allocate arrays required for RFM only
 
    IF (rivers_type == rivers_rfm .OR. rivers_type == rivers_umrfm) THEN

        ALLOCATE( rfm_flowobs1(nx_rivers, ny_rivers),  STAT=error )	    
        ALLOCATE( rfm_surfstore(nx_rivers, ny_rivers), STAT=error )
        ALLOCATE( rfm_substore(nx_rivers, ny_rivers),  STAT=error )
        ALLOCATE( rfm_flowin(nx_rivers, ny_rivers),    STAT=error )
        ALLOCATE( rfm_bflowin(nx_rivers, ny_rivers),   STAT=error )
        ALLOCATE( rfm_rivflow(nx_rivers, ny_rivers),   STAT=error )
        ALLOCATE( rfm_baseflow(nx_rivers, ny_rivers),  STAT=error )
        ALLOCATE( rfm_inext(nx_rivers, ny_rivers),     STAT=error )
        ALLOCATE( rfm_jnext(nx_rivers, ny_rivers),     STAT=error )
        ALLOCATE( rfm_land(nx_rivers, ny_rivers),      STAT=error )
        ALLOCATE( rfm_iarea(nx_rivers, ny_rivers),     STAT=error )
    
        rfm_inext(:,:)   = imdi
        rfm_jnext(:,:)   = imdi

        rfm_land(:,:) = 0       ! Initialise all points to sea
        rfm_iarea(:,:) = imdi

!! TO FIX - read initial flow from file if available. 
!!          For now, initialise variable to zero
    ! initialise observed flows
        rfm_flowobs1(:,:) = 0.0
    
    ! initialise surface and sub-surface stores and inflows                  
        rfm_surfstore(:,:) = 0.0
        rfm_substore(:,:) = 0.0

        rfm_flowin(:,:) = 0.0                                             
        rfm_bflowin(:,:) = 0.0                                            
        rfm_rivflow(:,:) = 0.0                                            
    
    END IF

  ELSE 

    CALL log_fatal("init_rivers_props",                                        &
                    "Invalid routing grid dimensions provided (nx_rivers = "// &
                    TRIM(to_string(nx_rivers)) // ", ny_rivers = " //          &
                    TRIM(to_string(ny_rivers)) // "). Check inputs. ")

  END IF
  
!-----------------------------------------------------------------------------
! Check which variables we will be using and partition them into variables
! set to constant values and variables set from file
!-----------------------------------------------------------------------------
  DO i = 1,nvars
     
!-----------------------------------------------------------------------------
! If the variable is one of the required vars, then we will be using it
!-----------------------------------------------------------------------------

    IF ( ANY(required_vars(1:nvars_required) == TRIM(var(i)))                 &
       .OR. (nvars_optional > 0 .AND.                                         &
       ANY(optional_vars(1:nvars_optional) == TRIM(var(i))) )) THEN
      IF ( use_file(i) ) THEN
        CALL log_info("init_rivers_props",                                    &
                      "'" // TRIM(var(i)) // "' will be read from file")

! If the variable will be filled from file, register it here
        nvars_file = nvars_file + 1
        file_var(nvars_file) = TRIM(var(i))
        file_var_name(nvars_file) = var_name(i)
        file_tpl_name(nvars_file) = TRIM(tpl_name(i))
      
      ELSE

! If the variable is being set as a constant, just populate it here
        CALL log_info("init_rivers_props",                                    &
                      "'" // TRIM(var(i)) // "' will be set to a " //         &
                      "constant = " // to_string(const_val(i)))

        CALL populate_var(get_var_id(var(i)), CONST_VAL=const_val(i))
      END IF
    
    ELSE

! If the variable is not a required variable, warn about not using it
      CALL log_warn("init_rivers_props",                                      &
                    "Provided variable '" // TRIM(var(i)) //                  &
                    "' is not required, so will be ignored")
    END IF
   
  END DO ! nvars
  
!-----------------------------------------------------------------------------
! Set variables from file
!-----------------------------------------------------------------------------

  IF ( nvars_file > 0 ) THEN
    
! Check that a file name was provided
    IF ( LEN_TRIM(file) == 0 )                                                &
       CALL log_fatal("init_rivers_props", "No file name provided")
    
    CALL log_info("init_rivers_props",                                        &
                "Reading routing information from file " // TRIM(file))
    
    IF ( tpl_has_var_name(file) ) THEN
! We are using a file name template, so loop through the variables setting
! one from each file
      DO i = 1,nvars_file
        CALL fill_variables_from_file(                                        &
          tpl_substitute_var(file, file_tpl_name(i)),                         &
          (/ file_var(i) /), (/ file_var_name(i) /)                           &
        )
      END DO

    ELSE
! We are not using a file name template, so set all variables from the same
! file
      CALL fill_variables_from_file(                                          &
        file,file_var(1:nvars_file), file_var_name(1:nvars_file))          
      
    END IF
  END IF

!-----------------------------------------------------------------------------
! Setup regular river routing grid (assuming same as model grid if not set 
! from namelist or model grid information)
!-----------------------------------------------------------------------------
  
  ! Wrap longitude inputs -180 to +180
  DO ix=1,nx
    IF (rivers_xgrid(ix) >= 180.0) rivers_xgrid(ix) = rivers_xgrid(ix)-360.0
  END DO
    
  ! Set regular Lat/Lon grid variables (S to N)
  rivers_lat1 = MINVAL(rivers_ygrid)
  rivers_lon1 = MINVAL(rivers_xgrid)
  rivers_dlat = ABS(rivers_ygrid(2) - rivers_ygrid(1))
  rivers_dlon = ABS(rivers_xgrid(2) - rivers_xgrid(1))

  IF (reg_lon1 <= -900.0 .OR. reg_lat1 <= -900.0 ) THEN   
      
    ! Infer model grid settings for 1-dimensional grid input
    IF ( model_grid%nx == 1 .OR. .NOT. rivers_reglatlon ) THEN   

       CALL log_info("init_rivers_props",                                     &
               "No regular model grid set - assuming same as routing grid")  
       nx_grid = nx_rivers
       ny_grid = ny_rivers
       reg_lat1 = rivers_lat1
       reg_lon1 = rivers_lon1
       reg_dlon = rivers_dlon
       reg_dlat = rivers_dlat
  
    ELSE   
    
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ntasks,error)
      IF (ntasks > 1 .AND.                                                    &
           (rivers_type == rivers_rfm .OR. rivers_type == rivers_trip) ) THEN   
          CALL log_fatal("init_rivers_props",                                 &
                         "This combination of settings in parallel mode " //  &
                         "requires reg_lat1 and reg_lon1 to be set in the "// &
                         "init_rivers_props namelist.")
      END IF

      CALL log_info("init_rivers_props",                                      &
               "Using 2D input model grid to define grid parameters")  
      nx_grid = SIZE(latitude,1)
      ny_grid = SIZE(latitude,2)
      reg_lat1 = MINVAL(latitude)
      reg_lon1 = MINVAL(longitude)
      reg_dlon = ABS(longitude(2,2) - longitude(1,1))
      reg_dlat = ABS(latitude(2,2) - latitude(1,1))

    END IF
  ELSE
    IF ( reg_dlat <= -900.0 .OR. reg_dlon <= -900.0 .OR.                          &
         ny_grid <= 0 .OR. nx_grid <= 0 ) THEN
      CALL log_info("init_rivers_props",                                          &
                    "If reg_lon1 and  reg_lon1 are specified in namelist, " //    &
                    "specify reg_dlat, reg_dlon, ny_grid, nx_grid too." )
    END IF
  END IF
 
  CALL log_info("init_rivers_props","Setting regular routing grid " //        &
            " Minimum " // TRIM(x_dim_name) // ": " //                        &
              TRIM(to_string(MINVAL(rivers_xgrid)))  //                       &
            " Minimum " // TRIM(y_dim_name) // ": " //                        &
              TRIM(to_string(MINVAL(rivers_ygrid))) //                        &
            " Maximum " // TRIM(x_dim_name) // ": " //                        &
              TRIM(to_string(MAXVAL(rivers_xgrid))) //                        &
            " Maximum " // TRIM(y_dim_name) // ": " //                        &
              TRIM(to_string(MAXVAL(rivers_ygrid)))) 
      
!-----------------------------------------------------------------------------
! Process the namelist values and set derived variables
!-----------------------------------------------------------------------------
 
! Consider routing points on valid domain only - otherwise unset rivers_dir
!-----------------------------------------------------------------------------

  reg_lat2 = reg_lat1 + ny_grid*reg_dlat
  reg_lon2 = reg_lon1 + nx_grid*reg_dlon

  DO iy=1,ny_rivers
    DO ix=1,nx_rivers
     IF (.NOT. ( (rivers_xgrid(ix) > reg_lon1) .AND.                           &
                  (rivers_ygrid(iy) > reg_lat1) .AND.                          &
                  (rivers_xgrid(ix) < reg_lon2) .AND.                          &
                  (rivers_ygrid(iy) < reg_lat2) ))THEN
       rivers_dir(ix,iy) = -1
     END IF
    END DO
  END DO

! Detect number of river points in river grid
!----------------------------------------------------------------------------- 
  
  np_rivers = COUNT(rivers_dir >= 1 .AND. rivers_dir <=10)

  CALL log_info("init_rivers_props",                                          &
                "River routing points = " // TRIM(to_string(np_rivers)))
  
  IF (np_rivers <= 0) THEN
     CALL log_fatal("init_rivers_props",                                      &
                    "Invalid number of valid routing grid points " //         & 
                    "(np_rivers = " // TRIM(to_string(np_rivers)) //          &
                    "). Check inputs. ")
  END IF

! Allocate land point arrays and initialise
!-----------------------------------------------------------------------------

  IF ( l_irrig_limit ) THEN
    ALLOCATE( rivers_sto_per_m2_on_landpts(land_pts), STAT=error )
    ALLOCATE( rivers_adj_on_landpts(land_pts), STAT=error )
    
    rivers_sto_per_m2_on_landpts(:) = 0.0
    rivers_adj_on_landpts(:) = 1.0
  END IF
  
! Allocate routing point arrays
!-----------------------------------------------------------------------------
  
  ALLOCATE(rivers_index_rp(np_rivers), STAT=error)
  ALLOCATE(rivers_next_rp(np_rivers), STAT=error)
  ALLOCATE(rivers_seq_rp(np_rivers), STAT=error)
  ALLOCATE(rivers_sto_rp(np_rivers), STAT=error)
  ALLOCATE(rivers_dir_rp(np_rivers), STAT=error)
  ALLOCATE(rivers_dra_rp(np_rivers), STAT=error)
  ALLOCATE(rivers_lat_rp(np_rivers), STAT=error)
  ALLOCATE(rivers_lon_rp(np_rivers), STAT=error)
  ALLOCATE(rivers_boxareas_rp(np_rivers), STAT=error)
 
  ! Initialise array values
 
  rivers_index_rp(:) = imdi
  rivers_next_rp(:) = 0
  rivers_seq_rp(:) = rmdi
  rivers_dir_rp(:) = rmdi
  rivers_sto_rp(:) = 0.0

  IF (rivers_type == rivers_rfm .OR. rivers_type == rivers_umrfm) THEN
    ALLOCATE( rfm_flowobs1_rp(np_rivers),  STAT=error )    
    ALLOCATE( rfm_surfstore_rp(np_rivers), STAT=error )
    ALLOCATE( rfm_substore_rp(np_rivers),  STAT=error )
    ALLOCATE( rfm_flowin_rp(np_rivers),    STAT=error )
    ALLOCATE( rfm_bflowin_rp(np_rivers),   STAT=error )
    ALLOCATE( rfm_iarea_rp(np_rivers),     STAT=error )

    ALLOCATE( rfm_land_rp(nx_rivers*ny_rivers),      STAT=error )
    ALLOCATE( rfm_rivflow_rp(nx_rivers*ny_rivers),   STAT=error )
    ALLOCATE( rfm_baseflow_rp(nx_rivers*ny_rivers),  STAT=error )

  ! Initialise array values

    rfm_land_rp(:) = 0       ! Initialise all points to sea
    rfm_iarea_rp(:) = imdi

    rfm_flowobs1_rp(:) = 0.0
    rfm_surfstore_rp(:) = 0.0
    rfm_substore_rp(:) = 0.0
    rfm_flowin_rp(:) = 0.0                                             
    rfm_bflowin_rp(:) = 0.0                                            
    rfm_rivflow_rp(:) = 0.0 
    rfm_baseflow_rp(:) = 0.0

  END IF

  inext = imdi
  jnext = imdi
  nseqmax = imdi
  ip = 0

  set_dra = .False.
  IF ( MAXVAL(rivers_dra) > 0.) set_dra=.True. ! Initialise iarea from drainage 

  min_lon = MINVAL(rivers_xgrid)

  ! Set up routing point arrays, and correct 2d lat/lon grid (S to N)

  DO ix=1,nx_rivers

    ilat = 0
    ilon = NINT( (rivers_xgrid(ix) - min_lon) / rivers_dlon ) + 1
    IF (ilon > nx_rivers) ilon = ilon-nx_rivers

    DO iy=i1,i2,step      
      
      ilat = ilat + 1

      IF ( rivers_reglatlon ) THEN
        ! Assume 1d x and y grid variables define lat / lon
        rivers_lat2d(ilon,ilat) = rivers_ygrid(iy)
        rivers_lon2d(ilon,ilat) = rivers_xgrid(ix)
      ENDIF
      
      ! Set appropriate routing flags etc
      IF ( rivers_dir(ix,iy) >= 1 .AND. rivers_dir(ix,iy) <=10) THEN

        ip = ip + 1        
        rivers_index_rp(ip) = (ilat-1)*nx_rivers + ilon
        rivers_lon_rp(ip) = rivers_lon2d(ilon,ilat) 
        rivers_lat_rp(ip) = rivers_lat2d(ilon,ilat)
        rivers_dir_rp(ip) = rivers_dir(ix,iy)
        rivers_seq_rp(ip) = rivers_seq(ix,iy)
        IF ( NINT(rivers_seq_rp(ip)) > nseqmax) nseqmax = NINT(rivers_seq_rp(ip))

        rivers_dra_rp(ip) = rivers_dra(ix,iy)        

! Set next variables from flow directions

        inext = ilon + flow_dir_delta(NINT(rivers_dir_rp(ip)),1)
        jnext = ilat + flow_dir_delta(NINT(rivers_dir_rp(ip)),2)

        rivers_next_rp(ip) = (jnext - 1)*nx_rivers + inext

! If the downstream point is the same as the current point, set to zero.
        IF (rivers_next_rp(ip) == rivers_index_rp(ip)) THEN
            rivers_next_rp(ip) = 0
        ENDIF

        IF (rivers_type == rivers_rfm .OR. rivers_type == rivers_umrfm) THEN

            rfm_inext(ilon,ilat) = inext
            rfm_jnext(ilon,ilat) = jnext

            rfm_land(ilon,ilat) = 2   ! Set to all routing grid points 
                                      ! to land point
    
            IF (set_dra) THEN       ! Initialise from iarea ancil 
                rfm_iarea(ilon,ilat) = rivers_dra(ix,iy)
            ELSE                    ! If not present, set all to river points
                rfm_iarea(ilon,ilat) = 2  !! potential to use riversSeq_rp(ip)
            ENDIF

!! HL: Note initial model conditions should be read via initial_conditions 
!      namelist in any case
!! TO FIX - initialise from ancil (n.b. units)
            rfm_flowobs1(ilon,ilat) = 0.0  !! potential to use riversSto_rp(ip) 

! 1D river points variables

           rfm_iarea_rp(ip) = rfm_iarea(ilon,ilat)
           rfm_flowobs1_rp(ip) = rfm_flowobs1(ilon,ilat)
           rfm_surfstore_rp(ip) = rfm_surfstore(ilon,ilat)
           rfm_substore_rp(ip) =  rfm_substore(ilon,ilat)
           rfm_flowin_rp(ip) = rfm_flowin(ilon,ilat)
           rfm_bflowin_rp(ip) = rfm_bflowin(ilon,ilat)

           rfm_land_rp(rivers_index_rp(ip)) = rfm_land(ilon,ilat)
           rfm_rivflow_rp(rivers_index_rp(ip)) = rfm_rivflow(ilon,ilat)
           rfm_baseflow_rp(rivers_index_rp(ip)) = rfm_baseflow(ilon,ilat)

        END IF        

      END IF

    END DO ! iy river rows

  END DO ! ix river columns

!-----------------------------------------------------------------------------
! Check routing and model grid settings are appropriate
!-----------------------------------------------------------------------------
  
  ! Check that full routing latitude/longitude grid within range
  IF ( ANY(rivers_lat2d < -90.0) .OR. ANY(rivers_lat2d > 90.0) )             &
    CALL log_fatal("init_rivers_props",                                      &
                   "Latitude is out of range - allowed range is " //         &
                   "-90.0 to 90.0, given range is " //                       &
                   TRIM(to_string(MINVAL(rivers_lat2d))) // " to " //        &
                   TRIM(to_string(MAXVAL(rivers_lat2d))))

  IF ( ANY(rivers_lon2d < -180.0) .OR. ANY(rivers_lon2d > 360.0) )           &
     CALL log_fatal("init_rivers_props",                                     &
                   "Longitude is out of range - allowed range is " //        &
                   "-180.0 to 360.0, given range is " //                     &
                   TRIM(to_string(MINVAL(rivers_lon2d))) // " to " //        &
                   TRIM(to_string(MAXVAL(rivers_lon2d))))

  ! Checks for regular latitude/longitude grid
  IF ( rivers_reglatlon ) THEN

    ! Check cells are 'square' if grid defined as regular
    IF ( abs(rivers_dlat) /= abs(rivers_dlon) ) THEN 
      CALL log_fatal("init_rivers_props",                                    &
         "Error in routing - latitude and longitude grid steps " //          &
         " must be equal if set as regular grid. rivers_dlat =" //           &
         TRIM(to_string(rivers_dlat)) // " , rivers_dlon =" //               &
         TRIM(to_string(rivers_dlon)))
    ENDIF
    
  ELSE
    
    ! Must be regular for any regridding methods
    IF ( rivers_regrid ) THEN
      WRITE(jules_message,*) "ERROR: init_rivers: routing codes with " //     &
                             "regridding only currently exist for " //        &
                             "regular lat/lon grids."
      CALL jules_print('init_rivers_props',jules_message)
      CALL log_fatal("init_rivers_props",                                     &
                     "Error in routing - non-regular lat/lon river grid " //  &
                     "defined and regridding requested. This is not a " //    &
                     "valid option currently.")
   
    END IF

    ! Check if rivers_dx set
    IF ( rivers_dx <= 0 ) THEN
        WRITE(jules_message,*) "ERROR: init_rivers: routing codes with " //   &
                             "non-regular lat/lon grid require setting " //   &
                             "rivers_dx in the namelist." 
        CALL jules_print('init_rivers_props',jules_message)
        CALL log_fatal("init_rivers_props",                                   &
                     "Error in routing - non-regular lat/lon river grid " //  &
                     "defined but no grid dimension size provided.")
    END IF
        
  END IF
   
!-----------------------------------------------------------------------------
! Check if routing and model grids are same (required for RFM implementation)
! Note some lat/lon rotation may still be required
!-----------------------------------------------------------------------------

  IF (rivers_type == rivers_umtrip) THEN

    IF (ABS(rivers_dlat) /= 1.0 .OR. ABS(rivers_dlon) /= 1.0) THEN ! NOT 1.deg
       WRITE(jules_message,*) &
             "ERROR: Dlat and Dlon must be 1degree if using UM TRIP " //      &
             "       implementation - assumed grid areas in use.    " //      &
             "       Try using rivers_type == rivers_trip instead?"
       CALL jules_print('init_rivers_props',jules_message)
       CALL log_fatal("init_rivers_props",                                    &
             "Error in routing - latitude and longitude grid steps " //       &
             TRIM(to_string(rivers_dlat)) // " , rivers_dlon =" //            &
             TRIM(to_string(rivers_dlon)))
     END IF
  END IF

  IF (rivers_type == rivers_umrfm .OR. rivers_type == rivers_umtrip) THEN   
  
     IF (rivers_reglatlon) THEN
     
       DO l=1,land_pts
         IF (MOD(latitude(l,1)-rivers_lat1, ABS(rivers_dlat)) /= 0 .OR.       &
              MOD(longitude(l,1)-rivers_lon1, ABS(rivers_dlon)) /= 0) THEN
            WRITE(jules_message,*)   &
                 "Routing grid points do not coincide with main model grid."
            CALL jules_print('init_rivers_props',jules_message)
            CALL log_fatal("init_rivers_props",                               &
             "Error in routing - latitude and longitude grid steps " //       &
             " must be equal if set as regular grid. rivers_dlat =" //        &
             TRIM(to_string(rivers_dlat)) // " , rivers_dlon =" //            &
             TRIM(to_string(rivers_dlon)))
  
         END IF

       END DO ! land_pts
     END IF

     rivers_regrid = .False.
     WRITE(jules_message,*)"UMRFM: setting riversRegrid options OFF -  " //   &
                            " assuming coincident grid"
     CALL jules_print('init_rivers_props',jules_message)

  END IF

  IF (rivers_reglatlon) THEN
    ! n.b. in parallel mode, this check is not comprehensive (there might be
    ! a different gap around the edges of the regions given to each task)
    DO l=1,land_pts
       IF (MOD(latitude(l,1)-reg_lat1, ABS(reg_dlat)) /= 0 .OR.               &
                & mod(longitude(l,1)-reg_lon1, ABS(reg_dlon)) /= 0) THEN
          WRITE(jules_message,*)"Land grid points do not coincide with " //   & 
                                "calculated main model grid."
          CALL jules_print('init_rivers_props',jules_message)
          CALL log_fatal("init_rivers_props",                                 &
             "Error in routing - model latitude and longitude grid steps " // &
             " must be equal if set as regular grid. reg_dlat =" //           &
             TRIM(to_string(reg_dlat)) // " , reg_dlon =" //                  &
             TRIM(to_string(reg_dlon)))
       END IF
    END DO
  END IF

  ! Checks for compatible dimension sizes if not regridding

  IF (.NOT. rivers_regrid) THEN
     
     IF (nx_rivers /= nx_grid .OR. ny_rivers /= ny_grid) THEN
          CALL jules_print('init_rivers_props',jules_message)
          CALL log_fatal("init_rivers_props",                                 &
             "Error in routing - model grid and river grid dimensions " //    &
             " must be equal size if no regridding attempted.")
     END IF
  END IF

!-----------------------------------------------------------------------------
! Initialise river gridbox areas (m2) and set grid spacing (m)
!-----------------------------------------------------------------------------

    IF ( rivers_dlat == 0 ) THEN
         rivers_dlat = rivers_lat2d(2,2)-rivers_lat2d(1,1)
    END IF
    IF ( rivers_dlon == 0 ) THEN
         rivers_dlon = rivers_lon2d(2,2)-rivers_lon2d(1,1) 
    END IF    

    IF ( rivers_reglatlon ) THEN
     
        IF ( rivers_dx < 0 ) THEN
            rivers_dx = planet_radius * (ABS(rivers_dlat)*pi_over_180) 
        END IF
 
        DO iy=1,ny_rivers
          DO ix=1,nx_rivers
             rivers_boxareas(ix,iy) = ABS(rivers_earth_area(                  &
                 rivers_lat2d(ix,iy)-0.5*rivers_dlat,                         &
                 rivers_lat2d(ix,iy)+0.5*rivers_dlat,                         &
                 rivers_lon2d(ix,iy)-0.5*rivers_dlon,                         &
                 rivers_lon2d(ix,iy)+0.5*rivers_dlon )) 

        END DO
      END DO

      DO ip=1,np_rivers
         CALL rivers_get_xy_pos(rivers_index_rp(ip),nx_rivers,ny_rivers,ix,iy)
         rivers_boxareas_rp(ip) = ABS(rivers_earth_area(                      &
                 rivers_lat2d(ix,iy)-0.5*rivers_dlat,                         &
                 rivers_lat2d(ix,iy)+0.5*rivers_dlat,                         &
                 rivers_lon2d(ix,iy)-0.5*rivers_dlon,                         &
                 rivers_lon2d(ix,iy)+0.5*rivers_dlon )) 
      END DO

    ELSE

      rivers_boxareas(:,:) = rivers_dx * rivers_dx
      rivers_boxareas_rp(:) = rivers_dx * rivers_dx

    ENDIF

!-----------------------------------------------------------------------------
! Reset saved JULES grid from river to full land model
!-----------------------------------------------------------------------------

   dummy_grid = grid_create(local_grid%is_1d,local_grid%dim_name,             &
                            local_grid%nx,local_grid%x_name,local_grid%nx,    &
                            local_grid%y_name,local_grid%ny)
   use_subgrid = use_sub_local

 END IF ! (check global_landpoint > 1)

 END SUBROUTINE init_rivers_props
#endif
