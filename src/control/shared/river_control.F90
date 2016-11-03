! *****************************COPYRIGHT****************************************
! (c) Crown copyright, Met Office 2014. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms and
! conditions set out therein.
!
! [Met Office Ref SC0237]
! *****************************COPYRIGHT****************************************
MODULE river_control_mod
  CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RIVER_CONTROL_MOD'

CONTAINS
SUBROUTINE river_control(                                                     &
#if !defined(UM_JULES)
   !INTEGER, INTENT(IN)
   land_pts,                                                                  &
   ! REAL, INTENT(IN)
   sub_surf_roff, surf_roff,                                                  &
   ! LOGICAL, INTENT(IN)
   SRFLOW, SRRUN,                                                             &
   ! REAL, INTENT (OUT)
   RFLOW, RRUN                                                                &
#else
  !LOGICAL, INTENT(IN)                                                        &
  invert_ocean,                                                               &
  !INTEGER, INTENT(IN)                                                        &
  n_proc, land_pts, row_length, rows, river_row_length, river_rows,           &
  land_index, ntype, i_river_vn, aocpl_row_length, aocpl_p_rows, g_p_field,   &
  g_r_field, mype, global_row_length, global_rows, global_river_row_length,   &
  global_river_rows, halo_i, halo_j, model_levels, nsurft,                    &
  !REAL, INTENT(IN)                                                           &
  fqw_surft, delta_lambda, delta_phi, xx_cos_theta_latitude,                  &
  xpa, xua, xva, ypa, yua, yva, flandg, river_vel, river_mcoef, trivdir,      &
  trivseq, r_area, slope, flowobs1, r_inext, r_jnext, r_land, substore,       &
  surfstore, flowin, bflowin, smvcst_gb, smvcwt_gb, surf_roff, sub_surf_roff, &
  frac_surft,                                                                 &
  !INTEGER, INTENT(INOUT)                                                     &
  rivers_count,                                                               &
  !REAL, INTENT(INOUT)                                                        &
  tot_surf_runoff_gb, tot_sub_runoff_gb, acc_lake_evap_gb, twatstor,          &
  soil_layer_moisture, sthu_gb,                                               &
  !LOGICAL, INTENT(OUT)                                                       &
  rivers_call,                                                                &
  !REAL, INTENT(OUT)                                                          &
  inlandout_atm, inlandout_atmos, inlandout_riv, riverout, box_outflow,       &
  box_inflow                                                                  &
#endif
  )

!Module imports

!Common modules
USE riv_intctl_mod_1A,  ONLY: riv_intctl_1A
USE riv_intctl_mod_2A,  ONLY: riv_intctl_2A
USE missing_data_mod, ONLY : rmdi
USE jules_surface_types_mod, ONLY : lake
USE jules_soil_mod, ONLY : sm_levels

!Module imports - Variables required only in UM-mode
#if defined(UM_JULES)
  USE river_inputs_mod,   ONLY: l_rivers, l_inland, river_step
  USE timestep_mod,       ONLY: timestep, timestep_number

  USE level_heights_mod,        ONLY:                                         &
    r_theta_levels

  USE atm_fields_bounds_mod,ONLY: tdims_s, pdims_s, pdims

  USE umPrintMgr
  USE ereport_mod, ONLY : ereport

!Variables required only in JULES standalone-mode
#else

! imported module routines
  USE jules_rivers_trip_mod, ONLY: rivers_drive_trip, regrid_routestore
  USE jules_rivers_rfm_mod, ONLY: rivers_drive_rfm

  USE jules_rivers_mod, ONLY :                                                &
!  imported scalar parameters
     rivers_trip, rivers_rfm, rivers_umtrip, rivers_umrfm                     &
!  imported scalars with intent(in)
    ,rivers_timestep,rivers_type                                              &
    ,np_rivers                                                                &
    ,river_row_length => nx_rivers                                            &
    ,river_rows => ny_rivers                                                  &
    ,rivers_dlat,rivers_dlon,reg_lat1,reg_lon1,reg_dlat,reg_dlon              &
    ,row_length => nx_grid                                                    &
    ,rows => ny_grid                                                          &
!  imported arrays with intent(in)
    ,rivers_lat2d, rivers_lon2d, rivers_boxareas                              &
    ,rivers_index_rp, rivers_next_rp                                          &
!  imported scalars with intent(inout)
    ,rivers_count,rivers_step,rivers_first                                    &
!  imported arrays with intent(inout)
    ,tot_surf_runoff_gb, tot_sub_runoff_gb, acc_lake_evap_gb                  &
    ,rivers_sto_rp,  rivers_dir_rp, rivers_seq_rp, rivers_dra_rp              &
    ,rfm_flowobs1, rfm_inext, rfm_jnext, rfm_substore, rfm_surfstore          &
    ,rfm_flowin, rfm_bflowin, rfm_iarea, rfm_land

  USE jules_riversparm, ONLY : rivers_meander, rivers_speed

  USE model_grid_mod, ONLY : latitude, longitude

  USE coastal, ONLY : flandg

  USE rivers_utils, ONLY :                                                    &
     rivers_map_to_model_grid,                                                &
     rivers_map_from_model_grid,                                              &
     rivers_get_xy_pos

  USE timestep_mod, ONLY : timestep

  USE ancil_info, ONLY : frac_surft, nx=>row_length, ny=>rows, land_index

  USE fluxes, ONLY :  fqw_surft

  USE p_s_parms, ONLY: smvcst_gb, smvcwt_gb, sthu_gb

  USE prognostics, ONLY: smcl_gb

  USE theta_field_sizes, ONLY : t_i_length, t_j_length

  USE c_pi, only: pi_over_180

  USE jules_print_mgr, ONLY :                                                 &
    jules_message,                                                            &
    jules_print

  USE model_grid_mod, ONLY : latitude_of_land_pts, longitude_of_land_pts,     &
    global_land_pts

  USE parallel_mod, ONLY : MASTER_TASK_ID, is_master_task,                    &
    gather_land_field, scatter_land_field

  USE jules_vegetation_mod, ONLY: l_irrig_limit

#endif

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Control routine for rivers
!
! Current Code Owner: Richard Gilham
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

! Subroutine arguments

!-----------------------------------------------------------------------------
! Jules-standalone only arguments
#if !defined(UM_JULES)

! Array arguments with intent(in)
  INTEGER, INTENT(IN) :: land_pts              !  number of land points

! Array arguments with intent(in)
  REAL, INTENT(IN) :: sub_surf_roff(land_pts)  ! Sub-surface runoff (kg m-2 s-1)
  REAL, INTENT(IN) :: surf_roff(land_pts)      ! Surface runoff (kg m-2 s-1)

! Array arguments with intent(inout)
  LOGICAL, INTENT(IN) :: SRFLOW         ! Flag for river flow diagnostic
  LOGICAL, INTENT(IN) :: SRRUN          ! Flag for runoff diagnostic

  REAL, INTENT(OUT) :: RFLOW(land_pts)  ! River flow diagnostic on land points
  REAL, INTENT(OUT) :: RRUN(land_pts)   ! Runoff diagnostic on land points

! Local array variables
  REAL :: surf_runoffin(row_length, rows) ! Total surface runoff on model grid
  REAL :: sub_runoffin(row_length, rows)  ! Total sub surf runoff on model grid

  REAL :: runoff_out(row_length, rows)   ! Total runoff on model grid
  REAL :: rivflow(row_length, rows)      ! River flow on model grid (kg m-2 s-1)

  REAL ::                                                                     &
     RIVEROUT_ATMOS(row_length,rows)                                          &
     ! river flow out from each  gridbox(KG/m2/S)
     ,BOX_OUTFLOW(RIVER_ROW_LENGTH, RIVER_ROWS)                               &
           ! gridbox outflow river grid (Kg/s)
     ,BOX_INFLOW(RIVER_ROW_LENGTH, RIVER_ROWS)                                &
           ! gridbox runoff river grid(Kg/s)
! Declare new variables for inland basin outflow
     ,INLANDOUT_RIV(RIVER_ROW_LENGTH,RIVER_ROWS)                              &
           !RFM OUTFLOW FROM INLAND BASINS ON RFM GRID Kg/s
     ,INLANDOUT_ATMOS(ROW_LENGTH,ROWS)
           ! RFM OUTFLOW FROM  INLAND BASINS ON atmos GRID Kg/m2/s

! Local scalar variables.
  REAL    :: delta_phi      ! RCM gridsize (radians)
  INTEGER :: ip, ix, iy     !  loop counters
  INTEGER :: maxLin         ! an array size

! Local array variables.
  INTEGER :: land_index_grid(land_pts)
            ! New land index array relative to routing grid

  REAL :: acc_evap_grid(row_length, rows)
            !  average lake evaporation since last routing call (kg m-2 s-1)
  REAL :: rivers_dir(river_row_length, river_rows)
            ! Values of river direction on routing grid
  REAL :: rivers_seq(river_row_length, river_rows)
            ! Values of river sequence on routing grid
  REAL :: rivers_sto(river_row_length, river_rows)
            ! Initial water storage
  REAL :: rivers_dra(river_row_length, river_rows)
            ! Values of river drainage area on routing grid

  REAL :: flandg_grid(row_length, rows) !  land fraction on global field
  REAL :: smcl_grid(land_pts)      !  smcl on DSM_LEVELS LAYER on global grid
  REAL :: sthu_grid(land_pts)      !  sthu on DSM_LEVELS LAYER on global grid
  REAL :: smvcst_grid(land_pts)    !  sthu on DSM_LEVELS LAYER on global grid
  REAL :: smvcwt_grid(land_pts)    !  sthu on DSM_LEVELS LAYER on global grid
  REAL :: dummy_xpa(row_length+1)  !  IN Atmosphere longitude coords
  REAL :: dummy_xua(0:row_length)  !  IN Atmosphere UV longitude coords
  REAL :: dummy_xva(row_length+1)  !  IN Atmosphere longitude coords
  REAL :: dummy_ypa(rows)          !  IN Atmosphere latitude coords
  REAL :: dummy_yua(rows)          !  IN Atmosphere latitude coords
  REAL :: dummy_yva(0:rows)        !  IN Atmosphere latitude coords
  REAL :: dummy_r(row_length, rows)
  REAL :: xua(0:row_length)        !  IN Atmosphere UV longitude coords
  REAL :: yva(0:rows)              !  IN Atmosphere latitude coords

  REAL :: dummy_real
  INTEGER :: dummy_int

  INTEGER ::                                                                  &
     global_row_length                                                        &
     , global_rows                                                            &
     , land_points                                                            &
     , global_river_row_length                                                &
     , global_river_rows                                                      &
     , dsm_levels               ! no. of deep soil moisture levels

  REAL ::                                                                     &
     RIV_STEP
                            ! IN river timestep (secs)

  REAL, ALLOCATABLE :: global_lat_of_land_pts(:)
  REAL, ALLOCATABLE :: global_lon_of_land_pts(:)
  REAL, ALLOCATABLE :: global_tot_sub_runoff(:)
  REAL, ALLOCATABLE :: global_tot_surf_runoff(:)
  REAL, ALLOCATABLE :: global_rrun(:)
  REAL, ALLOCATABLE :: global_rflow(:)

!! HL: Original expectation that RFLOW and RRUN would be output on rivers grid.
!!     Updates required to I/O to enable this. Stick with land_pts for now

  LOGICAL :: rivers_call

!-----------------------------------------------------------------------------
!UM-only arguments
#else

LOGICAL, INTENT(IN) ::                                                        &
  invert_ocean

INTEGER, INTENT(IN) ::                                                        &
  n_proc,                                                                     &
  land_pts,                                                                   &
  row_length,                                                                 &
  rows,                                                                       &
  river_row_length,                                                           &
  river_rows,                                                                 &
  land_index(land_pts),                                                       &
  ntype,                                                                      &
  i_river_vn,                                                                 &
  aocpl_row_length,                                                           &
  aocpl_p_rows,                                                               &
  g_p_field,                                                                  &
  g_r_field,                                                                  &
  mype,                                                                       &
  global_row_length,                                                          &
  global_rows,                                                                &
  global_river_row_length,                                                    &
  global_river_rows,                                                          &
  halo_i,                                                                     &
  halo_j,                                                                     &
  model_levels,                                                               &
  nsurft

REAL, INTENT(IN) ::                                                           &
  fqw_surft(land_pts,nsurft),                                                 &
  delta_lambda,                                                               &
  delta_phi,                                                                  &
!  Comment left to save rummaging around for the dims
!  r_theta_levels(1-halo_i:row_length+halo_i,                                 &
!                 1-halo_j:rows+halo_j, 0:model_levels),                      &
  xx_cos_theta_latitude(tdims_s%i_start:tdims_s%i_end,                        &
                        tdims_s%j_start:tdims_s%j_end),                       &
  xpa(aocpl_row_length+1),                                                    &
  xua(0:aocpl_row_length),                                                    &
  xva(aocpl_row_length+1),                                                    &
  ypa(aocpl_p_rows),                                                          &
  yua(aocpl_p_rows),                                                          &
  yva(0:aocpl_p_rows),                                                        &
  flandg(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),        &
  river_vel,                                                                  &
  river_mcoef,                                                                &
  trivdir(river_row_length, river_rows),                                      &
  trivseq(river_row_length, river_rows),                                      &
  r_area(row_length, rows),                                                   &
  slope(row_length, rows),                                                    &
  flowobs1(row_length, rows),                                                 &
  r_inext(row_length, rows),                                                  &
  r_jnext(row_length, rows),                                                  &
  r_land(row_length, rows),                                                   &
  substore(row_length, rows),                                                 &
  surfstore(row_length, rows),                                                &
  flowin(row_length, rows),                                                   &
  bflowin(row_length, rows),                                                  &
  smvcst_gb(land_pts),                                                        &
  smvcwt_gb(land_pts),                                                        &
  surf_roff(land_pts),                                                        &
  sub_surf_roff(land_pts),                                                    &
  frac_surft(land_pts,ntype)

INTEGER, INTENT(INOUT)  ::                                                    &
  rivers_count

REAL, INTENT(INOUT) ::                                                        &
  tot_surf_runoff_gb(land_pts),                                               &
  tot_sub_runoff_gb(land_pts),                                                &
  acc_lake_evap_gb(row_length,rows),                                          &
  twatstor(river_row_length, river_rows),                                     &
  soil_layer_moisture(land_pts,sm_levels),                                    &
  sthu_gb(land_pts,sm_levels)

LOGICAL, INTENT(OUT) ::                                                       &
    rivers_call

REAL, INTENT(OUT) ::                                                          &
  inlandout_atm(land_pts),                                                    &
  inlandout_atmos(row_length,rows),                                           &
  inlandout_riv(river_row_length,river_rows),                                 &
  riverout(row_length, rows),                                                 &
  box_outflow(river_row_length, river_rows),                                  &
  box_inflow(river_row_length, river_rows)

INTEGER, PARAMETER :: i_river_vn_1A = 1
INTEGER, PARAMETER :: i_river_vn_2A = 2
#endif

!Local variables
INTEGER ::                                                                    &
  nstep_rivers,                                                               &
  gather_pe_rivers,                                                           &
  l,i,j

REAL ::                                                                       &
  a_boxareas(row_length,rows)

LOGICAL ::                                                                    &
  first_routing,                                                              &
  invert_atmos

!Error reporting
CHARACTER(LEN=256)       :: message
INTEGER                  :: errorstatus
CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'RIVER_CONTROL'

!Dr Hook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------------
!end of header
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

#if defined(UM_JULES)
  !Set the river routing to run on the 'last' PE as PE0 is very busy
  gather_pe_rivers = n_proc-1

  !Initialise diagnostics on non-river routing timesteps
  riverout    = 0.0
  box_outflow = 0.0
  box_inflow  = 0.0
#endif

  !Initialise the accumulated surface and subsurface runoff to zero
  !at the beginning of river routing timestep
  IF ( rivers_count == 0 ) THEN
    tot_surf_runoff_gb = 0.0
    tot_sub_runoff_gb  = 0.0
    acc_lake_evap_gb   = 0.0
  END IF

  ! Increment counters.
!  rivers_step = rivers_step + 1
  rivers_count = rivers_count + 1
#if defined(UM_JULES)
  nstep_rivers = INT(river_step/timestep)
#else
  nstep_rivers = INT(rivers_timestep)
#endif

  IF (rivers_count == nstep_rivers) THEN
    rivers_call=.TRUE.
  ELSE
    rivers_call=.FALSE.
  END IF

  !Accumulate the runoff as Kg/m2/s over the River Routing period
  DO l = 1, land_pts
    IF(surf_roff(l) <  0.0)THEN
!      WRITE(umMessage,*)'surf_roff(',l,')= ',surf_roff(l)
!      CALL umPrint(umMessage,src='river_control')
    ELSE
      tot_surf_runoff_gb(l) = tot_surf_runoff_gb(l) +                         &
                           (surf_roff(l)/ REAL(nstep_rivers))
    END IF
  END DO
  DO l = 1, land_pts
    IF(sub_surf_roff(l) <  0.0)THEN
!      WRITE(umMessage,*)'sub_surf_roff(',l,')= ',sub_surf_roff(L)
!      CALL umPrint(umMessage,src='river_control')
    ELSE
      tot_sub_runoff_gb(l) = tot_sub_runoff_gb(l) +                           &
                          (sub_surf_roff(l)/ REAL(nstep_rivers))
    END IF
  END DO

# if defined(UM_JULES)
  DO l = 1, land_pts
    j = (land_index(l)-1)/row_length +1
    i = land_index(l) - (j-1)*row_length
    acc_lake_evap_gb(i,j) = acc_lake_evap_gb(i,j) +                           &
                         frac_surft(l,lake)*fqw_surft(l,lake)*timestep
  END DO

  !Detect first entry into river routing
  first_routing = .FALSE.
  IF (timestep_number == nstep_rivers) first_routing = .TRUE.

#else

! Assume inputs to UM routing code on regular grid, irrespective of input
! standalone data (e.g. if provided as vector). Need to rationalise this in
! time...

  IF (row_length /= nx .OR. rows /= ny ) THEN

    DO l=1,land_pts
      ip = land_index(l)
      ix = NINT( (longitude(ip,1) - reg_lon1) / reg_dlon ) + 1
      iy = NINT( (latitude(ip,1) - reg_lat1) / reg_dlat ) + 1

      land_index_grid(l) = (iy - 1) * row_length + ix
      acc_evap_grid(ix,iy) = acc_lake_evap_gb(ip)
      flandg_grid(ix,iy) = flandg(ip,1)

    ENDDO

  ELSE
    land_index_grid(:) = land_index(:)
    flandg_grid(:,:) = flandg(:,:)

    DO l=1,land_pts
      j = (land_index(l)-1)/t_i_length +1
      i = land_index(l) - (j-1)*t_j_length
      acc_evap_grid(i,j) = acc_lake_evap_gb(l)
    ENDDO

  ENDIF

! Map river ancillary information from river vector to river routing grid
  rivers_dir(:,:)=-1.
  rivers_seq(:,:)=-1.
  rivers_sto(:,:)=-1.

  DO ip=1,np_rivers

    CALL rivers_get_xy_pos(rivers_index_rp(ip), river_row_length, river_rows, &
                           ix, iy)
    rivers_dir(ix,iy) = rivers_dir_rp(ip)
    rivers_seq(ix,iy) = rivers_seq_rp(ip)
    rivers_sto(ix,iy) = rivers_sto_rp(ip)
    rivers_dra(ix,iy) = rivers_dra_rp(ip)

  ENDDO

#endif

  IF ( rivers_call ) THEN

    !If ATMOS fields are as Ocean (i.e. inverted NS) set invert_atmos
    invert_atmos = .FALSE.

# if defined(UM_JULES)
    IF (.NOT.invert_ocean) THEN
      invert_atmos = .TRUE.
    ENDIF

    !Calculate the Atmosphere gridbox areas
    DO j = 1, rows
      DO i = 1, row_length
        a_boxareas(i,j) = r_theta_levels(i,j,0)                               &
                      * r_theta_levels(i,j,0)                                 &
                      * delta_lambda * delta_phi                              &
                      * xx_cos_theta_latitude(i,j)
      END DO
    END DO

    SELECT CASE ( i_river_vn )
    CASE ( i_river_vn_1A )
      CALL RIV_INTCTL_1A(                                                     &
        xpa, xua, xva, ypa, yua, yva,                                         &
        g_p_field, g_r_field, n_proc, mype, rmdi,                             &
        gather_pe_rivers,land_pts,land_index,                                 &
        invert_atmos, row_length, rows,                                       &
        global_row_length, global_rows,                                       &
        river_row_length, river_rows,                                         &
        global_river_row_length, global_river_rows,                           &
        flandg(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
        river_step, river_vel, river_mcoef,                                   &
        trivdir, trivseq, twatstor, a_boxareas,                               &
        delta_phi,first_routing,                                              &
        r_area, slope, flowobs1,r_inext,r_jnext,r_land,                       &
        substore,surfstore,flowin,bflowin,                                    &
        !in/out accumulated runoff
        tot_surf_runoff_gb, tot_sub_runoff_gb,                                &
        !out
        box_outflow, box_inflow, riverout,                                    &
        !add inland basin arguments in call to rivctl
        inlandout_atmos,inlandout_riv,                                        &
        !required for soil moisture correction for water conservation
        sm_levels,acc_lake_evap_gb,smvcst_gb,smvcwt_gb,                       &
        soil_layer_moisture(1:,sm_levels),sthu_gb(1:,sm_levels)               &
        )

    CASE ( i_river_vn_2A )
      CALL RIV_INTCTL_2A(                                                     &
        xpa, xua, xva, ypa, yua, yva,                                         &
        g_p_field, g_r_field, n_proc, mype, rmdi,                             &
        gather_pe_rivers,land_pts,land_index,                                 &
        invert_atmos, row_length, rows,                                       &
        global_row_length, global_rows,                                       &
        river_row_length, river_rows,                                         &
        global_river_row_length, global_river_rows,                           &
        flandg(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
        river_step, river_vel, river_mcoef,                                   &
        trivdir, trivseq, twatstor, a_boxareas,                               &
        delta_phi,first_routing,                                              &
        r_area, slope, flowobs1,r_inext,r_jnext,r_land,                       &
        substore,surfstore,flowin,bflowin,                                    &
        !in/out accumulated runoff
        tot_surf_runoff_gb, tot_sub_runoff_gb,                                &
        !out
        box_outflow, box_inflow, riverout,                                    &
        !add inland basin arguments in call to rivctl
        inlandout_atmos,inlandout_riv,                                        &
        !required for soil moisture correction for water conservation
        sm_levels,acc_lake_evap_gb,smvcst_gb,smvcwt_gb,                       &
        soil_layer_moisture(1:,sm_levels),sthu_gb(1:,sm_levels)               &
        )

    CASE DEFAULT

      errorstatus = 10
      WRITE (message,'(A,I6,A)') 'River model type option ',                  &
                         i_river_vn,' not recognised.'

      CALL Ereport ( RoutineName, errorstatus, message)

    END SELECT

    !compress inland basin outputs to land points only
    IF (l_inland) THEN
      DO l = 1,land_pts
        j = (land_index(l)-1)/row_length +1
        i = land_index(l) - (j-1)*row_length
        inlandout_atm(l) = inlandout_atmos(i,j)
      END DO
    END IF

#else
!-------------------------------------------------------------------------------
! STANDALONE versions
!-------------------------------------------------------------------------------

! Translate JULES variables to UM-like structure
    runoff_out(:,:) = 0.0
    rivflow(:,:) = 0.0

    IF ( rivers_dlat == 0 ) rivers_dlat = rivers_lat2d(2,2)-rivers_lat2d(1,1)
    IF ( rivers_dlon == 0 ) rivers_dlon = rivers_lon2d(2,2)-rivers_lon2d(1,1)

    delta_phi = (ABS(rivers_dlat) * pi_over_180)
    riv_step = REAL(rivers_timestep) * timestep

    riverout_atmos(:,:) = 0.0
    box_outflow(:,:) = 0.0
    box_inflow(:,:) = 0.0

    dsm_levels = sm_levels
    smvcst_grid(:) = smvcst_gb(:,dsm_levels)
    smvcwt_grid(:) = smvcwt_gb(:,dsm_levels)
    smcl_grid(:) = smcl_gb(:,dsm_levels)
    sthu_grid(:) = sthu_gb(:,dsm_levels)

    global_row_length = row_length
    global_rows = rows

    global_river_row_length = river_row_length
    global_river_rows = river_rows

!-------------------------------------------------------------------------------
!   Select routing driver to be called
!-------------------------------------------------------------------------------

    SELECT CASE ( rivers_type )

!-------------------------------------------------------------------------------
!   Call TRIP or RFM routing driver
!-------------------------------------------------------------------------------
      CASE ( rivers_trip,  rivers_rfm )

        IF ( is_master_task() ) THEN
          ALLOCATE(global_lat_of_land_pts(global_land_pts))
          ALLOCATE(global_lon_of_land_pts(global_land_pts))
          ALLOCATE(global_tot_sub_runoff(global_land_pts))
          ALLOCATE(global_tot_surf_runoff(global_land_pts))
          ALLOCATE(global_rrun(global_land_pts))
          ALLOCATE(global_rflow(global_land_pts))
        ELSE
          ALLOCATE(global_lat_of_land_pts(1))
          ALLOCATE(global_lon_of_land_pts(1))
          ALLOCATE(global_tot_sub_runoff(1))
          ALLOCATE(global_tot_surf_runoff(1))
          ALLOCATE(global_rrun(1))
          ALLOCATE(global_rflow(1))
        END IF

        CALL gather_land_field(latitude_of_land_pts, global_lat_of_land_pts)
        CALL gather_land_field(longitude_of_land_pts, global_lon_of_land_pts)
        CALL gather_land_field(tot_sub_runoff_gb, global_tot_sub_runoff)
        CALL gather_land_field(tot_surf_runoff_gb, global_tot_surf_runoff)

        IF ( is_master_task() ) THEN
          !  Copy land point arrays to full model fields
          CALL rivers_map_to_model_grid( global_lat_of_land_pts,              &
                                         global_lon_of_land_pts,              &
                                         global_tot_sub_runoff,               &
                                         global_tot_surf_runoff,              &
                                         sub_runoffin, surf_runoffin )

          IF ( rivers_type == rivers_trip) THEN
            CALL rivers_drive_trip( sub_runoffin, surf_runoffin,              &
                                          runoff_out, rivflow )
          ELSE
            CALL rivers_drive_rfm( sub_runoffin, surf_runoffin,               &
                                       runoff_out, rivflow )
          END IF

        END IF

        IF (SRFLOW .OR. SRRUN) THEN
          IF ( is_master_task() ) THEN
            CALL rivers_map_from_model_grid( global_lat_of_land_pts,          &
                                             global_lon_of_land_pts,          &
                                             runoff_out, rivflow,             &
                                             global_rrun, global_rflow )
          END IF
          CALL scatter_land_field(global_rrun, rrun)
          CALL scatter_land_field(global_rflow, rflow)
        END IF

        IF ( l_irrig_limit ) THEN
          CALL regrid_routestore()
        END IF

        DEALLOCATE(global_rflow)
        DEALLOCATE(global_rrun)
        DEALLOCATE(global_tot_surf_runoff)
        DEALLOCATE(global_tot_sub_runoff)
        DEALLOCATE(global_lon_of_land_pts)
        DEALLOCATE(global_lat_of_land_pts)

!-------------------------------------------------------------------------------
!   Call UM-like TRIP routing driver
!-------------------------------------------------------------------------------
      CASE ( rivers_umtrip )

! Set coordinates of edges of model gridboxes.
        DO ix=0,global_row_length
          xua(ix) = reg_lon1 + (REAL(ix)-0.5) * reg_dlon
        END DO
        DO iy=0,global_rows
          yva(iy) = reg_lat1 + (REAL(iy)-0.5) * reg_dlat
        END DO

        CALL riv_intctl_1a ( dummy_xpa, xua, dummy_xva, dummy_ypa, dummy_yua, &
                             yva, dummy_int, dummy_int, dummy_int, dummy_int, &
                             rmdi, dummy_int, land_pts, land_index_grid,      &
                             INVERT_ATMOS, ROW_LENGTH, ROWS,                  &
                             GLOBAL_ROW_LENGTH, GLOBAL_ROWS,                  &
                             RIVER_ROW_LENGTH, RIVER_ROWS,                    &
                             GLOBAL_RIVER_ROW_LENGTH, GLOBAL_RIVER_ROWS,      &
                             flandg_grid, RIV_STEP, rivers_speed,             &
                             rivers_meander, rivers_dir, rivers_seq,          &
                             rivers_sto, rivers_boxareas,                     &
                             dummy_real, rivers_first,                        &
                             dummy_r, dummy_r, dummy_r, dummy_r, dummy_r,     &
                             dummy_r, dummy_r, dummy_r, dummy_r, dummy_r,     &
                             ! IN/OUT accumulated runoff
                             tot_surf_runoff_gb, tot_sub_runoff_gb,           &
                             ! OUT
                             BOX_OUTFLOW, BOX_INFLOW, RIVEROUT_ATMOS          &
                             ! Add new arguments for inland basin outflow
                             ! OUT INLAND BASINS
                             ,INLANDOUT_ATMOS,INLANDOUT_RIV                   &
                             ! Required for soil moisture correction for water
                             ! conservation
                             ,dsm_levels,acc_evap_grid,smvcst_grid,smvcwt_grid&
                             ,smcl_grid,sthu_grid)

! Tidy up standalone river routing variables
        rivers_first = .False.

        rivflow(:,:) = RIVEROUT_ATMOS(:,:)

        DO l=1,land_pts
          ip = land_index(l)
          iy = (land_index_grid(l)-1)/row_length +1
          ix = land_index_grid(l) - (iy-1)*row_length
          runoff_out(ix,iy) = tot_surf_runoff_gb(ip)+ tot_sub_runoff_gb(ip)
        END DO

        ! Update jules_routing module variables
        DO ip=1,np_rivers
          CALL rivers_get_xy_pos( rivers_index_rp(ip), river_row_length,      &
                               river_rows, ix, iy )
          rivers_sto_rp(ip) = rivers_sto(ix,iy)
        END DO

        IF (SRFLOW .OR. SRRUN) THEN
          CALL rivers_map_from_model_grid( latitude_of_land_pts,              &
                                           longitude_of_land_pts,             &
                                           runoff_out, rivflow,               &
                                           rrun, rflow )
        END IF

!-------------------------------------------------------------------------------
!   Call UM-like RFM routing driver
!-------------------------------------------------------------------------------

      CASE ( rivers_umrfm )

        CALL riv_intctl_2a ( dummy_xpa, dummy_xua, dummy_xva,                 &
                             dummy_ypa, dummy_yua, dummy_yva,                 &
                             dummy_int, dummy_int, dummy_int, 1, rmdi,        &
                             1, land_pts, land_index_grid,                    &
                             INVERT_ATMOS, ROW_LENGTH, ROWS,                  &
                             GLOBAL_ROW_LENGTH, GLOBAL_ROWS,                  &
                             RIVER_ROW_LENGTH, RIVER_ROWS,                    &
                             GLOBAL_RIVER_ROW_LENGTH, GLOBAL_RIVER_ROWS,      &
                             flandg_grid, RIV_STEP, rivers_speed,             &
                             rivers_meander, rivers_dir, rivers_seq,          &
                             rivers_sto, rivers_boxareas,                     &
                             delta_phi, rivers_first,                         &
                             rfm_iarea, dummy_r, rfm_flowobs1,                &
                             rfm_inext, rfm_jnext,                            &
                             rfm_land, rfm_substore, rfm_surfstore,           &
                             rfm_flowin,rfm_bflowin,                          &
                       ! IN/OUT accumulated runoff
                             tot_surf_runoff_gb, tot_sub_runoff_gb,           &
!                       ! OUT
                             BOX_OUTFLOW, BOX_INFLOW, RIVEROUT_ATMOS          &
                       ! Add new arguments for inland basin outflow
!                       ! OUT INLAND BASINS  !! NOT USED
                             ,INLANDOUT_ATMOS,INLANDOUT_RIV                   &
!                       ! Required for soil moisture correction for water
                        ! conservation !! NOT USED
                             ,dsm_levels,acc_evap_grid,smvcst_grid,smvcwt_grid&
                             ,smcl_grid,sthu_grid)

   ! Tidy up standalone river routing variables
        rivers_first = .False.

        rivflow(:,:) = RIVEROUT_ATMOS(:,:)

        DO l=1,land_pts
          ip = land_index(l)
          iy = (land_index_grid(l)-1)/row_length +1
          ix = land_index_grid(l) - (iy-1)*row_length
          runoff_out(ix,iy) = tot_surf_runoff_gb(ip)+ tot_sub_runoff_gb(ip)
        END DO

        ! Update jules_routing module variables
        DO ip=1,np_rivers
          CALL rivers_get_xy_pos( rivers_index_rp(ip), river_row_length,      &
                               river_rows, ix, iy )
          rivers_sto_rp(ip) = rivers_sto(ix,iy)
        END DO

        IF (SRFLOW .OR. SRRUN) THEN
          CALL rivers_map_from_model_grid( latitude_of_land_pts,              &
                                           longitude_of_land_pts,             &
                                           runoff_out, rivflow,               &
                                           rrun, rflow )
        END IF

!-------------------------------------------------------------------------------
!       Default case for rivers_type.
!-------------------------------------------------------------------------------
      CASE default

        WRITE(jules_message,*)'ERROR: rivers_drive: ' //                      &
           'do not recognise rivers_type=',TRIM(rivers_type)
        CALL jules_print('rivers_route_drive_standalone',jules_message)

    END SELECT

#endif

!-------------------------------------------------------------------------------
!   Reset counters after a call to routing.
!-------------------------------------------------------------------------------
    !Mult RIVEROUT by the number of physics timesteps per River routing
    !timestep as DAGHYD stores RIVEROUT every timestep. Non-routing
    !timestep vals are passed in as 0.0
    rivers_count = 0

  END IF ! rivers_call

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END SUBROUTINE river_control

END MODULE river_control_mod

