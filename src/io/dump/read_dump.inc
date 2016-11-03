#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE read_dump(file_name, identifiers)

  USE mpi

  USE io_constants, ONLY : MODE_READ, MAX_DIM_VAR

  USE parallel_mod, ONLY : MASTER_TASK_ID, is_master_task, scatter_land_field

  USE string_utils_mod, ONLY : to_string

  USE file_mod, ONLY : file_handle, file_open, file_introspect,               &
                       file_inquire_dim, file_inquire_var, file_read_var,     &
                       file_close

  USE model_grid_mod, ONLY : global_land_pts

  USE ancil_info, ONLY : land_pts, dim_cs1, frac_surft, nsurft

  USE jules_surface_types_mod, ONLY : npft, ntype, ncpft

  USE imogen_constants, ONLY : n_olevs, nfarray

  USE prognostics, ONLY : canht_pft, canopy_surft, cs_pool_gb, gs_gb, lai_pft,&
                          nsnow_surft, rgrain_surft, rgrainl_surft,           &
                          rho_snow_grnd_surft, sice_surft, sliq_surft,        &
                          snow_grnd_surft, snow_surft, snowdepth_surft,       &
                          t_soil_gb, tsnow_surft, tstar_surft, seed_rain,     &
                          ds_surft, tsoil_deep_gb, tsurf_elev_surft,          &
                          frac_agr_prev_gb,                                   &
                          wood_prod_fast_gb, wood_prod_med_gb,                &
                          wood_prod_slow_gb, n_inorg_gb, ns_pool_gb,          &
                          frac_past_prev_gb

  USE p_s_parms, ONLY : albsoil_gb, bexp_gb, sathh_gb, satcon_gb,             &
                        smvcst_gb, smvccl_gb, smvcwt_gb, hcap_gb,             & 
                        hcon_gb, sthu_gb, clay_gb

  USE crop_vars_mod, ONLY : dvi_cpft, rootc_cpft, harvc_cpft, reservec_cpft, croplai_cpft, cropcanht_cpft,  &
                            sthu_irr_gb, sow_date_cpft, tt_veg_cpft, tt_rep_cpft, irrtiles,     &
                            frac_irr_all_tiles, nirrtile, frac_irr_all

  USE imogen_progs, ONLY : co2_ppmv, co2_change_ppmv, dtemp_o, fa_ocean

  USE trifctl, ONLY : cv_gb, frac_agr_gb

  USE trif_vars_mod, ONLY: frac_past_gb

  USE jules_snow_mod, ONLY : nsmax

  USE jules_soil_mod, ONLY : sm_levels, ns_deep

  USE top_pdm, ONLY : sthzw_gb, zw_gb, fexp_gb, ti_mean_gb, ti_sig_gb

  USE fire_mod,            ONLY : fire_prog

  USE metstats_mod,        ONLY : metstats_prog

  USE jules_rivers_mod, ONLY : np_rivers, rivers_sto_rp

  USE aero,                ONLY : co2_mmr

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Check that the given file is a JULES dump compatible with the current
!   run, and read the given identifiers from it.
!   Note that the reading of the dump is done by the master task and the
!   results scattered to other tasks. This means that dumps written with
!   different amounts of tasks should be interchangable.
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  CHARACTER(len=*) :: file_name  ! The dump file
  CHARACTER(len=*) :: identifiers(:)  ! The model identifiers for the variables
                                      ! to define

! Work variables
  TYPE(file_handle) :: file  ! The opened file

  INTEGER :: nvars  ! The number of variables we are processing

  INTEGER :: dim_size_file  ! The size of the dimension currently being
                            ! processed in the file

! Used when defining dimensions and variables
  INTEGER :: ndims
                            ! The number of dimensions the current variable has
  CHARACTER(len=MAX_SDF_NAME_LEN) :: dim_names(MAX_DIM_VAR)
                            ! The dimension names the current variable should use
  INTEGER :: dim_sizes(MAX_DIM_VAR)
                            ! The dimension sizes for the current variable

  INTEGER :: dim_ids(MAX_DIM_VAR)  ! The dimension ids for the current variable as
                                   ! calculated from file_inquire_dim
  LOGICAL :: is_record_dim  ! Detects if the current dimension is a record dim

  INTEGER :: ndims_file  ! The number of dimensions the variable has in file
                         ! Compared to ndims above for each variable
  INTEGER :: dim_ids_file(MAX_DIM_VAR)
                         ! The ids of the dimensions the variable has in the file
                         ! Compared to dim_ids above to verify the variable has the
                         ! correct dimensions

  INTEGER :: var_ids(SIZE(identifiers))  ! The ids of the variables in the dump file

  LOGICAL :: is_record_var  ! Indicates if a variable uses the record dimension

  LOGICAL :: l_read !Used to bypass dimension check for ancil variables that are
                    !not read in

  REAL :: frac_irr_all_tiles_real
  REAL :: irrtiles_real(npft)
  REAL :: nirrtile_real

! Real versions of integer valued variables
  REAL :: nsnow_real(land_pts, nsurft)
  REAL :: seed_rain_real(SIZE(seed_rain))

  INTEGER :: i, j, m, n  ! Loop counters

  INTEGER :: error  ! Variable to collect MPI errors - most MPI implementations
                    ! bail on error, so this is not checked


! Arrays to hold global land points version of data read in master task
! before scattering
  REAL, ALLOCATABLE :: global_data(:)  ! For data with no vertical levels
  REAL, ALLOCATABLE :: global_data_levs(:,:)   ! With one vertical level
  REAL, ALLOCATABLE :: global_data_snow(:,:,:) ! With two "vertical" levels
                                               ! I.E. snow variables


!-----------------------------------------------------------------------------


  nvars = SIZE(identifiers)

!-----------------------------------------------------------------------------
! In the master task only, we open the file and check that the correct
! dimensions exist and are of a size compatible with this run
!-----------------------------------------------------------------------------
  IF ( is_master_task() ) THEN
!-----------------------------------------------------------------------------
! We use the lowest level file API here, as we don't want to impose the input
! grid
!-----------------------------------------------------------------------------
    file = file_open(file_name, MODE_READ)

! We want to auto-detect the dimensions and variables in the file
    CALL file_introspect(file)

    DO i = 1,nvars

      l_read = .TRUE.

!-----------------------------------------------------------------------------
! Get information about the dimensions used by the variable
!-----------------------------------------------------------------------------
      SELECT CASE ( identifiers(i) )
        CASE ( 'gs', 'sthzw', 'zw', 'cv', 'frac_agr_prev', 'frac_past_prev', &
               'wood_prod_fast', 'wood_prod_med', 'wood_prod_slow',          &
               'n_inorg',                                                    &
               'temp_max_00h_r', 'temp_ave_00h_r', 'prec_tot_00h_r',         &
               'prec_tot_12h_r', 'rhum_min_00h_r', 'dewp_ave_00h_r',         &
               'wind_ave_00h_r', 'temp_max_00h',   'temp_ave_00h',           &
               'temp_pnt_12h',   'prec_tot_00h',   'prec_tot_12h',           &
               'rhum_min_00h',   'rhum_pnt_12h',   'dewp_ave_00h',           &
               'wind_ave_00h',   'wind_pnt_12h',                             &
               'fire_mcarthur_r_dr', 'fire_mcarthur_n_dr',                   &
               'fire_canadian_ffmc', 'fire_canadian_ffmc_mois',              &
               'fire_canadian_dmc',  'fire_canadian_dc',                     &
               'fire_nesterov')
          ndims = 1
          dim_names(1) = LAND_DIM_NAME
          dim_sizes(1) = global_land_pts

        CASE ( 'canht', 'lai' )
          ndims = 2
          dim_names(1:2) = (/ LAND_DIM_NAME, PFT_DIM_NAME /)
          dim_sizes(1:2) = (/ global_land_pts, npft /)

        CASE ( 'cropdvi', 'croprootc', 'cropharvc', 'cropreservec',           &
               'croplai', 'cropcanht' )
          ndims = 2
          dim_names(1:2) = (/ LAND_DIM_NAME, CPFT_DIM_NAME /)
          dim_sizes(1:2) = (/ global_land_pts, ncpft /)

        CASE ( 'cs','ns' )
          ndims = 2
          dim_names(1:2) = (/ LAND_DIM_NAME, SC_POOL_DIM_NAME /)
          dim_sizes(1:2) = (/ global_land_pts, dim_cs1 /)

        CASE ( 'rgrainl', 'snow_ds', 'snow_ice', 'snow_liq', 'tsnow' )
          ndims = 3
          dim_names(1:3) = (/ LAND_DIM_NAME, TILE_DIM_NAME, SNOW_DIM_NAME /)
          dim_sizes(1:3) = (/ global_land_pts, nsurft, nsmax /)

        CASE ( 'sthuf', 't_soil', 'sthu_irr' )
          ndims = 2
          dim_names(1:2) = (/ LAND_DIM_NAME, SOIL_DIM_NAME /)
          dim_sizes(1:2) = (/ global_land_pts, sm_levels /)

        CASE ( 'tsoil_deep' )
          ndims = 2
          dim_names(1:2) = (/ LAND_DIM_NAME, BEDROCK_DIM_NAME /)
          dim_sizes(1:2) = (/ global_land_pts, ns_deep /)

        CASE ( 'canopy', 'nsnow', 'rgrain', 'rho_snow', 'snow_tile',          &
               'snow_depth', 'snow_grnd', 'tstar_tile', 'tsurf_elev_surft' )
          ndims = 2
          dim_names(1:2) = (/ LAND_DIM_NAME, TILE_DIM_NAME /)
          dim_sizes(1:2) = (/ global_land_pts, nsurft /)

! Cases for IMOGEN variables
        CASE ( 'co2_ppmv', 'co2_change_ppmv' )
! scalar variables are represented by an array of dimension 1
          ndims = 1
          dim_names(1) = SCALAR_DIM_NAME
          dim_sizes(1) = 1

        CASE ( 'dtemp_o' )
          ndims = 1
          dim_names(1) = NOLEVS_DIM_NAME
          dim_sizes(1) = n_olevs

        CASE ( 'fa_ocean' )
          ndims = 1
          dim_names(1) = NFARRAY_DIM_NAME
          dim_sizes(1) = nfarray

        CASE ( 'seed_rain' )
          ndims = 1
          dim_names(1) = SEED_DIM_NAME
          dim_sizes(1) = SIZE(seed_rain)

        CASE ( 'rivers_sto_rp' )
          ndims = 1
          dim_names(1) = P_RIVERS_DIM_NAME
          dim_sizes(1) = np_rivers

        !-----------------------------------------------------------------------
        ! Ancilliary variables
        !-----------------------------------------------------------------------

        !Frac ancil namelist
        CASE ( 'frac' )
          IF ( ancil_dump_read%frac ) THEN
            ndims = 2
            dim_names(1:2) = (/ LAND_DIM_NAME, TYPE_DIM_NAME /)
            dim_sizes(1:2) = (/ global_land_pts, ntype /)
          ELSE
            l_read = .FALSE.
          END IF

        !Soil properties ancil namelist
        CASE ( 'b      ', 'sathh  ', 'satcon ', 'sm_sat ', 'sm_crit',         &
               'sm_wilt', 'hcap   ', 'hcon   ' )
          IF ( ancil_dump_read%soil_props ) THEN
            ndims = 2
            dim_names(1:2) = (/ LAND_DIM_NAME, SOIL_DIM_NAME /)
            dim_sizes(1:2) = (/ global_land_pts, sm_levels /)
          ELSE
            l_read = .FALSE.
          END IF

        CASE ( 'albsoil', 'clay   ' )
          IF ( ancil_dump_read%soil_props ) THEN
            ndims = 1
            dim_names(1) = LAND_DIM_NAME
            dim_sizes(1) = global_land_pts
          ELSE
            l_read = .FALSE.
          END IF

        !Topmodel ancil namelist
        CASE ( 'fexp   ', 'ti_mean', 'ti_sig ' )
          IF ( ancil_dump_read%top ) THEN
            ndims = 1
            dim_names(1) = LAND_DIM_NAME
            dim_sizes(1) = global_land_pts
          ELSE
            l_read = .FALSE.
          END IF

        !Agric ancil namelist
        CASE ( 'frac_agr', 'frac_past' )
          IF ( ancil_dump_read%agric ) THEN
            ndims = 1
            dim_names(1) = LAND_DIM_NAME
            dim_sizes(1) = global_land_pts
          ELSE
            l_read = .FALSE.
          END IF

        !Crop props ancillaries namelist
        CASE ( 'cropsowdate', 'cropttveg  ', 'cropttrep  ')
          IF ( ancil_dump_read%crop_props ) THEN
            ndims = 1
            dim_names(1) = LAND_DIM_NAME
            dim_sizes(1) = global_land_pts
          ELSE
            l_read = .FALSE.
          END IF

        !Irrigation ancillaries namelist
        CASE ( 'frac_irrig' )
          IF ( ancil_dump_read%irrig ) THEN
            ndims = 1
            dim_names(1) = LAND_DIM_NAME
            dim_sizes(1) = global_land_pts
          ELSE
            l_read = .FALSE.
          END IF

        CASE ( 'frac_irr_all_tiles' )
          IF ( ancil_dump_read%irrig ) THEN
            !scalar variables are represented by an array of dimension 1
            ndims = 1
            dim_names(1) = SCALAR_DIM_NAME
            dim_sizes(1) = 1
          ELSE
            l_read = .FALSE.
          END IF

        CASE ( 'irrtiles' )
          IF ( ancil_dump_read%irrig ) THEN
            ndims = 1
            dim_names(1) = PFT_DIM_NAME
            dim_sizes(1) = npft
          ELSE
            l_read = .FALSE.
          END IF

        CASE ( 'nirrtile' )
          IF ( ancil_dump_read%irrig ) THEN
            !scalar variables are represented by an array of dimension 1
            ndims = 1
            dim_names(1) = SCALAR_DIM_NAME
            dim_sizes(1) = 1
          ELSE
            l_read = .FALSE.
          END IF

        !CO2 ancil namelist
        CASE ( 'co2_mmr')
          IF ( ancil_dump_read%co2 ) THEN
            !scalar variables are represented by an array of dimension 1
            ndims = 1
            dim_names(1) = SCALAR_DIM_NAME
            dim_sizes(1) = 1
          ELSE
            l_read = .FALSE.
          END IF

        CASE DEFAULT
          CALL log_fatal("read_dump",                                         &
                         "Unrecognised variable for dump - " // TRIM(identifiers(i)))
      END SELECT

!-----------------------------------------------------------------------------
! Check the dimensions exist and have the correct size
!-----------------------------------------------------------------------------
      IF ( l_read ) THEN
        DO j = 1,ndims
  ! Retrive information about the dimension from the file we store the id
  ! for use outside this loop
          CALL file_inquire_dim(                                                &
            file, dim_names(j), dim_ids(j), dim_size_file, is_record_dim        &
          )

  ! Check that we found a dimension
          IF ( dim_ids(j) < 0 )                                                 &
            CALL log_fatal("read_dump",                                         &
                           "Could not find expected dimension '" //             &
                           TRIM(dim_names(j)) // "' in dump file")

  ! Check that the dimension is not a record dimension (there shouldn't be one
  ! in dump files)
          IF ( is_record_dim )                                                  &
            CALL log_fatal("read_dump",                                         &
                           "Dimension '" // TRIM(dim_names(j)) // "' is a " //  &
                           "record dimension - should not exist in dump file")

  ! Check that the dimension has the correct size
          IF ( dim_size_file /= dim_sizes(j) )                                  &
            CALL log_fatal("read_dump",                                         &
                           "Dimension '" // TRIM(dim_names(j)) // "' has " //   &
                           "size incompatible with current run (required: " //  &
                           TRIM(to_string(dim_sizes(j))) // ", found: " //      &
                           TRIM(to_string(dim_size_file)) // ")")
        END DO  ! dims

!-----------------------------------------------------------------------------
! Check that the variable exists and has the correct dimensions
!-----------------------------------------------------------------------------
! Retrieve information about the variable from the file
      CALL file_inquire_var(                                                  &
        file, identifiers(i), var_ids(i), ndims_file, dim_ids_file,           &
        is_record_var                                                         &
      )

! Check that we found a variable
      IF ( var_ids(i) < 1 )                                                   &
        CALL log_fatal("read_dump",                                           &
                       "Failed to find requested variable '" //               &
                       TRIM(identifiers(i)) // "' in dump file")

! Check that the number of dimensions match
      IF ( ndims_file /= ndims )                                              &
        CALL log_fatal("read_dump",                                           &
                       "Variable '" // TRIM(identifiers(i)) // "' has " //    &
                       "incorrect number of dimensions in dump file (" //     &
                       "expected: " // TRIM(to_string(ndims)) // ", " //      &
                       "found: " // TRIM(to_string(ndims_file)) // ")")

! Check that the dimension ids match
      IF ( .NOT. ALL(dim_ids(1:ndims) == dim_ids_file(1:ndims)) )             &
        CALL log_fatal("read_dump",                                           &
                       "Variable '" // TRIM(identifiers(i)) // "' has " //    &
                       "incorrect dimensions in dump file")

      END IF !l_read

    END DO  ! vars

  END IF  ! MASTER TASK


!-----------------------------------------------------------------------------
! Set the requested variables from the file
!
! This is done by reading the value of the variable on global land points
! in the master task, then scattering it to the other tasks
!
! We assume that if the file passed all the checks on dimensions above, then
! it will be fine to fill variables here (i.e. we don't check the dimensions
! associated with the variables)
!-----------------------------------------------------------------------------
! Allocate the global data arrays
  IF ( is_master_task() ) THEN
    ALLOCATE(global_data(global_land_pts))
    ALLOCATE(global_data_levs(global_land_pts, MAX(npft, dim_cs1, sm_levels,  &
                                                   nsurft, ntype, ns_deep)))
    ALLOCATE(global_data_snow(global_land_pts, nsurft, nsmax))
  ELSE
    ALLOCATE(global_data(1))
    ALLOCATE(global_data_levs(1,1))
    ALLOCATE(global_data_snow(1,1,1))
  END IF

  DO i = 1,nvars
!-----------------------------------------------------------------------------
! In the master task, read the global data
!-----------------------------------------------------------------------------
    IF ( is_master_task() ) THEN
      SELECT CASE ( identifiers(i) )
! If it is a land_pts array with no levels associated, read into the
! global_data array
        CASE ( 'gs', 'sthzw', 'zw', 'cv', 'frac_agr_prev', 'frac_past_prev', &
               'wood_prod_fast', 'wood_prod_med', 'wood_prod_slow',          &
               'n_inorg',                                                    &
               'temp_max_00h_r', 'temp_ave_00h_r', 'prec_tot_00h_r',         &
               'prec_tot_12h_r', 'rhum_min_00h_r', 'dewp_ave_00h_r',         &
               'wind_ave_00h_r', 'temp_max_00h',   'temp_ave_00h',           &
               'temp_pnt_12h',   'prec_tot_00h',   'prec_tot_12h',           &
               'rhum_min_00h',   'rhum_pnt_12h',   'dewp_ave_00h',           &
               'wind_ave_00h',   'wind_pnt_12h',                             &
               'fire_mcarthur_r_dr', 'fire_mcarthur_n_dr',                   &
               'fire_canadian_ffmc', 'fire_canadian_ffmc_mois',              &
               'fire_canadian_dmc',  'fire_canadian_dc',                     &
               'fire_nesterov')
          CALL file_read_var(file, var_ids(i), global_data)

! If it is a variable with one vertical level, read the appropriate number of
! levels into global_data_levs
        CASE ( 'canht', 'lai' )
          CALL file_read_var(file, var_ids(i), global_data_levs(:,1:npft))

        CASE ( 'cropdvi', 'croprootc', 'cropharvc', 'cropreservec',           &
               'croplai', 'cropcanht' )
          CALL file_read_var(file, var_ids(i), global_data_levs(:,1:ncpft))

        CASE ( 'cs','ns' )
          CALL file_read_var(file, var_ids(i), global_data_levs(:,1:dim_cs1))

        CASE ( 'sthuf', 't_soil', 'sthu_irr' )
          CALL file_read_var(file, var_ids(i), global_data_levs(:,1:sm_levels))

        CASE ( 'tsoil_deep' )
          CALL file_read_var(file, var_ids(i), global_data_levs(:,1:ns_deep))

        CASE ( 'canopy', 'nsnow', 'rgrain', 'rho_snow', 'snow_tile',          &
               'snow_depth', 'snow_grnd', 'tstar_tile', 'tsurf_elev_surft' )
          CALL file_read_var(file, var_ids(i), global_data_levs(:,1:nsurft))

! If it is a snow variable, first allocate the array if required, then read it
        CASE ( 'rgrainl', 'snow_ds', 'snow_ice', 'snow_liq', 'tsnow' )
          CALL file_read_var(file, var_ids(i), global_data_snow)

! Cases for IMOGEN variables
! Each task runs its own version of IMOGEN - these variables are broadcast to
! all tasks below
        CASE ( 'co2_ppmv' )
          CALL file_read_var(file, var_ids(i), co2_ppmv)

        CASE ( 'co2_change_ppmv' )
          CALL file_read_var(file, var_ids(i), co2_change_ppmv)

        CASE ( 'dtemp_o' )
          CALL file_read_var(file, var_ids(i), dtemp_o)

        CASE ( 'fa_ocean' )
          CALL file_read_var(file, var_ids(i), fa_ocean)

        CASE ( 'seed_rain' )
          CALL file_read_var(file, var_ids(i), seed_rain_real)

        CASE ( 'rivers_sto_rp' )
          CALL file_read_var(file, var_ids(i), rivers_sto_rp)

!-------------------------------------------------------------------------------
! Ancilliary variables
!-------------------------------------------------------------------------------

        !Frac ancil namelist
        CASE ( 'frac' )
          IF ( ancil_dump_read%frac ) THEN
            CALL file_read_var(file, var_ids(i), global_data_levs(:,1:ntype))
          END IF

        !Soil properties ancil namelist
        CASE ( 'b      ', 'sathh  ', 'satcon ', 'sm_sat ', 'sm_crit',         &
               'sm_wilt', 'hcap   ', 'hcon   ' )
          IF ( ancil_dump_read%soil_props ) THEN
            CALL file_read_var(file, var_ids(i), global_data_levs(:,1:sm_levels))
          END IF

        CASE ( 'albsoil', 'clay   ' )
          IF ( ancil_dump_read%soil_props ) THEN
            CALL file_read_var(file, var_ids(i), global_data)
          END IF

        !Topmodel ancil namelist
        CASE ( 'fexp   ', 'ti_mean', 'ti_sig ' )
          IF ( ancil_dump_read%top ) THEN
            CALL file_read_var(file, var_ids(i), global_data)
          END IF

        !Agric ancil namelist
        CASE ( 'frac_agr', 'frac_past' )
          IF ( ancil_dump_read%agric ) THEN
            CALL file_read_var(file, var_ids(i), global_data)
          END IF

        !Crop props ancillaries namelist
        CASE ( 'cropsowdate', 'cropttveg  ', 'cropttrep  ')
          IF ( ancil_dump_read%crop_props ) THEN
            CALL file_read_var(file, var_ids(i), global_data)
          END IF

        !Irrigation ancillaries namelist
        CASE ( 'frac_irrig' )
          IF ( ancil_dump_read%irrig ) THEN
            CALL file_read_var(file, var_ids(i), global_data)
          END IF

        CASE ( 'frac_irr_all_tiles' )
          IF ( ancil_dump_read%irrig ) THEN
            CALL file_read_var(file, var_ids(i), frac_irr_all_tiles_real)
          END IF

        CASE ( 'irrtiles' )
          IF ( ancil_dump_read%irrig ) THEN
            CALL file_read_var(file, var_ids(i), irrtiles_real)
          END IF

        CASE ( 'nirrtile' )
          IF ( ancil_dump_read%irrig ) THEN
            CALL file_read_var(file, var_ids(i), nirrtile_real)
          END IF

        !CO2 ancil namelist
        CASE ( 'co2_mmr' )
          IF ( ancil_dump_read%co2 ) THEN
            CALL file_read_var(file, var_ids(i), co2_mmr)
          END IF

        CASE DEFAULT
          CALL log_fatal("read_dump",                                         &
                         "Unexpected variable in dump - " // TRIM(identifiers(i)))
      END SELECT
    END IF  ! MASTER TASK

! Now scatter the variables into their final destinations
! Note that scatter_land_field can only scatter one land_pts array at a time,
! so to scatter variables with multiple levels we must loop
    SELECT CASE ( identifiers(i) )
      CASE ( 'gs' )
        CALL scatter_land_field(global_data, gs_gb)

      CASE ( 'sthzw' )
        CALL scatter_land_field(global_data, sthzw_gb)

      CASE ( 'zw' )
        CALL scatter_land_field(global_data, zw_gb)

      CASE ( 'cv' )
        CALL scatter_land_field(global_data, cv_gb)

      CASE ( 'frac_agr_prev' )
        CALL scatter_land_field(global_data, frac_agr_prev_gb)

      CASE ( 'frac_past_prev' )
        CALL scatter_land_field(global_data, frac_past_prev_gb)

      CASE ( 'wood_prod_fast' )
        CALL scatter_land_field(global_data, wood_prod_fast_gb)

      CASE ( 'wood_prod_med' )
        CALL scatter_land_field(global_data, wood_prod_med_gb)

      CASE ( 'wood_prod_slow' )
        CALL scatter_land_field(global_data, wood_prod_slow_gb)

      CASE ( 'n_inorg' )
        CALL scatter_land_field(global_data, n_inorg_gb)

      CASE ( 'canht' )
        DO n = 1,npft
          CALL scatter_land_field(global_data_levs(:,n), canht_pft(:,n))
        END DO

      CASE ( 'lai' )
        DO n = 1,npft
          CALL scatter_land_field(global_data_levs(:,n), lai_pft(:,n))
        END DO

      CASE ( 'cropdvi' )
        DO n = 1,ncpft
          CALL scatter_land_field(global_data_levs(:,n), dvi_cpft(:,n))
        END DO

      CASE ( 'croprootc' )
        DO n = 1,ncpft
          CALL scatter_land_field(global_data_levs(:,n), rootc_cpft(:,n))
        END DO

      CASE ( 'cropharvc' )
        DO n = 1,ncpft
          CALL scatter_land_field(global_data_levs(:,n), harvc_cpft(:,n))
        END DO

      CASE ( 'cropreservec' )
        DO n = 1,ncpft
          CALL scatter_land_field(global_data_levs(:,n), reservec_cpft(:,n))
        END DO

      CASE ( 'croplai' )
        DO n = 1,ncpft
          CALL scatter_land_field(global_data_levs(:,n), croplai_cpft(:,n))
        END DO

      CASE ( 'cropcanht' )
        DO n = 1,ncpft
          CALL scatter_land_field(global_data_levs(:,n), cropcanht_cpft(:,n))
        END DO

      CASE ( 'cs' )
        DO n = 1,dim_cs1
          CALL scatter_land_field(global_data_levs(:,n), cs_pool_gb(:,n))
        END DO

      CASE ( 'ns' )
        DO n = 1,dim_cs1
          CALL scatter_land_field(global_data_levs(:,n), ns_pool_gb(:,n))
        END DO

      CASE ( 'sthuf' )
! sthuf is held in sthu until it is processed
        DO n = 1,sm_levels
          CALL scatter_land_field(global_data_levs(:,n), sthu_gb(:,n))
        END DO

      CASE ( 't_soil' )
        DO n = 1,sm_levels
          CALL scatter_land_field(global_data_levs(:,n), t_soil_gb(:,n))
        END DO

      CASE ( 'tsoil_deep' )
        DO n = 1,ns_deep
          CALL scatter_land_field(global_data_levs(:,n), tsoil_deep_gb(:,n))
        END DO

      CASE ( 'sthu_irr' )
        DO n = 1,sm_levels
          CALL scatter_land_field(global_data_levs(:,n), sthu_irr_gb(:,n))
        END DO

      CASE ( 'canopy' )
        DO n = 1,nsurft
          CALL scatter_land_field(global_data_levs(:,n), canopy_surft(:,n))
        END DO

      CASE ( 'nsnow' )
        DO n = 1,nsurft
          CALL scatter_land_field(global_data_levs(:,n), nsnow_real(:,n))
        END DO
        nsnow_surft = NINT(nsnow_real)

      CASE ( 'rgrain' )
        DO n = 1,nsurft
          CALL scatter_land_field(global_data_levs(:,n), rgrain_surft(:,n))
        END DO

      CASE ( 'rho_snow' )
        DO n = 1,nsurft
          CALL scatter_land_field(global_data_levs(:,n), rho_snow_grnd_surft(:,n))
        END DO

      CASE ( 'snow_tile' )
        DO n = 1,nsurft
          CALL scatter_land_field(global_data_levs(:,n), snow_surft(:,n))
        END DO

      CASE ( 'snow_depth' )
        DO n = 1,nsurft
          CALL scatter_land_field(global_data_levs(:,n), snowdepth_surft(:,n))
        END DO

      CASE ( 'snow_grnd' )
        DO n = 1,nsurft
          CALL scatter_land_field(global_data_levs(:,n), snow_grnd_surft(:,n))
        END DO

      CASE ( 'tstar_tile' )
        DO n = 1,nsurft
          CALL scatter_land_field(global_data_levs(:,n), tstar_surft(:,n))
        END DO

      CASE ( 'tsurf_elev_surft' )
        DO n = 1,nsurft
          CALL scatter_land_field(global_data_levs(:,n), tsurf_elev_surft(:,n))
        END DO

      CASE ( 'rgrainl' )
        DO n = 1,nsmax
          DO m = 1,nsurft
            CALL scatter_land_field(global_data_snow(:,m,n), rgrainl_surft(:,m,n))
          END DO
        END DO

      CASE ( 'snow_ds' )
        DO n = 1,nsmax
          DO m = 1,nsurft
            CALL scatter_land_field(global_data_snow(:,m,n), ds_surft(:,m,n))
          END DO
        END DO

      CASE ( 'snow_ice' )
        DO n = 1,nsmax
          DO m = 1,nsurft
            CALL scatter_land_field(global_data_snow(:,m,n), sice_surft(:,m,n))
          END DO
        END DO


      CASE ( 'snow_liq' )
        DO n = 1,nsmax
          DO m = 1,nsurft
            CALL scatter_land_field(global_data_snow(:,m,n), sliq_surft(:,m,n))
          END DO
        END DO


      CASE ( 'tsnow' )
        DO n = 1,nsmax
          DO m = 1,nsurft
            CALL scatter_land_field(global_data_snow(:,m,n), tsnow_surft(:,m,n))
          END DO
        END DO

! Fire and metstats
      CASE ( 'temp_max_00h_r' )
        CALL scatter_land_field(global_data, metstats_prog(:)%temp_max_00h%run)

      CASE ( 'temp_ave_00h_r' )
        CALL scatter_land_field(global_data, metstats_prog(:)%temp_ave_00h%run)

      CASE ( 'prec_tot_00h_r' )
        CALL scatter_land_field(global_data, metstats_prog(:)%prec_tot_00h%run)

      CASE ( 'prec_tot_12h_r' )
        CALL scatter_land_field(global_data, metstats_prog(:)%prec_tot_12h%run)

      CASE ( 'rhum_min_00h_r' )
        CALL scatter_land_field(global_data, metstats_prog(:)%rhum_min_00h%run)

      CASE ( 'dewp_ave_00h_r' )
        CALL scatter_land_field(global_data, metstats_prog(:)%dewp_ave_00h%run)

      CASE ( 'wind_ave_00h_r' )
        CALL scatter_land_field(global_data, metstats_prog(:)%wind_ave_00h%run)

      CASE ( 'temp_max_00h' )
        CALL scatter_land_field(global_data, metstats_prog(:)%temp_max_00h%fin)

      CASE ( 'temp_ave_00h' )
        CALL scatter_land_field(global_data, metstats_prog(:)%temp_ave_00h%fin)

      CASE ( 'temp_pnt_12h' )
        CALL scatter_land_field(global_data, metstats_prog(:)%temp_pnt_12h%fin)

      CASE ( 'prec_tot_00h' )
        CALL scatter_land_field(global_data, metstats_prog(:)%prec_tot_00h%fin)

      CASE ( 'prec_tot_12h' )
        CALL scatter_land_field(global_data, metstats_prog(:)%prec_tot_12h%fin)

      CASE ( 'rhum_min_00h' )
        CALL scatter_land_field(global_data, metstats_prog(:)%rhum_min_00h%fin)

      CASE ( 'rhum_pnt_12h' )
        CALL scatter_land_field(global_data, metstats_prog(:)%rhum_pnt_12h%fin)

      CASE ( 'dewp_ave_00h' )
        CALL scatter_land_field(global_data, metstats_prog(:)%dewp_ave_00h%fin)

      CASE ( 'wind_ave_00h' )
        CALL scatter_land_field(global_data, metstats_prog(:)%wind_ave_00h%fin)

      CASE ( 'wind_pnt_12h' )
        CALL scatter_land_field(global_data, metstats_prog(:)%wind_pnt_12h%fin)

! Fire module variables- land points only
      CASE ( 'fire_mcarthur_r_dr' )
        CALL scatter_land_field(global_data, fire_prog(:)%mcarthur%r_dr)

      CASE ( 'fire_mcarthur_n_dr' )
        CALL scatter_land_field(global_data, fire_prog(:)%mcarthur%n_dr)

      CASE ( 'fire_canadian_ffmc' )
        CALL scatter_land_field(global_data, fire_prog(:)%canadian%ffmc)

      CASE ( 'fire_canadian_ffmc_mois' )
        CALL scatter_land_field(global_data, fire_prog(:)%canadian%ffmc_mois)

      CASE ( 'fire_canadian_dmc' )
        CALL scatter_land_field(global_data, fire_prog(:)%canadian%dmc)

      CASE ( 'fire_canadian_dc' )
        CALL scatter_land_field(global_data, fire_prog(:)%canadian%dc)

      CASE ( 'fire_nesterov' )
        CALL scatter_land_field(global_data, fire_prog(:)%nesterov%findex)

! IMOGEN variables are just broadcast to all tasks
      CASE ( 'co2_ppmv' )
        CALL MPI_BCAST(co2_ppmv, 1, MPI_REAL,                                 &
                       MASTER_TASK_ID, MPI_COMM_WORLD, error)

      CASE ( 'co2_change_ppmv' )
        CALL MPI_BCAST(co2_change_ppmv, 1, MPI_REAL,                          &
                       MASTER_TASK_ID, MPI_COMM_WORLD, error)

      CASE ( 'dtemp_o' )
        CALL MPI_BCAST(dtemp_o, n_olevs, MPI_REAL,                            &
                       MASTER_TASK_ID, MPI_COMM_WORLD, error)

      CASE ( 'fa_ocean' )
        CALL MPI_BCAST(fa_ocean, nfarray, MPI_REAL,                           &
                       MASTER_TASK_ID, MPI_COMM_WORLD, error)

      CASE ( 'seed_rain' )
        CALL MPI_BCAST(seed_rain_real, SIZE(seed_rain_real), MPI_REAL,        &
                       MASTER_TASK_ID, MPI_COMM_WORLD, error)
        seed_rain = NINT(seed_rain_real)

      CASE ( 'rivers_sto_rp' )
        ! nothing to do

!-------------------------------------------------------------------------------
! Ancilliary variables
!-------------------------------------------------------------------------------

      !Frac ancil namelist
      CASE ( 'frac' )
        IF ( ancil_dump_read%frac ) THEN
          DO n = 1,ntype
            CALL scatter_land_field(global_data_levs(:,n), frac_surft(:,n))
          END DO
        END IF

      !Soil properties ancil namelist
      CASE ( 'b      ' )
        IF ( ancil_dump_read%soil_props ) THEN
          DO n = 1,sm_levels
            CALL scatter_land_field(global_data_levs(:,n), bexp_gb(:,n))
          END DO
        END IF

      CASE ( 'sathh  ' )
        IF ( ancil_dump_read%soil_props ) THEN
          DO n = 1,sm_levels
            CALL scatter_land_field(global_data_levs(:,n), sathh_gb(:,n))
          END DO
        END IF

      CASE ( 'satcon ' )
        IF ( ancil_dump_read%soil_props ) THEN
          DO n = 1,sm_levels
            CALL scatter_land_field(global_data_levs(:,n), satcon_gb(:,n))
          END DO
        END IF

      CASE ( 'sm_sat ' )
        IF ( ancil_dump_read%soil_props ) THEN
          DO n = 1,sm_levels
            CALL scatter_land_field(global_data_levs(:,n), smvcst_gb(:,n))
          END DO
        END IF

      CASE ( 'sm_crit' )
        IF ( ancil_dump_read%soil_props ) THEN
          DO n = 1,sm_levels
            CALL scatter_land_field(global_data_levs(:,n), smvccl_gb(:,n))
          END DO
        END IF

      CASE ( 'sm_wilt' )
        IF ( ancil_dump_read%soil_props ) THEN
          DO n = 1,sm_levels
            CALL scatter_land_field(global_data_levs(:,n), smvcwt_gb(:,n))
          END DO
        END IF

      CASE ( 'hcap   ' )
        IF ( ancil_dump_read%soil_props ) THEN
          DO n = 1,sm_levels
            CALL scatter_land_field(global_data_levs(:,n), hcap_gb(:,n))
          END DO
        END IF

      CASE ( 'hcon   ' )
        IF ( ancil_dump_read%soil_props ) THEN
          DO n = 1,sm_levels
            CALL scatter_land_field(global_data_levs(:,n), hcon_gb(:,n))
          END DO
        END IF

      CASE ( 'albsoil' )
        IF ( ancil_dump_read%soil_props ) THEN
          CALL scatter_land_field(global_data, albsoil_gb)
        END IF

      CASE ( 'clay   ' )
        IF ( ancil_dump_read%soil_props ) THEN
          CALL scatter_land_field(global_data, clay_gb)
        END IF

      !Topmodel ancil namelist
      CASE ( 'fexp   ' )
        IF ( ancil_dump_read%top ) THEN
          CALL scatter_land_field(global_data, fexp_gb)
        END IF

      CASE ( 'ti_mean' )
        IF ( ancil_dump_read%top ) THEN
          CALL scatter_land_field(global_data, ti_mean_gb)
        END IF

      CASE ( 'ti_sig ' )
        IF ( ancil_dump_read%top ) THEN
          CALL scatter_land_field(global_data, ti_sig_gb)
        END IF

      !Agric ancil namelist
      CASE ( 'frac_agr' )
        IF ( ancil_dump_read%agric ) THEN
          CALL scatter_land_field(global_data, frac_agr_gb)
        END IF

      CASE ( 'frac_past' )
        IF ( ancil_dump_read%agric ) THEN
          CALL scatter_land_field(global_data, frac_past_gb)
        END IF

      !Crop props ancillaries namelist
      CASE ( 'cropsowdate' )
        IF ( ancil_dump_read%crop_props ) THEN
          CALL scatter_land_field(global_data, sow_date_cpft)
        END IF

      CASE ( 'cropttveg  ' )
        IF ( ancil_dump_read%crop_props ) THEN
          CALL scatter_land_field(global_data, tt_veg_cpft)
        END IF

      CASE ( 'cropttrep  ' )
        IF ( ancil_dump_read%crop_props ) THEN
          CALL scatter_land_field(global_data, tt_rep_cpft)
        END IF

      !Irrigation ancillaries namelist
      CASE ( 'frac_irrig' )
        IF ( ancil_dump_read%irrig ) THEN
          CALL scatter_land_field(global_data, frac_irr_all(:,1))
        END IF

      CASE ( 'frac_irr_all_tiles' )
        IF ( ancil_dump_read%irrig ) THEN
          IF ( frac_irr_all_tiles_real > 0.5 ) THEN
            frac_irr_all_tiles = .TRUE.
          ELSE
            frac_irr_all_tiles = .FALSE.
          END IF
          CALL MPI_BCAST(frac_irr_all_tiles, 1, MPI_LOGICAL,                  &
               MASTER_TASK_ID, MPI_COMM_WORLD, error)
        END IF

      CASE ( 'irrtiles' )
        IF ( ancil_dump_read%irrig ) THEN
          CALL MPI_BCAST(NINT(irrtiles_real), npft, MPI_INTEGER,             &
               MASTER_TASK_ID, MPI_COMM_WORLD, error)
        END IF

      CASE ( 'nirrtile' )
        IF ( ancil_dump_read%irrig ) THEN
          CALL MPI_BCAST(NINT(nirrtile_real), 1, MPI_INTEGER,                 &
               MASTER_TASK_ID, MPI_COMM_WORLD, error)
        END IF

      !CO2 ancil namelist
      CASE ( 'co2_mmr' )
        IF ( ancil_dump_read%co2 ) THEN
          CALL MPI_BCAST(co2_mmr, 1, MPI_REAL,                                &
               MASTER_TASK_ID, MPI_COMM_WORLD, error)
        END IF

      CASE DEFAULT
        CALL log_fatal("read_dump",                                           &
                       "Unexpected variable in dump - " // TRIM(identifiers(i)))
    END SELECT
  END DO

! We are done with the file
  IF ( is_master_task() ) CALL file_close(file)

  IF ( ALLOCATED(global_data) ) DEALLOCATE(global_data)
  IF ( ALLOCATED(global_data_levs) ) DEALLOCATE(global_data_levs)
  IF ( ALLOCATED(global_data_snow) ) DEALLOCATE(global_data_snow)

  RETURN

END SUBROUTINE read_dump
#endif
