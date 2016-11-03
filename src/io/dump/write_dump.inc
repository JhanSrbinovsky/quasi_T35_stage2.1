#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE write_dump()

  USE io_constants, ONLY : MAX_FILE_NAME_LEN, MAX_DIM_VAR, MODE_WRITE

  USE parallel_mod, ONLY : MASTER_TASK_ID, is_master_task, gather_land_field

  USE model_interface_mod, ONLY : IDENTIFIER_LEN

  USE dictionary_mod, ONLY : dict, dict_create, dict_set, dict_get,           &
                             dict_has_key, dict_free

  USE string_utils_mod, ONLY : to_string

  USE file_mod, ONLY : file_handle, file_open, file_def_dim, file_def_var,    &
                       file_enddef, file_write_var, file_close

  USE output_mod, ONLY : output_dir, run_id

  USE model_time_mod, ONLY : current_time, is_spinup, spinup_cycle

  USE model_grid_mod, ONLY : global_land_pts, latitude, longitude

  USE ancil_info, ONLY : dim_cs1, frac_surft, land_pts, nsurft

  USE jules_surface_types_mod, ONLY : npft, ntype, ncpft

  USE imogen_constants, ONLY : n_olevs, nfarray

  USE prognostics, ONLY : canht_pft, canopy_surft, cs_pool_gb, gs_gb, lai_pft,&
                          nsnow_surft, rgrain_surft, rgrainl_surft,           &
                          rho_snow_grnd_surft, sice_surft, sliq_surft,        &
                          snow_grnd_surft, snow_surft, snowdepth_surft,       &
                          t_soil_gb, tsnow_surft, tstar_surft, seed_rain,     &
                          ds_surft, tsoil_deep_gb, tsurf_elev_surft,          &
                          frac_agr_prev_gb,                                   &
                          wood_prod_slow_gb, wood_prod_med_gb,                &
                          wood_prod_fast_gb, n_inorg_gb, ns_pool_gb,          &
                          frac_past_prev_gb

  USE crop_vars_mod, ONLY : dvi_cpft, rootc_cpft, harvc_cpft, reservec_cpft,  &
                            croplai_cpft, cropcanht_cpft,                     &
                            sthu_irr_gb, sow_date_cpft, tt_veg_cpft,          &
                            tt_rep_cpft, irrtiles,                            &
                            frac_irr_all_tiles, nirrtile, frac_irr_all

  USE imogen_progs, ONLY : co2_ppmv, co2_change_ppmv, dtemp_o, fa_ocean

  USE trifctl, ONLY : cv_gb, frac_agr_gb

  USE trif_vars_mod, ONLY: frac_past_gb

  USE jules_snow_mod, ONLY : nsmax

  USE jules_soil_mod, ONLY : sm_levels, ns_deep

  USE p_s_parms, ONLY : sthu_gb, sthf_gb, albsoil_gb, bexp_gb, sathh_gb,      &
                        satcon_gb, smvcst_gb, smvccl_gb,                      &
                        smvcwt_gb, hcap_gb, hcon_gb, clay_gb

  USE top_pdm, ONLY : sthzw_gb, zw_gb, fexp_gb, ti_mean_gb, ti_sig_gb

  USE fire_mod,            ONLY : fire_prog

  USE metstats_mod,        ONLY : metstats_prog

  USE aero,                ONLY : co2_mmr

  USE jules_rivers_mod, ONLY : np_rivers, rivers_sto_rp

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Writes a dump file for the current timestep
!   Note that the writing of the dump is done by the master task with the
!   values gathered from other tasks. This means that dumps written with
!   different amounts of tasks should be interchangable.
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
  CHARACTER(len=MAX_FILE_NAME_LEN) :: file_name
                                      ! The filename to use for the dump file
  CHARACTER(len=MAX_FILE_NAME_LEN) :: dt_string
                                      ! The datetime string to use in the file
                                      ! name
  CHARACTER(len=IDENTIFIER_LEN) :: identifiers(MAX_VAR_DUMP)
                                      ! The model identifiers for the variables
                                      ! to put in the dump
  CHARACTER(len=IDENTIFIER_LEN) :: vars_from_ancil(MAX_VAR_DUMP)
                               ! The variable identifiers of the ancil
                               ! variables (not used in this subroutine)

  TYPE(file_handle) :: file  ! The dump file

  INTEGER :: nvars  ! The number of variables we are processing
  INTEGER :: nvars_from_ancil

! Variables used when defining dimensions
  TYPE(dict) :: file_dim_ids  ! Dictionary of the dimensions that have been
                              ! defined
                              ! Maps dim_name => dim_id

  INTEGER :: ndims  ! The number of levels dims for the current variable
  CHARACTER(len=MAX_SDF_NAME_LEN) :: dim_names(MAX_DIM_VAR)
                    ! The levels dimension names for the current variable
  INTEGER :: dim_sizes(MAX_DIM_VAR)
                    ! The sizes of the levels dims for the current variable
  INTEGER :: dim_ids(MAX_DIM_VAR)
                   ! The ids in file of the levels dims for the current variable

  INTEGER :: var_ids(MAX_VAR_DUMP)
                        ! The ids of the variables in the dump file

  INTEGER :: i, j, m, n  ! Loop counters

  REAL :: frac_irr_all_tiles_real
  REAL :: irrtiles_real(npft)
  REAL :: nirrtile_real

! Arrays to hold global land points version of data gathered in master task
! before writing
  REAL, ALLOCATABLE :: global_data(:)  ! For data with no vertical levels
  REAL, ALLOCATABLE :: global_data_levs(:,:)   ! With one vertical level
  REAL, ALLOCATABLE :: global_data_snow(:,:,:) ! With two "vertical" levels
                                               ! I.E. snow variables


!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Get the list of identifiers that we are going to output
!-----------------------------------------------------------------------------
  CALL required_vars_for_configuration(nvars, identifiers,                    &
                                       nvars_from_ancil, vars_from_ancil)

! Add latitude and longitude to the list to help offline inspection of the
! dump file, ie they are diagnostics, not prognostics or ancillaries.
! Lat & lon will not be read in from the dump to prevent confusion with the
! model grid namelists

  nvars = nvars + 1
  identifiers(nvars) = 'latitude'
  nvars = nvars + 1
  identifiers(nvars) = 'longitude'

!-----------------------------------------------------------------------------
! In the master task only, we open a new file and define the required
! dimensions and variables
!-----------------------------------------------------------------------------
  IF ( is_master_task() ) THEN
!-----------------------------------------------------------------------------
! Generate the file name that we want to use and open the file
!-----------------------------------------------------------------------------
! File name starts with run id + indicator of a dump file
    file_name = TRIM(run_id) // ".dump."

! Include the current spinup cycle if there is one
    IF ( is_spinup )                                                          &
      file_name = TRIM(file_name) //                                          &
                  "spin" // TRIM(to_string(spinup_cycle)) // "."

! Then current date and time
    WRITE(dt_string, '(I4.4,I2.2,I2.2)') current_time%year,                   &
                                         current_time%month,                  &
                                         current_time%day
    dt_string = TRIM(dt_string) // "." // TRIM(to_string(current_time%time))
    file_name = TRIM(file_name) // TRIM(dt_string)

! Add the extension based on dump format
    SELECT CASE ( DUMP_FORMAT )
      CASE ( FORMAT_ASCII )
        file_name = TRIM(file_name) // ".asc"

      CASE ( FORMAT_NCDF )
        file_name = TRIM(file_name) // ".nc"

      CASE DEFAULT
        CALL log_fatal("write_dump",                                          &
                       "Unrecognised file format - " // TRIM(DUMP_FORMAT))
    END SELECT

! Prepend the output directory
    file_name = TRIM(output_dir) // "/" // TRIM(file_name)

! We use the lowest level file API here, as we don't want to impose a grid
    file = file_open(file_name, MODE_WRITE)

!-----------------------------------------------------------------------------
! Create the dimensions and variables
!-----------------------------------------------------------------------------
    file_dim_ids = dict_create(MAX_DIM_DUMP, INT(1))

    DO i = 1,nvars

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
               'fire_nesterov', 'latitude', 'longitude')
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
          ndims = 2
          dim_names(1:2) = (/ LAND_DIM_NAME, TYPE_DIM_NAME /)
          dim_sizes(1:2) = (/ global_land_pts, ntype /)

        !Soil properties ancil namelist
        CASE ( 'b      ', 'sathh  ', 'satcon ', 'sm_sat ', 'sm_crit',         &
               'sm_wilt', 'hcap   ', 'hcon   ' )
          ndims = 2
          dim_names(1:2) = (/ LAND_DIM_NAME, SOIL_DIM_NAME /)
          dim_sizes(1:2) = (/ global_land_pts, sm_levels /)

        CASE ( 'albsoil', 'clay   ' )
          ndims = 1
          dim_names(1) = LAND_DIM_NAME
          dim_sizes(1) = global_land_pts

        !Topmodel ancillaries namelist
        CASE ( 'fexp   ', 'ti_mean', 'ti_sig ')
          ndims = 1
          dim_names(1) = LAND_DIM_NAME
          dim_sizes(1) = global_land_pts

        !Agric ancillaries namelist
        CASE ( 'frac_agr')
          ndims = 1
          dim_names(1) = LAND_DIM_NAME
          dim_sizes(1) = global_land_pts

        !Agric ancillaries namelist
        CASE ( 'frac_past')
          ndims = 1
          dim_names(1) = LAND_DIM_NAME
          dim_sizes(1) = global_land_pts

        !Crop props ancillaries namelist
        CASE ( 'cropsowdate', 'cropttveg  ', 'cropttrep  ')
          ndims = 1
          dim_names(1) = LAND_DIM_NAME
          dim_sizes(1) = global_land_pts

        !Irrigation ancillaries namelist
        CASE ( 'frac_irrig' )
          ndims = 1
          dim_names(1) = LAND_DIM_NAME
          dim_sizes(1) = global_land_pts

        CASE ( 'frac_irr_all_tiles' )
          !scalar variables are represented by an array of dimension 1
          ndims = 1
          dim_names(1) = SCALAR_DIM_NAME
          dim_sizes(1) = 1

        CASE ( 'irrtiles' )
          ndims = 1
          dim_names(1) = PFT_DIM_NAME
          dim_sizes(1) = npft

        CASE ( 'nirrtile' )
          !scalar variables are represented by an array of dimension 1
          ndims = 1
          dim_names(1) = SCALAR_DIM_NAME
          dim_sizes(1) = 1

        !CO2 ancil namelist
        CASE ( 'co2_mmr' )
          !scalar variables are represented by an array of dimension 1
          ndims = 1
          dim_names(1) = SCALAR_DIM_NAME
          dim_sizes(1) = 1

        CASE DEFAULT
          CALL log_fatal("write_dump",                                        &
                         "Unrecognised variable for dump - " //               &
                         TRIM(identifiers(i)))
      END SELECT

!-----------------------------------------------------------------------------
! Define the dimensions if they have not already been defined
! We use a dictionary to keep track of defined dimension ids
!
! At the same time, gather up the dimension ids needed by the current variable
!-----------------------------------------------------------------------------
      DO j = 1,ndims
! If it has not yet been defined, define the dimension, storing its id
        IF ( .NOT. dict_has_key(file_dim_ids, dim_names(j)) )                 &
          CALL dict_set(                                                      &
            file_dim_ids, dim_names(j),                                       &
            file_def_dim(file, dim_names(j), dim_sizes(j))                    &
          )

! Get the dimension id from the dict and add it to the list for this variable
        CALL dict_get(file_dim_ids, dim_names(j), dim_ids(j))
      END DO

!-----------------------------------------------------------------------------
! Define the variable, saving the id in the file for later
!-----------------------------------------------------------------------------
      var_ids(i) = file_def_var(file, identifiers(i), dim_ids(1:ndims), .FALSE.)

    END DO

!-----------------------------------------------------------------------------
! We have finished defining things
!-----------------------------------------------------------------------------
    CALL file_enddef(file)
    CALL dict_free(file_dim_ids)

  END IF  ! MASTER TASK


!-----------------------------------------------------------------------------
! Gather data from other tasks and write it to file
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
! Gather the variables into a global array to write to file
! Note that gather_land_field can only gather one land_pts array at a time,
! so to gather variables with multiple levels we must loop

    CALL log_info("write_dump", identifiers(i))

    SELECT CASE ( identifiers(i) )
      CASE ( 'gs' )
        CALL gather_land_field(gs_gb, global_data)

      CASE ( 'sthzw' )
        CALL gather_land_field(sthzw_gb, global_data)

      CASE ( 'zw' )
        CALL gather_land_field(zw_gb, global_data)

      CASE ( 'cv' )
        CALL gather_land_field(cv_gb, global_data)

      CASE ( 'frac_agr_prev' )
        CALL gather_land_field(frac_agr_prev_gb, global_data)

      CASE ( 'frac_past_prev' )
        CALL gather_land_field(frac_past_prev_gb, global_data)

      CASE ( 'wood_prod_fast' )
        CALL gather_land_field(wood_prod_fast_gb, global_data)

      CASE ( 'wood_prod_med' )
        CALL gather_land_field(wood_prod_med_gb, global_data)

      CASE ( 'wood_prod_slow' )
        CALL gather_land_field(wood_prod_slow_gb, global_data)

      CASE ( 'latitude' )
        CALL gather_land_field(latitude, global_data)

      CASE ( 'longitude' )
        CALL gather_land_field(longitude, global_data)

      CASE ( 'n_inorg' )
        CALL gather_land_field(n_inorg_gb, global_data)

      CASE ( 'canht' )
        DO n = 1,npft
          CALL gather_land_field(canht_pft(:,n), global_data_levs(:,n))
        END DO

      CASE ( 'lai' )
        DO n = 1,npft
          CALL gather_land_field(lai_pft(:,n), global_data_levs(:,n))
        END DO

      CASE ( 'cs' )
        DO n = 1,dim_cs1
          CALL gather_land_field(cs_pool_gb(:,n), global_data_levs(:,n))
        END DO

      CASE ( 'cropdvi' )
        DO n = 1,ncpft
          CALL gather_land_field(dvi_cpft(:,n), global_data_levs(:,n))
        END DO

      CASE ( 'croprootc' )
        DO n = 1,ncpft
          CALL gather_land_field(rootc_cpft(:,n), global_data_levs(:,n))
        END DO

      CASE ( 'cropharvc' )
        DO n = 1,ncpft
          CALL gather_land_field(harvc_cpft(:,n), global_data_levs(:,n))
        END DO

      CASE ( 'cropreservec' )
        DO n = 1,ncpft
          CALL gather_land_field(reservec_cpft(:,n), global_data_levs(:,n))
        END DO

      CASE ( 'croplai' )
        DO n = 1,ncpft
          CALL gather_land_field(croplai_cpft(:,n), global_data_levs(:,n))
        END DO

      CASE ( 'cropcanht' )
        DO n = 1,ncpft
          CALL gather_land_field(cropcanht_cpft(:,n), global_data_levs(:,n))
        END DO

      CASE ( 'ns' )
        DO n = 1,dim_cs1
          CALL gather_land_field(ns_pool_gb(:,n), global_data_levs(:,n))
        END DO

      CASE ( 'sthuf' )
! sthuf is held in sthu until it is processed
        DO n = 1,sm_levels
          CALL gather_land_field(sthu_gb(:,n) + sthf_gb(:,n),                 &
                                 global_data_levs(:,n))
        END DO

      CASE ( 't_soil' )
        DO n = 1,sm_levels
          CALL gather_land_field(t_soil_gb(:,n), global_data_levs(:,n))
        END DO

      CASE ( 'tsoil_deep' )
        DO n = 1,ns_deep
          CALL gather_land_field(tsoil_deep_gb(:,n), global_data_levs(:,n))
        END DO

      CASE ( 'sthu_irr' )
        DO n = 1,sm_levels
          CALL gather_land_field(sthu_irr_gb(:,n), global_data_levs(:,n))
        END DO

      CASE ( 'canopy' )
        DO n = 1,nsurft
          CALL gather_land_field(canopy_surft(:,n), global_data_levs(:,n))
        END DO

      CASE ( 'nsnow' )
        DO n = 1,nsurft
          CALL gather_land_field(REAL(nsnow_surft(:,n)), global_data_levs(:,n))
        END DO

      CASE ( 'rgrain' )
        DO n = 1,nsurft
          CALL gather_land_field(rgrain_surft(:,n), global_data_levs(:,n))
        END DO

      CASE ( 'rho_snow' )
        DO n = 1,nsurft
          CALL gather_land_field(rho_snow_grnd_surft(:,n),                    &
                                 global_data_levs(:,n))
        END DO

      CASE ( 'snow_tile' )
        DO n = 1,nsurft
          CALL gather_land_field(snow_surft(:,n), global_data_levs(:,n))
        END DO

      CASE ( 'snow_depth' )
        DO n = 1,nsurft
          CALL gather_land_field(snowdepth_surft(:,n), global_data_levs(:,n))
        END DO

      CASE ( 'snow_grnd' )
        DO n = 1,nsurft
          CALL gather_land_field(snow_grnd_surft(:,n), global_data_levs(:,n))
        END DO

      CASE ( 'tstar_tile' )
        DO n = 1,nsurft
          CALL gather_land_field(tstar_surft(:,n), global_data_levs(:,n))
        END DO

      CASE ( 'tsurf_elev_surft' )
        DO n = 1,nsurft
          CALL gather_land_field(tsurf_elev_surft(:,n), global_data_levs(:,n))
        END DO

      CASE ( 'rgrainl' )
        DO n = 1,nsmax
          DO m = 1,nsurft
            CALL gather_land_field(rgrainl_surft(:,m,n),                       &
                                   global_data_snow(:,m,n))
          END DO
        END DO

      CASE ( 'snow_ds' )
        DO n = 1,nsmax
          DO m = 1,nsurft
            CALL gather_land_field(ds_surft(:,m,n), global_data_snow(:,m,n))
          END DO
        END DO

      CASE ( 'snow_ice' )
        DO n = 1,nsmax
          DO m = 1,nsurft
            CALL gather_land_field(sice_surft(:,m,n), global_data_snow(:,m,n))
          END DO
        END DO


      CASE ( 'snow_liq' )
        DO n = 1,nsmax
          DO m = 1,nsurft
            CALL gather_land_field(sliq_surft(:,m,n), global_data_snow(:,m,n))
          END DO
        END DO


      CASE ( 'tsnow' )
        DO n = 1,nsmax
          DO m = 1,nsurft
            CALL gather_land_field(tsnow_surft(:,m,n), global_data_snow(:,m,n))
          END DO
        END DO

      ! Cases for metstats variables
      CASE('temp_max_00h_r')
        CALL gather_land_field(metstats_prog(:)%temp_max_00h%run, global_data)

      CASE('temp_ave_00h_r')
        CALL gather_land_field(metstats_prog(:)%temp_ave_00h%run, global_data)

      CASE('prec_tot_00h_r')
        CALL gather_land_field(metstats_prog(:)%prec_tot_00h%run, global_data)

      CASE('prec_tot_12h_r')
        CALL gather_land_field(metstats_prog(:)%prec_tot_12h%run, global_data)

      CASE('rhum_min_00h_r')
        CALL gather_land_field(metstats_prog(:)%rhum_min_00h%run, global_data)

      CASE('dewp_ave_00h_r')
        CALL gather_land_field(metstats_prog(:)%dewp_ave_00h%run, global_data)

      CASE('wind_ave_00h_r')
        CALL gather_land_field(metstats_prog(:)%wind_ave_00h%run, global_data)

      CASE('temp_max_00h')
        CALL gather_land_field(metstats_prog(:)%temp_max_00h%fin, global_data)

      CASE('temp_ave_00h')
        CALL gather_land_field(metstats_prog(:)%temp_ave_00h%fin, global_data)

      CASE('temp_pnt_12h')
        CALL gather_land_field(metstats_prog(:)%temp_pnt_12h%fin, global_data)

      CASE('prec_tot_00h')
        CALL gather_land_field(metstats_prog(:)%prec_tot_00h%fin, global_data)

      CASE('prec_tot_12h')
        CALL gather_land_field(metstats_prog(:)%prec_tot_12h%fin, global_data)

      CASE('rhum_min_00h')
        CALL gather_land_field(metstats_prog(:)%rhum_min_00h%fin, global_data)

      CASE('rhum_pnt_12h')
        CALL gather_land_field(metstats_prog(:)%rhum_pnt_12h%fin, global_data)

      CASE('dewp_ave_00h')
        CALL gather_land_field(metstats_prog(:)%dewp_ave_00h%fin, global_data)

      CASE('wind_ave_00h')
        CALL gather_land_field(metstats_prog(:)%wind_ave_00h%fin, global_data)

      CASE('wind_pnt_12h')
        CALL gather_land_field(metstats_prog(:)%wind_pnt_12h%fin, global_data)

! Cases for Fire variables
      CASE( 'fire_mcarthur_r_dr' )
        CALL gather_land_field(fire_prog(:)%mcarthur%r_dr, global_data)

      CASE( 'fire_mcarthur_n_dr' )
        CALL gather_land_field(fire_prog(:)%mcarthur%n_dr, global_data)

      CASE( 'fire_canadian_ffmc' )
        CALL gather_land_field(fire_prog(:)%canadian%ffmc, global_data)

      CASE( 'fire_canadian_ffmc_mois' )
        CALL gather_land_field(fire_prog(:)%canadian%ffmc_mois, global_data)

      CASE( 'fire_canadian_dmc' )
        CALL gather_land_field(fire_prog(:)%canadian%dmc, global_data)

      CASE( 'fire_canadian_dc' )
        CALL gather_land_field(fire_prog(:)%canadian%dc, global_data)

      CASE( 'fire_nesterov' )
        CALL gather_land_field(fire_prog(:)%nesterov%findex, global_data)


! Since each task runs its own version of IMOGEN, we just use the values from
! the master task
      CASE ( 'co2_ppmv', 'co2_change_ppmv', 'dtemp_o', 'fa_ocean', 'seed_rain' )
! Nothing to do

      CASE ( 'rivers_sto_rp')
! Nothing to do

!-------------------------------------------------------------------------------
! Ancilliary variables
!-------------------------------------------------------------------------------

      !Frac ancil namelist
      CASE ( 'frac' )
        DO n = 1,ntype
          CALL gather_land_field(frac_surft(:,n), global_data_levs(:,n))
        END DO


      !Soil properties ancil namelist
      CASE ( 'b      ')
        DO n = 1,sm_levels
          CALL gather_land_field(bexp_gb(:,n), global_data_levs(:,n))
        END DO

      CASE ( 'sathh  ')
        DO n = 1,sm_levels
          CALL gather_land_field(sathh_gb(:,n), global_data_levs(:,n))
        END DO

      CASE ( 'satcon ')
        DO n = 1,sm_levels
          CALL gather_land_field(satcon_gb(:,n), global_data_levs(:,n))
        END DO

      CASE ( 'sm_sat ')
        DO n = 1,sm_levels
          CALL gather_land_field(smvcst_gb(:,n), global_data_levs(:,n))
        END DO

      CASE ( 'sm_crit')
        DO n = 1,sm_levels
          CALL gather_land_field(smvccl_gb(:,n), global_data_levs(:,n))
        END DO

      CASE ( 'sm_wilt')
        DO n = 1,sm_levels
          CALL gather_land_field(smvcwt_gb(:,n), global_data_levs(:,n))
        END DO

      CASE ( 'hcap   ')
        DO n = 1,sm_levels
          CALL gather_land_field(hcap_gb(:,n), global_data_levs(:,n))
        END DO

      CASE ( 'hcon   ' )
        DO n = 1,sm_levels
          CALL gather_land_field(hcon_gb(:,n), global_data_levs(:,n))
        END DO

      CASE ( 'albsoil' )
        CALL gather_land_field(albsoil_gb, global_data)

      CASE ( 'clay   ' )
        CALL gather_land_field(clay_gb, global_data)

      !Topmodel ancillaries namelist
      CASE ( 'fexp   ' )
        CALL gather_land_field(fexp_gb, global_data)

      CASE ( 'ti_mean' )
        CALL gather_land_field(ti_mean_gb, global_data)

      CASE ( 'ti_sig ' )
        CALL gather_land_field(ti_sig_gb, global_data)

      !Agric ancillaries namelist
      CASE ( 'frac_agr' )
        CALL gather_land_field(frac_agr_gb, global_data)

      CASE ( 'frac_past' )
        CALL gather_land_field(frac_past_gb, global_data)

      !Crop props ancillaries namelist
      CASE ( 'cropsowdate' )
        CALL gather_land_field(sow_date_cpft, global_data)

      CASE ( 'cropttveg  ' )
        CALL gather_land_field(tt_veg_cpft, global_data)

      CASE ( 'cropttrep  ')
        CALL gather_land_field(tt_rep_cpft, global_data)

      !Irrigation ancillaries namelist
        CASE ( 'frac_irrig' )
          CALL gather_land_field(frac_irr_all(:,1), global_data)

        CASE ( 'frac_irr_all_tiles' )
          ! Convert to a real
          IF ( frac_irr_all_tiles ) THEN
            frac_irr_all_tiles_real = 1.0
          ELSE
            frac_irr_all_tiles_real = 0.0
          END IF

        CASE ( 'irrtiles' )
          irrtiles_real = REAL(irrtiles(1:npft))

        CASE ( 'nirrtile' )
          nirrtile_real = REAL(nirrtile)

        !CO2 ancil namelist
        CASE ( 'co2_mmr' )
          !Nothing to do

      CASE DEFAULT
        CALL log_fatal("write_dump",                                          &
                       "No code to gather variable for dump - " //            &
                       TRIM(identifiers(i)))
    END SELECT

!-----------------------------------------------------------------------------
! In the master task, write the global data to file
!-----------------------------------------------------------------------------
    IF ( is_master_task() ) THEN

      SELECT CASE ( identifiers(i) )
! If it is a land_pts array with no levels associated, write the global_data array
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
               'fire_nesterov', 'latitude', 'longitude')
          CALL file_write_var(file, var_ids(i), global_data)

! If it is a variable with one vertical level, write the appropriate number of
! levels from global_data_levs
        CASE ( 'canht', 'lai' )
          CALL file_write_var(file, var_ids(i), global_data_levs(:,1:npft))

        CASE ( 'cropdvi', 'croprootc', 'cropharvc', 'cropreservec',           &
               'croplai', 'cropcanht' )
          CALL file_write_var(file, var_ids(i), global_data_levs(:,1:ncpft))

        CASE ( 'cs','ns' )
          CALL file_write_var(file, var_ids(i), global_data_levs(:,1:dim_cs1))

        CASE ( 'sthuf', 't_soil', 'sthu_irr' )
          CALL file_write_var(file, var_ids(i), global_data_levs(:,1:sm_levels))

        CASE ( 'tsoil_deep' )
          CALL file_write_var(file, var_ids(i), global_data_levs(:,1:ns_deep))

        CASE ( 'canopy', 'nsnow', 'rgrain', 'rho_snow', 'snow_tile',          &
               'snow_depth', 'snow_grnd', 'tstar_tile', 'tsurf_elev_surft' )
          CALL file_write_var(file, var_ids(i), global_data_levs(:,1:nsurft))

! If it is a snow variable, first allocate the array if required, then read it
        CASE ( 'rgrainl', 'snow_ds', 'snow_ice', 'snow_liq', 'tsnow' )
          CALL file_write_var(file, var_ids(i), global_data_snow)

! Cases for IMOGEN variables
! Each task runs its own version of IMOGEN - we just write the master task's
! versions
        CASE ( 'co2_ppmv' )
          CALL file_write_var(file, var_ids(i), co2_ppmv)

        CASE ( 'co2_change_ppmv' )
          CALL file_write_var(file, var_ids(i), co2_change_ppmv)

        CASE ( 'dtemp_o' )
          CALL file_write_var(file, var_ids(i), dtemp_o)

        CASE ( 'fa_ocean' )
          CALL file_write_var(file, var_ids(i), fa_ocean)

        CASE ( 'seed_rain' )
          CALL file_write_var(file, var_ids(i), REAL(seed_rain))

        CASE ( 'rivers_sto_rp' )
          CALL file_write_var(file, var_ids(i), rivers_sto_rp)

!-------------------------------------------------------------------------------
! Ancilliary variables
!-------------------------------------------------------------------------------

        !Frac ancil namelist
        CASE ( 'frac' )
          CALL file_write_var(file, var_ids(i), global_data_levs(:,1:ntype))

        !Soil properties ancil namelist
        CASE ( 'b      ', 'sathh  ', 'satcon ', 'sm_sat ', 'sm_crit',         &
               'sm_wilt', 'hcap   ', 'hcon   ' )
          CALL file_write_var(file, var_ids(i), global_data_levs(:,1:sm_levels))

        CASE ( 'albsoil', 'clay   ' )
          CALL file_write_var(file, var_ids(i), global_data)

        !Topmodel ancil namelist
        CASE ( 'fexp   ', 'ti_mean', 'ti_sig ' )
          CALL file_write_var(file, var_ids(i), global_data)

        !Agric ancil namelist
        CASE ( 'frac_agr' )
          CALL file_write_var(file, var_ids(i), global_data)

        !Crop props ancillaries namelist
        CASE ( 'cropsowdate', 'cropttveg  ', 'cropttrep  ')
          CALL file_write_var(file, var_ids(i), global_data)

        !Irrigation ancil namelist
        CASE ( 'frac_irrig' )
          CALL file_write_var(file, var_ids(i), global_data)

        CASE ( 'frac_irr_all_tiles' )
          CALL file_write_var(file, var_ids(i), frac_irr_all_tiles_real)

        CASE ( 'irrtiles' )
          CALL file_write_var(file, var_ids(i), irrtiles_real)

        CASE ( 'nirrtile' )
          CALL file_write_var(file, var_ids(i), nirrtile_real)

        !CO2 ancil namelist
        CASE ( 'co2_mmr' )
          CALL file_write_var(file, var_ids(i), co2_mmr)

        CASE DEFAULT
          CALL log_fatal("write_dump",                                        &
                         "Unrecognised variable for dump - " //               &
                         TRIM(identifiers(i)))
      END SELECT
    END IF  ! MASTER TASK
  END DO

! We are done with the file and dictionaries
  IF ( is_master_task() ) CALL file_close(file)

  IF ( ALLOCATED(global_data) ) DEALLOCATE(global_data)
  IF ( ALLOCATED(global_data_levs) ) DEALLOCATE(global_data_levs)
  IF ( ALLOCATED(global_data_snow) ) DEALLOCATE(global_data_snow)

  RETURN

END SUBROUTINE write_dump
#endif
