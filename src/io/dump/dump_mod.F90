#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


MODULE dump_mod

  USE io_constants, ONLY : FORMAT_LEN, FORMAT_ASCII, FORMAT_NCDF,             &
                           MAX_SDF_NAME_LEN

  USE logging_mod, ONLY : log_info, log_debug, log_warn, log_error, log_fatal

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Module constants
!-----------------------------------------------------------------------------

  INTEGER, PARAMETER :: MAX_DIM_DUMP = 13
  INTEGER, PARAMETER :: MAX_VAR_DUMP = 80

  CHARACTER(len=MAX_SDF_NAME_LEN), PARAMETER ::                               &
    LAND_DIM_NAME    = "land",                                                &
    PFT_DIM_NAME     = "pft",                                                 &
    CPFT_DIM_NAME    = "cpft",                                                &
    SC_POOL_DIM_NAME = "scpool",                                              &
    SNOW_DIM_NAME    = "snow",                                                &
    SOIL_DIM_NAME    = "soil",                                                &
    TILE_DIM_NAME    = "tile",                                                &
    TYPE_DIM_NAME    = "type",                                                &
    SCALAR_DIM_NAME  = "scalar",                                              &
    NOLEVS_DIM_NAME  = "olevs",                                               &
    NFARRAY_DIM_NAME = "nfarray",                                             &
    SEED_DIM_NAME    = "seed",                                                &
    BEDROCK_DIM_NAME = "bedrock",                                             &
    P_RIVERS_DIM_NAME = "p_rivers"

#if defined(NCDF_DUMMY)
  CHARACTER(len=FORMAT_LEN), PARAMETER :: DUMP_FORMAT = FORMAT_ASCII
#else
  CHARACTER(len=FORMAT_LEN), PARAMETER :: DUMP_FORMAT = FORMAT_NCDF
#endif

  !Create a defined type to store the flags for which ancils are read
  !from the dump file
  TYPE ancil_flags
    LOGICAL :: frac
    LOGICAL :: soil_props
    LOGICAL :: top
    LOGICAL :: agric
    LOGICAL :: crop_props
    LOGICAL :: irrig
    LOGICAL :: rivers_props
    LOGICAL :: co2
  END TYPE

  TYPE(ancil_flags) :: ancil_dump_read

!-----------------------------------------------------------------------------
! Visibility declarations
!-----------------------------------------------------------------------------
  PRIVATE
  PUBLIC  MAX_VAR_DUMP, required_vars_for_configuration,                      &
          read_dump, write_dump,                                              &
          ancil_dump_read

CONTAINS


! Fortran INCLUDE statements would be preferred, but (at least) the pgf90
! compiler objects to their use when the included file contains pre-processor
! directives. At present, such directives are used to exclude files from
! the UM build, so are required. This may change in the future.
#include "required_vars_for_configuration.inc"
#include "read_dump.inc"
#include "write_dump.inc"

END MODULE dump_mod
#endif
