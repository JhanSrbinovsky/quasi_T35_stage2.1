#if !defined(UM_JULES)
! Module containing logical switches for diagnostic outputs

  MODULE diag_swchs

  LOGICAL, PARAMETER :: STF_HF_SNOW_MELT=.TRUE.  ! Flag for snowmelt heat flux
  LOGICAL, PARAMETER :: STF_SUB_SURF_ROFF=.TRUE. ! Flag for sub-surface runoff
  LOGICAL, PARAMETER :: SRFLOW = .TRUE.          ! Flag for river outflow
  LOGICAL, PARAMETER :: SRRUN = .TRUE.           ! Flag for runoff after routing

  END MODULE diag_swchs
#endif
