module Machine

!!! define machine-related constants and parameters
!!! To define real data type precision, use "-DOUBLE_PREC" in CPPFLAG in user_build_options file
!!! By default, Noah-MP uses single precision

! ------------------------ Code history -----------------------------------
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! -------------------------------------------------------------------------
  use netcdf, only : NF90_FILL_FLOAT, NF90_FILL_DOUBLE, NF90_FILL_INT
  implicit none
  save
  private

#ifdef DOUBLE_PREC
  integer, public, parameter :: kind_noahmp = 8 ! double precision
  real(kind=kind_noahmp), public, parameter :: undefined_real = NF90_FILL_DOUBLE
  ! for variable initializatin
#else
  integer, public, parameter :: kind_noahmp = 4 ! single precision
  real(kind=kind_noahmp), public, parameter :: undefined_real = NF90_FILL_FLOAT
  ! for variable initializatin
#endif

  integer,                public, parameter :: undefined_int  = huge(1)    ! undefined integer for variable initialization
  integer,                public, parameter :: undefined_int_neg  = -999   ! undefined integer negative for variable initialization
  real(kind=kind_noahmp), public, parameter :: undefined_real_neg = -999.0 ! undefined real negative for variable initializatin

end module Machine
