module NoahmpMathHostMod

!!! Host-side counterparts to the names NoahmpAccDeviceMathShimMod provides.
!!!
!!! Only `pow` is needed. Every other name the shim exports (exp, log, log10,
!!! sqrt, tanh, tan, atan, acos, cos) is a Fortran intrinsic, so dropping the
!!! shim import on the host path resolves them automatically. `pow` is a
!!! C-ism with no Fortran equivalent, so without this module the ~12 files
!!! that call pow() fail to compile whenever NOAHMP_ACC_COLUMNS is undefined
!!! -- i.e. the plain CPU build of this submodule is broken without it.
!!!
!!! This is the exact operation, not an approximation: acc_powf computes
!!! exp(y*log(x)) through two more approximations and returns 0 for x <= 0,
!!! whereas this is the compiler's own `**`.

  use Machine, only : kind_noahmp

  implicit none
  private
  public :: pow

contains

  pure function pow(x, y) result(res)
    real(kind=kind_noahmp), intent(in) :: x, y
    real(kind=kind_noahmp)             :: res
    res = x ** y
  end function pow

end module NoahmpMathHostMod
