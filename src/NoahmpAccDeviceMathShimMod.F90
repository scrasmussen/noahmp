module NoahmpAccDeviceMathShimMod

!!! Small OpenACC device-side math shims for Cray CCE device helper symbols.
!!! This keeps the experimental NoahMP column offload linkable while the
!!! active physics path is being brought up.

  use, intrinsic :: iso_c_binding, only : c_float, c_int

  implicit none
  private
  public :: acc_expf, acc_logf, acc_log10f, acc_sqrtf, acc_atanf, acc_acosf
  public :: acc_tanf, acc_tanhf, acc_cosf, acc_powf

contains

  real(c_float) function acc_absf(x)
!$acc routine seq
    real(c_float), value :: x
    if (x < 0.0_c_float) then
       acc_absf = -x
    else
       acc_absf = x
    endif
  end function acc_absf

  real(c_float) function acc_sqrtf(x)
!$acc routine seq
    real(c_float), value :: x
    integer :: i
    real(c_float) :: y

    if (x <= 0.0_c_float) then
       acc_sqrtf = 0.0_c_float
       return
    endif

    if (x > 1.0_c_float) then
       y = x
    else
       y = 1.0_c_float
    endif
    do i = 1, 10
       y = 0.5_c_float * (y + x / y)
    enddo
    acc_sqrtf = y
  end function acc_sqrtf

  real(c_float) function acc_expf(x)
!$acc routine seq
    real(c_float), value :: x
    integer :: i, n
    real(c_float) :: r, term, y
    real(c_float), parameter :: ln2 = 0.6931471805599453_c_float

    r = x
    if (r > 80.0_c_float) r = 80.0_c_float
    if (r < -80.0_c_float) r = -80.0_c_float

    n = 0
    do while (r > 0.34657359027997264_c_float)
       r = r - ln2
       n = n + 1
    enddo
    do while (r < -0.34657359027997264_c_float)
       r = r + ln2
       n = n - 1
    enddo

    y = 1.0_c_float
    term = 1.0_c_float
    do i = 1, 8
       term = term * r / real(i, c_float)
       y = y + term
    enddo

    if (n > 0) then
       do i = 1, n
          y = y * 2.0_c_float
       enddo
    elseif (n < 0) then
       do i = 1, -n
          y = y * 0.5_c_float
       enddo
    endif
    acc_expf = y
  end function acc_expf

  real(c_float) function acc_logf(x)
!$acc routine seq
    real(c_float), value :: x
    integer :: i, n
    real(c_float) :: y, z, z2, term, sum
    real(c_float), parameter :: ln2 = 0.6931471805599453_c_float

    if (x <= 0.0_c_float) then
       acc_logf = -80.0_c_float
       return
    endif

    y = x
    n = 0
    do while (y > 1.5_c_float)
       y = y * 0.5_c_float
       n = n + 1
    enddo
    do while (y < 0.75_c_float)
       y = y * 2.0_c_float
       n = n - 1
    enddo

    z = (y - 1.0_c_float) / (y + 1.0_c_float)
    z2 = z * z
    term = z
    sum = term
    do i = 3, 15, 2
       term = term * z2
       sum = sum + term / real(i, c_float)
    enddo
    acc_logf = 2.0_c_float * sum + real(n, c_float) * ln2
  end function acc_logf

  real(c_float) function acc_log10f(x)
!$acc routine seq
    real(c_float), value :: x
    acc_log10f = acc_logf(x) * 0.4342944819032518_c_float
  end function acc_log10f

  real(c_float) function acc_atan_reduced_f(x)
!$acc routine seq
    real(c_float), value :: x
    real(c_float) :: z2, term, sum

    z2 = x * x
    term = x
    sum = term
    term = -term * z2
    sum = sum + term / 3.0_c_float
    term = -term * z2
    sum = sum + term / 5.0_c_float
    term = -term * z2
    sum = sum + term / 7.0_c_float
    term = -term * z2
    sum = sum + term / 9.0_c_float
    term = -term * z2
    sum = sum + term / 11.0_c_float
    acc_atan_reduced_f = sum
  end function acc_atan_reduced_f

  real(c_float) function acc_atanf(x)
!$acc routine seq
    real(c_float), value :: x
    logical :: neg
    real(c_float) :: y
    real(c_float), parameter :: pi2 = 1.5707963267948966_c_float
    real(c_float), parameter :: pi4 = 0.7853981633974483_c_float

    y = x
    neg = .false.
    if (y < 0.0_c_float) then
       y = -y
       neg = .true.
    endif

    if (y > 1.0_c_float) then
       y = pi2 - acc_atan_reduced_f(1.0_c_float / y)
    elseif (y > 0.5_c_float) then
       y = pi4 + acc_atan_reduced_f((y - 1.0_c_float) / (y + 1.0_c_float))
    else
       y = acc_atan_reduced_f(y)
    endif

    if (neg) y = -y
    acc_atanf = y
  end function acc_atanf

  real(c_float) function acc_sinf(x)
!$acc routine seq
    real(c_float), value :: x
    real(c_float) :: y, y2, term, sum
    real(c_float), parameter :: pi = 3.141592653589793_c_float
    real(c_float), parameter :: twopi = 6.283185307179586_c_float

    y = x
    do while (y > pi)
       y = y - twopi
    enddo
    do while (y < -pi)
       y = y + twopi
    enddo
    y2 = y * y
    term = y
    sum = term
    term = -term * y2 / 6.0_c_float
    sum = sum + term
    term = -term * y2 / 20.0_c_float
    sum = sum + term
    term = -term * y2 / 42.0_c_float
    sum = sum + term
    term = -term * y2 / 72.0_c_float
    sum = sum + term
    acc_sinf = sum
  end function acc_sinf

  real(c_float) function acc_cosf(x)
!$acc routine seq
    real(c_float), value :: x
    real(c_float) :: y, y2, term, sum
    real(c_float), parameter :: pi = 3.141592653589793_c_float
    real(c_float), parameter :: twopi = 6.283185307179586_c_float

    y = x
    do while (y > pi)
       y = y - twopi
    enddo
    do while (y < -pi)
       y = y + twopi
    enddo
    y2 = y * y
    term = 1.0_c_float
    sum = term
    term = -term * y2 / 2.0_c_float
    sum = sum + term
    term = -term * y2 / 12.0_c_float
    sum = sum + term
    term = -term * y2 / 30.0_c_float
    sum = sum + term
    term = -term * y2 / 56.0_c_float
    sum = sum + term
    acc_cosf = sum
  end function acc_cosf

  real(c_float) function acc_tanf(x)
!$acc routine seq
    real(c_float), value :: x
    real(c_float) :: c

    c = acc_cosf(x)
    if (acc_absf(c) < 1.0e-6_c_float) then
       if (c < 0.0_c_float) then
          c = -1.0e-6_c_float
       else
          c = 1.0e-6_c_float
       endif
    endif
    acc_tanf = acc_sinf(x) / c
  end function acc_tanf

  real(c_float) function acc_tanhf(x)
!$acc routine seq
    real(c_float), value :: x
    real(c_float) :: e

    if (x > 10.0_c_float) then
       acc_tanhf = 1.0_c_float
    elseif (x < -10.0_c_float) then
       acc_tanhf = -1.0_c_float
    else
       e = acc_expf(2.0_c_float * x)
       acc_tanhf = (e - 1.0_c_float) / (e + 1.0_c_float)
    endif
  end function acc_tanhf

  real(c_float) function acc_acosf(x)
!$acc routine seq
    real(c_float), value :: x
    real(c_float) :: y
    real(c_float), parameter :: pi = 3.141592653589793_c_float
    real(c_float), parameter :: pi2 = 1.5707963267948966_c_float

    y = x
    if (y >= 1.0_c_float) then
       acc_acosf = 0.0_c_float
    elseif (y <= -1.0_c_float) then
       acc_acosf = pi
    else
       acc_acosf = pi2 - acc_atanf(y / acc_sqrtf(max(1.0e-12_c_float, 1.0_c_float - y*y)))
    endif
  end function acc_acosf

  real(c_float) function acc_powif(x, n)
!$acc routine seq
    real(c_float), value :: x
    integer(c_int), value :: n
    integer :: i, m
    real(c_float) :: y

    m = n
    y = 1.0_c_float
    if (m < 0) m = -m
    do i = 1, m
       y = y * x
    enddo
    if (n < 0) y = 1.0_c_float / y
    acc_powif = y
  end function acc_powif

  real(c_float) function acc_powf(x, y)
!$acc routine seq
    real(c_float), value :: x, y
    if (x <= 0.0_c_float) then
       acc_powf = 0.0_c_float
    else
       acc_powf = acc_expf(y * acc_logf(x))
    endif
  end function acc_powf

end module NoahmpAccDeviceMathShimMod
