module NoahmpAccDeviceMathShimMod

!!! Small OpenACC device-side math shims for Cray CCE device helper symbols.
!!! This keeps the experimental NoahMP column offload linkable while the
!!! active physics path is being brought up.

  use, intrinsic :: iso_c_binding, only : c_float, c_int, c_double

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

  ! atan(t) for |t| <= tan(pi/12) = 0.26795, in double.
  ! Terms through t**21 leave a first-omitted term of t**23/23 ~ 3e-15 at the
  ! range edge, so the double result is good to ~1e-14 relative and the
  ! single-precision value it is rounded to is correctly rounded.
  real(c_double) function acc_atan_reduced_d(t)
!$acc routine seq
    real(c_double), value :: t
    real(c_double) :: t2, term, sum
    integer :: k

    t2   = t * t
    term = t
    sum  = t
    do k = 3, 21, 2
       term = -term * t2
       sum  = sum + term / real(k, c_double)
    end do
    acc_atan_reduced_d = sum
  end function acc_atan_reduced_d

  ! Two-stage argument reduction, which is what the previous version got
  ! wrong. It reduced |y|>1 via pi/2 - atan(1/y), but for y just above 1 that
  ! leaves 1/y ~ 1, where a 6-term Taylor series is badly unconverged --
  ! atan(1) came out 0.7440 instead of 0.7854, a 5.2% error. Here the first
  ! stage maps |y|>1 into (0,1], and the second maps (tan(pi/12),1] down to
  ! [0,tan(pi/12)] using
  !     atan(y) = pi/6 + atan( (y*sqrt(3) - 1) / (y + sqrt(3)) )
  ! so the series is only ever evaluated on a genuinely small argument.
  real(c_double) function acc_atan_d(x)
!$acc routine seq
    real(c_double), value :: x
    logical :: neg, inv
    real(c_double) :: y, r
    real(c_double), parameter :: pi2   = 1.5707963267948966_c_double
    real(c_double), parameter :: pi6   = 0.5235987755982988_c_double
    real(c_double), parameter :: sqrt3 = 1.7320508075688772_c_double
    real(c_double), parameter :: tan15 = 0.2679491924311227_c_double  ! tan(pi/12)

    y   = x
    neg = .false.
    if (y < 0.0_c_double) then
       y   = -y
       neg = .true.
    endif

    inv = .false.
    if (y > 1.0_c_double) then          ! stage 1: |y| > 1  ->  (0,1]
       y   = 1.0_c_double / y
       inv = .true.
    endif

    if (y > tan15) then                 ! stage 2: (tan(pi/12),1] -> [0,tan(pi/12)]
       r = pi6 + acc_atan_reduced_d((y*sqrt3 - 1.0_c_double) / (y + sqrt3))
    else
       r = acc_atan_reduced_d(y)
    endif

    if (inv) r = pi2 - r
    if (neg) r = -r
    acc_atan_d = r
  end function acc_atan_d

  real(c_float) function acc_atanf(x)
!$acc routine seq
    real(c_float), value :: x
    acc_atanf = real(acc_atan_d(real(x, c_double)), c_float)
  end function acc_atanf

  ! sin(r) and cos(r) for |r| <= pi/4, in double. At the range edge the first
  ! omitted terms are r**17/17! ~ 3e-17 and r**16/16! ~ 4e-16, so both are
  ! accurate to well under single-precision rounding.
  real(c_double) function acc_sin_reduced_d(r)
!$acc routine seq
    real(c_double), value :: r
    real(c_double) :: r2, term, sum
    integer :: k

    r2   = r * r
    term = r
    sum  = r
    do k = 3, 15, 2
       term = -term * r2 / real(k*(k-1), c_double)
       sum  = sum + term
    end do
    acc_sin_reduced_d = sum
  end function acc_sin_reduced_d

  real(c_double) function acc_cos_reduced_d(r)
!$acc routine seq
    real(c_double), value :: r
    real(c_double) :: r2, term, sum
    integer :: k

    r2   = r * r
    term = 1.0_c_double
    sum  = 1.0_c_double
    do k = 2, 16, 2
       term = -term * r2 / real(k*(k-1), c_double)
       sum  = sum + term
    end do
    acc_cos_reduced_d = sum
  end function acc_cos_reduced_d

  ! Quadrant reduction: k = nint(x/(pi/2)), r = x - k*pi/2, so |r| <= pi/4 and
  ! only a short series is needed. pi/2 is split hi+lo (Cody-Waite) so the
  ! subtraction stays accurate for arguments of moderate size.
  !
  ! The previous versions reduced with `do while (y > pi) y = y - twopi` and
  ! then used a 5-term Taylor series valid only near 0. Near the reduced-range
  ! edge (|y| ~ pi) that truncation left cos in error by up to 20% -- e.g. at
  ! x = -33 -- which is why this rewrite matters more than a precision tweak.
  ! The quadrant reduction is written inline in both functions rather than
  ! factored into a subroutine with intent(out) arguments. That shared-helper
  ! version computed correct values at -O0 and when called directly, but Cray
  ! -O2 produced a wrong result for acc_tanf inside a tight loop (sin and cos
  ! effectively disagreeing about the quadrant). Since the device build is
  ! also -O2, the construct is avoided altogether -- a little duplication is
  ! cheaper than a silent wrong answer.
  real(c_float) function acc_sinf(x)
!$acc routine seq
    real(c_float), value :: x
    real(c_double) :: xd, r, fk, s
    integer :: k, q
    real(c_double), parameter :: two_over_pi = 0.6366197723675814_c_double
    real(c_double), parameter :: pio2_hi     = 1.5707963267341256_c_double
    real(c_double), parameter :: pio2_lo     = 6.077100506506192e-11_c_double

    xd = real(x, c_double)
    k  = nint(xd * two_over_pi)
    fk = real(k, c_double)
    r  = (xd - fk*pio2_hi) - fk*pio2_lo
    q  = modulo(k, 4)

    if (q == 0) then
       s =  acc_sin_reduced_d(r)
    else if (q == 1) then
       s =  acc_cos_reduced_d(r)
    else if (q == 2) then
       s = -acc_sin_reduced_d(r)
    else
       s = -acc_cos_reduced_d(r)
    endif
    acc_sinf = real(s, c_float)
  end function acc_sinf

  real(c_float) function acc_cosf(x)
!$acc routine seq
    real(c_float), value :: x
    real(c_double) :: xd, r, fk, c
    integer :: k, q
    real(c_double), parameter :: two_over_pi = 0.6366197723675814_c_double
    real(c_double), parameter :: pio2_hi     = 1.5707963267341256_c_double
    real(c_double), parameter :: pio2_lo     = 6.077100506506192e-11_c_double

    xd = real(x, c_double)
    k  = nint(xd * two_over_pi)
    fk = real(k, c_double)
    r  = (xd - fk*pio2_hi) - fk*pio2_lo
    q  = modulo(k, 4)

    if (q == 0) then
       c =  acc_cos_reduced_d(r)
    else if (q == 1) then
       c = -acc_sin_reduced_d(r)
    else if (q == 2) then
       c = -acc_cos_reduced_d(r)
    else
       c =  acc_sin_reduced_d(r)
    endif
    acc_cosf = real(c, c_float)
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
    real(c_double) :: yd, sd, rd
    real(c_float),  parameter :: pi   = 3.141592653589793_c_float
    real(c_double), parameter :: pid  = 3.141592653589793_c_double
    real(c_double), parameter :: pi2d = 1.5707963267948966_c_double

    ! acos(x) = atan2(sqrt(1-x^2), x), evaluated in double.
    !
    ! The previous form was pi2 - atan(x/sqrt(1-x^2)). That inherited the old
    ! atan's 5% error, and is additionally ill-conditioned near |x| = 1, where
    ! the atan tends to +/-pi/2 and the subtraction cancels. Splitting on the
    ! sign of x keeps the atan argument finite and the result well scaled.
    y = x
    if (y >= 1.0_c_float) then
       acc_acosf = 0.0_c_float
    elseif (y <= -1.0_c_float) then
       acc_acosf = pi
    else
       yd = real(y, c_double)
       sd = sqrt(1.0_c_double - yd*yd)
       if (yd > 0.0_c_double) then
          rd = acc_atan_d(sd / yd)
       elseif (yd < 0.0_c_double) then
          rd = pid - acc_atan_d(sd / (-yd))
       else
          rd = pi2d
       endif
       acc_acosf = real(rd, c_float)
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
