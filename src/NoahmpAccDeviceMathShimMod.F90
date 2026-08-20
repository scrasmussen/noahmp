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

  ! exp(x) in DOUBLE. The previous single-precision version had two accuracy
  ! defects that this fixes:
  !   1. Range reduction subtracted a float-rounded ln2 repeatedly in a loop, so
  !      each subtraction injected ~1e-8 of error and they accumulated.
  !      Here n is chosen directly (n = nint(x/ln2)) and ln2 is subtracted as a
  !      Cody-Waite hi/lo pair, so n*ln2 is removed with no rounding error.
  !   2. Only 8 Taylor terms. At the reduced-range edge |r| <= ln2/2 = 0.3466
  !      the first omitted term was r**9/9! ~ 2.6e-10 -- visible in float.
  !      16 terms leave r**17/17! ~ 1e-21, far below double epsilon.
  real(c_double) function acc_exp_d(x)
!$acc routine seq
    real(c_double), value :: x
    integer :: n, k, m
    real(c_double) :: r, term, sum, p, b
    real(c_double), parameter :: invln2 = 1.4426950408889634_c_double
    real(c_double), parameter :: ln2_hi = 6.93147180369123816490e-01_c_double
    real(c_double), parameter :: ln2_lo = 1.90821492927058770002e-10_c_double

    ! Clamp to +/-150. Every caller (acc_expf, acc_powf, acc_tanhf) returns a
    ! FLOAT, and float exp overflows by x ~ 88.7 and underflows by x ~ -103, so
    ! this loses nothing. It also bounds |n| <= 217, which keeps the squaring
    ! loop below to 8 iterations with b peaking at 2**128 -- comfortably inside
    ! double range.
    r = x
    if (r >  150.0_c_double) r =  150.0_c_double
    if (r < -150.0_c_double) r = -150.0_c_double

    n = nint(r * invln2)
    ! Two-part subtraction: (x - n*ln2_hi) is exact-ish, then remove the tail.
    r = (r - real(n, c_double) * ln2_hi) - real(n, c_double) * ln2_lo

    sum  = 1.0_c_double
    term = 1.0_c_double
    do k = 1, 16
       term = term * r / real(k, c_double)
       sum  = sum + term
    end do

    ! Scale by 2**n WITHOUT using `2.0_c_double ** n`.
    !
    ! That expression looks like plain integer exponentiation, but Cray emits a
    ! call to the helper `_RTOI` (real-to-integer power), which nvlink cannot
    ! resolve on the device:
    !     nvlink error : Undefined reference to '_RTOI'
    ! Same family as the _HEXP/_HLOG failures for `x**y` with a runtime real
    ! exponent -- see the intrinsic-availability table in CLAUDE.md. Assuming
    ! the compiler would expand it inline was wrong.
    !
    ! Binary exponentiation instead: <= 8 iterations, all plain multiplies.
    p = 1.0_c_double
    if (n >= 0) then
       b = 2.0_c_double
       m = n
    else
       b = 0.5_c_double
       m = -n
    endif
    do while (m > 0)
       if (mod(m, 2) == 1) p = p * b
       b = b * b
       m = m / 2
    end do

    acc_exp_d = sum * p
  end function acc_exp_d

  real(c_float) function acc_expf(x)
!$acc routine seq
    real(c_float), value :: x
    acc_expf = real(acc_exp_d(real(x, c_double)), c_float)
  end function acc_expf

  ! log(x) in DOUBLE. Improvements over the single-precision version:
  !   1. Mantissa reduced to [1/sqrt2, sqrt2) instead of [0.75, 1.5). That caps
  !      |z| = |(y-1)/(y+1)| at 0.1716 instead of 0.2, and more importantly makes
  !      the interval symmetric in log space so the series converges uniformly.
  !   2. 11 series terms (through z**21) instead of 7, leaving z**23/23 ~ 1e-19.
  !   3. n*ln2 recombined from a Cody-Waite hi/lo pair rather than one rounded
  !      float constant -- that constant alone cost ~1e-7 relative for large |n|.
  real(c_double) function acc_log_d(x)
!$acc routine seq
    real(c_double), value :: x
    integer :: k, n
    real(c_double) :: y, z, z2, term, sum
    real(c_double), parameter :: ln2_hi = 6.93147180369123816490e-01_c_double
    real(c_double), parameter :: ln2_lo = 1.90821492927058770002e-10_c_double
    real(c_double), parameter :: sqrt2  = 1.4142135623730951_c_double
    real(c_double), parameter :: isqrt2 = 0.7071067811865476_c_double

    if (x <= 0.0_c_double) then
       acc_log_d = -700.0_c_double     ! preserves the old sentinel behaviour
       return
    endif

    y = x
    n = 0
    do while (y >= sqrt2)
       y = y * 0.5_c_double
       n = n + 1
    end do
    do while (y < isqrt2)
       y = y * 2.0_c_double
       n = n - 1
    end do

    z    = (y - 1.0_c_double) / (y + 1.0_c_double)
    z2   = z * z
    term = z
    sum  = z
    do k = 3, 21, 2
       term = term * z2
       sum  = sum + term / real(k, c_double)
    end do

    acc_log_d = 2.0_c_double * sum &
              + real(n, c_double) * ln2_hi + real(n, c_double) * ln2_lo
  end function acc_log_d

  real(c_float) function acc_logf(x)
!$acc routine seq
    real(c_float), value :: x
    if (x <= 0.0_c_float) then
       acc_logf = -80.0_c_float        ! unchanged sentinel for the float path
    else
       acc_logf = real(acc_log_d(real(x, c_double)), c_float)
    endif
  end function acc_logf

  real(c_float) function acc_log10f(x)
!$acc routine seq
    real(c_float), value :: x
    real(c_double), parameter :: invln10 = 0.43429448190325182765_c_double
    if (x <= 0.0_c_float) then
       acc_log10f = -80.0_c_float * 0.4342944819032518_c_float
    else
       ! Scale in double before rounding, rather than rounding log(x) to float
       ! first and then multiplying -- that double rounding cost ~1 ulp.
       acc_log10f = real(acc_log_d(real(x, c_double)) * invln10, c_float)
    endif
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

  ! tanh(x). The previous version was the WORST shim in the module at 1.0e-04
  ! worst-case relative error -- ~1000x float epsilon -- and its error peaked at
  ! SMALL |x|, which is precisely the shallow-snow regime FSNO is most sensitive
  ! to (SnowCoverGroundNiu07Mod is the only consumer).
  !
  ! The defect was catastrophic cancellation. It formed e = exp(2x) and then
  ! (e-1)/(e+1). For x = 1e-5, e = 1.00002, so "e - 1" subtracts two nearly
  ! equal floats: the 2e-5 result keeps only ~3 significant digits because
  ! everything above them cancelled.
  !
  ! Two fixes:
  !   * small |x|: use the Maclaurin series, which has no subtraction at all.
  !     tanh(x) = x - x^3/3 + 2x^5/15 - ...; at |x| <= 1e-3 the first omitted
  !     term (17x^7/315) is ~5e-23 relative, i.e. exact in double.
  !   * elsewhere: evaluate (e-1)/(e+1) in DOUBLE. The cancellation still
  !     happens but it removes ~16 digits, not ~7, so the float result is
  !     correctly rounded.
  real(c_float) function acc_tanhf(x)
!$acc routine seq
    real(c_float), value :: x
    real(c_double) :: xd, x2, e

    xd = real(x, c_double)

    if (xd >  20.0_c_double) then
       acc_tanhf =  1.0_c_float          ! 1 - tanh(20) ~ 8e-18, below float eps
    elseif (xd < -20.0_c_double) then
       acc_tanhf = -1.0_c_float
    elseif (abs(xd) <= 1.0e-3_c_double) then
       x2 = xd * xd
       acc_tanhf = real(xd * (1.0_c_double &
                   + x2 * (-1.0_c_double/3.0_c_double &
                   + x2 * ( 2.0_c_double/15.0_c_double))), c_float)
    else
       e = acc_exp_d(2.0_c_double * xd)
       acc_tanhf = real((e - 1.0_c_double) / (e + 1.0_c_double), c_float)
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

  ! pow(x,y) = exp(y*log(x)). This is the highest-leverage shim in the module:
  ! 13 consumers, including RunoffSurfaceXinAnJiangMod, which produces
  ! SFCRUNOFF -- the single largest GPU-vs-CPU difference (498,689 cells), and
  ! the head of the chain SFCRUNOFF -> infxsrt -> sfcheadsubrt (7.2M cells).
  !
  ! The old version compounded three float roundings: log rounded to float, the
  ! product y*log(x) rounded to float, then exp rounded again. The middle one
  ! dominated -- an absolute error d in the exponent becomes a RELATIVE error of
  ! ~d in the result, so for |y*log(x)| ~ 10 a float product carried ~1e-6
  ! relative error into every result. Measured worst case was 8.5e-07, ~7x
  ! float epsilon.
  !
  ! Keeping the entire chain in double leaves one rounding, at the end.
  real(c_float) function acc_powf(x, y)
!$acc routine seq
    real(c_float), value :: x, y
    if (x <= 0.0_c_float) then
       acc_powf = 0.0_c_float          ! unchanged: callers rely on 0**y = 0
    else
       acc_powf = real(acc_exp_d(real(y, c_double) * &
                                 acc_log_d(real(x, c_double))), c_float)
    endif
  end function acc_powf

end module NoahmpAccDeviceMathShimMod
