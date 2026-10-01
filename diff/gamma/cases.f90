! The input sweep for the gamma differential, shared by fdrv.f90 and mirrored
! literal for literal in odrv.ml.  Both drivers write every input out as a bit
! pattern as well as every result, so cmp proves the two sweeps were the same
! sweep rather than only that the answers agreed.
MODULE cases
  USE service
  IMPLICIT NONE

  REAL(DP), PARAMETER :: lgmcx(22) = [ 1.0E3_DP, 2.0E3_DP, 1.0E4_DP, 1.0E5_DP, &
    1.0E6_DP, 1.0E7_DP, 9.0E7_DP, 9.4E7_DP, 9.48E7_DP, 9.49E7_DP, 9.5E7_DP, &
    1.0E8_DP, 1.0E9_DP, 1.0E20_DP, 1.0E100_DP, 1.0E200_DP, 1.0E300_DP, &
    1.0E305_DP, 1.2E306_DP, 1.3E306_DP, 1.0E307_DP, 1.0E308_DP ]

  INTEGER, PARAMETER :: binn(20) = [ 41, 41, 45, 50, 60, 100, 170, 200, 300, &
    500, 1000, 1000, 1000, 1000, 1000, 1000, 100000, 100000, 100000, 100000 ]
  INTEGER, PARAMETER :: binm(20) = [ 20, 21, 22, 25, 30, 50, 85, 100, 150, &
    250, 0, 1, 25, 100, 250, 500, 2, 20, 21, 40 ]

  REAL(DP), PARAMETER :: pcha(27) = [ -3._DP, -5._DP, -5._DP, -5._DP, -10._DP, &
    -10._DP, -15._DP, -15._DP, -15._DP, -15._DP, -15._DP, -20._DP, -21._DP, &
    -30._DP, -30._DP, -50._DP, -50._DP, -100._DP, -170._DP, -5._DP, 0._DP, &
    0._DP, 0._DP, -1._DP, -1._DP, -2._DP, -170._DP ]
  REAL(DP), PARAMETER :: pchx(27) = [ 0._DP, -1._DP, -3._DP, -5._DP, -4._DP, &
    3._DP, -3._DP, -10._DP, 3._DP, 10._DP, 20._DP, -10._DP, -1._DP, &
    -5._DP, 5._DP, -20._DP, 10._DP, -30._DP, -1._DP, 10._DP, 0._DP, &
    -3._DP, 5._DP, -1._DP, 0._DP, -1._DP, 0._DP ]

  REAL(DP), PARAMETER :: poa(34) = [ -170.5_DP, -100.5_DP, -60.25_DP, -40.5_DP, &
    -30.5_DP, -25.75_DP, -21.5_DP, -20.5_DP, -16.5_DP, -15._DP, -13.5_DP, &
    -10.5_DP, -5.5_DP, -2.5_DP, -0.5_DP, -0.25_DP, -0.001_DP, 0._DP, 0.001_DP, &
    0.25_DP, 0.5_DP, 1._DP, 1.5_DP, 2._DP, 3._DP, 5.5_DP, 10.5_DP, 15._DP, &
    20.5_DP, 21.5_DP, 30.5_DP, 50.5_DP, 100.5_DP, 170.5_DP ]
  REAL(DP), PARAMETER :: pox(26) = [ -30.5_DP, -21._DP, -20.5_DP, -20._DP, &
    -10._DP, -7.5_DP, -5._DP, -3._DP, -1._DP, -0.5_DP, -0.25_DP, -0.001_DP, &
    0._DP, 0.001_DP, 0.25_DP, 0.5_DP, 1._DP, 2._DP, 3._DP, 5._DP, 7.5_DP, &
    10._DP, 20._DP, 20.5_DP, 21._DP, 30.5_DP ]

  REAL(DP), PARAMETER :: p1a(23) = [ -100.5_DP, -50.5_DP, -30.5_DP, -20.5_DP, &
    -10.5_DP, -5.5_DP, -2.5_DP, -1.5_DP, -0.75_DP, -0.5_DP, -0.25_DP, 0.25_DP, &
    0.5_DP, 1._DP, 2._DP, 5._DP, 10._DP, 11._DP, 15._DP, 20._DP, 50._DP, &
    100._DP, 1000._DP ]
  REAL(DP), PARAMETER :: p1x(25) = [ 0._DP, 1.0E-8_DP, -1.0E-8_DP, 1.0E-5_DP, &
    -1.0E-5_DP, 0.001_DP, -0.001_DP, 0.01_DP, -0.01_DP, 0.02_DP, -0.02_DP, &
    0.05_DP, -0.05_DP, 0.1_DP, -0.1_DP, 0.25_DP, -0.25_DP, 0.5_DP, -0.5_DP, &
    1._DP, -1._DP, 2._DP, -2._DP, 5._DP, -5._DP ]

  REAL(DP), PARAMETER :: btp(23) = [ 0.001_DP, 0.01_DP, 0.1_DP, 0.25_DP, &
    0.5_DP, 1._DP, 1.5_DP, 2._DP, 3._DP, 5._DP, 9._DP, 9.5_DP, 9.9_DP, &
    10._DP, 10.5_DP, 15._DP, 20._DP, 50._DP, 100._DP, 500._DP, 1000._DP, &
    1.0E4_DP, 1.0E5_DP ]

  REAL(DP), PARAMETER :: bix(19) = [ 0._DP, 1.0E-300_DP, 1.0E-20_DP, &
    1.0E-8_DP, 1.0E-4_DP, 0.001_DP, 0.01_DP, 0.1_DP, 0.19_DP, 0.2_DP, &
    0.3_DP, 0.5_DP, 0.7_DP, 0.79_DP, 0.8_DP, 0.9_DP, 0.99_DP, 0.999_DP, 1._DP ]
  REAL(DP), PARAMETER :: biq(14) = [ 0.001_DP, 0.1_DP, 0.5_DP, 1._DP, 1.5_DP, &
    2._DP, 3._DP, 5._DP, 9.5_DP, 10._DP, 20._DP, 50._DP, 200._DP, 1000._DP ]

  REAL(DP), PARAMETER :: gia(13) = [ 0.001_DP, 0.01_DP, 0.1_DP, 0.5_DP, 1._DP, &
    1.5_DP, 2._DP, 3._DP, 5._DP, 10._DP, 20._DP, 50._DP, 100._DP ]
  REAL(DP), PARAMETER :: gix(17) = [ 0._DP, 1.0E-10_DP, 1.0E-5_DP, 0.001_DP, &
    0.01_DP, 0.1_DP, 0.5_DP, 0.9_DP, 1._DP, 1.5_DP, 2._DP, 5._DP, 10._DP, &
    20._DP, 50._DP, 100._DP, 200._DP ]

  REAL(DP), PARAMETER :: gca(20) = [ -20.5_DP, -10._DP, -5.001_DP, -5._DP, &
    -4.999_DP, -2.5_DP, -1.001_DP, -1._DP, -0.999_DP, -0.5_DP, -0.001_DP, &
    0._DP, 0.001_DP, 0.5_DP, 1._DP, 2.5_DP, 5._DP, 10._DP, 20.5_DP, 50._DP ]
  REAL(DP), PARAMETER :: gcx(15) = [ 1.0E-8_DP, 1.0E-4_DP, 0.001_DP, 0.01_DP, &
    0.1_DP, 0.5_DP, 0.9_DP, 1._DP, 1.5_DP, 2._DP, 5._DP, 10._DP, 20._DP, &
    50._DP, 100._DP ]

  REAL(DP), PARAMETER :: gta(18) = [ -20.5_DP, -10._DP, -5._DP, -4.999_DP, &
    -2.5_DP, -1._DP, -0.999_DP, -0.5_DP, -0.001_DP, 0._DP, 0.001_DP, 0.5_DP, &
    1._DP, 2.5_DP, 5._DP, 10._DP, 20.5_DP, 50._DP ]
  REAL(DP), PARAMETER :: gtx(15) = [ 0._DP, 1.0E-8_DP, 0.001_DP, 0.01_DP, &
    0.1_DP, 0.5_DP, 0.9_DP, 1._DP, 1.5_DP, 2._DP, 5._DP, 10._DP, 20._DP, &
    50._DP, 100._DP ]

  REAL(DP), PARAMETER :: mica(13) = [ 0._DP, -0.1_DP, -0.4_DP, -0.5_DP, &
    -0.6_DP, -1._DP, -1.5_DP, -2._DP, -2.3_DP, -3._DP, -5._DP, -10._DP, &
    -20.5_DP ]
  REAL(DP), PARAMETER :: micx(8) = [ 1.0E-10_DP, 1.0E-5_DP, 0.001_DP, 0.01_DP, &
    0.1_DP, 0.5_DP, 0.9_DP, 1._DP ]

  ! mitg is log|Gamma(a+1)| and mits the sign of Gamma(a+1), handed to d9gmit
  ! as literals so that both sides get identical inputs rather than each
  ! side's own log gamma.
  REAL(DP), PARAMETER :: mita(15) = [ -10.5_DP, -5.5_DP, -2.5_DP, -1.5_DP, &
    -0.999_DP, -0.5_DP, -0.25_DP, 0._DP, 0.25_DP, 0.5_DP, 1._DP, 2.5_DP, &
    5._DP, 10._DP, 20.5_DP ]
  REAL(DP), PARAMETER :: mitg(15) = [ -12.795895333554364_DP, &
    -2.8130840817693166_DP, 0.86004701537648121_DP, 1.2655121234846449_DP, &
    6.9071788853838525_DP, 0.57236494292470042_DP, 0.20328095143129499_DP, &
    0._DP, -0.098271836421812697_DP, -0.12078223763524543_DP, 0._DP, &
    1.2009736023470738_DP, 4.7874917427820467_DP, 15.104412573075514_DP, &
    43.851925860675159_DP ]
  REAL(DP), PARAMETER :: mits(15) = [ 1._DP, -1._DP, 1._DP, -1._DP, 1._DP, &
    1._DP, 1._DP, 1._DP, 1._DP, 1._DP, 1._DP, 1._DP, 1._DP, 1._DP, 1._DP ]
  REAL(DP), PARAMETER :: mitx(9) = [ 1.0E-8_DP, 1.0E-4_DP, 0.001_DP, 0.01_DP, &
    0.1_DP, 0.25_DP, 0.5_DP, 0.9_DP, 1._DP ]

  REAL(DP), PARAMETER :: gica(12) = [ -20.5_DP, -10._DP, -5._DP, -1._DP, &
    0._DP, 0.5_DP, 1._DP, 2.5_DP, 5._DP, 10._DP, 20._DP, 50._DP ]
  REAL(DP), PARAMETER :: gicx(9) = [ 1._DP, 1.5_DP, 2._DP, 5._DP, 10._DP, &
    20._DP, 50._DP, 100._DP, 500._DP ]

  REAL(DP), PARAMETER :: gita(8) = [ 0.5_DP, 1._DP, 2._DP, 5._DP, 10._DP, &
    20._DP, 50._DP, 100._DP ]
  REAL(DP), PARAMETER :: gitg(8) = [ -0.12078223763524543_DP, 0._DP, &
    0.69314718055994495_DP, 4.7874917427820467_DP, 15.104412573075514_DP, &
    42.335616460753485_DP, 148.47776695177305_DP, 363.73937555556347_DP ]
  REAL(DP), PARAMETER :: gitf(10) = [ 0.001_DP, 0.01_DP, 0.05_DP, 0.1_DP, &
    0.25_DP, 0.5_DP, 0.75_DP, 0.9_DP, 0.99_DP, 1._DP ]

  REAL(DP), PARAMETER :: pfx(22) = [ 1.0E-20_DP, 1.0E-18_DP, 1.0E-17_DP, &
    1.0E-16_DP, 1.0E-10_DP, 1.0E-5_DP, 0.001_DP, 0.01_DP, 0.1_DP, 0.5_DP, &
    1._DP, 2._DP, 5._DP, 7._DP, 10._DP, 15._DP, 20._DP, 50._DP, 100._DP, &
    1000._DP, 1.0E5_DP, 1.0E10_DP ]
  INTEGER, PARAMETER :: pfn(9) = [ 0, 1, 2, 5, 10, 20, 50, 100, 200 ]
  INTEGER, PARAMETER :: pfm(4) = [ 1, 2, 3, 5 ]

  REAL(DP), PARAMETER :: chua(8) = [ 0.5_DP, 1._DP, 1.5_DP, 2._DP, 3._DP, &
    -0.5_DP, -1.5_DP, 0.25_DP ]
  REAL(DP), PARAMETER :: chub(11) = [ 0.3_DP, 0.5_DP, 0.7_DP, 1._DP, 1.2_DP, &
    1.5_DP, 2._DP, 2.5_DP, 3._DP, 4.4_DP, -0.5_DP ]
  REAL(DP), PARAMETER :: chuz(9) = [ 1._DP, 2._DP, 5._DP, 10._DP, 20._DP, &
    50._DP, 100._DP, 1000._DP, 1.0E4_DP ]
  REAL(DP), PARAMETER :: chux(12) = [ 0.01_DP, 0.1_DP, 0.5_DP, 0.9_DP, 1._DP, &
    1.1_DP, 2._DP, 3._DP, 5._DP, 10._DP, 50._DP, 100._DP ]

  REAL(DP), PARAMETER :: jjl2(27) = [ 0._DP, 1._DP, 1._DP, 1._DP, 2._DP, &
    2._DP, 3._DP, 3._DP, 5._DP, 10._DP, 10._DP, 10._DP, 20._DP, 50._DP, &
    100._DP, 200._DP, 400._DP, 0.5_DP, 0.5_DP, 1.5_DP, 1.5_DP, 10.5_DP, &
    10.5_DP, 20.5_DP, 50.5_DP, 100.5_DP, 200.5_DP ]
  REAL(DP), PARAMETER :: jjl3(27) = [ 0._DP, 0._DP, 1._DP, 1._DP, 2._DP, &
    3._DP, 3._DP, 4._DP, 5._DP, 10._DP, 10._DP, 10._DP, 30._DP, 50._DP, &
    100._DP, 150._DP, 400._DP, 0.5_DP, 1.5_DP, 2.5_DP, 1.5_DP, 20.5_DP, &
    10.5_DP, 20.5_DP, 50.5_DP, 100.5_DP, 200.5_DP ]
  REAL(DP), PARAMETER :: jjm2(27) = [ 0._DP, 0._DP, 0._DP, 1._DP, 2._DP, &
    1._DP, 0._DP, 1._DP, 2._DP, 0._DP, 5._DP, 3._DP, 5._DP, 10._DP, &
    0._DP, 20._DP, 1._DP, 0.5_DP, 0.5_DP, 0.5_DP, 1.5_DP, 0.5_DP, &
    10.5_DP, 0.5_DP, 10.5_DP, 0.5_DP, 0.5_DP ]
  REAL(DP), PARAMETER :: jjm3(27) = [ 0._DP, 0._DP, 0._DP, -1._DP, -2._DP, &
    -1._DP, 0._DP, 1._DP, -2._DP, 0._DP, -5._DP, 2._DP, -5._DP, -10._DP, &
    0._DP, -20._DP, -1._DP, -0.5_DP, -0.5_DP, -0.5_DP, -1.5_DP, -0.5_DP, &
    -10.5_DP, 0.5_DP, -10.5_DP, -0.5_DP, -0.5_DP ]

  REAL(DP), PARAMETER :: jml1(20) = [ 1._DP, 1._DP, 1._DP, 2._DP, 2._DP, &
    5._DP, 10._DP, 10._DP, 20._DP, 50._DP, 100._DP, 200._DP, 300._DP, &
    0.5_DP, 1.5_DP, 1.5_DP, 10.5_DP, 10.5_DP, 50.5_DP, 100.5_DP ]
  REAL(DP), PARAMETER :: jml2(20) = [ 1._DP, 1._DP, 0._DP, 3._DP, 2._DP, &
    5._DP, 10._DP, 10._DP, 30._DP, 50._DP, 100._DP, 150._DP, 300._DP, &
    0.5_DP, 2.5_DP, 1.5_DP, 10.5_DP, 20.5_DP, 50.5_DP, 100.5_DP ]
  REAL(DP), PARAMETER :: jml3(20) = [ 1._DP, 1._DP, 1._DP, 4._DP, 2._DP, &
    5._DP, 10._DP, 10._DP, 40._DP, 50._DP, 100._DP, 100._DP, 300._DP, &
    1._DP, 2._DP, 1._DP, 20._DP, 10._DP, 100._DP, 200._DP ]
  REAL(DP), PARAMETER :: jmm1(20) = [ 0._DP, 1._DP, 0._DP, 1._DP, 0._DP, &
    2._DP, 0._DP, 5._DP, 5._DP, 10._DP, 0._DP, 20._DP, 0._DP, &
    0.5_DP, 0.5_DP, 1.5_DP, 0.5_DP, 10.5_DP, 10.5_DP, 0.5_DP ]

  REAL(DP), PARAMETER :: sjl2(17) = [ 0._DP, 1._DP, 1._DP, 1._DP, 2._DP, &
    3._DP, 5._DP, 10._DP, 20._DP, 50._DP, 100._DP, 200._DP, 0.5_DP, 1.5_DP, &
    1.5_DP, 10.5_DP, 50.5_DP ]
  REAL(DP), PARAMETER :: sjl3(17) = [ 0._DP, 1._DP, 1._DP, 1._DP, 2._DP, &
    4._DP, 5._DP, 10._DP, 30._DP, 50._DP, 100._DP, 200._DP, 0.5_DP, 1.5_DP, &
    2.5_DP, 10.5_DP, 50.5_DP ]
  REAL(DP), PARAMETER :: sjl4(17) = [ 0._DP, 1._DP, 2._DP, 0._DP, 2._DP, &
    5._DP, 5._DP, 10._DP, 40._DP, 50._DP, 100._DP, 200._DP, 1._DP, 1._DP, &
    2._DP, 10._DP, 50._DP ]
  REAL(DP), PARAMETER :: sjl5(17) = [ 0._DP, 1._DP, 1._DP, 1._DP, 2._DP, &
    6._DP, 5._DP, 10._DP, 50._DP, 50._DP, 100._DP, 200._DP, 0.5_DP, 1.5_DP, &
    2.5_DP, 10.5_DP, 50.5_DP ]
  REAL(DP), PARAMETER :: sjl6(17) = [ 0._DP, 1._DP, 1._DP, 1._DP, 2._DP, &
    7._DP, 5._DP, 10._DP, 60._DP, 50._DP, 100._DP, 200._DP, 0.5_DP, 1.5_DP, &
    1.5_DP, 10.5_DP, 50.5_DP ]

CONTAINS

  ! dpoch stops when a+x is a non-positive integer and a is not.  It also
  ! reaches d9lgmc below ten, which stops in turn, when both a and a+x are
  ! non-positive integers, MIN(a+x,a) is under -20 and either of them is
  ! above -9, because the branch at dpoch.inc:84 asks for d9lgmc(1-a) and
  ! d9lgmc(1-a-x) without checking.  The port raises in the same place, so
  ! the pair is left out of the sweep on both sides either way.
  LOGICAL PURE FUNCTION pochbad(a, x)
    REAL(DP), INTENT(IN) :: a, x
    REAL(DP) :: ax
    LOGICAL :: axint, aint_
    ax = a + x
    axint = ax <= 0._DP .AND. AINT(ax) == ax
    aint_ = a <= 0._DP .AND. AINT(a) == a
    pochbad = axint .AND. .NOT. aint_
    IF (axint .AND. aint_ .AND. x /= 0._DP) THEN
      IF (MIN(ax, a) < -20._DP .AND. (a > -9._DP .OR. ax > -9._DP)) &
        pochbad = .TRUE.
    END IF
  END FUNCTION pochbad

  ! dpsifn reports Ierr = 2 where the port raises instead, so the overflow
  ! corner is left out and the underflow one, which only sets Nz, is kept.
  LOGICAL PURE FUNCTION psifnbad(x, n, m)
    REAL(DP), INTENT(IN) :: x
    INTEGER, INTENT(IN) :: n, m
    REAL(DP) :: t, elim
    elim = 2.302_DP*(MIN(-min_exp_dp, max_exp_dp)*log10_radix_dp - 3._DP)
    t = (n + m)*LOG(x)
    psifnbad = ABS(t) > elim .AND. t <= 0._DP
  END FUNCTION psifnbad

  ! d9chu's own guard at d9chu.inc:66 is live, because i there is the loop
  ! variable, and the fraction never settles in 300 terms when 1+a-b is a
  ! negative integer, at any z.  Both sides agree about that and both give up,
  ! the Fortran by stopping and the port by raising, so the pair is left out.
  LOGICAL PURE FUNCTION chuzbad(a, b)
    REAL(DP), INTENT(IN) :: a, b
    REAL(DP) :: bp
    bp = 1._DP + a - b
    chuzbad = bp < 0._DP .AND. AINT(bp) == bp
  END FUNCTION chuzbad

  ! dchu stops on 1+a-b near zero, and reaches dpoch and dpsi at arguments
  ! that stop in turn, all of them only on the small x path.
  LOGICAL PURE FUNCTION chubad(a, b, x)
    REAL(DP), INTENT(IN) :: a, b, x
    REAL(DP) :: aintb, beps, xi, s
    INTEGER :: n, istrt
    chubad = .FALSE.
    IF (MAX(ABS(a), 1._DP)*MAX(ABS(1._DP+a-b), 1._DP) < 0.99_DP*ABS(x)) THEN
      chubad = chuzbad(a, b)
      RETURN
    END IF
    IF (ABS(1._DP+a-b) < SQRT(eps_2_dp)) THEN
      chubad = .TRUE.
      RETURN
    END IF
    IF (b >= 0._DP) THEN
      aintb = AINT(b+0.5_DP)
    ELSE
      aintb = AINT(b-0.5_DP)
    END IF
    beps = b - aintb
    n = INT(aintb)
    istrt = 0
    IF (n < 1) istrt = 1 - n
    xi = istrt
    IF (n < 1) THEN
      s = 1._DP - b
      IF ((s <= 0._DP .AND. AINT(s) == s) .AND. &
          .NOT. (1._DP+a-b <= 0._DP .AND. AINT(1._DP+a-b) == 1._DP+a-b)) THEN
        chubad = .TRUE.
        RETURN
      END IF
    END IF
    s = a + xi - beps
    IF ((s <= 0._DP .AND. AINT(s) == s) .AND. &
        .NOT. (a <= 0._DP .AND. AINT(a) == a)) THEN
      chubad = .TRUE.
      RETURN
    END IF
    s = a + xi
    IF (beps == 0._DP .AND. s <= 0._DP .AND. AINT(s) == s) chubad = .TRUE.
  END FUNCTION chubad

END MODULE cases
