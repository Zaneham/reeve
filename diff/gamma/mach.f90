! D1MACH, I1MACH and XERMSG for the differential driver.  The values are the
! IEEE ones src/modern/service/service.f90 defines, taken from that module
! rather than typed again here.  The original-source dpsifn.f is the only
! caller, and it asks for D1MACH(4), D1MACH(5), I1MACH(14), I1MACH(15) and
! I1MACH(16).
DOUBLE PRECISION FUNCTION D1MACH(I)
  USE service, ONLY: tiny_dp, huge_dp, eps_2_dp, eps_dp, log10_radix_dp
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: I
  SELECT CASE (I)
  CASE (1)
    D1MACH = tiny_dp
  CASE (2)
    D1MACH = huge_dp
  CASE (3)
    D1MACH = eps_2_dp
  CASE (4)
    D1MACH = eps_dp
  CASE (5)
    D1MACH = log10_radix_dp
  CASE DEFAULT
    ERROR STOP 'D1MACH : I OUT OF BOUNDS'
  END SELECT
END FUNCTION D1MACH

INTEGER FUNCTION I1MACH(I)
  USE service, ONLY: radix_int, digits_int, huge_int, radix_fp, digits_sp, &
    min_exp_sp, max_exp_sp, digits_dp, min_exp_dp, max_exp_dp
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: I
  SELECT CASE (I)
  CASE (1)
    I1MACH = 5
  CASE (2)
    I1MACH = 6
  CASE (3)
    I1MACH = 7
  CASE (4)
    I1MACH = 0
  CASE (5)
    I1MACH = BIT_SIZE(1)
  CASE (6)
    I1MACH = BIT_SIZE(1)/8
  CASE (7)
    I1MACH = radix_int
  CASE (8)
    I1MACH = digits_int
  CASE (9)
    I1MACH = huge_int
  CASE (10)
    I1MACH = radix_fp
  CASE (11)
    I1MACH = digits_sp
  CASE (12)
    I1MACH = min_exp_sp
  CASE (13)
    I1MACH = max_exp_sp
  CASE (14)
    I1MACH = digits_dp
  CASE (15)
    I1MACH = min_exp_dp
  CASE (16)
    I1MACH = max_exp_dp
  CASE DEFAULT
    ERROR STOP 'I1MACH : I OUT OF BOUNDS'
  END SELECT
END FUNCTION I1MACH

SUBROUTINE XERMSG(Librar, Subrou, Messg, Nerr, Level)
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN) :: Librar, Subrou, Messg
  INTEGER, INTENT(IN) :: Nerr, Level
END SUBROUTINE XERMSG
