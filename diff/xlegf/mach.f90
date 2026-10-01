! Reeve differential, machine constants for the extended-range Legendre
! package.  The DATA blocks in the SLATEC originals of D1MACH and I1MACH are
! all commented out, so these are supplied here from the IEEE values that
! src/modern/service/service.f90 defines.
!
! I1MACH(8) is the one value that is not simply the IEEE double format: it is
! the number of magnitude bits in an integer, and the port's target is an
! OCaml int, which has 62.  The Fortran side is built with
! -fdefault-integer-8 so that a 62-bit index arithmetic has room, the same
! room OCaml's 63-bit int gives it.

FUNCTION D1MACH(I) RESULT(R)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: I
  DOUBLE PRECISION :: R
  SELECT CASE(I)
  CASE(1); R = TINY(1.0D0)
  CASE(2); R = HUGE(1.0D0)
  CASE(3); R = EPSILON(1.0D0)/RADIX(1.0D0)
  CASE(4); R = EPSILON(1.0D0)
  CASE(5); R = LOG10(REAL(RADIX(1.0D0), KIND(1.0D0)))
  CASE DEFAULT; R = 0.0D0
  END SELECT
END FUNCTION D1MACH

FUNCTION I1MACH(I) RESULT(K)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: I
  INTEGER :: K
  SELECT CASE(I)
  CASE(7);  K = 2
  CASE(8);  K = 62
  CASE(10); K = 2
  CASE(11); K = DIGITS(1.0E0)
  CASE(12); K = MINEXPONENT(1.0E0)
  CASE(13); K = MAXEXPONENT(1.0E0)
  CASE(14); K = DIGITS(1.0D0)
  CASE(15); K = MINEXPONENT(1.0D0)
  CASE(16); K = MAXEXPONENT(1.0D0)
  CASE DEFAULT; K = 0
  END SELECT
END FUNCTION I1MACH
