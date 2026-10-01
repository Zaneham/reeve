! Reeve differential, interpolation. The machine constant and error handler
! seam these routines need. DPLINT, DPOLCF and DPOLVL call nothing but XERMSG,
! so D1MACH and I1MACH are not referenced and are not supplied here.

MODULE xerstate
  IMPLICIT NONE
  INTEGER :: xercnt = 0
END MODULE xerstate

! XERMSG at level 1 is recoverable and returns to the caller, which is what the
! original DPLINT assumes when it falls through to a RETURN after the call. The
! counter is the only trace it leaves.
SUBROUTINE XERMSG(Librar, Subrou, Messg, Nerr, Level)
  USE xerstate, ONLY : xercnt
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN) :: Librar, Subrou, Messg
  INTEGER, INTENT(IN) :: Nerr, Level
  xercnt = xercnt + 1
END SUBROUTINE XERMSG
