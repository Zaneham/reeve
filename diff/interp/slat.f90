! Reeve differential, interpolation. The modernised SLATEC interpolation trio,
! included from C:\dev\numerical\slatec-modern\src\modern\interpolation.

MODULE slatec_interp
  USE service, ONLY : DP
  IMPLICIT NONE
CONTAINS
  include "dplint.inc"
  include "dpolcf.inc"
  include "dpolvl.inc"
END MODULE slatec_interp
