! The reference Fortran for the gamma differential, built straight out of
! src/modern/special_functions rather than vendored here.  DPSIFN is missing on
! purpose and comes from src/original/src/dpsifn.f instead, because the
! modernised dpsifn.inc declares X, Ans and every local as REAL(SP) while
! filling its table with _DP literals, which leaves the double precision
! routine working in single precision.
MODULE gmod
  USE service
  IMPLICIT NONE
CONTAINS
  INCLUDE "dcsevl.inc"
  INCLUDE "dlnrel.inc"
  INCLUDE "dcot.inc"
  INCLUDE "dexprl.inc"
  INCLUDE "dpsi.inc"
  INCLUDE "dpsixn.inc"
  INCLUDE "d9lgmc.inc"
  INCLUDE "dgamlm.inc"
  INCLUDE "dgamr.inc"
  INCLUDE "dlgams.inc"
  INCLUDE "dfac.inc"
  INCLUDE "dbinom.inc"
  INCLUDE "dpoch.inc"
  INCLUDE "dpoch1.inc"
  INCLUDE "dlbeta.inc"
  INCLUDE "dbeta.inc"
  INCLUDE "dbetai.inc"
  INCLUDE "d9gmic.inc"
  INCLUDE "d9gmit.inc"
  INCLUDE "d9lgic.inc"
  INCLUDE "d9lgit.inc"
  INCLUDE "dgami.inc"
  INCLUDE "dgamic.inc"
  INCLUDE "dgamit.inc"
  INCLUDE "d9chu.inc"
  INCLUDE "dchu.inc"
  INCLUDE "drc3jj.inc"
  INCLUDE "drc3jm.inc"
  INCLUDE "drc6j.inc"
END MODULE gmod
