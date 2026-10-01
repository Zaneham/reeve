MODULE service
  IMPLICIT NONE
  INTEGER, PARAMETER :: SP = SELECTED_REAL_KIND(6,37)
  INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(15,307)
  INTEGER, PARAMETER :: radix_fp = RADIX(1._SP), &
    min_exp_sp = MINEXPONENT(1._SP), max_exp_sp = MAXEXPONENT(1._SP), &
    min_exp_dp = MINEXPONENT(1._DP), max_exp_dp = MAXEXPONENT(1._DP)
  REAL(SP), PARAMETER :: tiny_sp = TINY(1._SP), huge_sp = HUGE(1._SP), &
    eps_2_sp = EPSILON(1._SP)/RADIX(1._SP), eps_sp = EPSILON(1._SP), &
    log10_radix_sp = LOG10( REAL( RADIX(1._SP), SP ) )
  REAL(DP), PARAMETER :: tiny_dp = TINY(1._DP), huge_dp = HUGE(1._DP), &
    eps_2_dp = EPSILON(1._DP)/RADIX(1._DP), eps_dp = EPSILON(1._DP), &
    log10_radix_dp = LOG10( REAL( RADIX(1._DP), DP ) )
END MODULE service

MODULE sf
  USE service
  IMPLICIT NONE
CONTAINS
  include "dcsevl.inc"
  include "csevl.inc"
  include "d9upak.inc"
  include "d9pak.inc"
  include "dcbrt.inc"
  include "dexprl.inc"
  include "dlnrel.inc"
  include "d9ln2r.inc"
  include "d9atn1.inc"
  include "dsindg.inc"
  include "dcosdg.inc"
  include "dspenc.inc"
  include "de1.inc"
  include "dei.inc"
  include "dli.inc"
  include "ddaws.inc"
  include "daws.inc"
  include "dpsixn.inc"
  include "dexint.inc"
END MODULE sf
