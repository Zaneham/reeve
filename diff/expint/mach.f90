FUNCTION D1MACH(I) RESULT(R)
  USE service, ONLY : DP, tiny_dp, huge_dp, eps_2_dp, eps_dp, log10_radix_dp
  INTEGER, INTENT(IN) :: I
  REAL(DP) :: R
  SELECT CASE(I)
  CASE(1); R = tiny_dp
  CASE(2); R = huge_dp
  CASE(3); R = eps_2_dp
  CASE(4); R = eps_dp
  CASE(5); R = log10_radix_dp
  CASE DEFAULT; R = 0._DP
  END SELECT
END FUNCTION D1MACH

FUNCTION I1MACH(I) RESULT(K)
  USE service, ONLY : DP, SP, min_exp_dp, max_exp_dp, min_exp_sp, max_exp_sp
  INTEGER, INTENT(IN) :: I
  INTEGER :: K
  SELECT CASE(I)
  CASE(10); K = 2
  CASE(11); K = DIGITS(1._SP)
  CASE(12); K = min_exp_sp
  CASE(13); K = max_exp_sp
  CASE(14); K = DIGITS(1._DP)
  CASE(15); K = min_exp_dp
  CASE(16); K = max_exp_dp
  CASE DEFAULT; K = 0
  END SELECT
END FUNCTION I1MACH
