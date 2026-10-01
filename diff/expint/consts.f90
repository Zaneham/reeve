PROGRAM consts
  USE service
  IMPLICIT NONE
  REAL(DP), PARAMETER :: alneps = LOG(eps_2_dp), xn = 3.72_DP - 0.3_DP*alneps, &
    xln = LOG((xn+1._DP)/1.36_DP)
  INTEGER, PARAMETER :: nterms = INT( xn - (xn*xln+alneps)/(xln+1.36_DP) + 1.5_DP )
  INTEGER, PARAMETER :: niter = INT( 1.443*LOG(-.106*LOG(0.1_SP*eps_2_dp)) )
  REAL(DP), PARAMETER :: a1n2b_inc = log10_radix_dp/LOG(2._DP)
  WRITE(*,'(A,I0)') 'dexprl nterms = ', nterms
  WRITE(*,'(A,I0)') 'dcbrt niter   = ', niter
  WRITE(*,'(A,I0,A,I0)') 'inc d9pak nmin/nmax = ', INT(a1n2b_inc*min_exp_dp), ' / ', INT(a1n2b_inc*max_exp_dp)
  WRITE(*,'(A,I0,A,I0)') 'orig d9pak nmin/nmax = ', min_exp_dp, ' / ', max_exp_dp
END PROGRAM consts
