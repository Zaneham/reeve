! Fortran side of the gamma differential.  It sweeps the inputs cases.f90
! holds and writes one line per value, a tag, a case index and the raw IEEE
! bit pattern.  odrv.ml runs the same sweep through lib/gamma.ml and the two
! outputs are compared with cmp.
!
! Every routine comes from src/modern/special_functions, except DPSIFN, which
! comes from src/original/src/dpsifn.f because the modernised include file
! declares X, Ans and every local as REAL(SP) and so is not the double
! precision routine at all.  DGAMLM follows the include file, where the
! original's 0.01 safety nudge has been swallowed into a comment and no longer
! happens, which is what the port follows too.
MODULE emit
  USE service, ONLY: DP
  IMPLICIT NONE
CONTAINS
  SUBROUTINE em(tag, idx, v)
    CHARACTER(*), INTENT(IN) :: tag
    INTEGER, INTENT(IN) :: idx
    REAL(DP), INTENT(IN) :: v
    CHARACTER(12) :: t
    INTEGER(8) :: b
    t = tag
    b = TRANSFER(v, b)
    WRITE (*, '(A,1X,I7,1X,Z16.16)') t, idx, b
  END SUBROUTINE em

  SUBROUTINE emi(tag, idx, k)
    CHARACTER(*), INTENT(IN) :: tag
    INTEGER, INTENT(IN) :: idx, k
    CHARACTER(12) :: t
    t = tag
    WRITE (*, '(A,1X,I7,1X,I16)') t, idx, k
  END SUBROUTINE emi
END MODULE emit

PROGRAM fdrv
  USE service
  USE cases
  USE gmod
  USE emit
  IMPLICIT NONE

  INTEGER :: i, j, k, jj, n, m, nc, nf, kode, nz, ierr
  REAL(DP) :: a, b, x, p, q, xmin, xmax, dl, sg, lmn, lmx
  REAL(DP) :: ans(8), wig(4100)

  ! ---- The derived constants, which gfortran folds at compile time while
  ! ---- OCaml works them out at run time from libm.
  CALL em('const', 1, LOG(tiny_dp))
  CALL em('const', 2, LOG(huge_dp))
  CALL em('const', 3, log10_radix_dp)
  CALL em('const', 4, eps_dp)
  CALL em('const', 5, eps_2_dp)
  CALL em('const', 6, tiny_dp)
  CALL em('const', 7, huge_dp)
  CALL em('const', 8, 1._DP/SQRT(eps_2_dp))
  CALL em('const', 9, EXP(MIN(LOG(huge_dp/12._DP), -LOG(12._DP*tiny_dp))))
  CALL em('const', 10, 1._DP/SQRT(24._DP*tiny_dp))
  CALL em('const', 11, LOG(eps_2_dp))
  CALL em('const', 12, LOG(huge_dp) - 0.0001_DP)
  CALL em('const', 13, 0.9_DP/eps_2_dp)
  CALL em('const', 14, -LOG(eps_2_dp))
  CALL em('const', 15, 0.5_DP*eps_2_dp)
  CALL em('const', 16, SQRT(huge_dp/20._DP))
  CALL em('const', 17, MAX(eps_dp, 1.E-18_DP))
  CALL em('const', 18, 2.302_DP*(MIN(-min_exp_dp, max_exp_dp)*log10_radix_dp - 3._DP))
  CALL em('const', 19, MAX(eps_dp*0.5_DP, 0.5E-18_DP))
  CALL em('const', 20, MIN(log10_radix_dp*digits_dp, 18.06_DP))

  ! ---- dgamlm
  CALL DGAMLM(xmin, xmax)
  CALL em('dgamlm:min', 1, xmin)
  CALL em('dgamlm:max', 1, xmax)

  ! ---- dgamr and dlgams over x = i/4 from -30 to 170
  nc = 0
  DO i = -120, 680
    x = i*0.25_DP
    IF (x <= 0._DP .AND. AINT(x) == x) CYCLE
    nc = nc + 1
    CALL em('dgamr:x', nc, x)
    CALL em('dgamr', nc, DGAMR(x))
    CALL DLGAMS(x, dl, sg)
    CALL em('dlgams:l', nc, dl)
    CALL em('dlgams:s', nc, sg)
  END DO

  ! ---- d9lgmc
  nc = 0
  DO i = 0, 400
    x = 10._DP + i*0.5_DP
    nc = nc + 1
    CALL em('d9lgmc:x', nc, x)
    CALL em('d9lgmc', nc, D9LGMC(x))
  END DO
  DO i = 1, SIZE(lgmcx)
    x = lgmcx(i)
    nc = nc + 1
    CALL em('d9lgmc:x', nc, x)
    CALL em('d9lgmc', nc, D9LGMC(x))
  END DO

  ! ---- dfac
  DO n = 0, 170
    CALL emi('dfac:n', n + 1, n)
    CALL em('dfac', n + 1, DFAC(n))
  END DO

  ! ---- dbinom
  nc = 0
  DO n = 0, 40
    DO m = 0, n
      nc = nc + 1
      CALL emi('dbinom:n', nc, n)
      CALL emi('dbinom:m', nc, m)
      CALL em('dbinom', nc, DBINOM(n, m))
    END DO
  END DO
  DO i = 1, SIZE(binn)
    nc = nc + 1
    CALL emi('dbinom:n', nc, binn(i))
    CALL emi('dbinom:m', nc, binm(i))
    CALL em('dbinom', nc, DBINOM(binn(i), binm(i)))
  END DO

  ! ---- dpoch, the a+x non-positive integer branch first
  nc = 0
  DO i = 1, SIZE(pcha)
    a = pcha(i)
    x = pchx(i)
    nc = nc + 1
    CALL em('dpoch:a', nc, a)
    CALL em('dpoch:x', nc, x)
    CALL em('dpoch', nc, DPOCH(a, x))
  END DO
  DO i = 1, SIZE(poa)
    DO j = 1, SIZE(pox)
      a = poa(i)
      x = pox(j)
      IF (pochbad(a, x)) CYCLE
      nc = nc + 1
      CALL em('dpoch:a', nc, a)
      CALL em('dpoch:x', nc, x)
      CALL em('dpoch', nc, DPOCH(a, x))
    END DO
  END DO

  ! ---- dpoch1
  nc = 0
  DO i = 1, SIZE(p1a)
    DO j = 1, SIZE(p1x)
      a = p1a(i)
      x = p1x(j)
      IF (pochbad(a, x)) CYCLE
      nc = nc + 1
      CALL em('dpoch1:a', nc, a)
      CALL em('dpoch1:x', nc, x)
      CALL em('dpoch1', nc, DPOCH1(a, x))
    END DO
  END DO

  ! ---- dlbeta and dbeta
  nc = 0
  DO i = 1, SIZE(btp)
    DO j = 1, SIZE(btp)
      p = btp(i)
      q = btp(j)
      nc = nc + 1
      CALL em('dlbeta:a', nc, p)
      CALL em('dlbeta:b', nc, q)
      CALL em('dlbeta', nc, DLBETA(p, q))
      CALL em('dbeta', nc, DBETA(p, q))
    END DO
  END DO

  ! ---- dbetai
  nc = 0
  DO i = 1, SIZE(bix)
    DO j = 1, SIZE(biq)
      DO k = 1, SIZE(biq)
        x = bix(i)
        p = biq(j)
        q = biq(k)
        nc = nc + 1
        CALL em('dbetai:x', nc, x)
        CALL em('dbetai:p', nc, p)
        CALL em('dbetai:q', nc, q)
        CALL em('dbetai', nc, DBETAI(x, p, q))
      END DO
    END DO
  END DO

  ! ---- dgami
  nc = 0
  DO i = 1, SIZE(gia)
    DO j = 1, SIZE(gix)
      a = gia(i)
      x = gix(j)
      nc = nc + 1
      CALL em('dgami:a', nc, a)
      CALL em('dgami:x', nc, x)
      CALL em('dgami', nc, DGAMI(a, x))
    END DO
  END DO

  ! ---- dgamic
  nc = 0
  DO i = 1, SIZE(gca)
    DO j = 1, SIZE(gcx)
      a = gca(i)
      x = gcx(j)
      nc = nc + 1
      CALL em('dgamic:a', nc, a)
      CALL em('dgamic:x', nc, x)
      CALL em('dgamic', nc, DGAMIC(a, x))
    END DO
  END DO
  DO i = 1, SIZE(gca)
    a = gca(i)
    IF (a <= 0._DP) CYCLE
    nc = nc + 1
    CALL em('dgamic:a', nc, a)
    CALL em('dgamic:x', nc, 0._DP)
    CALL em('dgamic', nc, DGAMIC(a, 0._DP))
  END DO

  ! ---- dgamit
  nc = 0
  DO i = 1, SIZE(gta)
    DO j = 1, SIZE(gtx)
      a = gta(i)
      x = gtx(j)
      nc = nc + 1
      CALL em('dgamit:a', nc, a)
      CALL em('dgamit:x', nc, x)
      CALL em('dgamit', nc, DGAMIT(a, x))
    END DO
  END DO

  ! ---- d9gmic, called directly
  nc = 0
  DO i = 1, SIZE(mica)
    DO j = 1, SIZE(micx)
      a = mica(i)
      x = micx(j)
      nc = nc + 1
      CALL em('d9gmic:a', nc, a)
      CALL em('d9gmic:x', nc, x)
      CALL em('d9gmic:l', nc, LOG(x))
      CALL em('d9gmic', nc, D9GMIC(a, x, LOG(x)))
    END DO
  END DO

  ! ---- d9gmit, called directly with log|Gamma(a+1)| and its sign handed in
  ! ---- as literals, so both sides get the same input
  nc = 0
  DO i = 1, SIZE(mita)
    DO j = 1, SIZE(mitx)
      a = mita(i)
      x = mitx(j)
      nc = nc + 1
      CALL em('d9gmit:a', nc, a)
      CALL em('d9gmit:x', nc, x)
      CALL em('d9gmit:g', nc, mitg(i))
      CALL em('d9gmit:s', nc, mits(i))
      CALL em('d9gmit', nc, D9GMIT(a, x, mitg(i), mits(i)))
    END DO
  END DO

  ! ---- d9lgic, called directly over the a < x region its callers use
  nc = 0
  DO i = 1, SIZE(gica)
    DO j = 1, SIZE(gicx)
      a = gica(i)
      x = gicx(j)
      IF (a >= x) CYCLE
      nc = nc + 1
      CALL em('d9lgic:a', nc, a)
      CALL em('d9lgic:x', nc, x)
      CALL em('d9lgic:l', nc, LOG(x))
      CALL em('d9lgic', nc, D9LGIC(a, x, LOG(x)))
    END DO
  END DO

  ! ---- d9lgit, called directly, needs 0 < x <= a
  nc = 0
  DO i = 1, SIZE(gita)
    DO j = 1, SIZE(gitf)
      a = gita(i)
      x = gitf(j)*a
      IF (x <= 0._DP .OR. x > a) CYCLE
      nc = nc + 1
      CALL em('d9lgit:a', nc, a)
      CALL em('d9lgit:x', nc, x)
      CALL em('d9lgit:g', nc, gitg(i))
      CALL em('d9lgit', nc, D9LGIT(a, x, gitg(i)))
    END DO
  END DO

  ! ---- dpsixn
  DO n = 1, 200
    CALL emi('dpsixn:n', n, n)
    CALL em('dpsixn', n, DPSIXN(n))
  END DO

  ! ---- dpsifn, from src/original/src/dpsifn.f
  nc = 0
  DO i = 1, SIZE(pfx)
    DO j = 1, SIZE(pfn)
      DO kode = 1, 2
        DO k = 1, SIZE(pfm)
          x = pfx(i)
          n = pfn(j)
          m = pfm(k)
          IF (psifnbad(x, n, m)) CYCLE
          nc = nc + 1
          ans = 0._DP
          CALL em('dpsifn:x', nc, x)
          CALL emi('dpsifn:n', nc, n)
          CALL emi('dpsifn:k', nc, kode)
          CALL emi('dpsifn:m', nc, m)
          CALL DPSIFN(x, n, kode, m, ans, nz, ierr)
          CALL emi('dpsifn:nz', nc, nz)
          CALL emi('dpsifn:ie', nc, ierr)
          DO jj = 1, m
            CALL em('dpsifn:a', nc*10 + jj, ans(jj))
          END DO
        END DO
      END DO
    END DO
  END DO

  ! ---- d9chu, called directly
  nc = 0
  DO i = 1, SIZE(chua)
    DO j = 1, SIZE(chub)
      DO k = 1, SIZE(chuz)
        a = chua(i)
        b = chub(j)
        x = chuz(k)
        IF (chuzbad(a, b)) CYCLE
        nc = nc + 1
        CALL em('d9chu:a', nc, a)
        CALL em('d9chu:b', nc, b)
        CALL em('d9chu:z', nc, x)
        CALL em('d9chu', nc, D9CHU(a, b, x))
      END DO
    END DO
  END DO

  ! ---- dchu
  nc = 0
  DO i = 1, SIZE(chua)
    DO j = 1, SIZE(chub)
      DO k = 1, SIZE(chux)
        a = chua(i)
        b = chub(j)
        x = chux(k)
        IF (chubad(a, b, x)) CYCLE
        nc = nc + 1
        CALL em('dchu:a', nc, a)
        CALL em('dchu:b', nc, b)
        CALL em('dchu:x', nc, x)
        CALL em('dchu', nc, DCHU(a, b, x))
      END DO
    END DO
  END DO

  ! ---- drc3jj
  DO i = 1, SIZE(jjl2)
    nf = INT(jjl2(i) + jjl3(i) &
             - MAX(ABS(jjl2(i)-jjl3(i)), ABS(-jjm2(i)-jjm3(i))) + 1._DP + 0.01_DP)
    CALL em('drc3jj:l2', i, jjl2(i))
    CALL em('drc3jj:l3', i, jjl3(i))
    CALL em('drc3jj:m2', i, jjm2(i))
    CALL em('drc3jj:m3', i, jjm3(i))
    CALL emi('drc3jj:nd', i, nf)
    wig = 0._DP
    CALL DRC3JJ(jjl2(i), jjl3(i), jjm2(i), jjm3(i), lmn, lmx, wig, nf, ierr)
    CALL em('drc3jj:mn', i, lmn)
    CALL em('drc3jj:mx', i, lmx)
    DO k = 1, nf
      CALL em('drc3jj', i*1000 + k, wig(k))
    END DO
  END DO

  ! ---- drc3jm
  DO i = 1, SIZE(jml1)
    nf = INT(MIN(jml2(i), jml3(i)-jmm1(i)) &
             - MAX(-jml2(i), -jml3(i)-jmm1(i)) + 1._DP + 0.01_DP)
    CALL em('drc3jm:l1', i, jml1(i))
    CALL em('drc3jm:l2', i, jml2(i))
    CALL em('drc3jm:l3', i, jml3(i))
    CALL em('drc3jm:m1', i, jmm1(i))
    CALL emi('drc3jm:nd', i, nf)
    wig = 0._DP
    CALL DRC3JM(jml1(i), jml2(i), jml3(i), jmm1(i), lmn, lmx, wig, nf, ierr)
    CALL em('drc3jm:mn', i, lmn)
    CALL em('drc3jm:mx', i, lmx)
    DO k = 1, nf
      CALL em('drc3jm', i*1000 + k, wig(k))
    END DO
  END DO

  ! ---- drc6j
  DO i = 1, SIZE(sjl2)
    nf = INT(MIN(sjl2(i)+sjl3(i), sjl5(i)+sjl6(i)) &
             - MAX(ABS(sjl2(i)-sjl3(i)), ABS(sjl5(i)-sjl6(i))) + 1._DP + 0.01_DP)
    CALL em('drc6j:l2', i, sjl2(i))
    CALL em('drc6j:l3', i, sjl3(i))
    CALL em('drc6j:l4', i, sjl4(i))
    CALL em('drc6j:l5', i, sjl5(i))
    CALL em('drc6j:l6', i, sjl6(i))
    CALL emi('drc6j:nd', i, nf)
    wig = 0._DP
    CALL DRC6J(sjl2(i), sjl3(i), sjl4(i), sjl5(i), sjl6(i), lmn, lmx, wig, nf, ierr)
    CALL em('drc6j:mn', i, lmn)
    CALL em('drc6j:mx', i, lmx)
    DO k = 1, nf
      CALL em('drc6j', i*1000 + k, wig(k))
    END DO
  END DO
END PROGRAM fdrv
