! Reeve differential, interpolation. Pass one writes every input of the sweep
! to cases.hex in hex, pass two reads that file back and runs
! the cases through the modernised DPLINT, DPOLCF and DPOLVL, dumping results,
! work arrays and the XERMSG call count in hex on stdout. The original
! F77 is run alongside on every case and the two are compared internally; the
! count of disagreements goes to stderr.

PROGRAM drv
  USE service, ONLY : DP
  USE xerstate, ONLY : xercnt
  USE slatec_interp, ONLY : m_dplint => DPLINT, m_dpolcf => DPOLCF, &
    m_dpolvl => DPOLVL
  IMPLICIT NONE

  EXTERNAL :: DPLINT, DPOLCF, DPOLVL

  INTEGER, PARAMETER :: NMAX = 64, NW = 2*NMAX + 8, NS = 15
  INTEGER, PARAMETER :: sizes(NS) = [1,2,3,4,5,6,7,8,10,12,16,20,24,32,40]
  REAL(DP), PARAMETER :: pi = 3.14159265358979323846264338327950288_DP
  INTEGER, PARAMETER :: CU = 11

  INTEGER :: ncase, nmis

  ncase = 0
  nmis = 0
  CALL gen()
  CALL sweep()
  WRITE(0,'(A,I8)') 'CASES ', ncase
  WRITE(0,'(A,I8)') 'MODVSORIG ', nmis

CONTAINS

  ! ---- Inputs ----

  SUBROUTINE mkx(n, xfam, x)
    INTEGER, INTENT(IN) :: n, xfam
    REAL(DP), INTENT(INOUT) :: x(NW)
    INTEGER :: i, j
    SELECT CASE( xfam )
    CASE( 1, 4, 5 )
      DO i = 1, n
        x(i) = -1._DP + (i-1)*0.25_DP
      END DO
      IF( xfam==4 .AND. n>=2 ) THEN
        j = MIN(3,n)
        x(j) = NEAREST(x(j-1), 1._DP)
      END IF
      IF( xfam==5 .AND. n>=2 ) THEN
        j = MIN(4,n)
        x(j) = x(MAX(j-2,1))
      END IF
    CASE( 2 )
      IF( n>=1 ) x(1) = -1._DP
      DO i = 2, n
        x(i) = x(i-1) + (1 + MOD(i*i,5))/8._DP
      END DO
    CASE( 3 )
      DO i = 1, n
        x(i) = COS(pi*(2*i-1)/(2._DP*n))
      END DO
    END SELECT
  END SUBROUTINE mkx

  SUBROUTINE mky(n, yfam, x, y)
    INTEGER, INTENT(IN) :: n, yfam
    REAL(DP), INTENT(IN) :: x(NW)
    REAL(DP), INTENT(INOUT) :: y(NW)
    INTEGER :: i
    INTEGER(8) :: s
    SELECT CASE( yfam )
    CASE( 1 )
      DO i = 1, n
        y(i) = ((((-0.25_DP*x(i) + 0.5_DP)*x(i) + 1._DP)*x(i) - 2._DP)*x(i) &
          + 3._DP)
      END DO
    CASE( 2 )
      s = 123456789_8 + 7_8*n
      DO i = 1, n
        s = MOD(1103515245_8*s + 12345_8, 2147483648_8)
        y(i) = (MOD(s/65536_8, 2048_8) - 1024_8)/16._DP
      END DO
    CASE( 3 )
      DO i = 1, n
        y(i) = (1 + MOD(i,3))*(1 - 2*MOD(i,2))
      END DO
    END SELECT
  END SUBROUTINE mky

  SUBROUTINE mkz(n, x, z)
    INTEGER, INTENT(IN) :: n
    REAL(DP), INTENT(IN) :: x(NW)
    REAL(DP), INTENT(OUT) :: z(3)
    IF( n>=1 ) THEN
      z(1) = x(1)
      z(2) = 0.5_DP*(x(1) + x(n))
      z(3) = x(n) + 0.75_DP
    ELSE
      z(1) = 0._DP
      z(2) = 0.5_DP
      z(3) = -0.75_DP
    END IF
  END SUBROUTINE mkz

  SUBROUTINE nders(n, nd, cnt)
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(OUT) :: nd(8), cnt
    INTEGER :: cand(8), i, j
    LOGICAL :: dup
    cand = [0, 1, 2, 3, 5, n-1, n, n+2]
    cnt = 0
    nd = 0
    DO i = 1, 8
      IF( cand(i)<0 ) CYCLE
      dup = .FALSE.
      DO j = 1, cnt
        IF( nd(j)==cand(i) ) dup = .TRUE.
      END DO
      IF( .NOT. dup ) THEN
        cnt = cnt + 1
        nd(cnt) = cand(i)
      END IF
    END DO
  END SUBROUTINE nders

  ! ---- cases.hex ----

  SUBROUTINE wrec(id, swp, n, nder, errc, xx, z, x, y)
    INTEGER, INTENT(IN) :: id, swp, n, nder, errc
    REAL(DP), INTENT(IN) :: xx, z(3), x(NW), y(NW)
    INTEGER :: i
    WRITE(CU,"(5I8)") id, swp, n, nder, errc
    WRITE(CU,'(Z16.16)') TRANSFER(xx, 0_8)
    DO i = 1, 3
      WRITE(CU,'(Z16.16)') TRANSFER(z(i), 0_8)
    END DO
    DO i = 1, MAX(n,0)
      WRITE(CU,'(Z16.16)') TRANSFER(x(i), 0_8)
    END DO
    DO i = 1, MAX(n,0)
      WRITE(CU,'(Z16.16)') TRANSFER(y(i), 0_8)
    END DO
  END SUBROUTINE wrec

  SUBROUTINE gen()
    INTEGER :: is, n, xf, yf, ixx, id, errc, nd(8), cnt, k, i, j
    REAL(DP) :: x(NW), y(NW), z(3)
    LOGICAL :: hasdup
    OPEN(CU, FILE='cases.hex', STATUS='REPLACE', ACTION='WRITE')
    id = 0

    ! dplint, every abscissa family including the repeated pair
    DO is = 1, NS
      n = sizes(is)
      DO xf = 1, 5
        DO yf = 1, 3
          x = 0._DP
          y = 0._DP
          CALL mkx(n, xf, x)
          CALL mky(n, yf, x, y)
          CALL mkz(n, x, z)
          hasdup = .FALSE.
          DO i = 1, n
            DO j = i+1, n
              IF( x(i)==x(j) ) hasdup = .TRUE.
            END DO
          END DO
          errc = 0
          IF( hasdup ) errc = 1
          id = id + 1
          CALL wrec(id, 1, n, 0, errc, z(1), z, x, y)
        END DO
      END DO
    END DO

    ! dplint with a bad N, which the original reports through XERMSG and the
    ! modernised .inc turns into ERROR STOP
    DO k = 1, 2
      n = -3*(k-1)
      x = 0._DP
      y = 0._DP
      CALL mkz(n, x, z)
      id = id + 1
      CALL wrec(id, 1, n, 0, 1, z(1), z, x, y)
    END DO

    ! dpolvl
    DO is = 1, NS
      n = sizes(is)
      CALL nders(n, nd, cnt)
      DO xf = 1, 4
        DO yf = 1, 2
          x = 0._DP
          y = 0._DP
          CALL mkx(n, xf, x)
          CALL mky(n, yf, x, y)
          CALL mkz(n, x, z)
          DO ixx = 1, 3
            DO k = 1, cnt
              id = id + 1
              CALL wrec(id, 2, n, nd(k), 0, z(ixx), z, x, y)
            END DO
          END DO
        END DO
      END DO
    END DO

    ! dpolcf
    DO is = 1, NS
      n = sizes(is)
      DO xf = 1, 4
        DO yf = 1, 2
          x = 0._DP
          y = 0._DP
          CALL mkx(n, xf, x)
          CALL mky(n, yf, x, y)
          CALL mkz(n, x, z)
          DO ixx = 1, 3
            id = id + 1
            CALL wrec(id, 3, n, 0, 0, z(ixx), z, x, y)
          END DO
        END DO
      END DO
    END DO

    ncase = id
    CLOSE(CU)
  END SUBROUTINE gen

  ! ---- Output ----

  SUBROUTINE dmp(tag, id, idx, v)
    CHARACTER(*), INTENT(IN) :: tag
    INTEGER, INTENT(IN) :: id, idx
    REAL(DP), INTENT(IN) :: v
    CHARACTER(8) :: t
    t = tag
    WRITE(*,'(A8,1X,I6,1X,I6,1X,Z16.16)') t, id, idx, TRANSFER(v, 0_8)
  END SUBROUTINE dmp

  SUBROUTINE dmpi(tag, id, idx, v)
    CHARACTER(*), INTENT(IN) :: tag
    INTEGER, INTENT(IN) :: id, idx, v
    CHARACTER(8) :: t
    t = tag
    WRITE(*,'(A8,1X,I6,1X,I6,1X,Z16.16)') t, id, idx, INT(v, 8)
  END SUBROUTINE dmpi

  SUBROUTINE fill(a)
    REAL(DP), INTENT(OUT) :: a(NW)
    INTEGER :: i
    DO i = 1, NW
      a(i) = -(1024 + i)/8._DP
    END DO
  END SUBROUTINE fill

  INTEGER FUNCTION cmpv(a, b, m)
    REAL(DP), INTENT(IN) :: a(NW), b(NW)
    INTEGER, INTENT(IN) :: m
    INTEGER :: i
    cmpv = 0
    DO i = 1, m
      IF( TRANSFER(a(i), 0_8)/=TRANSFER(b(i), 0_8) ) cmpv = cmpv + 1
    END DO
  END FUNCTION cmpv

  ! ---- The sweep ----

  SUBROUTINE sweep()
    INTEGER :: id, swp, n, nder, errc, i, k, ierr, ierro, nd1
    INTEGER(8) :: ib
    REAL(DP) :: xx, z(3), x(NW), y(NW), cm(NW), co(NW)
    REAL(DP) :: work(NW), worko(NW), yp(NW), ypo(NW), d(NW), dd(NW)
    REAL(DP) :: yfit, yfito, acc, w2(NW), ypz(NW)

    OPEN(CU, FILE='cases.hex', STATUS='OLD', ACTION='READ')
    DO
      READ(CU,"(5I8)",END=900) id, swp, n, nder, errc
      READ(CU,'(Z16)') ib
      xx = TRANSFER(ib, 0._DP)
      DO i = 1, 3
        READ(CU,'(Z16)') ib
        z(i) = TRANSFER(ib, 0._DP)
      END DO
      x = 0._DP
      y = 0._DP
      DO i = 1, MAX(n,0)
        READ(CU,'(Z16)') ib
        x(i) = TRANSFER(ib, 0._DP)
      END DO
      DO i = 1, MAX(n,0)
        READ(CU,'(Z16)') ib
        y(i) = TRANSFER(ib, 0._DP)
      END DO

      SELECT CASE( swp )

      CASE( 1 )
        CALL fill(cm)
        CALL fill(co)
        xercnt = 0
        IF( errc==0 ) THEN
          CALL m_dplint(n, x, y, cm)
          CALL DPLINT(n, x, y, co)
          nmis = nmis + cmpv(cm, co, n)
        ELSE
          CALL DPLINT(n, x, y, co)
          cm = co
        END IF
        DO i = 1, MAX(n,1)
          CALL dmp('a.c', id, i, cm(i))
        END DO
        CALL dmpi('a.xerr', id, 0, xercnt)

      CASE( 2 )
        CALL fill(cm)
        CALL m_dplint(n, x, y, cm)
        nd1 = MAX(nder,1)
        CALL fill(work)
        CALL fill(yp)
        CALL fill(worko)
        CALL fill(ypo)
        yfit = -7._DP
        yfito = -7._DP
        ierr = -7
        ierro = -7
        CALL m_dpolvl(nder, xx, yfit, yp, n, x, cm, work, ierr)
        CALL DPOLVL(nder, xx, yfito, ypo, n, x, cm, worko, ierro)
        IF( TRANSFER(yfit, 0_8)/=TRANSFER(yfito, 0_8) ) nmis = nmis + 1
        IF( ierr/=ierro ) nmis = nmis + 1
        nmis = nmis + cmpv(yp, ypo, nd1)
        nmis = nmis + cmpv(work, worko, 2*n)
        CALL dmp('b.yfit', id, 0, yfit)
        CALL dmpi('b.ierr', id, 0, ierr)
        DO i = 1, nd1
          CALL dmp('b.yp', id, i, yp(i))
        END DO
        DO i = 1, 2*n
          CALL dmp('b.work', id, i, work(i))
        END DO

      CASE( 3 )
        CALL fill(cm)
        CALL m_dplint(n, x, y, cm)
        CALL fill(work)
        CALL fill(d)
        CALL fill(worko)
        CALL fill(dd)
        CALL m_dpolcf(xx, n, x, cm, d, work)
        CALL DPOLCF(xx, n, x, cm, dd, worko)
        nmis = nmis + cmpv(d, dd, n)
        nmis = nmis + cmpv(work, worko, 2*n)
        DO i = 1, n
          CALL dmp('c.d', id, i, d(i))
        END DO
        DO i = 1, 2*n
          CALL dmp('c.work', id, i, work(i))
        END DO
        DO k = 1, 3
          acc = d(n)
          DO i = n-1, 1, -1
            acc = acc*(z(k) - xx) + d(i)
          END DO
          CALL dmp('c.horn', id, k, acc)
          CALL fill(w2)
          CALL fill(ypz)
          yfit = -7._DP
          ierr = -7
          CALL m_dpolvl(0, z(k), yfit, ypz, n, x, cm, w2, ierr)
          CALL dmp('c.yfit', id, k, yfit)
        END DO

      END SELECT
    END DO
900 CLOSE(CU)
  END SUBROUTINE sweep

END PROGRAM drv
