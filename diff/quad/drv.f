      PROGRAM DRV
      IMPLICIT NONE
      DOUBLE PRECISION X(401), Y(401), ANS, ERR, QUAD, TOL
      DOUBLE PRECISION C1(3,1), XI1(2), C2(2,2), XI2(3)
      INTEGER IERR, I, N, K
      DOUBLE PRECISION FXSQ, FSIN, FEXPN, FLOR, FX14, FX15, FONE, FID
      DOUBLE PRECISION FINVS, FDAMP, FEXP, FINVX, FSQRT, FX
      EXTERNAL FXSQ, FSIN, FEXPN, FLOR, FX14, FX15, FONE, FID
      EXTERNAL FINVS, FDAMP, FEXP, FINVX, FSQRT, FX

      ERR = 1.0D-12
      CALL DGAUS8(FXSQ, 0.0D0, 1.0D0, ERR, ANS, IERR)
      CALL P('g8 xsq', ANS, IERR)
      ERR = 1.0D-12
      CALL DGAUS8(FSIN, 0.0D0, 3.14159265358979323846D0, ERR, ANS, IERR)
      CALL P('g8 sin', ANS, IERR)
      ERR = 1.0D-12
      CALL DGAUS8(FEXPN, 0.0D0, 1.0D0, ERR, ANS, IERR)
      CALL P('g8 expn', ANS, IERR)
      ERR = 1.0D-12
      CALL DGAUS8(FLOR, 0.0D0, 1.0D0, ERR, ANS, IERR)
      CALL P('g8 lor', ANS, IERR)
      ERR = 1.0D-12
      CALL DGAUS8(FX14, -1.0D0, 1.0D0, ERR, ANS, IERR)
      CALL P('g8 x14', ANS, IERR)
      ERR = 1.0D-12
      CALL DGAUS8(FX15, -1.0D0, 1.0D0, ERR, ANS, IERR)
      CALL P('g8 x15', ANS, IERR)
      ERR = 1.0D-12
      CALL DGAUS8(FSIN, 3.14159265358979323846D0, 0.0D0, ERR, ANS, IERR)
      CALL P('g8 rev', ANS, IERR)
      ERR = 1.0D-14
      CALL DGAUS8(FINVS, 0.0D0, 1.0D0, ERR, ANS, IERR)
      CALL P('g8 invs', ANS, IERR)
      ERR = -1.0D-10
      CALL DGAUS8(FEXP, 0.0D0, 2.0D0, ERR, ANS, IERR)
      CALL P('g8 exp', ANS, IERR)
      CALL P('g8 expce', ERR, IERR)
      ERR = -1.0D-10
      CALL DGAUS8(FINVX, 1.0D0, 10.0D0, ERR, ANS, IERR)
      CALL P('g8 invx', ANS, IERR)
      CALL P('g8 invxce', ERR, IERR)
      ERR = -1.0D-10
      CALL DGAUS8(FSQRT, 0.0D0, 1.0D0, ERR, ANS, IERR)
      CALL P('g8 sqrt', ANS, IERR)
      CALL P('g8 sqrtce', ERR, IERR)
      ERR = -1.0D0
      CALL DGAUS8(FSIN, 1.0D0 - 1.0D-15, 1.0D0, ERR, ANS, IERR)
      CALL P('g8 near', ANS, IERR)
      ERR = 1.0D-4
      CALL DGAUS8(FINVX, 1.0D0, 10.0D0, ERR, ANS, IERR)
      CALL P('g8 coarse', ANS, IERR)
      ERR = 1.0D-14
      CALL DGAUS8(FINVX, 1.0D0, 10.0D0, ERR, ANS, IERR)
      CALL P('g8 fine', ANS, IERR)

      CALL DQNC79(FXSQ, 0.0D0, 1.0D0, 1.0D-12, ANS, IERR, K)
      CALL PK('q7 xsq', ANS, IERR, K)
      CALL DQNC79(FSIN, 0.0D0, 3.14159265358979323846D0, 1.0D-12, ANS,
     +   IERR, K)
      CALL PK('q7 sin', ANS, IERR, K)
      CALL DQNC79(FEXPN, 0.0D0, 1.0D0, 1.0D-12, ANS, IERR, K)
      CALL PK('q7 expn', ANS, IERR, K)
      CALL DQNC79(FLOR, 0.0D0, 1.0D0, 1.0D-12, ANS, IERR, K)
      CALL PK('q7 lor', ANS, IERR, K)
      CALL DQNC79(FX15, 0.0D0, 1.0D0, 1.0D-12, ANS, IERR, K)
      CALL PK('q7 x15', ANS, IERR, K)
      CALL DQNC79(FSIN, 3.14159265358979323846D0, 0.0D0, 1.0D-12, ANS,
     +   IERR, K)
      CALL PK('q7 rev', ANS, IERR, K)
      CALL DQNC79(FINVS, 0.0D0, 1.0D0, 1.0D-14, ANS, IERR, K)
      CALL PK('q7 invs', ANS, IERR, K)
      CALL DQNC79(FDAMP, 0.0D0, 3.0D0, 1.0D-13, ANS, IERR, K)
      CALL PK('q7 damp', ANS, IERR, K)
      CALL DQNC79(FDAMP, 0.0D0, 1.25D0, 1.0D-13, ANS, IERR, K)
      CALL PK('q7 dampl', ANS, IERR, K)
      CALL DQNC79(FDAMP, 1.25D0, 3.0D0, 1.0D-13, ANS, IERR, K)
      CALL PK('q7 dampr', ANS, IERR, K)
      CALL DQNC79(FSIN, 1.0D0, 1.0D0, 1.0D-12, ANS, IERR, K)
      CALL PK('q7 same', ANS, IERR, K)
      CALL DQNC79(FSIN, 1.0D0 - 1.0D-15, 1.0D0, 1.0D-12, ANS, IERR, K)
      CALL PK('q7 near', ANS, IERR, K)
      CALL DQNC79(FEXPN, 0.0D0, 2.0D0, 1.0D-12, ANS, IERR, K)
      CALL PK('q7 gauss', ANS, IERR, K)

      DO 10 I = 1, 5
        X(I) = (I-1)*0.25D0
        Y(I) = 1.0D0 + 2.0D0*X(I) + 3.0D0*X(I)*X(I)
   10 CONTINUE
      CALL DAVINT(X, Y, 5, 0.0D0, 1.0D0, ANS, IERR)
      CALL P('av par', ANS, IERR)
      CALL DAVINT(X, Y, 5, 0.5D0, 0.5D0, ANS, IERR)
      CALL P('av same', ANS, IERR)
      CALL DAVINT(X, Y, 5, 0.75D0, 0.25D0, ANS, IERR)
      CALL P('av bad', ANS, IERR)
      DO 20 I = 1, 9
        X(I) = (I-1)*0.125D0
        Y(I) = 1.0D0 + 2.0D0*X(I) + 3.0D0*X(I)*X(I)
   20 CONTINUE
      CALL DAVINT(X, Y, 9, 0.1D0, 0.9D0, ANS, IERR)
      CALL P('av int', ANS, IERR)
      X(1) = 0.0D0
      X(2) = 0.1D0
      X(3) = 0.35D0
      X(4) = 0.4D0
      X(5) = 0.72D0
      X(6) = 0.9D0
      X(7) = 1.0D0
      DO 30 I = 1, 7
        Y(I) = 1.0D0 + 2.0D0*X(I) + 3.0D0*X(I)*X(I)
   30 CONTINUE
      CALL DAVINT(X, Y, 7, 0.0D0, 1.0D0, ANS, IERR)
      CALL P('av une', ANS, IERR)
      X(1) = 0.0D0
      X(2) = 1.0D0
      Y(1) = 3.0D0
      Y(2) = 5.0D0
      CALL DAVINT(X, Y, 2, 0.0D0, 1.0D0, ANS, IERR)
      CALL P('av tr1', ANS, IERR)
      CALL DAVINT(X, Y, 2, 0.5D0, 2.0D0, ANS, IERR)
      CALL P('av tr2', ANS, IERR)
      X(1) = 0.0D0
      X(2) = 1.0D0
      X(3) = 2.0D0
      X(4) = 3.0D0
      DO 40 I = 1, 4
        Y(I) = 1.0D0 + 2.0D0*X(I) + 3.0D0*X(I)*X(I)
   40 CONTINUE
      CALL DAVINT(X, Y, 4, 2.5D0, 3.0D0, ANS, IERR)
      CALL P('av f3a', ANS, IERR)
      CALL DAVINT(X, Y, 4, 0.0D0, 0.5D0, ANS, IERR)
      CALL P('av f3b', ANS, IERR)
      X(3) = 1.0D0
      CALL DAVINT(X, Y, 4, 0.0D0, 3.0D0, ANS, IERR)
      CALL P('av rep', ANS, IERR)
      X(2) = 2.0D0
      X(3) = 1.0D0
      CALL DAVINT(X, Y, 4, 0.0D0, 3.0D0, ANS, IERR)
      CALL P('av dec', ANS, IERR)
      CALL DAVINT(X, Y, 1, 0.0D0, 1.0D0, ANS, IERR)
      CALL P('av n1', ANS, IERR)
      N = 401
      DO 50 I = 1, N
        X(I) = (I-1)*2.0D0/(N-1)
        Y(I) = EXP(-X(I))*COS(3.0D0*X(I))
   50 CONTINUE
      CALL DAVINT(X, Y, N, 0.0D0, 2.0D0, ANS, IERR)
      CALL P('av big', ANS, IERR)

      C1(1,1) = 1.0D0
      C1(2,1) = 2.0D0
      C1(3,1) = 6.0D0
      XI1(1) = 0.0D0
      XI1(2) = 1.0D0
      ERR = 1.0D-12
      CALL DPPGQ8(FONE, 3, C1, XI1, 1, 3, 0, 0.0D0, 1.0D0, 1, ERR,
     +   ANS, IERR)
      CALL P('pg id0', ANS, IERR)
      ERR = 1.0D-12
      CALL DPPGQ8(FONE, 3, C1, XI1, 1, 3, 1, 0.0D0, 1.0D0, 1, ERR,
     +   ANS, IERR)
      CALL P('pg id1', ANS, IERR)
      ERR = 1.0D-12
      CALL DPPGQ8(FONE, 3, C1, XI1, 1, 3, 2, 0.0D0, 1.0D0, 1, ERR,
     +   ANS, IERR)
      CALL P('pg id2', ANS, IERR)
      ERR = 1.0D-12
      CALL DPPGQ8(FID, 3, C1, XI1, 1, 3, 0, 0.0D0, 1.0D0, 1, ERR,
     +   ANS, IERR)
      CALL P('pg wt', ANS, IERR)
      ERR = 1.0D-12
      CALL DPPGQ8(FONE, 3, C1, XI1, 1, 3, 0, 1.0D0, 0.0D0, 1, ERR,
     +   ANS, IERR)
      CALL P('pg rev', ANS, IERR)
      ERR = -1.0D-10
      CALL DPPGQ8(FEXP, 3, C1, XI1, 1, 3, 0, 0.0D0, 1.0D0, 1, ERR,
     +   ANS, IERR)
      CALL P('pg exp', ANS, IERR)
      CALL P('pg expbe', ERR, IERR)
      ERR = 1.0D-12
      CALL DPPGQ8(FEXP, 3, C1, XI1, 1, 3, 0, 0.0D0, 1.0D0, 1, ERR,
     +   ANS, IERR)
      CALL P('pg exp2', ANS, IERR)

      C2(1,1) = 1.0D0
      C2(2,1) = 2.0D0
      C2(1,2) = 3.0D0
      C2(2,2) = 4.0D0
      XI2(1) = 0.0D0
      XI2(2) = 1.0D0
      XI2(3) = 2.0D0
      TOL = 1.0D-12
      CALL DPFQAD(FONE, 2, C2, XI2, 2, 2, 0, 0.0D0, 2.0D0, TOL, QUAD,
     +   IERR)
      CALL P('pf all', QUAD, IERR)
      TOL = 1.0D-12
      CALL DPFQAD(FONE, 2, C2, XI2, 2, 2, 0, 0.0D0, 0.5D0, TOL, QUAD,
     +   IERR)
      CALL P('pf one', QUAD, IERR)
      TOL = 1.0D-12
      CALL DPFQAD(FONE, 2, C2, XI2, 2, 2, 0, 0.5D0, 1.5D0, TOL, QUAD,
     +   IERR)
      CALL P('pf strd', QUAD, IERR)
      TOL = 1.0D-12
      CALL DPFQAD(FID, 2, C2, XI2, 2, 2, 0, 0.0D0, 2.0D0, TOL, QUAD,
     +   IERR)
      CALL P('pf wt', QUAD, IERR)
      TOL = 1.0D-12
      CALL DPFQAD(FONE, 2, C2, XI2, 2, 2, 0, 2.0D0, 0.0D0, TOL, QUAD,
     +   IERR)
      CALL P('pf rev', QUAD, IERR)
      TOL = 1.0D-12
      CALL DPFQAD(FONE, 2, C2, XI2, 2, 2, 0, 1.0D0, 1.0D0, TOL, QUAD,
     +   IERR)
      CALL P('pf same', QUAD, IERR)
      TOL = 1.0D-12
      CALL DPFQAD(FONE, 2, C2, XI2, 2, 2, 1, 0.0D0, 1.0D0, TOL, QUAD,
     +   IERR)
      CALL P('pf der', QUAD, IERR)
      TOL = 1.0D-12
      CALL DPFQAD(FEXP, 3, C1, XI1, 1, 3, 0, 0.0D0, 1.0D0, TOL, QUAD,
     +   IERR)
      CALL P('pf one1', QUAD, IERR)
      END

      SUBROUTINE P (TAG, V, IERR)
      CHARACTER*(*) TAG
      DOUBLE PRECISION V
      INTEGER IERR
      WRITE (*,'(A,1X,Z16,1X,I3)') TAG, V, IERR
      RETURN
      END

      SUBROUTINE PK (TAG, V, IERR, K)
      CHARACTER*(*) TAG
      DOUBLE PRECISION V
      INTEGER IERR, K
      WRITE (*,'(A,1X,Z16,1X,I3,1X,I6)') TAG, V, IERR, K
      RETURN
      END

      DOUBLE PRECISION FUNCTION FXSQ (X)
      DOUBLE PRECISION X
      FXSQ = X*X
      RETURN
      END
      DOUBLE PRECISION FUNCTION FSIN (X)
      DOUBLE PRECISION X
      FSIN = SIN(X)
      RETURN
      END
      DOUBLE PRECISION FUNCTION FEXPN (X)
      DOUBLE PRECISION X
      FEXPN = EXP(-X)
      RETURN
      END
      DOUBLE PRECISION FUNCTION FLOR (X)
      DOUBLE PRECISION X
      FLOR = 1.0D0/(1.0D0 + X*X)
      RETURN
      END
      DOUBLE PRECISION FUNCTION FX14 (X)
      DOUBLE PRECISION X
      FX14 = X**14.0D0
      RETURN
      END
      DOUBLE PRECISION FUNCTION FX15 (X)
      DOUBLE PRECISION X
      FX15 = X**15.0D0
      RETURN
      END
      DOUBLE PRECISION FUNCTION FONE (X)
      DOUBLE PRECISION X
      FONE = 1.0D0
      RETURN
      END
      DOUBLE PRECISION FUNCTION FID (X)
      DOUBLE PRECISION X
      FID = X
      RETURN
      END
      DOUBLE PRECISION FUNCTION FINVS (X)
      DOUBLE PRECISION X
      IF (X .GT. 0.0D0) THEN
        FINVS = 1.0D0/SQRT(X)
      ELSE
        FINVS = 0.0D0
      END IF
      RETURN
      END
      DOUBLE PRECISION FUNCTION FDAMP (X)
      DOUBLE PRECISION X
      FDAMP = EXP(-X)*COS(X)
      RETURN
      END
      DOUBLE PRECISION FUNCTION FEXP (X)
      DOUBLE PRECISION X
      FEXP = EXP(X)
      RETURN
      END
      DOUBLE PRECISION FUNCTION FINVX (X)
      DOUBLE PRECISION X
      FINVX = 1.0D0/X
      RETURN
      END
      DOUBLE PRECISION FUNCTION FSQRT (X)
      DOUBLE PRECISION X
      FSQRT = SQRT(X)
      RETURN
      END
      DOUBLE PRECISION FUNCTION FX (X)
      DOUBLE PRECISION X
      FX = X
      RETURN
      END
