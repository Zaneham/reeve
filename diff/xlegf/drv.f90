! Reeve differential driver for the extended-range Legendre package.
!
! One DXSET per process: the SLATEC DXSET keeps a saved IFLAG and refuses to
! run twice, and DXLEGF and DXNRMP both call DXSET(0,0,0,0) on entry, so a
! non-default configuration has to be installed before the first of those
! calls and cannot be changed afterwards.  The first input line carries the
! configuration, the rest carry one case each, and run.sh starts the program
! once per configuration.
!
! PQA is handed to DXLEGF offset by one slot, with a guard slot before it and
! five after, because DXQNU reads PQA(K-1) with K=1 and can write PQA(0).  The
! guard fill comes from the command line so the same sweep can be run with two
! different fills to show the printed window does not depend on it.
!
! A case counts as an error when XERMSG fired during it, whatever IERROR comes
! back as.  DXPMUP and DXPNRM put their IF (IERROR.NE.0) RETURN after the DO
! loop that calls DXADJ, and DXADJ sets IERROR=0 on entry, so an index
! overflow inside the loop is discarded and the routine returns success with a
! principal part it scaled and an index it did not.  The error did reach
! XERMSG, which is where the port raises, so the driver reads it there.

PROGRAM drv
  IMPLICIT NONE
  INTEGER, PARAMETER :: NMAX = 140
  INTEGER :: irad, nradpl, nbits, ierror, ios, i, nv, iguard
  INTEGER :: ix, iy, nudiff, mu1, mu2, id, nu, mode, isig
  INTEGER :: ipqa(0:NMAX+5)
  DOUBLE PRECISION :: pqa(0:NMAX+5)
  DOUBLE PRECISION :: dzero, x, y, z, dnu1, theta, darg, guard
  INTEGER(8) :: bits
  CHARACTER(LEN=400) :: line
  CHARACTER(LEN=16) :: cmd
  CHARACTER(LEN=16) :: h1, h2
  CHARACTER(LEN=64) :: arg
  INTEGER :: nxerr
  COMMON /XFLAG/ nxerr

  guard = 0.0D0
  iguard = 0
  IF( COMMAND_ARGUMENT_COUNT() >= 1 ) THEN
    CALL GET_COMMAND_ARGUMENT(1, arg)
    READ(arg,'(Z16)') bits
    guard = TRANSFER(bits, guard)
  END IF
  IF( COMMAND_ARGUMENT_COUNT() >= 2 ) THEN
    CALL GET_COMMAND_ARGUMENT(2, arg)
    READ(arg,*) iguard
  END IF

  nxerr = 0
  READ(*,'(A)',IOSTAT=ios) line
  IF( ios /= 0 ) STOP
  READ(line,*) irad, nradpl, h1, nbits
  dzero = hx(h1)
  ierror = 0
  CALL DXSET(irad, nradpl, dzero, nbits, ierror)
  IF( ierror /= 0 ) THEN
    WRITE(*,'(A)') 'set err'
    STOP
  END IF
  WRITE(*,'(A)') 'set ok'

  DO
    READ(*,'(A)',IOSTAT=ios) line
    IF( ios /= 0 ) EXIT
    READ(line,*) cmd
    nxerr = 0
    SELECT CASE( TRIM(cmd) )

    CASE('adj')
      READ(line,*) cmd, h1, ix
      x = hx(h1)
      ierror = 0
      CALL DXADJ(x, ix, ierror)
      CALL emit('adj', x, ix, ierror)

    CASE('add')
      READ(line,*) cmd, h1, ix, h2, iy
      x = hx(h1)
      y = hx(h2)
      ierror = 0
      CALL DXADD(x, ix, y, iy, z, i, ierror)
      CALL emit('add', z, i, ierror)

    CASE('rt')
      READ(line,*) cmd, h1, ix
      x = hx(h1)
      ierror = 0
      CALL DXADJ(x, ix, ierror)
      IF( ierror == 0 ) CALL DXRED(x, ix, ierror)
      CALL emit('rt', x, ix, ierror)

    CASE('red')
      READ(line,*) cmd, h1, ix
      x = hx(h1)
      ierror = 0
      CALL DXRED(x, ix, ierror)
      CALL emit('red', x, ix, ierror)

    CASE('con')
      READ(line,*) cmd, h1, ix
      x = hx(h1)
      ierror = 0
      CALL DXCON(x, ix, ierror)
      CALL emit('con', x, ix, ierror)

    CASE('legf')
      READ(line,*) cmd, h1, nudiff, mu1, mu2, h2, id
      dnu1 = hx(h1)
      theta = hx(h2)
      nv = (mu2 - mu1) + nudiff + 1
      CALL fill()
      ierror = 0
      CALL DXLEGF(dnu1, nudiff, mu1, mu2, theta, id, pqa(1), ipqa(1), ierror)
      CALL vec('legf', nv, ierror, 0, 0)

    CASE('nrmp')
      READ(line,*) cmd, nu, mu1, mu2, h1, mode
      darg = hx(h1)
      nv = mu2 - mu1 + 1
      CALL fill()
      ierror = 0
      isig = 0
      CALL DXNRMP(nu, mu1, mu2, darg, mode, pqa(1), ipqa(1), isig, ierror)
      CALL vec('nrmp', nv, ierror, 1, isig)

    CASE DEFAULT
      WRITE(*,'(A,A)') 'unknown ', TRIM(cmd)
    END SELECT
  END DO

CONTAINS

  DOUBLE PRECISION FUNCTION hx(s)
    CHARACTER(LEN=*), INTENT(IN) :: s
    INTEGER(8) :: b
    DOUBLE PRECISION :: t
    READ(s,'(Z16)') b
    t = TRANSFER(b, t)
    hx = t
  END FUNCTION hx

  SUBROUTINE emit(tag, v, k, e)
    CHARACTER(LEN=*), INTENT(IN) :: tag
    DOUBLE PRECISION, INTENT(IN) :: v
    INTEGER, INTENT(IN) :: k, e
    INTEGER(8) :: b
    IF( e /= 0 .OR. nxerr /= 0 ) THEN
      WRITE(*,'(A,1X,A)') tag, 'err'
    ELSE
      b = TRANSFER(v, b)
      WRITE(*,'(A,1X,Z16.16,1X,I0)') tag, b, k
    END IF
  END SUBROUTINE emit

  SUBROUTINE fill()
    INTEGER :: k
    DO k = 0, NMAX+5
      pqa(k) = guard
      ipqa(k) = iguard
    END DO
  END SUBROUTINE fill

  SUBROUTINE vec(tag, n, e, withsig, sg)
    CHARACTER(LEN=*), INTENT(IN) :: tag
    INTEGER, INTENT(IN) :: n, e, withsig, sg
    INTEGER(8) :: b
    INTEGER :: k
    IF( e /= 0 .OR. nxerr /= 0 ) THEN
      WRITE(*,'(A,1X,A)') tag, 'err'
      RETURN
    END IF
    IF( withsig == 1 ) THEN
      WRITE(*,'(A,1X,I0,1X,I0)') tag, n, sg
    ELSE
      WRITE(*,'(A,1X,I0)') tag, n
    END IF
    DO k = 1, n
      b = TRANSFER(pqa(k), b)
      WRITE(*,'(A,1X,I0,1X,Z16.16,1X,I0)') '.', k, b, ipqa(k)
    END DO
  END SUBROUTINE vec

END PROGRAM drv
