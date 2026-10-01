PROGRAM drv2
  USE service
  USE sf
  IMPLICIT NONE
  CHARACTER(LEN=32) :: name
  CHARACTER(LEN=32) :: hx
  INTEGER(8) :: bits
  REAL(DP) :: x, r
  INTEGER :: ios, n, i, nz, ierr
  REAL(DP) :: en(20)
  EXTERNAL :: DEXINTO
  DO
    READ(*,*,IOSTAT=ios) name, hx, n
    IF( ios/=0 ) EXIT
    READ(hx,'(Z16)') bits
    x = TRANSFER(bits,x)
    SELECT CASE( TRIM(name) )
    CASE('de1');    r = DE1(x)
    CASE('dei');    r = DEI(x)
    CASE('dli');    r = DLI(x)
    CASE('dspenc'); r = DSPENC(x)
    CASE('ddaws');  r = DDAWS(x)
    CASE('dexprl'); r = DEXPRL(x)
    CASE('dlnrel'); r = DLNREL(x)
    CASE('dcbrt');  r = DCBRT(x)
    CASE('dsindg'); r = DSINDG(x)
    CASE('dcosdg'); r = DCOSDG(x)
    CASE('d9atn1'); r = D9ATN1(x)
    CASE('d9ln2r'); r = D9LN2R(x)
    CASE('d9pak');  r = D9PAK(x,n)
    CASE('d9upak')
      CALL D9UPAK(x,r,i)
      WRITE(*,'(A,1X,Z16.16,1X,I0)') 'd9upak_n', TRANSFER(r,bits), i
    CASE('dexint')
      en = -1._DP
      CALL DEXINTO(x,n,1,4,1.0D-14,en,nz,ierr)
      WRITE(*,'(A,1X,I0,1X,I0)') 'dexint_flags', nz, ierr
      IF( ierr==0 ) THEN
        DO i = 1, 4
          WRITE(*,'(A,I0,1X,Z16.16)') 'dexint_', i, TRANSFER(en(i),bits)
        END DO
      END IF
      CYCLE
    CASE('dexint2')
      en = -1._DP
      CALL DEXINTO(x,n,2,4,1.0D-14,en,nz,ierr)
      WRITE(*,'(A,1X,I0,1X,I0)') 'dexint2_flags', nz, ierr
      IF( ierr==0 ) THEN
        DO i = 1, 4
          WRITE(*,'(A,I0,1X,Z16.16)') 'dexint2_', i, TRANSFER(en(i),bits)
        END DO
      END IF
      CYCLE
    CASE DEFAULT
      WRITE(*,*) 'unknown ', TRIM(name)
      CYCLE
    END SELECT
    WRITE(*,'(A,1X,Z16.16)') TRIM(name), TRANSFER(r,bits)
  END DO
END PROGRAM drv2
