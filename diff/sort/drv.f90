subroutine gen(n, sel, dx)
  implicit none
  integer :: n, sel, i
  double precision :: dx(*)
  integer(8) :: s
  s = 123456789_8 + int(sel,8)*7_8
  do i = 1, n
    s = mod(1103515245_8*s + 12345_8, 2147483648_8)
    if (sel == 1) then
      dx(i) = dble(mod(s/65536_8, 1000_8) - 500_8)/8.0d0
    else if (sel == 2) then
      dx(i) = dble(mod(s/65536_8, 7_8))
    else if (sel == 3) then
      dx(i) = 2.5d0
    else if (sel == 4) then
      dx(i) = dble(i)/4.0d0
    else if (sel == 5) then
      dx(i) = dble(n-i+1)/4.0d0
    else
      dx(i) = dble(min(i, n-i+1))/2.0d0
    end if
  end do
end subroutine gen

subroutine dumpd(tag, n, a)
  implicit none
  character(len=*) :: tag
  integer :: n, i
  double precision :: a(*)
  integer(8) :: b
  do i = 1, n
    b = transfer(a(i), b)
    write(*,'(A,1X,I6,1X,Z16.16)') tag, i, b
  end do
end subroutine dumpd

subroutine dumpi(tag, n, a)
  implicit none
  character(len=*) :: tag
  integer :: n, i
  integer :: a(*)
  do i = 1, n
    write(*,'(A,1X,I6,1X,I8)') tag, i, a(i)
  end do
end subroutine dumpi

program drv
  implicit none
  integer, parameter :: nmax = 2000
  double precision :: dx(nmax), dy(nmax), tc(3*nmax)
  integer :: iperm(nmax), ier, n, kf, i, m1, m2, c, sel, kflag
  integer, parameter :: sizes(12) = (/1,2,3,5,8,11,12,13,100,257,1000,2000/)
  integer, parameter :: flags(4) = (/1,2,-1,-2/)
  character(len=64) :: tag

  do sel = 1, 6
    do c = 1, 12
      n = sizes(c)
      do kf = 1, 4
        kflag = flags(kf)
        call gen(n, sel, dx)
        do i = 1, n
          dy(i) = dble(i)
        end do
        call dsort(dx, dy, n, kflag)
        write(tag,'(A,I1,A,I5,A,I2)') 'dsort s', sel, ' n', n, ' k', kflag
        call dumpd(trim(tag)//' x', n, dx)
        call dumpd(trim(tag)//' y', n, dy)

        call gen(n, sel, dx)
        call dpsort(dx, n, iperm, kflag, ier)
        write(tag,'(A,I1,A,I5,A,I2)') 'dpsort s', sel, ' n', n, ' k', kflag
        call dumpd(trim(tag)//' x', n, dx)
        call dumpi(trim(tag)//' p', n, iperm)
        write(*,'(A,1X,I3)') trim(tag)//' ier', ier
      end do

      call gen(n, sel, dx)
      call dpsort(dx, n, iperm, 1, ier)
      call gen(n, sel, dx)
      call dpperm(dx, n, iperm, ier)
      write(tag,'(A,I1,A,I5)') 'dpperm s', sel, ' n', n
      call dumpd(trim(tag)//' x', n, dx)
      call dumpi(trim(tag)//' p', n, iperm)
      write(*,'(A,1X,I3)') trim(tag)//' ier', ier

      m1 = n/2
      m2 = n - m1
      call gen(n, sel, tc)
      if (m1 > 0) call dsort(tc(1), dy, m1, 1)
      call dsort(tc(m1+1), dy, m2, 1)
      do i = n+1, 3*nmax
        tc(i) = -7.5d0
      end do
      call d1merg(tc, 0, m1, m1, m2, n)
      write(tag,'(A,I1,A,I5)') 'd1merg s', sel, ' n', n
      call dumpd(trim(tag)//' t', 2*n, tc)
    end do
  end do
end program drv
