module rndmod
  implicit none
  integer(8) :: sd = 0
contains
  subroutine rseed(v)
    integer, intent(in) :: v
    sd = int(v, 8)
  end subroutine rseed

  double precision function rnd()
    sd = mod(sd*1103515245_8 + 12345_8, 2147483648_8)
    rnd = dble(mod(sd, 20001_8) - 10000_8)/1000.0d0
    if (rnd == 0.0d0) rnd = 0.25d0
  end function rnd

  subroutine pr(x)
    double precision, intent(in) :: x
    integer(8) :: b
    b = transfer(x, b)
    write(*,'(Z16.16)') b
  end subroutine pr

  subroutine pin(i)
    integer, intent(in) :: i
    write(*,'(I0)') i
  end subroutine pin

  subroutine lbl(c, i1, i2, i3, i4, i5, i6)
    integer, intent(in) :: c, i1, i2, i3, i4, i5, i6
    write(*,'(I0,6(1X,I0))') c, i1, i2, i3, i4, i5, i6
  end subroutine lbl

  subroutine fillb(x, n)
    integer, intent(in) :: n
    double precision, intent(inout) :: x(*)
    integer :: i
    do i = 1, n
      x(i) = rnd()
    end do
  end subroutine fillb

  subroutine dumpb(x, n)
    integer, intent(in) :: n
    double precision, intent(in) :: x(*)
    integer :: i
    do i = 1, n
      call pr(x(i))
    end do
  end subroutine dumpb

  subroutine fillspd(x, nbt, off, ld, n)
    integer, intent(in) :: nbt, off, ld, n
    double precision, intent(inout) :: x(*)
    integer :: i, j
    double precision :: v
    call fillb(x, nbt)
    do j = 1, n
      do i = j + 1, n
        v = rnd()
        x(1 + off + (i-1) + (j-1)*ld) = v
        x(1 + off + (j-1) + (i-1)*ld) = v
      end do
    end do
    do j = 1, n
      x(1 + off + (j-1) + (j-1)*ld) = 10.0d0*dble(n) + dble(j)
    end do
  end subroutine fillspd
end module rndmod

subroutine tgemm
  use rndmod
  implicit none
  external dgemm
  integer, parameter :: ld = 8, nb = 64
  double precision :: a(nb), b(nb), c(nb)
  character :: tc(3) = (/ 'N', 'T', 'C' /)
  integer :: dm(3) = (/ 0, 1, 3 /)
  integer :: of(2) = (/ 0, 10 /)
  double precision :: al(3) = (/ 0.0d0, 1.0d0, 2.5d0 /)
  double precision :: be(3) = (/ 0.0d0, 1.0d0, 0.5d0 /)
  integer :: ita, itb, im, ix, ik, ia, ib, io, m, n, k
  call rseed(20261001)
  do ita = 1, 3
    do itb = 1, 3
      do im = 1, 3
        do ix = 1, 3
          do ik = 1, 3
            do ia = 1, 3
              do ib = 1, 3
                do io = 1, 2
                  m = dm(im)
                  n = dm(ix)
                  k = dm(ik)
                  call fillb(a, nb)
                  call fillb(b, nb)
                  call fillb(c, nb)
                  call lbl(1, ita, itb, m*100 + n*10 + k, ia, ib, io)
                  call dgemm(tc(ita), tc(itb), m, n, k, al(ia), a(1+of(io)), &
                             ld, b(1+of(io)), ld, be(ib), c(1+of(io)), ld)
                  call dumpb(c, nb)
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end do
end subroutine tgemm

subroutine ttri(which)
  use rndmod
  implicit none
  integer, intent(in) :: which
  external dtrsm, dtrmm
  integer, parameter :: ld = 8, nb = 64
  double precision :: a(nb), b(nb)
  character :: sc(2) = (/ 'L', 'R' /)
  character :: uc(2) = (/ 'U', 'L' /)
  character :: tc(3) = (/ 'N', 'T', 'C' /)
  character :: dc(2) = (/ 'U', 'N' /)
  integer :: dm(4) = (/ 0, 1, 3, 5 /)
  integer :: of(2) = (/ 0, 10 /)
  double precision :: al(3) = (/ 0.0d0, 1.0d0, 2.5d0 /)
  integer :: is, iu, it, id, im, ix, ia, io, m, n
  call rseed(20261002 + which)
  do is = 1, 2
    do iu = 1, 2
      do it = 1, 3
        do id = 1, 2
          do im = 1, 4
            do ix = 1, 4
              do ia = 1, 3
                do io = 1, 2
                  m = dm(im)
                  n = dm(ix)
                  call fillb(a, nb)
                  call fillb(b, nb)
                  call lbl(2 + which, is*1000 + iu*100 + it*10 + id, &
                           m*10 + n, ia, io, 0, 0)
                  if (which == 0) then
                    call dtrsm(sc(is), uc(iu), tc(it), dc(id), m, n, al(ia), &
                               a(1+of(io)), ld, b(1+of(io)), ld)
                  else
                    call dtrmm(sc(is), uc(iu), tc(it), dc(id), m, n, al(ia), &
                               a(1+of(io)), ld, b(1+of(io)), ld)
                  end if
                  call dumpb(b, nb)
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end do
end subroutine ttri

subroutine tsyrk
  use rndmod
  implicit none
  external dsyrk
  integer, parameter :: ld = 8, nb = 64
  double precision :: a(nb), c(nb)
  character :: uc(2) = (/ 'U', 'L' /)
  character :: tc(3) = (/ 'N', 'T', 'C' /)
  integer :: dm(4) = (/ 0, 1, 3, 5 /)
  integer :: of(2) = (/ 0, 10 /)
  double precision :: al(3) = (/ 0.0d0, 1.0d0, 2.5d0 /)
  double precision :: be(3) = (/ 0.0d0, 1.0d0, 0.5d0 /)
  integer :: iu, it, ix, ik, ia, ib, io, n, k
  call rseed(20261004)
  do iu = 1, 2
    do it = 1, 3
      do ix = 1, 4
        do ik = 1, 4
          do ia = 1, 3
            do ib = 1, 3
              do io = 1, 2
                n = dm(ix)
                k = dm(ik)
                call fillb(a, nb)
                call fillb(c, nb)
                call lbl(4, iu*10 + it, n*10 + k, ia, ib, io, 0)
                call dsyrk(uc(iu), tc(it), n, k, al(ia), a(1+of(io)), ld, &
                           be(ib), c(1+of(io)), ld)
                call dumpb(c, nb)
              end do
            end do
          end do
        end do
      end do
    end do
  end do
end subroutine tsyrk

subroutine tgemv
  use rndmod
  implicit none
  external dgemv
  integer, parameter :: ld = 8, nb = 64
  double precision :: a(nb), x(nb), y(nb)
  character :: tc(3) = (/ 'N', 'T', 'C' /)
  integer :: dm(4) = (/ 0, 1, 3, 5 /)
  integer :: ic(2) = (/ 1, -2 /)
  integer :: of(2) = (/ 0, 10 /)
  double precision :: al(3) = (/ 0.0d0, 1.0d0, 2.5d0 /)
  double precision :: be(3) = (/ 0.0d0, 1.0d0, 0.5d0 /)
  integer :: it, im, ix, jx, jy, ia, ib, io, m, n
  call rseed(20261005)
  do it = 1, 3
    do im = 1, 4
      do ix = 1, 4
        do jx = 1, 2
          do jy = 1, 2
            do ia = 1, 3
              do ib = 1, 3
                do io = 1, 2
                  m = dm(im)
                  n = dm(ix)
                  call fillb(a, nb)
                  call fillb(x, nb)
                  call fillb(y, nb)
                  call lbl(5, it, m*10 + n, jx*10 + jy, ia, ib, io)
                  call dgemv(tc(it), m, n, al(ia), a(1+of(io)), ld, &
                             x(1+of(io)), ic(jx), be(ib), y(1+of(io)), ic(jy))
                  call dumpb(y, nb)
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end do
end subroutine tgemv

subroutine tger
  use rndmod
  implicit none
  external dger
  integer, parameter :: ld = 8, nb = 64
  double precision :: a(nb), x(nb), y(nb)
  integer :: dm(4) = (/ 0, 1, 3, 5 /)
  integer :: ic(2) = (/ 1, -2 /)
  integer :: of(2) = (/ 0, 10 /)
  double precision :: al(3) = (/ 0.0d0, 1.0d0, 2.5d0 /)
  integer :: im, ix, jx, jy, ia, io, m, n
  call rseed(20261006)
  do im = 1, 4
    do ix = 1, 4
      do jx = 1, 2
        do jy = 1, 2
          do ia = 1, 3
            do io = 1, 2
              m = dm(im)
              n = dm(ix)
              call fillb(a, nb)
              call fillb(x, nb)
              call fillb(y, nb)
              call lbl(6, m*10 + n, jx*10 + jy, ia, io, 0, 0)
              call dger(m, n, al(ia), x(1+of(io)), ic(jx), y(1+of(io)), &
                        ic(jy), a(1+of(io)), ld)
              call dumpb(a, nb)
            end do
          end do
        end do
      end do
    end do
  end do
end subroutine tger

subroutine tvec
  use rndmod
  implicit none
  external dswap, dscal
  integer, parameter :: nb = 64
  integer :: idamax
  double precision :: ddot, dasum, dnrm2
  external idamax, ddot, dasum, dnrm2
  double precision :: x(nb), y(nb)
  integer :: dm(4) = (/ 0, 1, 3, 5 /)
  integer :: dn(6) = (/ 0, 1, 3, 5, 7, 13 /)
  integer :: ic(4) = (/ 1, 2, -1, -2 /)
  integer :: of(2) = (/ 0, 10 /)
  double precision :: al(3) = (/ 0.0d0, 1.0d0, 2.5d0 /)
  double precision :: sc(3) = (/ 1.0d0, 2.0d0**(-600), 2.0d0**600 /)
  integer :: i, is, ix, jx, jy, io, ia, n
  call rseed(20261007)
  do ix = 1, 4
    do jx = 1, 4
      do jy = 1, 4
        do io = 1, 2
          n = dm(ix)
          call fillb(x, nb)
          call fillb(y, nb)
          call lbl(7, n, jx*10 + jy, io, 0, 0, 0)
          call dswap(n, x(1+of(io)), ic(jx), y(1+of(io)), ic(jy))
          call dumpb(x, nb)
          call dumpb(y, nb)
        end do
      end do
    end do
  end do
  do ix = 1, 4
    do jx = 1, 4
      do ia = 1, 3
        do io = 1, 2
          n = dm(ix)
          call fillb(x, nb)
          call lbl(8, n, jx, ia, io, 0, 0)
          call dscal(n, al(ia), x(1+of(io)), ic(jx))
          call dumpb(x, nb)
        end do
      end do
    end do
  end do
  do ix = 1, 4
    do jx = 1, 4
      do io = 1, 2
        n = dm(ix)
        call fillb(x, nb)
        call lbl(9, n, jx, io, 0, 0, 0)
        call pin(idamax(n, x(1+of(io)), ic(jx)))
      end do
    end do
  end do
  do ix = 1, 6
    do jx = 1, 4
      do jy = 1, 4
        do io = 1, 2
          n = dn(ix)
          call fillb(x, nb)
          call fillb(y, nb)
          call lbl(50, n, jx*10 + jy, io, 0, 0, 0)
          call pr(ddot(n, x(1+of(io)), ic(jx), y(1+of(io)), ic(jy)))
        end do
      end do
    end do
  end do
  do ix = 1, 6
    do jx = 1, 4
      do io = 1, 2
        n = dn(ix)
        call fillb(x, nb)
        call lbl(51, n, jx, io, 0, 0, 0)
        call pr(dasum(n, x(1+of(io)), ic(jx)))
      end do
    end do
  end do
  do ix = 1, 6
    do jx = 1, 4
      do is = 1, 3
        do io = 1, 2
          n = dn(ix)
          call fillb(x, nb)
          do i = 2, nb, 2
            x(i) = x(i)*sc(is)
          end do
          call lbl(52, n, jx*10 + is, io, 0, 0, 0)
          call pr(dnrm2(n, x(1+of(io)), ic(jx)))
        end do
      end do
    end do
  end do
end subroutine tvec

subroutine tlu(which)
  use rndmod
  implicit none
  integer, intent(in) :: which
  external dgetf2, dgetrf2, dgetrf
  integer, parameter :: ld = 16, nb = 256
  double precision :: a(nb)
  integer :: ipiv(ld), info, i, j, ic, io, iz, m, n
  integer :: ms(11) = (/ 1, 1, 3, 2, 3, 5, 5, 3, 8, 13, 7 /)
  integer :: ns(11) = (/ 1, 3, 1, 3, 2, 5, 3, 5, 8, 7, 13 /)
  integer :: of(2) = (/ 0, 17 /)
  call rseed(20261010 + which)
  do ic = 1, 11
    do io = 1, 2
      do iz = 0, 1
        m = ms(ic)
        n = ns(ic)
        call fillb(a, nb)
        if (iz == 1) then
          j = (n + 1)/2
          do i = 1, m
            a(1 + of(io) + i - 1 + (j-1)*ld) = 0.0d0
          end do
        end if
        do i = 1, ld
          ipiv(i) = -7
        end do
        info = -99
        call lbl(10 + which, m, n, io, iz, 0, 0)
        if (which == 0) then
          call dgetf2(m, n, a(1+of(io)), ld, ipiv, info)
        else if (which == 1) then
          call dgetrf2(m, n, a(1+of(io)), ld, ipiv, info)
        else
          call dgetrf(m, n, a(1+of(io)), ld, ipiv, info)
        end if
        call pin(info)
        do i = 1, ld
          call pin(ipiv(i))
        end do
        call dumpb(a, nb)
      end do
    end do
  end do
end subroutine tlu

subroutine tchol(which)
  use rndmod
  implicit none
  integer, intent(in) :: which
  external dpotf2, dpotrf2, dpotrf
  integer, parameter :: ld = 16, nb = 256
  double precision :: a(nb)
  character :: uc(2) = (/ 'U', 'L' /)
  integer :: ns(6) = (/ 1, 2, 3, 5, 8, 13 /)
  integer :: of(2) = (/ 0, 17 /)
  integer :: ic, iu, io, iz, info, n
  call rseed(20261020 + which)
  do ic = 1, 6
    do iu = 1, 2
      do io = 1, 2
        do iz = 0, 1
          n = ns(ic)
          call fillspd(a, nb, of(io), ld, n)
          if (iz == 1) a(1 + of(io) + (n-1) + (n-1)*ld) = -1.0d0
          info = -99
          call lbl(20 + which, n, iu, io, iz, 0, 0)
          if (which == 0) then
            call dpotf2(uc(iu), n, a(1+of(io)), ld, info)
          else if (which == 1) then
            call dpotrf2(uc(iu), n, a(1+of(io)), ld, info)
          else
            call dpotrf(uc(iu), n, a(1+of(io)), ld, info)
          end if
          call pin(info)
          call dumpb(a, nb)
        end do
      end do
    end do
  end do
end subroutine tchol

subroutine tsolve
  use rndmod
  implicit none
  external dgesv, dgetrf, dgetrs, dposv, dpotrs, dpotrf
  integer, parameter :: ld = 16, nb = 256
  double precision :: a(nb), b(nb)
  character :: uc(2) = (/ 'U', 'L' /)
  character :: tc(3) = (/ 'N', 'T', 'C' /)
  integer :: ns(5) = (/ 1, 2, 3, 5, 8 /)
  integer :: rs(3) = (/ 1, 2, 3 /)
  integer :: of(2) = (/ 0, 17 /)
  integer :: ipiv(ld)
  integer :: ic, ir, io, it, iu, info, i, n, nrhs
  call rseed(20261030)
  do ic = 1, 5
    do ir = 1, 3
      do io = 1, 2
        n = ns(ic)
        nrhs = rs(ir)
        call fillb(a, nb)
        call fillb(b, nb)
        do i = 1, ld
          ipiv(i) = -7
        end do
        info = -99
        call lbl(30, n, nrhs, io, 0, 0, 0)
        call dgesv(n, nrhs, a(1+of(io)), ld, ipiv, b(1+of(io)), ld, info)
        call pin(info)
        do i = 1, ld
          call pin(ipiv(i))
        end do
        call dumpb(a, nb)
        call dumpb(b, nb)
      end do
    end do
  end do
  do ic = 1, 5
    do ir = 1, 3
      do it = 1, 3
        do io = 1, 2
          n = ns(ic)
          nrhs = rs(ir)
          call fillb(a, nb)
          call fillb(b, nb)
          do i = 1, ld
            ipiv(i) = -7
          end do
          info = -99
          call dgetrf(n, n, a(1+of(io)), ld, ipiv, info)
          call lbl(31, n, nrhs, it, io, info, 0)
          info = -99
          call dgetrs(tc(it), n, nrhs, a(1+of(io)), ld, ipiv, b(1+of(io)), &
                      ld, info)
          call pin(info)
          call dumpb(b, nb)
        end do
      end do
    end do
  end do
  do ic = 1, 5
    do ir = 1, 3
      do iu = 1, 2
        do io = 1, 2
          n = ns(ic)
          nrhs = rs(ir)
          call fillspd(a, nb, of(io), ld, n)
          call fillb(b, nb)
          info = -99
          call lbl(32, n, nrhs, iu, io, 0, 0)
          call dposv(uc(iu), n, nrhs, a(1+of(io)), ld, b(1+of(io)), ld, info)
          call pin(info)
          call dumpb(a, nb)
          call dumpb(b, nb)
        end do
      end do
    end do
  end do
  do ic = 1, 5
    do ir = 1, 3
      do iu = 1, 2
        do io = 1, 2
          n = ns(ic)
          nrhs = rs(ir)
          call fillspd(a, nb, of(io), ld, n)
          call fillb(b, nb)
          info = -99
          call dpotrf(uc(iu), n, a(1+of(io)), ld, info)
          call lbl(33, n, nrhs, iu, io, info, 0)
          info = -99
          call dpotrs(uc(iu), n, nrhs, a(1+of(io)), ld, b(1+of(io)), ld, info)
          call pin(info)
          call dumpb(b, nb)
        end do
      end do
    end do
  end do
end subroutine tsolve

subroutine tbig
  use rndmod
  implicit none
  external dgetrf, dpotrf, dgesv, dposv
  integer, parameter :: ld = 104, nb = 10816
  double precision, allocatable :: a(:), b(:)
  integer, allocatable :: ipiv(:)
  character :: uc(2) = (/ 'U', 'L' /)
  integer :: ns(5) = (/ 63, 64, 65, 70, 100 /)
  integer :: ic, iu, info, i, n
  allocate(a(nb), b(nb), ipiv(ld))
  call rseed(20261040)
  do ic = 1, 5
    n = ns(ic)
    call fillb(a, nb)
    do i = 1, ld
      ipiv(i) = -7
    end do
    info = -99
    call lbl(40, n, 0, 0, 0, 0, 0)
    call dgetrf(n, n, a, ld, ipiv, info)
    call pin(info)
    do i = 1, ld
      call pin(ipiv(i))
    end do
    call dumpb(a, nb)
  end do
  do ic = 1, 5
    do iu = 1, 2
      n = ns(ic)
      call fillspd(a, nb, 0, ld, n)
      info = -99
      call lbl(41, n, iu, 0, 0, 0, 0)
      call dpotrf(uc(iu), n, a, ld, info)
      call pin(info)
      call dumpb(a, nb)
    end do
  end do
  do ic = 1, 5
    n = ns(ic)
    call fillb(a, nb)
    call fillb(b, nb)
    do i = 1, ld
      ipiv(i) = -7
    end do
    info = -99
    call lbl(42, n, 0, 0, 0, 0, 0)
    call dgesv(n, 3, a, ld, ipiv, b, ld, info)
    call pin(info)
    call dumpb(b, nb)
  end do
  do ic = 1, 5
    do iu = 1, 2
      n = ns(ic)
      call fillspd(a, nb, 0, ld, n)
      call fillb(b, nb)
      info = -99
      call lbl(43, n, iu, 0, 0, 0, 0)
      call dposv(uc(iu), n, 3, a, ld, b, ld, info)
      call pin(info)
      call dumpb(b, nb)
    end do
  end do
  deallocate(a, b, ipiv)
end subroutine tbig

program fdiff
  implicit none
  external tgemm, ttri, tsyrk, tgemv, tger, tvec, tlu, tchol, tsolve, tbig
  call tgemm
  call ttri(0)
  call ttri(1)
  call tsyrk
  call tgemv
  call tger
  call tvec
  call tlu(0)
  call tlu(1)
  call tlu(2)
  call tchol(0)
  call tchol(1)
  call tchol(2)
  call tsolve
  call tbig
end program fdiff
