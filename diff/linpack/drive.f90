program drive
  use locs
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  call gb_tri(); call gb_penta(); call gb_diag(); call gb_lap(); call gb_cd()
  call gb_slice(); call gb_inv(); call pb_tri(); call pb_diag(); call pb_b2()
  call pb_fem(); call pb_mass(); call pb_string(); call pb_sym()
contains
  subroutine pv(tag, b, n)
    character(*) :: tag
    integer :: n, i
    real(dp) :: b(n)
    write(*,'(A,A)', advance='no') tag, ' '
    do i = 1, n
      write(*,'(ES24.16)', advance='no') b(i)
    end do
    write(*,*)
  end subroutine
  subroutine gb_tri()
    integer, parameter :: n=4, ml=1, mu=1, lda=2*ml+mu+1
    real(dp) :: abd(lda,n), b(n)
    integer :: ipvt(n), info, i
    abd = 0; abd(ml+mu+1,:) = 2
    do i = 1, n-1
      abd(ml+mu,i+1) = -1; abd(ml+mu+2,i) = -1
    end do
    b = [1.0_dp,0.0_dp,0.0_dp,1.0_dp]
    call dgbfa_local(abd,lda,n,ml,mu,ipvt,info); call dgbsl_local(abd,lda,n,ml,mu,ipvt,b,0)
    call pv('gb_tri', b, n)
  end subroutine
  subroutine gb_penta()
    integer, parameter :: n=5, ml=2, mu=2, lda=2*ml+mu+1
    real(dp) :: abd(lda,n), b(n)
    integer :: ipvt(n), info
    abd = 0
    abd(ml+mu+1,:) = 6; abd(ml+mu,2:n) = -2; abd(ml+mu-1,3:n) = -1
    abd(ml+mu+2,1:n-1) = -2; abd(ml+mu+3,1:n-2) = -1
    b = [3.0_dp,1.0_dp,0.0_dp,1.0_dp,3.0_dp]
    call dgbfa_local(abd,lda,n,ml,mu,ipvt,info); call dgbsl_local(abd,lda,n,ml,mu,ipvt,b,0)
    call pv('gb_penta', b, n)
  end subroutine
  subroutine gb_diag()
    integer, parameter :: n=3, ml=0, mu=0, lda=1
    real(dp) :: abd(lda,n), b(n)
    integer :: ipvt(n), info
    abd(1,:) = [2.0_dp,4.0_dp,8.0_dp]
    b = [4.0_dp,12.0_dp,32.0_dp]
    call dgbfa_local(abd,lda,n,ml,mu,ipvt,info); call dgbsl_local(abd,lda,n,ml,mu,ipvt,b,0)
    call pv('gb_diag', b, n)
  end subroutine
  subroutine gb_lap()
    integer, parameter :: n=4, ml=1, mu=1, lda=2*ml+mu+1
    real(dp) :: abd(lda,n), b(n), h
    integer :: ipvt(n), info
    h = 0.2_dp
    abd = 0; abd(ml+mu+1,:) = 2/h**2; abd(ml+mu,2:n) = -1/h**2; abd(ml+mu+2,1:n-1) = -1/h**2
    b = 2
    call dgbfa_local(abd,lda,n,ml,mu,ipvt,info); call dgbsl_local(abd,lda,n,ml,mu,ipvt,b,0)
    call pv('gb_lap', b, n)
  end subroutine
  subroutine gb_cd()
    integer, parameter :: n=4, ml=1, mu=1, lda=2*ml+mu+1
    real(dp) :: abd(lda,n), b(n), h, eps, diag, sub, sup
    integer :: ipvt(n), info
    h = 0.2_dp; eps = 0.1_dp
    sub = -eps/h**2 - 1/h; diag = 2*eps/h**2 + 1/h; sup = -eps/h**2
    abd = 0; abd(ml+mu+1,:) = diag; abd(ml+mu,2:n) = sup; abd(ml+mu+2,1:n-1) = sub
    b = 0; b(n) = -sup
    call dgbfa_local(abd,lda,n,ml,mu,ipvt,info); call dgbsl_local(abd,lda,n,ml,mu,ipvt,b,0)
    call pv('gb_cd', b, n)
  end subroutine
  subroutine gb_slice()
    integer, parameter :: n=3, ml=1, mu=1, lda=2*ml+mu+1
    real(dp) :: abd(lda,n), b(n)
    integer :: ipvt(n), info
    abd = 0; abd(ml+mu+1,:) = 4; abd(ml+mu,2:n) = -1; abd(ml+mu+2,1:n-1) = -1
    b = [3.0_dp,2.0_dp,3.0_dp]
    call dgbfa_local(abd,lda,n,ml,mu,ipvt,info); call dgbsl_local(abd,lda,n,ml,mu,ipvt,b,0)
    call pv('gb_slice', b, n)
  end subroutine
  subroutine gb_inv()
    integer, parameter :: n=3, ml=1, mu=1, lda=2*ml+mu+1
    real(dp) :: abd(lda,n), b(n)
    integer :: ipvt(n), info
    abd = 0; abd(ml+mu+1,:) = 4; abd(ml+mu,2:n) = -1; abd(ml+mu+2,1:n-1) = -1
    b = [1.0_dp,0.0_dp,0.0_dp]
    call dgbfa_local(abd,lda,n,ml,mu,ipvt,info); call dgbsl_local(abd,lda,n,ml,mu,ipvt,b,0)
    call pv('gb_inv', b, n)
  end subroutine
  subroutine pb_tri()
    integer, parameter :: n=4, m=1, lda=m+1
    real(dp) :: abd(lda,n), b(n)
    integer :: info
    abd = 0; abd(m+1,:) = 2; abd(m,2:n) = -1
    b = [1.0_dp,0.0_dp,0.0_dp,1.0_dp]
    call dpbfa_local(abd,lda,n,m,info); call dpbsl_local(abd,lda,n,m,b)
    call pv('pb_tri', b, n)
  end subroutine
  subroutine pb_diag()
    integer, parameter :: n=3, m=0, lda=1
    real(dp) :: abd(lda,n), b(n)
    integer :: info
    abd(1,:) = [4.0_dp,9.0_dp,16.0_dp]
    b = [8.0_dp,27.0_dp,64.0_dp]
    call dpbfa_local(abd,lda,n,m,info); call dpbsl_local(abd,lda,n,m,b)
    call pv('pb_diag', b, n)
  end subroutine
  subroutine pb_b2()
    integer, parameter :: n=4, m=2, lda=m+1
    real(dp) :: abd(lda,n), b(n)
    integer :: info
    abd = 0; abd(m+1,:) = 10; abd(m,2:n) = -3; abd(m-1,3:n) = 1
    b = [8.0_dp,5.0_dp,5.0_dp,8.0_dp]
    call dpbfa_local(abd,lda,n,m,info); call dpbsl_local(abd,lda,n,m,b)
    call pv('pb_b2', b, n)
  end subroutine
  subroutine pb_fem()
    integer, parameter :: n=4, m=1, lda=m+1
    real(dp) :: abd(lda,n), b(n)
    integer :: info
    abd = 0; abd(m+1,:) = 2; abd(m,2:n) = -1
    b = 1
    call dpbfa_local(abd,lda,n,m,info); call dpbsl_local(abd,lda,n,m,b)
    call pv('pb_fem', b, n)
  end subroutine
  subroutine pb_mass()
    integer, parameter :: n=3, m=1, lda=m+1
    real(dp) :: abd(lda,n), b(n)
    integer :: info
    abd = 0; abd(m+1,:) = 4.0_dp/6.0_dp; abd(m,2:n) = 1.0_dp/6.0_dp
    b = 1
    call dpbfa_local(abd,lda,n,m,info); call dpbsl_local(abd,lda,n,m,b)
    call pv('pb_mass', b, n)
  end subroutine
  subroutine pb_string()
    integer, parameter :: n=5, m=1, lda=m+1
    real(dp) :: abd(lda,n), b(n)
    integer :: info, i
    real(dp) :: pi
    pi = 4*atan(1.0_dp)
    abd = 0; abd(m+1,:) = 2; abd(m,2:n) = -1
    do i = 1, n
      b(i) = sin(pi*real(i,dp)/real(n+1,dp))
    end do
    call dpbfa_local(abd,lda,n,m,info); call dpbsl_local(abd,lda,n,m,b)
    call pv('pb_string', b, n)
  end subroutine
  subroutine pb_sym()
    integer, parameter :: n=3, m=1, lda=m+1
    real(dp) :: abd(lda,n), b(n), c(n)
    integer :: info
    abd = 0; abd(m+1,:) = 4; abd(m,2:n) = -1
    call dpbfa_local(abd,lda,n,m,info)
    b = [1.0_dp,0.0_dp,0.0_dp]; call dpbsl_local(abd,lda,n,m,b); call pv('pb_sym1', b, n)
    c = [0.0_dp,1.0_dp,0.0_dp]; call dpbsl_local(abd,lda,n,m,c); call pv('pb_sym2', c, n)
  end subroutine
end program
