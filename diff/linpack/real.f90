module service
  use, intrinsic :: iso_fortran_env, only: real32, real64
  implicit none
  integer, parameter :: SP = real32, DP = real64
end module service
module blas
  use service, only : DP
  implicit none
contains
  pure subroutine DAXPY(n, da, dx, incx, dy, incy)
    integer, intent(in) :: n, incx, incy
    real(DP), intent(in) :: da, dx(*)
    real(DP), intent(inout) :: dy(*)
    integer :: i, ix, iy
    if (n <= 0 .or. da == 0.0_DP) return
    ix = 1; iy = 1
    if (incx < 0) ix = (1-n)*incx + 1
    if (incy < 0) iy = (1-n)*incy + 1
    do i = 1, n
      dy(iy) = dy(iy) + da*dx(ix)
      ix = ix + incx; iy = iy + incy
    end do
  end subroutine DAXPY
end module blas
module linpack
  use service, only : SP, DP
  implicit none
contains
  include "dgefa.inc"
  include "dgesl.inc"
  include "dgeco.inc"
  include "dpofa.inc"
  include "dposl.inc"
  include "dgbfa.inc"
  include "dgbsl.inc"
end module linpack
