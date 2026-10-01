module locs
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
contains
  subroutine dgbfa_local(abd, lda, n, ml, mu, ipvt, info)
    integer, intent(in) :: lda, n, ml, mu
    real(dp), intent(inout) :: abd(lda, *)
    integer, intent(out) :: ipvt(*), info
    integer :: i, j, k, l, m, lm, mm, ju, jz
    real(dp) :: t

    m = ml + mu + 1
    info = 0

    ! Zero initial fill-in area
    jz = min(mu+1, n)
    if (jz > 1) then
      do j = 1, jz - 1
        do i = ml + 2 - j, ml
          abd(i, j) = 0.0_dp
        end do
      end do
    end if

    jz = jz + 1
    ju = 0

    do k = 1, n - 1
      ! Zero next fill-in column
      if (jz <= n) then
        do i = 1, ml
          abd(i, jz) = 0.0_dp
        end do
        jz = jz + 1
      end if

      ! Find pivot
      lm = min(ml, n - k)
      l = m
      do j = m + 1, m + lm
        if (abs(abd(j, k)) > abs(abd(l, k))) l = j
      end do
      ipvt(k) = l + k - m

      if (abd(l, k) == 0.0_dp) then
        info = k
        cycle
      end if

      ! Swap
      if (l /= m) then
        t = abd(l, k)
        abd(l, k) = abd(m, k)
        abd(m, k) = t
      end if

      ! Compute multipliers
      t = -1.0_dp / abd(m, k)
      do i = m + 1, m + lm
        abd(i, k) = abd(i, k) * t
      end do

      ! Row elimination
      ju = max(ju, min(mu + ipvt(k), n))
      mm = m
      do j = k + 1, ju
        l = l - 1
        mm = mm - 1
        t = abd(l, j)
        if (l /= mm) then
          abd(l, j) = abd(mm, j)
          abd(mm, j) = t
        end if
        do i = 1, lm
          abd(mm + i, j) = abd(mm + i, j) + t * abd(m + i, k)
        end do
      end do
    end do

    ipvt(n) = n
    if (abd(m, n) == 0.0_dp) info = n

  end subroutine dgbfa_local

  ! DGBSL: Solve using banded LU factorization
  subroutine dgbsl_local(abd, lda, n, ml, mu, ipvt, b, job)
    integer, intent(in) :: lda, n, ml, mu, job
    real(dp), intent(in) :: abd(lda, *)
    integer, intent(in) :: ipvt(*)
    real(dp), intent(inout) :: b(*)
    integer :: k, l, m, lm
    real(dp) :: t

    m = mu + ml + 1

    if (job == 0) then
      ! Solve A*x = b

      ! Forward substitution
      if (ml > 0) then
        do k = 1, n - 1
          lm = min(ml, n - k)
          l = ipvt(k)
          t = b(l)
          if (l /= k) then
            b(l) = b(k)
            b(k) = t
          end if
          do l = 1, lm
            b(k + l) = b(k + l) + t * abd(m + l, k)
          end do
        end do
      end if

      ! Back substitution
      do k = n, 1, -1
        b(k) = b(k) / abd(m, k)
        t = -b(k)
        lm = min(k - 1, m - 1)
        do l = 1, lm
          b(k - l) = b(k - l) + t * abd(m - l, k)
        end do
      end do
    end if

  end subroutine dgbsl_local
  subroutine dpbfa_local(abd, lda, n, m, info)
    integer, intent(in) :: lda, n, m
    real(dp), intent(inout) :: abd(lda, *)
    integer, intent(out) :: info
    integer :: j, k, ik, jk, mu
    real(dp) :: s, t

    do j = 1, n
      s = 0.0_dp
      ik = m + 1
      jk = max(j - m, 1)
      mu = max(m + 2 - j, 1)

      if (m >= mu) then
        do k = mu, m
          t = abd(k, j) - sum(abd(ik:m, jk) * abd(mu:k-1, j))
          t = t / abd(m + 1, jk)
          abd(k, j) = t
          s = s + t * t
          ik = ik - 1
          jk = jk + 1
        end do
      end if

      s = abd(m + 1, j) - s
      if (s <= 0.0_dp) then
        info = j
        return
      end if
      abd(m + 1, j) = sqrt(s)
    end do
    info = 0

  end subroutine dpbfa_local

  ! DPBSL: Solve using banded Cholesky factorization
  subroutine dpbsl_local(abd, lda, n, m, b)
    integer, intent(in) :: lda, n, m
    real(dp), intent(in) :: abd(lda, *)
    real(dp), intent(inout) :: b(*)
    integer :: k, l, lm, la, lb
    real(dp) :: t

    ! Solve R'*y = b (forward substitution)
    do k = 1, n
      lm = min(k - 1, m)
      la = m + 1 - lm
      lb = k - lm
      t = sum(abd(la:m, k) * b(lb:k-1))
      b(k) = (b(k) - t) / abd(m + 1, k)
    end do

    ! Solve R*x = y (back substitution)
    ! R(k, k+l) is stored at ABD(M+1-l, k+l)
    do k = n, 1, -1
      lm = min(m, n - k)
      t = 0.0_dp
      do l = 1, lm
        t = t + abd(m + 1 - l, k + l) * b(k + l)
      end do
      b(k) = (b(k) - t) / abd(m + 1, k)
    end do

  end subroutine dpbsl_local
end module locs
