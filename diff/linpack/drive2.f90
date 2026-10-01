program drive2
  use linpack
  use service, only : DP
  implicit none
  real(DP) :: A(3,3), A2(2,2), b(3), b2(2), z(3), rcond
  integer :: ipvt(3), ipvt2(2), info, i, j
  A2(1,:) = [2.0_dp,1.0_dp]; A2(2,:) = [1.0_dp,3.0_dp]
  b2 = [5.0_dp,7.0_dp]
  call DGEFA(A2,2,2,ipvt2,info); call DGESL(A2,2,2,ipvt2,b2,0); call pv('ge_2x2',b2,2,info)
  do i = 1, 3
    do j = 1, 3
      A(i,j) = 1.0_dp/real(i+j-1,dp)
    end do
  end do
  b = [A(1,1)+A(1,2)+A(1,3), A(2,1)+A(2,2)+A(2,3), A(3,1)+A(3,2)+A(3,3)]
  call DGEFA(A,3,3,ipvt,info); call DGESL(A,3,3,ipvt,b,0); call pv('ge_hilb',b,3,info)
  A(1,:) = [1.0_dp,2.0_dp,3.0_dp]; A(2,:) = [4.0_dp,5.0_dp,6.0_dp]; A(3,:) = [7.0_dp,8.0_dp,0.0_dp]
  b = [14.0_dp,32.0_dp,23.0_dp]
  call DGEFA(A,3,3,ipvt,info); call DGESL(A,3,3,ipvt,b,0); call pv('ge_doc',b,3,info)
  A(1,:) = [1.0_dp,2.0_dp,3.0_dp]; A(2,:) = [4.0_dp,5.0_dp,6.0_dp]; A(3,:) = [7.0_dp,8.0_dp,0.0_dp]
  b = [36.0_dp,42.0_dp,18.0_dp]
  call DGEFA(A,3,3,ipvt,info); call DGESL(A,3,3,ipvt,b,1); call pv('ge_trans',b,3,info)
  A(1,:) = [1.0_dp,1.0_dp,1.0_dp]; A(2,:) = [1.0_dp,2.0_dp,3.0_dp]; A(3,:) = [1.0_dp,3.0_dp,6.0_dp]
  b = [3.0_dp,6.0_dp,10.0_dp]
  call DGEFA(A,3,3,ipvt,info); call DGESL(A,3,3,ipvt,b,0); call pv('ge_pascal',b,3,info)
  A(1,:) = [1.0_dp,1.0_dp,1.0_dp]; A(2,:) = [1.0_dp,2.0_dp,4.0_dp]; A(3,:) = [1.0_dp,3.0_dp,9.0_dp]
  b = [6.0_dp,17.0_dp,34.0_dp]
  call DGEFA(A,3,3,ipvt,info); call DGESL(A,3,3,ipvt,b,0); call pv('ge_vander',b,3,info)
  A(1,:) = [0.0001_dp,1.0_dp,0.0_dp]; A(2,:) = [1.0_dp,1.0_dp,0.0_dp]; A(3,:) = [0.0_dp,0.0_dp,1.0_dp]
  b = [1.0001_dp,2.0_dp,1.0_dp]
  call DGEFA(A,3,3,ipvt,info); call DGESL(A,3,3,ipvt,b,0); call pv('ge_fm',b,3,info)
  A(1,:) = [1.0e-20_dp,1.0_dp,0.0_dp]; A(2,:) = [1.0_dp,1.0_dp,0.0_dp]; A(3,:) = [0.0_dp,0.0_dp,1.0_dp]
  b = [1.0_dp+1.0e-20_dp,2.0_dp,1.0_dp]
  call DGEFA(A,3,3,ipvt,info); call DGESL(A,3,3,ipvt,b,0); call pv('ge_wilk',b,3,info)
  A(1,:) = [1.0_dp,2.0_dp,3.0_dp]; A(2,:) = [0.0_dp,4.0_dp,5.0_dp]; A(3,:) = [0.0_dp,0.0_dp,6.0_dp]
  b = [14.0_dp,23.0_dp,18.0_dp]
  call DGEFA(A,3,3,ipvt,info); call DGESL(A,3,3,ipvt,b,0); call pv('ge_upper',b,3,info)
  A2(1,:) = [4.0_dp,2.0_dp]; A2(2,:) = [2.0_dp,5.0_dp]
  b2 = [6.0_dp,7.0_dp]
  call DPOFA(A2,2,2,info); call DPOSL(A2,2,2,b2); call pv('po_2x2',b2,2,info)
  A(1,:) = [1.0_dp,0.5_dp,1.0_dp/3.0_dp]; A(2,:) = [0.5_dp,1.0_dp,2.0_dp/3.0_dp]
  A(3,:) = [1.0_dp/3.0_dp,2.0_dp/3.0_dp,1.0_dp]
  b = [11.0_dp/6.0_dp,13.0_dp/6.0_dp,2.0_dp]
  call DPOFA(A,3,3,info); call DPOSL(A,3,3,b); call pv('po_lehmer',b,3,info)
  A(1,:) = [2.0_dp,-1.0_dp,0.0_dp]; A(2,:) = [-1.0_dp,2.0_dp,-1.0_dp]; A(3,:) = [0.0_dp,-1.0_dp,2.0_dp]
  b = [1.0_dp,0.0_dp,1.0_dp]
  call DPOFA(A,3,3,info); call DPOSL(A,3,3,b); call pv('po_tri',b,3,info)
  A(1,:) = [2.0_dp,1.0_dp,1.0_dp]; A(2,:) = [1.0_dp,2.0_dp,1.0_dp]; A(3,:) = [1.0_dp,1.0_dp,2.0_dp]
  b = [4.0_dp,4.0_dp,4.0_dp]
  call DPOFA(A,3,3,info); call DPOSL(A,3,3,b); call pv('po_ones',b,3,info)
  A2(1,:) = [1.0_dp,2.0_dp]; A2(2,:) = [2.0_dp,1.0_dp]
  call DPOFA(A2,2,2,info); print '(A,I3)', 'po_indef info', info
  A(1,:) = [1.0_dp,2.0_dp,3.0_dp]; A(2,:) = [2.0_dp,4.0_dp,6.0_dp]; A(3,:) = [1.0_dp,1.0_dp,1.0_dp]
  call DGEFA(A,3,3,ipvt,info); print '(A,I3)', 'ge_sing info', info
  do i = 1, 3
    do j = 1, 3
      A(i,j) = 1.0_dp/real(i+j-1,dp)
    end do
  end do
  call DGECO(A,3,3,ipvt,rcond,z)
  print '(A,ES24.16)', 'geco_hilb', rcond
  call pv('geco_hilb_z', z, 3, 0)
contains
  subroutine pv(tag, b, n, info)
    character(*) :: tag
    integer :: n, i, info
    real(DP) :: b(n)
    write(*,'(A,A,I2)', advance='no') tag, ' info', info
    do i = 1, n
      write(*,'(ES24.16)', advance='no') b(i)
    end do
    write(*,*)
  end subroutine
end program
