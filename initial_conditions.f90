subroutine initial_conditions()
use variables
implicit none

!---------------------------------------------------!
!         Initial conditions calculation            !
!---------------------------------------------------!

!$OMP PARALLEL DO PRIVATE(j,i)  
do k=-1,nz+2; do j=-1,ny+2; do i=-1,nx+2

   u(i,j,k) = 0.0
   v(i,j,k) = 0.0
   w(i,j,k) = 0.0
   u1(i,j,k) = 0.0
   v1(i,j,k) = 0.0
   w1(i,j,k) = 0.0
   u2(i,j,k) = 0.0
   v2(i,j,k) = 0.0
   w2(i,j,k) = 0.0
   p(i,j,k) = 0.0
   u_star(i,j,k) = 0.0
   v_star(i,j,k) = 0.0
   w_star(i,j,k) = 0.0
   ETA(i,j,k) = 0.0
   F_tavex(i,j,k) = 0.0
   F_tavey(i,j,k) = 0.0

end do; end do; end do
!$OMP END PARALLEL DO

end subroutine initial_conditions



