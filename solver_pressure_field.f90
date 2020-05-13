subroutine solver_pressure_field()
use variables
implicit none

!-------------------------------------------------!
!              Update pressure field              !
!-------------------------------------------------!

ik = 1 + myid * gcount(myid) * (nx) * (ny)
do k=istart,iend; do j=1,ny; do i=1,nx 

    pre(i,j,k)=x1(ik)
    p(i,j,k) = pre(i,j,k)
    ik=ik+1

enddo; enddo; enddo




end subroutine solver_pressure_field
