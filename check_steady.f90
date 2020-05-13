subroutine check_steady()
use variables
implicit none

    
    if(istep /= isto) then

        VelocityDifference = 0
        !error =  abs(  ( u(nx/2,ny/2,nz/2) - last_velocity ) / last_velocity  ) 

        do k=1,nz; do j=1,ny; do i=1,nx
            if(abs( wc(i,j,k)- last_velocity(i,j,k) ) > VelocityDifference ) then
                VelocityDifference =  abs( wc(i,j,k)- last_velocity(i,j,k) )
            end if
        end do; end do; end do

        if(myid==master)then
            open (16,file='velocity_difference.dat',position='append')
            write(16,'(E12.5,a,E12.5)') time,'      ',VelocityDifference
            close(16)
        end if

    end if

    do k=1,nz; do j=1,ny; do i=1,nx
        last_velocity(i,j,k) = wc(i,j,k)
    end do; end do; end do

end subroutine check_steady
