subroutine reading_data()
    use variables
    implicit none

    
    if(Gridder=='non-uniform')then
        !Small interval
        dxSml = lxSml/nxSml
        dySml = lySml/nySml
        dzSml = lzSml/nzSml

        !Middle intervel
        dyMid = (lyMid-lySml)/(nyMid-nySml)
        dzMid = (lzMid-lzSml)/(nzMid-nzSml)

        !Large interval
        dx = ( lx-lxSml ) / ( nx - nxSml ) 
        dy = ( ly-lyMid ) / ( ny - nyMid )
        dz = ( lz-lzMid ) / ( nz - nzMid )

        

        dx = lx / nx 
    

    else if(Gridder=='uniform')then

        dx = lx / nx 
        dy = ly / ny
        dz = lz / nz


    end if




    nu = 1./Re

end subroutine reading_data






subroutine reading_variables()
use variables
implicit none

open (18,file=inputfile,form='unformatted')
read(18) inblocks
read(18) inx, iny, inz
read(18) temp, temp, temp, temp
read(18) ( ( ( ( Qout(i,j,k,h), i = 1, nx), j = 1, ny), k = 1, nz), h = 1, 5 )
close(18)


do k=1,nz-1; do j=1,ny-1; do i=1,nx-1

u(i,j,k) = 0.5*( Qout(i,j,k,2)+Qout(i+1,j,k,2) )
v(i,j,k) = 0.5*( Qout(i,j,k,3)+Qout(i,j+1,k,3) )
w(i,j,k) = 0.5*( Qout(i,j,k,4)+Qout(i,j,k+1,4) )

enddo; enddo; enddo

do k=1,nz; do j=1,ny; do i=1,nx

   p(i,j,k) = Qout(i,j,k,1)
   pre(i,j,k) = p(i,j,k)

end do; end do; end do


do k=1,nz; do j=1,ny; do i=1,nx

uc(i,j,k) = 0.5*(u(i,j,k)+u(i-1,j,k))
vc(i,j,k) = 0.5*(v(i,j,k)+v(i,j-1,k))
wc(i,j,k) = 0.5*(w(i,j,k)+w(i,j,k-1))

end do; end do; end do

end subroutine reading_variables


