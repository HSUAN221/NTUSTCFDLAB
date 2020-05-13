subroutine gridder_Unequal()
use variables
implicit none

!------------Following are calculated using the above variables------------!
real*8 ,parameter :: xSBgn = GridderXc - lxSml/2.0
real*8 ,parameter :: xSEnd = GridderXc + lxSml/2.0
real*8 ,parameter :: ySBgn = GridderYc - lySml/2.0
real*8 ,parameter :: ySEnd = GridderYc + lySml/2.0
real*8 ,parameter :: zSBgn = GridderZc
real*8 ,parameter :: zSEnd = GridderZc + lzSml


real*8 ,parameter :: zMBgn = GridderZc
real*8 ,parameter :: zMEnd = GridderZc + lzMid

real*8 ,parameter :: yMBgn = GridderYc - lyMid/2.0
real*8 ,parameter :: yMEnd = GridderYc + lyMid/2.0

real*8 :: xNextLrgValue
real*8 :: xNextSmlValue

real*8 :: yNextLrgValue
real*8 :: yNextMidValue
real*8 :: yNextSmlValue

real*8 :: zNextLrgValue
real*8 :: zNextMidValue
real*8 :: zNextSmlValue
!------------Following are calculated using the above variables------------!


!----------------Unequal grid intervals----------------!
    do i=1,nx+3
        if(i == 1) then
            X(i) = 0.0
            X(i-1) = X(i) - dx
            X(i-2) = X(i-1) - dx
        else
            X(i) = X(i-1) + dx
        end if 
    end do


    do j=1,ny+1
        yNextLrgValue =  Y(j-1) + dy   
        yNextMidValue =  Y(j-1) + dyMid   
        yNextSmlValue =  Y(j-1) + dySml   

        if (j==1) then
            Y(j) = 0.0
        else if(yNextLrgValue > yMBgn) then
            if(yNextMidValue>yMEnd) then
                Y(j) = yNextLrgValue
            else if(yNextMidValue>ySBgn) then
                if(yNextSmlValue>ySEnd) then
                    Y(j) = yNextMidValue
                else
                    Y(j) = yNextSmlValue
                end if
            else
                Y(j) = yNextMidValue
            end if
        else
            Y(j) = yNextLrgValue
        end if
    end do

   

    do k=1,nz+1
        zNextLrgValue =  Z(k-1) + dz   
        zNextMidValue =  Z(k-1) + dzMid   
        zNextSmlValue =  Z(k-1) + dzSml   

        if (k==1) then
            Z(k) = 0.0
        else if(zNextLrgValue > zMBgn) then
            if(zNextMidValue>zMEnd) then
                Z(k) = zNextLrgValue
            else if(zNextMidValue>zSBgn) then
                if(zNextSmlValue>zSEnd) then
                    Z(k) = zNextMidValue
                else
                    Z(k) = zNextSmlValue
                end if
            else
                Z(k) = zNextMidValue
            end if
        else
            Z(k) = zNextLrgValue
        end if
    end do
    !----------------Unequal grid intervals----------------!

    !Define each of the directional grid lengths
    do i=1,nx-1

        iDx(i) = ( X(i+1) - X(i) )
        Dxs(i) = ( X(i+2) - X(i) ) / 2.0

    end do

    do j=1,ny-1

        iDy(j) = ( Y(j+1) - Y(j) )
        Dys(j) = ( Y(j+2) - Y(j) ) / 2.0

    end do

    do k=1,nz-1

        iDz(k) = ( Z(k+1) - Z(k) )
        Dzs(k) = ( Z(k+2) - Z(k) ) / 2.0

    end do



    !Ghost boundary grid lengths
    iDx(0) = iDx(1)
    iDx(-1) = iDx(1)
    iDx(nx) = X(nx+1) - X(nx)
    iDx(nx+1) = iDx(nx)
    iDx(nx+2) = iDx(nx)

    Dxs(0) = Dxs(1)
    Dxs(-1) = Dxs(1)
    Dxs(nx) = Dxs(nx-1)
    Dxs(nx+1) = Dxs(nx-1)
    Dxs(nx+2) = Dxs(nx-1)


    iDy(0) = iDy(1)
    iDy(-1) = iDy(1)
    iDy(ny) = Y(ny+1) - Y(ny)
    iDy(ny+1) = iDy(ny)
    iDy(ny+2) = iDy(ny)

    Dys(0) = Dys(1)
    Dys(-1) = Dys(1)
    Dys(ny) = Dys(ny-1)
    Dys(ny+1) = Dys(ny-1)
    Dys(ny+2) = Dys(ny-1)


    iDz(0) = iDz(1)
    iDz(-1) = iDz(1)
    iDz(nz) = Z(nz+1) - Z(nz)
    iDz(nz+1) = iDz(nz)
    iDz(nz+2) = iDz(nz)

    Dzs(0) = Dzs(1)
    Dzs(-1) = Dzs(1)
    Dzs(nz) = Dzs(nz-1)
    Dzs(nz+1) = Dzs(nz-1)
    Dzs(nz+2) = Dzs(nz-1)

    !Modifying the index of X, Y and Z arrays to represent the actual grid
    do i=1,nx+1
        Xa(i) = X(i)
    end do

    do j=1,ny+1
        Ya(j) = Y(j)
    end do

    do k=1,nz+1
        Za(k) = Z(k)
    end do

    !Defining the midpoint values of the grids
    do i=1,nx
        Xs(i) = ( Xa(i+1) + Xa(i) ) / 2.0
    end do

    do j=1,ny
        Ys(j) = ( Ya(j+1) + Ya(j) ) / 2.0
    end do

    do k=1,nz
        Zs(k) = ( Za(k+1) + Za(k) ) / 2.0
    end do



    !Output values of the grids
    do k=1,nz; do j=1,ny; do i=1,nx
        Xout(i,j,k) = Xs(i)
    enddo; enddo; enddo

    do k=1,nz; do j=1,ny; do i=1,nx
        Yout(i,j,k) = Ys(j)
    enddo; enddo; enddo

    do k=1,nz; do j=1,ny; do i=1,nx
        Zout(i,j,k) = Zs(k)
    enddo; enddo; enddo

    if(myid==master)then
        open (unit=1,form='unformatted',file='mesh.x')
        write(1) nblocks
        write(1) nx, ny, nz

        write(1)    (((Xout(i,j,k),i=1,nx),j=1,ny),k=1,nz), &
                    (((Yout(i,j,k),i=1,nx),j=1,ny),k=1,nz), &
                    (((Zout(i,j,k),i=1,nx),j=1,ny),k=1,nz)

        

        close(1)
    end if

    !if(myid==master)then
    !open (1,file='mesh.dat',position='append')
    !do j=1,ny
    !    write (1,'(E12.5)')  &
    !    Ys(j)
    !enddo
    !close(1)
    !end if

end subroutine gridder_Unequal




subroutine gridder_equal()
    use variables
    implicit none

    

    !-----------------unequal grid intervals----------------!
    do i=1,nx+3
        if(i == 1) then
            X(i) = 0.0
            X(i-1) = X(i) - dx
            X(i-2) = X(i-1) - dx
        else
            X(i) = X(i-1) + dx
        end if 
    end do

    do j=1,ny+3
        if(j == 1) then
            Y(j) = 0.0
            Y(j-1) = Y(j) - dy
            Y(j-2) = Y(j-1) - dy
        else
            Y(j) = Y(j-1) + dy
        end if 
    end do

    do k=1,nz+3
        if(k == 1) then
            Z(k) = 0.0
            Z(k-1) = Z(k) - dz
            Z(k-2) = Z(k-1) - dz
        else
            Z(k) = Z(k-1) + dz
        end if 
    end do
    !-----------------unequal grid intervals----------------!

    !Define each of the directional grid lengths
    do i=1,nx-1

        iDx(i) = ( X(i+1) - X(i) )
        Dxs(i) = ( X(i+2) - X(i) ) / 2.0

    end do

    do j=1,ny-1

        iDy(j) = ( Y(j+1) - Y(j) )
        Dys(j) = ( Y(j+2) - Y(j) ) / 2.0

    end do

    do k=1,nz-1

        iDz(k) = ( Z(k+1) - Z(k) )
        Dzs(k) = ( Z(k+2) - Z(k) ) / 2.0

    end do



    !Ghost boundary grid lengths
    iDx(0) = iDx(1)
    iDx(-1) = iDx(1)
    iDx(nx) = X(nx+1) - X(nx)
    iDx(nx+1) = iDx(nx)
    iDx(nx+2) = iDx(nx)

    Dxs(0) = Dxs(1)
    Dxs(-1) = Dxs(1)
    Dxs(nx) = Dxs(nx-1)
    Dxs(nx+1) = Dxs(nx-1)
    Dxs(nx+2) = Dxs(nx-1)


    iDy(0) = iDy(1)
    iDy(-1) = iDy(1)
    iDy(ny) = Y(ny+1) - Y(ny)
    iDy(ny+1) = iDy(ny)
    iDy(ny+2) = iDy(ny)

    Dys(0) = Dys(1)
    Dys(-1) = Dys(1)
    Dys(ny) = Dys(ny-1)
    Dys(ny+1) = Dys(ny-1)
    Dys(ny+2) = Dys(ny-1)


    iDz(0) = iDz(1)
    iDz(-1) = iDz(1)
    iDz(nz) = Z(nz+1) - Z(nz)
    iDz(nz+1) = iDz(nz)
    iDz(nz+2) = iDz(nz)

    Dzs(0) = Dzs(1)
    Dzs(-1) = Dzs(1)
    Dzs(nz) = Dzs(nz-1)
    Dzs(nz+1) = Dzs(nz-1)
    Dzs(nz+2) = Dzs(nz-1)

    !Modifying the index of X, Y and Z arrays to represent the actual grid
    do i=1,nx+1
        Xa(i) = X(i)
    end do

    do j=1,ny+1
        Ya(j) = Y(j)
    end do

    do k=1,nz+1
        Za(k) = Z(k)
    end do

    !Defining the midpoint values of the grids
    do i=1,nx
        Xs(i) = ( Xa(i+1) + Xa(i) ) / 2.0
    end do

    do j=1,ny
        Ys(j) = ( Ya(j+1) + Ya(j) ) / 2.0
    end do

    do k=1,nz
        Zs(k) = ( Za(k+1) + Za(k) ) / 2.0
    end do

    !Output values of the grids
    do k=1,nz; do j=1,ny; do i=1,nx
        Xout(i,j,k) = Xs(i)
    enddo; enddo; enddo

    do k=1,nz; do j=1,ny; do i=1,nx
        Yout(i,j,k) = Ys(j)
    enddo; enddo; enddo

    do k=1,nz; do j=1,ny; do i=1,nx
        Zout(i,j,k) = Zs(k)
    enddo; enddo; enddo

    if(myid==master)then
        open (unit=1,form='unformatted',file='mesh.x')
        write(1) nblocks
        write(1) nx, ny, nz

        write(1)  (((Xout(i,j,k),i=1,nx),j=1,ny),k=1,nz), &
                  (((Yout(i,j,k),i=1,nx),j=1,ny),k=1,nz), &
                  (((Zout(i,j,k),i=1,nx),j=1,ny),k=1,nz)

        

        close(1)
    end if


end subroutine gridder_equal

