subroutine gauss_seidel() 
    use variables
    implicit none



    ik=0
    pChangeMax = 1.0
    pChangeMax_= 1.0

    do while (pChangeMax_>zeta .AND. ik < itmax)

    ik=ik+1
    pChangeMax = 0.0
    pChangeMax_= 0.0
    mChangeMax = 0.0

    !----------data transformation among nodes----------!
       
    icount = (nx+4)*(ny+4)

    itag = 250
    call MPI_SENDRECV( p(-1,-1,istart), icount, MPI_REAL8, l_nbr, itag, &
                        p(-1,-1,iend+1), icount, MPI_REAL8, r_nbr, itag, MPI_COMM_WORLD, status, ierr )

    itag = 260
    call MPI_SENDRECV( p(-1,-1,iend), icount, MPI_REAL8, r_nbr, itag, &
                        p(-1,-1,istart-1), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )
                
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    
    !----------data transformation among nodes----------!

     
        do k=istart,iend
        !$OMP PARALLEL DO PRIVATE(i)  
        do j=1,ny; do i=1,nx

            mChange = ( u_star(i,j,k) - u_star(i-1,j,k) ) * iDy(j) * iDz(k) &
                    + ( v_star(i,j,k) - v_star(i,j-1,k) ) * iDx(i) * iDz(k) &
                    + ( w_star(i,j,k) - w_star(i,j,k-1) ) * iDx(i) * iDy(j)

            pNEW = (- p(i+1,j,k) * iDy(j) * iDz(k) / Dxs(i) &
                    - p(i-1,j,k) * iDy(j) * iDz(k) / Dxs(i-1) &
                    - p(i,j+1,k) * iDx(i) * iDz(k) / Dys(j) &
                    - p(i,j-1,k) * iDx(i) * iDz(k) / Dys(j-1) &
                    - p(i,j,k+1) * iDx(i) * iDy(j) / Dzs(k) &
                    - p(i,j,k-1) * iDx(i) * iDy(j) / Dzs(k-1) &
                    + mChange / dt) /  (- iDy(j) * iDz(k) / Dxs(i) - iDy(j) * iDz(k) / Dxs(i-1) &
                                        - iDx(i) * iDz(k) / Dys(j) - iDx(i) * iDz(k) / Dys(j-1) &
                                        - iDx(i) * iDy(j) / Dzs(k) - iDx(i) * iDy(j) / Dzs(k-1) )
            
            pChange = abs(pNew - P(i,j,k))

            P(i,j,k) = P(i,j,k) + ( omega * (pNew - P(i,j,k)) )

            if (abs(mChange) > mChangeMax) then
                mChangeMax = abs(mChange)
            end if

            if (pChange > pChangeMax) then
                pChangeMax = pChange
            end if

             
        enddo; enddo
        !$OMP END PARALLEL DO
        enddo

        call MPI_ALLREDUCE( pChangeMax, pChangeMax_, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr )

        !if(myid==master)then
        !    write(*,*) pChangeMax_
        !end if
    
    end do









    do k=istart,iend
    !$OMP PARALLEL DO PRIVATE(i)  
    do j=1,ny; do i=1,nx

        pre(i,j,k) = p(i,j,k)

    end do; end do
    !$OMP END PARALLEL DO
    end do

    if(myid==master)then

        
            print*, '          '
            print*, 'Iterations Gauss Seidel =',ik
        

    end if


end subroutine gauss_seidel
