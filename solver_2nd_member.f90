subroutine solver_2nd_member() ! Calculate the second membre of equation AX=B
use variables
implicit none

	
	do k=istart,iend;
	!$OMP PARALLEL DO PRIVATE(j,i)   
	do j=1,ny; do i=1,nx

		div(i,j,k)= ddt*(   (u_star(i,j,k)-u_star(i-1,j,k))/iDx(i)+ &
							(v_star(i,j,k)-v_star(i,j-1,k))/iDy(j)+ &
							(w_star(i,j,k)-w_star(i,j,k-1))/iDz(k)  )

		if(i.eq.1)then
                 div(i,j,k) = div(i,j,k) - p_west/iDx(i)/Dxs(i-1)
		elseif(i.eq.nx)then
				 div(i,j,k) = div(i,j,k) - 0
        endif
		
		if(j.eq.1)then
                 div(i,j,k) = div(i,j,k) - p_south/iDy(j)/Dys(j-1)
		elseif(j.eq.ny)then
				 div(i,j,k) = div(i,j,k) - p_north/iDy(j)/Dys(j)
        endif
		
		if(k.eq.1)then
                 div(i,j,k) = div(i,j,k) - p_back/iDz(k)/Dzs(k-1)
		elseif(k.eq.nz)then
				 div(i,j,k) = div(i,j,k) - p_front/iDz(k)/Dzs(k)
        endif

	enddo; enddo
	!$OMP END PARALLEL DO
	enddo

	
	


							




	!cpu0 	ik = 1 ～ ik = 10*(nx)*(ny)
	!cpu1	ik = 10*(nx)*(ny) + 1 ～ ik = 20*(nx)*(ny) 
	!cpu2	ik = 20*(nx)*(ny) + 1 ～ ik = 30*(nx)*(ny)
	!cpu3	ik = 30*(nx)*(ny) + 1 ～ ik = 40*(nx)*(ny)

	
	ik = 1 + myid * gcount(myid) * (nx) * (ny)
	do k=istart,iend; do j=1,ny; do i=1,nx
	
        div1(ik) = -div(i,j,k)
	    x1(ik)=0.
	    ik=ik+1

	enddo; enddo; enddo




end subroutine solver_2nd_member
