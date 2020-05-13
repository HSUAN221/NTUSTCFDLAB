	subroutine solver_1st_member(coef,jcoef,nx,ny,nz,ndim,mdim,iDx,Dxs,iDy,Dys,iDz,Dzs) 
	implicit none

	integer :: nx,ny,nz,i,j,k,ik,mdim,ndim
	integer ,dimension(1:mdim) :: jcoef 


	double precision ,dimension(1:ndim,1:mdim)::coef 
	double precision :: aa,bb,cc

	double precision, dimension (-1:nx+2) :: iDx
	double precision, dimension (-1:nx+2) :: Dxs
	double precision, dimension (-1:nx+2) :: iDy
	double precision, dimension (-1:nx+2) :: Dys
	double precision, dimension (-1:nx+2) :: iDz
	double precision, dimension (-1:nx+2) :: Dzs

	ik=1
	do k=1,nz; do j=1,ny; do i=1,nx
		coef(ik,1) =  1/iDx(i)/Dxs(i) + 1/iDx(i)/Dxs(i-1) + 1/iDy(j)/Dys(j) + 1/iDy(j)/Dys(j-1) + 1/iDz(k)/Dzs(k) + 1/iDz(k)/Dzs(k-1) 
		coef(ik,2) = -1/iDx(i)/Dxs(i) 		! AREA INTERPOLATION AT POINT CENTER (i+1,j,k)
		coef(ik,3) = -1/iDx(i)/Dxs(i-1) 	! AREA INTERPOLATION AT POINT CENTER (i-1,j,k)
		coef(ik,4) = -1/iDy(j)/Dys(j)		! AREA INTERPOLATION AT POINT CENTER (i,j+1,k)
		coef(ik,5) = -1/iDy(j)/Dys(j-1)  	! AREA INTERPOLATION AT POINT CENTER (i,j-1,k)
		coef(ik,6) = -1/iDz(k)/Dzs(k)		! AREA INTERPOLATION AT POINT CENTER (i,j,k+1)
		coef(ik,7) = -1/iDz(k)/Dzs(k-1)		! AREA INTERPOLATION AT POINT CENTER (i,j,k-1)
	    ik=ik+1
	enddo; enddo; enddo

	

	ik=1 
	do k=1,nz; do j=1,ny; do i=1,nx
		
              if(i.eq.1)then
                 aa=0
				 coef(ik,3)=0.0 ! AREA INTERPOLATION AT POINT CENTER (i-1,j,k)
              elseif(i.eq.nx)then
                 aa=-1/iDx(i)/Dxs(i)
				 coef(ik,2)=0.0 ! AREA INTERPOLATION AT POINT CENTER (i+1,j,k)
              else
                 aa=0.0
              endif

              if(j.eq.1)then
                 bb=0
				 coef(ik,5)=0.0 ! AREA INTERPOLATION AT POINT CENTER (i,j-1,k)
               elseif(j.eq.ny)then
                 bb=0	
				 coef(ik,4)=0.0 ! AREA INTERPOLATION AT POINT CENTER (i,j+1,k)
               else
                 bb=0.0
              endif

              if(k.eq.1)then
                 cc=0
				 coef(ik,7)=0.0 ! AREA INTERPOLATION AT POINT CENTER (i,j,k-1)
               elseif(k.eq.nz)then
                 cc=0
				 coef(ik,6)=0.0 ! AREA INTERPOLATION AT POINT CENTER (i,j,k+1)
               else
                 cc=0.0
              endif

	      	coef(ik,1)=coef(ik,1)+aa+bb+cc
	      	ik=ik+1

	enddo; enddo; enddo

	jcoef(1)=0
	jcoef(2)=1
	jcoef(3)=nx
	jcoef(4)=nx*ny 
	
	return
	
	
	
	end subroutine solver_1st_member
