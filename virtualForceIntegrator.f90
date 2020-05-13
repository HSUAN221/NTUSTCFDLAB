
 !    only cpu 0 answer correct
subroutine virtualForceIntegrator()
	use variables
	implicit none

    !---------------------------------------------------!
    !    LOCAL VARIABLES                                !
    !---------------------------------------------------!

	integer ,parameter :: iBgnVOS = 1
	integer ,parameter :: iEndVOS = nx
	integer ,parameter :: jBgnVOS = 1
	integer ,parameter :: jEndVOS = ny
	integer ,parameter :: kBgnVOS = 1
	integer ,parameter :: kEndVOS = nz




	!Integrating components in z-direction using Simpson's 1/3 rule
	do j=jBgnVOS,jEndVOS ; do i= iBgnVOS,iEndVOS

		FXz(i,j) = 0.0
		FYz(i,j) = 0.0
		
		do k= kBgnVOS,kEndVOS

			
				FXz(i,j) = FXz(i,j) + ( FZ(i,j,k-1)*Dzs(k-1) + 4.0*FZ(i,j,k)*Dzs(k) + FZ(i,j,k+1)*Dzs(k+1) ) / 6.0
				FYz(i,j) = FYz(i,j) + ( FY(i,j,k-1)*Dzs(k-1) + 4.0*FY(i,j,k)*Dzs(k) + FY(i,j,k+1)*Dzs(k+1) ) / 6.0
			
			
		end do

	end do; end do



	!Integrating components in y-direction using Simpson's 1/3 rule
	do i= iBgnVOS,iEndVOS

		FXy(i) = 0.0
		FYy(i) = 0.0
		
		do j=jBgnVOS,jEndVOS

			if (FXz(i,j) /= 0.0 ) then
				FXy(i) = FXy(i) + (FXz(i,j-1)*Dys(j-1) + 4.0*FXz(i,j)*Dys(j) + FXz(i,j+1)*Dys(j+1)) / 6.0
			end if
			
			if (FYz(i,j) /= 0.0 ) then
				FYy(i) = FYy(i) + (FYz(i,j-1)*Dys(j-1) + 4.0*FYz(i,j)*Dys(j) + FYz(i,j+1)*Dys(j+1)) / 6.0
			end if
			
		end do
		
	end do



	!Integrating components in x-direction using Simpson's 1/3 rule
	totalFX = 0.0
	totalFY = 0.0
	do i= iBgnVOS,iEndVOS

		if (FXy(i) /= 0.0 ) then
			totalFX = totalFX + (FXy(i-1)*Dxs(i-1) + 4.0*FXy(i)*Dxs(i) + FXy(i+1)*Dxs(i+1) ) / 6.0
		end if
		
		if (FYy(i) /= 0.0 ) then
			totalFY = totalFy + (FYy(i-1)*Dxs(i-1) + 4.0*FYy(i)*Dxs(i) + FYy(i+1)*Dxs(i+1) ) / 6.0
		end if

	end do




	cDrag = (-2.0) * totalFX / lx                 
	cLift = (-2.0) * totalFY / lx



end subroutine virtualForceIntegrator





subroutine virtualForceIntegrator_nima()
	use variables
	implicit none

	totalFX = 0.0
	totalFY = 0.0

    !$OMP PARALLEL DO PRIVATE(j,i)  
    do k=1,nz;do j=1,ny;do i=1,nx

	totalFX = totalFX + ( FZ(i,j,k) * ( iDx(i)*iDy(j)*iDz(k) ) )

    end do;end do;end do
    !$OMP END PARALLEL DO


    !$OMP PARALLEL DO PRIVATE(j,i)  
    do k=1,nz;do j=1,ny;do i=1,nx

	totalFY = totalFY + ( FY(i,j,k) * ( iDx(i)*iDy(j)*iDz(k) ) )

    end do;end do;end do
    !$OMP END PARALLEL DO

    
	cDrag = (-2.0) * totalFX / lx                  
	cLift = (-2.0) * totalFY / lx 

	!cDrag = (-2.0) * totalFX                 
	!cLift = (-2.0) * totalFY


end subroutine virtualForceIntegrator_nima




