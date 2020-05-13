subroutine DFIB_Sphere()
use variables
implicit none
real*4 :: DISTANCE,SDIST
real*4 :: DIAGONAL
integer :: l ,m ,n
real*4 :: xi
real*4 :: dxg, dyg, dzg
real*4 ,dimension(1:nSubGrids+1) :: SX
real*4 ,dimension(1:nSubGrids+1) :: SY
real*4 ,dimension(1:nSubGrids+1) :: SZ

!open(88,file='ETA.dat')

do k=kBgnVOS,kEndVOS 
  do j=jBgnVOS,jEndVOS
    do i=iBgnVOS,iEndVOS

         DIAGONAL = sqrt( iDx(i)*iDx(i) + iDy(j)*iDy(j) + iDz(k)*iDz(k) ) / 2.D0
         DISTANCE = sqrt((((X(i)+X(i+1))/2.D0)-xc)**2.D0+ &
                         (((Y(j)+Y(j+1))/2.D0)- ( yc + dis_Y ) )**2.D0+ &
                         (((Z(k)+Z(k+1))/2.D0)-zc)**2.D0  )

           if( abs(DISTANCE - r) < DIAGONAL ) then
              
              do l=1,nSubGrids+1,1
              do m=1,nSubGrids+1,1
              do n=1,nSubGrids+1,1
                dxg = iDx(i) / nSubGrids
                dyg = iDy(j) / nSubGrids
                dzg = iDz(k) / nSubGrids
                SX(n)=X(i)+(n-1)*dxg
                SY(m)=Y(j)+(m-1)*dyg
                SZ(l)=Z(k)+(l-1)*dzg
              end do
              end do
              end do

              xi = 0.0
              do l=1,nSubGrids
                do m=1,nSubGrids
                  do n=1,nSubGrids
                    
                    SDIST=SQRT((((SX(n)+SX(n+1))/2.D0)-xc)**2.D0+ &
                               (((SY(m)+SY(m+1))/2.D0)- ( yc + dis_Y ) )**2.D0+ &
                               (((SZ(l)+SZ(l+1))/2.D0)-zc)**2.D0  )
                    if(SDIST <= r) then
                      xi = xi + 1
                    end if

                  end do
                end do
              end do   
              ETA(i,j,k) = xi / (nSubGrids*nSubGrids*nSubGrids)
              
              !print*,'ETA(i,j,k) = ',ETA(i,j,k),'xi = ',xi,i,j,k
           
           else if (DISTANCE <= r) then
             ETA(i,j,k)=1.D0
           else
             ETA(i,j,k)=0.D0
           end if
                
    end do
  end do
end do



end subroutine DFIB_Sphere



subroutine DFIB_Cylinder()

use variables
implicit none
real*4 :: DISTANCE,SDIST
real*4 :: DIAGONAL
integer :: l ,m ,n
real*4 :: xi
real*4 :: dxg, dyg, dzg
real*4 ,dimension(1:nSubGrids+1) :: SX
real*4 ,dimension(1:nSubGrids+1) :: SY
real*4 ,dimension(1:nSubGrids+1) :: SZ




i=iBgnVOS
do k=kBgnVOS,kEndVOS
do j=jBgnVOS,jEndVOS

	 DIAGONAL = sqrt(  iDx(i)*iDx(i) + iDy(j)*iDy(j) + iDz(k)*iDz(k) ) / 2.D0

	 DISTANCE = sqrt( (( (Z(k)+Z(k+1)) / 2.D0 ) -   zc )           ** 2.D0 &
                 +  (( (Y(j)+Y(j+1)) / 2.D0 ) - ( yc + dis_Y ) ) ** 2.D0 )  

	   if( abs(DISTANCE - r) < DIAGONAL ) then
		  
		  xi = 0.0
		  
		 
			do m=1,nSubGrids+1; do n=1,nSubGrids+1
				dzg = iDz(k) / nSubGrids
				dyg = iDy(j) / nSubGrids
				SZ(n) = Z(k) + (n-1) * dzg
				SY(m) = Y(j) + (m-1) * dyg		
			end do; end do
		  
		  

		  
			do m=1,nSubGrids; do n=1,nSubGrids
				  SDIST = sqrt( (( (SZ(n)+SZ(n+1)) / 2.D0 ) -   zc )           ** 2.D0 &
                    +   (( (SY(m)+SY(m+1)) / 2.D0 ) - ( yc + dis_Y ) ) ** 2.D0 )
				  if(SDIST <= r) then
					xi = xi + 1
				  end if
			end do; end do
		   
		  
		  ETA(i,j,k) = xi / (nSubGrids*nSubGrids)
		  
		  !print*,'ETA(i,j,k) = ',ETA(i,j,k),'xi = ',xi,i,j,k
	   
	   else if (DISTANCE > r) then
		 ETA(i,j,k)=0.D0
	   else
		 ETA(i,j,k)=1.D0
	   end if
			
end do
end do


do k=kBgnVOS,kEndVOS
do j=jBgnVOS,jEndVOS
do i=iBgnVOS+1,iEndVOS
	ETA(i,j,k) = ETA(iBgnVOS,j,k)
end do
end do
end do    

end subroutine DFIB_Cylinder






