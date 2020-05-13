subroutine discretisation_upwind()
use variables 
implicit none


!!! u_star calculation (x component) !!!

!$OMP PARALLEL DO PRIVATE(i,j,k)
do k=1,nz
    do j=1,ny
        do i=1,nx-1

u_tilde_x1 = 0.5 * (u(i,j,k)+u(i+1,j,k)) !1 for +(i or j)
u_tilde_x2 = 0.5 * (u(i-1,j,k)+u(i,j,k))  !2 for -(i or j)
v_tilde_x1 = 0.5 * (v(i,j,k)+v(i+1,j,k))
v_tilde_x2 = 0.5 * (v(i,j-1,k)+v(i+1,j-1,k))
w_tilde_x1 = 0.5 * (w(i+1,j,k)+w(i,j,k))
w_tilde_x2 = 0.5 * (w(i,j,k-1)+w(i+1,j,k-1))

if (u_tilde_x1 >= 0) then 
        
        ue=u(i,j,k)
    else 
        ue=u(i+1,j,k)
end if


if (u_tilde_x2 >= 0) then

        uw=u(i-1,j,k)
    else 
        uw=u(i,j,k)
end if


if (v_tilde_x1 >= 0) then

        un=u(i,j,k)
    else 
        un=u(i,j+1,k)
end if


if (v_tilde_x2 >= 0) then 
        
        us=u(i,j-1,k)
    else 
        us=u(i,j,k)
end if

if (w_tilde_x1 >= 0) then

        uf=u(i,j,k)
    else 
        uf=u(i,j,k+1)
end if


if (w_tilde_x2 >= 0) then 
        
        ub=u(i,j,k-1)
    else 
        ub=u(i,j,k)
end if

		if (Q1 == 1) then 
            delta = ( Dxs(i)*iDy(j)*iDz(k) ) ** (1./3)
			nut = (Cs*delta)**2 * (   0.5*(1/Dxs(i)/Dxs(i))*((u(i+1,j,k)+u(i,j,k))-(u(i,j,k)+u(i-1,j,k)))**2 + &
									  0.5*(1/iDy(j)/iDy(j))*((v(i+1,j,k)+v(i,j,k))-(v(i+1,j-1,k)+v(i,j-1,k)))**2 + &
									  0.5*(1/iDz(k)/iDz(k))*((w(i+1,j,k)+w(i,j,k))-(w(i+1,j,k-1)+w(i,j,k-1)))**2 + &

									  ( 0.5*(1/iDy(j))*((u(i,j+1,k)+u(i,j,k))-(u(i,j,k)+u(i,j-1,k)))+ &
										0.5*(1/Dxs(i))*((v(i+1,j,k)+v(i+1,j-1,k))-(v(i,j,k)+v(i,j-1,k))) )**2 + &

									  ( 0.5*(1/iDz(k))*((u(i,j,k+1)+u(i,j,k))-(u(i,j,k)+u(i,j,k-1)))+ &
										0.5*(1/Dxs(i))*((w(i+1,j,k)+w(i+1,j,k-1))-(w(i,j,k)+w(i,j,k-1))) )**2 + &

									  ( 0.5*(1/iDz(k))*((v(i+1,j,k+1)+v(i,j-1,k))-(v(i+1,j,k-1)+v(i,j-1,k)))+ &
										0.5*(1/iDy(j))*((w(i+1,j+1,k)+w(i,j,k-1))-(w(i+1,j-1,k)+w(i,j,k-1))) )**2       )**0.5
		end if


		u_star(i,j,k) = u(i,j,k)-dt*(u_tilde_x1*ue-u_tilde_x2*uw) / Dxs(i) &
                                -dt*(v_tilde_x1*un-v_tilde_x2*us) / iDy(j) &
                                -dt*(w_tilde_x1*uf-w_tilde_x2*ub) / iDz(k) &

                        +(nu+nut*Q1)*dt*( (u(i+1,j,k)-u(i,j,k)) / iDx(i+1) - (u(i,j,k)-u(i-1,j,k)) / iDx(i) ) / Dxs(i) &
                        +(nu+nut*Q1)*dt*( (u(i,j+1,k)-u(i,j,k)) / Dys(j) - (u(i,j,k)-u(i,j-1,k)) / Dys(j-1) ) / iDy(j) &
                        +(nu+nut*Q1)*dt*( (u(i,j,k+1)-u(i,j,k)) / Dzs(k) - (u(i,j,k)-u(i,j,k-1)) / Dzs(k-1) ) / iDz(k) 

        end do
    end do
end do
!$OMP END PARALLEL DO


!!! v_star calculation (y component) !!!

!$OMP PARALLEL DO PRIVATE(i,j,k)
do k=1,nz
    do j=1,ny-1
        do i=1,nx

u_tilde_y1 = 0.5 * (u(i,j,k)+u(i,j+1,k))
u_tilde_y2 = 0.5 * (u(i-1,j,k)+u(i-1,j+1,k))
v_tilde_y1 = 0.5 * (v(i,j,k)+v(i,j+1,k))
v_tilde_y2 = 0.5 * (v(i,j-1,k)+v(i,j,k))
w_tilde_y1 = 0.5 * (w(i,j,k)+w(i,j+1,k))
w_tilde_y2 = 0.5 * (w(i,j,k-1)+w(i,j+1,k-1))

if (u_tilde_y1 >= 0) then 

        ve=v(i,j,k)
    else 
        ve=v(i+1,j,k)
end if


if (u_tilde_y2 >= 0) then 
        
        vw=v(i-1,j,k)
    else 
        vw=v(i,j,k)
end if


if (v_tilde_y1 >= 0) then 
        
        vn=v(i,j,k)
    else 
        vn=v(i,j+1,k)
end if


if (v_tilde_y2 >= 0) then 
        
        vs=v(i,j-1,k)
    else 
        vs=v(i,j,k)
end if

if (w_tilde_y1 >= 0) then 
        
        vf=v(i,j,k)
    else 
        vf=v(i,j,k+1)
end if

if (w_tilde_y2 >= 0) then 
        
        vb=v(i,j,k-1)
    else 
        vb=v(i,j,k)
end if

		if (Q1 == 1) then
            delta = ( iDx(i)*Dys(j)*iDz(k) ) ** (1./3)
			nut = (Cs*delta)**2 * (   0.5*(1/iDx(i)/iDx(i))*((u(i,j+1,k)+u(i,j,k))-(u(i-1,j+1,k)+u(i-1,j,k)))**2 + &
									  0.5*(1/Dys(j)/Dys(j))*((v(i,j+1,k)+v(i,j,k))-(v(i,j,k)+v(i,j-1,k)))**2 + &
									  0.5*(1/iDz(k)/iDz(k))*((w(i,j+1,k)+w(i,j,k))-(w(i,j+1,k-1)+w(i,j,k-1)))**2 + &

									  ( 0.5*(1/Dys(j))*((u(i,j+1,k)+u(i-1,j+1,k))-(u(i,j,k)+u(i-1,j,k)))+ &
										0.5*(1/iDx(i))*((v(i+1,j,k)+v(i,j,k))-(v(i,j,k)+v(i-1,j,k))) )**2 + &

									  ( 0.5*(1/iDz(k))*((u(i,j+1,k+1)+u(i-1,j,k))-(u(i,j+1,k-1)+u(i-1,j,k-1)))+ &
										0.5*(1/iDx(i))*((w(i+1,j+1,k)+w(i,j,k-1))-(w(i-1,j+1,k)+w(i,j,k-1))) )**2 + &

									  ( 0.5*(1/iDz(k))*((v(i,j,k+1)+v(i,j,k))-(v(i,j,k-1)+v(i,j,k)))+ &
										0.5*(1/Dys(j))*((w(i,j+1,k)+w(i,j+1,k-1))-(w(i,j,k)+w(i,j,k-1))) )**2       )**0.5
		end if


		v_star(i,j,k) = v(i,j,k)-dt*(u_tilde_y1*ve-u_tilde_y2*vw) / iDx(i) &
                                -dt*(v_tilde_y1*vn-v_tilde_y2*vs) / Dys(j) &
                                -dt*(w_tilde_y1*vf-w_tilde_y2*vb) / iDz(k) &

                        +(nu+nut*Q1)*dt*( (v(i+1,j,k)-v(i,j,k)) / Dxs(i) - (v(i,j,k) - v(i-1,j,k)) / Dxs(i-1) ) / iDx(i) &
                        +(nu+nut*Q1)*dt*( (v(i,j+1,k)-v(i,j,k)) / iDy(j+1) - (v(i,j,k) - v(i,j-1,k)) / iDy(j) ) / Dys(j) &
                        +(nu+nut*Q1)*dt*( (v(i,j,k+1)-v(i,j,k)) / Dzs(k) - (v(i,j,k) - v(i,j,k-1)) / Dzs(k-1) ) / iDz(k)

        end do
    end do
end do
!$OMP END PARALLEL DO

!!! w_star calculation (z component) !!!

!$OMP PARALLEL DO PRIVATE(i,j,k)
do k=1,nz-1
    do j=1,ny
        do i=1,nx

u_tilde_z1 = 0.5 * (u(i,j,k+1)+u(i,j,k)) !1 for +(i or j)
u_tilde_z2 = 0.5 * (u(i-1,j,k+1)+u(i-1,j,k))  !2 for -(i or j)
v_tilde_z1 = 0.5 * (v(i,j,k)+v(i+1,j,k+1))
v_tilde_z2 = 0.5 * (v(i,j-1,k+1)+v(i,j-1,k))
w_tilde_z1 = 0.5 * (w(i,j,k+1)+w(i,j,k))
w_tilde_z2 = 0.5 * (w(i,j,k)+w(i,j,k-1))

if (u_tilde_z1 >= 0) then 

        we=w(i,j,k)
    else 
        we=w(i+1,j,k)
end if


if (u_tilde_z2 >= 0) then 
        
        ww=w(i-1,j,k)
    else 
        ww=w(i,j,k)
end if


if (v_tilde_z1 >= 0) then 
        
        wn=w(i,j,k)
    else 
        wn=w(i,j+1,k)
end if


if (v_tilde_z2 >= 0) then 
        
        ws=w(i,j-1,k)
    else 
        ws=w(i,j,k)
end if

if (w_tilde_z1 >= 0) then 
        
        wf=w(i,j,k)
    else 
        wf=w(i,j,k+1)
end if

if (w_tilde_z2 >= 0) then 
        
        wb=w(i,j,k-1)
    else 
        wb=w(i,j,k)
end if

		if (Q1 == 1) then
            delta = ( iDx(i)*iDy(j)*Dzs(k) ) ** (1./3)
            nut = (Cs*delta)**2 * (   0.5*(1/iDx(i)/iDx(i))*((u(i,j,k+1)+u(i,j,k))-(u(i-1,j,k+1)+u(i-1,j,k)))**2 + &
                                      0.5*(1/iDy(j)/iDy(j))*((v(i,j,k+1)+v(i,j,k))-(v(i,j-1,k+1)+v(i,j-1,k)))**2 + &
                                      0.5*(1/Dzs(k)/Dzs(k))*((w(i,j,k+1)+w(i,j,k))-(w(i,j,k)+w(i,j,k-1)))**2 + &

                                    ( 0.5*(1/iDy(j))*((u(i,j+1,k)+u(i,j,k))-(u(i,j,k)+u(i,j-1,k)))+ &
                                      0.5*(1/iDx(i))*((v(i+1,j,k+1)+v(i,j-1,k))-(v(i-1,j,k+1)+v(i,j-1,k))) )**2 + &

                                    ( 0.5*(1/Dzs(k))*((u(i,j,k+1)+u(i-1,j,k+1))-(u(i,j,k)+u(i-1,j,k)))+ &
                                      0.5*(1/iDx(i))*((w(i+1,j,k)+w(i,j,k))-(w(i-1,j,k)+w(i,j,k))) )**2 + &

                                    ( 0.5*(1/Dzs(k))*((v(i,j,k+1)+v(i,j-1,k+1))-(v(i,j,k)+v(i,j-1,k)))+ &
                                      0.5*(1/iDy(j))*((w(i,j+1,k)+w(i,j,k))-(w(i,j-1,k)+w(i,j,k))) )**2       )**0.5
        end if


		w_star(i,j,k) = w(i,j,k)-dt*(u_tilde_z1*we-u_tilde_z2*ww) / iDx(i) &
                                -dt*(v_tilde_z1*wn-v_tilde_z2*ws) / iDy(j) &
                                -dt*(w_tilde_z1*wf-w_tilde_z2*wb) / Dzs(k) &
                                
                        +(nu+nut*Q1)*dt*( ( w(i+1,j,k)-w(i,j,k) ) / Dxs(i) -  ( w(i,j,k)-w(i-1,j,k) ) / Dxs(i-1)  ) / iDx(i) &
                        +(nu+nut*Q1)*dt*( ( w(i,j+1,k)-w(i,j,k) ) / Dys(j) -  ( w(i,j,k)-w(i,j-1,k) ) / Dys(j-1)  ) / iDy(j) &
                        +(nu+nut*Q1)*dt*( ( w(i,j,k+1)-w(i,j,k) ) / iDz(k+1) - ( w(i,j,k)-w(i,j,k-1) ) / iDz(k)   ) / Dzs(k)

        end do
    end do
end do
!$OMP END PARALLEL DO

end subroutine discretisation_upwind
