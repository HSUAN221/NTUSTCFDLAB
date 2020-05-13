subroutine discretisation_QUICK()
use variables 
implicit none

!!! u_star calculation (x component) !!!

!$OMP PARALLEL DO PRIVATE(i,j,k)   
      do k=1,nz
    do j=1,ny
  do i=1,nx-1  

     if (u(i,j,k) > 0) then 
    
        ue  = 0.5*(u(i,j,k)     + u(i+1,j,k)  ) -0.125*iDx(i+1)*iDx(i+1)/Dxs(i) &
               *(  (u(i+1,j,k)   - u(i,j,k)   ) / iDx(i+1) &
                - (u(i,j,k)     - u(i-1,j,k)  ) / iDx(i) )

        uw  = 0.5*(u(i-1,j,k)   + u(i,j,k)    ) -0.125*iDx(i)*iDx(i)/Dxs(i-1) &
               *(  (u(i,j,k)     - u(i-1,j,k) ) / iDx(i) &
                - (u(i-1,j,k)   - u(i-2,j,k)  ) / iDx(i-1) ) 

        vnu = 0.5*(v(i,j,k)     + v(i+1,j,k)  ) -0.125*Dxs(i)*Dxs(i)/iDx(i) &
              *(   (v(i+1,j,k)   - v(i,j,k)   ) / Dxs(i) &
                - (v(i,j,k)     - v(i-1,j,k)  ) / Dxs(i-1) )

        vsu = 0.5*(v(i,j-1,k)   + v(i+1,j-1,k)) -0.125*Dxs(i)*Dxs(i)/iDx(i) &
              *(   (v(i+1,j-1,k) - v(i,j-1,k) ) / Dxs(i) &
                - (v(i,j-1,k)   - v(i-1,j-1,k)) / Dxs(i-1) ) 

        wfu = 0.5*(w(i,j,k)     + w(i+1,j,k)  ) -0.125*Dxs(i)*Dxs(i)/iDx(i) &
              *(   (w(i+1,j,k)   - w(i,j,k)   ) / Dxs(i) &
                - (w(i,j,k)     - w(i-1,j,k)  ) / Dxs(i-1) )

        wbu = 0.5*(w(i,j,k-1)   + w(i+1,j,k-1)) -0.125*Dxs(i)*Dxs(i)/iDx(i) &
              *(   (w(i+1,j,k-1) - w(i,j,k-1) ) / Dxs(i) &
                - (w(i,j,k-1)   - w(i-1,j,k-1)) / Dxs(i-1) )

    else

        ue  = 0.5*(u(i,j,k)     + u(i+1,j,k)  ) -0.125*iDx(i+1)*iDx(i+1)/Dxs(i+1) &
              *(   (u(i+2,j,k)   - u(i+1,j,k) ) / iDx(i+2) &
                - (u(i+1,j,k)   - u(i,j,k)    ) / iDx(i+1) ) 

        uw  = 0.5*(u(i-1,j,k)   + u(i,j,k)    ) -0.125*iDx(i)*iDx(i)/Dxs(i) &
              *(   (u(i+1,j,k)   - u(i,j,k)   ) / iDx(i+1) &
                - (u(i,j,k)     - u(i-1,j,k)  ) / iDx(i) ) 

        vnu = 0.5*(v(i,j,k)     + v(i+1,j,k)  ) -0.125*Dxs(i)*Dxs(i)/iDx(i+1) &
              *(   (v(i+2,j,k)   - v(i+1,j,k) ) / Dxs(i+1) &
                - (v(i+1,j,k)   - v(i,j,k)    ) / Dxs(i) ) 

        vsu = 0.5*(v(i,j-1,k)   + v(i+1,j-1,k) ) -0.125*Dxs(i)*Dxs(i)/iDx(i+1) &
              *(   (v(i+2,j-1,k) - v(i+1,j-1,k)) / Dxs(i+1) &
                - (v(i+1,j-1,k) - v(i,j-1,k)   ) / Dxs(i) ) 

        wfu = 0.5*(w(i,j,k)     + w(i+1,j,k)  ) -0.125*Dxs(i)*Dxs(i)/iDx(i+1) &
              *(   (w(i+2,j,k)   - w(i+1,j,k) ) / Dxs(i+1) &
                - (w(i+1,j,k)   - w(i,j,k)    ) / Dxs(i)  )

        wbu = 0.5*(w(i,j,k-1)   + w(i+1,j,k-1) ) -0.125*Dxs(i)*Dxs(i)/iDx(i+1) &
              *(   (w(i+2,j,k-1) - w(i+1,j,k-1)) / Dxs(i+1) &
                - (w(i+1,j,k-1) - w(i,j,k-1)   ) / Dxs(i)  )

    end if


    if (v(i,j,k) > 0) then 

        un  = 0.5*(u(i,j,k)     + u(i,j+1,k)  ) -0.125*Dys(j)*Dys(j)/iDy(j) &
              *(   (u(i,j+1,k)   - u(i,j,k)    ) / Dys(j)  &
                - (u(i,j,k)     - u(i,j-1,k)  ) / Dys(j-1)) 

        us  = 0.5*(u(i,j-1,k)   + u(i,j,k)    ) -0.125*Dys(j-1)*Dys(j-1)/iDy(j-1) &
              *(   (u(i,j,k)     - u(i,j-1,k)  ) / Dys(j-1) &
                - (u(i,j-1,k)   - u(i,j-2,k)  ) / Dys(j-2) ) 
            
    else

        un  = 0.5*(u(i,j,k)     + u(i,j+1,k)  )-0.125*Dys(j)*Dys(j)/iDy(j+1) &
              *(   (u(i,j+2,k)   - u(i,j+1,k) ) / Dys(j+1) &
                - (u(i,j+1,k)   - u(i,j,k)    ) / Dys(j) )

        us  = 0.5*(u(i,j-1,k)   + u(i,j,k)    ) -0.125*Dys(j-1)*Dys(j-1)/iDy(j) &
              *(   (u(i,j+1,k)   - u(i,j,k)   ) / Dys(j) &
                - (u(i,j,k)     - u(i,j-1,k)  ) / Dys(j-1) )  

    end if


    if (w(i,j,k) > 0) then 

        uf  = 0.5*(u(i,j,k)     + u(i,j,k+1)  ) -0.125*Dzs(k)*Dzs(k)/iDz(k) &
              *(   (u(i,j,k+1)   - u(i,j,k)   ) / Dzs(k) &
                - (u(i,j,k)     - u(i,j,k-1)  ) / Dzs(k-1) ) !!!

        ub  = 0.5*(u(i,j,k-1)   + u(i,j,k)    ) -0.125*Dzs(k-1)*Dzs(k-1)/iDz(k-1) &
              *(   (u(i,j,k)     - u(i,j,k-1) ) / Dzs(k-1) &
                - (u(i,j,k-1)   - u(i,j,k-2)  ) / Dzs(k-2) ) 
        
    else

        uf  = 0.5*(u(i,j,k)     + u(i,j,k+1)  ) -0.125*Dzs(k)*Dzs(k)/iDz(k+1) &
              *(   (u(i,j,k+2)   - u(i,j,k+1) ) / Dzs(k+1) &
                - (u(i,j,k+1)   - u(i,j,k)    ) / Dzs(k) )

        ub  = 0.5*(u(i,j,k-1)   + u(i,j,k)    ) -0.125*Dzs(k-1)*Dzs(k-1)/iDz(k) &
              *(   (u(i,j,k+1)   - u(i,j,k)   ) / Dzs(k)  &
                - (u(i,j,k)     - u(i,j,k-1)  ) / Dzs(k-1))

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

    u_star(i,j,k) = u(i,j,k)-dt*(ue*ue-uw*uw)   / Dxs(i) &
                            -dt*(un*vnu-us*vsu) / iDy(j) &
                            -dt*(uf*wfu-ub*wbu) / iDz(k) &

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


    if (u(i,j,k) > 0) then 
    
      ve  = 0.5*(v(i,j,k)     + v(i+1,j,k)   ) -0.125*Dxs(i)*Dxs(i)/iDx(i) &
              *(   (v(i+1,j,k)   - v(i,j,k)  ) / Dxs(i) &
                - (v(i,j,k)     - v(i-1,j,k) ) / Dxs(i-1) )

      vw  = 0.5*(v(i-1,j,k)   + v(i,j,k)     ) -0.125*Dxs(i-1)*Dxs(i-1)/iDx(i-1) &
              *(   (v(i,j,k)     - v(i-1,j,k)) / Dxs(i-1) &
                - (v(i-1,j,k)   - v(i-2,j,k) ) / Dxs(i-2) ) 

    else

      ve  = 0.5*(v(i,j,k)     + v(i+1,j,k)    ) -0.125*Dxs(i)*Dxs(i)/iDx(i+1) &
              *(   (v(i+2,j,k)   - v(i+1,j,k) ) / Dxs(i+1) &
                - (v(i+1,j,k)   - v(i,j,k)    ) / Dxs(i) )  

      vw  = 0.5*(v(i-1,j,k)   + v(i,j,k)      ) -0.125*Dxs(i-1)*Dxs(i-1)/iDx(i) &
              *(   (v(i+1,j,k)   - v(i,j,k)   ) / Dxs(i) &
                - (v(i,j,k)     - v(i-1,j,k)  ) / Dxs(i-1) )  

    end if


    if (v(i,j,k) > 0) then 

      uev = 0.5*(u(i,j,k)     + u(i,j+1,k)   ) -0.125*Dys(j)*Dys(j)/iDy(j) &
              *(   (u(i,j+1,k)   - u(i,j,k)  ) / Dys(j) &
                - (u(i,j,k)     - u(i,j-1,k) ) / Dys(j-1) )

      uwv = 0.5*(u(i-1,j,k)   + u(i-1,j+1,k)  ) -0.125*Dys(j)*Dys(j)/iDy(j) &
              *(   (u(i-1,j+1,k) - u(i-1,j,k) ) / Dys(j) &
                - (u(i-1,j,k)   - u(i-1,j-1,k)) / Dys(j-1) )

      vn  = 0.5*(v(i,j,k)     + v(i,j+1,k)   ) -0.125*iDy(j+1)*iDy(j+1)/Dys(j) &
              *(   (v(i,j+1,k)   - v(i,j,k)  ) / iDy(j+1) &
                - (v(i,j,k)     - v(i,j-1,k) ) / iDy(j) ) 

      vs  = 0.5*(v(i,j-1,k)   + v(i,j,k)      ) -0.125*iDy(j)*iDy(j)/Dys(j-1) &
              *(   (v(i,j,k)     - v(i,j-1,k) ) / iDy(j) &
                - (v(i,j-1,k)   - v(i,j-2,k)  ) / iDy(j-1) ) 

      wfv = 0.5*(w(i,j,k)     + w(i,j+1,k)   ) -0.125*Dys(j)*Dys(j)/iDy(j) &
              *(   (w(i,j+1,k)   - w(i,j,k)  ) / Dys(j) &
                - (w(i,j,k)     - w(i,j-1,k) ) / Dys(j-1) )

      wbv = 0.5*(w(i,j,k-1)   + w(i,j+1,k-1)  ) -0.125*Dys(j)*Dys(j)/iDy(j) &
              *(   (w(i,j+1,k-1) - w(i,j,k-1) ) / Dys(j) &
                - (w(i,j,k-1)   - w(i,j-1,k-1)) / Dys(j-1))
            
    else

      uev = 0.5*(u(i,j,k)     + u(i,j+1,k)    ) -0.125*Dys(j)*Dys(j)/iDy(j+1) &
              *(   (u(i,j+2,k)   - u(i,j+1,k) ) / Dys(j+1) &
                - (u(i,j+1,k)   - u(i,j,k)    ) / Dys(j) )

      uwv = 0.5*(u(i-1,j,k)   + u(i-1,j+1,k)   ) -0.125*Dys(j)*Dys(j)/iDy(j+1) &
              *(   (u(i-1,j+2,k) - u(i-1,j+1,k)) / Dys(j+1) &
                - (u(i-1,j+1,k) - u(i-1,j,k)   ) / Dys(j) ) 

      vn  = 0.5*(v(i,j,k)     + v(i,j+1,k)    ) -0.125*iDy(j+1)*iDy(j+1)/Dys(j+1) &
              *(   (v(i,j+2,k)   - v(i,j+1,k) ) / iDy(j+2) &
                - (v(i,j+1,k)   - v(i,j,k)    ) / iDy(j+1) ) 

      vs  = 0.5*(v(i,j-1,k)   + v(i,j,k)      ) -0.125*iDy(j)*iDy(j)/Dys(j) &
              *(   (v(i,j+1,k)   - v(i,j,k)   ) / iDy(j+1) &
                - (v(i,j,k)     - v(i,j-1,k)  ) / iDy(j) ) 

      wfv = 0.5*(w(i,j,k)     + w(i,j+1,k)    ) -0.125*Dys(j)*Dys(j)/iDy(j+1) &
              *(   (w(i,j+2,k)   - w(i,j+1,k) ) / Dys(j+1) &
                - (w(i,j+1,k)   - w(i,j,k)    ) / Dys(j) )

      wbv = 0.5*(w(i,j,k-1)   + w(i,j+1,k-1)   ) -0.125*Dys(j)*Dys(j)/iDy(j+1) &
              *(   (w(i,j+2,k-1) - w(i,j+1,k-1)) / Dys(j+1) &
                - (w(i,j+1,k-1) - w(i,j,k-1)   ) / Dys(j) )

    end if


    if (w(i,j,k) > 0) then 

      vf  = 0.5*(v(i,j,k)     + v(i,j,k+1)    ) -0.125*Dzs(k)*Dzs(k)/iDz(k) &
              *(   (v(i,j,k+1)   - v(i,j,k)   ) / Dzs(k) &
                - (v(i,j,k)     - v(i,j,k-1)  ) / Dzs(k-1) )

      vb  = 0.5*(v(i,j,k-1)   + v(i,j,k)      ) -0.125*Dzs(k-1)*Dzs(k-1)/iDz(k-1) &
              *(   (v(i,j,k)     - v(i,j,k-1) ) / Dzs(k-1) &
                - (v(i,j,k-1)   - v(i,j,k-2)  ) / Dzs(k-2) )
        
    else

      vf  = 0.5*(v(i,j,k)     + v(i,j,k+1)    ) -0.125*Dzs(k)*Dzs(k)/iDz(k+1) &
              *(   (v(i,j,k+2)   - v(i,j,k+1) ) / Dzs(k+1) &
                - (v(i,j,k+1)   - v(i,j,k)    ) / Dzs(k)   )

      vb  = 0.5*(v(i,j,k-1)   + v(i,j,k)      ) -0.125*Dzs(k-1)*Dzs(k-1)/iDz(k) &
              *(   (v(i,j,k+1)   - v(i,j,k)   ) / Dzs(k) &
                - (v(i,j,k)     - v(i,j,k-1)  ) / Dzs(k-1) )

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

    v_star(i,j,k) = v(i,j,k)-dt*(uev*ve-uwv*vw) / iDx(i)  &
                            -dt*(vn*vn-vs*vs) / Dys(j)  &
                            -dt*(vf*wfv-vb*wbv) / iDz(k) &

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

        if (u(i,j,k) > 0) then 
    
          we  = 0.5*(w(i,j,k)     + w(i+1,j,k)  ) -0.125*Dxs(i)*Dxs(i)/iDx(i) &
                  *(   (w(i+1,j,k)   - w(i,j,k) ) / Dxs(i) &
                    - (w(i,j,k)     - w(i-1,j,k)) / Dxs(i-1) ) 

          ww  = 0.5*(w(i-1,j,k)   + w(i,j,k)      ) -0.125*Dxs(i-1)*Dxs(i-1)/iDx(i-1) &
                  *(   (w(i,j,k)     - w(i-1,j,k) ) / Dxs(i-1) &
                    - (w(i-1,j,k)   - w(i-2,j,k)  ) / Dxs(i-2) ) 

        else

          we  = 0.5*(w(i,j,k)     + w(i+1,j,k)    ) -0.125*Dxs(i)*Dxs(i)/iDx(i+1) &
                  *(   (w(i+2,j,k)   - w(i+1,j,k) ) / Dxs(i+1) &
                    - (w(i+1,j,k)   - w(i,j,k)    )  / Dxs(i) )

          ww  = 0.5*(w(i-1,j,k)   + w(i,j,k)     ) -0.125*Dxs(i-1)*Dxs(i-1)/iDx(i) &
                  *(   (w(i+1,j,k)   - w(i,j,k)  ) / Dxs(i) &
                    - (w(i,j,k)     - w(i-1,j,k) ) / Dxs(i-1) ) 

        end if


        if (v(i,j,k) > 0) then 

          wn  = 0.5*(w(i,j,k)     + w(i,j+1,k)   ) -0.125*Dys(j)*Dys(j)/iDy(j) &
                  *(   (w(i,j+1,k)   - w(i,j,k)  ) / Dys(j) &
                    - (w(i,j,k)     - w(i,j-1,k) ) / Dys(j-1) )

          ws  = 0.5*(w(i,j-1,k)   + w(i,j,k)      ) -0.125*Dys(j-1)*Dys(j-1)/iDy(j-1) &
                  *(   (w(i,j,k)     - w(i,j-1,k) ) / Dys(j-1) &
                    - (w(i,j-1,k)   - w(i,j-2,k)  ) / Dys(j-2) )
                
        else

          wn  = 0.5*(w(i,j,k)     + w(i,j+1,k)    ) -0.125*Dys(j)*Dys(j)/iDy(j+1) &
                  *(   (w(i,j+2,k)   - w(i,j+1,k) ) / Dys(j+1) &
                    - (w(i,j+1,k)   - w(i,j,k)    ) / Dys(j) ) 

          ws  = 0.5*(w(i,j-1,k)   + w(i,j,k)     ) -0.125*Dys(j-1)*Dys(j-1)/iDy(j) &
                  *(   (w(i,j+1,k)   - w(i,j,k)  ) / Dys(j) &
                    - (w(i,j,k)     - w(i,j-1,k) ) / Dys(j-1) )

        end if


        if (w(i,j,k) > 0) then 

          uew = 0.5*(u(i,j,k)     + u(i,j,k+1)    ) -0.125*Dzs(k)*Dzs(k)/iDz(k) &
                  *(   (u(i,j,k+1)   - u(i,j,k)   ) / Dzs(k) &
                    - (u(i,j,k)     - u(i,j,k-1)  ) / Dzs(k-1) ) 

          uww = 0.5*(u(i-1,j,k)   + u(i-1,j,k+1)  ) -0.125*Dzs(k)*Dzs(k)/iDz(k) &
                  *(   (u(i-1,j,k+1) - u(i-1,j,k) ) / Dzs(k) &
                    - (u(i-1,j,k)   - u(i-1,j,k-1)) / Dzs(k-1) )

          vnw = 0.5*(v(i,j,k)     + v(i,j,k+1)    ) -0.125*Dzs(k)*Dzs(k)/iDz(k) &
                  *(   (v(i,j,k+1)   - v(i,j,k)   ) / Dzs(k) &
                    - (v(i,j,k)     - v(i,j,k-1)  ) / Dzs(k-1) )

          vsw = 0.5*(v(i,j-1,k)   + v(i,j-1,k+1)  ) -0.125*Dzs(k)*Dzs(k)/iDz(k) &
                  *(   (v(i,j-1,k+1) - v(i,j-1,k) ) / Dzs(k) &
                    - (v(i,j-1,k)   - v(i,j-1,k-1)) / Dzs(k-1) ) 

          wf  = 0.5*(w(i,j,k)     + w(i,j,k+1)    ) -0.125*iDz(k+1)*iDz(k+1)/Dzs(k) &
                  *(   (w(i,j,k+1)   - w(i,j,k)   ) / iDz(k+1) &
                    - (w(i,j,k)     - w(i,j,k-1)  ) / iDz(k) )

           wb  = 0.5*(w(i,j,k-1)   + w(i,j,k)     ) -0.125*iDz(k)*iDz(k)/Dzs(k-1) &
                  *(   (w(i,j,k)     - w(i,j,k-1) ) / iDz(k) &
                    - (w(i,j,k-1)   - w(i,j,k-2)  ) / iDz(k-1) )
            
        else

          uew = 0.5*(u(i,j,k)     + u(i,j,k+1)    ) -0.125*Dzs(k)*Dzs(k)/iDz(k+1) &
                  *(   (u(i,j,k+2)   - u(i,j,k+1) ) / Dzs(k+1) &
                    - (u(i,j,k+1)   - u(i,j,k)    ) / Dzs(k) )

          uww = 0.5*(u(i-1,j,k)   + u(i-1,j,k+1)   ) -0.125*Dzs(k)*Dzs(k)/iDz(k+1) &
                  *(   (u(i-1,j,k+2) - u(i-1,j,k+1)) / Dzs(k+1) &
                    - (u(i-1,j,k+1) - u(i-1,j,k)   ) / Dzs(k) )

          vnw = 0.5*(v(i,j,k)     + v(i,j,k+1)    ) -0.125*Dzs(k)*Dzs(k)/iDz(k+1) &
                  *(   (v(i,j,k+2)   - v(i,j,k+1) ) / Dzs(k+1) &
                    - (v(i,j,k+1)   - v(i,j,k)    ) / Dzs(k) )

          vsw = 0.5*(v(i,j-1,k)   + v(i,j-1,k+1)   ) -0.125*Dzs(k)*Dzs(k)/iDz(k+1) &
                  *(   (v(i,j-1,k+2) - v(i,j-1,k+1)) / Dzs(k+1) &
                    - (v(i,j-1,k+1) - v(i,j-1,k)   ) / Dzs(k) )

          wf  = 0.5*(w(i,j,k)     + w(i,j,k+1)    ) -0.125*iDz(k+1)*iDz(k+1)/Dzs(k+1) &
                  *(   (w(i,j,k+2)   - w(i,j,k+1) ) / iDz(k+2) &
                    - (w(i,j,k+1)   - w(i,j,k)    ) / iDz(k+1)  )

          wb  = 0.5*(w(i,j,k-1)   + w(i,j,k)      ) -0.125*iDz(k)*iDz(k)/Dzs(k) &
                  *(  (w(i,j,k+1)   - w(i,j,k)    ) / iDz(k+1) &
                   - (   w(i,j,k)     - w(i,j,k-1)) / iDz(k)  )

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
      
        w_star(i,j,k) = w(i,j,k)-dt*(uew*we-uww*ww) / iDx(i) &
                                -dt*(vnw*wn-vsw*ws) / iDy(j) &
                                -dt*(wf*wf-wb*wb) / Dzs(k) &

                        +(nu+nut*Q1)*dt*( ( w(i+1,j,k)-w(i,j,k) ) / Dxs(i) -  ( w(i,j,k)-w(i-1,j,k) ) / Dxs(i-1)  ) / iDx(i) &
                        +(nu+nut*Q1)*dt*( ( w(i,j+1,k)-w(i,j,k) ) / Dys(j) -  ( w(i,j,k)-w(i,j-1,k) ) / Dys(j-1)  ) / iDy(j) &
                        +(nu+nut*Q1)*dt*( ( w(i,j,k+1)-w(i,j,k) ) / iDz(k+1) - ( w(i,j,k)-w(i,j,k-1) ) / iDz(k)   ) / Dzs(k)

        end do
    end do
end do
!$OMP END PARALLEL DO

end subroutine discretisation_QUICK