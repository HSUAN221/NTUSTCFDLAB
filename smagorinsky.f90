subroutine CalculateSmagorinskyViscosity()
    use variables 
    implicit none


    call GradientPhiGauss(u,dudx,Dxs,iDy,iDz)
    call GradientPhiGauss(v,dvdx,iDx,Dys,iDz)
    call GradientPhiGauss(w,dwdx,iDx,iDy,Dzs)

    do k=istart,iend
    !$OMP PARALLEL DO PRIVATE(i)  
    do j=1,ny; do i=1,nx

    mutsgs = abs(&
                 2*dudx(i,j,k,1)**2 + 2*dvdx(i,j,k,2)**2 + 2*dwdx(i,j,k,3)**2 + &
                 ( dudx(i,j,k,2)+dvdx(i,j,k,1) )**2 + &
                 ( dudx(i,j,k,3)+dwdx(i,j,k,1) )**2 + &
                 ( dvdx(i,j,k,3)+dwdx(i,j,k,2) )**2         )

    delta = ( iDx(i)*iDy(j)*iDz(k) ) ** (1./3.)             

    mutsgs = (Cs*delta)**2*sqrt(mutsgs)    

    Viseff(i,j,k) = mutsgs


    


    end do;end do
    !$OMP END PARALLEL DO
    end do
    



end subroutine CalculateSmagorinskyViscosity




subroutine GradientPhiGauss(Phi,dPhidX,Dxm,Dym,Dzm)
    use variables 
    implicit none

    real*8, dimension(nx,ny,nz) :: Phi
    real*8, dimension(nx,ny,nz,3) :: dPhidX

    !--------Local variables--------!
    real*8 :: fact
    real*8 ,dimension (-1:nx+2) :: Dxm, Dym, Dzm 
    real*8 :: Phiface_w, Phiface_e, Phiface_n, Phiface_s, Phiface_b, Phiface_f
    !--------Local variables--------!


    do k=istart,iend
    !$OMP PARALLEL DO PRIVATE(i)  
    do j=1,ny; do i=1,nx

        Phiface_e = 0.5 * ( Phi(i+1,j,k) + Phi(i,j,k) )
        Phiface_w = 0.5 * ( Phi(i-1,j,k) + Phi(i,j,k) )

        Phiface_n = 0.5 * ( Phi(i,j+1,k) + Phi(i,j,k) )
        Phiface_s = 0.5 * ( Phi(i,j-1,k) + Phi(i,j,k) )

        Phiface_f = 0.5 * ( Phi(i,j,k+1) + Phi(i,j,k) )
        Phiface_b = 0.5 * ( Phi(i,j,k-1) + Phi(i,j,k) )

        dPhidX(i,j,k,1) = Phiface_w*(-1)*(Dym(j)*Dzm(k)) + Phiface_e*(1)*(Dym(j)*Dzm(k)) 

        dPhidX(i,j,k,2) = Phiface_s*(-1)*(Dxm(i)*Dzm(k)) + Phiface_n*(1)*(Dxm(i)*Dzm(k)) 

        dPhidX(i,j,k,3) = Phiface_b*(-1)*(Dxm(i)*Dym(j)) + Phiface_f*(1)*(Dxm(i)*Dym(j))

        fact = 1.0/(Dxm(i)*Dym(j)*Dzm(k))

        dPhidX(i,j,k,1) =  fact * dPhidX(i,j,k,1) 
        dPhidX(i,j,k,2) =  fact * dPhidX(i,j,k,2) 
        dPhidX(i,j,k,3) =  fact * dPhidX(i,j,k,3) 


    end do;end do
    !$OMP END PARALLEL DO
    end do


end subroutine GradientPhiGauss