subroutine plasma(Zs, Ys, nx, ny, nz, F_tavex, F_tavey, edelta, E, myid, master, az, ay, poly, AOA)
    implicit none

!           Shyy et. al (2002)
!           y
            ! \
            !  
            !    a
            !  
!===========! /
            !------------------------------------------> x
!           \      =================/
!                        b

    !---------------------------------------------------!
    !    VARIABLES IN CALL                              !
    !---------------------------------------------------!
    integer, intent(in) :: poly
    real*8, intent(in) :: AOA
    real*8, dimension(1:nz), intent(in) :: Zs
    real*8, dimension(1:ny), intent(in) :: Ys
    real*8, dimension(1:poly), intent(in) :: az
    real*8, dimension(1:poly), intent(in) :: ay
    real*8, dimension(1:nx,1:ny,1:nz), intent(inout) :: F_tavex, F_tavey
    real*8, dimension(1:nx,1:ny,1:nz), intent(inout) :: edelta, E
    integer, intent(in) :: nx, ny, nz
    integer, intent(in) :: myid, master
    
    

    !---------------------------------------------------!
    !    LOCAL VARIABLES                                !
    !---------------------------------------------------!
    real*8 :: PlasmaZc
    real*8 :: PlasmaYc
    real*8 :: PlasmaZr
    real*8 :: PlasmaYr
    real*8 :: PlasmaZu
    real*8 :: PlasmaYu
    real*8 :: PI
    real*8 :: referpointZ
    real*8 :: referpointY
    real*8 :: ror

    real*8 :: E0, k1, k2
    real*8, dimension(1:nx,1:ny,1:nz) :: Ex, Ey, fx, fy
    real*8, dimension(1:nx,1:ny,1:nz) :: f_effx, f_effy
    integer :: i, j, k
 
    real*8 :: nu = 3500             !
    real*8 :: roc = 5e12
    real*8 :: poto = 7000.0           !
    real*8 :: Eb = 3000.0
    real*8 :: delta_t = 67.d0       !
    real*8 :: len_a = 0.025         !
    real*8 :: len_b = 0.05          !
    real*8 :: len_d = 0.00125       !
    real*8 :: alfa = 1.0
    real*8 :: ec = 1.6e-19


   

    PI=4.d0*ATAN(1.d0)

    !initialize zero point of Plasma
    PlasmaZc = az(591) ! U can change data point
    PlasmaYc = ay(591) ! U can change data point

    PlasmaZr = az(656)
    PlasmaYr = ay(656)


    ! shift point to zero point
    PlasmaZr = PlasmaZr - PlasmaZc
    PlasmaYr = PlasmaYr - PlasmaYc

    PlasmaZc = 0
    PlasmaYc = 0

    
   
    
    !initialize other points of Plasma
    PlasmaZu = sqrt( len_a**2 / ( 1 + (PlasmaZr**2/PlasmaYr**2) ) )
    
    if(PlasmaYr > PlasmaYc)then
        PlasmaZu = PlasmaZu * (-1)
    end if

    PlasmaYu = - ( PlasmaZr / PlasmaYr ) * PlasmaZu
    
    ! shift point go back
    PlasmaZc = az(591) ! U can change data point
    PlasmaYc = ay(591) ! U can change data point

    PlasmaZr = PlasmaZr + PlasmaZc
    PlasmaYr = PlasmaYr + PlasmaYc

    PlasmaZu = PlasmaZu + PlasmaZc
    PlasmaYu = PlasmaYu + PlasmaYc

    
   

    


    !initialize parameter of Plasma
    do k=1,nz; do j=1,ny; do i=1,nx
        E(i,j,k) = 1.D0
        Ex(i,j,k) = 0.D0
        Ey(i,j,k) = 0.D0
        fx(i,j,k) = 0.D0
        fy(i,j,k) = 0.D0
        edelta(i,j,k) = 0.D0
        f_effx(i,j,k) = 0.D0
        f_effy(i,j,k) = 0.D0
        F_tavex(i,j,k) = 0.D0
        F_tavey(i,j,k) = 0.D0
    end do; end do; end do

    

    E0 = poto / len_d

    k1 = ( E0 - Eb ) / len_b
    k2 = ( E0 - Eb ) / len_a

  
    i = 1


    do k=1,nz; do j=1,ny



        if(PlasmaYr > PlasmaYc)then

                    if( ( (Zs(k)-PlasmaZu)/(PlasmaZc-PlasmaZu) - (Ys(j)-PlasmaYu)/(PlasmaYc-PlasmaYu) ) >= 0 ) then  ! III   
                        if( (  (Zs(k)-PlasmaZc)/(PlasmaZr-PlasmaZc) - (Ys(j)-PlasmaYc)/(PlasmaYr-PlasmaYc)  ) <= 0  ) then   ! II 
                            if( (  (Zs(k)-PlasmaZu)/(PlasmaZr-PlasmaZu) - (Ys(j)-PlasmaYu)/(PlasmaYr-PlasmaYu)  ) <= 0 ) then ! I

                                E(i,j,k) = ABS( E0 - k1*(Zs(k)-PlasmaZc) - k2*(Ys(j)-PlasmaYc) )
                                edelta(i,j,k) = 1
                            
                            else

                                edelta(i,j,k) = 0
                        
                            end if
                        end if
                    end if

        else

                    if( ( (Zs(k)-PlasmaZu)/(PlasmaZc-PlasmaZu) - (Ys(j)-PlasmaYu)/(PlasmaYc-PlasmaYu) ) <= 0 ) then  ! III   
                        if( (  (Zs(k)-PlasmaZc)/(PlasmaZr-PlasmaZc) - (Ys(j)-PlasmaYc)/(PlasmaYr-PlasmaYc)  ) >= 0  ) then   ! II 
                            if( (  (Zs(k)-PlasmaZu)/(PlasmaZr-PlasmaZu) - (Ys(j)-PlasmaYu)/(PlasmaYr-PlasmaYu)  ) <= 0 ) then ! I

                                E(i,j,k) = ABS( E0 - k1*(Zs(k)-PlasmaZc) - k2*(Ys(j)-PlasmaYc) )
                                edelta(i,j,k) = 1
                            
                            else

                                edelta(i,j,k) = 0
                        
                            end if
                        end if
                    end if   

        end if

        

 

        Ex(i,j,k) = ( E(i,j,k)*k2 ) / ( k1**2 + k2**2 ) ** 0.5
        Ey(i,j,k) = ( E(i,j,k)*k1 ) / ( k1**2 + k2**2 ) ** 0.5

        fx(i,j,k) = Ex(i,j,k) * roc * ec
        fy(i,j,k) = Ey(i,j,k) * roc * ec

        f_effx(i,j,k) = alfa * fx(i,j,k) * edelta(i,j,k)
        f_effy(i,j,k) = alfa * fy(i,j,k) * edelta(i,j,k)

        F_tavex(i,j,k) = nu * f_effx(i,j,k) * delta_t
        F_tavey(i,j,k) = nu * f_effy(i,j,k) * delta_t



        !if(edelta(i,j,k)==1 .AND. myid==master)then
        !    write(*,*) F_tavex(i,j,k), F_tavey(i,j,k), E(i,j,k)
        !end if
        


    end do; end do

    !if(myid==0)then
    !    write(*,*) PlasmaZc, PlasmaYc
    !    write(*,*) PlasmaZr, PlasmaYr
    !    write(*,*) PlasmaZu, PlasmaYu
    !end if
    

    do k=1,nz ;do j=1,ny; do i=1,nx
        edelta(i,j,k)  = edelta(1,j,k)     
        F_tavex(i,j,k) = F_tavex(1,j,k)     
        F_tavey(i,j,k) = F_tavey(1,j,k)     
    end do; end do; end do


    
    return

end subroutine plasma