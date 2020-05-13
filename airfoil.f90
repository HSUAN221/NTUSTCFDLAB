subroutine MATRIX_PARSEC(myid,zc,yc,AOA, POINT)
 
    implicit none

    real*8  :: yc 
    real*8  :: zc
    integer, intent(in) :: myid, AOA, POINT

    !---------------------------------------------------!
    !    LOCAL VARIABLES                                !
    !---------------------------------------------------!

    INTEGER,PARAMETER:: N=6,M=6
    DOUBLE PRECISION :: PI   
    DOUBLE PRECISION :: A_UP(N,M),X_UP(N),B_UP(N)
    DOUBLE PRECISION :: A_LO(N,M),X_LO(N),B_LO(N)
    DOUBLE PRECISION :: P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11
    DOUBLE PRECISION :: TEMP_UP,TEMP_LO,ANG_UP,ANG_LO
    DOUBLE PRECISION :: PARTONE,PARTTWO,PIECES
    DOUBLE PRECISION :: Y(1:POINT),X(1:POINT),Z1,Z2
    INTEGER :: I,J 
    CHARACTER(LEN=20) :: POSITIONNAME

    
    


    !!$ Original value NACA 0012
    P1 = 0.0155
    P2 = 0.296632
    P3 = 0.060015
    P4 = -0.4515
    P5 = 0.296632
    P6 = -0.06055
    P7 = 0.453
    P8 = 0.00
    P9 = 0.00126
    P10 = 0.0
    P11 = 7.36




 


    POSITIONNAME='NACA0012.DAT'



    !------------------UPPER FUNCTION-------------------
        A_UP=RESHAPE((/DBLE(1.0),DBLE(1),DBLE(0.5),P2**(0.5),(0.5)*P2**(-0.5),(-0.25)*P2**(-1.5),&
                    DBLE(0.0),DBLE(1),DBLE(1.5),P2**(1.5),(1.5)*P2**(0.5),(0.75)*P2**(-0.5),&
                    DBLE(0.0),DBLE(1),DBLE(2.5),P2**(2.5),(2.5)*P2**(1.5),(3.75)*P2**(0.5),&
                    DBLE(0.0),DBLE(1),DBLE(3.5),P2**(3.5),(3.5)*P2**(2.5),(3.75)*P2**(0.5),&             !MODIFIED THE COEFFICIENT
                    DBLE(0.0),DBLE(1),DBLE(4.5),P2**(4.5),(4.5)*P2**(3.5),(15.75)*P2**(2.5),&
                    DBLE(0.0),DBLE(1),DBLE(5.5),P2**(5.5),(5.5)*P2**(4.5),(24.75)*P2**(3.5)/),(/N,M/))    ! RESHAPING 6*6

        
        PI=4*ATAN(1.0)                                                                              !DEFINE 'PI'
        ANG_UP=(P10-(P11/2.0))*PI/180                                                               !CONVERT THE ANGLE TO ARC
        
        B_UP=(/(2*P1)**0.5,(P8+(P9/2.0)),TAN(ANG_UP),P3,DBLE(0),P4/)     !CONSTANT VALUES           !TAN(A),'A' CANNOT BE CALCULATE IN THE MATRIX
        

    !!$     DO I=1,N
    !!$         !WRITE(*,*)'MATRIX_B_UP',B_UP(I) 
    !!$     END DO
    
    !------------------LOWER FUNCTION-------------------
    !
        A_LO=RESHAPE((/DBLE(1.0),DBLE(1),DBLE(0.5),P5**(0.5),(0.5)*P5**(-0.5),(-0.25)*P5**(-1.5),&
                    DBLE(0.0),DBLE(1),DBLE(1.5),P5**(1.5),(1.5)*P5**(0.5),(0.75)*P5**(-0.5),&
                    DBLE(0.0),DBLE(1),DBLE(2.5),P5**(2.5),(2.5)*P5**(1.5),(3.75)*P5**(0.5),&
                    DBLE(0.0),DBLE(1),DBLE(3.5),P5**(3.5),(3.5)*P5**(2.5),(3.75)*P5**(0.5),&
                    DBLE(0.0),DBLE(1),DBLE(4.5),P5**(4.5),(4.5)*P5**(3.5),(15.75)*P5**(2.5),&
                    DBLE(0.0),DBLE(1),DBLE(5.5),P5**(5.5),(5.5)*P5**(4.5),(24.75)*P5**(3.5)/),(/N,M/))     ! RESHAPING IN TO 6*6          

        
        ANG_LO=(P10+(P11/2.0))*PI/180                                                                !CONVERT THE ANGLE TO ARC
        
        B_LO=(/-(2*P1)**0.5,(P8-(P9/2.0)),TAN(ANG_LO),P6,DBLE(0),P7/)                               !CONSTANT VALUES AND TAN(A),'A' CANNOT BE CALCULATE IN THE MATRIX
        


    !------------------------------------GAUSS ELIMINATION 
    !-----------------------------PART ONE /2ND ROW  -------------------------------------------------------------------------------------    

        DO I=2,N
            TEMP_UP = -(A_UP(I,1)/A_UP(1,1))
            TEMP_LO = -(A_LO(I,1)/A_LO(1,1)) 
            B_UP(I) = B_UP(I)+TEMP_UP*B_UP(1)							!TOTAL MATRIX_B_UP COEFFICIENTS 
            B_LO(I)= B_LO(I)+TEMP_LO*B_LO(1)							!TOTAL MATRIX_B_LO COEFFICIENTS 
            DO J=1,M
                A_UP(I,J)=A_UP(I,J)+TEMP_UP*A_UP(1,J)	                ! FIRST	ROW AND SECOND ROW ADDTION TO CANCEL X1 FROM 2ND  ROW 				
                A_LO(I,J)=A_LO(I,J)+TEMP_LO*A_LO(1,J)                   
                !WRITE(*,*)'A_UP(I,J)',A_UP(I,J)
            END DO
        END DO    
        
        DO I=1,N
            DO J=1,M
                !WRITE(*,*)'MATRIX_A_UP',A_UP(I,J),J
                !WRITE(*,*)'MATRIX_A_LO',A_LO(I,J),J
            END DO 
        END DO

        DO I=1,N
            !WRITE(*,*)'MATRIX_B',B_UP(I) 
        END DO  
    !----------------------PART TWO /3RD ROW --------------------------------------------------------------
    
        DO I=3,N
            TEMP_UP=-(A_UP(I,2)/A_UP(2,2))
            TEMP_LO=-(A_LO(I,2)/A_LO(2,2))        
            B_UP(I)=B_UP(I)+TEMP_UP*B_UP(2)							!TOTAL MATRIX_B_UP COEFFICIENTS 
            B_LO(I)=B_LO(I)+TEMP_LO*B_LO(2)							!TOTAL MATRIX_B_LO COEFFICIENTS 
            DO J=2,M
                A_UP(I,J)=A_UP(I,J)+TEMP_UP*A_UP(2,J)               ! FIRST	ROW AND THIRD ROW ADDTION TO CANCEL X1 FROM 3RD   ROW 	
                A_LO(I,J)=A_LO(I,J)+TEMP_LO*A_LO(2,J)                    
                !WRITE(*,*)'A_UP(I,J)',A_UP(I,J)
            END DO
        END DO   
        
        DO I=1,N
            DO J=1,M
                !WRITE(*,*)'MATRIX_A_UP',A_UP(I,J),J
                !WRITE(*,*)'MATRIX_A_LO',A_LO(I,J),J
            END DO 
        END DO
    !---------------------PART THREE/ 4TH ROW ---------------------------------------------------------------

        DO I=4,N
            TEMP_UP=-(A_UP(I,3)/A_UP(3,3))
            TEMP_LO=-(A_LO(I,3)/A_LO(3,3))
            B_UP(I)=B_UP(I)+TEMP_UP*B_UP(3)							!TOTAL MATRIX_B_UP COEFFICIENTS   
            B_LO(I)=B_LO(I)+TEMP_LO*B_LO(3)							!TOTAL MATRIX_B_LO COEFFICIENTS   
            DO J=3,M
                A_UP(I,J)=A_UP(I,J)+TEMP_UP*A_UP(3,J)
                A_LO(I,J)=A_LO(I,J)+TEMP_LO*A_LO(3,J)                !FIRST	ROW AND FOURTH ROW ADDTION TO CANCEL X1 FROM 4TH ROW 	      
                !WRITE(*,*)'A_UP(I,J)',A_UP(I,J)
            END DO
        END DO
        
        DO I=1,N
            DO J=1,M
                !WRITE(*,*)'MATRIX_A_UP',A_UP(I,J),J
                !WRITE(*,*)'MATRIX_A_LO',A_LO(I,J),J
            END DO 
        END DO
    
        DO I=1,N
            !WRITE(*,*)'MATRIX_B',B_UP(I) 
        END DO     
    !---------------------PART FOUR/5TH ROW ---------------------------------------------------------------

        DO I=5,N
            TEMP_UP=-(A_UP(I,4)/A_UP(4,4))
            TEMP_LO=-(A_LO(I,4)/A_LO(4,4))
            B_UP(I)=B_UP(I)+TEMP_UP*B_UP(4)                           !TOTAL MATRIX_B_UP COEFFICIENTS 
            B_LO(I)=B_LO(I)+TEMP_LO*B_LO(4)                           !TOTAL MATRIX_B_LO COEFFICIENTS 
            DO J=4,M
                A_UP(I,J)=A_UP(I,J)+TEMP_UP*A_UP(4,J)          
                A_LO(I,J)=A_LO(I,J)+TEMP_LO*A_LO(4,J)                !FIRST	ROW AND FITH ROW ADDTION TO CANCEL X1 FROM 5TH ROW 	   
                !WRITE(*,*)'A_UP(I,J)',A_UP(I,J)
            END DO
        END DO
        
        DO I=1,N
            DO J=1,M
                !WRITE(*,*)'MATRIX_A_UP',A_UP(I,J),J
                !WRITE(*,*)'MATRIX_A_LO',A_LO(I,J),J
            END DO 
        END DO


        DO I=1,N
            !WRITE(*,*)'MATRIX_B',B_UP(I) 
        END DO
    !---------------------PART FIVE / 6TH  ROW ---------------------------------------------------------------

        DO I=6,N
            TEMP_UP=-(A_UP(I,5)/A_UP(5,5))
            TEMP_LO=-(A_LO(I,5)/A_LO(5,5))
            B_UP(I)=B_UP(I)+TEMP_UP*B_UP(5)                           !TOTAL MATRIX_B_UP COEFFICIENTS 
            B_LO(I)=B_LO(I)+TEMP_LO*B_LO(5)                           !TOTAL MATRIX_B_LO COEFFICIENTS 
            DO J=5,M
                A_UP(I,J)=A_UP(I,J)+TEMP_UP*A_UP(5,J)
                A_LO(I,J)=A_LO(I,J)+TEMP_LO*A_LO(5,J)                 !FIRST ROW AND SIXTH ROW ADDITION TO CANCEL X1 FROM 6TH ROW
                !WRITE(*,*)'A_UP(I,J)',A_UP(I,J)
            END DO
        END DO
        
        DO I=1,N
            DO J=1,M
                IF(ABS(A_UP(I,J))<0.00000000001 .OR. ABS(A_LO(I,J))<0.00000000001 )THEN               !MAKE SURE EVERY COMPONENT EQUAL ZERO
                    A_UP(I,J)=0.0
                    A_LO(I,J)=0.0
                END IF               
                !WRITE(*,*)'MATRIX_A_UP',A_UP(I,J),J
                !WRITE(*,*)'MATRIX_A_LO',A_LO(I,J),J
            END DO 
        END DO

    !!$     DO I=1,N
    !!$        !WRITE(*,*)'MATRIX_B_UP',B_UP(I)
    !!$        !WRITE(*,*)'MATRIX_B_LO',B_LO(I)
    !!$     END DO

    !--------------------CALCULATE THE ANSWER 'X'/COEFFICIENTS(A1,A2..A6)  FROM  LAST ROW TO FIRST ROW BY SUBSTITUTION -------------------
                                                                            
        X_UP=0.0																!MAKE SURE THE COEFFICIENT 'X'= 0
        X_LO=0.0
    
        DO I=N,1,-1
            X_UP(I)=0                           ! ALL VALUES OF COEFIIECIENT OF X_UP  ARE ZERO AFTER ELIMINATION 
            !WRITE(*,*)'A_UP',A_UP(I,I),'I',I  
            X_UP(I)=(B_UP(I)-A_UP(I,6)*X_UP(6)-A_UP(I,5)*X_UP(5)-A_UP(I,4)*X_UP(4)-A_UP(I,3)*X_UP(3)-&  ! X_UP(I) = a1,a2...a6)
            A_UP(I,2)*X_UP(2)-A_UP(I,1)*X_UP(1))/A_UP(I,I)
            !WRITE(*,*)'X_UP',X_UP(I),'I',I        

            X_LO(I)=0
            X_LO(I)=(B_LO(I)-A_LO(I,6)*X_LO(6)-A_LO(I,5)*X_LO(5)-A_LO(I,4)*X_LO(4)-A_LO(I,3)*X_LO(3)-& !  X_LO(I) = a1,a2...a6)
            A_LO(I,2)*X_LO(2)-A_LO(I,1)*X_LO(1))/A_LO(I,I)
            !WRITE(*,*)'X_LO',X_LO(I),'I',I   
        END DO


    !-------------------------------------------------------------------------------------------------------------
    !                       NORMALIZATION OF X  BETWEEN  0 AND  1 
    !
    !-------------------THE FUNCTION OF UPPER LINE----------------------------------------------------------------
        PARTONE=0.05                                                                        !  !TRAILING AND LEADING EDGEES 
        PARTTWO=0.9				       														! 2*PARTONE + PARTTWO =  1    TOTAL CURVE SIZE  
        !WRITE(*,*)'ONE',PARTONE,'TWO',PARTTWO,'TOTAL',PARTONE+PARTTWO
        OPEN(UNIT=32,FILE=POSITIONNAME)                                                     ! POSITION = 'REWIND')
            !REWIND(32)																		! IN CASE THE FILE LEAVE THE OLD POSITION OF AIRFOIL
            PIECES=POINT/2/3                                                                ! 90 
            Z1=PARTONE/PIECES                                                               ! 0.0005555555
            Z2=PARTTWO/PIECES                                                               ! 0.01
            !WRITE(*,*)'PIECES  ',PIECES,'   Z1',Z1,'   Z2',Z2
                DO J=1,POINT                                                                ! Condition ?
                
    ! FOR LOWER CURVES
                IF (J<=(POINT/2.0))THEN                                                 ! [1 270]
                    IF (J<=(POINT/6.0))THEN                                               ! [1 90]
                            X(J)=J*Z1                                                       !                                         
                            !WRITE(*,*)'X',X,'J',J
                            Y(J)= X_UP(1)*X(J)**(0.5)+X_UP(2)*X(J)**(1.5)+X_UP(3)*X(J)**(2.5) +  X_UP(4)*X(J)**(3.5) + &             ! COEFFICIENTS X_UP(I)
                                X_UP(5)*X(J)**(4.5)+X_UP(6)*X(J)**(5.5)
                                
                        ELSE IF (J>=(POINT/6) .AND. J <= (POINT/3)) THEN                  ![90, 180] 
                            X(J)=PARTONE+(J-(POINT/6))*Z2   
                            !WRITE(*,*)'X',X,'J',J,'  POINT',POINT/6                   
                            Y(J)=X_UP(1)*X(J)**(0.5)+X_UP(2)*X(J)**(1.5)+X_UP(3)*X(J)**(2.5)+X_UP(4)*X(J)**(3.5)+&
                                X_UP(5)*X(J)**(4.5)+X_UP(6)*X(J)**(5.5)
                        ELSE
                            X(J) = PARTONE + PARTTWO + (J-(POINT/3))*Z1                    ![180  270]
                            !WRITE(*,*)'X',X,'J',J,'  POINT',POINT/3                                      
                            Y(J)=X_UP(1)*X(J)**(0.5)+X_UP(2)*X(J)**(1.5)+ X_UP(3)*X(J)**(2.5)+X_UP(4)*X(J)**(3.5)+&
                                X_UP(5)*X(J)**(4.5)+X_UP(6)*X(J)**(5.5)
                        END IF         
    ! FOR LOWER CURVES              
                    ELSE                                                                              
                        X(J)=POINT-J                                                   ! J= [270 540] AND  FOR X(J) 0 ...270
                        IF (X(J)<=(POINT/6.0))THEN                                     ! [0 90]
                            X(J)=X(J)*Z1                  
                            !WRITE(*,*)'X',X,'J',J
                            Y(J)=X_LO(1)*X(J)**(0.5)+X_LO(2)*X(J)**(1.5)+X_LO(3)*X(J)**(2.5)+X_LO(4)*X(J)**(3.5)+&
                                X_LO(5)*X(J)**(4.5)+X_LO(6)*X(J)**(5.5) 
                                
                        ELSE IF (X(J)>(POINT/6) .AND. X(J)<=(POINT/3))THEN              ! [90 180]
                            X(J)=PARTONE+(X(J)-(POINT/6))*Z2                            
                            !WRITE(*,*)'X',X,'J',J,'  POINT',POINT/6,'  PART2'                   
                            Y(J)=X_LO(1)*X(J)**(0.5)+X_LO(2)*X(J)**(1.5)+X_LO(3)*X(J)**(2.5)+X_LO(4)*X(J)**(3.5)+&
                                X_LO(5)*X(J)**(4.5)+X_LO(6)*X(J)**(5.5) 
                        ELSE
                            X(J)=PARTONE+PARTTWO +(X(J)-(POINT/3))*Z1                     ! [180 270]
                            !WRITE(*,*)'X',X,'J',J,'  POINT',POINT/3,'  PART3'                                      
                            Y(J)=X_LO(1)*X(J)**(0.5)+X_LO(2)*X(J)**(1.5)+X_LO(3)*X(J)**(2.5)+X_LO(4)*X(J)**(3.5)+&
                                X_LO(5)*X(J)**(4.5)+X_LO(6)*X(J)**(5.5) 
                        END IF
                    END IF                 
                END DO


            
            if(myid==0)then
            DO  J= 1 , POINT-1
                
            
                !X(J)= COS((-AOA)*PI/180)*X(J) - SIN((-AOA)*PI/180)*Y(J)     ! ROTATION
                !Y(J)= SIN((-AOA)*PI/180)*X(J) + COS((-AOA)*PI/180)*Y(J) 
                !X(J)=X(J) + zc					                    
                !Y(J)=Y(J) + yc
                
                WRITE(32,*)  X(J), Y(J)
            END DO  

            write(32,*) X(1), Y(1)

            CLOSE(32)   
            end if

end subroutine  MATRIX_PARSEC





subroutine dynamic_AOA() 
    use variables
    implicit none
    real*8 ,parameter :: PI = 3.14159265359

    frequency = reduce_frequency / PI
    AOA = AOA1 + 10 * abs( sin( 2 * PI * frequency * ( time - dt*(StartDynamic-1) ) ) )

   

end subroutine dynamic_AOA



subroutine RayCasting_HSUAN()
    use variables
    implicit none
    integer :: l ,m ,n, intercount
    integer :: ll ,mm ,nn
    real*8,dimension(nSubGrids+1) :: ZZ, YY
    real*8 :: total
    real*8 :: diagonal
    DOUBLE PRECISION :: PI


    real*8 ,dimension(1:nSubGrids+1) :: SY
    real*8 ,dimension(1:nSubGrids+1) :: SZ
    real*8 :: dxg, dyg, dzg
    real*8 :: xi

    
    open( 24,file = NACA_filename, form = 'FORMATTED' )
        do i = 1, poly
            read(24,*) az(i), ay(i)
        end do
    close(24)

    PI=4.d0*ATAN(1.d0)

    
    !$OMP PARALLEL DO
    do i = 1, poly
        az(i) = az(i) - 0.25   

        az(i)= COS((-AOA)*PI/180)*az(i) - SIN((-AOA)*PI/180)*ay(i)     ! ROTATION
        ay(i)= SIN((-AOA)*PI/180)*az(i) + COS((-AOA)*PI/180)*ay(i) 

        az(i) = az(i) + zc
        ay(i) = ay(i) + yc

        az(i) = az(i) + 0.25
    end do

    !GET THE AIRFOIL BORDERS BY COMPARING  MESH POINTS WITHIN DOMAIN   
    do j=1,ny
        if( Y(j) > minval(ay) )then
            A = j-1
            !if(myid==0)then
            !    write(*,*) 'A', A, minval(ay)
            !end if
            exit
        end if
    end do

    do j=1,ny
        if( Y(j) > maxval(ay) )then
            B = j
            !if(myid==0)then
            !    write(*,*) 'B', B, maxval(ay)
            !end if
            exit
        end if
    end do

    do k=1,nz
        if( Z(k) > minval(az) )then
            C = k-1
            !if(myid==0)then
            !    write(*,*) 'C', C, minval(az)
            !end if
            exit
        end if
    end do


    do k=1,nz
        if( Z(k) > maxval(az) )then
            E = k
            !if(myid==0)then
            !    write(*,*) 'E', E, maxval(az)
            !end if
            exit
        end if
    end do    

   
    !How many points on edge 
    points = 0
    ik = 0
    do m =1,poly
    !-------------------------!
        do k=1,nz
            if( az(m) >= Z(k) .AND. az(m) < Z(k+1) )then
                do j=1,ny
                    if( ay(m) >=Y(j) .AND. ay(m) < Y(j+1) )then


                            points(K,j) = 1
                            ETA_1(k,j) = points(k,j)
                            ik = ik + 1


                    end if
                end do
            end if
        end do
    !-------------------------!
    end do 

    intersection = 0
    do k=C,E-1
        do j=A,B-1

            intercount = 0
            do l=j,B-1
                if(ETA_1(k,l) /= 0)then
                    intercount = intercount + 1
                end if
            end do


            if(intercount>0 .AND. ETA_1(k,j-1)==1 .AND. ETA_1(k,j)==0)then
                ETA_1(k,j) = 1
            end if


        end do
    end do

    !-----------------------------------------------!
    !              DEFINE THE SUBGRID               !
    !-----------------------------------------------!
    
    
    !-------------------!
    do k=C,E-1   
        do j=A,B-1      
    !-------------------!


            if( points(k,j) > 0 )then
                    sub_intersection = 0
                    position = 0
                    do l = 1,nSubGrids+1
                        position = position  + 1
                        ZZ(l) = Z(k) + (position-1) * ( Z(k+1)-Z(k) ) / DBLE(nSubGrids) !SUB_GRID SIZE IN X
                        YY(l) = Y(j) + (position-1) * ( Y(j+1)-Y(j) ) / DBLE(nSubGrids) !SUB_GRID SIZE IN Y
                    end do

                    !--------------CALCULATION THE SUBGRID-----------------!
                    !$OMP PARALLEL DO PRIVATE(l)  
                    !-------------------!
                    do m=1,poly-1       !
                        do l=1,nSubGrids!
                    !-------------------!

                        !---------- AY(I) > AY(I+1)  ----------!
                        if( ay(m) >=  (YY(l)+YY(l+1))/2.d0  )then
                            if( ay(m+1) <= (YY(l)+YY(l+1))/2.d0 )then
                                !-------------------------!
                                do n=1, nSubGrids
                                    if( (ZZ(n)+ZZ(n+1))/2.d0 > max(az(m),az(m+1)) )then
                                        sub_intersection(n,l) = sub_intersection(n,l) + 0
                                    else if( (ZZ(n)+ZZ(n+1))/2.d0 < min(az(m),az(m+1)) )then
                                        sub_intersection(n,l) = sub_intersection(n,l) + 1
                                    else
                                        if( YY(l)==ay(m+1) .OR. ZZ(n)==az(m+1) )then
                                            YY(l) = YY(l) + 0.d01
                                            ZZ(n) = ZZ(n) + 0.d01
                                        end if
                                        m_pa = ( YY(l)-ay(m+1) ) / ( ZZ(n)-az(m+1) ) 
                                        m_ab = ( ay(m)-ay(m+1) ) / ( az(m)-az(m+1) )

                                        if( m_pa > m_ab )then
                                            sub_intersection(n,l) = sub_intersection(n,l) + 1
                                        else
                                            sub_intersection(n,l) = sub_intersection(n,l) + 0
                                        end if
                                    end if
                                end do
                                !-------------------------!
                            end if
                        end if
                        !---------- AY(I) > AY(I+1)  ----------!


                        !---------- AY(I) < AY(I+1)  ----------!
                        if( ay(m) <=  (YY(l)+YY(l+1))/2.d0  )then
                            if( ay(m+1) >= (YY(l)+YY(l+1))/2.d0 )then
                                !-------------------------!
                                do n=1, nSubGrids
                                    if( (ZZ(n)+ZZ(n+1))/2.d0 > max(az(m),az(m+1)) )then
                                        sub_intersection(n,l) = sub_intersection(n,l) + 0
                                    else if( (ZZ(n)+ZZ(n+1))/2.d0 < min(az(m),az(m+1)) )then
                                        sub_intersection(n,l) = sub_intersection(n,l) + 1
                                    else
                                        if( YY(l)==ay(m+1) .OR. ZZ(n)==az(m+1) )then
                                            YY(l) = YY(l) + 0.d01
                                            ZZ(n) = ZZ(n) + 0.d01
                                        end if
                                        m_pa = ( YY(l)-ay(m) ) / ( ZZ(n)-az(m) ) 
                                        m_ab = ( ay(m+1)-ay(m) ) / ( az(m+1)-az(m) )

                                        if( m_pa > m_ab )then
                                            sub_intersection(n,l) = sub_intersection(n,l) + 1
                                        else
                                            sub_intersection(n,l) = sub_intersection(n,l) + 0
                                        end if
                                    end if
                                end do
                                !-------------------------!
                            end if
                        end if
                        !---------- AY(I) < AY(I+1)  ----------!

                    !---------!
                        end do!
                    end do    !
                    !---------!
                    !$OMP END PARALLEL DO
                    !--------------CALCULATION THE SUBGRID-----------------!

                total = 0
                
                do m=1,nSubGrids
                    do n=1,nSubGrids

                        if( mod(sub_intersection(m,n),2) == 0 )then
                            sub_intersection(m,n) = 0
                        else
                            sub_intersection(m,n) = 1
                        end if
                        total = total + sub_intersection(m,n)

                    end do
                end do
                 

                ETA_1(k,j) = DBLE(total)/DBLE(nSubGrids)/DBLE(nSubGrids)
        
            end if


    !---------!
        end do!
    end do    !
    !---------!
    

    !$OMP PARALLEL DO PRIVATE(i)  
    do k=1,nz ;do j=1,ny; do i=1,nx
        ETA(i,j,k) = ETA_1(k,j)      
    end do; end do; end do
    !$OMP END PARALLEL DO

endsubroutine RayCasting_HSUAN



subroutine RayCasting() 
    use variables
    implicit none
    integer :: l ,m ,n
    real*8,dimension(nSubGrids+1) :: ZZ, YY
    real*8 :: total
    DOUBLE PRECISION :: PI
    
    
    intersection = 0
    open( 24,file = NACA_filename, form = 'FORMATTED' )
        do i = 1, poly
            read(24,*) az(i), ay(i)
        end do
    close(24)


    PI=4.d0*ATAN(1.d0)

    
    !$OMP PARALLEL DO
    do i = 1, poly
        az(i) = az(i) - 0.25   

        az(i)= COS((-AOA)*PI/180)*az(i) - SIN((-AOA)*PI/180)*ay(i)     ! ROTATION
        ay(i)= SIN((-AOA)*PI/180)*az(i) + COS((-AOA)*PI/180)*ay(i) 

        az(i) = az(i) + zc
        ay(i) = ay(i) + yc

        az(i) = az(i) + 0.25
    end do
    !$OMP END PARALLEL DO

    !GET THE AIRFOIL BORDERS BY COMPARING  MESH POINTS WITHIN DOMAIN   
    do j=1,ny
        if( Y(j) > minval(ay) )then
            A = j-1
            !if(myid==0)then
            !    write(*,*) 'A', A, minval(ay)
            !end if
            exit
        end if
    end do

    do j=1,ny
        if( Y(j) > maxval(ay) )then
            B = j
            !if(myid==0)then
            !    write(*,*) 'B', B, maxval(ay)
            !end if
            exit
        end if
    end do

    do k=1,nz
        if( Z(k) > minval(az) )then
            C = k-1
            !if(myid==0)then
            !    write(*,*) 'C', C, minval(az)
            !end if
            exit
        end if
    end do


    do k=1,nz
        if( Z(k) > maxval(az) )then
            E = k
            !if(myid==0)then
            !    write(*,*) 'E', E, maxval(az)
            !end if
            exit
        end if
    end do    



    !How many points on edge 
    points = 0
    ik = 0
    do m =1,poly
    !-------------------------!
        do k=C,E-1
            if( az(m) >= Z(k) .AND. az(m) < Z(k+1) )then
                do j=A,B-1
                    if( ay(m) >=Y(j) .AND. ay(m) < Y(j+1) )then
                        points(K,j) = 1
                        ik = ik + 1
                    end if
                end do
            end if
        end do
    !-------------------------!
    end do 


    

    
    !-----------------------------------------------!
    !   Define ETA which point inside the airfoil   !
    !-----------------------------------------------!
    
    !--------------!
    do m=1,poly-1  !
        do j=1,ny  !
    !--------------!


            !---------- AY(I) > AY(I+1)  ----------!
            if( ay(m) >= (y(j)+y(j+1))/2.d0 )then
                if( ay(m+1) <= (y(j)+y(j+1))/2.d0 )then
                !-------------------------!
                    do k=1,nz
                        if( (Z(k)+Z(k+1))/2.d0 >= max( az(m),az(m+1) ) )then
                            intersection(k,j) = intersection(k,j) + 0
                        else if( (Z(k)+Z(k+1))/2.d0 < min( az(m),az(m+1) ) )then
                            intersection(k,j) = intersection(k,j) + 1
                        else
                            if( Y(j)==ay(m+1) .OR. Z(k)==az(m+1) )then
                                Y(j) = Y(j) + 0.d01
                                Z(k) = Z(k) + 0.d01
                            end if
                            m_pa = ( Y(j)-ay(m+1) ) / ( Z(k)-az(m+1) )
                            m_ab = ( ay(m)-ay(m+1) ) / ( az(m)-az(m+1) )

                            if( m_pa > m_ab )then
                                intersection(k,j) = intersection(k,j) + 1
                            else
                                intersection(k,j) = intersection(k,j) + 0
                            end if
                        end if
                    end do
                !-------------------------!
                end if
            end if
            !---------- AY(I) > AY(I+1)  ----------!
                


            !---------- AY(I) < AY(I+1)  ----------!
            if( ay(m) <= (y(j)+y(j+1))/2.d0 )then
                if( ay(m+1) >= (y(j)+y(j+1))/2.d0 )then
                !-------------------------!
                    do k=1,nz
                        if( (Z(k)+Z(k+1))/2.d0 >= max( az(m),az(m+1) ) )then
                            intersection(k,j) = intersection(k,j) + 0
                        else if( (Z(k)+Z(k+1))/2.d0 < min( az(m),az(m+1) ) )then
                            intersection(k,j) = intersection(k,j) + 1
                        else
                            if( Y(j)==ay(m) .OR. Z(k)==az(m) )then
                                Y(j) = Y(j) + 0.d01
                                Z(k) = Z(k) + 0.d01
                            end if
                            m_pa = ( Y(j)-ay(m) ) / ( Z(k)-az(m) )
                            m_ab = ( ay(m+1)-ay(m) ) / ( az(m+1)-az(m) )

                            if( m_pa > m_ab )then
                                intersection(k,j) = intersection(k,j) + 1
                            else
                                intersection(k,j) = intersection(k,j) + 0
                            end if
                        end if
                    end do
                !-------------------------!
                end if
            end if
            !---------- AY(I) < AY(I+1)  ----------!

    !---------!
        end do!
    end do    !
    !---------!

    !$OMP PARALLEL DO PRIVATE(j) 
    do k=1,nz;do j=1,ny

            if( mod( intersection(k,j),2 ) == 0 )then
                ETA_1(k,j) = 0
            else
                ETA_1(k,j) = 1
            end if
    
    end do;end do    
    !$OMP END PARALLEL DO

    
    !-----------------------------------------------!
    !              DEFINE THE SUBGRID               !
    !-----------------------------------------------!
    
    !-------------!
    do k=C,E-1    !
        do j=A,B-1!
    !-------------!


            if( points(k,j) > 0 )then
                    sub_intersection = 0
                    position = 0
                    do l = 1,nSubGrids+1
                        position = position  + 1
                        ZZ(l) = Z(k) + (position-1) * ( Z(k+1)-Z(k) ) / DBLE(nSubGrids) !SUB_GRID SIZE IN X
                        YY(l) = Y(j) + (position-1) * ( Y(j+1)-Y(j) ) / DBLE(nSubGrids) !SUB_GRID SIZE IN Y
                    end do

                    !--------------CALCULATION THE SUBGRID-----------------!
                    !$OMP PARALLEL DO PRIVATE(l)  
                    !-------------------!
                    do m=1,poly-1       !
                        do l=1,nSubGrids!
                    !-------------------!

                        !---------- AY(I) > AY(I+1)  ----------!
                        if( ay(m) >=  (YY(l)+YY(l+1))/2.d0  )then
                            if( ay(m+1) <= (YY(l)+YY(l+1))/2.d0 )then
                                !-------------------------!
                                do n=1, nSubGrids
                                    if( (ZZ(n)+ZZ(n+1))/2.d0 > max(az(m),az(m+1)) )then
                                        sub_intersection(n,l) = sub_intersection(n,l) + 0
                                    else if( (ZZ(n)+ZZ(n+1))/2.d0 < min(az(m),az(m+1)) )then
                                        sub_intersection(n,l) = sub_intersection(n,l) + 1
                                    else
                                        if( YY(l)==ay(m+1) .OR. ZZ(n)==az(m+1) )then
                                            YY(l) = YY(l) + 0.d01
                                            ZZ(n) = ZZ(n) + 0.d01
                                        end if
                                        m_pa = ( YY(l)-ay(m+1) ) / ( ZZ(n)-az(m+1) ) 
                                        m_ab = ( ay(m)-ay(m+1) ) / ( az(m)-az(m+1) )

                                        if( m_pa > m_ab )then
                                            sub_intersection(n,l) = sub_intersection(n,l) + 1
                                        else
                                            sub_intersection(n,l) = sub_intersection(n,l) + 0
                                        end if
                                    end if
                                end do
                                !-------------------------!
                            end if
                        end if
                        !---------- AY(I) > AY(I+1)  ----------!


                        !---------- AY(I) < AY(I+1)  ----------!
                        if( ay(m) <=  (YY(l)+YY(l+1))/2.d0  )then
                            if( ay(m+1) >= (YY(l)+YY(l+1))/2.d0 )then
                                !-------------------------!
                                do n=1, nSubGrids
                                    if( (ZZ(n)+ZZ(n+1))/2.d0 > max(az(m),az(m+1)) )then
                                        sub_intersection(n,l) = sub_intersection(n,l) + 0
                                    else if( (ZZ(n)+ZZ(n+1))/2.d0 < min(az(m),az(m+1)) )then
                                        sub_intersection(n,l) = sub_intersection(n,l) + 1
                                    else
                                        if( YY(l)==ay(m+1) .OR. ZZ(n)==az(m+1) )then
                                            YY(l) = YY(l) + 0.d01
                                            ZZ(n) = ZZ(n) + 0.d01
                                        end if
                                        m_pa = ( YY(l)-ay(m) ) / ( ZZ(n)-az(m) ) 
                                        m_ab = ( ay(m+1)-ay(m) ) / ( az(m+1)-az(m) )

                                        if( m_pa > m_ab )then
                                            sub_intersection(n,l) = sub_intersection(n,l) + 1
                                        else
                                            sub_intersection(n,l) = sub_intersection(n,l) + 0
                                        end if
                                    end if
                                end do
                                !-------------------------!
                            end if
                        end if
                        !---------- AY(I) < AY(I+1)  ----------!

                    !---------!
                        end do!
                    end do    !
                    !---------!
                    !$OMP END PARALLEL DO
                    !--------------CALCULATION THE SUBGRID-----------------!

                total = 0
                
                do m=1,nSubGrids
                    do n=1,nSubGrids

                        if( mod(sub_intersection(m,n),2) == 0 )then
                            sub_intersection(m,n) = 0
                        else
                            sub_intersection(m,n) = 1
                        end if
                        total = total + sub_intersection(m,n)

                    end do
                end do
                 

                ETA_1(k,j) = DBLE(total)/DBLE(nSubGrids)/DBLE(nSubGrids)
        
            end if


    !---------!
        end do!
    end do    !
    !---------!
    


    !$OMP PARALLEL DO PRIVATE(i)  
    do k=1,nz ;do j=1,ny; do i=1,nx
        ETA(i,j,k) = ETA_1(k,j)      
    end do; end do; end do
    !$OMP END PARALLEL DO


    

    !if(myid==0)then
    !    open(UNIT=18,FILE='test.dat')
    !    write(18,*)'VARIABLES=Z,Y,point,ETA,intersection'
    !    write(18,*)'ZONE k=',nz,' ,j=',ny
    ! 
    !    do  k=1,nz
    !        do  j=1,ny
    !  
    !            write(18,*) Zs(k), Ys(j), points(K,j), ETA_1(k,j), intersection(k,j)
    !        
    !        end do
    !    end do
    !    close(18)
    !   
    !end if
    

    return       
end subroutine RayCasting



