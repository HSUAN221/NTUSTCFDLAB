SUBROUTINE solver_BICG3(coef,jcoef,b,x,n,m,eps,p,r,r2,v,ss,t,itmax,myid,iDIMstart,iDIMend,DIMcount)

!-----------------------------------------------------------------------!
!                                                                       !
! 	  BICONJUGUATE GRADIENT FILL-IN 2                                   !
!                                                                       !
!     ATTENTION : coef, jcoef et b are modified in this routine         !
!     ---------                                                         !
!                                                                       !
!-----------------------------------------------------------------------!
	use omp_lib
	implicit none
	
!---------------------------------------------------!
!    VARIABLES IN CALL                              !
!---------------------------------------------------!

	integer, intent(in) :: n,m,itmax
	double precision, intent(in) :: eps
	double precision :: starttime

	integer, dimension(m), intent(inout) :: jcoef
	double precision, dimension(n,m), intent (inout) :: coef
	
	double precision, dimension(n), intent (inout) :: b,x
	double precision, dimension(n), intent(inout) :: p,r,r2 ,v,ss,t

!---------------------------------------------------!
!    VARIABLES IN CALL OF MPI                       !
!---------------------------------------------------!
	integer :: myid, iDIMstart, iDIMend
	integer ,dimension(1:40) :: DIMcount 

!---------------------------------------------------!
!    LOCAL VARIABLES                                !
!---------------------------------------------------!
 
	double precision :: alpha,beta,nu,mu,norm0,norm,sum,scal,norm1,norm2 ,omega,rho1,rho2
	integer :: i,j,col,  nx,ny,nz,ik,k
	


	
 	


		!iDIMstart, iDIMend


!---------------------------------------------------!
!    BI CONJUGUATE GRADIENT                         !
!---------------------------------------------------!
	norm0=norm2(b,n)


	call matmul_ell(p,coef,jcoef,x,n,m)

	!$OMP PARALLEL DO 
	do i=1,n; r(i)=b(i)-p(i); end do
	!$OMP END PARALLEL DO

	!$OMP PARALLEL DO 
	do i=1,n; r2(i) = r(i) ; end do
	!$OMP END PARALLEL DO

    rho1 = 1
	alpha = 1
	omega = 1

	!$OMP PARALLEL DO 
	do i=1,n; v(i) = 0; end do
	!$OMP END PARALLEL DO

	!$OMP PARALLEL DO 
    do i=1,n; p(i) = 0; end do
	!$OMP END PARALLEL DO
	
	norm=0.
	
	!$OMP PARALLEL DO 
	do i=1,n; norm=norm+r(i)*r(i); end do
	!$OMP END PARALLEL DO

	norm=sqrt(norm)/norm0

	j=0

	do while (norm>eps.and.j.lt.itmax)
		j=j+1

		rho2 = scal(r2,r,n)
		
		beta = (rho2/rho1) * (alpha/omega)

		!$OMP PARALLEL DO 
		do i=1,n; p(i) = r(i) + beta*(p(i)-omega*v(i)); end do		
		!$OMP END PARALLEL DO

		call matmul_ell(v,coef,jcoef,p,n,m)

		
		alpha = rho2 / scal(r2,v,n)

		!$OMP PARALLEL DO 
		do i=1,n; ss(i) = r(i) - alpha*v(i); end do
		!$OMP END PARALLEL DO
		
		call matmul_ell(t,coef,jcoef,ss,n,m)
		
		omega = scal(t,ss,n) / scal(t,t,n)

		!$OMP PARALLEL DO 
		do i=1,n; x(i)=x(i)+alpha*p(i)+omega*ss(i); end do
		!$OMP END PARALLEL DO

		!$OMP PARALLEL DO 
		do i=1,n; r(i)=ss(i)-omega*t(i); end do
		!$OMP END PARALLEL DO

		rho1 = rho2 

		norm=0.0

		!$OMP PARALLEL DO 
		do i=1,n; norm=norm+r(i)*r(i); end do 
		!$OMP END PARALLEL DO

		norm=sqrt(norm)/norm0
		
	end do

!---------------------------------------------------!
!    SOLUTION SCALING                               !
!---------------------------------------------------!
	if(myid==0)then


			if(j.ge.itmax) then
					print*, '      '
					print*, 'non convergence =', j, norm
				else
					print*, '      '
					print*, '  Iterations BICG3 =', j, norm
			endif

	endif


	return

END SUBROUTINE solver_BICG3

SUBROUTINE matmul_ell(x,coef,jcoef,y,n,m)

	implicit none

	integer, intent(in) :: n,m
	real*8, dimension(n,m), intent (in) :: coef
	integer, dimension(m), intent(in) :: jcoef
	real*8, dimension(n), intent(in) :: y
	real*8, dimension(n), intent(out) :: x
	
	integer :: i,j,col,a,b
	
	!$OMP PARALLEL DO 
	do i=1,n
		x(i)=coef(i,1)*y(i)
	end do
	!$OMP END PARALLEL DO

	a = 0
	b = 0
	do j=2,m

		col=jcoef(j)

		if(j.eq.2)then
			a = 0 
			b = 1
		elseif (j.eq.3) then
			a = 1 
			b = 2
		elseif (j.eq.4) then
			a = 2 
			b = 3
		end if

		!$OMP PARALLEL DO 
		do i=1,n-col
			x(i)=x(i)+coef(i,j+a)*y(i+col)
		end do
		!$OMP END PARALLEL DO

		!$OMP PARALLEL DO 
		do i=1+col,n
			x(i)=x(i)+coef(i,j+b)*y(i-col)
		end do
		!$OMP END PARALLEL DO

		
		!x(i+col)=x(i+col)+coef(i+col,j+b)*y(i)
		
	end do

	return

END SUBROUTINE matmul_ell


FUNCTION scal(x,y,n) result(res)
	use omp_lib
	implicit none

	integer, intent(in) :: n
	double precision, dimension(n), intent(in) :: x,y

	integer :: i

	double precision :: res

	res=0.
	!$OMP PARALLEL DO REDUCTION(+:res)
	do i=1,n; res=res+x(i)*y(i); end do
	!$OMP END PARALLEL DO

END FUNCTION


FUNCTION norm1(x,n) result(res)
	use omp_lib
	implicit none

	integer, intent(in) :: n
	double precision, dimension(n), intent(in) :: x

	integer :: i

	double precision :: res

	res=0.
	!$OMP PARALLEL DO REDUCTION(+:res)
	do i=1,n; res=res+x(i)*x(i); end do
	!$OMP END PARALLEL DO

END FUNCTION


FUNCTION norm2(x,n) result(res)
	use omp_lib
	implicit none

	integer, intent(in) :: n
	double precision, dimension(n), intent(in) :: x

	integer :: i

	double precision :: res

	res=0.
	!$OMP PARALLEL DO REDUCTION(+:res)
	do i=1,n; res=res+x(i)*x(i); end do
	!$OMP END PARALLEL DO
	res=sqrt(res)

END FUNCTION
