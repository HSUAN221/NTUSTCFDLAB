subroutine filer_final()
use variables
implicit none

   
   open (11,file='information.dat',position='append')
   write(11,*) ' ' 
   write(11,*) '======================================================='
   write(11,*) 'total cost time = ',totalcosttime,'sec'
   write(11,*) 'final time = ',time,'sec'
   write(11,*) '======================================================='
   write(11,*) ' ' 
   close(11)

   close(21)
   close(31)
  
end subroutine filer_final

subroutine filerInfo()
use variables
implicit none
   open (1,file='information.dat',position='append')

   write(1,*) '!-----------------------------------------------------------!'
   write(1,*) '!  Finite Volume Method by projection method with mpi       !'
   write(1,*) '!-----------------------------------------------------------!'
   write(1,*) ' '
   write(1,*) '======================================================='
   write(1,*) 'number of processor(MPI) = ',nproc
   write(1,*) 'number of threads(OpenMP) = ',nthreads
   write(1,*) 'number of grid = ',nx*ny*nz
   write(1,*) 'AOA = ',AOA
   write(1,*) '                '
   write(1,*) 'nx = ',nx
   write(1,*) 'ny = ',ny  
   write(1,*) 'nz = ',nz
   write(1,*) '                '
   write(1,*) 'dt = ',dt
   write(1,*) 'CFL = ', dt / dzSml + dt / dySml
   write(1,*) '                ' 
   write(1,*) 'DFIB xc = ',xc
   write(1,*) 'DFIB yc = ',yc
   write(1,*) 'DFIB zc = ',zc
   write(1,*) '                '
   write(1,*) 'Re = ',Re
   write(1,*) '                '
   if(Gridder=='non-uniform')then
   write(1,*) '====================non-uniform========================'
  !write(1,*) 'nxSml = ',nxSml
   write(1,*) 'nySml = ',nySml
   write(1,*) 'nzSml = ',nzSml
   write(1,*) '                '
   write(1,*) 'x-direction length = ',lx
   write(1,*) 'y-direction length = ',ly
   write(1,*) 'z-direction length = ',lz
   write(1,*) '                '
   write(1,*) 'x-direction small length = ',lxSml
   write(1,*) 'y-direction small length = ',lySml
   write(1,*) 'z-direction small length = ',lzSml
   write(1,*) '                '
   write(1,*) 'dx = ',dx
   write(1,*) 'Large dy = ',dy
   write(1,*) 'Large dz = ',dz
   write(1,*) '                '
   write(1,*) 'Middle dy = ',dyMid
   write(1,*) 'Middle dz = ',dzMid
   write(1,*) '                '
  !write(1,*) 'Small dx = ',dxSml
   write(1,*) 'Small dy = ',dySml
   write(1,*) 'Small dz = ',dzSml
   write(1,*) '                '
   write(1,*) 'Small grid xc = ',GridderXc
   write(1,*) 'Small grid yc = ',GridderYc
   write(1,*) 'Small grid zc = ',GridderZc
   write(1,*) '                '
   write(1,*) 'lxSml  = ',lxSml
   write(1,*) 'lySml  = ',lySml
   write(1,*) 'lzSml  = ',lzSml
   write(1,*) '====================non-uniform========================'
   write(1,*) '                '
   end if
   if(Gridder=='uniform')then
   write(1,*) '=====================uniform==========================='
   write(1,*) 'x-direction length = ',lx
   write(1,*) 'y-direction length = ',ly
   write(1,*) 'z-direction length = ',lz
   write(1,*) '                '
   write(1,*) 'dx = ',dx
   write(1,*) 'dy = ',dy
   write(1,*) 'dz = ',dz
   write(1,*) '=====================uniform==========================='
   write(1,*) '                '
   endif
   if (steadiness == 1) then; write(1,*) 'steadiness = steady'
   else if(steadiness == 2) then; write(1,*) 'steadiness = unsteady'
   end if
   write(1,*) '======================================================='
   write(1,*) ' '
   write(1,*) '======================================================='
   if (DBD == 1) then; write(1,*) 'DBD actuator on'
   else if(DBD == 0) then; write(1,*) 'DBD actuator off'
   end if
   write(1,*) '======================================================='
   write(1,*) ' '
   write(1,*) '======================================================='
   if (LES == 1) then; write(1,*) 'LES mode on'
   else if(LES == 0) then; write(1,*) 'LES mode off'
   end if
   write(1,*) 'LES Cs = ',Cs
   write(1,*) '======================================================='
   write(1,*) ' '
   write(1,*) '======================================================='
   write(1,*) 'dt = ',dt
   write(1,*) 'max time step = ',nstep*dt,'sec'
   write(1,*) 'each time step = ',isto*dt,'sec'
   write(1,*) '======================================================='
   write(1,*) ' '
   write(1,*) '======================================================='
   write(1,*) 'p Residual =',zeta
   write(1,*) 'velocity Residual =',zeta_vel
   write(1,*) '======================================================='
   

   close(1)



   open (21,file='CD_Time.dat',position='append')
   write(21,*) ' TITLE     = "" '
   write(21,*) ' VARIABLES = t(sec),C<sub>D</sub>,C<sub>L</sub> '
  

   open (31,file='CD_AOA.dat',position='append')
   write(31,*) ' TITLE     = "" '
   write(31,*) ' VARIABLES = AOA,C<sub>D</sub>,C<sub>L</sub> '


  


end subroutine filerInfo

subroutine filerProcess()
use variables
implicit none
   


   write(21,*) time, cDrag, cLift
   


   write(31,*) AOA, cDrag, cLift



   
end subroutine filerProcess

subroutine filereachtime()
use variables
implicit none

   do k=1,nz; do j=1,ny; do i=1,nx
      Qout(i,j,k,1)=pre(i,j,k)
      Qout(i,j,k,2)=uc(i,j,k)
      Qout(i,j,k,3)=vc(i,j,k)
      Qout(i,j,k,4)=wc(i,j,k)
      Qout(i,j,k,5)=ETA(i,j,k)
      !Qout(i,j,k,6)=F_tavex(i,j,k)
   enddo; enddo; enddo

	write(filename,'(I4.4)')num
 	fileformat = '.q'

   open (17,file=TRIM(filename)//fileformat,position='append',form='unformatted')
   write(17) nblocks
   write(17) nx, ny, nz
   write(17) temp, temp, temp, temp
   write(17) ( ( ( ( Qout(i,j,k,h), i = 1, nx), j = 1, ny), k = 1, nz), h = 1, 5 )
  
   close(17)
   
   num = num + 1
   
end subroutine filereachtime


