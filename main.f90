!    y=1 ______________                                                                                 
!       /             /|       |  Author  : Zi-Hsuan Wei                                                 
!      /             / |       |  Version : 1.8                                                          
!     /____________ /  |       |  Web     : http://smetana.me.ntust.edu.tw/                              
!     |  |         |   |                                                        
!     |  |         |   |                                          
!     |  | x=y=z=0 |   |                                           
!     |  |_________|___|x=1                                        
!     |  /         |  /                                         
!     | /          | /                                        
!     |/___________|/                                         
!    z=1                                                  
!                                                       
program main
   use variables
   use mpi
   use omp_lib
   implicit none

   

   !------------------- OPENMP ------------------------!
   nthreads = 4    
   call omp_set_num_threads(nthreads)
   !------------------- OPENMP ------------------------!

   !------------------------ MPI ------------------------!
   call MPI_INIT(ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
   call Mpi_division()
   !------------------------ MPI ------------------------!
   
   !--------for nz--------!
   istart = gstart(myid)  !
   iend = gend0(myid)     !
   igcount = gcount(myid) !
   !--------for nz--------!


   
   !-----------------Parameters for the simulation------------------!
   
   omega                          = 1.5               ! Set value for SOR method

   zeta                           = 1.e-4             ! zeta for solving pressure matrix

   itmax                          = 3000              ! maximum for zeta in Gauss Seidel subroutines

   zeta_vel                       = 1.e-10            ! convergence condition for velocity field

   time                           = 0.0               ! initialize time of simulation
   
   AOA1                           = 5.D0              ! initialize AOA of simulation
   
   nstep                          = 100000            ! number of timesteps for the simulation
   
   isto                           = 100               ! data stored every 'isto' steps
   
   LES                            = 0                 ! the LES mode. 1 : on ; 0 : off
   
   steadiness                     = 2                 ! steady : 1 ; unsteady : 2 
   
   totalcosttime                  = 0                 ! initialize wall time
   
   inputfile                      = '0000.Q'          !input first q file
   
   StartDynamic                   = 1000000

   NACA_filename                  = 'NACA0012.DAT'

   DBD                            = 0                ! DBD actuator on : 1 ; DBD actuator off : 0
   
   !-----------------Parameters for the simulation------------------!



   !--------------------------------------------------------------------!
   call reading_data()              ! define dx, dy, dz, nu = 1./Re
   !--------------------------------------------------------------------!

   
   !--------------------------------------------------------------------!
   if(Gridder=='non-uniform')then
      call gridder_unequal()           ! define unequal mesh
   else if(Gridder=='uniform')then
      call gridder_equal()             ! define equal mesh
   end if


   !--------------------------------------------------------------------!

   !--------------------------------------------------------------------!
   call initial_conditions()        !call initial conditions
   !--------------------------------------------------------------------!

   !--------------------------------------------------------------------!
   !call reading_variables()        !input first q file
   !--------------------------------------------------------------------!

   !--------------------------------------------------------------------!
   call boundary_conditions()       !call boundary conditions
   !--------------------------------------------------------------------!
   
   !--------------------------------------------------------------------!
   !call MATRIX_PARSEC(myid,zc,yc,AOA,poly)   ! Creat function of NACA0012
   !call MPI_BARRIER(MPI_COMM_WORLD, ierr)
   !--------------------------------------------------------------------!

   !--------------------------------------------------------------------!
   AOA = AOA1
   call RayCasting()                     ! Creat ETA of airfoil         
   !--------------------------------------------------------------------!

   !--------------------------------------------------------------------!
   if(DBD == 1)then
      AOA = AOA1
      call plasma(Zs, Ys, nx, ny, nz, F_tavex, F_tavey, edelta, &
                  EE, myid, master, az, ay, poly, AOA)
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      do k=1,nz; do j=1,ny
         if(edelta(i,j,k)==1 .AND. myid==master)then
            write(*,*) F_tavex(i,j,k), F_tavey(i,j,k), edelta(i,j,k)
         end if
      end do; end do
   end if
   !define plasma force field
   !--------------------------------------------------------------------!
   

   if(myid==master)then

      call filereachtime()     
      call filerInfo()

   end if

    

   !----------------------------for sendrecv----------------------------!
   l_nbr = myid - 1                                                     !
   r_nbr = myid + 1                                                     !         
   if(myid == 0) then; l_nbr=MPI_PROC_NULL; end if                      !
   if(myid == (nproc-1)) then; r_nbr=MPI_PROC_NULL; end if              !
   !----------------------------for sendrecv----------------------------!


   
   write(*,*) myid, 'istart = ', istart, 'iend = ', iend, 'gcount = ', igcount


   call MPI_BARRIER(MPI_COMM_WORLD, ierr)
   
!-------------------------main loop on the timesteps----------------------!
   do istep=1,nstep                                                                      

      totalstarttime = MPI_WTIME()

      
      if( istep >= StartDynamic )then

         !--------------------------------------------------------------------!
         call dynamic_AOA()                     ! Change AOA
         !--------------------------------------------------------------------!

         !--------------------------------------------------------------------!
         call RayCasting()                      ! Creat ETA of airfoil
         !--------------------------------------------------------------------!

         !--------------------------------------------------------------------!
         if(DBD == 1)then
            call plasma(Zs, Ys, nx, ny, nz, F_tavex, F_tavey, edelta, &
                        EE, myid, master, az, ay, poly, AOA)
            call MPI_BARRIER(MPI_COMM_WORLD, ierr)
         end if
         !define plasma force field
         !--------------------------------------------------------------------!

      end if
      
   
      



      !----------data transformation among nodes----------!
      icount = 2*(nx+4)*(ny+4)
      itag = 110
      call MPI_SENDRECV( u(-1,-1,istart), icount, MPI_REAL8, l_nbr, itag, &
                         u(-1,-1,iend+1), icount, MPI_REAL8, r_nbr, itag, MPI_COMM_WORLD, status, ierr )
      itag = 120
      call MPI_SENDRECV( v(-1,-1,istart), icount, MPI_REAL8, l_nbr, itag, &
                         v(-1,-1,iend+1), icount, MPI_REAL8, r_nbr, itag, MPI_COMM_WORLD, status, ierr )
      itag = 130
      call MPI_SENDRECV( w(-1,-1,istart), icount, MPI_REAL8, l_nbr, itag, &
                         w(-1,-1,iend+1), icount, MPI_REAL8, r_nbr, itag, MPI_COMM_WORLD, status, ierr )


      itag = 140
      call MPI_SENDRECV( u(-1,-1,iend-1), icount, MPI_REAL8, r_nbr, itag, &
                         u(-1,-1,istart-2), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )
      itag = 150
      call MPI_SENDRECV( v(-1,-1,iend-1), icount, MPI_REAL8, r_nbr, itag, &
                         v(-1,-1,istart-2), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )
      itag = 160
      call MPI_SENDRECV( w(-1,-1,iend-1), icount, MPI_REAL8, r_nbr, itag, &
                         w(-1,-1,istart-2), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )



      icount = 2
      itag = 170
      call MPI_SENDRECV( iDz(istart), icount, MPI_REAL8, l_nbr, itag, &
                         iDz(iend+1), icount, MPI_REAL8, r_nbr, itag, MPI_COMM_WORLD, status, ierr )
      itag = 180
      call MPI_SENDRECV( Dzs(istart), icount, MPI_REAL8, l_nbr, itag, &
                         Dzs(iend+1), icount, MPI_REAL8, r_nbr, itag, MPI_COMM_WORLD, status, ierr )


      itag = 190
      call MPI_SENDRECV( iDz(iend-1), icount, MPI_REAL8, r_nbr, itag, &
                         iDz(istart-2), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )
      itag = 200
      call MPI_SENDRECV( Dzs(iend-1), icount, MPI_REAL8, r_nbr, itag, &
                         Dzs(istart-2), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )
                         
      
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      !----------data transformation among nodes----------!


      !------------------------------------------------------------------------------------------------------!
      call CalculateSmagorinskyViscosity()    ! calculate Smagorinsky Viscosity
      !------------------------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------------------------!
      call discretisation_QUICK_centre()      ! calculate velocity field
      !------------------------------------------------------------------------------------------------------!



      
      
      !----------data transformation among nodes----------!
      icount = (nx+4)*(ny+4)
      itag = 240
      call MPI_SENDRECV( w_star(-1,-1,iend), icount, MPI_REAL8, r_nbr, itag, &
                         w_star(-1,-1,istart-1), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )
      
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      !----------data transformation among nodes----------!


      !------------------------------------------------------------------------------------------------------!
      call gauss_seidel()       ! calculate pressure field
      !------------------------------------------------------------------------------------------------------!


      !----------data transformation among nodes----------!
      icount = (nx+4)*(ny+4)

      itag = 270
      call MPI_SENDRECV( p(-1,-1,istart), icount, MPI_REAL8, l_nbr, itag, &
                         p(-1,-1,iend+1), icount, MPI_REAL8, r_nbr, itag, MPI_COMM_WORLD, status, ierr )
                    
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      !----------data transformation among nodes----------!


      !------------------------------------------------------------------------------------------------------!
      call calcul_new_velocity() ! update velocity field 
      !------------------------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------------------------!
      call virtualForceIntegrator() ! COEFFICIENTS: CD and Center_CL,Distance,Velocity for solid  
      !------------------------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------------------------!
      call Updating_velocity()   ! update velocity field 
      !------------------------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------------------------!
      call boundary_conditions() ! recall boundary conditions to update them
      !------------------------------------------------------------------------------------------------------!
      









      time = time + dt
      !----------data collect among nodes for filer----------!
      icount = igcount*(nx)*(ny)
      !Send my results back to the master
      if(myid>master)then
         itag = 10
         call MPI_SEND( pre(1,1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
         itag = 20
         call MPI_SEND( uc(1,1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
         itag = 30
         call MPI_SEND( vc(1,1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
         itag = 40
         call MPI_SEND( wc(1,1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
      end if
      !call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      !Wait to receive results from each task
      if(myid==master)then
         do i = 1, (nproc-1)
            itag = 10
            call MPI_RECV( pre(1,1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
            itag = 20
            call MPI_RECV( uc(1,1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
            itag = 30
            call MPI_RECV( vc(1,1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
            itag = 40
            call MPI_RECV( wc(1,1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
         end do
      end if
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      !----------data collect among nodes for filer----------!



      if(myid==master)then

         write(*,*) 'time = ',time  
         call filerProcess()
               
      end if
      

       



      VelocityDifference = 1024
      if (mod(istep,isto)==0) then ! write results if isto = istep
               
         if(myid==master)then

            call filereachtime()
            
         end if 

         call check_steady()

      end if

           








      !----------calculate wall time----------!
      totalfinaltime = MPI_WTIME()
      totalcosttime = totalcosttime + (totalfinaltime-totalstarttime)
      if(myid==master)then 
         write(*,*) 'each cost = ', totalfinaltime-totalstarttime, 'total cost = ', totalcosttime, &
                    'AOA = ', AOA
      end if
      !----------calculate wall time----------!




      
      


   end do
   istep = istep - 1
!-------------------------main loop on the timesteps----------------------!


   if(myid==master)then; call filer_final(); end if

   call MPI_FINALIZE(ierr) 

end program main