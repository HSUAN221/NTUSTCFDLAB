NTUSTCFD
--------
*****
<pre><code>
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
</code></pre>
This code is a Fortran implementation of a 3D flow using projection method with FVM. Navier Stokes equations are solved for velocity and pressure fields. The output data can be visualized with tecplot or Paraview.

How to use?
--------
The code being implemented in Fortran 90, and Fortran being a compiled language, it requires a compiler such as mpiifort or mpif90.

**File descriptions**
*****

1.**main.f90**

Simulation parameters to be defined in the main program 'main.f90'. Selecting the number of threads of OpenMP.



<pre><code>
!------------------- OPENMP ------------------------!
nthreads = 4    
call omp_set_num_threads(nthreads)
!------------------- OPENMP ------------------------!

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

DBD                            = 0                 ! DBD actuator on : 1 ; DBD actuator off : 0
   
!-----------------Parameters for the simulation------------------!
</code></pre>



2.**variables_module.f90**

Numerical values for the non-uniform mesh size, Re number and what kind of B.Cs to be defined in the program 'variables_module.f90'.


<pre><code>
real*8, parameter            :: Re           = 10.0
          
real*8, parameter            :: dt           = 1.0e-3
          
integer, parameter           :: nx           = 40
          
integer, parameter           :: ny           = 40
          
integer, parameter           :: nz           = 40 
          
real*8, parameter            :: lx           = 1.0
          
real*8, parameter            :: ly           = 1.0
          
real*8, parameter            :: lz           = 1.0

character(len=20)            :: Gridder      = 'uniform'  !non-uniform, uniform

</code></pre>

<pre><code>
!--------------------- Unequal grid ---------------------!

real*8, parameter            :: GridderXc        = 12.5
    
real*8, parameter            :: GridderYc        = 3.948
    
real*8, parameter            :: GridderZc        = 3.95
    
integer, parameter           :: nxSml            = 16
        
integer, parameter           :: nySml            = 100
    
integer, parameter           :: nyMid            = 150
    
integer, parameter           :: nzSml            = 200

integer, parameter           :: nzMid            = 300
    
real*8, parameter            :: lxSml            = 1.0
        
real*8, parameter            :: lySml            = 0.3
    
real*8, parameter            :: lyMid            = 0.6
    
real*8, parameter            :: lzSml            = 1.1

real*8, parameter            :: lzMid            = 2.0

real*8                       :: dySml, dyMid, dy

real*8                       :: dxSml, dx
    
real*8                       :: dzSml, dzMid, dz

!--------------------- Unequal grid ---------------------!

</code></pre>

<pre><code>
!----------------------------B.Cs---------------------------!

!    y=1 ______________                                                                                 
!       /             /|                                                     
!      /       N     / |                                                          
!     /____________ /  |                                
!     |  |         |   |                                                        
!     |  | B       |   |                                          
!   W |  | x=y=z=0 | E |                                           
!     |  |_________|___|x=1                                        
!     |  /         |  /                                         
!     | /     S    | /                                        
!     |/___________|/                                         
!    z=1 F     
!   Neumann     du/dn = 0
!   Dirichlet   u = 1
!   no-slip     u = 0

character(len=20)            :: WestWall_u         = 'no-slip'
character(len=20)            :: WestWall_v         = 'no-slip'
character(len=20)            :: WestWall_w         = 'no-slip'
    
character(len=20)            :: EastWall_u         = 'no-slip'
character(len=20)            :: EastWall_v         = 'no-slip'
character(len=20)            :: EastWall_w         = 'no-slip'
    
character(len=20)            :: SouthWall_u        = 'no-slip'
character(len=20)            :: SouthWall_v        = 'no-slip'
character(len=20)            :: SouthWall_w        = 'no-slip'
    
character(len=20)            :: NorthWall_u        = 'no-slip'
character(len=20)            :: NorthWall_v        = 'no-slip'
character(len=20)            :: NorthWall_w        = 'Dirichlet'
    
character(len=20)            :: BackWall_u         = 'no-slip'
character(len=20)            :: BackWall_v         = 'no-slip'
character(len=20)            :: BackWall_w         = 'no-slip'
    
character(len=20)            :: FrontWall_u        = 'no-slip'
character(len=20)            :: FrontWall_v        = 'no-slip'
character(len=20)            :: FrontWall_w        = 'no-slip'



!----------------------------B.Cs---------------------------!

</code></pre>

3.**Running the code**


1.Use the Makefile to compile all the files and create the Executive file (sol0).
2.Launch the Executive file (sol0) by command :
<pre><code>
mpirun -np A ./sol0
</code></pre>
where A is the number of processor

Mesh
--------
![Alt text](https://github.com/HSUAN221/NTUSTCFDLAB/blob/master/case/mesh.jpg)

Airfoil with plasma actuator
--------
![Alt text](https://github.com/HSUAN221/NTUSTCFDLAB/blob/master/case/plasma-on.gif)

The one-degree-of-freedom vortex-induced vibration (VIV) response of circular cylinder
--------
![Alt text](https://github.com/HSUAN221/NTUSTCFDLAB/blob/master/case/cylinder_FSI.gif)

Cavity at Re = 1000
![Alt text](https://github.com/HSUAN221/NTUSTCFDLAB/blob/master/case/cavity_Re1000.gif)
