MODULE variables
    use mpi
    implicit none 

    !This module gathers all structures and variable declarations. 
    !The different variables no need to be redefined and can this file can just be called in the subroutines.

    
    real*8, parameter            :: Re           = 10.0
          
    real*8, parameter            :: dt           = 1.0e-3
          
    integer, parameter           :: nx           = 40
          
    integer, parameter           :: ny           = 40
          
    integer, parameter           :: nz           = 40 
          
    real*8, parameter            :: lx           = 1.0
          
    real*8, parameter            :: ly           = 1.0
          
    real*8, parameter            :: lz           = 1.0

    character(len=20)            :: Gridder      = 'uniform'  !non-uniform, uniform


    !------------------- OPENMP ------------------------!

    integer                      :: nthreads

    !------------------- OPENMP ------------------------!



    !----------------------------MPI----------------------------!

    integer                      :: nproc, myid, ierr, dest

    integer                      :: status(MPI_STATUS_SIZE)

    integer, parameter           :: master=0

    integer                      :: Zdv, Zr

    integer, dimension(1:1024)   :: gstart ,gend, gend0, gcount

    integer                      :: l_nbr, r_nbr, icount, iend, istart, itag, igcount

    !----------------------------MPI----------------------------!

    

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



    !----------------Dynamic airfoil model funtion---------------------!

    real*8, parameter            :: xc               = 1.0

    real*8, parameter            :: yc               = 4.0

    real*8, parameter            :: zc               = 4.0

    real*8, parameter            :: reduce_frequency = 2.0

    real*8                       :: AOA, AOA1, frequency, StartDynamic


    !----------------Dynamic airfoil model funtion---------------------!



    !-----------------------RayCasting-------------------------!

    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)       :: ETA

    integer, parameter                               :: poly = 5000  ! equal POINT

    character(len=50)                                :: NACA_filename

    integer, parameter                               :: nSubGrids = 50

    integer                                          :: A, B, C, E

    real*8, dimension(1:poly)                        :: az, ay

    integer, dimension(1:nz,1:ny)                    :: points

    integer, dimension(1:nz,1:ny)                    :: intersection, sub_intersection

    real*8                                           :: m_pa, m_ab

    real*8, dimension(1:nz,1:ny)                     :: ETA_1

    integer                                          :: position

    !-----------------------RayCasting-------------------------!



    !-------------------------Plasma---------------------------!

    real*8 ,dimension(1:nx,1:ny,1:nz)                :: F_tavex, F_tavey

    real*8 ,dimension(1:nx,1:ny,1:nz)                :: edelta, EE

    real*8                                           :: PlasmaZc

    real*8                                           :: PlasmaYc

    !-------------------------Plasma---------------------------!



    !-------------------------iteration variable---------------------------!

    integer                                          :: ik, k, i, j, isto, istep, nstep

    real*8                                           :: time

    real*8                                           ::  VelocityDifference

    !-------------------------iteration variable---------------------------!
   
    

    !---------------------------Physical variable-----------------------!

    real*8                                              :: nu

    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)          :: p
    
    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)          :: u, v, w
    
    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)          :: u1, v1, w1
    
    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)          :: u_star, v_star, w_star
    
    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)          :: last_velocity

    real*8, dimension(1:nx,1:ny,1:nz)                   :: div,uc,vc,wc,pre

    !---------------------------Physical variable-----------------------!

  
    
    !---------------------------Gauss Seidel-----------------------!

    real*8                              :: pNew, pChange, mChange, omega, pChangeMax, mChangeMax, itmax

    real*8                              :: pChangeMax_

    !---------------------------Gauss Seidel-----------------------!



    !---------------------------QUICK---------------------------!

    real*8                              :: u_tilde_x1, u_tilde_x2, u_tilde_y1, u_tilde_y2, u_tilde_z1, u_tilde_z2 
    
    real*8                              :: v_tilde_x1, v_tilde_x2, v_tilde_y1, v_tilde_y2, v_tilde_z1, v_tilde_z2
    
    real*8                              :: w_tilde_x1, w_tilde_x2, w_tilde_y1, w_tilde_y2, w_tilde_z1, w_tilde_z2
    
    real*8                              :: ue, uw, un, us, uf, ub, vnu, vsu, wfu, wbu
    
    real*8                              :: ve, vw, vn, vs, vf, vb, uev, uwv, wfv, wbv
    
    real*8                              :: we, ww, wn, ws, wf, wb, uew, uww, vnw, vsw

    !---------------------------QUICK---------------------------!



    !---------------------------BICG----------------------------!

    integer, parameter                  :: ndim = nx * ny * nz, mdim = 4

    real*8, dimension(1:ndim,1:7)       :: coef

    integer, dimension(1:mdim)          :: jcoef  

    real*8, dimension(1:ndim)           :: div1, p_s, r_s, r2_s, v_s, ss_s, t_s

    real*8, dimension(1:ndim)           :: x1

    !---------------------------BICG----------------------------!



    !---------------------------LES-----------------------------!

    real*8, parameter                   :: Cs = 0.18

    real*8                              :: nut, delta, mutsgs

    real*8, dimension(nx,ny,nz,3)       :: dudx, dvdx, dwdx

    real*8, dimension(nx,ny,nz)         :: Viseff

    !---------------------------LES-----------------------------!



    !---------------------------Grid----------------------------!

    !Initial grid coordinates for evaluating grid lengths
    real*8, dimension (-1:nx+3)         :: X

    real*8, dimension (-1:ny+3)         :: Y

    real*8, dimension (-1:nz+3)         :: Z 

    !Actual grid cooridinates (with adjusted index)
    real*8, dimension (1:nx+1)          :: Xa

    real*8, dimension (1:ny+1)          :: Ya

    real*8, dimension (1:nz+1)          :: Za

    !Grid lengths
    real*8, dimension (-1:nx+2)         :: iDx

    real*8, dimension (-1:nx+2)         :: Dxs

    real*8, dimension (-1:ny+2)         :: iDy

    real*8, dimension (-1:ny+2)         :: Dys

    real*8, dimension (-1:nz+2)         :: iDz

    real*8, dimension (-1:nz+2)         :: Dzs

    !Midpoints of grid coordinates
    real*8, dimension (1:nx)            :: Xs

    real*8, dimension (1:ny)            :: Ys

    real*8, dimension (1:nz)            :: Zs

    !---------------------------Grid----------------------------!



    !------------------virtualForceIntegrator------------------!

    real*8                                         :: u_solid = 0
    
    real*8                                         :: v_solid = 0
    
    real*8                                         :: w_solid = 0
    
    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)     :: u2, v2, w2 
    
    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)     :: FX, FY, FZ
    
    real*8                                         :: totalFX, totalFY, totalFZ
    
    real*8                                         :: cDrag ,cLift

    real*8, dimension(1:nx,1:ny)                   :: FXz, FYz, FZz
    
    real*8, dimension(1:nx)                        :: FXy, FYy, FZy

    !------------------virtualForceIntegrator------------------!



    !--------------------- output ---------------------!
    
    character(len=20)                     :: filename, fileformat

    integer, parameter                    :: nblocks = 1
    
    real, dimension(1:nx,1:ny,1:nz)       :: Xout, Yout, Zout
    
    real, dimension(1:nx,1:ny,1:nz,5)     :: Qout
    
    real                                  :: temp = 1.0    ! mach, alpha, reyn, time 
    
    integer                               :: h,num
    
    !--------------------- output ---------------------!



    !---------------------- input ---------------------!
    
    integer                               :: inblocks
    
    integer                               :: inx
    
    integer                               :: iny
    
    integer                               :: inz 
    
    character(len=20)                     :: inputfile
    
    !---------------------- input ---------------------!



    !--------------------calculate wall time--------------------!
    
    real*8                   :: totalstarttime

    real*8                   :: totalfinaltime

    real*8                   :: totalcosttime

    !--------------------calculate wall time--------------------!



    !-------------------------chooser---------------------------!

    integer                  :: steadiness
    
    integer                  :: LES

    real*8                   :: zeta_vel

    real*8                   :: zeta

    integer                  :: DBD

    !-------------------------chooser---------------------------!



end module
