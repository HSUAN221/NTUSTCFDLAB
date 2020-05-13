subroutine boundary_conditions()
use variables
implicit none


!---------------------------------------------------!
!         Boundary conditions calculation           !
!---------------------------------------------------!


!!!!!!!!Velocity boundary conditions!!!!!!!!

!West vertical wall
  do k=-1,nz+2
do j=-1,ny+2

  !-------------------------------------!
  if(WestWall_u == 'no-slip') then
    u(0,j,k) = 0.0
  else if(WestWall_u == 'Neumann') then
    u(0,j,k) = u(1,j,k)
  else if(WestWall_u == 'Dirichlet') then
    u(0,j,k) = 1.0
  end if


  u(-1,j,k) = u(0,j,k)
  !-------------------------------------!
  
  !-------------------------------------!
  if(WestWall_v == 'no-slip') then
    v(0,j,k) = 2*0-v(1,j,k)
  else if(WestWall_v == 'Neumann') then
    v(0,j,k) = v(1,j,k)
  else if(WestWall_v == 'Dirichlet') then
    v(0,j,k) = 2*1-v(1,j,k)
  end if


  v(-1,j,k) = v(0,j,k)
  !-------------------------------------!

  !-------------------------------------!
  if(WestWall_w == 'no-slip') then
    w(0,j,k) = 2*0-w(1,j,k)
  else if(WestWall_w == 'Neumann') then
    w(0,j,k) = w(1,j,k)
  else if(WestWall_w == 'Dirichlet') then
    w(0,j,k) = 2*1-w(1,j,k)
  end if


  w(-1,j,k) = w(0,j,k)
  !-------------------------------------!

end do
  end do

!East vertical wall
  do k=-1,nz+2
do j=-1,ny+2

  !-------------------------------------!
  if(EastWall_u == 'no-slip') then
    u(nx,j,k) = 0.0
  else if(EastWall_u == 'Neumann') then
    u(nx,j,k) = u(nx-1,j,k)
  else if(EastWall_u == 'Dirichlet') then
    u(nx,j,k) = 1.0
  end if


  u(nx+1,j,k) = u(nx,j,k)
  u(nx+2,j,k) = u(nx+1,j,k)
  !-------------------------------------!


  !-------------------------------------!
  if(EastWall_v == 'no-slip') then
    v(nx+1,j,k) = 2*0-v(nx,j,k)
  else if(EastWall_v == 'Neumann') then
    v(nx+1,j,k) = v(nx,j,k)
  else if(EastWall_v == 'Dirichlet') then
    v(nx+1,j,k) = 2*1-v(nx,j,k)
  end if


  v(nx+2,j,k) = v(nx+1,j,k)
  !-------------------------------------!

  !-------------------------------------!
  if(EastWall_w == 'no-slip') then
    w(nx+1,j,k) = 2*0-w(nx,j,k)
  else if(EastWall_w == 'Neumann') then
    w(nx+1,j,k) = w(nx,j,k)
  else if(EastWall_w == 'Dirichlet') then
    w(nx+1,j,k) = 2*1-w(nx,j,k)
  end if


  w(nx+2,j,k) = w(nx+1,j,k)
  !-------------------------------------!
  
end do
  end do

!South horizontal wall
  do k=-1,nz+2
do i=-1,nx+2

  !-------------------------------------!
  if(SouthWall_u == 'no-slip') then
    u(i,0,k) = 2*0-u(i,1,k)
  else if(SouthWall_u == 'Neumann') then
    u(i,0,k) = u(i,1,k)
  else if(SouthWall_u == 'Dirichlet') then
    u(i,0,k) = 2*1-u(i,1,k)
  end if


  u(i,-1,k) = u(i,0,k)
  !-------------------------------------!
  

  !-------------------------------------!
  if(SouthWall_v == 'no-slip') then
    v(i,0,k) = 0.0
  else if(SouthWall_v == 'Neumann') then
    v(i,0,k) = v(i,1,k)
  else if(SouthWall_v == 'Dirichlet') then
    v(i,0,k) = 1.0
  end if


  v(i,-1,k) = v(i,0,k)
  !-------------------------------------!
  

  !-------------------------------------!
  if(SouthWall_w == 'no-slip') then
    w(i,0,k) = 2*0-w(i,1,k)
  else if(SouthWall_w == 'Neumann') then
    w(i,0,k) = w(i,1,k)
  else if(SouthWall_w == 'Dirichlet') then
    w(i,0,k) = 2*1-w(i,1,k)
  end if


  w(i,-1,k) = w(i,0,k)
  !-------------------------------------!
  
end do
  end do

!North horizontal wall
  do k=-1,nz+2
do i=-1,nx+2

  !-------------------------------------!
  if(NorthWall_u == 'no-slip') then
    u(i,ny+1,k) = 2*0-u(i,ny,k)
  else if(NorthWall_u == 'Neumann') then
    u(i,ny+1,k) = u(i,ny,k)
  else if(NorthWall_u == 'Dirichlet') then
    u(i,ny+1,k) = 2*1-u(i,ny,k)
  end if


  u(i,ny+2,k) = u(i,ny+1,k) 
  !-------------------------------------!
  

  !-------------------------------------!
  if(NorthWall_v == 'no-slip') then
    v(i,ny,k) = 0.0
  else if(NorthWall_v == 'Neumann') then
    v(i,ny,k) = v(i,ny-1,k)
  else if(NorthWall_v == 'Dirichlet') then
    v(i,ny,k) = 1.0
  end if


  v(i,ny+1,k) = v(i,ny,k)
  v(i,ny+2,k) = v(i,ny+1,k)
  !-------------------------------------!
  

  !-------------------------------------!
  if(NorthWall_w == 'no-slip') then
    w(i,ny+1,k) = 2*0-w(i,ny,k)
  else if(NorthWall_w == 'Neumann') then
    w(i,ny+1,k) = w(i,ny,k)
  else if(NorthWall_w == 'Dirichlet') then
    w(i,ny+1,k) = 2*1-w(i,ny,k)
  end if


  w(i,ny+2,k) = w(i,ny+1,k)
  !-------------------------------------!
  
end do
  end do


!inlet
!Back horizontal wall
  do j=-1,ny+2
do i=-1,nx+2

  !-------------------------------------!
  if(BackWall_u == 'no-slip') then
    u(i,j,0) =  2*0-u(i,j,1)
  else if(BackWall_u == 'Neumann') then
    u(i,j,0) = u(i,j,1)
  else if(BackWall_u == 'Dirichlet') then
    u(i,j,0) =  2*1-u(i,j,1)
  end if


  u(i,j,-1) = u(i,j,0)
  !-------------------------------------!
  

  !-------------------------------------!
  if(BackWall_v == 'no-slip') then
    v(i,j,0) =  2*0-v(i,j,1)
  else if(BackWall_v == 'Neumann') then
    v(i,j,0) = v(i,j,0)
  else if(BackWall_v == 'Dirichlet') then
    v(i,j,0) =  2*1-v(i,j,1)
  end if


  v(i,j,-1) = v(i,j,0)
  !-------------------------------------!
  

  !-------------------------------------!
  if(BackWall_w == 'no-slip') then
    w(i,j,0) =  0.0
  else if(BackWall_w == 'Neumann') then
    w(i,j,0) = w(i,j,1)
  else if(BackWall_w == 'Dirichlet') then
    w(i,j,0) =  1.0
  end if


  w(i,j,-1) = w(i,j,0)
  !-------------------------------------!
  
end do
  end do


!outlet
!Front horizontal wall
  do j=-1,ny+2
do i=-1,nx+2

  !-------------------------------------!
  if(FrontWall_u == 'no-slip') then
    u(i,j,nz+1) = 2*0-u(i,j,nz)
  else if(FrontWall_u == 'Neumann') then
    u(i,j,nz+1) = u(i,j,nz)
  else if(FrontWall_u == 'Dirichlet') then
    u(i,j,nz+1) = 2*1-u(i,j,nz)
  end if

  u(i,j,nz+2) = u(i,j,nz+1)
  !-------------------------------------!
  

  !-------------------------------------!
  if(FrontWall_v == 'no-slip') then
    v(i,j,nz+1) = 2*0-v(i,j,nz)
  else if(FrontWall_v == 'Neumann') then
    v(i,j,nz+1) = v(i,j,nz)
  else if(FrontWall_v == 'Dirichlet') then
    v(i,j,nz+1) = 2*1-v(i,j,nz)
  end if

  v(i,j,nz+2) = v(i,j,nz+1)
  !-------------------------------------!


  !-------------------------------------!
  if(FrontWall_w== 'no-slip') then
    w(i,j,nz) = 0.0
  else if(FrontWall_w == 'Neumann') then
    w(i,j,nz) = w(i,j,nz-1)
  else if(FrontWall_w == 'Dirichlet') then
    w(i,j,nz) = 1.0
  end if


  w(i,j,nz+1) = w(i,j,nz)
  w(i,j,nz+2) = w(i,j,nz+1)
  !-------------------------------------!

end do
  end do



   do k=-1,nz+2 
  do j=-1,ny+2
do i=-1,nx+2

  u1(i,j,k)=u(i,j,k)
  v1(i,j,k)=v(i,j,k)
  w1(i,j,k)=w(i,j,k)

  u2(i,j,k)=u(i,j,k)
  v2(i,j,k)=v(i,j,k)
  w2(i,j,k)=w(i,j,k)
  
  u_star(i,j,k)=u(i,j,k)
  v_star(i,j,k)=v(i,j,k)
  w_star(i,j,k)=w(i,j,k)


end do
  end do
    end do
    
!!!!!!!Pressure boundary conditions!!!!!!!

!West vertical wall
  do k=-1,nz+2 
do j=-1,ny+2
   p(0,j,k)=p(1,j,k)
   p(-1,j,k)=p(0,j,k)
end do
  end do

!East vertical wall
  do k=-1,nz+2
do j=-1,ny+2
   p(nx+1,j,k)=p(nx,j,k)
   p(nx+2,j,k)=p(nx+1,j,k)
end do
  end do

!South horizontal wall
  do k=-1,nz+2
do i=-1,nx+2
  p(i,0,k)=p(i,1,k)
  p(i,-1,k)=p(i,0,k)
end do
  end do
  
!North horizontal wall
  do k=-1,nz+2
do i=-1,nx+2
  p(i,ny+1,k)=p(i,ny,k)
  p(i,ny+2,k)=p(i,ny+1,k)
end do
  end do

!Back horizontal wall
  do j=-1,ny+2
do i=-1,nx+2
   p(i,j,0)=p(i,j,1)
   p(i,j,-1)=p(i,j,0)
end do
  end do

!Front horizontal wall
  do j=-1,ny+2
do i=-1,nx+2
   p(i,j,nz+1)=p(i,j,nz)
   p(i,j,nz+2)=p(i,j,nz+1)
end do
  end do

end subroutine boundary_conditions
