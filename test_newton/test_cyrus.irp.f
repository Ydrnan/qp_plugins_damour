subroutine test_cyrus(n,H,Hm1g,f_t)
  implicit none
  
  integer, intent(in) :: n
  double precision, intent(in) :: H(n,n), Hm1g(n)

  double precision, allocatable :: prev2_Hm1g(:),prev_Hm1g(:)
  double precision, allocatable :: v(:),w(:),hv(:),hw(:)
  double precision :: vhv,whw,vhw, cos_uw, epsilon, f_t
  double precision :: ddot
  integer :: i
 
  allocate(prev2_Hm1g(n),prev_Hm1g(n))
  allocate(v(n),w(n),hv(n),hw(n))
  
  ! Read the vectors Hm1g from the previous iteration 
  open(unit=10,file='prev2_Hm1g.dat')
  do i=1,n
    read(10,*) prev2_Hm1g(i)
  enddo
  close(10)

  open(unit=11,file='prev_Hm1g.dat')
  do i=1,n
    read(11,*) prev_Hm1g(i)
  enddo
  close(11)

  v = prev_Hm1g - prev2_Hm1g
  w = Hm1g - prev_Hm1g

  ! <v,v> 
  ! <Hm1g,Hm1g>  = Hm1g^T.H.Hm1g
  call dgemv('N',n,n,1d0,H,size(H,1),v,1,0d0,hv,1)
  vhv = ddot(n,v,1,hv,1)

  print*, 'test_cyrus part_1', vhv
 
  ! <w,w>  
  call dgemv('N',n,n,1d0,H,size(H,1),w,1,0d0,hw,1)
  whw = ddot(n,w,1,hw,1)

  print*, 'test_cyrus part_2', whw


  ! <v,w>
  call dgemv('N',n,n,1d0,H,size(H,1),w,1,0d0,hw,1)
  vhw = ddot(n,v,1,hw,1)

  print*, 'test_cyrus part_3', vhw 
 
  ! cos(u,w)
  if (vhv/=0 .and. whw /=0) then
    cos_uw = vhw / dsqrt(vhv*whw)
  else 
    cos_uw = 1d0
  endif

  print*, 'cos_uw', cos_uw

  epsilon = 0.01d0

  f_t = MIN(1d0/(2d0-cos_uw),1d0/epsilon)

  print*,'f_t', f_t

  ! Write the vectors for the next time
  open(unit=10,file='prev2_Hm1g.dat')
  do i=1,n
    write(10,*) prev_Hm1g(i)
  enddo
  close(10)

  ! Car étape modifiée avec facteur f_t
  open(unit=11,file='prev_Hm1g.dat')
  do i=1,n
    write(11,*) (Hm1g(i))
  enddo
  close(11)



end subroutine
 
  
