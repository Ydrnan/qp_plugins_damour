subroutine test_cyrus(n,H,Hm1g,f_t)
  
  include 'constants.h'
  
  implicit none

  !================================
  ! Compute the Umrigar factor f_t
  !================================
  
  !====
  ! in
  !====
  integer, intent(in)          :: n
  double precision, intent(in) :: H(n,n), Hm1g(n)
  ! n    : integer, n = mo_num*(mo_num-1)/2
  ! H    : n by n double precision matrix containing the hessian
  ! Hm1g : double precision vector of size n containing the next step

  !=====
  ! out
  !=====
  double precision, intent(out) :: f_t
  ! f_t : double precision, Umrigar factor  

  !==========
  ! Internal
  !==========
  double precision, allocatable :: prev2_Hm1g(:),prev_Hm1g(:)
  double precision, allocatable :: v(:),w(:),hv(:),hw(:)
  double precision              :: vhv,whw,vhw, cos_vw, epsilon
  integer :: i
  ! prev2_Hm1g : double precision vector of size n containing the previous step
  ! prev_Hm1g  : double precision vector of size n containing the actual step
  ! v          : double precision vector of size n containing the difference 
  !              between the actual and the previous step 
  ! w          : double precision vector of size n containing the difference 
  !              between the next and the actual step 
  ! hv         : double precision vector of size n containing h.v
  ! hw         : double precision vector of size n containing h.w
  ! vhv        : double precision, v^T.h.v
  ! whw        : double precision, w^T.h.w
  ! vhw        : double precision, v^T.h.w
  ! cos_uw     : double precision, cos(v,w)
  ! epsilon    : double precision, epsilon
  
  ! Function
  double precision :: ddot
  ! ddot : double precision Blas function, dot product

  !============
  ! Allocation 
  !============

  allocate(prev2_Hm1g(n),prev_Hm1g(n))
  allocate(v(n),w(n),hv(n),hw(n))

  !=============
  ! Calculation
  !=============
 
  if (debug) then
    print*,'Enter in test_cyrus'
  endif
 
  ! Read the vectors Hm1g from the previous iteration 
  open(unit=10,file='prev2_Hm1g.dat')
    do i=1,n
      read(10,*) prev2_Hm1g(i)
    enddo
  close(10)

  ! Read the vectors Hm1g from the actual iteration
  open(unit=11,file='prev_Hm1g.dat')
    do i=1,n
      read(11,*) prev_Hm1g(i)
    enddo
  close(11)

  ! Difference between the actual and previous step
  v = prev_Hm1g - prev2_Hm1g
  
  ! Difference between the next and actual step
  w = Hm1g - prev_Hm1g

  ! <v,v> 
  ! <Hm1g,Hm1g>  = Hm1g^T.H.Hm1g
  call dgemv('N',n,n,1d0,H,size(H,1),v,1,0d0,hv,1)
  vhv = ddot(n,v,1,hv,1)
 
  if (debug) then
    print*, 'test_cyrus part_1', vhv
  endif 

  ! <w,w>  
  call dgemv('N',n,n,1d0,H,size(H,1),w,1,0d0,hw,1)
  whw = ddot(n,w,1,hw,1)

  if (debug) then
    print*, 'test_cyrus part_2', whw
  endif

  ! <v,w>
  call dgemv('N',n,n,1d0,H,size(H,1),w,1,0d0,hw,1)
  vhw = ddot(n,v,1,hw,1)

  if (debug) then
    print*, 'test_cyrus part_3', vhw 
  endif 
 
  ! cos(u,w)
  if (vhv/=0 .and. whw /=0) then
    cos_vw = vhw / dsqrt(vhv*whw)
  else 
    cos_vw = 1d0
  endif

  print*, 'cos_vw', cos_vw

  epsilon = 0.01d0
  !if (cos_uw < 0d0) then 
  !  epsilon = epsilon**(0.8)
  !endif

  ! Compute the f_t factor
  f_t = MIN(1d0/(2d0-cos_vw),1d0/epsilon)

  print*,'Umrigar factor f_t :', f_t

  ! Write the vectors for the next time
  open(unit=10,file='prev2_Hm1g.dat')
    do i=1,n
      write(10,*) prev_Hm1g(i)
    enddo
  close(10)

  open(unit=11,file='prev_Hm1g.dat')
  do i=1,n
    write(11,*) (Hm1g(i))
  enddo
  close(11)

  if (debug) then
    print*,'Leaves test_cyrus'
  endif

end subroutine
