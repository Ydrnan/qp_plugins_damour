subroutine umrigar_acc_newton_v2(n,H,prev2_Hm1g,prev_Hm1g,Hm1g,prev_m,m,f_t)
  
  include 'constants.h'
  
  implicit none

  !================================
  ! Compute the Umrigar factor f_t
  !================================
  
  !====
  ! in
  !====
  integer, intent(in)          :: n
  double precision, intent(in) :: H(n,n), Hm1g(n),prev2_Hm1g(n),prev_Hm1g(n)
  double precision, intent(in) :: prev_m(n)
  ! n    : integer, n = mo_num*(mo_num-1)/2
  ! H    : n by n double precision matrix containing the hessian
  ! Hm1g : double precision vector of size n containing the next step
  ! prev2_Hm1g : double precision vector of size n containing the previous step
  ! prev_Hm1g  : double precision vector of size n containing the actual step

  !=====
  ! out
  !=====
  double precision, intent(out) :: f_t, m(n)
  ! f_t : double precision, Umrigar factor  

  !==========
  ! Internal
  !==========
  double precision, allocatable :: v(:),w(:),hv(:),hw(:)
  double precision, allocatable :: norm_v(:), accu_v(:), norm_prev_m(:), accu_prev_m(:)
  double precision              :: vhv,whw,vhw, cos_vw, epsilon, beta
  integer :: i
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

  ! Ref:
  ! Orbital Optimization in Selected Configuration Interaction Methods
  ! Yuan Yao, C. J. Umrigar
 

  !============
  ! Allocation 
  !============

  allocate(v(n),w(n),hv(n),hw(n))
  allocate(norm_v(n), accu_v(n), norm_prev_m(n), accu_prev_m(n))

  !=============
  ! Calculation
  !=============
 
  if (debug) then
    print*,'Enter in test_cyrus'
  endif

  beta = 0.5d0
 
  ! Difference between the actual and previous step
  ! x_t - x_{t-1}
  v = prev_Hm1g - prev2_Hm1g

  do i = 1, n
    norm_v(i) = dsqrt(v(i)**2)
  enddo

  do i = 1, n
    if (ABS(norm_v(i)) > 1d-6) then
      accu_v(i) = v(i)/norm_v(i)
    endif
  enddo
 
  do i = 1, n
    norm_prev_m(i) = dsqrt(prev_m(i)**2)
  enddo
  
  do i = 1, n
    if (ABS(norm_prev_m(i)) > 1d-6) then
      accu_prev_m(i) = prev_m(i)/norm_prev_m(i)
    endif
  enddo

  m = beta * accu_v + (1d0-beta)*accu_prev_m
  v = m
 
  ! Difference between the next and actual step
  ! x_{t+1} - x_t
  w = Hm1g - prev_Hm1g

  ! <v,v> 
  ! <Hm1g,Hm1g>  = Hm1g^T.H.Hm1g
  call dgemv('N',n,n,1d0,H,size(H,1),v,1,0d0,hv,1)
  vhv = ddot(n,v,1,hv,1)
 
  if (debug) then
    print*, 'Part_1', vhv
  endif 

  ! <w,w>  
  call dgemv('N',n,n,1d0,H,size(H,1),w,1,0d0,hw,1)
  whw = ddot(n,w,1,hw,1)

  if (debug) then
    print*, 'Part_2', whw
  endif

  ! <v,w>
  call dgemv('N',n,n,1d0,H,size(H,1),w,1,0d0,hw,1)
  vhw = ddot(n,v,1,hw,1)

  if (debug) then
    print*, 'Part_3', vhw 
  endif 
 
  ! cos(v,w)
  if (vhv/=0 .and. whw /=0) then
    cos_vw = vhw / dsqrt(vhv*whw)
  else 
    cos_vw = 1d0
  endif

  print*, 'cos_vw', cos_vw

  epsilon = 0.01d0
  if (cos_vw < 0d0) then 
    epsilon = epsilon**(0.8)
  endif

  ! Compute the f_t factor
  f_t = MIN(2d0/(1d0-cos_vw),1d0/epsilon)

  print*,'Umrigar factor f_t :', f_t

  if (debug) then
    print*,'Leave umrigar_acc_newton_v2'
  endif

end subroutine
