subroutine trust_newton(n,e_val,w,v_grad,delta,lambda)
  !===============================================
  ! Compute the lambda value for the trust region
  !===============================================

  implicit none

  integer, intent(in) :: n
  double precision, intent(in) :: e_val(n)
  double precision, intent(in) :: w(n,n)
  double precision, intent(in) :: v_grad(n)
  double precision, intent(in) :: delta 

  double precision, intent(out) :: lambda

  ! Functions
  double precision :: dN, d2N, fN  !derivees premiere et seconde par rapport à lambda de ||p||^2 - Delta 
  
  double precision :: d_N,d2_N,f_N
  double precision :: t1, t2
  integer :: i
  
  lambda = 0d0

  ! Debug
!   open(unit=10,file='norm_p.dat')
!   do i = 1, 10000
!   f_N = fN(n,e_val,w,v_grad,lambda)
!   d_N = dN(n,e_val,w,v_grad,lambda,delta)
!   d2_N = d2N(n,e_val,w,v_grad,lambda,delta)
!   write(10,*) lambda, (f_N-1d0)**2, d_N, d2_N
!   lambda = lambda + 0.0001
!   enddo
!   close(10)
 

  ! Iteration de la methode de Newton pour trouver lambda
  CALL CPU_TIME(t1)
  do i = 1, 100
    d_N = dN(n,e_val,w,v_grad,lambda,delta)
    d2_N = d2N(n,e_val,w,v_grad,lambda,delta)
    lambda = lambda - (1d0/ABS(d2_N))*d_N
    f_N = fN(n,e_val,w,v_grad,lambda)
    print*,i,d_N, lambda, f_N
  enddo
  CALL CPU_TIME(t2)
  t2=t2-t1
  print*,'Recherche du lambda optimal :', t2

end subroutine

function dN(n,e_val,w,v_grad,lambda,delta)
  implicit none
  
  integer, intent(in) :: n
  double precision, intent(in) :: e_val(n)
  double precision, intent(in) :: w(n,n)
  double precision, intent(in) :: v_grad(n)
  double precision, intent(in) :: lambda
  double precision, intent(in) :: delta
  
  double precision :: dN
  double precision :: ddot
  
  double precision :: wtg,accu1,accu2
  integer :: i

  accu1 = 0d0

  do i = 1, n
    wtg = ddot(n,w(:,i),1,v_grad,1)
    if (e_val(i)>1e-6) then
    accu1 = accu1 - 2d0 * wtg**2 / (e_val(i) + lambda)**3
    endif
  enddo 

  accu2 = 0d0

  do i = 1, n 
   wtg = ddot(n,w(:,i),1,v_grad,1)
    if (e_val(i)>1e-6) then
    accu2 = accu2 + wtg**2 / (e_val(i) + lambda)**2
    endif
  enddo

  accu2 = accu2 - delta**2 

  dN = 2d0 * accu1 * accu2 
end function

function d2N(n,e_val,w,v_grad,lambda,delta)
  implicit none

  integer, intent(in) :: n
  double precision, intent(in) :: e_val(n)
  double precision, intent(in) :: w(n,n)
  double precision, intent(in) :: v_grad(n)
  double precision, intent(in) :: lambda
  double precision, intent(in) :: delta

  double precision :: d2N
  double precision :: ddot

  double precision :: wtg,accu1,accu2,accu3
  integer :: i
 
  accu1 = 0d0
  
  do i = 1, n
    if (e_val(i)>1e-6) then
    wtg = ddot(n,w(:,i),1,v_grad,1)
    accu1 = accu1 + 6d0 * wtg**2 / (e_val(i) + lambda)**4
    endif
  enddo
  
  accu2 = 0d0
  
  do i = 1, n
    if (e_val(i)>1e-6) then
    wtg = ddot(n,w(:,i),1,v_grad,1)
    accu2 = accu2 + wtg**2 / (e_val(i) + lambda)**2
    endif
  enddo

  accu3 = 0d0
  
  do i = 1, n
    if (e_val(i)>1e-6) then
    wtg = ddot(n,w(:,i),1,v_grad,1)
    accu3 = accu3 -2d0* wtg**2 / (e_val(i) + lambda)**3
    endif
  enddo
  
  d2N = 2d0 * (accu1 * (- delta**2 + accu2) + accu3**2)

end function

function fN(n,e_val,w,v_grad,lambda)
  implicit none

  integer, intent(in) :: n
  double precision, intent(in) :: e_val(n)
  double precision, intent(in) :: w(n,n)
  double precision, intent(in) :: v_grad(n)
  double precision, intent(in) :: lambda

  double precision :: fN
  double precision :: ddot

  double precision :: wtg
  integer :: i

  fN = 0d0

  do i = 1, n
    if (e_val(i)>1e-6) then
    wtg = ddot(n,w(:,i),1,v_grad,1)
    fN = fN + wtg**2 / (e_val(i) + lambda)**2
    endif
  enddo

end function

