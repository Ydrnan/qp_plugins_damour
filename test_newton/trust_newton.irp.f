subroutine trust_newton(n,e_val,w,v_grad,delta,lambda)

  include 'constants.h'

  !===============================================
  ! Compute the lambda value for the trust region
  !===============================================

  implicit none

  !===========
  ! Variables
  !===========

  !====
  ! in
  !====
  integer, intent(in) :: n
  double precision, intent(in) :: e_val(n)
  double precision, intent(in) :: w(n,n)
  double precision, intent(in) :: v_grad(n)
  double precision, intent(in) :: delta 

  !=====
  ! out
  !=====
  double precision, intent(out) :: lambda

  !==========
  ! Internal
  !==========
  double precision :: d_N,d2_N,f_N
  double precision :: t1,t2,t3
  integer          :: i
  ! d_N      : double precision, value of dN
  ! d2_N     : double precision, value of d2N
  ! f_N      : double precision, value of fN
  ! t1,t2,t3 : double precision, t3 = t2 - t1, time to search the optimal
  !            lambda value
  ! i        : integer, index

  !===========
  ! Functions
  !===========
  double precision :: dN, d2N, fN 
  ! dN  : double precision function, first derivative with respect to lambda of ||p||^2 - Delta 
  ! d2N : double precision function, second derivative with respect to lambda of ||p||^2 - Delta 
  ! fN  : double precision function, value of ||p||^2 
 
  !=============
  ! Calculation
  !=============

  ! Initialization
  lambda = 0d0

  ! Debug
  if (debug) then
     open(unit=10,file='norm_x.dat')
     write(10,*) 'lambda, f_N, d_N, d2_N'
     do i = 1, 10000
       f_N = fN(n,e_val,w,v_grad,lambda)
       d_N = dN(n,e_val,w,v_grad,lambda,delta)
       d2_N = d2N(n,e_val,w,v_grad,lambda,delta)
       write(10,*) lambda, f_N, d_N, d2_N
       lambda = lambda + 0.0001
     enddo
     close(10)
  endif
 
  ! Iteration de la methode de Newton pour trouver lambda
  CALL wall_time(t1)
  ! Il y a beaucoup trop d'itérations 

  ! Diŝplay
  if (debug) then
      print*, 'Iteration   First derivative   lambda    ||x||^2'
  endif
  
  i = 1
  f_N = 0d0
  do while (i <= 100 .and. ABS(1d0-f_N/delta**2)>1d-6)
    d_N = dN(n,e_val,w,v_grad,lambda,delta)
    d2_N = d2N(n,e_val,w,v_grad,lambda,delta)
    lambda = lambda - (1d0/ABS(d2_N))*d_N
    f_N = fN(n,e_val,w,v_grad,lambda)
    
    ! Diŝplay
    if (debug) then
      print*, i, d_N, lambda, f_N, ABS(1d0-f_N/delta**2)
    endif  
    i = i+1
  enddo

  CALL wall_time(t2)

  t3 = t2 - t1
  print*,'Time to search the optimal lambda :', t3
  print*,'Error on the trust region :', 1d0-f_N/delta**2

end subroutine

function dN(n,e_val,w,v_grad,lambda,delta)

  !=============================================================
  ! Function to compute the first derivative of ||x||^2 - Delta
  !=============================================================

  implicit none
  
  !===========
  ! Variables
  !===========

  !====
  ! in
  !====
  integer, intent(in)          :: n
  double precision, intent(in) :: e_val(n)
  double precision, intent(in) :: w(n,n)
  double precision, intent(in) :: v_grad(n)
  double precision, intent(in) :: lambda
  double precision, intent(in) :: delta
  ! n      : integer, n = mo_num*(mo_num-1)/2
  ! e_val  : double precision vector of size n containing the eignevalues of the hessian
  ! w      : n by n double precision matrix containing the eigenvectors of the hessian
  ! v_grad : double precision vector of size n containing the gradient
  ! lambda : double precision, lagrange multiplier for the trust region
  ! delta  : double precision, trust radius 
   
  !==========
  ! Internal
  !==========
  double precision :: wtg,accu1,accu2
  integer          :: i
  ! wtg   : double precision, w_i^T.v_grad
  ! accu1 : double precision, temporary variable
  ! accu2 : double precision, temporary variable 
  ! i     : integer, index

  !===========
  ! Functions
  !===========
  double precision :: dN
  double precision :: ddot
  ! dN   : double precision function, first derivative with respect to lambda of ||x||^2 - Delta
  ! ddot : double precision Blas function, dot product  

  !=============
  ! Calculation
  !=============

  ! Initialization
  accu1 = 0d0
  accu2 = 0d0

  do i = 1, n
    wtg = ddot(n,w(:,i),1,v_grad,1)
    if (e_val(i)>1e-6) then
    accu1 = accu1 - 2d0 * wtg**2 / (e_val(i) + lambda)**3
    endif
  enddo 

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

  !=============================================================
  ! Function to compute the second derivative of ||x||^2 - Delta
  !=============================================================

  implicit none

  !===========
  ! Variables
  !===========

  !====
  ! in
  !====
  integer, intent(in) :: n
  double precision, intent(in) :: e_val(n)
  double precision, intent(in) :: w(n,n)
  double precision, intent(in) :: v_grad(n)
  double precision, intent(in) :: lambda
  double precision, intent(in) :: delta
  ! n      : integer, n = mo_num*(mo_num-1)/2
  ! e_val  : double precision vector of size n containing the eignevalues of the hessian
  ! w      : n by n double precision matrix containing the eigenvectors of the hessian
  ! v_grad : double precision vector of size n containing the gradient
  ! lambda : double precision, lagrange multiplier for the trust region
  ! delta  : double precision, trust radius

  !===========
  ! Functions
  !===========
  double precision :: d2N
  double precision :: ddot
  ! dN   : double precision function, first derivative with respect to lambda of ||x||^2 - Delta
  ! ddot : double precision Blas function, dot product 

  !==========
  ! Internal
  !==========
  double precision :: wtg,accu1,accu2,accu3
  integer :: i
  ! wtg   : double precision, w_i^T.v_grad
  ! accu1 : double precision, temporary variable
  ! accu2 : double precision, temporary variable 
  ! accu3 : double precision, temporary variable
  ! i     : integer, index

 
  !=============
  ! Calculation
  !=============
  
  ! Initialization
  accu1 = 0d0
  accu2 = 0d0
  accu3 = 0d0 
 
  do i = 1, n
    if (e_val(i)>1d-6) then
      wtg = ddot(n,w(:,i),1,v_grad,1)
      accu1 = accu1 + 6d0 * wtg**2 / (e_val(i) + lambda)**4
    endif
  enddo
  
  do i = 1, n
    if (e_val(i)>1d-6) then
      wtg = ddot(n,w(:,i),1,v_grad,1)
      accu2 = accu2 + wtg**2 / (e_val(i) + lambda)**2
    endif
  enddo

  do i = 1, n
    if (e_val(i)>1d-6) then
      wtg = ddot(n,w(:,i),1,v_grad,1)
      accu3 = accu3 -2d0* wtg**2 / (e_val(i) + lambda)**3
    endif
  enddo
  
  d2N = 2d0 * (accu1 * (- delta**2 + accu2) + accu3**2)

end function

function fN(n,e_val,w,v_grad,lambda)

  !==============================
  ! Compute the value of ||x||^2
  !==============================

  implicit none

  !===========
  ! Variables
  !===========

  integer, intent(in) :: n
  double precision, intent(in) :: e_val(n)
  double precision, intent(in) :: w(n,n)
  double precision, intent(in) :: v_grad(n)
  double precision, intent(in) :: lambda

  double precision :: fN
  double precision :: ddot

  double precision :: wtg
  integer :: i

  !=============
  ! Calculation
  !=============

  ! Initialization
  fN = 0d0

  do i = 1, n
    if (e_val(i)>1d-6) then
      wtg = ddot(n,w(:,i),1,v_grad,1)
      fN = fN + wtg**2 / (e_val(i) + lambda)**2
    endif
  enddo

end function

