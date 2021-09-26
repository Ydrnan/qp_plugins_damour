subroutine trust_newton_omp(n,e_val,w,v_grad,delta,lambda)

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
  double precision :: dN,dN_omp, d2N, d2N_omp, fN, fN_omp
  ! dN  : double precision function, first derivative with respect to lambda of ||p||^2 - Delta 
  ! d2N : double precision function, second derivative with respect to lambda of ||p||^2 - Delta 
  ! fN  : double precision function, value of ||p||^2 
 
  !=============
  ! Calculation
  !=============

  if (debug) then
    print*,'Enter in trust_newton'
  endif

  ! Initialization
  lambda = 0d0

  ! Debug
  if (debug) then
     open(unit=10,file='debug_trust_region.dat')
     write(10,*) 'lambda, f_N, d_N, d2_N'
     do i = 1, 10000
       f_N = fN(n,e_val,w,v_grad,lambda)
       d_N = dN(n,e_val,w,v_grad,lambda,delta)
       d2_N = d2N(n,e_val,w,v_grad,lambda,delta)
       write(10,*) lambda, f_N, d_N, d2_N
       f_N = fN_omp(n,e_val,w,v_grad,lambda)
       d_N = dN_omp(n,e_val,w,v_grad,lambda,delta)
       d2_N = d2N_omp(n,e_val,w,v_grad,lambda,delta)
       write(10,*) lambda, f_N, d_N, d2_N, '(OMP)'
       lambda = lambda + 0.0001
     enddo
     close(10)
  endif
 
  ! Iteration de la methode de Newton pour trouver lambda
  CALL wall_time(t1)

  ! Diŝplay
  if (debug) then
      print*, 'Iteration   First derivative   lambda    ||x||^2'
  endif

  ! Initialization  
  i = 1
  f_N = 0d0

  do while (i <= 100 .and. ABS(1d0-f_N/delta**2)>1d-6)
    d_N = dN_omp(n,e_val,w,v_grad,lambda,delta)
    d2_N = d2N_omp(n,e_val,w,v_grad,lambda,delta)
    lambda = lambda - (1d0/ABS(d2_N))*d_N
    f_N = fN_omp(n,e_val,w,v_grad,lambda)
    
    ! Display
    if (debug) then
      print*, i, d_N, lambda, f_N, ABS(1d0-f_N/delta**2)
    endif  

    i = i+1
  enddo

  CALL wall_time(t2)

  t3 = t2 - t1
  print*,'Time to search the optimal lambda :', t3
  print*,'Number of iterations :', i
  print*,'Error on the trust region :', 1d0-f_N/delta**2

  if (debug) then
    print*,'Leave trust_newton'
  endif

end subroutine

function dN_omp(n,e_val,w,v_grad,lambda,delta)

  use omp_lib

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
  double precision :: wtg,accu1,accu2, accu12, accu22
  integer          :: i,j
  double precision, allocatable :: tmp_accu1(:), tmp_accu2(:), tmp_wtg(:)
  ! wtg   : double precision, w_i^T.v_grad
  ! accu1 : double precision, temporary variable
  ! accu2 : double precision, temporary variable 
  ! i     : integer, index

  !===========
  ! Functions
  !===========
  double precision :: dN_omp
  double precision :: ddot
  ! dN   : double precision function, first derivative with respect to lambda of ||x||^2 - Delta
  ! ddot : double precision Blas function, dot product  

  !=============
  ! Calculation
  !=============

  allocate(tmp_accu1(n), tmp_accu2(n), tmp_wtg(n))

  call omp_set_max_active_levels(1)

  ! OMP 
  !$OMP PARALLEL                                                     &
      !$OMP PRIVATE(i,j,wtg)                                         &
      !$OMP SHARED(n,lambda,v_grad, w, e_val, &
      !$OMP tmp_accu1, tmp_accu2, tmp_wtg, accu1,accu2,accu12,accu22)&
      !$OMP DEFAULT(NONE)

  ! Initialization
  !$OMP MASTER
  accu1 = 0d0
  accu2 = 0d0
  accu12 = 0d0
  accu22 = 0d0
  !$OMP END MASTER

  !$OMP DO
  do i = 1, n
    tmp_accu1(i) = 0d0
  enddo
  !$OMP END DO

  !$OMP DO
  do i = 1, n
    tmp_accu2(i) = 0d0
  enddo
  !$OMP END DO

  !$OMP DO
  do i = 1, n
    tmp_wtg = 0d0
  enddo
  !$OMP END DO

  ! Part 1
  !$OMP DO
  do i = 1, n
    do j = 1, n
      tmp_wtg(i) = tmp_wtg(i) +  w(j,i) * v_grad(j)
    enddo
  enddo
  !$OMP END DO

  !$OMP DO
  do i = 1, n
    if (e_val(i)>1e-6) then
      tmp_accu1(i) = - 2d0 * tmp_wtg(i)**2 /  (e_val(i) + lambda)**3
    endif
  enddo
  !$OMP END DO
 
  !$OMP MASTER
  do i = 1, n 
    accu1 = accu1 + tmp_accu1(i)
  enddo
  !$OMP END MASTER

  ! Part2
  !$OMP DO
  do i = 1, n
    if (e_val(i)>1e-6) then
      tmp_accu2(i) =  tmp_wtg(i)**2 / (e_val(i) + lambda)**2
    endif
  enddo
  !$OMP END DO

  !$OMP MASTER
  do i = 1, n
    accu2 = accu2 + tmp_accu2(i)
  enddo
  !$OMP END MASTER

  !$OMP END PARALLEL

  call omp_set_max_active_levels(4)

  accu2 = accu2 - delta**2 

  dN_omp = 2d0 * accu1 * accu2

  deallocate(tmp_accu1, tmp_accu2, tmp_wtg)

end function

function d2N_omp(n,e_val,w,v_grad,lambda,delta)

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
  double precision :: d2N_omp
  double precision :: ddot
  ! dN   : double precision function, first derivative with respect to lambda of ||x||^2 - Delta
  ! ddot : double precision Blas function, dot product 

  !==========
  ! Internal
  !==========
  double precision :: wtg,accu1,accu2,accu3
  double precision, allocatable :: tmp_wtg(:), tmp_accu1(:), tmp_accu2(:), tmp_accu3(:)
  integer :: i, j
  ! wtg   : double precision, w_i^T.v_grad
  ! accu1 : double precision, temporary variable
  ! accu2 : double precision, temporary variable 
  ! accu3 : double precision, temporary variable
  ! i     : integer, index

 
  !=============
  ! Calculation
  !=============
 
  allocate(tmp_wtg(n), tmp_accu1(n), tmp_accu2(n), tmp_accu3(n))

  call omp_set_max_active_levels(1)

  ! OMP 
  !$OMP PARALLEL                                                     &
      !$OMP PRIVATE(i,j,wtg)                                         &
      !$OMP SHARED(n,lambda,v_grad, w, e_val, &
      !$OMP tmp_accu1, tmp_accu2, tmp_accu3, tmp_wtg, accu1,accu2,accu3)&
      !$OMP DEFAULT(NONE)
 
  ! Initialization
  !$OMP MASTER
  accu1 = 0d0
  accu2 = 0d0
  accu3 = 0d0 
  !$OMP END MASTER

  !$OMP DO
  do i = 1, n 
    tmp_wtg(i) = 0d0
  enddo
  !$OMP END DO
  !$OMP DO
  do i = 1, n 
    tmp_accu1(i) = 0d0
  enddo
  !$OMP END DO
  !$OMP DO
  do i = 1, n
    tmp_accu2(i) = 0d0
  enddo
  !$OMP END DO
  !$OMP DO
  do i = 1, n
    tmp_accu3(i) = 0d0
  enddo
  !$OMP END DO

  ! Part 1
  !$OMP DO
  do i = 1, n
    do j = 1, n
      tmp_wtg(i) = tmp_wtg(i) +  w(j,i) * v_grad(j)
    enddo
  enddo
  !$OMP END DO

  !$OMP DO
  do i = 1, n
    if (e_val(i)>1e-6) then
      tmp_accu1(i) = 6d0 * tmp_wtg(i)**2 /  (e_val(i) + lambda)**4
    endif
  enddo
  !$OMP END DO

  !$OMP MASTER
  do i = 1, n
    accu1 = accu1 + tmp_accu1(i)
  enddo
  !$OMP END MASTER
  
  !$OMP DO
  do i = 1, n
    if (e_val(i)>1e-6) then
      tmp_accu2(i) = tmp_wtg(i)**2 /  (e_val(i) + lambda)**2
    endif
  enddo
  !$OMP END DO
 
  !$OMP MASTER
  do i = 1, n
    accu2 = accu2 + tmp_accu2(i)
  enddo
  !$OMP END MASTER

  !$OMP DO
  do i = 1, n
    if (e_val(i)>1e-6) then
      tmp_accu3(i) = -2d0 * tmp_wtg(i)**2 /  (e_val(i) + lambda)**3
    endif
  enddo
  !$OMP END DO

  !$OMP MASTER
  do i = 1, n
    accu3 = accu3 + tmp_accu3(i)
  enddo
  !$OMP END MASTER

  !$OMP END PARALLEL

  d2N_omp = 2d0 * (accu1 * (- delta**2 + accu2) + accu3**2)

  deallocate(tmp_wtg, tmp_accu1, tmp_accu2, tmp_accu3)

end function

function fN_omp(n,e_val,w,v_grad,lambda)

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

  double precision :: fN_omp
  double precision :: ddot

  double precision, allocatable :: tmp_wtg(:), tmp_fN(:)

  double precision :: wtg
  integer :: i,j

  !=============
  ! Calculation
  !=============

  allocate(tmp_wtg(n), tmp_fN(n))

  call omp_set_max_active_levels(1)

  ! OMP 
  !$OMP PARALLEL                                                     &
      !$OMP PRIVATE(i,j,wtg)                                         &
      !$OMP SHARED(n,lambda,v_grad, w, e_val, &
      !$OMP tmp_fN, tmp_wtg, fN_omp)&
      !$OMP DEFAULT(NONE)

  ! Initialization
  !$OMP MASTER
  fN_omp = 0d0
  !$OMP END MASTER

  !$OMP DO
  do i = 1, n
    tmp_fN(i) = 0d0
  enddo
  !$OMP END DO

  !$OMP DO
  do i = 1, n
    tmp_wtg(i) = 0d0
  enddo
  !$OMP END DO

  !$OMP DO
  do i = 1, n
    do j = 1, n
      tmp_wtg(i) = tmp_wtg(i) +  w(j,i) * v_grad(j)
    enddo
  enddo
  !$OMP END DO

  !$OMP DO
  do i = 1, n
    if (e_val(i)>1d-6) then
      tmp_fN(i) = tmp_wtg(i)**2 / (e_val(i) + lambda)**2
    endif
  enddo
  !$OMP END DO
  
  !$OMP MASTER
  do i = 1, n
    fN_omp = fN_omp + tmp_fN(i)
  enddo
  !$OMP END MASTER

  !$OMP END PARALLEL

  deallocate(tmp_wtg, tmp_fN)

end function

