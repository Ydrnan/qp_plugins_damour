subroutine trust_region(n,method,nb_iter,H,v_grad,rho,e_val,w,x,m_x,delta)

  include 'constants.h'

  !=======================================================
  ! Compute the next step with the trust region algorithm
  !=======================================================

  implicit none

  !===========
  ! Variables
  !===========

  !====
  ! in
  !====
  integer, intent(in)          :: n
  integer, intent(in)          :: method ! pour la verif
  double precision, intent(in) :: H(n,n), v_grad(n), rho
  integer, intent(in)  :: nb_iter
  double precision, intent(in) :: e_val(n), w(n,n)
  ! n      : integer, n = mo_num*(mo_num-1)/2
  ! method : integer, method used to compute the hessian
  ! H      : n by n double precision matrix containing the hessian
  ! v_grad : double precision vector of size n containing the gradient
  ! rho    : double precision, represent the quality of the energy prediction
  !          with respect to the reality
  ! nb_iter : integer, number of iterations
  ! e_val : double precision vector of size n containing the eigenvalues of the hessian H
  ! w     : n by n double precision matrix containing the eigenvectors of the hessian H 

  !=======
  ! inout
  !=======
  double precision, intent(inout) :: delta

  !=====
  ! out
  !=====
  double precision, intent(out) :: m_x(mo_num,mo_num), x(n)
  ! m_x : mo_num by mo_num double precision matrix containing the next step
  ! x   : double precision vector of size n containing the next step

  !==========
  ! Internal
  !==========
  double precision, allocatable :: diff(:)
  double precision, allocatable :: Hm1(:,:), Hm1g(:)
  double precision              :: accu, lambda, trust_radius
  double precision              :: norm_x, norm_g
  integer                       :: i,j,k
  ! x            : double precision vector of size n containing the next step
  ! W            : double precision matrix containing the eigenvectors of the hessian matrix
  ! Hm1g         : double precision vector of size n containing the next step
  ! Hm1          : double precision matrix containing the inverse of the hessian matrix
  ! accu         : double precision, temporary variable
  ! lambda       : double precision, lagrange multiplier to put the trust region contraint
  ! trust_radius : double precision, trust region radius
  ! norm_x       : double precision, norm^2 of the vector x
  ! norm_g       : double precision, norm^2 of the gradient
  ! delta        : double precision, sqrt(trust_radius)
  ! rho          : double precision, ratio for the trust region
  ! e_val        : double precision vector of size n containing the eigenvalues of the hessian
  ! work         : lwork by n double precision matrix, working array for the diagonalization
  ! nb_iter      : integer, number of iteration
  ! info         : integer, info for the diagonalization
  ! lwork        : integer, for the diagonalization
  ! i,j,k        : integer, indexes

  !===========
  ! Functions
  !===========
  double precision :: ddot, dnrm2
  double precision :: fN
  ! ddot  : double precision Blas function, dot product
  ! dnrm2 : double precision Blas function, norm
  ! fN    : double precision function, compute ||p||^2

  !============
  ! Allocation
  !============

  allocate(diff(n))
  allocate(Hm1(n,n),Hm1g(n))

  !=============
  ! Calculation
  !=============

  print*,''
  print*,'==========================='
  print*,'---Enter in trust_region---'
  print*,'==========================='

  ! Initialization of the Lagrange multiplier
  lambda =0d0

  ! Norm^2 of p
  print*,'||x||^2 :'
  norm_x = fN(n,e_val,W,v_grad,0d0)
  print*, norm_x

  ! Norm^2 of v_grad
  norm_g = (dnrm2(n,v_grad,1))**2
  print*,'norm grad^2 :'
  print*, norm_g

  ! trust radius
  if (nb_iter == 0) then
    trust_radius = norm_x 
    ! Compute delta, delta = sqrt(trust_radius)
    delta = dsqrt(trust_radius)
  endif

  ! Modification of the trust radius in function of rho
  if (rho >= 0.75d0) then
    delta = 2d0 * delta
  elseif (rho >= 0.5d0) then
    delta = delta
  elseif (rho >= 0.25d0) then
    delta = 0.5d0 * delta
  else
    delta = 0.25d0 * delta
  endif
 
  print*, 'Delta :', delta

  trust_radius = delta**2
  print*, 'trust_radius :', trust_radius
  
  ! En donnant delta, on cherche ||x||^2 - delta^2 = 0
  ! et non ||x||^2 - delta = 0

  ! Newton method to find lambda such as: ||x(lambda)|| = Delta
  if (trust_radius < norm_x ) then
    print*,'Computation of the optimal lambda for the next step...'
    call trust_newton_omp(n,e_val,W,v_grad,delta,lambda)
  else
    ! Unconstraint solution, lambda = 0
    print*,'Step in the trust region, no lambda optimization'
    lambda = 0d0
  endif

  ! Initialisation
  x = 0d0

  ! Calculation of the step x
  do i = 1, n
    if (e_val(i) > 1d-4) then
      accu = 0d0
      do j = 1, n 
        accu = accu + W(j,i) * v_grad(j) ! j !!!!!
      enddo 
      !accu = ddot(n,W(:,i),1,v_grad,1)
      do j = 1, n
        x(j) = x(j) - accu * W(j,i) / (e_val(i) + lambda)
      enddo 
      !x = x - accu * W(:,i) / (e_val(i) + lambda)
    endif
  enddo

  ! pour avoir la meme chose que gHm1
  x = -x

  ! Step transformation vector -> matrix
  ! Vector with n element -> mo_num by mo_num matrix
  do j = 1, mo_num
    do i = 1, mo_num
      if (i>j) then
        call mat_to_vec_index(i,j,k)
        m_x(i,j) = x(k)
      else
        m_x(i,j)=0d0
      endif
    enddo
  enddo

  ! Antisymmetrization of the previous matrix
  do j = 1, mo_num
    do i = 1, mo_num
      if (i<j) then
        m_x(i,j) = - m_x(j,i)
      endif
    enddo
  enddo

  ! Debug
  if (debug) then
    integer :: nb_error
    double precision :: max_error

    print*,'x'
    write(*,'(100(F10.5))') x(:)

    ! Verification
    call matrix_inversion(method,n,H,Hm1)

    print*,''
    call dgemv('T',n,n,1d0,Hm1,size(Hm1,1),v_grad,1,0d0,Hm1g,1)

    print*,'vector Hm1.g :'
    write(*,'(100(F10.5))') Hm1g(:)

    ! Calculation of the error
    diff = x - Hm1g

    nb_error = 0
    max_error = 0d0
   
    print*,'diff'
    do i = 1, n
      if (ABS(x(i)) > 1e-12) then
        print*,i, diff(i)
        nb_error = nb_error + 1
        if (ABS(x(i)) > max_error) then
          max_error =  x(i)
        endif
      endif
    enddo

    print*, 'Number of errors :', nb_error
    print*, 'Max error :', max_error    

  endif

  !==============
  ! Deallocation
  !==============

  deallocate(Hm1,Hm1g)

  if (debug) then
    print*,'========================'
    print*,'---Leave trust_region---'
    print*,'========================'
    print*,''
  endif

end subroutine
