subroutine trust_region(n,method,H,v_grad,m_Hm1g, prev_energy,nb_iter,trust_radius,e_model,cancel_step,prev_mos)

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
  double precision, intent(in) :: H(n,n), v_grad(n)
  double precision, intent(inout) :: prev_energy,trust_radius
  integer, intent(inout)  :: nb_iter
  double precision, intent(inout) :: e_model
  logical, intent(inout) :: cancel_step
  double precision, intent(in) :: prev_mos(ao_num,mo_num)
  ! n      : integer, n = mo_num*(mo_num-1)/2
  ! method : integer, method used to compute the hessian
  ! H      : n by n double precision matrix containing the hessian
  ! v_grad : double precision vector of size n containing the gradient

  !=====
  ! out
  !=====
  double precision, intent(out) :: m_Hm1g(mo_num,mo_num)
  ! m_Hm1g : mo_num by mo_num double precision matrix containing the next step

  !==========
  ! Internal
  !==========
  double precision, allocatable :: x(:),W(:,:)
  double precision, allocatable :: Hm1(:,:), Hm1g(:)
  double precision              :: accu, lambda!, trust_radius
  double precision              :: norm_x, norm_g, delta, rho!,e_model
  double precision, allocatable :: e_val(:),work(:,:)
  integer                       :: info,lwork!, nb_iter
  integer                       :: i,j,k
  integer                       :: nb_negative_vp  
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

  lwork=3*n-1

  allocate(x(n),W(n,n))
  allocate(Hm1(n,n),Hm1g(n))
  allocate(e_val(n),work(lwork,n))

  !=============
  ! Calculation
  !=============

  print*,''
  print*,'==========================='
  print*,'---Enter in trust_region---'
  print*,'==========================='

  ! Copy the hessian matrix, the eigenvectors will be store in W
  W=H

  ! Diagonalization of the hessian
  call dsyev('V','U',n,W,size(W,1),e_val,work,lwork,info)

  if (info /= 0) then
      print*, 'Error diagonalization : trust_region'
      call ABORT
  endif

  if (debug) then
    print *, 'vp Hess:'
    write(*,'(100(F10.5))')  real(e_val(:))
  endif

  ! Number of negative eigenvalues
  nb_negative_vp = 0
  do i = 1, n
    if (e_val(i) < -1d-12) then 
      nb_negative_vp = nb_negative_vp + 1
      print*,'e_val < 0 :', e_val(i)
    endif
  enddo 
  print*,'Number of negative eigenvalues :', nb_negative_vp

  ! Initialization of the Lagrange multiplier
  lambda =0d0

  ! Norm^2 of p
  print*,'||p||^2 :'
  norm_x = fN(n,e_val,W,v_grad,0d0)
  print*, norm_x

  ! Norm^2 of v_grad
  !norm_g = (dnrm2(n,v_grad,1))**2
  norm_g = (dnrm2(n,v_grad,1))**2
  print*,'norm grad^2 :'
  print*, norm_g

  ! Compute rho <=> the quality of the model
  call rho_model(rho,nb_iter,prev_energy,e_model,cancel_step)

  ! trust radius
  ! For the first iteration trust_radius = norm_p
  ! else read the trust_radius
  if (nb_iter == 0) then

    trust_radius = norm_x 

  endif

  ! Compute delta, delta = sqrt(trust_radius)
  delta = dsqrt(trust_radius)
  print*, 'delta', delta

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
  
  trust_radius = delta**2

  ! En donnant delta, on cherche ||x||^2 - delta^2 = 0
  ! et non ||x||^2 - delta = 0
  print*,'trust_radius :',trust_radius
  print*,'delta * coef :',delta

  if (rho >= 0.1d0) then ! minimal rho for step acceptance

    print*,'!!!previous step accepted !!!'
    cancel_step = .False.

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

  else
    ! If rho < 0.1
    print*,'!!! previous step rejected !!!'

    ! Cancellation of the previous step
    x = 0d0
    cancel_step = .True.

  endif

  ! Compute the predicted energy for the next step
  if (.not. cancel_step) then
    call trust_e_model(n,v_grad,H,x,prev_energy,e_model)
  else
    print*,'Cancellation of the previous step no energy prediction'
  endif
 
  ! Step transformation vector -> matrix
  ! Vector with n element -> mo_num by mo_num matrix
  do j = 1, mo_num
    do i = 1, mo_num
      if (i>j) then
        call mat_to_vec_index(i,j,k)
        m_Hm1g(i,j) = x(k)
      else
        m_Hm1g(i,j)=0d0
      endif
    enddo
  enddo

  ! Antisymmetrization of the previous matrix
  do j = 1, mo_num
    do i = 1, mo_num
      if (i<j) then
        m_Hm1g(i,j) = - m_Hm1g(j,i)
      endif
    enddo
  enddo

  ! Debug
  if (debug) then

    print*,'x'
    write(*,'(100(F10.5))') x(:)

    ! Verification
    call matrix_inversion(method,n,H,Hm1)

    print*,''
    call dgemv('T',n,n,1d0,Hm1,size(Hm1,1),v_grad,1,0d0,Hm1g,1)

    print*,'vector Hm1.g :'
    write(*,'(100(F10.5))') Hm1g(:)
    !print*, gHm1

    ! Calculation of the error
    x = x - Hm1g

    print*,'diff'
    do i = 1, n
      if (ABS(x(i)) > 1e-12) then
        print*,i,x(i)
      endif
    enddo

  endif

  !==============
  ! Deallocation
  !==============

  deallocate(x,W)
  deallocate(Hm1,Hm1g)
  deallocate(e_val,work)

  if (debug) then
    print*,'========================='
    print*,'---Leaves trust_region---'
    print*,'========================='
    print*,''
  endif

end subroutine
