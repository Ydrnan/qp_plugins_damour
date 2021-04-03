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
  double precision, allocatable :: p(:),W(:,:)
  double precision, allocatable :: Hm1(:,:), Hm1g(:)
  double precision              :: accu, lambda!, trust_radius
  double precision              :: norm_p, norm_g, delta, rho!,e_model
  double precision, allocatable :: e_val(:),work(:,:)
  integer                       :: info,lwork!, nb_iter
  integer                       :: i,j,k
  ! p            : double precision vector of size n containing the next step
  ! W            : double precision matrix containing the eigenvectors of the hessian matrix
  ! Hm1g         : double precision vector of size n containing the next step
  ! Hm1          : double precision matrix containing the inverse of the hessian matrix
  ! accu         : double precision, temporary variable
  ! lambda       : double precision, lagrange multiplier to put the trust region contraint
  ! trust_radius : double precision, trust region radius
  ! norm_p       : double precision, norm^2 of the vector p
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

  allocate(p(n),W(n,n))
  allocate(Hm1(n,n),Hm1g(n))
  allocate(e_val(n),work(lwork,n))

  !=============
  ! Calculation
  !=============

 ! if (debug) then  
    print*,''
    print*,'==========================='
    print*,'---Enter in trust_region---'
    print*,'==========================='
 ! endif

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

  ! Initialization
  lambda =0d0

  ! Norm^2 of p
  print*,'||p||^2 :'
  norm_p = fN(n,e_val,W,v_grad,0d0)
  print*, norm_p

  ! Norm^2 of v_grad
  !norm_g = (dnrm2(n,v_grad,1))**2
  norm_g = (dnrm2(n,v_grad,1))**2
  print*,'norm grad^2 :'
  print*, norm_g

!  ! Read the iteration number
!  open(unit=10,file='nb_iteration.dat')
!    read(10,*) nb_iter
!  close(10)
!
!  ! Add 1
!  ! + 2 in the case of the cancellation of the previous step
!  open(unit=10,file='nb_iteration.dat')
!    if (nb_iter == -1) then
!      write(10,*) nb_iter + 2
!    else
!      write(10,*) nb_iter + 1
!    endif
!  close(10)

  ! Compute rho <=> the quality of the model
  call dn_rho_model(rho,nb_iter,prev_energy,e_model,cancel_step)

  ! trust radius
  ! For the first iteration trust_radius = norm_p
  ! else read the trust_radius
  if (nb_iter ==0) then

    trust_radius = norm_p !MIN(norm_g,norm_p)

!    open(unit=10,file='trust_radius.dat')
!      write(10,*) trust_radius
!    close(10)
!
!  else
!
!    open(unit=10,file='trust_radius.dat')
!      read(10,*) trust_radius
!    close(10)

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

!  ! Replacement of the trust radius for the next step
!  open(unit=10,file='trust_radius.dat')
!    write(10,*) delta**2 ! Delta
!  close(10)

  ! En donnant delta, on cherche ||p||^2 - delta^2 = 0
  ! et non ||p||^2 - delta = 0
  print*,'trust_radius :',trust_radius
  print*,'delta * coef :',delta

  if (rho >= 0.1d0) then ! minimal rho for step acceptance
    print*,'!!!previous step accepted !!!'
    cancel_step = .False.
    ! Newton method to find ||p(lambda)|| = Delta
    if (trust_radius < norm_p ) then

      print*,'Computation of the optimal lambda for the next step...'
      call trust_newton(n,e_val,W,v_grad,delta,lambda)

    else

      ! Unconstraint solution, lambda = 0
      print*,'Step in the trust region, no lambda optimization'
      lambda = 0d0

    endif

    ! Initialisation
    p = 0d0

    ! Calculation of the step p
    do i = 1, n
      if (e_val(i) > 1d-4) then
        accu = 0d0
        accu = ddot(n,W(:,i),1,v_grad,1)
        p = p - accu * W(:,i) / (e_val(i) + lambda)
      endif
    enddo

    ! pour avoir la meme chose que gHm1
    p = -p

    ! Storage of the step (in order to cancel it if rho < 0.1
    ! for the next step
    !open(unit=10,file='Hm1g.dat')
    !  do i = 1, n
    !    write(10,*) p(i)
    !  enddo
    !close(10)

  else
    ! If rho < 0.1
    print*,'!!! previous step rejected !!!'

    ! Cancellation of the previous step by applying
    ! step = - previous step
    !open(unit=10,file='Hm1g.dat')
    !  do i = 1, n
    !    read(10,*) p(i)
    !  enddo
    !close(10)

    ! Cancellation of the previous step
    !p = -p
    p=0d0
    cancel_step = .True.

    ! Replacement of e_model and prev_energy by simulating
    ! the next step as a first step
!    open(unit=10,file='nb_iteration.dat')
!      write(10,*) -1
!    close(10)

  endif

  ! Compute the predicted energy for the next step
  if (.not. cancel_step) then
    call dn_e_model(n,v_grad,H,p,prev_energy,e_model)
  else
    print*,'Cancellation of the previous step no energy prediction'
  endif
 
  ! Step transformation vector -> matrix
  ! Vector with n element -> mo_num by mo_num matrix
  do i=1,mo_num
    do j=1,mo_num
      if (i>j) then
        call in_vec_to_mat(i,j,k)
        m_Hm1g(i,j) = p(k)
      else
        m_Hm1g(i,j)=0d0
      endif
    enddo
  enddo

  ! Antisymmetrization of the previous matrix
  do i=1,mo_num
    do j=1,mo_num
      if (i<j) then
        m_Hm1g(i,j) = - m_Hm1g(j,i)
      endif
    enddo
  enddo

  ! Debug
  if (debug) then

    print*,'p'
    write(*,'(100(F10.5))') p(:)

    ! Verification
    call dm_inversion(method,n,H,Hm1)

    print*,''
    call dgemv('T',n,n,1d0,Hm1,size(Hm1,1),v_grad,1,0d0,Hm1g,1)

    print*,'vector g^T.Hm1 :'
    write(*,'(100(F10.5))') Hm1g(:)
    !print*, gHm1

    ! Calculation of the error
    p = p - Hm1g

    print*,'diff'
    do i = 1, n
      if (ABS(p(i)) > 1e-12) then
        print*,i,p(i)
      endif
    enddo

  endif

  !==============
  ! Deallocation
  !==============

  deallocate(p,W)
  deallocate(Hm1,Hm1g)
  deallocate(e_val,work)

  if (debug) then
    print*,'========================='
    print*,'---Leaves trust_region---'
    print*,'========================='
    print*,''
  endif

end subroutine
