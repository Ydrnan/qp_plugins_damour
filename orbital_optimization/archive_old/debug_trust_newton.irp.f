program debug_trust_newton
  implicit none

  !===================================
  ! Main program to optimize orbitals
  !===================================

  !===========
  ! Variables
  !===========

  double precision, allocatable :: H(:,:),h_f(:,:,:,:)
  double precision, allocatable :: v_grad(:),W(:,:),e_val(:)
  integer                       :: info,method
  integer                       :: n
  integer                       :: i,j,p,q,k
  double precision              :: delta,lambda1,lambda2
  integer :: nb_error
  double precision :: max_error, threshold, max_elem
  ! grad   : mo_num by mo_num double precision matrix, the gradient for the gradient method
  ! R      : mo_num by mo_num double precision matrix, rotation matrix to change the MOs
  ! H      : n by n double precision matrix, Hessian matrix
  ! Hm1    : n by n double precision matrix, inverse of the Hessian matrix
  ! v_grad : double precision vector of length n, gradient
  ! Hm1g   : double precision vector of length n, result of the product Hm1.v_grad
  ! m_Hm1g : double precision matrix builds from Hm1g
  ! info   : integer, if 0 ok, else problem in a Lapack routine
  ! method : - 0 : gradient
  !          - 1 : Full hessian
  !          - 2 : Diagonal hessian
  ! n      :  integer, n = mo_num*(mo_num-1)/2, number of orbital pairs (p,q) with p < q
  ! i,j,p,q,k : integer, indexes
  ! max_elem : double precision, maximum element value in the gradient
  ! converged : logical, if the algorithm is converged
  ! nb_iter : integer, number of iteration
  ! prev_mos : ao_num by mo_num double precision matrix containing the previous mos
  ! new_mos : ao_num by mo_num double precision matrix containing the new mos 
 
  PROVIDE mo_two_e_integrals_in_map ci_energy

  ! Choice of the method
  method = 2  ! 1 -> full h, 2 -> diag_h

  ! Display the method
  print*, 'method :', method

  ! Definition of n
  n = mo_num*(mo_num-1)/2

  !============
  ! Allocation
  !============

  allocate(v_grad(n))
  allocate(H(n,n))
  allocate(h_f(mo_num,mo_num,mo_num,mo_num))
  allocate(W(n,n),e_val(n))
  !=============
  ! Calculation
  !=============

  delta = 1d0
  
! Gradient and norm
  call gradient(n,v_grad,max_elem)
  
  ! Hessian and norm
  if (method == 1) then
   call hess(n,H,h_f)
  else
    call diag_hess(n,H,h_f)
  endif

  ! Diagonalization
  call diagonalization_hessian(n,H,e_val,w)

  ! Trust newton
  call trust_newton(n,e_val,w,v_grad,delta,lambda1)
  call trust_newton_omp(n,e_val,w,v_grad,delta,lambda2)

  print*,'Without omp:', lambda1
  print*,'Without:', lambda1

  lambda1 = 1d0 - lambda1/lambda2

  print*,'Error', lambda1

  deallocate(v_grad,H)
  deallocate(h_f)
  deallocate(W)
  deallocate(e_val)

end
