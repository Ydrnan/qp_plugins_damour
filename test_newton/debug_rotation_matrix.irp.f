program debug_rotation_matrix
  implicit none

  !===================================
  ! Main program to optimize orbitals
  !===================================

  !===========
  ! Variables
  !===========

  double precision, allocatable :: R(:,:),R2(:,:)
  double precision, allocatable :: H(:,:),h_f(:,:,:,:)
  double precision, allocatable :: v_grad(:),m_x(:,:),x(:),W(:,:),e_val(:)
  integer                       :: info,method
  integer                       :: n
  integer                       :: i,j,p,q,k
  integer                       :: nb_iter
  double precision              :: delta,rho,energy,e_model
  logical :: cancel_step
  integer :: nb_error
  double precision :: max_error, threshold, max_elem
  double precision :: prev_energy
  logical :: converged
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

  allocate(v_grad(n),R(mo_num,mo_num),R2(mo_num,mo_num))
  allocate(H(n,n),m_x(mo_num,mo_num),x(n))
  allocate(h_f(mo_num,mo_num,mo_num,mo_num))
  allocate(W(n,n),e_val(n))
  !=============
  ! Calculation
  !=============

  converged = .False.
  prev_energy = 0d0
  rho = 0.5d0
  delta = 0d0
  cancel_step = .False.

  nb_iter = 0
  
! Gradient and norm
  call gradient(n,v_grad,max_elem)
  
  ! Hessian and norm
  if (method == 1) then
   call hess(n,H,h_f) !h_f -> debug
  else
    call diag_hess(n,H,h_f) !h_f -> debug
  endif

  call diagonalization_hessian(n,H,e_val,w)

  ! Inversion of the hessian
  call trust_region(n,method,nb_iter,H,v_grad,rho,e_val,W,x,m_x,delta)

  ! Rotation matrix
  call rotation_matrix(m_x,mo_num,R,mo_num,mo_num,info) 

  call org_rotation_matrix(m_x,mo_num,R2,mo_num,mo_num,info)

  R = R2 - R
  
  nb_error = 0
  max_error = 0d0
  threshold = 1d-12
  
  do j = 1, mo_num
    do i = 1, ao_num
      if (ABS(R(i,j)) > threshold) then
        print*,R(i,j)
        nb_error = nb_error + 1
        if (ABS(R(i,j)) > max_error) then
          max_error = R(i,j)
        endif
      endif
    enddo
  enddo

  print*,''
  print*,'Threshold', threshold
  print*,'Nb error', nb_error
  print*,'max_error', max_error

  deallocate(v_grad,H,m_x,x,R,R2)
  deallocate(h_f)
  deallocate(W)
  deallocate(e_val)

end
