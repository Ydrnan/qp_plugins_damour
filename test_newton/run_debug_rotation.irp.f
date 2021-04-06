program run_debug_rotation
  implicit none

  !===================================
  ! Main program to optimize orbitals
  !===================================

  !===========
  ! Variables
  !===========

  double precision, allocatable :: grad(:,:),R(:,:),R1(:,:)
  double precision, allocatable :: H(:,:),H1(:,:),h_f(:,:,:,:)
  double precision, allocatable :: Hm1(:,:),v_grad(:),m_Hm1g(:,:),Hm1g(:)
  double precision, allocatable :: prev_mos(:,:),new_mos(:,:),new_mos1(:,:)
  integer                       :: info,method
  integer                       :: n
  integer                       :: i,j,p,q,k
  integer                       :: nb_iter_trust
  double precision              :: trust_radius,rho,energy,e_model
  logical :: cancel_step
  integer :: nb_error
  double precision :: max_error, threshold
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
  ! rho    : double precision : test
  ! f_t    : double precision : test

  double precision ::  norm

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

  allocate(v_grad(n),R(mo_num,mo_num),R1(mo_num,mo_num))
  allocate(H(n,n),H1(n,n),Hm1(n,n),m_Hm1g(mo_num,mo_num),Hm1g(n))
  allocate(h_f(mo_num,mo_num,mo_num,mo_num))
  allocate(prev_mos(ao_num,mo_num))

  !=============
  ! Calculation
  !=============

  logical :: converged

  converged = .False.

  integer :: nb_iter
  double precision :: prev_energy

  trust_radius = 0d0
  prev_energy = 0.d0
  cancel_step = .False.

  nb_iter_trust = 0
  ! Gradient and norm
  call gradient(n,v_grad)
  
  ! Hessian and norm
  if (method == 1) then
    print*,'Use the full hessian matrix'
   !call first_hess(n,H)
   call hess(n,H,h_f) !h_f -> debug
  else
    print*, 'Use the diagonal hessian matrix'
    !call first_diag_hess(n,H)
    call diag_hess(n,H,h_f) !h_f -> debug
  endif

  ! Inversion of the hessian
  call trust_region(n,method,H,v_grad,m_Hm1g,prev_energy,nb_iter,trust_radius,e_model,cancel_step,prev_mos)

  ! Rotation matrix
  call dm_rotation(m_Hm1g,mo_num,R,mo_num,mo_num,info)

  call omp_dm_rotation(m_Hm1g,mo_num,R1,mo_num,mo_num,info) 

  R = R1 - R
  
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

  deallocate(v_grad,H,Hm1,m_Hm1g,R,Hm1g)
  deallocate(h_f)

end
