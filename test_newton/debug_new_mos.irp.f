program debug_new_mos
  implicit none

  !===================================
  ! Main program to optimize orbitals
  !===================================

  !===========
  ! Variables
  !===========

  double precision, allocatable :: grad(:,:),R(:,:)
  double precision, allocatable :: H(:,:),H1(:,:),h_f(:,:,:,:)
  double precision, allocatable :: Hm1(:,:),v_grad(:),m_Hm1g(:,:),Hm1g(:)
  double precision, allocatable :: prev_mos(:,:),new_mos(:,:),new_mos1(:,:)
  integer                       :: info,method
  integer                       :: n
  integer                       :: i,j,p,q,k
  integer                       :: nb_iter_trust
  double precision              :: trust_radius,e_model
  logical :: cancel_step
  integer :: nb_error
  double precision :: max_error, threshold
  integer :: nb_iter
  double precision :: prev_energy
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
  ! trust_radius : double precision, radius of the trust region
  ! e_model : double precision, predicted energy after the orbital rotation
  ! max_elem : double precision, maximum element value in the gradient
  ! cancel_step : logical, if the previous step must be cancel
  ! converged : logical, if the algorithm is converged
  ! nb_iter : integer, number of iteration
  ! prev_energy : double precision, previous energy
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

  allocate(v_grad(n),R(mo_num,mo_num))
  allocate(H(n,n),H1(n,n),Hm1(n,n),m_Hm1g(mo_num,mo_num),Hm1g(n))
  allocate(h_f(mo_num,mo_num,mo_num,mo_num))
  allocate(prev_mos(ao_num,mo_num),new_mos(ao_num,mo_num),new_mos1(ao_num,mo_num))

  !=============
  ! Calculation
  !=============

  ! Initialization
  trust_radius = 0d0
  prev_energy = 0.d0
  cancel_step = .False.

  nb_iter_trust = 0
  ! Gradient 
  call gradient(n,v_grad)
  
  ! Hessian and norm
  if (method == 1) then
   call hess(n,H,h_f) !h_f -> debug
  else
    call diag_hess(n,H,h_f) !h_f -> debug
  endif

  ! Inversion of the hessian
  call trust_region(n,method,H,v_grad,m_Hm1g,prev_energy,nb_iter,trust_radius,e_model,cancel_step,prev_mos)

  ! Rotation matrix
  call rotation_matrix(m_Hm1g,mo_num,R,mo_num,mo_num,info)

  ! Orbital optimization
  call apply_mo_rotation(R,prev_mos,new_mos)
  
  mo_coef = prev_mos
  call save_mos
  call apply_mo_rotation_omp(R,prev_mos,new_mos1)
 
  new_mos = new_mos -new_mos1
  
  nb_error = 0
  max_error = 0d0
  threshold = 1d-12
  
  do j = 1, mo_num
    do i = 1, ao_num
      if (ABS(new_mos(i,j)) > threshold) then
        print*,new_mos(i,j)
        nb_error = nb_error + 1
        if (ABS(new_mos(i,j)) > max_error) then
          max_error = new_mos(i,j)
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
