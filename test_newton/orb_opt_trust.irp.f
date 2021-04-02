program orb_opt
  read_wf = .True.
  TOUCH read_wf
  call run
end

subroutine run
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
  double precision, allocatable :: prev_mos(:,:)
  integer                       :: info,method
  integer                       :: n
  integer                       :: i,j,p,q,k
  integer                       :: nb_iter_trust
  double precision              :: trust_radius,rho,energy,e_model
  logical :: cancel_step
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

  PROVIDE mo_two_e_integrals_in_map

  ! Choice of the method
  method = 2  ! 1 -> full h, 2 -> diag_h

  ! Display the method
  print*, 'method :', method

  ! Definition of n
  n = mo_num*(mo_num-1)/2

  print *, 'CI energy : ', ci_energy

  !============
  ! Allocation
  !============

  allocate(v_grad(n),R(mo_num,mo_num))
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

  nb_iter = 0
  do while (.not.converged)
    if (.not.cancel_step) then
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
    endif

    ! Inversion of the hessian
    call trust_region(n,method,H,v_grad,m_Hm1g,prev_energy,nb_iter,trust_radius,e_model,cancel_step,prev_mos)

    if (cancel_step) then
      
      mo_coef = prev_mos
      call save_mos
      nb_iter_trust += 1
      print*,'step cancel :',nb_iter_trust
      
    else
      ! Rotation matrix
      call dm_rotation(m_Hm1g,mo_num,R,mo_num,mo_num,info)

      ! Orbital optimization
      call dm_newton_test(R,prev_mos)

      call clear_mo_map
      TOUCH mo_coef

      call diagonalize_ci
      nb_iter += 1
    endif

    if (nb_iter == 100) then
      converged = .True.
    endif
    !nb_iter += 1
  end do

  deallocate(v_grad,H,Hm1,m_Hm1g,R,Hm1g)
  deallocate(h_f)

end
