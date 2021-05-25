program orb_opt_trust
  read_wf = .True.
  TOUCH read_wf
  call run_orb_opt_trust
end

subroutine run_orb_opt_trust
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
  double precision, allocatable :: prev_mos(:,:),new_mos(:,:)
  integer                       :: info,method
  integer                       :: n
  integer                       :: i,j,p,q,k
  double precision              :: trust_radius,e_model,max_elem
  logical :: cancel_step
  logical :: converged
  integer :: nb_iter
  double precision :: prev_energy
  double precision :: t1, t2, t3
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
  
  PROVIDE mo_two_e_integrals_in_map ci_energy psi_det psi_coef

  ! Choice of the method
  method = 2  ! 1 -> full h, 2 -> diag_h

  ! Display the method
  print*, 'Method :', method

  ! Definition of n
  n = mo_num*(mo_num-1)/2

  print *, 'CI energy : ', ci_energy

  !============
  ! Allocation
  !============

  allocate(v_grad(n),R(mo_num,mo_num))
  allocate(H(n,n),H1(n,n),Hm1(n,n),m_Hm1g(mo_num,mo_num),Hm1g(n))
  allocate(h_f(mo_num,mo_num,mo_num,mo_num))
  allocate(prev_mos(ao_num,mo_num),new_mos(ao_num,mo_num))

  !=============
  ! Calculation
  !=============

  ! Initialization
  converged = .False.
  trust_radius = 0d0
  prev_energy = 0.d0
  cancel_step = .False.

  call diagonalize_ci

  nb_iter = 0
  do while (.not.converged)

  print*,'*********************'
  print*,'Iteration :',nb_iter
  print*,'*********************'


    if (.not.cancel_step) then ! Car inutile de recalculer le gardient et l'hessien si on annule l'étape
      
      ! Gradient
      call gradient(n,v_grad,max_elem)
      
      ! Hessian
      if (method == 1) then
        call hess(n,H,h_f) !h_f -> debug
      else
        call diag_hess(n,H,h_f) !h_f -> debug
      endif
      call wall_time(t2)

    endif

    ! Step in the trust region
    call wall_time(t1)
    call trust_region(n,method,H,v_grad,m_Hm1g,prev_energy,nb_iter,trust_radius,e_model,cancel_step,prev_mos)
    call wall_time(t2)
    t3 = t2 - t1
    print*, 'Time for trust region :', t3

    if (cancel_step) then
      print*,'Cancellation of the previous step :', t3
      mo_coef = prev_mos
      call save_mos

      print*,''
      print*,'========================================'
      print*,'---trust region with a smaller radius---'
      print*,'========================================'
    
    else
      ! Rotation matrix
      call wall_time(t1)
      call rotation_matrix(m_Hm1g,mo_num,R,mo_num,mo_num,info)
      call wall_time(t2)
      t3 = t2 - t1
      print*, 'Time to compute the rotation matrix :', t3

      ! Orbital optimization
      call wall_time(t1)
      call apply_mo_rotation(R,prev_mos,new_mos)
      call wall_time(t2)
      t3 = t2 - t1
      print*, 'Time to apply MO rotations :', t3
      nb_iter += 1
    endif
   
    call wall_time(t1)
    call clear_mo_map
    TOUCH mo_coef psi_det psi_coef
    call diagonalize_ci
    call save_wavefunction_unsorted
    call wall_time(t2)
    t3 = t2 - t1
    print*, 'Time to diagonalize ci :', t3

    if (nb_iter == 40 .or. ABS(max_elem) <= 1d-5) then
      converged = .True.
    endif

  enddo

  deallocate(v_grad,H,Hm1,m_Hm1g,R,Hm1g)
  deallocate(h_f,prev_mos,new_mos)

end program
