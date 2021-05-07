program orb_opt_trust_test2
  read_wf = .True.
  TOUCH read_wf
  call run_opt_trust_test2
end

subroutine run_opt_trust_test2
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
  integer                       :: nb_iter_trust
  double precision              :: trust_radius,e_model,max_elem
  logical :: cancel_step
  logical :: converged
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

  !============
  ! Allocation
  !============
 
  ! Definition of n
  n = mo_num*(mo_num-1)/2

  allocate(v_grad(n),R(mo_num,mo_num))
  allocate(H(n,n),H1(n,n),Hm1(n,n),m_Hm1g(mo_num,mo_num),Hm1g(n))
  allocate(h_f(mo_num,mo_num,mo_num,mo_num))
  allocate(prev_mos(ao_num,mo_num),new_mos(ao_num,mo_num))

  !=============
  ! Calculation
  !=============

  ! Choice of the method
  method = 2  ! 1 -> full h, 2 -> diag_h

  ! Display the method
  print*, 'method :', method

  print *, 'CI energy : ', ci_energy

  converged = .False.
  trust_radius = 0d0
  prev_energy = 0.d0
  cancel_step = .False.
  nb_iter = 0

  call diagonalize_ci

  do while (.not.converged)

  print*,'*********************'
  print*,'Iteration :',nb_iter
  print*,'*********************'


    if (.not.cancel_step) then ! Car inutile de recalculer le gardient et l'hessien si on annule l'étape
      
      ! Gradient
      call gradient(n,v_grad,max_elem)
  
      ! Hessian 
      if (method == 1 .or. nb_iter >= 50) then
       call hess(n,H,h_f) !h_f -> debug
      else
        call diag_hess(n,H,h_f) !h_f -> debug
      endif
    endif

    ! Step in the trust region
    call trust_region(n,method,H,v_grad,m_Hm1g,prev_energy,nb_iter,trust_radius,e_model,cancel_step,prev_mos)

    if (cancel_step) then
      print*,'Cancellation of the previous step' 
      mo_coef = prev_mos
      call save_mos

      print*,''
      print*,'========================================'
      print*,'---trust region with a smaller radius---'
      print*,'========================================'
    
    else
      ! Rotation matrix
      call rotation_matrix(m_Hm1g,mo_num,R,mo_num,mo_num,info)

      ! Orbital optimization
      call apply_mo_rotation(R,prev_mos,new_mos)
  
      nb_iter += 1
    endif

    call clear_mo_map
    TOUCH mo_coef
    call diagonalize_ci
    call save_wavefunction_unsorted

    if (nb_iter == 100 .or. ABS(max_elem) <= 1d-5) then
      converged = .True.
    endif
  end do

  deallocate(v_grad,H,Hm1,m_Hm1g,R,Hm1g)
  deallocate(h_f,new_mos,prev_mos)

end
