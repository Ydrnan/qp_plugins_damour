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
  double precision, allocatable :: H(:,:),h_f(:,:,:,:)
  double precision, allocatable :: v_grad(:),m_x(:,:),x(:)
  double precision, allocatable :: e_val(:), w(:,:)
  double precision, allocatable :: prev_mos(:,:),new_mos(:,:)
  integer                       :: info,method
  integer                       :: n
  integer                       :: i,j,p,q,k
  double precision              :: max_elem, delta, rho
  logical :: converged, cancel_step
  integer :: nb_iter, nb_diag, nb_cancel, nb_cancel_tot
  double precision :: prev_energy, e_model
  double precision :: t1, t2, t3
  ! grad   : mo_num by mo_num double precision matrix, the gradient for the gradient method
  ! R      : mo_num by mo_num double precision matrix, rotation matrix to change the MOs
  ! H      : n by n double precision matrix, Hessian matrix
  ! Hm1    : n by n double precision matrix, inverse of the Hessian matrix
  ! v_grad : double precision vector of length n, gradient
  ! Hm1g   : double precision vector of length n, result of the product Hm1.v_grad
  ! m_Hm1g : double precision matrix builds from Hm1g
  ! info   : integer, if 0 ok, else problem in a Lapack routine
  ! method : - 1 : Full hessian
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
  method = 2 ! 1 -> full h, 2 -> diag_h

  ! Display the method
  print*, 'Method :', method

  ! Definition of n
  n = mo_num*(mo_num-1)/2

  !============
  ! Allocation
  !============

  allocate(v_grad(n),R(mo_num,mo_num))
  allocate(H(n,n),m_x(mo_num,mo_num),x(n))
  allocate(e_val(n), w(n,n))
  allocate(h_f(mo_num,mo_num,mo_num,mo_num))
  allocate(prev_mos(ao_num,mo_num),new_mos(ao_num,mo_num))

  !=============
  ! Calculation
  !=============

  ! Initialization
  converged = .False.
  cancel_step = .False.
  delta = 0d0
  rho = 0.5d0

  call diagonalize_ci

  print *, 'CI energy : ', ci_energy
  prev_energy = 0d0
  do i = 1, N_states
    prev_energy = prev_energy + ci_energy(i) * state_average_weight(i)
  enddo
  prev_energy = prev_energy / DBLE(N_states)
  print*, 'State av energy :', prev_energy

  nb_iter = 0
  nb_cancel = 0
  nb_diag = 0
  nb_cancel_tot = 0

  do while (.not.converged)

    print*,'*********************'
    print*,'Iteration :', nb_iter
    print*,'*********************'

    print *, 'CI energy : ', ci_energy

    ! Gradient
    call gradient(n,v_grad,max_elem)
    
    ! Hessian
    if (method == 1) then
      call hess(n,H,h_f) !h_f -> debug
    else
      call diag_hess(n,H,h_f) !h_f -> debug
    endif

    call diagonalization_hessian(n,H,e_val,w)

    cancel_step = .True.
    nb_cancel = 0

    do while ( cancel_step )
      
      call trust_region(n,method,nb_iter,H,v_grad,rho,e_val,w,x,m_x,delta)

      call trust_e_model(n,v_grad,H,x,prev_energy,e_model)
 
      call rotation_matrix(m_x,mo_num,R,mo_num,mo_num,info)

      call apply_mo_rotation(R,prev_mos,new_mos)

      call clear_mo_map
      TOUCH mo_coef psi_det psi_coef
      call diagonalize_ci
      call save_wavefunction_unsorted

      call rho_model(prev_energy,e_model,rho)
     
      if (nb_iter == 0) then
        nb_iter = 1 ! in order to enable the change of delta if the first iteration is cancelled  
      endif

      if (rho >= 0.1d0) then
        cancel_step = .False.
      else
        mo_coef = prev_mos
        call save_mos
        nb_cancel = nb_cancel + 1
        nb_cancel_tot = nb_cancel_tot + 1
        print*, '***********************'
        print*, 'Step cancel : rho < 0.1'
        print*, '***********************'
      endif
      nb_diag = nb_diag + 1

      print*, 'nb_diag :', nb_diag
      print*, 'nb_cancel :', nb_cancel
      print*, 'nb_cancel_tot :', nb_cancel_tot
     
      ! exit  
      if (nb_diag >= 100) then
        print*,'nb_diag >= 100 : end'
        return
      endif

    enddo

    nb_iter = nb_iter + 1
   
    if (nb_diag >= 100 .or. ABS(max_elem) <= 1d-5) then
      converged = .True.
    endif

  enddo    

  deallocate(v_grad,H,m_x,R,x,e_val,w)
  deallocate(h_f,prev_mos,new_mos)

end program
