program orb_opt
  read_wf = .True.
  TOUCH read_wf
  call run_orb_opt
end

subroutine run_orb_opt
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
  integer                       :: info,method
  integer                       :: n
  integer                       :: i,j,p,q,k
  double precision, allocatable :: prev_mos(:,:), new_mos(:,:)
  logical :: converged
  integer :: nb_iter
  double precision :: max_elem
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
 
  PROVIDE mo_two_e_integrals_in_map

  !============
  ! Allocation
  !============

  ! Definition of n  
  n = mo_num*(mo_num-1)/2

  allocate(prev_mos(ao_num,mo_num),new_mos(ao_num,mo_num)) 
  allocate(v_grad(n),R(mo_num,mo_num))
  allocate(H(n,n),H1(n,n),Hm1(n,n),m_Hm1g(mo_num,mo_num),Hm1g(n))  
  allocate(h_f(mo_num,mo_num,mo_num,mo_num))

  !=============
  ! Calculation
  !=============
 
  ! Choice of the method 
  method = 2  ! 1 -> full h, 2 -> diag_h

  ! Display the method
  print*, 'method :', method

  print *, 'CI energy : ', ci_energy

  call diagonalize_ci

  ! Initialization
  converged = .False.
  nb_iter = 0
  max_elem = 1d0
  do while (.not.converged)

    print*,'*********************'
    print*,'Iteration :',nb_iter
    print*,'*********************'
   
  
    if (method == 0) then
      
      call gradient(n,v_grad,max_elem)
     
      allocate(grad(mo_num,mo_num))
  
      call vec_to_mat(v_grad,size(v_grad,1),grad,size(grad,1),info)  
      
      call matrix_antisym(grad,mo_num,mo_num,info)
      
      call rotation_matrix(grad,mo_num,R,mo_num,mo_num,info)
  
      call apply_mo_rotation(R) 
  
      deallocate(grad,R)
    
    else  ! Full or diagonal hessian 
    
      ! Gradient 
      call gradient(n,v_grad,max_elem)
    
      ! Hessian
      if (method == 1) then 
       call hess(n,H,h_f) !h_f -> debug
      else
        call diag_hess(n,H,h_f) !h_f -> debug
      endif
   
      ! Inversion of the hessian
      call matrix_inversion(method,n,H,Hm1)
  
      ! Hm1.g product
      call Hm1g_product(n,Hm1,v_grad,m_Hm1g,Hm1g)     
  
      ! Rotation matrix
      call rotation_matrix(m_Hm1g,mo_num,R,mo_num,mo_num,info)
     
      ! Orbital optimization
      call apply_mo_rotation(R,prev_mos,new_mos)
    
    endif
  
    call clear_mo_map
    TOUCH mo_coef
    call diagonalize_ci
    call save_wavefunction_unsorted
  
    nb_iter = nb_iter+1
    if (nb_iter == 40 .or. ABS(max_elem) <= 1d-5) then
        converged = .True.
    endif

  enddo

  deallocate(v_grad,H,Hm1,m_Hm1g,R,Hm1g)
  deallocate(h_f,prev_mos,new_mos)
end 
