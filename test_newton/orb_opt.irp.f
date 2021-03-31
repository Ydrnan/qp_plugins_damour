program orb_opt
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
  integer                       :: info,method, trust_method, cyrus
  integer                       :: n
  integer                       :: i,j,p,q,k
  double precision              :: rho,f_t
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
 
  ! Choice of the method 
  method = 2  ! 1 -> full h, 2 -> diag_h
  trust_method = 1 ! 0 -> without trust region, 1 -> with trust region
  cyrus = 0 ! 0 -> no cyrus, 1 -> cyrus

  ! Display the method
  print*, 'method :', method
  print*, 'trust_method :', trust_method
  print*, 'cyrus :', cyrus 
 
  ! Definition of n  
  n = mo_num*(mo_num-1)/2

  !============
  ! Allocation
  !============
 
  allocate(v_grad(n),R(mo_num,mo_num))
  allocate(H(n,n),H1(n,n),Hm1(n,n),m_Hm1g(mo_num,mo_num),Hm1g(n))  
  allocate(h_f(mo_num,mo_num,mo_num,mo_num))

  !=============
  ! Calculation
  !=============

  if (method == 0) then
    
    call gradient(n,v_grad)
   
    allocate(grad(mo_num,mo_num))

    call dm_vec_to_mat(v_grad,size(v_grad,1),grad,size(grad,1),info)  
    
    call dm_antisym(grad,mo_num,mo_num,info)
    
    call dm_rotation(grad,mo_num,R,mo_num,mo_num,info)

    call dm_newton_test(R) 

    deallocate(grad,R)
  
  else  ! Full or diagonal hessian 
  
    ! Gradient and norm 
    call gradient(n,v_grad)
  
    ! Hessian and norm
    if (method == 1) then 
      print*,'Use the full hessian matrix'
     !call first_hess(n,H)
     call hess(n,H,h_f) !h_f -> debug
     deallocate(h_f)
    else
      print*, 'Use the diagonal hessian matrix'
      !call first_diag_hess(n,H)
      call diag_hess(n,H,h_f) !h_f -> debug
      deallocate(h_f)
    endif
 
    ! Inversion of the hessian
    if (trust_method == 0) then
      call dm_inversion(method,n,H,Hm1)
    endif

    ! Hm1.g product
    if (trust_method == 0) then
      call dm_Hm1g(n,Hm1,v_grad,m_Hm1g,Hm1g)     
    else
      call trust_region(n,method,H,v_grad,m_Hm1g)
    endif

    if (cyrus == 1) then
      !Test cyrus : f_t
      call test_cyrus(n,H,Hm1g,f_t)
      m_Hm1g = f_t * m_Hm1g
    endif

    ! Rotation matrix
    call dm_rotation(m_Hm1g,mo_num,R,mo_num,mo_num,info)
   
    ! send MOs to ocaml 
    call ocaml_debug

    ! Orbital optimization
    call dm_newton_test(R)
  
    deallocate(v_grad,H,Hm1,m_Hm1g,R,Hm1g)
 endif

end program
