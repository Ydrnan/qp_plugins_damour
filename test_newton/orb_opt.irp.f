program orb_opt
  implicit none

  !===================================
  ! Main program to optimize orbitals
  !===================================
  
  !===========
  ! Variables
  !===========

  double precision, allocatable :: grad(:,:),R(:,:)
  double precision, allocatable :: H(:,:)
  double precision, allocatable :: Hm1(:,:),v_grad(:),m_Hm1g(:,:),Hm1g(:)
  integer                       :: info,method
  integer                       :: n
  integer                       :: i,j,p,q,k
  double precision              :: rho,f_t
  integer, allocatable          :: ipiv(:) 
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

  ! Choice of the method 
  method = 2
 
  ! Def of n  
  n = mo_num*(mo_num-1)/2

  ! Debug
  !mo_coef(:,1) = -mo_coef(:,1)
  !mo_coef(:,2) = -mo_coef(:,2)
  !mo_coef(:,3) = -mo_coef(:,3)
  !mo_coef(:,4) = -mo_coef(:,4)
  !mo_coef(:,6) = -mo_coef(:,6)
  !mo_coef(:,7) = -mo_coef(:,7)
  !call save_mos

  !============
  ! Allocation
  !============
 
  allocate(v_grad(n),R(mo_num,mo_num))
  allocate(H(n,n),Hm1(n,n),m_Hm1g(mo_num,mo_num),Hm1g(n))  
  allocate(ipiv(n))

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
   
    ! test : Rho 
    !call dn_rho_model(rho)  
 
    ! Gradient and norm 
    call gradient(n,v_grad)
    
    ! Hessian and norm
    if (method == 1) then 
      print*,'Use the full hessian matrix'
      call hess(n,H)
    else
      print*, 'Use the diagonal hessian matrix'
      call diag_hess(n,H)
    endif
 
    ! Inversion of the hessian
    call dm_inversion(method,n,H,Hm1)
   
!    Hm1 = H 
!    gHm1 = -v_grad
!    call dgesv (n , 1, H , size(H,1), ipiv , gHm1 , size(gHm1,1), info)
!    if (info/=0) then
!    call abort
!    endif       
!    call dgemv('T',n,n,1d0,Hm1,size(Hm1,1),gHm1,1,0d0,v_grad,1)

    ! Hm1.g product
    call dm_Hm1g(n,Hm1,v_grad,m_Hm1g,Hm1g)     

    ! test : Rho
    ! Model energy calculation
    !call dn_e_model(n,v_grad,H,Hm1g)   

    !Test cyrus : f_t
    !call test_cyrus(n,H,Hm1g,f_t)
    !m_Hm1g = f_t * m_Hm1g

    ! Rotation matrix
    call dm_rotation(m_Hm1g,mo_num,R,mo_num,mo_num,info)
   
    ! send MOs to ocaml 
    call ocaml_debug

    ! Orbital optimization
    call dm_newton_test(R)
  
    deallocate(v_grad,H,Hm1,m_Hm1g,R,Hm1g)
 endif

end program
