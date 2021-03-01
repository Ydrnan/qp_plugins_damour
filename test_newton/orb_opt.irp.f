program orb_opt
  implicit none
  
  double precision, allocatable :: grad(:,:),R(:,:)
  double precision, allocatable :: H(:,:),e_val(:),work(:,:)
  double precision, allocatable :: Hm1(:,:),v_grad(:),m_Hm1g(:,:),vec(:),Hm1_tmpr(:,:)
  integer :: info,method,n,i,j,lwork
  integer :: p,q,k
  double precision :: angle
  integer, allocatable :: ipiv(:) 
  method = 1
  
  n = mo_num*(mo_num-1)/2
  !n = mo_num**2

  ! Pour avoir les mêmes OMs que Titou
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
  allocate(H(n,n),Hm1(n,n),m_Hm1g(mo_num,mo_num))  
  allocate(Hm1_tmpr(n,n))
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
  
  else 
   
    ! Gradient and norm 
    call gradient(n,v_grad)
    
    ! Hessian and norm
    call hess(n,H)
   
    !normH = norm2(H)
    !print*, 'NormH : ', normH

    ! Inversion of the hessian
    call dm_inversion(n,H,Hm1)
     
!    Hm1 = H 
!    gHm1 = -v_grad
!    call dgesv (n , 1, H , size(H,1), ipiv , gHm1 , size(gHm1,1), info)
!    if (info/=0) then
!    call abort
!    endif       
!    call dgemv('T',n,n,1d0,Hm1,size(Hm1,1),gHm1,1,0d0,v_grad,1)

    ! Hm1.g product
    call dm_Hm1g(n,Hm1,v_grad,m_Hm1g)     

    ! Rotation matrix
    call dm_rotation(m_Hm1g,mo_num,R,mo_num,mo_num,info)
    
    ! send MOs to ocaml 
    call ocaml_debug

    ! Orbital optimization
    call dm_newton_test(R)
  
    deallocate(v_grad,H,Hm1,Hm1_tmpr,m_Hm1g,R)
 endif

end program
