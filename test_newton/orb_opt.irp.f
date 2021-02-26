program orb_opt
  implicit none
  
  double precision, allocatable :: grad(:,:),R(:,:)
  double precision, allocatable :: H(:,:),e_val(:),work(:,:),c_grad(:)
  double precision, allocatable :: Hm1(:,:),v_grad(:),gHm1(:),A(:,:),vec(:),Hm1_tmpr(:,:)
  integer :: info,method,n,i,j,lwork
  integer :: p,q 
  double precision :: angle, norm, normH
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
  allocate(H(n,n),Hm1(n,n),gHm1(n),A(mo_num,mo_num))  
  allocate(Hm1_tmpr(n,n))
  allocate(c_grad(n),ipiv(n))
  !=============
  ! Calculation
  !=============
  
  if (method == 0) then
    
    call gradient(n,v_grad)
    
    print*,'v_grad'
    write(*,'(100(F10.5))') v_grad(:)
 
    norm = norm2(v_grad)
    print*, 'Norm : ', norm 

    v_grad=0.01d0*v_grad

    allocate(grad(mo_num,mo_num))

    call dm_vec_to_mat(v_grad,size(v_grad,1),grad,size(grad,1),info)  
    
    print*,'grad'
    do i = 1, mo_num
      write(*,'(100(F10.5))') grad(i,:)
    enddo 

    call dm_antisym(grad,mo_num,mo_num,info)
  
    print*,'Antisymmetrization'
    do i = 1, mo_num
      write(*,'(100(F10.5))') grad(i,:)
    enddo

    call dm_rotation(grad,mo_num,R,mo_num,mo_num,info)
   
    print*,'Rotation matrix'
    do i = 1, mo_num
      write(*,'(100(F10.5))') R(i,:)
    enddo


    call dm_newton_test(R) 

    deallocate(grad,R)
  
  else 
   
    ! Gradient and norm 
    call gradient(n,v_grad)
    
    norm = norm2(v_grad)
    print*, 'Norm : ', norm
    
    ! Hessian and norm
    call hess(n,H)
   
    !normH = norm2(H)
    !print*, 'NormH : ', normH
    
   print*,'H'
    do i=1,n
    write(*,'(100(F10.5))') H(i,:)
    enddo   
    Hm1 = H 
    gHm1 = -v_grad
    ! Inversion of the hessian
    !call dm_inversion(n,H,Hm1)
    call dgesv (n , 1, H , size(H,1), ipiv , gHm1 , size(gHm1,1), info)
    if (info/=0) then
    call abort
    endif    
   
    ! Product Hm1^T.g
    !call dgemv('T',n,n,1d0,Hm1,size(Hm1,1),v_grad,1,0d0,gHm1,1)
     call dgemv('T',n,n,1d0,Hm1,size(Hm1,1),gHm1,1,0d0,v_grad,1)
    
    print*,'vector g^T.Hm1 :'
    write(*,'(100(F10.5))') gHm1(:) 
    write(*,'(100(F10.5))') v_grad(:)

    ! Vector with n element -> mo_num by mo_num matrix
    call dm_vec_to_mat(gHm1,n,A,mo_num,info)
    ! 
    print*,'matrix g^T.Hm1 :'  
    do i=1,mo_num
       write(*,'(100(F10.5))') A(i,:)
    enddo

!test 
!   do i=1,mo_num
!      do j=1,mo_num
!        A(i,j) = 0.001d0!angle
!      enddo
!    enddo
    
    ! Inutile d'antisym car fait dans vec to mat
    !call dm_antisym(A,mo_num,mo_num,info)
    ! Rotation matrix
    call dm_rotation(A,mo_num,R,mo_num,mo_num,info)

    print*,'Rotation matrix :'
    do i=1,mo_num
       write(*,'(100(F10.5))') R(i,:)
    enddo   

    print*,'mo_coef before rotation : '
    do i=1,mo_num
       write(*,'(100(F10.5))') mo_coef(i,:)
    enddo

   ! Orbital optimization
   ! call dm_newton_test(R)
  
    print*,'New mo_coef : '
    do i=1,mo_num
       write(*,'(100(F10.5))') mo_coef(i,:)
    enddo

!    do i = 1,15 
!    call in_mat_vec_index(i,p,q)   
!    print*,i,'p,q',p,q
!    enddo

    deallocate(v_grad,H,Hm1,Hm1_tmpr,gHm1,A,R)
 endif

end program
