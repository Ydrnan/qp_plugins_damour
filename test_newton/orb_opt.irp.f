program orb_opt
  implicit none
  
  double precision, allocatable :: grad(:,:),R(:,:)
  double precision, allocatable :: H(:,:),e_val(:),work(:,:)
  double precision, allocatable :: Hm1(:,:),v_grad(:),gHm1(:),A(:,:),vec(:),Hm1_tmpr(:,:)
  integer :: info,method,n,i,j,lwork
  integer :: p,q 
  double precision :: angle, norm, normH
  
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
    
    call gradient(n,v_grad)
    
    norm = norm2(v_grad)
    print*, 'Norm : ', norm
    
    call hess(n,H)
   
    !normH = norm2(H)
    !print*, 'NormH : ', normH

    call dm_inversion(n,H,Hm1)
   
!    ! Inversion de H 
!    lwork=3*n-1
!    allocate(work(lwork,n),e_val(n))
!    
!    call dsyev('V','U',n,H,size(H,1),e_val,work,lwork,info)
!    if (info /= 0) then
!        call ABORT
!    endif
!   
!    Hm1=0d0 
!    !H_diag=0d0
!    print *, 'Hess:', real(e_val(:))
!    do i=1,n
!        !print*,'H_val',e_val(i)
!        !H_diag(i,i)=e_val(i)
!        if ( (e_val(i)>0.d0)) then
!              Hm1(i,i)=1d0/max(1.d-4,e_val(i))
!        else
!              Hm1(i,i)=0d0  !1d0/min(-1.d-4,e_val(i))
!        endif
!    enddo
!    
!    deallocate(work,e_val)
!    
!    call dgemm('N','T',n,n,n,1d0,Hm1,size(Hm1,1),H,size(H,1),0d0,Hm1_tmpr,size(Hm1_tmpr,1))
!    call dgemm('N','N',n,n,n,1d0,H,size(H,1),Hm1_tmpr,size(Hm1_tmpr,1),0d0,Hm1,size(Hm1,1))
!    ! Fin de l'inversion

    v_grad = 0.5d0 * v_grad
    
    call dgemv('T',n,n,1d0,Hm1,size(Hm1,1),v_grad,1,0d0,gHm1,1)
    
    print*,'vecteur gHm1'
    write(*,'(100(F10.5))') gHm1(:) 

    call dm_vec_to_mat(gHm1,n,A,mo_num,info)
    
    ! print matrice a exponentialiser
    print*,'gHm1 en matrice'  
    do i=1,mo_num
       write(*,'(100(F10.5))') A(i,:)
    enddo

!test 
!   do i=1,mo_num
!      do j=1,mo_num
!        A(i,j) = 0.001d0!angle
!      enddo
!    enddo
    
!     A=0d0
!    A(1,3) = -0.0054456d0
!    A(2,4) = 1.87241d0
!    A(3,1) = 0.0054456d0
!    A(4,2) = -1.87241d0

!    print*,'A'
!    do i=1,mo_num
!      print*,A(i,:)
!    enddo

    call dm_antisym(A,mo_num,mo_num,info)
    call dm_rotation(A,mo_num,R,mo_num,mo_num,info)

    print*,'Rotation matrix : '
    do i=1,mo_num
       write(*,'(100(F10.5))') R(i,:)
    enddo   

    print*,'mo_coef before rotation : '
    do i=1,mo_num
       write(*,'(100(F10.5))') mo_coef(i,:)
    enddo


   ! call dm_newton_test(R)
  
    print*,'new_mo_coef : '
    do i=1,mo_num
       write(*,'(100(F10.5))') mo_coef(i,:)
    enddo


    !do i = 1,15 
    !call in_mat_vec_index(i,p,q)   
    !print*,i,'p,q',p,q
    !enddo

    deallocate(v_grad,H,Hm1,Hm1_tmpr,gHm1,A,R)
 endif

end program
