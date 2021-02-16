program orb_opt
  implicit none
  
  double precision, allocatable :: grad(:,:),R(:,:)
  double precision, allocatable :: H(:,:),e_val(:),work(:,:)
  double precision, allocatable :: Hm1(:,:),v_grad(:),gHm1(:),A(:,:),vec(:),Hm1_tmpr(:,:)
  integer :: info,method,n,i,j,lwork
  
  double precision :: angle
  
  read(*,*) angle
  print*,"angle", angle
  method = 1 
  
  n = mo_num*(mo_num-1)/2
 
  !============
  ! Allocation
  !============
 
  allocate(grad(mo_num,mo_num),R(mo_num,mo_num))
  allocate(H(n,n),Hm1(n,n),v_grad(n),gHm1(n),A(mo_num,mo_num),vec(n))
  allocate(Hm1_tmpr(n,n))

  !=============
  ! Calculation
  !=============
  
  if (method == 0) then
    call gradient(grad)
    
    !call dm_antisym(grad,mo_num,mo_num,info)
    call dm_rotation(grad,mo_num,R,mo_num,mo_num,info)

    call dm_newton_test(R) 
    print*,'Grad : ', grad(1,2)
    deallocate(grad,R)
  
  else 
    
    call gradient(grad)
    call hess(H,n)
    
    lwork=3*n-1
    allocate(work(lwork,n),e_val(n))
    
    call dsyev('V','U',n,H,size(H,1),e_val,work,lwork,info)
    if (info /= 0) then
            call ABORT
    endif
    
    do j=1,n
            do i=1,n
                    if ((i==j) .and. (ABS(e_val(i))>10d0**(-10))) then
                           Hm1(i,j)=1d0/e_val(i)
                    else
                           Hm1(i,j)=0d0
                    endif
            enddo
    enddo
    
    deallocate(work,e_val)
    
    call dgemm('N','T',n,n,n,1d0,Hm1,size(Hm1,1),H,size(H,1),0d0,Hm1_tmpr,size(Hm1_tmpr,1))
    !print*,'Hm1', Hm1_tmpr(:,:) 
    call dgemm('N','N',n,n,n,1d0,H,size(H,1),Hm1_tmpr,size(Hm1_tmpr,1),0d0,Hm1,size(Hm1,1))
    
    !print*,'grad'
    !print*, grad(:,:)
    !print*,'Hm1'
    !print*, Hm1(:,:)

    call dv_mat_to_vec(grad,size(grad,1),mo_num,v_grad,n,info)
    !print*,'vgrad',v_grad(:) 
    call dgemv('T',n,n,1d0,Hm1,size(Hm1,1),v_grad,1,0d0,gHm1,1)
    
    call dm_vec_to_mat(gHm1,n,A,mo_num,info)
    
    do i=1,mo_num
      do j=1,mo_num
        A(i,j) = angle
      enddo
    enddo
    call dm_antisym(A,mo_num,mo_num,info)
    call dm_rotation(A,mo_num,R,mo_num,mo_num,info)
    
    print*,'R'
    print*,R(:,:)
   
    call dm_newton_test(R)
    
    deallocate(grad,H,Hm1,Hm1_tmpr,gHm1,A,R)
 endif

end program
