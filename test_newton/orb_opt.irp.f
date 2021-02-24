program orb_opt
  implicit none
  
  double precision, allocatable :: grad(:,:),R(:,:)
  double precision, allocatable :: H(:,:),e_val(:),work(:,:)
  double precision, allocatable :: Hm1(:,:),v_grad(:),gHm1(:),A(:,:),vec(:),Hm1_tmpr(:,:)
  integer :: info,method,n,i,j,lwork
  !double precision, allocatable ::  Aexp(:,:),R_tmpr(:,:), H_diag(:,:),H_tmpr(:,:)
  
  double precision :: angle, norm, normH
  
  method = 1 
  
  n = mo_num*(mo_num-1)/2
  !n = mo_num**2

  !============
  ! Allocation
  !============
 
  !allocate(H_diag(n,n),H_tmpr(n,n))
  !allocate(Aexp(mo_num,mo_num),R_tmpr(mo_num,mo_num))
  allocate(v_grad(n),R(mo_num,mo_num))
  allocate(H(n,n),Hm1(n,n),gHm1(n),A(mo_num,mo_num))  
  allocate(Hm1_tmpr(n,n))

  !=============
  ! Calculation
  !=============
  
  if (method == 0) then
    
    call gradient(n,v_grad)
    
    print*, 'Norm : ', norm 
    !print*, 'v_grad' 
    !print*, v_grad(:)

    allocate(grad(mo_num,mo_num))

    call dm_vec_to_mat(v_grad,size(v_grad,1),grad,size(grad,1),info)  
  
    call dm_antisym(grad,mo_num,mo_num,info)
    call dm_rotation(grad,mo_num,R,mo_num,mo_num,info)

    !call dm_newton_test(R) 

    deallocate(grad,R)
  
  else 
    
    call gradient(n,v_grad)
    
    !norm = norm2(v_grad)
    !print*, 'Norm : ', norm
    
    call hess(n,H)
   
    !normH = norm2(H)
    !print*, 'NormH : ', normH

    ! Inversion de H 
    lwork=3*n-1
    allocate(work(lwork,n),e_val(n))
    
    call dsyev('V','U',n,H,size(H,1),e_val,work,lwork,info)
    if (info /= 0) then
        call ABORT
    endif
   
    Hm1=0d0 
    !H_diag=0d0
    do i=1,n
        !print*,'H_val',e_val(i)
        !H_diag(i,i)=e_val(i)
        if ( (ABS(e_val(i))>1.d-7)) then
            Hm1(i,i)=1d0/e_val(i)
        else
            Hm1(i,i)=0d0
        endif
    enddo
    
    deallocate(work,e_val)
    
    call dgemm('N','T',n,n,n,1d0,Hm1,size(Hm1,1),H,size(H,1),0d0,Hm1_tmpr,size(Hm1_tmpr,1))
    call dgemm('N','N',n,n,n,1d0,H,size(H,1),Hm1_tmpr,size(Hm1_tmpr,1),0d0,Hm1,size(Hm1,1))
    ! Fin de l'inversion

!test     
!    call dgemm('N','T',n,n,n,1d0,H_diag,size(H_diag,1),H,size(H,1),0d0,H_tmpr,size(H_tmpr,1))
!    call dgemm('N','N',n,n,n,1d0,H,size(H,1),H_tmpr,size(H_tmpr,1),0d0,H_diag,size(H_diag,1))
!
!    print*,'H test'
!    do i=1,n
!      print*,H_diag(i,:)
!    enddo

    v_grad = 0.5d0 * v_grad
    
    call dgemv('T',n,n,1d0,Hm1,size(Hm1,1),v_grad,1,0d0,gHm1,1)
    
    call dm_vec_to_mat(gHm1,n,A,mo_num,info)
    
    ! print matrice a exponentialiser
    print*,'matrice A'  
    do i=1,mo_num
      print*,A(i,:)
    enddo

!test 
    !do i=1,mo_num
    !  do j=1,mo_num
    !    A(i,j) = -0.0017d0!angle
    !  enddo
    !enddo
    
     A=0d0
    A(1,3) = -0.0054456d0
    A(2,4) = 1.87241d0
    A(3,1) = 0.0054456d0
    A(4,2) = -1.87241d0

    print*,'A'
    do i=1,mo_num
      print*,A(i,:)
    enddo


    call dm_antisym(A,mo_num,mo_num,info)
    call dm_rotation(A,mo_num,R,mo_num,mo_num,info)
 
!test
!    R=transpose(R)
!    ! exp 
!    lwork=3*mo_num-1
!    allocate(work(lwork,mo_num),e_val(mo_num))
!
!    call dsyev('V','U',mo_num,A,size(A,1),e_val,work,lwork,info)
!    if (info /= 0) then
!        call ABORT
!    endif
!
!    Aexp=0d0
!    do i=1,mo_num
!      Aexp(i,i)=exp(e_val(i))
!    enddo
!
!    deallocate(work,e_val)
!
!    call dgemm('N','T',mo_num,mo_num,mo_num,1d0,Aexp,size(Aexp,1),A,size(A,1),0d0,R_tmpr,size(R_tmpr,1))
!    call dgemm('N','N',mo_num,mo_num,mo_num,1d0,A,size(A,1),R_tmpr,size(R_tmpr,1),0d0,R,size(R,1))
!    ! Fin de l'exp
!   
!    print*,'R'
!    do i=1,mo_num
!    print*,R(i,:)
!    enddo   

!    call dm_newton_test(R)
    
    deallocate(v_grad,H,Hm1,Hm1_tmpr,gHm1,A,R)
 endif

end program
