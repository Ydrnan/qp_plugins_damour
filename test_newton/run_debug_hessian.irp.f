program run_debug_hessian
  implicit none

  !======================================
  ! Program to verify the hessian matrix
  !======================================
  
  !===========
  ! Variables
  !===========

  double precision, allocatable :: H(:,:),H1(:,:), h_f(:,:,:,:), h_f2(:,:,:,:)
  integer                       :: method
  integer                       :: n
  integer                       :: i,j,k,l
  double precision :: max_error,threshold
  integer :: nb_error
  ! H      : n by n double precision matrix, Hessian matrix
  ! method : 
  !          - 1 : Full hessian
  !          - 2 : Diagonal hessian
  ! n      :  integer, n = mo_num*(mo_num-1)/2, number of orbital pairs (p,q) with p < q
  ! i,j,p,q,k : integer, indexes
  
  ! Choice of the method 
  method = 1 ! 1 -> full hessian, 2 -> diagonal hessian
 
  ! Def of n  
  n = mo_num*(mo_num-1)/2

  !============
  ! Allocation
  !============
 
  allocate(H(n,n),H1(n,n))  
  allocate(h_f(mo_num,mo_num,mo_num,mo_num),h_f2(mo_num,mo_num,mo_num,mo_num))

  !=============
  ! Calculation
  !=============

  
  ! Hessian and norm
  if (method == 1) then 

    print*,'Use the full hessian matrix'
    call first_hess(n,H,h_f)
    !call hess(n,H1,h_f2)
    call old_hess(n,H1,h_f2) 
 
    h_f = h_f - h_f2
    max_error = 0d0
    nb_error = 0    
    threshold = 1d-12

    do i = 1, mo_num
      do j= 1, mo_num
        do k = 1, mo_num
          do l = 1, mo_num

            if (h_f(i,j,k,l) > threshold) then

              print*,h_f(i,j,k,l)
              nb_error = nb_error + 1

              if (ABS(h_f(i,j,k,l)) > ABS(max_error)) then
                max_error = h_f(i,j,k,l)
              endif

            endif

          enddo
        enddo
      enddo
    enddo

  else

    print*, 'Use the diagonal hessian matrix'
    call first_diag_hess(n,H,h_f)
    call diag_hess(n,H1,h_f2)
    
    h_f = h_f - h_f2
    max_error = 0d0
    nb_error = 0
    threshold = 1d-12

    do i = 1, mo_num
      do j= 1, mo_num
        do k = 1, mo_num
          do l = 1, mo_num

            if (h_f(i,j,k,l) > threshold) then

              print*,h_f(i,j,k,l)
              nb_error = nb_error + 1

              if (ABS(h_f(i,j,k,l)) > ABS(max_error)) then
                max_error = h_f(i,j,k,l)
              endif

            endif

          enddo
        enddo
      enddo
    enddo

  endif
  
  print*,''
  if (method == 1) then
    print*,'Check the full hessian'
  else
    print*,'Check the diagonal hessian'
  endif
   
  print*,'Threshold :', threshold
  print*,'Nb error :', nb_error
  print*,'Max error :', max_error
 
  deallocate(H,H1,h_f,h_f2)

end program
