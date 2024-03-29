program debug_hessian
  implicit none

  !======================================
  ! Program to verify the hessian matrix
  !======================================
  
  !===========
  ! Variables
  !===========

  double precision, allocatable :: H(:,:),H2(:,:), h_f(:,:,:,:), h_f2(:,:,:,:)
  integer                       :: method
  integer                       :: n
  integer                       :: i,j,k,l
  double precision :: max_error,threshold,max_error_H
  integer :: nb_error, nb_error_H
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

  PROVIDE mo_two_e_integrals_in_map

  !============
  ! Allocation
  !============
 
  allocate(H(n,n),H2(n,n))  
  allocate(h_f(mo_num,mo_num,mo_num,mo_num),h_f2(mo_num,mo_num,mo_num,mo_num))

  !=============
  ! Calculation
  !=============

  
  ! Hessian and norm
  if (method == 1) then 

    print*,'Use the full hessian matrix'
    call org_first_hess(n,H,h_f)
    call hess(n,H2,h_f2)

    !do i = 1, n
    !  print*,i,H(i,i)
    !enddo

    h_f =h_f - h_f2
    H = H - H2
    max_error = 0d0
    nb_error = 0    
    threshold = 1d-12

    do l = 1, mo_num
      do k= 1, mo_num
        do j = 1, mo_num
          do i = 1, mo_num

            if (ABS(h_f(i,j,k,l)) > threshold) then

              print*,i,j,k,l,h_f(i,j,k,l),h_f2(i,j,k,l)
              nb_error = nb_error + 1

              if (ABS(h_f(i,j,k,l)) > ABS(max_error)) then
                max_error = h_f(i,j,k,l)
              endif

            endif

          enddo
        enddo
      enddo
    enddo

   max_error_H = 0d0
   nb_error_H = 0

   do j = 1, n
     do i = 1, n
       if (ABS(H(i,j)) > threshold) then
         print*, H(i,j)
         nb_error_H = nb_error_H + 1

         if (ABS(H(i,j)) > ABS(max_error_H)) then
           max_error_H = H(i,j)
         endif

       endif
     enddo
   enddo 

  else

    print*, 'Use the diagonal hessian matrix'
    call first_diag_hess(n,H,h_f)
    call diag_hess(n,H2,h_f2)
    
    do i = 1, n
     print*,i,H(i,i)
    enddo

    h_f = h_f - h_f2
    max_error = 0d0
    nb_error = 0
    threshold = 1d-12

    do l = 1, mo_num
      do k = 1, mo_num
        do j = 1, mo_num
          do i = 1, mo_num

            if (ABS(h_f(i,j,k,l)) > threshold) then

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

    h=H-H2
  
    max_error_H = 0d0
    nb_error_H = 0
 
    do j = 1, n
      do i = 1, n
        if (ABS(H(i,j)) > threshold) then
          print*, H(i,j)
          nb_error_H = nb_error_H + 1
 
          if (ABS(H(i,j)) > ABS(max_error_H)) then
            max_error_H = H(i,j)
          endif
 
        endif
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
  print*,''
  print*,'Nb error_H :', nb_error_H
  print*,'Max error_H :', max_error_H
 
  deallocate(H,H2,h_f,h_f2)

end program
