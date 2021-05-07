program debug_grad
  implicit none

  !================================
  ! Program to verify the gradient
  !================================
  
  !===========
  ! Variables
  !===========

  double precision, allocatable :: v_grad(:),v_grad2(:)
  integer                       :: n
  integer                       :: i
  double precision :: max_error, threshold,max_elem
  integer :: nb_error
  ! v_grad : double precision vector of length n, gradient
  ! n      :  integer, n = mo_num*(mo_num-1)/2, number of orbital pairs (p,q) with p < q
  ! i,j,p,q,k : integer, indexes
  
  ! Def of n  
  n = mo_num*(mo_num-1)/2

  PROVIDE mo_two_e_integrals_in_map

  !============
  ! Allocation
  !============
 
  allocate(v_grad(n))
  allocate(v_grad2(n))

  !=============
  ! Calculation
  !=============

  call diagonalize_ci

  ! Gradient  
  call first_gradient(n,v_grad,max_elem)
  call gradient(n,v_grad2,max_elem)
  
  v_grad = v_grad - v_grad2
  nb_error = 0
  max_error = 0d0 
  threshold = 1d-12 

  do i = 1,n
    if (ABS(v_grad(i)) > threshold) then
       print*,v_grad(i)
       nb_error = nb_error + 1

       if (ABS(v_grad(i)) > max_error) then
         max_error = v_grad(i)
       endif

    endif
  enddo
 
  print*,''
  print*,'Check the gradient' 
  print*,'Threshold :', threshold
  print*,'Nb error :', nb_error
  print*,'Max error :', max_error

  deallocate(v_grad,v_grad2)

end program
