subroutine gradient(n,v_grad,max_elem)
  use omp_lib
  include 'constants.h'

  implicit none

  !===================================================================
  ! Compute the gradient of energy with respects to orbital rotations
  !===================================================================

  ! Check if read_wf = true, else :
  ! qp set determinant read_wf true

  END_DOC

  !====
  ! in
  !====
  integer, intent(in) :: n
  ! n : integer, n = mo_num*(mo_num-1)/2
  
  !=====
  ! out
  !=====
  double precision, intent(out) :: v_grad(n), max_elem
  ! v_grad : double precision vector of length n containeing the gradient

  !==========
  ! internal
  !==========
  double precision, allocatable :: grad(:,:),A(:,:)
  double precision              :: norm
  integer                       :: i,p,q,r,s,t
  double precision              :: t1,t2,t3,t4,t5,t6
  ! grad : double precision matrix containing the gradient before the permutation
  ! A : double precision matrix containing the gradient after the permutation
  ! norm : double precision number, the norm of the vector gradient
  ! i,p,q,r,s,t : integer, indexes 
  ! istate : integer, the electronic state
  ! t1,t2,t3 : t3 = t2 - t1, time to compute the gradient
  ! t4,t5,t6 : t6 = t5 - t4, time to compute each element

  double precision, allocatable :: tmp_accu(:,:)
  double precision, allocatable :: tmp_bi_int_3(:,:,:), tmp_2rdm_3(:,:,:)
  ! tmp_bi_int_3 : mo_num 3D double precision matrix containinig the bi electronic
  !                integrals with 1 fix index
  ! tmp_2rdm_3   : mo_num 3D double precision matrix containinig the 2 body reduce
  !                density matrix with 1 fix index
  ! tmp_accu     : mo_num 2D double precision matrix, temporary matrix

  !===========
  ! Functions
  !===========
  double precision :: get_two_e_integral, dnrm2
  ! get_two_e_integral :  double precision function that gives the two e integrals
  ! norm2 : double precision function that gives the norm of a vector
 
  ! Provided :
  ! mo_one_e_integrals : mono e- integrals
  ! get_two_e_integral : two e- integrals
  ! one_e_dm_mo_alpha, one_e_dm_mo_beta : one body density matrix
  ! two_e_dm_mo : two body density matrix

  !============
  ! Allocation
  !============

  allocate(grad(mo_num,mo_num),A(mo_num,mo_num))
  
  ! From Anderson et. al. (2014) 
  ! The Journal of Chemical Physics 141, 244104 (2014); doi: 10.1063/1.4904384 

  ! LaTeX formula

  ! \begin{align*}
  ! G_{pq}&= \dfrac{\partial E(x)}{\partial x_{pq}} \\ 
  ! &= \mathcal{P}_{pq} \left[ \sum_r ( h_p^r \gamma_r^q - h_r^q \gamma_p^r) 
  !+ \sum_{rst} (v_{pt}^{rs} \Gamma_{rs}^{qt} - v_{rs}^{qt} \Gamma_{pt}^{rs}) \right]
  ! \end{align*} 

 ! Initialization
  call omp_set_max_active_levels(1)

  !$OMP PARALLEL                                                 &
      !$OMP PRIVATE(                                             &
      !$OMP   p,q,r,s,t,                                         &
      !$OMP   tmp_accu, tmp_bi_int_3, tmp_2rdm_3)                &
      !$OMP SHARED(grad, one_e_dm_mo, mo_num,mo_one_e_integrals, &
      !$OMP mo_integrals_map,t4,t5,t6)                           &
      !$OMP DEFAULT(SHARED)
 
  !==============================
  ! Allocation of private arrays
  !==============================

  allocate(tmp_accu(mo_num,mo_num))
  allocate(tmp_bi_int_3(mo_num,mo_num,mo_num))
  allocate(tmp_2rdm_3(mo_num,mo_num,mo_num))
 
  !=============
  ! Calculation
  !============= 

  if (debug) then
    print*,'Enter in gradient'
  endif

  ! Initialization

  !$OMP DO
  do q = 1, mo_num
    do p = 1,mo_num
      grad(p,q) = 0d0
    enddo
  enddo
  !$OMP END DO
  
  !========
  ! Term 1
  !======== 

  !do p = 1, mo_num
  !    do q = 1, mo_num
  !       do r = 1, mo_num
  !        
  !          grad(p,q) = grad(p,q) &
  !                + mo_one_e_integrals(p,r) * one_e_rdm_mo_y(r,q) &
  !                - mo_one_e_integrals(r,q) * one_e_rdm_mo_y(p,r)
  !
  !      enddo
  !   enddo
  !enddo

  !****************
  ! Opt first term
  !****************

  !$OMP MASTER
  CALL wall_TIME(t4)
  !$OMP END MASTER

  call dgemm('N','N',mo_num,mo_num,mo_num,1d0,mo_one_e_integrals,&
  mo_num,one_e_dm_mo,mo_num,0d0,tmp_accu,mo_num)
  
  !$OMP DO
  do q = 1, mo_num
    do p = 1, mo_num

      grad(p,q) = grad(p,q) + (tmp_accu(p,q) - tmp_accu(q,p))

    enddo
  enddo 
  !$OMP END DO
  
  !$OMP MASTER
  CALL wall_TIME(t5)
  t6 = t5-t4
  print*,'Gradient, first term (s) :', t6 
  !$OMP END MASTER

  !========
  ! Term 2
  !========
  
  !do p = 1, mo_num
  !  do q = 1, mo_num 
  !    do r = 1, mo_num
  !      do s = 1, mo_num
  !        do t= 1, mo_num
  !
  !        grad(p,q) = grad(p,q) &
  !                + get_two_e_integral(p,t,r,s,mo_integrals_map) * two_e_dm_mo(r,s,q,t) &
  !                - get_two_e_integral(r,s,q,t,mo_integrals_map) * two_e_dm_mo(p,t,r,s)
  !       enddo
  !      enddo
  !    enddo
  !  enddo
  !enddo

  !*****************
  ! Opt second term  
  !*****************

  !$OMP MASTER
  CALL wall_TIME(t4)
  !$OMP END MASTER 
 
  !$OMP DO
  do t = 1, mo_num
    
    do p = 1, mo_num
      do s = 1, mo_num
        do r = 1, mo_num
            
          tmp_bi_int_3(r,s,p) = get_two_e_integral(r,s,p,t,mo_integrals_map)
         
        enddo
      enddo
    enddo

    do q = 1, mo_num
      do s = 1, mo_num
        do r = 1, mo_num
             
           tmp_2rdm_3(r,s,q) = two_e_dm_mo(r,s,q,t)
  
        enddo
      enddo
    enddo

    call dgemm('T','N',mo_num,mo_num,mo_num*mo_num,1d0,tmp_bi_int_3,&
      mo_num*mo_num,tmp_2rdm_3,mo_num*mo_num,0d0,tmp_accu,mo_num)

    !$OMP CRITICAL   
    do q = 1, mo_num
      do p = 1, mo_num

        grad(p,q) = grad(p,q) + tmp_accu(p,q) - tmp_accu(q,p)

      enddo
    enddo
    !$OMP END CRITICAL

  enddo
  !$OMP END DO

  !$OMP MASTER
  CALL wall_TIME(t5)
  t6 = t5-t4
  print*,'Gradient second term (s) : ', t6
  !$OMP END MASTER  

  !================================
  ! Deallocation of private arrays
  !================================

  deallocate(tmp_bi_int_3,tmp_2rdm_3,tmp_accu)

  !$OMP END PARALLEL

  call omp_set_max_active_levels(4)

  !=====================
  ! 2D matrix to vector
  !=====================

  ! Conversion mo_num*mo_num matrix to mo_num(mo_num-1)/2 vector
  do i=1,n
    call vec_to_mat_index(i,p,q)
    v_grad(i)=(grad(p,q) - grad(q,p))
  enddo  

  ! Display, vector containing the gradient elements 
  if (debug) then  
    print*,'Vector containing the gradient :'
    write(*,'(100(F10.5))') v_grad(1:n)
  endif  

  !======
  ! Norm
  !======

  ! Norm of the vector
  norm = dnrm2(n,v_grad,1)
  print*, 'Gradient norm : ', norm

  !=============
  ! Max element
  !=============

  ! Max element of the gradient
  max_elem = 0d0
  do i = 1, n
    if (ABS(v_grad(i)) > ABS(max_elem)) then
      max_elem = v_grad(i)
    endif
  enddo

  print*,'Max element in gardient :', max_elem  

  ! Display, matrix containting the gradient elements
  if (debug) then
    ! Matrix gradient
    A = 0d0
    do q=1,mo_num
      do p=1,mo_num
        A(p,q) = grad(p,q) - grad(q,p)
      enddo
    enddo

    print*,'Matrix containing the gradient :'
    do i = 1, mo_num
      write(*,'(100(F10.5))') A(i,1:mo_num)
    enddo
  endif

  !==============
  ! Deallocation
  !==============

  deallocate(grad,A)

  if (debug) then
    print*,'Leave gradient'
  endif

end subroutine
