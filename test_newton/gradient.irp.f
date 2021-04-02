subroutine gradient(n,v_grad)

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
  double precision, intent(out) :: v_grad(n)
  ! v_grad : double precision vector of length n containeing the gradient

  !==========
  ! internal
  !==========
  double precision, allocatable :: grad(:,:),A(:,:)
  double precision              :: norm
  integer                       :: i,p,q,r,s,t
  integer                       :: istate
  double precision              :: t1,t2,t3,t4,t5,t6
  ! grad : double precision matrix containing the gradient before the permutation
  ! A : double precision matrix containing the gradient after the permutation
  ! norm : double precision number, the norm of the vector gradient
  ! i,p,q,r,s,t : integer, indexes 
  ! istate : integer, the electronic state
  ! t1,t2,t3 : t3 = t2 - t1, time to compute the gradient
  ! t4,t5,t6 : t6 = t5 - t4, time to compute each element

  !double precision, allocatable :: one_e_rdm_mo_y(:,:)
  ! one_e_rdm_mo_y : mo_num 2D double precision matrix containing the one e density matrix,
  !                  compute as the sum of one_e_dm_mo_alpha and one_e_dm_mo_beta 

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
  !allocate(one_e_rdm_mo_y(mo_num,mo_num))
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
  v_grad = 0d0
  grad = 0d0
 
  ! Electronic state
  !!do istate = 1, N_states
  istate = 1
 
   ! do q = 1, mo_num
   !   do p = 1, mo_num

   !      one_e_rdm_mo_y(p,q) = one_e_dm_mo_alpha(p,q,istate) + one_e_dm_mo_beta(p,q,istate)

   !   enddo
   ! enddo

    ! From Anderson et. al. (2014) 
    ! The Journal of Chemical Physics 141, 244104 (2014); doi: 10.1063/1.4904384

    !!! first term 

    !do p = 1, mo_num
    !    do q = 1, mo_num
    !       grad(p,q) = 0d0
    !       do r = 1, mo_num
    !        
    !          grad(p,q) = grad(p,q) &
    !                +2d0 * mo_one_e_integrals(p,r) * one_e_rdm_mo_y(r,q) !&
    !                !- mo_one_e_integrals(r,q) * one_e_rdm_mo_y(p,r)
    !
    !      enddo
    !   enddo
    !enddo

    !!! Opt first term
  
    CALL CPU_TIME(t4)
    tmp_accu = 0d0
  
    call dgemm('N','N',mo_num,mo_num,mo_num,2d0,mo_one_e_integrals,&
    mo_num,one_e_dm_mo,mo_num,0d0,tmp_accu,mo_num)
  
    do q = 1, mo_num
      do p = 1, mo_num

        grad(p,q) = grad(p,q) + tmp_accu(p,q)

      enddo
    enddo 
    
    CALL CPU_TIME(t5)
    t6 = t5-t4
    print*,'Gradient, first term (s) :', t6 
  
    !!!!! Second term
    
    !do p = 1, mo_num
    !  do q = 1, mo_num 
    !    do r = 1, mo_num
    !      do s = 1, mo_num
    !        do t= 1, mo_num
    !
    !        grad(p,q) = grad(p,q) &
    !                + get_two_e_integral(p,t,r,s,mo_integrals_map) * two_e_dm_mo(r,s,q,t,1) &
    !                - get_two_e_integral(r,s,q,t,mo_integrals_map) * two_e_dm_mo(p,t,r,s,1)
    !       enddo
    !      enddo
    !    enddo
    !  enddo
    !enddo
  
    !!! Opt second term  
    
    CALL CPU_TIME(t4)
  
    tmp_accu = 0d0
  
    do t = 1, mo_num
     
      do p = 1, mo_num
        do s = 1, mo_num
          do r = 1, mo_num
              
            tmp_bi_int_3(r,s,p) = 2d0 * get_two_e_integral(r,s,p,t,mo_integrals_map)
           
          enddo
        enddo
      enddo
      
      do q = 1, mo_num
        do s = 1, mo_num
          do r = 1, mo_num
               
             tmp_2rdm_3(r,s,q) = two_e_dm_mo(r,s,q,t,1)
    
          enddo
        enddo
      enddo
  
  !     tmp_accu = 0d0
  !     do p = 1, mo_num
  !      do q = 1, mo_num
  !        do s = 1, mo_num
  !          do r = 1, mo_num
  !   
  !          tmp_accu(p,q) = tmp_accu(p,q) + tmp_bi_int_3(r,s,p) * tmp_2rdm_3(r,s,q)
  !         !grad(p,q) = grad(p,q) + tmp_bi_int_3(r,s,p) * tmp_2rdm_3(r,s,q)   
  !         
  !          enddo
  !        enddo
  !      enddo
  !    enddo
     
      call dgemm('T','N',mo_num,mo_num,mo_num*mo_num,1d0,tmp_bi_int_3,&
        mo_num*mo_num,tmp_2rdm_3,mo_num*mo_num,0d0,tmp_accu,mo_num)
     
      do q = 1, mo_num
        do p = 1, mo_num
  
          grad(p,q) = grad(p,q) + tmp_accu(p,q)
  
        enddo
      enddo
  
    enddo
    
    CALL CPU_TIME(t5)
    t6 = t5-t4
    print*,'Gradient second term (s) : ', t6

 ! Debug ocaml
!  print*,'two e rdm'
!  do p=1,mo_num
!    do q=1,mo_num
!      do r=1,mo_num
!        do s=1,mo_num
!          print*,p,q,r,s,two_e_dm_mo(p,q,r,s,1)
!        enddo
!      enddo
!    enddo
!  enddo
!
!  print*,'bi int'
!  do p=1,mo_num
!    do q=1,mo_num
!      do r=1,mo_num
!        do s=1,mo_num
!          if (ABS(get_two_e_integral(p,q,r,s,mo_integrals_map))>1.d-14) then
!            print*,p,q,r,s, get_two_e_integral(p,q,r,s,mo_integrals_map)
!          else
!            print*,p,q,r,s,0d0
!          endif
!        enddo
!      enddo
!    enddo
!  enddo
!
!  print*,'mono int'
!  do p=1,mo_num
!    print*,mo_one_e_integrals(p,:)
!  enddo
!
!  print*, 'one e rdm'
!  do p=1,mo_num
!    print*, (one_e_dm_mo_alpha(p,:,istate) + one_e_dm_mo_beta(p,:,istate))
!  enddo

 !Ecriture des intégrales dans un fichier pour le lire avec OCaml
 ! Debug pour moi 

  if (ocaml) then
    open(unit=10,file='../../../../../../App_y/miniconda3/Work_yann/one_e_dm.dat')
    do p = 1, mo_num
      do q = 1, mo_num
        write(10,*) p, q, (one_e_dm_mo_alpha(p,q,istate) + one_e_dm_mo_beta(p,q,istate))
      enddo
    enddo
    close(10)
  
    open(unit=11,file='../../../../../../App_y/miniconda3/Work_yann/one_e_int.dat')
    do p = 1, mo_num
      do q = 1, mo_num
        write(11,*) p, q, mo_one_e_integrals(p,q)
      enddo
    enddo
    close(11)
    
    open(unit=12,file='../../../../../../App_y/miniconda3/Work_yann/two_e_int.dat')
    do p = 1, mo_num
      do q = 1, mo_num
        do r = 1, mo_num
          do s = 1, mo_num
            write(12,*) p, q, r, s, get_two_e_integral(p,q,r,s,mo_integrals_map)
          enddo
        enddo
      enddo
    enddo
    close(12)
    
    open(unit=13,file='../../../../../../App_y/miniconda3/Work_yann/two_e_dm.dat')
    do p = 1, mo_num
      do q = 1, mo_num
        do r = 1, mo_num
          do s = 1, mo_num
            write(13,*) p, q, r, s, two_e_dm_mo(p,q,r,s,1)
          enddo
        enddo
      enddo
    enddo
    close(13)
  endif

 ! Debug ocaml
!  double precision, allocatable :: two_e_integrals(:,:,:,:)
!  allocate(two_e_integrals(mo_num,mo_num,mo_num,mo_num))
!
!  do p=1,mo_num
!    do q=1,mo_num
!      do r=1,mo_num
!        do s=1,mo_num
!          two_e_integrals(p,q,r,s) = get_two_e_integral(p,q,r,s,mo_integrals_map)
!        enddo
!      enddo
!    enddo
!  enddo
!
!
!  print*,'two e int'
!  do p = 1,mo_num
!    do q= 1, mo_num
!      write(*,'(100(F10.5))') two_e_integrals(p,q,:,:)
!    enddo
!  enddo
!
!  print*,'two e rdm'
!  do p = 1,mo_num
!    do q= 1, mo_num
!      write(*,'(100(F10.5))') two_e_dm_mo(p,q,:,:,1)
!    enddo
!  enddo

  ! Conversion mo_num*mo_num matrix to mo_num(mo_num-1)/2 vector

  do i=1,n
    call in_mat_vec_index(i,p,q)
    v_grad(i)=(grad(p,q) - grad(q,p))
  enddo  

  ! Display, vector containing the gradient elements 

  if (debug) then  
    print*,'Vector containing the gradient :'
    write(*,'(100(F10.5))') v_grad(1:n)
  endif  

  ! Norm of the vector

  norm = dnrm2(v_grad)
  print*, 'Gradient norm : ', norm

  ! Matrix gradient

  A = 0d0
  do q=1,mo_num
    do p=1,mo_num
      A(p,q) = grad(p,q) - grad(q,p)
    enddo
  enddo

  ! Display, matrix containting the gradient elements

  if (debug) then
    print*,'Matrix containing the gradient :'
    do i = 1, mo_num
      write(*,'(100(F10.5))') A(i,1:mo_num)
    enddo
  endif

  !==============
  ! Deallocation
  !==============

  deallocate(grad,A,tmp_bi_int_3,tmp_2rdm_3)
  deallocate(tmp_accu)!,one_e_rdm_mo_y)

  if (debug) then
    print*,'Leaves gradient'
  endif

end subroutine
