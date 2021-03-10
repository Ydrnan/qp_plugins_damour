subroutine hess(n,H)

  include 'constants.h' 

  implicit none

  !==================================================================
  ! Compute the hessian of energy with respects to orbital rotations
  !==================================================================

  !===========
  ! Variables  
  !===========

  ! in
  integer, intent(in)           :: n 
  !n         : integer, n = mo_num*(mo_num-1)/2
  
  ! out
  double precision, intent(out) :: H(n,n)
  ! H        : n by n double precision matrix containing the 2D hessian
 
  ! internal
  double precision, allocatable :: hessian(:,:,:,:), h_tmpr(:,:,:,:)
  double precision, allocatable :: H_test(:,:)
  integer                       :: p,q
  integer                       :: r,s,t,u,v
  integer                       :: pq,rs
  integer                       :: istate
  double precision              :: t1,t2,t3,t4,t5,t6
  ! hessian  : mo_num 4D double precision matrix containing the hessian before the permutations
  ! h_tmpr   : mo_num 4D double precision matrix containing the hessian after the permutations
  ! H_test   : monum**2 by mo_num**2 double precision matrix to debug the H matrix
  ! p,q,r,s  : integer, indexes of the 4D hessian matrix
  ! t,u,v    : integer, indexes to compute hessian elements
  ! pq,rs    : integer, indexes for the conversion from 4D to 2D hessian matrix
  ! istate   : integer, electronic state
  ! t1,t2,t3 : double precision, t3 = t2 - t1, time to compute the hessian 

  double precision, allocatable :: tmp_bi_int_3(:,:,:), tmp_2rdm_3(:,:,:)
  double precision, allocatable :: tmp_accu(:,:), tmp_accu_sym(:,:)

  ! Funtion 
  double precision              :: get_two_e_integral
  ! get_two_e_integral :  double precision function, two e integrals 

  double precision :: ddot
  
  ! Provided :
  ! mo_one_e_integrals : mono e- integrals
  ! get_two_e_integral : two e- integrals
  ! one_e_dm_mo_alpha, one_e_dm_mo_beta : one body density matrix
  ! two_e_dm_mo : two body density matrix

  !============
  ! Allocation
  !============

  allocate(hessian(mo_num,mo_num,mo_num,mo_num),h_tmpr(mo_num,mo_num,mo_num,mo_num))
  allocate(H_test(mo_num**2,mo_num**2))

  !=============
  ! Calculation
  !=============

  ! Initialization
  hessian = 0d0

  !do istate = 1, N_states
  istate = 1

  ! From Anderson et. al. (2014) 
  ! The Journal of Chemical Physics 141, 244104 (2014); doi: 10.1063/1.4904384

  CALL CPU_TIME(t1)

  ! First line, first term
  !  do p = 1, mo_num
  !    do q = 1, mo_num
  !      do r = 1, mo_num
  !        do s = 1, mo_num

  !            if (q==r) then
  !              do u = 1, mo_num

  !                hessian(p,q,r,s) = hessian(p,q,r,s) + 0.5d0 * (  &
  !                  mo_one_e_integrals(u,p) * (one_e_dm_mo_alpha(u,s,istate) + one_e_dm_mo_beta(u,s,istate)) &
  !                + mo_one_e_integrals(s,u) * (one_e_dm_mo_alpha(p,u,istate) + one_e_dm_mo_beta(p,u,istate)))

  !              enddo
  !            endif

  !        enddo
  !      enddo
  !    enddo
  !  enddo
 
 ! Opt First line, first term
 double precision, allocatable :: one_e_rdm_mo_y(:,:)
 allocate(one_e_rdm_mo_y(mo_num,mo_num))
 allocate(tmp_accu(mo_num,mo_num), tmp_accu_sym(mo_num,mo_num))

 do s = 1, mo_num
   do p = 1, mo_num
     one_e_rdm_mo_y(p,s) = one_e_dm_mo_alpha(p,s,istate) + one_e_dm_mo_beta(p,s,istate)
   enddo
 enddo

!   do s = 1, mo_num
!     do p = 1, mo_num
!       tmp_accu(p,s) = 0d0
!       do u = 1, mo_num
!
!          tmp_accu(p,s) = tmp_accu(p,s) +  &
!            mo_one_e_integrals(u,p) * one_e_dm_mo(u,s) 
!
!        enddo
!     enddo
!   enddo
 
  call dgemm('T','N',mo_num,mo_num,mo_num,1d0,mo_one_e_integrals,&
  size(mo_one_e_integrals,1), one_e_rdm_mo_y, size(one_e_rdm_mo_y,1),&
  0d0, tmp_accu, size(tmp_accu,1))
 
   do s = 1, mo_num
     do p = 1, mo_num
       tmp_accu_sym(p,s) = 0.5d0 * (tmp_accu(p,s) + tmp_accu(s,p))
     enddo
   enddo 

   do s = 1, mo_num
      do p = 1, mo_num
         do r = 1, mo_num
           hessian(p,r,r,s) = hessian(p,r,r,s) + tmp_accu_sym(p,s)
         enddo
      enddo
   enddo

 deallocate(tmp_accu,tmp_accu_sym)

 ! First line, second term
 !    do p = 1, mo_num
 !      do q = 1, mo_num
 !        do r = 1, mo_num
 !          do s = 1, mo_num

 !            if (p==s) then
 !              do u = 1, mo_num

 !                    hessian(p,q,r,s) = hessian(p,q,r,s) + 0.5d0 * ( &
 !                      mo_one_e_integrals(u,r) * (one_e_dm_mo_alpha(u,q,istate) + one_e_dm_mo_beta(u,q,istate)) &
 !                    + mo_one_e_integrals(q,u) * (one_e_dm_mo_alpha(r,u,istate) + one_e_dm_mo_beta(r,u,istate)))
 !              enddo
 !            endif

 !          enddo
 !        enddo
 !      enddo
 !    enddo
 
 ! Opt First line, second term
  allocate(tmp_accu(mo_num,mo_num), tmp_accu_sym(mo_num,mo_num))


     !  do q = 1, mo_num
     !    do r = 1, mo_num
     !      tmp_accu(q,r) = 0d0

     !          do u = 1, mo_num

     !                tmp_accu(q,r) = tmp_accu(q,r) + 0.5d0 * ( &
     !                  mo_one_e_integrals(u,r) * one_e_rdm_mo_y(u,q)  &
     !                + mo_one_e_integrals(q,u) * one_e_rdm_mo_y(r,u))
     !          enddo
     !    enddo
     !  enddo
  call dgemm('T','N',mo_num,mo_num,mo_num,1d0,mo_one_e_integrals,&
  size(mo_one_e_integrals,1), one_e_rdm_mo_y, size(one_e_rdm_mo_y,1),&
  0d0, tmp_accu, size(tmp_accu,1))

   do r = 1, mo_num
     do q = 1, mo_num
       tmp_accu_sym(r,q) = 0.5d0 * (tmp_accu(q,r) + tmp_accu(r,q))
     enddo
   enddo

   do r = 1, mo_num
      do q = 1, mo_num
         do s = 1, mo_num
           hessian(s,q,r,s) = hessian(s,q,r,s) + tmp_accu_sym(q,r)
         enddo
      enddo
   enddo

 deallocate(tmp_accu,tmp_accu_sym)

 


 ! First line, third term
 !   do p = 1, mo_num
 !     do q = 1, mo_num
 !       do r = 1, mo_num
 !         do s = 1, mo_num

 !           hessian(p,q,r,s) = hessian(p,q,r,s) &
 !           - mo_one_e_integrals(s,p) * (one_e_dm_mo_alpha(r,q,istate) + one_e_dm_mo_beta(r,q,istate)) &
 !           - mo_one_e_integrals(q,r) * (one_e_dm_mo_alpha(p,s,istate) + one_e_dm_mo_beta(p,s,istate))

 !         enddo
 !       enddo
 !     enddo
 !   enddo

 ! Opt First line, third term
    do p = 1, mo_num
      do q = 1, mo_num
        do r = 1, mo_num
          do s = 1, mo_num

            hessian(p,q,r,s) = hessian(p,q,r,s) &
            - mo_one_e_integrals(s,p) * one_e_rdm_mo_y(r,q) &
            - mo_one_e_integrals(r,q) * one_e_rdm_mo_y(s,p)

          enddo
        enddo
      enddo
    enddo


 ! Second line, first term
 !    do p = 1, mo_num
 !      do q = 1, mo_num
 !        do r = 1, mo_num
 !          do s = 1, mo_num

 !             if (q==r) then
 !               do t = 1, mo_num
 !                 do u = 1, mo_num
 !                   do v = 1, mo_num

 !                     hessian(p,q,r,s) = hessian(p,q,r,s) + 0.5d0 * (  &
 !                       get_two_e_integral(u,v,p,t,mo_integrals_map) * two_e_dm_mo(u,v,s,t,istate) &
 !                     + get_two_e_integral(s,t,u,v,mo_integrals_map) * two_e_dm_mo(p,t,u,v,istate))

 !                   enddo
 !                 enddo
 !               enddo
 !             endif

 !          enddo
 !        enddo
 !      enddo
 !    enddo

 ! Opt Second line, first term
  allocate(tmp_bi_int_3(mo_num,mo_num,mo_num))
  allocate(tmp_2rdm_3(mo_num,mo_num,mo_num))
  allocate(tmp_accu(mo_num,mo_num),tmp_accu_sym(mo_num,mo_num))
 call CPU_TIME(t4)
    do t = 1, mo_num
 
      do p = 1, mo_num
        do v = 1, mo_num
          do u = 1, mo_num
   
            tmp_bi_int_3(u,v,p) = get_two_e_integral(u,v,p,t,mo_integrals_map)              
 
          enddo
        enddo
      enddo
  
      do p = 1, mo_num
        do v = 1, mo_num
          do u = 1, mo_num
   
            tmp_2rdm_3(u,v,p) = two_e_dm_mo(u,v,p,t,istate)  
 
          enddo
        enddo
      enddo
 
      call dgemm('T','N', mo_num, mo_num, mo_num*mo_num, 1.d0, &
      tmp_bi_int_3   , mo_num*mo_num, tmp_2rdm_3, mo_num*mo_num, &
      0.d0, tmp_accu, size(tmp_accu,1))
 
      do p = 1, mo_num
        do s = 1, mo_num
          tmp_accu_sym(p,s) = 0.5d0 * (tmp_accu(p,s)+tmp_accu(s,p))
        enddo
      enddo
 
      do s = 1, mo_num
        do r = 1, mo_num
          do p = 1, mo_num
 
            hessian(p,r,r,s) = hessian(p,r,r,s) + tmp_accu_sym(p,s) 
 
          enddo
        enddo
      enddo
 
    enddo
 call CPU_TIME(t5)
 t6=t5-t4
 print*,'t6', t6
 deallocate(tmp_bi_int_3,tmp_2rdm_3) 
 deallocate(tmp_accu,tmp_accu_sym) 

 ! Second line, second term
!     do p = 1, mo_num
!       do q = 1, mo_num
!         do r = 1, mo_num
!           do s = 1, mo_num
!
!             if (p==s) then
!               do t = 1, mo_num
!                 do u = 1, mo_num
!                   do v = 1, mo_num
!
!                     hessian(p,q,r,s) = hessian(p,q,r,s) + 0.5d0 * ( &
!                       get_two_e_integral(q,t,u,v,mo_integrals_map) * two_e_dm_mo(r,t,u,v,istate) &
!                     + get_two_e_integral(u,v,r,t,mo_integrals_map) * two_e_dm_mo(u,v,q,t,istate))
!
!                   enddo
!                 enddo
!               enddo
!             endif
!
!           enddo
!         enddo
!       enddo
!     enddo

 ! Opt Second line, second term
 allocate(tmp_bi_int_3(mo_num,mo_num,mo_num), tmp_2rdm_3(mo_num,mo_num,mo_num))
 allocate(tmp_accu(mo_num,mo_num),tmp_accu_sym(mo_num,mo_num))

     do t = 1, mo_num
       
       do r = 1, mo_num
         do v = 1, mo_num
           do u = 1, mo_num
              tmp_bi_int_3(u,v,r) = get_two_e_integral(u,v,r,t,mo_integrals_map)
           enddo
         enddo
       enddo
  
       do r = 1, mo_num
         do v = 1, mo_num
           do u = 1, mo_num
              tmp_2rdm_3(u,v,r) = two_e_dm_mo(u,v,r,t,istate)
           enddo
         enddo
       enddo

       call dgemm('T','N', mo_num, mo_num, mo_num*mo_num, 1.d0, &
       tmp_bi_int_3 , mo_num*mo_num, tmp_2rdm_3, mo_num*mo_num, &
       0.d0, tmp_accu, size(tmp_accu,1))


       do r = 1, mo_num
         do q = 1, mo_num
           tmp_accu_sym(q,r) = 0.5d0 * (tmp_accu(q,r) + tmp_accu(r,q))
         enddo
       enddo
       
       do r = 1, mo_num
         do q = 1, mo_num
           do s = 1, mo_num
             hessian(s,q,r,s) = hessian(s,q,r,s) + tmp_accu_sym(q,r)
           enddo
         enddo
       enddo
     enddo

  deallocate(tmp_bi_int_3,tmp_2rdm_3)
  deallocate(tmp_accu,tmp_accu_sym)

 ! Third line, first term
 !   do p = 1, mo_num
 !     do q = 1, mo_num
 !       do r = 1, mo_num
 !         do s = 1, mo_num

 !           do u = 1, mo_num
 !             do v = 1, mo_num

 !               hessian(p,q,r,s) = hessian(p,q,r,s) &
 !                + get_two_e_integral(u,v,p,r,mo_integrals_map) * two_e_dm_mo(u,v,q,s,istate) &
 !                + get_two_e_integral(q,s,u,v,mo_integrals_map) * two_e_dm_mo(p,r,u,v,istate)

 !             enddo
 !           enddo

 !         enddo
 !       enddo
 !     enddo
 !   enddo

  ! Opt Third line, first term
!  double precision, allocatable :: tmp_bi_int_2(:,:), tmp_2rdm_2(:,:)
!  allocate(tmp_bi_int_2(mo_num,mo_num), tmp_2rdm_2(mo_num,mo_num))
!  call CPU_TIME(t4)
!
!  do v = 1, mo_num
!    do u = 1, mo_num
!
!      do r = 1, mo_num
!        do p = 1, mo_num  
!          tmp_bi_int_2(p,r) = get_two_e_integral(p,r,u,v,mo_integrals_map)
!        enddo
!      enddo
!
!      do r = 1, mo_num
!        do p = 1, mo_num 
!          tmp_2rdm_2(p,r) = two_e_dm_mo(p,r,u,v,istate)
!        enddo
!      enddo
!
!    do s = 1, mo_num
!      do r = 1, mo_num
!        do q = 1, mo_num
!          do p = 1, mo_num
!
!              hessian(p,q,r,s) = hessian(p,q,r,s) &
!               + tmp_bi_int_2(p,r) * tmp_2rdm_2(q,s) &
!               + tmp_bi_int_2(q,s) * tmp_2rdm_2(p,r)
!         
!          enddo
!        enddo
!      enddo
!    enddo
!
!    enddo
!  enddo
!  call CPU_TIME(t5)
!  t6 = t5 - t4
!  print*,'t',t6
!  deallocate(tmp_bi_int_2,tmp_2rdm_2)
 
 allocate(tmp_bi_int_3(mo_num,mo_num,mo_num), tmp_2rdm_3(mo_num,mo_num,mo_num))
 call CPU_TIME(t4)
 do v = 1, mo_num
  
  do u = 1, mo_num
    do r = 1, mo_num
     do p = 1, mo_num
       tmp_bi_int_3(u,r,p) = get_two_e_integral(p,r,u,v,mo_integrals_map)
     enddo
   enddo
  enddo

  do u = 1, mo_num
    do r = 1, mo_num
     do p = 1, mo_num
       tmp_2rdm_3(p,r,u) = two_e_dm_mo(p,r,u,v,istate)
     enddo
   enddo
  enddo

   do u = 1, mo_num
     do s = 1, mo_num
       do r = 1, mo_num
         do q = 1, mo_num
           do p = 1, mo_num

                 hessian(p,q,r,s) = hessian(p,q,r,s) &
                  + tmp_bi_int_3(p,r,u) * tmp_2rdm_3(q,s,u) &
                  + tmp_bi_int_3(q,s,u) * tmp_2rdm_3(p,r,u)

           enddo
         enddo
       enddo
     enddo
   enddo
 enddo
 call CPU_TIME(t5)
 t6 = t5 - t4
 print*,'t', t6
 deallocate(tmp_bi_int_3,tmp_2rdm_3)

 ! Third line, second term
 !   do p = 1, mo_num
 !     do q = 1, mo_num
 !       do r = 1, mo_num
 !         do s = 1, mo_num

 !           do t = 1, mo_num
 !             do u = 1, mo_num

 !               hessian(p,q,r,s) = hessian(p,q,r,s) &
 !                - get_two_e_integral(s,t,p,u,mo_integrals_map) * two_e_dm_mo(r,t,q,u,istate) &
 !                - get_two_e_integral(t,s,p,u,mo_integrals_map) * two_e_dm_mo(t,r,q,u,istate) &
 !                - get_two_e_integral(q,u,r,t,mo_integrals_map) * two_e_dm_mo(p,u,s,t,istate) &
 !                - get_two_e_integral(q,u,t,r,mo_integrals_map) * two_e_dm_mo(p,u,t,s,istate)

 !             enddo
 !           enddo

 !         enddo
 !       enddo
 !     enddo
 !   enddo
   
! Opt Third line, second term
  allocate(tmp_bi_int_3(mo_num,mo_num,mo_num), tmp_2rdm_3(mo_num,mo_num,mo_num))
  call CPU_TIME(t4)

  do u = 1, mo_num

    do p = 1, mo_num
      do t = 1, mo_num
        do s = 1, mo_num
          tmp_bi_int_3(s,t,p) = get_two_e_integral(s,t,p,u,mo_integrals_map)
        enddo
      enddo
    enddo

    do p = 1, mo_num
      do t = 1, mo_num
        do s = 1, mo_num
          tmp_2rdm_3(s,t,p) = two_e_dm_mo(s,t,p,u,istate) 
        enddo
      enddo
    enddo

    do p = 1, mo_num
      do q = 1, mo_num
        do r = 1, mo_num
          do s = 1, mo_num

            do t = 1, mo_num

                hessian(p,q,r,s) = hessian(p,q,r,s) &
                 - tmp_bi_int_3(s,t,p) * tmp_2rdm_3(r,t,q) &
                 - tmp_bi_int_3(t,s,p) * tmp_2rdm_3(t,r,q) &
                 - tmp_bi_int_3(r,t,q) * tmp_2rdm_3(s,t,p) &
                 - tmp_bi_int_3(t,r,q) * tmp_2rdm_3(t,s,p)

            enddo

          enddo
        enddo
      enddo
    enddo
  enddo
  call CPU_TIME(t5)
  t6 = t5-t4
  print*,'t',t6
  deallocate(tmp_bi_int_3,tmp_2rdm_3)

  !enddo

  CALL CPU_TIME(t2)
  t3 = t2 -t1
  print*,'Time to compute the hessian : ', t3

  !===========
  ! 2D matrix
  !===========

  ! Convert the hessian mo_num * mo_num * mo_num * mo_num matrix in a
  ! 2D n * n matrix (n = mo_num*(mo_num-1)/2)
  ! H(pq,rs) : p<q and r<s
  ! Hessian(p,q,r,s) = P_pq P_rs [ ...]
  ! => Hessian(p,q,r,s) = (p,q,r,s) - (q,p,r,s) - (p,q,s,r) + (q,p,s,r)

  ! Permutations 
  do r = 1, mo_num
    do s = 1, mo_num
      do q = 1, mo_num
        do p = 1, mo_num

          h_tmpr(p,q,r,s) = (hessian(p,q,r,s) - hessian(q,p,r,s) - hessian(p,q,s,r) + hessian(q,p,s,r))

        enddo
      enddo
    enddo
  enddo

  ! Debug, all the elements of the 4D hessian in a mo_num**2 by mo_num**2 matrix
  if (debug) then
    pq=0
    rs=0
    do p=1,mo_num
      do q = 1,mo_num
      pq=pq+1
      rs=0
        do r = 1, mo_num
          do s = 1, mo_num
          rs = rs+1 
           H_test(pq,rs) = h_tmpr(p,q,r,s)
          enddo
        enddo
      enddo
    enddo
  
    print*,'mo_num**2 by mo_num**2 hessian matrix'
    do pq=1,mo_num**2
      write(*,'(100(F10.5))') H_test(pq,:)
    enddo
  endif

  ! Debug, eigenvalues of the 4D hessian to compare with an other code
  !double precision :: e_val(mo_num**2),H_v(mo_num**2,mo_num**2), H_u(mo_num,mo_num,mo_num,mo_num)
  !call lapack_diag(e_val,H_v,H_test,mo_num**2,mo_num**2)

  !print*,'Eigenvalues of the 4D hessian as a mo_num**2 by mo_num**2 matrix :'
  !write(*,'(100(F10.5))') e_val(:) 
  
  !deallocate(e_val,H_v,H_u)

  ! 4D mo_num matrix to 2D n matrix
  do pq = 1, n
    call in_mat_vec_index(pq,p,q)
    do rs = 1, n
      call in_mat_vec_index(rs,r,s)
      H(pq,rs) = h_tmpr(p,q,r,s)   
    enddo
  enddo


  ! Display
  if (debug) then 
    print*,'2D Hessian matrix'
    do pq = 1, n
      write(*,'(100(F10.5))') H(pq,:)
    enddo 
  endif

  !==============
  ! Deallocation
  !==============

  deallocate(hessian,h_tmpr,H_test)

end subroutine
