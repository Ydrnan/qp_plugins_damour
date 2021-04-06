subroutine diag_hess(n,H, h_tmpr)
 
  use omp_lib
  
  include 'constants.h' 

  implicit none

  !===========================================================================
  ! Compute the diagonal hessian of energy with respects to orbital rotations
  !===========================================================================
  
  ! By diagonal hessian we mean, diagonal elements of the hessian  

  !===========
  ! Variables 
  !===========
 
  !====
  ! in
  !====
  integer, intent(in)           :: n 
  ! n        : integer, n = mo_num*(mo_num-1)/2
 
  !=====
  ! out
  !=====
  double precision, intent(out) :: H(n,n),  h_tmpr(mo_num,mo_num,mo_num,mo_num)
  ! H        : n by n double precision matrix containing the 2D hessian
  
  !==========
  ! internal
  !==========
  double precision, allocatable :: hessian(:,:,:,:)!, h_tmpr(:,:,:,:)
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
  ! t4,t5,t6 : double precision, t6 = t5 - t4, time to compute each term

  !===========
  ! temporary
  !===========
  double precision :: accu
  ! accu : tmp double precision variable to sum elements

  double precision, allocatable :: tmp_bi_int_3(:,:,:),tmp_bi_int_3_shared(:,:,:)
  double precision, allocatable :: tmp_2rdm_3(:,:,:),tmp_2rdm_3_shared(:,:,:)
  double precision, allocatable :: tmp_accu(:,:), tmp_accu_1(:), tmp_accu_shared(:,:)
  ! tmp_bi_int_2 : mo_num 2D double precision matrix containing the bi electronic integrals
  !                with 2 fix indexes
  ! tmp_bi_int_3 : mo_num 3D double precision matrix containing the bi electronic integrals
  !                with 1 fix index
  ! tmp_2rdm_2   : mo_num 2D double precision matrix containing the 2 body reduce density
  !                matrix with 2 fix indexes
  ! tmp_2rdm_3   : mo_num 3D double precision matrix containing the 2 body reduce density
  !                matrix with 1 fix index
  ! tmp_accu     : mo_num 2D double precision matrix temporary matrix
  ! tmp_accu_1   : mo_num 1D double precision matrix temporary matrix
  
  double precision, allocatable :: tmp_h_pppp(:), tmp_h_pqpq(:,:), tmp_h_pqqp(:,:)
  ! tmp_h_pppp : mo_num 1D double precision matrix containing the hessien elements hessian(p,p,p,p)
  ! tmp_h_pqpq : mo_num 2D double precision matrix containing the hessien elements hessian(p,q,p,q)
  ! tmp_h_pqqp : mo_num 2D double precision matrix containing the hessien elements hessian(p,q,q,p)
  

  !double precision, allocatable :: one_e_rdm_mo_y(:,:)
  ! one_e_rdm_mo_y : mo_num 2D double precision matrix containing the one e density matrix,
  !                  compute as the sum of one_e_dm_mo_alpha and one_e_dm_mo_beta
 
  ! Function
  double precision :: get_two_e_integral
  ! get_two_e_integral : double precision function, two e integrals
   
  double precision :: ddot
  ! ddot :: double precision Blas function, dot product
 
  ! Provided :
  ! mo_one_e_integrals : mono e- integrals
  ! get_two_e_integral : two e- integrals
  ! one_e_dm_mo_alpha, one_e_dm_mo_beta : one body density matrix
  ! two_e_dm_mo : two body density matrix

  !============
  ! Allocation
  !============

  allocate(hessian(mo_num,mo_num,mo_num,mo_num))!,h_tmpr(mo_num,mo_num,mo_num,mo_num))
  allocate(H_test(mo_num**2,mo_num**2))
  allocate(tmp_h_pppp(mo_num),tmp_h_pqpq(mo_num,mo_num),tmp_h_pqqp(mo_num,mo_num))
  allocate(tmp_accu(mo_num,mo_num),tmp_accu_shared(mo_num,mo_num))
  allocate(tmp_2rdm_3(mo_num,mo_num,mo_num),tmp_2rdm_3_shared(mo_num,mo_num,mo_num))
  allocate(tmp_bi_int_3(mo_num,mo_num,mo_num),tmp_bi_int_3_shared(mo_num,mo_num,mo_num))
  allocate(tmp_accu_1(mo_num))

  !=============
  ! Calculation
  !=============

  if (debug) then
    print*,'Enter in diag_hess'
  endif

  ! Initialization

   call omp_set_max_active_levels(1)

  !$OMP PARALLEL                                                     &
      !$OMP PRIVATE(                                                 &
      !$OMP   p,q,r,s, tmp_accu, tmp_accu_1,                         &
      !$OMP   u,v,t, tmp_bi_int_3, tmp_2rdm_3)                       &
      !$OMP SHARED(hessian, tmp_h_pppp, tmp_h_pqpq, tmp_h_pqqp,      &
      !$OMP  mo_num, mo_one_e_integrals, one_e_dm_mo,    &
      !$OMP  tmp_bi_int_3_shared, tmp_2rdm_3_shared,tmp_accu_shared,      &
      !$OMP two_e_dm_mo,mo_integrals_map,t1,t2,t3,t4,t5,t6) &
      !$OMP DEFAULT(NONE)

  !$OMP DO
  do s = 1,mo_num
    do r = 1, mo_num
      do q = 1, mo_num
        do p = 1, mo_num
          hessian(p,q,r,s) = 0d0
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO
  
  !$OMP MASTER
  do p = 1, mo_num
    tmp_h_pppp(p) = 0d0
  enddo
  !$OMP END MASTER

  !$OMP DO
  do q = 1, mo_num
    do p = 1, mo_num
      tmp_h_pqpq(p,q) = 0d0
    enddo
  enddo
  !$OMP END DO
   
  !$OMP DO
  do q = 1, mo_num
    do p = 1, mo_num
      tmp_h_pqqp(p,q) = 0d0
    enddo
  enddo
  !$OMP END DO
 
    ! From Anderson et. al. (2014) 
    ! The Journal of Chemical Physics 141, 244104 (2014); doi: 10.1063/1.4904384

    !$OMP MASTER
    CALL wall_TIME(t1)
    !$OMP END MASTER

    !!!!! First line, first term
    ! do p = 1, mo_num
    !   do q = 1, mo_num
    !     do r = 1, mo_num
    !       do s = 1, mo_num

    !       ! Et oui il y a des permutations après, donc il faut prendre aussi les permutations !!!! 
    !       if (((p==r) .and. (q==s)) .or. ((q==r) .and. (p==s)) &
    !          .or. ((p==s) .and. (q==r))) then
    !        
    !           if (q==r) then
    !             do u = 1, mo_num

    !               hessian(p,q,r,s) = hessian(p,q,r,s) + 0.5d0 * (  &
    !                 mo_one_e_integrals(u,p) * one_e_rdm_mo_y(u,s) &
    !               !mo_one_e_integrals(u,p) * (one_e_dm_mo_alpha(u,s,istate) + one_e_dm_mo_beta(u,s,istate)) &
    !               + mo_one_e_integrals(s,u) * one_e_rdm_mo_y(p,u))
    !               !+ mo_one_e_integrals(s,u) * (one_e_dm_mo_alpha(p,u,istate) + one_e_dm_mo_beta(p,u,istate)))

    !             enddo
    !           endif
    !         endif
    !       enddo
    !     enddo
    !   enddo
    ! enddo
  
    !!!!! Opt First line, first term

    !!! (p==r) .and. (q==s) .and. (q==r)
    !
    ! hessian(p,q,r,s) -> hessian(p,p,p,p)
    !
    !  0.5d0 * (  &
    !  mo_one_e_integrals(u,p) * one_e_rdm_mo_y(u,s) &
    !+ mo_one_e_integrals(s,u) * one_e_rdm_mo_y(p,u))
    !  =  
    !  0.5d0 * (  &
    !  mo_one_e_integrals(u,p) * one_e_rdm_mo_y(u,p) &
    !+ mo_one_e_integrals(p,u) * one_e_rdm_mo_y(p,u))
    ! =  
    !  mo_one_e_integrals(u,p) * one_e_rdm_mo_y(u,p) 

    !$OMP MASTER
    CALL wall_TIME(t4) 
    !$OMP END MASTER
  
    !$OMP MASTER
    do p = 1, mo_num
      tmp_accu_1(p) = 0d0
    enddo
    !$OMP END MASTER
  
    !$OMP MASTER
    do p = 1, mo_num
      do u = 1, mo_num
  
        tmp_accu_1(p) = tmp_accu_1(p) + mo_one_e_integrals(u,p) * one_e_dm_mo(u,p)
  
      enddo
    enddo
    !$OMP END MASTER
    
    !$OMP MASTER
    do p = 1, mo_num
      tmp_h_pppp(p) = tmp_h_pppp(p) + tmp_accu_1(p)
    enddo
    !$OMP END MASTER
 
    !!! (q==r) .and. (p==s) .and. (q==r)
    !
    ! hessian(p,q,r,s) -> hessian(p,q,q,p)
    !   
    !  0.5d0 * (  &
    !  mo_one_e_integrals(u,p) * one_e_rdm_mo_y(u,s) &
    !+ mo_one_e_integrals(s,u) * one_e_rdm_mo_y(p,u))
    !  =  
    !  0.5d0 * (  &
    !  mo_one_e_integrals(u,p) * one_e_rdm_mo_y(u,p) &
    !+ mo_one_e_integrals(p,u) * one_e_rdm_mo_y(p,u))
    ! =  
    !  mo_one_e_integrals(u,p) * one_e_rdm_mo_y(u,p)    

    !$OMP MASTER
    do p = 1, mo_num
      tmp_accu_1(p) = 0d0
    enddo
    !$OMP END MASTER

    !$OMP MASTER 
    do p = 1, mo_num
      do u = 1, mo_num

        tmp_accu_1(p) = tmp_accu_1(p) + mo_one_e_integrals(u,p) * one_e_dm_mo(u,p)

      enddo
    enddo
    !$OMP END MASTER
    
    !$OMP MASTER
    do q = 1, mo_num
      do p = 1, mo_num

        tmp_h_pqqp(p,q) = tmp_h_pqqp(p,q) + tmp_accu_1(p)

      enddo
    enddo
    !$OMP END MASTER

    !$OMP MASTER
    CALL wall_TIME(t5)
    t6= t5-t4
    print*,'l1 1',t6
    !$OMP END MASTER

    !!!!! First line, second term
    ! do p = 1, mo_num
    !   do q = 1, mo_num
    !     do r = 1, mo_num
    !       do s = 1, mo_num

    !       ! Et oui il y a des permutations après, donc il faut prendre aussi les permutations !!!! 
    !      if (((p==r) .and. (q==s)) .or. ((q==r) .and. (p==s)) &
    !         .or. ((p==s) .and. (q==r))) then
    !        
    !         if (p==s) then
    !           do u = 1, mo_num

    !                 hessian(p,q,r,s) = hessian(p,q,r,s) + 0.5d0 * ( &
    !                   mo_one_e_integrals(u,r) * one_e_rdm_mo_y(u,q) &
    !                 !mo_one_e_integrals(u,r) * (one_e_dm_mo_alpha(u,q,istate) + one_e_dm_mo_beta(u,q,istate)) &
    !                 + mo_one_e_integrals(q,u) * one_e_rdm_mo_y(r,u))
    !                 !+ mo_one_e_integrals(q,u) * (one_e_dm_mo_alpha(r,u,istate) + one_e_dm_mo_beta(r,u,istate)))
    !           enddo
    !         endif
    !         endif
    !       enddo
    !     enddo
    !   enddo
    ! enddo

    !!!!! Opt First line, second term

    !!! (p==r) .and. (q==s) .and. (p==s)
    !
    ! hessian(p,q,r,s) -> hessian(p,p,p,p)
    !
    ! 0.5d0 * (&
    !  mo_one_e_integrals(u,r) * one_e_rdm_mo_y(u,q) &
    !+ mo_one_e_integrals(q,u) * one_e_rdm_mo_y(r,u))
    ! =
    ! 0.5d0 * ( &
    ! mo_one_e_integrals(u,p) * one_e_rdm_mo_y(u,p) &
    !+ mo_one_e_integrals(p,u) * one_e_rdm_mo_y(p,u))
    ! = 
    ! mo_one_e_integrals(u,p) * one_e_rdm_mo_y(u,p)

    !$OMP MASTER
    CALL wall_TIME(t4)
    !$OMP END MASTER  

    !$OMP MASTER
    do p = 1, mo_num
      tmp_accu_1(p) = 0d0 
    enddo
    !$OMP END MASTER  

    !$OMP MASTER
    do p = 1, mo_num
      do u = 1, mo_num
 
        tmp_accu_1(p) = tmp_accu_1(p) +  mo_one_e_integrals(u,p) * one_e_dm_mo(u,p) 
 
      enddo
    enddo
    !$OMP END MASTER
    
    !$OMP MASTER
    do p = 1, mo_num

      tmp_h_pppp(p) = tmp_h_pppp(p) + tmp_accu_1(p)

    enddo
    !$OMP END MASTER

    !!! (q==r) .and. (p==s) .and. (p==s)
    !
    ! hessian(p,q,r,s) -> hessian(p,q,q,p)
    !
    ! 0.5d0 * (&
    !  mo_one_e_integrals(u,r) * one_e_rdm_mo_y(u,q) &
    !+ mo_one_e_integrals(q,u) * one_e_rdm_mo_y(r,u))
    ! =
    ! 0.5d0 * ( &
    ! mo_one_e_integrals(u,q) * one_e_rdm_mo_y(u,q) &
    !+ mo_one_e_integrals(q,u) * one_e_rdm_mo_y(q,u))
    ! = 
    ! mo_one_e_integrals(u,q) * one_e_rdm_mo_y(u,q)

    !$OMP MASTER
    do p = 1, mo_num
      tmp_accu_1(p) = 0d0
    enddo
    !$OMP END MASTER
 
    !$OMP MASTER
    do q = 1, mo_num
      do u = 1, mo_num
 
        tmp_accu_1(q) = tmp_accu_1(q) + mo_one_e_integrals(u,q) * one_e_dm_mo(u,q)
 
      enddo
    enddo
    !$OMP END MASTER

    !$OMP MASTER
    do q = 1, mo_num
      do p = 1, mo_num
 
        tmp_h_pqqp(p,q) = tmp_h_pqqp(p,q) + tmp_accu_1(q)
 
      enddo
    enddo
    !$OMP END MASTER

    !$OMP MASTER
    CALL wall_TIME(t5)
    t6= t5-t4
    print*,'l1 2',t6
    !$OMP END MASTER

    !!!!! First line, third term
    ! do p = 1, mo_num
    !   do q = 1, mo_num
    !     do r = 1, mo_num
    !       do s = 1, mo_num
 
    !       ! Et oui il y a des permutations après, donc il faut prendre aussi les permutations !!!! 
    !       if (((p==r) .and. (q==s)) .or. ((q==r) .and. (p==s)) &
    !          .or. ((p==s) .and. (q==r))) then
 
    !         hessian(p,q,r,s) = hessian(p,q,r,s) &
    !         !- mo_one_e_integrals(s,p) * (one_e_dm_mo_alpha(r,q,istate) + one_e_dm_mo_beta(r,q,istate)) &
    !         - mo_one_e_integrals(s,p) * one_e_rdm_mo_y(r,q) &
    !         !- mo_one_e_integrals(q,r) * (one_e_dm_mo_alpha(p,s,istate) + one_e_dm_mo_beta(p,s,istate))
    !         - mo_one_e_integrals(q,r) * one_e_rdm_mo_y(p,s)
 
    !        endif
    !       enddo
    !     enddo
    !   enddo
    ! enddo
 
    !!!!! Opt First line, third term
 
    !$OMP MASTER
    CALL wall_TIME(t4)
    !$OMP END MASTER 

    !!! (p==r) .and. (q==s)
    ! 
    ! hessian(p,q,r,s) -> hessian(p,q,p,q)
    ! 
    ! - mo_one_e_integrals(s,p) * one_e_rdm_mo_y(r,q) &
    ! - mo_one_e_integrals(q,r) * one_e_rdm_mo_y(p,s)
    ! =
    ! - mo_one_e_integrals(q,p) * one_e_rdm_mo_y(p,q) &
    ! - mo_one_e_integrals(q,p) * one_e_rdm_mo_y(p,q) 
    ! = 
    ! - 2d0 mo_one_e_integrals(q,p) * one_e_rdm_mo_y(p,q)
 
    !$OMP DO
    do q = 1, mo_num
      do p = 1, mo_num

        tmp_h_pqpq(p,q) = tmp_h_pqpq(p,q) &
        - 2d0 * mo_one_e_integrals(q,p) * one_e_dm_mo(p,q)

      enddo
    enddo
    !$OMP END DO
 
    !!! (q==r) .and. (p==s)
    ! 
    ! hessian(p,q,r,s) -> hessian(p,q,p,q)
    ! 
    ! - mo_one_e_integrals(s,p) * one_e_rdm_mo_y(r,q) &
    ! - mo_one_e_integrals(q,r) * one_e_rdm_mo_y(p,s)
    ! =
    ! - mo_one_e_integrals(q,p) * one_e_rdm_mo_y(p,q) &
    ! - mo_one_e_integrals(q,p) * one_e_rdm_mo_y(p,q) 
    ! = 
    ! - 2d0 mo_one_e_integrals(q,p) * one_e_rdm_mo_y(p,q)

    !$OMP DO
    do q = 1, mo_num
      do p = 1, mo_num

        tmp_h_pqqp(p,q) = tmp_h_pqqp(p,q) &
        - 2d0 * mo_one_e_integrals(p,p) * one_e_dm_mo(q,q)

      enddo
    enddo
    !$OMP END DO

    !$OMP MASTER
    CALL wall_TIME(t5)
    t6= t5-t4
    print*,'l1 3',t6
    !$OMP END MASTER

    !!!!! Second line, first term
    ! do p = 1, mo_num
    !   do q = 1, mo_num
    !      do r = 1, mo_num
    !        do s = 1, mo_num
 
    !        ! Et oui il y a des permutations après, donc il faut prendre aussi les permutations !!!! 
    !       if (((p==r) .and. (q==s)) .or. ((q==r) .and. (p==s)) &
    !          .or. ((p==s) .and. (q==r))) then
 
    !           if (q==r) then
    !             do t = 1, mo_num
    !               do u = 1, mo_num
    !                 do v = 1, mo_num
 
    !                   hessian(p,q,r,s) = hessian(p,q,r,s) + 0.5d0 * (  &
    !                     get_two_e_integral(u,v,p,t,mo_integrals_map) * two_e_dm_mo(u,v,s,t,1) &
    !                   + get_two_e_integral(s,t,u,v,mo_integrals_map) * two_e_dm_mo(p,t,u,v,1))
 
    !                 enddo
    !               enddo
    !             enddo
    !           endif
    !         endif
    !        enddo
    !      enddo
    !    enddo
    !  enddo

    !!!!! Opt Second line, first term
    
    ! (p==r) .and. (q==s) .and. (q==r)
    !
    ! hessian(p,q,r,s) -> hessian(p,q,p,q)
    ! 
    ! 0.5d0 * (  &
    !  get_two_e_integral(u,v,p,t,mo_integrals_map) * two_e_dm_mo(u,v,s,t,1) &
    !+ get_two_e_integral(s,t,u,v,mo_integrals_map) * two_e_dm_mo(p,t,u,v,1))
    ! =
    ! 0.5d0 * (  &
    !  get_two_e_integral(u,v,p,t,mo_integrals_map) * two_e_dm_mo(u,v,q,t,&) &
    !+ get_two_e_integral(q,t,u,v,mo_integrals_map) * two_e_dm_mo(p,t,u,v,&))
    ! = 
    ! get_two_e_integral(u,v,p,t,mo_integrals_map) * two_e_dm_mo(u,v,q,t,&)

    !$OMP MASTER
    CALL wall_TIME(t4)
    !$OMP END MASTER

    !$OMP MASTER
    do p = 1, mo_num
      tmp_accu_1(p) = 0d0
    enddo
    !$OMP END MASTER

    !$OMP MASTER
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

             tmp_2rdm_3(u,v,p) = two_e_dm_mo(u,v,p,t,1)

          enddo
        enddo
      enddo
   
      do p = 1, mo_num
        do v = 1, mo_num
          do u = 1, mo_num
  
            tmp_accu_1(p) = tmp_accu_1(p) &
            + tmp_bi_int_3(u,v,p) * tmp_2rdm_3(u,v,p) 

          enddo
        enddo
      enddo
    enddo
    !$OMP END MASTER
  
    !$OMP MASTER
    do p =1, mo_num
  
      tmp_h_pppp(p) = tmp_h_pppp(p) + tmp_accu_1(p)
  
    enddo
    !$OMP END MASTER
   
  !!! (q==r) .and. (p==s) .and. (q=r)
  !
  ! hessian(p,q,r,s) -> hessian(p,q,q,p)
  ! 
  ! 0.5d0 * (  &
  !  get_two_e_integral(u,v,p,t,mo_integrals_map) * two_e_dm_mo(u,v,s,t,1) &
  !+ get_two_e_integral(s,t,u,v,mo_integrals_map) * two_e_dm_mo(p,t,u,v,1))
  ! =
  ! 0.5d0 * (  &
  !  get_two_e_integral(u,v,p,t,mo_integrals_map) * two_e_dm_mo(u,v,p,t,1) &
  !+ get_two_e_integral(p,t,u,v,mo_integrals_map) * two_e_dm_mo(p,t,u,v,1))
  ! = 
  ! get_two_e_integral(u,v,p,t,mo_integrals_map) * two_e_dm_mo(u,v,p,t,1)

    !$OMP MASTER
    do p = 1, mo_num
      tmp_accu_1(p) = 0d0
    enddo
    !$OMP END MASTER
  
    !$OMP MASTER
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

             tmp_2rdm_3(u,v,p) = two_e_dm_mo(u,v,p,t,1)

          enddo
        enddo
      enddo

      do p = 1, mo_num
        do v = 1, mo_num
          do u = 1, mo_num
  
            tmp_accu_1(p) = tmp_accu_1(p) + & 
              tmp_bi_int_3(u,v,p) * tmp_2rdm_3(u,v,p)
 
          enddo
        enddo
      enddo
    enddo
    !$OMP END MASTER
  
    !$OMP MASTER
    do q = 1, mo_num
      do p = 1, mo_num
  
        tmp_h_pqqp(p,q) = tmp_h_pqqp(p,q) + tmp_accu_1(p) 
  
      enddo
    enddo
    !$OMP END MASTER
  
    !$OMP MASTER
    CALL wall_TIME(t5)
    t6 = t5-t4
    print*,'l2 1',t6
    !$OMP END MASTER
 
 !!!!! Second line, second term
    ! do p = 1, mo_num
    !   do q = 1, mo_num
    !     do r = 1, mo_num
    !       do s = 1, mo_num

    !       ! Et oui il y a des permutations après, donc il faut prendre aussi les permutations !!!! 
    !      if (((p==r) .and. (q==s)) .or. ((q==r) .and. (p==s)) &
    !         .or. ((p==s) .and. (q==r))) then 

    !        if (p==s) then
    !           do t = 1, mo_num
    !             do u = 1, mo_num
    !               do v = 1, mo_num

    !                 hessian(p,q,r,s) = hessian(p,q,r,s) + 0.5d0 * ( &
    !                   get_two_e_integral(q,t,u,v,mo_integrals_map) * two_e_dm_mo(r,t,u,v,1) &
    !                 + get_two_e_integral(u,v,r,t,mo_integrals_map) * two_e_dm_mo(u,v,q,t,1))

    !               enddo
    !             enddo
    !           enddo
    !         endif
    !        endif
    !       enddo
    !     enddo
    !   enddo
    ! enddo

  !!!!! Opt Second line, second term
  
  !!! (p==r) .and. (q==s) .and. (p==s)
  !
  ! hessian(p,q,r,s) -> hessian(p,p,p,p)
  !
  ! 0.5d0 * ( &
  !  get_two_e_integral(q,t,u,v,mo_integrals_map) * two_e_dm_mo(r,t,u,v,1) &
  !+ get_two_e_integral(u,v,r,t,mo_integrals_map) * two_e_dm_mo(u,v,q,t,1)) 
  ! = 
  ! 0.5d0 * ( &
  !  get_two_e_integral(p,t,u,v,mo_integrals_map) * two_e_dm_mo(p,t,u,v,1) &
  !+ get_two_e_integral(u,v,p,t,mo_integrals_map) * two_e_dm_mo(u,v,p,t,1))
  ! =
  ! get_two_e_integral(p,t,u,v,mo_integrals_map) * two_e_dm_mo(p,t,u,v,1)

      !$OMP MASTER
      CALL wall_TIME(t4)
      !$OMP END MASTER
   
      !$OMP MASTER
      do p = 1, mo_num
        tmp_accu_1(p) = 0d0
      enddo
      !$OMP END MASTER   

      !$OMP MASTER
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

               tmp_2rdm_3(u,v,p) = two_e_dm_mo(u,v,p,t,1)

            enddo
          enddo
        enddo
  
        do p = 1, mo_num
          do v = 1, mo_num
            do u = 1, mo_num
   
              tmp_accu_1(p) = tmp_accu_1(p) +&
                tmp_bi_int_3(u,v,p) * tmp_2rdm_3(u,v,p)
  
            enddo
          enddo
        enddo
      enddo
      !$OMP END MASTER
   
      !$OMP MASTER
      do p = 1, mo_num
   
        tmp_h_pppp(p) = tmp_h_pppp(p) + tmp_accu_1(p)
   
      enddo
      !$OMP END MASTER
   
    ! (q==r) .and. (p==s) .and. (p==s)
    !
    ! hessian(p,q,r,s) -> hessian(p,q,q,p)
    !
    ! 0.5d0 * ( &
    !  get_two_e_integral(q,t,u,v,mo_integrals_map) * two_e_dm_mo(r,t,u,v,1) &
    !+ get_two_e_integral(u,v,r,t,mo_integrals_map) * two_e_dm_mo(u,v,q,t,1)) 
    ! = 
    ! 0.5d0 * ( &
    !  get_two_e_integral(p,t,u,v,mo_integrals_map) * two_e_dm_mo(q,t,u,v,1) &
    !+ get_two_e_integral(u,v,q,t,mo_integrals_map) * two_e_dm_mo(u,v,p,t,1))
    ! =
    ! get_two_e_integral(q,t,u,v,mo_integrals_map) * two_e_dm_mo(q,t,u,v,1)

    !$OMP MASTER
    do p = 1,mo_num
      tmp_accu_1(p) = 0d0
    enddo
    !$OMP END MASTER
  
    !$OMP MASTER
    do t = 1, mo_num
  
      do q = 1, mo_num
        do v = 1, mo_num
          do u = 1, mo_num
  
             tmp_bi_int_3(u,v,q) = get_two_e_integral(u,v,q,t,mo_integrals_map)
  
          enddo
        enddo
      enddo
  
      do q = 1, mo_num
        do v = 1, mo_num
          do u = 1, mo_num
  
             tmp_2rdm_3(u,v,q) = two_e_dm_mo(u,v,q,t,1)
  
          enddo
        enddo
      enddo
      
      do q = 1, mo_num
        do v = 1, mo_num
          do u = 1, mo_num
   
            tmp_accu_1(q) = tmp_accu_1(q) +&
             tmp_bi_int_3(u,v,q) * tmp_2rdm_3(u,v,q)
  
          enddo
        enddo
      enddo
    enddo
    !$OMP END MASTER
    
    !$OMP MASTER
    do q = 1, mo_num
      do p = 1, mo_num
   
        tmp_h_pqqp(p,q) = tmp_h_pqqp(p,q) + tmp_accu_1(p)
   
      enddo
    enddo
    !$OMP END MASTER   
 
    !$OMP MASTER
    CALL wall_TIME(t5)
    t6 = t5-t4
    print*,'l2 2',t6
    !$OMP END MASTER

    !!!!! Third line, first term
    !   do p = 1, mo_num
    !     do q = 1, mo_num
    !      do r = 1, mo_num
    !         do s = 1, mo_num
   
    !         ! Et oui il y a des permutations après, donc il faut prendre aussi les permutations !!!! 
    !         if (((p==r) .and. (q==s)) .or. ((q==r) .and. (p==s)) &
    !            .or. ((p==s) .and. (q==r))) then
   
    !           do u = 1, mo_num
    !             do v = 1, mo_num
   
    !               hessian(p,q,r,s) = hessian(p,q,r,s) &
    !                + get_two_e_integral(u,v,p,r,mo_integrals_map) * two_e_dm_mo(u,v,q,s,1) &
    !                + get_two_e_integral(q,s,u,v,mo_integrals_map) * two_e_dm_mo(p,r,u,v,1)
   
    !             enddo
    !           enddo
    !          endif
   
    !         enddo
    !       enddo
    !     enddo
    !   enddo

    !!!!! Opt Third line, first term
    
    !!!(p==r) .and. (q==s)
    !
    ! hessian(p,q,r,s) -> hessian(p,q,p,q)
    ! 
    !  get_two_e_integral(u,v,p,r,mo_integrals_map) * two_e_dm_mo(u,v,q,s,1) &
    !+ get_two_e_integral(q,s,u,v,mo_integrals_map) * two_e_dm_mo(p,r,u,v,1) 
    ! = 
    !  get_two_e_integral(u,v,p,p,mo_integrals_map) * two_e_dm_mo(u,v,q,q,1) &
    !+ get_two_e_integral(q,q,u,v,mo_integrals_map) * two_e_dm_mo(p,p,u,v,1)
    ! = 
    ! 2d0 * get_two_e_integral(u,v,p,p,mo_integrals_map) * two_e_dm_mo(u,v,q,q,1)

    !$OMP MASTER
    CALL wall_TIME(t4)
    !$OMP END MASTER

    !$OMP DO
    do q = 1, mo_num
      do v = 1, mo_num
        do u = 1, mo_num

          tmp_2rdm_3_shared(u,v,q) = two_e_dm_mo(u,v,q,q,1)

        enddo
      enddo
    enddo
    !$OMP END DO 

    !$OMP DO
    do p = 1, mo_num
      do v = 1, mo_num
        do u = 1, mo_num

          tmp_bi_int_3_shared(u,v,p) = 2d0* get_two_e_integral(u,v,p,p,mo_integrals_map)

        enddo
      enddo
    enddo
    !$OMP END DO

  !  do p = 1, mo_num
  !    do q = 1, mo_num
  !      tmp_accu(p,q) = 0d0
  !      do v = 1, mo_num
  !        do u = 1, mo_num

  !          tmp_accu(p,q) = tmp_accu(p,q) + tmp_bi_int_3(u,v,p) * tmp_2rdm_3(u,v,q)

  !        enddo
  !      enddo
  !    enddo
  !  enddo
   
    call dgemm('T','N',mo_num,mo_num,mo_num*mo_num,1d0,tmp_bi_int_3_shared,&
    mo_num*mo_num,tmp_2rdm_3_shared,mo_num*mo_num,0d0,tmp_accu,mo_num)
 
    !$OMP DO
    do q = 1, mo_num
      do p = 1, mo_num

        tmp_h_pqpq(p,q) = tmp_h_pqpq(p,q) + tmp_accu(p,q)

      enddo
    enddo
    !$OMP END DO
  
    !(q==r) .and. (p==s)
    !
    ! hessian(p,q,r,s) -> hessian(p,q,q,p)
    ! 
    !  get_two_e_integral(u,v,p,r,mo_integrals_map) * two_e_dm_mo(u,v,q,s,1) &
    !+ get_two_e_integral(q,s,u,v,mo_integrals_map) * two_e_dm_mo(p,r,u,v,1) 
    ! = 
    !  get_two_e_integral(u,v,p,q,mo_integrals_map) * two_e_dm_mo(u,v,q,p,1) &
    !+ get_two_e_integral(q,p,u,v,mo_integrals_map) * two_e_dm_mo(p,q,u,v,1)
    ! = 
    ! 2d0 * get_two_e_integral(u,v,p,q,mo_integrals_map) * two_e_dm_mo(u,v,q,p,1)

    !$OMP DO
    do q = 1, mo_num
      
      do p = 1, mo_num
        do v = 1, mo_num
          do u = 1, mo_num

            tmp_bi_int_3(u,v,p) = 2d0* get_two_e_integral(u,v,q,p,mo_integrals_map)

          enddo
        enddo 
      enddo
      
      do p = 1, mo_num
        do v = 1, mo_num
          do u = 1, mo_num

            tmp_2rdm_3(u,v,p) = two_e_dm_mo(u,v,p,q,1)

          enddo
        enddo
      enddo
      
      do p = 1, mo_num
        do v = 1, mo_num
          do u = 1, mo_num
 
            tmp_h_pqqp(p,q) = tmp_h_pqqp(p,q) &
            + tmp_bi_int_3(u,v,p) * tmp_2rdm_3(u,v,p)

          enddo
        enddo
      enddo

    enddo 
    !$OMP END DO

    !$OMP MASTER
    CALL wall_TIME(t5)
    t6= t5-t4
    print*,'l3 1',t6
    !$OMP END MASTER

    !!!!! Third line, second term
    ! do p = 1, mo_num
    !   do q = 1, mo_num
    !     do r = 1, mo_num
    !       do s = 1, mo_num
    !       ! Et oui il y a des permutations après, donc il faut prendre aussi les permutations !!!! 
    !       if (((p==r) .and. (q==s)) .or. ((q==r) .and. (p==s)) &
    !          .or. ((p==s) .and. (q==r))) then
 
    !         do t = 1, mo_num
    !           do u = 1, mo_num
 
    !             hessian(p,q,r,s) = hessian(p,q,r,s) &
    !              - get_two_e_integral(s,t,p,u,mo_integrals_map) * two_e_dm_mo(r,t,q,u,1) &
    !              - get_two_e_integral(t,s,p,u,mo_integrals_map) * two_e_dm_mo(t,r,q,u,1) &
    !              - get_two_e_integral(q,u,r,t,mo_integrals_map) * two_e_dm_mo(p,u,s,t,1) &
    !              - get_two_e_integral(q,u,t,r,mo_integrals_map) * two_e_dm_mo(p,u,t,s,1)
 
    !           enddo
    !         enddo
 
    !       endif     
   
    !       enddo
    !     enddo
    !   enddo
    ! enddo

    !!!!! Opt Third line, second term

    !!! (p==r) .and. (q==s)
    ! 
    ! hessian(p,q,r,s) -> hessian(p,q,p,q)
    !
    ! - get_two_e_integral(s,t,p,u,mo_integrals_map) * two_e_dm_mo(r,t,q,u,1) &
    ! - get_two_e_integral(t,s,p,u,mo_integrals_map) * two_e_dm_mo(t,r,q,u,1) &
    ! - get_two_e_integral(q,u,r,t,mo_integrals_map) * two_e_dm_mo(p,u,s,t,1) &
    ! - get_two_e_integral(q,u,t,r,mo_integrals_map) * two_e_dm_mo(p,u,t,s,1)
    ! =
    ! - get_two_e_integral(q,t,p,u,mo_integrals_map) * two_e_dm_mo(p,t,q,u,1) &
    ! - get_two_e_integral(t,q,p,u,mo_integrals_map) * two_e_dm_mo(t,p,q,u,1) &
    ! - get_two_e_integral(q,u,p,t,mo_integrals_map) * two_e_dm_mo(p,u,q,t,1) &
    ! - get_two_e_integral(q,u,t,p,mo_integrals_map) * two_e_dm_mo(p,u,t,q,1)
    ! =
    ! - 2d0 * get_two_e_integral(q,t,p,u,mo_integrals_map) * two_e_dm_mo(p,t,q,u,1) &
    ! - 2d0 * get_two_e_integral(t,q,p,u,mo_integrals_map) * two_e_dm_mo(t,p,q,u,1) 
    
    !$OMP MASTER
    CALL wall_TIME(t4)
    !$OMP END MASTER

    !  part 1
    !$OMP MASTER
    do q = 1, mo_num
      do p = 1, mo_num
        tmp_accu(p,q) = 0d0
      enddo
    enddo 
    !$OMP END MASTER

    !$OMP MASTER
    do t = 1, mo_num

      do p = 1, mo_num
        do u = 1, mo_num
          do q = 1, mo_num

            tmp_bi_int_3(q,u,p) = 2d0 * get_two_e_integral(q,u,p,t,mo_integrals_map)

          enddo
        enddo
      enddo

      do p = 1, mo_num
        do u = 1, mo_num
          do q = 1, mo_num

             tmp_2rdm_3(q,u,p) = two_e_dm_mo(q,u,p,t,1)

          enddo
        enddo
      enddo

      do p = 1, mo_num
        do u = 1, mo_num
          do q = 1, mo_num
  
             tmp_accu(p,q) = tmp_accu(p,q) &
             - tmp_bi_int_3(q,u,p) * tmp_2rdm_3(q,u,p) 

          enddo
        enddo
      enddo

    enddo
    !$OMP END MASTER
    
    !$OMP MASTER
    do q = 1, mo_num
      do p = 1, mo_num
 
        tmp_h_pqpq(p,q) = tmp_h_pqpq(p,q) + tmp_accu(p,q)
 
      enddo
    enddo
    !$OMP END MASTER

   ! part 2 of the term
    
    !$OMP MASTER
    do q = 1, mo_num
      do p = 1, mo_num
        tmp_accu(p,q) = 0d0
      enddo
    enddo
    !$OMP END MASTER

    !$OMP MASTER
    do u = 1, mo_num
      
      do p = 1, mo_num
        do q = 1, mo_num
           do t = 1, mo_num

             tmp_bi_int_3(t,q,p) = 2d0*get_two_e_integral(t,q,p,u,mo_integrals_map)

           enddo
        enddo
      enddo

      do p= 1, mo_num
        do q = 1, mo_num
           do t = 1, mo_num
 
              tmp_2rdm_3(t,q,p) = two_e_dm_mo(t,p,q,u,1)
              
           enddo
        enddo
      enddo

   
      do q = 1, mo_num
        do p = 1, mo_num
          do t = 1, mo_num

             tmp_accu(p,q) = tmp_accu(p,q) &
             - tmp_bi_int_3(t,q,p) * tmp_2rdm_3(t,q,p)

          enddo
        enddo
      enddo

    enddo
    !$OMP END MASTER

    !$OMP MASTER
    do q = 1, mo_num
      do p = 1, mo_num

        tmp_h_pqpq(p,q) = tmp_h_pqpq(p,q) + tmp_accu(p,q)

      enddo
    enddo
    !$OMP END MASTER

    ! (q==r) .and. (p==s)
    !
    ! hessian(p,q,r,s) -> hessian(p,q,p,q)
    !
    ! - get_two_e_integral(s,t,p,u,mo_integrals_map) * two_e_dm_mo(r,t,q,u,1) &
    ! - get_two_e_integral(t,s,p,u,mo_integrals_map) * two_e_dm_mo(t,r,q,u,1) &
    ! - get_two_e_integral(q,u,r,t,mo_integrals_map) * two_e_dm_mo(p,u,s,t,1) &
    ! - get_two_e_integral(q,u,t,r,mo_integrals_map) * two_e_dm_mo(p,u,t,s,1)
    ! =
    ! - get_two_e_integral(p,t,p,u,mo_integrals_map) * two_e_dm_mo(q,t,q,u,1) &
    ! - get_two_e_integral(t,p,p,u,mo_integrals_map) * two_e_dm_mo(t,q,q,u,1) &
    ! - get_two_e_integral(q,u,q,t,mo_integrals_map) * two_e_dm_mo(p,u,p,t,1) &
    ! - get_two_e_integral(q,u,t,q,mo_integrals_map) * two_e_dm_mo(p,u,t,p,1)

!!     do q = 1, mo_num
!!       do u = 1, mo_num
!!         do t = 1, mo_num
!!     
!!           tmp_2rdm_3(t,u,q) = two_e_dm_mo(t,q,u,q,1)
!!
!!         enddo
!!       enddo
!!     enddo
!!
!!    do p = 1, mo_num
!!        
!!      do u = 1, mo_num
!!        do t = 1, mo_num
!!  
!!          tmp_bi_int_2(t,u) = 2d0* get_two_e_integral(t,p,u,p,mo_integrals_map)
!!
!!        enddo
!!      enddo
!!  
!!      do q = 1, mo_num 
!!        tmp_accu(q,p) = 0d0
!!        do u = 1, mo_num
!!          do t = 1, mo_num
!!  
!!            tmp_accu(q,p)  = tmp_accu(q,p) &
!!            - tmp_bi_int_2(t,u) * tmp_2rdm_3(t,u,q) 
!!
!!          enddo
!!        enddo
!!      enddo
!!    enddo
!!    
!!    do p = 1, mo_num
!!      do q = 1, mo_num
!!  
!!        tmp_h_pqqp(q,p) = tmp_h_pqqp(q,p) + tmp_accu(q,p)
!!  
!!      enddo
!!    enddo
  
    !!!!!!!!!!!!!!!!!!
    
    !$OMP DO
    do q = 1, mo_num
      do u = 1, mo_num
        do t = 1, mo_num

          tmp_2rdm_3_shared(t,u,q) = two_e_dm_mo(t,q,u,q,1)

        enddo
      enddo
    enddo
    !$OMP END DO
  
    !$OMP DO
    do p = 1, mo_num
      do u = 1, mo_num
        do t = 1, mo_num

          tmp_bi_int_3_shared(t,u,p) = 2d0* get_two_e_integral(t,p,u,p,mo_integrals_map)

        enddo
      enddo
    enddo
    !$OMP END DO

   !   do p = 1, mo_num
   !     do q = 1, mo_num
   !       tmp_accu(q,p) = 0d0
   !       do u = 1, mo_num
   !         do t = 1, mo_num

   !           tmp_accu(q,p)  = tmp_accu(q,p) &
   !           - tmp_bi_int_3(t,u,p) * tmp_2rdm_3(t,u,q)

   !         enddo
   !       enddo
   !      enddo
   !   enddo

   call dgemm('T','N',mo_num,mo_num,mo_num*mo_num,1d0,tmp_bi_int_3_shared,&
    mo_num*mo_num,tmp_2rdm_3_shared,mo_num*mo_num,0d0,tmp_accu,mo_num)

    !$OMP MASTER
    do p = 1, mo_num
      do q = 1, mo_num

        tmp_h_pqqp(q,p) = tmp_h_pqqp(q,p) - tmp_accu(q,p)

      enddo
    enddo
    !$OMP END MASTER
   
  !!! second part of the term 
  
!!    do q = 1, mo_num
!!      do t = 1, mo_num
!!        do u = 1, mo_num
!!
!!          tmp_2rdm_3(u,t,q) = two_e_dm_mo(q,u,t,q,1)
!!
!!        enddo
!!      enddo 
!!    enddo
!! 
!!    do p = 1, mo_num 
!!       
!!        do u = 1, mo_num
!!          do t = 1, mo_num
!!  
!!             tmp_bi_int_3(t,u,p) = 2d0 * get_two_e_integral(t,u,p,p,mo_integrals_map)
!!
!!          enddo
!!        enddo
!!      enddo
!!   
!!     do p = 1, mo_num    
!!       do q = 1, mo_num
!!        tmp_accu(q,p) = 0d0
!!        do u = 1, mo_num
!!          do t = 1, mo_num
!!  
!!            tmp_accu(q,p) = tmp_accu(q,p)& 
!!            - tmp_bi_int_3(t,u,p) * tmp_2rdm_3(t,u,q)
!!
!!          enddo
!!        enddo
!!      enddo
!!    enddo
!!  
!!    do p = 1, mo_num
!!      do q = 1, mo_num
!!  
!!        tmp_h_pqqp(q,p) = tmp_h_pqqp(q,p) + tmp_accu(q,p)
!!  
!!      enddo
!!    enddo

    !!!!!!!!!!!!!!!!!
  
    !$OMP DO
    do p = 1, mo_num
      do u = 1, mo_num
        do t = 1, mo_num
 
          tmp_bi_int_3_shared(t,u,p) = 2d0* get_two_e_integral(t,u,p,p,mo_integrals_map)
 
        enddo
      enddo
    enddo
    !$OMP END DO

    !$OMP DO
    do q = 1, mo_num
      do t = 1, mo_num
        do u = 1, mo_num

          tmp_2rdm_3_shared(u,t,q) = two_e_dm_mo(q,u,t,q,1)

        enddo
      enddo
    enddo
    !$OMP END DO
 
    call dgemm('T','N',mo_num,mo_num,mo_num*mo_num,1d0,tmp_2rdm_3_shared,&
      mo_num*mo_num,tmp_bi_int_3_shared,mo_num*mo_num,0d0,tmp_accu,mo_num)

    !$OMP MASTER
    do p = 1, mo_num
      do q = 1, mo_num

        tmp_h_pqqp(p,q) = tmp_h_pqqp(p,q) - tmp_accu(p,q)

       enddo
    enddo
    !$OMP END MASTER
 
    !$OMP MASTER
    CALL wall_TIME(t5)
    t6= t5-t4
    print*,'l3 2',t6
    !$OMP END MASTER  

    !$OMP MASTER
    CALL wall_TIME(t2)
    t2 = t2 - t1
    print*, 'Time to compute the hessian :', t2
    !$OMP END MASTER

    !$OMP END PARALLEL
    call omp_set_max_active_levels(4)

  deallocate(tmp_2rdm_3,tmp_bi_int_3,tmp_accu,tmp_accu_1)

  !===========
  ! 2D matrix
  !===========

  ! Convert the hessian mo_num * mo_num * mo_num * mo_num matrix in a
  ! 2D n * n matrix (n = mo_num*(mo_num-1)/2)
  ! H(pq,rs) : p<q and r<s
  ! Hessian(p,q,r,s) = P_pq P_rs [ ...]
  ! => Hessian(p,q,r,s) = (p,q,r,s) - (q,p,r,s) - (p,q,s,r) + (q,p,s,r)

  ! Permutations  
  
  do p = 1, mo_num
    hessian(p,p,p,p) = hessian(p,p,p,p) + tmp_h_pppp(p)
  enddo

  do q = 1, mo_num
    do p = 1, mo_num
      hessian(p,q,p,q) = hessian(p,q,p,q) + tmp_h_pqpq(p,q)
    enddo
  enddo
  
  do q = 1, mo_num
    do p = 1, mo_num
      hessian(p,q,q,p) = hessian(p,q,q,p) + tmp_h_pqqp(p,q)
    enddo
  enddo

 deallocate(tmp_h_pppp,tmp_h_pqpq,tmp_h_pqqp)

!  do q = 1, mo_num
!    do p = 1, mo_num
!      h_tmpr(p,q,p,q) = (hessian(p,q,p,q) - hessian(q,p,p,q) - hessian(p,q,q,p) + hessian(q,p,q,p))
!    enddo
!  enddo
 
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
!  if (debug) then
!    pq=0
!    rs=0
!    do p=1,mo_num
!      do q = 1,mo_num
!      pq=pq+1
!      rs=0
!        do r = 1, mo_num
!          do s = 1, mo_num
!          rs = rs+1 
!           H_test(pq,rs) = h_tmpr(p,q,r,s)
!          enddo
!        enddo
!      enddo
!    enddo
!  
!    print*,'mo_num**2 by mo_num**2 hessian matrix'
!    do pq=1,mo_num**2
!      write(*,'(100(F10.5))') H_test(pq,:)
!    enddo
!  endif

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
    print*,'2D diag Hessian matrix'
    do pq = 1, n
      write(*,'(100(F10.5))') H(pq,:)
    enddo 
  endif

  !==============
  ! Deallocation
  !==============

  deallocate(hessian)!,h_tmpr,H_test)

  if (debug) then
    print*,'Leaves diag_hess'
  endif

end subroutine
