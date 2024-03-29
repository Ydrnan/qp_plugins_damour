subroutine hess(n,H,h_tmpr)
  use omp_lib
  include 'constants.h' 

  implicit none

  !==================================================================
  ! Compute the hessian of energy with respects to orbital rotations
  !==================================================================

  !===========
  ! Variables  
  !===========

  !====
  ! in
  !====
  integer, intent(in)           :: n 
  !n         : integer, n = mo_num*(mo_num-1)/2
  
  !=====
  ! out
  !=====
  double precision, intent(out) :: H(n,n),h_tmpr(mo_num,mo_num,mo_num,mo_num)
  ! H        : n by n double precision matrix containing the 2D hessian
 
  !==========
  ! internal
  !==========
  double precision, allocatable :: hessian(:,:,:,:)!, h_tmpr(:,:,:,:)
  double precision, allocatable :: H_test(:,:)
  integer                       :: p,q
  integer                       :: r,s,t,u,v,k
  integer                       :: pq,rs
  double precision              :: t1,t2,t3,t4,t5,t6
  ! hessian  : mo_num 4D double precision matrix containing the hessian before the permutations
  ! h_tmpr   : mo_num 4D double precision matrix containing the hessian after the permutations
  ! H_test   : monum**2 by mo_num**2 double precision matrix to debug the H matrix
  ! p,q,r,s  : integer, indexes of the 4D hessian matrix
  ! t,u,v    : integer, indexes to compute hessian elements
  ! pq,rs    : integer, indexes for the conversion from 4D to 2D hessian matrix
  ! t1,t2,t3 : double precision, t3 = t2 - t1, time to compute the hessian 
  ! t4,t5,t6 : double precision, t6 = t5 - t4, time to compute each element

  double precision, allocatable :: tmp_bi_int_3(:,:,:), tmp_2rdm_3(:,:,:), ind_3(:,:,:)
  double precision, allocatable :: tmp_accu(:,:), tmp_accu_sym(:,:)
  ! tmp_bi_int_3 : mo_num 3D double precision matrix containinig the bi electronic
  !                integrals with 1 fix index
  ! tmp_2rdm_3   : mo_num 3D double precision matrix containinig the 2 body reduce
  !                density matrix with 1 fix index
  ! ind_3        : mo_num 3D double precision matrix, temporary matrix
  ! tmp_accu     : mo_num 2D double precision matrix, temporary matrix
  ! tmp_accu_sym : mo_num 2D double precision matrix, temporary matrix

  ! Function 
  double precision              :: get_two_e_integral
  ! get_two_e_integral :  double precision function, two e integrals 

  ! Provided :
  ! mo_one_e_integrals : mono e- integrals
  ! get_two_e_integral : two e- integrals
  ! one_e_dm_mo : one body density matrix (state average)
  ! two_e_dm_mo : two body density matrix (state average)

  !============
  ! Allocation
  !============

  allocate(hessian(mo_num,mo_num,mo_num,mo_num))!,h_tmpr(mo_num,mo_num,mo_num,mo_num))

  !=============
  ! Calculation
  !=============

  print*,'Use the full hessian'

  if (debug) then
    print*,'Enter in hessien.irp.f'
  endif

  ! From Anderson et. al. (2014) 
  ! The Journal of Chemical Physics 141, 244104 (2014); doi: 10.1063/1.4904384

  ! LaTeX formula :

  !\begin{align*}
  !H_{pq,rs} &= \dfrac{\partial^2 E(x)}{\partial x_{pq}^2} \\
  !&= \mathcal{P}_{pq} \mathcal{P}_{rs} [ \frac{1}{2} \sum_u [\delta_{qr}(h_p^u \gamma_u^s + h_u^s \gamma_p^u) 
  !+ \delta_{ps}(h_r^u \gamma_u^q + h_u^q \gamma_u^r)]
  !-(h_p^s \gamma_r^q + h_r^q \gamma_p^s) \\
  !&+ \frac{1}{2} \sum_{tuv} [\delta_{qr}(v_{pt}^{uv} \Gamma_{uv}^{st} +v_{uv}^{st} \Gamma_{pt}^{uv}) 
  !+ \delta_{ps}(v_{uv}^{qt} \Gamma_{rt}^{uv} + v_{rt}^{uv}\Gamma_{uv}^{qt})] \\
  !&+ \sum_{uv} (v_{pr}^{uv} \Gamma_{uv}^{qs} + v_{uv}^{qs}  \Gamma_{ps}^{uv}) \\
  !&- \sum_{tu} (v_{pu}^{st} \Gamma_{rt}^{qu}+v_{pu}^{tr} \Gamma_{tr}^{qu}+v_{rt}^{qu}\Gamma_{pu}^{st} + v_{tr}^{qu}\Gamma_{pu}^{ts}) 
  !\end{align*}

  call omp_set_max_active_levels(1)

  ! OMP 
  !$OMP PARALLEL                                                     &
      !$OMP PRIVATE(                                                 &
      !$OMP   p,q,r,s, tmp_accu, tmp_accu_sym,                       &
      !$OMP   u,v,t, tmp_bi_int_3, tmp_2rdm_3, ind_3)                       &
      !$OMP SHARED(hessian,h_tmpr,H, mo_num,n, & 
      !$OMP mo_one_e_integrals, one_e_dm_mo, &
      !$OMP two_e_dm_mo,mo_integrals_map, &
      !$OMP t1,t2,t3,t4,t5,t6)&
      !$OMP DEFAULT(NONE)
 
  ! Allocation of private arrays 
  allocate(tmp_bi_int_3(mo_num,mo_num,mo_num))
  allocate(tmp_2rdm_3(mo_num,mo_num,mo_num), ind_3(mo_num,mo_num,mo_num))
  allocate(tmp_accu(mo_num,mo_num), tmp_accu_sym(mo_num,mo_num))

  !$OMP MASTER
  do q = 1, mo_num
    do p = 1, mo_num
      tmp_accu(p,q) = 0d0
    enddo
  enddo
  !$OMP END MASTER

  !$OMP MASTER
  do q = 1, mo_num
    do p = 1, mo_num
      tmp_accu_sym(p,q) = 0d0
    enddo
  enddo
  !$OMP END MASTER

  !$OMP DO
  do s=1,mo_num
    do r=1,mo_num
      do q=1,mo_num
        do p=1,mo_num
          hessian(p,q,r,s) = 0d0
        enddo
      enddo
    enddo
  enddo
  !$OMP ENDDO
 
  !$OMP MASTER
  CALL wall_TIME(t1)
  !$OMP END MASTER

  !========================
  ! First line, first term
  !========================

  !do p = 1, mo_num
  !  do q = 1, mo_num
  !    do r = 1, mo_num
  !      do s = 1, mo_num

  !          if (q==r) then
  !            do u = 1, mo_num

  !              hessian(p,q,r,s) = hessian(p,q,r,s) + 0.5d0 * (  &
  !                mo_one_e_integrals(u,p) * one_e_dm_mo(u,s) &
  !              + mo_one_e_integrals(s,u) * one_e_dm_mo(p,u))

  !            enddo
  !          endif

  !      enddo
  !    enddo
  !  enddo
  !enddo

  !**************************** 
  ! Opt First line, first term
  !****************************

  !$OMP MASTER    
  CALL wall_TIME(t4)
  !$OMP END MASTER

  call dgemm('T','N', mo_num, mo_num, mo_num, 1d0, mo_one_e_integrals,&
             size(mo_one_e_integrals,1), one_e_dm_mo, size(one_e_dm_mo,1),&
             0d0, tmp_accu, size(tmp_accu,1))

  !$OMP DO
  do s = 1, mo_num
    do p = 1, mo_num

      tmp_accu_sym(p,s) = 0.5d0 * (tmp_accu(p,s) + tmp_accu(s,p))

    enddo
  enddo 
  !$OMP END DO

  !$OMP DO
  do s = 1, mo_num
    do p = 1, mo_num
      do r = 1, mo_num

        hessian(p,r,r,s) = hessian(p,r,r,s) + tmp_accu_sym(p,s)

      enddo
    enddo
  enddo
  !$OMP END DO
 
  !$OMP MASTER
  CALL wall_TIME(t5)
  t6=t5-t4
  print*,'l1 1',t6
  !$OMP END MASTER
  
  !=========================
  ! First line, second term
  !=========================

  !do p = 1, mo_num
  !  do q = 1, mo_num
  !    do r = 1, mo_num
  !      do s = 1, mo_num
 
  !        if (p==s) then
  !          do u = 1, mo_num
 
  !                hessian(p,q,r,s) = hessian(p,q,r,s) + 0.5d0 * ( &
  !                  mo_one_e_integrals(u,r) * (one_e_dm_mo(u,q) &
  !                + mo_one_e_integrals(q,u) * (one_e_dm_mo(r,u))
  !          enddo
  !        endif
 
  !      enddo
  !    enddo
  !  enddo
  !enddo

  !*****************************
  ! Opt First line, second term
  !*****************************

  !$OMP MASTER
  CALL wall_TIME(t4)
  !$OMP END MASTER

  call dgemm('T','N', mo_num, mo_num, mo_num, 1d0, mo_one_e_integrals,&
             size(mo_one_e_integrals,1), one_e_dm_mo, size(one_e_dm_mo,1),&
             0d0, tmp_accu, size(tmp_accu,1))

  !$OMP DO
  do r = 1, mo_num
    do q = 1, mo_num

      tmp_accu_sym(q,r) = 0.5d0 * (tmp_accu(q,r) + tmp_accu(r,q))

    enddo
  enddo
  !OMP END DO

  !$OMP DO
  do r = 1, mo_num
    do q = 1, mo_num
      do s = 1, mo_num

        hessian(s,q,r,s) = hessian(s,q,r,s) + tmp_accu_sym(q,r)

      enddo
    enddo
  enddo
  !OMP END DO

  !$OMP MASTER
  CALL wall_TIME(t5)
  t6=t5-t4
  print*,'l1 2',t6
  !$OMP END MASTER
 
  !========================
  ! First line, third term
  !========================

  !do p = 1, mo_num
  !  do q = 1, mo_num
  !    do r = 1, mo_num
  !      do s = 1, mo_num
 
  !        hessian(p,q,r,s) = hessian(p,q,r,s) &
  !        - mo_one_e_integrals(s,p) * one_e_dm_mo(r,q) &
  !        - mo_one_e_integrals(q,r) * one_e_dm_mo(p,s))
 
  !      enddo
  !    enddo
  !  enddo
  !enddo
 
  !****************************
  ! Opt First line, third term
  !****************************

  !$OMP MASTER
  CALL wall_TIME(t4)
  !$OMP END MASTER

  !$OMP DO
  do s = 1, mo_num
    do r = 1, mo_num
      do q = 1, mo_num
        do p = 1, mo_num

          hessian(p,q,r,s) = hessian(p,q,r,s) &
            - mo_one_e_integrals(s,p) * one_e_dm_mo(r,q)&
            - mo_one_e_integrals(q,r) * one_e_dm_mo(p,s)

        enddo
      enddo
    enddo
  enddo
  !$OMP END DO

  !$OMP MASTER
  CALL wall_TIME(t5)
  t6=t5-t4
  print*,'l1 3',t6
  !$OMP END MASTER
  
  !=========================
  ! Second line, first term
  !=========================

  !do p = 1, mo_num
  !  do q = 1, mo_num
  !    do r = 1, mo_num
  !      do s = 1, mo_num
 
  !         if (q==r) then
  !           do t = 1, mo_num
  !             do u = 1, mo_num
  !               do v = 1, mo_num
 
  !                 hessian(p,q,r,s) = hessian(p,q,r,s) + 0.5d0 * (  &
  !                   get_two_e_integral(u,v,p,t,mo_integrals_map) * two_e_dm_mo(u,v,s,t) &
  !                 + get_two_e_integral(s,t,u,v,mo_integrals_map) * two_e_dm_mo(p,t,u,v))
 
  !               enddo
  !             enddo
  !           enddo
  !         endif
 
  !      enddo
  !    enddo
  !  enddo
  !enddo
 
  !*****************************
  ! Opt Second line, first term
  !*****************************

  !$OMP MASTER 
  call wall_TIME(t4)
  !$OMP END MASTER

  !$OMP DO
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
 
          tmp_2rdm_3(u,v,p) = two_e_dm_mo(u,v,p,t)  

        enddo
      enddo
    enddo

    call dgemm('T','N', mo_num, mo_num, mo_num*mo_num, 1.d0, &
               tmp_bi_int_3, mo_num*mo_num, tmp_2rdm_3, mo_num*mo_num, &
               0.d0, tmp_accu, size(tmp_accu,1))

    do p = 1, mo_num
      do s = 1, mo_num

        tmp_accu_sym(s,p) = 0.5d0 * (tmp_accu(p,s)+tmp_accu(s,p))

      enddo
    enddo

    !$OMP CRITICAL 
    do s = 1, mo_num
      do r = 1, mo_num
        do p = 1, mo_num

          hessian(p,r,r,s) = hessian(p,r,r,s) + tmp_accu_sym(p,s) 

        enddo
      enddo
    enddo
    !$OMP END CRITICAL

  enddo
  !$OMP END DO

  !$OMP MASTER
  call wall_TIME(t5)
  t6=t5-t4
  print*,'l2 1', t6 
  !$OMP END MASTER

  !==========================
  ! Second line, second term
  !==========================

  !do p = 1, mo_num
  !  do q = 1, mo_num
  !    do r = 1, mo_num
  !      do s = 1, mo_num
  !
  !        if (p==s) then
  !          do t = 1, mo_num
  !            do u = 1, mo_num
  !              do v = 1, mo_num
  !
  !                hessian(p,q,r,s) = hessian(p,q,r,s) + 0.5d0 * ( &
  !                  get_two_e_integral(q,t,u,v,mo_integrals_map) * two_e_dm_mo(r,t,u,v) &
  !                + get_two_e_integral(u,v,r,t,mo_integrals_map) * two_e_dm_mo(u,v,q,t))
  !
  !              enddo
  !            enddo
  !          enddo
  !        endif
  !
  !      enddo
  !    enddo
  !  enddo
  !enddo

  !******************************
  ! Opt Second line, second term
  !******************************

  !$OMP MASTER 
  CALL wall_TIME(t4)
  !$OMP END MASTER

  !$OMP DO
  do t = 1, mo_num

    do q = 1, mo_num
      do v = 1, mo_num
        do u = 1, mo_num

          tmp_bi_int_3(u,v,q) = get_two_e_integral(u,v,q,t,mo_integrals_map)

        enddo
      enddo
    enddo

    do r = 1, mo_num
      do v = 1, mo_num
        do u = 1, mo_num

           tmp_2rdm_3(u,v,r) = two_e_dm_mo(u,v,r,t)

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

    !$OMP CRITICAL
    do r = 1, mo_num
      do q = 1, mo_num
        do s = 1, mo_num

          hessian(s,q,r,s) = hessian(s,q,r,s) + tmp_accu_sym(q,r)

        enddo
      enddo
    enddo
    !$OMP END CRITICAL

  enddo
  !$OMP END DO

  !$OMP MASTER
  CALL wall_TIME(t5)
  t6=t5-t4
  print*,'l2 2',t6
  !$OMP END MASTER

  !========================
  ! Third line, first term
  !========================

  !do p = 1, mo_num
  !  do q = 1, mo_num
  !    do r = 1, mo_num
  !      do s = 1, mo_num
 
  !        do u = 1, mo_num
  !          do v = 1, mo_num
 
  !            hessian(p,q,r,s) = hessian(p,q,r,s) &
  !             + get_two_e_integral(u,v,p,r,mo_integrals_map) * two_e_dm_mo(u,v,q,s) &
  !             + get_two_e_integral(q,s,u,v,mo_integrals_map) * two_e_dm_mo(p,r,u,v)
 
  !          enddo
  !        enddo
 
  !      enddo
  !    enddo
  !  enddo
  !enddo

  !****************************
  ! Opt Third line, first term
  !****************************

  !$OMP MASTER 
  call wall_TIME(t4)
  !$OMP END MASTER 

  !--------
  ! part 1
  ! get_two_e_integral(u,v,p,r,mo_integrals_map) * two_e_dm_mo(u,v,q,s)
  !--------

  !$OMP DO
  do v = 1, mo_num

    do u = 1, mo_num 
      do r = 1, mo_num
        do p = 1, mo_num

            tmp_bi_int_3(p,r,u) = get_two_e_integral(p,r,u,v,mo_integrals_map)

        enddo
      enddo
    enddo

    do s = 1, mo_num
      do q = 1, mo_num
        do u = 1, mo_num

          tmp_2rdm_3(u,q,s) = two_e_dm_mo(q,s,u,v)

        enddo
      enddo
    enddo
    
    do s = 1, mo_num
                
      call dgemm('N','N',mo_num*mo_num, mo_num, mo_num, 1d0, tmp_bi_int_3,&
                 size(tmp_bi_int_3,1)*size(tmp_bi_int_3,2), tmp_2rdm_3(1,1,s),&
                 size(tmp_2rdm_3,1), 0d0, ind_3, size(ind_3,1) * size(ind_3,2))

      !$OMP CRITICAL
      do r = 1, mo_num
        do q = 1, mo_num
          do p = 1, mo_num
            hessian(p,q,r,s) = hessian(p,q,r,s) + ind_3(p,r,q)
          enddo
        enddo
      enddo
      !$OMP END CRITICAL

    enddo

  enddo
  !$OMP END DO

  !--------
  ! part 2
  ! get_two_e_integral(q,s,u,v,mo_integrals_map) * two_e_dm_mo(p,r,u,v)
  !--------

  !$OMP DO
  do v = 1, mo_num

    do u = 1, mo_num
      do s = 1, mo_num
        do q = 1, mo_num

            tmp_bi_int_3(q,s,u) = get_two_e_integral(q,s,u,v,mo_integrals_map)

        enddo
      enddo
    enddo

    do r = 1, mo_num
      do p = 1, mo_num
        do u = 1, mo_num

          tmp_2rdm_3(u,p,r) = two_e_dm_mo(p,r,u,v)

        enddo
      enddo
    enddo

    do r = 1, mo_num
      call dgemm('N','N', mo_num*mo_num, mo_num, mo_num, 1d0, tmp_bi_int_3,& 
                 size(tmp_bi_int_3,1)*size(tmp_bi_int_3,2), tmp_2rdm_3(1,1,r),&
                 size(tmp_2rdm_3,1), 0d0, ind_3, size(ind_3,1) * size(ind_3,2))

      !$OMP CRITICAL
      do s = 1, mo_num
        do q = 1, mo_num
          do p = 1, mo_num
            hessian(p,q,r,s) = hessian(p,q,r,s) + ind_3(q,s,p)
          enddo
        enddo
      enddo
      !$OMP END CRITICAL

    enddo

  enddo
  !$OMP END DO
   
  !$OMP MASTER
  call wall_TIME(t5)
  t6 = t5 - t4
  print*,'l3 1', t6
  !$OMP END MASTER 

  !=========================
  ! Third line, second term
  !=========================

  !do p = 1, mo_num
  !  do q = 1, mo_num
  !    do r = 1, mo_num
  !      do s = 1, mo_num
 
  !        do t = 1, mo_num
  !          do u = 1, mo_num
 
  !            hessian(p,q,r,s) = hessian(p,q,r,s) &
  !             - get_two_e_integral(s,t,p,u,mo_integrals_map) * two_e_dm_mo(r,t,q,u) &
  !             - get_two_e_integral(t,s,p,u,mo_integrals_map) * two_e_dm_mo(t,r,q,u) &
  !             - get_two_e_integral(q,u,r,t,mo_integrals_map) * two_e_dm_mo(p,u,s,t) &
  !             - get_two_e_integral(q,u,t,r,mo_integrals_map) * two_e_dm_mo(p,u,t,s)
 
  !          enddo
  !        enddo
 
  !      enddo
  !    enddo
  !  enddo
  !enddo
 
  !*****************************
  ! Opt Third line, second term
  !*****************************

  !--------
  ! Part 1 
  ! - get_two_e_integral(s,t,p,u,mo_integrals_map) * two_e_dm_mo(r,t,q,u)
  !--------

  !$OMP MASTER
  CALL wall_TIME(t4) 
  !$OMP END MASTER

  !$OMP DO
  do q = 1, mo_num

    do r = 1, mo_num
      do t = 1, mo_num
        do u = 1, mo_num

          tmp_2rdm_3(u,t,r) = two_e_dm_mo(q,u,r,t)

        enddo
      enddo
    enddo

    do p = 1, mo_num

      do s = 1, mo_num
        do t = 1, mo_num
          do u = 1, mo_num

            tmp_bi_int_3(u,t,s) = - get_two_e_integral(u,s,t,p,mo_integrals_map)

          enddo
        enddo
      enddo

      call dgemm('T','N', mo_num, mo_num, mo_num*mo_num, 1d0, tmp_bi_int_3,&
                 mo_num*mo_num, tmp_2rdm_3, mo_num*mo_num, 0d0, tmp_accu, mo_num)

      !$OMP CRITICAL
      do s = 1, mo_num
        do r = 1, mo_num

           hessian(p,q,r,s) = hessian(p,q,r,s) + tmp_accu(s,r)

        enddo
      enddo
      !$OMP END CRITICAL

    enddo

  enddo
  !$OMP END DO

  !--------
  ! Part 2
  !- get_two_e_integral(t,s,p,u,mo_integrals_map) * two_e_dm_mo(t,r,q,u)
  !--------
  
  !$OMP DO 
  do q = 1, mo_num

    do r = 1, mo_num
      do t = 1, mo_num
        do u = 1, mo_num

          tmp_2rdm_3(u,t,r) = two_e_dm_mo(q,u,t,r)

        enddo
      enddo
    enddo
    
    do p = 1, mo_num

      do s = 1, mo_num
        do t = 1, mo_num
          do u = 1, mo_num

            tmp_bi_int_3(u,t,s) = - get_two_e_integral(u,t,s,p,mo_integrals_map)

          enddo
        enddo
      enddo   

      call dgemm('T','N', mo_num, mo_num, mo_num*mo_num, 1d0, tmp_bi_int_3,&
                 mo_num*mo_num, tmp_2rdm_3, mo_num*mo_num, 0d0, tmp_accu, mo_num)

      !$OMP CRITICAL
      do s = 1, mo_num
        do r = 1, mo_num

          hessian(p,q,r,s) = hessian(p,q,r,s) + tmp_accu(s,r)

        enddo
      enddo
      !$OMP END CRITICAL

    enddo

  enddo
  !$OMP END DO

  !--------
  ! Part 3
  !- get_two_e_integral(q,u,r,t,mo_integrals_map) * two_e_dm_mo(p,u,s,t) 
  !--------
 
  !$OMP DO 
  do q = 1, mo_num

    do r = 1, mo_num
      do t = 1, mo_num
        do u = 1, mo_num

          tmp_bi_int_3(u,t,r) = - get_two_e_integral(u,q,t,r,mo_integrals_map)

        enddo
      enddo
    enddo

    do p = 1, mo_num

      do s = 1, mo_num
        do t = 1, mo_num
          do u = 1, mo_num

            tmp_2rdm_3(u,t,s) = two_e_dm_mo(p,u,s,t)

          enddo
        enddo
      enddo

      call dgemm('T','N', mo_num, mo_num, mo_num*mo_num, 1d0, tmp_2rdm_3,&
                 mo_num*mo_num, tmp_bi_int_3, mo_num*mo_num, 0d0, tmp_accu, mo_num)

      !$OMP CRITICAL
      do s = 1, mo_num
        do r = 1, mo_num

          hessian(p,q,r,s) = hessian(p,q,r,s) + tmp_accu(s,r)

        enddo
      enddo
      !$OMP END CRITICAL

    enddo

  enddo
  !$OMP END DO

  !--------
  ! Part 4
  ! - get_two_e_integral(q,u,t,r,mo_integrals_map) * two_e_dm_mo(p,u,t,s)
  !--------
  
  !$OMP DO
  do q = 1, mo_num

    do r = 1, mo_num
      do t = 1, mo_num
        do u = 1, mo_num

          tmp_bi_int_3(u,t,r) = - get_two_e_integral(u,t,r,q,mo_integrals_map)

        enddo
      enddo
    enddo

    do p = 1, mo_num

      do s = 1, mo_num
        do t = 1, mo_num
          do u = 1, mo_num

            tmp_2rdm_3(u,t,s) = two_e_dm_mo(p,u,t,s)

          enddo
        enddo
      enddo

      call dgemm('T','N', mo_num, mo_num, mo_num*mo_num, 1d0, tmp_2rdm_3,&
                 mo_num*mo_num, tmp_bi_int_3, mo_num*mo_num, 0d0, tmp_accu, mo_num)

      !$OMP CRITICAL
      do s = 1, mo_num
        do r = 1, mo_num

          hessian(p,q,r,s) = hessian(p,q,r,s) + tmp_accu(s,r)

        enddo
      enddo
      !$OMP END CRITICAL

    enddo

  enddo
  !$OMP END DO  

  !$OMP MASTER
  call wall_TIME(t5)
  t6 = t5-t4
  print*,'l3 2',t6
  !$OMP END MASTER

  !$OMP MASTER
  CALL wall_TIME(t2)
  t3 = t2 -t1
  print*,'Time to compute the hessian : ', t3
  !$OMP END MASTER

  !================================ 
  ! Deallocation of private arrays
  !================================

  deallocate(tmp_bi_int_3, tmp_2rdm_3, tmp_accu, tmp_accu_sym, ind_3)

  !==============
  ! Permutations 
  !==============

  ! Hessian(p,q,r,s) = P_pq P_rs [ ...]
  ! => Hessian(p,q,r,s) = (p,q,r,s) - (q,p,r,s) - (p,q,s,r) + (q,p,s,r)

  !$OMP MASTER
  CALL wall_TIME(t4)
  !$OMP END MASTER

  !$OMP DO
  do s = 1, mo_num
    do r = 1, mo_num
      do q = 1, mo_num
        do p = 1, mo_num

          h_tmpr(p,q,r,s) = (hessian(p,q,r,s) - hessian(q,p,r,s) - hessian(p,q,s,r) + hessian(q,p,s,r))

        enddo
      enddo
    enddo
  enddo
  !$OMP END DO

  !$OMP MASTER
  call wall_TIME(t5)
  t6 = t5-t4
  print*,'Time for permutations :',t6
  !$OMP END MASTER

  !===================
  ! 2D hessian matrix
  !===================

  ! Convert the hessian mo_num * mo_num * mo_num * mo_num matrix in a
  ! 2D n * n matrix (n = mo_num*(mo_num-1)/2)
  ! H(pq,rs) : p<q and r<s

  ! 4D mo_num matrix to 2D n matrix

  !$OMP MASTER
  CALL wall_TIME(t4)
  !$OMP END MASTER

  !$OMP DO
  do rs = 1, n
    call vec_to_mat_index(rs,r,s)
    do pq = 1, n
      call vec_to_mat_index(pq,p,q)
      H(pq,rs) = h_tmpr(p,q,r,s)   
    enddo
  enddo
  !$OMP END DO 

  !$OMP MASTER
  call wall_TIME(t5)
  t6 = t5-t4
  print*,'4D -> 2D :',t6
  !$OMP END MASTER

  !$OMP END PARALLEL
  call omp_set_max_active_levels(4)

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

  deallocate(hessian)!,h_tmpr)

  if (debug) then
    print*,'Leave hess'
  endif

end subroutine
