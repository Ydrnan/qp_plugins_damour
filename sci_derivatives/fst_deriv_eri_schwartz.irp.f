
! ---

double precision function fst_deriv_eri_schwartz(i_axis, i_nucl, i, j, k, l)

  BEGIN_DOC
  ! 
  ! \partial_a ( Gi(1) Gj(1) | 1 / r_12 | Gk(2) Gl(2) )
  ! where a = Xa, Ya or Za 
  !
  END_DOC

  implicit none
  include 'utils/constants.include.F'

  integer, intent(in)           :: i_axis ! index of axis; 1 : \partial_x, 2 : \partial_y, 3 : \partial_z
  integer, intent(in)           :: i_nucl ! index of nuclueus
  integer, intent(in)           :: i, j, k, l
                                
  integer                       :: p, q, r, s, ii
  integer                       :: num_i, num_j, num_k, num_l, dim1, I_power(3), J_power(3), K_power(3), L_power(3)
  integer                       :: iorder_p(3) , iorder_q(3)
  double precision              :: I_center(3), J_center(3), K_center(3), L_center(3)
  double precision              :: P_new(0:max_dim,3) , P_center(3) , fact_p , pp , p_inv
  double precision              :: Q_new(0:max_dim,3) , Q_center(3) , fact_q , qq , q_inv
  double precision              :: fact_int, thr, tmp_ij, tmp_kl, integral
  double precision              :: coef1, coef2, coef3, coef4
  double precision              :: expo1, expo2, expo3, expo4
  double precision, allocatable :: schwartz_ij(:,:), schwartz_kl(:,:)

  double precision              :: general_primitive_integral


  fst_deriv_eri_schwartz = 0.d0

  dim1 = n_pt_max_integrals

  num_i = ao_nucl(i)
  num_j = ao_nucl(j)
  num_k = ao_nucl(k)
  num_l = ao_nucl(l)

  do ii = 1, 3
    I_power(ii)  = ao_power(i,ii)
    J_power(ii)  = ao_power(j,ii)
    K_power(ii)  = ao_power(k,ii)
    L_power(ii)  = ao_power(l,ii)

    I_center(ii) = nucl_coord(num_i,ii)
    J_center(ii) = nucl_coord(num_j,ii)
    K_center(ii) = nucl_coord(num_k,ii)
    L_center(ii) = nucl_coord(num_l,ii)
  enddo

  if((num_i .eq. num_j) .and. (num_k .eq. num_l) .and. (num_i .eq. num_k)) then
    return
  endif

  thr = ao_integrals_threshold * ao_integrals_threshold

  ! ---

  if( (i_nucl .eq. num_i) .or. (i_nucl .eq. num_j) ) then  

    allocate( schwartz_kl(0:ao_prim_num(l),0:ao_prim_num(k)) )
    
    schwartz_kl(0,0) = 0.d0
    do r = 1, ao_prim_num(k)
      coef1 = ao_coef_normalized_ordered_transp(r,k) * ao_coef_normalized_ordered_transp(r,k)
      expo1 = ao_expo_ordered_transp(r,k)
    
      schwartz_kl(0,r) = 0.d0
      do s = 1, ao_prim_num(l)
        coef2 = coef1 * ao_coef_normalized_ordered_transp(s,l) * ao_coef_normalized_ordered_transp(s,l)
        expo2 = ao_expo_ordered_transp(s,l)
    
        call give_explicit_poly_and_gaussian( Q_new, Q_center, qq, fact_q, iorder_q                    &
                                            , expo1, expo2, K_power, L_power, K_center, L_center, dim1 )
        q_inv = 1.d0 / qq
    
        schwartz_kl(s,r) = coef2 * general_primitive_integral( dim1, Q_new, Q_center, fact_q, qq, q_inv, iorder_q &
                                                                   , Q_new, Q_center, fact_q, qq, q_inv, iorder_q )
    
        schwartz_kl(0,r) = max(schwartz_kl(0,r), schwartz_kl(s,r))
      enddo
    
      schwartz_kl(0,0) = max(schwartz_kl(0,r), schwartz_kl(0,0))
    enddo

  endif

  ! ---

  if( (i_nucl .eq. num_k) .or. (i_nucl .eq. num_l) ) then  

    allocate( schwartz_ij(0:ao_prim_num(j),0:ao_prim_num(i)) )
    
    schwartz_ij(0,0) = 0.d0
    do p = 1, ao_prim_num(i)
      coef1 = ao_coef_normalized_ordered_transp(p,i) * ao_coef_normalized_ordered_transp(p,i)
      expo1 = ao_expo_ordered_transp(p,i)

      schwartz_ij(0,p) = 0.d0
      do q = 1, ao_prim_num(j)
        coef2 = coef1 * ao_coef_normalized_ordered_transp(q,j) * ao_coef_normalized_ordered_transp(q,j)
        expo2 = ao_expo_ordered_transp(q,j)
    
        call give_explicit_poly_and_gaussian( P_new, P_center, pp, fact_p, iorder_p                    &
                                            , expo1, expo2, I_power, J_power, I_center, J_center, dim1 )
        p_inv = 1.d0 / pp
    
        schwartz_ij(q,p) = general_primitive_integral( dim1, P_new, P_center, fact_p, pp, p_inv, iorder_p &
                                                           , P_new, P_center, fact_p, pp, p_inv, iorder_p )
    
        schwartz_ij(0,p) = max(schwartz_ij(0,p), schwartz_ij(q,p))
      enddo
    
      schwartz_ij(0,0) = max(schwartz_ij(0,p), schwartz_ij(0,0))
    enddo

  endif

  ! ---

  ! -----------------------------------------------------------------------------------------------

  if(i_nucl .eq. num_i) then  

    ! ---

    fact_int        = -1.d0 * dble(I_power(i_axis))
    I_power(i_axis) = I_power(i_axis) - 1
    if(I_power(i_axis) .gt. -1) then

      do p = 1, ao_prim_num(i)
        coef1 = ao_coef_normalized_ordered_transp(p,i)
        expo1 = ao_expo_ordered_transp(p,i)

        do q = 1, ao_prim_num(j)
          coef2 = coef1 * ao_coef_normalized_ordered_transp(q,j)
          expo2 = ao_expo_ordered_transp(q,j)

          call give_explicit_poly_and_gaussian( P_new, P_center, pp, fact_p, iorder_p                    &
                                              , expo1, expo2, I_power, J_power, I_center, J_center, dim1 )
          p_inv = 1.d0 / pp

          tmp_ij = coef2 * coef2 * general_primitive_integral( dim1, P_new, P_center, fact_p, pp, p_inv, iorder_p &
                                                                   , P_new, P_center, fact_p, pp, p_inv, iorder_p )

          if(schwartz_kl(0,0) * tmp_ij < thr) cycle

          do r = 1, ao_prim_num(k)
            if(schwartz_kl(0,r) * tmp_ij < thr) cycle
            coef3 = coef2 * ao_coef_normalized_ordered_transp(r,k)
            expo3 = ao_expo_ordered_transp(r,k)

            do s = 1, ao_prim_num(l)
              if(schwartz_kl(s,r) * tmp_ij < thr) cycle
              coef4 = coef3 * ao_coef_normalized_ordered_transp(s,l)
              expo4 = ao_expo_ordered_transp(s,l)

              call give_explicit_poly_and_gaussian( Q_new, Q_center, qq, fact_q, iorder_q                    &
                                                  , expo3, expo4, K_power, L_power, K_center, L_center, dim1 )
              q_inv = 1.d0 / qq

              integral = general_primitive_integral( dim1, P_new, P_center, fact_p, pp, p_inv, iorder_p &
                                                         , Q_new, Q_center, fact_q, qq, q_inv, iorder_q )

              fst_deriv_eri_schwartz = fst_deriv_eri_schwartz + coef4 * fact_int * integral
            enddo  ! s
          enddo  ! r
        enddo  ! q
      enddo  ! p

    endif 

    I_power(i_axis) = I_power(i_axis) + 1

    ! ---

    I_power(i_axis) = I_power(i_axis) + 1

    do p = 1, ao_prim_num(i)
      coef1 = ao_coef_normalized_ordered_transp(p,i)
      expo1 = ao_expo_ordered_transp(p,i)

      fact_int = 2.d0 * expo1

      do q = 1, ao_prim_num(j)
        coef2 = coef1 * ao_coef_normalized_ordered_transp(q,j)
        expo2 = ao_expo_ordered_transp(q,j)

        call give_explicit_poly_and_gaussian( P_new, P_center, pp, fact_p, iorder_p                    &
                                            , expo1, expo2, I_power, J_power, I_center, J_center, dim1 )
        p_inv = 1.d0 / pp

        tmp_ij = coef2 * coef2 * general_primitive_integral( dim1, P_new, P_center, fact_p, pp, p_inv, iorder_p &
                                                                 , P_new, P_center, fact_p, pp, p_inv, iorder_p )

        if(schwartz_kl(0,0) * tmp_ij < thr) cycle

        do r = 1, ao_prim_num(k)
          if(schwartz_kl(0,r) * tmp_ij < thr) cycle
          coef3 = coef2 * ao_coef_normalized_ordered_transp(r,k)
          expo3 = ao_expo_ordered_transp(r,k)

          do s = 1, ao_prim_num(l)
            if(schwartz_kl(s,r) * tmp_ij < thr) cycle
            coef4 = coef3 * ao_coef_normalized_ordered_transp(s,l)
            expo4 = ao_expo_ordered_transp(s,l)

            call give_explicit_poly_and_gaussian( Q_new, Q_center, qq, fact_q, iorder_q                    &
                                                , expo3, expo4, K_power, L_power, K_center, L_center, dim1 )
            q_inv = 1.d0 / qq

            integral = general_primitive_integral( dim1, P_new, P_center, fact_p, pp, p_inv, iorder_p &
                                                       , Q_new, Q_center, fact_q, qq, q_inv, iorder_q )

            fst_deriv_eri_schwartz = fst_deriv_eri_schwartz + coef4 * fact_int * integral
          enddo  ! s
        enddo  ! r
      enddo  ! q
    enddo  ! p

    I_power(i_axis) = I_power(i_axis) - 1

    ! ---

  endif

  ! -----------------------------------------------------------------------------------------------

  ! ---

  ! -----------------------------------------------------------------------------------------------

  if(i_nucl .eq. num_j) then  

    ! ---

    fact_int        = -1.d0 * dble(J_power(i_axis))
    J_power(i_axis) = J_power(i_axis) - 1
    if(J_power(i_axis) .gt. -1) then

      do p = 1, ao_prim_num(i)
        coef1 = ao_coef_normalized_ordered_transp(p,i)
        expo1 = ao_expo_ordered_transp(p,i)

        do q = 1, ao_prim_num(j)
          coef2 = coef1 * ao_coef_normalized_ordered_transp(q,j)
          expo2 = ao_expo_ordered_transp(q,j)

          call give_explicit_poly_and_gaussian( P_new, P_center, pp, fact_p, iorder_p                    &
                                              , expo1, expo2, I_power, J_power, I_center, J_center, dim1 )
          p_inv = 1.d0 / pp

          tmp_ij = coef2 * coef2 * general_primitive_integral( dim1, P_new, P_center, fact_p, pp, p_inv, iorder_p &
                                                                   , P_new, P_center, fact_p, pp, p_inv, iorder_p )

          if(schwartz_kl(0,0) * tmp_ij < thr) cycle

          do r = 1, ao_prim_num(k)
            if(schwartz_kl(0,r) * tmp_ij < thr) cycle
            coef3 = coef2 * ao_coef_normalized_ordered_transp(r,k)
            expo3 = ao_expo_ordered_transp(r,k)

            do s = 1, ao_prim_num(l)
              if(schwartz_kl(s,r) * tmp_ij < thr) cycle
              coef4 = coef3 * ao_coef_normalized_ordered_transp(s,l)
              expo4 = ao_expo_ordered_transp(s,l)

              call give_explicit_poly_and_gaussian( Q_new, Q_center, qq, fact_q, iorder_q                    &
                                                  , expo3, expo4, K_power, L_power, K_center, L_center, dim1 )
              q_inv = 1.d0 / qq

              integral = general_primitive_integral( dim1, P_new, P_center, fact_p, pp, p_inv, iorder_p &
                                                         , Q_new, Q_center, fact_q, qq, q_inv, iorder_q )

              fst_deriv_eri_schwartz = fst_deriv_eri_schwartz + coef4 * fact_int * integral
            enddo  ! s
          enddo  ! r
        enddo  ! q
      enddo  ! p

    endif 

    J_power(i_axis) = J_power(i_axis) + 1

    ! ---

    J_power(i_axis) = J_power(i_axis) + 1

    do q = 1, ao_prim_num(j)
      coef2 = ao_coef_normalized_ordered_transp(q,j)
      expo2 = ao_expo_ordered_transp(q,j)

      fact_int = 2.d0 * expo2

      do p = 1, ao_prim_num(i)
        coef1 = coef2 * ao_coef_normalized_ordered_transp(p,i)
        expo1 = ao_expo_ordered_transp(p,i)
     
        call give_explicit_poly_and_gaussian( P_new, P_center, pp, fact_p, iorder_p                    &
                                            , expo1, expo2, I_power, J_power, I_center, J_center, dim1 )
        p_inv = 1.d0 / pp

        tmp_ij = coef2 * coef2 * general_primitive_integral( dim1, P_new, P_center, fact_p, pp, p_inv, iorder_p &
                                                                 , P_new, P_center, fact_p, pp, p_inv, iorder_p )

        if(schwartz_kl(0,0) * tmp_ij < thr) cycle

        do r = 1, ao_prim_num(k)
          if(schwartz_kl(0,r) * tmp_ij < thr) cycle
          coef3 = coef1 * ao_coef_normalized_ordered_transp(r,k)
          expo3 = ao_expo_ordered_transp(r,k)

          do s = 1, ao_prim_num(l)
            if(schwartz_kl(s,r) * tmp_ij < thr) cycle
            coef4 = coef3 * ao_coef_normalized_ordered_transp(s,l)
            expo4 = ao_expo_ordered_transp(s,l)

            call give_explicit_poly_and_gaussian( Q_new, Q_center, qq, fact_q, iorder_q                    &
                                                , expo3, expo4, K_power, L_power, K_center, L_center, dim1 )
            q_inv = 1.d0 / qq

            integral = general_primitive_integral( dim1, P_new, P_center, fact_p, pp, p_inv, iorder_p &
                                                       , Q_new, Q_center, fact_q, qq, q_inv, iorder_q )

            fst_deriv_eri_schwartz = fst_deriv_eri_schwartz + coef4 * fact_int * integral
          enddo  ! s
        enddo  ! r
      enddo  ! q
    enddo  ! p

    J_power(i_axis) = J_power(i_axis) - 1

    ! ---

  endif

  ! -----------------------------------------------------------------------------------------------

  ! ---

  ! -----------------------------------------------------------------------------------------------

  if(i_nucl .eq. num_k) then  

    ! ---

    fact_int        = -1.d0 * dble(K_power(i_axis))
    K_power(i_axis) = K_power(i_axis) - 1
    if(K_power(i_axis) .gt. -1) then

      do r = 1, ao_prim_num(k)
        coef3 = ao_coef_normalized_ordered_transp(r,k)
        expo3 = ao_expo_ordered_transp(r,k)

        do s = 1, ao_prim_num(l)
          coef4 = coef3 * ao_coef_normalized_ordered_transp(s,l)
          expo4 = ao_expo_ordered_transp(s,l)

          call give_explicit_poly_and_gaussian( Q_new, Q_center, qq, fact_q, iorder_q                    &
                                              , expo3, expo4, K_power, L_power, K_center, L_center, dim1 )
          q_inv = 1.d0 / qq

          tmp_kl = coef4 * coef4 * general_primitive_integral( dim1, Q_new, Q_center, fact_q, qq, q_inv, iorder_q &
                                                                   , Q_new, Q_center, fact_q, qq, q_inv, iorder_q )

          if(schwartz_ij(0,0) * tmp_kl < thr) cycle

          do p = 1, ao_prim_num(i)
            if(schwartz_ij(0,p) * tmp_kl < thr) cycle
            coef1 = coef4 * ao_coef_normalized_ordered_transp(p,i)
            expo1 = ao_expo_ordered_transp(p,i)

            do q = 1, ao_prim_num(j)
              if(schwartz_ij(q,p) * tmp_kl < thr) cycle
              coef2 = coef1 * ao_coef_normalized_ordered_transp(q,j)
              expo2 = ao_expo_ordered_transp(q,j)

              call give_explicit_poly_and_gaussian( P_new, P_center, pp, fact_p, iorder_p                    &
                                                  , expo1, expo2, I_power, J_power, I_center, J_center, dim1 )
              p_inv = 1.d0 / pp

              integral = general_primitive_integral( dim1, P_new, P_center, fact_p, pp, p_inv, iorder_p &
                                                         , Q_new, Q_center, fact_q, qq, q_inv, iorder_q )

              fst_deriv_eri_schwartz = fst_deriv_eri_schwartz + coef2 * fact_int * integral
            enddo  ! s
          enddo  ! r
        enddo  ! q
      enddo  ! p

    endif 

    K_power(i_axis) = K_power(i_axis) + 1

    ! ---

    K_power(i_axis) = K_power(i_axis) + 1

    do r = 1, ao_prim_num(k)
      coef3 = ao_coef_normalized_ordered_transp(r,k)
      expo3 = ao_expo_ordered_transp(r,k)

      fact_int = 2.d0 * expo3

      do s = 1, ao_prim_num(l)
        coef4 = coef3 * ao_coef_normalized_ordered_transp(s,l)
        expo4 = ao_expo_ordered_transp(s,l)

        call give_explicit_poly_and_gaussian( Q_new, Q_center, qq, fact_q, iorder_q                    &
                                            , expo3, expo4, K_power, L_power, K_center, L_center, dim1 )
        q_inv = 1.d0 / qq

        tmp_kl = coef4 * coef4 * general_primitive_integral( dim1, Q_new, Q_center, fact_q, qq, q_inv, iorder_q &
                                                                 , Q_new, Q_center, fact_q, qq, q_inv, iorder_q )

        if(schwartz_ij(0,0) * tmp_kl < thr) cycle

        do p = 1, ao_prim_num(i)
          if(schwartz_ij(0,p) * tmp_kl < thr) cycle
          coef1 = coef4 * ao_coef_normalized_ordered_transp(p,i)
          expo1 = ao_expo_ordered_transp(p,i)

          do q = 1, ao_prim_num(j)
            if(schwartz_ij(q,p) * tmp_kl < thr) cycle
            coef2 = coef1 * ao_coef_normalized_ordered_transp(q,j)
            expo2 = ao_expo_ordered_transp(q,j)

            call give_explicit_poly_and_gaussian( P_new, P_center, pp, fact_p, iorder_p                    &
                                                , expo1, expo2, I_power, J_power, I_center, J_center, dim1 )
            p_inv = 1.d0 / pp

            integral = general_primitive_integral( dim1, P_new, P_center, fact_p, pp, p_inv, iorder_p &
                                                       , Q_new, Q_center, fact_q, qq, q_inv, iorder_q )

            fst_deriv_eri_schwartz = fst_deriv_eri_schwartz + coef2 * fact_int * integral
          enddo  ! s
        enddo  ! r
      enddo  ! q
    enddo  ! p

    K_power(i_axis) = K_power(i_axis) - 1

    ! ---

  endif

  ! -----------------------------------------------------------------------------------------------

  ! ---

  ! -----------------------------------------------------------------------------------------------

  if(i_nucl .eq. num_l) then  

    ! ---

    fact_int        = -1.d0 * dble(L_power(i_axis))
    L_power(i_axis) = L_power(i_axis) - 1
    if(L_power(i_axis) .gt. -1) then

      do r = 1, ao_prim_num(k)
        coef3 = ao_coef_normalized_ordered_transp(r,k)
        expo3 = ao_expo_ordered_transp(r,k)

        do s = 1, ao_prim_num(l)
          coef4 = coef3 * ao_coef_normalized_ordered_transp(s,l)
          expo4 = ao_expo_ordered_transp(s,l)

          call give_explicit_poly_and_gaussian( Q_new, Q_center, qq, fact_q, iorder_q                    &
                                              , expo3, expo4, K_power, L_power, K_center, L_center, dim1 )
          q_inv = 1.d0 / qq

          tmp_kl = coef4 * coef4 * general_primitive_integral( dim1, Q_new, Q_center, fact_q, qq, q_inv, iorder_q &
                                                                   , Q_new, Q_center, fact_q, qq, q_inv, iorder_q )

          if(schwartz_ij(0,0) * tmp_kl < thr) cycle

          do p = 1, ao_prim_num(i)
            if(schwartz_ij(0,p) * tmp_kl < thr) cycle
            coef1 = coef4 * ao_coef_normalized_ordered_transp(p,i)
            expo1 = ao_expo_ordered_transp(p,i)

            do q = 1, ao_prim_num(j)
              if(schwartz_ij(q,p) * tmp_kl < thr) cycle
              coef2 = coef1 * ao_coef_normalized_ordered_transp(q,j)
              expo2 = ao_expo_ordered_transp(q,j)

              call give_explicit_poly_and_gaussian( P_new, P_center, pp, fact_p, iorder_p                    &
                                                  , expo1, expo2, I_power, J_power, I_center, J_center, dim1 )
              p_inv = 1.d0 / pp

              integral = general_primitive_integral( dim1, P_new, P_center, fact_p, pp, p_inv, iorder_p &
                                                         , Q_new, Q_center, fact_q, qq, q_inv, iorder_q )

              fst_deriv_eri_schwartz = fst_deriv_eri_schwartz + coef2 * fact_int * integral
            enddo  ! s
          enddo  ! r
        enddo  ! q
      enddo  ! p

    endif 

    L_power(i_axis) = L_power(i_axis) + 1

    ! ---

    L_power(i_axis) = L_power(i_axis) + 1

    do s = 1, ao_prim_num(l)
      coef4 = ao_coef_normalized_ordered_transp(s,l)
      expo4 = ao_expo_ordered_transp(s,l)

      fact_int = 2.d0 * expo4

      do r = 1, ao_prim_num(k)
        coef3 = coef4 * ao_coef_normalized_ordered_transp(r,k)
        expo3 = ao_expo_ordered_transp(r,k)

        call give_explicit_poly_and_gaussian( Q_new, Q_center, qq, fact_q, iorder_q                    &
                                            , expo3, expo4, K_power, L_power, K_center, L_center, dim1 )
        q_inv = 1.d0 / qq

        tmp_kl = coef4 * coef4 * general_primitive_integral( dim1, Q_new, Q_center, fact_q, qq, q_inv, iorder_q &
                                                                 , Q_new, Q_center, fact_q, qq, q_inv, iorder_q )

        if(schwartz_ij(0,0) * tmp_kl < thr) cycle

        do p = 1, ao_prim_num(i)
          if(schwartz_ij(0,p) * tmp_kl < thr) cycle
          coef1 = coef3 * ao_coef_normalized_ordered_transp(p,i)
          expo1 = ao_expo_ordered_transp(p,i)

          do q = 1, ao_prim_num(j)
            if(schwartz_ij(q,p) * tmp_kl < thr) cycle
            coef2 = coef1 * ao_coef_normalized_ordered_transp(q,j)
            expo2 = ao_expo_ordered_transp(q,j)

            call give_explicit_poly_and_gaussian( P_new, P_center, pp, fact_p, iorder_p                    &
                                                , expo1, expo2, I_power, J_power, I_center, J_center, dim1 )
            p_inv = 1.d0 / pp

            integral = general_primitive_integral( dim1, P_new, P_center, fact_p, pp, p_inv, iorder_p &
                                                       , Q_new, Q_center, fact_q, qq, q_inv, iorder_q )

            fst_deriv_eri_schwartz = fst_deriv_eri_schwartz + coef2 * fact_int * integral
          enddo  ! s
        enddo  ! r
      enddo  ! q
    enddo  ! p

    L_power(i_axis) = L_power(i_axis) - 1

    ! ---

  endif

  ! -----------------------------------------------------------------------------------------------

  if( (i_nucl .eq. num_i) .or. (i_nucl .eq. num_j) ) deallocate( schwartz_kl )
  if( (i_nucl .eq. num_k) .or. (i_nucl .eq. num_l) ) deallocate( schwartz_ij )

end function fst_deriv_eri_schwartz

! ---

