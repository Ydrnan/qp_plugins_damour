
! ---

double precision function fst_deriv_eri(i_axis, i_nucl, i, j, k, l)

  BEGIN_DOC
  ! 
  ! \partial_a ( Gi(1) Gj(1) | 1 / r_12 | Gk(2) Gl(2) )
  ! where a = Xa, Ya or Za 
  !
  END_DOC

  implicit none
  include 'utils/constants.include.F'

  integer, intent(in) :: i_axis ! index of axis; 1 : \partial_x, 2 : \partial_y, 3 : \partial_z
  integer, intent(in) :: i_nucl ! index of nuclueus
  integer, intent(in) :: i, j, k, l

  integer             :: p, q, r, s, ii
  integer             :: num_i, num_j, num_k, num_l, dim1, I_power(3), J_power(3), K_power(3), L_power(3)
  integer             :: iorder_p(3) , iorder_q(3)
  integer             :: iorder_pm(3), iorder_qm(3)
  integer             :: iorder_pp(3), iorder_qp(3)
  double precision    :: I_center(3), J_center(3), K_center(3), L_center(3)
  double precision    :: P_new(0:max_dim,3) , P_center(3) , fact_p , pp , p_inv
  double precision    :: Q_new(0:max_dim,3) , Q_center(3) , fact_q , qq , q_inv
  double precision    :: Pm_new(0:max_dim,3), Pm_center(3), fact_pm, ppm, pm_inv
  double precision    :: Qm_new(0:max_dim,3), Qm_center(3), fact_qm, qqm, qm_inv
  double precision    :: Pp_new(0:max_dim,3), Pp_center(3), fact_pp, ppp, pp_inv
  double precision    :: Qp_new(0:max_dim,3), Qp_center(3), fact_qp, qqp, qp_inv
  double precision    :: fact_intm, fact_intp, integral_m, integral_p
  double precision    :: coef1, coef2, coef3, coef4
  double precision    :: expo1, expo2, expo3, expo4

  double precision    :: general_primitive_integral
  double precision    :: fst_deriv_eri_schwartz

  if(ao_prim_num(i) * ao_prim_num(j) * ao_prim_num(k) * ao_prim_num(l) > 1024) then
    fst_deriv_eri = fst_deriv_eri_schwartz(i_axis, i_nucl, i, j, k, l)
    return
  endif


  fst_deriv_eri = 0.d0

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

  ! -----------------------------------------------------------------------------------------------

  if(i_nucl .eq. num_i) then  

    do p = 1, ao_prim_num(i)
      coef1 = ao_coef_normalized_ordered_transp(p,i)
      expo1 = ao_expo_ordered_transp(p,i)

      do q = 1, ao_prim_num(j)
        coef2 = coef1 * ao_coef_normalized_ordered_transp(q,j)
        expo2 = ao_expo_ordered_transp(q,j)

        ! ---

        fact_intm       = -1.d0 * dble(I_power(i_axis))
        I_power(i_axis) = I_power(i_axis) - 1
        if(I_power(i_axis) .gt. -1) then
      
          call give_explicit_poly_and_gaussian( Pm_new, Pm_center, ppm, fact_pm, iorder_pm               &
                                              , expo1, expo2, I_power, J_power, I_center, J_center, dim1 )
          pm_inv = 1.d0 / ppm

        endif
        I_power(i_axis) = I_power(i_axis) + 1

        ! ---

        fact_intp  = 2.d0 * expo1
        I_power(i_axis) = I_power(i_axis) + 1
        call give_explicit_poly_and_gaussian( Pp_new, Pp_center, ppp, fact_pp, iorder_pp               &
                                            , expo1, expo2, I_power, J_power, I_center, J_center, dim1 )
        pp_inv = 1.d0 / ppp
        I_power(i_axis) = I_power(i_axis) - 1

        ! ---

        do r = 1, ao_prim_num(k)
          coef3 = coef2 * ao_coef_normalized_ordered_transp(r,k)
          expo3 = ao_expo_ordered_transp(r,k)

          do s = 1, ao_prim_num(l)
            coef4 = coef3 * ao_coef_normalized_ordered_transp(s,l)
            expo4 = ao_expo_ordered_transp(s,l)

            call give_explicit_poly_and_gaussian( Q_new, Q_center, qq, fact_q, iorder_q                    &
                                                , expo3, expo4, K_power, L_power, K_center, L_center, dim1 )
            q_inv = 1.d0 / qq

            ! ---

            if(I_power(i_axis) .gt. 0) then
              integral_m = fact_intm * general_primitive_integral( dim1, Pm_new, Pm_center, fact_pm, ppm, pm_inv, iorder_pm &
                                                                       , Q_new , Q_center , fact_q , qq , q_inv , iorder_q  )
            else
              integral_m = 0.d0
            endif

            ! ---

            integral_p = fact_intp * general_primitive_integral( dim1, Pp_new, Pp_center, fact_pp, ppp, pp_inv, iorder_pp &
                                                                     , Q_new , Q_center , fact_q , qq , q_inv , iorder_q  )

            ! ---

            fst_deriv_eri = fst_deriv_eri + coef4 * (integral_m + integral_p)
            
          enddo ! s
        enddo  ! r
      enddo   ! q
    enddo    ! p

  endif

  ! -----------------------------------------------------------------------------------------------

  ! ---

  ! -----------------------------------------------------------------------------------------------

  if(i_nucl .eq. num_j) then  

    do p = 1, ao_prim_num(i)
      coef1 = ao_coef_normalized_ordered_transp(p,i)
      expo1 = ao_expo_ordered_transp(p,i)

      do q = 1, ao_prim_num(j)
        coef2 = coef1 * ao_coef_normalized_ordered_transp(q,j)
        expo2 = ao_expo_ordered_transp(q,j)

        ! ---

        fact_intm       = -1.d0 * dble(J_power(i_axis))
        J_power(i_axis) = J_power(i_axis) - 1
        if(J_power(i_axis) .gt. -1) then
      
          call give_explicit_poly_and_gaussian( Pm_new, Pm_center, ppm, fact_pm, iorder_pm               &
                                              , expo1, expo2, I_power, J_power, I_center, J_center, dim1 )
          pm_inv = 1.d0 / ppm

        endif
        J_power(i_axis) = J_power(i_axis) + 1

        ! ---

        fact_intp  = 2.d0 * expo2
        J_power(i_axis) = J_power(i_axis) + 1
        call give_explicit_poly_and_gaussian( Pp_new, Pp_center, ppp, fact_pp, iorder_pp               &
                                            , expo1, expo2, I_power, J_power, I_center, J_center, dim1 )
        pp_inv = 1.d0 / ppp
        J_power(i_axis) = J_power(i_axis) - 1

        ! ---

        do r = 1, ao_prim_num(k)
          coef3 = coef2 * ao_coef_normalized_ordered_transp(r,k)
          expo3 = ao_expo_ordered_transp(r,k)

          do s = 1, ao_prim_num(l)
            coef4 = coef3 * ao_coef_normalized_ordered_transp(s,l)
            expo4 = ao_expo_ordered_transp(s,l)

            call give_explicit_poly_and_gaussian( Q_new, Q_center, qq, fact_q, iorder_q                    &
                                                , expo3, expo4, K_power, L_power, K_center, L_center, dim1 )
            q_inv = 1.d0 / qq

            ! ---

            if(J_power(i_axis) .gt. 0) then
              integral_m = fact_intm * general_primitive_integral( dim1, Pm_new, Pm_center, fact_pm, ppm, pm_inv, iorder_pm &
                                                                       , Q_new , Q_center , fact_q , qq , q_inv , iorder_q  )
            else
              integral_m = 0.d0
            endif

            ! ---

            integral_p = fact_intp * general_primitive_integral( dim1, Pp_new, Pp_center, fact_pp, ppp, pp_inv, iorder_pp &
                                                                     , Q_new , Q_center , fact_q , qq , q_inv , iorder_q  )

            ! ---

            fst_deriv_eri = fst_deriv_eri + coef4 * (integral_m + integral_p)
          enddo ! s
        enddo  ! r
      enddo   ! q
    enddo    ! p

  endif

  ! -----------------------------------------------------------------------------------------------

  ! ---

  ! -----------------------------------------------------------------------------------------------

  if(i_nucl .eq. num_k) then  

    do p = 1, ao_prim_num(i)
      coef1 = ao_coef_normalized_ordered_transp(p,i)
      expo1 = ao_expo_ordered_transp(p,i)

      do q = 1, ao_prim_num(j)
        coef2 = coef1 * ao_coef_normalized_ordered_transp(q,j)
        expo2 = ao_expo_ordered_transp(q,j)

        call give_explicit_poly_and_gaussian( P_new, P_center, pp, fact_p, iorder_p                    &
                                            , expo1, expo2, I_power, J_power, I_center, J_center, dim1 )
        p_inv = 1.d0 / pp

        do r = 1, ao_prim_num(k)
          coef3 = coef2 * ao_coef_normalized_ordered_transp(r,k)
          expo3 = ao_expo_ordered_transp(r,k)

          do s = 1, ao_prim_num(l)
            coef4 = coef3 * ao_coef_normalized_ordered_transp(s,l)
            expo4 = ao_expo_ordered_transp(s,l)

            ! ---

            fact_intm       = -1.d0 * dble(K_power(i_axis))
            K_power(i_axis) = K_power(i_axis) - 1
            if(K_power(i_axis) .gt. -1) then
      
              call give_explicit_poly_and_gaussian( Qm_new, Qm_center, qqm, fact_qm, iorder_qm               &
                                                  , expo3, expo4, K_power, L_power, K_center, L_center, dim1 )
              qm_inv = 1.d0 / qqm

            endif
            K_power(i_axis) = K_power(i_axis) + 1

            ! ---

            fact_intp  = 2.d0 * expo3
            K_power(i_axis) = K_power(i_axis) + 1
            call give_explicit_poly_and_gaussian( Qp_new, Qp_center, qqp, fact_qp, iorder_qp               &
                                                , expo3, expo4, K_power, L_power, K_center, L_center, dim1 )
            qp_inv = 1.d0 / qqp
            K_power(i_axis) = K_power(i_axis) - 1

            ! ---

            if(K_power(i_axis) .gt. 0) then
              integral_m = fact_intm * general_primitive_integral( dim1, P_new , P_center , fact_p , pp , p_inv , iorder_p  &
                                                                       , Qm_new, Qm_center, fact_qm, qqm, qm_inv, iorder_qm )
            else
              integral_m = 0.d0
            endif

            ! ---

            integral_p = fact_intp * general_primitive_integral( dim1, P_new , P_center , fact_p , pp , p_inv , iorder_p  &
                                                                     , Qp_new, Qp_center, fact_qp, qqp, qp_inv, iorder_qp )

            ! ---

            fst_deriv_eri = fst_deriv_eri + coef4 * (integral_m + integral_p)
          enddo ! s
        enddo  ! r
      enddo   ! q
    enddo    ! p

  endif

  ! -----------------------------------------------------------------------------------------------

  ! ---

  ! -----------------------------------------------------------------------------------------------

  if(i_nucl .eq. num_l) then  

    do p = 1, ao_prim_num(i)
      coef1 = ao_coef_normalized_ordered_transp(p,i)
      expo1 = ao_expo_ordered_transp(p,i)

      do q = 1, ao_prim_num(j)
        coef2 = coef1 * ao_coef_normalized_ordered_transp(q,j)
        expo2 = ao_expo_ordered_transp(q,j)

        call give_explicit_poly_and_gaussian( P_new, P_center, pp, fact_p, iorder_p                    &
                                            , expo1, expo2, I_power, J_power, I_center, J_center, dim1 )
        p_inv = 1.d0 / pp

        do r = 1, ao_prim_num(k)
          coef3 = coef2 * ao_coef_normalized_ordered_transp(r,k)
          expo3 = ao_expo_ordered_transp(r,k)

          do s = 1, ao_prim_num(l)
            coef4 = coef3 * ao_coef_normalized_ordered_transp(s,l)
            expo4 = ao_expo_ordered_transp(s,l)

            ! ---

            fact_intm       = -1.d0 * dble(L_power(i_axis))
            L_power(i_axis) = L_power(i_axis) - 1
            if(L_power(i_axis) .gt. -1) then
      
              call give_explicit_poly_and_gaussian( Qm_new, Qm_center, qqm, fact_qm, iorder_qm               &
                                                  , expo3, expo4, K_power, L_power, K_center, L_center, dim1 )
              qm_inv = 1.d0 / qqm

            endif
            L_power(i_axis) = L_power(i_axis) + 1

            ! ---

            fact_intp  = 2.d0 * expo4
            L_power(i_axis) = L_power(i_axis) + 1
            call give_explicit_poly_and_gaussian( Qp_new, Qp_center, qqp, fact_qp, iorder_qp               &
                                                , expo3, expo4, K_power, L_power, K_center, L_center, dim1 )
            qp_inv = 1.d0 / qqp
            L_power(i_axis) = L_power(i_axis) - 1

            ! ---

            if(L_power(i_axis) .gt. 0) then
              integral_m = fact_intm * general_primitive_integral( dim1, P_new , P_center , fact_p , pp , p_inv , iorder_p  &
                                                                       , Qm_new, Qm_center, fact_qm, qqm, qm_inv, iorder_qm )
            else
              integral_m = 0.d0
            endif

            ! ---

            integral_p = fact_intp * general_primitive_integral( dim1, P_new , P_center , fact_p , pp , p_inv , iorder_p  &
                                                                     , Qp_new, Qp_center, fact_qp, qqp, qp_inv, iorder_qp )

            ! ---

            fst_deriv_eri = fst_deriv_eri + coef4 * (integral_m + integral_p)
          enddo ! s
        enddo  ! r
      enddo   ! q
    enddo    ! p

  endif


end function fst_deriv_eri

! ---

