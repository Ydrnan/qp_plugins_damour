
! ---

double precision function fst_deriv_nai(i_axis, i_nucl, i, j)

  BEGIN_DOC
  ! 
  ! \partial_a < Gi | -\sum_n Zn / |r - Rn| | Gj >
  ! where a = xa, ya or za 
  !
  END_DOC

  implicit none
  
  integer, intent(in) :: i_axis ! index of axis; 1 : \partial_x, 2 : \partial_y, 3 : \partial_z
  integer, intent(in) :: i_nucl ! index of nuclueus
  integer, intent(in) :: i, j

  integer             :: ii, k, l, m, n, power_A(3), power_B(3), n_pt_in
  integer             :: prim_ao_i, prim_ao_j
  double precision    :: alpha, beta, c, A_center(3), B_center(3), C_center(3)
  double precision    :: power_tmp, deriv_tmp, Z, p
  double precision    :: A_center_x, B_center_x, C_center_x, P_center_x
  double precision    :: integ_m, integ_p, integ_sum
  double precision    :: integ_0, integ_1, integ_2

  double precision    :: NAI_standard, NAI_shift_t2


  fst_deriv_nai = 0.d0

  do ii = 1, 3

    A_center(ii) = nucl_coord(ao_nucl(j),ii)
    B_center(ii) = nucl_coord(ao_nucl(i),ii)

    power_A(ii)  = ao_power(j,ii)
    power_B(ii)  = ao_power(i,ii)

  enddo

  n_pt_in   = n_pt_max_integrals
  prim_ao_j = ao_prim_num(j)
  prim_ao_i = ao_prim_num(i)

  ! ---

  if(i_nucl .eq. ao_nucl(j)) then


! !$OMP PARALLEL                                                                                    &
! !$OMP DEFAULT (NONE)                                                                              &
! !$OMP PRIVATE ( i, j, k, l, m, i_axis, alpha, beta, c, A_center, B_center, C_center               &
! !$OMP         , power_A, power_B, power_tmp, Z, integ_sum, integ_m, integ_p, deriv_tmp, n_pt_in   &
! !$OMP         , prim_ao_i, prim_ao_j )                                                            &
! !$OMP SHARED ( ao_expo_ordered_transp, ao_coef_normalized_ordered_transp    &
! !$OMP        , nucl_coord, nucl_num, nucl_charge )
!
! !$OMP DO SCHEDULE (dynamic)

    deriv_tmp = 0.d0

    do l = 1, prim_ao_j
      alpha = ao_expo_ordered_transp(l,j)

      do m = 1, prim_ao_i
        beta = ao_expo_ordered_transp(m,i)
        c    = ao_coef_normalized_ordered_transp(l,j) * ao_coef_normalized_ordered_transp(m,i)

        ! --- 

        power_tmp       = -1.d0 * dble(power_A(i_axis))
        power_A(i_axis) = power_A(i_axis) - 1
        if(power_A(i_axis) > -1) then

          integ_sum = 0.d0
          do k = 1, nucl_num
            Z             = nucl_charge(k)
            C_center(1:3) = nucl_coord(k,1:3)

            integ_sum = integ_sum &
                      - Z * NAI_standard(A_center, B_center, power_A, power_B, alpha, beta, C_center, n_pt_in)
          enddo

          integ_m = power_tmp * integ_sum

        else

          integ_m = 0.d0

        endif

        ! --- 

        power_A(i_axis) = power_A(i_axis) + 2
        integ_sum = 0.d0
        do k = 1, nucl_num
          Z             = nucl_charge(k)
          C_center(1:3) = nucl_coord(k,1:3)

          integ_sum = integ_sum &
                    - Z * NAI_standard(A_center, B_center, power_A, power_B, alpha, beta, C_center, n_pt_in)
        enddo

        integ_p = 2.d0 * alpha * integ_sum

        ! --- 

        power_A(i_axis) = power_A(i_axis) - 1

        deriv_tmp = deriv_tmp + c * (integ_m + integ_p)
      enddo
    enddo

! !$OMP END DO
! !$OMP END PARALLEL

    fst_deriv_nai = fst_deriv_nai + deriv_tmp

  endif

  ! ---

  if(i_nucl .eq. ao_nucl(i)) then

! !$OMP PARALLEL       &
! !$OMP DEFAULT (NONE) &
! !$OMP PRIVATE ( i, j, k, l, m, i_axis, alpha, beta, c, A_center, B_center, C_center            &
! !$OMP         , power_A, power_B, power_tmp, Z, integ_sum, integ_m, integ_p, deriv_tmp )       &
! !$OMP SHARED ( prim_ao_j, prim_ao_i, ao_expo_ordered_transp, ao_coef_normalized_ordered_transp &
! !$OMP        , nucl_coord, nucl_num, nucl_charge, n_pt_in )
!
! !$OMP DO SCHEDULE (dynamic)

    deriv_tmp = 0.d0

    do l = 1, prim_ao_j
      alpha = ao_expo_ordered_transp(l,j)

      do m = 1, prim_ao_i
        beta = ao_expo_ordered_transp(m,i)
        c    = ao_coef_normalized_ordered_transp(l,j) * ao_coef_normalized_ordered_transp(m,i)

        ! --- 

        power_tmp       = -1.d0 * dble(power_B(i_axis))
        power_B(i_axis) = power_B(i_axis) - 1
        if(power_B(i_axis) > -1) then

          integ_sum = 0.d0
          do k = 1, nucl_num
            Z             = nucl_charge(k)
            C_center(1:3) = nucl_coord(k,1:3)

            integ_sum = integ_sum &
                      - Z * NAI_standard(A_center, B_center, power_A, power_B, alpha, beta, C_center, n_pt_in)
          enddo

          integ_m = power_tmp * integ_sum

        else

          integ_m = 0.d0

        endif

        ! --- 

        power_B(i_axis) = power_B(i_axis) + 2
        integ_sum = 0.d0
        do k = 1, nucl_num
          Z             = nucl_charge(k)
          C_center(1:3) = nucl_coord(k,1:3)

          integ_sum = integ_sum &
                    - Z * NAI_standard(A_center, B_center, power_A, power_B, alpha, beta, C_center, n_pt_in)
        enddo

        integ_p = 2.d0 * beta * integ_sum

        ! --- 

        power_B(i_axis) = power_B(i_axis) - 1

        deriv_tmp = deriv_tmp + c * (integ_m + integ_p)
      enddo
    enddo

! !$OMP END DO
! !$OMP END PARALLEL

    fst_deriv_nai = fst_deriv_nai + deriv_tmp

  endif

  ! ---

  k             = i_nucl
  Z             = nucl_charge(k)
  C_center(1:3) = nucl_coord(k,1:3)

  A_center_x = A_center(i_axis)
  B_center_x = B_center(i_axis)
  C_center_x = C_center(i_axis)

! !$OMP PARALLEL                                                                                     &
! !$OMP DEFAULT (NONE)                                                                               &
! !$OMP PRIVATE ( i, j, k, l, m, i_axis, alpha, beta, p, c, A_center, B_center, C_center, P_center_x &
! !$OMP         , power_A, power_B, power_tmp, integ_0, integ_1, integ_2, deriv_tmp, n_pt_in )       &
! !$OMP SHARED ( prim_ao_j, prim_ao_i, ao_expo_ordered_transp, ao_coef_normalized_ordered_transp     &
! !$OMP        , A_center_x, B_center_x, C_center_x, Z, nucl_coord, nucl_num, nucl_charge )
!
! !$OMP DO SCHEDULE (dynamic)

  deriv_tmp = 0.d0

  do l = 1, prim_ao_j
    alpha = ao_expo_ordered_transp(l,j)

    do m = 1, prim_ao_i
      beta = ao_expo_ordered_transp(m,i)
      c    = ao_coef_normalized_ordered_transp(l,j) * ao_coef_normalized_ordered_transp(m,i)

      ! --- 

      p          = alpha + beta
      P_center_x = (alpha * A_center_x + beta * B_center_x) / p
      integ_0    = 2.d0 * p * (P_center_x - C_center_x) &
                 * NAI_shift_t2(A_center, B_center, power_A, power_B, alpha, beta, C_center, n_pt_in)

      ! --- 

      power_tmp       = dble(power_A(i_axis))
      power_A(i_axis) = power_A(i_axis) - 1
      if(power_A(i_axis) > -1) then

        integ_1 = power_tmp &
                * NAI_shift_t2(A_center, B_center, power_A, power_B, alpha, beta, C_center, n_pt_in)

      else

        integ_1 = 0.d0

      endif

      power_A(i_axis) = power_A(i_axis) + 1

      ! --- 

      power_tmp       = dble(power_B(i_axis))
      power_B(i_axis) = power_B(i_axis) - 1
      if(power_B(i_axis) > -1) then

        integ_2 = power_tmp &
                * NAI_shift_t2(A_center, B_center, power_A, power_B, alpha, beta, C_center, n_pt_in)

      else

        integ_2 = 0.d0

      endif

      power_B(i_axis) = power_B(i_axis) + 1

      ! --- 

      deriv_tmp = deriv_tmp - Z * c * (integ_0 + integ_1 + integ_2)
    enddo
  enddo

! !$OMP END DO
! !$OMP END PARALLEL

  fst_deriv_nai = fst_deriv_nai + deriv_tmp

end function fst_deriv_nai

! ---

double precision function NAI_standard(A_center, B_center, power_A, power_B, alpha, beta, C_center, n_pt_in)

  implicit none
  include 'utils/constants.include.F'

  integer,          intent(in) :: n_pt_in, power_A(3), power_B(3)
  double precision, intent(in) :: C_center(3), A_center(3), B_center(3), alpha, beta

  integer                      :: i, j, k, l, n_pt, n_pt_out
  double precision             :: P_center(3)
  double precision             :: d(0:n_pt_in), coeff, rho, const, p, p_inv, factor
  double precision             :: const_factor, dist_integral, tmp, accu, dist_AB, dist_AC

  double precision             :: rint


  dist_AB = (A_center(1) - B_center(1)) * (A_center(1) - B_center(1)) &
          + (A_center(2) - B_center(2)) * (A_center(2) - B_center(2)) &
          + (A_center(3) - B_center(3)) * (A_center(3) - B_center(3)) 

  dist_AC = (A_center(1) - C_center(1)) * (A_center(1) - C_center(1)) &
          + (A_center(2) - C_center(2)) * (A_center(2) - C_center(2)) &
          + (A_center(3) - C_center(3)) * (A_center(3) - C_center(3)) 
  
  if((dist_AB .lt. 1d-12) .and. (dist_AC .lt. 1d-12)) then
    NAI_standard = 0.d0
    return
  endif

  p     = alpha + beta
  p_inv = 1.d0 / p
  rho   = alpha * beta * p_inv

  dist_integral = 0.d0
  do i = 1, 3
    P_center(i)    = (alpha * A_center(i) + beta * B_center(i)) * p_inv
    dist_integral += (P_center(i) - C_center(i)) * (P_center(i) - C_center(i))
  enddo
  const_factor = dist_AB * rho
  const        = p * dist_integral
  if(const_factor > 80.d0) then
    NAI_standard = 0.d0
    return
  endif
  factor = dexp(-const_factor)
  coeff  = dtwo_pi * factor * p_inv

  do i = 0, n_pt_in
    d(i) = 0.d0
  enddo

  n_pt = 2 * ( (power_A(1) + power_B(1)) + (power_A(2) + power_B(2)) + (power_A(3) + power_B(3)) )
  if(n_pt == 0) then
    NAI_standard = coeff * rint(0, const)
    return
  endif

  call give_polynomial_mult_center_one_e(A_center, B_center, alpha, beta, power_A, power_B, C_center, n_pt_in, d, n_pt_out)
  if(n_pt_out < 0) then
    NAI_standard = 0.d0
    return
  endif

  accu = 0.d0
  do i = 0, n_pt_out, 2
    accu +=  d(i) * rint(shiftr(i, 1), const)
  enddo
  NAI_standard = accu * coeff

end function NAI_standard

! ---

double precision function NAI_shift_t2(A_center, B_center, power_A, power_B, alpha, beta, C_center, n_pt_in)

  implicit none
  include 'utils/constants.include.F'

  integer,          intent(in) :: n_pt_in, power_A(3), power_B(3)
  double precision, intent(in) :: C_center(3), A_center(3), B_center(3), alpha, beta

  integer                      :: i, j, k, l, n_pt, n_pt_out
  double precision             :: P_center(3)
  double precision             :: d(0:n_pt_in), coeff, rho, const, p, p_inv, factor
  double precision             :: const_factor, dist_integral, tmp
  double precision             :: accu, dist_AB, dist_AC

  double precision             :: rint


  dist_AB = (A_center(1) - B_center(1)) * (A_center(1) - B_center(1)) &
          + (A_center(2) - B_center(2)) * (A_center(2) - B_center(2)) &
          + (A_center(3) - B_center(3)) * (A_center(3) - B_center(3)) 

  dist_AC = (A_center(1) - C_center(1)) * (A_center(1) - C_center(1)) &
          + (A_center(2) - C_center(2)) * (A_center(2) - C_center(2)) &
          + (A_center(3) - C_center(3)) * (A_center(3) - C_center(3)) 
  
  if((dist_AB .lt. 1d-12) .and. (dist_AC .lt. 1d-12)) then
    NAI_shift_t2 = 0.d0
    return
  endif

  p     = alpha + beta
  p_inv = 1.d0 / p
  rho   = alpha * beta * p_inv

  dist_integral = 0.d0
  do i = 1, 3
    P_center(i)    = (alpha * A_center(i) + beta * B_center(i)) * p_inv
    dist_integral += (P_center(i) - C_center(i)) * (P_center(i) - C_center(i))
  enddo
  const_factor = dist_AB * rho
  const        = p * dist_integral
  if(const_factor > 80.d0) then
    NAI_shift_t2 = 0.d0
    return
  endif
  factor = dexp(-const_factor)
  coeff  = dtwo_pi * factor * p_inv

  do i = 0, n_pt_in
    d(i) = 0.d0
  enddo

  n_pt = 2 * ( (power_A(1) + power_B(1)) + (power_A(2) + power_B(2)) + (power_A(3) + power_B(3)) )
  if(n_pt == 0) then
    NAI_shift_t2 = coeff * rint(1, const)
    return
  endif

  call give_polynomial_mult_center_one_e(A_center, B_center, alpha, beta, power_A, power_B, C_center, n_pt_in, d, n_pt_out)
  if(n_pt_out < 0) then
    NAI_shift_t2 = 0.d0
    return
  endif

  accu = 0.d0
  do i = 0, n_pt_out, 2
    accu += d(i) * rint(shiftr(i, 1) + 1, const)
  enddo
  NAI_shift_t2 = accu * coeff

end function NAI_shift_t2

! ---

!subroutine give_polynomial_mult_center_one_e(A_center, B_center, alpha, beta, power_A, power_B, C_center, n_pt_in, d, n_pt_out)
!
!  BEGIN_DOC
!  ! Returns the explicit polynomial in terms of the "t" variable of the following
!  !
!  ! $I_{x1}(a_x, d_x, p, q) \times I_{x1}(a_y, d_y, p, q) \times I_{x1}(a_z, d_z, p, q)$.
!  END_DOC
!
!  implicit none
!
!  integer,          intent(in)  :: n_pt_in, power_A(3), power_B(3)
!  double precision, intent(in)  :: alpha, beta, A_center(3), B_center(3), C_center(3)
!  integer,          intent(out) :: n_pt_out
!  double precision, intent(out) :: d(0:n_pt_in)
!
!  integer                       :: a_x, b_x, a_y, b_y, a_z, b_z
!  integer                       :: n_pt1, n_pt2, n_pt3, i, n_pt_tmp
!  double precision              :: p, P_center(3), rho, p_inv, p_inv_2
!  double precision              :: d1(0:n_pt_in), d2(0:n_pt_in), d3(0:n_pt_in)
!  double precision              :: R1x(0:2), B01(0:2), R1xp(0:2), R2x(0:2)
!
!  ASSERT (n_pt_in > 1)
!
!  p       = alpha + beta
!  p_inv   = 1.d0 / p
!  p_inv_2 = 0.5d0 * p_inv
!
!  do i = 1, 3
!    P_center(i) = (alpha * A_center(i) + beta * B_center(i)) * p_inv
!  enddo
!
!  R1x(0) =  (P_center(1) - A_center(1))
!  R1x(1) =  0.d0
!  R1x(2) = -(P_center(1) - C_center(1))
!
!  R1xp(0) =  (P_center(1) - B_center(1))
!  R1xp(1) =  0.d0
!  R1xp(2) = -(P_center(1) - C_center(1))
!
!  R2x(0) =  p_inv_2
!  R2x(1) =  0.d0
!  R2x(2) = -p_inv_2
!
!  do i = 0, n_pt_in
!    d(i)  = 0.d0
!    d1(i) = 0.d0
!    d2(i) = 0.d0
!    d3(i) = 0.d0
!  enddo
!
!  n_pt1 = n_pt_in
!  n_pt2 = n_pt_in
!  n_pt3 = n_pt_in
!
!  a_x = power_A(1)
!  b_x = power_B(1)
!  call I_x1_pol_mult_one_e(a_x, b_x, R1x, R1xp, R2x, d1, n_pt1, n_pt_in)
!
!  if(n_pt1 < 0) then
!    n_pt_out = -1
!    do i = 0, n_pt_in
!      d(i) = 0.d0
!    enddo
!    return
!  endif
!
!  R1x(0) =  (P_center(2) - A_center(2))
!  R1x(1) =  0.d0
!  R1x(2) = -(P_center(2) - C_center(2))
!
!  R1xp(0) =  (P_center(2) - B_center(2))
!  R1xp(1) =  0.d0
!  R1xp(2) = -(P_center(2) - C_center(2))
!
!  a_y = power_A(2)
!  b_y = power_B(2)
!
!  call I_x1_pol_mult_one_e(a_y, b_y, R1x, R1xp, R2x, d2, n_pt2, n_pt_in)
!  if(n_pt2 < 0) then
!    n_pt_out = -1
!    do i = 0, n_pt_in
!      d(i) = 0.d0
!    enddo
!    return
!  endif
!
!  R1x(0) =  (P_center(3) - A_center(3))
!  R1x(1) =  0.d0
!  R1x(2) = -(P_center(3) - C_center(3))
!
!  R1xp(0) =  (P_center(3) - B_center(3))
!  R1xp(1) =  0.d0
!  R1xp(2) = -(P_center(3) - C_center(3))
!
!  a_z = power_A(3)
!  b_z = power_B(3)
!
!  call I_x1_pol_mult_one_e(a_z, b_z, R1x, R1xp, R2x, d3, n_pt3, n_pt_in)
!  if(n_pt3 < 0) then
!    n_pt_out = -1
!    do i = 0, n_pt_in
!      d(i) = 0.d0
!    enddo
!    return
!  endif
!
!  n_pt_tmp = 0
!  call multiply_poly(d1, n_pt1, d2, n_pt2, d, n_pt_tmp)
!  do i = 0, n_pt_tmp
!    d1(i) = 0.d0
!  enddo
!  n_pt_out = 0
!  call multiply_poly(d, n_pt_tmp, d3, n_pt3, d1, n_pt_out)
!  do i = 0, n_pt_out
!    d(i) = d1(i)
!  enddo
!
!end subroutine give_polynomial_mult_center_one_e

! ---

!recursive subroutine I_x1_pol_mult_one_e(a, c, R1x, R1xp, R2x, d, nd, n_pt_in)
!
!  BEGIN_DOC
!  !  Recursive routine involved in the electron-nucleus potential
!  END_DOC
!
!  implicit none
!  include 'utils/constants.include.F'
!
!  integer,          intent(in)    :: n_pt_in, a, c
!  double precision, intent(in)    :: R1x(0:2), R1xp(0:2), R2x(0:2)
!  integer,          intent(inout) :: nd
!  double precision, intent(inout) :: d(0:n_pt_in)
!
!  integer                         :: nx, ix, dim, ny
!  double precision                :: X(0:max_dim), Y(0:max_dim)
!  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: X, Y
!
!  dim = n_pt_in
!
!  if( (a==0) .and. (c==0)) then
!
!    nd = 0
!    d(0) = 1.d0
!    return
!
!  elseif( (c<0).or.(nd<0) )then
!
!    nd = -1
!    return
!
!  else if ((a==0).and.(c.ne.0)) then
!
!    call I_x2_pol_mult_one_e(c, R1x, R1xp, R2x, d, nd, n_pt_in)
!
!  else if (a==1) then
!
!    nx = nd
!    do ix = 0, n_pt_in
!      X(ix) = 0.d0
!      Y(ix) = 0.d0
!    enddo
!    call I_x2_pol_mult_one_e(c-1, R1x, R1xp, R2x, X, nx, n_pt_in)
!    do ix = 0, nx
!      X(ix) *= dble(c)
!    enddo
!    call multiply_poly(X, nx, R2x, 2, d, nd)
!    ny = 0
!    call I_x2_pol_mult_one_e(c, R1x, R1xp, R2x, Y, ny, n_pt_in)
!    call multiply_poly(Y, ny, R1x, 2, d, nd)
!
!  else
!
!    do ix = 0, n_pt_in
!      X(ix) = 0.d0
!      Y(ix) = 0.d0
!    enddo
!    nx = 0
!    call I_x1_pol_mult_one_e(a-2, c, R1x, R1xp, R2x, X, nx, n_pt_in)
!    do ix = 0, nx
!      X(ix) *= dble(a-1)
!    enddo
!    call multiply_poly(X, nx, R2x, 2, d, nd)
!
!    nx = nd
!    do ix = 0, n_pt_in
!      X(ix) = 0.d0
!    enddo
!    call I_x1_pol_mult_one_e(a-1, c-1, R1x, R1xp, R2x, X, nx, n_pt_in)
!    do ix = 0, nx
!      X(ix) *= dble(c)
!    enddo
!    call multiply_poly(X, nx, R2x, 2, d, nd)
!    ny = 0
!    call I_x1_pol_mult_one_e(a-1, c, R1x, R1xp, R2x, Y, ny, n_pt_in)
!    call multiply_poly(Y, ny, R1x, 2, d, nd)
!
!  endif
!
!end subroutine I_x1_pol_mult_one_e

! ---

!recursive subroutine I_x2_pol_mult_one_e(c, R1x, R1xp, R2x, d, nd, dim)
!
!  BEGIN_DOC
!  !  Recursive routine involved in the electron-nucleus potential
!  END_DOC
!
!  implicit none
!  include 'utils/constants.include.F'
!
!  integer,          intent(in)    :: c, dim
!  double precision, intent(in)    :: R1x(0:2), R1xp(0:2), R2x(0:2)
!  integer,          intent(inout) :: nd
!  double precision, intent(inout) :: d(0:max_dim)
!
!  integer                         :: i, nx, ix, ny
!  double precision                :: X(0:max_dim), Y(0:max_dim)
!  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: X, Y
!
!  if(c==0) then
!    nd = 0
!    d(0) = 1.d0
!    return
!  elseif ((nd<0).or.(c<0))then
!    nd = -1
!    return
!  else
!    do ix = 0, dim
!      X(ix) = 0.d0
!      Y(ix) = 0.d0
!    enddo
!    nx = 0
!    call I_x1_pol_mult_one_e(0, c-2, R1x, R1xp, R2x, X, nx, dim)
!    do ix = 0, nx
!      X(ix) *= dble(c-1)
!    enddo
!    call multiply_poly(X, nx, R2x, 2, d, nd)
!    ny = 0
!    do ix = 0, dim
!      Y(ix) = 0.d0
!    enddo
!
!    call I_x1_pol_mult_one_e(0, c-1, R1x, R1xp, R2x, Y, ny, dim)
!    if(ny .ge. 0)then
!      call multiply_poly(Y, ny, R1xp, 2, d, nd)
!    endif
!  endif
!
!end subroutine I_x2_pol_mult_one_e

! ---

