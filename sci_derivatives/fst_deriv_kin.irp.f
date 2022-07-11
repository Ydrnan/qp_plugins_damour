
! ---

double precision function fst_deriv_kin(i_axis, i_nucl, i, j)

  BEGIN_DOC
  ! 
  ! -0.5 x \partial_n < G_i | \partial_x^2 G_j + \partial_y^2 G_j + \partial_z^2 G_j >
  ! where n = xn, yn or zn 
  !
  END_DOC

  implicit none
  
  integer, intent(in) :: i_axis ! index of axis; 1 : \partial_x, 2 : \partial_y, 3 : \partial_z
  integer, intent(in) :: i_nucl ! index of nuclueus
  integer, intent(in) :: i, j

  integer             :: ii, n, l, dim1, power_A(3), power_B(3)
  integer             :: prim_ao_i, prim_ao_j
  double precision    :: overlap_x, overlap_y, overlap_z, overlap
  double precision    :: alpha, beta, c, A_center(3), B_center(3)
  double precision    :: power_tmp, dist_AB, fact_tmp
  double precision    :: integ1, integ2, integ_tmp, Ix, Iy, Iz

  fst_deriv_kin = 0.d0

  if( (i_nucl .ne. ao_nucl(i)) .and. (i_nucl .ne. ao_nucl(j)) ) then
    return
  endif

  do ii = 1, 3

    A_center(ii) = nucl_coord(ao_nucl(i),ii)
    B_center(ii) = nucl_coord(ao_nucl(j),ii)

    power_A(ii)  = ao_power(i,ii)
    power_B(ii)  = ao_power(j,ii)

  enddo

  dist_AB = (A_center(1) - B_center(1)) * (A_center(1) - B_center(1)) &
          + (A_center(2) - B_center(2)) * (A_center(2) - B_center(2)) &
          + (A_center(3) - B_center(3)) * (A_center(3) - B_center(3)) 
  if(dist_AB .lt. 1d-12) return

  
  prim_ao_i = ao_prim_num(i)
  prim_ao_j = ao_prim_num(j)

  dim1 = 100

  ! ---
  ! < \partial_n G_i | \partial_x^2 G_j + \partial_y^2 G_j + \partial_z^2 G_j >

  if(i_nucl .eq. ao_nucl(i)) then

! !$OMP PARALLEL DO SCHEDULE(GUIDED) &
! !$OMP DEFAULT(NONE)                &
! !$OMP PRIVATE( A_center, B_center, power_A, power_B   &
! !$OMP        , Ix, Iy, Iz, integ1, integ2, integ_tmp  &
! !$OMP        , alpha, beta, c, i, j, n, l, fact_tmp ) &
! !$OMP SHARED( prim_ao_i, prim_ao_j, dim1, i_axis      &
! !$OMP       , ao_coef_normalized_ordered_transp, ao_expo_ordered_transp )
    
    integ_tmp = 0.d0

    do n = 1, prim_ao_i
      alpha = ao_expo_ordered_transp(n,i)

      do l = 1, prim_ao_j
        beta = ao_expo_ordered_transp(l,j)
        c = ao_coef_normalized_ordered_transp(n,i) * ao_coef_normalized_ordered_transp(l,j)

        ! ---
  
        if(power_A(i_axis) .gt. 0) then

          fact_tmp = -dble(power_A(i_axis))

          power_A(i_axis) = power_A(i_axis) - 1
          call overlap_3d_G_Gxx(A_center, B_center, alpha, beta, power_A, power_B, dim1, Ix, Iy, Iz)
          power_A(i_axis) = power_A(i_axis) + 1

          integ1 = fact_tmp * (Ix + Iy + Iz)

        else

          integ1 = 0.d0

        endif

        ! ---

        fact_tmp = 2.d0 * alpha

        power_A(i_axis) = power_A(i_axis) + 1
        call overlap_3d_G_Gxx(A_center, B_center, alpha, beta, power_A, power_B, dim1, Ix, Iy, Iz)
        power_A(i_axis) = power_A(i_axis) - 1

        integ2 = fact_tmp * (Ix + Iy + Iz)

        ! ---

        integ_tmp = integ_tmp + c * (integ1 + integ2)

      enddo
    enddo
! !$OMP END PARALLEL DO

    fst_deriv_kin = fst_deriv_kin -0.5d0 * integ_tmp

  endif

  ! ---
  ! < G_i | \partial_n [ \partial_x^2 G_j + \partial_y^2 G_j + \partial_z^2 G_j ] >

  if(i_nucl .eq. ao_nucl(j)) then

! !$OMP PARALLEL DO SCHEDULE(GUIDED) &
! !$OMP DEFAULT(NONE)                &
! !$OMP PRIVATE( A_center, B_center, power_A, power_B   &
! !$OMP        , Ix, Iy, Iz, integ_tmp                  &
! !$OMP        , alpha, beta, c, i, j, n, l )           &
! !$OMP SHARED( prim_ao_i, prim_ao_j, dim1, i_axis      &
! !$OMP       , ao_coef_normalized_ordered_transp, ao_expo_ordered_transp )
    
    integ_tmp = 0.d0

    do n = 1, prim_ao_i
      alpha = ao_expo_ordered_transp(n,i)

      do l = 1, prim_ao_j
        beta = ao_expo_ordered_transp(l,j)
        c = ao_coef_normalized_ordered_transp(n,i) * ao_coef_normalized_ordered_transp(l,j)

        ! ---
  
        if(i_axis .eq. 1) then
          call overlap_3d_G_GXbxx(A_center, B_center, alpha, beta, power_A, power_B, dim1, Ix, Iy, Iz)
        elseif(i_axis .eq. 2) then
          call overlap_3d_G_GYbxx(A_center, B_center, alpha, beta, power_A, power_B, dim1, Ix, Iy, Iz)
        elseif(i_axis .eq. 3) then
          call overlap_3d_G_GZbxx(A_center, B_center, alpha, beta, power_A, power_B, dim1, Ix, Iy, Iz)
        else
          print*, ' i_axis = 1, 2 or 3.', i_axis
        endif

        integ_tmp = integ_tmp + c * (Ix + Iy + Iz)

      enddo
    enddo
! !$OMP END PARALLEL DO

    fst_deriv_kin = fst_deriv_kin -0.5d0 * integ_tmp

  endif


end function fst_deriv_kin

! ---

