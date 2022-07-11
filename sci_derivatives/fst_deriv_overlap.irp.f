
! ---

double precision function fst_deriv_overlap(i_axis, i_nucl, i, j)

  BEGIN_DOC
  ! 
  ! \partial_a S_{ij}
  ! where a = xa, ya or za 
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
  double precision    :: power_tmp, dist_AB, deriv_tmp
  double precision    :: integ_m, integ_p

  fst_deriv_overlap = 0.d0

  if( (i_nucl .ne. ao_nucl(i)) .and. (i_nucl .ne. ao_nucl(j)) ) then
    return
  endif

  do ii = 1, 3

    A_center(ii) = nucl_coord(ao_nucl(j),ii)
    B_center(ii) = nucl_coord(ao_nucl(i),ii)

    power_A(ii)  = ao_power(j,ii)
    power_B(ii)  = ao_power(i,ii)

  enddo

  dist_AB = (A_center(1) - B_center(1)) * (A_center(1) - B_center(1)) &
          + (A_center(2) - B_center(2)) * (A_center(2) - B_center(2)) &
          + (A_center(3) - B_center(3)) * (A_center(3) - B_center(3)) 
  if(dist_AB .lt. 1d-12) return

  
  prim_ao_j = ao_prim_num(j)
  prim_ao_i = ao_prim_num(i)

  dim1 = 100

  ! ---

  if(i_nucl .eq. ao_nucl(j)) then

! !$OMP PARALLEL DO SCHEDULE(GUIDED) &
! !$OMP DEFAULT(NONE)                &
! !$OMP PRIVATE( A_center, B_center, power_A, power_B, i_axis, overlap_x, overlap_y, overlap_z, overlap &
! !$OMP        , alpha, beta, i, j, n, l, c, power_tmp, integ_m, integ_p, deriv_tmp )                   &
! !$OMP SHARED(prim_ao_i, prim_ao_j, ao_coef_normalized_ordered_transp, ao_expo_ordered_transp, dim1)

    deriv_tmp = 0.d0

    do n = 1, prim_ao_j
      alpha = ao_expo_ordered_transp(n,j)

      do l = 1, prim_ao_i
        beta = ao_expo_ordered_transp(l,i)
        c    = ao_coef_normalized_ordered_transp(n,j) * ao_coef_normalized_ordered_transp(l,i)

        ! ---

        power_tmp       = -1.d0 * dble(power_A(i_axis))
        power_A(i_axis) = power_A(i_axis) - 1
        if(power_A(i_axis) > -1) then
          call overlap_gaussian_xyz(A_center, B_center, alpha, beta, power_A, power_B, overlap_x, overlap_y, overlap_z, overlap, dim1)
          integ_m = power_tmp * overlap
        else
          integ_m = 0.d0
        endif

        ! ---

        power_A(i_axis) = power_A(i_axis) + 2
        call overlap_gaussian_xyz(A_center, B_center, alpha, beta, power_A, power_B, overlap_x, overlap_y, overlap_z, overlap, dim1)
        integ_p = 2.d0 * alpha * overlap

        ! ---

        power_A(i_axis) = power_A(i_axis) - 1

        deriv_tmp = deriv_tmp + c * (integ_m + integ_p)
      enddo
    enddo
! !$OMP END PARALLEL DO

    fst_deriv_overlap = fst_deriv_overlap + deriv_tmp

  endif

  ! ---

  if(i_nucl .eq. ao_nucl(i)) then

! !$OMP PARALLEL DO SCHEDULE(GUIDED) &
! !$OMP DEFAULT(NONE)                &
! !$OMP PRIVATE( A_center, B_center, power_A, power_B, i_axis, overlap_x, overlap_y, overlap_z, overlap &
! !$OMP        , alpha, beta, i, j, n, l, c, power_tmp, integ_m, integ_p, deriv_tmp )                   &
! !$OMP SHARED(prim_ao_i, prim_ao_j, ao_coef_normalized_ordered_transp, ao_expo_ordered_transp, dim1)

    deriv_tmp = 0.d0

    do n = 1, prim_ao_j
      alpha = ao_expo_ordered_transp(n,j)

      do l = 1, prim_ao_i
        beta = ao_expo_ordered_transp(l,i)
        c    = ao_coef_normalized_ordered_transp(n,j) * ao_coef_normalized_ordered_transp(l,i)

        ! ---

        power_tmp       = -1.d0 * dble(power_B(i_axis))
        power_B(i_axis) = power_B(i_axis) - 1
        if(power_B(i_axis) > -1) then
          call overlap_gaussian_xyz(A_center, B_center, alpha, beta, power_A, power_B, overlap_x, overlap_y, overlap_z, overlap, dim1)
          integ_m = power_tmp * overlap
        else
          integ_m = 0.d0
        endif

        ! ---

        power_B(i_axis) = power_B(i_axis) + 2
        call overlap_gaussian_xyz(A_center, B_center, alpha, beta, power_A, power_B, overlap_x, overlap_y, overlap_z, overlap, dim1)
        integ_p = 2.d0 * beta * overlap

        ! ---

        power_B(i_axis) = power_B(i_axis) - 1

        deriv_tmp = deriv_tmp + c * (integ_m + integ_p)
      enddo
    enddo
! !$OMP END PARALLEL DO

    fst_deriv_overlap = fst_deriv_overlap + deriv_tmp

  endif

end function fst_deriv_overlap

! ---

