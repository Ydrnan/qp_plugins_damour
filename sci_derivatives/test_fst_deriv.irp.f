! ---

program test_fst_deriv

  implicit none

  !call test_fst_deriv_overlap()
  !call test_fst_deriv_kin()
  !call test_fst_deriv_nai()
  !call test_fst_deriv_eri()
  call test_fst_deriv_eri_schwartz()

end

! ---

subroutine test_fst_deriv_overlap()

  implicit none
  integer                       :: i, j, i_axis, i_nucl
  double precision              :: deriv_ex, deriv_nm
  double precision              :: delta, accu_abs, eps, inv_eps
  double precision, allocatable :: Int_1(:,:), Int_2(:,:)

  double precision              :: fst_deriv_overlap

  print*, ' check fst_deriv_overlap:'

  !print *, nucl_num
  !print *, ao_nucl

  eps     = 1d-8
  inv_eps = 1.d0 / eps

  allocate( Int_1(ao_num,ao_num), Int_2(ao_num,ao_num) )

  accu_abs = 0.d0
  accu_abs = 0.d0

  do i_nucl = 1, nucl_num
    do i_axis = 1, 3

      !print*, nucl_coord
      nucl_coord(i_nucl,i_axis) = nucl_coord(i_nucl,i_axis) + eps
      TOUCH nucl_coord
      !print*, nucl_coord
      Int_2(1:ao_num,1:ao_num) = ao_overlap(1:ao_num,1:ao_num)

      nucl_coord(i_nucl,i_axis) = nucl_coord(i_nucl,i_axis) - 2.d0 * eps
      TOUCH nucl_coord
      Int_1(1:ao_num,1:ao_num) = ao_overlap(1:ao_num,1:ao_num)

      nucl_coord(i_nucl,i_axis) = nucl_coord(i_nucl,i_axis) + eps
      TOUCH nucl_coord

      do j = 1, ao_num
        do i = 1, ao_num

          deriv_nm = 0.5d0 * inv_eps * (Int_2(i,j) - Int_1(i,j))
          deriv_ex = fst_deriv_overlap(i_axis, i_nucl, i, j)

          delta = dabs(deriv_nm - deriv_ex)
          accu_abs += delta

          if(delta .gt. 1.d-6) then
            print*, ' problem on: '
            print*, i_nucl, i_axis, i, j
            print*, deriv_nm, deriv_ex, delta
            !stop
          endif

        enddo
      enddo
    enddo
  enddo

  deallocate( Int_1, Int_2 )

  print*,'accu_abs = ', accu_abs

end subroutine test_fst_deriv_overlap

! ---

subroutine test_fst_deriv_kin()

  implicit none
  integer                       :: i, j, i_axis, i_nucl
  double precision              :: deriv_ex, deriv_nm
  double precision              :: delta, accu_abs, eps, inv_eps
  double precision, allocatable :: Int_1(:,:), Int_2(:,:)

  double precision              :: fst_deriv_kin

  print*, ' check fst_deriv_kin:'

  eps     = 1d-8
  inv_eps = 1.d0 / eps

  allocate( Int_1(ao_num,ao_num), Int_2(ao_num,ao_num) )

  accu_abs = 0.d0
  accu_abs = 0.d0

  do i_nucl = 1, nucl_num
    do i_axis = 1, 3

      nucl_coord(i_nucl,i_axis) = nucl_coord(i_nucl,i_axis) + eps
      TOUCH nucl_coord
      Int_2(1:ao_num,1:ao_num) = ao_kinetic_integrals(1:ao_num,1:ao_num)

      nucl_coord(i_nucl,i_axis) = nucl_coord(i_nucl,i_axis) - 2.d0 * eps
      TOUCH nucl_coord
      Int_1(1:ao_num,1:ao_num) = ao_kinetic_integrals(1:ao_num,1:ao_num)

      nucl_coord(i_nucl,i_axis) = nucl_coord(i_nucl,i_axis) + eps
      TOUCH nucl_coord

      do j = 1, ao_num
        do i = 1, ao_num

          deriv_nm = 0.5d0 * inv_eps * (Int_2(i,j) - Int_1(i,j))
          deriv_ex = fst_deriv_kin(i_axis, i_nucl, i, j)

          delta = dabs(deriv_nm - deriv_ex)
          accu_abs += delta

          if(delta .gt. 1.d-6) then
            print*, ' problem on: '
            print*, i_nucl, i_axis, i, j
            print*, deriv_nm, deriv_ex, delta
            !stop
          endif

        enddo
      enddo
    enddo
  enddo

  deallocate( Int_1, Int_2 )

  print*,'accu_abs = ', accu_abs

end subroutine test_fst_deriv_kin

! ---

subroutine test_fst_deriv_nai()

  implicit none
  integer                       :: i, j, i_axis, i_nucl
  double precision              :: deriv_ex, deriv_nm
  double precision              :: delta, accu_abs, eps, inv_eps
  double precision, allocatable :: Int_1(:,:), Int_2(:,:)

  double precision              :: fst_deriv_nai


  print*, ' check fst_deriv_nai:'

  eps     = 1d-8
  inv_eps = 1.d0 / eps

  allocate( Int_1(ao_num,ao_num), Int_2(ao_num,ao_num) )

  accu_abs = 0.d0
  accu_abs = 0.d0

  do i_nucl = 1, nucl_num
    do i_axis = 1, 3

      nucl_coord(i_nucl,i_axis) = nucl_coord(i_nucl,i_axis) + eps
      TOUCH nucl_coord
      Int_2(1:ao_num,1:ao_num) = ao_integrals_n_e(1:ao_num,1:ao_num)

      nucl_coord(i_nucl,i_axis) = nucl_coord(i_nucl,i_axis) - 2.d0 * eps
      TOUCH nucl_coord
      Int_1(1:ao_num,1:ao_num) = ao_integrals_n_e(1:ao_num,1:ao_num)

      nucl_coord(i_nucl,i_axis) = nucl_coord(i_nucl,i_axis) + eps
      TOUCH nucl_coord

      do j = 1, ao_num
        do i = 1, ao_num

          deriv_nm = 0.5d0 * inv_eps * (Int_2(i,j) - Int_1(i,j))
          deriv_ex = fst_deriv_nai(i_axis, i_nucl, i, j)

          delta = dabs(deriv_nm - deriv_ex)
          accu_abs += delta

          if(delta .gt. 1.d-6) then
            print*, ' problem on: '
            print*, i_nucl, i_axis, i, j
            print*, deriv_nm, deriv_ex, delta
            !stop
          endif

        enddo
      enddo
    enddo
  enddo

  deallocate( Int_1, Int_2 )

  print*,'accu_abs = ', accu_abs

end subroutine test_fst_deriv_nai

! ---

subroutine test_fst_deriv_eri()

  implicit none
  integer                       :: i, j, k, l, i_axis, i_nucl
  double precision              :: deriv_ex, deriv_nm
  double precision              :: delta, accu_abs, eps, inv_eps
  double precision, allocatable :: Int_1(:,:,:,:), Int_2(:,:,:,:)

  double precision              :: fst_deriv_eri, ao_two_e_integral


  print*, ' check fst_deriv_eri:'

  allocate( Int_1(ao_num,ao_num,ao_num,ao_num), Int_2(ao_num,ao_num,ao_num,ao_num) )

  eps     = 1d-8
  inv_eps = 1.d0 / eps

  accu_abs = 0.d0
  accu_abs = 0.d0

  do i_nucl = 1, nucl_num
    do i_axis = 1, 3

      nucl_coord(i_nucl,i_axis) = nucl_coord(i_nucl,i_axis) + eps
      TOUCH nucl_coord
      do i = 1, ao_num
        do j = 1, ao_num
          do k = 1, ao_num
            do l = 1, ao_num
              Int_2(i,j,k,l) = ao_two_e_integral(i, j, k, l)
            enddo
          enddo
        enddo
      enddo

      nucl_coord(i_nucl,i_axis) = nucl_coord(i_nucl,i_axis) - 2.d0 * eps
      TOUCH nucl_coord
      do i = 1, ao_num
        do j = 1, ao_num
          do k = 1, ao_num
            do l = 1, ao_num
              Int_1(i,j,k,l) = ao_two_e_integral(i, j, k, l)
            enddo
          enddo
        enddo
      enddo

      nucl_coord(i_nucl,i_axis) = nucl_coord(i_nucl,i_axis) + eps
      TOUCH nucl_coord

      do i = 1, ao_num
        do j = 1, ao_num
          do k = 1, ao_num
            do l = 1, ao_num
              deriv_nm = 0.5d0 * inv_eps * (Int_2(i,j,k,l) - Int_1(i,j,k,l))
              deriv_ex = fst_deriv_eri(i_axis, i_nucl, i, j, k, l)

              delta     = dabs(deriv_nm - deriv_ex)
              accu_abs += delta

              if(delta .gt. 1.d-6) then
                print*, ' problem on: '
                print*, i_nucl, i_axis, i, j, k, l
                print*, deriv_nm, deriv_ex, delta
                !stop
              endif

            enddo
          enddo
        enddo
      enddo

    enddo
  enddo

!  do i_nucl = 1, nucl_num
!    do i_axis = 1, 3
!
!      do i = 1, ao_num
!        do j = 1, ao_num
!          do k = 1, ao_num
!            do l = 1, ao_num
!
!              nucl_coord(i_nucl,i_axis) = nucl_coord(i_nucl,i_axis) + eps
!              TOUCH nucl_coord
!              Int_2 = ao_two_e_integral(i, j, k, l)
!              
!              nucl_coord(i_nucl,i_axis) = nucl_coord(i_nucl,i_axis) - 2.d0 * eps
!              TOUCH nucl_coord
!              Int_1 = ao_two_e_integral(i, j, k, l)
!              
!              nucl_coord(i_nucl,i_axis) = nucl_coord(i_nucl,i_axis) + eps
!              TOUCH nucl_coord
!
!              deriv_nm = 0.5d0 * inv_eps * (Int_2 - Int_1)
!              deriv_ex = fst_deriv_eri(i_axis, i_nucl, i, j, k, l)
!
!              delta     = dabs(deriv_nm - deriv_ex)
!              accu_abs += delta
!
!              if(delta .gt. 1.d-6) then
!                print*, ' problem on: '
!                print*, i_nucl, i_axis, i, j, k, l
!                print*, deriv_nm, deriv_ex, delta
!                !stop
!              endif
!
!            enddo
!          enddo
!        enddo
!      enddo
!    enddo
!  enddo

  deallocate( Int_1, Int_2 )
  print*,'accu_abs = ', accu_abs

end subroutine test_fst_deriv_eri

! ---

subroutine test_fst_deriv_eri_schwartz()

  implicit none
  integer                       :: i, j, k, l, i_axis, i_nucl
  double precision              :: deriv_ex, deriv_nm
  double precision              :: delta, accu_abs, eps, inv_eps
  double precision, allocatable :: Int_1(:,:,:,:), Int_2(:,:,:,:)

  double precision              :: fst_deriv_eri_schwartz, ao_two_e_integral


  print*, ' check fst_deriv_eri_schwartz:'

  allocate( Int_1(ao_num,ao_num,ao_num,ao_num), Int_2(ao_num,ao_num,ao_num,ao_num) )

  eps     = 1d-8
  inv_eps = 1.d0 / eps

  accu_abs = 0.d0
  accu_abs = 0.d0

  do i_nucl = 1, nucl_num
    do i_axis = 1, 3

      nucl_coord(i_nucl,i_axis) = nucl_coord(i_nucl,i_axis) + eps
      TOUCH nucl_coord
      do i = 1, ao_num
        do j = 1, ao_num
          do k = 1, ao_num
            do l = 1, ao_num
              Int_2(i,j,k,l) = ao_two_e_integral(i, j, k, l)
            enddo
          enddo
        enddo
      enddo

      nucl_coord(i_nucl,i_axis) = nucl_coord(i_nucl,i_axis) - 2.d0 * eps
      TOUCH nucl_coord
      do i = 1, ao_num
        do j = 1, ao_num
          do k = 1, ao_num
            do l = 1, ao_num
              Int_1(i,j,k,l) = ao_two_e_integral(i, j, k, l)
            enddo
          enddo
        enddo
      enddo

      nucl_coord(i_nucl,i_axis) = nucl_coord(i_nucl,i_axis) + eps
      TOUCH nucl_coord

      do i = 1, ao_num
        do j = 1, ao_num
          do k = 1, ao_num
            do l = 1, ao_num
              deriv_nm = 0.5d0 * inv_eps * (Int_2(i,j,k,l) - Int_1(i,j,k,l))
              deriv_ex = fst_deriv_eri_schwartz(i_axis, i_nucl, i, j, k, l)

              delta     = dabs(deriv_nm - deriv_ex)
              accu_abs += delta

              if(delta .gt. 1.d-6) then
                print*, ' problem on: '
                print*, i_nucl, i_axis, i, j, k, l
                print*, deriv_nm, deriv_ex, delta
                !stop
              endif

            enddo
          enddo
        enddo
      enddo

    enddo
  enddo

  deallocate( Int_1, Int_2 )
  print*,'accu_abs = ', accu_abs

end subroutine test_fst_deriv_eri_schwartz

! ---


